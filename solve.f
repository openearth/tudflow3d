!    Part of Dflow3d a 3D Navier Stokes solver with variable density for 
!    simulations of near field dredge plume mixing

!    Copyright (C) 2012  Lynyrd de Wit

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.


      subroutine adamsb2(ib,ie,jb,je,kb,ke)
c
c    Performes time integration with second order 
c    Adams-Bashforth scheme, i.e
c
c
c             n+1     n
c     dU     U     - U                               n
c    ---- = ------------ = 1.5*( -ADV + DIFF + Force)     -
c     dt        dt                                   n-1
c                          0.5*( -ADV + DIFF + Force)
c
c   This scheme is weakly instabel for pure advection,
c   and therefore a very small amount of physical diffusion
c   is necessary.
c   The timestep is limited with CFL=0.3 (see routine chkdt)
c
c 
      USE nlist
      USE sediment

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer  ib,ie,jb,je,kb,ke,n,t,me
      real     doldc(nfrac,0:i1,0:j1,0:k1),dnewc(nfrac,0:i1,0:j1,0:k1)
      real     dold(0:i1,0:j1,0:k1),dnew(0:i1,0:j1,0:k1)
      real     wsed(nfrac,0:i1,0:j1,0:k1),sumWkm(0:i1,0:j1,0:k1)
	real*8 pplus(imax,kmax)
!	real dUdt1(0:i1,0:j1,0:k1),dVdt1(0:i1,0:j1,0:k1),dWdt1(0:i1,0:j1,0:k1)
!	
!	dUdt1=dudt
!	dVdt1=dvdt
!	dWdt1=dwdt
	

!	real     Diffcof(0:i1,0:j1,0:k1)

!	Diffcof=ekm/Sc/Rnew

	me = momentum_exchange_obstacles !100 or 101 or 111 means no advection momentum exchange with bed-cells with HYB6

c********************************************************************
c     CALCULATE slipvelocity
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(cnew,Wnew,wsed,rnew,1,kmax,sumWkm,dt,dz)
	       !! Set boundary conditions jet in:
	       do t=1,tmax_inPpunt
	 	i=i_inPpunt(t)
	 	j=j_inPpunt(t)
		do n=1,nfrac
	    	  do k=kmax,kmax !kmax-kjet,kmax
		 	Wsed(n,i,j,k)=Wnew(i,j,k)
		  enddo
		enddo
		if (LOA>0.and.outflow_overflow_down.eq.1) then
		 do n=1,nfrac
	    	  do k=kmax-kjet+1,kmax 
		 	Wsed(n,i,j,k)=Wnew(i,j,k)
		  enddo
		 enddo	
		endif
	       enddo
	  else 
	    do n=1,nfrac
		wsed(n,:,:,:)=wnew
	    enddo
	  endif

c********************************************************************
c     CALCULATE advection, diffusion Concentration
c********************************************************************
	  do n=1,nfrac
	      call advecc_VLE(dnewc(n,:,:,:),cnew(n,:,:,:),Unew,Vnew,Wsed(n,:,:,:),rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)
	      call diffc_CDS2 (dnewc(n,:,:,:),Cnew(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)

		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(1.5*dnewc(n,:,:,:)-0.5*cc(n,:,:,:)) !Adams-Bashford
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n,dt) ! bc for intermediate dcdt AB2
		cc(n,:,:,:)=dnewc(n,:,:,:)
	  enddo

	  call state(dcdt,drdt) ! determine drdt with intermediate dcdt
	endif
c********************************************************************
c     CALCULATE advection, diffusion and Force U-velocity
c********************************************************************
	if (diffusion.eq.'COM4') then
		call stress_terms(Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	endif

      	if (convection.eq.'CDS2') then
      call advecu_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecu_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2)	 
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd,numdiff2)
	elseif(convection.eq.'uTVD') then
      call advecu_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif

	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(dnew,Srr,Spr,Szr,Spp,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif


      do k=1,kmax
         do j=1,jmax
            do i=1,imax
		
            dUdt(i,j,k)=Unew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i+1,j,k))+
     1     dt*(     1.5*dnew(i,j,k)-
     1              0.5*wx(i,j,k)
     1        )
     1     -(gx*cos_u(j)+gy*sin_u(j))*dt*(0.75*(rnew(i,j,k)+rnew(i+1,j,k))-0.25*(rold(i,j,k)+rold(i+1,j,k))-rho_b)  
     1     +Ppropx(i,j,k)*dt
              wx(i,j,k)=dnew(i,j,k)
            enddo
         enddo
      enddo

c********************************************************************
c     CALCULATE advection, diffusion and Force V-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecv_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecv_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecv_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif

	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(dnew,Spr,Spp,Spz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
            dVdt(i,j,k)=Vnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j+1,k))+
     1     dt*(      1.5*dnew(i,j,k)-
     1               0.5*wy(i,j,k) 
     1        )
     1     -(-gx*sin_v(j)+gy*cos_v(j))*dt*(0.75*(rnew(i,j,k)+rnew(i,j+1,k))-0.25*(rold(i,j,k)+rold(i,j+1,k))-rho_b)
!     1     +Ppropy(i,j,k)/(MAX(1.e-2,ABS(1.5*Vnew(i,j,k)-0.5*Vold(i,j,k))))*dt 
     1     +Ppropy(i,j,k)*dt

             wy(i,j,k)=dnew(i,j,k)
            enddo
         enddo
      enddo
c********************************************************************
c     CALCULATE advection, diffusion and Force W-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecw_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecw_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecw_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif

	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(dnew,Szr,Spz,Szz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 call advecw_driftfluxCDS2(dnew,driftfluxforce_calfac*sumWkm,
     &   rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
	endif
		

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
            dWdt(i,j,k)=Wnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j,k+1))+
     1          dt*(1.5*dnew(i,j,k)-0.5*wz(i,j,k)) -gz*dt*(0.75*(rnew(i,j,k)+rnew(i,j,k+1))-0.25*(rold(i,j,k)+rold(i,j,k+1))-rho_b)
     1     +Ppropz(i,j,k)*dt
		    wz(i,j,k)=dnew(i,j,k)
            enddo
         enddo
      enddo
	  pold2=pold !save pressure previous timestep [pold2 = timestep before pold]
	pold=p+pold    !what is called p here was dp in reality, now p is 
      call pshiftb(pold,pplus) !,rank,imax,jmax,kmax,px)
c********************************************************************
c     CALCULATE pressure-gradient with old pressure:
c********************************************************************
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -dt * ( pold(i+1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
	!!! do nothing on imax and i=1, so p(imax+1)=p(imax) --> dp/dn=0 at outflow and inflow 
	if (periodicx.eq.1) then
	!!! !periodic bc in x: 
	      i = imax
	      do k=1,kmax
		do j=1,jmax
		  dUdt(i,j,k) = dUdt(i,j,k) -dt * ( pold(1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		 enddo
	      enddo
	endif
	
      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * (phip(j+1)-phip(j)) )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
      	  endif
	endif

      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          dWdt(i,j,k) = dWdt(i,j,k) -dt * ( pold(i,j,k+1) - pold(i,j,k) ) / dz
          enddo
        enddo
      enddo

      return
      end

      subroutine adamsb3(ib,ie,jb,je,kb,ke)
c
c    Performes time integration with second order 
c    Adams-Bashforth-3 scheme, i.e
c
c
c             n+1     n
c     dU     U     - U                                 n
c    ---- = ------------ = 23/12*( -ADV + DIFF + Force)     
c     dt        dt                                     n-1
c                          -4/3 *( -ADV + DIFF + Force)
c                                                      n-2
c                          +5/12*( -ADV + DIFF + Force)
c
c   The timestep is limited with CFL=0.3-0.6 (see routine chkdt)
c
c 
      USE nlist
      USE sediment

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer  ib,ie,jb,je,kb,ke,n,t,me
      real     doldc(nfrac,0:i1,0:j1,0:k1),dnewc(nfrac,0:i1,0:j1,0:k1)
      real     dold(0:i1,0:j1,0:k1),dnew(0:i1,0:j1,0:k1)
      real     wsed(nfrac,0:i1,0:j1,0:k1),sumWkm(0:i1,0:j1,0:k1)
	real*8 pplus(imax,kmax)
	real dnewcbot(nfrac,0:i1,0:j1)
!	real dUdt1(0:i1,0:j1,0:k1),dVdt1(0:i1,0:j1,0:k1),dWdt1(0:i1,0:j1,0:k1)
!	
!	dUdt1=dudt
!	dVdt1=dvdt
!	dWdt1=dwdt
	
!	real     Diffcof(0:i1,0:j1,0:k1)

!	Diffcof=ekm/Sc/Rnew 


	me = momentum_exchange_obstacles !100 or 101 or 111 means no advection momentum exchange with bed-cells with HYB6
c********************************************************************
c     CALCULATE slipvelocity
c********************************************************************
c********************************************************************
c     CALCULATE advection, diffusion Concentration
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(cnew,Wnew,wsed,rnew,1,kmax,sumWkm,dt,dz)
	       !! Set boundary conditions jet in:
	       do t=1,tmax_inPpunt
	 	i=i_inPpunt(t)
	 	j=j_inPpunt(t)
		do n=1,nfrac
	    	  do k=kmax,kmax !-kjet,kmax
		 	Wsed(n,i,j,k)=Wnew(i,j,k)
		  enddo
		enddo
		if (LOA>0.and.outflow_overflow_down.eq.1) then
		 do n=1,nfrac
	    	  do k=kmax-kjet+1,kmax 
		 	Wsed(n,i,j,k)=Wnew(i,j,k)
		  enddo
		 enddo	
		endif	
	       enddo
	  else
	    do n=1,nfrac
		wsed(n,:,:,:)=wnew
	    enddo
	  endif
	  do n=1,nfrac
	      call advecc_VLE(dnewc(n,:,:,:),cnew(n,:,:,:),Unew,Vnew,Wsed(n,:,:,:),rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)
		
	      call diffc_CDS2 (dnewc(n,:,:,:),Cnew(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)
		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(1.5*dnewc(n,:,:,:)-0.5*cc(n,:,:,:)) !AB2
		cc(n,:,:,:)=dnewc(n,:,:,:)
	  enddo
	  call slipvelocity_bed(cnew,wnew,wsed,rnew,sumWkm,dt,dz) !driftflux_force must be calculated with settling velocity at bed
      	  if (interaction_bed>0) then
	    CALL erosion_deposition(dcdt,dcdtbot,unew,vnew,wnew,rnew,cnew,cnewbot,dt,dz) !first two vars are adjusted
		do n=1,nfrac 
		  call bound_cbot(dcdtbot(n,:,:))
		enddo			
		!! advec concentration in bed with velocity TSHD:
		IF (ABS(U_TSHD)>1e-6) THEN 
		  do n=1,nfrac 
			call adveccbot_TVD(dnewcbot(n,:,:),dcdtbot(n,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,phiv,phipt,dz,
     +            	i1,j1,ib,ie,jb,je,dt,rank,px,periodicx,periodicy)
			dcdtbot(n,:,:)= dcdtbot(n,:,:) + dt*dnewcbot(n,:,:) ! time update	
			call bound_cbot(dcdtbot(n,:,:))			
		  enddo 
		ENDIF 		
	    !! for kn1: erosion_deposition must be after advecc_TVD and after dcdt update, because cnew and dcdt are two separate vars
		if (interaction_bed.eq.4.or.interaction_bed.eq.6) then
			call advec_update_Clivebed(dcdt,dcdtbot,dt) 
		endif	
	    do n=1,nfrac 
	        call bound_cbot(dcdtbot(n,:,:))
	    enddo
      	  endif
	    do n=1,nfrac 
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n,dt) ! bc after erosion_deposition AB3
	    enddo
	call state(dcdt,drdt) ! determine drdt with intermediate dcdt
	endif


	if (diffusion.eq.'COM4') then
		call stress_terms(Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	endif
c********************************************************************
c     CALCULATE advection, diffusion and Force U-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecu_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff
     & ,periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecu_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff
     & ,periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2)	 
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd,numdiff2)	  
	elseif(convection.eq.'uTVD') then
      call advecu_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif

	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(dnew,Srr,Spr,Szr,Spp,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif


      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-(gx*cos_u(j)+gy*sin_u(j))*(0.5*(rnew(i,j,k)+rnew(i+1,j,k))-rho_b)
!     1     +Ppropx(i,j,k)/(MAX(1.e-2,(ABS(Unew(i,j,k)))))*dt 
     1     +Ppropx(i,j,k)

            dUdt(i,j,k)=Unew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i+1,j,k))+
     1     dt*(     23./12.*dnew(i,j,k)-
     1               4./3. *wx(i,j,k)+
     1               5./12.*wxold(i,j,k))
!     1     -(gx*cos_u(j)+gy*sin_u(j))*dt*(0.75*(rnew(i,j,k)+rnew(i+1,j,k))-0.25*(rold(i,j,k)+rold(i+1,j,k))-rho_b)  
              wxold(i,j,k) = wx(i,j,k)
              wx(i,j,k)    = dnew(i,j,k)
            enddo
         enddo
      enddo

c********************************************************************
c     CALCULATE advection, diffusion and Force V-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecv_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecv_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecv_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif

	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(dnew,Spr,Spp,Spz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-(-gx*sin_v(j)+gy*cos_v(j))*(0.5*(rnew(i,j,k)+rnew(i,j+1,k))-rho_b)
     1     +Ppropy(i,j,k)
            dVdt(i,j,k)=Vnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j+1,k))+
     1     dt*(      23./12.*dnew(i,j,k)-
     1                4./3. *wy(i,j,k)+
     1                5./12.*wyold(i,j,k))
!     1     -(-gx*sin_v(j)+gy*cos_v(j))*dt*(0.75*(rnew(i,j,k)+rnew(i,j+1,k))-0.25*(rold(i,j,k)+rold(i,j+1,k))-rho_b)
             wyold(i,j,k) = wy(i,j,k)
             wy(i,j,k)    = dnew(i,j,k)
            enddo
         enddo
      enddo
c********************************************************************
c     CALCULATE advection, diffusion and Force W-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecw_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecw_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecw_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif

	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(dnew,Szr,Spz,Szz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 call advecw_driftfluxCDS2(dnew,driftfluxforce_calfac*sumWkm,
     &   rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
	endif
	


      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-gz*(0.5*(rnew(i,j,k)+rnew(i,j,k+1))-rho_b)
     1     +Ppropz(i,j,k)
            dWdt(i,j,k)=Wnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j,k+1))+
     1     dt*(      23./12.*dnew(i,j,k)-
     1                4./3. *wz(i,j,k)+
     1                5./12.*wzold(i,j,k))
!     1          dt*(1.5*dnew(i,j,k)-0.5*wz(i,j,k)) -gz*dt*(0.75*(rnew(i,j,k)+rnew(i,j,k+1))-0.25*(rold(i,j,k)+rold(i,j,k+1))-rho_b)
		    wzold(i,j,k) = wz(i,j,k)
		    wz(i,j,k)    = dnew(i,j,k)
            enddo
         enddo
      enddo
	pold2=pold !save pressure previous timestep [pold2 = timestep before pold]
	pold=p+pold    !what is called p here was dp in reality, now p is 
      call pshiftb(pold,pplus) !,rank,imax,jmax,kmax,px)
c********************************************************************
c     CALCULATE pressure-gradient with old pressure:
c********************************************************************
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -dt * ( pold(i+1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
	  !!! do nothing on imax and i=1, so p(imax+1)=p(imax) --> dp/dn=0 at outflow and inflow 
	if (periodicx.eq.1) then
	!!! !periodic bc in x: 
	      i = imax
	      do k=1,kmax
		do j=1,jmax
		  dUdt(i,j,k) = dUdt(i,j,k) -dt * ( pold(1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		 enddo
	      enddo
	endif
      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * (phip(j+1)-phip(j)) )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
      	  endif
	endif
      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          dWdt(i,j,k) = dWdt(i,j,k) -dt * ( pold(i,j,k+1) - pold(i,j,k) ) / dz
          enddo
        enddo
      enddo

      return
      end

      subroutine adamsbv(ib,ie,jb,je,kb,ke)
c
c    Performes time integration with second order 
c    Adams-Bashforth-3 scheme with variable timestep, i.e
c
c
c             n+1     n
c     dU     U     - U                                 n
c    ---- = ------------ = 23/12*( -ADV + DIFF + Force)     
c     dt        dt                                     n-1
c                          -4/3 *( -ADV + DIFF + Force)
c                                                      n-2
c                          +5/12*( -ADV + DIFF + Force)
c
c   The timestep is limited with CFL=0.7 (see routine chkdt)
c
c 
      USE nlist
      USE sediment

      implicit none
	  
	  include 'mpif.h'
!       include 'param.txt'
!       include 'common.txt'
      integer  ib,ie,jb,je,kb,ke,n,t,kpp,kp,km,ip,im,jp,jm,me
      real     doldc(nfrac,0:i1,0:j1,0:k1),dnewc(nfrac,0:i1,0:j1,0:k1)
      real     dold(0:i1,0:j1,0:k1),dnew(0:i1,0:j1,0:k1),dnew2(0:i1,0:j1,0:k1)
      real     wsed(nfrac,0:i1,0:j1,0:k1),sumWkm(0:i1,0:j1,0:k1)
	real*8 pplus(imax,kmax)
	real fAB3_1,fAB3_2,fAB3_3,fAB3_4,facAB3_1,facAB3_2,facAB3_3,facAB3_4
	real dist,timeAB_desired(1:4)
	real dnewcbot(nfrac,0:i1,0:j1)	
	real aaa(0:k1),bbb(0:k1),ccc(0:k1),ekm_min,ekm_plus,rhss(0:k1),rho_min,rho_plus
	real aaax(0:i1),bbbx(0:i1),cccx(0:i1),rhssx(0:i1)
	real aaay(0:jmax*px+1),bbby(0:jmax*px+1),cccy(0:jmax*px+1),rhssy(0:jmax*px+1)
	real ans_T(1:i1,0:jmax*px+1,1:kmax/px+1),ekm_T(1:i1,0:jmax*px+1,1:kmax/px+1)!,ekm2(0:i1,0:j1,0:k1) 
	real uu_T(1:imax,0:jmax*px+1,1:kmax/px),vv_T(1:imax,0:jmax*px+1,1:kmax/px),ww_T(1:imax,0:jmax*px+1,1:kmax/px)
	real rhU_T(1:imax,0:jmax*px+1,1:kmax/px),rhV_T(1:imax,0:jmax*px+1,1:kmax/px),rhW_T(1:imax,0:jmax*px+1,1:kmax/px)
	real interface_T(1:i1,0:jmax*px+1)
	real utr(0:i1,0:j1,0:k1),vtr(0:i1,0:j1,0:k1),wtr(0:i1,0:j1,0:k1)	
	real ws(0:i1,0:j1,0:k1),gvector,CNz,dir_order_real
	integer ileng,ierr,itag,status(MPI_STATUS_SIZE),perx,pery,dir_order
	integer clock
	INTEGER, DIMENSION(:), ALLOCATABLE :: seed
	!real rhU2(0:i1,0:j1,0:k1),rhV2(0:i1,0:j1,0:k1),rhW2(0:i1,0:j1,0:k1)
!!!	real dUdt2(0:i1,0:j1,0:k1),dVdt2(0:i1,0:j1,0:k1),dWdt2(0:i1,0:j1,0:k1)
	real ddd1(0:i1,0:j1,0:k1),ddd2(0:i1,0:j1,0:k1),ddd3(0:i1,0:j1,0:k1)
!	real eee1(0:i1,0:j1,0:k1),eee2(0:i1,0:j1,0:k1),eee3(0:i1,0:j1,0:k1) !als deze twee regels uitgezet worden dan komen er andere uitkomsten bij CN31?
!	real fff1(0:i1,0:j1,0:k1),fff2(0:i1,0:j1,0:k1),fff3(0:i1,0:j1,0:k1)

	me = momentum_exchange_obstacles !100 or 101 or 111 means no advection momentum exchange with bed-cells with HYB6

c********************************************************************
c     CALCULATE coefficients for AB3 interpolation with variable dt
c********************************************************************

	do i=4,1,-1
		if (i.eq.1) then
		  timeAB_real(i)=time_n
		else
		  timeAB_real(i)=timeAB_real(i-1)
		endif
		timeAB_desired(i)=time_n-(i-1)*dt
	enddo
	
	!IF (CNdiffz.eq.11.or.CNdiffz.eq.12) THEN !switch to AB2 with variable time step:
	IF (.false.) THEN !always use AB3 with variable time step
!	!	!EE1:
!	!   facAB3_1=1. !1.5 !23./12.
!	!   facAB3_2=0. !-0.5 !-4./3.
!	!   facAB3_3=0. !5./12.
!	!   facAB3_4=0.
		if (istep.lt.5) then !start with standard AB2
		   facAB3_1=factorsABv3(1) !1.5 !23./12.
		   facAB3_2=factorsABv3(2) !-0.5 !-4./3.
		   facAB3_3=factorsABv3(3) !0. !5./12.
		   facAB3_4=0.
		else

			dist=timeAB_desired(2)-timeAB_real(2)
			if (dist.lt.0) then
				fAB3_2=(timeAB_real(2)-timeAB_real(3)+dist)/(timeAB_real(2)-timeAB_real(3))
				fAB3_3=1.-(timeAB_real(2)-timeAB_real(3)+dist)/(timeAB_real(2)-timeAB_real(3))
				fAB3_1=0.
			else
				fAB3_2=(timeAB_real(1)-timeAB_real(2)-dist)/(timeAB_real(1)-timeAB_real(2))
				fAB3_1=1.-(timeAB_real(1)-timeAB_real(2)-dist)/(timeAB_real(1)-timeAB_real(2))
				fAB3_3=0.
			endif

			facAB3_1=factorsABv3(1)+fAB3_1*factorsABv3(2) !23./12.-fAB3_1*4./3.
			facAB3_2=fAB3_2*factorsABv3(2)!-fAB3_2*4./3.
			facAB3_3=fAB3_3*factorsABv3(2) !-fAB3_3*4./3.
			facAB3_4=0.

		endif
	ELSE !AB3
		if (istep.lt.5) then !start with standard AB3
	!	   facAB3_1=1.
	!	   facAB3_2=0.
	!	   facAB3_3=0.
	!	   facAB3_4=0.
		   facAB3_1=factorsABv3(1) !23./12.
		   facAB3_2=factorsABv3(2) !-4./3.
		   facAB3_3=factorsABv3(3) !5./12.
		   facAB3_4=0.
		else

			dist=timeAB_desired(2)-timeAB_real(2)
			if (dist.lt.0) then
				fAB3_2=(timeAB_real(2)-timeAB_real(3)+dist)/(timeAB_real(2)-timeAB_real(3))
				fAB3_3=1.-(timeAB_real(2)-timeAB_real(3)+dist)/(timeAB_real(2)-timeAB_real(3))
				fAB3_1=0.
			else
				fAB3_2=(timeAB_real(1)-timeAB_real(2)-dist)/(timeAB_real(1)-timeAB_real(2))
				fAB3_1=1.-(timeAB_real(1)-timeAB_real(2)-dist)/(timeAB_real(1)-timeAB_real(2))
				fAB3_3=0.
			endif

			facAB3_1=factorsABv3(1)+factorsABv3(2)*fAB3_1 !23./12.-fAB3_1*4./3.
			facAB3_2=factorsABv3(2)*fAB3_2 !-fAB3_2*4./3.
			facAB3_3=factorsABv3(2)*fAB3_3 !-fAB3_3*4./3.

			dist=timeAB_desired(3)-timeAB_real(3)
			if (dist.lt.0) then
				fAB3_3=(timeAB_real(3)-timeAB_real(4)+dist)/(timeAB_real(3)-timeAB_real(4))
				fAB3_4=1.-(timeAB_real(3)-timeAB_real(4)+dist)/(timeAB_real(3)-timeAB_real(4))
				fAB3_2=0.
			else
				fAB3_3=(timeAB_real(2)-timeAB_real(3)-dist)/(timeAB_real(2)-timeAB_real(3))
				fAB3_2=1.-(timeAB_real(2)-timeAB_real(3)-dist)/(timeAB_real(2)-timeAB_real(3))
				fAB3_4=0.
			endif

			facAB3_2=facAB3_2+fAB3_2*factorsABv3(3) !facAB3_2+fAB3_2*5./12.
			facAB3_3=facAB3_3+fAB3_3*factorsABv3(3) !facAB3_3+fAB3_3*5./12.
			facAB3_4=fAB3_4*factorsABv3(3) !fAB3_4*5./12.

		endif	
	ENDIF 

	   if (applyVOF.eq.1) then ! VOF fraction apply real density for gravity and sediment settling related things
		 CALL state(cnew,rnew)
	   endif 
c********************************************************************
c     CALCULATE tau walls based on divergence free velocity field previous time step
c********************************************************************	
	dUdt = Unew*rhU ! for some continuity_solver dUdt-end is without rho 
	dVdt = Vnew*rhV	
	dWdt = Wnew*rhW
	call bound_rhoU(dUdt,dVdt,dWdt,rnew,MIN(3,slip_bot),monopile,time_n,
     & 		Ub1old,Vb1old,Wb1old,Ub2old,Vb2old,Wb2old,Ub3old,Vb3old,Wb3old) ! apply wall shear stress based on UVWnew (valid at time_n and UB1old...)

	dcdt = 0.
	sumWkm = 0.
	dnewc=0.
	dnew2=0. !initialize as zero
	
c********************************************************************
c     CALCULATE slipvelocity
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(cnew,Wnew,wsed,rnew,1,kmax,sumWkm,dt,dz)
	       !! Set boundary conditions jet in:
	       do t=1,tmax_inPpunt
	 	i=i_inPpunt(t)
	 	j=j_inPpunt(t)
		do n=1,nfrac
	    	  do k=kmax,kmax !-kjet,kmax
		 	Wsed(n,i,j,k)=Wnew(i,j,k)
		  enddo
		enddo
		if (LOA>0.and.outflow_overflow_down.eq.1) then
		 do n=1,nfrac
	    	  do k=kmax-kjet+1,kmax 
			Wsed(n,i,j,k)=Wnew(i,j,k)
		  enddo
		 enddo	
		endif
	       enddo
	  else
	    do n=1,nfrac
		wsed(n,:,:,:)=wnew
	    enddo
	  endif

c********************************************************************
c     CALCULATE advection, diffusion Concentration
c********************************************************************

	  do n=1,nfrac
		if (transporteq_fracs.eq.'massfrac') then
			utr=rhoU
			vtr=rhoV		
			wtr(:,:,0:kmax)=wsed(n,:,:,0:kmax)*0.5*(rnew(:,:,0:kmax)+rnew(:,:,1:k1))
			cnew(n,:,:,:)=cnew(n,:,:,:)*frac(n)%rho/rnew(:,:,:)  ! go from volume fraction to mass fraction
			Diffcof=Diffcof*rnew
		elseif (transporteq_fracs.eq.'massfra2') then
			cnew(n,:,:,:)=cnew(n,:,:,:)*frac(n)%rho  ! go from volume fraction to (mass-fraction*rho_mix)
			utr=Unew
			vtr=Vnew	
			wtr=Wsed(n,:,:,:)		
		else
			utr=Unew
			vtr=Vnew	
			wtr=Wsed(n,:,:,:)
		endif	
		if (settling_along_gvector.eq.1) then
			gvector=sqrt(gx**2+gy**2+gz**2)
			ws = Wsed(n,:,:,:)-Wnew 
			utr=utr+gx/gvector*ws 
			vtr=vtr+gy/gvector*ws
			wtr=wtr-ws+gz/gvector*ws
		endif 
	    call make_UtransC_zeroIBM(utr,vtr,wtr)
		if (frac(n)%type.eq.-1) then ! VOF fraction
	         call advecc_VOF(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
!	       call advecc_NVD_VOF(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
!     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)	 
			 dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(dnewc(n,:,:,:)) !update in time with EE1 for TVD, no diffusion
			 do i=1,imax
			   do j=1,jmax
			     do k=1,kmax
					dcdt(n,i,j,k)=MAX(0.,dcdt(n,i,j,k))
					dcdt(n,i,j,k)=MIN(1.,dcdt(n,i,j,k))
				 enddo
			   enddo
			 enddo 
		else
		  if (advec_conc.eq.'NVD') then
	       call advecc_NVD(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)	 
          elseif (advec_conc.eq.'VLE') then	 
	       call advecc_VLE(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)
          elseif (advec_conc.eq.'SBE') then	 
	       call advecc_SBE(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)	 
          elseif (advec_conc.eq.'VL2') then	 
	       call advecc_VL2_nocfl(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
          elseif (advec_conc.eq.'SB2') then	 
	       call advecc_SB2_nocfl(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)	 
		  elseif (advec_conc.eq.'ARO') then	 
	       call advecc_ARO(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)	 
          endif
		
		  if (Sc<1.e18) then
	        call diffc_CDS2 (dnewc(n,:,:,:),Cnew(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)
		  endif 

		  dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(dnewc(n,:,:,:)) !update in time with EE1 for TVD
		endif
		if (transporteq_fracs.eq.'massfrac') then
			dcdt(n,:,:,:)=dcdt(n,:,:,:)+cnew(n,:,:,:)*rnew-cnew(n,:,:,:) 
			dcdt(n,:,:,:)=dcdt(n,:,:,:)/frac(n)%rho ! from dcdt (mass-fraction*rho_mix) back to volume-fraction
			cnew(n,:,:,:)=cnew(n,:,:,:)*rnew/frac(n)%rho ! from mass fraction back to volume fraction
			Diffcof=Diffcof/rnew
		elseif (transporteq_fracs.eq.'massfra2') then
			dcdt(n,:,:,:)=dcdt(n,:,:,:)/frac(n)%rho ! from dcdt (mass-fraction*rho_mix) back to volume-fraction
			cnew(n,:,:,:)=cnew(n,:,:,:)/frac(n)%rho ! from (mass-fraction*rho_mix) back to volume fraction
		endif
		
			IF ((CNdiffz.eq.1.or.CNdiffz.eq.2.or.CNdiffz.eq.11.or.CNdiffz.eq.12).and.Sc<1.e18) THEN !CN semi implicit treatment diff-z:
			 IF (CNdiffz.eq.1.or.CNdiffz.eq.11) THEN 
				CNz=0.5
			 ELSE 
				CNz=1.
			 ENDIF 
			 call bound_c(dcdt(n,:,:,:),frac(n)%c,n,0.) ! bc start CN-diffz ABv
			 do k=k1,k1 !-kjet,k1
			   do t=1,tmax_inPpunt
				 i=i_inPpunt(t)
				 j=j_inPpunt(t)
				 dcdt(n,i,j,k) = dcdt(n,i,j,kmax) ! No diffusion over vertical inflow boundary, this line is needed for exact influx
			   enddo
			 enddo			 
!!!			 dUdt2=dcdt(n,:,:,:) !for non-updated RHS 
			 IF (CNdiffz.eq.11.or.CNdiffz.eq.12) THEN 
				 perx=0
				 pery=0 
				 IF (periodicx.eq.1) perx=1
				 IF (periodicy.eq.1) pery=1
				 call t2np(Diffcof(1:imax,1:jmax,1:kmax),ekm_T(1:imax,1:jmax*px,1:kmax/px))
				 call t2np(dcdt(n,1:imax,1:jmax,1:kmax),ans_T(1:imax,1:jmax*px,1:kmax/px))
				 !! also pass over boundaries at j=0 :
				 IF (rank.eq.0) THEN
				  do i=1,px-1
				  call mpi_send(Diffcof(1:imax,0,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i,MPI_COMM_WORLD,status,ierr)
				  call mpi_send(dcdt(n,1:imax,0,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+100000,MPI_COMM_WORLD,status,ierr)
				  enddo
				  ekm_T(1:imax,0,1:kmax/px)=Diffcof(1:imax,0,1:kmax/px)
				  ans_T(1:imax,0,1:kmax/px)=dcdt(n,1:imax,0,1:kmax/px)			
				 ELSE
					call mpi_recv(ekm_T(1:imax,0,1:kmax/px),imax*kmax/px,MPI_REAL8,0,rank,MPI_COMM_WORLD,status,ierr)
					call mpi_recv(ans_T(1:imax,0,1:kmax/px),imax*kmax/px,MPI_REAL8,0,rank+100000,MPI_COMM_WORLD,status,ierr)
				 ENDIF
				 !! also pass over boundaries at j=jmax+1 :
				 IF (rank.eq.px-1) THEN
				  do i=0,px-2
					call mpi_send(Diffcof(1:imax,jmax+1,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+200000,MPI_COMM_WORLD
     &					,status,ierr)
					call mpi_send(dcdt(n,1:imax,jmax+1,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+300000,MPI_COMM_WORLD
     &					,status,ierr)	 
				  enddo
				  ekm_T(1:imax,jmax*px+1,1:kmax/px)=Diffcof(1:imax,jmax+1,rank*kmax/px+1:(rank+1)*kmax/px)
				  ans_T(1:imax,jmax*px+1,1:kmax/px)=dcdt(n,1:imax,jmax+1,rank*kmax/px+1:(rank+1)*kmax/px)
				 ELSE
				  call mpi_recv(ekm_T(1:imax,jmax*px+1,1:kmax/px),imax*kmax/px,MPI_REAL8,px-1,rank+200000,MPI_COMM_WORLD,status,ierr)
				  call mpi_recv(ans_T(1:imax,jmax*px+1,1:kmax/px),imax*kmax/px,MPI_REAL8,px-1,rank+300000,MPI_COMM_WORLD,status,ierr)
				 ENDIF

				IF (periodicy.eq.0) THEN !bc defined at both lateral boundaries 
				  do k=1,kmax/px 
				   do i=1,imax
					 do j=1,px*jmax !0,px*jmax+1 
					   jm=j-1 !MAX(0,j-1)
					   jp=j+1 !MIN(px*jmax+1,j+1)
					   ekm_min=0.5* (ekm_T(i,jm,k)+ekm_T(i,j,k))*MIN(fc_global(i,jm,k+rank*kmax/px),fc_global(i,j,k+rank*kmax/px))
					   ekm_plus=0.5*(ekm_T(i,jp,k)+ekm_T(i,j,k))*MIN(fc_global(i,jp,k+rank*kmax/px),fc_global(i,j,k+rank*kmax/px))
					   aaay(j)=-CNz*ekm_min*dt/((Rp(i)*(phipt(j)-phipt(jm)))*Rp(i)*dphi2t(j))
					   cccy(j)=-CNz*ekm_plus*dt/((Rp(i)*(phipt(jp)-phipt(j)))*Rp(i)*dphi2t(j))
					   bbby(j)=1.-aaay(j)-cccy(j) 				   
					 enddo
					 aaay(0)=0.
					 aaay(px*jmax+1)=0.
					 cccy(0)=0.
					 cccy(px*jmax+1)=0.
					 bbby(0)=1.
					 bbby(px*jmax+1)=1.
					 rhssy=ans_T(i,0:px*jmax+1,k)
					 CALL solve_tridiag_switchperiodic2(ans_T(i,0:px*jmax+1,k),aaay,bbby,cccy,rhssy,px*jmax+2,pery)					
				   enddo
				  enddo
				ELSEIF (periodicy.eq.1) THEN !periodic lateral boundaries 
				  do k=1,kmax/px 
				   do i=1,imax
					 do j=1,px*jmax !0,px*jmax+1 
					   jm=j-1 !MAX(0,j-1)
					   jp=j+1 !MIN(px*jmax+1,j+1)
					   ekm_min=0.5* (ekm_T(i,jm,k)+ekm_T(i,j,k))*MIN(fc_global(i,jm,k+rank*kmax/px),fc_global(i,j,k+rank*kmax/px))
					   ekm_plus=0.5*(ekm_T(i,jp,k)+ekm_T(i,j,k))*MIN(fc_global(i,jp,k+rank*kmax/px),fc_global(i,j,k+rank*kmax/px))
					   aaay(j)=-CNz*ekm_min*dt/((Rp(i)*(phipt(j)-phipt(jm)))*Rp(i)*dphi2t(j))
					   cccy(j)=-CNz*ekm_plus*dt/((Rp(i)*(phipt(jp)-phipt(j)))*Rp(i)*dphi2t(j))
					   bbby(j)=1.-aaay(j)-cccy(j) 				   
					 enddo
					 rhssy(1:px*jmax)=ans_T(i,1:px*jmax,k)
					 CALL solve_tridiag_switchperiodic2(ans_T(i,1:px*jmax,k),aaay(1:px*jmax),bbby(1:px*jmax),cccy(1:px*jmax),
     &			rhssy(1:px*jmax),px*jmax,pery)
				   enddo
				  enddo
				ELSEIF (periodicy.eq.2) THEN !free slip lateral boundaries dUWdn=0 at both ends and V=0
				  do k=1,kmax/px 
				   do i=1,imax
					 do j=1,px*jmax !0,px*jmax+1 
					   jm=j-1 !MAX(0,j-1)
					   jp=j+1 !MIN(px*jmax+1,j+1)
					   ekm_min=0.5* (ekm_T(i,jm,k)+ekm_T(i,j,k))*MIN(fc_global(i,jm,k+rank*kmax/px),fc_global(i,j,k+rank*kmax/px))
					   ekm_plus=0.5*(ekm_T(i,jp,k)+ekm_T(i,j,k))*MIN(fc_global(i,jp,k+rank*kmax/px),fc_global(i,j,k+rank*kmax/px))
					   aaay(j)=-CNz*ekm_min*dt/((Rp(i)*(phipt(j)-phipt(jm)))*Rp(i)*dphi2t(j))
					   cccy(j)=-CNz*ekm_plus*dt/((Rp(i)*(phipt(jp)-phipt(j)))*Rp(i)*dphi2t(j))
					   bbby(j)=1.-aaay(j)-cccy(j) 				   
					 enddo
					bbby(px*jmax)=bbby(px*jmax)+cccy(px*jmax) !d.dn=0
					cccy(px*jmax)=0.
					bbby(1)=bbby(1)+aaay(1) !d.dn=0
					aaay(1)=0.
					rhssy(1:px*jmax)=ans_T(i,1:px*jmax,k)
					CALL solve_tridiag_switchperiodic2(ans_T(i,1:px*jmax,k),aaay(1:px*jmax),bbby(1:px*jmax),
     &			cccy(1:px*jmax),rhssy(1:px*jmax),px*jmax,pery)	
				   enddo
				  enddo
				ENDIF 
				call t2fp(ans_T(1:imax,1:jmax*px,1:kmax/px),dcdt(n,1:imax,1:jmax,1:kmax))	
				 
				IF (periodicx.eq.0) THEN !Dirichlet inflow and Neumann outflow 
				 do k=1,kmax 
				   do j=1,jmax 
					 do i=1,imax !0,i1
					   im=i-1 !MAX(0,i-1)
					   ip=i+1 !MIN(i1,i+1)
					   ekm_min=0.5*(Diffcof(im,j,k)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,k),fc_global(im,j+rank*jmax,k))
					   ekm_plus=0.5*(Diffcof(ip,j,k)+Diffcof(i,j,k))*MIN(fc_global(ip,j+rank*jmax,k),fc_global(i,j+rank*jmax,k))
					   aaax(i)=-CNz*ekm_min*Ru(im)*dt/((Rp(i)-Rp(im))*dr(i)*Rp(i))
					   cccx(i)=-CNz*ekm_plus*Ru(i)*dt/((Rp(ip)-Rp(i))*dr(i)*Rp(i))
					   bbbx(i)=1.-aaax(i)-cccx(i) 
					 enddo
					 aaax(0)=0.
					 cccx(0)=0.
					 bbbx(0)=1.
					 bbbx(imax)=bbbx(imax)+cccx(imax) !Neumann outflow 
					 cccx(imax)=0.
					 rhssx(0:imax)=dCdt(n,0:imax,j,k) 
					 CALL solve_tridiag_switchperiodic2(dcdt(n,0:imax,j,k),aaax(0:imax),bbbx(0:imax),cccx(0:imax),rhssx(0:imax),i1,perx) 
				   enddo
				 enddo	
				ELSEIF (periodicx.eq.2) THEN !Neumann inflow and outflow 
				 do k=1,kmax 
				   do j=1,jmax 
					 do i=1,imax !0,i1
					   im=i-1 !MAX(0,i-1)
					   ip=i+1 !MIN(i1,i+1)
					   ekm_min=0.5*(Diffcof(im,j,k)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,k),fc_global(im,j+rank*jmax,k))
					   ekm_plus=0.5*(Diffcof(ip,j,k)+Diffcof(i,j,k))*MIN(fc_global(ip,j+rank*jmax,k),fc_global(i,j+rank*jmax,k))
					   aaax(i)=-CNz*ekm_min*Ru(im)*dt/((Rp(i)-Rp(im))*dr(i)*Rp(i))
					   cccx(i)=-CNz*ekm_plus*Ru(i)*dt/((Rp(ip)-Rp(i))*dr(i)*Rp(i))
					   bbbx(i)=1.-aaax(i)-cccx(i) 
					 enddo
					 bbbx(1)=bbbx(1)+aaax(1)
					 aaax(1)=0.
					 bbbx(imax)=bbbx(imax)+cccx(imax)
					 cccx(imax)=0.
					 rhssx(1:imax)=dCdt(n,1:imax,j,k) 
					 CALL solve_tridiag_switchperiodic2(dcdt(n,1:imax,j,k),aaax(1:imax),bbbx(1:imax),cccx(1:imax),rhssx(1:imax),imax,perx) 
				   enddo
				 enddo					
				ELSEIF (periodicx.eq.1) THEN !periodic bc
				 do k=1,kmax 
				   do j=1,jmax 
					 do i=1,imax !0,i1
					   im=i-1 !MAX(0,i-1)
					   ip=i+1 !MIN(i1,i+1)
					   ekm_min=0.5*(Diffcof(im,j,k)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,k),fc_global(im,j+rank*jmax,k))
					   ekm_plus=0.5*(Diffcof(ip,j,k)+Diffcof(i,j,k))*MIN(fc_global(ip,j+rank*jmax,k),fc_global(i,j+rank*jmax,k))
					   aaax(i)=-CNz*ekm_min*Ru(im)*dt/((Rp(i)-Rp(im))*dr(i)*Rp(i))
					   cccx(i)=-CNz*ekm_plus*Ru(i)*dt/((Rp(ip)-Rp(i))*dr(i)*Rp(i))
					   bbbx(i)=1.-aaax(i)-cccx(i) 
					 enddo
					 rhssx(1:imax)=dCdt(n,1:imax,j,k) 
					 CALL solve_tridiag_switchperiodic2(dcdt(n,1:imax,j,k),aaax(1:imax),bbbx(1:imax),cccx(1:imax),rhssx(1:imax),imax,perx) 
				   enddo
				 enddo
				ENDIF 
			 ENDIF

             do j=1,jmax
               do i=1,imax
                 do k=1,kmax !0,k1
				   km=k-1 !MAX(0,k-1)
				   kp=k+1 !MIN(k1,k+1)
				   ekm_min=0.5*(Diffcof(i,j,km)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km))
				   ekm_plus=0.5*(Diffcof(i,j,kp)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,kp),fc_global(i,j+rank*jmax,k))
				   aaa(k)=-CNz*ekm_min*dt/dz**2
				   ccc(k)=-CNz*ekm_plus*dt/dz**2
				   bbb(k)=1.-aaa(k)-ccc(k) 
                 enddo
				 bbb(1)=bbb(1)+aaa(1)
				 aaa(1)=0.
				 bbb(kmax)=bbb(kmax)+ccc(kmax)
				 ccc(kmax)=0.
				 rhss(1:kmax)=dCdt(n,i,j,1:kmax) 
				 CALL solve_tridiag(dcdt(n,i,j,1:kmax),aaa(1:kmax),bbb(1:kmax),ccc(1:kmax),rhss(1:kmax),kmax) 
               enddo
             enddo
			ENDIF
          enddo !end nfrac
		  call slipvelocity_bed(cnew,wnew,wsed,rnew,sumWkm,dt,dz) !driftflux_force must be calculated with settling velocity at bed
      	  if (interaction_bed>0) then
	        CALL erosion_deposition(dcdt,dcdtbot,unew,vnew,wnew,rnew,cnew,cnewbot,dt,dz) !first two vars are adjusted
	        do n=1,nfrac 
	          call bound_cbot(dcdtbot(n,:,:))
	        enddo			
		    !! advec concentration in bed with velocity TSHD:
			IF (ABS(U_TSHD)>1e-6) THEN 
			  do n=1,nfrac 
				call adveccbot_TVD(dnewcbot(n,:,:),dcdtbot(n,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,phiv,phipt,dz,
     +            	i1,j1,ib,ie,jb,je,dt,rank,px,periodicx,periodicy)
				dcdtbot(n,:,:)= dcdtbot(n,:,:) + dt*dnewcbot(n,:,:) ! time update	
				call bound_cbot(dcdtbot(n,:,:))			
			  enddo 
			ENDIF 
	        !! erosion_deposition must be after advecc_TVD and after dcdt update, because cnew and dcdt are two separate vars
		    if (interaction_bed.eq.4.or.interaction_bed.eq.6) then
			  call advec_update_Clivebed(dcdt,dcdtbot,dt) 
		    endif		
	        do n=1,nfrac 
	          call bound_cbot(dcdtbot(n,:,:))
	        enddo
      	  endif
	    do n=1,nfrac 
		  call bound_c(dcdt(n,:,:,:),frac(n)%c,n,dt) ! bc after erosion_deposition ABv
	    enddo
	    call state(dcdt,drdt) ! determine drdt with intermediate dcdt

	   if (applyVOF.eq.1.or.applyVOF.eq.2) then ! VOF fraction, do not solve fully variable density NavierStokes but NS with constant density in advection and diffusion; only apply variable density in gravity term
		 drdt=rho_b
		 rhU=rho_b
		 rhV=rho_b
		 rhW=rho_b
	   endif
	endif

	if (diffusion.eq.'COM4') then
		call stress_terms(Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	endif
c********************************************************************
c     CALCULATE advection, diffusion and Force U-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecu_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff
     & ,periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecu_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff
     & ,periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2)	 
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd,numdiff2)	  
	elseif(convection.eq.'uTVD') then
      call advecu_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)		
	elseif(convection.eq.'C2Bl') then
      call advecu_C2Blend(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,
     & rank,px,dt,periodicx,periodicy,numdiff)
	elseif(convection.eq.'C4Bl') then
      call advecu_C4Blend(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,
     & rank,px,dt,periodicx,periodicy,numdiff)	 
	endif

	if (CNdiffz.eq.1.or.CNdiffz.eq.2.or.CNdiffz.eq.11.or.CNdiffz.eq.12.or.CNdiffz.eq.31) then 
		!!Hirsch 2007 2nd edition Eq. 9.4.21 p.478 mentions that also explicit part of CN should build on for 3 subsequent directions, hence here the explicit part is skipped
		if (CNdiffz.eq.11.or.CNdiffz.eq.12.or.CNdiffz.eq.31) then
		!if (CNdiffz.eq.12) then
			dnew2=0. 
		else 
			if (diffusion.eq.'CDS2') then
			  !call diffu_CDS2_CNexpl (dnew2,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			  dnew2=0.
			  call diffu_CDS2 (dnew2,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			elseif(diffusion.eq.'COM4') then
			  call diffu_com4(dnew2,Srr,Spr,Szr,Spp,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
			endif
		endif 
	else
		if (diffusion.eq.'CDS2') then
		  call diffu_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
		elseif(diffusion.eq.'COM4') then
		  call diffu_com4(dnew,Srr,Spr,Szr,Spp,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
		endif
	endif 
	
	if (applyVOF.eq.1) then ! VOF fraction apply real density only here for correct determination gravity
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-(gx*cos_u(j)+gy*sin_u(j))*(0.5*(rnew(i,j,k)+rnew(i+1,j,k))-rho_b)/
     1 (0.5*(rnew(i,j,k)+rnew(i+1,j,k))/rho_b)
     1     +Ppropx(i,j,k) 
            dUdt(i,j,k)=dUdt(i,j,k)+ !0.5*(rnew(i,j,k)+rnew(i+1,j,k))+
     1     dt*(     facAB3_1*dnew(i,j,k)+
     1              facAB3_2*wx(i,j,k)+
     1              facAB3_3*wxold(i,j,k)+
     1              facAB3_4*wxolder(i,j,k)+  
     1              dnew2(i,j,k)	 ) !AB3v for advection and body forces, CN for diffusion when CNdiffz=1 or 2
              wxolder(i,j,k) = wxold(i,j,k)
              wxold(i,j,k)   = wx(i,j,k)
              wx(i,j,k)      = dnew(i,j,k)
            enddo
         enddo
      enddo	
	else 
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-(gx*cos_u(j)+gy*sin_u(j))*(0.5*(rnew(i,j,k)+rnew(i+1,j,k))-rho_b)
     1     +Ppropx(i,j,k) 
            dUdt(i,j,k)=dUdt(i,j,k)+ !0.5*(rnew(i,j,k)+rnew(i+1,j,k))+
     1     dt*(     facAB3_1*dnew(i,j,k)+
     1              facAB3_2*wx(i,j,k)+
     1              facAB3_3*wxold(i,j,k)+
     1              facAB3_4*wxolder(i,j,k)+
     1              dnew2(i,j,k)	 ) !AB3v for advection and body forces, CN for diffusion when CNdiffz=1 or 2	 
              wxolder(i,j,k) = wxold(i,j,k)
              wxold(i,j,k)   = wx(i,j,k)
              wx(i,j,k)      = dnew(i,j,k)
            enddo
         enddo
      enddo
	endif 

c********************************************************************
c     CALCULATE advection, diffusion and Force V-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecv_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff
     & ,periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecv_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff
     & ,periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecv_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)
	elseif(convection.eq.'C2Bl') then
      call advecv_C2Blend(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,
     & rank,px,dt,periodicx,periodicy,numdiff)	
	elseif(convection.eq.'C4Bl') then
      call advecv_C4Blend(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,
     & rank,px,dt,periodicx,periodicy,numdiff)	 
	endif

	if (CNdiffz.eq.1.or.CNdiffz.eq.2.or.CNdiffz.eq.11.or.CNdiffz.eq.12.or.CNdiffz.eq.31) then 
		!!Hirsch 2007 2nd edition Eq. 9.4.21 p.478 mentions that also explicit part of CN should build on for 3 subsequent directions, hence here the explicit part is skipped
		if (CNdiffz.eq.11.or.CNdiffz.eq.12.or.CNdiffz.eq.31) then
		!if (CNdiffz.eq.12) then
			dnew2=0. 
		else 
			if (diffusion.eq.'CDS2') then
			  !call diffv_CDS2_CNexpl (dnew2,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			  dnew2=0.
			  call diffv_CDS2 (dnew2,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			elseif(diffusion.eq.'COM4') then
			  call diffv_com4(dnew2,Spr,Spp,Spz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
			endif
		endif 
	else 
		if (diffusion.eq.'CDS2') then
		  call diffv_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
		elseif(diffusion.eq.'COM4') then
		  call diffv_com4(dnew,Spr,Spp,Spz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
		endif
	endif 
	
	if (applyVOF.eq.1) then ! VOF fraction apply real density only here for correct determination gravity
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-(-gx*sin_v(j)+gy*cos_v(j))*(0.5*(rnew(i,j,k)+rnew(i,j+1,k))-rho_b)/
     1 (0.5*(rnew(i,j,k)+rnew(i,j+1,k))/rho_b)
     1     +Ppropy(i,j,k) 
            dVdt(i,j,k)=dVdt(i,j,k)+ !*0.5*(rnew(i,j,k)+rnew(i,j+1,k))+
     1     dt*(      facAB3_1*dnew(i,j,k)+
     1               facAB3_2*wy(i,j,k)+
     1               facAB3_3*wyold(i,j,k)+
     1               facAB3_4*wyolder(i,j,k)+
     1               dnew2(i,j,k)	 ) !AB3v for advection and body forces, CN for diffusion when CNdiffz=1 or 2	 
             wyolder(i,j,k) = wyold(i,j,k)
             wyold(i,j,k)   = wy(i,j,k)
             wy(i,j,k)      = dnew(i,j,k)
            enddo
         enddo
      enddo	
	else 
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-(-gx*sin_v(j)+gy*cos_v(j))*(0.5*(rnew(i,j,k)+rnew(i,j+1,k))-rho_b)
     1     +Ppropy(i,j,k) 
            dVdt(i,j,k)=dVdt(i,j,k)+ !*0.5*(rnew(i,j,k)+rnew(i,j+1,k))+
     1     dt*(      facAB3_1*dnew(i,j,k)+
     1               facAB3_2*wy(i,j,k)+
     1               facAB3_3*wyold(i,j,k)+
     1               facAB3_4*wyolder(i,j,k)+
     1               dnew2(i,j,k)	 ) !AB3v for advection and body forces, CN for diffusion when CNdiffz=1 or 2	 
             wyolder(i,j,k) = wyold(i,j,k)
             wyold(i,j,k)   = wy(i,j,k)
             wy(i,j,k)      = dnew(i,j,k)
            enddo
         enddo
      enddo
	endif 
c********************************************************************
c     CALCULATE advection, diffusion and Force W-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecw_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecw_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecw_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	
	elseif(convection.eq.'C2Bl') then
      call advecw_C2Blend(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,
     & rank,px,dt,periodicx,periodicy,numdiff)
	elseif(convection.eq.'C4Bl') then
      call advecw_C4Blend(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,
     & rank,px,dt,periodicx,periodicy,numdiff)
	endif

	if (CNdiffz.eq.1.or.CNdiffz.eq.2.or.CNdiffz.eq.11.or.CNdiffz.eq.12.or.CNdiffz.eq.31) then 
		!Hirsch 2007 2nd edition Eq. 9.4.21 p.478 mentions that also explicit part of CN should build on for 3 subsequent directions, hence here the explicit part is skipped
		if (CNdiffz.eq.11.or.CNdiffz.eq.12.or.CNdiffz.eq.31) then 
		!if (CNdiffz.eq.12) then
			dnew2=0. 
		else 
			if (diffusion.eq.'CDS2') then
			  !call diffw_CDS2_CNexpl (dnew2,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			  dnew2=0.
			  call diffw_CDS2 (dnew2,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)			  
			elseif(diffusion.eq.'COM4') then
			  call diffw_com4(dnew2,Szr,Spz,Szz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
			endif
		endif 
	else
		if (diffusion.eq.'CDS2') then
		  call diffw_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
		elseif(diffusion.eq.'COM4') then
		  call diffw_com4(dnew,Szr,Spz,Szz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
		endif
	endif 
		
		
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 call advecw_driftfluxCDS2(dnew,driftfluxforce_calfac*sumWkm,
     &   rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
	endif
	
	   if (applyVOF.eq.1) then ! VOF fraction apply real density only here for correct determination gravity
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-gz*(0.5*(rnew(i,j,k)+rnew(i,j,k+1))-rho_b)/(0.5*(rnew(i,j,k)+rnew(i,j,k+1))/rho_b)
     1     +Ppropz(i,j,k)
            dWdt(i,j,k)= dWdt(i,j,k)+ !0.5*(rnew(i,j,k)+rnew(i,j,k+1))+
     1     dt*(      facAB3_1*dnew(i,j,k)+
     1               facAB3_2*wz(i,j,k)+
     1               facAB3_3*wzold(i,j,k)+
     1               facAB3_4*wzolder(i,j,k)+
     1               dnew2(i,j,k)	 ) !AB3v for advection and body forces, CN for diffusion when CNdiffz=1 or 2	 
!     1          dt*(1.5*dnew(i,j,k)-0.5*wz(i,j,k)) -gz*dt*(0.75*(rnew(i,j,k)+rnew(i,j,k+1))-0.25*(rold(i,j,k)+rold(i,j,k+1))-rho_b)
		    wzolder(i,j,k) = wzold(i,j,k)
		    wzold(i,j,k)   = wz(i,j,k)
		    wz(i,j,k)      = dnew(i,j,k)
            enddo
         enddo
      enddo	
	   else 
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-gz*(0.5*(rnew(i,j,k)+rnew(i,j,k+1))-rho_b)
     1     +Ppropz(i,j,k)
            dWdt(i,j,k)= dWdt(i,j,k)+ !0.5*(rnew(i,j,k)+rnew(i,j,k+1))+
     1     dt*(      facAB3_1*dnew(i,j,k)+
     1               facAB3_2*wz(i,j,k)+
     1               facAB3_3*wzold(i,j,k)+
     1               facAB3_4*wzolder(i,j,k)+
     1              dnew2(i,j,k)	 ) !AB3v for advection and body forces, CN for diffusion when CNdiffz=1 or 2	 
!     1          dt*(1.5*dnew(i,j,k)-0.5*wz(i,j,k)) -gz*dt*(0.75*(rnew(i,j,k)+rnew(i,j,k+1))-0.25*(rold(i,j,k)+rold(i,j,k+1))-rho_b)
		    wzolder(i,j,k) = wzold(i,j,k)
		    wzold(i,j,k)   = wz(i,j,k)
		    wz(i,j,k)      = dnew(i,j,k)
            enddo
         enddo
      enddo		
	   endif
  
	  pold=p+pold    !what is called p here was dp in reality, now p is 			--> P^n 

		  !IF (continuity_solver.eq.35) THEN
			!pold1=pold2 !-->P^n-2
			!pold2=pold3 !-->P^n-1
			!pold3=pold !--> P^n		  
		  IF (continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
			pold1=pold2 !-->P^n-2
			pold2=pold3 !-->P^n-1
			pold3=pold !--> P^n
		  ENDIF 
		  IF (pres_in_predictor_step.eq.0) THEN 
		    pold=0.
		  ENDIF
		  IF (applyVOF.eq.1) THEN 
			pold=rho_b*pold/rnew(1:imax,1:jmax,1:kmax)
		  ENDIF 
		!ENDIF 
	  
      call pshiftb(pold,pplus) !,rank,imax,jmax,kmax,px)
c********************************************************************
c     CALCULATE pressure-gradient with old pressure:
c********************************************************************
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -dt * ( pold(i+1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
	  !!! do nothing on imax and i=1, so p(imax+1)=p(imax) --> dp/dn=0 at outflow and inflow 
	  if (periodicx.eq.1) then
	  !!! !periodic bc in x: 
	    i = imax
	      do k=1,kmax
		    do j=1,jmax
		      dUdt(i,j,k) = dUdt(i,j,k) -dt * ( pold(1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		    enddo
	      enddo
	  endif
	
      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * (phip(j+1)-phip(j)) )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
      	  endif
	endif
      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          dWdt(i,j,k) = dWdt(i,j,k) -dt * ( pold(i,j,k+1) - pold(i,j,k) ) / dz
          enddo
        enddo
      enddo


      IF (CNdiffz.eq.1.or.CNdiffz.eq.2.or.CNdiffz.eq.11.or.CNdiffz.eq.12.or.CNdiffz.eq.31) THEN !CN semi implicit treatment diff-z:

		  ! add source terms for cylindrical coordinate system (following PhD dissertation Pourqeui 1994 TU p. 82)
		  ! --> LdW 24-5-2021 choice to do it explicitly at start to avoid issues with boundary conditions at 0,i1/0,j1
		  ! --> in diff_compact.f the diffusion terms are discretised a bit differently (via the stress tensor), but the below source terms are included there as well
		  IF (CNdiffz.eq.12.or.CNdiffz.eq.11.or.CNdiffz.eq.31) THEN 
		 ! IF (.false.) THEN 
		  !IF (CNdiffz.eq.12) THEN !for now 11 uses 0.5*DIFF(U^n+1-U^n)+1*DIFF(U^n) and therefore no need for this explicit source term
		  do k=1,kmax ! add source term [-u/r^2 - 2/r^2*dvdphi] to U and [-v/r^2 + 2/r^2*dudphi] to V 
           do i=1,imax
			ip=i+1 
            do j=1,jmax 
				jm=j-1
				jp=j+1 !MIN(px*jmax+1,j+1)
				dUdt(i,j,k)=dUdt(i,j,k)-CNdiff_factor*dt*(0.5*(ekm(i,j,k)+ekm(ip,j,k)))*
     &          (Unew(i,j,k)+2.*(Vnew(i,j,k)-Vnew(i,jm,k))/(phiv(j)-phiv(jm)))/(Ru(i)*Ru(i)) 				
				dVdt(i,j,k)=dVdt(i,j,k)-CNdiff_factor*dt*(0.5*(ekm(i,j,k)+ekm(i,jp,k)))*
     &          (Vnew(i,j,k)-2.*(Unew(i,jp,k)-Unew(i,j,k))/(phip(jp)-phip(j)))/(Rp(i)*Rp(i)) 	 
			enddo
           enddo
		  enddo	
		  ENDIF 		
		
		call bound_rhoU(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
     & 	Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)

! apply diffusion on dUdt/rhU with ekm instead on dUdt with ekm/drdt which because of staggering is not the same and in explicit diff the former is done:
		dUdt = dUdt/rhU 
		dVdt = dVdt/rhV 
		dWdt = dWdt/rhW 
		
		
		IF (CNdiffz.eq.11) THEN
		  IF (rank.eq.0) THEN 
			CALL SYSTEM_CLOCK(COUNT=clock)
			CALL RANDOM_SEED(size = n)
			ALLOCATE(seed(n))
			CALL SYSTEM_CLOCK(COUNT=clock)
			seed = clock + 37 * (/ (i - 1, i = 1, n) /)
			CALL RANDOM_SEED(PUT = seed)		  
			call random_number(dir_order_real) ! uniform distribution 0,1
			dir_order_real=dir_order_real*6. ! uniform distribution 0,6
			dir_order=CEILING(dir_order_real) ! either 1,2,3,4,5,6
		  ENDIF 
		  call mpi_bcast(dir_order,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

		  !randomizer in order of implicit steps to make sure that this will not impact quality of outcome after many time steps
!			IF (dir_order.eq.1) THEN 
				! CNdiff_factor acts as theta in CN --> 0.5 gives 2nd order dt Douglas-Gunn ADI and 1 gives 1st order Douglass-Rachford; factor 0.5-1 gives blend
				! Douglas-Gunn and Douglass-Rachford ADI first take 100% explicit diffusion of 2 dirs not implicit
				call diffuvw_ydir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + (1.-CNdiff_factor)*dt*dnew/rhU  
				dVdt = dVdt + (1.-CNdiff_factor)*dt*dnew2/rhV
				dWdt = dWdt + (1.-CNdiff_factor)*dt*dold/rhW				
				call diffuvw_zdir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + dt*dnew/rhU  
				dVdt = dVdt + dt*dnew2/rhV
				dWdt = dWdt + dt*dold/rhW
				call diffuvw_xdir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + dt*dnew/rhU  
				dVdt = dVdt + dt*dnew2/rhV
				dWdt = dWdt + dt*dold/rhW				
!!!				if (CNdiff_ho>0.and.CNdiff_ho<1.01.and.CNdiff_factor>0.501) then ! high order correction to arrive almost at fully 3D theta-CN based on unew-uold from Douglas-Kim 2001 :
!!!					call diffuvw_xdir_CDS2_CNexpl(ddd1,eee1,fff1,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
!!!					call diffuvw_ydir_CDS2_CNexpl(ddd2,eee2,fff2,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
!!!					call diffuvw_zdir_CDS2_CNexpl(ddd3,eee3,fff3,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
!!!					ddd1=-ddd1/rhU 
!!!					ddd2=-ddd2/rhU
!!!					ddd3=-ddd3/rhU
!!!					eee1=-eee1/rhV 
!!!					eee2=-eee2/rhV
!!!					eee3=-eee3/rhV
!!!					fff1=-fff1/rhW 
!!!					fff2=-fff2/rhW
!!!					fff3=-fff3/rhW					
!!!					! add correction to RHS of first step Douglas-Gunn according to Douglas-Kim 2001 to obtain >2nd order in time and get close to accuracy full 3D CN
!!!					! in my own derivation last term AxAyAz should be 1/8dt^3 instead of 2/8 as mentioned in Douglas-Kim 2001 
!!!					dUdt = dUdt+CNdiff_factor**2*dt*dt*(ddd1*ddd2+ddd1*ddd3+ddd2*ddd3)+CNdiff_factor**3*dt*dt*dt*(ddd1*ddd2*ddd3)  
!!!					dVdt = dVdt+CNdiff_factor**2*dt*dt*(eee1*eee2+eee1*eee3+eee2*eee3)+CNdiff_factor**3*dt*dt*dt*(eee1*eee2*eee3) 
!!!					dWdt = dWdt+CNdiff_factor**2*dt*dt*(fff1*fff2+fff1*fff3+fff2*fff3)+CNdiff_factor**3*dt*dt*dt*(fff1*fff2*fff3) 
!!!					call diffuvw_xdir_CDS2_CNexpl(ddd1,eee1,fff1,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
!!!					call diffuvw_ydir_CDS2_CNexpl(ddd2,eee2,fff2,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
!!!					call diffuvw_zdir_CDS2_CNexpl(ddd3,eee3,fff3,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
!!!					ddd1=ddd1/rhU 
!!!					ddd2=ddd2/rhU
!!!					ddd3=ddd3/rhU
!!!					eee1=eee1/rhV 
!!!					eee2=eee2/rhV
!!!					eee3=eee3/rhV
!!!					fff1=fff1/rhW 
!!!					fff2=fff2/rhW
!!!					fff3=fff3/rhW					
!!!					! add correction to RHS of first step Douglas-Gunn according to Douglas-Kim 2001 to obtain >2nd order in time and get close to accuracy full 3D CN
!!!					! in my own derivation last term AxAyAz should be 1/8dt^3 instead of 2/8 as mentioned in Douglas-Kim 2001 
!!!					dUdt = dUdt+(1.-CNdiff_factor)**2*dt*dt*(ddd1*ddd2+ddd1*ddd3+ddd2*ddd3)+(1.-CNdiff_factor)**3*dt*dt*dt*(ddd1*ddd2*ddd3)  
!!!					dVdt = dVdt+(1.-CNdiff_factor)**2*dt*dt*(eee1*eee2+eee1*eee3+eee2*eee3)+(1.-CNdiff_factor)**3*dt*dt*dt*(eee1*eee2*eee3) 
!!!					dWdt = dWdt+(1.-CNdiff_factor)**2*dt*dt*(fff1*fff2+fff1*fff3+fff2*fff3)+(1.-CNdiff_factor)**3*dt*dt*dt*(fff1*fff2*fff3) 					
!!!				elseif (CNdiff_ho>0.and.CNdiff_ho<1.01) then !CNdiff_factor=0.5
!!!					call diffuvw_xdir_CDS2_CNexpl(ddd1,eee1,fff1,Unew-Uold,Vnew-Vold,Wnew-Wold,ib,ie,jb,je,kb,ke)
!!!					call diffuvw_ydir_CDS2_CNexpl(ddd2,eee2,fff2,Unew-Uold,Vnew-Vold,Wnew-Wold,ib,ie,jb,je,kb,ke)
!!!					call diffuvw_zdir_CDS2_CNexpl(ddd3,eee3,fff3,Unew-Uold,Vnew-Vold,Wnew-Wold,ib,ie,jb,je,kb,ke)
!!!					ddd1=-ddd1/rhU 
!!!					ddd2=-ddd2/rhU
!!!					ddd3=-ddd3/rhU
!!!					eee1=-eee1/rhV 
!!!					eee2=-eee2/rhV
!!!					eee3=-eee3/rhV
!!!					fff1=-fff1/rhW 
!!!					fff2=-fff2/rhW
!!!					fff3=-fff3/rhW					
!!!					! add correction to RHS of first step Douglas-Gunn according to Douglas-Kim 2001 to obtain >2nd order in time and get close to accuracy full 3D CN
!!!					! in my own derivation last term AxAyAz should be 1/8dt^3 instead of 2/8 as mentioned in Douglas-Kim 2001 
!!!					dUdt = dUdt+0.25*dt*dt*(ddd1*ddd2+ddd1*ddd3+ddd2*ddd3)+0.125*dt*dt*dt*(ddd1*ddd2*ddd3)  
!!!					dVdt = dVdt+0.25*dt*dt*(eee1*eee2+eee1*eee3+eee2*eee3)+0.125*dt*dt*dt*(eee1*eee2*eee3) 
!!!					dWdt = dWdt+0.25*dt*dt*(fff1*fff2+fff1*fff3+fff2*fff3)+0.125*dt*dt*dt*(fff1*fff2*fff3) 				
!!!				
!!!				elseif (CNdiff_ho>1.01) then ! high order correction based on dudt-Unew with dudt determined with expl diffusion (my own idea)
!!!					dUdt = dUdt + CNdiff_factor*dt*dnew/rhU  
!!!					dVdt = dVdt + CNdiff_factor*dt*dnew2/rhV
!!!					dWdt = dWdt + CNdiff_factor*dt*dold/rhW	
!!!					call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!!!     & 				Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)					
!!!					call diffuvw_xdir_CDS2_CNexpl(ddd1,eee1,fff1,dUdt-Unew,dVdt-Vnew,dWdt-Wnew,ib,ie,jb,je,kb,ke)
!!!					call diffuvw_ydir_CDS2_CNexpl(ddd2,eee2,fff2,dUdt-Unew,dVdt-Vnew,dWdt-Wnew,ib,ie,jb,je,kb,ke)
!!!					call diffuvw_zdir_CDS2_CNexpl(ddd3,eee3,fff3,dUdt-Unew,dVdt-Vnew,dWdt-Wnew,ib,ie,jb,je,kb,ke)
!!!					ddd1=-ddd1/rhU 
!!!					ddd2=-ddd2/rhU
!!!					ddd3=-ddd3/rhU
!!!					eee1=-eee1/rhV 
!!!					eee2=-eee2/rhV
!!!					eee3=-eee3/rhV
!!!					fff1=-fff1/rhW 
!!!					fff2=-fff2/rhW
!!!					fff3=-fff3/rhW					
!!!					! add correction to RHS of first step Douglas-Gunn according to Douglas-Kim 2001 to obtain >2nd order in time and get close to accuracy full 3D CN
!!!					! in my own derivation last term AxAyAz should be 1/8dt^3 instead of 2/8 as mentioned in Douglas-Kim 2001 
!!!					dUdt=dUdt+0.25*dt*dt*(ddd1*ddd2+ddd1*ddd3+ddd2*ddd3)+0.125*dt*dt*dt*(ddd1*ddd2*ddd3)-CNdiff_factor*dt*dnew/rhU    
!!!					dVdt=dVdt+0.25*dt*dt*(eee1*eee2+eee1*eee3+eee2*eee3)+0.125*dt*dt*dt*(eee1*eee2*eee3)-CNdiff_factor*dt*dnew2/rhV 
!!!					dWdt=dWdt+0.25*dt*dt*(fff1*fff2+fff1*fff3+fff2*fff3)+0.125*dt*dt*dt*(fff1*fff2*fff3)-CNdiff_factor*dt*dold/rhW	 				
!!!				endif 
				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)	
				CALL diffuvw_ydir_CDS2_CNimpl
				
				! Douglas-Gunn ADI after each direction 50% implicit is finished subtract 50% explicit 
				! diffuvw_xdir_CDS2_CNexpl was executed as last direction; therefore dnew,dnew2,dold are already correct for xdir_expl:
				!call diffuvw_xdir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke) 
				dUdt = dUdt - CNdiff_factor*dt*dnew/rhU  
				dVdt = dVdt - CNdiff_factor*dt*dnew2/rhV
				dWdt = dWdt - CNdiff_factor*dt*dold/rhW
				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)		 
				CALL diffuvw_xdir_CDS2_CNimpl 
				
				! Douglas-Gunn ADI after each direction 50% implicit is finished subtract 50% explicit 
				call diffuvw_zdir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
				dUdt = dUdt - CNdiff_factor*dt*dnew/rhU  
				dVdt = dVdt - CNdiff_factor*dt*dnew2/rhV
				dWdt = dWdt - CNdiff_factor*dt*dold/rhW			
				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
				CALL diffuvw_zdir_CDS2_CNimpl
!			ELSEIF (dir_order.eq.2) THEN  
!				CALL diffuvw_ydir_CDS2_CNimpl
!				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
!				CALL diffuvw_zdir_CDS2_CNimpl
!				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
!				CALL diffuvw_xdir_CDS2_CNimpl 
!			ELSEIF (dir_order.eq.3) THEN 
!				CALL diffuvw_xdir_CDS2_CNimpl
!				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
!				CALL diffuvw_ydir_CDS2_CNimpl
!				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
!				CALL diffuvw_zdir_CDS2_CNimpl			
!			ELSEIF (dir_order.eq.4) THEN 
!				CALL diffuvw_xdir_CDS2_CNimpl
!				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
!				CALL diffuvw_zdir_CDS2_CNimpl			
!				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
!				CALL diffuvw_ydir_CDS2_CNimpl
!			ELSEIF (dir_order.eq.5) THEN 
!				CALL diffuvw_zdir_CDS2_CNimpl
!				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
!				CALL diffuvw_ydir_CDS2_CNimpl
!				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
!				CALL diffuvw_xdir_CDS2_CNimpl						
!			ELSE 
!				CALL diffuvw_zdir_CDS2_CNimpl
!				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
!				CALL diffuvw_xdir_CDS2_CNimpl									
!				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
!				CALL diffuvw_ydir_CDS2_CNimpl				
!			ENDIF 
		ELSEIF (CNdiffz.eq.12) THEN
				! Chang 1991 ADI method 
				CNdiffz=0 !to get 100% explicit diff from diffuvw_.dir_CDS2_CNexpl
				call diffuvw_xdir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + dt*CNdiff_factor/3.*dnew/rhU  
				dVdt = dVdt + dt*CNdiff_factor/3.*dnew2/rhV
				dWdt = dWdt + dt*CNdiff_factor/3.*dold/rhW
				call diffuvw_zdir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + dt*CNdiff_factor/3.*dnew/rhU  
				dVdt = dVdt + dt*CNdiff_factor/3.*dnew2/rhV
				dWdt = dWdt + dt*CNdiff_factor/3.*dold/rhW			
				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)	
				CNdiffz=12 !to get correct implicit diff from diffuvw_.dir_CDS2_CNimpl
				CALL diffuvw_ydir_CDS2_CNimpl
				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)				
				! Chang 1991 ADI method 
				CNdiffz=0 !to get 100% explicit diff from diffuvw_.dir_CDS2_CNexpl
				call diffuvw_ydir_CDS2_CNexpl(dnew,dnew2,dold,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
				call diffuvw_zdir_CDS2_CNexpl(ddd1,ddd2,ddd3,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + dt*CNdiff_factor/3.*dnew/rhU  
				dVdt = dVdt + dt*CNdiff_factor/3.*dnew2/rhV
				dWdt = dWdt + dt*CNdiff_factor/3.*dold/rhW
				dUdt = dUdt + dt*CNdiff_factor/3.*ddd1/rhU  
				dVdt = dVdt + dt*CNdiff_factor/3.*ddd2/rhV
				dWdt = dWdt + dt*CNdiff_factor/3.*ddd3/rhW	
				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
				CNdiffz=12 !to get correct implicit diff from diffuvw_.dir_CDS2_CNimpl	 
				CALL diffuvw_xdir_CDS2_CNimpl 
				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)				
				! Chang 1991 ADI method 
				CNdiffz=0 !to get 100% explicit diff from diffuvw_.dir_CDS2_CNexpl
				call diffuvw_xdir_CDS2_CNexpl(dnew,dnew2,dold,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
				call diffuvw_ydir_CDS2_CNexpl(ddd1,ddd2,ddd3,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + dt*CNdiff_factor/3.*dnew/rhU  
				dVdt = dVdt + dt*CNdiff_factor/3.*dnew2/rhV
				dWdt = dWdt + dt*CNdiff_factor/3.*dold/rhW
				dUdt = dUdt + dt*CNdiff_factor/3.*ddd1/rhU  
				dVdt = dVdt + dt*CNdiff_factor/3.*ddd2/rhV
				dWdt = dWdt + dt*CNdiff_factor/3.*ddd3/rhW	
				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
				CNdiffz=12 !to get correct implicit diff from diffuvw_.dir_CDS2_CNimpl
				CALL diffuvw_zdir_CDS2_CNimpl
		ENDIF 
		
		IF (CNdiffz.eq.31) THEN
			IF (CNdiff_factor<1.) THEN 
				call diffuvw_ydir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + (1.-CNdiff_factor)*dt*dnew/rhU  
				dVdt = dVdt + (1.-CNdiff_factor)*dt*dnew2/rhV
				dWdt = dWdt + (1.-CNdiff_factor)*dt*dold/rhW				
				call diffuvw_zdir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + (1.-CNdiff_factor)*dt*dnew/rhU  
				dVdt = dVdt + (1.-CNdiff_factor)*dt*dnew2/rhV
				dWdt = dWdt + (1.-CNdiff_factor)*dt*dold/rhW
				call diffuvw_xdir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + (1.-CNdiff_factor)*dt*dnew/rhU  
				dVdt = dVdt + (1.-CNdiff_factor)*dt*dnew2/rhV
				dWdt = dWdt + (1.-CNdiff_factor)*dt*dold/rhW		
			ENDIF 
			call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 		Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
			CALL diffuvw_CDS2_3DCNimpl 
			
		ENDIF 
		
		
		IF (CNdiffz.eq.1.or.CNdiffz.eq.2) THEN 
			CALL diffuvw_zdir_CDS2_CNimpl
		ENDIF 
		
		dUdt = dUdt*rhU 
		dVdt = dVdt*rhV
		dWdt = dWdt*rhW
      ENDIF


		!IF (applyVOF.eq.1) THEN 
		! IF (continuity_solver.eq.35) THEN
		!  pold=pold*rnew(1:imax,1:jmax,1:kmax)
		!  rnew=1.
		! ENDIF
		!ENDIF 
		IF (applyVOF.eq.1) THEN 
			pold=pold/rho_b*rnew(1:imax,1:jmax,1:kmax)  
			rnew=rho_b
		ENDIF 
		

      return
      end

      subroutine RK3(ib,ie,jb,je,kb,ke)
c
c    Performes time integration with three step RK (RungeKutta) scheme
c    dependent on the used constants the time stepping can be up to third order accurate
c    the implementation of the pressure term is only second order accurate in time, 
c    so overall accuracy is restricted to second order 
c    
c    The RK3 constants can be written in Butcher notation like:
c
c    0  0    |
c    1  c2   | a21
c    2  c3   | a31   a32
c      ------|------------------
c            | b1    b2    b3    
c
c   Here the constants derived by Wray are used (Spalart, Moser and Rogers, JCP 96, 1991) 
c   It is 3th order accurate, CFL max = 1.73  (see routine chkdt)
c    0  0    |
c    1  8/15 | 8/15
c    2  2/3  | 1/4   5/12
c      ------|------------------
c            | 1/4   0     3/4   
c
c                    time_n
c    k1 = ADV+DIFF(U) 
c
c                        time_n+8/15dt
c    k2 = ADV+DIFF(pred1) 
c
c                        time_n+2/3dt
c    k3 = ADV+DIFF(pred2)
c
c     n+1    n
c    U    = U  +1/4*dt*k1 +3/4*dt*k3
c
c    with predictors:
c                n
c    pred1    = U  +8/15*dt*k1
c                n
c    pred2    = U  +1/4*dt*k1 +5/12*dt*k2
c

      USE nlist
      USE sediment

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer  ib,ie,jb,je,kb,ke,n,t,kpp,km,kp,me
      real     wsed(nfrac,0:i1,0:j1,0:k1),sumWkm(0:i1,0:j1,0:k1)
      real     c2,c3,a21,a31,a32,b1,b2,b3,cn1,cn2,cn3
	real   b1tvd,b2tvd,b3tvd
!	real     Diffcof(0:i1,0:j1,0:k1)
	real	 invrho_b
	real*8 pplus(imax,kmax),pplus2(imax,kmax)
!	real aaa(1:kmax),bbb(1:kmax),ccc(1:kmax),ekm_min,ekm_plus,rhss(1:kmax)
	real aaa(0:k1),bbb(0:k1),ccc(0:k1),ekm_min,ekm_plus,rhss(0:k1)
	real utr(0:i1,0:j1,0:k1),vtr(0:i1,0:j1,0:k1),wtr(0:i1,0:j1,0:k1)
!	real dUdt1(0:i1,0:j1,0:k1),dVdt1(0:i1,0:j1,0:k1),dWdt1(0:i1,0:j1,0:k1)
	real ws(0:i1,0:j1,0:k1),gvector 

	real   k1c(nfrac,0:i1,0:j1,0:k1)
	real   k1cbot(nfrac,0:i1,0:j1)
	real   k1u(0:i1,0:j1,0:k1)
	real   k1v(0:i1,0:j1,0:k1)
	real   k1w(0:i1,0:j1,0:k1)
	real   k2c(nfrac,0:i1,0:j1,0:k1)
	real   k2cbot(nfrac,0:i1,0:j1)
	real   k2u(0:i1,0:j1,0:k1)
	real   k2v(0:i1,0:j1,0:k1)
	real   k2w(0:i1,0:j1,0:k1)
	real   k3c(nfrac,0:i1,0:j1,0:k1)
	real   k3cbot(nfrac,0:i1,0:j1)
	real   k3u(0:i1,0:j1,0:k1)
	real   k3v(0:i1,0:j1,0:k1)
	real   k3w(0:i1,0:j1,0:k1)

!! predictors:	
	real   dcdt1bot(nfrac,0:i1,0:j1)
	real   dcdt2bot(nfrac,0:i1,0:j1)
	
!	dUdt1=dudt
!	dVdt1=dvdt
!	dWdt1=dwdt

	me = momentum_exchange_obstacles !100 or 101 or 111 means no advection momentum exchange with bed-cells with HYB6

!	Diffcof=ekm/Sc/Rnew 

	dcdt=0.
	sumWkm=0.

!!	Wray RK3 constants:
	c2=8./15.
	c3=2./3.
	a21=8./15.
	a31=1./4.
	a32=5./12.
	b1=1./4.
	b2=0.
	b3=3./4.
!!!!	Wray Crack Nicolson constants for pressure:
	cn1=8./15.
	cn2=2./15.
	cn3=1./3.
!!!	TVD RK3 constants for adv-diff concentration belonging to Wray (found by accident, are they correct?):
	b1tvd=1./4.
	b2tvd=0.
	b3tvd=3./4.


!!!	TVD RK3 constants (TVD with CFLmax=1)
!	c2 = 1.
!	c3=1./2.
!	a21=1.
!	a31=1./4.
!	a32=1./4.
!	b1=1./6.
!	b2=1./6.
!	b3=2./3.
!!!	Wray Crack Nicolson constants for pressure:
!	cn1=1.
!	cn2=-0.5
!	cn3=0.5

!!!!	TVD RK3 constants (TVD with CFLmax=1.5, in practice CFLmax=1 for correct results)
!	c2 = 1./5.
!	c3=1./2.
!	a21=1./5.
!	a31=0.
!	a32=1./2. 
!	b1=0.
!	b2=0.
!	b3=1.
!!!!	Crack Nicolson constants for pressure:
!	cn1=1./5.
!	cn2=3./10.
!	cn3=1./2.
	pold=p+pold    !what is called p here was dp in reality, now p is 			--> P^n 
	!IF (applyVOF.eq.1) THEN 
	  !IF (continuity_solver.eq.35) THEN
		!pold1=pold2 !-->P^n-2
		!pold2=pold3 !-->P^n-1
		!pold3=pold !--> P^n		  
	  !  CALL state(cnew,rnew)		  
	  !  pold=pold/rnew(1:imax,1:jmax,1:kmax)
	  IF (continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
		pold1=pold2 !-->P^n-2
		pold2=pold3 !-->P^n-1
		pold3=pold !--> P^n
	  ENDIF 
	  IF (pres_in_predictor_step.eq.0) THEN 
		pold=0.
	  ENDIF
		  
        call pshiftb(pold,pplus) !,rank,imax,jmax,kmax,px)
c********************************************************************
c     CALCULATE k1 and predictor 1
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(cnew,Wnew,wsed,rnew,1,kmax,sumWkm,cn1*dt,dz)
	       !! Set boundary conditions jet in:
	       do t=1,tmax_inPpunt
	 	i=i_inPpunt(t)
	 	j=j_inPpunt(t)
		do n=1,nfrac
	    	  do k=kmax,kmax !-kjet,kmax
		 	Wsed(n,i,j,k)=Wnew(i,j,k)
		  enddo
		enddo
		if (LOA>0.and.outflow_overflow_down.eq.1) then
		 do n=1,nfrac
	    	  do k=kmax-kjet+1,kmax 
		 	Wsed(n,i,j,k)=Wnew(i,j,k)
		  enddo
		 enddo	
		endif	
	       enddo
	  else
	    do n=1,nfrac
		wsed(n,:,:,:)=wnew
	    enddo
	  endif
	  do n=1,nfrac
	      utr=Unew
	      vtr=Vnew
	      wtr=Wsed(n,:,:,:)
			if (settling_along_gvector.eq.1) then
				gvector=sqrt(gx**2+gy**2+gz**2)
				ws = Wsed(n,:,:,:)-Wnew 
				utr=utr+gx/gvector*ws 
				vtr=vtr+gy/gvector*ws
				wtr=wtr-ws+gz/gvector*ws
			endif 		  
	      call make_UtransC_zeroIBM(utr,vtr,wtr)
		  if (advec_conc.eq.'NVD') then
	      call advecc_NVD(k1c(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn1*dt,rank,px,periodicx,periodicy)	 
          elseif (advec_conc.eq.'VLE') then	 
	      call advecc_VLE(k1c(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn1*dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)
          elseif (advec_conc.eq.'SBE') then	 
	      call advecc_SBE(k1c(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn1*dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)	
          elseif (advec_conc.eq.'VL2') then	 
	      call advecc_VL2_nocfl(k1c(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
          elseif (advec_conc.eq.'SB2') then	 
	      call advecc_SB2_nocfl(k1c(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)	 
		  elseif (advec_conc.eq.'ARO') then	
	      call advecc_ARO(k1c(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn1*dt,rank,px,periodicx,periodicy)	 
          endif
		if (Sc<1.e18) then
	      call diffc_CDS2 (k1c(n,:,:,:),Cnew(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)
	    endif
!		dcdt(n,:,:,:) =cnew(n,:,:,:) + a21*dt*(k1c(n,:,:,:)) !pred1
		dcdt(n,:,:,:) =cnew(n,:,:,:) + cn1*dt*(k1c(n,:,:,:)) !pred1
           IF (CNdiffz.eq.1) THEN !CN semi implicit treatment diff-z:
	     call bound_c(dcdt(n,:,:,:),frac(n)%c,n,0.) ! bc start CN-diffz k1 RK3
	     do k=k1,k1 !-kjet,k1
	       do t=1,tmax_inPpunt
	         i=i_inPpunt(t)
	         j=j_inPpunt(t)
	         dcdt(n,i,j,k) = dcdt(n,i,j,kmax) ! No diffusion over vertical inflow boundary, this line is needed for exact influx
	       enddo
	     enddo
             do j=1,jmax
               do i=1,imax
                 do k=0,k1
	           km=MAX(0,k-1)
	           kp=MIN(k1,k+1)
	           kpp=MIN(k1,k+2)
	           ekm_min=0.5*(Diffcof(i,j,km)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,km),fc_global(i,j+rank*jmax,k))
	           ekm_plus=0.5*(Diffcof(i,j,kp)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,kp),fc_global(i,j+rank*jmax,k))
	           aaa(k)=-0.5*ekm_min*cn1*dt/dz**2
	           bbb(k)=1.+0.5*(ekm_min+ekm_plus)*cn1*dt/dz**2
	           ccc(k)=-0.5*ekm_plus*cn1*dt/dz**2
                 enddo
		 aaa(0)=0.
		 aaa(k1)=0.
		 ccc(0)=0.
		 ccc(k1)=0.
		 bbb(0)=1.
		 bbb(k1)=1.
	         rhss=dcdt(n,i,j,0:k1)
	         CALL solve_tridiag(dcdt(n,i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
               enddo
             enddo
	   ENDIF !endif semi-implicit CN
          enddo
		  call slipvelocity_bed(cnew,wnew,wsed,rnew,sumWkm,cn1*dt,dz) !driftflux_force must be calculated with settling velocity at bed
      	  if (interaction_bed>0) then
	    CALL erosion_deposition(dcdt,dcdt1bot,unew,vnew,wnew,rnew,cnew,cnewbot,cn1*dt,dz) !first two vars are adjusted
		do n=1,nfrac 
		  call bound_cbot(dcdt1bot(n,:,:))
		enddo			
		!! advec concentration in bed with velocity TSHD:
		IF (ABS(U_TSHD)>1e-6) THEN 
		  do n=1,nfrac 
	      	    call adveccbot_TVD(k1cbot(n,:,:),dcdt1bot(n,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,phiv,phipt,dz,
     +            	i1,j1,ib,ie,jb,je,cn1*dt,rank,px,periodicx,periodicy)
		    dcdt1bot(n,:,:)= dcdt1bot(n,:,:) + cn1*dt*k1cbot(n,:,:) ! pred1	
			call bound_cbot(dcdt1bot(n,:,:))			
		  enddo 
		ENDIF 		
	    !! for kn1: erosion_deposition must be after advecc_TVD and after dcdt update, because cnew and dcdt are two separate vars
		if (interaction_bed.eq.4.or.interaction_bed.eq.6) then
			call advec_update_Clivebed(dcdt,dcdt1bot,cn1*dt) 	
		endif			
	    do n=1,nfrac 
	        call bound_cbot(dcdt1bot(n,:,:))
	    enddo
      	  endif
	    do n=1,nfrac 
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n,cn1*dt) ! bc after erosion_deposition k1 RK3
	    enddo
	endif

	if (diffusion.eq.'COM4') then
		call stress_terms(Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	endif

      	if (convection.eq.'CDS2') then
      call advecu_CDS2(k1u,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(k1u,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(k1u,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecu_CDS4(k1u,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(k1u,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(k1u,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2)	 
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(k1u,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd,numdiff2)	  
	elseif(convection.eq.'uTVD') then
      call advecu_TVD(k1u,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif
	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (k1u,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(k1u,Srr,Spr,Szr,Spp,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

	if(convection.eq.'CDS2') then
      call advecv_CDS2(k1v,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(k1v,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(k1v,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecv_CDS4(k1v,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(k1v,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(k1v,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(k1v,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecv_TVD(k1v,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)		  
	endif
	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (k1v,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(k1v,Spr,Spp,Spz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      	if (convection.eq.'CDS2') then
      call advecw_CDS2(k1w,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(k1w,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(k1w,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecw_CDS4(k1w,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(k1w,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(k1w,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(k1w,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecw_TVD(k1w,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)		  
	endif
	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (k1w,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(k1w,Szr,Spz,Szz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 call advecw_driftfluxCDS2(k1w,driftfluxforce_calfac*sumWkm,
     &   rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
	endif
	
	

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    k1u(i,j,k) = k1u(i,j,k)-(gx*cos_u(j)+gy*sin_u(j))*(0.5*(rnew(i,j,k)+rnew(i+1,j,k))-rho_b)
            dUdt(i,j,k)=Unew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i+1,j,k))+(a21*k1u(i,j,k)	! pred1
     1     +(Ppropx(i,j,k))*c2)*dt 
            enddo
         enddo
      enddo
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    k1v(i,j,k) = k1v(i,j,k)-(-gx*sin_v(j)+gy*cos_v(j))*(0.5*(rnew(i,j,k)+rnew(i,j+1,k))-rho_b)
            dVdt(i,j,k)=Vnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j+1,k))+(a21*k1v(i,j,k) 	! pred1
     1     +(Ppropy(i,j,k))*c2)*dt 
            enddo
         enddo
      enddo
       do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    k1w(i,j,k)=k1w(i,j,k)-gz*(0.5*(rnew(i,j,k)+rnew(i,j,k+1))-rho_b)
            dWdt(i,j,k)=Wnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j,k+1))+(a21*k1w(i,j,k) 	! pred1
     1     +(Ppropz(i,j,k))*c2)*dt 
            enddo
         enddo
       enddo

      IF (CNdiffz.eq.1) THEN !CN semi implicit treatment diff-z:
	if (nfrac>0) then
	  call state(dcdt,drdt) ! determine drdt1 with pred1
	endif
      call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),0,time_n+c2*dt,(1.-c2)*Ub1old+c2*Ub1new,(1.-c2)*Vb1old+c2*Vb1new,
     & (1.-c2)*Wb1old+c2*Wb1new,(1.-c2)*Ub2old+c2*Ub2new,(1.-c2)*Vb2old+c2*Vb2new,(1.-c2)*Wb2old+c2*Wb2new,
     $ (1.-c2)*Ub3old+c2*Ub3new,(1.-c2)*Vb3old+c2*Vb3new,(1.-c2)*Wb3old+c2*Wb3new)
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))		
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))		
	    aaa(k)=-0.5*ekm_min*c2*dt/dz**2/(0.5*(rnew(i,j,km)+rnew(i+1,j,km)))
	    bbb(k)=1.+0.5*(ekm_min+ekm_plus)*c2*dt/dz**2/(0.5*(rnew(i,j,k)+rnew(i+1,j,k)))
	    ccc(k)=-0.5*ekm_plus*c2*dt/dz**2/(0.5*(rnew(i,j,kp)+rnew(i+1,j,kp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dUdt(i,j,0:k1)
	    CALL solve_tridiag(dUdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
      enddo
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))		
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))		
	    aaa(k)=-0.5*ekm_min*c2*dt/dz**2/(0.5*(rnew(i,j,km)+rnew(i,j+1,km)))
	    bbb(k)=1.+0.5*(ekm_min+ekm_plus)*c2*dt/dz**2/(0.5*(rnew(i,j,k)+rnew(i,j+1,k)))
	    ccc(k)=-0.5*ekm_plus*c2*dt/dz**2/(0.5*(rnew(i,j,kp)+rnew(i,j+1,kp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dVdt(i,j,0:k1)
	    CALL solve_tridiag(dVdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
      enddo

       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=ekm(i,j,k)*fc_global(i,j+rank*jmax,k)
	    ekm_plus=ekm(i,j,kp)*fc_global(i,j+rank*jmax,kp)
	    aaa(k)=-0.5*ekm_min*c2*dt/dz**2/(0.5*(rnew(i,j,km)+rnew(i+1,j,k)))
	    bbb(k)=1.+0.5*(ekm_min+ekm_plus)*c2*dt/dz**2/(0.5*(rnew(i,j,k)+rnew(i+1,j,kp)))
	    ccc(k)=-0.5*ekm_plus*c2*dt/dz**2/(0.5*(rnew(i,j,kp)+rnew(i+1,j,kpp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dWdt(i,j,0:k1)
	    CALL solve_tridiag(dWdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
       enddo
      ENDIF
	if (nfrac>0) then
	  call state(dcdt,drdt) ! determine drdt1 with pred1
	endif

c********************************************************************
c     CALCULATE CN pressure-gradient on predictor 1
c********************************************************************
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -c2*dt * ( pold(i+1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
	  if (periodicx.eq.1) then
	  !!! !periodic bc in x: 
	    i = imax
	      do k=1,kmax
		    do j=1,jmax
		      dUdt(i,j,k) = dUdt(i,j,k) -c2*dt * ( pold(1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		    enddo
	      enddo
	  endif  
      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -c2*dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * (phip(j+1)-phip(j)) )   
          enddo
        enddo
      enddo

	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c2*dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c2*dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
      	  endif
	endif

      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          dWdt(i,j,k) = dWdt(i,j,k) -c2*dt * ( pold(i,j,k+1) - pold(i,j,k) ) / dz
          enddo
        enddo
      enddo

!      DO n=1,npresIBM+1	!default npresIBM=0, then only one pres-cor, with npresIBM>0 extra pres-cor for better impermeable IBM boundaries
      call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),0,time_n+c2*dt,(1.-c2)*Ub1old+c2*Ub1new,(1.-c2)*Vb1old+c2*Vb1new,
     & (1.-c2)*Wb1old+c2*Wb1new,(1.-c2)*Ub2old+c2*Ub2new,(1.-c2)*Vb2old+c2*Vb2new,(1.-c2)*Wb2old+c2*Wb2new,
     & (1.-c2)*Ub3old+c2*Ub3new,(1.-c2)*Vb3old+c2*Vb3new,(1.-c2)*Wb3old+c2*Wb3new)

      call fillps2(dUdt,dVdt,dWdt,dRdt,time_n+c2*dt,c2*dt) 
!	  IF (poissolver.eq.1) THEN
!	   CALL SOLVEpois_vg(p)
!	  IF (poissolver.eq.2) THEN
!	   CALL SOLVEpois_vg_mumps(p)
	  IF (poissolver.eq.3) THEN
	   CALL SOLVEpois_vg_pardiso(p)	   
	  ELSE
	   CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
	  ENDIF

      call pshiftb(p,pplus2) !,rank,imax,jmax,kmax,px)
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -c2*dt * ( p(i+1,j,k) - p(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
	  !!! do nothing on imax and i=1, so p(imax+1)=p(imax) --> dp/dn=0 at outflow and inflow 
	  if (periodicx.eq.1) then
	  !!! !periodic bc in x: 
	    i = imax
	      do k=1,kmax
		    do j=1,jmax
		      dUdt(i,j,k) = dUdt(i,j,k) -c2*dt * ( p(1,j,k) - p(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		    enddo
	      enddo
	  endif	
      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -c2*dt * ( p(i,j+1,k) - p(i,j,k) ) /( Rp(i) * (phip(j+1)-phip(j)) )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c2*dt * ( pplus2(i,k) - p(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c2*dt * ( pplus2(i,k) - p(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
      	  endif
	endif

      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          dWdt(i,j,k) = dWdt(i,j,k) -c2*dt * ( p(i,j,k+1) - p(i,j,k) ) / dz
          enddo
        enddo
      enddo
!      ENDDO



c********************************************************************
c     split off the density and apply boundary conditions on predictor 1
c********************************************************************	
!	  dUdt1=dudt
!	  dVdt1=dvdt
!	  dWdt1=dWdt
      call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),0,time_n+c2*dt,(1.-c2)*Ub1old+c2*Ub1new,(1.-c2)*Vb1old+c2*Vb1new,
     & (1.-c2)*Wb1old+c2*Wb1new,(1.-c2)*Ub2old+c2*Ub2new,(1.-c2)*Vb2old+c2*Vb2new,(1.-c2)*Wb2old+c2*Wb2new,
     & (1.-c2)*Ub3old+c2*Ub3new,(1.-c2)*Vb3old+c2*Vb3new,(1.-c2)*Wb3old+c2*Wb3new)

	 dUdt=dUdt/rhU 
	 dVdt=dVdt/rhV
	 dWdt=dWdt/rhW
!      do k=1,kmax
!        do j=1,jmax
!          do i=1,imax
!           dUdt(i,j,k)=dUdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i+1,j,k)))
!           dVdt(i,j,k)=dVdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i,j+1,k)))
!           dWdt(i,j,k)=dWdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i,j,k+1)))
!          enddo
!        enddo
!      enddo   

      call bound_incljet(dUdt,dVdt,dWdt,drdt,0,time_n+c2*dt,(1.-c2)*Ub1old+c2*Ub1new,(1.-c2)*Vb1old+c2*Vb1new,
     & (1.-c2)*Wb1old+c2*Wb1new,(1.-c2)*Ub2old+c2*Ub2new,(1.-c2)*Vb2old+c2*Vb2new,(1.-c2)*Wb2old+c2*Wb2new,
     & (1.-c2)*Ub3old+c2*Ub3new,(1.-c2)*Vb3old+c2*Vb3new,(1.-c2)*Wb3old+c2*Wb3new)
!      call bound_incljet(dUdt1,dVdt1,dWdt1,drdt1,slip_bot,time_n+c2*dt,(1.-c2)*Ub1old+c2*Ub1new,(1.-c2)*Vb1old+c2*Vb1new,
!     & (1.-c2)*Wb1old+c2*Wb1new,(1.-c2)*Ub2old+c2*Ub2new,(1.-c2)*Vb2old+c2*Vb2new,(1.-c2)*Wb2old+c2*Wb2new)

c********************************************************************
c     CALCULATE k2 and predictor 2
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(dcdt,dWdt,wsed,drdt,1,kmax,sumWkm,cn2*dt,dz)
	       !! Set boundary conditions jet in:
	       do t=1,tmax_inPpunt
	 	i=i_inPpunt(t)
	 	j=j_inPpunt(t)
		do n=1,nfrac
	    	  do k=kmax,kmax !-kjet,kmax
		 	Wsed(n,i,j,k)=dWdt(i,j,k)
		  enddo
		enddo
		if (LOA>0.and.outflow_overflow_down.eq.1) then
		 do n=1,nfrac
	    	  do k=kmax-kjet+1,kmax 
		 	Wsed(n,i,j,k)=dWdt(i,j,k)
		  enddo
		 enddo	
		endif	
	       enddo
	  else
	    do n=1,nfrac
		wsed(n,:,:,:)=dWdt
	    enddo
	  endif

	  do n=1,nfrac
	      utr=dUdt
	      vtr=dVdt
	      wtr=Wsed(n,:,:,:)
		if (settling_along_gvector.eq.1) then
			gvector=sqrt(gx**2+gy**2+gz**2)
			ws = Wsed(n,:,:,:)-Wnew 
			utr=utr+gx/gvector*ws 
			vtr=vtr+gy/gvector*ws
			wtr=wtr-ws+gz/gvector*ws
		endif 		  
	      call make_UtransC_zeroIBM(utr,vtr,wtr)
		  if (advec_conc.eq.'NVD') then		  
	      call advecc_NVD(k2c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn2*dt,rank,px,periodicx,periodicy)
          elseif (advec_conc.eq.'VLE') then	 
	      call advecc_VLE(k2c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn2*dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)
          elseif (advec_conc.eq.'SBE') then	 
	      call advecc_SBE(k2c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn2*dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)	 
          elseif (advec_conc.eq.'VL2') then	 
	      call advecc_VL2_nocfl(k2c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn2*dt,rank,px,periodicx,periodicy)	 
          elseif (advec_conc.eq.'SB2') then	 
	      call advecc_SB2_nocfl(k2c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn2*dt,rank,px,periodicx,periodicy)		 
          elseif (advec_conc.eq.'ARO') then	 
		  call advecc_ARO(k2c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn2*dt,rank,px,periodicx,periodicy)
		  endif
		 if (Sc<1.e18) then
	      call diffc_CDS2 (k2c(n,:,:,:),dCdt(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)
	     endif
	  enddo
	  dcdt2=dcdt
	  do n=1,nfrac
		dcdt(n,:,:,:) =dcdt(n,:,:,:) + cn2*dt*k2c(n,:,:,:) !pred2
!		dcdt(n,:,:,:) =cnew(n,:,:,:) + a31*dt*k1c(n,:,:,:) + a32*dt*k2c(n,:,:,:) !pred2
!		call bound_c(dcdt(n,:,:,:),frac(n)%c,n,0.) ! switched off 6-7-16 RK3 k2 because bound_c after dep-ero
           IF (CNdiffz.eq.1) THEN !CN semi implicit treatment diff-z:
			call bound_c(dcdt(n,:,:,:),frac(n)%c,n,0.) ! bc start CN-diffz k2 RK3 
	     do k=k1,k1 !-kjet,k1
	       do t=1,tmax_inPpunt
	         i=i_inPpunt(t)
	         j=j_inPpunt(t)
	         dcdt(n,i,j,k) = dcdt(n,i,j,kmax) ! No diffusion over vertical inflow boundary, this line is needed for exact influx
	       enddo
	     enddo
             do j=1,jmax
               do i=1,imax
                 do k=0,k1
	           km=MAX(0,k-1)
	           kp=MIN(k1,k+1)
	           kpp=MIN(k1,k+2)
	           ekm_min=0.5*(Diffcof(i,j,km)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km))
	           ekm_plus=0.5*(Diffcof(i,j,kp)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp))
	           aaa(k)=-0.5*ekm_min*cn2*dt/dz**2
	           bbb(k)=1.+0.5*(ekm_min+ekm_plus)*cn2*dt/dz**2
	           ccc(k)=-0.5*ekm_plus*cn2*dt/dz**2
                 enddo
		    aaa(0)=0.
		    aaa(k1)=0.
		    ccc(0)=0.
		    ccc(k1)=0.
		    bbb(0)=1.
		    bbb(k1)=1.
	         rhss=dcdt(n,i,j,0:k1)
	         CALL solve_tridiag(dcdt(n,i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
               enddo
             enddo
	     !call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc after CN-diffz switched off 6-7-16 
	   ENDIF
	  enddo
	  call slipvelocity_bed(dcdt2,dwdt,wsed,drdt,sumWkm,cn2*dt,dz) !get sumWkm and wsed at bed for driftflux_force 
	  if (interaction_bed>0) then
	    CALL erosion_deposition(dcdt,dcdt2bot,dudt,dvdt,dwdt,drdt,dcdt2,dcdt1bot,cn2*dt,dz) !first two vars are adjusted 
		do n=1,nfrac 
		  call bound_cbot(dcdt2bot(n,:,:))
		enddo			
		!! advec concentration in bed with velocity TSHD:
		IF (ABS(U_TSHD)>1e-6) THEN 
		  do n=1,nfrac 
	            call adveccbot_TVD(k2cbot(n,:,:),dcdt2bot(n,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,phiv,phipt,dz,
     +              i1,j1,ib,ie,jb,je,cn2*dt,rank,px,periodicx,periodicy)
		    dcdt2bot(n,:,:)= dcdt2bot(n,:,:) + cn2*dt*k2cbot(n,:,:) ! pred2	
			call bound_cbot(dcdt2bot(n,:,:))		
		  enddo 
		ENDIF 		
 	    !! for kn2: variable dcdt2 made
		if (interaction_bed.eq.4.or.interaction_bed.eq.6) then
			call advec_update_Clivebed(dcdt,dcdt2bot,cn2*dt) 
		endif			
	    do n=1,nfrac
		call bound_cbot(dcdt2bot(n,:,:))
	    enddo
	  endif
	    do n=1,nfrac 
			call bound_c(dcdt(n,:,:,:),frac(n)%c,n,cn2*dt) ! bc after erosion_deposition k2 RK3
	    enddo	  
	endif

	if (diffusion.eq.'COM4') then
		call stress_terms(dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	endif

      	if (convection.eq.'CDS2') then
      call advecu_CDS2(k2u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(k2u,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(k2u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecu_CDS4(k2u,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(k2u,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(k2u,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2)	 
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(k2u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd,numdiff2)	  
	elseif(convection.eq.'uTVD') then
      call advecu_TVD(k2u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif
	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (k2u,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(k2u,Srr,Spr,Szr,Spp,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      	if (convection.eq.'CDS2') then
      call advecv_CDS2(k2v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(k2v,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(k2v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecv_CDS4(k2v,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(k2v,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(k2v,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(k2v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecv_TVD(k2v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif
	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (k2v,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(k2v,Spr,Spp,Spz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      	if (convection.eq.'CDS2') then
      call advecw_CDS2(k2w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(k2w,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(k2w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecw_CDS4(k2w,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(k2w,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(k2w,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(k2w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecw_TVD(k2w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif
	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (k2w,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(k2w,Szr,Spz,Szz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 call advecw_driftfluxCDS2(k2w,driftfluxforce_calfac*sumWkm,
     &   drdt,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
	endif
	
	

        do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    k2u(i,j,k) = k2u(i,j,k)-(gx*cos_u(j)+gy*sin_u(j))*(0.5*(drdt(i,j,k)+drdt(i+1,j,k))-rho_b)
            dUdt(i,j,k)=Unew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i+1,j,k))+(a31*k1u(i,j,k)+a32*k2u(i,j,k)	! pred2
     1     +(Ppropx(i,j,k))*c3)*dt 
            enddo
         enddo
        enddo
        do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    k2v(i,j,k) = k2v(i,j,k)-(-gx*sin_v(j)+gy*cos_v(j))*(0.5*(drdt(i,j,k)+drdt(i,j+1,k))-rho_b)
            dVdt(i,j,k)=Vnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j+1,k))+(a31*k1v(i,j,k)+a32*k2v(i,j,k) ! pred2
     1     +(Ppropy(i,j,k))*c3)*dt 
            enddo
         enddo
        enddo
        do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    k2w(i,j,k)=k2w(i,j,k)-gz*(0.5*(drdt(i,j,k)+drdt(i,j,k+1))-rho_b)
            dWdt(i,j,k)=Wnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j,k+1))+(a31*k1w(i,j,k)+a32*k2w(i,j,k) ! pred2
     1     +(Ppropz(i,j,k))*c3)*dt 
            enddo
         enddo
        enddo

      IF (CNdiffz.eq.1) THEN !CN semi implicit treatment diff-z:
	if (nfrac>0) then
	  call state(dcdt,drdt) ! determine drdt1 with pred1
	endif
      call bound_rhoU(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new,
     & (1.-c3)*Ub3old+c3*Ub3new,(1.-c3)*Vb3old+c3*Vb3new,(1.-c3)*Wb3old+c2*Wb3new)
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))		
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))			
	    aaa(k)=-0.5*ekm_min*c3*dt/dz**2/(0.5*(drdt(i,j,km)+drdt(i+1,j,km)))
	    bbb(k)=1.+0.5*(ekm_min+ekm_plus)*c3*dt/dz**2/(0.5*(drdt(i,j,k)+drdt(i+1,j,k)))
	    ccc(k)=-0.5*ekm_plus*c3*dt/dz**2/(0.5*(drdt(i,j,kp)+drdt(i+1,j,kp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dUdt(i,j,0:k1)
	    CALL solve_tridiag(dUdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
      enddo
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))			
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))			
	    aaa(k)=-0.5*ekm_min*c3*dt/dz**2/(0.5*(drdt(i,j,km)+drdt(i,j+1,km)))
	    bbb(k)=1.+0.5*(ekm_min+ekm_plus)*c3*dt/dz**2/(0.5*(drdt(i,j,k)+drdt(i,j+1,k)))
	    ccc(k)=-0.5*ekm_plus*c3*dt/dz**2/(0.5*(drdt(i,j,kp)+drdt(i,j+1,kp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dVdt(i,j,0:k1)
	    CALL solve_tridiag(dVdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
      enddo

       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=ekm(i,j,k)*fc_global(i,j+rank*jmax,k)
	    ekm_plus=ekm(i,j,kp)*fc_global(i,j+rank*jmax,kp)
	    aaa(k)=-0.5*ekm_min*c3*dt/dz**2/(0.5*(drdt(i,j,km)+drdt(i+1,j,k)))
	    bbb(k)=1.+0.5*(ekm_min+ekm_plus)*c3*dt/dz**2/(0.5*(drdt(i,j,k)+drdt(i+1,j,kp)))
	    ccc(k)=-0.5*ekm_plus*c3*dt/dz**2/(0.5*(drdt(i,j,kp)+drdt(i+1,j,kpp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dWdt(i,j,0:k1)
	    CALL solve_tridiag(dWdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
       enddo
      ENDIF
	if (nfrac>0) then
	  call state(dcdt,drdt) ! determine drdt1 with pred1
	endif

c********************************************************************
c     CALCULATE CN pressure-gradient on predictor 2
c********************************************************************
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -c3*dt * ( pold(i+1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
	  if (periodicx.eq.1) then
	  !!! !periodic bc in x: 
	    i = imax
	      do k=1,kmax
		    do j=1,jmax
		      dUdt(i,j,k) = dUdt(i,j,k) -c3*dt * ( pold(1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		    enddo
	      enddo
	  endif
      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -c3*dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * (phip(j+1)-phip(j)) )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c3*dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c3*dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
      	  endif
	endif
      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          dWdt(i,j,k) = dWdt(i,j,k) -c3*dt * ( pold(i,j,k+1) - pold(i,j,k) ) / dz
          enddo
        enddo
      enddo  

!      DO n=1,npresIBM+1	!default npresIBM=0, then only one pres-cor, with npresIBM>0 extra pres-cor for better impermeable IBM boundaries
      call bound_rhoU(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new,
     & (1.-c3)*Ub3old+c3*Ub3new,(1.-c3)*Vb3old+c3*Vb3new,(1.-c3)*Wb3old+c2*Wb3new) 
!      call bound_rhoU(dUdt2,dVdt2,dWdt2,dRdt2,slip_bot,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
!     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new)
      call fillps2(dUdt,dVdt,dWdt,dRdt,time_n+c3*dt,c3*dt) 
!	  IF (poissolver.eq.1) THEN
!	   CALL SOLVEpois_vg(p)
!	  IF (poissolver.eq.2) THEN
!	   CALL SOLVEpois_vg_mumps(p)
	  IF (poissolver.eq.3) THEN
	   CALL SOLVEpois_vg_pardiso(p)	   
	  ELSE
	   CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
	  ENDIF
      call pshiftb(p,pplus2) !,rank,imax,jmax,kmax,px)
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -c3*dt * ( p(i+1,j,k) - p(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
	  !!! do nothing on imax and i=1, so p(imax+1)=p(imax) --> dp/dn=0 at outflow and inflow 
	  if (periodicx.eq.1) then
	  !!! !periodic bc in x: 
	    i = imax
	      do k=1,kmax
		    do j=1,jmax
		      dUdt(i,j,k) = dUdt(i,j,k) -c3*dt * ( p(1,j,k) - p(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		    enddo
	      enddo
	  endif 
      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -c3*dt * ( p(i,j+1,k) - p(i,j,k) ) /( Rp(i) * (phip(j+1)-phip(j)) )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c3*dt * ( pplus2(i,k) - p(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c3*dt * ( pplus2(i,k) - p(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
      	  endif
	endif
      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          dWdt(i,j,k) = dWdt(i,j,k) -c3*dt * ( p(i,j,k+1) - p(i,j,k) ) / dz
          enddo
        enddo
      enddo   
!      ENDDO
c********************************************************************
c     split off the density and apply boundary conditions to predictor 2
c********************************************************************	
!		dUdt1=dudt
!		dVdt1=dvdt
!		dWdt1=dwdt
      call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),0,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new,
     & (1.-c3)*Ub3old+c3*Ub3new,(1.-c3)*Vb3old+c3*Vb3new,(1.-c3)*Wb3old+c2*Wb3new)

	 dUdt=dUdt/rhU 
	 dVdt=dVdt/rhV
	 dWdt=dWdt/rhW
!      do k=1,kmax
!        do j=1,jmax
!          do i=1,imax
!           dUdt(i,j,k)=dUdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i+1,j,k)))
!           dVdt(i,j,k)=dVdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i,j+1,k)))
!           dWdt(i,j,k)=dWdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i,j,k+1)))
!          enddo
!        enddo
!      enddo       

      call bound_incljet(dUdt,dVdt,dWdt,dRdt,0,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new,
     & (1.-c3)*Ub3old+c3*Ub3new,(1.-c3)*Vb3old+c3*Vb3new,(1.-c3)*Wb3old+c2*Wb3new)
!      call bound_incljet(dUdt2,dVdt2,dWdt2,dRdt2,slip_bot,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
!     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new)	

!      call LES_smagorinsky(dUdt2,dVdt2,dWdt2,drdt2)
!	Diffcof=ekm/Sc/dRdt2+ekm_mol*invrho_b 
!	ekm = ekm + ekm_mol
c********************************************************************
c     CALCULATE k3 and n+1 values with RK3
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(dcdt,dWdt,wsed,drdt,1,kmax,sumWkm,cn3*dt,dz)
	       !! Set boundary conditions jet in:
	       do t=1,tmax_inPpunt
	 	i=i_inPpunt(t)
	 	j=j_inPpunt(t)
		do n=1,nfrac
	    	  do k=kmax,kmax !-kjet,kmax
		 	Wsed(n,i,j,k)=dWdt(i,j,k)
		  enddo		
		enddo
		if (LOA>0.and.outflow_overflow_down.eq.1) then
		 do n=1,nfrac
	    	  do k=kmax-kjet+1,kmax 
		 	Wsed(n,i,j,k)=dWdt(i,j,k)
		  enddo
		 enddo	
		endif	
	       enddo
	  else
	    do n=1,nfrac
		wsed(n,:,:,:)=dWdt
	    enddo
	  endif
	  	
	  do n=1,nfrac
	      utr=dUdt
	      vtr=dVdt
	      wtr=Wsed(n,:,:,:)
		if (settling_along_gvector.eq.1) then
			gvector=sqrt(gx**2+gy**2+gz**2)
			ws = Wsed(n,:,:,:)-Wnew 
			utr=utr+gx/gvector*ws 
			vtr=vtr+gy/gvector*ws
			wtr=wtr-ws+gz/gvector*ws
		endif 		  
	      call make_UtransC_zeroIBM(utr,vtr,wtr)
		  if (advec_conc.eq.'NVD') then		  
	      call advecc_NVD(k3c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn3*dt,rank,px,periodicx,periodicy)
          elseif (advec_conc.eq.'VLE') then	 	 
	      call advecc_VLE(k3c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn3*dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)
          elseif (advec_conc.eq.'SBE') then	 	 
	      call advecc_SBE(k3c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn3*dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)	 
          elseif (advec_conc.eq.'VL2') then	 	 
	      call advecc_VL2_nocfl(k3c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn3*dt,rank,px,periodicx,periodicy)
	      elseif (advec_conc.eq.'SB2') then	 	 
	      call advecc_SB2_nocfl(k3c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn3*dt,rank,px,periodicx,periodicy)
		  elseif (advec_conc.eq.'ARO') then	
	      call advecc_ARO(k3c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn3*dt,rank,px,periodicx,periodicy)
		  endif
		  if (Sc<1.e18) then
	      call diffc_CDS2 (k3c(n,:,:,:),dCdt(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)
		  endif
	  enddo
	  dcdt2=dcdt 
	  do n=1,nfrac
!		k3c(n,:,:,:) = dnewc(n,:,:,:)
!		dcdt(n,:,:,:) =cnew(n,:,:,:) + b1*dt*k1c(n,:,:,:) + b2*dt*k2c(n,:,:,:)+ b3*dt*k3c(n,:,:,:) !n+1
		dcdt(n,:,:,:) =dcdt(n,:,:,:) + cn3*dt*k3c(n,:,:,:) !n+1
!		CALL air_bubbles_free_surface ! added air_bubbles_free_surface here, just before bc are applied
!		9-10:2012: air_bubbles_free_surface not called because wsed is calculated on top row now, so air is dissappearing anyway
!		call bound_c(dcdt(n,:,:,:),frac(n)%c,n,0.) ! switched off 6-7-16 RK3 k3 because bound_c after dep-ero
           IF (CNdiffz.eq.1) THEN !CN semi implicit treatment diff-z:
		   call bound_c(dcdt(n,:,:,:),frac(n)%c,n,0.) ! start CN-diffz k3 RK3 
	     do k=k1,k1 !-kjet,k1
	       do t=1,tmax_inPpunt
	         i=i_inPpunt(t)
	         j=j_inPpunt(t)
	         dcdt(n,i,j,k) = dcdt(n,i,j,kmax) ! No diffusion over vertical inflow boundary, this line is needed for exact influx
	       enddo
	     enddo
             do j=1,jmax
               do i=1,imax
                 do k=0,k1
	           km=MAX(0,k-1)
	           kp=MIN(k1,k+1)
	           kpp=MIN(k1,k+2)
	           ekm_min=0.5*(Diffcof(i,j,km)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km))
	           ekm_plus=0.5*(Diffcof(i,j,kp)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp))
	           aaa(k)=-0.5*ekm_min*cn3*dt/dz**2
	           bbb(k)=1.+0.5*(ekm_min+ekm_plus)*cn3*dt/dz**2
	           ccc(k)=-0.5*ekm_plus*cn3*dt/dz**2
                 enddo
		    aaa(0)=0.
		    aaa(k1)=0.
		    ccc(0)=0.
		    ccc(k1)=0.
		    bbb(0)=1.
		    bbb(k1)=1.
	         rhss=dcdt(n,i,j,0:k1)
	         CALL solve_tridiag(dcdt(n,i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
               enddo
             enddo
	      
	   ENDIF
	  enddo
	    call slipvelocity_bed(dcdt2,dwdt,wsed,drdt,sumWkm,cn3*dt,dz) !get sumWkm and wsed at bed for driftflux_force 
		if (interaction_bed>0) then
	    CALL erosion_deposition(dcdt,dcdtbot,dudt,dvdt,dwdt,drdt,dcdt2,dcdt2bot,cn3*dt,dz) !first two vars are adjusted
		do n=1,nfrac 
		  call bound_cbot(dcdtbot(n,:,:))
		enddo			
		!! advec concentration in bed with velocity TSHD:
		IF (ABS(U_TSHD)>1e-6) THEN 
		  do n=1,nfrac 
	            call adveccbot_TVD(k3cbot(n,:,:),dcdtbot(n,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,phiv,phipt,dz,
     +              i1,j1,ib,ie,jb,je,cn3*dt,rank,px,periodicx,periodicy)
		    dcdtbot(n,:,:)= dcdtbot(n,:,:) + cn3*dt*k3cbot(n,:,:) ! n+1	
			call bound_cbot(dcdtbot(n,:,:))			
		  enddo 
		ENDIF 		
		    !! for kn3: extra var dcdt2 is used
	    ! should correct dwdt(n,:,:,0) with (1.-tau/frac(n)%tau_d), but this correction almost always is zero, for now leave
		if (interaction_bed.eq.4.or.interaction_bed.eq.6) then
			call advec_update_Clivebed(dcdt,dcdtbot,cn3*dt) 
		endif		
	    do n=1,nfrac
		call bound_cbot(dcdtbot(n,:,:))
	    enddo
	  endif
	    do n=1,nfrac 
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n,cn3*dt) ! bc after erosion_deposition k3 RK3
	    enddo	  
	endif

	if (diffusion.eq.'COM4') then
		call stress_terms(dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	endif

      	if (convection.eq.'CDS2') then
      call advecu_CDS2(k3u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(k3u,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(k3u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecu_CDS4(k3u,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(k3u,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(k3u,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2)	 
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(k3u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd,numdiff2)	  
	elseif(convection.eq.'uTVD') then
      call advecu_TVD(k3u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif

	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (k3u,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(k3u,Srr,Spr,Szr,Spp,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      	if (convection.eq.'CDS2') then
      call advecv_CDS2(k3v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(k3v,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(k3v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecv_CDS4(k3v,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(k3v,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(k3v,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(k3v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecv_TVD(k3v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif
	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (k3v,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(k3v,Spr,Spp,Spz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      	if (convection.eq.'CDS2') then
      call advecw_CDS2(k3w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(k3w,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(k3w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecw_CDS4(k3w,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(k3w,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(k3w,dUdt,dVdt,dWdt,dRdt,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(k3w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecw_TVD(k3w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif
	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (k3w,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(k3w,Szr,Spz,Szz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 call advecw_driftfluxCDS2(k3w,driftfluxforce_calfac*sumWkm,
     &   drdt,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
	endif
	
	
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    k3u(i,j,k) = k3u(i,j,k)-(gx*cos_u(j)+gy*sin_u(j))*(0.5*(drdt(i,j,k)+drdt(i+1,j,k))-rho_b)
            dUdt(i,j,k)=Unew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i+1,j,k))+(b1*k1u(i,j,k)+b2*k2u(i,j,k)+b3*k3u(i,j,k)	!n+1
     1     +(Ppropx(i,j,k)))*dt 
            enddo
         enddo
      enddo
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    k3v(i,j,k) = k3v(i,j,k)-(-gx*sin_v(j)+gy*cos_v(j))*(0.5*(drdt(i,j,k)+drdt(i,j+1,k))-rho_b)
            dVdt(i,j,k)=Vnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j+1,k))+(b1*k1v(i,j,k)+b2*k2v(i,j,k)+b3*k3v(i,j,k) !n+1
     1     +(Ppropy(i,j,k)))*dt 
            enddo
         enddo
      enddo
       do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    k3w(i,j,k)=k3w(i,j,k)-gz*(0.5*(drdt(i,j,k)+drdt(i,j,k+1))-rho_b)
            dWdt(i,j,k)=Wnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j,k+1))+(b1*k1w(i,j,k)+b2*k2w(i,j,k)+b3*k3w(i,j,k) !n+1
     1     +(Ppropz(i,j,k)))*dt 
            enddo
         enddo
       enddo

      IF (CNdiffz.eq.1) THEN !CN semi implicit treatment diff-z:
	if (nfrac>0) then
	  call state(dcdt,drdt) ! determine drdt1 with pred1
	endif
      call bound_rhoU(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
     & Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))			
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))		
	    aaa(k)=-0.5*ekm_min*dt/dz**2/(0.5*(drdt(i,j,km)+drdt(i+1,j,km)))
	    bbb(k)=1.+0.5*(ekm_min+ekm_plus)*dt/dz**2/(0.5*(drdt(i,j,k)+drdt(i+1,j,k)))
	    ccc(k)=-0.5*ekm_plus*dt/dz**2/(0.5*(drdt(i,j,kp)+drdt(i+1,j,kp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dUdt(i,j,0:k1)
	    CALL solve_tridiag(dUdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
      enddo
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))			
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))			
	    aaa(k)=-0.5*ekm_min*dt/dz**2/(0.5*(drdt(i,j,km)+drdt(i,j+1,km)))
	    bbb(k)=1.+0.5*(ekm_min+ekm_plus)*dt/dz**2/(0.5*(drdt(i,j,k)+drdt(i,j+1,k)))
	    ccc(k)=-0.5*ekm_plus*dt/dz**2/(0.5*(drdt(i,j,kp)+drdt(i,j+1,kp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dVdt(i,j,0:k1)
	    CALL solve_tridiag(dVdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
      enddo

       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=ekm(i,j,k)*fc_global(i,j+rank*jmax,k)
	    ekm_plus=ekm(i,j,kp)*fc_global(i,j+rank*jmax,kp)
	    aaa(k)=-0.5*ekm_min*dt/dz**2/(0.5*(drdt(i,j,km)+drdt(i+1,j,k)))
	    bbb(k)=1.+0.5*(ekm_min+ekm_plus)*dt/dz**2/(0.5*(drdt(i,j,k)+drdt(i+1,j,kp)))
	    ccc(k)=-0.5*ekm_plus*dt/dz**2/(0.5*(drdt(i,j,kp)+drdt(i+1,j,kpp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dWdt(i,j,0:k1)
	    CALL solve_tridiag(dWdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
       enddo
       ENDIF
	if (nfrac>0) then
	  call state(dcdt,drdt) ! determine drdt1 with pred1
	endif

c********************************************************************
c     ready now, split off the density and apply boundary conditions 
c     to n+1 values takes place in main.f after correction with pressure-gradient
c********************************************************************	

!      call fillps2(dUdt,dVdt,dWdt,dRdt,time_n+dt,dt) 
!      CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
!      call pshiftb(p,pplus,rank,imax,jmax,kmax,px)


c********************************************************************
c     CALCULATE CN pressure-gradient on predictor 3
c********************************************************************
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -dt * ( pold(i+1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo
      enddo
	  if (periodicx.eq.1) then
	  !!! !periodic bc in x: 
	    i = imax
	      do k=1,kmax
		    do j=1,jmax
		      dUdt(i,j,k) = dUdt(i,j,k) -dt * ( pold(1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		    enddo
	      enddo
	  endif
      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * (phip(j+1)-phip(j)) )
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
      	  endif
	endif
      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          dWdt(i,j,k) = dWdt(i,j,k) -dt * ( pold(i,j,k+1) - pold(i,j,k) ) / dz
          enddo
        enddo
      enddo

      return
      end


      subroutine euler_expl(ib,ie,jb,je,kb,ke)
c
c    Performes time integration with first order 
c    Euler-explicit scheme, i.e
c
c
c             n+1     n
c     dU     U     - U                            n
c    ---- = ------------ = ( -ADV + DIFF + Force)     
c     dt        dt                                   

c
c   The timestep is limited with CFL=1 (see routine chkdt)
c
c 
      USE nlist
      USE sediment

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer  ib,ie,jb,je,kb,ke,n,t,kpp,kp,km,me
      real     doldc(nfrac,0:i1,0:j1,0:k1),dnewc(nfrac,0:i1,0:j1,0:k1)
      real     dold(0:i1,0:j1,0:k1),dnew(0:i1,0:j1,0:k1)
	  real     wsed(nfrac,0:i1,0:j1,0:k1),sumWkm(0:i1,0:j1,0:k1)
	real*8 pplus(imax,kmax)
	real dnewcbot(nfrac,0:i1,0:j1)	
	real aaa(0:k1),bbb(0:k1),ccc(0:k1),ekm_min,ekm_plus,rhss(0:k1)
	real utr(0:i1,0:j1,0:k1),vtr(0:i1,0:j1,0:k1),wtr(0:i1,0:j1,0:k1)	
!		real dUdt1(0:i1,0:j1,0:k1),dVdt1(0:i1,0:j1,0:k1),dWdt1(0:i1,0:j1,0:k1)
	real ws(0:i1,0:j1,0:k1),gvector,CNz


	me = momentum_exchange_obstacles !100 or 101 or 111 means no advection momentum exchange with bed-cells with HYB6


	dcdt = 0.
	sumWkm = 0.
	dnewc=0.
c********************************************************************
c     CALCULATE slipvelocity
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(cnew,Wnew,wsed,rnew,1,kmax,sumWkm,dt,dz)
	       !! Set boundary conditions jet in:
	       do t=1,tmax_inPpunt
	 	i=i_inPpunt(t)
	 	j=j_inPpunt(t)
		do n=1,nfrac
	    	  do k=kmax,kmax !-kjet,kmax
		 	Wsed(n,i,j,k)=Wnew(i,j,k)
		  enddo
		enddo
		if (LOA>0.and.outflow_overflow_down.eq.1) then
		 do n=1,nfrac
	    	  do k=kmax-kjet+1,kmax 
			Wsed(n,i,j,k)=Wnew(i,j,k)
		  enddo
		 enddo	
		endif
	       enddo
	  else
	    do n=1,nfrac
		wsed(n,:,:,:)=wnew
	    enddo
	  endif

c********************************************************************
c     CALCULATE advection, diffusion Concentration
c********************************************************************
!	  do n=1,nfrac
!	      call advecc_TVD2(dnewc(n,:,:,:),cnew(n,:,:,:),Unew,Vnew,Wsed(n,:,:,:),rnew,Ru,Rp,dr,phiv,phipt,dz,
!     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
!	      call diffc_CDS2 (dnewc(n,:,:,:),Cnew(n,:,:,:),Diffcof,
!     +            ib,ie,jb,je,kb,ke)

!!		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(1.5*dnewc(n,:,:,:)-0.5*cc(n,:,:,:)) !AB2
!		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*dnewc(n,:,:,:)                       !EE1
!		call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc for intermediate dcdt
!		cc(n,:,:,:)=dnewc(n,:,:,:)
!	  enddo

	  do n=1,nfrac
	      utr=Unew
	      vtr=Vnew
	      wtr=Wsed(n,:,:,:)
		if (settling_along_gvector.eq.1) then
			gvector=sqrt(gx**2+gy**2+gz**2)
			ws = Wsed(n,:,:,:)-Wnew 
			utr=utr+gx/gvector*ws 
			vtr=vtr+gy/gvector*ws
			wtr=wtr-ws+gz/gvector*ws
		endif 		  
	      call make_UtransC_zeroIBM(utr,vtr,wtr)
		if (frac(n)%type.eq.-1) then ! VOF fraction
	         call advecc_VOF(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
!	      call advecc_NVD_VOF(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
!     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)	 
			 dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(dnewc(n,:,:,:)) !update in time with EE1 for TVD, no diffusion
			 do i=1,imax
			   do j=1,jmax
			     do k=1,kmax
					dcdt(n,i,j,k)=MAX(0.,dcdt(n,i,j,k))
					dcdt(n,i,j,k)=MIN(1.,dcdt(n,i,j,k))
				 enddo
			   enddo
			 enddo 
		else		  
		  if (advec_conc.eq.'NVD') then
	      call advecc_NVD(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)	 
          elseif (advec_conc.eq.'VLE') then	 
	      call advecc_VLE(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)
          elseif (advec_conc.eq.'SBE') then	 
	      call advecc_SBE(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)	 
          elseif (advec_conc.eq.'VL2') then	 
	      call advecc_VL2_nocfl(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
          elseif (advec_conc.eq.'SB2') then	 
	      call advecc_SB2_nocfl(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)	 
		  elseif (advec_conc.eq.'ARO') then	 
	      call advecc_ARO(dnewc(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)	 
          endif
		
		if (Sc<1.e18) then
	      call diffc_CDS2 (dnewc(n,:,:,:),Cnew(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)
		endif

		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(dnewc(n,:,:,:)) !update in time with EE1 for TVD
	   endif 
           IF ((CNdiffz.eq.1.or.CNdiffz.eq.2).and.Sc<1.e18) THEN !CN semi implicit treatment diff-z:
			IF (CNdiffz.eq.1) THEN 
				CNz=0.5 
			ELSE 
				CNz=1.
			ENDIF 		   
	     call bound_c(dcdt(n,:,:,:),frac(n)%c,n,0.) ! bc start CN-diffz EE1
	     do k=k1,k1 !-kjet,k1
	       do t=1,tmax_inPpunt
	         i=i_inPpunt(t)
	         j=j_inPpunt(t)
	         dcdt(n,i,j,k) = dcdt(n,i,j,kmax) ! No diffusion over vertical inflow boundary, this line is needed for exact influx
	       enddo
	     enddo
             do j=1,jmax
               do i=1,imax
                 do k=0,k1
	           km=MAX(0,k-1)
	           kp=MIN(k1,k+1)
	           kpp=MIN(k1,k+2)
	           ekm_min=0.5*(Diffcof(i,j,km)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km))
	           ekm_plus=0.5*(Diffcof(i,j,kp)+Diffcof(i,j,k))*MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp))
	           aaa(k)=-CNz*ekm_min*dt/dz**2
	           bbb(k)=1.+CNz*(ekm_min+ekm_plus)*dt/dz**2
	           ccc(k)=-CNz*ekm_plus*dt/dz**2
                 enddo
		 aaa(0)=0.
		 aaa(k1)=0.
		 ccc(0)=0.
		 ccc(k1)=0.
		 bbb(0)=1.
		 bbb(k1)=1.
	         rhss=dcdt(n,i,j,0:k1)
	         CALL solve_tridiag(dcdt(n,i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
               enddo
             enddo
	   ENDIF
          enddo
		  call slipvelocity_bed(cnew,wnew,wsed,rnew,sumWkm,dt,dz) !driftflux_force must be calculated with settling velocity at bed
      	  if (interaction_bed>0) then
	    CALL erosion_deposition(dcdt,dcdtbot,unew,vnew,wnew,rnew,cnew,cnewbot,dt,dz) !first two vars are adjusted
		do n=1,nfrac 
		  call bound_cbot(dcdtbot(n,:,:))
		enddo			
		!! advec concentration in bed with velocity TSHD:
		IF (ABS(U_TSHD)>1e-6) THEN 
		  do n=1,nfrac 
			call adveccbot_TVD(dnewcbot(n,:,:),dcdtbot(n,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,phiv,phipt,dz,
     +            	i1,j1,ib,ie,jb,je,dt,rank,px,periodicx,periodicy)
			dcdtbot(n,:,:)= dcdtbot(n,:,:) + dt*dnewcbot(n,:,:) ! time update	
			call bound_cbot(dcdtbot(n,:,:))			
		  enddo 
		ENDIF 		
	    !! erosion_deposition must be after advecc_TVD and after dcdt update, because cnew and dcdt are two separate vars
		if (interaction_bed.eq.4.or.interaction_bed.eq.6) then
			call advec_update_Clivebed(dcdt,dcdtbot,dt) 
		endif		
	    do n=1,nfrac 
	        call bound_cbot(dcdtbot(n,:,:))
	    enddo
      	  endif !end (interaction_bed>0) 
	    do n=1,nfrac 
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n,dt) ! bc after erosion_deposition EE1
	    enddo
	  call state(dcdt,drdt) ! determine drdt with intermediate dcdt
	   if (applyVOF.eq.1.or.applyVOF.eq.2) then ! VOF fraction, do not solve fully variable density NavierStokes but NS with constant density in advection and diffusion; only apply variable density in gravity term
		 drdt=rho_b
		 rhU=rho_b
		 rhV=rho_b
		 rhW=rho_b 
	   endif		  
	endif


	if (diffusion.eq.'COM4') then
		call stress_terms(Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	endif
c********************************************************************
c     CALCULATE advection, diffusion and Force U-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecu_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecu_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     & periodicx,periodicy,wf,wd,numdiff2)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,numdiff2)	 
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd,numdiff2)
	elseif(convection.eq.'uTVD') then
      call advecu_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)		  
	endif

	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(dnew,Srr,Spr,Szr,Spp,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif


      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-(gx*cos_u(j)+gy*sin_u(j))*(0.5*(rnew(i,j,k)+rnew(i+1,j,k))-rho_b)
     1     +Ppropx(i,j,k) 
            dUdt(i,j,k)=Unew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i+1,j,k))+
     1     dt*dnew(i,j,k)

            enddo
         enddo
      enddo

c********************************************************************
c     CALCULATE advection, diffusion and Force V-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecv_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecv_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,
     & rank,px,numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,
     & rank,px,numdiff,periodicx,periodicy,kbed,wf,wd)	 
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecv_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif

	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(dnew,Spr,Spp,Spz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-(-gx*sin_v(j)+gy*cos_v(j))*(0.5*(rnew(i,j,k)+rnew(i,j+1,k))-rho_b)
     1     +Ppropy(i,j,k) 
            dVdt(i,j,k)=Vnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j+1,k))+
     1     dt*dnew(i,j,k)
            enddo
         enddo
      enddo
c********************************************************************
c     CALCULATE advection, diffusion and Force W-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecw_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS4') then
      call advecw_CDS4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,
     &  periodicx,periodicy,wf,wd)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd,fc_global,me)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(dnew,Unew,Vnew,Wnew,Rnew,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,
     & numdiff,periodicx,periodicy,kbed,wf,wd)
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy
     &,phiv,wf,wd)	  
	elseif(convection.eq.'uTVD') then
      call advecw_TVD(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt
     & ,periodicx,periodicy,istep)	  
	endif

	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(dnew,Szr,Spz,Szz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 call advecw_driftfluxCDS2(dnew,driftfluxforce_calfac*sumWkm,
     &   rnew,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
	endif
	
	
	   if (applyVOF.eq.1) then ! VOF fraction apply real density only here for correct determination gravity
		 CALL state(cnew,rnew)
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-gz*(0.5*(rnew(i,j,k)+rnew(i,j,k+1))-rho_b)/(0.5*(rnew(i,j,k)+rnew(i,j,k+1)))
     1     +Ppropz(i,j,k)
            dWdt(i,j,k)=Wnew(i,j,k)*rhW(i,j,k)+ ! *0.5*(rnew(i,j,k)+rnew(i,j,k+1))+
     1     dt*dnew(i,j,k)
            enddo
         enddo
      enddo	
	    rnew=rho_b ! make density 1 again
	   else 
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-gz*(0.5*(rnew(i,j,k)+rnew(i,j,k+1))-rho_b)
     1     +Ppropz(i,j,k)
            dWdt(i,j,k)=Wnew(i,j,k)*rhW(i,j,k)+
     1     dt*dnew(i,j,k)
            enddo
         enddo
      enddo
	 endif 

      IF (CNdiffz.eq.1.or.CNdiffz.eq.2) THEN !CN semi implicit treatment diff-z:
			IF (CNdiffz.eq.1) THEN 
				CNz=0.5 !0.55 
			ELSE 
				CNz=1.
			ENDIF 	  
      call bound_rhoU(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
     & Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))		
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))		
	    aaa(k)=-CNz*ekm_min*dt/dz**2/(0.5*(drdt(i,j,km)+drdt(i+1,j,km)))
	    bbb(k)=1.+CNz*(ekm_min+ekm_plus)*dt/dz**2/(0.5*(drdt(i,j,k)+drdt(i+1,j,k)))
	    ccc(k)=-CNz*ekm_plus*dt/dz**2/(0.5*(drdt(i,j,kp)+drdt(i+1,j,kp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dUdt(i,j,0:k1)
	    CALL solve_tridiag(dUdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
      enddo
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))		
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))
     & *MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))			
	    aaa(k)=-CNz*ekm_min*dt/dz**2/(0.5*(drdt(i,j,km)+drdt(i,j+1,km)))
	    bbb(k)=1.+CNz*(ekm_min+ekm_plus)*dt/dz**2/(0.5*(drdt(i,j,k)+drdt(i,j+1,k)))
	    ccc(k)=-CNz*ekm_plus*dt/dz**2/(0.5*(drdt(i,j,kp)+drdt(i,j+1,kp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dVdt(i,j,0:k1)
	    CALL solve_tridiag(dVdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
      enddo

       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=ekm(i,j,k)*fc_global(i,j+rank*jmax,k)
	    ekm_plus=ekm(i,j,kp)*fc_global(i,j+rank*jmax,kp)
	    aaa(k)=-CNz*ekm_min*dt/dz**2/(0.5*(drdt(i,j,km)+drdt(i+1,j,k)))
	    bbb(k)=1.+CNz*(ekm_min+ekm_plus)*dt/dz**2/(0.5*(drdt(i,j,k)+drdt(i+1,j,kp)))
	    ccc(k)=-CNz*ekm_plus*dt/dz**2/(0.5*(drdt(i,j,kp)+drdt(i+1,j,kpp)))
            enddo
	    aaa(0)=0.
	    aaa(k1)=0.
	    ccc(0)=0.
	    ccc(k1)=0.
	    bbb(0)=1.
	    bbb(k1)=1.
	    rhss=dWdt(i,j,0:k1)
	    CALL solve_tridiag(dWdt(i,j,0:k1),aaa,bbb,ccc,rhss,k1+1) 
         enddo
       enddo
       ENDIF
		pold=p+pold    !what is called p here was dp in reality, now p is 			--> P^n 
		!IF (applyVOF.eq.1) THEN 
		  !IF (continuity_solver.eq.35) THEN
		  !  CALL state(cnew,rnew)		  
		  !  pold=pold/rnew(1:imax,1:jmax,1:kmax)
		  IF (continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
		    pold1=pold2 !-->P^n-2
			pold2=pold3 !-->P^n-1
			pold3=pold !--> P^n
		  ENDIF 
		  IF (pres_in_predictor_step.eq.0) THEN 
		    pold=0.
		  ENDIF		  
		!ENDIF 
		
      call pshiftb(pold,pplus) !,rank,imax,jmax,kmax,px)
c********************************************************************
c     CALCULATE pressure-gradient with old pressure:
c********************************************************************
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -dt * ( pold(i+1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
	  !!! do nothing on imax and i=1, so p(imax+1)=p(imax) --> dp/dn=0 at outflow and inflow 
	  if (periodicx.eq.1) then
	  !!! !periodic bc in x: 
	    i = imax
	      do k=1,kmax
		    do j=1,jmax
		      dUdt(i,j,k) = dUdt(i,j,k) -dt * ( pold(1,j,k) - pold(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		    enddo
	      enddo
	  endif
      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * (phip(j+1)-phip(j)) )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * (phip(jmax+1)-phip(jmax)) )
                enddo
	      enddo
      	  endif
	endif
      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          dWdt(i,j,k) = dWdt(i,j,k) -dt * ( pold(i,j,k+1) - pold(i,j,k) ) / dz
          enddo
        enddo
      enddo

		!IF (applyVOF.eq.1) THEN 
		! IF (continuity_solver.eq.35) THEN
		!  pold=pold*rnew(1:imax,1:jmax,1:kmax)
		!  rnew=1.
		! ENDIF
		!ENDIF 
		
      return
      end

      subroutine fillps
      USE nlist
      implicit none
	  include 'mpif.h'

	real xx,yy,r_orifice2,twodtdt,Qsource,rrri,rrrj,rrrk,rrrim,rrrjm,rrrkm
	real rrr2i,rrr2j,rrr2k,rrr2im,rrr2jm,rrr2km,ph_ref
	integer t,n,ierr
	real xTSHD(1:4),yTSHD(1:4),phi,ddrr,cin,s_in,dpdt(1:imax,1:jmax,1:kmax)
	real sum_c_ws,ctot,ws(nfrac),ccc
	integer inout,n2
	real upred(0:i1,0:j1,0:k1),vpred(0:i1,0:j1,0:k1),wpred(0:i1,0:j1,0:k1)
c
!       include 'param.txt'
!       include 'common.txt'
c
c**************************************************************
c
c       *** Fill the right hand for the poisson solver. ***
c
c**************************************************************
c
	IF (continuity_solver.eq.2) THEN ! Optional: 2 (neglect drdt)
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
			
      p(i,j,k)  =(
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz )
     +              ) /dt 
          enddo
        enddo
      enddo
	ELSEIF (continuity_solver.eq.3) THEN ! Optional: 3 (dudx=0) with input d(ur^*/r^n+1)dx

!! rhU,rhV,rhW are filled in bound_rhoU
	 
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =drdt(i,j,k)*(
     1  ( Ru(i)*dUdt(i,j,k)/rhU(i,j,k) - Ru(i-1)*dUdt(i-1,j,k)/rhU(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k)/rhV(i,j,k) -         dVdt(i,j-1,k)/rhV(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k)/rhW(i,j,k) -         dWdt(i,j,k-1)/rhW(i,j,k-1) ) / ( dz )
     +              ) /dt 
          enddo
        enddo
      enddo	
	ELSEIF (continuity_solver.eq.33) THEN ! Optional: 33 (dudx=0) improved version of 3 which gives dudx=e-14
	  dUdt(1:imax,1:jmax,1:kmax)=dUdt(1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)
	  dVdt(1:imax,1:jmax,1:kmax)=dVdt(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
	  dWdt(1:imax,1:jmax,1:kmax)=dWdt(1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)		
	  call bound_incljet(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz ) ) / dt 
          enddo
        enddo
      enddo	
	ELSEIF (continuity_solver.eq.34) THEN ! Optional: 34 (dudx=0) improved version of 3 which gives dudx=e-14 and with correction to account for difference Um and Uv
	 dUdt(1:imax,1:jmax,1:kmax)=dUdt(1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)
	 dVdt(1:imax,1:jmax,1:kmax)=dVdt(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
	 dWdt(1:imax,1:jmax,1:kmax)=dWdt(1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)		  
	 IF (slipvel.eq.2) THEN
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 	  
		    sum_c_ws=0.
			do n=1,nfrac
				ws(n)=-frac(n)%ws !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ELSE
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt  
		    sum_c_ws=0.
		    ctot=0.
		    do n=1,nfrac
				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
			enddo
			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
			do n=1,nfrac
				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ENDIF
	  call bound_incljet(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
!	 IF (slipvel.eq.2) THEN  ! make Wmix at interface kbed zero --> Wvol=-sum_c_ws 
!	  do j=0,j1
!	    do i=0,i1 
!		    k=kbed(i,j) 	  
!		    sum_c_ws=0.
!			do n=1,nfrac
!				ws(n)=-frac(n)%ws !ws defined positive downwards
!				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
!			enddo
!			dWdt(i,j,k)=-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
!        enddo
!      enddo
!	 ELSE
!	  do j=0,j1
!	    do i=0,i1 
!		    k=kbed(i,j)  
!		    sum_c_ws=0.
!		    ctot=0.
!		    do n=1,nfrac
!				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
!			enddo
!			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
!			do n=1,nfrac
!				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
!				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
!			enddo
!			dWdt(i,j,k)=-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
!        enddo
!      enddo
!	 ENDIF 
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz ) ) / dt 
          enddo
        enddo
      enddo	
	ELSEIF (continuity_solver.eq.35) THEN ! Optional: 35 (dudx=0) improved version of 33 (better treatment variable density in source Poisson) which gives dudx=e-14
	  dUdt(1:imax,1:jmax,1:kmax)=dUdt(1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)
	  dVdt(1:imax,1:jmax,1:kmax)=dVdt(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
	  dWdt(1:imax,1:jmax,1:kmax)=dWdt(1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)		
	  call bound_incljet(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
	  !dp=0.
	  IF (pres_in_predictor_step.eq.0) THEN 
		IF (oPRHO.eq.3) THEN 
!			phdt(:,:,kmax)=0.
!			do k=kmax-1,1,-1
!				phdt(:,:,k)=phdt(:,:,k+1)+(drdt(:,:,k)-rho_b)*ABS(gz)*dz
!			enddo
!			if (rank.eq.0) then
!			  ph_ref=phdt(imax,1,1)
!			endif
!			call mpi_bcast(ph_ref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!			phdt=phdt-ph_ref	
!		    dp(1:imax,1:jmax,1:kmax)= pold3+phdt-phnew								!dp contains predicted pressure 
!			phnew = phdt 
			dp(1:imax,1:jmax,1:kmax) = (pold3+pold2)*0.5+dt*(pold3-pold1)/(time_n-time_nm2) ! seems to be a slight flip-flop in time which is tried to be prevented 
		ELSEIF (oPRHO.eq.4) THEN 
			dp(1:imax,1:jmax,1:kmax)= dt*(3.*pold3-4.*pold2+pold1)/(3.*time_n-4.*time_nm+time_nm2)+pold3
		ELSEIF (oPRHO.eq.1) THEN 
			dp(1:imax,1:jmax,1:kmax)= pold3  										!dp contains predicted pressure 
		ELSE 
			dp(1:imax,1:jmax,1:kmax)= dt/dt_old*(pold3-pold2)+pold3  !pold3 		!dp contains predicted pressure
		ENDIF
	  ELSE 
		IF (oPRHO.eq.3) THEN 
!			phdt(:,:,kmax)=0.
!			do k=kmax-1,1,-1
!				phdt(:,:,k)=phdt(:,:,k+1)+(drdt(:,:,k)-rho_b)*ABS(gz)*dz
!			enddo
!			if (rank.eq.0) then
!			  ph_ref=phdt(imax,1,1)
!			endif
!			call mpi_bcast(ph_ref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!			phdt=phdt-ph_ref	
!		    dp(1:imax,1:jmax,1:kmax)= phdt-phnew									!dp contains predicted increment of pressure  
			dp(1:imax,1:jmax,1:kmax) = dt*(pold3-pold1)/(time_n-time_nm2) ! seems to be a slight flip-flop in time which is tried to be prevented 
		ELSEIF (oPRHO.eq.4) THEN 
			dp(1:imax,1:jmax,1:kmax)= dt*(3.*pold3-4.*pold2+pold1)/(3.*time_n-4.*time_nm+time_nm2)
		ELSEIF (oPRHO.eq.1) THEN 
			dp(1:imax,1:jmax,1:kmax)= 0.  											!dp contains predicted increment of pressure 
		ELSE 
			dp(1:imax,1:jmax,1:kmax)= dt/dt_old*(pold3-pold2) 						!dp contains predicted increment of pressure 
		ENDIF	
	  ENDIF	  
	  call bound_3D(dp)
	  !write(*,*),'rank,pold,p,dp',rank,pold(45,1,1),dp(45,1,1),p(45,1,1)
	  if (applyVOF.eq.1) then 
		call state_edges(cU,rhU)
		call state_edges(cV,rhV)
		call state_edges(cW,rhW)
	  endif
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz ) ) / dt  + 
     4  ((Ru(i  )*(1./rho_b2-1./rhU(i  ,j,k))*(dp(i+1,j,k)-dp(i  ,j,k))/(Rp(i+1)-Rp(i  )) - 
     4    Ru(i-1)*(1./rho_b2-1./rhU(i-1,j,k))*(dp(i  ,j,k)-dp(i-1,j,k))/(Rp(i  )-Rp(i-1)) ) / ( Rp(i)*dr(i) )
     +              +
     5  ( (1./rho_b2-1./rhV(i,j  ,k))*(dp(i,j+1,k)-dp(i,j  ,k))/ (Rp(i)*(phip(j+1)-phip(j  )))-         
     5    (1./rho_b2-1./rhV(i,j-1,k))*(dp(i,j  ,k)-dp(i,j-1,k))/ (Rp(i)*(phip(j  )-phip(j-1)))) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     6  ( (1./rho_b2-1./rhW(i,j,k  ))*(dp(i,j,k+1)-dp(i,j,k  ))/ dz -
     6    (1./rho_b2-1./rhW(i,j,k-1))*(dp(i,j,k  )-dp(i,j,k-1))/ dz	) / dz ) 	 
          enddo
        enddo
      enddo	
	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
		!  if (continuity_solver.eq.33.or.continuity_solver.eq.34.or.continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
			p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				  			! div(u)=0 --> total volume flux in 
		!  else
		!	p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				! div(u)=0 --> total volume flux in 
		!  endif
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop	  
	 do n=1,npresPRHO
	  ! 2nd iteration loop with updated P-predictor:
	  IF (poissolver.eq.3) THEN
	   CALL SOLVEpois_vg_pardiso(p)		  
	  ELSE
	   CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
	  ENDIF
	  dp(1:imax,1:jmax,1:kmax)=p*rho_b2 !0.8*p*rho_b+0.2*dp(1:imax,1:jmax,1:kmax) !use first iteration loop to obtain more accurate p-predictor 
	  call bound_3D(dp)	  
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz ) ) / dt  + 
     4  ((Ru(i  )*(1./rho_b2-1./rhU(i  ,j,k))*(dp(i+1,j,k)-dp(i  ,j,k))/(Rp(i+1)-Rp(i  )) - 
     4    Ru(i-1)*(1./rho_b2-1./rhU(i-1,j,k))*(dp(i  ,j,k)-dp(i-1,j,k))/(Rp(i  )-Rp(i-1)) ) / ( Rp(i)*dr(i) )
     +              +
     5  ( (1./rho_b2-1./rhV(i,j  ,k))*(dp(i,j+1,k)-dp(i,j  ,k))/ (Rp(i)*(phip(j+1)-phip(j  )))-         
     5    (1./rho_b2-1./rhV(i,j-1,k))*(dp(i,j  ,k)-dp(i,j-1,k))/ (Rp(i)*(phip(j  )-phip(j-1)))) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     6  ( (1./rho_b2-1./rhW(i,j,k  ))*(dp(i,j,k+1)-dp(i,j,k  ))/ dz -
     6    (1./rho_b2-1./rhW(i,j,k-1))*(dp(i,j,k  )-dp(i,j,k-1))/ dz	) / dz ) 	 
          enddo
        enddo
      enddo
	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
		!  if (continuity_solver.eq.33.or.continuity_solver.eq.34.or.continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
			p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				  			! div(u)=0 --> total volume flux in 
		!  else
		!	p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				! div(u)=0 --> total volume flux in 
		!  endif
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop	  
	 enddo
	 	   
	ELSEIF (continuity_solver.eq.36) THEN ! Optional: 36 (dudx=0) improved version of 34 (better treatment variable density in source Poisson) which gives dudx=e-14 and with correction to account for difference Um and Uv
	 dUdt(1:imax,1:jmax,1:kmax)=dUdt(1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)
	 dVdt(1:imax,1:jmax,1:kmax)=dVdt(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
	 dWdt(1:imax,1:jmax,1:kmax)=dWdt(1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)	
	  if (applyVOF.eq.1) then 
		call state_edges(cU,rhU)
		call state_edges(cV,rhV)
		call state_edges(cW,rhW)
	  endif	 
	 IF (slipvel.eq.2) THEN
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 	  
		    sum_c_ws=0.
			do n=1,nfrac
				ws(n)=-frac(n)%ws !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ELSE
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt  
		    sum_c_ws=0.
		    ctot=0.
		    do n=1,nfrac
				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
			enddo
			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
			do n=1,nfrac
				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ENDIF
	  call bound_incljet(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
!	 IF (slipvel.eq.2) THEN  ! make Wmix at interface kbed zero --> Wvol=-sum_c_ws 
!	  do j=0,j1
!	    do i=0,i1 
!		    k=kbed(i,j) 	  
!		    sum_c_ws=0.
!			do n=1,nfrac
!				ws(n)=-frac(n)%ws !ws defined positive downwards
!				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
!			enddo
!			dWdt(i,j,k)=-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
!        enddo
!      enddo
!	 ELSE
!	  do j=0,j1
!	    do i=0,i1 
!		    k=kbed(i,j)  
!		    sum_c_ws=0.
!		    ctot=0.
!		    do n=1,nfrac
!				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
!			enddo
!			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
!			do n=1,nfrac
!				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
!				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
!			enddo
!			dWdt(i,j,k)=-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
!        enddo
!      enddo
!	 ENDIF 
	  IF (pres_in_predictor_step.eq.0) THEN 
		IF (oPRHO.eq.3) THEN 
!			phdt(:,:,kmax)=0.
!			do k=kmax-1,1,-1
!				phdt(:,:,k)=phdt(:,:,k+1)+(drdt(:,:,k)-rho_b)*ABS(gz)*dz
!			enddo
!			if (rank.eq.0) then
!			  ph_ref=phdt(imax,1,1)
!			endif
!			call mpi_bcast(ph_ref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!			phdt=phdt-ph_ref	
!		    dp(1:imax,1:jmax,1:kmax)= pold3+phdt-phnew								!dp contains predicted pressure 
!			phnew = phdt 
			dp(1:imax,1:jmax,1:kmax) = (pold3+pold2)*0.5+dt*(pold3-pold1)/(time_n-time_nm2) ! seems to be a slight flip-flop in time which is tried to be prevented 
		ELSEIF (oPRHO.eq.4) THEN 
			dp(1:imax,1:jmax,1:kmax)= dt*(3.*pold3-4.*pold2+pold1)/(3.*time_n-4.*time_nm+time_nm2)+pold3
		ELSEIF (oPRHO.eq.1) THEN 
			dp(1:imax,1:jmax,1:kmax)= pold3  										!dp contains predicted pressure 
		ELSE 
			dp(1:imax,1:jmax,1:kmax)= dt/dt_old*(pold3-pold2)+pold3  !pold3 		!dp contains predicted pressure
		ENDIF
	  ELSE 
		IF (oPRHO.eq.3) THEN 
!			phdt(:,:,kmax)=0.
!			do k=kmax-1,1,-1
!				phdt(:,:,k)=phdt(:,:,k+1)+(drdt(:,:,k)-rho_b)*ABS(gz)*dz
!			enddo
!			if (rank.eq.0) then
!			  ph_ref=phdt(imax,1,1)
!			endif
!			call mpi_bcast(ph_ref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!			phdt=phdt-ph_ref	
!		    dp(1:imax,1:jmax,1:kmax)= phdt-phnew									!dp contains predicted increment of pressure  
!			phnew = phdt 
			dp(1:imax,1:jmax,1:kmax) = dt*(pold3-pold1)/(time_n-time_nm2) ! seems to be a slight flip-flop in time which is tried to be prevented 
		ELSEIF (oPRHO.eq.4) THEN 
			dp(1:imax,1:jmax,1:kmax)= dt*(3.*pold3-4.*pold2+pold1)/(3.*time_n-4.*time_nm+time_nm2)
		ELSEIF (oPRHO.eq.1) THEN 
			dp(1:imax,1:jmax,1:kmax)= 0.  											!dp contains predicted increment of pressure 
		ELSE 
			dp(1:imax,1:jmax,1:kmax)= dt/dt_old*(pold3-pold2) 						!dp contains predicted increment of pressure 
		ENDIF		   
	  ENDIF
	  call bound_3D(dp)
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz ) ) / dt +
     4  ((Ru(i  )*(1./rho_b2-1./rhU(i  ,j,k))*(dp(i+1,j,k)-dp(i  ,j,k))/(Rp(i+1)-Rp(i  )) - 
     4    Ru(i-1)*(1./rho_b2-1./rhU(i-1,j,k))*(dp(i  ,j,k)-dp(i-1,j,k))/(Rp(i  )-Rp(i-1)) ) / ( Rp(i)*dr(i) )
     +              +
     5  ( (1./rho_b2-1./rhV(i,j  ,k))*(dp(i,j+1,k)-dp(i,j  ,k))/ (Rp(i)*(phip(j+1)-phip(j  )))-         
     5    (1./rho_b2-1./rhV(i,j-1,k))*(dp(i,j  ,k)-dp(i,j-1,k))/ (Rp(i)*(phip(j  )-phip(j-1)))) / (Rp(i)*(phiv(j)-phiv(j-1)))
     +              +
     6  ( (1./rho_b2-1./rhW(i,j,k  ))*(dp(i,j,k+1)-dp(i,j,k  )) / dz - 
     6    (1./rho_b2-1./rhW(i,j,k-1))*(dp(i,j,k  )-dp(i,j,k-1)) / dz	) / dz ) 
          enddo
        enddo
      enddo	
	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
		!  if (continuity_solver.eq.33.or.continuity_solver.eq.34.or.continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
			p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				  			! div(u)=0 --> total volume flux in 
		!  else
		!	p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				! div(u)=0 --> total volume flux in 
		!  endif
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop
	
	 do n=1,npresPRHO
	  ! 2nd iteration loop with updated P-predictor:
	  IF (poissolver.eq.3) THEN
	   CALL SOLVEpois_vg_pardiso(p)		  
	  ELSE
	   CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
	  ENDIF
	  dp(1:imax,1:jmax,1:kmax)=p*rho_b2 !use first iteration loop to obtain more accurate p-predictor 
	  call bound_3D(dp)	   
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz ) ) / dt  + 
     4  ((Ru(i  )*(1./rho_b2-1./rhU(i  ,j,k))*(dp(i+1,j,k)-dp(i  ,j,k))/(Rp(i+1)-Rp(i  )) - 
     4    Ru(i-1)*(1./rho_b2-1./rhU(i-1,j,k))*(dp(i  ,j,k)-dp(i-1,j,k))/(Rp(i  )-Rp(i-1)) ) / ( Rp(i)*dr(i) )
     +              +
     5  ( (1./rho_b2-1./rhV(i,j  ,k))*(dp(i,j+1,k)-dp(i,j  ,k))/ (Rp(i)*(phip(j+1)-phip(j  )))-         
     5    (1./rho_b2-1./rhV(i,j-1,k))*(dp(i,j  ,k)-dp(i,j-1,k))/ (Rp(i)*(phip(j  )-phip(j-1)))) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     6  ( (1./rho_b2-1./rhW(i,j,k  ))*(dp(i,j,k+1)-dp(i,j,k  ))/ dz -
     6    (1./rho_b2-1./rhW(i,j,k-1))*(dp(i,j,k  )-dp(i,j,k-1))/ dz	) / dz ) 	 
          enddo
        enddo
      enddo	
	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
		!  if (continuity_solver.eq.33.or.continuity_solver.eq.34.or.continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
			p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				  			! div(u)=0 --> total volume flux in 
		!  else
		!	p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				! div(u)=0 --> total volume flux in 
		!  endif
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop
	  
	 enddo	

	ELSEIF (continuity_solver.eq.38) THEN ! Optional: 38 (dudx=0), alternative to 36 (also better treatment variable density in source Poisson) which gives dudx=e-14 and with correction to account for difference Um and Uv
	  if (applyVOF.eq.1) then 
		call state_edges(cU,rhU)
		call state_edges(cV,rhV)
		call state_edges(cW,rhW)
	  endif	 
	 IF (slipvel.eq.2) THEN
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 	  
		    sum_c_ws=0.
			do n=1,nfrac
				ws(n)=-frac(n)%ws !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)-sum_c_ws*rhW(i,j,k)  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ELSE
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt  
		    sum_c_ws=0.
		    ctot=0.
		    do n=1,nfrac
				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
			enddo
			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
			do n=1,nfrac
				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)-sum_c_ws*rhW(i,j,k)  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ENDIF
	 call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),monopile,time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,
     &	 Ub3new,Vb3new,Wb3new)
	 IF (oPRHO.eq.1) THEN !use Unew for contribution of u*div(rho)
			upred = Unew
			vpred = Vnew
			wpred = Wnew
			 IF (slipvel.eq.2) THEN
			  do j=1,jmax
				do i=1,imax   
				  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 	  
					sum_c_ws=0.
					do n=1,nfrac
						ws(n)=-frac(n)%ws !ws defined positive downwards
						ccc=0.5*(cnew(n,i,j,k)+cnew(n,i,j,k+1))
						sum_c_ws=sum_c_ws+ws(n)*(ccc*frac(n)%rho/(0.5*(rnew(i,j,k)+rnew(i,j,k+1)))-ccc)
					enddo
					wpred(i,j,k)=wpred(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
				  enddo
				enddo
			  enddo
			 ELSE
			  do j=1,jmax
				do i=1,imax   
				  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt  
					sum_c_ws=0.
					ctot=0.
					do n=1,nfrac
						ccc=0.5*(cnew(n,i,j,k)+cnew(n,i,j,k+1))
						ctot=ccc*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
					enddo
					ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
					do n=1,nfrac
						ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
						ccc=0.5*(cnew(n,i,j,k)+cnew(n,i,j,k+1))
						sum_c_ws=sum_c_ws+ws(n)*(ccc*frac(n)%rho/(0.5*(rnew(i,j,k)+rnew(i,j,k+1)))-ccc)
					enddo
					wpred(i,j,k)=wpred(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
				  enddo
				enddo
			  enddo
			 ENDIF			
	 ELSE  !use U* for contribution of u*div(rho)
			upred = dUdt/rhU 
			vpred = dVdt/rhV
			wpred = dWdt/rhW
	 ENDIF
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz ) ) / dt - (
     1  0.5*(upred(i,j,k)+upred(i-1,j,k))*( Ru(i)*rhU(i,j,k) - Ru(i-1)*rhU(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  0.5*(vpred(i,j,k)+vpred(i,j-1,k))*(       rhV(i,j,k) -         rhV(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  0.5*(wpred(i,j,k)+wpred(i,j,k-1))*(       rhW(i,j,k) -         rhW(i,j,k-1) ) / ( dz ) ) / dt 	 
          enddo
        enddo
      enddo	
	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
		!  if (continuity_solver.eq.33.or.continuity_solver.eq.34.or.continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
		!	p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				  			! div(u)=0 --> total volume flux in 
		!  else
			p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				! div(u)=0 --> total volume flux in 
		!  endif
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop
	 
	ELSE !default is 1 (drdt+drudx=0)
	twodtdt=MAX((3.*time_np-4.*time_n+time_nm)*dt,1.e-12)

      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
c
      p(i,j,k)  =(
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz )
     +              ) /dt  +  (3.*drdt(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/twodtdt ! normally this is used before 25-10-2018
!     +              ) /dt  +  (drdt(i,j,k)-rnew(i,j,k))/((time_np-time_n)*dt) !test  
          enddo
        enddo
      enddo
!	DO n2=1,nbedplume
!    	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end)) THEN
!    	! rotation ship for ambient side current
!    	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
!    	  phi=atan2(V_b,1.e-12)
!    	else
!    	  phi=atan2(V_b,(U_TSHD-U_b))
!    	endif
!      do k=1,kmax
!       do i=1,imax  
!         do j=jmax,1,-1 
!	  xx=Rp(i)*cos_u(j)-schuif_x
!	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
!		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
!		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
!	  ELSE 
!	 	inout=0
!	  ENDIF
!	  if (inout.eq.1) then
!		ddrr=0.
!		do n=1,nfrac
!		  ddrr = ddrr + bp(n2)%sedflux(n)*dt/bp(n2)%volncells +
!     &		  dcdt(n,i,j,k)*MIN(0.,bp(n2)%Q*bp(n2)%changesedsuction)*frac(n)%rho*dt/bp(n2)%volncells  
!			! change of rho in cell due to sedflux or suction Q --> correct determination drhodt in continuity equation for this
!			! IMPLICIT: c^n+1-c^n=-Qout_cel/Vol_cel*dt*c^n+1 --> c^n+1 = c^n/(1+Qout_cel/Vol_cel*dt)
!			!when Q negative, remove sediment from cell as well   IMPLICIT 	 
!			p(i,j,k)=p(i,j,k)-bp(n2)%sedflux(n)/bp(n2)%volncells/dt  ! total mass flux in 
!		enddo
!		!p(i,j,k) = p(i,j,k) - ddrr/dt/dt
!	 
!	   endif
!	  enddo
!	 enddo
!	enddo
!	ENDIF	
!	ENDDO ! bedplume loop	

	ENDIF ! continuity_solver loop

	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
!		cin=0.
!		s_in=0.
!		ddrr=0.		
!		do n=1,nfrac
!		  cin = cin + bp(n2)%sedflux(n)/(frac(n)%rho*bp(n2)%Q)
!		  s_in = s_in + bp(n2)%sedflux(n) 
!		  ddrr = ddrr + bp(n2)%sedflux(n)*dt/bp(n2)%volncells +
!     &		  dcdt(n,i,j,k)*MIN(0.,bp(n2)%Q*bp(n2)%changesedsuction)*frac(n)%rho*dt/bp(n2)%volncells  
!			! change of rho in cell due to sedflux or suction Q --> correct determination drhodt in continuity equation for this
!			! IMPLICIT: c^n+1-c^n=-Qout_cel/Vol_cel*dt*c^n+1 --> c^n+1 = c^n/(1+Qout_cel/Vol_cel*dt)
!			!when Q negative, remove sediment from cell as well   IMPLICIT 		  
!		enddo	
		  if (continuity_solver.eq.33.or.continuity_solver.eq.34) THEN
			p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				  			! div(u)=0 --> total volume flux in 
		  elseif (continuity_solver.eq.35.or.continuity_solver.eq.36) THEN 
			! bp%Q already added to p input in lines above 
		  else
			p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/dt  				! div(u)=0 --> total volume flux in 
		  endif
!!			p(i,j,k)=p(i,j,k)-(s_in+(1.-MIN(cin,1.))*drdt(i,j,k)*bp(n2)%Q)/bp(n2)%volncells/dt  ! total mass flux in 
!			p(i,j,k)=p(i,j,k)-((1.-MIN(cin,1.))*drdt(i,j,k)*bp(n2)%Q)/bp(n2)%volncells/dt  ! total mass flux in 
!!			p(i,j,k)=p(i,j,k)-rho_b*bp(n2)%Q/bp(n2)%volncells/dt
!		  endif
!		  !bp(n2)%Q positive means influx
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop


      return
      end

      subroutine fillps2(uu,vv,ww,rr,tt,ddtt)
      USE nlist
      implicit none

	real xx,yy,r_orifice2,twodtdt,Qsource,rrri,rrrj,rrrk,rrrim,rrrjm,rrrkm
	real tt,ddtt,uu(0:i1,0:j1,0:k1),vv(0:i1,0:j1,0:k1),ww(0:i1,0:j1,0:k1),rr(0:i1,0:j1,0:k1)
	integer t,n
	real xTSHD(1:4),yTSHD(1:4),phi
	integer inout,n2	
	real sum_c_ws,ws(1:nfrac),ctot 
c
!       include 'param.txt'
!       include 'common.txt'
c
c**************************************************************
c
c       *** Fill the right hand for the poisson solver. ***
c
c**************************************************************
c

	IF (continuity_solver.eq.2) THEN ! Optional: 2 (neglect drdt)
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
c
      p(i,j,k)  =(
     1  ( Ru(i)*uu(i,j,k) - Ru(i-1)*uu(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       vv(i,j,k) -         vv(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       ww(i,j,k) -         ww(i,j,k-1) ) / ( dz )
     +              ) /ddtt  
          enddo
        enddo
      enddo
	ELSEIF (continuity_solver.eq.3) THEN ! Optional: 3 (dudx=0)
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
		   rrri=0.5*(rr(i,j,k)+rr(i+1,j,k))
		   rrrj=0.5*(rr(i,j,k)+rr(i,j+1,k))
		   rrrk=0.5*(rr(i,j,k)+rr(i,j,k+1))
		   rrrim=0.5*(rr(i,j,k)+rr(i-1,j,k))
		   rrrjm=0.5*(rr(i,j,k)+rr(i,j-1,k))
		   rrrkm=0.5*(rr(i,j,k)+rr(i,j,k-1))
      p(i,j,k)  =rr(i,j,k)*(
     1  ( Ru(i)*uu(i,j,k)/rrri - Ru(i-1)*uu(i-1,j,k)/rrrim ) / ( Rp(i)*dr(i) )
     +              +
     2  (       vv(i,j,k)/rrrj -         vv(i,j-1,k)/rrrjm ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       ww(i,j,k)/rrrk -         ww(i,j,k-1)/rrrkm ) / ( dz )
     +              ) /ddtt  
          enddo
        enddo
      enddo
	ELSEIF (continuity_solver.eq.33) THEN ! Optional: 33 (dudx=0) improved version of 3 which gives dudx=e-14
	  uu(1:imax,1:jmax,1:kmax)=uu(1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)
	  vv(1:imax,1:jmax,1:kmax)=vv(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
	  ww(1:imax,1:jmax,1:kmax)=ww(1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)		
	  call bound_incljet(uu,vv,ww,drdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)	
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*uu(i,j,k) - Ru(i-1)*uu(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       vv(i,j,k) -         vv(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       ww(i,j,k) -         ww(i,j,k-1) ) / ( dz ) ) / ddtt 
          enddo
        enddo
      enddo	
	ELSEIF (continuity_solver.eq.34) THEN ! Optional: 34 (dudx=0) improved version of 3 which gives dudx=e-14 and with correction to account for difference Um and Uv
	 uu(1:imax,1:jmax,1:kmax)=uu(1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)
	 vv(1:imax,1:jmax,1:kmax)=vv(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
	 ww(1:imax,1:jmax,1:kmax)=ww(1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)		  
	 IF (slipvel.eq.2) THEN
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 	  
		    sum_c_ws=0.
			do n=1,nfrac
				ws(n)=-frac(n)%ws !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			ww(i,j,k)=ww(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ELSE
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt  
		    sum_c_ws=0.
		    ctot=0.
		    do n=1,nfrac
				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
			enddo
			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
			do n=1,nfrac
				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			ww(i,j,k)=ww(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ENDIF
	  call bound_incljet(uu,vv,ww,drdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)	
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*uu(i,j,k) - Ru(i-1)*uu(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       vv(i,j,k) -         vv(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       ww(i,j,k) -         ww(i,j,k-1) ) / ( dz ) ) / ddtt
          enddo
        enddo
      enddo	
	ELSEIF (continuity_solver.eq.35) THEN ! Optional: 35 (dudx=0) improved version of 33 (better treatment variable density in source Poisson) which gives dudx=e-14
	  uu(1:imax,1:jmax,1:kmax)=uu(1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)
	  vv(1:imax,1:jmax,1:kmax)=vv(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
	  ww(1:imax,1:jmax,1:kmax)=ww(1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)		
	  call bound_incljet(uu,vv,ww,drdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)	
	  IF (pres_in_predictor_step.eq.0) THEN 
		   dp(1:imax,1:jmax,1:kmax)=ddtt/dt_old*(pold3-pold2)+pold3 		!dp contains predicted pressure 
	  ELSE 
		   dp(1:imax,1:jmax,1:kmax)=ddtt/dt_old*(pold3-pold2)		!dp contains predicted increment of pressure 
	  ENDIF	  
	  call bound_3D(dp)
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*uu(i,j,k) - Ru(i-1)*uu(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       vv(i,j,k) -         vv(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       ww(i,j,k) -         ww(i,j,k-1) ) / ( dz ) ) / ddtt  + 
     4  ((Ru(i  )*(1./rho_b2-1./rhU(i  ,j,k))*(dp(i+1,j,k)-dp(i  ,j,k))/(Rp(i+1)-Rp(i  )) - 
     4    Ru(i-1)*(1./rho_b2-1./rhU(i-1,j,k))*(dp(i  ,j,k)-dp(i-1,j,k))/(Rp(i  )-Rp(i-1)) ) / ( Rp(i)*dr(i) )
     +              +
     5  ( (1./rho_b2-1./rhV(i,j  ,k))*(dp(i,j+1,k)-dp(i,j  ,k))/ (Rp(i)*(phip(j+1)-phip(j  )))-         
     5    (1./rho_b2-1./rhV(i,j-1,k))*(dp(i,j  ,k)-dp(i,j-1,k))/ (Rp(i)*(phip(j  )-phip(j-1)))) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     6  ( (1./rho_b2-1./rhW(i,j,k  ))*(dp(i,j,k+1)-dp(i,j,k  ))/ dz -
     6    (1./rho_b2-1./rhW(i,j,k-1))*(dp(i,j,k  )-dp(i,j,k-1))/ dz	) / dz ) 	 
          enddo
        enddo
      enddo	
	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
		!  if (continuity_solver.eq.33.or.continuity_solver.eq.34.or.continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
			p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/ddtt  				  			! div(u)=0 --> total volume flux in 
		!  else
		!	p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/ddtt  				! div(u)=0 --> total volume flux in 
		!  endif
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop	  
	 do n=1,npresPRHO
	  ! 2nd iteration loop with updated P-predictor:
	  IF (poissolver.eq.3) THEN
	   CALL SOLVEpois_vg_pardiso(p)		  
	  ELSE
	   CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
	  ENDIF
	  dp(1:imax,1:jmax,1:kmax)=p*rho_b2 !0.8*p*rho_b+0.2*dp(1:imax,1:jmax,1:kmax) !use first iteration loop to obtain more accurate p-predictor 
	  call bound_3D(dp)	   
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*uu(i,j,k) - Ru(i-1)*uu(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       vv(i,j,k) -         vv(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       ww(i,j,k) -         ww(i,j,k-1) ) / ( dz ) ) / ddtt  + 
     4  ((Ru(i  )*(1./rho_b2-1./rhU(i  ,j,k))*(dp(i+1,j,k)-dp(i  ,j,k))/(Rp(i+1)-Rp(i  )) - 
     4    Ru(i-1)*(1./rho_b2-1./rhU(i-1,j,k))*(dp(i  ,j,k)-dp(i-1,j,k))/(Rp(i  )-Rp(i-1)) ) / ( Rp(i)*dr(i) )
     +              +
     5  ( (1./rho_b2-1./rhV(i,j  ,k))*(dp(i,j+1,k)-dp(i,j  ,k))/ (Rp(i)*(phip(j+1)-phip(j  )))-         
     5    (1./rho_b2-1./rhV(i,j-1,k))*(dp(i,j  ,k)-dp(i,j-1,k))/ (Rp(i)*(phip(j  )-phip(j-1)))) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     6  ( (1./rho_b2-1./rhW(i,j,k  ))*(dp(i,j,k+1)-dp(i,j,k  ))/ dz -
     6    (1./rho_b2-1./rhW(i,j,k-1))*(dp(i,j,k  )-dp(i,j,k-1))/ dz	) / dz ) 	 
          enddo
        enddo
      enddo
	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
		!  if (continuity_solver.eq.33.or.continuity_solver.eq.34.or.continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
			p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/ddtt  				  			! div(u)=0 --> total volume flux in 
		!  else
		!	p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/ddtt  				! div(u)=0 --> total volume flux in 
		!  endif
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop	  
	 enddo
	 	   
	ELSEIF (continuity_solver.eq.36) THEN ! Optional: 36 (dudx=0) improved version of 34 (better treatment variable density in source Poisson) which gives dudx=e-14 and with correction to account for difference Um and Uv
	 uu(1:imax,1:jmax,1:kmax)=uu(1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)
	 vv(1:imax,1:jmax,1:kmax)=vv(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
	 ww(1:imax,1:jmax,1:kmax)=ww(1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)		  
	 IF (slipvel.eq.2) THEN
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 	  
		    sum_c_ws=0.
			do n=1,nfrac
				ws(n)=-frac(n)%ws !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			ww(i,j,k)=ww(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ELSE
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt  
		    sum_c_ws=0.
		    ctot=0.
		    do n=1,nfrac
				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
			enddo
			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
			do n=1,nfrac
				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			ww(i,j,k)=ww(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ENDIF
	  call bound_incljet(uu,vv,ww,drdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)	
	  IF (pres_in_predictor_step.eq.0) THEN 
		   dp(1:imax,1:jmax,1:kmax)=ddtt/dt_old*(pold3-pold2)+pold3 		!dp contains predicted pressure 
	  ELSE 
		   dp(1:imax,1:jmax,1:kmax)=ddtt/dt_old*(pold3-pold2)		!dp contains predicted increment of pressure 
	  ENDIF
	  call bound_3D(dp)
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*uu(i,j,k) - Ru(i-1)*uu(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       vv(i,j,k) -         vv(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       ww(i,j,k) -         ww(i,j,k-1) ) / ( dz ) ) / ddtt +
     4  ((Ru(i  )*(1./rho_b2-1./rhU(i  ,j,k))*(dp(i+1,j,k)-dp(i  ,j,k))/(Rp(i+1)-Rp(i  )) - 
     4    Ru(i-1)*(1./rho_b2-1./rhU(i-1,j,k))*(dp(i  ,j,k)-dp(i-1,j,k))/(Rp(i  )-Rp(i-1)) ) / ( Rp(i)*dr(i) )
     +              +
     5  ( (1./rho_b2-1./rhV(i,j  ,k))*(dp(i,j+1,k)-dp(i,j  ,k))/ (Rp(i)*(phip(j+1)-phip(j  )))-         
     5    (1./rho_b2-1./rhV(i,j-1,k))*(dp(i,j  ,k)-dp(i,j-1,k))/ (Rp(i)*(phip(j  )-phip(j-1)))) / (Rp(i)*(phiv(j)-phiv(j-1)))
     +              +
     6  ( (1./rho_b2-1./rhW(i,j,k  ))*(dp(i,j,k+1)-dp(i,j,k  )) / dz - 
     6    (1./rho_b2-1./rhW(i,j,k-1))*(dp(i,j,k  )-dp(i,j,k-1)) / dz	) / dz ) 
          enddo
        enddo
      enddo	
	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
		!  if (continuity_solver.eq.33.or.continuity_solver.eq.34.or.continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
			p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/ddtt  				  			! div(u)=0 --> total volume flux in 
		!  else
		!	p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/ddtt  				! div(u)=0 --> total volume flux in 
		!  endif
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop
	
	 do n=1,npresPRHO
	  ! 2nd iteration loop with updated P-predictor:
	  IF (poissolver.eq.3) THEN
	   CALL SOLVEpois_vg_pardiso(p)		  
	  ELSE
	   CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
	  ENDIF
	  dp(1:imax,1:jmax,1:kmax)=p*rho_b2 !use first iteration loop to obtain more accurate p-predictor 
	  call bound_3D(dp)	   
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
      p(i,j,k)  =(
     1  ( Ru(i)*uu(i,j,k) - Ru(i-1)*uu(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       vv(i,j,k) -         vv(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       ww(i,j,k) -         ww(i,j,k-1) ) / ( dz ) ) / ddtt  + 
     4  ((Ru(i  )*(1./rho_b2-1./rhU(i  ,j,k))*(dp(i+1,j,k)-dp(i  ,j,k))/(Rp(i+1)-Rp(i  )) - 
     4    Ru(i-1)*(1./rho_b2-1./rhU(i-1,j,k))*(dp(i  ,j,k)-dp(i-1,j,k))/(Rp(i  )-Rp(i-1)) ) / ( Rp(i)*dr(i) )
     +              +
     5  ( (1./rho_b2-1./rhV(i,j  ,k))*(dp(i,j+1,k)-dp(i,j  ,k))/ (Rp(i)*(phip(j+1)-phip(j  )))-         
     5    (1./rho_b2-1./rhV(i,j-1,k))*(dp(i,j  ,k)-dp(i,j-1,k))/ (Rp(i)*(phip(j  )-phip(j-1)))) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     6  ( (1./rho_b2-1./rhW(i,j,k  ))*(dp(i,j,k+1)-dp(i,j,k  ))/ dz -
     6    (1./rho_b2-1./rhW(i,j,k-1))*(dp(i,j,k  )-dp(i,j,k-1))/ dz	) / dz ) 	 
          enddo
        enddo
      enddo	
	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
		!  if (continuity_solver.eq.33.or.continuity_solver.eq.34.or.continuity_solver.eq.35.or.continuity_solver.eq.36) THEN
			p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/ddtt  				  			! div(u)=0 --> total volume flux in 
		!  else
		!	p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/ddtt  				! div(u)=0 --> total volume flux in 
		!  endif
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop	  
	 enddo	  
	ELSE !default is 1 (drdt+drudx=0)
	twodtdt=MAX((3.*tt-4.*time_n+time_nm)*ddtt,1.e-12)
	
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
c
      p(i,j,k)  =(
     1  ( Ru(i)*uu(i,j,k) - Ru(i-1)*uu(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       vv(i,j,k) -         vv(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       ww(i,j,k) -         ww(i,j,k-1) ) / ( dz )
!     +              ) /dt  +  (3.*drdt(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/(2*dt*dt)
     +              ) /ddtt  +  (3.*rr(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/twodtdt
          enddo
        enddo
      enddo
	ENDIF
	
	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
       do i=1,imax  
         do j=jmax,1,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		!! do not use bp(n2)%i,j indices as they are defined from 0,j1 instead of 1,jmax needed for pressure !!
		if (bp(n2)%radius.gt.0.) then 
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
!	  ELSE 
!	 	inout=0
!	  ENDIF
	  if (inout.eq.1) then
!		cin=0.
!		s_in=0.
!		ddrr=0.		
!		do n=1,nfrac
!		  cin = cin + bp(n2)%sedflux(n)/(frac(n)%rho*bp(n2)%Q)
!		  s_in = s_in + bp(n2)%sedflux(n) 
!		  ddrr = ddrr + bp(n2)%sedflux(n)*dt/bp(n2)%volncells +
!     &		  dcdt(n,i,j,k)*MIN(0.,bp(n2)%Q*bp(n2)%changesedsuction)*frac(n)%rho*dt/bp(n2)%volncells  
!			! change of rho in cell due to sedflux or suction Q --> correct determination drhodt in continuity equation for this
!			! IMPLICIT: c^n+1-c^n=-Qout_cel/Vol_cel*dt*c^n+1 --> c^n+1 = c^n/(1+Qout_cel/Vol_cel*dt)
!			!when Q negative, remove sediment from cell as well   IMPLICIT 		  
!		enddo	
		  if (continuity_solver.eq.33.or.continuity_solver.eq.34) THEN
			p(i,j,k)=p(i,j,k)-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/ddtt  				  			! div(u)=0 --> total volume flux in 
		  elseif (continuity_solver.eq.35.or.continuity_solver.eq.36) THEN 
			! bp%Q already added to p input in lines above 
		  else
			p(i,j,k)=p(i,j,k)-drdt(i,j,k)*bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/ddtt  				! div(u)=0 --> total volume flux in 
		  endif
!!			p(i,j,k)=p(i,j,k)-(s_in+(1.-MIN(cin,1.))*drdt(i,j,k)*bp(n2)%Q)/bp(n2)%volncells/dt  ! total mass flux in 
!			p(i,j,k)=p(i,j,k)-((1.-MIN(cin,1.))*drdt(i,j,k)*bp(n2)%Q)/bp(n2)%volncells/dt  ! total mass flux in 
!!			p(i,j,k)=p(i,j,k)-rho_b*bp(n2)%Q/bp(n2)%volncells/dt
!		  endif
!		  !bp(n2)%Q positive means influx
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop
 
      return
      end



      subroutine correc
      USE nlist
      implicit none
c
!       include 'param.txt'
!       include 'common.txt'
      real*8 pplus(imax,kmax)
	  real sum_c_ws,ctot,ws(nfrac)
      integer t,n


      call pshiftb(p,pplus) !,rank,imax,jmax,kmax,px)


      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -
     +      dt * ( p(i+1,j,k) - p(i,j,k) ) /
     +                ( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
	  !!! do nothing on imax and i=1, so p(imax+1)=p(imax) --> dp/dn=0 at outflow and inflow 
	  if (periodicx.eq.1) then
	  !!! !periodic bc in x: 
	    i = imax
	      do k=1,kmax
		    do j=1,jmax
		      dUdt(i,j,k) = dUdt(i,j,k) -dt * ( p(1,j,k) - p(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		    enddo
	      enddo
	  endif

      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -
     +     dt * ( p(i,j+1,k) - p(i,j,k) ) /
     +                ( Rp(i) * (phip(j+1)-phip(j)) )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
	      do k=1,kmax
		do i=1,imax
		dVdt(i,jmax,k) = dVdt(i,jmax,k) -
     +   dt * ( pplus(i,k) - p(i,jmax,k) ) /
     +              ( Rp(i) * (phip(jmax+1)-phip(jmax)) )
		enddo
	      enddo
	else
	      if (rank.ne.px-1) then
 	       do k=1,kmax
	        do i=1,imax
	        dVdt(i,jmax,k) = dVdt(i,jmax,k) -
     +   dt * ( pplus(i,k) - p(i,jmax,k) ) /
     +              ( Rp(i) * (phip(jmax+1)-phip(jmax)) )
	        enddo
	       enddo
      	      endif
	endif

      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          dWdt(i,j,k) = dWdt(i,j,k) -
     +     dt * ( p(i,j,k+1) - p(i,j,k) ) / dz
          enddo
        enddo
      enddo

!	  if (periodicx.eq.1.or.periodicy.eq.1) then
!      do k=1,kmax
!         do j=1,jmax
!            do i=1,imax
!	    dUdt(i,j,k) = dUdt(i,j,k)+Ppropx(i,j,k) 
!	    dVdt(i,j,k) = dVdt(i,j,k)+Ppropy(i,j,k) 	 !apply periodic driving forces after pressure correction step otherwise pressure gives exactly reaction force leading to same total flux as at t0
!            enddo
!         enddo
!      enddo
!	  endif

		! apply upper limit per fraction if defined:
	  DO n=1,nfrac
		IF (frac(n)%cmax<99.99) THEN 
			DO i=0,i1 
			  DO j=0,j1 
			    DO k=0,k1 
					dcdt(n,i,j,k)=MIN(dcdt(n,i,j,k),frac(n)%cmax)
				ENDDO 
			  ENDDO 
			ENDDO 
		ENDIF 
	  ENDDO 
	  
	  
      do k=0,k1
        do j=0,j1
          do i=0,i1
           Uold(i,j,k)=Unew(i,j,k)  !! ALS test even alle 4 de old=new statements uitgezet omdat volgens mij deze overbodig zijn
           Vold(i,j,k)=Vnew(i,j,k)
           Wold(i,j,k)=Wnew(i,j,k) 
	   do n=1,nfrac
             cold(n,i,j,k)=cnew(n,i,j,k)
             cnew(n,i,j,k)=dcdt(n,i,j,k)
	   enddo
          enddo
        enddo
      enddo 
        do j=0,j1
          do i=0,i1
	   do n=1,nfrac
             coldbot(n,i,j)=cnewbot(n,i,j)
             cnewbot(n,i,j)=dcdtbot(n,i,j)
	   enddo
          enddo
        enddo

	IF (continuity_solver.eq.33) THEN 
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
           Unew(i,j,k)=dUdt(i,j,k)
           Vnew(i,j,k)=dVdt(i,j,k)
           Wnew(i,j,k)=dWdt(i,j,k)
          enddo
        enddo
      enddo 	 
		p=p(1:imax,1:jmax,1:kmax)*drdt(1:imax,1:jmax,1:kmax) !scale P back with rho
	ELSEIF (continuity_solver.eq.34) THEN 
	 IF (slipvel.eq.2) THEN
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j),kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt
		    sum_c_ws=0.
			do n=1,nfrac
				ws(n)=-frac(n)%ws !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)+sum_c_ws  !go from velocity of volume centre to mixture velocity (centre of mass velocity)
          enddo
        enddo
      enddo
	 ELSE
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j),kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 
		    sum_c_ws=0.
		    ctot=0.
		    do n=1,nfrac
				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
			enddo
			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
			do n=1,nfrac
				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)+sum_c_ws  !go from velocity of volume centre to mixture velocity (centre of mass velocity)
          enddo
        enddo
      enddo
	 ENDIF	
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
           Unew(i,j,k)=dUdt(i,j,k)
           Vnew(i,j,k)=dVdt(i,j,k)
           Wnew(i,j,k)=dWdt(i,j,k)
          enddo
        enddo
      enddo 	 
		p=p(1:imax,1:jmax,1:kmax)*drdt(1:imax,1:jmax,1:kmax) !scale P back with rho	 
	ELSEIF (continuity_solver.eq.35) THEN 
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
		   dUdt(i,j,k)=dUdt(i,j,k)-dt*(1./rhU(i,j,k)-1./rho_b2)*(dp(i+1,j,k)-dp(i,j,k))/(Rp(i+1)-Rp(i))
		   dVdt(i,j,k)=dVdt(i,j,k)-dt*(1./rhV(i,j,k)-1./rho_b2)*(dp(i,j+1,k)-dp(i,j,k))/(Rp(i)*(phip(j+1)-phip(j)) ) 
		   dWdt(i,j,k)=dWdt(i,j,k)-dt*(1./rhW(i,j,k)-1./rho_b2)*(dp(i,j,k+1)-dp(i,j,k))/ dz		   
           Unew(i,j,k)=dUdt(i,j,k)
           Vnew(i,j,k)=dVdt(i,j,k)
           Wnew(i,j,k)=dWdt(i,j,k)
          enddo
        enddo
      enddo 	 
		p=p(1:imax,1:jmax,1:kmax)*rho_b2  !drdt(1:imax,1:jmax,1:kmax) !scale P back with rho_b not with drdt because in source pressure Poisson eq. already extra source-term 1/rho_b-1/drdt is included 
		!p=(p(1:imax,1:jmax,1:kmax)+(1./drdt(1:imax,1:jmax,1:kmax)+1./rho_b)*dp(1:imax,1:jmax,1:kmax))*rho_b 
		
	  if (applyVOF.eq.1.or.applyVOF.eq.2) then !make rhU,rhV,rhW 1 again 
		rhU=rho_b
		rhV=rho_b
		rhW=rho_b
		drdt=rho_b
	  endif		
	ELSEIF (continuity_solver.eq.36) THEN 
	 IF (slipvel.eq.2) THEN
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j),kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt
		    sum_c_ws=0.
			do n=1,nfrac
				ws(n)=-frac(n)%ws !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)+sum_c_ws  !go from velocity of volume centre to mixture velocity (centre of mass velocity)
          enddo
        enddo
      enddo	  
	 ELSE
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j),kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 
		    sum_c_ws=0.
		    ctot=0.
		    do n=1,nfrac
				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
			enddo
			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
			do n=1,nfrac
				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)+sum_c_ws  !go from velocity of volume centre to mixture velocity (centre of mass velocity)
          enddo
        enddo
      enddo
	 ENDIF	
	 
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
		   dUdt(i,j,k)=dUdt(i,j,k)-dt*(1./rhU(i,j,k)-1./rho_b2)*(dp(i+1,j,k)-dp(i,j,k))/(Rp(i+1)-Rp(i))
		   dVdt(i,j,k)=dVdt(i,j,k)-dt*(1./rhV(i,j,k)-1./rho_b2)*(dp(i,j+1,k)-dp(i,j,k))/(Rp(i)*(phip(j+1)-phip(j)) ) 
		   dWdt(i,j,k)=dWdt(i,j,k)-dt*(1./rhW(i,j,k)-1./rho_b2)*(dp(i,j,k+1)-dp(i,j,k))/ dz		  
           Unew(i,j,k)=dUdt(i,j,k)
           Vnew(i,j,k)=dVdt(i,j,k)
           Wnew(i,j,k)=dWdt(i,j,k)
          enddo
        enddo
      enddo 	
		p=p(1:imax,1:jmax,1:kmax)*rho_b2 !drdt(1:imax,1:jmax,1:kmax)  !scale P back with rho	  
	  if (applyVOF.eq.1.or.applyVOF.eq.2) then !make rhU,rhV,rhW 1 again 
		rhU=rho_b
		rhV=rho_b
		rhW=rho_b
		drdt=rho_b
	  endif		  

	ELSEIF (continuity_solver.eq.38) THEN 
	 IF (slipvel.eq.2) THEN
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j),kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt
		    sum_c_ws=0.
			do n=1,nfrac
				ws(n)=-frac(n)%ws !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)+sum_c_ws*rhW(i,j,k)  !go from velocity of volume centre to mixture velocity (centre of mass velocity)
          enddo
        enddo
      enddo	  
	 ELSE
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j),kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 
		    sum_c_ws=0.
		    ctot=0.
		    do n=1,nfrac
				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
			enddo
			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
			do n=1,nfrac
				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)+sum_c_ws*rhW(i,j,k)  !go from velocity of volume centre to mixture velocity (centre of mass velocity)
          enddo
        enddo
      enddo
	 ENDIF	
	 
		Unew(1:imax,1:jmax,1:kmax)=dUdt(1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)
		Vnew(1:imax,1:jmax,1:kmax)=dVdt(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
		Wnew(1:imax,1:jmax,1:kmax)=dWdt(1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)  
	  if (applyVOF.eq.1.or.applyVOF.eq.2) then !make rhU,rhV,rhW 1 again 
		rhU=rho_b
		rhV=rho_b
		rhW=rho_b
		drdt=rho_b
	  endif		  


	ELSE
		Unew(1:imax,1:jmax,1:kmax)=dUdt(1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)
		Vnew(1:imax,1:jmax,1:kmax)=dVdt(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
		Wnew(1:imax,1:jmax,1:kmax)=dWdt(1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)
	ENDIF
	  
	  
	  if (transporteq_fracs.eq.'massfrac') then
	    rhoU=dudt
		rhoV=dvdt
		rhoW=dwdt
	  endif
	  
      end


	  subroutine correc2(uu,vv,ww,ddtt)
      USE nlist
      implicit none
c
!       include 'param.txt'
!       include 'common.txt'
      real*8 pplus(imax,kmax)
      integer t,n
	  real uu(0:i1,0:j1,0:k1),vv(0:i1,0:j1,0:k1),ww(0:i1,0:j1,0:k1)
	  real ddtt
	  real sum_c_ws,ws(1:nfrac),ctot 


      call pshiftb(p,pplus) !,rank,imax,jmax,kmax,px)

      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          uu(i,j,k) = uu(i,j,k) -
     +      ddtt * ( p(i+1,j,k) - p(i,j,k) ) /
     +                ( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
	  !!! do nothing on imax and i=1, so p(imax+1)=p(imax) --> dp/dn=0 at outflow and inflow 
	  if (periodicx.eq.1) then
	  !!! !periodic bc in x: 
	    i = imax
	      do k=1,kmax
		    do j=1,jmax
		      uu(i,j,k) = uu(i,j,k) -ddtt * ( p(1,j,k) - p(i,j,k) ) /( Rp(i+1) - Rp(i) ) 
		    enddo
	      enddo
	  endif

      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          vv(i,j,k) = vv(i,j,k) -
     +     ddtt * ( p(i,j+1,k) - p(i,j,k) ) /
     +                ( Rp(i) * (phip(j+1)-phip(j)) )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
	      do k=1,kmax
		do i=1,imax
		vv(i,jmax,k) = vv(i,jmax,k) -
     +   ddtt * ( pplus(i,k) - p(i,jmax,k) ) /
     +              ( Rp(i) * (phip(jmax+1)-phip(jmax)) )
		enddo
	      enddo
	else
	      if (rank.ne.px-1) then
 	       do k=1,kmax
	        do i=1,imax
	        vv(i,jmax,k) = vv(i,jmax,k) -
     +   ddtt * ( pplus(i,k) - p(i,jmax,k) ) /
     +              ( Rp(i) * (phip(jmax+1)-phip(jmax)) )
	        enddo
	       enddo
      	      endif
	endif

      do k=1,kmax-1 ! rij kmax wordt niet meegenomen, dus is p(k+1)=p(k), dus dp/dn=0 bij wateroppervlak
        do j=1,jmax
          do i=1,imax
          ww(i,j,k) = ww(i,j,k) -
     +     ddtt * ( p(i,j,k+1) - p(i,j,k) ) / dz
          enddo
        enddo
      enddo
	  
	IF (continuity_solver.eq.34.or.continuity_solver.eq.36) THEN 
	 IF (slipvel.eq.2) THEN
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j),kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt
		    sum_c_ws=0.
			do n=1,nfrac
				ws(n)=-frac(n)%ws !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)+sum_c_ws  !go from velocity of volume centre to mixture velocity (centre of mass velocity)
          enddo
        enddo
      enddo
	 ELSE
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j),kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 
		    sum_c_ws=0.
		    ctot=0.
		    do n=1,nfrac
				ctot=cW(n,i,j,k)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,k+1))+ctot
			enddo
			ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
			do n=1,nfrac
				ws(n)=-frac(n)%ws*(1.-ctot)**(frac(n)%n) !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			dWdt(i,j,k)=dWdt(i,j,k)+sum_c_ws  !go from velocity of volume centre to mixture velocity (centre of mass velocity)
          enddo
        enddo
      enddo
	 ENDIF
	ENDIF 
	 IF (continuity_solver.eq.35.or.continuity_solver.eq.36) THEN 
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
		   dUdt(i,j,k)=dUdt(i,j,k)-dt*(1./rhU(i,j,k)-1./rho_b2)*(dp(i+1,j,k)-dp(i,j,k))/(Rp(i+1)-Rp(i))
		   dVdt(i,j,k)=dVdt(i,j,k)-dt*(1./rhV(i,j,k)-1./rho_b2)*(dp(i,j+1,k)-dp(i,j,k))/(Rp(i)*(phip(j+1)-phip(j)) ) 
		   dWdt(i,j,k)=dWdt(i,j,k)-dt*(1./rhW(i,j,k)-1./rho_b2)*(dp(i,j,k+1)-dp(i,j,k))/ dz		
		   dUdt(i,j,k)=dUdt(i,j,k)*rhU(i,j,k)
		   dVdt(i,j,k)=dVdt(i,j,k)*rhV(i,j,k)
		   dWdt(i,j,k)=dWdt(i,j,k)*rhW(i,j,k)
          enddo
        enddo
      enddo	 
	 ENDIF 
	 IF (continuity_solver.eq.33.or.continuity_solver.eq.34) THEN 
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
		   dUdt(i,j,k)=dUdt(i,j,k)*rhU(i,j,k)
		   dVdt(i,j,k)=dVdt(i,j,k)*rhV(i,j,k)
		   dWdt(i,j,k)=dWdt(i,j,k)*rhW(i,j,k)
          enddo
        enddo
      enddo	 
	 ENDIF 	 
	 
      end


 
!      SUBROUTINE SOLVEpois_vg_init_mumps
!!		initialises 2D poisson solver matrix structure A in cylindrical coordinates 
!!		for grid with variable grid size in two directions and fixed grid size (FFT) in one direction  
!!		Lynyrd de Wit September 2015
!	  
!      USE nlist
!      USE sparse_matrix
!	  
!      implicit none
!
!      include 'mpif.h'
!! not needed, because USE sparse_matrix      INCLUDE 'dmumps_struc.h'	  
!
!      ! Local performance variables
!      INTEGER          :: dmumps_flag,dmumps_niter
!      DOUBLE PRECISION :: dmumps_runtime,dmumps_residual,dmumps_res0,dmumps_res,t1,t2,t3
!	  
!	  
!      real      ar(imax),br(imax),cr(imax)
!      real      zrt(kmax),yrt(px*jmax)
!      real      aphi(jmax*px),bphi(jmax*px),cphi(jmax*px)
!      integer   tel,tel2,tel3,tel4,ind
!      real      dzi
!      integer   knd,nnz
!	  integer   newcom,newgroup,group,ierr
!	  
!	  
!      pi = 4.*atan(1.)
!	  
!	  !   tridiagonal system in r-direction:
!      do i=1,imax
!        ar(i)= Ru(I-1)/((Rp(I)-Rp(I-1))*Rp(I)*(Ru(I)-Ru(I-1)))
!	    cr(i)= Ru(I) /((Rp(I+1)-Rp(I))*Rp(I)*(Ru(I)-Ru(I-1)))
!!        br(i)=-(Ru(I)/(Rp(I+1)-Rp(I))+Ru(I-1)/(Rp(I)-Rp(I-1)))/
!!     $      (Rp(I)*(Ru(I)-Ru(I-1)))
!	    br(i)=-ar(i)-cr(i)
!      end do
!! DON'T SWITCH ON:		br(imax)=br(imax)+cr(imax) !! p(imax+1)= p(imax); dpdn=0
!		br(imax)=br(imax)-cr(imax) !! p(imax+1)=-p(imax); p=0 at full bc; don't do this, make p=0 at one location at outflow bc
!		br(1)=br(1)+ar(1) !dpdn=0
!		cr(imax)=0. !no interaction
!		ar(1)=0. !no interaction	  
!!       if (periodicx.eq.0) then !--> periodic need to be added later
!	  
!!  tridiagonal system in phi-direction
!      do j=1,jmax*px
!        aphi(j)=1. / ( (phipt(j  )-phipt(j-1))*(phivt(j  )-phivt(j-1)) )
!        cphi(j)=1. / ( (phipt(j+1)-phipt(j  ))*(phivt(j  )-phivt(j-1)) )
!        bphi(j)=-aphi(j)-cphi(j)
!      enddo
!      bphi(jmax*px)=bphi(jmax*px)+cphi(jmax*px) !dpdn=0
!      bphi(1)=bphi(1)+aphi(1) !dpdn=0
!      cphi(jmax*px)=0. !no interaction
!      aphi(1)=0. !no interaction
!	!!! to be added: periodic boundaries in j-dir for variable grid in j-dir
!	
!	
!      !  K --> direction      (zrt)
!      dzi = 1./dz
!      do k=1,kmax
!      zrt(k)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax)))**2
!      enddo
!	 ! build sparse matrix for 2D poisson r-phi plane
!      do k=1,kmax/px
!       knd = (rank)*kmax/px+k
!       tel=0
!       beg(1)=1
!       tel2=1
!       tel3=0	  
!	   tel4=0
!       do j=1,jmax*px
!        do i=1,imax
!         ind=(j-1)*imax+i !diagonal location
!         tel3=tel3+1
!         if ((ind-imax).gt.0.and.aphi(j).gt.0.) then
!            tel=tel+1
!            lhs(tel,k)=aphi(j)*1./(Rp(i)**2) !imax left from diagonal 
!            jco(tel)=ind-imax
!            iro(tel)=tel3
!         endif
!         if ((ind-1).gt.0.and.ar(i).gt.0.) then
!            tel=tel+1
!            lhs(tel,k)=ar(i) !one left from diagonal 
!            jco(tel)=ind-1
!            iro(tel)=tel3
!         endif
!         tel=tel+1
!         lhs(tel,k)=(br(i)+bphi(j)*1./(Rp(i)**2)+zrt(knd)) ! diagonal entry
!         jco(tel)=ind
!		 tel4=tel4+1
!		 di(tel4)=tel
!         iro(tel)=tel3
!         if ((ind+1).le.imax*jmax*px.and.cr(i).gt.0.) then
!            tel=tel+1
!            lhs(tel,k)=cr(i) !one right from diagonal 
!            jco(tel)=ind+1
!            iro(tel)=tel3
!         endif
!         if ((ind+imax).le.imax*jmax*px.and.cphi(j).gt.0.) then
!            tel=tel+1
!            lhs(tel,k)=cphi(j)*1./(Rp(i)**2) !imax right from diagonal 
!            jco(tel)=ind+imax
!            iro(tel)=tel3
!         endif
!         tel2=tel2+1
!         beg(tel2)=tel+1
!        enddo
!       enddo
!      enddo
!
!      ALLOCATE(id2(kmax/px))
!      ! Define communicator
!      CALL CPU_TIME(t1)
!      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,group,ierr) !gives group containing all ranks in mpi_comm_world
!      CALL MPI_GROUP_INCL(group,1,rank,newgroup,ierr) !gives newgroup containing only present rank
!      CALL MPI_COMM_CREATE(MPI_COMM_WORLD,newgroup,newcom,ierr) !gives communicator belonging to newgroup is only present rank
!!      id%COMM = newcom
!      ! MUMPS symmetric matrix
!!      id%SYM=0
!      ! Host working
!!      id%PAR = 1
!      ! Initialize an instance of the package
!!      id%JOB = -1
!!      CALL DMUMPS(id)
!      ! Initialize matrix structure
!!!      id%ICNTL(18)=3 !3=local input of matrix
!!      id%N=imax*jmax*px !nnodeg
!!      id%NZ=imax*jmax*px*5-2*imax-2*jmax*px !nmat
!!		ALLOCATE( id%IRN ( id%NZ ) )
!!		ALLOCATE( id%JCN ( id%NZ ) )
!!		ALLOCATE( id%A( id%NZ ) )
!!		ALLOCATE( id%RHS ( id%N ) )
!!      id%IRN=iro
!!      id%JCN=jco
!!      id%JOB=1
!
!      DO k=1,kmax/px
!	    id2(k)%COMM = newcom
!        ! MUMPS unsymmetric matrix
!        id2(k)%SYM=0 !0 !-->with SYM=2 then off-diags of A must be *0.5!
!        ! Host working
!        id2(k)%PAR = 1
!        ! Initialize an instance of the package
!        id2(k)%JOB = -1
!        CALL DMUMPS(id2(k))
!        ! Initialize matrix structure
!!       id%ICNTL(18)=3 !3=local input of matrix
!        id2(k)%N=imax*jmax*px !nnodeg
!        nnz=imax*jmax*px*5-2*imax-2*jmax*px !imax*jmax*px*3-1*imax-1*jmax*px
!        id2(k)%NZ=nnz  !nmat
!!        ALLOCATE( id2(k)%IRN ( id%NZ ) )
!!        ALLOCATE( id2(k)%JCN ( id%NZ ) )
!!        ALLOCATE( id2(k)%A( id%NZ ) )
!!        ALLOCATE( id2(k)%RHS ( id%N ) )
!        id2(k)%IRN=>iro(1:nnz)
!        id2(k)%JCN=>jco(1:nnz)
!        !To suppress a lot of output information
!        id2(k)%ICNTL(3)=0
!        id2(k)%ICNTL(4)=1
!        id2(k)%ICNTL(14)=40
!		!! perfomance tweakers:
!		!id2(k)%ICNTL(6)=7 !7 is standard
!		
!		! assign lhs:
!        id2(k)%A=>lhs(1:nnz,k)		
!        id2(k)%JOB = 4
!        CALL DMUMPS(id2(k))			  
!      enddo
!      !To suppress a lot of output information
!!      id%ICNTL(3)=0
!!      id%ICNTL(4)=1
!!      id%ICNTL(14)=40
!!      CALL CPU_TIME(t2)
!!      CALL DMUMPS(id)
!!      CALL CPU_TIME(t3)
!!      write(*,'(A,2F18.8,4I)') 't2-t1, t3-t2 = ',t2-t1,t3-t2,rank,id%MYID,id%INFOG(32),id%INFOG(7)
!	  
!
!
!      END SUBROUTINE SOLVEpois_vg_init_mumps
	  
 
!      SUBROUTINE SOLVEpois_vg_mumps(rhs)
!!		poisson solver in cylindrical coordinates 
!!		for grid with variable grid size in two directions and fixed grid size (FFT) in one direction
!!		Lynyrd de Wit, September 2015
!
!      USE nlist
!      USE sparse_matrix
!	  
!      implicit none
!
!      include 'mpif.h'
!
!!       include   'param.txt'
!      REAL      RHS(IMAX,JMAX,KMAX)
!      real      dzi
!      real      vfftk(imax*jmax,kmax),vfftj(imax*kmax/px,jmax*px)
!      real wj(4*px*jmax+15),wk(4*kmax+15),bb(imax),rtmp(imax,jmax*px,kmax/px)
!      integer   ipos,ierr,ind
!      real	rhs_ref
!      REAL, TARGET :: rhs_vec(imax*jmax*px,kmax/px)
!!	  real dumm(1,kmax),dumm2(1,kmax),dumm3(1,jmax*px),dumm4(1,jmax*px)
!!	  real p0_vec(imax*jmax*px,kmax/px),pprev(imax,jmax,kmax),pprev2(imax,jmax*px,kmax/px),vfftkp(imax*jmax,kmax)
!	  real time0,time1
!		
!	!   set up lookup tables
!      call vcosqi(kmax,wk)
!      do j=1,jmax
!         do i=1,imax
!         ipos = (j-1)*imax + i
!            do k=1,kmax
!               vfftk(ipos,k)=rhs(i,j,k)
!            enddo
!         enddo
!      enddo
!
!      call vcosqb(imax*jmax,kmax,vfftk,rhs,imax*jmax,wk)
!
!      do j=1,jmax
!         do i=1,imax
!         ipos = (j-1)*imax + i
!            do k=1,kmax
!               rhs(i,j,k)=vfftk(ipos,k)
!            enddo
!         enddo
!      enddo
!
!      call t2np(rhs,rtmp) !rtmp contains rhs shifted parallel in k-dir, size (i=1:imax,j=1:jmax*px,k=1:kmax/px)
!	  
!      do j=1,jmax*px
!       do i=1,imax
!        do k=1,kmax/px
!         ind = (j-1)*imax+i
!         rhs_vec(ind,k) = rtmp(i,j,k)
!        enddo
!       enddo
!      enddo
!
!      do k=1,kmax/px
!	  
!!	    rhs_vec(:,k)=0.
!!		do i=1,1
!!		  do j=1,jmax*px
!!		    ind = (j-1)*imax+i
!!			rhs_vec(ind,k)=1.15
!!		  enddo
!!		enddo
!		
!		! solve system lhs(:,k)*P=rhs_vec(:,k)
!      ! Assign left-hand side
!      !CALL CPU_TIME(time0)	  
!!      id%A=lhs(:,k)		
!!      id%RHS=rhs_vec(:,k)
!!      id%JOB = 5
!!      CALL DMUMPS(id)	
!      id2(k)%RHS=>rhs_vec(:,k)
!      id2(k)%JOB = 3 !5
!      CALL DMUMPS(id2(k))		  
!!	  	   write(*,*),'rank,k,sol',rank,k,id2(k).RHS
!      !CALL CPU_TIME(time1)
!       do j=1,jmax*px
!        do i=1,imax
!         ind = (j-1)*imax+i
!!         rtmp(i,j,k)=id%RHS(ind) ! sol(ind)
!		 rtmp(i,j,k)=id2(k)%RHS(ind) ! sol(ind)
!        enddo
!       enddo
!      enddo
!	  
!      call t2fp(rtmp,rhs)  !rhs (1:imax,1:jmax,1:kmax) contains solution P
!
!      do j=1,jmax
!         do i=1,imax
!         ipos = (j-1)*imax + i
!            do k=1,kmax
!               vfftk(ipos,k)=rhs(i,j,k)
!            enddo
!         enddo
!      enddo
!      call vcosqf(imax*jmax,kmax,vfftk,rhs,imax*jmax,wk)
!      do j=1,jmax
!         do i=1,imax
!         ipos = (j-1)*imax + i
!            do k=1,kmax
!               rhs(i,j,k)=vfftk(ipos,k)
!            enddo
!         enddo
!      enddo
!
!	!!    Make pressure zero at one point in outflow:
!	!!    Only b(imax)-c(imax) in matrix is not sufficient to make pressure exactly zero in (imax,1,1)
!	!!    Some small drift is occuring without following lines:
!		if (rank.eq.0) then
!		  rhs_ref=rhs(imax,1,1)
!		endif
!		call mpi_bcast(rhs_ref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!		rhs=rhs-rhs_ref
!
!      return
!      end
	  
	  

	  
      SUBROUTINE SOLVEpois_vg_init_pardiso
	  
      USE nlist	  
	  
!      USE mkl_pardiso
	  
        IMPLICIT NONE
!        include 'mkl_pardiso.f77'
C.. Internal solver memory pointer for 64-bit architectures
C.. INTEGER*8 pt(64)
C.. Internal solver memory pointer for 32-bit architectures
C.. INTEGER*4 pt(64)
C.. This is OK in both cases

!      INCLUDE 'mkl_pardiso.f90'
!      INCLUDE 'mkl_pardiso.h'	  

!        INTEGER*8 pt(64)
C.. All other variables
        INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
        !INTEGER iparm(64)
!        INTEGER ia(imax*jmax*px+1)
!        INTEGER ja(imax*jmax*px*3-1*imax-1*jmax*px)
!        REAL*8 a(imax*jmax*px*3-1*imax-1*jmax*px)
!!        REAL*8 b(8)
!!        REAL*8 x(8)

		
        INTEGER idum(1)
        REAL*8 ddum(1)

      real      ar(imax),br(imax),cr(imax),aar(imax),ccr(imax)
      real      zrt(kmax),yrt(px*jmax)
      real      aphi(jmax*px),bphi(jmax*px),cphi(jmax*px),aaphi(jmax*px),ccphi(jmax*px)
      real      time1,time0
      integer   tel,tel2,tel3,tel4,ind
      real      dzi
      integer   knd,nnz,rank_target
	  real      bbrr(imax,jmax*px,kmax/px)
	  
	  
      pi = 4.*atan(1.)
	  
	  !! TEST 9-5-2018 FOUND THAT PERIODIC (X AND Y) SIMS DO RUN WITH PARDISO AND RESULTS LOOK FINE, HOWEVER THE CONTINUITY ERROR IS E-3 INSTEAD OF E-10 WITH FFT-P-SOLVER
	  !! HAVE TO LOOK INTO THIS FURTHER IN CASE PERIODIC SIMS WITH PARDISO ARE NEEDED
	  
	  !   tridiagonal system in r-direction:
      do i=1,imax
        ar(i)= Ru(I-1)/((Rp(I)-Rp(I-1))*Rp(I)*(Ru(I)-Ru(I-1)))
	    cr(i)= Ru(I) /((Rp(I+1)-Rp(I))*Rp(I)*(Ru(I)-Ru(I-1)))
!        br(i)=-(Ru(I)/(Rp(I+1)-Rp(I))+Ru(I-1)/(Rp(I)-Rp(I-1)))/
!     $      (Rp(I)*(Ru(I)-Ru(I-1)))
	    br(i)=-ar(i)-cr(i)
      end do
! DON'T SWITCH ON:		br(imax)=br(imax)+cr(imax) !! p(imax+1)= p(imax); dpdn=0
		aar(1:imax)=ar(1:imax)
		ccr(1:imax)=cr(1:imax)
		if((U_b.ge.0.and.LOA<0).or.((U_TSHD-U_b).ge.0.and.LOA>0.)) THEN 
			if (periodicx.eq.0.or.periodicx.eq.2) then
			  if (Poutflow.eq.1) then !Poutflow zero at just 1 location
				br(imax)=br(imax)+cr(imax) !! p(imax+1)=p(imax); dpdn=0 at full bc; 
			  else 
				br(imax)=br(imax)-cr(imax) !! p(imax+1)=-p(imax); p=0 at full bc; 
			  endif 
				br(1)=br(1)+ar(1) !dpdn=0
				cr(imax)=0. !no interaction
				ar(1)=0. !no interaction	
			else ! periodicx=1
				cr(imax)=0. !no interaction
				ar(1)=0. !no interaction				
			endif	
			do i=1,imax
			  do j=1,jmax*px 
				do k=1,kmax/px 
				  bbrr(i,j,k)=br(i)
				enddo 
			  enddo 
			enddo 
			if (Poutflow.eq.1) then !Poutflow zero at just 1 location
			 rank_target = FLOOR(REAL(k_pzero)/REAL(kmax/px))
			 if (rank.eq.rank_target.and.(periodicx.eq.0.or.periodicx.eq.2)) then 
			  bbrr(imax,1,k_pzero-rank_target*kmax/px)=bbrr(imax,1,k_pzero-rank_target*kmax/px)-2*ccr(imax) !get p=0 at only one outflow gridpoint (first b+c; now -2c to arrive at b-c)
			 endif 
			endif
			rank_target = FLOOR(REAL(k_pzero)/REAL(kmax/px))
			if (rank.eq.rank_target.and.periodicx.eq.1) then !make P zero at one location (for solve_tridiag_periodic2 this is necessary otherwise that crashes, so would be natural to also for pardiso make P zero at one location) !17-10-2021 added
			  bbrr(imax,1,k_pzero-rank_target*kmax/px)=bbrr(imax,1,k_pzero-rank_target*kmax/px)-ccr(imax) !get p=0 at only one outflow gridpoint (b-c to give P=0)
			endif
		else !negative flow direction with inflow at i=imax and outflow at i=0
			if (periodicx.eq.0.or.periodicx.eq.2) then
			  if (Poutflow.eq.1) then !Poutflow zero at just 1 location
				br(1)=br(1)+ar(1) !! p(0)=p(1); dpdn=0 at full bc; 
			  else 
				br(1)=br(1)-ar(1) !! p(0)=-p(1); p=0 at full bc; 
			  endif 
				br(imax)=br(imax)+cr(imax) !dpdn=0
				cr(imax)=0. !no interaction
				ar(1)=0. !no interaction	
			else ! periodicx=1
				cr(imax)=0. !no interaction
				ar(1)=0. !no interaction				
			endif	
			do i=1,imax
			  do j=1,jmax*px 
				do k=1,kmax/px 
				  bbrr(i,j,k)=br(i)
				enddo 
			  enddo 
			enddo 
			if (Poutflow.eq.1) then !Poutflow zero at just 1 location
			 rank_target = FLOOR(REAL(k_pzero)/REAL(kmax/px))
			 if (rank.eq.rank_target.and.(periodicx.eq.0.or.periodicx.eq.2)) then 
			  bbrr(1,1,k_pzero-rank_target*kmax/px)=bbrr(1,1,k_pzero-rank_target*kmax/px)-2*aar(1) !get p=0 at only one outflow gridpoint (first b+a; now -2a to arrive at b-a)
			 endif 
			endif
			rank_target = FLOOR(REAL(k_pzero)/REAL(kmax/px))
			if (rank.eq.rank_target.and.periodicx.eq.1) then !make P zero at one location (for solve_tridiag_periodic2 this is necessary otherwise that crashes, so would be natural to also for pardiso make P zero at one location) !17-10-2021 added
			  bbrr(1,1,k_pzero-rank_target*kmax/px)=bbrr(1,1,k_pzero-rank_target*kmax/px)-aar(1) !get p=0 at only one outflow gridpoint (b-a to give P=0)
			endif 		
		endif 
	  
!  tridiagonal system in phi-direction
      do j=1,jmax*px
        aphi(j)=1. / ( (phipt(j  )-phipt(j-1))*(phivt(j  )-phivt(j-1)) )
        cphi(j)=1. / ( (phipt(j+1)-phipt(j  ))*(phivt(j  )-phivt(j-1)) )
        bphi(j)=-aphi(j)-cphi(j)
      enddo
	  aaphi(1:jmax*px)=aphi(1:jmax*px)
	  ccphi(1:jmax*px)=cphi(1:jmax*px)
	  if (periodicy.eq.0.or.periodicy.eq.2) then
		bphi(jmax*px)=bphi(jmax*px)+cphi(jmax*px) !dpdn=0
		bphi(1)=bphi(1)+aphi(1) !dpdn=0
		cphi(jmax*px)=0. !no interaction
		aphi(1)=0. !no interaction
	  else ! periodicy=1 
		cphi(jmax*px)=0. !no interaction
		aphi(1)=0. !no interaction	  
	  endif ! periodicy=1 
	
      !  K --> direction      (zrt)
      dzi = 1./dz
      do k=1,kmax
      zrt(k)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax)))**2
      enddo
	 ! build sparse matrix for 2D poisson r-phi plane
      do k=1,kmax/px
       knd = (rank)*kmax/px+k
       tel=0
       beg(1)=1
       tel2=1
       tel3=0	  
	   tel4=0
       do j=1,jmax*px
        do i=1,imax
         ind=(j-1)*imax+i !diagonal location
         tel3=tel3+1
         if (periodicy.eq.1.and.(ind-(jmax*px-1)*imax).gt.0.and.ccphi(j).gt.1.e-12) then
            tel=tel+1
            lhs(tel,k)=ccphi(j)*1./(Rp(i)**2) !(jmax*px-1)*imax left from diagonal
            jco(tel)=ind-(jmax*px-1)*imax
            iro(tel)=tel3
         endif		 
         if ((ind-imax).gt.0.and.aphi(j).gt.1.e-12) then
            tel=tel+1
            lhs(tel,k)=aphi(j)*1./(Rp(i)**2) !imax left from diagonal
            jco(tel)=ind-imax
            iro(tel)=tel3
         endif
         if (periodicx.eq.1.and.i.eq.imax.and.(ind+1-imax).gt.0.and.ccr(i).gt.1.e-12) then
            tel=tel+1
            lhs(tel,k)=ccr(i) !imax minus one left from diagonal
            jco(tel)=ind+1-imax
            iro(tel)=tel3
         endif		 
         if ((ind-1).gt.0.and.ar(i).gt.1.e-12) then
            tel=tel+1
            lhs(tel,k)=ar(i) !one left from diagonal
            jco(tel)=ind-1
            iro(tel)=tel3
         endif
         tel=tel+1
         lhs(tel,k)=(bbrr(i,j,k)+bphi(j)*1./(Rp(i)**2)+zrt(knd)) ! diagonal entry for x,y,z component (z is from fft)
         jco(tel)=ind
		 tel4=tel4+1
		 di(tel4)=tel
         iro(tel)=tel3
         if ((ind+1).le.imax*jmax*px.and.cr(i).gt.1.e-12) then
            tel=tel+1
            lhs(tel,k)=cr(i) !one right from diagonal
            jco(tel)=ind+1
            iro(tel)=tel3
         endif
         if (periodicx.eq.1.and.i.eq.1.and.(ind-1+imax).le.imax*jmax*px.and.aar(i).gt.1.e-12) then
            tel=tel+1
            lhs(tel,k)=aar(i) !imax minus one right from diagonal
            jco(tel)=ind-1+imax
            iro(tel)=tel3
         endif
         if ((ind+imax).le.imax*jmax*px.and.cphi(j).gt.1.e-12) then
            tel=tel+1
            lhs(tel,k)=cphi(j)*1./(Rp(i)**2) !imax right from diagonal
            jco(tel)=ind+imax
            iro(tel)=tel3
         endif
         if (periodicy.eq.1.and.(ind+(jmax*px-1)*imax).le.imax*jmax*px.and.aaphi(j).gt.1.e-12) then
            tel=tel+1
            lhs(tel,k)=aaphi(j)*1./(Rp(i)**2) !(jmax*px-1)*imax right from diagonal
            jco(tel)=ind+(jmax*px-1)*imax
            iro(tel)=tel3
         endif		 		 
         tel2=tel2+1
         beg(tel2)=tel+1
        enddo
       enddo
      enddo

	 
C.. Fill all arrays containing matrix data.
        n = imax*jmax*px
        nrhs =1
		maxfct=1
		mnum = 1
		if (periodicx.eq.1.and.periodicy.eq.1) then
			nnz = imax*jmax*px*5
!			if (rank.eq.0) then
!				lhs(1,1)=(1.+bphi(1)*1./(Rp(1)**2)+zrt(1))
!				lhs(2,1)=0. !only diagonal entry for first cell must be non-zero to enforce p=0 at one grid point in mesh
!				lhs(3,1)=0. ! make 4 off-diags zero
!				lhs(4,1)=0.						
!				lhs(5,1)=0.			
!			endif
		elseif (periodicx.eq.1) then
			nnz = imax*jmax*px*5-2*imax
!			if (rank.eq.0) then
!				lhs(1,1)=(1.+bphi(1)*1./(Rp(1)**2)+zrt(1))
!				lhs(2,1)=0. !only diagonal entry for first cell must be non-zero to enforce p=0 at one grid point in mesh
!				lhs(3,1)=0. ! make 3 off-diags zero
!				lhs(4,1)=0.
!			endif
		elseif (periodicy.eq.1) then
			nnz = imax*jmax*px*5-2*jmax*px 
		else
			nnz = imax*jmax*px*5-2*imax-2*jmax*px !imax*jmax*px*3-1*imax-1*jmax*px
		endif

!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.
!		ALLOCATE (pt(64))
!		DO i = 1, 64
!		   pt(i)%DUMMY =  0 
!		END DO
		if (.not. allocated(pt)) allocate(pt(64*kmax/px))
!        ALLOCATE(pt(64*kmax/px))
        DO i = 1, 64*kmax/px
            pt(i) = 0
        END DO
	 
      CALL mkl_set_num_threads(1) 
	  !CALL mkl_disable_fast_mm() !does not help with memory issues, but does make pardiso slower 
	  DO k=1,kmax/px
	    call pardisoinit (pt((k-1)*64+1:(k-1)*64+64), mtype, iparm)
	  ENDDO 
	 
C..
C.. Set up PARDISO control parameter
C..
        DO i = 1, 64
            iparm(i) = 0
        END DO
        iparm(1) = 0 ! solver defaults iparm(2)-iparm(64)
        iparm(2) = 2 ! fill-in reordering from METIS
        iparm(4) = 0 ! no iterative-direct algorithm
        iparm(5) = 0 ! no user fill-in reducing permutation
        iparm(6) = 0 ! =0 solution on the first n components of x
        iparm(8) = 2 ! numbers of iterative refinement steps
        iparm(10) = 13 !13 ! perturb the pivot elements with 1E-13; 8 default for symmetric matrix according to https://software.intel.com/en-us/node/470298 !3-9-2020 switch to 13 for k_pzero>1
        iparm(11) = 1 !0 !1 ! use nonsymmetric permutation and scaling MPS; 0 default for symmetric matrix according to https://software.intel.com/en-us/node/470298 !3-9-2020 switch to 1 for k_pzero>1
        iparm(13) = 1 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
        iparm(14) = 0 ! Output: number of perturbed pivots
        iparm(18) = -1 ! Output: number of nonzeros in the factor LU
        iparm(19) = -1 ! Output: Mflops for LU factorization
        iparm(20) = 0 ! Output: Numbers of CG Iterations
        error = 0 ! initialize error flag
        msglvl = 0 !1 ! print statistical information
        mtype = 11 !1 !11 !1=Real and structurally symmetric; 11=Real and nonsymmetric matrix both work; 1 uses 10% less memory and is 10% faster for large problems !3-9-2020 switch to 11 nonsymmetric for k_pzero>1

		

		
		
        DO k=1,kmax/px
C.. Reordering and Symbolic Factorization, This step also allocates
C all memory that is necessary for the factorization
		CALL CPU_TIME(time0)
        phase = 11 ! only reordering and symbolic factorization
        CALL pardiso (pt((k-1)*64+1:(k-1)*64+64), maxfct, mnum, mtype, phase, n, lhs(1:nnz,k), beg, jco(1:nnz),
     &  idum, nrhs, iparm, msglvl, ddum, ddum, error)
        !WRITE(*,*) 'Reordering completed ... '
        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF
        !WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
        !WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
C.. Factorization.
        phase = 22 ! only factorization
        CALL pardiso (pt((k-1)*64+1:(k-1)*64+64), maxfct, mnum, mtype, phase, n, lhs(1:nnz,k), beg, jco(1:nnz),
     &  idum, nrhs, iparm, msglvl, ddum, ddum, error)
        !WRITE(*,*) 'Factorization completed ... '
        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF
		CALL CPU_TIME(time1)
		!write(*,*),'rank,CPU wall clock time reorder and factorize=',time1-time0
C.. Termination and release of memory
!        phase = -1 ! release internal memory
!        CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, 
!     &  idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)
	 
        ENDDO
      END SUBROUTINE SOLVEpois_vg_init_pardiso
	  

      SUBROUTINE SOLVEpois_vg_pardiso(rhs)
!		poisson solver in cylindrical coordinates 
!		for grid with variable grid size in two directions and fixed grid size (FFT) in one direction
!		Lynyrd de Wit, September 2015

      USE nlist
!	  USE fftpk90
      !USE sparse_matrix
	  
      implicit none

      include 'mpif.h'

!       include   'param.txt'
      REAL      RHS(IMAX,JMAX,KMAX)
      real      dzi
      real      vfftk(imax*jmax,kmax),vfftj(imax*kmax/px,jmax*px)
      real wj(4*px*jmax+15),wk(4*kmax+15),bb(imax),rtmp(imax,jmax*px,kmax/px)
      integer   ipos,ierr,ind,i_target 
      real	rhs_ref
      REAL, TARGET :: rhs_vec(imax*jmax*px,kmax/px)
      REAL sol(imax*jmax*px)
      INTEGER nnz
!	  real dumm(1,kmax),dumm2(1,kmax),dumm3(1,jmax*px),dumm4(1,jmax*px)
!	  real p0_vec(imax*jmax*px,kmax/px),pprev(imax,jmax,kmax),pprev2(imax,jmax*px,kmax/px),vfftkp(imax*jmax,kmax)
	  real time0,time1
!        INTEGER*8 pt(64)
C.. All other variables
        INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
        !INTEGER iparm(64)
        INTEGER idum(1)
!		real, dimension(:), allocatable :: wk2
		
		
	 
C.. Fill all arrays containing matrix data.
        n = imax*jmax*px
        nrhs =1
		maxfct=1
		mnum = 1
		if (periodicx.eq.1.and.periodicy.eq.1) then
			nnz = imax*jmax*px*5
!			if (rank.eq.0) then !enforce zero pressure in one grid cell
!			  rhs(1,1,1)=0.
!			endif
		elseif (periodicx.eq.1) then
			nnz = imax*jmax*px*5-2*imax
!			if (rank.eq.0) then !enforce zero pressure in one grid cell
!			  rhs(1,1,1)=0.
!			endif			
		elseif (periodicy.eq.1) then
			nnz = imax*jmax*px*5-2*jmax*px 
		else
			nnz = imax*jmax*px*5-2*imax-2*jmax*px !imax*jmax*px*3-1*imax-1*jmax*px
		endif

		
	 
	
!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.
!		ALLOCATE (pt(64))
!		DO i = 1, 64
!		   pt(i)%DUMMY =  0 
!		END DO
!        DO i = 1, 64
!            pt(i) = 0
!        END DO

        error = 0 ! initialize error flag
        msglvl = 0 !1 ! print statistical information
        mtype = 11 !11 !1 !11 !-2 ! symmetric, indefinite !3-9-2020 switch to 11 nonsymmetric for k_pzero>1
		
		
      do j=1,jmax
         do i=1,imax
         ipos = (j-1)*imax + i
            do k=1,kmax
               vfftk(ipos,k)=rhs(i,j,k)
            enddo
         enddo
      enddo
	  
!	  if (fft_routines.eq.1) then !good-old vfft.f 
		call vcosqi(kmax,wk) !   set up lookup tables
		call vcosqb(imax*jmax,kmax,vfftk,rhs,imax*jmax,wk)
!	  elseif (fft_routines.eq.2) then !vfftpack_netlib_Clarke 
!		allocate(wk2(ftqiq(kmax)))
!		call ftqi (wk2) !   set up lookup tables
!		call ftqcb (vfftk,wk2) !  ftqcb is similar to vcosqb
!	  endif 
		

      do j=1,jmax
         do i=1,imax
         ipos = (j-1)*imax + i
            do k=1,kmax
               rhs(i,j,k)=vfftk(ipos,k)
            enddo
         enddo
      enddo

      call t2np(rhs,rtmp) !rtmp contains rhs shifted parallel in k-dir, size (i=1:imax,j=1:jmax*px,k=1:kmax/px)
	  
      do j=1,jmax*px
       do i=1,imax
        do k=1,kmax/px
         ind = (j-1)*imax+i
         rhs_vec(ind,k) = rtmp(i,j,k)
        enddo
       enddo
      enddo

	  n = imax*jmax*px
      do k=1,kmax/px
	  
!	    rhs_vec(:,k)=0.
!		do i=1,1
!		  do j=1,jmax*px
!		    ind = (j-1)*imax+i
!			rhs_vec(ind,k)=1.15
!		  enddo
!		enddo
			
		
	        !CALL CPU_TIME(time1)
		! solve system lhs(:,k)*P=rhs_vec(:,k)
		!.. Back substitution and iterative refinement
        phase = 33 !33 !13 !analyze,factorize,solve !33 ! only solution
        CALL pardiso (pt((k-1)*64+1:(k-1)*64+64), maxfct, mnum, mtype, phase, n, lhs(1:nnz,k), beg, jco(1:nnz),
     &  idum, nrhs, iparm, msglvl, rhs_vec(:,k), sol, error)
!	   write(*,*),'rank,k,sol',rank,k,sol
	 
      !CALL CPU_TIME(time2)
       do j=1,jmax*px
        do i=1,imax
         ind = (j-1)*imax+i
		 rtmp(i,j,k)=sol(ind)
        enddo
       enddo
      enddo
	  
      call t2fp(rtmp,rhs)  !rhs (1:imax,1:jmax,1:kmax) contains solution P

      do j=1,jmax
         do i=1,imax
         ipos = (j-1)*imax + i
            do k=1,kmax
               vfftk(ipos,k)=rhs(i,j,k)
            enddo
         enddo
      enddo
!	  if (fft_routines.eq.1) then !good-old vfft.f 
		call vcosqf(imax*jmax,kmax,vfftk,rhs,imax*jmax,wk)
!	  elseif (fft_routines.eq.2) then !vfftpack_netlib_Clarke 
!		call ftqcf (vfftk,wk2) !  ftqcb is similar to vcosqb
!	  endif 	  
      
      do j=1,jmax
         do i=1,imax
         ipos = (j-1)*imax + i
            do k=1,kmax
               rhs(i,j,k)=vfftk(ipos,k)
            enddo
         enddo
      enddo

!!!	!!    Make pressure zero at one point in outflow:
!!!	!!    Only b(imax)-c(imax) in matrix is not sufficient to make pressure exactly zero in (imax,1,1)
!!!	!!    Some small drift is occuring without following lines:
!!!	    if((U_b.ge.0.and.LOA<0).or.((U_TSHD-U_b).ge.0.and.LOA>0.)) THEN 
!!!	      i_target = imax 
!!!	    else 
!!!	      i_target = 1 
!!!	    endif 
!!!		if (rank.eq.0) then
!!!		  rhs_ref=rhs(i_target,1,k_pzero)
!!!		endif
!!!		call mpi_bcast(rhs_ref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!!!		rhs=rhs-rhs_ref
!!!		! switched off 22-2-2023; should not be needed and this correction makes dpdx at imax contribution incorrect 
      return
      end
	  

	  
	  
      SUBROUTINE SOLVEpois_vg_init_pardiso3D !only rank=0 calls this subroutine
	  
      USE nlist	  
	  
!      USE mkl_pardiso
	  
        IMPLICIT NONE
!        include 'mkl_pardiso.f77'
C.. Internal solver memory pointer for 64-bit architectures
C.. INTEGER*8 pt(64)
C.. Internal solver memory pointer for 32-bit architectures
C.. INTEGER*4 pt(64)
C.. This is OK in both cases

!      INCLUDE 'mkl_pardiso.f90'

!        INTEGER*8 pt(64)
C.. All other variables
        INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
        !INTEGER iparm(64)
!        INTEGER ia(imax*jmax*px+1)
!        INTEGER ja(imax*jmax*px*3-1*imax-1*jmax*px)
!        REAL*8 a(imax*jmax*px*3-1*imax-1*jmax*px)
!!        REAL*8 b(8)
!!        REAL*8 x(8)

		
        INTEGER idum(1)
        REAL*8 ddum(1)

      real      ar(imax),br(imax),cr(imax)
      real      aphi(jmax*px),bphi(jmax*px),cphi(jmax*px)
      real      time1,time0
      real      az(kmax),bz(kmax),cz(kmax)
      integer   tel,tel2,tel3,tel4,ind
      real      dzi
      integer   knd,nnz
	  
	  
	! only rank=0 runs this subroutine once, therefore the following allocation:  
	nnz=imax*jmax*px*kmax*7-2*imax-2*jmax*px-2*kmax 
	ALLOCATE(jco3(nnz)) !col nr CSR 
	ALLOCATE(iro3(nnz)) !row nr CSR 
	ALLOCATE(di3(imax*jmax*px*kmax)) !diag nr CSR
	ALLOCATE(beg3(imax*jmax*px*kmax+1)) !begin new line CSR	  
	ALLOCATE(lhs3(nnz))
	
!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.
!		ALLOCATE (pt(64))
!		DO i = 1, 64
!		   pt(i)%DUMMY =  0 
!		END DO
        ALLOCATE(pt3(64))
        DO i = 1, 64
            pt3(i) = 0
        END DO	
	  
      pi = 4.*atan(1.)
	  
	  !   tridiagonal system in r-direction:
      do i=1,imax
        ar(i)= Ru(I-1)/((Rp(I)-Rp(I-1))*Rp(I)*(Ru(I)-Ru(I-1)))
	    cr(i)= Ru(I) /((Rp(I+1)-Rp(I))*Rp(I)*(Ru(I)-Ru(I-1)))
!        br(i)=-(Ru(I)/(Rp(I+1)-Rp(I))+Ru(I-1)/(Rp(I)-Rp(I-1)))/
!     $      (Rp(I)*(Ru(I)-Ru(I-1)))
	    br(i)=-ar(i)-cr(i)
      end do
! DON'T SWITCH ON:		br(imax)=br(imax)+cr(imax) !! p(imax+1)= p(imax); dpdn=0
		br(imax)=br(imax)-cr(imax) !! p(imax+1)=-p(imax); p=0 at full bc; don't do this, make p=0 at one location at outflow bc
		br(1)=br(1)+ar(1) !dpdn=0
		cr(imax)=0. !no interaction
		ar(1)=0. !no interaction	  
!       if (periodicx.eq.0) then !--> periodic need to be added later
	  

!  tridiagonal system in phi-direction
      do j=1,jmax*px
        aphi(j)=1. / ( (phipt(j  )-phipt(j-1))*(phivt(j  )-phivt(j-1)) )
        cphi(j)=1. / ( (phipt(j+1)-phipt(j  ))*(phivt(j  )-phivt(j-1)) )
        bphi(j)=-aphi(j)-cphi(j)
      enddo
      bphi(jmax*px)=bphi(jmax*px)+cphi(jmax*px) !dpdn=0
      bphi(1)=bphi(1)+aphi(1) !dpdn=0
      cphi(jmax*px)=0. !no interaction
      aphi(1)=0. !no interaction
	!!! to be added: periodic boundaries in j-dir for variable grid in j-dir

!  tridiagonal system in z-direction
      do k=1,kmax
        az(k)=1. / (dz**2) 
        cz(k)=1. / (dz**2) 
        bz(k)=-az(k)-cz(k)
      enddo
      bz(kmax)=bz(kmax)+cz(kmax) !dpdn=0
      bz(1)=bz(1)+az(1) !dpdn=0
      cz(kmax)=0. !no interaction
      az(1)=0. !no interaction
	
	
	 ! build sparse matrix for 3D poisson r-phi-z cube
       tel=0
       beg3(1)=1
       tel2=1
       tel3=0	  
	   tel4=0
      do k=1,kmax	   
       do j=1,jmax*px
        do i=1,imax
         ind=(k-1)*imax*jmax*px+(j-1)*imax+i !diagonal location
         tel3=tel3+1
         if ((ind-imax*jmax*px).gt.0.and.az(k).gt.0.) then
            tel=tel+1
            lhs3(tel)=az(k) !imax*jmax*px left from diagonal
            jco3(tel)=ind-imax*jmax*px
            iro3(tel)=tel3
         endif 
         if ((ind-imax).gt.0.and.aphi(j).gt.0.) then
            tel=tel+1
            lhs3(tel)=aphi(j)*1./(Rp(i)**2) !imax left from diagonal
            jco3(tel)=ind-imax
            iro3(tel)=tel3
         endif
         if ((ind-1).gt.0.and.ar(i).gt.0.) then
            tel=tel+1
            lhs3(tel)=ar(i) !one left from diagonal
            jco3(tel)=ind-1
            iro3(tel)=tel3
         endif
         tel=tel+1
         lhs3(tel)=(br(i)+bphi(j)*1./(Rp(i)**2)+bz(k)) ! diagonal entry
         jco3(tel)=ind
		 tel4=tel4+1
		 di3(tel4)=tel
         iro3(tel)=tel3
         if ((ind+1).le.imax*jmax*px*kmax.and.cr(i).gt.0.) then
            tel=tel+1
            lhs3(tel)=cr(i) !one right from diagonal
            jco3(tel)=ind+1
            iro3(tel)=tel3
         endif
         if ((ind+imax).le.imax*jmax*px*kmax.and.cphi(j).gt.0.) then
            tel=tel+1
            lhs3(tel)=cphi(j)*1./(Rp(i)**2) !imax right from diagonal
            jco3(tel)=ind+imax
            iro3(tel)=tel3
         endif
         if ((ind+imax*jmax*px).le.imax*jmax*px*kmax.and.cz(k).gt.0.) then
            tel=tel+1
            lhs3(tel)=cz(k) !imax*jmax*px right from diagonal
            jco3(tel)=ind+imax*jmax*px
            iro3(tel)=tel3
         endif 
         tel2=tel2+1
         beg3(tel2)=tel+1
        enddo
       enddo
      enddo

	 
C.. Fill all arrays containing matrix data.
        n = imax*jmax*px*kmax
        nrhs =1
		maxfct=1
		mnum = 1

      CALL mkl_set_num_threads(px) 
	 
C..
C.. Set up PARDISO control parameter
C..
        DO i = 1, 64
            iparm(i) = 0
        END DO
        iparm(1) = 0 !1 ! no solver default
        iparm(2) = 2 !2 ! fill-in reordering from METIS
        iparm(4) = 0 ! no iterative-direct algorithm
        iparm(5) = 0 ! no user fill-in reducing permutation
        iparm(6) = 0 ! =0 solution on the first n components of x
        iparm(8) = 2 ! numbers of iterative refinement steps
        iparm(10) = 8 !13 ! perturb the pivot elements with 1E-13; 8 default for symmetric matrix according to https://software.intel.com/en-us/node/470298
        iparm(11) = 1 !1 ! use nonsymmetric permutation and scaling MPS; 0 default for symmetric matrix according to https://software.intel.com/en-us/node/470298
        iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
        iparm(14) = 0 ! Output: number of perturbed pivots
        iparm(18) = -1 ! Output: number of nonzeros in the factor LU
        iparm(19) = -1 ! Output: Mflops for LU factorization
        iparm(20) = 0 ! Output: Numbers of CG Iterations
        error = 0 ! initialize error flag
        msglvl = 0 !1 ! print statistical information
        mtype = 11 !unsymmetric
		
		
		
C.. Reordering and Symbolic Factorization, This step also allocates
C all memory that is necessary for the factorization
		CALL CPU_TIME(time0)
        phase = 11 ! only reordering and symbolic factorization
        CALL pardiso (pt3, maxfct, mnum, mtype, phase, n, lhs3(1:nnz), beg3, jco3(1:nnz),
     &  idum, nrhs, iparm, msglvl, ddum, ddum, error)
        !WRITE(*,*) 'Reordering completed ... '
        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF
        !WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
        !WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
C.. Factorization.
        phase = 22 ! only factorization
        CALL pardiso (pt3, maxfct, mnum, mtype, phase, n, lhs3(1:nnz), beg3, jco3(1:nnz),
     &  idum, nrhs, iparm, msglvl, ddum, ddum, error)
        !WRITE(*,*) 'Factorization completed ... '
        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF
		CALL CPU_TIME(time1)
		write(*,*),'rank,CPU wall clock time poisson 3D reorder and factorize=',rank,time1-time0
C.. Termination and release of memory
!        phase = -1 ! release internal memory
!        CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, 
!     &  idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)
 

      END SUBROUTINE SOLVEpois_vg_init_pardiso3D
	  
	  
      SUBROUTINE SOLVEpois_vg_pardiso3D(rhs) !only rank=0 calls this subroutine; pardiso uses openmp
!		poisson solver in cylindrical coordinates 
!		for grid with variable grid size in 3 directions
!		Lynyrd de Wit, October 2015

      USE nlist
      !USE sparse_matrix
	  
      implicit none

      include 'mpif.h'

!       include   'param.txt'
      REAL      RHS(IMAX,JMAX*px,KMAX)
      real      dzi
      integer   ipos,ierr,ind
      REAL, TARGET :: rhs_vec(imax*jmax*px*kmax)
      REAL sol(imax*jmax*px*kmax),rhs_ref

      INTEGER nnz
!	  real dumm(1,kmax),dumm2(1,kmax),dumm3(1,jmax*px),dumm4(1,jmax*px)
!	  real p0_vec(imax*jmax*px,kmax/px),pprev(imax,jmax,kmax),pprev2(imax,jmax*px,kmax/px),vfftkp(imax*jmax,kmax)
	  real time0,time1
!        INTEGER*8 pt(64)
C.. All other variables
        INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
        !INTEGER iparm(64)
        INTEGER idum(1)
		
		
	 
C.. Fill all arrays containing matrix data.
        n = imax*jmax*px*kmax
        nrhs =1
		maxfct=1
		mnum = 1
		nnz = imax*jmax*px*7-2*imax-2*jmax*px-2*kmax
		n = imax*jmax*px*kmax

        error = 0 ! initialize error flag
        msglvl = 0 !1 ! print statistical information
        mtype = 11 ! unsymmetric, indefinite

		DO k=1,kmax
		  DO j=1,jmax*px
		    DO i=1,imax
			  ind=(k-1)*imax*jmax*px+(j-1)*imax+i
			  rhs_vec(ind)=rhs(i,j,k)
		    ENDDO
		  ENDDO
		ENDDO
		
	        !CALL CPU_TIME(time1)
		! solve system lhs*P=rhs_vec
        phase = 33 ! only solution
        CALL pardiso (pt3, maxfct, mnum, mtype, phase, n, lhs3(1:nnz), beg3, jco3(1:nnz),
     &  idum, nrhs, iparm, msglvl, rhs_vec, sol, error)
!	   write(*,*),'rank,k,sol',rank,k,sol

		DO k=1,kmax
		  DO j=1,jmax*px
		    DO i=1,imax
			  ind=(k-1)*imax*jmax*px+(j-1)*imax+i
			  rhs(i,j,k)=sol(ind)
		    ENDDO
		  ENDDO
		ENDDO
		
!!!	!!    Make pressure zero at one point in outflow:
!!!	!!    Only b(imax)-c(imax) in matrix is not sufficient to make pressure exactly zero in (imax,1,1)
!!!	!!    Some small drift is occuring without following lines:
!!!        rhs_ref=rhs(imax,1,kmax)
!!!		rhs=rhs-rhs_ref
!!! switched off 22-2-2023; should not be needed and this correction makes dpdx at imax contribution incorrect

      return
      end
	  

	  
      SUBROUTINE SOLVEpois(rhs) !,Ru,Rp,dphi,dz,rank,imax,jmax,kmax,px)
c
c  FAST POISSON SOLVER IN CYLINDRICAL COORDINATES
c  BENDIKS JAN BOERSMA
C  LABORATORIUM VOOR AERO EN HYDRODYNAMICA
C  ROTTERDAMSEWEG 145
C  2628 AL DELFT
C  email ::::   b.j.boersma@wbmt.tudelft.nl
*

      USE nlist
      implicit none

      include 'mpif.h'

!       include   'param.txt'
      REAL      RHS(IMAX,JMAX,KMAX)
!       REAL      RHS(IMAX,JMAX,KMAX)
!       real      dz,dzi,pi,d(IMAX,JMAX,kmax),bbb,z
      real      dzi,d(IMAX,JMAX,kmax),bbb,z,bbbb(IMAX,JMAX,kmax)
      real      a(imax),b(imax),c(imax),dd(imax)
      real      zrt(kmax),yrt(px*jmax)
      real      vfftk(imax*jmax,kmax),vfftj(imax*kmax/px,jmax*px)
      real wj(4*px*jmax+15),wk(4*kmax+15),bb(imax),rtmp(imax,jmax*px,kmax/px)
      integer   ipos,ierr,i_target 
      real	rhs_ref,dphi3
	real bbbbb(IMAX,JMAX,kmax)
      
!	  rhs=0.
!	  IF (rank.eq.px-1) THEN
!	    do i=1,imax
!		  do k=1,1 !,kmax
!	        rhs(i,jmax,k)=-1. !./(Rp(i)*dphi)
!		  enddo
!		enddo
!		endif
!		rhs(1,1:jmax,1)=1.

      dphi3=(phiv(2)-phiv(1)) !This subroutine with FFT in y-dir can only be applied with equidistant dphi
	
	  
c   generate tridiagonal systems
      pi = 4.*atan(1.)
      do i=1,imax
      a(i)= Ru(I-1)/((Rp(I)-Rp(I-1))*Rp(I)*(Ru(I)-Ru(I-1)))
      b(i)=-(Ru(I)/(Rp(I+1)-Rp(I))+Ru(I-1)/(Rp(I)-Rp(I-1)))/
     $      (Rp(I)*(Ru(I)-Ru(I-1)))
      c(i)= Ru(I) /((Rp(I+1)-Rp(I))*Rp(I)*(Ru(I)-Ru(I-1)))
      end do
	
	if((U_b.ge.0.and.LOA<0).or.((U_TSHD-U_b).ge.0.and.LOA>0.)) THEN 
		if (periodicx.eq.0.or.periodicx.eq.2) then
			  b(1)=   -(Ru(1)/(Rp(2)-Rp(1)))/(Rp(1)*(Ru(1)-Ru(0))) !dpdn=0 
			  if (Poutflow.eq.1) then !Poutflow zero at just 1 location
				b(imax)=b(imax)+c(imax)      !! p(imax+1)= p(imax) --> bc dpdn=0
			  else !default P=0 complete outflow 
				b(imax)=b(imax)-c(imax)     !! p(imax+1)=-p(imax) --> bc p=0
			  endif 
	!			  i=imax !! p(imax+1)= p(imax) --> bc dpdn=0
	 !     	      b(i)=-(Ru(I-1)/(Rp(I)-Rp(I-1)))/   !! this is equivalent to b(imax)+c(imax) --> dpdn=0
	!     $        (Rp(I)*(Ru(I)-Ru(I-1)))
		!		  b(imax)=b(imax)+c(imax)
	!	else !do nothing
		endif
		! make bbbb three dimensional to get p=0 only at one loc of outflow (not everywhere at outflow):
		do i=1,imax
		  do j=1,jmax
			do k=1,kmax
			  bbbb(i,j,k)=b(i)
			enddo
		  enddo
		enddo
		 if (periodicx.eq.0.or.periodicx.eq.2) then
		  if (Poutflow.eq.1) then !Poutflow zero at just 1 location 
		   if (rank==0) then
			 bbbb(imax,1,k_pzero)=b(imax)-2.*c(imax)     !! p(imax+1)=-p(imax) --> bc p=0 at 1 loc in mesh (first +c(imax) so now -2c(imax) to arrive at b-c)
		   endif
		  endif	   
		   c(imax)=0.
		   a(1)=0.
		 else !make P zero at one location for periodicx.eq.1 otherwise P is undetermined and solve_tridiag_periodic2 doesn't give answer
		  if (rank==0) then
			bbbb(imax,1,k_pzero)=b(imax)-c(imax)     !! p(imax+1)=-p(imax) --> bc p=0 at 1 loc in mesh 
		  endif
		 endif
	else !negative flow direction with inflow at i=imax and outflow at i=0
		if (periodicx.eq.0.or.periodicx.eq.2) then
			  b(imax) = b(imax)+c(imax) !dpdn = 0
			  if (Poutflow.eq.1) then !Poutflow zero at just 1 location
				b(1)=b(1)+a(1)      !! p(0)= p(1) --> bc dpdn=0
			  else !default P=0 complete outflow 
				b(1)=b(1)-a(1)     !! p(0)=-p(1) --> bc p=0
			  endif 
		endif
		! make bbbb three dimensional to get p=0 only at one loc of outflow (not everywhere at outflow):
		do i=1,imax
		  do j=1,jmax
			do k=1,kmax
			  bbbb(i,j,k)=b(i)
			enddo
		  enddo
		enddo
		 if (periodicx.eq.0.or.periodicx.eq.2) then
		  if (Poutflow.eq.1) then !Poutflow zero at just 1 location 
		   if (rank==0) then
			 bbbb(1,1,k_pzero)=b(1)-2.*a(1)     !! p(0)=-p(1) --> bc p=0 at 1 loc in mesh (first +a(1) so now -2a(1) to arrive at b-a)
		   endif
		  endif	   
		   c(imax)=0.
		   a(1)=0.
		 else !make P zero at one location for periodicx.eq.1 otherwise P is undetermined and solve_tridiag_periodic2 doesn't give answer
		  if (rank==0) then
			bbbb(1,1,k_pzero)=b(1)-a(1)     !! p(0)=-p(1) --> bc p=0 at 1 loc in mesh 
		  endif
		 endif
	endif 

      do i=1,imax
         dd(i) = 1/(Rp(i)**2)
      enddo
c Generate Eigenvalues
c  K --> direction      (zrt)
      dzi = 1./dz
      do k=1,kmax
      zrt(k)=-4.*dzi*dzi*(sin(float((k-1))*pi/(2.*kmax)))**2
      enddo
c  J --> direction      (yrt)
	if (periodicy.eq.0.or.periodicy.eq.2) then
	      do j=1,jmax*px
	        yrt(j)=(-4./(dphi3*dphi3))*(sin(float((j-1))*pi/(2.*jmax*px)))**2
      	      enddo 
		!   set up lookup tables
      	      call vcosqi(px*jmax,wj)
	elseif (periodicy.eq.1) then
	      do j=3,jmax*px,2
	        yrt(j-1)=(-4./(dphi3*dphi3))*(sin(float((j-1))*pi/(2.*jmax*px)))**2
		yrt(j)=yrt(j-1)
      	      enddo 
	      yrt(1)=0.
		if (mod(jmax*px,2).eq.0) then
		  yrt(jmax*px)=(-4./(dphi3*dphi3))
		endif		
		!   set up lookup tables
      	      call vrffti(px*jmax,wj)	
	endif	
	
	!   set up lookup tables
      call vcosqi(kmax,wk)
      do j=1,jmax
         do i=1,imax
         ipos = (j-1)*imax + i
            do k=1,kmax
               vfftk(ipos,k)=rhs(i,j,k)
            enddo
         enddo
      enddo
      call vcosqb(imax*jmax,kmax,vfftk,rhs,imax*jmax,wk)
      do j=1,jmax
         do i=1,imax
         ipos = (j-1)*imax + i
            do k=1,kmax
               rhs(i,j,k)=vfftk(ipos,k)
            enddo
         enddo
      enddo
  
      call t2np(rhs,rtmp)
      do k=1,kmax/px
         do i=1,imax
         ipos = (k-1)*imax+i
           do j=1,jmax*px
             vfftj(ipos,j)=rtmp(i,j,k)
           enddo
         enddo
      enddo
      if (periodicy.eq.0.or.periodicy.eq.2) then
	      call vcosqb(imax*kmax/px,jmax*px,vfftj,rhs,imax*kmax/px,wj)
      else 
	      call vrfftf(imax*kmax/px,jmax*px,vfftj,rhs,imax*kmax/px,wj) !not a typo that forward is used instead of backward for periodicy.eq.0, test 23-8-2024 shows that this is right 
      endif
      do k=1,kmax/px
         do i=1,imax
         ipos = (k-1)*imax+i
            do j=1,jmax*px
               rtmp(i,j,k) =vfftj(ipos,j)
            enddo
         enddo
      enddo
      call t2fp(rtmp,rhs)

      if (periodicx.eq.0.or.periodicx.eq.2) then !  solve tridiagonal systems with Gaussian elemination
	      do k=1,kmax
		 do j=1,jmax
	!            bbb        = b(1)+yrt(j+rank*jmax)*dd(1)+zrt(k)
		    bbb        = bbbb(1,j,k)+yrt(j+rank*jmax)*dd(1)+zrt(k)
		    z          = 1./bbb 
		    d(1,j,k)   = c(1)*z
		    rhs(1,j,k) = rhs(1,j,k)*z
		 enddo
	      enddo
	      do k=1,kmax
		do j=1,jmax
		   do i=2,imax  - 1 
	!            bb(i)      = b(i)+yrt(j+rank*jmax)*dd(i)+zrt(k)
		    bb(i)      = bbbb(i,j,k)+yrt(j+rank*jmax)*dd(i)+zrt(k)
		    z          = 1./(bb(i)-a(i)*d(i-1,j,k))
		    d(i,j,k)   = c(i)*z
		    rhs(i,j,k) = (rhs(i,j,k)-a(i)*rhs(i-1,j,k))*z
		   enddo
		enddo
	       enddo
	      do k=1,kmax
		 do j=1,jmax
	!          z = b(imax)+yrt(j+rank*jmax)*dd(imax)+zrt(k)-a(imax)*d(imax-1,j,k)
		  z = bbbb(imax,j,k)+yrt(j+rank*jmax)*dd(imax)+zrt(k)-a(imax)*d(imax-1,j,k)
		  if (z .ne. 0.) then
		  rhs(imax,j,k) = (rhs(imax,j,k)-a(imax)*rhs(imax-1,j,k))/z
		  else 
		  write(*,*),'TDMA trouble z=0,rank,j,k',rank,j,k
		  rhs(imax,j,k) = 0.
		  endif
		  enddo
		enddo
	      do k=1,kmax
		  do j=1,jmax
		     do  i=imax-1,1,-1
		          rhs(i,j,k) = rhs(i,j,k)-d(i,j,k)*rhs(i+1,j,k)
		     enddo
		   enddo
	      enddo
      else ! periodic solver x: solve tridiagonal systems with Gaussian elemination
	      do k=1,kmax
		do j=1,jmax
		  do i=1,imax
		    bbbbb(i,j,k)=bbbb(i,j,k)+yrt(j+rank*jmax)*dd(i)+zrt(k)
		  enddo
		enddo
	      enddo
	      do k=1,kmax
		 do j=1,jmax
		    CALL solve_tridiag_periodic2(rhs(1:imax,j,k),a,bbbbb(1:imax,j,k),c,rhs(1:imax,j,k),imax) 
		 enddo
	      enddo
      endif


      do j=1,jmax
         do i=1,imax
         ipos = (j-1)*imax + i
            do k=1,kmax
               vfftk(ipos,k)=rhs(i,j,k)
            enddo
         enddo
      enddo
      call vcosqf(imax*jmax,kmax,vfftk,rhs,imax*jmax,wk)
      do j=1,jmax
         do i=1,imax
         ipos = (j-1)*imax + i
            do k=1,kmax
               rhs(i,j,k)=vfftk(ipos,k)
            enddo
         enddo
      enddo
      call t2np(rhs,rtmp)
      do k=1,kmax/px
         do i=1,imax
         ipos = (k-1)*imax+i
           do j=1,jmax*px
             vfftj(ipos,j)=rtmp(i,j,k)
           enddo
         enddo
      enddo
      if (periodicy.eq.0.or.periodicy.eq.2) then
	      call vcosqf(imax*kmax/px,jmax*px,vfftj,rhs,imax*kmax/px,wj)
      else 
	      call vrfftb(imax*kmax/px,jmax*px,vfftj,rhs,imax*kmax/px,wj) !not a typo that forward is used instead of backward for periodicy.eq.0, test 23-8-2024 shows that this is right 
      endif

      do k=1,kmax/px
         do i=1,imax
         ipos = (k-1)*imax+i
            do j=1,jmax*px
               rtmp(i,j,k) =vfftj(ipos,j)
            enddo
         enddo
      enddo
      call t2fp(rtmp,rhs)
	

!!!	!!    Make pressure zero at one point in outflow:
!!!	!!    Only b(imax)-c(imax) in matrix is not sufficient to make pressure exactly zero in (imax,1,1)
!!!	!!    Some small drift is occuring without following lines:
!!!	    if((U_b.ge.0.and.LOA<0).or.((U_TSHD-U_b).ge.0.and.LOA>0.)) THEN 
!!!	      i_target = imax 
!!!	    else 
!!!	      i_target = 1 
!!!	    endif 	
!!!		if (rank.eq.0) then
!!!		  rhs_ref=rhs(i_target,1,k_pzero)
!!!		endif
!!!		call mpi_bcast(rhs_ref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
!!!		rhs=rhs-rhs_ref
!!!! switched off 22-2-2023; should not be needed and this correction makes dpdx at imax contribution incorrect


      return
      end
      
      
      subroutine state(c,rho)

      USE nlist
! 	include 'param.txt'
      real c(nfrac,0:i1,0:j1,0:k1)
      real rho(0:i1,0:j1,0:k1)
      integer n,n2,inout,ierr,tel
	real xTSHD(4),yTSHD(4),phi,xx,yy
      
      do k=0,k1
	do j=0,j1
	 do i=0,i1
	   rho(i,j,k)=rho_b
	   do n=1,nfrac

!	     rho(i,j,k)=rho(i,j,k)-c(n,i,j,k)*rho_b+c(n,i,j,k)*frac(n)%rho
	     rho(i,j,k)=rho(i,j,k)+c(n,i,j,k)*(frac(n)%rho-rho_b)

!	    rho(i,j,k)=rho_b+(rho_j-rho_b)*c(i,j,k) !rho=(1-c)*rho_w+c*rho_b=rho_w+c*(rho_b-rho_w)
	   enddo
	   do n=1,nair !correction for compressible air fraction; compressible air fills different volume not occupied by water
	     rho(i,j,k)=rho(i,j,k)-(rhocorr_air_z(nfrac_air(n),k)*c(nfrac_air(n),i,j,k)-c(nfrac_air(n),i,j,k))*rho_b
	   enddo
	 enddo
        enddo
      enddo
	  
	DO n2=1,nbedplume
	IF (bp(n2)%forever.eq.0.and.time_n.lt.bp(n2)%t0.and.time_np.gt.bp(n2)%t0) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      !! Search for P,V:
      do k=CEILING(bp(n2)%zbottom/dz),FLOOR(bp(n2)%height/dz) !do k=0,k1
	   do tel=1,bp(n2)%tmax 
		 i=bp(n2)%i(tel) 
		 j=bp(n2)%j(tel) 	  
!!       do i=0,i1  
!!         do j=j1,0,-1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
!!	  xx=Rp(i)*cos_u(j)-schuif_x
!!	  yy=Rp(i)*sin_u(j)
!!!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
!!		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
!!		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
!!!	  ELSE 
!!!	 	inout=0
!!!	  ENDIF
!!	  if (inout.eq.1) then
		rnew(i,j,k)=rho(i,j,k)
		rold(i,j,k)=rho(i,j,k)
		! prevent large source in pres-corr by sudden increase in density with bp%c 
!!	   endif
!!	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop

     
	 
	  
      end
	  
      subroutine state_edges(c,rho)

      USE nlist
      real c(nfrac,0:i1,0:j1,0:k1)
      real rho(0:i1,0:j1,0:k1)
      integer n
      
		do k=0,k1
			do j=0,j1
				do i=0,i1
					rho(i,j,k)=rho_b
					do n=1,nfrac
						rho(i,j,k)=rho(i,j,k)+c(n,i,j,k)*(frac(n)%rho-rho_b)
					enddo
					do n=1,nair !correction for compressible air fraction; compressible air fills different volume not occupied by water
						rho(i,j,k)=rho(i,j,k)-(rhocorr_air_z(nfrac_air(n),k)*c(nfrac_air(n),i,j,k)-c(nfrac_air(n),i,j,k))*rho_b
					enddo
				enddo
			enddo
		enddo
    
	end
	  

	        subroutine state_massfrac(c,rho)

      USE nlist
! 	include 'param.txt'
      real c(nfrac,0:i1,0:j1,0:k1)
      real rho(0:i1,0:j1,0:k1)
      integer n,n2,inout,ierr,tel
	real xTSHD(4),yTSHD(4),phi,xx,yy,summ
      
      do k=0,k1
	do j=0,j1
	 do i=0,i1
	   rho(i,j,k)=rho_b
	   summ=0.
	   do n=1,nfrac
		summ=summ+c(n,i,j,k)*(1.-rho_b/(frac(n)%rho/rhocorr_air_z(n,k)))
	   enddo
		rho(i,j,k)=rho_b/summ
	 enddo
        enddo
      enddo
	  
	DO n2=1,nbedplume
	IF (bp(n2)%forever.eq.0.and.time_n.lt.bp(n2)%t0.and.time_np.gt.bp(n2)%t0) THEN
	! rotation ship for ambient side current
	if (LOA<0.) then 
	  phi=0. !don't rotate grid
	else 
		if ((U_TSHD-U_b).eq.0) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
	endif 
      !! Search for P,V:
      do k=CEILING(bp(n2)%zbottom/dz),FLOOR(bp(n2)%height/dz) !do k=0,k1
	   do tel=1,bp(n2)%tmax 
		 i=bp(n2)%i(tel) 
		 j=bp(n2)%j(tel) 	  
!!       do i=0,i1  
!!         do j=j1,0,-1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
!!	  xx=Rp(i)*cos_u(j)-schuif_x
!!	  yy=Rp(i)*sin_u(j)
!!!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
!!		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
!!		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
!!!	  ELSE 
!!!	 	inout=0
!!!	  ENDIF
!!	  if (inout.eq.1) then
		rnew(i,j,k)=rho(i,j,k)
		rold(i,j,k)=rho(i,j,k)
		! prevent large source in pres-corr by sudden increase in density with bp%c 
!!	   endif
!!	  enddo
	 enddo
	enddo
	ENDIF
	ENDDO ! bedplume loop

      end
	  

c  All the MPI - shit
c

      subroutine pshiftb(UT,UP) !(UT,UP,rank,imax,jmax,kmax,px)

       USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf,ierr
!       parameter (ileng= (kmax)*(imax) )
      include 'mpif.h'
      integer itag,status(MPI_STATUS_SIZE),l  !,rank,imax,jmax,kmax,px
      real ut(imax,jmax,kmax)
      real*8 up(imax,kmax),UTMP(imax,kmax)
!      integer i,j,k
!	real up(imax,kmax),UTMP(imax,kmax)
      do i=1,imax
	 do k=1,kmax
	  utmp(i,k) =UT(i,1,k)
          enddo
      enddo
      itag = 10
      ileng= (kmax)*(imax)
      rankf=rank+1
      rankb=rank-1
	if (periodicy.eq.0.or.periodicy.eq.2) then
           if(rank.eq.px-1)rankf=MPI_PROC_NULL
           if(rank.eq.   0)rankb=MPI_PROC_NULL
	else 
      	   if(rank.eq.px-1)rankf=0 ! MPI_PROC_NULL
	   if(rank.eq.   0)rankb=px-1 !MPI_PROC_NULL
	endif
      call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankb,itag,
     $                  up ,ileng,MPI_REAL8,rankf,itag, MPI_COMM_WORLD,status,ierr)
!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankb,itag,
!     $                  up ,ileng,MPI_REAL,rankf,itag, MPI_COMM_WORLD,status,ierr)
c      if (rank.eq.   0) then
c         call MPI_SEND(utmp  ,ileng,MPI_REAL8,px-1,itag,MPI_COMM_WORLD,ierr)
c        endif
c       if (rank.eq.px-1) then
c           call MPI_RECV(up,ileng,MPI_REAL8,0 ,itag,MPI_COMM_WORLD,status,ierr)
c          endif

      end

   


C ... Perfroms data transform across the nodes.  The data
C     starts out Q(j,k,i) node distributed in x and is transformed
C     to Qt(i,j,k) node distributed in theta.

      SUBROUTINE t2fp(Q,Qt)

      USE nlist
!         INCLUDE 'param.txt'
        INCLUDE 'mpif.h'
	integer N
!         parameter (nr=imax,mt=jmax,nx=kmax,mx=kmax/px,NT=jmax*px)
	integer nr,mt,nx,mx,nt,ierr
! 	real*8  Q(Nr,Nt,Mx),Qt(Nr,Mt,Nx)
	real*8  Q(imax,jmax*px,kmax/px),Qt(imax,jmax,kmax)

	!real  Q(Nr,Nt,Mx),Qt(Nr,Mt,Nx)
		
C ...  Locals
	integer ii,kk,l
	
! 	real*8  W1(Mx,Nr,Nt),W2(Mx,Nr,Nt)
	real*8  W1(kmax/px,imax,jmax*px),W2(kmax/px,imax,jmax*px)
	!real  W1(Mx,Nr,Nt),W2(Mx,Nr,Nt)
! 	common  W1,W2

	nr=imax
	mt=jmax
	nx=kmax	
	mx=kmax/px
	NT=jmax*px
	
	  
	  do i = 1,Mx
	    do k = 1,Nt
	      do j = 1,Nr
		W1(i,j,k) = Q(j,k,i)
	      end do
	    end do
	  end do
          
	  call MPI_ALLTOALL(W1,Nr*Mt*Mx,MPI_REAL8,W2,Nr*Mt*Mx,
     &         MPI_REAL8,MPI_COMM_WORLD,ierr)
!	  call MPI_ALLTOALL(W1,Nr*Mt*Mx,MPI_REAL,W2,Nr*Mt*Mx,
!     &         MPI_REAL,MPI_COMM_WORLD,ierr)	  



	  do k = 1,Mt
	    do j = 1,Nr
	      do i = 1,Nx
		Qt(j,k,i) = W2(Xii(i),j,Xkk(i,k))
	      end do
	    end do
	  end do
	  
		
      RETURN
      END
      
C ... Perfroms data transform across the nodes.  The data
C     starts out Q(j,k,i) node distributed in x and is transformed
C     to Qt(i,j,k) node distributed in theta.

      SUBROUTINE t2fp_0ie(Q,Qt)

      USE nlist
!         INCLUDE 'param.txt'
        INCLUDE 'mpif.h'
	integer N
!         parameter (nr=imax,mt=jmax,nx=kmax,mx=kmax/px,NT=jmax*px)
	integer nr,mt,nx,mx,nt,ierr
! 	real*8  Q(Nr,Nt,Mx),Qt(Nr,Mt,Nx)
	real*8  Q(imax+1,jmax*px,kmax/px),Qt(imax+1,jmax,kmax)

	!real  Q(Nr,Nt,Mx),Qt(Nr,Mt,Nx)
		
C ...  Locals
	integer ii,kk,l
	
! 	real*8  W1(Mx,Nr,Nt),W2(Mx,Nr,Nt)
	real*8  W1(kmax/px,imax+1,jmax*px),W2(kmax/px,imax+1,jmax*px)
	!real  W1(Mx,Nr,Nt),W2(Mx,Nr,Nt)
! 	common  W1,W2

	nr=imax+1
	mt=jmax
	nx=kmax	
	mx=kmax/px
	NT=jmax*px
	
	  
	  do i = 1,Mx
	    do k = 1,Nt
	      do j = 1,Nr
		W1(i,j,k) = Q(j,k,i)
	      end do
	    end do
	  end do
          
	  call MPI_ALLTOALL(W1,Nr*Mt*Mx,MPI_REAL8,W2,Nr*Mt*Mx,
     &         MPI_REAL8,MPI_COMM_WORLD,ierr)
!	  call MPI_ALLTOALL(W1,Nr*Mt*Mx,MPI_REAL,W2,Nr*Mt*Mx,
!     &         MPI_REAL,MPI_COMM_WORLD,ierr)	  



	  do k = 1,Mt
	    do j = 1,Nr
	      do i = 1,Nx
		Qt(j,k,i) = W2(Xii(i),j,Xkk(i,k))
	      end do
	    end do
	  end do
	  
		
      RETURN
      END
      
C ... Performs data transform across the nodes.  The data
C     starts out (i,j,k) node distributed in theta and is transformed
C     to (j,k,i) node distributed in x.

      SUBROUTINE t2np(Q,Qt)
      USE nlist
!         INCLUDE 'param.txt'
        INCLUDE 'mpif.h'
!         parameter (nr=imax,mt=jmax,nx=kmax,mx=kmax/px,NT=jmax*px)
	integer nr,mt,nx,mx,nt,ierr
	integer N
!         real*8 Q(Nr,Mt,Nx),Qt(Nr,Nt,Mx)
	real*8 Q(imax,jmax,kmax),Qt(imax,jmax*px,kmax/px)
	!real Q(Nr,Mt,Nx),Qt(Nr,Nt,Mx)

C ...  Locals
	integer ii,kk,l

C ...  Globals
! 	integer Tkk(Nt),Tii(Nt,Mx)
! 	common /TPOSE/ Tkk,Tii
	  	  
! 	real*8  W1(Nr,Mt,Nx),W2(Nr,Mt,Nx)
	real*8  W1(imax,jmax,kmax),W2(imax,jmax,kmax)
	!real  W1(Nr,Mt,Nx),W2(Nr,Mt,Nx)
! 	common  W1,W2

	nr=imax
	mt=jmax
	nx=kmax	
	mx=kmax/px
	NT=jmax*px		  
	  do k = 1,Mt
	    do j = 1,Nr
	      do i = 1,Nx
		W1(j,k,i) = Q(j,k,i)
	      end do
	    end do
	  end do
	  
	  call MPI_ALLTOALL(W1,Nr*Mt*Mx,MPI_REAL8,W2,Nr*Mt*Mx,
     &         MPI_REAL8,MPI_COMM_WORLD,ierr)
!	  call MPI_ALLTOALL(W1,Nr*Mt*Mx,MPI_REAL,W2,Nr*Mt*Mx,
!     &         MPI_REAL,MPI_COMM_WORLD,ierr)    

 
	  do i = 1,Mx
	    do k = 1,Nt
	      do j = 1,Nr
		Qt(j,k,i) = W2(j,Tkk(k),Tii(k,i))
	      end do
	    end do
	  end do
	  

      RETURN
      END

C ... Performs data transform across the nodes.  The data
C     starts out (i,j,k) node distributed in theta and is transformed
C     to (j,k,i) node distributed in x.

      SUBROUTINE t2np_i1(Q,Qt)
      USE nlist
!         INCLUDE 'param.txt'
        INCLUDE 'mpif.h'
!         parameter (nr=imax,mt=jmax,nx=kmax,mx=kmax/px,NT=jmax*px)
	integer nr,mt,nx,mx,nt,ierr
	integer N
!         real*8 Q(Nr,Mt,Nx),Qt(Nr,Nt,Mx)
	real*8 Q(imax+1,jmax,kmax),Qt(imax+1,jmax*px,kmax/px)
	!real Q(Nr,Mt,Nx),Qt(Nr,Nt,Mx)

C ...  Locals
	integer ii,kk,l

C ...  Globals
! 	integer Tkk(Nt),Tii(Nt,Mx)
! 	common /TPOSE/ Tkk,Tii
	  	  
! 	real*8  W1(Nr,Mt,Nx),W2(Nr,Mt,Nx)
	real*8  W1(imax+1,jmax,kmax),W2(imax+1,jmax,kmax)
	!real  W1(Nr,Mt,Nx),W2(Nr,Mt,Nx)
! 	common  W1,W2

	nr=imax+1
	mt=jmax
	nx=kmax	
	mx=kmax/px
	NT=jmax*px		  
	  do k = 1,Mt
	    do j = 1,Nr
	      do i = 1,Nx
		W1(j,k,i) = Q(j,k,i)
	      end do
	    end do
	  end do
	  
	  call MPI_ALLTOALL(W1,Nr*Mt*Mx,MPI_REAL8,W2,Nr*Mt*Mx,
     &         MPI_REAL8,MPI_COMM_WORLD,ierr)
!	  call MPI_ALLTOALL(W1,Nr*Mt*Mx,MPI_REAL,W2,Nr*Mt*Mx,
!     &         MPI_REAL,MPI_COMM_WORLD,ierr)    

 
	  do i = 1,Mx
	    do k = 1,Nt
	      do j = 1,Nr
		Qt(j,k,i) = W2(j,Tkk(k),Tii(k,i))
	      end do
	    end do
	  end do
	  

      RETURN
      END
  
