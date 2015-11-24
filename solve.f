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
      integer  ib,ie,jb,je,kb,ke,n,t
      real     doldc(nfrac,0:i1,0:j1,0:k1),dnewc(nfrac,0:i1,0:j1,0:k1)
      real     dold(0:i1,0:j1,0:k1),dnew(0:i1,0:j1,0:k1)
      real     wsed(nfrac,0:i1,0:j1,0:k1),wfluid(0:i1,0:j1,0:k1)
	real*8 pplus(imax,kmax)
	

!	real     Diffcof(0:i1,0:j1,0:k1)

!	Diffcof=ekm/Sc/Rnew

c********************************************************************
c     CALCULATE slipvelocity
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(cnew,Wnew,wsed,rnew,1,kmax,wfluid,dt,dz)
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
	      call advecc_TVD(dnewc(n,:,:,:),cnew(n,:,:,:),Unew,Vnew,Wsed(n,:,:,:),rnew,Ru,Rp,dr,dphi,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
	      call diffc_CDS2 (dnewc(n,:,:,:),Cnew(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)

		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(1.5*dnewc(n,:,:,:)-0.5*cc(n,:,:,:)) !Adams-Bashford
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc for intermediate dcdt
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
      call advecu_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(dnew,Srr,Spr,Szr,Spp,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
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
      call advecv_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(dnew,Spr,Spp,Spz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
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
      call advecw_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(dnew,Szr,Spz,Szz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	if (slipvel.eq.1.or.slipvel.eq.2) then
	do n=1,nfrac
	  call advecw_driftfluxCDS2(dnew,0.,0.,wsed(n,:,:,:)-Wnew,
     &  cnew(n,:,:,:)*frac(n)%rho,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
	enddo
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

!      i = imax  !! at imax do nothing -> dpdn=0 at outflow
!      do k=1,kmax
!	do j=1,jmax
!		dUdt(i,j,k) = dUdt(i,j,k) -dt * ( -2* pold(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
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
          dVdt(i,j,k) = dVdt(i,j,k) -dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * dphi )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
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
      integer  ib,ie,jb,je,kb,ke,n,t
      real     doldc(nfrac,0:i1,0:j1,0:k1),dnewc(nfrac,0:i1,0:j1,0:k1)
      real     dold(0:i1,0:j1,0:k1),dnew(0:i1,0:j1,0:k1)
      real     wsed(nfrac,0:i1,0:j1,0:k1),wfluid(0:i1,0:j1,0:k1),W_km_sum(0:i1,0:j1,0:k1)
	real*8 pplus(imax,kmax)
	real dnewcbot(nfrac,0:i1,0:j1)
	

!	real     Diffcof(0:i1,0:j1,0:k1)

!	Diffcof=ekm/Sc/Rnew 

c********************************************************************
c     CALCULATE slipvelocity
c********************************************************************
c********************************************************************
c     CALCULATE advection, diffusion Concentration
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(cnew,Wnew,wsed,rnew,1,kmax,wfluid,dt,dz)
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
	      call advecc_TVD(dnewc(n,:,:,:),cnew(n,:,:,:),Unew,Vnew,Wsed(n,:,:,:),rnew,Ru,Rp,dr,dphi,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
		
	      call diffc_CDS2 (dnewc(n,:,:,:),Cnew(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)
		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(1.5*dnewc(n,:,:,:)-0.5*cc(n,:,:,:)) !AB2
		cc(n,:,:,:)=dnewc(n,:,:,:)
	      	  if (interaction_bed>0) then
		    !! advec concentration in bed with velocity TSHD:
	      	    call adveccbot_TVD(dnewcbot(n,:,:),cnewbot(n,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,dphi,dz,
     +            	i1,j1,ib,ie,jb,je,dt,rank,px,periodicx,periodicy)
		    dcdtbot(n,:,:) =cnewbot(n,:,:) + dt*(1.5*dnewcbot(n,:,:)-0.5*ccbot(n,:,:)) !AB2
		    ccbot(n,:,:)=dnewcbot(n,:,:)
	      	  endif
	  enddo
      	  if (interaction_bed>0) then
	    call slipvelocity(cnew,wnew,wsed,rnew,0,0,wfluid,dt,dz) !driftflux_force must be calculated with settling velocity at k=0
	    CALL erosion_deposition(dcdt,dcdtbot,unew,vnew,wnew,rnew,cnew,cnewbot,dt,dz) !first two vars are adjusted
	    !! for kn1: erosion_deposition must be after advecc_TVD and after dcdt update, because cnew and dcdt are two separate vars
	    do n=1,nfrac 
	        call bound_cbot(dcdtbot(n,:,:))
	    enddo
      	  endif
	    do n=1,nfrac 
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc after erosion_deposition 
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
      call advecu_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(dnew,Srr,Spr,Szr,Spp,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
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
      call advecv_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(dnew,Spr,Spp,Spz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
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
      call advecw_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(dnew,Szr,Spz,Szz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 W_km_sum=0.
	 do n=1,nfrac
	  W_km_sum=W_km_sum+wsed(n,:,:,:)-Wnew
	 enddo 
	 W_km_sum=W_km_sum+Wfluid !Wfluid is difference with Ucfd 
	 call advecw_driftfluxCDS2(dnew,0.,0.,W_km_sum,
     &   rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
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
!      i = imax  !! at imax do nothing -> dpdn=0 at outflow
!      do k=1,kmax
!	do j=1,jmax
!		dUdt(i,j,k) = dUdt(i,j,k) -dt * ( -2* pold(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
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
          dVdt(i,j,k) = dVdt(i,j,k) -dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * dphi )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
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
!       include 'param.txt'
!       include 'common.txt'
      integer  ib,ie,jb,je,kb,ke,n,t
      real     doldc(nfrac,0:i1,0:j1,0:k1),dnewc(nfrac,0:i1,0:j1,0:k1)
      real     dold(0:i1,0:j1,0:k1),dnew(0:i1,0:j1,0:k1)
      real     wsed(nfrac,0:i1,0:j1,0:k1),wfluid(0:i1,0:j1,0:k1),W_km_sum(0:i1,0:j1,0:k1)
	real*8 pplus(imax,kmax)
	real fAB3_1,fAB3_2,fAB3_3,fAB3_4,facAB3_1,facAB3_2,facAB3_3,facAB3_4
	real dist,timeAB_desired(1:4)
	

!	real     Diffcof(0:i1,0:j1,0:k1)

!	Diffcof=ekm/Sc/Rnew 

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

	if (istep.lt.5) then
	   facAB3_1=1.
	   facAB3_2=0.
	   facAB3_3=0.
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

		facAB3_1=23./12.-fAB3_1*4./3.
		facAB3_2=-fAB3_2*4./3.
		facAB3_3=-fAB3_3*4./3.

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

		facAB3_2=facAB3_2+fAB3_2*5./12.
		facAB3_3=facAB3_3+fAB3_3*5./12.
		facAB3_4=fAB3_4*5./12.


!		write(*,*) 'rank,istep,f1,f2,f3,f4',rank,istep,facAB3_1,facAB3_2,facAB3_3,facAB3_4
	endif


c********************************************************************
c     CALCULATE slipvelocity
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(cnew,Wnew,wsed,rnew,1,kmax,wfluid,dt,dz)
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
	      call advecc_TVD(dnewc(n,:,:,:),cnew(n,:,:,:),Unew,Vnew,Wsed(n,:,:,:),rnew,Ru,Rp,dr,dphi,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
	      call diffc_CDS2 (dnewc(n,:,:,:),Cnew(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)

!		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(1.5*dnewc(n,:,:,:)-0.5*cc(n,:,:,:)) !AB2
		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*dnewc(n,:,:,:)                       !EE1
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc for intermediate dcdt
		cc(n,:,:,:)=dnewc(n,:,:,:)
	  enddo

!	  do n=1,nfrac
!	      call advecc_TVD(dnewc(n,:,:,:),cnew(n,:,:,:),Unew,Vnew,Wsed(n,:,:,:),rnew,Ru,Rp,dr,dphi,dz,
!     +            i1,j1,k1,ib,ie,jb,je,kb,ke)
!	      call diffc (dnewc(n,:,:,:),Cnew(n,:,:,:),Diffcof,
!     +            ib,ie,jb,je,kb,ke)

!		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(23./12.*dnewc(n,:,:,:)-4./3.*cc(n,:,:,:)+5./12.*ccold(n,:,:,:)) !AB3
!		call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc for intermediate dcdt
!		ccold(n,:,:,:)=cc(n,:,:,:)
!		cc(n,:,:,:)=dnewc(n,:,:,:)
!	  enddo

	  call state(dcdt,drdt) ! determine drdt with intermediate dcdt
	endif

	if (diffusion.eq.'COM4') then
		call stress_terms(Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	endif
c********************************************************************
c     CALCULATE advection, diffusion and Force U-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecu_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(dnew,Srr,Spr,Szr,Spp,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif


      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-(gx*cos_u(j)+gy*sin_u(j))*(0.5*(rnew(i,j,k)+rnew(i+1,j,k))-rho_b)
     1     +Ppropx(i,j,k) 

            dUdt(i,j,k)=Unew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i+1,j,k))+
     1     dt*(     facAB3_1*dnew(i,j,k)+
     1              facAB3_2*wx(i,j,k)+
     1              facAB3_3*wxold(i,j,k)+
     1              facAB3_4*wxolder(i,j,k))
!     1     -(gx*cos_u(j)+gy*sin_u(j))*dt*(0.75*(rnew(i,j,k)+rnew(i+1,j,k))-0.25*(rold(i,j,k)+rold(i+1,j,k))-rho_b)  
              wxolder(i,j,k) = wxold(i,j,k)
              wxold(i,j,k)   = wx(i,j,k)
              wx(i,j,k)      = dnew(i,j,k)
            enddo
         enddo
      enddo

c********************************************************************
c     CALCULATE advection, diffusion and Force V-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecv_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(dnew,Spr,Spp,Spz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-(-gx*sin_v(j)+gy*cos_v(j))*(0.5*(rnew(i,j,k)+rnew(i,j+1,k))-rho_b)
     1     +Ppropy(i,j,k) 
            dVdt(i,j,k)=Vnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j+1,k))+
     1     dt*(      facAB3_1*dnew(i,j,k)+
     1               facAB3_2*wy(i,j,k)+
     1               facAB3_3*wyold(i,j,k)+
     1               facAB3_4*wyolder(i,j,k))
!     1     -(-gx*sin_v(j)+gy*cos_v(j))*dt*(0.75*(rnew(i,j,k)+rnew(i,j+1,k))-0.25*(rold(i,j,k)+rold(i,j+1,k))-rho_b)
             wyolder(i,j,k) = wyold(i,j,k)
             wyold(i,j,k)   = wy(i,j,k)
             wy(i,j,k)      = dnew(i,j,k)
            enddo
         enddo
      enddo
c********************************************************************
c     CALCULATE advection, diffusion and Force W-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecw_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(dnew,Szr,Spz,Szz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

	if (slipvel.eq.1.or.slipvel.eq.2) then
	 if (interaction_bed.gt.0) then
	  call slipvelocity(cnew,wnew,wsed,rnew,0,0,wfluid,dt,dz) !driftflux_force must be calculated with settling velocity at k=0
	! should correct Wsed(n,:,:,0) with (1.-tau/frac(n)%tau_d), but this correction almost always is zero, for now leave it like this
	 endif
	 W_km_sum=0.
	 do n=1,nfrac
	  W_km_sum=W_km_sum+wsed(n,:,:,:)-Wnew
	 enddo 
	 W_km_sum=W_km_sum+Wfluid !Wfluid is difference with Ucfd 
	 call advecw_driftfluxCDS2(dnew,0.,0.,W_km_sum,
     &   rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
	endif

      do k=1,kmax
         do j=1,jmax
            do i=1,imax
	    dnew(i,j,k)=dnew(i,j,k)-gz*(0.5*(rnew(i,j,k)+rnew(i,j,k+1))-rho_b)
     1     +Ppropz(i,j,k)
            dWdt(i,j,k)=Wnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j,k+1))+
     1     dt*(      facAB3_1*dnew(i,j,k)+
     1               facAB3_2*wz(i,j,k)+
     1               facAB3_3*wzold(i,j,k)+
     1               facAB3_4*wzolder(i,j,k))
!     1          dt*(1.5*dnew(i,j,k)-0.5*wz(i,j,k)) -gz*dt*(0.75*(rnew(i,j,k)+rnew(i,j,k+1))-0.25*(rold(i,j,k)+rold(i,j,k+1))-rho_b)
		    wzolder(i,j,k) = wzold(i,j,k)
		    wzold(i,j,k)   = wz(i,j,k)
		    wz(i,j,k)      = dnew(i,j,k)
            enddo
         enddo
      enddo


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
!      i = imax  !! at imax do nothing -> dpdn=0 at outflow
!      do k=1,kmax
!	do j=1,jmax
!		dUdt(i,j,k) = dUdt(i,j,k) -dt * ( -2* pold(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
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
          dVdt(i,j,k) = dVdt(i,j,k) -dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * dphi )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
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
      integer  ib,ie,jb,je,kb,ke,n,t,kpp,km,kp
      real     wsed(nfrac,0:i1,0:j1,0:k1),wfluid(0:i1,0:j1,0:k1)
      real     c2,c3,a21,a31,a32,b1,b2,b3,cn1,cn2,cn3
	real   b1tvd,b2tvd,b3tvd
!	real     Diffcof(0:i1,0:j1,0:k1)
	real	 invrho_b
	real*8 pplus(imax,kmax),pplus2(imax,kmax)
	real   W_km_sum(0:i1,0:j1,0:k1),t_output,dWdt_old(0:i1,0:j1,0:k1),c_sum(0:i1,0:j1,0:k1)
!	real aaa(1:kmax),bbb(1:kmax),ccc(1:kmax),ekm_min,ekm_plus,rhss(1:kmax)
	real aaa(0:k1),bbb(0:k1),ccc(0:k1),ekm_min,ekm_plus,rhss(0:k1)
	real utr(0:i1,0:j1,0:k1),vtr(0:i1,0:j1,0:k1),wtr(0:i1,0:j1,0:k1)

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


!	Diffcof=ekm/Sc/Rnew 

	dcdt=0.
	Wfluid=0.

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

	pold=p+pold    !what is called p here was dp in reality, now p is 
        call pshiftb(pold,pplus) !,rank,imax,jmax,kmax,px)
c********************************************************************
c     CALCULATE k1 and predictor 1
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(cnew,Wnew,wsed,rnew,1,kmax,wfluid,cn1*dt,dz)
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
	      call make_UtransC_zeroIBM(utr,vtr,wtr)
	      call advecc_TVD(k1c(n,:,:,:),cnew(n,:,:,:),utr,vtr,wtr,rnew,Ru,Rp,dr,dphi,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn1*dt,rank,px,periodicx,periodicy)
!	      call advecc_TVD(k1c(n,:,:,:),cnew(n,:,:,:),Unew,Vnew,Wsed(n,:,:,:),rnew,Ru,Rp,dr,dphi,dz,
!     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn1*dt,rank,px,periodicx,periodicy)
		
	      call diffc_CDS2 (k1c(n,:,:,:),Cnew(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)
!		dcdt(n,:,:,:) =cnew(n,:,:,:) + a21*dt*(k1c(n,:,:,:)) !pred1
		dcdt(n,:,:,:) =cnew(n,:,:,:) + cn1*dt*(k1c(n,:,:,:)) !pred1
	      	  if (interaction_bed>0) then
		    !! advec concentration in bed with velocity TSHD:
	      	    call adveccbot_TVD(k1cbot(n,:,:),cnewbot(n,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,dphi,dz,
     +            	i1,j1,ib,ie,jb,je,cn1*dt,rank,px,periodicx,periodicy)
		    dcdt1bot(n,:,:)= cnewbot(n,:,:) + cn1*dt*k1cbot(n,:,:) ! pred1
	      	  endif
           IF (CNdiffz.eq.1) THEN !CN semi implicit treatment diff-z:
	     call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc after erosion_deposition 
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
	           ekm_min=0.5*(Diffcof(i,j,km)+Diffcof(i,j,k))
	           ekm_plus=0.5*(Diffcof(i,j,kp)+Diffcof(i,j,k))
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
	   ENDIF
          enddo
      	  if (interaction_bed>0) then
	    call slipvelocity(cnew,wnew,wsed,rnew,0,0,wfluid,cn1*dt,dz) !driftflux_force must be calculated with settling velocity at k=0
	    CALL erosion_deposition(dcdt,dcdt1bot,unew,vnew,wnew,rnew,cnew,cnewbot,cn1*dt,dz) !first two vars are adjusted
	    !! for kn1: erosion_deposition must be after advecc_TVD and after dcdt update, because cnew and dcdt are two separate vars
	    do n=1,nfrac 
	        call bound_cbot(dcdt1bot(n,:,:))
	    enddo
      	  endif
	    do n=1,nfrac 
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc after erosion_deposition 
	    enddo
	endif

	if (diffusion.eq.'COM4') then
		call stress_terms(Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	endif

      	if (convection.eq.'CDS2') then
      call advecu_CDS2(k1u,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(k1u,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(k1u,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(k1u,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(k1u,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(k1u,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif
	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (k1u,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(k1u,Srr,Spr,Szr,Spp,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

	if(convection.eq.'CDS2') then
      call advecv_CDS2(k1v,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(k1v,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(k1v,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(k1v,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(k1v,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(k1v,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif
	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (k1v,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(k1v,Spr,Spp,Spz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      	if (convection.eq.'CDS2') then
      call advecw_CDS2(k1w,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(k1w,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(k1w,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(k1w,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(k1w,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(k1w,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif
	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (k1w,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(k1w,Szr,Spz,Szz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 W_km_sum=0.
	 c_sum=0.
	 do n=1,nfrac
	  W_km_sum=W_km_sum+cnew(n,:,:,:)*frac(n)%rho*(wsed(n,:,:,:)-Wnew)
	  c_sum=c_sum+cnew(n,:,:,:)
	 enddo 
	 W_km_sum=W_km_sum+(1.-c_sum)*rho_b*Wfluid !Wfluid is difference with Ucfd 
	 W_km_sum=W_km_sum/rnew
	 call advecw_driftfluxCDS2(k1w,0.,0.,W_km_sum,
     &   rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
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
      call bound_rhoU(dUdt,dVdt,dWdt,drdt,0,time_n+c2*dt,(1.-c2)*Ub1old+c2*Ub1new,(1.-c2)*Vb1old+c2*Vb1new,
     & (1.-c2)*Wb1old+c2*Wb1new,(1.-c2)*Ub2old+c2*Ub2new,(1.-c2)*Vb2old+c2*Vb2new,(1.-c2)*Wb2old+c2*Wb2new)
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))
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
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))
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
	    ekm_min=ekm(i,j,k)
	    ekm_plus=ekm(i,j,kp)
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
!      i = imax
!      do k=1,kmax
!	do j=1,jmax
!		dUdt(i,j,k) = dUdt(i,j,k) -c2*dt * ( -2* p(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
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
          dVdt(i,j,k) = dVdt(i,j,k) -c2*dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * dphi )   
          enddo
        enddo
      enddo

	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c2*dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c2*dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
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
      call bound_rhoU(dUdt,dVdt,dWdt,drdt,0,time_n+c2*dt,(1.-c2)*Ub1old+c2*Ub1new,(1.-c2)*Vb1old+c2*Vb1new,
     & (1.-c2)*Wb1old+c2*Wb1new,(1.-c2)*Ub2old+c2*Ub2new,(1.-c2)*Vb2old+c2*Vb2new,(1.-c2)*Wb2old+c2*Wb2new)

      call fillps2(dUdt,dVdt,dWdt,dRdt,time_n+c2*dt,c2*dt) 
      CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
      call pshiftb(p,pplus2) !,rank,imax,jmax,kmax,px)
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -c2*dt * ( p(i+1,j,k) - p(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
!      i = imax
!      do k=1,kmax
!	do j=1,jmax
!		dUdt(i,j,k) = dUdt(i,j,k) -c2*dt * ( -2* p(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
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
          dVdt(i,j,k) = dVdt(i,j,k) -c2*dt * ( p(i,j+1,k) - p(i,j,k) ) /( Rp(i) * dphi )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c2*dt * ( pplus2(i,k) - p(i,jmax,k) ) /( Rp(i) * dphi )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c2*dt * ( pplus2(i,k) - p(i,jmax,k) ) /( Rp(i) * dphi )
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
!      call bound_rhoU(dUdt1,dVdt1,dWdt1,drdt1,slip_bot,time_n+c2*dt,(1.-c2)*Ub1old+c2*Ub1new,(1.-c2)*Vb1old+c2*Vb1new,
!     & (1.-c2)*Wb1old+c2*Wb1new,(1.-c2)*Ub2old+c2*Ub2new,(1.-c2)*Vb2old+c2*Vb2new,(1.-c2)*Wb2old+c2*Wb2new)

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
           dUdt(i,j,k)=dUdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i+1,j,k)))
           dVdt(i,j,k)=dVdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i,j+1,k)))
           dWdt(i,j,k)=dWdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i,j,k+1)))
          enddo
        enddo
      enddo   

      call bound_incljet(dUdt,dVdt,dWdt,drdt,0,time_n+c2*dt,(1.-c2)*Ub1old+c2*Ub1new,(1.-c2)*Vb1old+c2*Vb1new,
     & (1.-c2)*Wb1old+c2*Wb1new,(1.-c2)*Ub2old+c2*Ub2new,(1.-c2)*Vb2old+c2*Vb2new,(1.-c2)*Wb2old+c2*Wb2new)
!      call bound_incljet(dUdt1,dVdt1,dWdt1,drdt1,slip_bot,time_n+c2*dt,(1.-c2)*Ub1old+c2*Ub1new,(1.-c2)*Vb1old+c2*Vb1new,
!     & (1.-c2)*Wb1old+c2*Wb1new,(1.-c2)*Ub2old+c2*Ub2new,(1.-c2)*Vb2old+c2*Vb2new,(1.-c2)*Wb2old+c2*Wb2new)

c********************************************************************
c     CALCULATE k2 and predictor 2
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(dcdt,dWdt,wsed,drdt,1,kmax,wfluid,cn2*dt,dz)
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
	      call make_UtransC_zeroIBM(utr,vtr,wtr)
	      call advecc_TVD(k2c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,dphi,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn2*dt,rank,px,periodicx,periodicy)
!	      call advecc_TVD(k2c(n,:,:,:),dcdt(n,:,:,:),dUdt,dVdt,Wsed(n,:,:,:),drdt,Ru,Rp,dr,dphi,dz,
!     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn2*dt,rank,px,periodicx,periodicy)
	      call diffc_CDS2 (k2c(n,:,:,:),dCdt(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)
	          if (interaction_bed>0) then
		    !! advec concentration in bed with velocity TSHD:
	            call adveccbot_TVD(k2cbot(n,:,:),dcdt1bot(n,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,dphi,dz,
     +              i1,j1,ib,ie,jb,je,cn2*dt,rank,px,periodicx,periodicy)
		    dcdt2bot(n,:,:)= dcdt1bot(n,:,:) + cn2*dt*k2cbot(n,:,:) ! pred2
	          endif
	  enddo
	  dcdt2=dcdt
	  do n=1,nfrac
		dcdt(n,:,:,:) =dcdt(n,:,:,:) + cn2*dt*k2c(n,:,:,:) !pred2
!		dcdt(n,:,:,:) =cnew(n,:,:,:) + a31*dt*k1c(n,:,:,:) + a32*dt*k2c(n,:,:,:) !pred2
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc for intermediate pred1
           IF (CNdiffz.eq.1) THEN !CN semi implicit treatment diff-z:
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
	           ekm_min=0.5*(Diffcof(i,j,km)+Diffcof(i,j,k))
	           ekm_plus=0.5*(Diffcof(i,j,kp)+Diffcof(i,j,k))
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
	     call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc after CN-diffz 
	   ENDIF
	  enddo
	  if (interaction_bed>0) then
	    call slipvelocity(dcdt,dwdt,wsed,drdt,0,0,wfluid,cn2*dt,dz) !get wfluid and wsed at k=0 for driftflux_force 
	    CALL erosion_deposition(dcdt,dcdt2bot,dudt,dvdt,dwdt,drdt,dcdt2,dcdt1bot,cn2*dt,dz) !first two vars are adjusted 
 	    !! for kn2: variable dcdt2 made
	    do n=1,nfrac
		call bound_cbot(dcdt2bot(n,:,:))
	    enddo
	  endif
	endif

	if (diffusion.eq.'COM4') then
		call stress_terms(dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	endif

      	if (convection.eq.'CDS2') then
      call advecu_CDS2(k2u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(k2u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(k2u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(k2u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(k2u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(k2u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif
	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (k2u,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(k2u,Srr,Spr,Szr,Spp,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      	if (convection.eq.'CDS2') then
      call advecv_CDS2(k2v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(k2v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(k2v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(k2v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(k2v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(k2v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif
	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (k2v,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(k2v,Spr,Spp,Spz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      	if (convection.eq.'CDS2') then
      call advecw_CDS2(k2w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(k2w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(k2w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(k2w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(k2w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(k2w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif
	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (k2w,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(k2w,Szr,Spz,Szz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 W_km_sum=0.
	 c_sum=0.
	 do n=1,nfrac
	  W_km_sum=W_km_sum+dcdt(n,:,:,:)*frac(n)%rho*(wsed(n,:,:,:)-dWdt)
	  c_sum=c_sum+dcdt(n,:,:,:)
	 enddo
	 W_km_sum=W_km_sum+(1.-c_sum)*rho_b*Wfluid  !Wfluid is difference with Ucfd
	 W_km_sum=W_km_sum/drdt
	 call advecw_driftfluxCDS2(k2w,0.,0.,W_km_sum,
     &   drdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
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
      call bound_rhoU(dUdt,dVdt,dWdt,dRdt,0,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new)
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))
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
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))
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
	    ekm_min=ekm(i,j,k)
	    ekm_plus=ekm(i,j,kp)
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
!      i = imax
!      do k=1,kmax
!	do j=1,jmax
!		dUdt(i,j,k) = dUdt(i,j,k) -c3*dt * ( -2* p(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
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
          dVdt(i,j,k) = dVdt(i,j,k) -c3*dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * dphi )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c3*dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c3*dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
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
      call bound_rhoU(dUdt,dVdt,dWdt,dRdt,0,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new)
!      call bound_rhoU(dUdt2,dVdt2,dWdt2,dRdt2,slip_bot,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
!     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new)
      call fillps2(dUdt,dVdt,dWdt,dRdt,time_n+c3*dt,c3*dt) 
      CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
      call pshiftb(p,pplus2) !,rank,imax,jmax,kmax,px)
      do k=1,kmax
        do j=1,jmax
          do i=1,imax-1
          dUdt(i,j,k) = dUdt(i,j,k) -c3*dt * ( p(i+1,j,k) - p(i,j,k) ) /( Rp(i+1) - Rp(i) )
          enddo
        enddo	
      enddo
!      i = imax
!      do k=1,kmax
!	do j=1,jmax
!		dUdt(i,j,k) = dUdt(i,j,k) -c3*dt * ( -2* p(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
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
          dVdt(i,j,k) = dVdt(i,j,k) -c3*dt * ( p(i,j+1,k) - p(i,j,k) ) /( Rp(i) * dphi )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c3*dt * ( pplus2(i,k) - p(i,jmax,k) ) /( Rp(i) * dphi )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -c3*dt * ( pplus2(i,k) - p(i,jmax,k) ) /( Rp(i) * dphi )
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
!      call bound_rhoU(dUdt2,dVdt2,dWdt2,drdt2,slip_bot,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
!     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new)	

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
           dUdt(i,j,k)=dUdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i+1,j,k)))
           dVdt(i,j,k)=dVdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i,j+1,k)))
           dWdt(i,j,k)=dWdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i,j,k+1)))
          enddo
        enddo
      enddo       

      call bound_incljet(dUdt,dVdt,dWdt,dRdt,0,time_n+c3*dt,(1.-c3)*Ub1old+c3*Ub1new,(1.-c3)*Vb1old+c3*Vb1new,
     & (1.-c3)*Wb1old+c3*Wb1new,(1.-c3)*Ub2old+c3*Ub2new,(1.-c3)*Vb2old+c3*Vb2new,(1.-c3)*Wb2old+c3*Wb2new)
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
	      call slipvelocity(dcdt,dWdt,wsed,drdt,1,kmax,wfluid,cn3*dt,dz)
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
	  	
!		dwdt_old=dwdt
!		dWdt(:,:,1:kmax)=0.5*(dcdt(1,:,:,1:kmax)*frac(1)%rho/drdt(:,:,1:kmax)
!     & +dcdt(1,:,:,2:kmax+1)*frac(1)%rho/drdt(:,:,2:kmax+1))*Wsed(1,:,:,1:kmax)
!     & +0.5*((1.-dcdt(1,:,:,1:kmax))*rho_b/drdt(:,:,1:kmax)+(1.-dcdt(1,:,:,2:kmax+1))*rho_b/drdt(:,:,2:kmax+1))
!     & *(Wfluid(:,:,1:kmax)+dWdt_old(:,:,1:kmax))

!	  if (mod(time_np,5.).le.1.e-1) then
!		write(*,*),'dwdt,wfluid,wsed,sigma(w_k):',rank,dwdt(3,3,1),wfluid(3,3,1),wsed(1,3,3,1),
!     &		0.5*(dcdt(1,3,3,1)*frac(1)%rho/drdt(3,3,1)+dcdt(1,3,3,2)*frac(1)%rho/drdt(3,3,2))*Wsed(1,3,3,1)
!     &          +0.5*((1.-dcdt(1,3,3,1))*rho_b/drdt(3,3,1)+(1.-dcdt(1,3,3,2))*rho_b/drdt(3,3,2))*(Wfluid(3,3,1)+dWdt(3,3,1))
!	  endif
	  do n=1,nfrac
	      utr=dUdt
	      vtr=dVdt
	      wtr=Wsed(n,:,:,:)
	      call make_UtransC_zeroIBM(utr,vtr,wtr)
	      call advecc_TVD(k3c(n,:,:,:),dcdt(n,:,:,:),utr,vtr,wtr,drdt,Ru,Rp,dr,dphi,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn3*dt,rank,px,periodicx,periodicy)

!	      call advecc_TVD(k3c(n,:,:,:),dcdt(n,:,:,:),dUdt,dVdt,Wsed(n,:,:,:),drdt,Ru,Rp,dr,dphi,dz,
!     +            i1,j1,k1,ib,ie,jb,je,kb,ke,cn3*dt,rank,px,periodicx,periodicy)
	      call diffc_CDS2 (k3c(n,:,:,:),dCdt(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)

	          if (interaction_bed>0) then
		    !! advec concentration in bed with velocity TSHD:
	            call adveccbot_TVD(k3cbot(n,:,:),dcdt2bot(n,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,dphi,dz,
     +              i1,j1,ib,ie,jb,je,cn3*dt,rank,px,periodicx,periodicy)
		    dcdtbot(n,:,:)= dcdt2bot(n,:,:) + cn3*dt*k3cbot(n,:,:) ! n+1
	          endif
	  enddo
	  dcdt2=dcdt 
	  do n=1,nfrac
!		k3c(n,:,:,:) = dnewc(n,:,:,:)
!		dcdt(n,:,:,:) =cnew(n,:,:,:) + b1*dt*k1c(n,:,:,:) + b2*dt*k2c(n,:,:,:)+ b3*dt*k3c(n,:,:,:) !n+1
		dcdt(n,:,:,:) =dcdt(n,:,:,:) + cn3*dt*k3c(n,:,:,:) !n+1
!		CALL air_bubbles_free_surface ! added air_bubbles_free_surface here, just before bc are applied
!		9-10:2012: air_bubbles_free_surface not called because wsed is calculated on top row now, so air is dissappearing anyway
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc 
           IF (CNdiffz.eq.1) THEN !CN semi implicit treatment diff-z:
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
	           ekm_min=0.5*(Diffcof(i,j,km)+Diffcof(i,j,k))
	           ekm_plus=0.5*(Diffcof(i,j,kp)+Diffcof(i,j,k))
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
	     call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc after CN-diffz 
	   ENDIF
	  enddo
	  if (interaction_bed>0) then
	    call slipvelocity(dcdt,dwdt,wsed,drdt,0,0,wfluid,cn3*dt,dz) !get wfluid and wsed at k=0 for driftflux_force 
	    CALL erosion_deposition(dcdt,dcdtbot,dudt,dvdt,dwdt,drdt,dcdt2,dcdt2bot,cn3*dt,dz) !first two vars are adjusted
		    !! for kn3: extra var dcdt2 is used
	    ! should correct dwdt(n,:,:,0) with (1.-tau/frac(n)%tau_d), but this correction almost always is zero, for now leave
	    do n=1,nfrac
		call bound_cbot(dcdtbot(n,:,:))
	    enddo
	  endif
	endif

	if (diffusion.eq.'COM4') then
		call stress_terms(dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	endif

      	if (convection.eq.'CDS2') then
      call advecu_CDS2(k3u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(k3u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(k3u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(k3u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(k3u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(k3u,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (k3u,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(k3u,Srr,Spr,Szr,Spp,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      	if (convection.eq.'CDS2') then
      call advecv_CDS2(k3v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(k3v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(k3v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(k3v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(k3v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(k3v,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif
	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (k3v,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(k3v,Spr,Spp,Spz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif

      	if (convection.eq.'CDS2') then
      call advecw_CDS2(k3w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(k3w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(k3w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(k3w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(k3w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(k3w,dUdt,dVdt,dWdt,dRdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif
	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (k3w,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(k3w,Szr,Spz,Szz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	if ((slipvel.eq.1.or.slipvel.eq.2).and.nfrac>0) then
	 W_km_sum=0.
	 c_sum=0.
	 do n=1,nfrac
	  W_km_sum=W_km_sum+dcdt(n,:,:,:)*frac(n)%rho*(wsed(n,:,:,:)-dWdt)
	  c_sum=c_sum+dcdt(n,:,:,:)
	 enddo
	 W_km_sum=W_km_sum+(1.-c_sum)*rho_b*Wfluid  !Wfluid is difference with Ucfd
	 W_km_sum=W_km_sum/drdt
	 call advecw_driftfluxCDS2(k3w,0.,0.,W_km_sum,
     &   drdt,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
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
      call bound_rhoU(dUdt,dVdt,dWdt,dRdt,0,time_np,Ub1new,Vb1new,
     & Wb1new,Ub2new,Vb2new,Wb2new)
       do j=1,jmax
         do i=1,imax
            do k=0,k1
	    km=MAX(0,k-1)
	    kp=MIN(k1,k+1)
	    kpp=MIN(k1,k+2)
	    ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))
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
	    ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))
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
	    ekm_min=ekm(i,j,k)
	    ekm_plus=ekm(i,j,kp)
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
!      i = imax
!      do k=1,kmax
!       do j=1,jmax
!               dUdt1(i,j,k) = dUdt1(i,j,k) -c2*dt * ( -2* p(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
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
          dVdt(i,j,k) = dVdt(i,j,k) -dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * dphi )
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
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
      integer  ib,ie,jb,je,kb,ke,n,t
      real     doldc(nfrac,0:i1,0:j1,0:k1),dnewc(nfrac,0:i1,0:j1,0:k1)
      real     dold(0:i1,0:j1,0:k1),dnew(0:i1,0:j1,0:k1)
      real     wsed(nfrac,0:i1,0:j1,0:k1),wfluid(0:i1,0:j1,0:k1)
	real*8 pplus(imax,kmax)
	

!	real     Diffcof(0:i1,0:j1,0:k1)

!	Diffcof=ekm/Sc/Rnew 
c********************************************************************
c     CALCULATE slipvelocity
c********************************************************************
	if (nfrac>0) then
	  if (slipvel.eq.1.or.slipvel.eq.2) then
	      call slipvelocity(cnew,Wnew,wsed,rnew,1,kmax,wfluid,dt,dz)
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
	      call advecc_TVD(dnewc(n,:,:,:),cnew(n,:,:,:),Unew,Vnew,Wsed(n,:,:,:),rnew,Ru,Rp,dr,dphi,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
	      call diffc_CDS2 (dnewc(n,:,:,:),Cnew(n,:,:,:),Diffcof,
     +            ib,ie,jb,je,kb,ke)
		dcdt(n,:,:,:) =cnew(n,:,:,:) + dt*(dnewc(n,:,:,:)) !euler_expl
		call bound_c(dcdt(n,:,:,:),frac(n)%c,n) ! bc for intermediate dcdt
	  enddo
	  call state(dcdt,drdt) ! determine drdt
	endif

	if (diffusion.eq.'COM4') then
		call stress_terms(Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	endif

c********************************************************************
c     CALCULATE advection, diffusion and Force U-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecu_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecu_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecu_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecu_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecu_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecu_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffu_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffu_com4(dnew,Srr,Spr,Szr,Spp,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
            dUdt(i,j,k)=Unew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i+1,j,k))+
     1      dt*dnew(i,j,k)
     1      -(gx*cos_u(j)+gy*sin_u(j))*dt*(0.5*(rnew(i,j,k)+rnew(i+1,j,k))-rho_b)  
!     1     +Ppropx(i,j,k)/(MAX(1.e-2,(ABS(Unew(i,j,k)))))*dt 
     1     +Ppropx(i,j,k)*dt
            enddo
         enddo
      enddo
c********************************************************************
c     CALCULATE advection, diffusion and Force V-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecv_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecv_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecv_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecv_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecv_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecv_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffv_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffv_com4(dnew,Spr,Spp,Spz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
            dVdt(i,j,k)=Vnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j+1,k))+
     1      dt*dnew(i,j,k) 
     1      -(-gx*sin_v(j)+gy*cos_v(j))*dt*(0.5*(rnew(i,j,k)+rnew(i,j+1,k))-rho_b)
!     1     +Ppropy(i,j,k)/(MAX(1.e-2,(ABS(Vnew(i,j,k)))))*dt 
     1     +Ppropy(i,j,k)*dt
            enddo
         enddo
      enddo
c********************************************************************
c     CALCULATE advection, diffusion and Force W-velocity
c********************************************************************
      	if (convection.eq.'CDS2') then
      call advecw_CDS2(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'CDS6') then
      call advecw_CDS6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'COM4') then
      call advecw_COM4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	elseif(convection.eq.'HYB4') then
      call advecw_HYB4(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff)
	elseif(convection.eq.'HYB6') then
      call advecw_HYB6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	elseif(convection.eq.'C4A6') then
      call advecw_C4A6(dnew,Unew,Vnew,Wnew,Rnew,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdiff,periodicx,periodicy)
	endif

	if (diffusion.eq.'CDS2') then
      call diffw_CDS2 (dnew,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
	elseif(diffusion.eq.'COM4') then
      call diffw_com4(dnew,Szr,Spz,Szz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
	endif
	if (slipvel.eq.1.or.slipvel.eq.2) then
	do n=1,nfrac
	  call advecw_driftfluxCDS2(dnew,0.,0.,wsed(n,:,:,:)-Wnew,
     &  cnew(n,:,:,:)*frac(n)%rho,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke)
	enddo
	endif
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
            dWdt(i,j,k)=Wnew(i,j,k)*0.5*(rnew(i,j,k)+rnew(i,j,k+1))+
     1      dt*dnew(i,j,k) 
     1      -gz*dt*(0.5*(rnew(i,j,k)+rnew(i,j,k+1))-rho_b)
            enddo
         enddo
      enddo

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
!      i = imax  !! at imax do nothing -> dpdn=0 at outflow
!      do k=1,kmax
!	do j=1,jmax
!		dUdt(i,j,k) = dUdt(i,j,k) -dt * ( -2* pold(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
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
          dVdt(i,j,k) = dVdt(i,j,k) -dt * ( pold(i,j+1,k) - pold(i,j,k) ) /( Rp(i) * dphi )   
          enddo
        enddo
      enddo

	if (periodicy.eq.1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
                enddo
	      enddo
	else
      	  if (rank.ne.px-1) then
              do k=1,kmax
        	do i=1,imax
        	  dVdt(i,jmax,k) = dVdt(i,jmax,k) -dt * ( pplus(i,k) - pold(i,jmax,k) ) /( Rp(i) * dphi )
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

      subroutine fillps
      USE nlist
      implicit none

	real theta_U,xx,yy,r_orifice2,twodtdt
	integer t
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

	twodtdt=MAX((3.*time_np-4.*time_n+time_nm)*dt,1.e-12)

      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
c
      p(i,j,k)  =(
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz )
!     +              ) /dt  +  (3.*drdt(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/(2*dt*dt)
     +              ) /dt  +  (3.*drdt(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/twodtdt
          enddo
        enddo
      enddo
      return
      end

      subroutine fillps2(uu,vv,ww,rr,tt,ddtt)
      USE nlist
      implicit none

	real theta_U,xx,yy,r_orifice2,twodtdt
	real tt,ddtt,uu(0:i1,0:j1,0:k1),vv(0:i1,0:j1,0:k1),ww(0:i1,0:j1,0:k1),rr(0:i1,0:j1,0:k1)
	integer t
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

	twodtdt=MAX((3.*tt-4.*time_n+time_nm)*ddtt,1.e-12)
	
      do  k=1,kmax
        do j=1,jmax
          do i=1,imax
c
      p(i,j,k)  =(
     1  ( Ru(i)*uu(i,j,k) - Ru(i-1)*uu(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       vv(i,j,k) -         vv(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       ww(i,j,k) -         ww(i,j,k-1) ) / ( dz )
!     +              ) /dt  +  (3.*drdt(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/(2*dt*dt)
     +              ) /ddtt  +  (3.*rr(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/twodtdt
          enddo
        enddo
      enddo
      return
      end



      subroutine correc
      USE nlist
      implicit none
c
!       include 'param.txt'
!       include 'common.txt'
      real*8 pplus(imax,kmax)
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
!!! Do nothing on imax, so p(imax+1)=p(imax) --> dp/dn=0 at outflow
!      i = imax
!      do k=1,kmax
!	do j=1,jmax
!		dUdt(i,j,k) = dUdt(i,j,k) -  dt * ( -2* p(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
	if (periodicx.eq.1) then
	!!! !periodic bc in x: 
	      i = imax
	      do k=1,kmax
		do j=1,jmax
		  dUdt(i,j,k) = dUdt(i,j,k) -
     +             dt * ( p(1,j,k) - p(i,j,k) ) /
     +                ( Rp(i+1) - Rp(i) ) 
		 enddo
	      enddo
	endif



      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          dVdt(i,j,k) = dVdt(i,j,k) -
     +     dt * ( p(i,j+1,k) - p(i,j,k) ) /
     +                ( Rp(i) * dphi )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
	      do k=1,kmax
		do i=1,imax
		dVdt(i,jmax,k) = dVdt(i,jmax,k) -
     +   dt * ( pplus(i,k) - p(i,jmax,k) ) /
     +              ( Rp(i) * dphi )
		enddo
	      enddo
	else
	      if (rank.ne.px-1) then
 	       do k=1,kmax
	        do i=1,imax
	        dVdt(i,jmax,k) = dVdt(i,jmax,k) -
     +   dt * ( pplus(i,k) - p(i,jmax,k) ) /
     +              ( Rp(i) * dphi )
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

      do k=0,k1
        do j=0,j1
          do i=0,i1
           Uold(i,j,k)=Unew(i,j,k)
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

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
           Unew(i,j,k)=dUdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i+1,j,k)))
           Vnew(i,j,k)=dVdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i,j+1,k)))
           Wnew(i,j,k)=dWdt(i,j,k)/(0.5*(drdt(i,j,k)+drdt(i,j,k+1)))
          enddo
        enddo
      enddo 
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
!!! Do nothing on imax, so p(imax+1)=p(imax) --> dp/dn=0 at outflow
!      i = imax
!      do k=1,kmax
!	do j=1,jmax
!		uu(i,j,k) = uu(i,j,k) -  ddtt * ( -2* p(i,j,k) ) /    ( Rp(i+1) - Rp(i) )
!         enddo
!      enddo
	if (periodicx.eq.1) then
	!!! !periodic bc in x: 
	      i = imax
	      do k=1,kmax
		do j=1,jmax
		  uu(i,j,k) = uu(i,j,k) -
     +             ddtt * ( p(1,j,k) - p(i,j,k) ) /
     +                ( Rp(i+1) - Rp(i) ) 
		 enddo
	      enddo
	endif


      do k=1,kmax
        do j=1,jmax-1
          do i=1,imax
          vv(i,j,k) = vv(i,j,k) -
     +     ddtt * ( p(i,j+1,k) - p(i,j,k) ) /
     +                ( Rp(i) * dphi )   
          enddo
        enddo
      enddo
	if (periodicy.eq.1) then
	      do k=1,kmax
		do i=1,imax
		vv(i,jmax,k) = vv(i,jmax,k) -
     +   ddtt * ( pplus(i,k) - p(i,jmax,k) ) /
     +              ( Rp(i) * dphi )
		enddo
	      enddo
	else
	      if (rank.ne.px-1) then
 	       do k=1,kmax
	        do i=1,imax
	        vv(i,jmax,k) = vv(i,jmax,k) -
     +   ddtt * ( pplus(i,k) - p(i,jmax,k) ) /
     +              ( Rp(i) * dphi )
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
      integer   ipos,ierr
      real	rhs_ref
	real bbbbb(IMAX,JMAX,kmax)
      
c   generate tridiagonal systems
      pi = 4.*atan(1.)
      do i=1,imax
      a(i)= Ru(I-1)/((Rp(I)-Rp(I-1))*Rp(I)*(Ru(I)-Ru(I-1)))
      b(i)=-(Ru(I)/(Rp(I+1)-Rp(I))+Ru(I-1)/(Rp(I)-Rp(I-1)))/
     $      (Rp(I)*(Ru(I)-Ru(I-1)))
      c(i)= Ru(I) /((Rp(I+1)-Rp(I))*Rp(I)*(Ru(I)-Ru(I-1)))
      end do

	if (periodicx.eq.0) then
	      b(1)=   -(Ru(1)/(Rp(2)-Rp(1)))/(Rp(1)*(Ru(1)-Ru(0)))
	!      b(imax)=b(imax)-c(imax)     !! p(imax+1)=-p(imax) --> bc p=0
	!      b(imax)=b(imax)+c(imax)      !! p(imax+1)= p(imax) --> bc dpdn=0
      	      i=imax !! p(imax+1)= p(imax) --> bc dpdn=0
      	      b(i)=-(Ru(I-1)/(Rp(I)-Rp(I-1)))/   !! this is equivalent to b(imax)+c(imax) --> dpdn=0
     $        (Rp(I)*(Ru(I)-Ru(I-1)))
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

	if (periodicx.eq.0) then
 	  if (rank==0) then
	    bbbb(imax,1,1)=b(imax)-2.*c(imax)     !! p(imax+1)=-p(imax) --> bc p=0 at 1 loc in mesh (first +c(imax) so now -2c(imax))
	  endif
	  c(imax)=0.
	  a(1)=0.
	else 
 	  if (rank==0) then
	    bbbb(imax,1,1)=b(imax)-c(imax)     !! p(imax+1)=-p(imax) --> bc p=0 at 1 loc in mesh (first +c(imax) so now -2c(imax))
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
	        yrt(j)=(-4./(dphi*dphi))*(sin(float((j-1))*pi/(2.*jmax*px)))**2
      	      enddo 
		!   set up lookup tables
      	      call vcosqi(px*jmax,wj)
	elseif (periodicy.eq.1) then
	      do j=3,jmax*px,2
	        yrt(j-1)=(-4./(dphi*dphi))*(sin(float((j-1))*pi/(2.*jmax*px)))**2
		yrt(j)=yrt(j-1)
      	      enddo 
	      yrt(1)=0.
		if (mod(jmax*px,2).eq.0) then
		  yrt(jmax*px)=(-4./(dphi*dphi))
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
	      call vrfftf(imax*kmax/px,jmax*px,vfftj,rhs,imax*kmax/px,wj)
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

      if (periodicx.eq.0) then !  solve tridiagonal systems with Gaussian elemination
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
		    CALL solve_tridiag_periodic(rhs(1:imax,j,k),a,bbbbb(1:imax,j,k),c,rhs(1:imax,j,k),imax)
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
	      call vrfftb(imax*kmax/px,jmax*px,vfftj,rhs,imax*kmax/px,wj)
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
	

	!!    Make pressure zero at one point in outflow:
	!!    Only b(imax)-c(imax) in matrix is not sufficient to make pressure exactly zero in (imax,1,1)
	!!    Some small drift is occuring without following lines:
		if (rank.eq.0) then
		  rhs_ref=rhs(imax,1,1)
		endif
		call mpi_bcast(rhs_ref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
		rhs=rhs-rhs_ref

      return
      end
      
      
      subroutine state(c,rho)

      USE nlist
! 	include 'param.txt'
      real c(nfrac,0:i1,0:j1,0:k1)
      real rho(0:i1,0:j1,0:k1)
      integer n
	real xTSHD(4),yTSHD(4),phi	  
      
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
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
      !! Search for P,V:
      do k=0,k1
       do i=0,i1  
         do j=j1,0,-1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  IF (k.le.FLOOR(bp(n2)%height/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
		rnew(i,j,k)=rho(i,j,k)
		rold(i,j,k)=rho(i,j,k)
		! prevent large source in pres-corr by sudden increase in density
	   endif
	  enddo
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
	integer nr,mt,nx,mx,nt
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
	integer nr,mt,nx,mx,nt
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
	integer nr,mt,nx,mx,nt
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

      SUBROUTINE t2np_0ie(Q,Qt)
      USE nlist
!         INCLUDE 'param.txt'
        INCLUDE 'mpif.h'
!         parameter (nr=imax,mt=jmax,nx=kmax,mx=kmax/px,NT=jmax*px)
	integer nr,mt,nx,mx,nt
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
