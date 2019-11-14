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

      subroutine chkdt

      USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      include 'mpif.h'
      integer ierr
      real  dr2,dz2,df,df2,kcoeff,tmp1,tmp2,tmp3,Courant,dtmp,dtold
	  double precision Tadapt,Uav,Vav,localsum,globalsum,du,duu,wsettlingmax,wsettlingmin

	  
!		IF (nobst>0.and.bp(1)%forever.eq.0.and.time_np+dt.gt.bp(1)%t0.and.counter<10) THEN
!			counter=counter+1
!			Courant=0.1
!			IF (rank.eq.0) THEN
!				write(*,*),'Placement bedplume, CFL=0.1 for 10 time steps ',counter,dt
!			ENDIF
!		ELSE 
			Courant = CFL
!		ENDIF
		
		dt_old=dt
		dtold=dt
      IF (istep.le.10) THEN
	dt = dt_ini
      ELSE
	dt = dt_max
      ENDIF
	  
!		write(*,*),'dt,rank voor',dt,rank
      IF (CNdiffz.eq.1) THEN
	dz2 = 1.e9*dz*dz ! no dt restriction for vertical diff with CN implicit scheme
      ELSE
        dz2 = dz    * dz
      ENDIF
	  IF (nfrac>0) THEN
		wsettlingmax=MAXVAL(frac(1:nfrac)%ws) !positive downward
		wsettlingmin=MIN(0.,MINVAL(frac(1:nfrac)%ws)) !find air rise velocity (otherwise zero)
	  ELSE
		wsettlingmax=0.
		wsettlingmin=0.
	  ENDIF
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
            df = Rp(i)*(phiv(j)-phiv(j-1))
      	    df2 = df * df
            dr2 = dr(i) * dr(i)
            kcoeff = ekm(i,j,k) /rnew(i,j,k)
            tmp1 = ( abs(Unew(i,j,k)) / ( Rp(i+1)-Rp(i) ) ) +
     &             ( abs(Vnew(i,j,k)) /         df        ) +
     &             ( MAX(abs(Wnew(i,j,k)-wsettlingmax),abs(Wnew(i,j,k)-wsettlingmin)) /         dz        )             
!			tmp1 = ( abs(Unew(i,j,k)) / ( Rp(i+1)-Rp(i) ) )              !2D TVD tests show that taking MAX of 3 dirs is not sufficient for stable results without overshoot/undershoot
!			tmp1 = MAX(tmp1,( abs(Vnew(i,j,k)) /         df        ))
!			tmp1 = MAX(tmp1,( abs(Wnew(i,j,k)) /         dz        ))
            tmp2 = ( 1.0/dr2 + 1.0/df2 + 1.0/dz2 )
            tmp3 = 1.0 / ( 1.0 * tmp2 * kcoeff + tmp1 + 1.e-12 )
            tmp3 =Courant *tmp3 
            dt = min( dt , tmp3 )
            dtmp = dt

            enddo
         enddo
      enddo

!	call mpi_allreduce(dtmp,dt,1,mpi_real8,mpi_min,mpi_comm_world,ierr)
	!call mpi_allreduce(dtmp,dt,1,mpi_real,mpi_min,mpi_comm_world,ierr)
	call mpi_allreduce(dtmp,dt,1,mpi_double_precision,mpi_min,mpi_comm_world,ierr)	

	  dt=MIN(dt,dtold*1.1) !change in time step never more than +10% in one timestep for extra stability

	  
	if (isnan(dt)) stop 'ERROR, QUITING TUDFLOW3D: "dt" is a NaN'
	if (dt<1.e-12) stop 'ERROR, QUITING TUDFLOW3D: "dt" is smaller than 1e-12'

	if (time_int.eq.'AB2'.or.time_int.eq.'AB3') then 
	  if (dt<dt_max) then
		stop 'ERROR, QUITING TUDFLOW3D: "dt" is adjusted in AB2 or AB3, sign of instability try again with smaller dt'
	  endif
	endif
!! 	if AB2 or AB3 timestep must be fixed, if unfortunately initially a too large timestep has been chosen 
!! 	then a 25% reduction must lead to stable simulation from this point on (one change in timestep can be 
!!	stable, but subsequent changes in timestep are not stable for AB2 or AB3)
!	if (time_int.eq.'AB2'.or.time_int.eq.'AB3') then 
!	  if (dt<dt_max) then  
!		WRITE(*,'(a,a,a,f12.9,a,f12.9,a)') ' WARNING: Using ',time_int,' with dt too large, dt changed from : '
!     &           ,dt_max,' s, to : ',0.75*dt_max,' s'  
!	    dt_max=0.75*dt_max
!	    dt=dt_max
!	  endif
!	endif
	  if (periodicx.eq.1.and.ABS(dpdx).lt.1.e-12) THEN ! test determination correct dpdx and dpdy based on wanted u or v
	    localsum=0.
        do k=1,kmax
          do j=1,jmax
             do i=1,imax	  
		       localsum=localsum+Unew(i,j,k)*(Rp(i+1)-Rp(i))*(phiv(j)-phiv(j-1))*Ru(i)*dz
			 enddo
		  enddo
	    enddo
		call mpi_allreduce(localsum,globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		
		Uav = globalsum/(SUM(Vol_V)*REAL(kmax)) 
		if (istep>1) THEN
			du = Uav-Uavold
		else
			du = 0.
		endif
		Uavold = Uav
		Tadapt = 100.*depth/sqrt(U_b**2+V_b**2+1.e-18) !for test just take a value
		duu = Tadapt*du/dt !30 seconds worked
		Uav = Uav + duu 
		dpdx1 = dpdx1 + ABS(U_b)*(U_b-Uav)/depth*dt/Tadapt
		!dpdx1 = dpdx1 + MIN(ABS(ABS(U_b)*(U_b-Uav)/depth*dt/Tadapt),0.1*ABS(dpdx1))*SIGN(1.,ABS(U_b)*(U_b-Uav)/depth*dt/Tadapt)
		Ppropx = dpdx1*rhU !rnew !with variable density important to use rnew and not rho_b 
		Uav = Uav - duu
		if (rank.eq.0) then
			if (mod(istep,10) .eq.0) then   
				write(*,'(a,f8.2,a,f6.2,a,f6.2,a)') ' # dpdx: ',dpdx1*rho_b,' Pa/m; Uav: ',Uav,'; U_b:',U_b,' m/s #'
     		    !write(*,*),du,duu,ABS(U_b)*(U_b-(Uav+duu))/depth*dt/Tadapt
			endif  
		endif			
	  endif
	  if (periodicy.eq.1.and.ABS(dpdy).lt.1.e-12) THEN ! test determination correct dpdx and dpdy based on wanted u or v
	    localsum=0.
        do k=1,kmax
          do j=1,jmax
             do i=1,imax	  
		       localsum=localsum+Vnew(i,j,k)*(Ru(i)-Ru(i-1))*(phip(j+1)-phip(j))*Rp(i)*dz
			 enddo
		  enddo
	    enddo
		call mpi_allreduce(localsum,globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)	
		Vav = globalsum/(SUM(Vol_V)*REAL(kmax)) 
		if (istep>1) THEN
			du = Vav-Vavold
		else
			du = 0.
		endif
		Vavold = Vav
		Tadapt = 100.*depth/sqrt(U_b**2+V_b**2+1.e-18) !for test just take a value
		duu = Tadapt*du/dt !30 seconds worked
		Vav = Vav + duu 
		dpdy1 = dpdy1 + ABS(V_b)*(V_b-Vav)/depth*dt/Tadapt 
		!dpdy1 = dpdy1 + MIN(ABS(ABS(V_b)*(V_b-Vav)/depth*dt/Tadapt),0.1*ABS(dpdy1))*SIGN(1.,ABS(V_b)*(V_b-Vav)/depth*dt/Tadapt)
		Ppropy = dpdy1*rhV !rnew !with variable density important to use rnew and not rho_b
		Vav = Vav - duu
		if (rank.eq.0) then
			if (mod(istep,10) .eq.0) then   
				write(*,'(a,f8.2,a,f6.2,a,f6.2,a)') ' # dpdy: ',dpdy1*rho_b,' Pa/m; Vav: ',Vav,'; V_b:',V_b,' m/s #'		
			endif  
		endif	
	   endif		

      end

      subroutine chkdiv
      USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      include 'mpif.h'
      real   div1,divmax1,divbar1,divmax_tot1,divbar_tot1,rhoip,rhoim,rhojp,rhojm,rhokp,rhokm
      real   div2,divmax2,divbar2,divmax_tot2,divbar_tot2
	  real   div3,divmax3,divbar3,divmax_tot3,divbar_tot3,rrri,rrrj,rrrk,rrrim,rrrjm,rrrkm
      integer ierr,n,inout,n2
	  real xTSHD(1:4),yTSHD(1:4),phi,xx,yy,kg_sus,kg_sed,kg_sus_tot,kg_sed_tot
	  real cU2(nfrac,0:i1,0:j1,0:k1),cV2(nfrac,0:i1,0:j1,0:k1),cW2(nfrac,0:i1,0:j1,0:k1)
	  real rhU2(0:i1,0:j1,0:k1),rhV2(0:i1,0:j1,0:k1),rhW2(0:i1,0:j1,0:k1)
	  
      divbar1 = 0.0
      divmax1 = 0.0
      divbar2 = 0.0
      divmax2 = 0.0
      divbar3 = 0.0
      divmax3 = 0.0	  
	  
	  
      do k=2,kmax-1
         do j=2,jmax-1
            do i=2,imax-1		
		div1 =
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz )
!     +              +  (3.*drdt(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/(2.*dt)
     +              +  (3.*drdt(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/((3.*time_np-4.*time_n+time_nm))

	    div2= ( Ru(i)*Unew(i,j,k) - Ru(i-1)*Unew(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vnew(i,j,k) -         Vnew(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       Wnew(i,j,k) -         Wnew(i,j,k-1) ) / ( dz )  
	 
	DO n2=1,nbedplume !correct div2 for adding or subtracting volume:
		IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end.and.bp(n2)%Q.ne.0.)) THEN
		! rotation ship for ambient side current
		if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
		  xx=Rp(i)*cos_u(j)-schuif_x
		  yy=Rp(i)*sin_u(j)
		  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
			xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
			yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
		    if (bp(n2)%radius.gt.0.) then
			  inout=0
		      IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			    inout=1
			  ENDIF
			else 
			  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
			endif 			
		  ELSE 
			inout=0
		  ENDIF
		  if (inout.eq.1) then
			div2=div2-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells
		   endif
		ENDIF
	ENDDO ! bedplume loop
	

		IF (split_rho_cont.eq.'VL2'.or.split_rho_cont.eq.'SB2') then
		! U,V,W,C has been updated already, but rho not!
		  IF (split_rho_cont.eq.'VL2') then
			do n=1,nfrac
			  call c_edges_VL2_nocfl(cU2(n,:,:,:),cV2(n,:,:,:),cW2(n,:,:,:),cold(n,:,:,:),Uold,Vold,Wold,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,1,imax,1,jmax,1,kmax,dt,rank,px,periodicx,periodicy)
			enddo
		  ENDIF
		  IF (split_rho_cont.eq.'SB2') then
			do n=1,nfrac
			  call c_edges_SB2_nocfl(cU2(n,:,:,:),cV2(n,:,:,:),cW2(n,:,:,:),cold(n,:,:,:),Uold,Vold,Wold,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,1,imax,1,jmax,1,kmax,dt,rank,px,periodicx,periodicy)
			enddo
		  ENDIF		  
			call state_edges(cU2,rhU2)
			call state_edges(cV2,rhV2)
			call state_edges(cW2,rhW2)		  
		   rrri= rhU2(i,j,k)
		   rrrj= rhV2(i,j,k)
		   rrrk= rhW2(i,j,k)
		   rrrim=rhU2(i-1,j,k)
		   rrrjm=rhV2(i,j-1,k)
		   rrrkm=rhW2(i,j,k-1)
		ELSE
		   rrri= 0.5*(rnew(i,j,k)+rnew(i+1,j,k))
		   rrrj= 0.5*(rnew(i,j,k)+rnew(i,j+1,k))
		   rrrk= 0.5*(rnew(i,j,k)+rnew(i,j,k+1))
		   rrrim=0.5*(rnew(i,j,k)+rnew(i-1,j,k))
		   rrrjm=0.5*(rnew(i,j,k)+rnew(i,j-1,k))
		   rrrkm=0.5*(rnew(i,j,k)+rnew(i,j,k-1))		
		ENDIF
	 
			div3 =
     1  ( Ru(i)*Uold(i,j,k)*rrri - Ru(i-1)*Uold(i-1,j,k)*rrrim ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vold(i,j,k)*rrrj -         Vold(i,j-1,k)*rrrjm ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       Wold(i,j,k)*rrrk -         Wold(i,j,k-1)*rrrkm ) / ( dz )
     +              +  (drdt(i,j,k)-rnew(i,j,k))/(dt)
	
      divbar1 = divbar1 + div1
      div1    = abs(div1)
      divmax1 = max( divmax1 , div1 )
      divbar2 = divbar2 + div2
      div2    = abs(div2)
      divmax2 = max( divmax2 , div2 )
      divbar3 = divbar3 + div3
      div3    = abs(div3)
      divmax3 = max( divmax3 , div3 )
           enddo
         enddo
      enddo

!		kg_sed=0.
!		kg_sus=0.
!      do  i=1,imax
!        do j=1,jmax
!		  do n=1,nfrac
!            do k=1,kmax
!				kg_sus=kg_sus+Cnew(n,i,j,k)*dr(i)*Rp(i)*(phiv(j)-phiv(j-1))*dz*frac(n)%rho
!				kg_sed=kg_sed+Clivebed(n,i,j,k)*dr(i)*Rp(i)*(phiv(j)-phiv(j-1))*dz*frac(n)%rho
!            enddo
!			kg_sed=kg_sed+cnewbot(n,i,j)*dr(i)*Rp(i)*(phiv(j)-phiv(j-1))*dz*frac(n)%rho
!		  enddo
!        enddo
!      enddo		  
!		call mpi_allreduce(kg_sus,kg_sus_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
!		call mpi_allreduce(kg_sed,kg_sed_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
!      if (rank.eq.0) write(6,104)kg_sus_tot,kg_sed_tot,
!     & (kg_sus_tot+kg_sed_tot)/(time_np*(bp(1)%sedflux(1)+bp(1)%sedflux(2)+bp(1)%sedflux(3)+bp(1)%sedflux(4)))  
!104   format('kg_sus = ',e13.6,'kg_sed = ',e13.6,'(kg_sus+kg_sed)/influx=',e13.6)
	 
	  
      call mpi_allreduce(divbar1,divbar_tot1,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(divmax1,divmax_tot1,1,mpi_real8,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(divbar2,divbar_tot2,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(divmax2,divmax_tot2,1,mpi_real8,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(divbar3,divbar_tot3,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(divmax3,divmax_tot3,1,mpi_real8,mpi_max,mpi_comm_world,ierr)	  
      !call mpi_allreduce(divbar,divbar_tot,1,mpi_real,mpi_max,mpi_comm_world,ierr)
      !call mpi_allreduce(divmax,divmax_tot,1,mpi_real,mpi_sum,mpi_comm_world,ierr)
      
	if (rank.eq.0) write(6,100)divbar_tot1,divmax_tot1   
100   format('Mass loss/gain A : Tot = ',e13.6,
     +                      '  Max = ',e13.6)
	if (rank.eq.0) write(6,102)divbar_tot3,divmax_tot3 
102   format('Mass loss/gain B : Tot = ',e13.6,
     +                      '  Max = ',e13.6)	
	 if (rank.eq.0) write(6,101)divbar_tot2,divmax_tot2  
101   format('Div(u)           : Tot = ',e13.6,
     +                      '  Max = ',e13.6)
 
      end
      
      
!      subroutine calc_div
!      USE nlist
!      implicit none

!      real dz_i,Rpdr_i,Rpdphi_i

!      dz_i=1./dz
!      do i=1,imax
!	Rpdr_i=1./(Rp(i)*dr(i))
!	Rpdphi_i=1./(Rp(i)*dphi)
!	do j=1,jmax
!	  do k=1,kmax

!	    div(i,j,k)= ( Ru(i)*Unew(i,j,k) - Ru(i-1)*Unew(i-1,j,k) ) *Rpdr_i ! / ( Rp(i)*dr(i) )
!     +              +
!     2  (       Vnew(i,j,k) -         Vnew(i,j-1,k) ) * Rpdphi_i !/ ( Rp(i)*dphi )
!     +              +
!     3  (       Wnew(i,j,k) -         Wnew(i,j,k-1) ) * dz_i !/ ( dz )  ) 
!	  enddo
!	enddo
!      enddo
!      end

!	MODULE work_array
!        REAL, DIMENSION(:,:,:), ALLOCATABLE :: Uavg,Vavg,Wavg,Ravg,Pavg,muavg
!	  REAL, DIMENSION(:,:,:), ALLOCATABLE :: sigU2,sigV2,sigW2,sigR2,sigUV,sigUW,sigVW
!	  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: sigC2,Cavg,sigUC,sigVC,sigWC,Cmax,Cmin
!	  INTEGER stat_count
!	REAL stat_time_count
!!        REAL, DIMENSION(:,:,:), ALLOCATABLE :: fUavg,sigfU2
!      END MODULE
      


	
	subroutine statistics
!	USE work_array
	USE nlist
	implicit none
! 	include 'param.txt'
! 	include 'common.txt'

	!REAL uu,vv,ww,uudt,vvdt,wwdt
	REAL uu(1:imax,1:jmax,1:kmax),vv(1:imax,1:jmax,1:kmax),ww(1:imax,1:jmax,1:kmax)
	REAL uudt(1:imax,1:jmax,1:kmax),vvdt(1:imax,1:jmax,1:kmax),wwdt(1:imax,1:jmax,1:kmax)
	INTEGER n

	stat_count=stat_count+1
	stat_time_count=stat_time_count+dt
	
	! below is 3x faster than elementwise within filling AVG and sig2 arrays within do-loop: now RMS takes 10% instead of 30% CPU time!
		do i=1,imax
		  do j=1,jmax
			do k=1,kmax
			  uu(i,j,k)=0.5*(Unew(i,j,k)+Unew(i-1,j,k))*cos_u(j)-0.5*(Vnew(i,j,k)+Vnew(i,j-1,k))*sin_u(j)
			  vv(i,j,k)=0.5*(Vnew(i,j,k)+Vnew(i,j-1,k))*cos_u(j)+0.5*(Unew(i,j,k)+Unew(i-1,j,k))*sin_u(j)
			  ww(i,j,k)=0.5*(Wnew(i,j,k)+Wnew(i,j,k-1))
			  do n=1,nfrac
				Cmax(n,i,j,k) = MAX(Cmax(n,i,j,k),Cnew(n,i,j,k))
				Cmin(n,i,j,k) = MIN(Cmin(n,i,j,k),Cnew(n,i,j,k))
			  enddo		  
			enddo
		  enddo
		enddo		  
		  uudt=uu*dt
		  vvdt=vv*dt
		  wwdt=ww*dt

		  sigU2 = sigU2 + uu*uudt
		  sigV2 = sigV2 + vv*vvdt
		  sigW2 = sigW2 + ww*wwdt
		  sigR2 = sigR2 + Rnew(1:imax,1:jmax,1:kmax)*Rnew(1:imax,1:jmax,1:kmax)*dt
		  sigUV = sigUV + uu*vvdt
		  sigUW = sigUW + uu*wwdt
		  sigVW = sigVW + vv*wwdt

		  Uavg  = Uavg + uudt
		  Vavg  = Vavg + vvdt
		  Wavg  = Wavg + wwdt
		  Umax  = MAX(Umax,ABS(uu))
		  Vmax  = MAX(Vmax,ABS(vv))
		  Wmax  = MAX(Wmax,ABS(ww))
!		  Umin  = MIN(Umax,ABS(uu))
!		  Vmin  = MIN(Vmax,ABS(vv))
!		  Wmin  = MIN(Wmax,ABS(ww))	
		  Uhormax  = MAX(Uhormax,SQRT(uu**2+vv**2))
		  U3dmax  = MAX(U3dmax,SQRT(uu**2+vv**2+ww**2))
		  Ravg  = Ravg + Rnew(1:imax,1:jmax,1:kmax)*dt
		  Pavg  = Pavg + (p(1:imax,1:jmax,1:kmax)+pold(1:imax,1:jmax,1:kmax))*dt
		  muavg = muavg + ekm(1:imax,1:jmax,1:kmax)*dt
		  Cavg  = Cavg + Cnew(1:nfrac,1:imax,1:jmax,1:kmax)*dt
		  do n=1,nfrac
			  sigC2(n,:,:,:) = sigC2(n,:,:,:) + Cnew(n,1:imax,1:jmax,1:kmax)*Cnew(n,1:imax,1:jmax,1:kmax)*dt
			  sigUC(n,:,:,:) = sigUC(n,:,:,:) + Cnew(n,1:imax,1:jmax,1:kmax)*uudt
			  sigVC(n,:,:,:) = sigVC(n,:,:,:) + Cnew(n,1:imax,1:jmax,1:kmax)*vvdt
			  sigWC(n,:,:,:) = sigWC(n,:,:,:) + Cnew(n,1:imax,1:jmax,1:kmax)*wwdt
		  enddo
		
!!! Favre average results in differences of <1 promille for Umean and <1% for U' between Favre and Reynolds avg
!!! for CO2 plume simulation (Wang et al. 2008), this is not important for dredge plumes, thus only Reynolds avg is calculated from now on
!		  fUavg(i,j,k) = fUavg(i,j,k) + uu*R(i,j,k)
!		  sigfU2(i,j,k) = sigfU2(i,j,k) + uu*uu*R(i,j,k)
	end

	MODULE history_array
	use nlist
	IMPLICIT NONE

	INTEGER,PARAMETER :: tdim=200000 ! length of history array in the time
	INTEGER :: nhispoint
	
	type history
	    integer ::	i,j,k
	    REAL ::	x,y,z
	    REAL*8, DIMENSION(tdim) :: U,V,W,P,RHO,t
	    REAL*8, DIMENSION(20,tdim) :: C ! 20 fractions in history, must be enough..
	end type history

	TYPE(history), DIMENSION(:), ALLOCATABLE :: his
      END MODULE	
	
	subroutine init_his
	USE history_array
	USE error_functions
	USE nlist
	implicit none
	
!       include 'param.txt'
!       include 'common.txt'

	integer ios,n

	type history_init
	  integer :: i,j,k
	end type history_init

	TYPE(history_init), DIMENSION(100) :: hist

	NAMELIST /histories/hist
	
	hist(:)%i=0
	OPEN(1,FILE=hisfile,IOSTAT=ios,ACTION='read')
	IF (ios/=0) CALL writeerror(100)

	READ (UNIT=1,NML=histories,IOSTAT=ios)

	nhispoint=0
	DO WHILE (hist(nhispoint+1)%i.NE.0)
	  nhispoint=nhispoint+1
	END DO

	!write(*,*)'hist1,2,3,4',hist(1)%i,hist(2)%i,hist(3)%i,hist(4)%i
	!write(*,*)'nhispoint',nhispoint
	CLOSE(1)

	ALLOCATE(his(nhispoint))
	!write(*,*)'his,2,3,4',his(1)%i,his(2)%i,his(3)%i,his(4)%i
	
	DO n=1,nhispoint
	  his(n)%i=hist(n)%i
	  his(n)%j=hist(n)%j
	  his(n)%k=hist(n)%k

	  IF (his(n)%i>imax) CALL writeerror(110)
	  IF (his(n)%i<0) CALL writeerror(111)
	  IF (his(n)%j>jmax*px) CALL writeerror(120)
	  IF (his(n)%j<0) CALL writeerror(121)
	  IF (his(n)%k>kmax) CALL writeerror(130)
	  IF (his(n)%k<0) CALL writeerror(131)

	  his(n)%x=Rp(his(n)%i)*cos_u(his(n)%j-rank*jmax)-schuif_x
	  his(n)%y=Rp(his(n)%i)*sin_u(his(n)%j-rank*jmax)
	  his(n)%z=his(n)%k*dz-0.5*dz
	ENDDO
	    

	END SUBROUTINE init_his

! 	subroutine appendhis(rank,istep,time,Unew,Vnew,Wnew,P,Cnew,Rnew)
	subroutine appendhis
	USE history_array
	USE nlist

	implicit none


!       include 'param.txt'
      !include 'common.txt'


! 	REAL Unew(0:i1,0:j1,0:k1),Vnew(0:i1,0:j1,0:k1),Wnew(0:i1,0:j1,0:k1)
! 	REAL Rnew(0:i1,0:j1,0:k1),Cnew(0:i1,0:j1,0:k1),P(0:i1,0:j1,0:k1)
! 	REAL time

	integer n,r,nf
	
	IF (istep.lt.200000) THEN
	his(1)%t(istep)=time_np !Unew is already updated, therefore it is output at time_np
	DO n=1,nhispoint
	  r=INT((his(n)%j-1)/jmax)
	  IF (rank.eq.r) THEN
	    i=his(n)%i
	    j=his(n)%j-rank*jmax
	    k=his(n)%k
	    his(n)%U(istep)=Unew(i,j,k)*cos_u(j)-Vnew(i,j,k)*sin_v(j)
	    his(n)%V(istep)=Vnew(i,j,k)*cos_v(j)+Unew(i,j,k)*sin_u(j)
	    his(n)%W(istep)=Wnew(i,j,k)
	    his(n)%P(istep)=Pold(i,j,k)+P(i,j,k)
	    his(n)%RHO(istep)=Rnew(i,j,k)
	    DO nf=1,MIN(20,nfrac)
	    	his(n)%C(nf,istep)=Cnew(nf,i,j,k)
	    ENDDO
	  ENDIF
	ENDDO
	ENDIF
	end SUBROUTINE

! 	subroutine finalize_his(rank,istep)
	subroutine finalize_his
	USE history_array
	USE nlist
        USE netcdf
	implicit none

!       include 'param.txt'
      include 'mpif.h'
!	real UU(1:200000)
	integer ierr,n,r,tag,status(MPI_STATUS_SIZE),nf
	INTEGER (8) :: ps
	REAL :: Uhis(1:nhispoint,1:istep),Vhis(1:nhispoint,1:istep),Whis(1:nhispoint,1:istep)
	REAL :: Phis(1:nhispoint,1:istep),Chis(1:nfrac,1:nhispoint,1:istep),RHOhis(1:nhispoint,1:istep)
	REAL :: xhis(1:nhispoint),yhis(1:nhispoint),zhis(1:nhispoint)
	INTEGER :: ihis(1:nhispoint),jhis(1:nhispoint),khis(1:nhispoint)
	REAL :: this(1:istep)

       ! We are writing 3D, 2D and 1D data, a nhispoint x istep grid.
       integer, parameter :: NDIMS1 = 2
       integer, parameter :: NDIMS2 = 1
       integer, parameter :: NDIMS3 = 3
       integer, parameter :: NDIMS4 = 1
!       integer, parameter :: NX = imax, NY = jmax, NZ = kmax
     
       ! When we create netCDF files, variables and dimensions, we get back
       ! an ID for each one.
       integer :: ncid, varid1,varid2,varid3, varid4, varid5, varid6, varid7, varid8 
       integer :: varid9,varid10,varid11,varid12,varid13
       integer :: dimids1(NDIMS1), dimids2(NDIMS2),dimids3(NDIMS3),dimids4(NDIMS4)
       integer :: nhis_dimid,time_dimid,nfrac_dimid,istep2
	character(1024) :: svnversion
	character(1024) :: svnurl
      include 'version.inc'


	IF (rank>0) THEN
	  DO n=1,nhispoint
	    r=INT((his(n)%j-1)/jmax)
	    IF (rank.eq.r) THEN
	      tag=INT(r*n)
	      call mpi_sendrecv_replace(his(n)%U ,200000,MPI_REAL8,0,10,0,10, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%V ,200000,MPI_REAL8,0,11,0,11, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%W ,200000,MPI_REAL8,0,12,0,12, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%P ,200000,MPI_REAL8,0,13,0,13, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%C,20*200000,MPI_REAL8,0,14,0,14, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%RHO ,200000,MPI_REAL8,0,15,0,15, MPI_COMM_WORLD,status,ierr)
 	      call mpi_sendrecv_replace(his(n)%x ,1,MPI_REAL8,0,16,0,16, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%y ,1,MPI_REAL8,0,17,0,17, MPI_COMM_WORLD,status,ierr)	      
	      call mpi_sendrecv_replace(his(n)%z ,1,MPI_REAL8,0,18,0,18, MPI_COMM_WORLD,status,ierr)
	    ENDIF
	  ENDDO
	ELSE
	  DO n=1,nhispoint
	    r=INT((his(n)%j-1)/jmax)
	    IF (rank.ne.r) THEN
 	      call mpi_sendrecv_replace(his(n)%U ,200000,MPI_REAL8,r,10,r,10, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%V ,200000,MPI_REAL8,r,11,r,11, MPI_COMM_WORLD,status,ierr)	      
	      call mpi_sendrecv_replace(his(n)%W ,200000,MPI_REAL8,r,12,r,12, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%P ,200000,MPI_REAL8,r,13,r,13, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%C ,20*200000,MPI_REAL8,r,14,r,14, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%RHO ,200000,MPI_REAL8,r,15,r,15, MPI_COMM_WORLD,status,ierr)
 	      call mpi_sendrecv_replace(his(n)%x ,1,MPI_REAL8,r,16,r,16, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%y ,1,MPI_REAL8,r,17,r,17, MPI_COMM_WORLD,status,ierr)	      
	      call mpi_sendrecv_replace(his(n)%z ,1,MPI_REAL8,r,18,r,18, MPI_COMM_WORLD,status,ierr)
	    ENDIF
	  ENDDO
	ENDIF

	istep2=MIN(istep,200000) ! array size is max 200000
	IF (rank.eq.0) THEN
	  DO n=1,nhispoint
	    Uhis(n,1:istep2)=his(n)%U(1:istep2)
	    Vhis(n,1:istep2)=his(n)%V(1:istep2)
	    Whis(n,1:istep2)=his(n)%W(1:istep2)
	    DO k=1,MIN(nfrac,20)
	      Chis(k,n,1:istep2)=his(n)%C(k,1:istep2)
	    ENDDO
	    Phis(n,1:istep2)=his(n)%P(1:istep2)
	    RHOhis(n,1:istep2)=his(n)%RHO(1:istep2)
		
	    xhis(n)=his(n)%x
	    yhis(n)=his(n)%y
	    zhis(n)=his(n)%z
	    ihis(n)=his(n)%i
	    jhis(n)=his(n)%j
	    khis(n)=his(n)%k

	  ENDDO
	this(1:istep)=his(1)%t(1:istep)

       ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
       ! overwrite this file, if it already exists.
       call check( nf90_create('history.nc', NF90_CLOBBER, ncid) )
     
       ! Define the dimensions. NetCDF will hand back an ID for each.
       call check( nf90_def_dim(ncid, "nhis_dim", nhispoint, nhis_dimid) )
       call check( nf90_def_dim(ncid, "time_dim", istep, time_dimid) )
	if (nfrac>0) then
       call check( nf90_def_dim(ncid, "nfrac_dim", nfrac, nfrac_dimid) )
	endif
     
       ! The dimids array is used to pass the IDs of the dimensions of
       ! the variables. Note that in fortran arrays are stored in
       ! column-major format.
       dimids1 =  (/ nhis_dimid, time_dimid/)
       dimids2 =  (/ nhis_dimid /)
	if (nfrac>0) then
       dimids3 =  (/ nfrac_dimid, nhis_dimid, time_dimid/)
	endif
       dimids4 =  (/ time_dimid/)
    
       ! Define the variable. The type of the variable in this case is
       ! NF90_DOUBLE (4-byte double).
       call check( nf90_def_var(ncid, "Uhis", NF90_DOUBLE, dimids1, varid1) )
       call check( nf90_put_att(ncid, varid1, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid1, 'long_name', 'Time history of U velocity at specific history location') )

       call check( nf90_def_var(ncid, "Vhis", NF90_DOUBLE, dimids1, varid2) )
       call check( nf90_put_att(ncid, varid2, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid2, 'long_name', 'Time history of V velocity at specific history location') )

       call check( nf90_def_var(ncid, "Whis", NF90_DOUBLE, dimids1, varid3) )
       call check( nf90_put_att(ncid, varid3, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid3, 'long_name', 'Time history of W velocity at specific history location') )

       call check( nf90_def_var(ncid, "RHOhis", NF90_DOUBLE, dimids1, varid4) )
       call check( nf90_put_att(ncid, varid4, 'units', 'kg/m3') )
       call check( nf90_put_att(ncid, varid4, 'long_name', 'Time history of density at specific history location') )

       call check( nf90_def_var(ncid, "Phis", NF90_DOUBLE, dimids1, varid5) )
       call check( nf90_put_att(ncid, varid5, 'units', 'Pa') )
       call check( nf90_put_att(ncid, varid5, 'long_name', 'Time history of pressure at specific history location') )

       if (nfrac>0) then
       call check( nf90_def_var(ncid, "Chis", NF90_DOUBLE, dimids3, varid6) )
       call check( nf90_put_att(ncid, varid6, 'units', '-') )
       call check( nf90_put_att(ncid, varid6, 'long_name', 'Time history of volume concentration 
     & of each fraction at specific history location') )
	endif

       call check( nf90_def_var(ncid, "xhis", NF90_DOUBLE, dimids2, varid7) )
       call check( nf90_put_att(ncid, varid7, 'units', 'm') )
       call check( nf90_put_att(ncid, varid7, 'long_name', 'X coordinate of specific history location') )

       call check( nf90_def_var(ncid, "yhis", NF90_DOUBLE, dimids2, varid8) )
       call check( nf90_put_att(ncid, varid8, 'units', 'm') )
       call check( nf90_put_att(ncid, varid8, 'long_name', 'Y coordinate of specific history location') )

       call check( nf90_def_var(ncid, "zhis", NF90_DOUBLE, dimids2, varid9) )
       call check( nf90_put_att(ncid, varid9, 'units', 'm') )
       call check( nf90_put_att(ncid, varid9, 'long_name', 'Z coordinate of specific history location') )

       call check( nf90_def_var(ncid, "ihis", NF90_INT, dimids2, varid10) )
       call check( nf90_put_att(ncid, varid10, 'units', '-') )
       call check( nf90_put_att(ncid, varid10, 'long_name', 'i coordinate of specific history location') )

       call check( nf90_def_var(ncid, "jhis", NF90_INT, dimids2, varid11) )
       call check( nf90_put_att(ncid, varid11, 'units', '-') )
       call check( nf90_put_att(ncid, varid11, 'long_name', 'j coordinate of specific history location') )

       call check( nf90_def_var(ncid, "khis", NF90_INT, dimids2, varid12) )
       call check( nf90_put_att(ncid, varid12, 'units', '-') )
       call check( nf90_put_att(ncid, varid12, 'long_name', 'k coordinate of specific history location') )

       call check( nf90_def_var(ncid, "time_his", NF90_DOUBLE, dimids4, varid13) )
       call check( nf90_put_att(ncid, varid13, 'units', 's') )
       call check( nf90_put_att(ncid, varid13, 'long_name', 'Time series of history') )

	! also add svn info in output files:
       CALL check( nf90_put_att(ncid,nf90_global, "svnversion", trim(svnversion)))
       CALL check( nf90_put_att(ncid,nf90_global, "svnurl", trim(svnurl)))

       ! End define mode. This tells netCDF we are done defining metadata.
       call check( nf90_enddef(ncid) )
     
       ! Write the pretend data to the file. Although netCDF supports
       ! reading and writing subsets of data, in this case we write all the
       ! data in one operation.
       call check( nf90_put_var(ncid, varid1, Uhis) )
       call check( nf90_put_var(ncid, varid2, Vhis) )
       call check( nf90_put_var(ncid, varid3, Whis) )
       call check( nf90_put_var(ncid, varid4, RHOhis) )
       call check( nf90_put_var(ncid, varid5, Phis) )
	if (nfrac>0) then
       call check( nf90_put_var(ncid, varid6, Chis) )
	endif
       call check( nf90_put_var(ncid, varid7, xhis) )
       call check( nf90_put_var(ncid, varid8, yhis) )
       call check( nf90_put_var(ncid, varid9, zhis) )
       call check( nf90_put_var(ncid, varid10, ihis) )
       call check( nf90_put_var(ncid, varid11, jhis) )
       call check( nf90_put_var(ncid, varid12, khis) )
       call check( nf90_put_var(ncid, varid13, this) )

       call check( nf90_close(ncid) )

	ENDIF

	END SUBROUTINE
    
	SUBROUTINE writeprogress(time_np,t_end,istep,dt,cput1,cput10a,cput10b,cput11a,cput11b,trestart)
	
	implicit none
	
	real time_np,t_end,dt,cput1,cput10a,cput10b,cput11a,cput11b,secrunt,hrrunt,trestart
	integer istep 
	integer,dimension(8) :: ttvalues
	CHARACTER(len=3) :: mons(12)
	INTEGER(2) :: monlen(12)
	INTEGER(2) :: tmpyear, tmpmonth, tmpday, tmphour, tmpminute, tmpsecond
	CHARACTER(len=2) :: dd, hh, mm
	CHARACTER(len=4) :: yyyy
	CHARACTER(len=17) :: datumtemp
	
		mons = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
		monlen = [31, 28,   31,   30,   31,    30,   31,  31,   30,   31,  30,    31]
		
		!secrunt= ((1.-(time_np/t_end))/(time_np/t_end))*(cput10b-cput1) ! remaining run time in seconds
		secrunt= ((1.-((time_np-trestart)/(t_end-trestart)))/((time_np-trestart)/(t_end-trestart)))*(cput10b-cput1) ! remaining run time in seconds
		hrrunt=secrunt/3600.
		
		tmpyear = floor(secrunt/(365*24*3600))
		IF (MOD(tmpyear, 4)==0) THEN
			monlen(2) = 29
		ENDIF			
		secrunt = secrunt - tmpyear * 365*24*3600
		tmpmonth = floor(secrunt/(30.25*24*3600))
		secrunt = secrunt - tmpmonth * 30.25*24*3600
		tmpday = floor(secrunt/(24*3600))
		secrunt = secrunt - tmpday*24*3600
		tmphour = floor(secrunt/3600)
		secrunt = secrunt - tmphour*3600
		tmpminute = floor(secrunt/60)
		secrunt = secrunt - tmpminute*60
		tmpsecond = secrunt
		
		CALL DATE_AND_TIME(VALUES=ttvalues)	
		
		tmpsecond = tmpsecond+ttvalues(7)
		tmpminute = tmpminute+ttvalues(6)
		tmphour = tmphour+ttvalues(5)
		tmpday = tmpday+ttvalues(3)
		tmpmonth = tmpmonth+ttvalues(2)
		tmpyear = tmpyear+ttvalues(1)
				
		IF (tmpsecond > 59) THEN
			tmpsecond = tmpsecond - 60
			tmpminute = tmpminute+1
		ENDIF
		
		IF (tmpminute > 59) THEN
			tmpminute = tmpminute - 60
			tmphour = tmphour+1
		ENDIF
		
		IF (tmphour > 23) THEN
			tmphour = tmphour - 24
			tmpday = tmpday+1
		ENDIF
		
		IF (tmpday > monlen(tmpmonth)) THEN
			tmpday = tmpday - monlen(tmpmonth)
			tmpmonth = tmpmonth+1
		ENDIF
		
		IF (tmpmonth > 12) THEN
			tmpmonth = tmpmonth - 12
			tmpyear = tmpyear+1
		ENDIF
		
		WRITE(yyyy,'(i4)') tmpyear
		WRITE(dd,'(i2)') tmpday
		WRITE(hh,'(i2.2)') tmphour
		WRITE(mm,'(i2.2)') tmpminute
		
		datumtemp = dd//'-'//mons(tmpmonth)//'-'//yyyy//' '//hh//':'//mm
	
		write(*,'(a,f10.4,a,f10.4,a,f8.2,a)') ' # Time: ',time_np,' s. of ',t_end,' s ',100.*time_np/t_end,' %   #'
		WRITE(*,'(a,i10.0,a,f9.6,a)') ' # Timestep: ',istep, ' dt : ',dt,' s            #' 			
		write(*,'(a,f6.3,a,f6.3,a,f5.2,a)'),' # CPU t=',NINT((cput10b-cput10a)*1000.)/1000.,'s, 1x pois=',
     &   NINT((cput11b-cput11a)*1000.)/1000.,'s = ',NINT(10000.*(cput11b-cput11a)/(cput10b-cput10a))/100.,
     &'%          #'	
		WRITE(*,'(a,f6.1,a,a,a)') ' # Remaining CPU: ',NINT(hrrunt*10.)/10., ' hr, ETA: ',datumtemp,' #' 			
		
	 END SUBROUTINE
			
			