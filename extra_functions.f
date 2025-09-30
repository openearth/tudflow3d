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
	  real Wmaxx,tmp11,tmp12,tmp13,tmp4,tmp11global,tmp12global,tmp13global
	  real dtnew,tmp10,tmp10global,Umaxx,Vmaxx
	  double precision localsum_Vol_V,globalsum_Vol_V

	  
!		IF (nobst>0.and.bp(1)%forever.eq.0.and.time_np+dt.gt.bp(1)%t0.and.counter<10) THEN
!			counter=counter+1
!			Courant=0.1
!			IF (rank.eq.0) THEN
!				write(*,*),'Placement bedplume, CFL=0.1 for 10 time steps ',counter,dt
!			ENDIF
!		ELSE 
			Courant = CFL
!		ENDIF
		IF (oPRHO.eq.1) THEN 
			dt_old = 1.e18 !1st order estimate P^n+1=P^n
		ELSE 
			dt_old=dt		!2nd order estimate P^n+1=P^n+dt/dt_old*(P^n - P^n-1)
		ENDIF 
		dtold=dt
      IF (istep.le.10) THEN
		  dt = dt_ini
		  IF (n_dtavg>0) THEN  
			dt_series=dt_max 
		  ENDIF 	
      ELSE
	dt = dt_max
      ENDIF
      IF (istep.le.1000) THEN
		IF (Uoutflow.eq.2) Uoutflow=999 !first 1000 time steps don't use convective outflow 
      ELSE
		IF (Uoutflow.eq.999) Uoutflow=2 !after 10 time steps use convective outflow 
      ENDIF	  
!		write(*,*),'dt,rank voor',dt,rank
      IF (CNdiffz.eq.1.or.CNdiffz.eq.2) THEN
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
!!!!!	  IF (CNdiffz.eq.11) THEN !AB2-CN, use stability limits viscosity from Wesseling for AB2-CN + standard CFL advection
!!!!!		  do k=1,kmax
!!!!!			 do j=1,jmax
!!!!!				do i=1,imax
!!!!!				df = Rp(i)*(phiv(j)-phiv(j-1))
!!!!!				df2 = df * df 
!!!!!				dr2 = dr(i) * dr(i) 
!!!!!				kcoeff = ekm(i,j,k) /rnew(i,j,k)
!!!!!				Wmaxx = MAX(abs(Wnew(i,j,k)-wsettlingmax),abs(Wnew(i,j,k)-wsettlingmin),
!!!!!     &				abs(Wnew(i,j,k-1)-wsettlingmax),abs(Wnew(i,j,k-1)-wsettlingmin))
!!!!!				Umaxx = MAX(abs(Unew(i,j,k)),abs(Unew(i-1,j,k))) 
!!!!!				Vmaxx = MAX(abs(Vnew(i,j,k)),abs(Vnew(i,j-1,k)))
!!!!!!!!				tmp1 = ( Umaxx / dr(i) ) +
!!!!!!!!     &             ( Vmaxx /         df        ) +
!!!!!!!!     &             ( Wmaxx /         dz        )    
!!!!!				tmp2 = ( 1.0/dr2 + 1.0/df2 + 1.0/dz2 )
!!!!!!				tmp3 = 1.0 / ( 0.001 * tmp2 * kcoeff + 3.*tmp1 + 1.e-12 ) !19-11-2021: included dt<1000*dx^2/(kcoeff) as additional time step restriction for CN-implicit as tests for very high viscosity and no flow gave unstable results; theoretically it is not needed for stability + factor 3 for CFL advection criterium as AB2 needs <0.3 and CFL user-input can be as high as 0.9	
!!!!!
!!!!!!!!				! switched off courant advection criterium 17-1-2022 as Wesseling doesn't have this 
!!!!!!!!				tmp3 = Courant / ( 3.*tmp1 + 1.e-12 )
!!!!!				tmp11 = Courant*1.333333*kcoeff/(Umaxx**2+Vmaxx**2+Wmaxx**2+1.e-12)
!!!!!				tmp12 = Courant*0.25/(kcoeff*(1./dr2+1./df2+1./dz2))
!!!!!				tmp13 = Courant*(3.*kcoeff)**0.333333*(Umaxx**4/dr2+Vmaxx**4/df2+Wmaxx**4/dz2+1.e-12)**-0.333333
!!!!!				
!!!!!!!!				dt = min( dt , tmp3 , max(tmp11,min(tmp12,tmp13)) )
!!!!!				dt = min( dt , max(tmp11,min(tmp12,tmp13)) )
!!!!!				enddo
!!!!!			 enddo
!!!!!		  enddo	  
!!!!!		dtmp=dt
!!!!!		call mpi_allreduce(dtmp,dt,1,mpi_double_precision,mpi_min,mpi_comm_world,ierr)	
!!!!!	  ELSEIF (CNdiffz.eq.12) THEN !AB2-CN, use stability limits viscosity from Wesseling for AB2-CN + standard CFL advection + diffusion limit from Chang 1991
!!!!!		  do k=1,kmax
!!!!!			 do j=1,jmax
!!!!!				do i=1,imax
!!!!!				df = Rp(i)*(phiv(j)-phiv(j-1))
!!!!!				df2 = df * df 
!!!!!				dr2 = dr(i) * dr(i) 
!!!!!				kcoeff = ekm(i,j,k) /rnew(i,j,k)
!!!!!				Wmaxx = MAX(abs(Wnew(i,j,k)-wsettlingmax),abs(Wnew(i,j,k)-wsettlingmin),
!!!!!     &				abs(Wnew(i,j,k-1)-wsettlingmax),abs(Wnew(i,j,k-1)-wsettlingmin))
!!!!!				Umaxx = MAX(abs(Unew(i,j,k)),abs(Unew(i-1,j,k))) 
!!!!!				Vmaxx = MAX(abs(Vnew(i,j,k)),abs(Vnew(i,j-1,k)))
!!!!!!!!				tmp1 = ( Umaxx / dr(i) ) +
!!!!!!!!     &             ( Vmaxx /         df        ) +
!!!!!!!!     &             ( Wmaxx /         dz        )    
!!!!!				tmp2 = ( 1.0/dr2 + 1.0/df2 + 1.0/dz2 )
!!!!!!!!				tmp3 = Courant / ( CNdiff_factor * 0.75* tmp2 * kcoeff + 3.*tmp1 + 1.e-12 ) !19-11-2021: included dt<0.75*dx^2/(kcoeff)/CNdiff_factor from Chang 1991 + factor 3 for CFL advection criterium as AB2 needs <0.3 and CFL user-input can be as high as 0.9	 
!!!!!				tmp3 = Courant / ( CNdiff_factor * 0.75* tmp2 * kcoeff + 1.e-12 ) !19-11-2021: included dt<0.75*dx^2/(kcoeff)/CNdiff_factor from Chang 1991; 17-1-2022 removed advection tmp1 criterium	 
!!!!!!				tmp3 = Courant / ( 3.*tmp1 + 1.e-12 )
!!!!!				tmp11 = Courant*1.333333*kcoeff/(Umaxx**2+Vmaxx**2+Wmaxx**2+1.e-12)
!!!!!				tmp12 = Courant*0.25/(kcoeff*(1./dr2+1./df2+1./dz2))
!!!!!				tmp13 = Courant*(3.*kcoeff)**0.333333*(Umaxx**4/dr2+Vmaxx**4/df2+Wmaxx**4/dz2+1.e-12)**-0.333333
!!!!!				
!!!!!				dt = min( dt , tmp3 , max(tmp11,min(tmp12,tmp13)) )
!!!!!				enddo
!!!!!			 enddo
!!!!!		  enddo	  
!!!!!		dtmp=dt
!!!!!		call mpi_allreduce(dtmp,dt,1,mpi_double_precision,mpi_min,mpi_comm_world,ierr)	
	IF (CNdiffz.eq.11) THEN  !CN-diff 3D implicit; theoretically unconditionally stable 
			! advection CFL dt criterium and diffusion dt criterium relaxed by factor CNdiff_dtfactor
			! when CNdiff_dtfactor is chosen very large then no dt limit for diffusion (theoretically unconditionally stable)
			! don't use Wesseling AB2-CN dt limits, because analyses of Wesseling AB2-CN dt limits reveal some problems: 
			! * at low viscosity, high velocity Wesseling AB2-CN dt gives much lower dt as advection CFL; why, that is not logical because even explicit diffusion scheme would give higher dt in this range of conditions and implicit diffusion is more stable 
			! * at high viscosity, all velocities Wesseling AB2-CN dt gives much higher dt as advection CFL; --> this seems dangerous
!		  do k=1,kmax
!			 do j=1,jmax
!				do i=1,imax
!				df = Rp(i)*(phiv(j)-phiv(j-1))
!				tmp1 = ( abs(Unew(i,j,k)) / ( Rp(i+1)-Rp(i) ) ) +
!     &             ( abs(Vnew(i,j,k)) /         df        ) +
!     &             ( MAX(abs(Wnew(i,j,k)-wsettlingmax),abs(Wnew(i,j,k)-wsettlingmin)) /         dz        )             
!				tmp3 = 1.0 / ( tmp1 + 1.e-12 )
!				tmp3 =Courant *tmp3 
!				dt = min( dt , tmp3 )
!				!dtmp = dt
!				enddo
!			 enddo
!		  enddo
!		dtmp=dt
!		call mpi_allreduce(dtmp,dt,1,mpi_double_precision,mpi_min,mpi_comm_world,ierr)		
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
				tmp3 = 1.0 / ( 2./CNdiff_dtfactor * tmp2 * kcoeff + tmp1 + 1.e-12 ) !included dt<CNdiff_dtfactor*dx^2/(2*kcoeff) --> Courant is extra safety factor for diffusion dt
				tmp3 =Courant *tmp3 
				dt = min( dt , tmp3 )
				!dtmp = dt

				enddo
			 enddo
		  enddo
		dtmp=dt
		call mpi_allreduce(dtmp,dt,1,mpi_double_precision,mpi_min,mpi_comm_world,ierr)	
	  ELSEIF (CNdiffz.eq.12) THEN !dt based on advection CFL and diffusion criterium relaxed by factor CNdiff_factor
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
				tmp3 = 1.0 / ( CNdiff_factor * 0.75 * tmp2 * kcoeff + tmp1 + 1.e-12 ) !included dt<0.75*dx^2/(kcoeff)/CNdiff_factor from Chang 1991 --> Courant is extra safety factor for diffusion dt
				tmp3 =Courant *tmp3 
				dt = min( dt , tmp3 )
				!dtmp = dt

				enddo
			 enddo
		  enddo
		dtmp=dt
		call mpi_allreduce(dtmp,dt,1,mpi_double_precision,mpi_min,mpi_comm_world,ierr)
	ELSEIF (CNdiffz.eq.31) THEN  !fully 3D CN-diff implicit; theoretically unconditionally stable 
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
				tmp3 = 1.0 / ( 2./CNdiff_dtfactor * tmp2 * kcoeff + tmp1 + 1.e-12 ) !included dt<CNdiff_dtfactor*dx^2/(2*kcoeff) --> Courant is extra safety factor for diffusion dt
				tmp3 =Courant *tmp3 
				dt = min( dt , tmp3 )
				enddo
			 enddo
		  enddo
		dtmp=dt
		call mpi_allreduce(dtmp,dt,1,mpi_double_precision,mpi_min,mpi_comm_world,ierr)	
	  ELSE 
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
				tmp3 = 1.0 / ( 2.0 * tmp2 * kcoeff + tmp1 + 1.e-12 ) !19-11-2021: adjusted with factor 2 on tmp2, because stability explicit diffusion operator is dt<dx^2/(2*kcoeff)
				tmp3 =Courant *tmp3 
				dt = min( dt , tmp3 )
				!dtmp = dt

				enddo
			 enddo
		  enddo
		dtmp=dt
		call mpi_allreduce(dtmp,dt,1,mpi_double_precision,mpi_min,mpi_comm_world,ierr)			  
	  ENDIF 
	
	  IF (n_dtavg>0) THEN  
		tel_dt=tel_dt+1 
		!IF (tel_dt>100) tel_dt=1
		!dt_series(tel_dt)=dt 
		!IF (dt.lt.dt_factor_avg * SUM(dt_series)/100.) THEN
		!	dt_series=dt/dt_factor_avg 
		!ENDIF 
		!dt = MIN(dt,dt_factor_avg * SUM(dt_series)/100.) !use minimum of avg_dt * factor or instantaneous dt
		IF (tel_dt>n_dtavg) tel_dt=1
		dt_series(tel_dt)=dt
		dtnew = MINVAL(dt_series)
		IF (dtnew/dtold>1.01) THEN !dtnew is growing
			dt=MIN(0.5*(dtnew+dtold),dtold*1.1,dt) !dampen growth, but should always obey timestep restrictions present time step
			dt_series(tel_dt)=dt
		ELSE 
			dt=dtnew
		ENDIF 
	  ENDIF 
		
	  dt=MIN(dt,dtold*1.1) !change in time step never more than +10% in one timestep for extra stability

	if (isnan(dt).or.dt<1.e-12) then
		drdt=rnew
		dcdt=cnew
		dudt=unew
		dvdt=vnew
		dwdt=wnew		
		rnew=rold 
		cnew=cold 
		unew=uold
		vnew=vold
		wnew=wold
		call output_nc('flow3D_',999999998,time_np) !output file with rcuvw from previous timestep (hopefully not crashed)
		rnew=drdt 
		cnew=dcdt 
		unew=dudt
		vnew=dvdt
		wnew=dwdt		
		call output_nc('flow3D_',999999999,time_np)
	endif 
	  
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
	if (periodicx.eq.1.and.ABS(dpdx).lt.1.e-12) THEN ! determination correct dpdx and dpdy based on wanted u or v
	   if (surf_layer>1.e-12) then !two different driving forces in the vertical are used 
	    localsum=0.
		localsum_Vol_V=0.
        do k=1,ksurf_bc !kmax
          do j=1,jmax
             do i=1,imax	  
		       localsum=localsum+Unew(i,j,k)*(Rp(i+1)-Rp(i))*(phiv(j)-phiv(j-1))*Ru(i)*dz*fc_global(i,j+jmax*rank,k)	
			   localsum_Vol_V=localsum_Vol_V+Vol_V(i,j+rank*jmax)*fc_global(i,j+jmax*rank,k)
			 enddo
		  enddo
	    enddo
		call mpi_allreduce(localsum,globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(localsum_Vol_V,globalsum_Vol_V,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		Uav = globalsum/(globalsum_Vol_V+1.e-12)
		if (istep>1) THEN
			du = Uav-Uavold
		else
			du = 0.
		endif
		Uavold = Uav
!		IF ((MAX(U_b,U_b3)**2+MAX(V_b,V_b3)**2)>1.e-18) THEN 
!			Tadapt = 100.*depth/sqrt(MAX(U_b,U_b3)**2+MAX(V_b,V_b3)**2+1.e-18) 
!		ELSE 
			Tadapt = 1000.*dt 
!		ENDIF
		duu = Tadapt*du/dt !30 seconds worked
		Uav = Uav + duu 
		dpdx1 = dpdx1 + ABS(U_b-Uav)*(U_b-Uav)/depth*dt/Tadapt
		!dpdx1 = dpdx1 + MIN(ABS(ABS(U_b)*(U_b-Uav)/depth*dt/Tadapt),0.1*ABS(dpdx1))*SIGN(1.,ABS(U_b)*(U_b-Uav)/depth*dt/Tadapt)
		do i=1,imax
		  do j=1,jmax
		    do k=1,ksurf_bc 
			  Ppropx(i,j,k) = dpdx1*rhU(i,j,k)*fc_global(i,j+jmax*rank,k) !rnew !with variable density important to use rnew and not rho_b 
		    enddo
		  enddo
		enddo
!		Ppropx(:,:,1:ksurf_bc) = dpdx1*rhU(:,:,1:ksurf_bc) !rnew !with variable density important to use rnew and not rho_b 
		Uav = Uav - duu
		if (rank.eq.0) then
			if (mod(istep,10) .eq.0) then   
				write(*,'(a,e11.4,a,f6.2,a,f6.2,a)') ' # bottom-layer dpdx: ',dpdx1,' Pa/m/(kg/m3); Uav: ',Uav,'; U_b:',U_b,' m/s #'
			endif  
		endif
	    localsum=0.
		localsum_Vol_V=0.
        do k=ksurf_bc+1,kmax !kmax
          do j=1,jmax
             do i=1,imax	  
		       localsum=localsum+Unew(i,j,k)*(Rp(i+1)-Rp(i))*(phiv(j)-phiv(j-1))*Ru(i)*dz*fc_global(i,j+jmax*rank,k)
			   localsum_Vol_V=localsum_Vol_V+Vol_V(i,j+rank*jmax)*fc_global(i,j+jmax*rank,k)
			 enddo
		  enddo
	    enddo
		call mpi_allreduce(localsum,globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(localsum_Vol_V,globalsum_Vol_V,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		Uav = globalsum/(globalsum_Vol_V+1.e-12)
		if (istep>1) THEN
			du = Uav-U3avold
		else
			du = 0.
		endif
		U3avold = Uav
!		IF ((MAX(U_b,U_b3)**2+MAX(V_b,V_b3)**2)>1.e-18) THEN 
!			Tadapt = 100.*depth/sqrt(MAX(U_b,U_b3)**2+MAX(V_b,V_b3)**2+1.e-18) 
!		ELSE 
			Tadapt = 1000.*dt 
!		ENDIF		
		duu = Tadapt*du/dt !30 seconds worked
		Uav = Uav + duu 
		dpdx3 = dpdx3 + ABS(U_b3-Uav)*(U_b3-Uav)/depth*dt/Tadapt
		!dpdx1 = dpdx1 + MIN(ABS(ABS(U_b)*(U_b-Uav)/depth*dt/Tadapt),0.1*ABS(dpdx1))*SIGN(1.,ABS(U_b)*(U_b-Uav)/depth*dt/Tadapt)
		do i=1,imax
		  do j=1,jmax
		    do k=ksurf_bc+1,kmax 
			  Ppropx(i,j,k) = dpdx3*rhU(i,j,k)*fc_global(i,j+jmax*rank,k) !rnew !with variable density important to use rnew and not rho_b 
		    enddo
		  enddo
		enddo		
!		Ppropx(:,:,ksurf_bc+1:kmax) = dpdx3*rhU(:,:,ksurf_bc+1:kmax) !rnew !with variable density important to use rnew and not rho_b 
		Uav = Uav - duu
		if (rank.eq.0) then
			if (mod(istep,10) .eq.0) then   
				write(*,'(a,e11.4,a,f6.2,a,f6.2,a)') ' # top-layer dpdx: ',dpdx3,' Pa/m/(kg/m3); Uav: ',Uav,'; U_b3:',U_b3,' m/s #'
			endif  
		endif		
	   else ! one driving force over whole vertical
	    localsum=0.
		localsum_Vol_V=0.
        do k=1,kmax
          do j=1,jmax
             do i=1,imax	  
		       localsum=localsum+Unew(i,j,k)*(Rp(i+1)-Rp(i))*(phiv(j)-phiv(j-1))*Ru(i)*dz*fc_global(i,j+jmax*rank,k)
			   localsum_Vol_V=localsum_Vol_V+Vol_V(i,j+rank*jmax)*fc_global(i,j+jmax*rank,k)
			 enddo
		  enddo
	    enddo
		call mpi_allreduce(localsum,globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(localsum_Vol_V,globalsum_Vol_V,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		Uav = globalsum/(globalsum_Vol_V+1.e-12)
		if (istep>1) THEN
			du = Uav-Uavold
		else
			du = 0.
		endif
		Uavold = Uav
!		IF ((U_b**2+V_b**2)>1.e-18) THEN 
!			Tadapt = 100.*depth/sqrt(U_b**2+V_b**2+1.e-18)
!		ELSE 
			Tadapt = 1000.*dt 
!		ENDIF		
		duu = Tadapt*du/dt !30 seconds worked
		Uav = Uav + duu 
		dpdx1 = dpdx1 + ABS(U_b-Uav)*(U_b-Uav)/depth*dt/Tadapt
		do i=1,imax
		  do j=1,jmax
		    do k=1,kmax 
			  Ppropx(i,j,k) = dpdx1*rhU(i,j,k)*fc_global(i,j+jmax*rank,k) !rnew !with variable density important to use rnew and not rho_b 
		    enddo
		  enddo
		enddo		
!		Ppropx = dpdx1*rhU !rnew !with variable density important to use rnew and not rho_b 
		
		Uav = Uav - duu
		if (rank.eq.0) then
			if (mod(istep,10) .eq.0) then   
				write(*,'(a,e11.4,a,f6.2,a,f6.2,a)') ' # dpdx: ',dpdx1,' Pa/m/(kg/m3); Uav: ',Uav,'; U_b:',U_b,' m/s #'
     		    !write(*,*),du,duu,ABS(U_b)*(U_b-(Uav+duu))/depth*dt/Tadapt
			endif  
		endif			
	   endif
	endif   
	if (periodicy.eq.1.and.ABS(dpdy).lt.1.e-12) THEN ! test determination correct dpdx and dpdy based on wanted u or v
	   if (surf_layer>1.e-12) then !two different driving forces in the vertical are used 
	    localsum=0.
		localsum_Vol_V=0.
        do k=1,ksurf_bc !kmax
          do j=1,jmax
             do i=1,imax	  
		       localsum=localsum+Vnew(i,j,k)*(Ru(i)-Ru(i-1))*(phip(j+1)-phip(j))*Rp(i)*dz*fc_global(i,j+jmax*rank,k)
			   localsum_Vol_V=localsum_Vol_V+Vol_V(i,j+rank*jmax)*fc_global(i,j+jmax*rank,k)
			 enddo
		  enddo
	    enddo
		call mpi_allreduce(localsum,globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)	
		call mpi_allreduce(localsum_Vol_V,globalsum_Vol_V,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		Vav = globalsum/(globalsum_Vol_V+1.e-12)
		if (istep>1) THEN
			du = Vav-Vavold
		else
			du = 0.
		endif
		Vavold = Vav
!		IF ((MAX(U_b,U_b3)**2+MAX(V_b,V_b3)**2)>1.e-18) THEN 
!			Tadapt = 100.*depth/sqrt(MAX(U_b,U_b3)**2+MAX(V_b,V_b3)**2+1.e-18) 
!		ELSE 
			Tadapt = 1000.*dt 
!		ENDIF
		duu = Tadapt*du/dt !30 seconds worked
		Vav = Vav + duu 
		dpdy1 = dpdy1 + ABS(V_b-Vav)*(V_b-Vav)/depth*dt/Tadapt 
		do i=1,imax
		  do j=1,jmax
		    do k=1,ksurf_bc 
			  Ppropy(i,j,k) = dpdy1*rhV(i,j,k)*fc_global(i,j+jmax*rank,k) !rnew !with variable density important to use rnew and not rho_b 
		    enddo
		  enddo
		enddo
		!Ppropy(:,:,1:ksurf_bc) = dpdy1*rhV(:,:,1:ksurf_bc) !rnew !with variable density important to use rnew and not rho_b
		Vav = Vav - duu
		if (rank.eq.0) then
			if (mod(istep,10) .eq.0) then   
				write(*,'(a,e11.4,a,f6.2,a,f6.2,a)') ' # bottom-layer dpdy: ',dpdy1,' Pa/m/(kg/m3); Vav: ',Vav,'; V_b:',V_b,' m/s #'
			endif  
		endif	
	    localsum=0.
		localsum_Vol_V=0.
        do k=ksurf_bc+1,kmax !kmax
          do j=1,jmax
             do i=1,imax	  
		       localsum=localsum+Vnew(i,j,k)*(Ru(i)-Ru(i-1))*(phip(j+1)-phip(j))*Rp(i)*dz*fc_global(i,j+jmax*rank,k)
			   localsum_Vol_V=localsum_Vol_V+Vol_V(i,j+rank*jmax)*fc_global(i,j+jmax*rank,k)
			 enddo
		  enddo
	    enddo
		call mpi_allreduce(localsum,globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)	
		call mpi_allreduce(localsum_Vol_V,globalsum_Vol_V,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		Vav = globalsum/(globalsum_Vol_V+1.e-12)
		if (istep>1) THEN
			du = Vav-V3avold
		else
			du = 0.
		endif
		V3avold = Vav
!		IF ((MAX(U_b,U_b3)**2+MAX(V_b,V_b3)**2)>1.e-18) THEN 
!			Tadapt = 100.*depth/sqrt(MAX(U_b,U_b3)**2+MAX(V_b,V_b3)**2+1.e-18) 
!		ELSE 
			Tadapt = 1000.*dt 
!		ENDIF
		duu = Tadapt*du/dt !30 seconds worked
		Vav = Vav + duu 
		dpdy3 = dpdy3 + ABS(V_b3-Vav)*(V_b3-Vav)/depth*dt/Tadapt 
		do i=1,imax
		  do j=1,jmax
		    do k=ksurf_bc+1,kmax
			  Ppropy(i,j,k) = dpdy3*rhV(i,j,k)*fc_global(i,j+jmax*rank,k) !rnew !with variable density important to use rnew and not rho_b 
		    enddo
		  enddo
		enddo
		!Ppropy(:,:,ksurf_bc+1:kmax) = dpdy3*rhV(:,:,ksurf_bc+1:kmax) !rnew !with variable density important to use rnew and not rho_b
		Vav = Vav - duu
		if (rank.eq.0) then
			if (mod(istep,10) .eq.0) then   
				write(*,'(a,e11.4,a,f6.2,a,f6.2,a)') ' # top-layer dpdy: ',dpdy3,' Pa/m/(kg/m3); Vav: ',Vav,'; V_b3:',V_b3,' m/s #'
			endif  
		endif		
	   else ! one driving force over whole vertical
	    localsum=0.
		localsum_Vol_V=0.
        do k=1,kmax
          do j=1,jmax
             do i=1,imax	  
		       localsum=localsum+Vnew(i,j,k)*(Ru(i)-Ru(i-1))*(phip(j+1)-phip(j))*Rp(i)*dz*fc_global(i,j+jmax*rank,k)
			   localsum_Vol_V=localsum_Vol_V+Vol_V(i,j+rank*jmax)*fc_global(i,j+jmax*rank,k)
			 enddo
		  enddo
	    enddo
		call mpi_allreduce(localsum,globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)	
		call mpi_allreduce(localsum_Vol_V,globalsum_Vol_V,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		Vav = globalsum/(globalsum_Vol_V+1.e-12)
		if (istep>1) THEN
			du = Vav-Vavold
		else
			du = 0.
		endif
		Vavold = Vav
!		IF ((U_b**2+V_b**2)>1.e-18) THEN 
!			Tadapt = 100.*depth/sqrt(U_b**2+V_b**2+1.e-18)
!		ELSE 
			Tadapt = 1000.*dt 
!		ENDIF
		duu = Tadapt*du/dt !30 seconds worked
		Vav = Vav + duu 
		dpdy1 = dpdy1 + ABS(V_b-Vav)*(V_b-Vav)/depth*dt/Tadapt 
		do i=1,imax
		  do j=1,jmax
		    do k=1,kmax 
			  Ppropy(i,j,k) = dpdy1*rhV(i,j,k)*fc_global(i,j+jmax*rank,k) !rnew !with variable density important to use rnew and not rho_b 
		    enddo
		  enddo
		enddo		
		!Ppropy = dpdy1*rhV !rnew !with variable density important to use rnew and not rho_b
		Vav = Vav - duu
		if (rank.eq.0) then
			if (mod(istep,10) .eq.0) then   
				write(*,'(a,e11.4,a,f6.2,a,f6.2,a)') ' # dpdy: ',dpdy1,' Pa/m/(kg/m3); Vav: ',Vav,'; V_b:',V_b,' m/s #'
			endif  
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
	  real wnew3(0:i1,0:j1,0:k1),sum_c_ws,ctot,ws(nfrac)
	  
      divbar1 = 0.0
      divmax1 = 0.0
      divbar2 = 0.0
      divmax2 = 0.0
      divbar3 = 0.0
      divmax3 = 0.0	 

	 wnew3=wnew	  
	  
	 IF (slipvel.eq.2) THEN
	  do j=1,jmax
	    do i=1,imax   
		  do k=kbed(i,j)+1,kmax !1,kmax ! all cells below kbed no correction needed as dwdt=dwdt 	  
		    sum_c_ws=0.
			do n=1,nfrac
				ws(n)=-frac(n)%ws !ws defined positive downwards
				sum_c_ws=sum_c_ws+ws(n)*(cW(n,i,j,k)*frac(n)%rho/rhW(i,j,k)-cW(n,i,j,k))
			enddo
			wnew3(i,j,k)=wnew3(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
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
			wnew3(i,j,k)=wnew3(i,j,k)-sum_c_ws  !go from mixture velocity (centre of mass velocity) to velocity of volume centre
          enddo
        enddo
      enddo
	 ENDIF
	 
	  
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
	 
	    div3= ( Ru(i)*Unew(i,j,k) - Ru(i-1)*Unew(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vnew(i,j,k) -         Vnew(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       Wnew3(i,j,k) -         Wnew3(i,j,k-1) ) / ( dz )  	 
	 
	DO n2=1,nbedplume !correct div2 for adding or subtracting volume:
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
			div3=div3-bp(n2)%Q*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells
		   endif
		ENDIF
	ENDDO ! bedplume loop
	
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

  
      call mpi_allreduce(divbar1,divbar_tot1,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(divmax1,divmax_tot1,1,mpi_real8,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(divbar2,divbar_tot2,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(divmax2,divmax_tot2,1,mpi_real8,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(divbar3,divbar_tot3,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(divmax3,divmax_tot3,1,mpi_real8,mpi_max,mpi_comm_world,ierr)	  
     
	if (jmax.gt.2) then 
	  if (rank.eq.0) write(6,100)divbar_tot1,divmax_tot1  
	 
100   format('Mass loss/gain dr/dt+div(ru)=0 (solved in continuity_solver=1) : Tot = ',e13.6,
     +                      '  Max = ',e13.6)
	  if (rank.eq.0) write(6,101)divbar_tot3,divmax_tot3  
101   format('Div(u) volume centre velocity (solved in continuity_solver=36) : Tot = ',e13.6,
     +                      '  Max = ',e13.6)
	  if (rank.eq.0) write(6,102)divbar_tot2,divmax_tot2 
102   format('Div(u) mass centre velocity   (solved in continuity_solver=35) : Tot = ',e13.6,
     +                      '  Max = ',e13.6)	 
	else 
	  if (rank.eq.0) then 
	    write(*,*),'No check on continuity because < 3 grid cells per CPU in lateral direction'
	  endif 
	endif 
    
	end
      

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
	REAL tau_flow_temp(1:imax,1:jmax),tau_frac_temp(1:nfrac,1:imax,1:jmax)
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
		tau_flow_temp=sqrt((0.5*(tau_fl_Unew(1:imax,1:jmax)+tau_fl_Unew(0:imax-1,1:jmax)))**2
     &   +(0.5*(tau_fl_Vnew(1:imax,1:jmax)+tau_fl_Vnew(1:imax,0:jmax-1)))**2)		
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
		  sig_tau_flow2 = sig_tau_flow2+tau_flow_temp*tau_flow_temp*dt
		  
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
		  tau_flow_avg = tau_flow_avg+tau_flow_temp*dt
		  
		  do n=1,nfrac
			  sigC2(n,:,:,:) = sigC2(n,:,:,:) + Cnew(n,1:imax,1:jmax,1:kmax)*Cnew(n,1:imax,1:jmax,1:kmax)*dt
			  sigUC(n,:,:,:) = sigUC(n,:,:,:) + Cnew(n,1:imax,1:jmax,1:kmax)*uudt
			  sigVC(n,:,:,:) = sigVC(n,:,:,:) + Cnew(n,1:imax,1:jmax,1:kmax)*vvdt
			  sigWC(n,:,:,:) = sigWC(n,:,:,:) + Cnew(n,1:imax,1:jmax,1:kmax)*wwdt
		  enddo
		  if (nfrac>0) then
			tau_flow_temp=ust_sl_new(1:imax,1:jmax)*ust_sl_new(1:imax,1:jmax)*rho_b
			sig_tau_sl2 = sig_tau_sl2+tau_flow_temp*tau_flow_temp*dt
			tau_sl_avg = tau_sl_avg+tau_flow_temp*dt
			tau_flow_temp=ust_bl_new(1:imax,1:jmax)*ust_bl_new(1:imax,1:jmax)*rho_b
			sig_tau_bl2 = sig_tau_bl2+tau_flow_temp*tau_flow_temp*dt
			tau_bl_avg = tau_bl_avg+tau_flow_temp*dt
			tau_flow_temp=ust_mud_new(1:imax,1:jmax)*ust_mud_new(1:imax,1:jmax)*rho_b
			sig_tau_mud2 = sig_tau_mud2+tau_flow_temp*tau_flow_temp*dt
			tau_mud_avg = tau_mud_avg+tau_flow_temp*dt
			tau_frac_temp=ust_frac_new(1:nfrac,1:imax,1:jmax)*ust_frac_new(1:nfrac,1:imax,1:jmax)*rho_b
			sig_tau_frac2 = sig_tau_frac2+tau_frac_temp*tau_frac_temp*dt
			tau_frac_avg = tau_frac_avg+tau_frac_temp*dt
		  endif 
		
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
	character(1024) :: gitversion
	character(1024) :: url
	character(1024) :: date_make
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
       CALL check( nf90_put_att(ncid,nf90_global, "gitversion", trim(gitversion)))
       CALL check( nf90_put_att(ncid,nf90_global, "url", trim(url)))
	   CALL check( nf90_put_att(ncid,nf90_global, "date make", trim(date_make)))

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
			
		SUBROUTINE CN3Dpcg_d(xnew,x0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,maxiter,tol,rank) 
		! Pre-conditioned ConjungateGradient implicit solver for 3D Laplacian Ax=RHS 
		! Matrix A is a 3D Laplacian matrix containing 7 diagonals:
		! Main diagonal values of matrix A in DD3D(0:i1,0:j1,0:k1)
		! Off-diagonal values left and right in x-dir, y-dir and z-dir: Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D with size (0:i1,0:j1,0:k1) 
		! Righthandside in RHS(0:i1,0:j1,0:k1) 
		! solution output in xnew(0:i1,0:j1,0:k1) 
		! starting condition in x0(0:i1,0:j1,0:k1)
		
		implicit none
		
        include 'mpif.h'
        integer ierr,rank
		integer  i,j,k,i1,j1,k1,im,ip,jm,jp,km,kp,tel,maxiter !,imax,jmax,kmax
		integer ib,ie,jb,je,kb,ke
		real Ax3D(0:i1,0:j1,0:k1),Ay3D(0:i1,0:j1,0:k1),Az3D(0:i1,0:j1,0:k1)
		real Cx3D(0:i1,0:j1,0:k1),Cy3D(0:i1,0:j1,0:k1),Cz3D(0:i1,0:j1,0:k1)
		real DD3D(0:i1,0:j1,0:k1),RHS3D(0:i1,0:j1,0:k1)
		real x(0:i1,0:j1,0:k1),xnew(0:i1,0:j1,0:k1),p(0:i1,0:j1,0:k1),r(0:i1,0:j1,0:k1),rnew(0:i1,0:j1,0:k1) !,pnew(0:i1,0:j1,0:k1)
		real ap(0:i1,0:j1,0:k1),x0(0:i1,0:j1,0:k1)
		real tol,alpha,beta,rnew2,r2,pap,rnew2_global,r2_global,pap_global,eps,b2,b2_global
		real znew(0:i1,0:j1,0:k1),rz2_global,rz2,rznew2,rznew2_global		

		!x = RHS3D !start with RHS as initial guess
		x = x0 !starting condition
		p=0. 
		b2=0.
		r2=0.
		rz2=0.		
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					! residual = RHS - A*x 
!					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k)*DD3D(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
     &						    Cx3D(i  ,j  ,k  ) * x(i+1,j  ,k  ) - Ax3D(i  ,j  ,k  ) * x(i-1,j  ,k  ) -
     &							Cy3D(i  ,j  ,k  ) * x(i  ,j+1,k  ) - Ay3D(i  ,j  ,k  ) * x(i  ,j-1,k  ) -
     &							Cz3D(i  ,j  ,k  ) * x(i  ,j  ,k+1) - Az3D(i  ,j  ,k  ) * x(i  ,j  ,k-1) 
					znew(i,j,k) = r(i,j,k)/DD3D(i,j,k) !preconditioning
					p(i,j,k) = znew(i,j,k)
					b2 = b2 + RHS3D(i,j,k)*RHS3D(i,j,k) 					
					r2 = r2 + r(i,j,k)*r(i,j,k) 
					rz2 = rz2 + r(i,j,k)*znew(i,j,k)
				enddo
			enddo 
		enddo 
		call mpi_allreduce(b2,b2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)		
		call mpi_allreduce(rz2,rz2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)	
		IF (sqrt(r2_global).le.tol.or.sqrt(b2_global).le.tol) THEN
			!no need to iterate
			xnew=x0
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
			return 
		ENDIF
		call bound_intern_and_periodic(p)
			
		
				
		do tel=1,maxiter 
			rnew2=0.
			rznew2=0.
			pap=0.
			rnew2_global=0.
			rznew2_global=0.
			pap_global=0. 
			eps=0. 
			alpha=0. 
			beta=0.

			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						!ap(i,j,k) = p(i,j,k) + !*DIAG/DIAG
						ap(i,j,k) = p(i,j,k)*DD3D(i,j,k) + !*DIAG/DIAG
     &								Cx3D(i  ,j  ,k  ) * p(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * p(i-1,j  ,k  ) +
     &								Cy3D(i  ,j  ,k  ) * p(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * p(i  ,j-1,k  ) +
     &								Cz3D(i  ,j  ,k  ) * p(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * p(i  ,j  ,k-1) 
!						ap(i,j,k) = ap(i,j,k)/DD3D(i,j,k) !preconditioning
						!r2  = r2  + r(i,j,k)* r(i,j,k) 
						pap = pap + p(i,j,k)*ap(i,j,k) 	 
					enddo
				enddo 
			enddo 	
			
!			do i=1,imax 
!				do j=1,jmax 
!					do k=1,kmax 
!						!r2  = r2  + r(i,j,k)* r(i,j,k) 
!						pap = pap + p(i,j,k)*ap(i,j,k) 
!					enddo
!				enddo 
!			enddo 			
			!call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			call mpi_allreduce(pap,pap_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			alpha = rz2_global/pap_global !with pre-conditioner !r2_global/pap_global
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						xnew(i,j,k)  = x(i,j,k)  + alpha* p(i,j,k)
						rnew(i,j,k)  = r(i,j,k)  - alpha*ap(i,j,k)
						rnew2 = rnew2 + rnew(i,j,k)*rnew(i,j,k)
						znew(i,j,k) = rnew(i,j,k)/DD3D(i,j,k) !pre-conditioning
						rznew2 = rznew2 + rnew(i,j,k)*znew(i,j,k) 
					enddo
				enddo 
			enddo 			
			call mpi_allreduce(rnew2,rnew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			call mpi_allreduce(rznew2,rznew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			eps = sqrt(rnew2_global)/sqrt(b2_global) 
!				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dpcg iteration working; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!				ENDIF			
			IF (eps<=tol) THEN 
				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dcg iteration finished; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',
!     & 					tel,tol,eps,b2_global,MINVAL(DD3D(ib:ie,jb:je,kb:ke)),MAXVAL(DD3D(ib:ie,jb:je,kb:ke))
					write(*,*),'CN3Dpcg_d finished; tel,tol,eps,b2_global:',tel,tol,eps,b2_global	 
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF			
				EXIT 
			ELSEIF (tel.eq.maxiter) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration did not reach desired tolerance, TUDflow3D will continue with end result'
					write(*,*)'tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',
     &					tel,tol,eps,b2_global,MINVAL(DD3D(ib:ie,jb:je,kb:ke)),MAXVAL(DD3D(ib:ie,jb:je,kb:ke))
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF
				EXIT 
			ENDIF 
			beta = rznew2_global/rz2_global 
			IF (isnan(alpha).or.isnan(beta)) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration crashed; alpha,beta,tel,tol,eps=',alpha,beta,tel,tol,eps
				ENDIF
				EXIT 
			ENDIF

			! not finished, so continue:
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax
						p(i,j,k) = znew(i,j,k)  + beta*p(i,j,k)
						x(i,j,k) = xnew(i,j,k)
						r(i,j,k) = rnew(i,j,k)
						!p(i,j,k) = pnew(i,j,k)
					enddo
				enddo 
			enddo
			r2_global = rnew2_global 
			rz2_global = rznew2_global 
			call bound_intern_and_periodic(p)

		enddo ! interation loop 
		
		END SUBROUTINE
		
		SUBROUTINE CN3Dpcg_d2(xnew,x0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,maxiter,tol,rank) 
		! Pre-conditioned ConjungateGradient implicit solver for 3D Laplacian Ax=RHS 
		! Matrix A is a 3D Laplacian matrix containing 7 diagonals:
		! Main diagonal values of matrix A in DD3D(0:i1,0:j1,0:k1)
		! Off-diagonal values left and right in x-dir, y-dir and z-dir: Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D with size (0:i1,0:j1,0:k1) 
		! Righthandside in RHS(0:i1,0:j1,0:k1) 
		! solution output in xnew(0:i1,0:j1,0:k1) 
		! starting condition in x0(0:i1,0:j1,0:k1)
		
		implicit none
		
        include 'mpif.h'
        integer ierr,rank
		integer  i,j,k,i1,j1,k1,im,ip,jm,jp,km,kp,tel,maxiter !,imax,jmax,kmax
		integer ib,ie,jb,je,kb,ke
		real Ax3D(0:i1,0:j1,0:k1),Ay3D(0:i1,0:j1,0:k1),Az3D(0:i1,0:j1,0:k1)
		real Cx3D(0:i1,0:j1,0:k1),Cy3D(0:i1,0:j1,0:k1),Cz3D(0:i1,0:j1,0:k1)
		real DD3D(0:i1,0:j1,0:k1),RHS3D(0:i1,0:j1,0:k1)
		real x(0:i1,0:j1,0:k1),xnew(0:i1,0:j1,0:k1),p(0:i1,0:j1,0:k1),r(0:i1,0:j1,0:k1),rnew(0:i1,0:j1,0:k1) !,pnew(0:i1,0:j1,0:k1)
		real ap(0:i1,0:j1,0:k1),x0(0:i1,0:j1,0:k1)
		real tol,alpha,beta,rnew2,r2,pap,rnew2_global,r2_global,pap_global,eps,b2,b2_global
		
		! Diagonal pre-conditioner
		! Appy regular CG to solve AAxx=bb with AA=P^-1AP^-T and x=P^-T*xx and bb=P^-1*b
		! With diagonal pre-conditioner p_ii=sqrt(a_ii)

!!		!x = RHS3D/sqrt(DD3D) !start condition equal to RHS 
!!		x = x0 ! !starting condition 	
!!		!! left preconditioning M^-1Ax = M^-1b with M=D -->D^-1Ax = b/D
!!		! works for sdc 
!!		Ax3D = Ax3D/DD3D
!!		Cx3D = Cx3D/DD3D
!!		Ay3D = Ay3D/DD3D
!!		Cy3D = Cy3D/DD3D
!!		Az3D = Az3D/DD3D
!!		Cz3D = Cz3D/DD3D
!!		RHS3D = RHS3D/DD3D
		
		
		!! central preconditioning leading to symmetric matrix:
		!A~x~=b~ with A~ = P^-1AP^-T x=P^-Tx~ and b~ = P^-1b 
		!P = sqrt(D)
		!at end x=P^-Tx~
		call bound_3D(DD3D) !--> DD3D (ip im jp jm kp km) are needed  
		x = x0*sqrt(DD3D) !starting condition for x~
		Ax3D = Ax3D/sqrt(DD3D)
		Cx3D = Cx3D/sqrt(DD3D)
		Ay3D = Ay3D/sqrt(DD3D)
		Cy3D = Cy3D/sqrt(DD3D)
		Az3D = Az3D/sqrt(DD3D)
		Cz3D = Cz3D/sqrt(DD3D)
		RHS3D = RHS3D/sqrt(DD3D)
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					im=MAX(ib,i-1)
					jm=MAX(jb,j-1)
					km=MAX(kb,k-1)
					ip=MIN(ie,i+1) 
					jp=MIN(je,j+1)
					kp=MIN(ke,k+1) 			
					!if (im>=0) then 
						Ax3D(i,j,k)=Ax3D(i,j,k)/sqrt(DD3D(im,j,k))
					!endif 
					!if (ip<=i1) then 
						Cx3D(i,j,k)=Cx3D(i,j,k)/sqrt(DD3D(ip,j,k))
					!endif 					
					!if (jm>=0) then  
						Ay3D(i,j,k)=Ay3D(i,j,k)/sqrt(DD3D(i,jm,k))
					!endif
					!if (jp<=j1) then 
						Cy3D(i,j,k)=Cy3D(i,j,k)/sqrt(DD3D(i,jp,k))
					!endif 						
					!if (km>=0) then 
						Az3D(i,j,k)=Az3D(i,j,k)/sqrt(DD3D(i,j,km))
					!endif 
					!if (kp<=k1) then 
						Cz3D(i,j,k)=Cz3D(i,j,k)/sqrt(DD3D(i,j,kp))
					!endif 				
				enddo 
			enddo 
		enddo
		p=0.
		b2=0.
		r2=0. 
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					! residual = RHS - A*x 
					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
					!r(i,j,k) = 	RHS3D(i,j,k)/sqrt(DD3D(i,j,k))-x(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
!					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k)*DD3D(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
     &						    Cx3D(i  ,j  ,k  ) * x(i+1,j  ,k  ) - Ax3D(i  ,j  ,k  ) * x(i-1,j  ,k  ) -
     &							Cy3D(i  ,j  ,k  ) * x(i  ,j+1,k  ) - Ay3D(i  ,j  ,k  ) * x(i  ,j-1,k  ) -
     &							Cz3D(i  ,j  ,k  ) * x(i  ,j  ,k+1) - Az3D(i  ,j  ,k  ) * x(i  ,j  ,k-1) 
					p(i,j,k) = r(i,j,k)
					b2 = b2 + RHS3D(i,j,k)*RHS3D(i,j,k) 					
					r2 = r2 + r(i,j,k)*r(i,j,k) 
				enddo
			enddo 
		enddo 
		call mpi_allreduce(b2,b2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)		
		IF (sqrt(r2_global).le.tol.or.sqrt(b2_global).le.tol) THEN
			!no need to iterate
			xnew=x0
			return 
		ENDIF
		call bound_intern_and_periodic(p)
			
		
				
		do tel=1,maxiter 
			rnew2=0.
			pap=0.
			rnew2_global=0.
			pap_global=0. 
			eps=0. 
			alpha=0. 
			beta=0.

			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						ap(i,j,k) = p(i,j,k) + !*DIAG/DIAG
						!ap(i,j,k) = p(i,j,k)*DD3D(i,j,k) + !*DIAG/DIAG
     &								Cx3D(i  ,j  ,k  ) * p(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * p(i-1,j  ,k  ) +
     &								Cy3D(i  ,j  ,k  ) * p(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * p(i  ,j-1,k  ) +
     &								Cz3D(i  ,j  ,k  ) * p(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * p(i  ,j  ,k-1) 
						!r2  = r2  + r(i,j,k)* r(i,j,k) 
						pap = pap + p(i,j,k)*ap(i,j,k) 	 
					enddo
				enddo 
			enddo 	
		
			!call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			call mpi_allreduce(pap,pap_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			alpha = r2_global/pap_global
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						xnew(i,j,k)  = x(i,j,k)  + alpha* p(i,j,k)
						rnew(i,j,k)  = r(i,j,k)  - alpha*ap(i,j,k)
						rnew2 = rnew2 + rnew(i,j,k)*rnew(i,j,k)
					enddo
				enddo 
			enddo 			
			call mpi_allreduce(rnew2,rnew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			eps = sqrt(rnew2_global)/sqrt(b2_global) 
!				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dpcg iteration working; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!				ENDIF			
			IF (eps<=tol) THEN 
				xnew = xnew/sqrt(DD3D)
				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dcg iteration finished; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',
!     & 					tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
					write(*,*),'CN3Dpcg_d2 finished; tel,tol,eps,b2_global:',tel,tol,eps,b2_global
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF			
				EXIT 
			ELSEIF (tel.eq.maxiter) THEN 
				xnew = xnew/sqrt(DD3D)
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dcg iteration did not reach desired tolerance, TUDflow3D will continue with end result'
					write(*,*)'tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF
				EXIT 
			ENDIF 
			beta = rnew2_global/r2_global 
			IF (isnan(alpha).or.isnan(beta)) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration crashed; alpha,beta,tel,tol,eps=',alpha,beta,tel,tol,eps
				ENDIF
				EXIT 
			ENDIF

			! not finished, so continue:
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax
						p(i,j,k) = rnew(i,j,k)  + beta*p(i,j,k)
						x(i,j,k) = xnew(i,j,k)
						r(i,j,k) = rnew(i,j,k)
						!p(i,j,k) = pnew(i,j,k)
					enddo
				enddo 
			enddo
			r2_global = rnew2_global 
			call bound_intern_and_periodic(p)
		enddo ! interation loop 	
		
		END SUBROUTINE		
		
		SUBROUTINE CN3Dcg(xnew,x0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,maxiter,tol,rank) 
		! ConjungateGradient implicit solver for 3D Laplacian Ax=RHS without pre-conditioner
		! Matrix A is a 3D Laplacian matrix containing 7 diagonals:
		! Main diagonal values of matrix A in DD3D(0:i1,0:j1,0:k1)
		! Off-diagonal values left and right in x-dir, y-dir and z-dir: Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D with size (0:i1,0:j1,0:k1) 
		! Righthandside in RHS(0:i1,0:j1,0:k1) 
		! solution output in xnew(0:i1,0:j1,0:k1) 
		! starting condition in x0(0:i1,0:j1,0:k1) 
		
		implicit none
		
        include 'mpif.h'
        integer ierr,rank
		integer  i,j,k,i1,j1,k1,im,ip,jm,jp,km,kp,tel,maxiter !,imax,jmax,kmax
		integer ib,ie,jb,je,kb,ke
		real Ax3D(0:i1,0:j1,0:k1),Ay3D(0:i1,0:j1,0:k1),Az3D(0:i1,0:j1,0:k1)
		real Cx3D(0:i1,0:j1,0:k1),Cy3D(0:i1,0:j1,0:k1),Cz3D(0:i1,0:j1,0:k1)
		real DD3D(0:i1,0:j1,0:k1),RHS3D(0:i1,0:j1,0:k1)
		real x(0:i1,0:j1,0:k1),xnew(0:i1,0:j1,0:k1),p(0:i1,0:j1,0:k1),r(0:i1,0:j1,0:k1),rnew(0:i1,0:j1,0:k1) !,pnew(0:i1,0:j1,0:k1)
		real ap(0:i1,0:j1,0:k1),x0(0:i1,0:j1,0:k1)
		real tol,alpha,beta,rnew2,r2,pap,rnew2_global,r2_global,pap_global,eps,b2,b2_global

		!x = RHS3D !start condition equal to RHS 
		x = x0 !starting condition 
		p=0. 
		b2=0.
		r2=0. 
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					! residual = RHS - A*x 
!					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k)*DD3D(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
     &						    Cx3D(i  ,j  ,k  ) * x(i+1,j  ,k  ) - Ax3D(i  ,j  ,k  ) * x(i-1,j  ,k  ) -
     &							Cy3D(i  ,j  ,k  ) * x(i  ,j+1,k  ) - Ay3D(i  ,j  ,k  ) * x(i  ,j-1,k  ) -
     &							Cz3D(i  ,j  ,k  ) * x(i  ,j  ,k+1) - Az3D(i  ,j  ,k  ) * x(i  ,j  ,k-1) 
!					r(i,j,k) = r(i,j,k)/DD3D(i,j,k) !preconditioning
					p(i,j,k) = r(i,j,k)
					b2 = b2 + RHS3D(i,j,k)*RHS3D(i,j,k) 					
					r2 = r2 + r(i,j,k)*r(i,j,k) 
				enddo
			enddo 
		enddo 
		call mpi_allreduce(b2,b2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)		
		IF (sqrt(r2_global).le.tol.or.sqrt(b2_global).le.tol) THEN
			!no need to iterate
			xnew=x0
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
			return 
		ENDIF
		call bound_intern_and_periodic(p)
			
		
				
		do tel=1,maxiter 
			rnew2=0.
			pap=0.
			rnew2_global=0.
			pap_global=0. 
			eps=0. 
			alpha=0. 
			beta=0.

			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						!ap(i,j,k) = p(i,j,k) + !*DIAG/DIAG
						ap(i,j,k) = p(i,j,k)*DD3D(i,j,k) + !*DIAG/DIAG
     &								Cx3D(i  ,j  ,k  ) * p(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * p(i-1,j  ,k  ) +
     &								Cy3D(i  ,j  ,k  ) * p(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * p(i  ,j-1,k  ) +
     &								Cz3D(i  ,j  ,k  ) * p(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * p(i  ,j  ,k-1) 
!						ap(i,j,k) = ap(i,j,k)/DD3D(i,j,k) !preconditioning
						!r2  = r2  + r(i,j,k)* r(i,j,k) 
						pap = pap + p(i,j,k)*ap(i,j,k) 	 
					enddo
				enddo 
			enddo 	
			
!			do i=1,imax 
!				do j=1,jmax 
!					do k=1,kmax 
!						!r2  = r2  + r(i,j,k)* r(i,j,k) 
!						pap = pap + p(i,j,k)*ap(i,j,k) 
!					enddo
!				enddo 
!			enddo 			
			!call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			call mpi_allreduce(pap,pap_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			alpha = r2_global/pap_global
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						xnew(i,j,k)  = x(i,j,k)  + alpha* p(i,j,k)
						rnew(i,j,k)  = r(i,j,k)  - alpha*ap(i,j,k)
						rnew2 = rnew2 + rnew(i,j,k)*rnew(i,j,k)
					enddo
				enddo 
			enddo 			
			call mpi_allreduce(rnew2,rnew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			eps = sqrt(rnew2_global)/sqrt(b2_global) 
!				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dpcg iteration working; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!				ENDIF			
			IF (eps<=tol) THEN 
				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dcg iteration finished; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',
!     & 					tel,tol,eps,b2_global,MINVAL(DD3D(ib:ie,jb:je,kb:ke)),MAXVAL(DD3D(ib:ie,jb:je,kb:ke))
					write(*,*),'CN3Dcg finished; tel,tol,eps,b2_global:',tel,tol,eps,b2_global	 
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF			
				EXIT 
			ELSEIF (tel.eq.maxiter) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dcg iteration did not reach desired tolerance, TUDflow3D will continue with end result'
					write(*,*)'tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',
     &					tel,tol,eps,b2_global,MINVAL(DD3D(ib:ie,jb:je,kb:ke)),MAXVAL(DD3D(ib:ie,jb:je,kb:ke))
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF
				EXIT 
			ENDIF 
			beta = rnew2_global/r2_global 
			IF (isnan(alpha).or.isnan(beta)) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration crashed; alpha,beta,tel,tol,eps=',alpha,beta,tel,tol,eps
				ENDIF
				EXIT 
			ENDIF

			! not finished, so continue:
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax
						p(i,j,k) = rnew(i,j,k)  + beta*p(i,j,k)
						x(i,j,k) = xnew(i,j,k)
						r(i,j,k) = rnew(i,j,k)
						!p(i,j,k) = pnew(i,j,k)
					enddo
				enddo 
			enddo
			r2_global = rnew2_global 
			call bound_intern_and_periodic(p)

		enddo ! interation loop 		
		END SUBROUTINE
		
		SUBROUTINE CN3Dcg2(xnew,x0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,ib,ie,jb,je,kb,ke,maxiter,tol) 
		! ConjungateGradient implicit solver for 3D Laplacian Ax=RHS without pre-conditioner
		! Matrix A is a 3D Laplacian matrix containing 7 diagonals:
		! Main diagonal values of matrix A in DD3D(0:i1,0:j1,0:k1)
		! Off-diagonal values left and right in x-dir, y-dir and z-dir: Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D with size (0:i1,0:j1,0:k1) 
		! Righthandside in RHS(0:i1,0:j1,0:k1) 
		! solution output in xnew(0:i1,0:j1,0:k1) 
		! starting condition in x0(0:i1,0:j1,0:k1) 
		
		USE nlist 
		
		implicit none
				
        include 'mpif.h'
        integer ierr
		integer im,ip,jm,jp,km,kp,tel,maxiter !,imax,jmax,kmax
		integer ib,ie,jb,je,kb,ke
		real Ax3D(0:i1,0:j1,0:k1),Ay3D(0:i1,0:j1,0:k1),Az3D(0:i1,0:j1,0:k1)
		real Cx3D(0:i1,0:j1,0:k1),Cy3D(0:i1,0:j1,0:k1),Cz3D(0:i1,0:j1,0:k1)
		real DD3D(0:i1,0:j1,0:k1),RHS3D(0:i1,0:j1,0:k1)
		real x(0:i1,0:j1,0:k1),xnew(0:i1,0:j1,0:k1),pp(0:i1,0:j1,0:k1),rr(0:i1,0:j1,0:k1),rrnew(0:i1,0:j1,0:k1) !,pnew(0:i1,0:j1,0:k1)
		real ap(0:i1,0:j1,0:k1),x0(0:i1,0:j1,0:k1)
		real tol,alpha,beta,rnew2,r2,pap,rnew2_global,r2_global,pap_global,eps,b2,b2_global

		!x = RHS3D !start condition equal to RHS 
		x = x0 !starting condition 
		
		pp=0.
		b2=0.
		r2=0. 
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					! residual = RHS - A*x 
!					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
					rr(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k)*DD3D(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
     &						    Cx3D(i  ,j  ,k  ) * x(i+1,j  ,k  ) - Ax3D(i  ,j  ,k  ) * x(i-1,j  ,k  ) -
     &							Cy3D(i  ,j  ,k  ) * x(i  ,j+1,k  ) - Ay3D(i  ,j  ,k  ) * x(i  ,j-1,k  ) -
     &							Cz3D(i  ,j  ,k  ) * x(i  ,j  ,k+1) - Az3D(i  ,j  ,k  ) * x(i  ,j  ,k-1) 
!					r(i,j,k) = r(i,j,k)/DD3D(i,j,k) !preconditioning
					pp(i,j,k) = rr(i,j,k)
					b2 = b2 + RHS3D(i,j,k)*RHS3D(i,j,k) 					
					r2 = r2 + rr(i,j,k)*rr(i,j,k) 
				enddo
			enddo 
		enddo 
		call mpi_allreduce(b2,b2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)		
		IF (sqrt(r2_global).le.tol.or.sqrt(b2_global).le.tol) THEN
			!no need to iterate
			xnew=x0
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
			return 
		ENDIF
		call bound_intern_and_periodic(pp)
			
		
				
		do tel=1,maxiter 
			rnew2=0.
			pap=0.
			rnew2_global=0.
			pap_global=0. 
			eps=0. 
			alpha=0. 
			beta=0.

			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						!ap(i,j,k) = p(i,j,k) + !*DIAG/DIAG
						ap(i,j,k) = pp(i,j,k)*DD3D(i,j,k) + !*DIAG/DIAG
     &								Cx3D(i  ,j  ,k  ) * pp(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * pp(i-1,j  ,k  ) +
     &								Cy3D(i  ,j  ,k  ) * pp(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * pp(i  ,j-1,k  ) +
     &								Cz3D(i  ,j  ,k  ) * pp(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * pp(i  ,j  ,k-1) 
!						ap(i,j,k) = ap(i,j,k)/DD3D(i,j,k) !preconditioning
						!r2  = r2  + r(i,j,k)* r(i,j,k) 
						pap = pap + pp(i,j,k)*ap(i,j,k) 	 
					enddo
				enddo 
			enddo 	
			
!			do i=1,imax 
!				do j=1,jmax 
!					do k=1,kmax 
!						!r2  = r2  + r(i,j,k)* r(i,j,k) 
!						pap = pap + p(i,j,k)*ap(i,j,k) 
!					enddo
!				enddo 
!			enddo 			
			!call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			call mpi_allreduce(pap,pap_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			alpha = r2_global/pap_global
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						xnew(i,j,k)  = x(i,j,k)  + alpha* pp(i,j,k)
						rrnew(i,j,k)  = rr(i,j,k)  - alpha*ap(i,j,k)
						rnew2 = rnew2 + rrnew(i,j,k)*rrnew(i,j,k)
					enddo
				enddo 
			enddo 			
			call mpi_allreduce(rnew2,rnew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			eps = sqrt(rnew2_global)/sqrt(b2_global) 
!				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dpcg iteration working; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!				ENDIF			
			IF (eps<=tol) THEN 
				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dcg iteration finished; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',
!     & 					tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
					write(*,*),'CN3Dcg2 finished; tel,tol,eps,b2_global:',tel,tol,eps,b2_global	 
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF			
				EXIT 
			ELSEIF (tel.eq.maxiter) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dcg iteration did not reach desired tolerance, TUDflow3D will continue with end result'
					write(*,*)'tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF
				EXIT 
			ENDIF 
			beta = rnew2_global/r2_global 
			IF (isnan(alpha).or.isnan(beta)) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration crashed; alpha,beta,tel,tol,eps=',alpha,beta,tel,tol,eps
				ENDIF
				EXIT 
			ENDIF

			! not finished, so continue:
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax
						pp(i,j,k) = rrnew(i,j,k)  + beta*pp(i,j,k)
						x(i,j,k) = xnew(i,j,k)
						rr(i,j,k) = rrnew(i,j,k)
						!p(i,j,k) = pnew(i,j,k)
					enddo
				enddo 
			enddo
			r2_global = rnew2_global 
			call bound_intern_and_periodic(pp)

		enddo ! interation loop 		
		END SUBROUTINE		
		
		SUBROUTINE CN3Dpcg_ic(xnew,x0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,maxiter,tol,rank) 
		! Pre-conditioned ConjungateGradient implicit solver for 3D Laplacian Ax=RHS 
		! Matrix A is a 3D Laplacian matrix containing 7 diagonals:
		! Main diagonal values of matrix A in DD3D(0:i1,0:j1,0:k1)
		! Off-diagonal values left and right in x-dir, y-dir and z-dir: Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D with size (0:i1,0:j1,0:k1) 
		! Righthandside in RHS(0:i1,0:j1,0:k1) 
		! solution output in xnew(0:i1,0:j1,0:k1) 
		! starting condition in x0(0:i1,0:j1,0:k1)
		
		implicit none
		
        include 'mpif.h'
        integer ierr,rank
		integer  i,j,k,i1,j1,k1,im,ip,jm,jp,km,kp,tel,maxiter !,imax,jmax,kmax
		integer ib,ie,jb,je,kb,ke,tel_end,imax2,jmax2,kmax2 
		real Ax3D(0:i1,0:j1,0:k1),Ay3D(0:i1,0:j1,0:k1),Az3D(0:i1,0:j1,0:k1)
		real Cx3D(0:i1,0:j1,0:k1),Cy3D(0:i1,0:j1,0:k1),Cz3D(0:i1,0:j1,0:k1)
		real DD3D(0:i1,0:j1,0:k1),RHS3D(0:i1,0:j1,0:k1)
		real x(0:i1,0:j1,0:k1),xnew(0:i1,0:j1,0:k1),p(0:i1,0:j1,0:k1),r(0:i1,0:j1,0:k1),rnew(0:i1,0:j1,0:k1) !,pnew(0:i1,0:j1,0:k1)
		real ap(0:i1,0:j1,0:k1),x0(0:i1,0:j1,0:k1)
		real tol,alpha,beta,rnew2,r2,pap,rnew2_global,r2_global,pap_global,eps,b2,b2_global
		real znew(0:i1,0:j1,0:k1),rz2_global,rz2,rznew2,rznew2_global
		integer*2 irks((i1+1)*(j1+1)*(k1+1)),jrks((i1+1)*(j1+1)*(k1+1)),krks((i1+1)*(j1+1)*(k1+1))
		real Ax3D_c(0:i1,0:j1,0:k1),Ay3D_c(0:i1,0:j1,0:k1),Az3D_c(0:i1,0:j1,0:k1)
		!real Cx3D_c(0:i1,0:j1,0:k1),Cy3D_c(0:i1,0:j1,0:k1),Cz3D_c(0:i1,0:j1,0:k1)
		real DD3D_c(0:i1,0:j1,0:k1)
		

		
		! make local (on this partition) incomplete Cholesky decomposition:
		call bound_3D(DD3D) 
		Ax3D_c=0. 
		Ay3D_c=0.
		Az3D_c=0.
		!Cx3D_c=0.
		!Cy3D_c=0.
		!Cz3D_c=0.
		DD3D_c=0.
		tel=0 
		do k=kb,ke 
			do j=jb,je
				do i=ib,ie
					tel=tel+1;
					!tel2=i+(j-1)*nx+(k-1)*ny*nx;
					irks(tel)=i;
					jrks(tel)=j;
					krks(tel)=k;		
				enddo
			enddo 
		enddo
		tel_end=tel 
		
		imax2=(ie-ib+1)
		jmax2=(je-jb+1)
		kmax2=(ke-kb+1) 
		tel=1;
		DD3D_c(irks(tel),jrks(tel),krks(tel))=DD3D(irks(tel),jrks(tel),krks(tel));
		do tel = 2,imax2
			DD3D_c(irks(tel),jrks(tel),krks(tel))=DD3D(irks(tel),jrks(tel),krks(tel))
     &			-Ax3D(irks(tel),jrks(tel),krks(tel))**2/DD3D_c(irks(tel-1),jrks(tel-1),krks(tel-1))
		enddo
		do tel = imax2+1,imax2*jmax2
			DD3D_c(irks(tel),jrks(tel),krks(tel))=DD3D(irks(tel),jrks(tel),krks(tel)) 
     &			-Ax3D(irks(tel),jrks(tel),krks(tel))**2/DD3D_c(irks(tel-1),jrks(tel-1),krks(tel-1)) 
     &			-Ay3D(irks(tel),jrks(tel),krks(tel))**2/DD3D_c(irks(tel-imax2),jrks(tel-imax2),krks(tel-imax2))
		enddo
		do tel = imax2*jmax2+1,imax2*jmax2*kmax2
			DD3D_c(irks(tel),jrks(tel),krks(tel))=DD3D(irks(tel),jrks(tel),krks(tel)) 
     &			-Ax3D(irks(tel),jrks(tel),krks(tel))**2/DD3D_c(irks(tel-1),jrks(tel-1),krks(tel-1)) 
     &			-Ay3D(irks(tel),jrks(tel),krks(tel))**2/DD3D_c(irks(tel-imax2),jrks(tel-imax2),krks(tel-imax2)) 
     &			-Az3D(irks(tel),jrks(tel),krks(tel))**2/DD3D_c(irks(tel-imax2*jmax2),jrks(tel-imax2*jmax2),krks(tel-imax2*jmax2))
		enddo 
		
		call bound_3D(DD3D_c) 
		do k=kb,ke 
			do j=jb,je
				do i=ib,ie	
					im=i-1 
					jm=j-1
					km=k-1
					ip=i+1 
					jp=j+1
					kp=k+1 
					if (im>=0) then 
						Ax3D_c(i,j,k)=Ax3D(i,j,k)/sqrt(DD3D_c(im,j,k))
					endif 
!					if (ip<=i1) then 
!						Cx3D_c(i,j,k)=Cx3D(i,j,k)/sqrt(DD3D_c(ip,j,k))
!					endif 					
					if (jm>=0) then  
						Ay3D_c(i,j,k)=Ay3D(i,j,k)/sqrt(DD3D_c(i,jm,k)) 
					endif
!					if (jp<=j1) then 
!						Cy3D_c(i,j,k)=Cy3D(i,j,k)/sqrt(DD3D_c(i,jp,k))
!					endif 									
					if (km>=0) then 
						Az3D_c(i,j,k)=Az3D(i,j,k)/sqrt(DD3D_c(i,j,km))
					endif 
!					if (kp<=k1) then 
!						Cz3D_c(i,j,k)=Cz3D(i,j,k)/sqrt(DD3D_c(i,j,kp))
!					endif 									
				enddo 
			enddo 
		enddo
		do k=kb,ke 
			do j=jb,je
				do i=ib,ie	
					DD3D_c(i,j,k)=DD3D_c(i,j,k)/sqrt(DD3D_c(i,j,k))
				enddo 
			enddo 
		enddo		
		call bound_3D(Ax3D_c)
		call bound_3D(Ay3D_c)
		call bound_3D(Az3D_c)
		call bound_3D(DD3D_c)
		! Local incomplete Cholesky matrix created, symmetric therefore both upper and lower triangle are filled 


		p=0. 
		znew=0. 
		
		!x = RHS3D !start with RHS as initial guess
		x = x0 !starting condition
		
		b2=0.
		r2=0.
		rz2=0.		
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					! residual = RHS - A*x 
!					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k)*DD3D(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
     &						    Cx3D(i  ,j  ,k  ) * x(i+1,j  ,k  ) - Ax3D(i  ,j  ,k  ) * x(i-1,j  ,k  ) -
     &							Cy3D(i  ,j  ,k  ) * x(i  ,j+1,k  ) - Ay3D(i  ,j  ,k  ) * x(i  ,j-1,k  ) -
     &							Cz3D(i  ,j  ,k  ) * x(i  ,j  ,k+1) - Az3D(i  ,j  ,k  ) * x(i  ,j  ,k-1) 
					! preconditioning in following two loops 
					b2 = b2 + RHS3D(i,j,k)*RHS3D(i,j,k) 					
					r2 = r2 + r(i,j,k)*r(i,j,k) 
				enddo
			enddo 
		enddo 
		
		! preconditioning: z = M^-1r -> Mz = r --> U^TDUz = r --> U^TDv=r solved by forward substitution followed by DUz = Dv solved by backward substitution 
		! preconditioning: Mz = r --> L^TLz = r --> L^Tv = r solved by forward substitution followed by Lz = v solved by backward substitution 
!		znew = r/DD3D_c 						!initialize all znew of ghost boundary cells 
!		call bound_intern_and_periodic(znew) 		!initialize all znew of ghost boundary cells 
		!call bound_intern_and_periodic(r) 			!initialize all znew of ghost boundary cells 
		znew = 0. 
		do k=kb,ke 
			do j=jb,je
				do i=ib,ie	!start
					im=MAX(i-1,0)
					jm=MAX(j-1,0)
					km=MAX(k-1,0) 				
					! preconditioning: z = M^-1r -> Mz = r with M=LL^T --> LL^Tz = r --> Lv=r solved by forward substitution followed by Lz = v solved by backward substitution 
					!forward substitution:
					znew(i,j,k) = (r(i,j,k)-Ax3D_c(i,j,k)*znew(im,j,k)-Ay3D_c(i,j,k)*znew(i,jm,k)-Az3D_c(i,j,k)*znew(i,j,km))/DD3D_c(i,j,k)
				enddo 
			enddo 
		enddo 
!		znew = znew / DD3D_c 						!initialize all znew of ghost boundary cells
		!call bound_intern_and_periodic(znew) 		!initialize all znew of ghost boundary cells 		
		do k=ke,kb,-1 
			do j=je,jb,-1 
				do i=ie,ib,-1 
					ip=MIN(i+1,i1)
					jp=MIN(j+1,j1)
					kp=MIN(k+1,k1) 				
					! preconditioning: z = M^-1r -> Mz = r with M=LL^T --> LL^Tz = r --> Lv=r solved by forward substitution followed by Lz = v solved by backward substitution  
					!backward substitution:
					!in Cholesky matrix Ax3D_c(i+1,j,k) is not equal to Cx3D_c(i,j,k) (due to scaling with column-diagonal) and the former needs to be used (tested in matlab)
					znew(i,j,k) = (znew(i,j,k)-Ax3D_c(ip,j,k)*znew(ip,j,k)-Ay3D_c(i,jp,k)*znew(i,jp,k)
     &					-Az3D_c(i,j,kp)*znew(i,j,kp))/DD3D_c(i,j,k)	 
					!p(i,j,k) = znew(i,j,k)
					!rz2 = rz2 + r(i,j,k)*znew(i,j,k)
				enddo 
			enddo 
		enddo		
		do k=ke,kb,-1 
			do j=je,jb,-1 
				do i=ie,ib,-1 
					p(i,j,k) = znew(i,j,k)
					rz2 = rz2 + r(i,j,k)*znew(i,j,k)
				enddo 
			enddo 
		enddo	
		
		
		
		call mpi_allreduce(b2,b2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)		
		call mpi_allreduce(rz2,rz2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)	
		IF (sqrt(r2_global).le.tol.or.sqrt(b2_global).le.tol) THEN
			!no need to iterate
			xnew=x0
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
			return 
		ENDIF
		call bound_intern_and_periodic(p)
			
		
				
		do tel=1,maxiter 
			rnew2=0.
			rznew2=0.
			pap=0.
			rnew2_global=0.
			rznew2_global=0.
			pap_global=0. 
			eps=0. 
			alpha=0. 
			beta=0.

			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						ap(i,j,k) = p(i,j,k)*DD3D(i,j,k) + 
     &								Cx3D(i  ,j  ,k  ) * p(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * p(i-1,j  ,k  ) +
     &								Cy3D(i  ,j  ,k  ) * p(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * p(i  ,j-1,k  ) +
     &								Cz3D(i  ,j  ,k  ) * p(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * p(i  ,j  ,k-1) 
						!r2  = r2  + r(i,j,k)* r(i,j,k) 
						pap = pap + p(i,j,k)*ap(i,j,k) 	 
					enddo
				enddo 
			enddo 	
			call mpi_allreduce(pap,pap_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			alpha = rz2_global/pap_global !with pre-conditioner !r2_global/pap_global
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						xnew(i,j,k)  = x(i,j,k)  + alpha* p(i,j,k)
						rnew(i,j,k)  = r(i,j,k)  - alpha*ap(i,j,k)
						rnew2 = rnew2 + rnew(i,j,k)*rnew(i,j,k)
					enddo
				enddo 
			enddo 
		znew = 0. 
		!call bound_intern_and_periodic(rnew) 
		do k=kb,ke 
			do j=jb,je
				do i=ib,ie	
					im=MAX(i-1,0)
					jm=MAX(j-1,0)
					km=MAX(k-1,0) 				
					! preconditioning: z = M^-1r -> Mz = r with M=LL^T --> LL^Tz = r --> Lv=r solved by forward substitution followed by Lz = v solved by backward substitution 
					!forward substitution:
					znew(i,j,k) = (rnew(i,j,k)-Ax3D_c(i,j,k)*znew(im,j,k)-Ay3D_c(i,j,k)*znew(i,jm,k)-Az3D_c(i,j,k)*znew(i,j,km))/DD3D_c(i,j,k)			
				enddo 
			enddo 
		enddo 
		!call bound_intern_and_periodic(znew) 
		do k=ke,kb,-1 
			do j=je,jb,-1 
				do i=ie,ib,-1 
					ip=MIN(i+1,i1)
					jp=MIN(j+1,j1)
					kp=MIN(k+1,k1) 				
					! preconditioning: z = M^-1r -> Mz = r with M=LL^T --> LL^Tz = r --> Lv=r solved by forward substitution followed by Lz = v solved by backward substitution  
					!backward substitution:
					!in Cholesky matrix Ax3D_c(i+1,j,k) is not equal to Cx3D_c(i,j,k) (due to scaling with column-diagonal) and the former needs to be used (tested in matlab)
					znew(i,j,k) = (znew(i,j,k)-Ax3D_c(ip,j,k)*znew(ip,j,k)-Ay3D_c(i,jp,k)*znew(i,jp,k)
     &					-Az3D_c(i,j,kp)*znew(i,j,kp))/DD3D_c(i,j,k)

					!rznew2 = rznew2 + rnew(i,j,k)*znew(i,j,k) 
				enddo 
			enddo 
		enddo
		do k=ke,kb,-1 
			do j=je,jb,-1 
				do i=ie,ib,-1 
					rznew2 = rznew2 + rnew(i,j,k)*znew(i,j,k) 
				enddo 
			enddo 
		enddo
		

			
			call mpi_allreduce(rnew2,rnew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			call mpi_allreduce(rznew2,rznew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			eps = sqrt(rnew2_global)/sqrt(b2_global) 
!				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dpcg iteration working; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!				ENDIF			
			IF (eps<=tol) THEN 
				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dcg iteration finished; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',
!     & 					tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
					write(*,*),'CN3Dpcg_ic finished; tel,tol,eps,b2_global:',tel,tol,eps,b2_global
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF			
				EXIT 
			ELSEIF (tel.eq.maxiter) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration did not reach desired tolerance, TUDflow3D will continue with end result'
					write(*,*)'tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF
				EXIT 
			ENDIF 
			beta = rznew2_global/rz2_global 
			IF (isnan(alpha).or.isnan(beta)) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration crashed; alpha,beta,tel,tol,eps=',alpha,beta,tel,tol,eps
				ENDIF
				EXIT 
			ENDIF

			! not finished, so continue:
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax
						p(i,j,k) = znew(i,j,k)  + beta*p(i,j,k)
						x(i,j,k) = xnew(i,j,k)
						r(i,j,k) = rnew(i,j,k)
						!p(i,j,k) = pnew(i,j,k)
					enddo
				enddo 
			enddo
			r2_global = rnew2_global 
			rz2_global = rznew2_global 
			call bound_intern_and_periodic(p)

		enddo ! interation loop 
		
		END SUBROUTINE	

		SUBROUTINE CN3Dpcg_ilu(xnew,x0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,maxiter,tol,rank) 
		! Pre-conditioned ConjungateGradient implicit solver for 3D Laplacian Ax=RHS 
		! Matrix A is a 3D Laplacian matrix containing 7 diagonals:
		! Main diagonal values of matrix A in DD3D(0:i1,0:j1,0:k1)
		! Off-diagonal values left and right in x-dir, y-dir and z-dir: Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D with size (0:i1,0:j1,0:k1) 
		! Righthandside in RHS(0:i1,0:j1,0:k1) 
		! solution output in xnew(0:i1,0:j1,0:k1) 
		! starting condition in x0(0:i1,0:j1,0:k1)
		
		implicit none
		
        include 'mpif.h'
        integer ierr,rank
		integer  i,j,k,i1,j1,k1,im,ip,jm,jp,km,kp,tel,maxiter !,imax,jmax,kmax
		integer ib,ie,jb,je,kb,ke,tel_end,imax2,jmax2,kmax2 
		real Ax3D(0:i1,0:j1,0:k1),Ay3D(0:i1,0:j1,0:k1),Az3D(0:i1,0:j1,0:k1)
		real Cx3D(0:i1,0:j1,0:k1),Cy3D(0:i1,0:j1,0:k1),Cz3D(0:i1,0:j1,0:k1)
		real DD3D(0:i1,0:j1,0:k1),RHS3D(0:i1,0:j1,0:k1)
		real x(0:i1,0:j1,0:k1),xnew(0:i1,0:j1,0:k1),p(0:i1,0:j1,0:k1),r(0:i1,0:j1,0:k1),rnew(0:i1,0:j1,0:k1) !,pnew(0:i1,0:j1,0:k1)
		real ap(0:i1,0:j1,0:k1),x0(0:i1,0:j1,0:k1)
		real tol,alpha,beta,rnew2,r2,pap,rnew2_global,r2_global,pap_global,eps,b2,b2_global
		real znew(0:i1,0:j1,0:k1),rz2_global,rz2,rznew2,rznew2_global
		integer*2 irks((i1+1)*(j1+1)*(k1+1)),jrks((i1+1)*(j1+1)*(k1+1)),krks((i1+1)*(j1+1)*(k1+1))
		real Ax3D_c(0:i1,0:j1,0:k1),Ay3D_c(0:i1,0:j1,0:k1),Az3D_c(0:i1,0:j1,0:k1)
		real Cx3D_c(0:i1,0:j1,0:k1),Cy3D_c(0:i1,0:j1,0:k1),Cz3D_c(0:i1,0:j1,0:k1)
		real DD3D_c(0:i1,0:j1,0:k1),AAA(1:(ie-ib+1)*(je-jb+1)*(ke-kb+1),1:(ie-ib+1)*(je-jb+1)*(ke-kb+1))
		

		
		! make local (on this partition) ILU(0) decomposition:
		tel=0 
		do k=kb,ke 
			do j=jb,je
				do i=ib,ie
					tel=tel+1;
					!tel2=i+(j-1)*nx+(k-1)*ny*nx;
					irks(tel)=i;
					jrks(tel)=j;
					krks(tel)=k;		
				enddo
			enddo 
		enddo
		tel_end=tel 
		
		imax2=(ie-ib+1)
		jmax2=(je-jb+1)
		kmax2=(ke-kb+1) 
		tel=1;
		AAA(1:imax2*jmax2*kmax2,1:imax2*jmax2*kmax2)=0.
		do tel = 1,imax2*jmax2*kmax2
			i=irks(tel)
			j=jrks(tel)
			k=krks(tel)
			im=i-1 
			jm=j-1
			km=k-1
			ip=i+1 
			jp=j+1
			kp=k+1		
			AAA(tel,tel)=DD3D(i,j,k) 
			!if (im>=0) then 
				AAA(tel-1,tel)=Ax3D(i,j,k)
			!endif 
			!if (ip<=i1) then 
				AAA(tel+1,tel)=Cx3D(i,j,k)
			!endif 					
			!if (jm>=0) then  
				AAA(tel-imax2,tel)=Ay3D(i,j,k) 
			!endif
			!if (jp<=j1) then 
				AAA(tel+imax2,tel)=Cy3D(i,j,k)
			!endif 						
			!if (km>=0) then 
				AAA(tel+imax2*jmax2,tel)=Az3D(i,j,k)
			!endif 
			!if (kp<=k1) then 
				AAA(tel-imax2*jmax2,tel)=Cz3D(i,j,k)
			!endif 
		enddo 
		
		write(*,*),'1 building ILU AAA matrix',rank
		
		! DD3D_c contains diagonal belonging to U; L has a diagonal with ones which is not stored 
		! first entry DD3D_c(1,1) = DD3D(1,1) 
		do i = 2,imax2*jmax2*kmax2 !do i=2,n 
			do k = 1,i-1 
				! a_ik = a_ik/a_kk
				AAA(i,k) = AAA(i,k)/AAA(k,k)
				do j=k+1,imax2*jmax2*kmax2
					! a_ij = a_ij - a_ik*a_kj
					AAA(i,j) = AAA(i,j) - AAA(i,k)*AAA(k,j) 
				enddo
			enddo 
		enddo 
		
		write(*,*),'2 building ILU AAA matrix',rank
		
		do tel = 1,imax2*jmax2*kmax2
			i=irks(tel)
			j=jrks(tel)
			k=krks(tel)
			im=i-1 
			jm=j-1
			km=k-1
			ip=i+1 
			jp=j+1
			kp=k+1			
			DD3D_c(i,j,k)=AAA(tel,tel)
			!if (im>=0) then 
				Ax3D_c(i,j,k)=AAA(tel-1,tel)
			!endif 
			!if (ip<=i1) then 
				Cx3D_c(i,j,k)=AAA(tel+1,tel)
			!endif 					
			!if (jm>=0) then  
				Ay3D_c(i,j,k)=AAA(tel-imax2,tel)
			!endif
			!if (jp<=j1) then 
				Cy3D_c(i,j,k)=AAA(tel+imax2,tel)
			!endif 						
			!if (km>=0) then 
				Az3D_c(i,j,k)=AAA(tel+imax2*jmax2,tel)
			!endif 
			!if (kp<=k1) then 
				Cz3D_c(i,j,k)=AAA(tel-imax2*jmax2,tel)
			!endif 
		enddo 
		! Local ILU(0) matrix created, DD3D_c contains diagonal belonging to U; L has a diagonal with ones which is not stored 
		write(*,*),'3 building ILU AAA matrix',rank
		
		!x = RHS3D !start with RHS as initial guess
		x = x0 !starting condition
		
		b2=0.
		r2=0.
		rz2=0.		
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					! residual = RHS - A*x 
!					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k)*DD3D(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
     &						    Cx3D(i  ,j  ,k  ) * x(i+1,j  ,k  ) - Ax3D(i  ,j  ,k  ) * x(i-1,j  ,k  ) -
     &							Cy3D(i  ,j  ,k  ) * x(i  ,j+1,k  ) - Ay3D(i  ,j  ,k  ) * x(i  ,j-1,k  ) -
     &							Cz3D(i  ,j  ,k  ) * x(i  ,j  ,k+1) - Az3D(i  ,j  ,k  ) * x(i  ,j  ,k-1) 
					! preconditioning in following two loops 
					b2 = b2 + RHS3D(i,j,k)*RHS3D(i,j,k) 					
					r2 = r2 + r(i,j,k)*r(i,j,k) 
				enddo
			enddo 
		enddo 
		
		! preconditioning: z = M^-1r -> Mz = r --> LUz = r --> Lv=r solved by forward substitution followed by Uz = v solved by backward substitution 
		znew = 0. 
		do k=kb,ke 
			do j=jb,je
				do i=ib,ie	
					!forward substitution:
					!znew(i,j,k) = (r(i,j,k)-Ax3D_c(i,j,k)*znew(i-1,j,k)-Ay3D_c(i,j,k)*znew(i,j-1,k)-Az3D_c(i,j,k)*znew(i,j,k-1)) !divided by diag which is one !/DD3D_c(i,j,k)
					znew(i,j,k) = (r(i,j,k)-Ax3D_c(i,j,k)*znew(i-1,j,k)-Ay3D_c(i,j,k)*znew(i,j-1,k)-Az3D_c(i,j,k)*znew(i,j,k-1))/DD3D_c(i,j,k)
				enddo 
			enddo 
		enddo 
		do k=ke,kb,-1 
			do j=je,jb,-1 
				do i=ie,ib,-1 
					!backward substitution:
!					znew(i,j,k) = (znew(i,j,k)-Cx3D_c(i,j,k)*znew(i+1,j,k)-Cy3D_c(i,j,k)*znew(i,j+1,k)
!     &					-Cz3D_c(i,j,k)*znew(i,j,k+1))/DD3D_c(i,j,k)
					znew(i,j,k) = (znew(i,j,k)*DD3D_c(i,j,k)-Cx3D_c(i,j,k)*znew(i+1,j,k)-Cy3D_c(i,j,k)*znew(i,j+1,k)
     &					-Cz3D_c(i,j,k)*znew(i,j,k+1))/DD3D_c(i,j,k)

					p(i,j,k) = znew(i,j,k)
					rz2 = rz2 + r(i,j,k)*znew(i,j,k)
				enddo 
			enddo 
		enddo		
		write(*,*),'4 preconditioning with ILU AAA matrix',rank
		
		
		call mpi_allreduce(b2,b2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)		
		call mpi_allreduce(rz2,rz2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)	
		IF (sqrt(r2_global).le.tol.or.sqrt(b2_global).le.tol) THEN
			!no need to iterate
			xnew=x0
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
			return 
		ENDIF
		call bound_intern_and_periodic(p)
			
		
				
		do tel=1,maxiter 
			rnew2=0.
			rznew2=0.
			pap=0.
			rnew2_global=0.
			rznew2_global=0.
			pap_global=0. 
			eps=0. 
			alpha=0. 
			beta=0.

			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						ap(i,j,k) = p(i,j,k)*DD3D(i,j,k) + 
     &								Cx3D(i  ,j  ,k  ) * p(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * p(i-1,j  ,k  ) +
     &								Cy3D(i  ,j  ,k  ) * p(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * p(i  ,j-1,k  ) +
     &								Cz3D(i  ,j  ,k  ) * p(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * p(i  ,j  ,k-1) 
						!r2  = r2  + r(i,j,k)* r(i,j,k) 
						pap = pap + p(i,j,k)*ap(i,j,k) 	 
					enddo
				enddo 
			enddo 	
			call mpi_allreduce(pap,pap_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			alpha = rz2_global/pap_global !with pre-conditioner !r2_global/pap_global
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						xnew(i,j,k)  = x(i,j,k)  + alpha* p(i,j,k)
						rnew(i,j,k)  = r(i,j,k)  - alpha*ap(i,j,k)
						rnew2 = rnew2 + rnew(i,j,k)*rnew(i,j,k)
					enddo
				enddo 
			enddo 

		znew = 0. 
		do k=kb,ke 
			do j=jb,je
				do i=ib,ie	
					! preconditioning: z = M^-1r 
					!forward substitution:
					!znew(i,j,k) = (rnew(i,j,k)-Ax3D_c(i,j,k)*znew(i-1,j,k)-Ay3D_c(i,j,k)*znew(i,j-1,k)-Az3D_c(i,j,k)*znew(i,j,k-1)) !divided by diag which is one !/DD3D_c(i,j,k)
					znew(i,j,k) = (rnew(i,j,k)-Ax3D_c(i,j,k)*znew(i-1,j,k)-Ay3D_c(i,j,k)*znew(i,j-1,k)-Az3D_c(i,j,k)*znew(i,j,k-1))/DD3D_c(i,j,k)
				enddo 
			enddo 
		enddo 
		do k=ke,kb,-1 
			do j=je,jb,-1 
				do i=ie,ib,-1 
					! preconditioning: z = M^-1r -> Mz = r --> U^TDUz = r --> U^TDv=r solved by forward substitution followed by DUz = Dv solved by backward substitution 
					
					!backward substitution:
!					znew(i,j,k) = (znew(i,j,k)-Cx3D_c(i,j,k)*znew(i+1,j,k)-Cy3D_c(i,j,k)*znew(i,j+1,k)
!     &					-Cz3D_c(i,j,k)*znew(i,j,k+1))/DD3D_c(i,j,k)	 
					znew(i,j,k) = (znew(i,j,k)*DD3D_c(i,j,k)-Cx3D_c(i,j,k)*znew(i+1,j,k)-Cy3D_c(i,j,k)*znew(i,j+1,k)
     &					-Cz3D_c(i,j,k)*znew(i,j,k+1))/DD3D_c(i,j,k)	 

					rznew2 = rznew2 + rnew(i,j,k)*znew(i,j,k) 
				enddo 
			enddo 
		enddo


			
			call mpi_allreduce(rnew2,rnew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			call mpi_allreduce(rznew2,rznew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			eps = sqrt(rnew2_global)/sqrt(b2_global) 
!				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dpcg iteration working; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!				ENDIF			
			IF (eps<=tol) THEN 
				IF (rank.eq.0) THEN 
!					write(*,*),'CN3Dcg iteration finished; tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',
!     & 					tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
					write(*,*),'CN3Dpcg_ilu finished; tel,tol,eps,b2_global:',tel,tol,eps,b2_global	 
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF			
				EXIT 
			ELSEIF (tel.eq.maxiter) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration did not reach desired tolerance, TUDflow3D will continue with end result'
					write(*,*)'tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF
				EXIT 
			ENDIF 
			beta = rznew2_global/rz2_global 
			IF (isnan(alpha).or.isnan(beta)) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration crashed; alpha,beta,tel,tol,eps=',alpha,beta,tel,tol,eps
				ENDIF
				EXIT 
			ENDIF

			! not finished, so continue:
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax
						p(i,j,k) = znew(i,j,k)  + beta*p(i,j,k)
						x(i,j,k) = xnew(i,j,k)
						r(i,j,k) = rnew(i,j,k)
						!p(i,j,k) = pnew(i,j,k)
					enddo
				enddo 
			enddo
			r2_global = rnew2_global 
			rz2_global = rznew2_global 
			call bound_intern_and_periodic(p)

		enddo ! interation loop 
		
		END SUBROUTINE		

		SUBROUTINE CN3Dpcg_pol(xnew,x0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,maxiter,tol,rank) 
		! Pre-conditioned ConjungateGradient implicit solver for 3D Laplacian Ax=RHS 
		! Matrix A is a 3D Laplacian matrix containing 7 diagonals:
		! Main diagonal values of matrix A in DD3D(0:i1,0:j1,0:k1)
		! Off-diagonal values left and right in x-dir, y-dir and z-dir: Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D with size (0:i1,0:j1,0:k1) 
		! Righthandside in RHS(0:i1,0:j1,0:k1) 
		! solution output in xnew(0:i1,0:j1,0:k1) 
		! starting condition in x0(0:i1,0:j1,0:k1)
		
		implicit none
		
        include 'mpif.h'
        integer ierr,rank
		integer  i,j,k,i1,j1,k1,im,ip,jm,jp,km,kp,tel,maxiter !,imax,jmax,kmax
		integer ib,ie,jb,je,kb,ke,contin
		real Ax3D(0:i1,0:j1,0:k1),Ay3D(0:i1,0:j1,0:k1),Az3D(0:i1,0:j1,0:k1)
		real Cx3D(0:i1,0:j1,0:k1),Cy3D(0:i1,0:j1,0:k1),Cz3D(0:i1,0:j1,0:k1)
		real DD3D(0:i1,0:j1,0:k1),RHS3D(0:i1,0:j1,0:k1)
		real x(0:i1,0:j1,0:k1),xnew(0:i1,0:j1,0:k1),p(0:i1,0:j1,0:k1),r(0:i1,0:j1,0:k1),rnew(0:i1,0:j1,0:k1) !,pnew(0:i1,0:j1,0:k1)
		real ap(0:i1,0:j1,0:k1),x0(0:i1,0:j1,0:k1)
		real tol,alpha,beta,rnew2,r2,pap,rnew2_global,r2_global,pap_global,eps,b2,b2_global
		real znew(0:i1,0:j1,0:k1),rz2_global,rz2,rznew2,rznew2_global	
		real con(3),g,z1(0:i1,0:j1,0:k1),z2(0:i1,0:j1,0:k1),g_local 

!!		! First apply diagonal pre-conditioner
!!		!! left preconditioning M^-1Ax = M^-1b with M=D -->D^-1Ax = b/D
!!		!! works for sdc, but not really lower number of iterations needed, except at end when no iterations are required  
!!		x = x0 
!!		Ax3D = Ax3D/DD3D
!!		Cx3D = Cx3D/DD3D
!!		Ay3D = Ay3D/DD3D
!!		Cy3D = Cy3D/DD3D
!!		Az3D = Az3D/DD3D
!!		Cz3D = Cz3D/DD3D
!!		RHS3D = RHS3D/DD3D

		!! central preconditioning leading to symmetric matrix:
		!A~x~=b~ with A~ = P^-1AP^-T x=P^-Tx~ and b~ = P^-1b 
		!P = sqrt(D)
		!at end x=P^-Tx~
		call bound_3D(DD3D) !--> DD3D (ip im jp jm kp km) are needed  
		x = x0*sqrt(DD3D) !starting condition for x~
		Ax3D = Ax3D/sqrt(DD3D)
		Cx3D = Cx3D/sqrt(DD3D)
		Ay3D = Ay3D/sqrt(DD3D)
		Cy3D = Cy3D/sqrt(DD3D)
		Az3D = Az3D/sqrt(DD3D)
		Cz3D = Cz3D/sqrt(DD3D)
		RHS3D = RHS3D/sqrt(DD3D)
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					im=MAX(ib,i-1)
					jm=MAX(jb,j-1)
					km=MAX(kb,k-1)
					ip=MIN(ie,i+1) 
					jp=MIN(je,j+1)
					kp=MIN(ke,k+1) 			
					!if (im>=0) then 
						Ax3D(i,j,k)=Ax3D(i,j,k)/sqrt(DD3D(im,j,k))
					!endif 
					!if (ip<=i1) then 
						Cx3D(i,j,k)=Cx3D(i,j,k)/sqrt(DD3D(ip,j,k))
					!endif 					
					!if (jm>=0) then  
						Ay3D(i,j,k)=Ay3D(i,j,k)/sqrt(DD3D(i,jm,k))
					!endif
					!if (jp<=j1) then 
						Cy3D(i,j,k)=Cy3D(i,j,k)/sqrt(DD3D(i,jp,k))
					!endif 						
					!if (km>=0) then 
						Az3D(i,j,k)=Az3D(i,j,k)/sqrt(DD3D(i,j,km))
					!endif 
					!if (kp<=k1) then 
						Cz3D(i,j,k)=Cz3D(i,j,k)/sqrt(DD3D(i,j,kp))
					!endif 				
				enddo 
			enddo 
		enddo


		r=0. 
		p=0. 
		z1=0.
		z2=0. 
		rnew=0. 
		xnew=0. 
		
		b2=0.
		r2=0. 
		rz2=0.			
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					! residual = RHS - A*x 
					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
!					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k)*DD3D(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
     &						    Cx3D(i  ,j  ,k  ) * x(i+1,j  ,k  ) - Ax3D(i  ,j  ,k  ) * x(i-1,j  ,k  ) -
     &							Cy3D(i  ,j  ,k  ) * x(i  ,j+1,k  ) - Ay3D(i  ,j  ,k  ) * x(i  ,j-1,k  ) -
     &							Cz3D(i  ,j  ,k  ) * x(i  ,j  ,k+1) - Az3D(i  ,j  ,k  ) * x(i  ,j  ,k-1) 
					!p(i,j,k) = r(i,j,k)
					b2 = b2 + RHS3D(i,j,k)*RHS3D(i,j,k) 					
					r2 = r2 + r(i,j,k)*r(i,j,k) 
				enddo
			enddo 
		enddo
		call mpi_allreduce(b2,b2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)	
		IF (sqrt(r2_global).le.tol.or.sqrt(b2_global).le.tol) THEN
			!no need to iterate
			xnew=x0
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
			return 
		ENDIF
		
		!Polynomial preconditioning 
		g_local=0.
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					g_local  = MAX(g_local,ABS(Ax3D(i,j,k))+ABS(Cx3D(i,j,k))+ABS(Ay3D(i,j,k))+ABS(Cy3D(i,j,k))
     & 						+ABS(Az3D(i,j,k))+ABS(Cz3D(i,j,k))+1.) !last value is diagonal which after scaling is one		
				enddo
			enddo 
		enddo 	
		call mpi_allreduce(g_local,g,1,mpi_real8,mpi_max,mpi_comm_world,ierr)	
		!g=2.
		con(1)=-9./4.*g 			
		con(2)=27./16.*g**2
		con(3)=-15./32.*g**3	
		call bound_intern_and_periodic(r)
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					z1(i,j,k) =	con(1)*r(i,j,k)+r(i,j,k) + 
     &						    Cx3D(i  ,j  ,k  ) * r(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * r(i-1,j  ,k  ) +
     &							Cy3D(i  ,j  ,k  ) * r(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * r(i  ,j-1,k  ) +
     &							Cz3D(i  ,j  ,k  ) * r(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * r(i  ,j  ,k-1) 
				enddo
			enddo 
		enddo
		call bound_intern_and_periodic(z1)
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					z2(i,j,k) = con(2)*r(i,j,k)+z1(i,j,k) + 
     &						    Cx3D(i  ,j  ,k  ) * z1(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * z1(i-1,j  ,k  ) +
     &							Cy3D(i  ,j  ,k  ) * z1(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * z1(i  ,j-1,k  ) +
     &							Cz3D(i  ,j  ,k  ) * z1(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * z1(i  ,j  ,k-1) 
				enddo
			enddo 
		enddo
		call bound_intern_and_periodic(z2)
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					znew(i,j,k) = con(3)*r(i,j,k)+z2(i,j,k) + 
     &						    Cx3D(i  ,j  ,k  ) * z2(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * z2(i-1,j  ,k  ) +
     &							Cy3D(i  ,j  ,k  ) * z2(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * z2(i  ,j-1,k  ) +
     &							Cz3D(i  ,j  ,k  ) * z2(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * z2(i  ,j  ,k-1) 
	 
					p(i,j,k) = znew(i,j,k)
					rz2 = rz2 + r(i,j,k)*znew(i,j,k)
				enddo
			enddo 
		enddo		
		
		call mpi_allreduce(rz2,rz2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)	
		call bound_intern_and_periodic(p)
	
				
		do tel=1,maxiter 
			rnew2=0.
			rznew2=0.
			pap=0.
			rnew2_global=0.
			rznew2_global=0.
			pap_global=0. 
			eps=0. 
			alpha=0. 
			beta=0.
			contin=0

			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						ap(i,j,k) = p(i,j,k) + !*DIAG/DIAG
						!ap(i,j,k) = p(i,j,k)*DD3D(i,j,k) + !*DIAG/DIAG
     &								Cx3D(i  ,j  ,k  ) * p(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * p(i-1,j  ,k  ) +
     &								Cy3D(i  ,j  ,k  ) * p(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * p(i  ,j-1,k  ) +
     &								Cz3D(i  ,j  ,k  ) * p(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * p(i  ,j  ,k-1) 
!						ap(i,j,k) = ap(i,j,k)/DD3D(i,j,k) !preconditioning
						!r2  = r2  + r(i,j,k)* r(i,j,k) 
						pap = pap + p(i,j,k)*ap(i,j,k) 	 
					enddo
				enddo 
			enddo 	
			
		
			!call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			call mpi_allreduce(pap,pap_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			alpha = rz2_global/pap_global !with pre-conditioner !r2_global/pap_global
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						xnew(i,j,k)  = x(i,j,k)  + alpha* p(i,j,k)
						rnew(i,j,k)  = r(i,j,k)  - alpha*ap(i,j,k)
						rnew2 = rnew2 + rnew(i,j,k)*rnew(i,j,k)
						!znew(i,j,k) = rnew(i,j,k)/DD3D(i,j,k) !pre-conditioning
						!rznew2 = rznew2 + rnew(i,j,k)*znew(i,j,k) 
					enddo
				enddo 
			enddo 	
			call mpi_allreduce(rnew2,rnew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			eps = sqrt(rnew2_global)/sqrt(b2_global) 
			IF (eps<=tol) THEN
				!check whether convergence has been achieved, due to round-off errors this could not be the case, then iteration should continue
				rnew2=0.
				rnew2_global=0. 
				call bound_intern_and_periodic(xnew)
				do i=ib,ie !1,imax 
					do j=jb,je !1,jmax 
						do k=kb,ke !1,kmax 
							! residual = RHS - A*x 
							rnew(i,j,k) = 	RHS3D(i,j,k)-xnew(i,j,k) - 
     &						    Cx3D(i  ,j  ,k  ) * xnew(i+1,j  ,k  ) - Ax3D(i  ,j  ,k  ) * xnew(i-1,j  ,k  ) -
     &							Cy3D(i  ,j  ,k  ) * xnew(i  ,j+1,k  ) - Ay3D(i  ,j  ,k  ) * xnew(i  ,j-1,k  ) -
     &							Cz3D(i  ,j  ,k  ) * xnew(i  ,j  ,k+1) - Az3D(i  ,j  ,k  ) * xnew(i  ,j  ,k-1) 
							rnew2 = rnew2 + rnew(i,j,k)*rnew(i,j,k) 
						enddo
					enddo 
				enddo			
				call mpi_allreduce(rnew2,rnew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
				eps = sqrt(rnew2_global)/sqrt(b2_global) 
				IF (eps<=tol) THEN
					IF (rank.eq.0) THEN 
!							write(*,*),'CN3Dcg iteration finished; tel,tol,eps,g,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',
!     & 						tel,tol,eps,g,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
					write(*,*),'CN3Dpcg_pol finished; tel,tol,eps,b2_global:',tel,tol,eps,b2_global	 
					ENDIF	
					xnew = xnew/sqrt(DD3D)
					EXIT 
				ELSE 
					IF (rank.eq.0) THEN 
							write(*,*),'CN3Dcg iteration continued after false stop; tel,tol,eps:',tel,tol,eps
					ENDIF					
					!continue with iteration with new rnew and rnew2_global 
					contin = 1
				ENDIF 
			ELSEIF (tel.eq.maxiter) THEN 
				xnew = xnew/sqrt(DD3D)
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration did not reach desired tolerance, TUDflow3D will continue with end result'
					write(*,*)'tel,tol,eps,g,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,g,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF
				EXIT 
			ENDIF 

			!Polynomial preconditioning 
			call bound_intern_and_periodic(rnew)
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						z1(i,j,k) =	con(1)*rnew(i,j,k)+rnew(i,j,k) + 
     &					    Cx3D(i  ,j  ,k  ) * rnew(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * rnew(i-1,j  ,k  ) +
     &						Cy3D(i  ,j  ,k  ) * rnew(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * rnew(i  ,j-1,k  ) +
     &						Cz3D(i  ,j  ,k  ) * rnew(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * rnew(i  ,j  ,k-1) 
					enddo
				enddo 
			enddo
			call bound_intern_and_periodic(z1)
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						z2(i,j,k) = con(2)*rnew(i,j,k)+z1(i,j,k) + 
     &					    Cx3D(i  ,j  ,k  ) * z1(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * z1(i-1,j  ,k  ) +
     &						Cy3D(i  ,j  ,k  ) * z1(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * z1(i  ,j-1,k  ) +
     &						Cz3D(i  ,j  ,k  ) * z1(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * z1(i  ,j  ,k-1) 
					enddo
				enddo 
			enddo
			call bound_intern_and_periodic(z2)
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						znew(i,j,k) = con(3)*rnew(i,j,k)+z2(i,j,k) + 
     &				    	Cx3D(i  ,j  ,k  ) * z2(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * z2(i-1,j  ,k  ) +
     &						Cy3D(i  ,j  ,k  ) * z2(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * z2(i  ,j-1,k  ) +
     &						Cz3D(i  ,j  ,k  ) * z2(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * z2(i  ,j  ,k-1) 

						rznew2 = rznew2 + rnew(i,j,k)*znew(i,j,k) 
					enddo
				enddo 
			enddo			
			call mpi_allreduce(rznew2,rznew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			if (contin.eq.1) then 
				beta = 0.
			else 
				beta = rznew2_global/rz2_global 
			endif 
			IF (isnan(alpha).or.isnan(beta)) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration crashed; alpha,beta,tel,tol,eps=',alpha,beta,tel,tol,eps
				ENDIF
				EXIT 
			ENDIF

			! not finished, so continue:
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax
						p(i,j,k) = znew(i,j,k)  + beta*p(i,j,k)
						x(i,j,k) = xnew(i,j,k) !should be removable; keep only x
						r(i,j,k) = rnew(i,j,k) !should be removable; keep only r
						!p(i,j,k) = pnew(i,j,k)
					enddo
				enddo 
			enddo
			r2_global = rnew2_global 
			rz2_global = rznew2_global 
			call bound_intern_and_periodic(p)

		enddo ! interation loop 
		
		END SUBROUTINE		
	
	
		SUBROUTINE CN3Dpcg_pol2(xnew,x0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,maxiter,tol,rank) 
		! Pre-conditioned ConjungateGradient implicit solver for 3D Laplacian Ax=RHS 
		! Matrix A is a 3D Laplacian matrix containing 7 diagonals:
		! Main diagonal values of matrix A in DD3D(0:i1,0:j1,0:k1)
		! Off-diagonal values left and right in x-dir, y-dir and z-dir: Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D with size (0:i1,0:j1,0:k1) 
		! Righthandside in RHS(0:i1,0:j1,0:k1) 
		! solution output in xnew(0:i1,0:j1,0:k1) 
		! starting condition in x0(0:i1,0:j1,0:k1)
		
		implicit none
		
        include 'mpif.h'
        integer ierr,rank
		integer  i,j,k,i1,j1,k1,im,ip,jm,jp,km,kp,tel,maxiter !,imax,jmax,kmax
		integer ib,ie,jb,je,kb,ke,contin
		real Ax3D(0:i1,0:j1,0:k1),Ay3D(0:i1,0:j1,0:k1),Az3D(0:i1,0:j1,0:k1)
		real Cx3D(0:i1,0:j1,0:k1),Cy3D(0:i1,0:j1,0:k1),Cz3D(0:i1,0:j1,0:k1)
		real DD3D(0:i1,0:j1,0:k1),RHS3D(0:i1,0:j1,0:k1)
		real x(0:i1,0:j1,0:k1),xnew(0:i1,0:j1,0:k1),p(0:i1,0:j1,0:k1),r(0:i1,0:j1,0:k1),rnew(0:i1,0:j1,0:k1) !,pnew(0:i1,0:j1,0:k1)
		real ap(0:i1,0:j1,0:k1),x0(0:i1,0:j1,0:k1)
		real tol,alpha,beta,rnew2,r2,pap,rnew2_global,r2_global,pap_global,eps,b2,b2_global
		real znew(0:i1,0:j1,0:k1),rz2_global,rz2,rznew2,rznew2_global	
		real con(3),g,z1(0:i1,0:j1,0:k1),z2(0:i1,0:j1,0:k1),g_local 


		r=0. 
		p=0. 
		z1=0.
		z2=0. 
		rnew=0. 
		xnew=0. 
!!		! No diagonal pre-conditioner
		b2=0.
		r2=0. 
		rz2=0.			
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					! residual = RHS - A*x 
!					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
					r(i,j,k) = 	RHS3D(i,j,k)-x(i,j,k)*DD3D(i,j,k) - !x(i,j,k)-x(i,j,k)*DD3D(i,j,k)/DD3D(i,j,k) !first x(i,j,k) is RHS, second is DIAG*1/DIAG*x(i,j,k)
     &						    Cx3D(i  ,j  ,k  ) * x(i+1,j  ,k  ) - Ax3D(i  ,j  ,k  ) * x(i-1,j  ,k  ) -
     &							Cy3D(i  ,j  ,k  ) * x(i  ,j+1,k  ) - Ay3D(i  ,j  ,k  ) * x(i  ,j-1,k  ) -
     &							Cz3D(i  ,j  ,k  ) * x(i  ,j  ,k+1) - Az3D(i  ,j  ,k  ) * x(i  ,j  ,k-1) 
					!p(i,j,k) = r(i,j,k)
					b2 = b2 + RHS3D(i,j,k)*RHS3D(i,j,k) 					
					r2 = r2 + r(i,j,k)*r(i,j,k) 
				enddo
			enddo 
		enddo
		call mpi_allreduce(b2,b2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
		call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)	
		IF (sqrt(r2_global).le.tol.or.sqrt(b2_global).le.tol) THEN
			!no need to iterate
			xnew=x0
			return 
		ENDIF
		
		!Polynomial preconditioning 
		g_local=0.
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					g_local  = MAX(g_local,ABS(Ax3D(i,j,k))+ABS(Cx3D(i,j,k))+ABS(Ay3D(i,j,k))+ABS(Cy3D(i,j,k))
     & 						+ABS(Az3D(i,j,k))+ABS(Cz3D(i,j,k))+ABS(DD3D(i,j,k))) 
				enddo
			enddo 
		enddo 	
		call mpi_allreduce(g_local,g,1,mpi_real8,mpi_max,mpi_comm_world,ierr)	
		!g=2.
		con(1)=-9./4.*g 			
		con(2)=27./16.*g**2
		con(3)=-15./32.*g**3	
		call bound_intern_and_periodic(r)
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					z1(i,j,k) =	con(1)*r(i,j,k)+r(i,j,k)*DD3D(i,j,k) + 
     &						    Cx3D(i  ,j  ,k  ) * r(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * r(i-1,j  ,k  ) +
     &							Cy3D(i  ,j  ,k  ) * r(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * r(i  ,j-1,k  ) +
     &							Cz3D(i  ,j  ,k  ) * r(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * r(i  ,j  ,k-1) 
				enddo
			enddo 
		enddo
		call bound_intern_and_periodic(z1)
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					z2(i,j,k) = con(2)*r(i,j,k)+z1(i,j,k)*DD3D(i,j,k) + 
     &						    Cx3D(i  ,j  ,k  ) * z1(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * z1(i-1,j  ,k  ) +
     &							Cy3D(i  ,j  ,k  ) * z1(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * z1(i  ,j-1,k  ) +
     &							Cz3D(i  ,j  ,k  ) * z1(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * z1(i  ,j  ,k-1) 
				enddo
			enddo 
		enddo
		call bound_intern_and_periodic(z2)
		do i=ib,ie !1,imax 
			do j=jb,je !1,jmax 
				do k=kb,ke !1,kmax 
					znew(i,j,k) = con(3)*r(i,j,k)+z2(i,j,k)*DD3D(i,j,k) + 
     &						    Cx3D(i  ,j  ,k  ) * z2(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * z2(i-1,j  ,k  ) +
     &							Cy3D(i  ,j  ,k  ) * z2(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * z2(i  ,j-1,k  ) +
     &							Cz3D(i  ,j  ,k  ) * z2(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * z2(i  ,j  ,k-1) 
	 
					p(i,j,k) = znew(i,j,k)
					rz2 = rz2 + r(i,j,k)*znew(i,j,k)
				enddo
			enddo 
		enddo		
		
		call mpi_allreduce(rz2,rz2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)	
		call bound_intern_and_periodic(p)
	
				
		do tel=1,maxiter 
			rnew2=0.
			rznew2=0.
			pap=0.
			rnew2_global=0.
			rznew2_global=0.
			pap_global=0. 
			eps=0. 
			alpha=0. 
			beta=0.
			contin=0

			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						!ap(i,j,k) = p(i,j,k) + !*DIAG/DIAG
						ap(i,j,k) = p(i,j,k)*DD3D(i,j,k) + !*DIAG/DIAG
     &								Cx3D(i  ,j  ,k  ) * p(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * p(i-1,j  ,k  ) +
     &								Cy3D(i  ,j  ,k  ) * p(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * p(i  ,j-1,k  ) +
     &								Cz3D(i  ,j  ,k  ) * p(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * p(i  ,j  ,k-1) 
!						ap(i,j,k) = ap(i,j,k)/DD3D(i,j,k) !preconditioning
						!r2  = r2  + r(i,j,k)* r(i,j,k) 
						pap = pap + p(i,j,k)*ap(i,j,k) 	 
					enddo
				enddo 
			enddo 	
			
		
			!call mpi_allreduce(r2,r2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			call mpi_allreduce(pap,pap_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			alpha = rz2_global/pap_global !with pre-conditioner !r2_global/pap_global
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						xnew(i,j,k)  = x(i,j,k)  + alpha* p(i,j,k)
						rnew(i,j,k)  = r(i,j,k)  - alpha*ap(i,j,k)
						rnew2 = rnew2 + rnew(i,j,k)*rnew(i,j,k)
						!znew(i,j,k) = rnew(i,j,k)/DD3D(i,j,k) !pre-conditioning
						!rznew2 = rznew2 + rnew(i,j,k)*znew(i,j,k) 
					enddo
				enddo 
			enddo 	
			call mpi_allreduce(rnew2,rnew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			eps = sqrt(rnew2_global)/sqrt(b2_global) 
			IF (eps<=tol) THEN
				!check whether convergence has been achieved, due to round-off errors this could not be the case, then iteration should continue
				rnew2=0.
				rnew2_global=0. 
				call bound_intern_and_periodic(xnew)
				do i=ib,ie !1,imax 
					do j=jb,je !1,jmax 
						do k=kb,ke !1,kmax 
							! residual = RHS - A*x 
							rnew(i,j,k) = 	RHS3D(i,j,k)-xnew(i,j,k)*DD3D(i,j,k) - 
     &						    Cx3D(i  ,j  ,k  ) * xnew(i+1,j  ,k  ) - Ax3D(i  ,j  ,k  ) * xnew(i-1,j  ,k  ) -
     &							Cy3D(i  ,j  ,k  ) * xnew(i  ,j+1,k  ) - Ay3D(i  ,j  ,k  ) * xnew(i  ,j-1,k  ) -
     &							Cz3D(i  ,j  ,k  ) * xnew(i  ,j  ,k+1) - Az3D(i  ,j  ,k  ) * xnew(i  ,j  ,k-1) 
							rnew2 = rnew2 + rnew(i,j,k)*rnew(i,j,k) 
						enddo
					enddo 
				enddo			
				call mpi_allreduce(rnew2,rnew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
				eps = sqrt(rnew2_global)/sqrt(b2_global) 
				IF (eps<=tol) THEN
					IF (rank.eq.0) THEN 
!							write(*,*),'CN3Dcg iteration finished; tel,tol,eps,g,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',
!     & 						tel,tol,eps,g,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
					write(*,*),'CN3Dpcg_pol2 finished; tel,tol,eps,b2_global:',tel,tol,eps,b2_global	 
					ENDIF	
					!xnew = xnew/sqrt(DD3D)
					EXIT 
				ELSE 
					IF (rank.eq.0) THEN 
							write(*,*),'CN3Dcg iteration continued after false stop; tel,tol,eps:',tel,tol,eps
					ENDIF					
					!continue with iteration with new rnew and rnew2_global 
					contin = 1
				ENDIF 
			ELSEIF (tel.eq.maxiter) THEN 
				!xnew = xnew/sqrt(DD3D)
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration did not reach desired tolerance, TUDflow3D will continue with end result'
					write(*,*)'tel,tol,eps,g,b2_global,MINVAL(DD3D),MAXVAL(DD3D):',tel,tol,eps,g,b2_global,MINVAL(DD3D),MAXVAL(DD3D)
!					write(*,*),'xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)',xnew(imax,1,1),xnew(imax,1,32),xnew(imax,1,kmax)
!					write(*,*),'RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)',RHS3D(imax,1,1),RHS3D(imax,1,32),RHS3D(imax,1,kmax)
!					write(*,*),'DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)',DD3D(imax,1,1),DD3D(imax,1,32),DD3D(imax,1,kmax)
				ENDIF
				EXIT 
			ENDIF 

			!Polynomial preconditioning 
			call bound_intern_and_periodic(rnew)
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						z1(i,j,k) =	con(1)*rnew(i,j,k)+rnew(i,j,k)*DD3D(i,j,k) + 
     &					    Cx3D(i  ,j  ,k  ) * rnew(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * rnew(i-1,j  ,k  ) +
     &						Cy3D(i  ,j  ,k  ) * rnew(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * rnew(i  ,j-1,k  ) +
     &						Cz3D(i  ,j  ,k  ) * rnew(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * rnew(i  ,j  ,k-1) 
					enddo
				enddo 
			enddo
			call bound_intern_and_periodic(z1)
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						z2(i,j,k) = con(2)*rnew(i,j,k)+z1(i,j,k)*DD3D(i,j,k) + 
     &					    Cx3D(i  ,j  ,k  ) * z1(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * z1(i-1,j  ,k  ) +
     &						Cy3D(i  ,j  ,k  ) * z1(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * z1(i  ,j-1,k  ) +
     &						Cz3D(i  ,j  ,k  ) * z1(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * z1(i  ,j  ,k-1) 
					enddo
				enddo 
			enddo
			call bound_intern_and_periodic(z2)
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax 
						znew(i,j,k) = con(3)*rnew(i,j,k)+z2(i,j,k)*DD3D(i,j,k) + 
     &				    	Cx3D(i  ,j  ,k  ) * z2(i+1,j  ,k  ) + Ax3D(i  ,j  ,k  ) * z2(i-1,j  ,k  ) +
     &						Cy3D(i  ,j  ,k  ) * z2(i  ,j+1,k  ) + Ay3D(i  ,j  ,k  ) * z2(i  ,j-1,k  ) +
     &						Cz3D(i  ,j  ,k  ) * z2(i  ,j  ,k+1) + Az3D(i  ,j  ,k  ) * z2(i  ,j  ,k-1) 

						rznew2 = rznew2 + rnew(i,j,k)*znew(i,j,k) 
					enddo
				enddo 
			enddo			
			call mpi_allreduce(rznew2,rznew2_global,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			if (contin.eq.1) then 
				beta = 0.
			else 
				beta = rznew2_global/rz2_global 
			endif 
			IF (isnan(alpha).or.isnan(beta)) THEN 
				IF (rank.eq.0) THEN 
					write(*,*),'WARNING : CN3Dpcg iteration crashed; alpha,beta,tel,tol,eps=',alpha,beta,tel,tol,eps
				ENDIF
				EXIT 
			ENDIF

			! not finished, so continue:
			do i=ib,ie !1,imax 
				do j=jb,je !1,jmax 
					do k=kb,ke !1,kmax
						p(i,j,k) = znew(i,j,k)  + beta*p(i,j,k)
						x(i,j,k) = xnew(i,j,k) !should be removable; keep only x
						r(i,j,k) = rnew(i,j,k) !should be removable; keep only r
						!p(i,j,k) = pnew(i,j,k)
					enddo
				enddo 
			enddo
			r2_global = rnew2_global 
			rz2_global = rznew2_global 
			call bound_intern_and_periodic(p)

		enddo ! interation loop 
		
		END SUBROUTINE		
	
	