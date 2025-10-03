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

      PROGRAM main

      USE nlist
      USE sediment
	
      implicit none
      include 'mpif.h'

      integer  ib,ie,jb,je,kb,ke,ploc,ierr,itag,n,t,status,kbp,i_target 
      real     cput1,cput2,tm_temp,interpseries
      real     cput10a,cput10b,cput11a,cput11b
	  !real     cput1b,cput2b,t_output,t_output_movie,crate
	  !integer cput1,cput2,cr
      real     bulk ,pold_ref,p_ref
      real A(3,3)
	character(1024) :: gitversion
	character(1024) :: url
	character(1024) :: date_make
		CHARACTER*20 username,hostname
		INTEGER(4) istat,istat2,GETPID,HOSTNAM
	
      include 'version.inc'

      !CHARACTER (LEN=20):: SVN_REV

      call mpi_init(ierr)
      call mpi_comm_rank (MPI_COMM_WORLD,rank,ierr)
      call mpi_comm_size (MPI_COMM_WORLD,ploc,ierr)
	  
	  !write(*,*),'rank,MPI_COMM_WORLD',rank,MPI_COMM_WORLD
	  
    
      call read_namelist
      call allocate_global_vars
      ib=1 
      ie=imax
      jb=1 
      je=jmax
      kb=1 
      ke=kmax
      
	t_output=dt_output
	istep=0
      istep_output=0      
      istep_output_movie=0
	  istep_output_bpmove = 0
	t_output_movie=t0_output_movie
	t_output=t0_output

      call init_transpose	
      call mkgrid
      call determine_indices_jet_in
      !call fkdat 
      call init_sediment
      call init_propeller
      call determine_indices_ship_in
      call inflow_profiles_front_sides(0.)
	  call init_location_bedplume
	  call update_fc_local
	  call fkdat
	  call determine_kn_flow	  
	  

      dt    = MIN(dt_ini,dt_max) 
      ekm   =0. 
      if (SEM.eq.1) then
	call create_init_eddies_SEM
! 	call SEM_turb_bc(Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub1old,Vb1old,Wb1old,Ub2old,Vb2old,Wb2old)
	call SEM_turb_bc
      endif
	if (hisfile.ne.'') then
	  call init_his
	endif
	if (bcfile.ne.'') then
	  call read_bc_from_coarse_sim(bcfile)
	endif      
!	if (poissolver.eq.1) THEN
!	   CALL SOLVEpois_vg_init
!	ENDIF
!	if (poissolver.eq.2) THEN
!	   CALL SOLVEpois_vg_init_mumps
!	ENDIF
      IF (poissolver.eq.3) THEN
	    CALL SOLVEpois_vg_init_pardiso
      ENDIF
!      IF (poissolver.eq.4) THEN
!	    CALL SOLVEpois_vg_init_AMG
!		CALL SOLVEpois_vg_AMG
!      ENDIF
      IF (poissolver.eq.5) THEN
        IF (rank.eq.0) THEN
	      CALL SOLVEpois_vg_init_pardiso3D
		  ALLOCATE(rhs3(imax,jmax*px,kmax)) ! allocate rhs3 with real size only at rank=0
		ELSE
		  ALLOCATE(rhs3(1,1,1)) ! allocate rhs3 with dummy size on other ranks
        ENDIF
      ENDIF
	  
      do n=1,nfrac
	   call bound_c(cold(n,:,:,:),frac(n)%c,n,0.)
 	   call bound_c(cnew(n,:,:,:),frac(n)%c,n,0.)
 	   call bound_c(dcdt(n,:,:,:),frac(n)%c,n,0.)
      enddo 
      call state(cold,rold)
      call state(cnew,rnew)
      call state(dcdt,drdt)
      call bound(Uold,Vold,Wold,rold,MIN(0,slip_bot),0.,Ub1old,Vb1old,Wb1old,Ub2old,Vb2old,Wb2old,Ub3old,Vb3old,Wb3old)
      call bound(Unew,Vnew,Wnew,rnew,MIN(0,slip_bot),0.,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
      call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(3,slip_bot),monopile,0.,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,
     & Vb3new,Wb3new)
	 
!	call stress_terms(Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
!	call diffu_com4(wx,Srr,Spr,Szr,Spp,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)

	pold(:,:,kmax)=0.
	IF (applyVOF.eq.1.or.initPhydrostatic.eq.1) THEN 
		do k=kmax-1,1,-1
			pold(:,:,k)=pold(:,:,k+1)+(rnew(:,:,k)-rho_b)*ABS(gz)*dz
		enddo
	    if((U_b.ge.0.and.LOA<0).or.((U_TSHD-U_b).ge.0.and.LOA>0.)) THEN 
	      i_target = imax 
	    else 
	      i_target = 1 
	    endif 		
		if (rank.eq.0) then
		  pold_ref=pold(i_target ,1,k_pzero)
		endif
		call mpi_bcast(pold_ref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
		pold=pold-pold_ref	
	endif 
	pold1=pold 
	pold2=pold 
	pold3=pold
	phdt=pold
	phnew=pold 
	kbedold=kbed

	if (Cs>0.or.sgs_model.eq.'MixLe'.or.sgs_model.eq.'ReaKE'.or.sgs_model.eq.'DSmag') then
          if (sgs_model.eq.'SSmag') then
            call LES_smagorinsky(Unew,Vnew,Wnew,rnew)
          elseif (sgs_model.eq.'DSmag') then
            call LES_DSmag(Unew,Vnew,Wnew,rnew)			
          elseif (sgs_model.eq.'FSmag') then
            call LES_filteredSmagorinsky(Unew,Vnew,Wnew,rnew)
          elseif (sgs_model.eq.'SWALE') then
            call LES_WALE(Unew,Vnew,Wnew,rnew)
	  elseif (sgs_model.eq.'Sigma') then
	    call LES_Sigma(Unew,Vnew,Wnew,rnew)
	  elseif (sgs_model.eq.'MixLe') then
	    call LES_mixinglengthdamped(Unew,Vnew,Wnew,rnew)
		elseif (sgs_model.eq.'ReaKE') then
            call RealizibleKEps(Unew,Vnew,Wnew,rnew)
          endif
	else 
	  ekm(:,:,:)=ekm_mol
	endif
	

	if (Non_Newtonian.eq.1) then
		!Erwin ten Brummelhuis corrected for ekm_mol by following line, Lynyrd de Wit 6-10-2021 switched off:
		!ekm=ekm-ekm_mol				!ekm_mol is replaced by the apparent viscosity,so correcting for ekm_mol from turbulence models
		if (Rheological_model.eq.'SIMPLE') then
			call Simple_Bingham
		elseif (Rheological_model.eq.'JACOBS') then
			call Rheo_Jacobs_and_vKesteren
		elseif (Rheological_model.eq.'WINTER') then
			call Rheo_Winterwerp_and_Kranenburg
		elseif (Rheological_model.eq.'THOMAS') then
			call Rheo_Thomas
		endif
		if (PAPANASTASIOUS_m>0.) then 
		  call Bingham_Papanastasiou
		else
		  call Bingham_Fluent_manner
		endif 
	elseif (Non_Newtonian.eq.2) then
		!Erwin ten Brummelhuis corrected for ekm_mol by following line, Lynyrd de Wit 6-10-2021 switched off:
		!ekm=ekm-ekm_mol				!ekm_mol is replaced by the apparent viscosity,so correcting for ekm_mol from turbulence models
		lambda_old(:,:,:)=Lambda_init
		lambda_new(:,:,:)=0.
		call Houska_Papanastasiou(ib,ie,jb,je,kb,ke)
	endif

	IF (time_int.eq.'AB2'.or.time_int.eq.'AB3'.or.time_int.eq.'ABv') THEN
      call advecu_CDS2(wx,Uold,Vold,Wold,Rold,Ru,Rp,dr,phiv,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
      call diffu_CDS2 (wx,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
      call advecv_CDS2(wy,Uold,Vold,Wold,Rold,Ru,Rp,dr,phiv,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
      call diffv_CDS2 (wy,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
      call advecw_CDS2(wz,Uold,Vold,Wold,Rold,Ru,Rp,dr,phiv,dz,
     +            i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
      call diffw_CDS2 (wz,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
	  wx=wx+Ppropx
	  wy=wy+Ppropy
	  wz=wz+Ppropz
	wxold=wx
	wyold=wy
	wzold=wz
	ENDIF
	

	
      call chkdt
	  if (restart_dir.eq.'') then
	    time_n=0.
		mortime=0.
	  endif
	  time_nm2=time_n-2*dt 
      time_nm=time_n-dt
	  time_np=time_n+dt
	  dt_old=dt
	  
      
	  call output_init_nc(time_np)
	  call update_nvol_bedplume(time_n)
	  call update_QSc_bedplume(time_n)
	  call update_Qc_plume(time_n)

!      call fillps 
!      CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
	  if (nfrac>0) then 
		do i=1,imax 
		  do j=1,jmax 
			!zbed(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+ SUM(Clivebed(1:nfrac,i,j,kbed(i,j)))/cfixedbed*dz
			zbed_old(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
			zbed_init(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
		  enddo 
		enddo 
	  endif 
		
	  if (rank.eq.0) then		
		call cpu_time(cput1)
		!call SYSTEM_CLOCK(cput1) !SYSTEM_CLOCK works correctly in case a processor is overloaded with more partitions than physical cores
		!CALL system_clock(count_rate=cr)
		!crate = REAL(cr)	
	  endif
	  CALL write_inputtxtfile	
	  ! write machine file already at start
		CALL GETLOG (username)
		istat = GETPID()
		istat2= HOSTNAM (hostname)	  
		WRITE(*,'(a,a,a)') ' # Machine               : ',hostname,'                           ' 
	  
      do while (time_n.le.t_end)
	  if (rank.eq.0) then		
		call cpu_time(cput10a)
	  endif	  
	!do istep=1,nstep
		istep=istep   + 1
!		call calc_div 
		if (U_b_tseriesfile.ne.''.or.V_b_tseriesfile.ne.''.or.W_b_tseriesfile.ne.'') then 
			call inflow_profiles_front_sides(time_n)
		endif 
		if (tmorf_tseriesfile.ne.'') then 
		  tm_temp=interpseries(tmorf_tseries,tmorf_series,tmorf_seriesloc,time_n)
		  if (tm_temp<0.5) then !when interpolated tm_temp is < 0.5
		    tstart_morf=t_end 	!giving present time is before tstart_morf which allows for NO exchange of sediment between bed and fluid
		  else 
			tstart_morf=0. 		!giving present time is after tstart_morf which allows for exchange of sediment between bed and fluid
		  endif
		endif
		if (tmorf2_tseriesfile.ne.'') then 
		  tm_temp=interpseries(tmorf2_tseries,tmorf2_series,tmorf2_seriesloc,time_n)
		  if (tm_temp<0.5) then !when interpolated tm_temp is < 0.5
		    tstart_morf2=t_end 	!giving present time is before tstart_morf2 which allows for exchange of sediment between bed and fluid but no bedupdate 
		  else 
			tstart_morf2=0. 		!giving present time is after tstart_morf2 which allows for exchange of sediment between bed and fluid and bedupdate
		  endif
		endif	
		if (morfac_tseriesfile.ne.'') then
			morfac=interpseries(morfac_tseries,morfac_series,morfac_seriesloc,time_n)
			IF (morfac<0.) CALL writeerror(101)
		endif 
		if (morfac2_tseriesfile.ne.'') then
			morfac2=interpseries(morfac2_tseries,morfac2_series,morfac2_seriesloc,time_n)
			IF (morfac2<1.) CALL writeerror(102)
		endif
		if (time_n>tstart_morf2) then 
  		  mortime = mortime + dt*morfac2*morfac 
		endif 
		
		if (Cs>0.or.sgs_model.eq.'MixLe'.or.sgs_model.eq.'ReaKE'.or.sgs_model.eq.'DSmag') then
		  if (sgs_model.eq.'SSmag') then
			call LES_smagorinsky(Unew,Vnew,Wnew,rnew)
          elseif (sgs_model.eq.'DSmag') then
            call LES_DSmag(Unew,Vnew,Wnew,rnew)						
		  elseif (sgs_model.eq.'FSmag') then
			call LES_filteredSmagorinsky(Unew,Vnew,Wnew,rnew)
          elseif (sgs_model.eq.'SWALE') then
            call LES_WALE(Unew,Vnew,Wnew,rnew)
		  elseif (sgs_model.eq.'Sigma') then
		    call LES_Sigma(Unew,Vnew,Wnew,rnew)
		  elseif (sgs_model.eq.'MixLe') then
		    call LES_mixinglengthdamped(Unew,Vnew,Wnew,rnew)
		elseif (sgs_model.eq.'ReaKE') then
            call RealizibleKEps(Unew,Vnew,Wnew,rnew)			
		  endif
		else
		      ekm(:,:,:)=ekm_mol 
		endif
		!add rheological apparent viscosity to turbulent viscosity calculated above
		if (Non_Newtonian.eq.1) then 
		!Erwin ten Brummelhuis corrected for ekm_mol by following line, Lynyrd de Wit 6-10-2021 switched off:
		!ekm=ekm-ekm_mol				!ekm_mol is replaced by the apparent viscosity,so correcting for ekm_mol from turbulence models
			if (Rheological_model.eq.'JACOBS') then
				call Rheo_Jacobs_and_vKesteren
			elseif (Rheological_model.eq.'WINTER') then
				call Rheo_Winterwerp_and_Kranenburg
			elseif (Rheological_model.eq.'THOMAS') then
				call Rheo_Thomas
			elseif (Rheological_model.eq.'SIMPLE'.and.MAXVAL(SIMPLE_climit)>1.e-12) then
				call Simple_Bingham				
			endif
			if (PAPANASTASIOUS_m>0.) then 
			  call Bingham_Papanastasiou
			else
			  call Bingham_Fluent_manner
			endif 
		elseif (Non_Newtonian.eq.2) then
		!Erwin ten Brummelhuis corrected for ekm_mol by following line, Lynyrd de Wit 6-10-2021 switched off:
		!ekm=ekm-ekm_mol				!ekm_mol is replaced by the apparent viscosity,so correcting for ekm_mol from turbulence models
			call Houska_Papanastasiou(ib,ie,jb,je,kb,ke)
		endif
		
		if ((interaction_bed.eq.4.or.interaction_bed.eq.6).and.time_n.ge.tstart_morf2) then 
			call update_fc_local
		endif 		
		if (kn_flow_d50_multiplier>0.) call determine_kn_flow
		if (SEM.eq.1) then
		  call SEM_turb_bc
		  call move_eddies_SEM
		endif
		kbedold=kbed 
		if (time_int.eq.'AB2') then
			call adamsb2(ib,ie,jb,je,kb,ke)
		elseif (time_int.eq.'AB3') then
			call adamsb3(ib,ie,jb,je,kb,ke)
		elseif (time_int.eq.'ABv') then
			call adamsbv(ib,ie,jb,je,kb,ke)
		elseif (time_int.eq.'EE1') then
			call euler_expl(ib,ie,jb,je,kb,ke)
		elseif (time_int.eq.'RK3') then
			call RK3(ib,ie,jb,je,kb,ke)
		endif
		IF (vel_start_after_ero.eq.1) THEN 
		do i=1,imax 
			do j=1,jmax 
				do k=kbed(i,j),kbedold(i,j)-1 ! do loop only if kbed(i,j)<kbedold(i,j), so erosion has occured
				    kbp=MIN(k1,kbed(i,j)+1)
					! initialise UV in just eroded cell with UV old of cell above. This is done at both sides of eroded cell and in bound_rhoU IBM makes UV zero when it fall below zbed (f.i. when at edge of erosion hole)
					!in case bed is lowered >1 cell in one timestep (which should not occur) UV gets initialized with zero --> this is not incorrect but will not help either
					dUdt(i,j,kbed(i,j))=Unew(i,j,kbp)*rhU(i,j,kbp)
					dUdt(i-1,j,kbed(i,j))=Unew(i-1,j,kbp)*rhU(i-1,j,kbp)
					dVdt(i,j,kbed(i,j))=Vnew(i,j,kbp)*rhV(i,j,kbp)
					dVdt(i,j-1,kbed(i,j))=Vnew(i,j-1,kbp)*rhV(i,j-1,kbp)	
				enddo 
			enddo 
		enddo 
		ENDIF 
		IF (time_int.eq.'ABv') THEN 
			call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),monopile,time_np,
     & 		Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new) !bound_rhoU on rhou^* without wall shear stress, which in ABv is done as first step 
		ELSE
			call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(3,slip_bot),monopile,time_np,
     & 		Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new) !bound_rhoU on rhou^* with wall shear stress	 
		ENDIF 
		call update_QSc_bedplume(time_np)
		call update_Qc_plume(time_np)
		call update_location_bedplume
		call update_nvol_bedplume(time_np)
		IF (MAXVAL(bp(1:nbedplume)%dt_history)>0.) THEN 
		   call update_his_bedplume(time_n)
		ENDIF
		if (time_np.ge.obst_starttimes(MAX(1,nobst_file)).and.obstfile.ne.'') then
		   CALL read_obstacle(tmax_inPpuntTSHDini,tmax_inUpuntTSHDini,tmax_inVpuntTSHDini,tmax_inWpuntTSHDini) 
		endif		
		
		if (comp_filter_n>0) then
		  if (mod(istep,comp_filter_n).eq.0) then
			  call compact_filter(dUdt,dVdt,dWdt,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,comp_filter_a,
     &   tmax_inPpunt,i_inPpunt,j_inPpunt,tmax_inUpunt,i_inUpunt,j_inUpunt,tmax_inVpunt,i_inVpunt,j_inVpunt,kjet)
			call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
		  endif
		endif
		
		
		!IF (Hs>0.and.time_int.eq.'RK3') THEN
		 DO n=1,npresIBM
			! extra pres-corr loops with IBM boundary for better impermeable IBM objects 
			! npresIBM=0 -> order 1/1000 Umax in IBM objects (default)
			! npresIBM=10 -> order 1/100000 Umax in IBM objects
			IF (npresIBM_viscupdate.eq.1) THEN 
				Uold=dUdt/rhU 
				Vold=dVdt/rhV 
				Wold=dWdt/rhW
			ENDIF 
			call fillps2(dudt,dvdt,dwdt,drdt,time_np,dt)
!			  IF (poissolver.eq.2) THEN
!			   CALL SOLVEpois_vg_mumps(p)
			IF (poissolver.eq.3) THEN
			   CALL SOLVEpois_vg_pardiso(p)	   
			ELSE
			   CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
			  ENDIF
	  
			call correc2(dudt,dvdt,dwdt,dt)
			IF (npresIBM_viscupdate.eq.1) THEN 
				call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
				Uold = dUdt/rhU - Uold !fill Uold with difference present iteration and previous iteration 
				Vold = dVdt/rhV - Vold 
				Wold = dWdt/rhW - Wold 
				viscf = 0. 
				call diffu_CDS2(viscf,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
				dUdt = dUdt + dt*viscf
				viscf = 0. 
				call diffv_CDS2(viscf,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
				dVdt = dVdt + dt*viscf
				viscf = 0. 
				call diffw_CDS2(viscf,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
				dWdt = dWdt + dt*viscf
			ELSEIF (npresIBM_viscupdate.eq.2) THEN 
				call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
				Uold = dUdt/rhU!fill Uold with present iteration 
				Vold = dVdt/rhV 
				Wold = dWdt/rhW 
				! first subtract diff-contribution of present iteration, then update diff-contribution implicit 
				viscf = 0. 
				call diffu_CDS2(viscf,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
				dUdt = dUdt - dt*CNdiff_factor*viscf
				viscf = 0. 
				call diffv_CDS2(viscf,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
				dVdt = dVdt - dt*CNdiff_factor*viscf
				viscf = 0. 
				call diffw_CDS2(viscf,Uold,Vold,Wold,ib,ie,jb,je,kb,ke)
				dWdt = dWdt - dt*CNdiff_factor*viscf
				dUdt = dUdt/rhU
				dVdt = dVdt/rhV 
				dWdt = dWdt/rhW
				call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
				CALL diffuvw_CDS2_3DCNimpl 		
				dUdt = dUdt*rhU 
				dVdt = dVdt*rhV 
				dWdt = dWdt*rhW 
			ENDIF 			
			 if (continuity_solver.eq.33.or.continuity_solver.eq.34) then 
			   pold=p(1:imax,1:jmax,1:kmax)*drdt(1:imax,1:jmax,1:kmax)+pold !scale P back with rho	
			 elseif  (continuity_solver.eq.35.or.continuity_solver.eq.36) then 
			   pold=p(1:imax,1:jmax,1:kmax)*rho_b2+pold !scale P back with rho	
			 else 
			   pold=p+pold    !what is called p here was dp in reality, now p is 
			 endif
			 call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)			 
		 ENDDO
		!ENDIF

		call fillps
	  if (rank.eq.0) then		
		call cpu_time(cput11a)
	  endif	  		
	  
	    
!	  IF (poissolver.eq.1) THEN
!	   CALL SOLVEpois_vg(p)
!	  IF (poissolver.eq.2) THEN
!	   CALL SOLVEpois_vg_mumps(p)
	  IF (poissolver.eq.3) THEN
	   CALL SOLVEpois_vg_pardiso(p)	   
	  ELSEIF (poissolver.eq.5) THEN
	  
	   IF (rank>0) THEN
        call mpi_send(p,imax*jmax*kmax,MPI_REAL8,0,rank+400,MPI_COMM_WORLD,status,ierr)    !! Gather p for all ranks to rank=0
	   ENDIF
	   IF (rank.eq.0) THEN
	     rhs3(:,1:jmax,:)=p
	     DO i=1,px-1
	       call mpi_recv(rhs3(:,i*jmax+1:i*jmax+jmax,:),imax*jmax*kmax,MPI_REAL8,i,i+400,MPI_COMM_WORLD,status,ierr) !Gather p for all ranks to rank=0
         ENDDO
	     CALL SOLVEpois_vg_pardiso3D(rhs3)	
		 p=rhs3(:,1:jmax,:)
	     DO i=1,px-1
	       call mpi_send(rhs3(:,i*jmax+1:i*jmax+jmax,:),imax*jmax*kmax,MPI_REAL8,i,i+600,MPI_COMM_WORLD,status,ierr) !! Seed rhs from rank=0 to all ranks
         ENDDO		 
	   ELSE
	     call mpi_recv(p,imax*jmax*kmax,MPI_REAL8,0,rank+600,MPI_COMM_WORLD,status,ierr)    !! Seed rhs from rank=0 to all ranks
	   ENDIF
	   
	  ELSE
	   CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
	  ENDIF

	  if (rank.eq.0) then		
		call cpu_time(cput11b)
	  endif	  		
		call correc

!		call bound(Uold,Vold,Wold,Rold,MIN(0,slip_bot),time_n,Ub1old,Vb1old,Wb1old,Ub2old,Vb2old,Wb2old)
		call bound(Unew,Vnew,Wnew,Rnew,MIN(0,slip_bot),time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
		!call bound_rhoU(dUdt,dVdt,dWdt,drdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new) !bound_rhoU on rhoU^n+1 
		!extra call bound_rhoU only needed for determination rhU,cU with split_rho_cont.eq.'VL2' based on direction of rhoU^n+1 instead of rhoU^* 
		!nov-2018 not used because 1) consistent with splitting rho of from rhoU^* 2) faster
		!only drawback is that it is not fully consistent with update C as this is based on c_edge based on direction U^n+1
		if (time_np.ge.tstart_morf2) then 
			if (bedupdatefile.ne.'') then ! start using b_update from bedupdatefile again 
				b_update = b_update_bu
			else  
				b_update=0. 	
				b_update(istart_morf2(1):i1,0:j1)=1.
				if (istart_morf1(2)>0) b_update(istart_morf1(2):i1,0:j1)=0.
				if ((istart_morf2(1)-istart_morf1(1))>0) then 
					do i=istart_morf1(1),istart_morf2(1)
					  b_update(i,0:j1)=DBLE(i-istart_morf1(1))/DBLE(istart_morf2(1)-istart_morf1(1)) !linear grow 0 -> 1
					enddo	
				endif 
				if ((istart_morf1(2)-istart_morf2(2))>0) then 		
					do i=istart_morf2(2),istart_morf1(2)
					  b_update(i,0:j1)=1.-DBLE(i-istart_morf2(2))/DBLE(istart_morf1(2)-istart_morf2(2)) !linear decrease 1 -> 0
					enddo
				endif
				b_update(0:1,0:j1)=0. ! no bed-update at inflow
				b_update(imax:i1,0:j1)=0. ! no bed-update at outflow
			endif 
		else 
			b_update=0. !force b_update to 0 because time_np is before tstart_morf2 
		endif 
		if (mod(istep,100).eq.0) then
			call chkdiv
			if (pres_in_predictor_step_internal.eq.2) then 
				pres_in_predictor_step = 0 ! use no Pold in predictor step every 100 time steps to remove spurious strange pressure relicts in immersed boundary zero flow zones
			endif 
		else
		  if (pres_in_predictor_step_internal.ne.3.or.pres_in_predictor_step_internal.ne.4.or.pres_in_predictor_step_internal.ne.5
     &		  .or.pres_in_predictor_step_internal.ne.6.or.pres_in_predictor_step_internal.ne.7) then 
			pres_in_predictor_step = pres_in_predictor_step_internal
		  endif 
		endif
		if (pres_in_predictor_step_internal.eq.3) then !spatial filter for Pold to remove spikes especially round obstacles/bed every timestep; use pres_in_predictor_step=1 in the rest of the code
		  dp(1:imax,1:jmax,1:kmax)=pold+p 
		  call bound_3D(dp)
		  do i=1,imax
			do j=1,jmax
			  do k=1,kmax !MIN(kmax,kbed(i,j))
			    pold(i,j,k) = (dp(i-1,j,k)+dp(i+1,j,k)+dp(i,j-1,k)+dp(i,j+1,k)+dp(i,j,k)+dp(i,j,k-1)+dp(i,j,k+1))/7.-p(i,j,k) ! filter pold 
			  enddo
			enddo	
		  enddo	
		elseif (pres_in_predictor_step_internal.eq.4) then !spatial filter for Pold only in bed to remove spikes especially; use pres_in_predictor_step=1 in the rest of the code
		  dp(1:imax,1:jmax,1:kmax)=pold+p 
		  call bound_3D(dp)
		  do i=1,imax
			do j=1,jmax
			  do k=1,MIN(kmax,kbed(i,j))
			    pold(i,j,k) = (dp(i-1,j,k)+dp(i+1,j,k)+dp(i,j-1,k)+dp(i,j+1,k)+dp(i,j,k)+dp(i,j,k-1)+dp(i,j,k+1))/7.-p(i,j,k) ! filter pold 
			  enddo
			enddo	
		  enddo
		elseif (pres_in_predictor_step_internal.eq.5) then !
		  if (mod(istep,1000).eq.0) then
		    pold=-p !forcing pold=0 in in predictor step in solve.f
		  else
			  dp(1:imax,1:jmax,1:kmax)=pold+p 
			  call bound_3D(dp)
			  do i=1,imax
				do j=1,jmax
				  do k=1,kmax !MIN(kmax,kbed(i,j))
					pold(i,j,k) = (dp(i-1,j,k)+dp(i+1,j,k)+dp(i,j-1,k)+dp(i,j+1,k)+dp(i,j,k)+dp(i,j,k-1)+dp(i,j,k+1))/7.-p(i,j,k) ! filter pold 
				  enddo
				enddo	
			  enddo
		  endif 
		elseif (pres_in_predictor_step_internal.eq.6) then !
			  dp(1:imax,1:jmax,1:kmax)=pold+p 
			  call bound_3D(dp)
			  do i=1,imax
				do j=1,jmax
				  do k=1,kmax !MIN(kmax,kbed(i,j))
					pold(i,j,k) = (dp(i-1,j,k)+dp(i+1,j,k)+2.*dp(i,j,k))/4.-p(i,j,k)  ! i-dir
				  enddo
				enddo	
			  enddo
			  dp(1:imax,1:jmax,1:kmax)=pold+p
			  call bound_3D(dp)
			  do i=1,imax
				do j=1,jmax
				  do k=1,kmax !MIN(kmax,kbed(i,j))
					pold(i,j,k) = (dp(i,j,k-1)+dp(i,j,k+1)+2.*dp(i,j,k))/4.-p(i,j,k)   !k-dir 
				  enddo
				enddo	
			  enddo	
			  dp(1:imax,1:jmax,1:kmax)=pold+p
			  call bound_3D(dp)
			  do i=1,imax
				do j=1,jmax
				  do k=1,kmax !MIN(kmax,kbed(i,j))
					pold(i,j,k) = (dp(i,j-1,k)+dp(i,j+1,k)+2.*dp(i,j,k))/4.-p(i,j,k)  !j-dir 
				  enddo
				enddo	
			  enddo	
		elseif (pres_in_predictor_step_internal.eq.7) then !
		  if (mod(istep,1000).eq.0) then
		    pold=-p !forcing pold=0 in in predictor step in solve.f
		  else
			  dp(1:imax,1:jmax,1:kmax)=pold+p 
			  call bound_3D(dp)
			  do i=1,imax
				do j=1,jmax
				  do k=1,kmax !MIN(kmax,kbed(i,j))
					pold(i,j,k) = (dp(i-1,j,k)+dp(i+1,j,k)+2.*dp(i,j,k))/4.-p(i,j,k)  ! i-dir
				  enddo
				enddo	
			  enddo
			  dp(1:imax,1:jmax,1:kmax)=pold+p
			  call bound_3D(dp)
			  do i=1,imax
				do j=1,jmax
				  do k=1,kmax !MIN(kmax,kbed(i,j))
					pold(i,j,k) = (dp(i,j,k-1)+dp(i,j,k+1)+2.*dp(i,j,k))/4.-p(i,j,k)   !k-dir 
				  enddo
				enddo	
			  enddo	
			  dp(1:imax,1:jmax,1:kmax)=pold+p
			  call bound_3D(dp)
			  do i=1,imax
				do j=1,jmax
				  do k=1,kmax !MIN(kmax,kbed(i,j))
					pold(i,j,k) = (dp(i,j-1,k)+dp(i,j+1,k)+2.*dp(i,j,k))/4.-p(i,j,k)  !j-dir 
				  enddo
				enddo	
			  enddo	
		  endif 		  
		endif 
		  !!    Make pressure zero at one point in outflow:
		  if((U_b.ge.0.and.LOA<0).or.((U_TSHD-U_b).ge.0.and.LOA>0.)) THEN 
			i_target = imax 
		  else 
			i_target = 1 
		  endif 		
		  if (rank.eq.0) then
			pold_ref=pold(i_target ,1,k_pzero)
		  endif
		  call mpi_bcast(pold_ref,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
		  pold=pold-pold_ref			
		if (Apvisc_force_eq.eq.1) then 
			Ppropx=0. !initialize zero every time step, in chkdt driving force periodic sims is determined and put in Ppropx,Ppropy
			Ppropy=0. 	
		endif 
		call chkdt

	   	rold=rnew
		rnew=drdt
		if (time_np.ge.tstart_rms.and.time_np.le.te_rms) then
			call statistics
		endif
		if (time_np.ge.te_rms.and.time_n.lt.te_rms) then ! write output statistics at te_rms to have them available also when simulation crashes after te_rms
			call output_stat_nc(time_n)		
			write(*,*),'RMS and AVG output written to file'			
		endif		

!		if (time_np.ge.t_output) then
		if (output_times_morphology.eq.1) then 
			if (mortime.ge.t_output.and.t_output.le.te_output+1.e-12) then
			   istep_output=istep_output+1
			   t_output=t_output+dt_output		   	
			   call output_nc('flow3D_',istep_output,time_np)
			endif
		else
			if (time_np.ge.t_output.and.t_output.le.te_output+1.e-12) then
			   istep_output=istep_output+1
			   t_output=t_output+dt_output		   	
			   call output_nc('flow3D_',istep_output,time_np)
			endif
		endif 
		if (output_times_morphology.eq.1) then
			if (mortime.ge.t_output_movie.and.t_output_movie.le.te_output_movie+1.e-12) then
			   istep_output_movie=istep_output_movie+1
			   t_output_movie=t_output_movie+dt_output_movie		   	
			   call output_nc_movie(istep_output_movie,time_np)
			   if (rank.eq.0) then
				 WRITE(*,'(a,i10.0,a)') ' # Movie file: ',istep_output_movie,'                   #' 	
			   endif
			endif
		else 
			if (time_np.ge.t_output_movie.and.t_output_movie.le.te_output_movie+1.e-12) then
			   istep_output_movie=istep_output_movie+1
			   t_output_movie=t_output_movie+dt_output_movie		   	
			   call output_nc_movie(istep_output_movie,time_np)
			   if (rank.eq.0) then
				 WRITE(*,'(a,i10.0,a)') ' # Movie file: ',istep_output_movie,'                   #' 	
			   endif
			endif
		endif
		if (hisfile.ne.'') then
		  call appendhis !(rank,istep,time_np,Unew,Vnew,Wnew,P,Cnew,Rnew)
		endif
		if (rank.eq.0) then
			if (mod(istep,10) .eq.0) then  
				call cpu_time(cput10b)			
				CALL writeprogress(time_np,t_end,istep,dt,cput1,cput10a,cput10b,cput11a,cput11b,trestart)
			endif
		endif
		time_nm2 = time_nm 
		time_nm = time_n
		time_n = time_np
		time_np  = time_np  + dt
	  if (rank.eq.0) then		

	  endif	  
      enddo
!	if (time_n.ge.t_output) then
!	   istep_output=istep_output+1
!	   t_output=t_output+dt_output		   	
!	   call output(istep_output,time_n)
!	endif

	if (time_n.ge.tstart_rms.and.te_rms.gt.t_end) then ! write output statistics only when not already written
		call output_stat_nc(time_n)
		write(*,*),'RMS and AVG output written to file'
	endif
		if (hisfile.ne.'') then
		  call finalize_his!(rank,istep)
		  write(*,*),'HIS output written to file'
		endif

      call mpi_finalize(ierr,istep)
	
	if (rank.eq.0) then
		IF (MAXVAL(bp(1:nbedplume)%dt_history)>0.) THEN 
		   call output_his_bedplume
		ENDIF	
	
	call cpu_time(cput2)
	!call SYSTEM_CLOCK(cput2) !SYSTEM_CLOCK works correctly in case a processor is overloaded with more partitions than physical cores, but strangely sometimes results in incorrect negative runtime...
	
!		cput2b=REAL(cput2)/crate
!		cput1b=REAL(cput1)/crate
		
		CALL GETLOG (username)
		istat = GETPID()
		istat2= HOSTNAM (hostname)	


		WRITE(*,*) '####################################################'
		WRITE(*,'(a,i12.0,a)') ' # Timesteps : ',istep,'                         #' 
		WRITE(*,'(a,f9.6,a,f9.6,a)') ' # dt end : ',dt,' s, dt avg :',time_n/istep,' s        #' 
		WRITE(*,'(a,a,a)') ' # Time integration : ',time_int,'                           #' 
		WRITE(*,'(a,a,a)') ' # Convection integration : ',convection,'                    #' 
                if (convection.eq.'CDS6'.or.convection.eq.'CDS4'.or.convection.eq.'HYB6'.or.convection.eq.'C4A6') then
                        WRITE(*,'(a,f9.6,a)') ' # Blend factor                  : ',numdiff/2.,'        #'
		endif
		WRITE(*,'(a,a,a)') ' # Diffusion integration : ',diffusion,'                     #' 
!		WRITE(*,'(a,f9.6,a)') ' # Compact filter alpha : ',comp_filter_a,'                 #' 
!		WRITE(*,'(a,i5.1,a)') ' # Compact filter every # steps : ',comp_filter_n,'             #' 
                WRITE(*,'(a,a,a)') ' # SGS model : ',sgs_model,'                                #'
                WRITE(*,'(a,f9.6,a)') ' # Cs : ',Cs,'                                   #'
		WRITE(*,'(a,i12.0,a,i4.0,a)') ' # Grid : ',imax*jmax*kmax*px,' cells on ',px,' cpu            #'
		WRITE(*,'(a,f12.2,a)') ' # Simulation time : ',time_n,' s                 #'				
		WRITE(*,'(a,i12,a,f12.2,a)') ' # Run time : ',floor((cput2-cput1)/3600),' h ',
     &          (cput2-cput1-floor((cput2-cput1)/3600)*3600)/60,' m         #' 				
		WRITE(*,*) '####################################################'
		WRITE(*,'(a,a)') ' # TUDflow3d git version : ',TRIM(gitversion)
		WRITE(*,'(a,a)') ' # TUDflow3d url : ',TRIM(url)
		WRITE(*,'(a,a)') ' # TUDflow3d date make : ',TRIM(date_make)		
		WRITE(*,'(a,a,a)') ' # User                  : ',username,'                           ' 
		WRITE(*,'(a,a,a)') ' # Machine               : ',hostname,'                           ' 
		WRITE(*,*) '####################################################'		
	endif

!       stop
      end program main
