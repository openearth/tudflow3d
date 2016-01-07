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

      integer  ib,ie,jb,je,kb,ke,ploc,ierr,istep_output,itag,n,t,istep_output_movie,status
      real     cput1,cput2,t_output,t_output_movie
      real     cput10a,cput10b,cput11a,cput11b
	  !real     cput1b,cput2b,t_output,t_output_movie,crate
	  !integer cput1,cput2,cr
      real     bulk ,stress
      real A(3,3)
	character(1024) :: svnversion
	character(1024) :: svnurl
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
	t_output_movie=t0_output_movie
	t_output=t0_output

      call init_transpose	
      call mkgrid
      call determine_indices_jet_in
      call fkdat
      call init_sediment
      call init_propeller
      call determine_indices_ship_in
      call bedroughness_init

      dt    = 0.5*dt_max !first timestep 50% for stable startup
      ekm   =0. 
      if (SEM.eq.1) then
	call create_init_eddies_SEM
! 	call SEM_turb_bc(Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub1old,Vb1old,Wb1old,Ub2old,Vb2old,Wb2old)
	call SEM_turb_bc
      endif
	if (restart_dir.ne.'') then
		write(*,*) 'Restart option obsolete, Dflow3d is stopped...'
		STOP
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
	if (poissolver.eq.2) THEN
	   CALL SOLVEpois_vg_init_mumps
	ENDIF
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
	 call bound_c(cold(n,:,:,:),frac(n)%c,n)
 	 call bound_c(cnew(n,:,:,:),frac(n)%c,n)
 	 call bound_c(dcdt(n,:,:,:),frac(n)%c,n)
      enddo 
      call state(cold,rold)
      call state(cnew,rnew)
      call state(dcdt,drdt)
      call bound(Uold,Vold,Wold,rold,0,0.,Ub1old,Vb1old,Wb1old,Ub2old,Vb2old,Wb2old,Ub3old,Vb3old,Wb3old)
      call bound(Unew,Vnew,Wnew,rnew,0,0.,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
      call bound_rhoU(dUdt,dVdt,dWdt,drdt,slip_bot,0.,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)

!	call stress_terms(Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
!	call diffu_com4(wx,Srr,Spr,Szr,Spp,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)

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
	wxold=wx
	wyold=wy
	wzold=wz
	ENDIF

	if (Cs>0.or.sgs_model.eq.'MixLe') then
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
          endif
	else 
	  ekm(:,:,:)=ekm_mol
	endif
      call chkdt
      time_nm=-dt
      time_n=0.
      time_np=dt
	

!      call fillps 
!      CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
	   call output_init_nc(time_np)

	  if (rank.eq.0) then		
		call cpu_time(cput1)
		!call SYSTEM_CLOCK(cput1) !SYSTEM_CLOCK works correctly in case a processor is overloaded with more partitions than physical cores
		!CALL system_clock(count_rate=cr)
		!crate = REAL(cr)	
	  endif
	

      do while (time_n.le.t_end)
	  if (rank.eq.0) then		
		call cpu_time(cput10a)
	  endif	  
	!do istep=1,nstep
		istep=istep   + 1
!		call calc_div 
		if (Cs>0.or.sgs_model.eq.'MixLe' ) then
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
		  endif
		else
		      ekm(:,:,:)=ekm_mol 
		endif
		if (SEM.eq.1) then
		  call SEM_turb_bc
		  call move_eddies_SEM
		endif
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
		call bound_rhoU(dUdt,dVdt,dWdt,drdt,slip_bot,time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)

		if (comp_filter_n>0) then
			if (mod(istep,comp_filter_n).eq.0) then
			  call compact_filter(dUdt,dVdt,dWdt,Ru,Rp,dr,phiv,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,comp_filter_a,
     &   tmax_inPpunt,i_inPpunt,j_inPpunt,tmax_inUpunt,i_inUpunt,j_inUpunt,tmax_inVpunt,i_inVpunt,j_inVpunt,kjet)
			  call bound_rhoU(dUdt,dVdt,dWdt,drdt,0,time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
			endif
		endif
		
		!IF (Hs>0.and.time_int.eq.'RK3') THEN
		 DO n=1,npresIBM
			! extra pres-corr loops with IBM boundary for better impermeable IBM objects 
			! npresIBM=0 -> order 1/1000 Umax in IBM objects (default)
			! npresIBM=10 -> order 1/100000 Umax in IBM objects
			call fillps2(dudt,dvdt,dwdt,drdt,time_np,dt)
			  IF (poissolver.eq.2) THEN
			   CALL SOLVEpois_vg_mumps(p)
			  ELSEIF (poissolver.eq.3) THEN
			   CALL SOLVEpois_vg_pardiso(p)	   
			ELSE
			   CALL SOLVEpois(p) !,Ru,Rp,DPHI,dz,rank,imax,jmax,kmax,px)
			  ENDIF
	  
			call correc2(dudt,dvdt,dwdt,dt)
			pold=p+pold    !what is called p here was dp in reality, now p is 
			call bound_rhoU(dUdt,dVdt,dWdt,drdt,0,time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
		 ENDDO
		!ENDIF

		call fillps
	  if (rank.eq.0) then		
		call cpu_time(cput11a)
	  endif	  		
	  
	    
!	  IF (poissolver.eq.1) THEN
!	   CALL SOLVEpois_vg(p)
	  IF (poissolver.eq.2) THEN
	   CALL SOLVEpois_vg_mumps(p)
	  ELSEIF (poissolver.eq.3) THEN
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

!		call bound(Uold,Vold,Wold,Rold,0,time_n,Ub1old,Vb1old,Wb1old,Ub2old,Vb2old,Wb2old)
		call bound(Unew,Vnew,Wnew,Rnew,0,time_np,Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)

		if (mod(istep,100).eq.0) then
			call chkdiv
		endif
		call chkdt

	   	rold=rnew
		rnew=drdt
		if (time_np.ge.tstart_rms) then
			call statistics(Unew,Vnew,Wnew,Cnew,Rnew)
		endif

!		if (time_np.ge.t_output) then
		if (time_np.ge.t_output.and.t_output.le.te_output+1.e-12) then
		   istep_output=istep_output+1
		   t_output=t_output+dt_output		   	
		   call output_nc(istep_output,time_np)
		endif
		if (time_np.ge.t_output_movie.and.t_output_movie.le.te_output_movie+1.e-12) then
		   istep_output_movie=istep_output_movie+1
		   t_output_movie=t_output_movie+dt_output_movie		   	
		   call output_nc_movie(istep_output_movie,time_np)
		   if (rank.eq.0) then
		     WRITE(*,'(a,i10.0,a)') ' # Movie file: ',istep_output_movie,'                   #' 	
		   endif
		endif
		if (hisfile.ne.'') then
		  call appendhis !(rank,istep,time_np,Unew,Vnew,Wnew,P,Cnew,Rnew)
		endif
		if (rank.eq.0) then
			if (mod(istep,10) .eq.0) then   
				write(*,'(a,f10.4,a,f10.4,a,f8.2,a)') ' # Time: ',time_np,' s. of ',t_end,' s ',100.*time_np/t_end,' %   #'
				WRITE(*,'(a,i10.0,a,f9.6,a)') ' # Timestep: ',istep, ' dt : ',dt,' s,           #' 			
				call cpu_time(cput10b)
		write(*,'(a,f6.3,a,f6.3,a,f5.2,a)'),' # CPU t=',NINT((cput10b-cput10a)*1000.)/1000.,'s, 1x pois=',
     &   NINT((cput11b-cput11a)*1000.)/1000.,'s = ',NINT(10000.*(cput11b-cput11a)/(cput10b-cput10a))/100.,
     &'%          #'				
			endif
		endif
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

	if (time_n.ge.tstart_rms) then
		call output_stat_nc(time_n)
	endif
		if (hisfile.ne.'') then
		  call finalize_his!(rank,istep)
		endif

      call mpi_finalize(ierr,istep)
	
	if (rank.eq.0) then
	call cpu_time(cput2)
	!call SYSTEM_CLOCK(cput2) !SYSTEM_CLOCK works correctly in case a processor is overloaded with more partitions than physical cores, but strangely sometimes results in incorrect negative runtime...
	
	CALL write_inputtxtfile
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
                if (convection.eq.'CDS6'.or.convection.eq.'HYB4'.or.convection.eq.'HYB6'.or.convection.eq.'C4A6') then
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
		WRITE(*,'(a,a,a)') ' # TUDflow3d svn revision:',TRIM(svnversion),'                  '
!		WRITE(*,'(a,a,a)') ' # TUDflow3d svn ',TRIM(svnurl),''
		WRITE(*,'(a,a,a)') ' # User                  : ',username,'                           ' 
		WRITE(*,'(a,a,a)') ' # Machine               : ',hostname,'                           ' 
		WRITE(*,*) '####################################################'		
	endif

!       stop
      end program main
