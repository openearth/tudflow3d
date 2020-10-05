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

      subroutine output_init_nc(tt)
      USE nlist
      USE netcdf

	implicit none


!       include 'param.txt'
!       include 'common.txt'
      include 'mpif.h'
	real x(1:imax,1:jmax,1:kmax)
	real y(1:imax,1:jmax,1:kmax)
	real z(1:imax,1:jmax,1:kmax)
	real uu(1:imax,1:jmax,1:kmax)
	real vv(1:imax,1:jmax,1:kmax),obstacle(1:imax,1:jmax,1:kmax)
	real mass_bed(1:nfrac,1:imax,1:jmax)
	real tt,ddxx(1:imax,1:jmax),ddyy(1:imax,1:jmax)
	integer tel,n,ios,j2,r,status(MPI_STATUS_SIZE),ierr
	character*60 FILE_NAME
	character*60 strng
	character*3 varname
	character*200 tline
	logical(4) res
	!real*8 fc_local(1:imax,1:jmax,1:kmax),fc_local_vec(imax*jmax*kmax),fc_global_vec(imax*jmax*px*kmax)
!	integer(4) ps

       ! We are writing 3D data, a nx x ny x nz grid.
       integer, parameter :: NDIMS = 3
       integer, parameter :: NDIMS2 = 4
       integer, parameter :: NDIMS3 = 1
	   integer, parameter :: NDIMS5 = 2
!       integer, parameter :: NX = imax, NY = jmax, NZ = kmax
     
       ! When we create netCDF files, variables and dimensions, we get back
       ! an ID for each one.
       integer :: ncid, varid1,varid2,varid3, varid4, varid5, varid6, varid7, varid8, dimids(NDIMS), dimids2(NDIMS2)
       integer :: x_dimid, y_dimid, z_dimid, nfrac_dimid, par_dimid,dimids3(NDIMS3),dimids5(NDIMS5)
       integer :: dimids4(NDIMS),varid20,varid21,t
	   integer :: varid22,varid23	   
	character(1024) :: svnversion
	character(1024) :: svnurl
      include 'version.inc'
     
!       ! This is the data array we will write. It will just be filled with
!       ! a progression of integers for this example.
!       integer, dimension(:,:), allocatable :: data_out
     
	tel=0
      do k=1,kmax
		do j=1,jmax
			do i=1,imax
				tel=tel+1
				x(i,j,k)=Rp(i)*cos_u(j)-schuif_x !pie-piece is shifted left to make (x,y)_jet = (0,0)
				y(i,j,k)=Rp(i)*sin_u(j)
				z(i,j,k)=(k-0.5)*dz-bc_obst_h
			enddo
		enddo
	enddo
	do j=1,jmax
		do i=1,imax
			ddxx(i,j)=dr(i)
			ddyy(i,j)=Rp(i)*(phiv(j)-phiv(j-1))
		enddo
	enddo
	
	WRITE(strng,'(a,a)')'mkdir ',TRIM(inpfile)

	CALL SYSTEM(strng)
	CALL Chdir(TRIM(inpfile))
	
	WRITE(*,'(a,i4.4,a)')'mesh3D_',INT(rank),'.nc'
	WRITE(FILE_NAME,'(a,i4.4,a)')'mesh3D_',INT(rank),'.nc'


       ! Always check the return code of every netCDF function call. In
       ! this example program, wrapping netCDF calls with "call check()"
       ! makes sure that any return which is not equal to nf90_noerr (0)
       ! will print a netCDF error message and exit.
     
       ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
       ! overwrite this file, if it already exists.
       call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )
     
       ! Define the dimensions. NetCDF will hand back an ID for each.
       call check( nf90_def_dim(ncid, "xdim", imax, x_dimid) )
       call check( nf90_def_dim(ncid, "ydim", jmax, y_dimid) )
       call check( nf90_def_dim(ncid, "zdim", kmax, z_dimid) )
     
       ! The dimids array is used to pass the IDs of the dimensions of
       ! the variables. Note that in fortran arrays are stored in
       ! column-major format.
       dimids =  (/ x_dimid, y_dimid, z_dimid /)
	   dimids5 =  (/ x_dimid, y_dimid /)
     
       ! Define the variable. The type of the variable in this case is
       ! NF90_DOUBLE (4-byte double).
       call check( nf90_def_var(ncid, "x", NF90_REAL, dimids, varid1) )
       call check( nf90_put_att(ncid, varid1, 'units', 'm') )
       call check( nf90_put_att(ncid, varid1, 'long_name', 'local x coordinate C-point') )
!       call check( nf90_put_att(ncid, varid1, 'axis', 'x') )
       call check( nf90_def_var(ncid, "y", NF90_REAL, dimids, varid2) )
       call check( nf90_put_att(ncid, varid2, 'units', 'm') )
       call check( nf90_put_att(ncid, varid2, 'long_name', 'local y coordinate C-point') )
!       call check( nf90_put_att(ncid, varid2, 'axis', 'y') )
       call check( nf90_def_var(ncid, "z", NF90_REAL, dimids, varid3) )
       call check( nf90_put_att(ncid, varid3, 'units', 'm') )
       call check( nf90_put_att(ncid, varid3, 'long_name', 'local z coordinate C-point') )
!       call check( nf90_put_att(ncid, varid3, 'axis', 'z') )

       call check( nf90_def_var(ncid, "dx", NF90_REAL, dimids5, varid22) )
       call check( nf90_put_att(ncid, varid22, 'units', 'm') )
       call check( nf90_put_att(ncid, varid22, 'long_name', 'grid size dx C-grid') )
	   call check( nf90_def_var(ncid, "dy", NF90_REAL, dimids5, varid23) )
       call check( nf90_put_att(ncid, varid23, 'units', 'm') )
       call check( nf90_put_att(ncid, varid23, 'long_name', 'grid size dy C-grid') )
	   
	   obstacle=1.-fc_global(1:imax,1+rank*jmax:jmax+rank*jmax,1:kmax)
       call check( nf90_def_var(ncid, "obstacle", NF90_REAL, dimids, varid21) )
       call check( nf90_put_att(ncid, varid21, 'units', '-') )
       call check( nf90_put_att(ncid, varid21, 'long_name', 'If cell is in obstacle 1 else 0') )
	   

	! also add svn info in output files:
       CALL check( nf90_put_att(ncid,nf90_global, "svnversion", trim(svnversion)))
       CALL check( nf90_put_att(ncid,nf90_global, "svnurl", trim(svnurl)))
     
       ! End define mode. This tells netCDF we are done defining metadata.
       call check( nf90_enddef(ncid) )
     
       ! Write the pretend data to the file. Although netCDF supports
       ! reading and writing subsets of data, in this case we write all the
       ! data in one operation.
       call check( nf90_put_var(ncid, varid1, x) )
       call check( nf90_put_var(ncid, varid2, y) )
       call check( nf90_put_var(ncid, varid3, z) )
       call check( nf90_put_var(ncid, varid22, ddxx) )
       call check( nf90_put_var(ncid, varid23, ddyy) )
       call check( nf90_put_var(ncid, varid21, obstacle) )

	

	 
	 
       ! Close the file. This frees up any internal netCDF resources
       ! associated with the file, and flushes any buffers.
       call check( nf90_close(ncid) )
     

!	tel=0
!      do k=1,kmax
!		do j=1,jmax
!			do i=1,imax
!				tel=tel+1
!				uu(i,j,k)=unew(i,j,k)*cos_u(j)-vnew(i,j,k)*sin_v(j) !unew(i,j,k)
!				vv(i,j,k)=vnew(i,j,k)*cos_v(j)+unew(i,j,k)*sin_u(j) !vnew(i,j,k)
!			enddo
!		enddo
!	enddo
!	do i=1,imax
!	   do j=1,jmax
!	     	do n=1,nfrac
!			mass_bed(n,i,j) = Cnewbot(n,i,j)*dz*frac(n)%rho ! Cnewbot(n,i,j)*dz*dr(i)*Rp(i)*dphi*frac(n)%rho
!		enddo		
!	   enddo
!	enddo
!	
!	WRITE(*,'(a,i9.9,a,i4.4,a)')'flow3D_',INT(0),'_',INT(rank),'.nc'
!	WRITE(FILE_NAME,'(a,i9.9,a,i4.4,a)')'flow3D_',INT(0),'_',INT(rank),'.nc'

!       ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
!       ! overwrite this file, if it already exists.
!       call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )
!     
!       ! Define the dimensions. NetCDF will hand back an ID for each.
!       call check( nf90_def_dim(ncid, "xdim", imax, x_dimid) )
!       call check( nf90_def_dim(ncid, "ydim", jmax, y_dimid) )
!       call check( nf90_def_dim(ncid, "zdim", kmax, z_dimid) )
!       if (nfrac>0) then
!	call check( nf90_def_dim(ncid, "nfracdim", nfrac, nfrac_dimid) )
!	endif
!       call check( nf90_def_dim(ncid, "pardim", 1, par_dimid) )
!     
!       ! The dimids array is used to pass the IDs of the dimensions of
!       ! the variables. Note that in fortran arrays are stored in
!       ! column-major format.
!       dimids =  (/ x_dimid, y_dimid, z_dimid /)
!       if (nfrac>0) then
!       	dimids2 =  (/ nfrac_dimid, x_dimid, y_dimid, z_dimid /)
!	dimids4 =  (/ nfrac_dimid, x_dimid, y_dimid /)
!	endif
!       dimids3 =  (/ par_dimid  /) 
!    
!       ! Define the variable. The type of the variable in this case is
!       ! NF90_DOUBLE (8-byte double). NF90_REAL (4-byte real)
!       call check( nf90_def_var(ncid, "U", NF90_REAL, dimids, varid1) )
!       call check( nf90_put_att(ncid, varid1, 'units', 'm/s') )
!       call check( nf90_put_att(ncid, varid1, 'long_name', 'U velocity') )

!       call check( nf90_def_var(ncid, "V", NF90_REAL, dimids, varid2) )
!       call check( nf90_put_att(ncid, varid2, 'units', 'm/s') )
!       call check( nf90_put_att(ncid, varid2, 'long_name', 'V velocity') )

!       call check( nf90_def_var(ncid, "W", NF90_REAL, dimids, varid3) )
!       call check( nf90_put_att(ncid, varid3, 'units', 'm/s') )
!       call check( nf90_put_att(ncid, varid3, 'long_name', 'W velocity') )

!       if (nfrac>0) then
!       call check( nf90_def_var(ncid, "C", NF90_REAL, dimids2, varid4) )
!       call check( nf90_put_att(ncid, varid4, 'units', '-') )
!       call check( nf90_put_att(ncid, varid4, 'long_name', 'Volume concentration for each fraction') )
!       call check( nf90_def_var(ncid, "mass_bed", NF90_REAL, dimids4, varid20) )
!       call check( nf90_put_att(ncid, varid4, 'units', 'kg/m2') )
!       call check( nf90_put_att(ncid, varid4, 'long_name', 'Mass per m2 sediment fractions in bed') )
!	endif
!       call check( nf90_def_var(ncid, "rho", NF90_REAL, dimids, varid5) )
!       call check( nf90_put_att(ncid, varid5, 'units', 'kg/m3') )
!       call check( nf90_put_att(ncid, varid5, 'long_name', 'Mixture density') )
!     
!       call check( nf90_def_var(ncid, "mu_t", NF90_REAL, dimids, varid6) )
!       call check( nf90_put_att(ncid, varid6, 'units', 'kg/(sm)') )
!       call check( nf90_put_att(ncid, varid6, 'long_name', 'Dynamic eddy viscosity') )

!       call check( nf90_def_var(ncid, "P", NF90_REAL, dimids, varid7) )
!       call check( nf90_put_att(ncid, varid7, 'units', 'Pa') )
!       call check( nf90_put_att(ncid, varid7, 'long_name', 'Pressure') )

!       call check( nf90_def_var(ncid, "time", NF90_REAL, dimids3, varid8) )
!       call check( nf90_put_att(ncid, varid8, 'units', 's') )
!       call check( nf90_put_att(ncid, varid8, 'long_name', 'Time from start simulation') )

!	! also add svn info in output files:
!       CALL check( nf90_put_att(ncid,nf90_global, "svnversion", trim(svnversion)))
!       CALL check( nf90_put_att(ncid,nf90_global, "svnurl", trim(svnurl)))

!       ! End define mode. This tells netCDF we are done defining metadata.
!       call check( nf90_enddef(ncid) )
!     
!       ! Write the pretend data to the file. Although netCDF supports
!       ! reading and writing subsets of data, in this case we write all the
!       ! data in one operation.
!       call check( nf90_put_var(ncid, varid1, uu(1:imax,1:jmax,1:kmax)) )
!       call check( nf90_put_var(ncid, varid2, vv(1:imax,1:jmax,1:kmax)) )
!       call check( nf90_put_var(ncid, varid3, Wnew(1:imax,1:jmax,1:kmax)) )
!       if (nfrac>0) then
!        call check( nf90_put_var(ncid, varid4, Cnew(1:nfrac,1:imax,1:jmax,1:kmax)) )
!	call check( nf90_put_var(ncid, varid20, mass_bed(1:nfrac,1:imax,1:jmax)) )
!	endif
!       call check( nf90_put_var(ncid, varid5, rnew(1:imax,1:jmax,1:kmax)) )
!       call check( nf90_put_var(ncid, varid6, ekm(1:imax,1:jmax,1:kmax)) )
!       call check( nf90_put_var(ncid, varid7, p(1:imax,1:jmax,1:kmax)+pold(1:imax,1:jmax,1:kmax)) )
!       call check( nf90_put_var(ncid, varid8, tt) )
!     
!       ! Close the file. This frees up any internal netCDF resources
!       ! associated with the file, and flushes any buffers.
!       call check( nf90_close(ncid) )

	end

      subroutine write_inputtxtfile
      USE nlist
      USE netcdf

	implicit none

      include 'mpif.h'

	integer ios
	character*120 FILE_NAME
	character*200 tline
	character(1024) :: svnversion
	character(1024) :: svnurl
      include 'version.inc'

	!! read input file line by line and place in txt file in output dir to always be able to find input:
	IF (rank.eq.0) THEN
	WRITE(FILE_NAME,'(a,a)')'../Dflow3d.',TRIM(inpfile)
	OPEN(1,FILE=FILE_NAME,IOSTAT=ios,ACTION='read')
	IF (ios/=0) THEN
	  write(*,*) 'input file:',FILE_NAME,' to write full input into outputdir does not exist'
	  CALL writeerror(200)
	ENDIF
	OPEN(2,FILE='Dflow3d.input_of_this_run',IOSTAT=ios,ACTION='write')	


	WRITE(2,'(a,a)'),'Dflow3d svn version:',TRIM(svnversion)
	WRITE(2,'(a,a)'),'Dflow3d svn ',TRIM(svnurl)
	ios=0
	DO WHILE (ios.eq.0)
		READ(1,'(a)',iostat=ios) tline  ! read full inputfile 
		write(2,'(a)'),TRIM(tline)      ! write into input txt file
	ENDDO
	CLOSE(1)
	CLOSE(2)

	ENDIF

	END 

      subroutine output_nc(fname_basis,istap,tt)
      USE nlist
      USE netcdf

	implicit none


!       include 'param.txt'
!       include 'common.txt'
      include 'mpif.h'
!	real x(1:imax,1:jmax,1:kmax)
!	real y(1:imax,1:jmax,1:kmax)
!	real z(1:imax,1:jmax,1:kmax)
	real uu(1:imax,1:jmax,1:kmax),obstacle(1:imax,1:jmax,1:kmax)
	real vv(1:imax,1:jmax,1:kmax),zzbed(1:imax,1:jmax)
	real mass_bed(1:nfrac,1:imax,1:jmax)
	real tt
	integer tel,istap,n
	character*60 FILE_NAME
	character*60 strng
	character*3 varname
	character*7 fname_basis
	logical(4) res
!	integer(4) ps

       ! We are writing 3D data, a nx x ny x nz grid.
       integer, parameter :: NDIMS = 3
       integer, parameter :: NDIMS2 = 4
       integer, parameter :: NDIMS3 = 1
	   integer, parameter :: NDIMS5 = 2
!       integer, parameter :: NX = imax, NY = jmax, NZ = kmax
     
       ! When we create netCDF files, variables and dimensions, we get back
       ! an ID for each one.
       integer :: ncid, varid1,varid2,varid3, varid4, varid5, varid6, varid7, varid8, varid9, varid10,varid11
	   integer :: dimids(NDIMS), dimids2(NDIMS2),dimids3(NDIMS3),dimids5(NDIMS5)
       integer :: x_dimid, y_dimid, z_dimid, nfrac_dimid, par_dimid
	integer :: dimids4(NDIMS),varid20,varid21,varid22,varid12,varid13,varid14,varid15,varid16,varid23
	character(1024) :: svnversion
	character(1024) :: svnurl
      include 'version.inc'

      do k=1,kmax
		do j=1,jmax
			do i=1,imax
				uu(i,j,k)=unew(i,j,k)*cos_u(j)-vnew(i,j,k)*sin_v(j) !unew(i,j,k)
				vv(i,j,k)=vnew(i,j,k)*cos_v(j)+unew(i,j,k)*sin_u(j) !vnew(i,j,k)
			enddo
		enddo
	enddo
	do i=1,imax
	   do j=1,jmax
	     	do n=1,nfrac
			mass_bed(n,i,j) = Cnewbot(n,i,j)*dz*frac(n)%rho ! Cnewbot(n,i,j)*dz*dr(i)*Rp(i)*dphi*frac(n)%rho
			enddo		
			IF (interaction_bed.ge.4) THEN
			  zzbed(i,j) = REAL(MAX(kbed(i,j)-1,0))*dz+ SUM(Clivebed(1:nfrac,i,j,kbed(i,j)))/cfixedbed*dz 
			ENDIF			
	   enddo
	enddo
	
	  if (applyVOF.eq.1) then 
		call state(cnew,rnew)
	  endif

	!WRITE(FILE_NAME,'(a,i9.9,a,i4.4,a)')'flow3D_',INT(istap),'_',INT(rank),'.nc'
	!WRITE(*,'(a,i9.9,a,i4.4,a)')'flow3D_',INT(istap),'_',INT(rank),'.nc'
	WRITE(FILE_NAME,'(a,i9.9,a,i4.4,a)'),fname_basis,INT(istap),'_',INT(rank),'.nc'
	WRITE(*,'(a,i9.9,a,i4.4,a)'),fname_basis,INT(istap),'_',INT(rank),'.nc'
	
       ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
       ! overwrite this file, if it already exists.
       call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )
	   
       ! Define the dimensions. NetCDF will hand back an ID for each.
       call check( nf90_def_dim(ncid, "xdim", imax, x_dimid) )
       call check( nf90_def_dim(ncid, "ydim", jmax, y_dimid) )
       call check( nf90_def_dim(ncid, "zdim", kmax, z_dimid) )
       if (nfrac>0) then
       call check( nf90_def_dim(ncid, "nfracdim", nfrac, nfrac_dimid) )
	endif
       call check( nf90_def_dim(ncid, "pardim", 1, par_dimid) )
     
       ! The dimids array is used to pass the IDs of the dimensions of
       ! the variables. Note that in fortran arrays are stored in
       ! column-major format.
       dimids =  (/ x_dimid, y_dimid, z_dimid /)
       if (nfrac>0) then
       	dimids2 =  (/ nfrac_dimid, x_dimid, y_dimid, z_dimid /)
	dimids4 =  (/ nfrac_dimid, x_dimid, y_dimid /)
	   dimids5 =  (/ x_dimid, y_dimid  /) 
		endif
       dimids3 =  (/ par_dimid  /) 


       ! Define the variable. The type of the variable in this case is
       ! NF90_DOUBLE (4-byte double).
       call check( nf90_def_var(ncid, "U", NF90_REAL, dimids, varid1) )
       call check( nf90_put_att(ncid, varid1, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid1, 'long_name', 'U velocity') )

       call check( nf90_def_var(ncid, "V", NF90_REAL, dimids, varid2) )
       call check( nf90_put_att(ncid, varid2, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid2, 'long_name', 'V velocity') )

       call check( nf90_def_var(ncid, "W", NF90_REAL, dimids, varid3) )
       call check( nf90_put_att(ncid, varid3, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid3, 'long_name', 'W velocity') )

       if (nfrac>0) then
       call check( nf90_def_var(ncid, "C", NF90_REAL, dimids2, varid4) )
       call check( nf90_put_att(ncid, varid4, 'units', '-') )
       call check( nf90_put_att(ncid, varid4, 'long_name', 'Volume concentration for each fraction') )
       call check( nf90_def_var(ncid, "mass_bed", NF90_REAL, dimids4, varid20) )
       call check( nf90_put_att(ncid, varid20, 'units', 'kg/m2') )
       call check( nf90_put_att(ncid, varid20, 'long_name', 'Mass per m2 sediment fractions in bed') )
		if (interaction_bed.ge.4) then
			call check( nf90_def_var(ncid, "zbed", NF90_REAL, dimids5, varid21) )
			call check( nf90_put_att(ncid, varid21, 'units', 'm') )
			call check( nf90_put_att(ncid, varid21, 'long_name', 'Bed level height of dynamic morphological bed') )
		   call check( nf90_def_var(ncid, "Cbed", NF90_REAL, dimids2, varid22) )
		   call check( nf90_put_att(ncid, varid22, 'units', '-') )
		   call check( nf90_put_att(ncid, varid22, 'long_name', 'Volume concentration for each fraction inside bed') )		
       call check( nf90_def_var(ncid, "kbed", NF90_SHORT, dimids5, varid11) )
       call check( nf90_put_att(ncid, varid11, 'units', '-') )
       call check( nf90_put_att(ncid, varid11, 'long_name', '2D index of highest bed cell') )

       call check( nf90_def_var(ncid, "Sbedload_u", NF90_REAL, dimids4, varid12) )
       call check( nf90_put_att(ncid, varid12, 'units', 'kg/m/s') )
       call check( nf90_put_att(ncid, varid12, 'long_name', 'U-bedload flux for each fraction (multiplied by morfac)') )	
       call check( nf90_def_var(ncid, "Sbedload_v", NF90_REAL, dimids4, varid13) )
       call check( nf90_put_att(ncid, varid13, 'units', 'kg/m/s') )
       call check( nf90_put_att(ncid, varid13, 'long_name', 'V-bedload flux for each fraction (multiplied by morfac)') )	   
		endif
	endif
       call check( nf90_def_var(ncid, "rho", NF90_REAL, dimids, varid5) )
       call check( nf90_put_att(ncid, varid5, 'units', 'kg/m3') )
       call check( nf90_put_att(ncid, varid5, 'long_name', 'Mixture density') )
     
       call check( nf90_def_var(ncid, "mu_t", NF90_REAL, dimids, varid6) )
       call check( nf90_put_att(ncid, varid6, 'units', 'kg/(sm)') )
       call check( nf90_put_att(ncid, varid6, 'long_name', 'Dynamic eddy viscosity') )

       call check( nf90_def_var(ncid, "P", NF90_REAL, dimids, varid7) )
       call check( nf90_put_att(ncid, varid7, 'units', 'Pa') )
       call check( nf90_put_att(ncid, varid7, 'long_name', 'Pressure') )
	   
	   if (sgs_model.eq.'DSmag') then
         call check( nf90_def_var(ncid, "Cs", NF90_REAL, dimids, varid9) )
         call check( nf90_put_att(ncid, varid9, 'units', '-') )
         call check( nf90_put_att(ncid, varid9, 'long_name', 'Smagorinsky constant from dynamic Germano-Lilly sgs model') )	   
       endif
	   if (sgs_model.eq.'ReaKE') then
         call check( nf90_def_var(ncid, "TKE", NF90_REAL, dimids, varid14) )
         call check( nf90_put_att(ncid, varid14, 'units', 'm2/s2') )
         call check( nf90_put_att(ncid, varid14, 'long_name', 'Simulated TKE from K-Epsilon model') )	  
         call check( nf90_def_var(ncid, "EPS", NF90_REAL, dimids, varid15) )
         call check( nf90_put_att(ncid, varid15, 'units', 'm2/s3') )
         call check( nf90_put_att(ncid, varid15, 'long_name', 'Simulated Epsilon (TKE dissipation) from K-Epsilon model') )		
         call check( nf90_def_var(ncid, "Cmu", NF90_REAL, dimids, varid16) )
         call check( nf90_put_att(ncid, varid16, 'units', '-') )
         call check( nf90_put_att(ncid, varid16, 'long_name', 'Simulated Cmu from realizible K-Epsilon model') )			 
       endif
	   
       call check( nf90_def_var(ncid, "wiggle_factor", NF90_REAL, dimids, varid10) )
       call check( nf90_put_att(ncid, varid10, 'units', '-') )
       call check( nf90_put_att(ncid, varid10, 'long_name', 'wiggle factor blend: 0=no wiggles, 1=wiggles') )
	   
       call check( nf90_def_var(ncid, "time", NF90_REAL, dimids3, varid8) )
       call check( nf90_put_att(ncid, varid8, 'units', 's') )
       call check( nf90_put_att(ncid, varid8, 'long_name', 'Time from start simulation') )

	   obstacle=1.-fc_global(1:imax,1+rank*jmax:jmax+rank*jmax,1:kmax)
       call check( nf90_def_var(ncid, "obstacle", NF90_REAL, dimids, varid23) )
       call check( nf90_put_att(ncid, varid23, 'units', '-') )
       call check( nf90_put_att(ncid, varid23, 'long_name', 'If cell is in obstacle 1 else 0') )
	   
	! also add svn info in output files:
       CALL check( nf90_put_att(ncid,nf90_global, "svnversion", trim(svnversion)))
       CALL check( nf90_put_att(ncid,nf90_global, "svnurl", trim(svnurl)))

       ! End define mode. This tells netCDF we are done defining metadata.
       call check( nf90_enddef(ncid) )
   
       ! Write the pretend data to the file. Although netCDF supports
       ! reading and writing subsets of data, in this case we write all the
       ! data in one operation.
       call check( nf90_put_var(ncid, varid1, uu(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid2, vv(1:imax,1:jmax,1:kmax)) )
	   !call check( nf90_put_var(ncid, varid3, Wnew(1:imax,1:jmax,1:kmax)) )
	   ! Strange error for very large sims: line above goes wrong, 2 lines below (which basically do the same) goes OK....
	   uu=Wnew(1:imax,1:jmax,1:kmax)
	   call check( nf90_put_var(ncid, varid3, uu(1:imax,1:jmax,1:kmax)) )
       

	   
       if (nfrac>0) then
       	call check( nf90_put_var(ncid, varid4, Cnew(1:nfrac,1:imax,1:jmax,1:kmax)) ) 
		call check( nf90_put_var(ncid, varid20, mass_bed(1:nfrac,1:imax,1:jmax)) ) 
		if (interaction_bed.ge.4) then
		  call check( nf90_put_var(ncid, varid21, zzbed(1:imax,1:jmax) ))
		  call check( nf90_put_var(ncid, varid22, Clivebed(1:nfrac,1:imax,1:jmax,1:kmax)) )
		  call check( nf90_put_var(ncid, varid11, kbed(1:imax,1:jmax) ))
		  call check( nf90_put_var(ncid, varid12, qbU(1:nfrac,1:imax,1:jmax) ))
		  call check( nf90_put_var(ncid, varid13, qbV(1:nfrac,1:imax,1:jmax) ))
		endif
		  
	endif
		uu=rnew(1:imax,1:jmax,1:kmax)
		call check( nf90_put_var(ncid, varid5, uu(1:imax,1:jmax,1:kmax)) )
		uu=ekm(1:imax,1:jmax,1:kmax)
		call check( nf90_put_var(ncid, varid6, uu(1:imax,1:jmax,1:kmax)) )
		uu=p(1:imax,1:jmax,1:kmax)+pold(1:imax,1:jmax,1:kmax)
		call check( nf90_put_var(ncid, varid7, uu(1:imax,1:jmax,1:kmax)) )
       !call check( nf90_put_var(ncid, varid5, rnew(1:imax,1:jmax,1:kmax)) )
       !call check( nf90_put_var(ncid, varid6, ekm(1:imax,1:jmax,1:kmax)) )
       !call check( nf90_put_var(ncid, varid7, p(1:imax,1:jmax,1:kmax)+pold(1:imax,1:jmax,1:kmax)) )
	   !!!call check( nf90_put_var(ncid, varid7, p(1:imax,1:jmax,1:kmax)) ) !! changed for exact solver output,
	   if (sgs_model.eq.'DSmag') then
	     call check( nf90_put_var(ncid, varid9, Csgrid(1:imax,1:jmax,1:kmax)) )
	   endif
	   if (sgs_model.eq.'ReaKE') then
	     call check( nf90_put_var(ncid, varid14, TKE(1:imax,1:jmax,1:kmax)) )
		 call check( nf90_put_var(ncid, varid15, EEE(1:imax,1:jmax,1:kmax)) )
		 call check( nf90_put_var(ncid, varid16, Cmu(1:imax,1:jmax,1:kmax)) )
	   endif	   
	   call check( nf90_put_var(ncid, varid10, wf(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid8, tt) )
	   call check( nf90_put_var(ncid, varid23, obstacle) )
    
       ! Close the file. This frees up any internal netCDF resources
       ! associated with the file, and flushes any buffers.
       call check( nf90_close(ncid) )
	   
	   
	   if (applyVOF.eq.1) then 
		 rnew=rho_b
	   endif

	end


      subroutine output_nc_movie(istap,tt)
      USE nlist
      USE netcdf

	implicit none


!       include 'param.txt'
!       include 'common.txt'
      include 'mpif.h'
	real tt,add_offset,scale_factor,data_range
	integer tel,istap
	character*60 FILE_NAME
	character*60 strng
	character*3 varname
	logical(4) res
	real uu(1:imax,1:jmax,1:kmax)
	real vv(1:imax,1:jmax,1:kmax)
!	integer(4) ps

       ! We are writing 3D data, a nx x ny x nz grid.
       integer, parameter :: NDIMS = 3
       integer, parameter :: NDIMS2 = 4
       integer, parameter :: NDIMS3 = 1
!       integer, parameter :: NX = imax, NY = jmax, NZ = kmax
     
       ! When we create netCDF files, variables and dimensions, we get back
       ! an ID for each one.
       integer :: ncid, varid9,varid10,varid4,varid8, dimids(NDIMS), dimids2(NDIMS2),dimids3(NDIMS3)
       integer :: varid11,varid12,varid13,varid14,varid15,varid16,varid17,varid18,varid19
       integer :: x_dimid, y_dimid, z_dimid, nfrac_dimid, par_dimid
	character(1024) :: svnversion
	character(1024) :: svnurl
      include 'version.inc'       

	WRITE(FILE_NAME,'(a,i9.9,a,i4.4,a)')'movie3D_',INT(istap),'_',INT(rank),'.nc'
	!WRITE(*,'(a,i9.9,a,i4.4,a)')'movie3D_',INT(istap),'_',INT(rank),'.nc'
	
       ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
       ! overwrite this file, if it already exists.
       call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )
     
       ! Define the dimensions. NetCDF will hand back an ID for each.
       call check( nf90_def_dim(ncid, "xdim", imax, x_dimid) )
       call check( nf90_def_dim(ncid, "ydim", jmax, y_dimid) )
       call check( nf90_def_dim(ncid, "zdim", kmax, z_dimid) )
       if (nfrac>0) then
       call check( nf90_def_dim(ncid, "nfracdim", nfrac, nfrac_dimid) )
	endif
       call check( nf90_def_dim(ncid, "pardim", 1, par_dimid) )
     
       ! The dimids array is used to pass the IDs of the dimensions of
       ! the variables. Note that in fortran arrays are stored in
       ! column-major format.
       dimids =  (/ x_dimid, y_dimid, z_dimid /)
       if (nfrac>0) then
       dimids2 =  (/ nfrac_dimid, x_dimid, y_dimid, z_dimid /)
	endif
       dimids3 =  (/ par_dimid  /) 
    
       ! Define the variable. The type of the variable in this case is
       ! NF90_SHORT (2-byte 16bit var).
       if (nfrac>0) then
       call check( nf90_def_var(ncid, "C", NF90_SHORT, dimids2, varid4) )
       call check( nf90_put_att(ncid, varid4, 'units', '-') )
       call check( nf90_put_att(ncid, varid4, 'long_name', 'Volume concentration for each fraction') )
       call check( nf90_def_var(ncid, "scale_factor_c", NF90_REAL, dimids3, varid9) )
       call check( nf90_put_att(ncid, varid9, 'units', '-') )
       call check( nf90_put_att(ncid, varid9, 'long_name', 'unpacked_value = packed_value * scale_factor + add_offset') )
       call check( nf90_def_var(ncid, "add_offset_c", NF90_REAL, dimids3, varid10) )
       call check( nf90_put_att(ncid, varid10, 'units', '-') )
       call check( nf90_put_att(ncid, varid10, 'long_name', 'unpacked_value = packed_value * scale_factor + add_offset') )
	endif

       call check( nf90_def_var(ncid, "U", NF90_SHORT, dimids, varid11) )
       call check( nf90_put_att(ncid, varid11, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid11, 'long_name', 'U velocity') )
       call check( nf90_def_var(ncid, "scale_factor_u", NF90_REAL, dimids3, varid12) )
       call check( nf90_put_att(ncid, varid12, 'units', '-') )
       call check( nf90_put_att(ncid, varid12, 'long_name', 'unpacked_value = packed_value * scale_factor + add_offset') )
       call check( nf90_def_var(ncid, "add_offset_u", NF90_REAL, dimids3, varid13) )
       call check( nf90_put_att(ncid, varid13, 'units', '-') )
       call check( nf90_put_att(ncid, varid13, 'long_name', 'unpacked_value = packed_value * scale_factor + add_offset') )
       call check( nf90_def_var(ncid, "V", NF90_SHORT, dimids, varid14) )
       call check( nf90_put_att(ncid, varid14, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid14, 'long_name', 'V velocity') )
       call check( nf90_def_var(ncid, "scale_factor_v", NF90_REAL, dimids3, varid15) )
       call check( nf90_put_att(ncid, varid15, 'units', '-') )
       call check( nf90_put_att(ncid, varid15, 'long_name', 'unpacked_value = packed_value * scale_factor + add_offset') )
       call check( nf90_def_var(ncid, "add_offset_v", NF90_REAL, dimids3, varid16) )
       call check( nf90_put_att(ncid, varid16, 'units', '-') )
       call check( nf90_put_att(ncid, varid16, 'long_name', 'unpacked_value = packed_value * scale_factor + add_offset') )
       call check( nf90_def_var(ncid, "W", NF90_SHORT, dimids, varid17) )
       call check( nf90_put_att(ncid, varid17, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid17, 'long_name', 'W velocity') )
       call check( nf90_def_var(ncid, "scale_factor_w", NF90_REAL, dimids3, varid18) )
       call check( nf90_put_att(ncid, varid18, 'units', '-') )
       call check( nf90_put_att(ncid, varid18, 'long_name', 'unpacked_value = packed_value * scale_factor + add_offset') )
       call check( nf90_def_var(ncid, "add_offset_w", NF90_REAL, dimids3, varid19) )
       call check( nf90_put_att(ncid, varid19, 'units', '-') )
       call check( nf90_put_att(ncid, varid19, 'long_name', 'unpacked_value = packed_value * scale_factor + add_offset') )

       call check( nf90_def_var(ncid, "time", NF90_REAL, dimids3, varid8) )
       call check( nf90_put_att(ncid, varid8, 'units', 's') )
       call check( nf90_put_att(ncid, varid8, 'long_name', 'Time from start simulation') )

	! also add svn info in output files:
       CALL check( nf90_put_att(ncid,nf90_global, "svnversion", trim(svnversion)))
       CALL check( nf90_put_att(ncid,nf90_global, "svnurl", trim(svnurl)))

       ! End define mode. This tells netCDF we are done defining metadata.
       call check( nf90_enddef(ncid) )
     
       ! Write the pretend data to the file. Although netCDF supports
       ! reading and writing subsets of data, in this case we write all the
       ! data in one operation.

       if (nfrac>0) then
	add_offset = MINVAL(Cnew(1:nfrac,1:imax,1:jmax,1:kmax))
	data_range = MAXVAL(Cnew(1:nfrac,1:imax,1:jmax,1:kmax))-MINVAL(Cnew(1:nfrac,1:imax,1:jmax,1:kmax))
	if (data_range > 0.) then
	  scale_factor = data_range/(2.**15-1.)
          call check( nf90_put_var(ncid, varid4, nint((Cnew(1:nfrac,1:imax,1:jmax,1:kmax)-add_offset)/scale_factor) ))
          call check( nf90_put_var(ncid, varid9, scale_factor) )
          call check( nf90_put_var(ncid, varid10, add_offset) )
	else
	  scale_factor = data_range/(2.**15-1.) ! scale_factor is zero, prevent division by zero:
          call check( nf90_put_var(ncid, varid4, nint(0.*Cnew(1:nfrac,1:imax,1:jmax,1:kmax)) ))
          call check( nf90_put_var(ncid, varid9, scale_factor) )
          call check( nf90_put_var(ncid, varid10, add_offset) )
        endif
       endif

        do k=1,kmax
		do j=1,jmax
			do i=1,imax
				uu(i,j,k)=unew(i,j,k)*cos_u(j)-vnew(i,j,k)*sin_v(j) !unew(i,j,k)
				vv(i,j,k)=vnew(i,j,k)*cos_v(j)+unew(i,j,k)*sin_u(j) !vnew(i,j,k)
			enddo
		enddo
	enddo

	add_offset = MINVAL(uu(1:imax,1:jmax,1:kmax))
	data_range = MAXVAL(uu(1:imax,1:jmax,1:kmax))-MINVAL(uu(1:imax,1:jmax,1:kmax))
	if (data_range > 0.) then
	  scale_factor = data_range/(2.**15-1.)
          call check( nf90_put_var(ncid, varid11, nint((uu(1:imax,1:jmax,1:kmax)-add_offset)/scale_factor) ))
          call check( nf90_put_var(ncid, varid12, scale_factor) )
          call check( nf90_put_var(ncid, varid13, add_offset) )
	else
	  scale_factor = data_range/(2.**15-1.) ! scale_factor is zero, prevent division by zero:
          call check( nf90_put_var(ncid, varid11, nint(0.*uu(1:imax,1:jmax,1:kmax)) ))
          call check( nf90_put_var(ncid, varid12, scale_factor) )
          call check( nf90_put_var(ncid, varid13, add_offset) )
        endif

	add_offset = MINVAL(vv(1:imax,1:jmax,1:kmax))
	data_range = MAXVAL(vv(1:imax,1:jmax,1:kmax))-MINVAL(vv(1:imax,1:jmax,1:kmax))
	if (data_range > 0.) then
	  scale_factor = data_range/(2.**15-1.)
          call check( nf90_put_var(ncid, varid14, nint((vv(1:imax,1:jmax,1:kmax)-add_offset)/scale_factor) ))
          call check( nf90_put_var(ncid, varid15, scale_factor) )
          call check( nf90_put_var(ncid, varid16, add_offset) )
	else
	  scale_factor = data_range/(2.**15-1.) ! scale_factor is zero, prevent division by zero:
          call check( nf90_put_var(ncid, varid14, nint(0.*vv(1:imax,1:jmax,1:kmax)) ))
          call check( nf90_put_var(ncid, varid15, scale_factor) )
          call check( nf90_put_var(ncid, varid16, add_offset) )
        endif

	add_offset = MINVAL(wnew(1:imax,1:jmax,1:kmax))
	data_range = MAXVAL(wnew(1:imax,1:jmax,1:kmax))-MINVAL(wnew(1:imax,1:jmax,1:kmax))
	if (data_range > 0.) then
	  scale_factor = data_range/(2.**15-1.)
          call check( nf90_put_var(ncid, varid17, nint((wnew(1:imax,1:jmax,1:kmax)-add_offset)/scale_factor) ))
          call check( nf90_put_var(ncid, varid18, scale_factor) )
          call check( nf90_put_var(ncid, varid19, add_offset) )
	else
	  scale_factor = data_range/(2.**15-1.) ! scale_factor is zero, prevent division by zero:
          call check( nf90_put_var(ncid, varid17, nint(0.*wnew(1:imax,1:jmax,1:kmax)) ))
          call check( nf90_put_var(ncid, varid18, scale_factor) )
          call check( nf90_put_var(ncid, varid19, add_offset) )
        endif

       call check( nf90_put_var(ncid, varid8, tt) )
     
       ! Close the file. This frees up any internal netCDF resources
       ! associated with the file, and flushes any buffers.
       call check( nf90_close(ncid) )

	end


      subroutine output_stat_nc(tt)
      USE nlist
      USE netcdf
!      USE work_array


	implicit none


!       include 'param.txt'
!       include 'common.txt'
      include 'mpif.h'

	real Urms(1:imax,1:jmax,1:kmax)
	real Vrms(1:imax,1:jmax,1:kmax)
	real Wrms(1:imax,1:jmax,1:kmax)
	real Rrms(1:imax,1:jmax,1:kmax)
	real Crms(nfrac,1:imax,1:jmax,1:kmax)
	real uv_shear(1:imax,1:jmax,1:kmax)
	real vw_shear(1:imax,1:jmax,1:kmax)
	real uw_shear(1:imax,1:jmax,1:kmax)
	real uc_avg(nfrac,1:imax,1:jmax,1:kmax),vc_avg(nfrac,1:imax,1:jmax,1:kmax),wc_avg(nfrac,1:imax,1:jmax,1:kmax)
	real uc_rms(nfrac,1:imax,1:jmax,1:kmax),vc_rms(nfrac,1:imax,1:jmax,1:kmax),wc_rms(nfrac,1:imax,1:jmax,1:kmax)
	real tt
	integer tel,istap,n
	character*60 FILE_NAME
	character*60 strng
	character*3 varname
	logical(4) res
!	integer(4) ps

       ! We are writing 3D data, a nx x ny x nz grid.
       integer, parameter :: NDIMS = 3
       integer, parameter :: NDIMS2 = 4
       integer, parameter :: NDIMS3 = 1
!       integer, parameter :: NX = imax, NY = jmax, NZ = kmax
     
       ! When we create netCDF files, variables and dimensions, we get back
       ! an ID for each one.
       integer :: ncid, varid1,varid2,varid3, varid4, varid5, varid6, varid7, varid8 
       integer :: varid9,varid10,varid11, varid12, varid13, varid14, varid15, varid16, varid17
       integer :: varid18,varid19,varid20,varid21,varid22,varid23,varid24,varid25
	   integer :: varid26,varid27,varid28,varid29,varid30,varid31,varid32,varid33
       integer :: dimids(NDIMS), dimids2(NDIMS2),dimids3(NDIMS3)
       integer :: x_dimid, y_dimid, z_dimid, nfrac_dimid, par_dimid
	character(1024) :: svnversion
	character(1024) :: svnurl
      include 'version.inc'

	do i=1,imax
	  do j=1,jmax
	    do k=1,kmax
		  Uavg(i,j,k)  = Uavg(i,j,k)/stat_time_count
		  Vavg(i,j,k)  = Vavg(i,j,k)/stat_time_count 
		  Wavg(i,j,k)  = Wavg(i,j,k)/stat_time_count 
		  Ravg(i,j,k)  = Ravg(i,j,k)/stat_time_count 
		  Pavg(i,j,k)  = Pavg(i,j,k)/stat_time_count
		  muavg(i,j,k)  = muavg(i,j,k)/stat_time_count

		  Urms(i,j,k)  = sqrt(ABS(sigU2(i,j,k)/stat_time_count-Uavg(i,j,k)*Uavg(i,j,k) ))
		  Vrms(i,j,k)  = sqrt(ABS(sigV2(i,j,k)/stat_time_count-Vavg(i,j,k)*Vavg(i,j,k) ))
		  Wrms(i,j,k)  = sqrt(ABS(sigW2(i,j,k)/stat_time_count-Wavg(i,j,k)*Wavg(i,j,k) ))
		  Rrms(i,j,k)  = sqrt(ABS(sigR2(i,j,k)/stat_time_count-Ravg(i,j,k)*Ravg(i,j,k) ))

		do n=1,nfrac
		  Cavg(n,i,j,k)  = Cavg(n,i,j,k)/stat_time_count 
		  Crms(n,i,j,k)  = sqrt(ABS(sigC2(n,i,j,k)/stat_time_count-Cavg(n,i,j,k)*Cavg(n,i,j,k) ))

		  uc_avg(n,i,j,k) = sigUC(n,i,j,k)/stat_time_count
		  vc_avg(n,i,j,k) = sigVC(n,i,j,k)/stat_time_count
		  wc_avg(n,i,j,k) = sigWC(n,i,j,k)/stat_time_count
		  uc_rms(n,i,j,k) = sigUC(n,i,j,k)/stat_time_count-Uavg(i,j,k)*Cavg(n,i,j,k)
		  vc_rms(n,i,j,k) = sigVC(n,i,j,k)/stat_time_count-Vavg(i,j,k)*Cavg(n,i,j,k)
		  wc_rms(n,i,j,k) = sigWC(n,i,j,k)/stat_time_count-Wavg(i,j,k)*Cavg(n,i,j,k)
		enddo

		  uv_shear(i,j,k) = sigUV(i,j,k)/stat_time_count-Uavg(i,j,k)*Vavg(i,j,k) 
		  uw_shear(i,j,k) = sigUW(i,j,k)/stat_time_count-Uavg(i,j,k)*Wavg(i,j,k)
		  vw_shear(i,j,k) = sigVW(i,j,k)/stat_time_count-Wavg(i,j,k)*Vavg(i,j,k)



!	   	  fUavg(i,j,k) = fUavg(i,j,k)/Ravg(i,j,k)/stat_count
!		  fUrms(i,j,k) = sqrt(ABS(sigfU2(i,j,k)/(stat_count*Ravg(i,j,k))-fUavg(i,j,k)**2 ))
	    enddo
	  enddo
	enddo

	WRITE(FILE_NAME,'(a,i4.4,a)')'Stat_output_',INT(rank),'.nc'
	
       ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
       ! overwrite this file, if it already exists.
       call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )
     
       ! Define the dimensions. NetCDF will hand back an ID for each.
       call check( nf90_def_dim(ncid, "xdim", imax, x_dimid) )
       call check( nf90_def_dim(ncid, "ydim", jmax, y_dimid) )
       call check( nf90_def_dim(ncid, "zdim", kmax, z_dimid) )
       if (nfrac>0) then
       call check( nf90_def_dim(ncid, "nfracdim", nfrac, nfrac_dimid) )
	endif
       call check( nf90_def_dim(ncid, "pardim", 1, par_dimid) )
     
       ! The dimids array is used to pass the IDs of the dimensions of
       ! the variables. Note that in fortran arrays are stored in
       ! column-major format.
       dimids =  (/ x_dimid, y_dimid, z_dimid /)
       if (nfrac>0) then
       dimids2 =  (/ nfrac_dimid, x_dimid, y_dimid, z_dimid /)
	endif
       dimids3 =  (/ par_dimid  /) 
    
       ! Define the variable. The type of the variable in this case is
       ! NF90_DOUBLE (4-byte double).
       call check( nf90_def_var(ncid, "Urms", NF90_REAL, dimids, varid1) )
       call check( nf90_put_att(ncid, varid1, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid1, 'long_name', 'Urms velocity') )

       call check( nf90_def_var(ncid, "Vrms", NF90_REAL, dimids, varid2) )
       call check( nf90_put_att(ncid, varid2, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid2, 'long_name', 'Vrms velocity') )

       call check( nf90_def_var(ncid, "Wrms", NF90_REAL, dimids, varid3) )
       call check( nf90_put_att(ncid, varid3, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid3, 'long_name', 'Wrms velocity') )

       call check( nf90_def_var(ncid, "Uavg", NF90_REAL, dimids, varid4) )
       call check( nf90_put_att(ncid, varid4, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid4, 'long_name', 'Uavg velocity') )

       call check( nf90_def_var(ncid, "Vavg", NF90_REAL, dimids, varid5) )
       call check( nf90_put_att(ncid, varid5, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid5, 'long_name', 'Vavg velocity') )

       call check( nf90_def_var(ncid, "Wavg", NF90_REAL, dimids, varid6) )
       call check( nf90_put_att(ncid, varid6, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid6, 'long_name', 'Wavg velocity') )
	   
       call check( nf90_def_var(ncid, "Umax", NF90_REAL, dimids, varid26) )
       call check( nf90_put_att(ncid, varid26, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid26, 'long_name', 'Umax velocity') )

       call check( nf90_def_var(ncid, "Vmax", NF90_REAL, dimids, varid27) )
       call check( nf90_put_att(ncid, varid27, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid27, 'long_name', 'Vmax velocity') )

       call check( nf90_def_var(ncid, "Wmax", NF90_REAL, dimids, varid28) )
       call check( nf90_put_att(ncid, varid28, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid28, 'long_name', 'Wmax velocity') )	
	   
!       call check( nf90_def_var(ncid, "Umin", NF90_REAL, dimids, varid29) )
!       call check( nf90_put_att(ncid, varid29, 'units', 'm/s') )
!       call check( nf90_put_att(ncid, varid29, 'long_name', 'Umin velocity') )
!
!       call check( nf90_def_var(ncid, "Vmin", NF90_REAL, dimids, varid30) )
!       call check( nf90_put_att(ncid, varid30, 'units', 'm/s') )
!       call check( nf90_put_att(ncid, varid30, 'long_name', 'Vmin velocity') )
!
!       call check( nf90_def_var(ncid, "Wmin", NF90_REAL, dimids, varid31) )
!       call check( nf90_put_att(ncid, varid31, 'units', 'm/s') )
!       call check( nf90_put_att(ncid, varid31, 'long_name', 'Wmin velocity') )		

       call check( nf90_def_var(ncid, "Uhormax", NF90_REAL, dimids, varid32) )
       call check( nf90_put_att(ncid, varid32, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid32, 'long_name', 'Uhor max velocity MAX of sqrt(uu^2+vv^2)') )
	   call check( nf90_def_var(ncid, "U3dmax", NF90_REAL, dimids, varid33) )
       call check( nf90_put_att(ncid, varid33, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid33, 'long_name', 'U3d max velocity MAX of sqrt(uu^2+vv^2+ww^2)') )
	   
       call check( nf90_def_var(ncid, "rho_rms", NF90_REAL, dimids, varid7) )
       call check( nf90_put_att(ncid, varid7, 'units', 'kg/m3') )
       call check( nf90_put_att(ncid, varid7, 'long_name', 'RMS Mixture density') )

       call check( nf90_def_var(ncid, "rho_avg", NF90_REAL, dimids, varid8) )
       call check( nf90_put_att(ncid, varid8, 'units', 'kg/m3') )
       call check( nf90_put_att(ncid, varid8, 'long_name', 'AVG Mixture density') )

       call check( nf90_def_var(ncid, "UV_shear", NF90_REAL, dimids, varid9) )
       call check( nf90_put_att(ncid, varid9, 'units', 'm2/s2') )
       call check( nf90_put_att(ncid, varid9, 'long_name', 'UV shear stress') )

       call check( nf90_def_var(ncid, "VW_shear", NF90_REAL, dimids, varid10) )
       call check( nf90_put_att(ncid, varid10, 'units', 'm2/s2') )
       call check( nf90_put_att(ncid, varid10, 'long_name', 'VW shear stress') )

       call check( nf90_def_var(ncid, "UW_shear", NF90_REAL, dimids, varid11) )
       call check( nf90_put_att(ncid, varid11, 'units', 'm2/s2') )
       call check( nf90_put_att(ncid, varid11, 'long_name', 'UW shear stress') )

       call check( nf90_def_var(ncid, "Pavg", NF90_REAL, dimids, varid12) )
       call check( nf90_put_att(ncid, varid12, 'units', 'Pa') )
       call check( nf90_put_att(ncid, varid12, 'long_name', 'AVG Pressure') )

       call check( nf90_def_var(ncid, "mu_avg", NF90_REAL, dimids, varid13) )
       call check( nf90_put_att(ncid, varid13, 'units', 'kg/(sm)') )
       call check( nf90_put_att(ncid, varid13, 'long_name', 'AVG Dynamic eddy viscosity') )

       if (nfrac>0) then
       call check( nf90_def_var(ncid, "Cavg", NF90_REAL, dimids2, varid14) )
       call check( nf90_put_att(ncid, varid14, 'units', '-') )
       call check( nf90_put_att(ncid, varid14, 'long_name', 'AVG Volume concentration for each fraction') )

       call check( nf90_def_var(ncid, "Crms", NF90_REAL, dimids2, varid15) )
       call check( nf90_put_att(ncid, varid15, 'units', '-') )
       call check( nf90_put_att(ncid, varid15, 'long_name', 'RMS Volume concentration for each fraction') )

       call check( nf90_def_var(ncid, "Cmax", NF90_REAL, dimids2, varid24) )
       call check( nf90_put_att(ncid, varid24, 'units', '-') )
       call check( nf90_put_att(ncid, varid24, 'long_name', 'MAX Volume concentration for each fraction') )

       call check( nf90_def_var(ncid, "Cmin", NF90_REAL, dimids2, varid25) )
       call check( nf90_put_att(ncid, varid25, 'units', '-') )
       call check( nf90_put_att(ncid, varid25, 'long_name', 'MIN Volume concentration for each fraction') )

       call check( nf90_def_var(ncid, "UC_rms", NF90_REAL, dimids2, varid18) )
       call check( nf90_put_att(ncid, varid18, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid18, 'long_name', 'Turbulent volume conc. flux UCrms u_accent*c_accent') )

       call check( nf90_def_var(ncid, "VC_rms", NF90_REAL, dimids2, varid19) )
       call check( nf90_put_att(ncid, varid19, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid19, 'long_name', 'Turbulent volume conc. flux VCrms v_accent*c_accent') )

       call check( nf90_def_var(ncid, "WC_rms", NF90_REAL, dimids2, varid20) )
       call check( nf90_put_att(ncid, varid20, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid20, 'long_name', 'Turbulent volume conc. flux WCrms w_accent*c_accent') )

       call check( nf90_def_var(ncid, "UC_avg", NF90_REAL, dimids2, varid21) )
       call check( nf90_put_att(ncid, varid21, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid21, 'long_name', 'Avg. volume conc. flux UCavg avg_of(u*c)') )

       call check( nf90_def_var(ncid, "VC_avg", NF90_REAL, dimids2, varid22) )
       call check( nf90_put_att(ncid, varid22, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid22, 'long_name', 'Avg. volume conc. flux VCavg avg_of(v*c)') )

       call check( nf90_def_var(ncid, "WC_avg", NF90_REAL, dimids2, varid23) )
       call check( nf90_put_att(ncid, varid23, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid23, 'long_name', 'Avg. volume conc. flux WCavg avg_of(w*c)') )

	endif     
       call check( nf90_def_var(ncid, "time", NF90_REAL, dimids3, varid16) )
       call check( nf90_put_att(ncid, varid16, 'units', 's') )
       call check( nf90_put_att(ncid, varid16, 'long_name', 'Time from start simulation') )

       call check( nf90_def_var(ncid, "avg_time", NF90_REAL, dimids3, varid17) )
       call check( nf90_put_att(ncid, varid17, 'units', 's') )
       call check( nf90_put_att(ncid, varid17, 'long_name', 'Time over which average is determined') )

	! also add svn info in output files:
       CALL check( nf90_put_att(ncid,nf90_global, "svnversion", trim(svnversion)))
       CALL check( nf90_put_att(ncid,nf90_global, "svnurl", trim(svnurl)))

       ! End define mode. This tells netCDF we are done defining metadata.
       call check( nf90_enddef(ncid) )
     
       ! Write the pretend data to the file. Although netCDF supports
       ! reading and writing subsets of data, in this case we write all the
       ! data in one operation.
       call check( nf90_put_var(ncid, varid1, Urms(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid2, Vrms(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid3, Wrms(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid4, Uavg(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid5, Vavg(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid6, Wavg(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid26, Umax(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid27, Vmax(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid28, Wmax(1:imax,1:jmax,1:kmax)) )	   
	   call check( nf90_put_var(ncid, varid32, Uhormax(1:imax,1:jmax,1:kmax)) )
	   call check( nf90_put_var(ncid, varid33, U3dmax(1:imax,1:jmax,1:kmax)) )
       !call check( nf90_put_var(ncid, varid29, Umin(1:imax,1:jmax,1:kmax)) )
       !call check( nf90_put_var(ncid, varid30, Vmin(1:imax,1:jmax,1:kmax)) )
       !call check( nf90_put_var(ncid, varid31, Wmin(1:imax,1:jmax,1:kmax)) )	   	   
       call check( nf90_put_var(ncid, varid7, Rrms(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid8, Ravg(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid9, uv_shear(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid10, vw_shear(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid11, uw_shear(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid12, Pavg(1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid13, muavg(1:imax,1:jmax,1:kmax)) )

       if (nfrac>0) then
       call check( nf90_put_var(ncid, varid14, Cavg(1:nfrac,1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid15, Crms(1:nfrac,1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid24, Cmax(1:nfrac,1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid25, Cmin(1:nfrac,1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid18, uc_rms(1:nfrac,1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid19, vc_rms(1:nfrac,1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid20, wc_rms(1:nfrac,1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid21, uc_avg(1:nfrac,1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid22, vc_avg(1:nfrac,1:imax,1:jmax,1:kmax)) )
       call check( nf90_put_var(ncid, varid23, wc_avg(1:nfrac,1:imax,1:jmax,1:kmax)) )
	endif

       call check( nf90_put_var(ncid, varid16, tt) )
       call check( nf90_put_var(ncid, varid17, stat_time_count) )


!	CALL wvar2matf(Urms(1:imax,1:jmax,1:kmax),fmatname,'Urms')
!	CALL wvar2matf(Vrms(1:imax,1:jmax,1:kmax),fmatname,'Vrms')
!	CALL wvar2matf(Wrms(1:imax,1:jmax,1:kmax),fmatname,'Wrms')
!	CALL wvar2matf(Rrms(1:imax,1:jmax,1:kmax),fmatname,'rho_rms')
!	CALL wvar2matf(Uavg(1:imax,1:jmax,1:kmax),fmatname,'Uavg')
!	CALL wvar2matf(Vavg(1:imax,1:jmax,1:kmax),fmatname,'Vavg')
!	CALL wvar2matf(Wavg(1:imax,1:jmax,1:kmax),fmatname,'Wavg')
!	CALL wvar2matf(Ravg(1:imax,1:jmax,1:kmax),fmatname,'rho_avg')
!	CALL wvar2matf(uv_shear(1:imax,1:jmax,1:kmax),fmatname,'uv_shear')
!	CALL wvar2matf(vw_shear(1:imax,1:jmax,1:kmax),fmatname,'vw_shear')
!	CALL wvar2matf(uw_shear(1:imax,1:jmax,1:kmax),fmatname,'uw_shear')
!        CALL wvar2matf(Pavg(1:imax,1:jmax,1:kmax),fmatname,'p_avg')
!        CALL wvar2matf(muavg(1:imax,1:jmax,1:kmax),fmatname,'mu_avg')
!	do n=1,nfrac
!	  	WRITE(varname,'(a,i2.2)')'Cavg',INT(n)
!		CALL wvar2matf(Cavg(n,1:imax,1:jmax,1:kmax),fmatname,varname)
!	  	WRITE(varname,'(a,i2.2)')'Crms',INT(n)
!		CALL wvar2matf(Crms(n,1:imax,1:jmax,1:kmax),fmatname,varname)
!	enddo
!	CALL wvar2matf(tt,fmatname,'time')
!	CALL wvar2matf(stat_count,fmatname,'stat_count')
     
       ! Close the file. This frees up any internal netCDF resources
       ! associated with the file, and flushes any buffers.
       call check( nf90_close(ncid) )

	end

	
		subroutine output_his_bedplume
		USE nlist
        USE netcdf
	implicit none

!       include 'param.txt'
      include 'mpif.h'
!	real UU(1:200000)
	integer ierr,n,r,tag,status(MPI_STATUS_SIZE),nf
	INTEGER (8) :: ps

       ! We are writing 3D, 2D data data
       integer, parameter :: NDIMS2 = 2
       integer, parameter :: NDIMS3 = 3
     
       ! When we create netCDF files, variables and dimensions, we get back
       ! an ID for each one.
       integer :: ncid, varid1,varid2,varid3,varid4,varid5,varid6
       integer :: dimids2(NDIMS2),dimids3(NDIMS3)
       integer :: nhis_dimid,time_dimid,nfrac_dimid,istep2
	character(1024) :: svnversion
	character(1024) :: svnurl
      include 'version.inc'

       ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
       ! overwrite this file, if it already exists.
       call check( nf90_create('history_bedplume.nc', NF90_CLOBBER, ncid) )
     
       ! Define the dimensions. NetCDF will hand back an ID for each.
       call check( nf90_def_dim(ncid, "nhis_bp_dim", nbedplume, nhis_dimid) )
       call check( nf90_def_dim(ncid, "time_dim", 20000, time_dimid) )
		if (nfrac>0) then
		   call check( nf90_def_dim(ncid, "nfrac_dim", nfrac, nfrac_dimid) )
		endif
     
       ! The dimids array is used to pass the IDs of the dimensions of
       ! the variables. Note that in fortran arrays are stored in
       ! column-major format.
       dimids2 =  (/ nhis_dimid, time_dimid/)
	if (nfrac>0) then
       dimids3 =  (/ nfrac_dimid, nhis_dimid, time_dimid/)
	endif
    
       ! Define the variable. The type of the variable in this case is
       ! NF90_DOUBLE (4-byte double).
       call check( nf90_def_var(ncid, "this_bp", NF90_DOUBLE, dimids2, varid1) )
       call check( nf90_put_att(ncid, varid1, 'units', 's') )
       call check( nf90_put_att(ncid, varid1, 'long_name', 'Time series at specific bedplume location') )

       call check( nf90_def_var(ncid, "zhis_bp", NF90_DOUBLE, dimids2, varid2) )
       call check( nf90_put_att(ncid, varid2, 'units', 'm') )
       call check( nf90_put_att(ncid, varid2, 'long_name', 'Z coordinate history at specific bedplume location') )

       if (nfrac>0) then
       call check( nf90_def_var(ncid, "Chis_bp", NF90_DOUBLE, dimids3, varid3) )
       call check( nf90_put_att(ncid, varid3, 'units', '-') )
       call check( nf90_put_att(ncid, varid3, 'long_name', 'Time history of volume concentration 
     & of each fraction at specific bedplume location') )
	 
       call check( nf90_def_var(ncid, "Uhis_bp", NF90_DOUBLE, dimids2, varid4) )
       call check( nf90_put_att(ncid, varid4, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid4, 'long_name', 'U velocity history at specific bedplume location') )	 
       call check( nf90_def_var(ncid, "Vhis_bp", NF90_DOUBLE, dimids2, varid5) )
       call check( nf90_put_att(ncid, varid5, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid5, 'long_name', 'V velocity history at specific bedplume location') )	 
	   call check( nf90_def_var(ncid, "Whis_bp", NF90_DOUBLE, dimids2, varid6) )
       call check( nf90_put_att(ncid, varid6, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid6, 'long_name', 'W velocity history at specific bedplume location') )	 
	endif


	! also add svn info in output files:
       CALL check( nf90_put_att(ncid,nf90_global, "svnversion", trim(svnversion)))
       CALL check( nf90_put_att(ncid,nf90_global, "svnurl", trim(svnurl)))

       ! End define mode. This tells netCDF we are done defining metadata.
       call check( nf90_enddef(ncid) )
     
       ! Write the pretend data to the file. Although netCDF supports
       ! reading and writing subsets of data, in this case we write all the
       ! data in one operation.
       call check( nf90_put_var(ncid, varid1, thisbp) )
       call check( nf90_put_var(ncid, varid2, zhisbp) )
	if (nfrac>0) then
       call check( nf90_put_var(ncid, varid3, Chisbp) )
	endif
	   call check( nf90_put_var(ncid, varid4, Uhisbp) )
	   call check( nf90_put_var(ncid, varid5, Vhisbp) )
	   call check( nf90_put_var(ncid, varid6, Whisbp) )
	   
       call check( nf90_close(ncid) )

	END SUBROUTINE
	
	


	  ! Error handling of netcdf errors 
	  subroutine check(status) 
	    use netcdf

	    integer, intent ( in) :: status
	    integer :: status2

	    if(status /= nf90_noerr) then
	       !UNIT=6 for stdout and UNIT=0 for stderr.
	       write(0,*) trim(nf90_strerror(status))
	       write(0,*) 'closing file'
	       status2 = nf90_close(ncid)
	       if (status2 /= nf90_noerr) then
		  write(0,*) trim(nf90_strerror(status2))
	       end if
	       stop 1
	    end if
	  end 

