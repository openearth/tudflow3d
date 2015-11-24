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


      MODULE nlist
      USE error_functions	

      IMPLICIT NONE
      SAVE
      
      INTEGER i,j,k,imax,jmax,kmax,i1,j1,k1,px,rank,kjet,nmax1,nmax2,istep,CNdiffz,npresIBM,counter
      INTEGER Lmix_type,slip_bot,SEM,azi_n,outflow_overflow_down,azi_n2
      REAL ekm_mol,nu_mol,pi,kappa,gx,gy,gz,Cs,Sc
      REAL dt,time_nm,time_n,time_np,t_end,t0_output,dt_output,te_output,dt_max,tstart_rms,CFL,dt_ini
      REAL dt_output_movie,t0_output_movie,te_output_movie
      REAL U_b,V_b,W_b,rho_b,W_j,Awjet,Aujet,Avjet,Strouhal,radius_j,kn,W_ox,U_bSEM,V_bSEM,U_w,V_w
      REAL U_j2,Awjet2,Aujet2,Avjet2,Strouhal2,radius_j2,zjet2
      REAL xj(4),yj(4),radius_inner_j,W_j_powerlaw,plume_z_outflow_belowsurf
      REAL dphi,dy,dz,schuif_x,depth,lm_min,Rmin,bc_obst_h
      REAL dr_grid(1:100)
      INTEGER imax_grid(1:100)
      INTEGER tmax_inPpunt,tmax_inUpunt,tmax_inVpunt,tmax_inPpuntrand
      INTEGER tmax_inPpuntTSHD,tmax_inUpuntTSHD,tmax_inVpuntTSHD,tmax_inWpuntTSHD
      INTEGER tmax_inUpunt_tauTSHD,tmax_inVpunt_tauTSHD,tmax_inVpunt_rudder
      INTEGER tmax_inWpunt2,tmax_inVpunt2,tmax_inPpunt2,tmax_inWpunt_suction
      INTEGER nfrac,slipvel,interaction_bed,nobst,kbed_bc,nbedplume
      CHARACTER*256 hisfile,restart_dir,inpfile,plumetseriesfile,bcfile,plumetseriesfile2,bedlevelfile
      CHARACTER*3 time_int
      CHARACTER*5 sgs_model
      CHARACTER*4 damping_drho_dz
      REAL damping_a1,damping_b1,damping_a2,damping_b2
      REAL plumetseries(1:100000) 
      REAL plumeUseries(1:100000)
      REAL plumetseries2(1:100000) 
      REAL plumeUseries2(1:100000)
      INTEGER plumeseriesloc,plumeseriesloc2
      INTEGER nr_HPfilter
      REAL timeAB_real(1:4),dpdx,dpdy
      INTEGER periodicx,periodicy,wallup
      REAL U_b3,V_b3,surf_layer
      INTEGER ksurf_bc,kmaxTSHD_ind,nair

      CHARACTER*4 convection,diffusion
      REAL numdiff,comp_filter_a
      INTEGER comp_filter_n

      INTEGER nprop,rudder,softnose
      REAL U_TSHD,LOA,Lfront,Breadth,Draught,Lback,Hback,xfront,yfront,Hfront
      REAL Dprop,xprop,yprop,zprop,Pprop,kn_TSHD,rot_prop
      REAL signU_b,signV_b,signU_bSEM,signV_bSEM
      REAL Dsp,xdh,perc_dh_suction
      CHARACTER*4 draghead

      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inPpunt,j_inPpunt
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inPpuntrand,j_inPpuntrand
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inUpunt,j_inUpunt
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inVpunt,j_inVpunt

      INTEGER*2, DIMENSION(:),ALLOCATABLE :: k_inVpunt2,j_inVpunt2
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: k_inWpunt2,j_inWpunt2
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: k_inPpunt2,j_inPpunt2

      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inPpuntTSHD,j_inPpuntTSHD,k_inPpuntTSHD
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inUpuntTSHD,j_inUpuntTSHD,k_inUpuntTSHD
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inVpuntTSHD,j_inVpuntTSHD,k_inVpuntTSHD
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inWpuntTSHD,j_inWpuntTSHD,k_inWpuntTSHD

      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inUpunt_tauTSHD,j_inUpunt_tauTSHD,k_inUpunt_tauTSHD
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inVpunt_tauTSHD,j_inVpunt_tauTSHD,k_inVpunt_tauTSHD
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inVpunt_rudder,j_inVpunt_rudder,k_inVpunt_rudder
      INTEGER*2, DIMENSION(:),ALLOCATABLE :: i_inWpunt_suction,j_inWpunt_suction,k_inWpunt_suction
      REAL, DIMENSION(:),ALLOCATABLE :: Ubot_TSHD,Vbot_TSHD

      INTEGER*2, DIMENSION(:,:,:),ALLOCATABLE :: llist1,llist2
      INTEGER*2, DIMENSION(:,:),ALLOCATABLE :: llmax1,llmax2,kbed,kbedt,kbed2
      INTEGER, DIMENSION(:,:),ALLOCATABLE :: Xkk,Tii
      INTEGER, DIMENSION(:),ALLOCATABLE :: Xii,Tkk,nfrac_air
      REAL, DIMENSION(:),ALLOCATABLE :: cos_u,cos_v,sin_u,sin_v
      REAL, DIMENSION(:),ALLOCATABLE :: Ru,Rp,dr,Lmix,vol_U,vol_V
      REAL*8, DIMENSION(:),ALLOCATABLE :: xSEM1,ySEM1,zSEM1,uSEM1
      REAL*8, DIMENSION(:),ALLOCATABLE :: lmxSEM1,lmySEM1,lmzSEM1
      REAL*8, DIMENSION(:),ALLOCATABLE :: xSEM2,ySEM2,zSEM2,uSEM2 
      REAL*8, DIMENSION(:,:),ALLOCATABLE :: epsSEM1,epsSEM2
      REAL*8, DIMENSION(:),ALLOCATABLE :: lmxSEM2,lmySEM2,lmzSEM2
      REAL, DIMENSION(:,:),ALLOCATABLE :: azi_angle_p,azi_angle_u,azi_angle_v,zbed,Ubc1,Vbc1,Ubc2,Vbc2,rhocorr_air_z
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: ekm,AA,Diffcof
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: Uold,Vold,Wold,Rold
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: Unew,Vnew,Wnew,Rnew
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: dUdt,dVdt,dWdt,drdt
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: Srr,Spr,Szr,Spp,Spz,Szz
	REAL, DIMENSION(:,:,:),ALLOCATABLE :: Ppropx_dummy,Ppropy_dummy,Ppropz_dummy

      REAL Hs,Tp,Lw,nx_w,ny_w,kabs_w,kx_w,ky_w,om_w

 !     REAL, DIMENSION(:,:,:),ALLOCATABLE :: Uf,Vf,Wf

!      REAL, DIMENSION(:,:,:),ALLOCATABLE :: div
!       REAL, DIMENSION(:,:,:),ALLOCATABLE :: Uavg,Vavg,Wavg,Cavg,Ravg
!       REAL, DIMENSION(:,:,:),ALLOCATABLE :: Urms,Vrms,Wrms,Crms,Rrms
!       REAL, DIMENSION(:,:,:),ALLOCATABLE :: sigU2,sigV2,sigW2,sigC2,sigR2
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: p,pold
      
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: wx,wy,wz,wxold,wyold,wzold,Ppropx,Ppropy,Ppropz
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: wxolder,wyolder,wzolder
      REAL, DIMENSION(:,:),ALLOCATABLE :: Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new
      REAL, DIMENSION(:,:),ALLOCATABLE :: Ub1old,Vb1old,Wb1old,Ub2old,Vb2old,Wb2old
      REAL, DIMENSION(:,:),ALLOCATABLE :: Ubcoarse1,Vbcoarse1,Wbcoarse1,Ubcoarse2,Vbcoarse2,Wbcoarse2
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: Cbcoarse1,Cbcoarse2

	REAL, DIMENSION(:,:,:,:),ALLOCATABLE :: Cold,Cnew,dcdt,cc,ccold,dcdt2
	REAL, DIMENSION(:,:,:),ALLOCATABLE :: Coldbot,Cnewbot,dcdtbot,ccbot,ccoldbot

	type fractions
	    REAL ::	ws,rho,c,dpart,dfloc,n,tau_d,tau_e,M,kn_sed,ws_dep,zair_ref_belowsurf
	end type fractions
	type bed_obstacles
	    REAL ::	x(4),y(4),height,zbottom
	end type bed_obstacles
	type bed_plumes
	    REAL ::	x(4),y(4),height,u,v,w,c(100),t0,zbottom !c(100) matches with size frac_init
	    INTEGER :: forever
	end type bed_plumes
	TYPE(fractions), DIMENSION(:), ALLOCATABLE :: frac
	TYPE(bed_obstacles), DIMENSION(:), ALLOCATABLE :: ob
	TYPE(bed_obstacles), DIMENSION(1000) :: obst !temporary array to read namelist with unknown size
	TYPE(bed_plumes), DIMENSION(:), ALLOCATABLE :: bp
	TYPE(bed_plumes), DIMENSION(25) :: bedplume !temporary array to read namelist with unknown size
      
!       REAL*8, DIMENSION(:,:,:),ALLOCATABLE :: UT
!       REAL*8, DIMENSION(:,:),ALLOCATABLE :: UP,UTMP

      CONTAINS

      subroutine read_namelist

	implicit none

!      include 'param.txt'
!      include 'common.txt'

	integer ios,n
	type frac_init
	  real :: ws,c,rho,dpart,dfloc,tau_d,tau_e,M,kn_sed,ws_dep,zair_ref_belowsurf
	end type frac_init
	TYPE(frac_init), DIMENSION(100) :: fract

! 	TYPE(gridtype), DIMENSION(1) :: grid

	NAMELIST /simulation/px,imax,jmax,kmax,imax_grid,dr_grid,Rmin,schuif_x,dy,depth,hisfile,restart_dir
	NAMELIST /times/t_end,t0_output,dt_output,te_output,tstart_rms,dt_max,dt_ini,time_int,CFL,
     & t0_output_movie,dt_output_movie,te_output_movie
	NAMELIST /num_scheme/convection,numdiff,diffusion,comp_filter_a,comp_filter_n,CNdiffz,npresIBM
	NAMELIST /ambient/U_b,V_b,W_b,bcfile,rho_b,SEM,nmax2,nmax1,lm_min,slip_bot,kn,interaction_bed,periodicx,periodicy,
     & dpdx,dpdy,W_ox,Hs,Tp,nx_w,ny_w,obst,bc_obst_h,U_b3,V_b3,surf_layer,wallup,bedlevelfile,U_bSEM,V_bSEM,U_w,V_w
	NAMELIST /plume/W_j,plumetseriesfile,Awjet,Aujet,Avjet,Strouhal,azi_n,kjet,radius_j,Sc,slipvel,outflow_overflow_down,
     & U_j2,plumetseriesfile2,Awjet2,Aujet2,Avjet2,Strouhal2,azi_n2,radius_j2,zjet2,bedplume,radius_inner_j,xj,yj,W_j_powerlaw,
     & plume_z_outflow_belowsurf
	NAMELIST /LESmodel/sgs_model,Cs,Lmix_type,nr_HPfilter,damping_drho_dz,damping_a1,damping_b1,damping_a2,damping_b2
	NAMELIST /constants/kappa,gx,gy,gz,ekm_mol
	NAMELIST /fractions_in_plume/fract
	NAMELIST /ship/U_TSHD,LOA,Lfront,Breadth,Draught,Lback,Hback,xfront,yfront,kn_TSHD,nprop,Dprop,xprop,yprop,zprop,
     &   Pprop,rudder,rot_prop,draghead,Dsp,xdh,perc_dh_suction,softnose,Hfront

	!! initialise:
	!! simulation:
	px = -999
	imax = -999
	jmax = -999
	kmax = -999
	imax_grid=0
	dr_grid=0.
	Rmin = -999.
	schuif_x = -999.
	dy = -999.
	depth = -999.
	hisfile = ''
	restart_dir = ''
	!! times
	t_end = -999.
	t0_output = 0. ! if not defined than zero
	dt_output = -999.
	te_output = 9.e18 ! if not defined than inf --> continue output towards end simulation
	tstart_rms = 999.
	dt_max = -999.
	dt_ini = -999.
	time_int=''
	CFL = -999.
	t0_output_movie = 9.e18
	dt_output_movie = 9.e18
	te_output_movie = 9.e18
	!! num_scheme
	convection = 'ARGH'
	numdiff = 0.
	diffusion = 'ARGH'
	comp_filter_a = 0.5
	comp_filter_n = 0
	CNdiffz = 0
	npresIBM = 0
	!! ambient:
	U_b = -999.
	V_b = -999.
	W_b = -999.
	bcfile = ''
	rho_b = -999.  
	SEM = -999
	U_bSEM = -999.
	V_bSEM = -999.
	nmax2 = -999
	nmax1 = -999
	lm_min = -999.
	slip_bot = -999
	kn = -999.
	interaction_bed=-999
	periodicx=0
	periodicy=0
	dpdx=0.
	dpdy=0.
	W_ox=0.
	Hs=0. 
	Tp=999.
	nx_w=0.
	ny_w=0.
	U_w=-9999.
	V_w=-9999.	
	bedlevelfile = ''
	DO i=1,4
		obst(:)%x(i) = -99999.
		obst(:)%y(i) = -99999. 
	ENDDO
	obst(:)%height = -99999. 
	obst(:)%zbottom = -99999.
	bc_obst_h = 0.
	U_b3=-9999.
	V_b3=-9999.
	surf_layer=0.
	wallup=0
	!! plume
	W_j = -999.
	plumetseriesfile=''
	Awjet = -999.
	Aujet = -999.
	Avjet = -999.
	Strouhal = -999.
	azi_n = -999
	kjet = -999
	radius_j = -999.
	W_j_powerlaw = 7. ! default 1/7 powerlaw, if 1.e12 block velocity profile
	radius_inner_j = 0. ! default no inner radius where jet inflow is zero
  	xj(1) = -1.e12 ! bounding box for jet --> and inside radius_j and inside bounding box --> default whole area
	xj(2) = 1.e12
  	xj(3) = 1.e12
	xj(4) = -1.e12
	yj(1) = -1.e12
	yj(2) = -1.e12	  
	yj(3) = 1.e12
	yj(4) = 1.e12
	Sc = -999.
	outflow_overflow_down = 0
	U_j2 = -999.
	plumetseriesfile2=''
	Awjet2 = -999.
	Aujet2 = -999.
	Avjet2 = -999.
	Strouhal2 = -999.
	azi_n2 = -999
	zjet2 = -999.
	radius_j2 = -999.
	DO i=1,4
		bedplume(:)%x(i) = -99999.
		bedplume(:)%y(i) = -99999. 
	ENDDO
	bedplume(:)%height = -99999. 
	bedplume(:)%u = -99999. 
	bedplume(:)%v = -99999. 
	bedplume(:)%w = -99999. 
	bedplume(:)%forever = 0
	bedplume(:)%t0 = 0.
	bedplume(:)%zbottom = -99999.
	DO i=1,100
		bedplume(:)%c(i) = 0.
	ENDDO
	plume_z_outflow_belowsurf=-999.
	!! Fractions_in_plume
	fract(:)%ws=-999.
	fract(:)%rho=-999.
	fract(:)%c=-999.
	fract(:)%dpart=-999.
	fract(:)%dfloc=-999.
	fract(:)%tau_e=-999.
	fract(:)%tau_d=-999.
	fract(:)%M=-999.
	fract(:)%kn_sed=-999.
	fract(:)%ws_dep=-999.
	fract(:)%zair_ref_belowsurf=-999.
	nfrac=0
	!! LESmodel
	sgs_model = 'ARGHH'
	Cs = -999.
	Lmix_type = -999
	nr_HPfilter = -999	
	damping_drho_dz = 'none'
	damping_a1 = -999.
	damping_b1 = -999.
	damping_a2 = -999.
	damping_b2 = -999.
	!!constants
	kappa=-999.
	gx = -999.
	gy = -999.
	gz = -999.
	ekm_mol = -999.
	time_nm = 0.
	time_n=0.
	time_np=0.
	!! ship
	U_TSHD=-999.
	LOA=-999.
	Lfront=-999.
	Breadth=-999.
	Draught=-999.
	Lback=-999.
	Hback=-999.
	xfront=-999.
	yfront=-999.
	kn_TSHD=-999.
	nprop=0
	Dprop=-999.
	xprop=-999.
	yprop=-999.
	zprop=-999.
	Pprop=-999.
	rot_prop=-9999.
	rudder=-9
	draghead='none'
	Dsp=0.
	xdh=0.
	perc_dh_suction=0.
	softnose=0
	Hfront=-999.

	CALL GETARG(1,inpfile)
	OPEN(1,FILE=inpfile,IOSTAT=ios,ACTION='read')
	
	IF (ios/=0) THEN
	  write(*,*) 'input file:',inpfile
	  CALL writeerror(200)
	ENDIF
	inpfile=inpfile(9:LEN_TRIM(inpfile))                             
	READ (UNIT=1,NML=simulation,IOSTAT=ios)
	!! checks on input simulation:      
	IF (px<0) CALL writeerror(1)
	IF (imax<0) CALL writeerror(2)
	IF (jmax<0) CALL writeerror(3)
	IF (kmax<0) CALL writeerror(4)
	IF (imax_grid(1).eq.0) CALL writeerror(5)
	IF (dr_grid(1).eq.0.) CALL writeerror(6)
	IF (Rmin<0.) CALL writeerror(7)
	IF (schuif_x<0.) CALL writeerror(8)
	IF (dy<0.) CALL writeerror(9)
	IF (depth<0.) CALL writeerror(10)
	IF (mod(jmax,px).ne.0) CALL writeerror(11)
	IF (mod(kmax,px).ne.0) CALL writeerror(12)
	jmax=jmax/px
	READ (UNIT=1,NML=times,IOSTAT=ios)
	!! check input times  
	IF (t_end<0.) CALL writeerror(30)
	IF (t0_output<0.) CALL writeerror(31)
	IF (dt_output<0.) CALL writeerror(31)
	IF (te_output<0.) CALL writeerror(31)
	IF (tstart_rms<0.) CALL writeerror(32)
	IF (dt_max<0.) CALL writeerror(33) 
	IF (dt_ini<0.) THEN 
	  dt_ini = dt_max
	ELSE
	  dt_ini = MIN(dt_max,dt_ini)
	ENDIF
	IF (time_int.ne.'EE1'.AND.time_int.ne.'RK3'.AND.time_int.ne.'AB2'.AND.time_int.ne.'AB3'.AND.time_int.ne.'ABv') 
     &     CALL writeerror(34) 	 
	IF (time_int.eq.'EE1'.or.time_int.eq.'AB2'.or.time_int.eq.'AB3'.or.time_int.eq.'ABv') THEN
		write(*,*),' WARNING: Your time integration scheme: ',time_int
		write(*,*),' is a testing option and not all functionalities of Dflow3d are working,'
		write(*,*),' use RK3 for a fully supported time integration scheme.'
	ENDIF
	IF (CFL<0.) CALL writeerror(35)

	READ (UNIT=1,NML=num_scheme,IOSTAT=ios)
	!! check input num_scheme
	IF (convection.ne.'CDS2'.AND.convection.ne.'CDS6'.AND.convection.ne.'COM4'.AND.convection.ne.'HYB4'
     &      .AND.convection.ne.'HYB6'.AND.convection.ne.'C4A6' ) CALL writeerror(401) 
	IF (numdiff<0.or.numdiff>1.) CALL writeerror(402)
	IF (diffusion.ne.'CDS2'.AND.diffusion.ne.'COM4') CALL writeerror(403) 
	IF (comp_filter_a<0.or.comp_filter_a>0.5) CALL writeerror(404)
	IF (comp_filter_n<0) CALL writeerror(405)
	numdiff=numdiff*2.  !needed to get correct value (in advec is a 'hidden' factor 2/4)
	IF (CNdiffz.ne.0.and.CNdiffz.ne.1) CALL writeerror(406)
	IF (CNdiffz<0) CALL writeerror(407)

	READ (UNIT=1,NML=ambient,IOSTAT=ios)
	!! check input ambient
	IF (U_b<-998.) CALL writeerror(40)
	IF(U_b<0.) THEN
	  signU_b=-1.
	ELSE 
	  signU_b=1.	
	ENDIF
	IF (V_b<-998.) CALL writeerror(41)
	IF(V_b<0.) THEN
	  signV_b=-1.
	ELSE 
	  signV_b=1.	
	ENDIF
	IF (W_b<0.) CALL writeerror(42)
	IF (rho_b<0.) CALL writeerror(43)  
	IF (SEM<0) CALL writeerror(44)
	IF (nmax2<0) CALL writeerror(45)
	IF (nmax1<0) CALL writeerror(46)
	IF (lm_min<0.) CALL writeerror(47)
	IF (slip_bot<0) CALL writeerror(48)
	IF (kn<0.) CALL writeerror(49)
	IF (interaction_bed<0) CALL writeerror(50)
	IF (periodicx.ne.0.and.periodicx.ne.1.and.periodicx.ne.2) CALL writeerror(52)
	IF (periodicy.ne.0.and.periodicy.ne.1.and.periodicy.ne.2) CALL writeerror(53)
	IF (periodicx.eq.1.and.dpdx.eq.0.) CALL writeerror(54)
	IF (Hs>0.and.Hs>depth) CALL writeerror(55)
	IF (Hs>0.and.Tp.le.0) CALL writeerror(56)
	IF (Hs>0.and.(nx_w>1.or.ny_w>1.or.nx_w<-1.or.ny_w<-1.or.(nx_w**2+ny_w**2)>1.01.or.(nx_w**2+ny_w**2)<0.99)) THEN
	   CALL writeerror(57)	
	ENDIF
	IF (Hs>0.and.U_w<-998.) CALL writeerror(611)
	IF (Hs>0.and.V_w<-998.) CALL writeerror(612)
	
	IF (surf_layer>0.and.U_b3<-999.) CALL writeerror(602)
	IF (surf_layer>0.and.V_b3<-999.) CALL writeerror(603)
	IF (surf_layer<0.or.surf_layer>depth) CALL writeerror(604)
	IF (U_bSEM<-998.) U_bSEM=U_b ! only if defined then different U_bSEM is used else equal to U_b
	IF (V_bSEM<-998.) V_bSEM=V_b
	IF(U_bSEM<0.) THEN
	  signU_bSEM=-1.
	ELSE 
	  signU_bSEM=1.	
	ENDIF
	IF(V_bSEM<0.) THEN
	  signV_bSEM=-1.
	ELSE 
	  signV_bSEM=1.	
	ENDIF

	nobst=0
	DO WHILE (obst(nobst+1)%height.NE.-99999.)
	  nobst=nobst+1
	END DO
	ALLOCATE(ob(nobst))
	DO n=1,nobst
	  ob(n)%x=obst(n)%x
	  ob(n)%y=obst(n)%y
	  ob(n)%height=obst(n)%height
	  ob(n)%zbottom=obst(n)%zbottom
	  DO i=1,4
	    IF (ob(n)%x(i).eq.-99999.or.ob(n)%y(i).eq.-99999) THEN
		write(*,*),' Obstacle:',n
		CALL writeerror(58)
	    ENDIF
	  ENDDO
	  IF (ob(n)%height<0.) THEN
	    write(*,*),' Obstacle:',n
	    CALL writeerror(59)
	  ENDIF
	ENDDO
	IF (bc_obst_h<0.or.bc_obst_h>depth) CALL writeerror(601)
	IF (wallup.ne.0.and.wallup.ne.1) CALL writeerror(605)

	READ (UNIT=1,NML=plume,IOSTAT=ios)
	!! check input plume
	IF (W_j.eq.-999.) CALL writeerror(60)
	IF (Awjet<0.) CALL writeerror(61)
	IF (Aujet<0.) CALL writeerror(62)
	IF (Avjet<0.) CALL writeerror(63)  
	IF (Strouhal<0.) CALL writeerror(69)  
	IF (azi_n<0) CALL writeerror(70)  
	IF (kjet<0) CALL writeerror(64)
	IF (radius_j<0.) CALL writeerror(66)
	IF (Sc<0.) CALL writeerror(67)
	IF (slipvel<0.) CALL writeerror(68)

	IF (U_j2.eq.999.and.radius_j2>0.) CALL writeerror(260)
	IF (U_j2>-999.and.Awjet2<0.) CALL writeerror(261)
	IF (U_j2>-999.and.Aujet2<0.) CALL writeerror(262)
	IF (U_j2>-999.and.Avjet2<0.) CALL writeerror(263)  
	IF (U_j2>-999.and.Strouhal2<0.) CALL writeerror(264)  
	IF (U_j2>-999.and.azi_n2<0) CALL writeerror(265)  
	IF (U_j2>-999.and.radius_j2<0.) CALL writeerror(266)
	IF (U_j2>-999.and.zjet2<0.) CALL writeerror(267)
	IF (U_j2>-999.and.zjet2>depth) CALL writeerror(268)
	IF (radius_inner_j<0.or.(radius_inner_j>radius_j.and.radius_j>0.)) CALL writeerror(275)
	IF (W_j_powerlaw<0.) CALL writeerror(276)

	READ (UNIT=1,NML=fractions_in_plume,IOSTAT=ios)
	nfrac=0
	DO WHILE (fract(nfrac+1)%ws.NE.-999.)
	  nfrac=nfrac+1
	END DO
	ALLOCATE(frac(nfrac))
	ALLOCATE(nfrac_air(nfrac))
	  i=0
	DO n=1,nfrac
	  frac(n)%ws=fract(n)%ws
	  frac(n)%rho=fract(n)%rho
	  frac(n)%c=fract(n)%c
	  frac(n)%dpart=fract(n)%dpart/1000000. !!convert from 10^-6m to m 
	  frac(n)%dfloc=fract(n)%dfloc/1000000. !!convert from 10^-6m to m 
	  frac(n)%tau_e=fract(n)%tau_e
	  frac(n)%tau_d=fract(n)%tau_d
	  frac(n)%M=fract(n)%M
	  frac(n)%kn_sed=fract(n)%kn_sed/1000000. !!convert from 10^-6m to m 
	  frac(n)%ws_dep=fract(n)%ws_dep
	  frac(n)%zair_ref_belowsurf=fract(n)%zair_ref_belowsurf	  
	!! check input fractions_in_plume
	  IF (frac(n)%ws.eq.-999.) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(71)
	  ENDIF
	  IF (frac(n)%rho<0.) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(72)
	  ENDIF		
	  IF (frac(n)%c<0.) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(73)
	  ENDIF
	  IF (frac(n)%c>1.) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(74)
	  ENDIF
	  IF (frac(n)%dpart<0.) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(75)
	  ENDIF
	  IF (frac(n)%dfloc<0.or.(frac(n)%dfloc.ne.0.and.frac(n)%dfloc<frac(n)%dpart)) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(79)
	  ENDIF
	  IF (frac(n)%tau_e<0.) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(76)
	  ENDIF
	  IF (frac(n)%tau_d<0.) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(77)
	  ENDIF
	  IF (frac(n)%M<0.) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(78)
	  ENDIF
	  IF (frac(n)%kn_sed<0.) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(84)
	  ENDIF
	  IF (frac(n)%ws_dep<0.) THEN
		frac(n)%ws_dep=frac(n)%ws  ! if ws_dep is not defined, use ws for deposition
 	  ENDIF
	  IF (frac(n)%zair_ref_belowsurf>-998.) THEN !a reference level for air is defined
		i=i+1
		nfrac_air(i)=n
		write(*,*),'air fraction found: ',n
	  ENDIF
	ENDDO
	  nair=i
	  
	! check op input bedplume after nfrac is known
	nbedplume=0
	DO WHILE (bedplume(nbedplume+1)%height.NE.-99999.)
	  nbedplume=nbedplume+1
	END DO
	ALLOCATE(bp(nbedplume))
	DO n=1,nbedplume
	  bp(n)%x=bedplume(n)%x
	  bp(n)%y=bedplume(n)%y
	  bp(n)%height=bedplume(n)%height
	  bp(n)%u=bedplume(n)%u
	  bp(n)%v=bedplume(n)%v
	  bp(n)%w=bedplume(n)%w
	  bp(n)%c=bedplume(n)%c
	  bp(n)%forever=bedplume(n)%forever
	  bp(n)%t0=bedplume(n)%t0
	  bp(n)%zbottom=bedplume(n)%zbottom
	  DO i=1,4
	    IF (bp(n)%x(i).eq.-99999.or.bp(n)%y(i).eq.-99999) THEN
		write(*,*),' Bedplume:',n
		CALL writeerror(269)
	    ENDIF
	  ENDDO
	  IF (bp(n)%height<0.) THEN
	    write(*,*),' Bedplume:',n
	    CALL writeerror(270)
	  ENDIF
	  DO i=1,nfrac
	    IF (bp(n)%c(i)<0.) THEN
	      write(*,*),' Bedplume:',n
	      write(*,*),' Sediment fraction:',i
	      CALL writeerror(271)
	    ENDIF
	  ENDDO
	  IF (bp(n)%forever.ne.0.and.bp(n)%forever.ne.1) THEN
	    write(*,*),' Bedplume:',n
	    CALL writeerror(272)
	  ENDIF
	  IF (bp(n)%t0.lt.0.) THEN
	    write(*,*),' Bedplume:',n
	    CALL writeerror(273)
	  ENDIF
	  IF (bp(n)%zbottom>bp(n)%height.or.bp(n)%zbottom>depth) THEN
	    write(*,*),' Bedplume:',n
	    CALL writeerror(274)
	  ENDIF
	ENDDO

		
	READ (UNIT=1,NML=LESmodel,IOSTAT=ios)
	!! check input LESmodel
	IF (sgs_model.ne.'SSmag'.and.sgs_model.ne.'FSmag'.and.sgs_model.ne.'SWALE'.and.
     &  sgs_model.ne.'Sigma'.and.sgs_model.ne.'MixLe') CALL writeerror(82)
	IF (Cs<0.) CALL writeerror(80)
	IF (Lmix_type<0) CALL writeerror(81)
	IF (nr_HPfilter<0) CALL writeerror(83)
	IF (damping_drho_dz.ne.'none'.and.damping_drho_dz.ne.'MuAn') CALL writeerror(85)
	IF ((damping_drho_dz.eq.'MuAn'.and.damping_a1.lt.0.).or.(damping_drho_dz.eq.'MuAn'.and.damping_a2.lt.0.)) 
     &  CALL writeerror(86)
	IF ((damping_drho_dz.eq.'MuAn'.and.damping_b1.lt.-990.).or.(damping_drho_dz.eq.'MuAn'.and.damping_b2.lt.-990.)) 
     &  CALL writeerror(87)

	READ (UNIT=1,NML=constants,IOSTAT=ios)
	!! check input constants
	IF (kappa<0.) CALL writeerror(90)
	IF (gx<-800.) CALL writeerror(91)
	IF (gy<-800.) CALL writeerror(92)
	IF (gz<-80000000.) CALL writeerror(93)
	IF (ekm_mol<0.) CALL writeerror(94)

	READ (UNIT=1,NML=ship,IOSTAT=ios)
	!! check input constants
	IF(LOA>0.and.U_TSHD<-998.) CALL writeerror(301)
	IF(LOA>0.and.(Lfront<0.or.Lfront>LOA)) CALL writeerror(302)
	IF(LOA>0.and.Breadth<0) CALL writeerror(303)
	IF(LOA>0.and.(Draught<0.)) CALL writeerror(304)
	IF(LOA>0.and.(Lback<0.or.Lback>LOA)) CALL writeerror(305)
	IF(LOA>0.and.(Hback<0.)) CALL writeerror(306)
	!IF(LOA>0.and.xfront>0) CALL writeerror(307) !switched off 10-4-15 to be able to place TSHD completely out of centre grid
	!IF(LOA>0.and.ABS(yfront)>0.5*Breadth) CALL writeerror(308) !switched off 10-4-15 to be able to place TSHD completely out of centre grid
	IF (nprop<0.or.nprop>2) CALL writeerror(310)
	IF (nprop>0.and.Dprop<0.) CALL writeerror(311)
	IF (nprop>0.and.xprop<0.) CALL writeerror(312)
	IF (nprop.eq.2.and.yprop<0.) CALL writeerror(313)
	IF (nprop>0.and.(zprop<0.or.zprop>depth)) CALL writeerror(314)
	IF (nprop>0.and.Pprop<0.) CALL writeerror(315)
	IF (kn_TSHD.eq.-999.) CALL writeerror(316)
	IF (rot_prop<-9990.) CALL writeerror(317)	
	IF (rudder<0) CALL writeerror(318)
	IF (draghead.ne.'star'.and.draghead.ne.'port'.and.draghead.ne.'both'.and.draghead.ne.'none') CALL writeerror(319)
	IF ((Dsp<0.and.draghead.eq.'star').or.(Dsp<0.and.draghead.eq.'port').or.(Dsp<0.and.draghead.eq.'both')) CALL writeerror(320)
	IF (draghead.eq.'both') THEN
	  perc_dh_suction=perc_dh_suction*0.5 !correction for 2 pipes
	ENDIF
	IF (softnose.ne.0.and.softnose.ne.1) CALL writeerror(321)
	IF (softnose.eq.1.and.Lfront.le.0.) CALL writeerror(322)
	IF (LOA<0.) U_TSHD=0.

	CLOSE(1)

	IF (plumetseriesfile.eq.'') THEN
	ELSE
	   call readtseries(plumetseriesfile,plumetseries,plumeUseries)
	   plumeseriesloc=1
	ENDIF
	IF (plumetseriesfile2.eq.'') THEN
	ELSE
	   call readtseries(plumetseriesfile2,plumetseries2,plumeUseries2)
	   plumeseriesloc2=1
	ENDIF

	END SUBROUTINE read_namelist

	SUBROUTINE allocate_global_vars

	integer ii

	i1=imax+1
	j1=jmax+1
	k1=kmax+1 

	ALLOCATE(i_inPpunt(i1*j1))
	ALLOCATE(j_inPpunt(i1*j1))
        ALLOCATE(i_inPpuntrand(i1*j1))
        ALLOCATE(j_inPpuntrand(i1*j1))
	ALLOCATE(i_inUpunt(i1*j1))
	ALLOCATE(j_inUpunt(i1*j1))
	ALLOCATE(i_inVpunt(i1*j1))
	ALLOCATE(j_inVpunt(i1*j1))
	ALLOCATE(k_inVpunt2(k1*j1))
	ALLOCATE(j_inVpunt2(k1*j1))
	ALLOCATE(k_inWpunt2(k1*j1))
	ALLOCATE(j_inWpunt2(k1*j1))
	ALLOCATE(k_inPpunt2(k1*j1))
	ALLOCATE(j_inPpunt2(k1*j1))

!	ii=(i1+1)*(j1+1)*(k1+1)
!	ALLOCATE(i_inPpuntTSHD(ii))
!	ALLOCATE(j_inPpuntTSHD(ii))
!	ALLOCATE(k_inPpuntTSHD(ii))
!	ALLOCATE(i_inUpuntTSHD(ii))
!	ALLOCATE(j_inUpuntTSHD(ii))
!	ALLOCATE(k_inUpuntTSHD(ii))
!	ALLOCATE(i_inVpuntTSHD(ii))
!	ALLOCATE(j_inVpuntTSHD(ii))
!	ALLOCATE(k_inVpuntTSHD(ii))
!	ALLOCATE(i_inWpuntTSHD(ii))
!	ALLOCATE(j_inWpuntTSHD(ii))
!	ALLOCATE(k_inWpuntTSHD(ii))
	ALLOCATE(Ubot_TSHD(0:j1))
	ALLOCATE(Vbot_TSHD(0:j1))
	ALLOCATE(i_inVpunt_rudder(k1*i1))
	ALLOCATE(j_inVpunt_rudder(k1*i1))
	ALLOCATE(k_inVpunt_rudder(k1*i1))
	ALLOCATE(i_inWpunt_suction(i1*j1))
	ALLOCATE(j_inWpunt_suction(i1*j1))
	ALLOCATE(k_inWpunt_suction(i1*j1))


	ALLOCATE(i_inUpunt_tauTSHD((i1+1)*(j1+1)))
	ALLOCATE(j_inUpunt_tauTSHD((i1+1)*(j1+1)))
	ALLOCATE(k_inUpunt_tauTSHD((i1+1)*(j1+1)))
	ALLOCATE(i_inVpunt_tauTSHD((i1+1)*(j1+1)))
	ALLOCATE(j_inVpunt_tauTSHD((i1+1)*(j1+1)))
	ALLOCATE(k_inVpunt_tauTSHD((i1+1)*(j1+1)))


	Ubot_TSHD=0.
	Vbot_TSHD=0.


	ALLOCATE(cos_u(0:j1))
	ALLOCATE(cos_v(0:j1))
	ALLOCATE(sin_u(0:j1))
	ALLOCATE(sin_v(0:j1))
	ALLOCATE(Ru(0:i1))
	ALLOCATE(Rp(0:i1))	
	ALLOCATE(dr(0:i1))
	ALLOCATE(Lmix(1:imax))
	ALLOCATE(vol_U(1:imax))
	ALLOCATE(vol_V(1:imax))
	ALLOCATE(kbed(0:i1,0:j1))
	ALLOCATE(kbed2(0:i1,0:px*jmax+1))
	ALLOCATE(kbedt(0:i1,0:j1))
	ALLOCATE(zbed(0:i1,0:j1))
	ALLOCATE(rhocorr_air_z(1:nfrac,0:k1))

	IF (SEM.eq.1) THEN
	IF (rank.eq.0.or.rank.eq.px-1) THEN
		ALLOCATE(llist1(0:i1,1:kmax,1:2000))
	ENDIF
	ALLOCATE(llist2(0:j1,1:kmax,1:1000))
	ALLOCATE(AA(3,3,1:kmax))
	ALLOCATE(llmax1(0:i1,1:kmax))
	ALLOCATE(llmax2(0:j1,1:kmax))
	ALLOCATE(xSEM1(nmax1))
	ALLOCATE(ySEM1(nmax1))
	ALLOCATE(zSEM1(nmax1))
	ALLOCATE(epsSEM1(3,nmax1))
	ALLOCATE(uSEM1(nmax1))
	ALLOCATE(lmxSEM1(nmax1))
	ALLOCATE(lmySEM1(nmax1))
	ALLOCATE(lmzSEM1(nmax1))
	ALLOCATE(xSEM2(nmax2))
	ALLOCATE(ySEM2(nmax2))
	ALLOCATE(zSEM2(nmax2))
	ALLOCATE(epsSEM2(3,nmax2))
	ALLOCATE(uSEM2(nmax2))
	ALLOCATE(lmxSEM2(nmax2))
	ALLOCATE(lmySEM2(nmax2))
	ALLOCATE(lmzSEM2(nmax2))
	ENDIF

	ALLOCATE(azi_angle_p(0:i1,0:j1))
	ALLOCATE(azi_angle_u(0:i1,0:j1))
	ALLOCATE(azi_angle_v(0:i1,0:j1))
	ALLOCATE(ekm(0:i1,0:j1,0:k1))
        ALLOCATE(Diffcof(0:i1,0:j1,0:k1))
	ALLOCATE(Uold(0:i1,0:j1,0:k1))  
	ALLOCATE(Vold(0:i1,0:j1,0:k1))
	ALLOCATE(Wold(0:i1,0:j1,0:k1))
	ALLOCATE(Rold(0:i1,0:j1,0:k1))
	ALLOCATE(Unew(0:i1,0:j1,0:k1))  
	ALLOCATE(Vnew(0:i1,0:j1,0:k1))
	ALLOCATE(Wnew(0:i1,0:j1,0:k1))
	ALLOCATE(Rnew(0:i1,0:j1,0:k1))
	ALLOCATE(dUdt(0:i1,0:j1,0:k1))  
	ALLOCATE(dVdt(0:i1,0:j1,0:k1))
	ALLOCATE(dWdt(0:i1,0:j1,0:k1))
	ALLOCATE(dRdt(0:i1,0:j1,0:k1))

!        ALLOCATE(Uf(0:i1,0:j1,0:k1))
!        ALLOCATE(Vf(0:i1,0:j1,0:k1))
!        ALLOCATE(Wf(0:i1,0:j1,0:k1))


	ALLOCATE(Cold(nfrac,0:i1,0:j1,0:k1))
	ALLOCATE(Cnew(nfrac,0:i1,0:j1,0:k1))
	ALLOCATE(dCdt(nfrac,0:i1,0:j1,0:k1))
	ALLOCATE(dCdt2(nfrac,0:i1,0:j1,0:k1))
	ALLOCATE(Coldbot(nfrac,0:i1,0:j1))
	ALLOCATE(Cnewbot(nfrac,0:i1,0:j1))
	ALLOCATE(dCdtbot(nfrac,0:i1,0:j1))

	ALLOCATE(Ppropx(0:i1,0:j1,0:k1))
	ALLOCATE(Ppropy(0:i1,0:j1,0:k1))
	ALLOCATE(Ppropz(0:i1,0:j1,0:k1))

	IF (diffusion.eq.'COM4') THEN
		ALLOCATE(Srr(1:imax,1:jmax,1:kmax))
		ALLOCATE(Spp(1:imax,1:jmax,1:kmax))
		ALLOCATE(Szz(1:imax,1:jmax,1:kmax))
		ALLOCATE(Spr(0:imax,0:jmax,1:kmax))
		ALLOCATE(Szr(0:imax,1:jmax,0:kmax))
		ALLOCATE(Spz(1:imax,0:jmax,0:kmax))
	ENDIF
!	IF(nfrac.eq.1) THEN
!		ALLOCATE(Cold1(0:i1,0:j1,0:k1))
!		ALLOCATE(Cnew1(0:i1,0:j1,0:k1))
!		ALLOCATE(dCdt1(0:i1,0:j1,0:k1))
!		ALLOCATE(cc1(0:i1,0:j1,0:k1))
!	ELSEIF(nfrac.eq.2) THEN
!		ALLOCATE(Cold2(0:i1,0:j1,0:k1))
!		ALLOCATE(Cnew2(0:i1,0:j1,0:k1))
!		ALLOCATE(dCdt2(0:i1,0:j1,0:k1))
!		ALLOCATE(cc2(0:i1,0:j1,0:k1))
!	ENDIF

!	ALLOCATE(div(1:imax,1:jmax,1:kmax))
!	div=0.
! 	ALLOCATE(Uavg(0:i1,0:j1,0:k1))  
! 	ALLOCATE(Vavg(0:i1,0:j1,0:k1))
! 	ALLOCATE(Wavg(0:i1,0:j1,0:k1))
! 	ALLOCATE(Cavg(0:i1,0:j1,0:k1))
! 	ALLOCATE(Ravg(0:i1,0:j1,0:k1))
! 	ALLOCATE(Urms(0:i1,0:j1,0:k1))  
! 	ALLOCATE(Vrms(0:i1,0:j1,0:k1))
! 	ALLOCATE(Wrms(0:i1,0:j1,0:k1))
! 	ALLOCATE(Crms(0:i1,0:j1,0:k1))
! 	ALLOCATE(Rrms(0:i1,0:j1,0:k1))
! 	ALLOCATE(sigU2(0:i1,0:j1,0:k1))  
! 	ALLOCATE(sigV2(0:i1,0:j1,0:k1))
! 	ALLOCATE(sigW2(0:i1,0:j1,0:k1))
! 	ALLOCATE(sigC2(0:i1,0:j1,0:k1))
! 	ALLOCATE(sigR2(0:i1,0:j1,0:k1))
	ALLOCATE(p(1:imax,1:jmax,1:kmax))
	ALLOCATE(pold(1:imax,1:jmax,1:kmax))
	p=0.
	pold=0.
	IF (time_int.eq.'AB2'.or.time_int.eq.'AB3'.or.time_int.eq.'ABv') THEN
		ALLOCATE(wx(0:i1,0:j1,0:k1))  
		ALLOCATE(wy(0:i1,0:j1,0:k1))
		ALLOCATE(wz(0:i1,0:j1,0:k1))
		ALLOCATE(wxold(0:i1,0:j1,0:k1))  
		ALLOCATE(wyold(0:i1,0:j1,0:k1))
		ALLOCATE(wzold(0:i1,0:j1,0:k1))
		ALLOCATE(wxolder(0:i1,0:j1,0:k1))  
		ALLOCATE(wyolder(0:i1,0:j1,0:k1))
		ALLOCATE(wzolder(0:i1,0:j1,0:k1))
		ALLOCATE(cc(nfrac,0:i1,0:j1,0:k1))
		ALLOCATE(ccold(nfrac,0:i1,0:j1,0:k1))
		ALLOCATE(ccbot(nfrac,0:i1,0:j1))
		ALLOCATE(ccoldbot(nfrac,0:i1,0:j1))
	ENDIF

	ALLOCATE(Ub1new(0:i1,0:k1))
	ALLOCATE(Vb1new(0:i1,0:k1))
	ALLOCATE(Wb1new(0:i1,0:k1))
	ALLOCATE(Ub2new(0:j1,0:k1))
	ALLOCATE(Vb2new(0:j1+1,0:k1))
	ALLOCATE(Wb2new(0:j1,0:k1))	
	ALLOCATE(Ub1old(0:i1,0:k1))
	ALLOCATE(Vb1old(0:i1,0:k1))
	ALLOCATE(Wb1old(0:i1,0:k1))
	ALLOCATE(Ub2old(0:j1,0:k1))
	ALLOCATE(Vb2old(0:j1+1,0:k1))
	ALLOCATE(Wb2old(0:j1,0:k1))    
	ALLOCATE(Ubcoarse1(1:imax,1:kmax))
	ALLOCATE(Vbcoarse1(1:imax,1:kmax))
	ALLOCATE(Wbcoarse1(1:imax,1:kmax))
	ALLOCATE(Ubcoarse2(0:j1,1:kmax))
	ALLOCATE(Vbcoarse2(0:j1,1:kmax))
	ALLOCATE(Wbcoarse2(0:j1,1:kmax)) 
	ALLOCATE(Cbcoarse1(nfrac,1:imax,1:kmax))
	ALLOCATE(Cbcoarse2(nfrac,0:j1,1:kmax)) 
	ALLOCATE(Ubc1(1:imax,1:kmax))
	ALLOCATE(Vbc1(1:imax,1:kmax))
	ALLOCATE(Ubc2(0:j1,1:kmax))
	ALLOCATE(Vbc2(0:j1,1:kmax))
	! init at zero (when no bcfile it stays zero)
	Ubcoarse1=0.
	Ubcoarse2=0.
	Vbcoarse1=0.
	Vbcoarse2=0.
	Wbcoarse1=0.
	Wbcoarse2=0.
	Cbcoarse1=0.
	Cbcoarse2=0.

	ALLOCATE(Xii(1:kmax))
	ALLOCATE(Tkk(1:jmax*px))
	ALLOCATE(Xkk(kmax,jmax))
	ALLOCATE(Tii(jmax*px,kmax/px))

! 	ALLOCATE(UT(0:i1,0:j1,0:k1))
! 	ALLOCATE(UP(0:i1,0:k1))
! 	ALLOCATE(UTMP(0:i1,0:k1))  

	Ub1new=0.
	Vb1new=0.
	Wb1new=0.
	Ub2new=0.
	Vb2new=0.
	Wb2new=0.
	Ub1old=0.
	Vb1old=0.
	Wb1old=0.
	Ub2old=0.
	Vb2old=0.
	Wb2old=0.

	tmax_inUpuntTSHD=0
	tmax_inVpuntTSHD=0
	tmax_inPpuntTSHD=0

	nu_mol=ekm_mol/rho_b
  

	END SUBROUTINE allocate_global_vars

	SUBROUTINE readtseries(seriesfile,tseries,series)
!	SUBROUTINE readtseries(seriesfile,tseries,series)
!	reads a time series from seriesfile
!	Lynyrd de Wit, October 2010

	CHARACTER*256 seriesfile
	REAL tseries(1:100000),series(1:100000)
	REAL begintime,timestep,endtime
	INTEGER i,k,p,ios

	NAMELIST /timeseries/begintime,timestep,endtime,series
	series=-9999.
	OPEN(10,FILE=seriesfile,IOSTAT=ios,ACTION='read')

		write(*,*) 'File :', seriesfile,'inlezen'

	IF (ios/=0) THEN
		write(*,*) 'File :', seriesfile
		CALL writeerror(1001)
	END IF
	READ (UNIT=10,NML=timeseries,IOSTAT=ios)
	IF (ios/=0) THEN
		write(*,*) 'File :', seriesfile
		CALL writeerror(1002)
	END IF
	CLOSE(10)
	k=1
	DO WHILE (series(k).NE.-9999)
		tseries(k)=begintime+REAL(k)*timestep-timestep
		k=k+1
	END DO
!		write(*,*) 'File :', seriesfile,'ingelezen'
!		write(*,*) 'tseries:', tseries(1:100)
!		write(*,*) 'series:', series(1:100)

	END SUBROUTINE readtseries


	END MODULE nlist
