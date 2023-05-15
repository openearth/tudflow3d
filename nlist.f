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
      
      INTEGER i,j,k,imax,jmax,kmax,i1,j1,k1,px,rank,kjet,nmax1,nmax2,nmax3,istep,CNdiffz,npresIBM,counter,npresPRHO,oPRHO,k_pzero
      INTEGER Lmix_type,slip_bot,SEM,azi_n,outflow_overflow_down,azi_n2,wiggle_detector,wd,applyVOF,Poutflow,k_ust_tau,Uoutflow
	  INTEGER k_ust_tau_flow
      REAL ekm_mol,nu_mol,pi,kappa,gx,gy,gz,Cs,Sc,calibfac_sand_pickup,calibfac_Shields_cr,morfac,morfac2,calibfac_sand_bedload
      REAL dt,time_nm,time_n,time_np,t_end,t0_output,dt_output,te_output,dt_max,tstart_rms,CFL,dt_ini,tstart_morf,trestart,dt_old
      REAL dt_output_movie,t0_output_movie,te_output_movie,te_rms,time_nm2,tstart_morf2,fcor
      REAL U_b,V_b,W_b,rho_b,W_j,Awjet,Aujet,Avjet,Strouhal,radius_j,kn,W_ox,U_bSEM,V_bSEM,U_w,V_w,U_init,V_init
      REAL U_j2,Awjet2,Aujet2,Avjet2,Strouhal2,radius_j2,zjet2,rho_b2
      REAL xj(4),yj(4),radius_inner_j,W_j_powerlaw,plume_z_outflow_belowsurf
      REAL dy,dz,schuif_x,depth,lm_min,lm_min3,Rmin,bc_obst_h
      REAL dr_grid(1:100),dy_grid(1:100),fac_y_grid(1:100),lim_y_grid(1:100),fac_r_grid(1:100),lim_r_grid(1:100)
      INTEGER imax_grid(1:100),jmax_grid(1:100),sym_grid_y
      INTEGER tmax_inPpunt,tmax_inUpunt,tmax_inVpunt,tmax_inPpuntrand
      INTEGER tmax_inPpuntTSHD,tmax_inUpuntTSHD,tmax_inVpuntTSHD,tmax_inWpuntTSHD
      INTEGER tmax_inUpunt_tauTSHD,tmax_inVpunt_tauTSHD,tmax_inVpunt_rudder
      INTEGER tmax_inWpunt2,tmax_inVpunt2,tmax_inPpunt2,tmax_inWpunt_suction
      INTEGER nfrac,slipvel,interaction_bed,nobst,nbedplume,continuity_solver,hindered_settling,settling_along_gvector
      INTEGER nfr_silt,nfr_sand,nfr_air,istart_morf2,i_periodicx,hindered_settling_c
      CHARACTER*256 hisfile,restart_dir,inpfile,plumetseriesfile,bcfile,plumetseriesfile2,bedlevelfile,initconditionsfile
	  CHARACTER*256 U_b_tseriesfile,V_b_tseriesfile,W_b_tseriesfile
      CHARACTER*3 time_int,advec_conc,cutter,split_rho_cont
      CHARACTER*5 sgs_model
      CHARACTER*4 damping_drho_dz,extra_mix_visc
	  CHARACTER*11 pickup_formula,bedload_formula
	  CHARACTER*8 transporteq_fracs
	  CHARACTER*20 pickup_correction
      REAL damping_a1,damping_b1,damping_a2,damping_b2,cfixedbed
      REAL plumetseries(1:10000) 
      REAL plumeUseries(1:10000)
      REAL plumetseries2(1:10000) 
      REAL plumeUseries2(1:10000),c_bed(100)
	  REAL U_b_series(1:10000),V_b_series(1:10000),W_b_series(1:10000)
	  REAL U_b_tseries(1:10000),V_b_tseries(1:10000),W_b_tseries(1:10000)
      INTEGER plumeseriesloc,plumeseriesloc2,plumeQseriesloc,plumecseriesloc
	  INTEGER U_b_seriesloc,V_b_seriesloc,W_b_seriesloc
      INTEGER nr_HPfilter,depo_implicit,depo_cbed_option,monopile
      REAL timeAB_real(1:4),dpdx,dpdy,kn_d50_multiplier,avalanche_slope(100),av_slope_z(100)
	  REAL dpdx1,dpdy1,Uavold,Vavold,U3avold,V3avold,dpdx3,dpdy3
      INTEGER periodicx,periodicy,wallup,dUVdn_IBMbed
      REAL U_b3,V_b3,surf_layer,reduction_sedimentation_shields,kn_mp,kn_sidewalls
      INTEGER ksurf_bc,kmaxTSHD_ind,nair
      INTEGER poissolver,nm1,istep_output_bpmove,avalanche_until_done,IBMorder
      INTEGER iparm(64)
      CHARACTER*256 plumeQtseriesfile,plumectseriesfile,avfile,obstfile,kn_flow_file      
      REAL Q_j,plumeQseries(1:10000),plumeQtseries(1:10000),plumectseries(1:10000),plumecseries(30,1:10000) !c(30) matches with size frac_init
      REAL Aplume,driftfluxforce_calfac,kn_flow_d50_multiplier
	  REAL vwal,vwal2,delta_nsed,nl,permeability_kl,pickup_fluctuations_ampl,z_tau_sed,kn_d50_multiplier_bl,bl_relax,power_VR2019
	  INTEGER pickup_fluctuations,cbed_method,k_layer_pickup,nu_minimum_wall,pickup_bedslope_geo,wbed_correction,bedslope_effect
	  INTEGER wallmodel_tau_sed,ndtbed,nrmsbed,telUVWbed,tel_dt
	  REAL Const1eps,Const2,Sc_k,Sc_eps,Cal_buoyancy_k,Cal_buoyancy_eps,Cs_relax
	  REAL bedslope_mu_s,alfabs_bl,alfabn_bl,phi_sediment
	  REAL dt_factor_avg,CNdiff_factor,CNdiff_ho,CNdiff_dtfactor,CNdiff_tol,Apvisc_shear_relax
	  INTEGER n_dtavg,CNdiff_pc,CNdiff_maxi,CNdiff_ini,initPhydrostatic,npresIBM_viscupdate,rheo_shear_method,Apvisc_force_eq
	  INTEGER tmax_inPpuntTSHDini,tmax_inUpuntTSHDini,tmax_inVpuntTSHDini,tmax_inWpuntTSHDini
	  INTEGER nobst_files,nobst_file,momentum_exchange_obstacles,erosion_cbed_start,movebed_absorb_cfluid
	  CHARACTER(len=256) :: obst_file_series(5000)
	  REAL :: obst_starttimes(5000)
	  INTEGER cbc_perx_j(2)
	  REAL cbc_relax
	  
	  
	  !new variables rheology
	  INTEGER Non_Newtonian, Apvisc_interp
	  REAL SIMPLE_tauy,SIMPLE_muB,SIMPLE_climit(30) 
	  REAL JACOBS_Ky,JACOBS_Kmu,JACOBS_Aclay,JACOBS_By,JACOBS_Bmu,JACOBS_muw
	  REAL WINTER_Ay,WINTER_Amu,WINTER_nf,WINTER_af,WINTER_muw
	  REAL THOMAS_Cy,THOMAS_Cmu,THOMAS_ky,THOMAS_kmu,THOMAS_Py,THOMAS_Pmu,THOMAS_phi_sand_max
	  REAL BAGNOLD_beta,BAGNOLD_phi_max,PAPANASTASIOUS_m,shear0limit
	  REAL Lambda_init,Kin_eq_a,Kin_eq_b,Kin_eq_lambda_0
	  REAL HOUSKA_n,HOUSKA_eta_0,HOUSKA_eta_inf,HOUSKA_tauy_0,HOUSKA_tauy_inf
	  

	  REAL, DIMENSION(:,:,:),ALLOCATABLE :: tauY,muB
	  REAL, DIMENSION(:,:,:),ALLOCATABLE :: stress,strain,muA
	  REAL, DIMENSION(:,:,:),ALLOCATABLE :: lambda_old,lambda_new
	  CHARACTER*6 Rheological_model	  
	  
	  
      CHARACTER*4 convection,diffusion
      REAL numdiff,comp_filter_a,numdiff2
      INTEGER comp_filter_n,pres_in_predictor_step,pres_in_predictor_step_internal

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
      REAL, DIMENSION(:),ALLOCATABLE :: Ubot_TSHD,Vbot_TSHD,dt_series !,facIBMu,facIBMv,facIBMw,
      INTEGER*8, DIMENSION(:),ALLOCATABLE :: pt,pt3

      INTEGER*2, DIMENSION(:,:,:),ALLOCATABLE :: llist1,llist2,llist3
      INTEGER*2, DIMENSION(:,:),ALLOCATABLE :: llmax1,llmax2,llmax3 !,kbed,kbedt,kbed2
      INTEGER*8, DIMENSION(:,:),ALLOCATABLE :: kbed,kbedt,kbed0,kbedold !,kbed2,kbed22
      INTEGER, DIMENSION(:,:),ALLOCATABLE :: Xkk,Tii
      INTEGER, DIMENSION(:),ALLOCATABLE :: Xii,Tkk,nfrac_air,nfrac_silt,nfrac_sand,nfrac_air2
	  INTEGER*8, DIMENSION(:),ALLOCATABLE :: kbedin
      REAL, DIMENSION(:),ALLOCATABLE :: cos_u,cos_v,sin_u,sin_v,cos_ut,sin_ut,cos_vt,sin_vt
      REAL, DIMENSION(:),ALLOCATABLE :: Ru,Rp,dr,phivt,phipt,dphi2t,phiv,phip,dphi2,b_update
      REAL*8, DIMENSION(:),ALLOCATABLE :: xSEM1,ySEM1,zSEM1,uSEM1
      REAL*8, DIMENSION(:),ALLOCATABLE :: lmxSEM1,lmySEM1,lmzSEM1
      REAL*8, DIMENSION(:),ALLOCATABLE :: xSEM2,ySEM2,zSEM2,uSEM2 
      REAL*8, DIMENSION(:,:),ALLOCATABLE :: epsSEM1,epsSEM2,epsSEM3,Lmix2,Lmix2hat,vol_V,vol_Vp
      REAL*8, DIMENSION(:),ALLOCATABLE :: lmxSEM2,lmySEM2,lmzSEM2
	  REAL*8, DIMENSION(:),ALLOCATABLE :: rSEM3,thetaSEM3,zSEM3,wSEM3,xSEM3,ySEM3
	  REAL, DIMENSION(:,:,:,:),ALLOCATABLE :: AA1,AA2,AA3
      REAL*8, DIMENSION(:),ALLOCATABLE :: lmrSEM3,lmzSEM3
      REAL, DIMENSION(:,:),ALLOCATABLE :: azi_angle_p,azi_angle_u,azi_angle_v,zbed,Ubc1,Vbc1,Ubc2,Vbc2,rhocorr_air_z,Wbed,wscorr_z
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: ekm,Diffcof,bednotfixed,bednotfixed_depo
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: Uold,Vold,Wold,Rold
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: Unew,Vnew,Wnew,Rnew
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: dUdt,dVdt,dWdt,drdt
	  REAL, DIMENSION(:,:,:),ALLOCATABLE :: rhoU,rhoV,rhoW
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: Srr,Spr,Szr,Spp,Spz,Szz,TKE,EEE,Cmu
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: Ppropx_dummy,Ppropy_dummy,Ppropz_dummy,wf,av_slope 
      INTEGER, DIMENSION(:),ALLOCATABLE, TARGET :: jco,iro,beg,di,di2
      REAL Hs,Tp,Lw,nx_w,ny_w,kabs_w,kx_w,ky_w,om_w
      REAL, DIMENSION(:),ALLOCATABLE :: LUB,LHS2 
      REAL, DIMENSION(:,:),ALLOCATABLE, TARGET :: lhs,LUBs,ubot
      INTEGER, DIMENSION(:),ALLOCATABLE, TARGET :: jco3,iro3,beg3,di3 
      REAL, DIMENSION(:),ALLOCATABLE, TARGET :: lhs3,sigtbed 
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: rhs3,d_cbotdelay
	  REAL, DIMENSION(:,:,:),ALLOCATABLE :: rhU,rhV,rhW
	  REAL, DIMENSION(:,:,:,:),ALLOCATABLE :: cU,cV,cW
	  REAL, DIMENSION(:,:),ALLOCATABLE :: thisbp,zhisbp
	  REAL, DIMENSION(:,:,:),ALLOCATABLE :: Chisbp
	  REAL, DIMENSION(:,:),ALLOCATABLE :: Uhisbp,Vhisbp,Whisbp
	  REAL*8, DIMENSION(:,:,:),ALLOCATABLE :: fc_global
	  REAL, DIMENSION(:,:),ALLOCATABLE :: uuR_relax,uuL_relax,vvR_relax,vvL_relax,tau2Vold,tau2Vnew,tau2Wold,tau2Wnew,qb_relax
	  REAL, DIMENSION(:,:),ALLOCATABLE :: tau_fl_Uold,tau_fl_Vold,tau_fl_Unew,tau_fl_Vnew,ust_sl_new,ust_sl_old,ust_bl_new,ust_bl_old
	  REAL, DIMENSION(:,:),ALLOCATABLE :: ust_mud_new,ust_mud_old,kn_flow
	  REAL, DIMENSION(:,:,:),ALLOCATABLE :: qbU,qbV,ust_frac_old,ust_frac_new,sigUWbed,sigVWbed,sigUbed,sigVbed,sigWbed

 !     REAL, DIMENSION(:,:,:),ALLOCATABLE :: Uf,Vf,Wf

!      REAL, DIMENSION(:,:,:),ALLOCATABLE :: div
!       REAL, DIMENSION(:,:,:),ALLOCATABLE :: Uavg,Vavg,Wavg,Cavg,Ravg
!       REAL, DIMENSION(:,:,:),ALLOCATABLE :: Urms,Vrms,Wrms,Crms,Rrms
!       REAL, DIMENSION(:,:,:),ALLOCATABLE :: sigU2,sigV2,sigW2,sigC2,sigR2
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: p,pold,dp,pold1,pold2,pold3,Csgrid,phdt,phnew,viscf
      
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: wx,wy,wz,wxold,wyold,wzold,Ppropx,Ppropy,Ppropz
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: wxolder,wyolder,wzolder
      REAL, DIMENSION(:,:),ALLOCATABLE :: Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new
      REAL, DIMENSION(:,:),ALLOCATABLE :: Ub1old,Vb1old,Wb1old,Ub2old,Vb2old,Wb2old,Ub3old,Vb3old,Wb3old
      REAL, DIMENSION(:,:),ALLOCATABLE :: Ubcoarse1,Vbcoarse1,Wbcoarse1,Ubcoarse2,Vbcoarse2,Wbcoarse2
      REAL, DIMENSION(:,:,:),ALLOCATABLE :: Cbcoarse1,Cbcoarse2

	REAL, DIMENSION(:,:,:,:),ALLOCATABLE :: Cold,Cnew,dcdt,cc,ccold,dcdt2,Clivebed
	REAL, DIMENSION(:,:,:),ALLOCATABLE :: Coldbot,Cnewbot,dcdtbot,ccbot,ccoldbot
	
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: Uavg,Vavg,Wavg,Ravg,Pavg,muavg
		REAL, DIMENSION(:,:,:), ALLOCATABLE :: Umax,Vmax,Wmax,Uhormax,U3dmax !,Umin,Vmin,Wmin
	  REAL, DIMENSION(:,:,:), ALLOCATABLE :: sigU2,sigV2,sigW2,sigR2,sigUV,sigUW,sigVW
	  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: sigC2,Cavg,sigUC,sigVC,sigWC,Cmax,Cmin
	  REAL, DIMENSION(:,:), ALLOCATABLE :: sig_tau_flow2,tau_flow_avg
	  REAL, DIMENSION(:,:), ALLOCATABLE :: tau_sl_avg,tau_bl_avg,tau_mud_avg,sig_tau_sl2,sig_tau_bl2,sig_tau_mud2
	  REAL, DIMENSION(:,:,:), ALLOCATABLE :: tau_frac_avg,sig_tau_frac2
	  
	  INTEGER stat_count
	REAL stat_time_count	

	type fractions
	    REAL ::	ws,rho,c,dpart,dfloc,n,tau_d,tau_e,M,kn_sed,ws_dep,zair_ref_belowsurf,CD,cmax
	    INTEGER :: type
	end type fractions
	type bed_obstacles
	    REAL ::	x(4),y(4),height,zbottom,ero,depo
	end type bed_obstacles
	type bed_plumes
	    REAL ::	x(4),y(4),height,u,v,w,c(30),t0,t_end,zbottom,Q,sedflux(30),volncells,changesedsuction !c(30) matches with size frac_init
		REAL :: move_zbed_criterium(100000),move_dx_series(100000),move_dy_series(100000),move_dz_series(100000)
		INTEGER*2 :: move_zbed_type(100000)
		REAL :: move_nx_series(100000),move_ny_series(100000),x2(4),y2(4),move_dx2_series(100000),move_dy2_series(100000)
		REAL :: move_outputfile_series(100000),uinput,dt_history,t_bphis_output,move_dz_height_factor,move_dz_zbottom_factor
	    INTEGER :: forever,h_seriesloc,zb_seriesloc,Q_seriesloc,S_seriesloc,c_seriesloc,velocity_force
		INTEGER :: u_seriesloc,v_seriesloc,w_seriesloc,istep_bphis_output
		CHARACTER*256 :: h_tseriesfile,zb_tseriesfile,Q_tseriesfile,S_tseriesfile,c_tseriesfile
		CHARACTER*256 :: u_tseriesfile,v_tseriesfile,w_tseriesfile
		REAL :: move_u,move_v,move_w,radius,radius2
	end type bed_plumes
	type bed_plumes2
	    REAL ::	x(4),y(4),height,u,v,w,c(30),t0,t_end,zbottom,Q,sedflux(30),volncells,changesedsuction !c(30) matches with size frac_init
		REAL :: h_tseries(10000),h_series(10000),zb_tseries(10000),zb_series(10000)
		REAL :: Q_tseries(10000),c_tseries(10000),S_tseries(10000),Q_series(10000),c_series(30,10000),S_series(30,10000)
		REAL :: move_zbed_criterium(100000),move_dx_series(100000),move_dy_series(100000),move_dz_series(100000)
		INTEGER*2 :: move_zbed_type(100000)
		REAL :: move_nx_series(100000),move_ny_series(100000),x2(4),y2(4),move_dx2_series(100000),move_dy2_series(100000)
		REAL :: move_outputfile_series(100000),uinput,dt_history,t_bphis_output,move_dz_height_factor,move_dz_zbottom_factor
		REAL :: u_tseries(10000),v_tseries(10000),w_tseries(10000)
		REAL :: u_series(10000),v_series(10000),w_series(10000)
	    INTEGER :: forever,h_seriesloc,zb_seriesloc,Q_seriesloc,S_seriesloc,c_seriesloc,nmove,nmove_present,velocity_force,tmax
		INTEGER :: u_seriesloc,v_seriesloc,w_seriesloc,istep_bphis_output
		CHARACTER*256 :: h_tseriesfile,zb_tseriesfile,Q_tseriesfile,S_tseriesfile,c_tseriesfile
		CHARACTER*256 :: u_tseriesfile,v_tseriesfile,w_tseriesfile
!		INTEGER :: tmax_iP,tmax_iU,iP_inbp(10000),jP_inbp(10000),kP_inbp(10000),iU_inbp(10000),jU_inbp(10000),kU_inbp(10000)
		REAL :: move_u,move_v,move_w,radius,radius2
		INTEGER*2 :: i(100000),j(100000) 
	end type bed_plumes2
	
	TYPE(fractions), DIMENSION(:), ALLOCATABLE :: frac
	TYPE(bed_obstacles), DIMENSION(:), ALLOCATABLE :: ob,obst
	TYPE(bed_plumes), DIMENSION(:), ALLOCATABLE :: bedplume
	TYPE(bed_plumes2), DIMENSION(:), ALLOCATABLE :: bp

      
!       REAL*8, DIMENSION(:,:,:),ALLOCATABLE :: UT
!       REAL*8, DIMENSION(:,:),ALLOCATABLE :: UP,UTMP

      CONTAINS

      subroutine read_namelist

	implicit none

!      include 'param.txt'
!      include 'common.txt'

	integer ios,n,n1,n2,n3,n4
	type frac_init
	  real :: ws,c,rho,dpart,dfloc,tau_d,tau_e,M,kn_sed,ws_dep,zair_ref_belowsurf,CD,cmax
	  integer :: type
	end type frac_init
	TYPE(frac_init), DIMENSION(30) :: fract

	ALLOCATE(bedplume(400)) !temporary array to read namelist with unknown size
	ALLOCATE(obst(100000)) !temporary array to read namelist with unknown size
! 	TYPE(gridtype), DIMENSION(1) :: grid

	NAMELIST /simulation/px,imax,jmax,kmax,imax_grid,dr_grid,Rmin,schuif_x,dy,depth,hisfile,restart_dir
     & ,lim_r_grid,fac_r_grid,jmax_grid,lim_y_grid,fac_y_grid,sym_grid_y,dy_grid
	NAMELIST /times/t_end,t0_output,dt_output,te_output,tstart_rms,dt_max,dt_ini,time_int,CFL,
     & t0_output_movie,dt_output_movie,te_output_movie,tstart_morf,te_rms,tstart_morf2,n_dtavg
	NAMELIST /num_scheme/convection,numdiff,wiggle_detector,diffusion,comp_filter_a,comp_filter_n,CNdiffz,npresIBM,advec_conc,
     & continuity_solver,transporteq_fracs,split_rho_cont,driftfluxforce_calfac,depo_implicit,IBMorder,npresPRHO,
     & pres_in_predictor_step,Poutflow,oPRHO,applyVOF,k_ust_tau,Uoutflow,dUVdn_IBMbed,k_pzero,numdiff2,CNdiff_factor,CNdiff_ho,
     & CNdiff_dtfactor,CNdiff_pc,CNdiff_maxi,CNdiff_tol,CNdiff_ini,initPhydrostatic,npresIBM_viscupdate,momentum_exchange_obstacles
     & ,k_ust_tau_flow	 
	NAMELIST /ambient/U_b,V_b,W_b,bcfile,rho_b,SEM,nmax2,nmax1,nmax3,lm_min,lm_min3,slip_bot,kn,interaction_bed,
     & periodicx,periodicy,dpdx,dpdy,W_ox,Hs,Tp,nx_w,ny_w,obst,bc_obst_h,U_b3,V_b3,surf_layer,wallup,bedlevelfile,
     & U_bSEM,V_bSEM,U_w,V_w,c_bed,cfixedbed,U_init,V_init,initconditionsfile,rho_b2,monopile,kn_mp,kn_sidewalls,obstfile
     & ,istart_morf2,i_periodicx,kn_flow_file,kn_flow_d50_multiplier,U_b_tseriesfile,V_b_tseriesfile,W_b_tseriesfile
     & ,cbc_perx_j,cbc_relax
	NAMELIST /plume/W_j,plumetseriesfile,Awjet,Aujet,Avjet,Strouhal,azi_n,kjet,radius_j,Sc,slipvel,outflow_overflow_down,
     & U_j2,plumetseriesfile2,Awjet2,Aujet2,Avjet2,Strouhal2,azi_n2,radius_j2,zjet2,bedplume,radius_inner_j,xj,yj,W_j_powerlaw,
     & plume_z_outflow_belowsurf,hindered_settling,hindered_settling_c,Q_j,plumeQtseriesfile,plumectseriesfile
	NAMELIST /LESmodel/sgs_model,Cs,Lmix_type,nr_HPfilter,damping_drho_dz,damping_a1,damping_b1,damping_a2,damping_b2,
     & extra_mix_visc,nu_minimum_wall,Const1eps,Const2,Sc_k,Sc_eps,Cal_buoyancy_k,Cal_buoyancy_eps,Cs_relax
	NAMELIST /constants/kappa,gx,gy,gz,ekm_mol,calibfac_sand_pickup,pickup_formula,kn_d50_multiplier,avalanche_slope,
     &	av_slope_z,calibfac_Shields_cr,reduction_sedimentation_shields,morfac,morfac2,avalanche_until_done,avfile,
     & settling_along_gvector,vwal,nl,permeability_kl,pickup_fluctuations_ampl,pickup_fluctuations,pickup_correction,cbed_method,
     & z_tau_sed,k_layer_pickup,pickup_bedslope_geo,bedload_formula,kn_d50_multiplier_bl,calibfac_sand_bedload,bl_relax,fcor,
     & wbed_correction,bedslope_effect,bedslope_mu_s,alfabs_bl,alfabn_bl,phi_sediment,wallmodel_tau_sed,nrmsbed,ndtbed,
     & erosion_cbed_start,movebed_absorb_cfluid,power_VR2019	 
	NAMELIST /fractions_in_plume/fract
	NAMELIST /ship/U_TSHD,LOA,Lfront,Breadth,Draught,Lback,Hback,xfront,yfront,kn_TSHD,nprop,Dprop,xprop,yprop,zprop,
     &   Pprop,rudder,rot_prop,draghead,Dsp,xdh,perc_dh_suction,softnose,Hfront,cutter
	NAMELIST /rheology/Non_Newtonian,Rheological_model,PAPANASTASIOUS_m,shear0limit,Apvisc_interp,SIMPLE_tauy,SIMPLE_muB,SIMPLE_climit,
     & JACOBS_Ky,JACOBS_Kmu,JACOBS_Aclay,JACOBS_By,JACOBS_Bmu,JACOBS_muw,
     & WINTER_Ay,WINTER_Amu,WINTER_nf,WINTER_af,WINTER_muw,
     & THOMAS_Cy,THOMAS_Cmu,THOMAS_ky,THOMAS_kmu,THOMAS_Py,THOMAS_Pmu,THOMAS_phi_sand_max,
     & Lambda_init,Kin_eq_a,Kin_eq_b,Kin_eq_lambda_0,HOUSKA_n,HOUSKA_eta_0,HOUSKA_eta_inf,HOUSKA_tauy_0,HOUSKA_tauy_inf,
     & BAGNOLD_beta,BAGNOLD_phi_max,rheo_shear_method,Apvisc_shear_relax,Apvisc_force_eq

	!! initialise:
	!! simulation:
	px = -999
	imax = -999
	jmax = -999
	kmax = -999
	imax_grid=0
	fac_r_grid=-999.
	lim_r_grid=-999.
	dr_grid=0.		
	fac_y_grid=-999.
	lim_y_grid=-999.
	jmax_grid=0
	jmax_grid(1)=1
	dy_grid=0.
	sym_grid_y=0

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
	tstart_rms = 9.e18
	te_rms = 1.e19
	dt_max = -999.
	dt_ini = -999.
	time_int=''
	CFL = -999.
	n_dtavg = -1
	t0_output_movie = 9.e18
	dt_output_movie = 9.e18
	te_output_movie = 9.e18
	tstart_morf=0.
	tstart_morf2=0.
	trestart=0.
	!! num_scheme
	convection = 'ARGH'
	numdiff = 0.
	numdiff2 = 0.
	wiggle_detector = 0
	diffusion = 'ARGH'
	comp_filter_a = 0.5
	comp_filter_n = 0
	CNdiffz = 0
	CNdiff_factor = 0.5
	CNdiff_dtfactor = 100.
	CNdiff_ho = 0.
	CNdiff_pc = 0 !Choice in pre-conditioner of implicit CG 3D CN diffusion solver; 0 no preconditioner (default), 1 DIAG pc, 2 IC(0) pc
	CNdiff_maxi = 1000 !Maximum number of iterations in implicit CG 3D diffusion solver, default 1000
	CNdiff_tol = 1.e-12 !Tolerance in implicit CG 3D diffusion solver, default 1.e-12  
	CNdiff_ini = 1 !1 (default) = starting condition predictor U* just before applying implicit diff; 2 = starting condition U from previous timestep 
	npresIBM = 0
	npresIBM_viscupdate = 0
	npresPRHO = 0
	oPRHO = 2
	pres_in_predictor_step = 1
	advec_conc='VLE' !optional advection scheme concentration, options: 'VLE' (default) 'VL2' 'ARO' 'SBE' SB2' 'NVD' :'VLE' VanLeer(via LW)(default) 'ARO' Arora(via LW) 'SBE' Superbee(via LW) 'VL2' VanLeer(via CDS2) or 'SB2' Superbee(via CDS2) TVD schemes or 'NVD' for NVD scheme 
	continuity_solver = 1 !nerd option, default is 1 (drdt+drudx=0). Optional: 2 (neglect drdt), 3 (dudx almost 0 U-mix), 33 (dudx=0 U-mix), 34 (dudx=0 U-vol with proper U-mix)
	transporteq_fracs = 'volufrac' !nerd option default volume fractions, but as option also 'massfrac' can be used internally (input/output still is volume frac!)
	split_rho_cont='CDS'  ! optional 'VL2' or 'SB2' TVD scheme (via CDS2 without 1-cfl term) to split off rho from rho*U, default 'CDS' scheme
	driftfluxforce_calfac=1.
	depo_implicit=0
	depo_cbed_option=0 !2-3-2020 obsolete; user input not read anymore from num_scheme always option 0 is used in code 
	IBMorder=0 
	Poutflow=0 ! 0 (default, most robust) Poutflow is zero for complete outflow crosssection at rmax; 1 (optional) Poutflow is zero at just one grid-location at rmax --> sometimes better but less robust
	k_pzero=1 ! default 1, defines vertical level where p=0 when Poutflow=1 is used
	Uoutflow=0 !0 (default) Neumann outflow dUdn=0 (for U,V,W); 2 means Convective outflow condition dUdt+U_normal*dUdn=0 (for U,V,W)
	applyVOF=0
	k_ust_tau=1
	k_ust_tau_flow=1
	dUVdn_IBMbed=-1 !default no correction, but with 0 then dUdn and dVdn is made zero over immersed bed and with -2  to apply UV(1:kbed)=0 also for IBM2
	initPhydrostatic = 0 !default 0 no initial hydrostatic pressure; 1 = initialize hydrostatic pressure before first timestep 
	momentum_exchange_obstacles = 1 !default there is momentum exchange between flow and obstacles, optional 0 makes momenum terms zero in case the advective velocity is inside or at edge obstacle 
	!! ambient:
	U_b = -999.
	V_b = -999.
	W_b = -999.
	U_init = -999.
	V_init = -999.
	bcfile = ''
	U_b_tseriesfile=''
	V_b_tseriesfile=''
	W_b_tseriesfile=''
	U_b_tseries=-99999.
	U_b_series=-99999.	
	V_b_tseries=-99999.
	V_b_series=-99999.	
	W_b_tseries=-99999.
	W_b_series=-99999.		
	rho_b = -999.  
	rho_b2 = -999.  
	SEM = -999
	U_bSEM = -999.
	V_bSEM = -999.
	nmax2 = -999
	nmax1 = -999
	nmax3 = -999
	lm_min = -999.
	lm_min3 = -999.
	slip_bot = -999
	kn = -999.
	kn_flow_file=''
	kn_flow_d50_multiplier = -2.
	interaction_bed=-999
	periodicx=0
	periodicy=0
	monopile=-1
	kn_mp = -999.
	kn_sidewalls = -999.
	dpdx=0.
	dpdy=0.
	dpdx1=0.
	dpdy1=0.
	dpdx3=0.
	dpdy3=0.	
	cbc_perx_j(1:2)=0 
	cbc_relax = 1. 
	Uavold=0.
	Vavold=0.
	U3avold=0. 
	V3avold=0.
	W_ox=0.
	Hs=0. 
	Tp=999.
	nx_w=0.
	ny_w=0.
	U_w=-9999.
	V_w=-9999.	
	bedlevelfile = ''
	initconditionsfile = ''
	cfixedbed=0.6
	c_bed(:)=-1.
	DO i=1,4
		obst(:)%x(i) = -99999.
		obst(:)%y(i) = -99999. 
	ENDDO
	obst(:)%height = -99999. 
	obst(:)%zbottom = -99999.
	obst(:)%ero = 0.  !default no erosion on top of obstacle with interaction_bed=4
	obst(:)%depo = 0. !default no deposition on top of obstacle with interaction_bed=4
	bc_obst_h = 0.
	U_b3=-9999.
	V_b3=-9999.
	surf_layer=0.
	wallup=0
	obstfile = ''
	nobst_files = 0
	i_periodicx=0
	istart_morf2=0
	
	!! plume
	W_j = -999.
	plumetseriesfile=''
	plumetseries=-99999.
	plumeQtseries=-99999.
	plumectseries=-99999.
	Q_j = -999.
	plumeQtseriesfile=''
	plumectseriesfile=''
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
	slipvel = -999
	hindered_settling = 1 !! hindered_settling = 1	!Hindered settling formula [-] 1=Rowe (1987) (smooth Ri-Za org) (default); 2=Garside (1977); 3=Di Felice (1999)
	hindered_settling_c = 0 
	outflow_overflow_down = 0
	U_j2 = -999.
	plumetseriesfile2=''
	plumetseries2=-99999.
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
		bedplume(:)%x2(i) = -99999.
		bedplume(:)%y2(i) = -99999. 		
	ENDDO
	bedplume(:)%height = -99999. 
	bedplume(:)%u = -99999. 
	bedplume(:)%v = -99999. 
	bedplume(:)%w = -99999. 
	bedplume(:)%forever = 0
	bedplume(:)%t0 = 0.
	bedplume(:)%t_end = 9.e18
	bedplume(:)%Q = 0.
	bedplume(:)%volncells = 0.
	bedplume(:)%zbottom = -99999.
	bedplume(:)%changesedsuction = 1.	
	bedplume(:)%h_tseriesfile=''
	bedplume(:)%h_seriesloc=1	
	bedplume(:)%zb_tseriesfile=''
	bedplume(:)%zb_seriesloc=1
	bedplume(:)%Q_tseriesfile=''
	bedplume(:)%Q_seriesloc=1
	bedplume(:)%S_tseriesfile=''
	bedplume(:)%S_seriesloc=1
	bedplume(:)%c_tseriesfile=''
	bedplume(:)%c_seriesloc=1
	bedplume(:)%velocity_force=1 
	bedplume(:)%u_tseriesfile=''
	bedplume(:)%v_tseriesfile=''
	bedplume(:)%w_tseriesfile=''
	bedplume(:)%u_seriesloc=1
	bedplume(:)%v_seriesloc=1
	bedplume(:)%w_seriesloc=1
	bedplume(:)%dt_history=0.
	bedplume(:)%move_dz_height_factor=1. 
	bedplume(:)%move_dz_zbottom_factor=1.
	bedplume(:)%move_u=0.
	bedplume(:)%move_v=0.
	bedplume(:)%move_w=0.
	bedplume(:)%radius=-9.
	bedplume(:)%radius2=-9.
	
	
	DO i=1,30
		bedplume(:)%c(i) = 0.
		bedplume(:)%sedflux(i) = 0. 
	ENDDO
	DO i=1,100000
		bedplume(:)%move_dx_series(i)=-99999999.
		bedplume(:)%move_dy_series(i)=-99999999.
		bedplume(:)%move_dz_series(i)=-99999999.
		bedplume(:)%move_zbed_criterium(i)=999999999.
		bedplume(:)%move_zbed_type(i)=-9
		bedplume(:)%move_outputfile_series(i)=-1
		bedplume(:)%move_nx_series(i)=-99999999.
		bedplume(:)%move_ny_series(i)=-99999999.
		bedplume(:)%move_dx2_series(i)=-99999999.
		bedplume(:)%move_dy2_series(i)=-99999999.		
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
	fract(:)%type=1
	fract(:)%CD=0.
	fract(:)%cmax=100.
	nfrac=0
	nfr_silt=0
	nfr_sand=0
	nfr_air=0
	!! LESmodel
	sgs_model = 'ARGHH'
	Cs = -999.
	Cs_relax = 0.01 
	Lmix_type = -999
	nr_HPfilter = 0
	damping_drho_dz = 'none'
	damping_a1 = -999.
	damping_b1 = -999.
	damping_a2 = -999.
	damping_b2 = -999.
	extra_mix_visc='none'
	nu_minimum_wall = 0 
	! four constants default given value for default Realizible K-Epsilon turbulence model:
	Const1eps = 1.44 
	Const2 = 1.9 
	Sc_k = 1.0
	Sc_eps = 1.2
	Cal_buoyancy_k = 1.
	Cal_buoyancy_eps = 1.
	
	
	!!constants
	kappa=-999.
	gx = -999.
	gy = -999.
	gz = -999.
	settling_along_gvector = 0
	ekm_mol = -999.
	time_nm = 0.
	time_n=0.
	time_np=0.
	calibfac_sand_pickup = 1.
	calibfac_sand_bedload = 1.
	calibfac_Shields_cr = 1.
	pickup_formula = 'vanrijn1984' !default
	bedload_formula ='nonenon0000' !default no bedload taken into account 
	kn_d50_multiplier = 2. !default, kn=2*d50 defined in paper Van Rijn 1984
	kn_d50_multiplier_bl = 2. !default, kn=d90~2*d50 
	avalanche_slope = -99. 
	avalanche_slope(1)=0. !default vertical slopes are allowed
	av_slope_z= -1.
	avfile = ''
	reduction_sedimentation_shields = 0. ! default no reduction in sedimentation by shear stresss (shields) ! PhD thesis vRhee p146 eq 7.74
	morfac = 1.
	morfac2 = 1.
	avalanche_until_done=0
	pickup_correction='nonenonenonenonenone'
	vwal=-999.
	wbed_correction = 0
	nl=-999.
	fcor =1.
	permeability_kl=-999.
	pickup_fluctuations_ampl=0.
	pickup_fluctuations=0
	cbed_method = 1
	k_layer_pickup = 1
	wallmodel_tau_sed = 1
	ndtbed = 10
	nrmsbed = 100 
	z_tau_sed = -999.
	pickup_bedslope_geo=0
	bl_relax=0.01 
	bedslope_effect=0
	bedslope_mu_s=0.63 
	alfabs_bl=1.0 
	alfabn_bl=1.5 
	phi_sediment=30.
	erosion_cbed_start = 1 !
	movebed_absorb_cfluid = 1 !default 1 absorb cfluid in bed when bed moves and new bed is higher than previous timestep (like in avalanche and regular bed-update), 0 = add cfluid buried in new bed to first fluid cell above new bed
	power_VR2019 = 1.
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
	cutter='not'
	
	!! rheology
	Non_Newtonian = 0			!Default is 0, Newtonian treatment
	Rheological_model ='AAARGH'
	PAPANASTASIOUS_m=9.e18		!if not defined than inf --> model reduces to standard Bingham
	Apvisc_interp = 1 ! 1 (default) is linear interpolation neighbouring cells for apparent viscosity; 2 is maximum of neighbouring cells 
	shear0limit = 0. !1/s if not defined then zero and no influence 
	rheo_shear_method = 1
	Apvisc_shear_relax=1. !default no relaxation 
	Apvisc_force_eq = 0 !default no additional direct force terms from rheology 
	SIMPLE_tauy=0.2
	SIMPLE_muB=0.1
	SIMPLE_climit(:) = 0.
	JACOBS_Ky=6.72e4
	JACOBS_Kmu=251
	JACOBS_Aclay=1.0
	JACOBS_By=-4.75
	JACOBS_Bmu=-2.64
	JACOBS_muw=0.001
	WINTER_Ay=7.3e5
	WINTER_Amu=932
	WINTER_nf=2.64
	WINTER_af=3.65
	WINTER_muw=0.001
	THOMAS_Cy=7.45e5
	THOMAS_Cmu=1e-3
	THOMAS_ky=1.5
	THOMAS_kmu=1.25
	THOMAS_Py=5.61
	THOMAS_Pmu=-3.03
	THOMAS_phi_sand_max=0.6		!Maximal solids fraction of sand
	Lambda_init=0.
	Kin_eq_a=1.
	Kin_eq_b=0.02
	Kin_eq_lambda_0=1.
	HOUSKA_n=1.
	HOUSKA_eta_0=0.3
	HOUSKA_eta_inf=	20.
	HOUSKA_tauy_0=400.
	HOUSKA_tauy_inf=20.
	BAGNOLD_beta=0.275
	BAGNOLD_phi_max=0.6			!Usually 0.6	

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
	IF (dy<0.and.dy_grid(1).eq.0.) CALL writeerror(9)
	IF (depth<0.) CALL writeerror(10)
	IF (mod(jmax,px).ne.0) CALL writeerror(11)
	IF (mod(kmax,px).ne.0) CALL writeerror(12)
	IF (sym_grid_y.ne.0.and.sym_grid_y.ne.1.and.sym_grid_y.ne.-1) CALL writeerror(13)
	jmax=jmax/px
	READ (UNIT=1,NML=times,IOSTAT=ios)
	!! check input times  
	IF (t_end<0.) CALL writeerror(30)
	IF (t0_output<0.) CALL writeerror(31)
	IF (dt_output<0.) CALL writeerror(31)
	IF (te_output<0.) CALL writeerror(31)
	IF (tstart_rms<0.) CALL writeerror(32)
	IF (te_rms<tstart_rms) CALL writeerror(36)
	IF (dt_max<0.) CALL writeerror(33) 
	IF (dt_ini<0.) THEN 
	  dt_ini = dt_max
	ELSE
	  dt_ini = MIN(dt_max,dt_ini)
	ENDIF
	IF (time_int.ne.'EE1'.AND.time_int.ne.'RK3'.AND.time_int.ne.'AB2'.AND.time_int.ne.'AB3'.AND.time_int.ne.'ABv') 
     &     CALL writeerror(34) 	 
	IF (time_int.eq.'AB2'.or.time_int.eq.'AB3') THEN
		write(*,*),' WARNING: Your time integration scheme: ',time_int
		write(*,*),' is a testing option and not all functionalities of Dflow3d are working,'
		write(*,*),' use ABv for a fully supported time integration scheme.'
	ENDIF
	IF (CFL<0.) CALL writeerror(35)
	IF (n_dtavg>0) THEN 
		ALLOCATE(dt_series(n_dtavg))
	ENDIF 
	READ (UNIT=1,NML=num_scheme,IOSTAT=ios)
	!! check input num_scheme
	IF (convection.ne.'CDS2'.AND.convection.ne.'CDS6'.AND.convection.ne.'COM4'.AND.convection.ne.'CDS4'
     &  .AND.convection.ne.'HYB6'.AND.convection.ne.'HYB4'.AND.convection.ne.'C4A6'.AND.convection.ne.'uTVD'
     &  .AND.convection.ne.'C2Bl'.AND.convection.ne.'C4Bl') CALL writeerror(401)
	IF (numdiff<0.or.numdiff>1.) CALL writeerror(402)
	IF (wiggle_detector.ne.0.and.wiggle_detector.ne.1.and.wiggle_detector.ne.2) CALL writeerror(408) 
	wd = wiggle_detector 
	IF (diffusion.ne.'CDS2'.AND.diffusion.ne.'COM4') CALL writeerror(403) 
	IF (comp_filter_a<0.or.comp_filter_a>0.5) CALL writeerror(404)
	IF (comp_filter_n<0) CALL writeerror(405)
	if (numdiff>0.and.numdiff2<numdiff) then 
		numdiff2=numdiff2/numdiff !used as minimum for wf
	else 
		numdiff2=0.
	endif 
	if (convection.eq.'C2Bl'.or.convection.eq.'C4Bl') then 
	  numdiff = numdiff   !no hidden factor 
	else
	  numdiff=numdiff*2.  !needed to get correct value (in advec is a 'hidden' factor 2/4)
	endif
	IF (CNdiffz.ne.0.and.CNdiffz.ne.1.and.CNdiffz.ne.2.and.CNdiffz.ne.11.and.CNdiffz.ne.12.and.CNdiffz.ne.31) CALL writeerror(406)
	IF (npresIBM<0) CALL writeerror(407)
	pres_in_predictor_step_internal = pres_in_predictor_step
	IF (pres_in_predictor_step_internal.eq.3.or.pres_in_predictor_step_internal.eq.4) THEN !make Pold in bed zero every timestep; use pres_in_predictor_step=1 in the rest of the code
		pres_in_predictor_step = 1 
	ENDIF 
	IF (k_ust_tau<1.or.k_ust_tau>kmax) CALL writeerror(409)
	IF (k_ust_tau_flow<1.or.k_ust_tau_flow>kmax) CALL writeerror(409)
	IF (k_pzero<1.or.k_pzero>kmax) CALL writeerror(410)
	IF (CNdiff_dtfactor.lt.0.9999) CALL writeerror(411) 

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
	IF (rho_b2<0.) rho_b2=rho_b
	IF (SEM<0) CALL writeerror(44)
	IF (nmax2<0) CALL writeerror(45)
	IF (nmax1<0) CALL writeerror(46)
	IF (nmax3<0) CALL writeerror(620)
	IF (lm_min<0.) CALL writeerror(47)
	IF (lm_min3<0) CALL writeerror(621)
	IF (slip_bot<-2) CALL writeerror(48)
	IF (kn<0.) CALL writeerror(49)
	IF (interaction_bed<0) CALL writeerror(50)
	IF (periodicx.ne.0.and.periodicx.ne.1.and.periodicx.ne.2) CALL writeerror(52)
	IF (periodicy.ne.0.and.periodicy.ne.1.and.periodicy.ne.2) CALL writeerror(53)
	IF (cbc_perx_j(1).ne.0.and.cbc_perx_j(2).ne.0) THEN 
	  IF(cbc_perx_j(2).lt.cbc_perx_j(1).or.cbc_perx_j(1).lt.1.or.cbc_perx_j(2).gt.jmax) CALL writeerror(622)
	ENDIF 
	IF (cbc_relax.gt.1.or.cbc_relax.lt.0) CALL writeerror(623)
	!IF (periodicx.eq.1.and.dpdx.eq.0.) CALL writeerror(54)
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
	IF (U_init<-998.) U_init=U_b ! only if defined then different U_init is used else equal to U_b
	IF (V_init<-998.) V_init=V_b
	IF (kn_mp<0.and.monopile.eq.1) CALL writeerror(610)
	
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
	  ob(n)%ero=obst(n)%ero
	  ob(n)%depo=obst(n)%depo
	    IF (ob(n)%ero.lt.0.or.ob(n)%depo.lt.0.or.ob(n)%ero.gt.1.or.ob(n)%depo.gt.1.) THEN
			write(*,*),' Obstacle:',n
			CALL writeerror(1051)
	    ENDIF	  
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
	IF (wallup.ne.0.and.wallup.ne.1.and.wallup.ne.2) CALL writeerror(605)
	IF (cfixedbed<0.or.cfixedbed>1.) CALL writeerror(607)
	IF (i_periodicx<0.or.istart_morf2<0.or.i_periodicx>imax.or.istart_morf2>imax) CALL writeerror(613)
	DEALLOCATE(obst)

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
	IF (slipvel.ne.0.and.slipvel.ne.1.and.slipvel.ne.2) CALL writeerror(68)
	IF (hindered_settling.ne.1.and.hindered_settling.ne.2.and.hindered_settling.ne.3) CALL writeerror(280)
	IF (hindered_settling_c.ne.0.and.hindered_settling_c.ne.1) CALL writeerror(282)
	
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
	ALLOCATE(nfrac_air2(nfrac))
	ALLOCATE(nfrac_silt(nfrac))
	ALLOCATE(nfrac_sand(nfrac))
	  i=0
	  n1=0
	  n2=0
	  n3=0	  
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
	  frac(n)%type = fract(n)%type
	  frac(n)%CD=fract(n)%CD
	  frac(n)%cmax=fract(n)%cmax	  
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
	  IF (frac(n)%CD<0.) THEN
		write(*,*)'Fraction:',n
		CALL writeerror(175)
	  ENDIF	  
	  IF (frac(n)%type.eq.1) THEN
		n1=n1+1
		nfrac_silt(n1)=n
		IF (n1>1) THEN
			IF (frac(n)%dpart.le.frac(nfrac_silt(n1-1))%dpart) THEN	
				write(*,*)'Fraction:',n
				CALL writeerror(172)
			ENDIF				
		ENDIF
	  ELSEIF (frac(n)%type.eq.2) THEN
		n2=n2+1
		nfrac_sand(n2)=n
		IF (n2>1) THEN
			IF (frac(n)%dpart.le.frac(nfrac_sand(n2-1))%dpart) THEN	
				write(*,*)'Fraction:',n
				CALL writeerror(172)
			ENDIF				
		ENDIF		
	  ELSEIF (frac(n)%type.eq.3) THEN 
		n3=n3+1
		nfrac_air2(n3)=n
	  ELSEIF (frac(n)%type.eq.-1) THEN
		! do nothing
		applyVOF=1
	  ELSE 
		write(*,*),'fraction nr:',n
	    CALL writeerror(171)
	  ENDIF
	ENDDO
	  nair=i !nr of compressible air fractions
	  nfr_silt=n1
	  nfr_sand=n2
	  nfr_air=n3
	  
	  !check op c_bed after nfrac
	DO n=1,nfrac
		IF(c_bed(n)<0.or.c_bed(n)>1.) THEN
		write(*,*),'fraction :',n		
			CALL writeerror(608)
		ENDIF
	ENDDO
	!IF (SUM(c_bed(1:nfrac))<cfixedbed.or.SUM(c_bed(1:nfrac))>1.) CALL writeerror(609)	  
	  
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
	  bp(n)%uinput=bedplume(n)%u
	  bp(n)%v=bedplume(n)%v
	  bp(n)%w=bedplume(n)%w
	  bp(n)%c=bedplume(n)%c
	  bp(n)%forever=bedplume(n)%forever
	  bp(n)%velocity_force=bedplume(n)%velocity_force
	  bp(n)%t0=bedplume(n)%t0
	  bp(n)%t_end=bedplume(n)%t_end
	  bp(n)%zbottom=bedplume(n)%zbottom
	  bp(n)%Q=bedplume(n)%Q
	  bp(n)%sedflux=bedplume(n)%sedflux
	  bp(n)%volncells = bedplume(n)%volncells
	  bp(n)%changesedsuction=bedplume(n)%changesedsuction
	  bp(n)%h_tseriesfile=bedplume(n)%h_tseriesfile
	  bp(n)%zb_tseriesfile=bedplume(n)%zb_tseriesfile
	  bp(n)%h_seriesloc=bedplume(n)%h_seriesloc
	  bp(n)%zb_seriesloc=bedplume(n)%zb_seriesloc
	  bp(n)%u_tseriesfile=bedplume(n)%u_tseriesfile
	  bp(n)%v_tseriesfile=bedplume(n)%v_tseriesfile
	  bp(n)%w_tseriesfile=bedplume(n)%w_tseriesfile
	  bp(n)%u_seriesloc=bedplume(n)%u_seriesloc
	  bp(n)%v_seriesloc=bedplume(n)%v_seriesloc
	  bp(n)%w_seriesloc=bedplume(n)%w_seriesloc
	  bp(n)%dt_history=bedplume(n)%dt_history
	  bp(n)%move_dz_height_factor=bedplume(n)%move_dz_height_factor
	  bp(n)%move_dz_zbottom_factor=bedplume(n)%move_dz_zbottom_factor
	  
	  bp(n)%move_zbed_criterium=bedplume(n)%move_zbed_criterium
	  bp(n)%move_zbed_type=bedplume(n)%move_zbed_type
	  bp(n)%move_dx_series=bedplume(n)%move_dx_series
	  bp(n)%move_dy_series=bedplume(n)%move_dy_series
	  bp(n)%move_dz_series=bedplume(n)%move_dz_series
	  bp(n)%move_nx_series=bedplume(n)%move_nx_series
	  bp(n)%move_ny_series=bedplume(n)%move_ny_series
	  bp(n)%x2=bedplume(n)%x2
	  bp(n)%y2=bedplume(n)%y2	  
	  bp(n)%move_dx2_series=bedplume(n)%move_dx2_series
	  bp(n)%move_dy2_series=bedplume(n)%move_dy2_series	  
	  bp(n)%move_u=bedplume(n)%move_u
	  bp(n)%move_v=bedplume(n)%move_v
	  bp(n)%move_w=bedplume(n)%move_w
	  bp(n)%radius=bedplume(n)%radius
	  IF (bp(n)%radius>0.) THEN
		bp(n)%x(2:4)=0.
		bp(n)%y(2:4)=0.
	  ENDIF
	  
	  bp(n)%move_outputfile_series=bedplume(n)%move_outputfile_series
	  bp(n)%nmove_present=0
	  bp(n)%nmove=0
	  DO WHILE (bp(n)%move_dx_series(bp(n)%nmove+1).NE.-99999999.)
		bp(n)%nmove=bp(n)%nmove+1
	  END DO
	  n3=0
	  DO WHILE (bp(n)%move_dy_series(n3+1).NE.-99999999.)
		n3=n3+1
	  END DO
	  n2=0
	  DO WHILE (bp(n)%move_dz_series(n2+1).NE.-99999999.)
		n2=n2+1
	  END DO
	  n1=0
	  DO WHILE (bp(n)%move_zbed_criterium(n1+1).NE.999999999.)
		n1=n1+1
	  END DO
	  IF (n1.eq.1) THEN 
	    bp(n)%move_zbed_criterium(1:100000)=bp(n)%move_zbed_criterium(1)
		n1=n3
	  ENDIF 
	  IF (n3.ne.bp(n)%nmove.or.n2.ne.bp(n)%nmove.or.n1.ne.bp(n)%nmove) THEN
		write(*,*),' Bedplume : ',n
		CALL writeerror(173)
	  ENDIF
	  n2=0
	  DO WHILE (bp(n)%move_nx_series(n2+1).NE.-99999999.)
		n2=n2+1
	  END DO
	  IF (n2.eq.1) THEN 
	    bp(n)%move_nx_series(1:100000)=bp(n)%move_nx_series(1)
		n2=bp(n)%nmove
	  ENDIF 	  
	  n1=0
	  DO WHILE (bp(n)%move_ny_series(n1+1).NE.-99999999.)
		n1=n1+1
	  END DO
	  IF (n1.eq.1) THEN 
	    bp(n)%move_ny_series(1:100000)=bp(n)%move_ny_series(1)
		n1=bp(n)%nmove
	  ENDIF 
	  IF (n1.ne.bp(n)%nmove.or.n2.ne.bp(n)%nmove) THEN
		write(*,*),' Bedplume : ',n
		CALL writeerror(173)
	  ENDIF	  
	  
	  n1=0
	  DO WHILE (bp(n)%move_outputfile_series(n1+1).NE.-1)
		n1=n1+1
	  END DO
	  IF (n1.eq.0) THEN 
	    bp(n)%move_outputfile_series(1:100000)=1
		n1=n3
	  ENDIF 
	  IF (n1.ne.bp(n)%nmove) THEN
		write(*,*),' Bedplume : ',n
		CALL writeerror(174)
	  ENDIF
	  n3=0
	  DO WHILE (bp(n)%move_dx2_series(n3+1).NE.-99999999.)
		n3=n3+1
	  END DO
	  IF (n3.eq.0) THEN
	    bp(n)%move_dx2_series=bp(n)%move_dx_series
      ELSEIF (n3.eq.1) THEN 
	    bp(n)%move_dx2_series(1:100000)=bp(n)%move_dx2_series(1)
      ELSEIF (n3.ne.bp(n)%nmove) THEN
		write(*,*),' Bedplume : ',n
		write(*,*),'Error in move_dx2_series'
		CALL writeerror(174)	    
	  ENDIF 
	  n3=0
	  DO WHILE (bp(n)%move_dy2_series(n3+1).NE.-99999999.)
		n3=n3+1
	  END DO
	  IF (n3.eq.0) THEN
	    bp(n)%move_dy2_series=bp(n)%move_dy_series
      ELSEIF (n3.eq.1) THEN 
	    bp(n)%move_dy2_series(1:100000)=bp(n)%move_dy2_series(1)
      ELSEIF (n3.ne.bp(n)%nmove) THEN
		write(*,*),' Bedplume : ',n
		write(*,*),'Error in move_dy2_series'
		CALL writeerror(174)	    
	  ENDIF 
	  IF (bp(n)%x2(1).eq.-99999.or.bp(n)%x2(2).eq.-99999.or.bp(n)%x2(3).eq.-99999.or.bp(n)%x2(4).eq.-99999.) THEN 
	    bp(n)%x2=bp(n)%x
	  ENDIF
	  IF (bp(n)%y2(1).eq.-99999.or.bp(n)%y2(2).eq.-99999.or.bp(n)%y2(3).eq.-99999.or.bp(n)%y2(4).eq.-99999.) THEN 
	    bp(n)%y2=bp(n)%y
	  ENDIF	 
	  IF (bp(n)%radius2<0.) THEN 
		bp(n)%radius2=bp(n)%radius
	  ENDIF 
	  n4=0
	  DO WHILE (bp(n)%move_zbed_type(n4+1).NE.-9)
		n4=n4+1
	  END DO
	  IF (n4.eq.1) THEN 
	    bp(n)%move_zbed_type(1:100000)=bp(n)%move_zbed_type(1)
	  ENDIF 	  
	  
	  DO i=1,10000
	  	bp(n)%h_tseries(i)=-99999.
	  	bp(n)%h_series(i)=-99999.
	  	bp(n)%zb_tseries(i)=-99999.
	  	bp(n)%zb_series(i)=-99999.
	  	bp(n)%Q_tseries(i)=-99999.
	  	bp(n)%Q_series(i)=-99999.
	  	bp(n)%S_tseries(i)=-99999.
	  	bp(n)%c_tseries(i)=-99999.
	  	bp(n)%u_tseries(i)=-99999.
	  	bp(n)%u_series(i)=-99999.
	  	bp(n)%v_tseries(i)=-99999.
	  	bp(n)%v_series(i)=-99999.		
	  	bp(n)%w_tseries(i)=-99999.
	  	bp(n)%w_series(i)=-99999.		
	  	DO j=1,30
	  		bp(n)%S_series(j,i)=-99999.
	  		bp(n)%c_series(j,i)=-99999.		
	  	ENDDO
	  ENDDO
	
	  
	  bp(n)%Q_tseriesfile=bedplume(n)%Q_tseriesfile
	  bp(n)%S_tseriesfile=bedplume(n)%S_tseriesfile
	  bp(n)%c_tseriesfile=bedplume(n)%c_tseriesfile
	  bp(n)%u_tseriesfile=bedplume(n)%u_tseriesfile
	  bp(n)%v_tseriesfile=bedplume(n)%v_tseriesfile
	  bp(n)%w_tseriesfile=bedplume(n)%w_tseriesfile
	  bp(n)%Q_seriesloc=bedplume(n)%Q_seriesloc
	  bp(n)%S_seriesloc=bedplume(n)%S_seriesloc
	  bp(n)%c_seriesloc=bedplume(n)%c_seriesloc


!	  bp(n)%h_tseries=bedplume(n)%h_tseries
!	  bp(n)%zb_tseries=bedplume(n)%zb_tseries
!	  bp(n)%h_series=bedplume(n)%h_series
!	  bp(n)%zb_series=bedplume(n)%zb_series
!	  
!	  bp(n)%Q_tseries=bedplume(n)%Q_tseries
!	  bp(n)%S_tseries=bedplume(n)%S_tseries
!	  bp(n)%c_tseries=bedplume(n)%c_tseries
!	  bp(n)%Q_series=bedplume(n)%Q_series
!	  bp(n)%S_series=bedplume(n)%S_series
!	  bp(n)%c_series=bedplume(n)%c_series
	  


	IF (bp(n)%Q_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries2(bp(n)%Q_tseriesfile,bp(n)%Q_tseries,bp(n)%Q_series)
	   bp(n)%Q_seriesloc=1
	   n3=0
	   DO WHILE (bp(n)%Q_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (bp(n)%Q_tseries(n3).lt.t_end) THEN
				write(*,*),' Bedplume : ',n
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',bp(n)%Q_tseriesfile
				CALL writeerror(279)
		ENDIF	   
	ENDIF	  
	IF (bp(n)%u_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries2(bp(n)%u_tseriesfile,bp(n)%u_tseries,bp(n)%u_series)
	   bp(n)%u_seriesloc=1
	   n3=0
	   DO WHILE (bp(n)%u_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (bp(n)%u_tseries(n3).lt.t_end) THEN
				write(*,*),' Bedplume : ',n
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',bp(n)%u_tseriesfile
				CALL writeerror(279)
		ENDIF	   
	ENDIF
	IF (bp(n)%v_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries2(bp(n)%v_tseriesfile,bp(n)%v_tseries,bp(n)%v_series)
	   bp(n)%v_seriesloc=1
	   n3=0
	   DO WHILE (bp(n)%v_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (bp(n)%v_tseries(n3).lt.t_end) THEN
				write(*,*),' Bedplume : ',n
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',bp(n)%v_tseriesfile
				CALL writeerror(279)
		ENDIF	   
	ENDIF
	IF (bp(n)%w_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries2(bp(n)%w_tseriesfile,bp(n)%w_tseries,bp(n)%w_series)
	   bp(n)%u_seriesloc=1
	   n3=0
	   DO WHILE (bp(n)%w_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (bp(n)%w_tseries(n3).lt.t_end) THEN
				write(*,*),' Bedplume : ',n
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',bp(n)%w_tseriesfile
				CALL writeerror(279)
		ENDIF	   
	ENDIF
	
	IF (bp(n)%S_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries3(bp(n)%S_tseriesfile,bp(n)%S_tseries,bp(n)%S_series(1:nfrac,:),nfrac)
	   bp(n)%S_seriesloc=1
	   n3=0
	   DO WHILE (bp(n)%S_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (bp(n)%S_tseries(n3).lt.t_end) THEN
				write(*,*),' Bedplume : ',n
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',bp(n)%S_tseriesfile
				CALL writeerror(279)
		ENDIF
	   
	   DO i=1,nfrac
		 DO j=1,n3
			IF (bp(n)%S_series(i,j)<0) THEN
				write(*,*),' Bedplume : ',n
				write(*,*),' Sediment fraction,# in serie:',i,j
				write(*,*),' value:',bp(n)%S_series(i,j)
				write(*,*),' series file:',bp(n)%S_tseriesfile
				CALL writeerror(271)
		    ENDIF
		  ENDDO
		ENDDO	   
	ENDIF
	IF (bp(n)%c_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries3(bp(n)%c_tseriesfile,bp(n)%c_tseries,bp(n)%c_series(1:nfrac,:),nfrac)
	   bp(n)%c_seriesloc=1
	   n3=0
	   DO WHILE (bp(n)%c_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (bp(n)%c_tseries(n3).lt.t_end) THEN
				write(*,*),' Bedplume : ',n
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',bp(n)%c_tseriesfile
				CALL writeerror(279)
		ENDIF	   
	   !write(*,*) 'bp_cseries 8', bp(n)%c_series(8,1:10)
	   DO i=1,nfrac
		 DO j=1,n3
			IF (bp(n)%c_series(i,j)<0) THEN
				write(*,*),' Bedplume : ',n
				write(*,*),' Sediment fraction, # n series:',i,j
				write(*,*),' value:',bp(n)%c_series(i,j)
				write(*,*),' series file:',bp(n)%c_tseriesfile
				CALL writeerror(271)
		    ENDIF
		  ENDDO
		ENDDO
	ENDIF	
	IF (bp(n)%c_tseriesfile.ne.''.and.bp(n)%S_tseriesfile.ne.'') THEN
		write(*,*),' Bedplume : ',n
		write(*,*),' series file:',bp(n)%c_tseriesfile
		write(*,*),' series file:',bp(n)%S_tseriesfile
		write(*,*),' should not be both defined'
		CALL writeerror(277)
	ENDIF
	IF (bp(n)%c_tseriesfile.ne.''.and.MAXVAL(bp(n)%sedflux(1:nfrac))>0.) THEN
		write(*,*),' Bedplume : ',n
		write(*,*),' series file:',bp(n)%c_tseriesfile
		write(*,*),' sedflux:',bp(n)%sedflux
		write(*,*),' should not be both defined'
		CALL writeerror(277)
	ENDIF
	IF (MAXVAL(bp(n)%c(1:nfrac))>0.and.bp(n)%S_tseriesfile.ne.'') THEN
		write(*,*),' Bedplume : ',n
		write(*,*),' series file:',bp(n)%S_tseriesfile
		write(*,*),' c:',bp(n)%c
		write(*,*),' should not be both defined'
		CALL writeerror(277)
	ENDIF
		
	IF (bp(n)%h_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries2(bp(n)%h_tseriesfile,bp(n)%h_tseries,bp(n)%h_series)
	   bp(n)%h_seriesloc=1
	   n3=0
	   DO WHILE (bp(n)%h_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (bp(n)%h_tseries(n3).lt.t_end) THEN
				write(*,*),' Bedplume : ',n
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',bp(n)%h_tseriesfile
				CALL writeerror(279)
		ENDIF	   
	ENDIF	  
	IF (bp(n)%zb_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries2(bp(n)%zb_tseriesfile,bp(n)%zb_tseries,bp(n)%zb_series)
	   bp(n)%zb_seriesloc=1
	   n3=0
	   DO WHILE (bp(n)%zb_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (bp(n)%zb_tseries(n3).lt.t_end) THEN
				write(*,*),' Bedplume : ',n
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',bp(n)%zb_tseriesfile
				CALL writeerror(279)
		ENDIF	   
	ENDIF
	  DO i=1,4
	    IF (bp(n)%x(i).eq.-99999.or.bp(n)%y(i).eq.-99999) THEN
		write(*,*),' Bedplume : ',n
		CALL writeerror(269)
	    ENDIF
	  ENDDO
	  IF (bp(n)%height<0.) THEN
	    write(*,*),' Bedplume : ',n
	    CALL writeerror(270)
	  ENDIF
	  DO i=1,nfrac
	    IF (bp(n)%c(i)<0.or.bp(n)%sedflux(i)<0.) THEN
	      write(*,*),' Bedplume:',n
	      write(*,*),' Sediment fraction:',i
	      CALL writeerror(271)
	    ENDIF
	  ENDDO
	  DO i=1,nfrac
	    IF (bp(n)%c(i)>0.and.bp(n)%sedflux(i)>0.) THEN
	      write(*,*),' Bedplume:',n
	      write(*,*),' Sediment fraction:',i
	      CALL writeerror(277)
	    ENDIF
	  ENDDO
	  IF (bp(n)%forever.ne.0.and.bp(n)%forever.ne.1) THEN
	    write(*,*),' Bedplume:',n
	    CALL writeerror(272)
	  ENDIF
	  IF (bp(n)%t0.lt.0.or.bp(n)%t_end.lt.0.) THEN
	    write(*,*),' Bedplume:',n
	    CALL writeerror(273)
	  ENDIF
	  IF (bp(n)%zbottom>bp(n)%height.or.bp(n)%zbottom>depth) THEN
	    write(*,*),' Bedplume:',n
	    CALL writeerror(274)
	  ENDIF
	  IF (bp(n)%changesedsuction.lt.0.or.bp(n)%changesedsuction.gt.1.) THEN
	    write(*,*),' Bedplume:',n
	    CALL writeerror(278)
	  ENDIF
	ENDDO

	DEALLOCATE(bedplume)
	
	IF (MAXVAL(bp(1:nbedplume)%dt_history)>0) THEN
	  ALLOCATE(thisbp(nbedplume,20000))
	  ALLOCATE(zhisbp(nbedplume,20000))
	  ALLOCATE(Chisbp(nfrac,nbedplume,20000))
	  ALLOCATE(Uhisbp(nbedplume,20000))
	  ALLOCATE(Vhisbp(nbedplume,20000))
	  ALLOCATE(Whisbp(nbedplume,20000))
	  DO n=1,nbedplume 
		bp(n)%istep_bphis_output=0
		bp(n)%t_bphis_output=bp(n)%t0
	  ENDDO
	  thisbp=0.
	  zhisbp=0.
	  Chisbp=0.
	  Uhisbp=0.
	  Vhisbp=0.
	  Whisbp=0.
	ENDIF
	
	  
	READ (UNIT=1,NML=LESmodel,IOSTAT=ios)
	!! check input LESmodel
	IF (sgs_model.ne.'SSmag'.and.sgs_model.ne.'FSmag'.and.sgs_model.ne.'SWALE'.and.
     &  sgs_model.ne.'Sigma'.and.sgs_model.ne.'MixLe'.and.sgs_model.ne.'DSmag'.and.sgs_model.ne.'ReaKE') CALL writeerror(82)
	IF (Cs<0.) CALL writeerror(80)
	IF (Lmix_type<0) CALL writeerror(81)
	IF (nr_HPfilter<0) CALL writeerror(83)
	IF (damping_drho_dz.ne.'none'.and.damping_drho_dz.ne.'MuAn') CALL writeerror(85)
	IF ((damping_drho_dz.eq.'MuAn'.and.damping_a1.lt.0.).or.(damping_drho_dz.eq.'MuAn'.and.damping_a2.lt.0.)) 
     &  CALL writeerror(86)
	IF ((damping_drho_dz.eq.'MuAn'.and.damping_b1.lt.-990.).or.(damping_drho_dz.eq.'MuAn'.and.damping_b2.lt.-990.)) 
     &  CALL writeerror(87)
	IF (extra_mix_visc.ne.'none'.and.extra_mix_visc.ne.'Krie') CALL writeerror(88) 
	IF (Cs_relax>1.or.Cs_relax<0.) CALL writeerror(089)

	READ (UNIT=1,NML=constants,IOSTAT=ios)
	!! check input constants
	IF (kappa<0.) CALL writeerror(90)
	IF (gx<-800.) CALL writeerror(91)
	IF (gy<-800.) CALL writeerror(92)
	IF (gz<-80000000.) CALL writeerror(93)
	IF (ekm_mol<0.) CALL writeerror(94)
	IF (pickup_formula.ne.'vanrijn1984'.and.pickup_formula.ne.'nielsen1992'.and.pickup_formula.ne.'okayasu2010'
     & .and.pickup_formula.ne.'vanrijn2019'.and.pickup_formula.ne.'VR2019_Cbed'
     & .and.pickup_formula.ne.'VR1984_Cbed')   CALL writeerror(95)
	IF (kn_d50_multiplier<0.or.kn_d50_multiplier_bl<0.) CALL writeerror(96)
	IF (MINVAL(avalanche_slope)<0.and.MINVAL(avalanche_slope).ne.-99.) CALL writeerror(97)
	IF (calibfac_sand_pickup<0.or.calibfac_sand_bedload<0.) CALL writeerror(98)
	IF (calibfac_Shields_cr<0.) CALL writeerror(99)
	IF (morfac<0.) CALL writeerror(101)
	IF (morfac2<1.) CALL writeerror(102)
	IF (pickup_correction.eq.'MastBergenvdBerg2003'.or.pickup_correction.eq.'MBvdBerg2003_vecheck'
     &	.or.pickup_correction.eq.'breach_via_avalanche') THEN 
		vwal2=vwal 
		IF (vwal<0) CALL writeerror(103)
		IF (nl<0.or.nl<(1.-cfixedbed)) CALL writeerror(104)
		delta_nsed=(nl-(1.-cfixedbed))/(cfixedbed) !delta_nsed=(nl-n0)/(1-n0)
		IF (permeability_kl<0) CALL writeerror(105)
		IF (fcor<0.) CALL writeerror(108)
	ENDIF
	IF (k_layer_pickup<1) CALL writeerror(106)
	IF (bl_relax>1.or.bl_relax<0.) CALL writeerror(107)
	IF (bedload_formula.ne.'vanrijn2007'.and.bedload_formula.ne.'vanrijn2003'.and.bedload_formula.ne.'MeyPeMu1947'
     & .and.bedload_formula.ne.'nonenon0000')   CALL writeerror(109)
	IF (bedslope_effect.ne.0.and.bedslope_effect.ne.1.and.bedslope_effect.ne.2.and.bedslope_effect.ne.3)
     &	CALL writeerror(143) 
	IF (phi_sediment.lt.(4.0 * atan(1.0))) CALL writeerror(144)
	phi_sediment = phi_sediment/180.*4.0 * atan(1.0) !change from degrees to rad 
	IF (wallmodel_tau_sed.ne.1.and.wallmodel_tau_sed.ne.3.and.wallmodel_tau_sed.ne.4.and.wallmodel_tau_sed.ne.5.
     &	and.wallmodel_tau_sed.ne.11) THEN 
		CALL writeerror(145)
	IF (power_VR2019<0.) CALL writeerror(146)
	ENDIF 
	 
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

	READ (UNIT=1,NML=rheology,IOSTAT=ios)
	!! check input rheology
	IF (Non_Newtonian.ne.0.and.Non_Newtonian.ne.1.and.Non_Newtonian.ne.2) CALL writeerror(801)
	IF (Rheological_model.ne.'SIMPLE'.AND.Rheological_model.ne.'JACOBS'.AND.Rheological_model.ne.'WINTER'
     &  .AND.Rheological_model.ne.'THOMAS') CALL writeerror(802)
	 IF (Apvisc_shear_relax>1.or.Apvisc_shear_relax<0.) CALL writeerror(803)

	CLOSE(1)

	IF (plumetseriesfile.eq.'') THEN
	ELSE
	   call readtseries(plumetseriesfile,plumetseries,plumeUseries)
	   plumeseriesloc=1
	   n3=0
	   DO WHILE (plumetseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (plumetseries(n3).lt.t_end) THEN
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',plumetseriesfile
				CALL writeerror(279)
		ENDIF	   
	ENDIF
	IF (plumeQtseriesfile.eq.'') THEN
	ELSE
	   call readtseries(plumeQtseriesfile,plumeQtseries,plumeQseries)
	   plumeQseriesloc=1
	   n3=0
	   DO WHILE (plumeQtseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (plumeQtseries(n3).lt.t_end) THEN
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',plumetseriesfile
				CALL writeerror(279)
		ENDIF	   
	ENDIF
	IF (plumectseriesfile.eq.'') THEN
	ELSE
	   call readtseries3(plumectseriesfile,plumectseries,plumecseries(1:nfrac,:),nfrac)
	   plumecseriesloc=1
	   n3=0
	   DO WHILE (plumectseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (plumectseries(n3).lt.t_end) THEN
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',plumectseriesfile
				CALL writeerror(279)
		ENDIF	   
	   DO i=1,nfrac
		 DO j=1,n3
			IF (plumecseries(i,j)<0) THEN
				write(*,*),' Sediment fraction, plume inflow c # n series:',i,j
				write(*,*),' value:',plumecseries(i,j)
				write(*,*),' series file:',plumectseriesfile
				CALL writeerror(271)
		    ENDIF
		  ENDDO
		ENDDO
	ENDIF	
	
	IF (plumetseriesfile2.eq.'') THEN
	ELSE
	   call readtseries(plumetseriesfile2,plumetseries2,plumeUseries2)
	   plumeseriesloc2=1
	   n3=0
	   DO WHILE (plumetseries2(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (plumetseries2(n3).lt.t_end) THEN
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',plumetseriesfile2
				CALL writeerror(279)
		ENDIF	   
	ENDIF
	
	IF (U_b_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries(U_b_tseriesfile,U_b_tseries,U_b_series)
	   U_b_seriesloc=1
	   n3=0
	   DO WHILE (U_b_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (U_b_tseries(n3).lt.t_end) THEN
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',U_b_tseriesfile
				CALL writeerror(039)
		ENDIF	   
	ENDIF	
	IF (V_b_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries(V_b_tseriesfile,V_b_tseries,V_b_series)
	   V_b_seriesloc=1
	   n3=0
	   DO WHILE (V_b_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (V_b_tseries(n3).lt.t_end) THEN
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',V_b_tseriesfile
				CALL writeerror(039)
		ENDIF	   
	ENDIF
	IF (W_b_tseriesfile.eq.'') THEN
	ELSE
	   call readtseries(W_b_tseriesfile,W_b_tseries,W_b_series)
	   W_b_seriesloc=1
	   n3=0
	   DO WHILE (W_b_tseries(n3+1).NE.-99999.)
		n3=n3+1
	   END DO	  
	   IF (W_b_tseries(n3).lt.t_end) THEN
				write(*,*),' time series shorter than t_end'
				write(*,*),' series file:',W_b_tseriesfile
				CALL writeerror(039)
		ENDIF	   
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

	ALLOCATE(Ubot_TSHD(0:j1))
	ALLOCATE(Vbot_TSHD(0:j1))
	IF (LOA>0) THEN
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
	ENDIF


	Ubot_TSHD=0.
	Vbot_TSHD=0.


	ALLOCATE(cos_u(0:j1))
	ALLOCATE(cos_v(0:j1))
	ALLOCATE(sin_u(0:j1))
	ALLOCATE(sin_v(0:j1))
	ALLOCATE(cos_ut(0:jmax*px+1))
	ALLOCATE(sin_ut(0:jmax*px+1))
	ALLOCATE(cos_vt(0:jmax*px+1))
	ALLOCATE(sin_vt(0:jmax*px+1))	
	ALLOCATE(Ru(0:i1))
	ALLOCATE(Rp(0:i1))	
	ALLOCATE(dr(0:i1))
	ALLOCATE(phivt(0:jmax*px+1))
	ALLOCATE(phipt(0:jmax*px+1))
	ALLOCATE(dphi2t(0:jmax*px+1))
	ALLOCATE(phiv(0:j1))
	ALLOCATE(phip(0:j1))
	ALLOCATE(dphi2(0:j1))	
	ALLOCATE(Lmix2(1:imax,1:jmax))
	ALLOCATE(Lmix2hat(1:imax,1:jmax))
	ALLOCATE(vol_V(1:imax,1:jmax*px))
	ALLOCATE(vol_Vp(0:i1,0:j1))
	ALLOCATE(kbed(0:i1,0:j1))
	ALLOCATE(kbed0(0:i1,0:j1))
	ALLOCATE(kbedin(0:j1))
	ALLOCATE(kbedold(0:i1,0:j1))
!!!	ALLOCATE(kbed2(0:i1,0:px*jmax+1))
	ALLOCATE(kbedt(0:i1,0:j1))
!!!	ALLOCATE(kbed22(0:i1,0:j1))
	ALLOCATE(zbed(0:i1,0:j1))
	ALLOCATE(rhocorr_air_z(1:nfrac,0:k1))
	ALLOCATE(av_slope(1:imax,1:jmax,0:k1))
	ALLOCATE(wscorr_z(1:nfrac,0:k1))
	
	IF (SEM.eq.1) THEN
	IF (nmax1.gt.0.or.nmax2.gt.0) THEN
	  IF (rank.eq.0.or.rank.eq.px-1) THEN
		ALLOCATE(llist1(0:i1,1:kmax,1:2000))
	  ENDIF
	  ALLOCATE(llist2(0:j1,1:kmax,1:1000))
	  ALLOCATE(llmax1(0:i1,1:kmax))
	  ALLOCATE(llmax2(0:j1,1:kmax))	  	  
	ENDIF	
	ALLOCATE(AA1(3,3,0:i1,1:kmax))
	ALLOCATE(AA2(3,3,0:j1,1:kmax))
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
	IF (nmax3.gt.0) THEN
	  ALLOCATE(llist3(0:i1,0:j1,1:1000))
	  ALLOCATE(llmax3(0:i1,0:j1))
	ENDIF
	ALLOCATE(AA3(3,3,0:i1,0:j1))	
	ALLOCATE(rSEM3(nmax3))
	ALLOCATE(thetaSEM3(nmax3))
	ALLOCATE(zSEM3(nmax3))
	ALLOCATE(xSEM3(nmax3))
	ALLOCATE(ySEM3(nmax3))
	ALLOCATE(epsSEM3(3,nmax3))
	ALLOCATE(wSEM3(nmax3))
	ALLOCATE(lmrSEM3(nmax3))
	ALLOCATE(lmzSEM3(nmax3))	
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
	
	
	ALLOCATE(rhU(0:i1,0:j1,0:k1))
	ALLOCATE(rhV(0:i1,0:j1,0:k1))
	ALLOCATE(rhW(0:i1,0:j1,0:k1))
	ALLOCATE(cU(nfrac,0:i1,0:j1,0:k1))
	ALLOCATE(cV(nfrac,0:i1,0:j1,0:k1))
	ALLOCATE(cW(nfrac,0:i1,0:j1,0:k1))
	ALLOCATE(fc_global(0:i1,0:jmax*px+1,0:k1))
	
	if (transporteq_fracs.eq.'massfrac') then
		ALLOCATE(rhoU(0:i1,0:j1,0:k1))  
		ALLOCATE(rhoV(0:i1,0:j1,0:k1))  
		ALLOCATE(rhoW(0:i1,0:j1,0:k1))  
	endif
	if (sgs_model.eq.'ReaKE') then
		ALLOCATE(TKE(0:i1,0:j1,0:k1))  
		ALLOCATE(EEE(0:i1,0:j1,0:k1))  
		ALLOCATE(Cmu(1:imax,1:jmax,1:kmax)) 
		TKE=1.e-12
		EEE=1.e-12
	endif

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
	IF (interaction_bed.ge.4) THEN
	  ALLOCATE(bednotfixed(0:i1,0:j1,0:k1))	
	  ALLOCATE(bednotfixed_depo(0:i1,0:j1,0:k1))
	  bednotfixed=1. !default avalanche or erosion is allowed everywhere, only in obstacles connected to bed not allowed, see init.f
	  bednotfixed_depo=1. !default deposition is allowed everywhere, only in obstacles connected to bed not allowed, see init.f
	  ALLOCATE(d_cbotdelay(nfrac,0:i1,0:j1))
	  d_cbotdelay = 0.
	ENDIF
	IF (interaction_bed.ge.4) THEN
	  ALLOCATE(Clivebed(nfrac,0:i1,0:j1,0:k1))
	  Clivebed=0.
	ENDIF
			Coldbot=0.
			Cnewbot=0.
			dCdtbot=0.
	ALLOCATE(b_update(0:i1))
	IF (tstart_morf2.gt.1e-6) THEN 
		b_update=0. 
	ELSE 
		b_update(istart_morf2:i1)=1.
	ENDIF 	
	if (cbc_perx_j(1)>0) then  !use quasi periodic bc for concentration:
	  b_update(0:1)=0. ! no bed-update at inflow
	  b_update(imax:i1)=0. ! no bed-update at outflow
	endif 
	
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
	ALLOCATE(pold1(1:imax,1:jmax,1:kmax))
	ALLOCATE(pold2(1:imax,1:jmax,1:kmax))
	ALLOCATE(pold3(1:imax,1:jmax,1:kmax))
	ALLOCATE(phnew(1:imax,1:jmax,1:kmax))
	ALLOCATE(phdt(1:imax,1:jmax,1:kmax))
	ALLOCATE(dp(0:i1,0:j1,0:k1))
	ALLOCATE(viscf(0:i1,0:j1,0:k1))

	if (sgs_model.eq.'DSmag') then
	  ALLOCATE(Csgrid(1:imax,1:jmax,1:kmax))
	  Csgrid=0.1 !start with default Cs=0.1
	endif
	p=0.
	pold=0.
	pold1=0.
	pold2=0.
	pold3=0.
	phdt=0.
	phnew=0.
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
		
		wxolder=0.
		wyolder=0.
		wzolder=0.
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
	ALLOCATE(Ub3new(0:i1,0:j1))	
	ALLOCATE(Vb3new(0:i1,0:j1))	
	ALLOCATE(Wb3new(0:i1,0:j1))	
	ALLOCATE(Ub3old(0:i1,0:j1))	
	ALLOCATE(Vb3old(0:i1,0:j1))	
	ALLOCATE(Wb3old(0:i1,0:j1))	
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
	ALLOCATE(ubot(0:i1,0:j1))
	ALLOCATE(Wbed(0:i1,0:j1))
	Wbed = 0.
	
	ALLOCATE(qbU(1:nfrac,0:i1,0:j1))
	ALLOCATE(qbV(1:nfrac,0:i1,0:j1))
	ALLOCATE(uuR_relax(0:i1,0:j1))
	ALLOCATE(uuL_relax(0:i1,0:j1))
	ALLOCATE(vvR_relax(0:i1,0:j1))
	ALLOCATE(vvL_relax(0:i1,0:j1))	
	ALLOCATE(qb_relax(0:i1,0:j1))
	uuR_relax = 0. 
	uuL_relax = 0. 
	vvR_relax = 0. 
	vvL_relax = 0. 
	qb_relax = 0. 
	qbU = 0.
	qbV = 0.
	
	if (monopile.eq.3) then 
		ALLOCATE(tau2Vold(0:j1,0:k1))
		ALLOCATE(tau2Wold(0:j1,0:k1))
		ALLOCATE(tau2Vnew(0:j1,0:k1))
		ALLOCATE(tau2Wnew(0:j1,0:k1))
		tau2Vold = 0.
		tau2Wold = 0.
		tau2Vnew = 0.
		tau2Wnew = 0.
	endif 
	ALLOCATE(tau_fl_Vold(0:i1,0:j1))
	ALLOCATE(tau_fl_Uold(0:i1,0:j1))
	ALLOCATE(tau_fl_Vnew(0:i1,0:j1))
	ALLOCATE(tau_fl_Unew(0:i1,0:j1))
	tau_fl_Vold = 0.
	tau_fl_Uold = 0.
	tau_fl_Vnew = 0.
	tau_fl_Unew = 0.
	ALLOCATE(ust_sl_new(0:i1,0:j1))
	ALLOCATE(ust_sl_old(0:i1,0:j1))
	ALLOCATE(ust_bl_new(0:i1,0:j1))
	ALLOCATE(ust_bl_old(0:i1,0:j1))
	ALLOCATE(ust_mud_new(0:i1,0:j1))
	ALLOCATE(ust_mud_old(0:i1,0:j1))
	ALLOCATE(ust_frac_new(1:nfrac,0:i1,0:j1))
	ALLOCATE(ust_frac_old(1:nfrac,0:i1,0:j1))	
	ust_sl_new=0.
	ust_sl_old=0.
	ust_bl_new=0.
	ust_bl_old=0.
	ust_mud_new=0.
	ust_mud_old=0.
	ust_frac_new=0.
	ust_frac_old=0.
	ALLOCATE(kn_flow(0:i1,0:j1))	
	kn_flow(0:i1,0:j1) = kn !default use kn defined in input file; in main.f subroutine determine_kn_flow is called once when kn_flow_file is defined or every timestep when kn_flow_d50_multiplies is defined
	
	!new variables rheology
	IF (Non_Newtonian.eq.1.or.Non_Newtonian.eq.2) THEN
		ALLOCATE(tauY(0:i1,0:j1,0:k1))
		ALLOCATE(muB(0:i1,0:j1,0:k1))
		ALLOCATE(stress(0:i1,0:j1,0:k1))
		ALLOCATE(strain(0:i1,0:j1,0:k1))
		ALLOCATE(muA(0:i1,0:j1,0:k1))
	ENDIF
	IF (Non_Newtonian.eq.2) THEN
		ALLOCATE(lambda_old(0:i1,0:j1,0:k1))
		ALLOCATE(lambda_new(0:i1,0:j1,0:k1))
	ENDIF	
	
	IF (wallmodel_tau_sed.eq.11) THEN 
		ALLOCATE(sigUWbed(1:nrmsbed,0:i1,0:j1))
		ALLOCATE(sigVWbed(1:nrmsbed,0:i1,0:j1))
		ALLOCATE(sigUbed(1:nrmsbed,0:i1,0:j1))
		ALLOCATE(sigVbed(1:nrmsbed,0:i1,0:j1))
		ALLOCATE(sigWbed(1:nrmsbed,0:i1,0:j1))
		ALLOCATE(sigtbed(1:nrmsbed))
		sigUWbed=0.
		sigVWbed=0.
		sigUbed=0.
		sigVbed=0.
		sigWbed=0.
		sigtbed=1e-6 !start with not zero, otherwise maybe divide by zero
	ENDIF
	
	! init at zero (when no bcfile it stays zero)
	Ubcoarse1=0.
	Ubcoarse2=0.
	Vbcoarse1=0.
	Vbcoarse2=0.
	Wbcoarse1=0.
	Wbcoarse2=0.
	Cbcoarse1=0.
	Cbcoarse2=0.
	ubot=0.

	ALLOCATE(Xii(1:kmax))
	ALLOCATE(Tkk(1:jmax*px))
	ALLOCATE(Xkk(kmax,jmax))
	ALLOCATE(Tii(jmax*px,kmax/px))
	nm1=imax*jmax*px*5 !allocate for full matrix periodicx and periodicy   !imax*jmax*px*5-2*imax-2*jmax*px 
	ALLOCATE(jco(nm1)) !col nr CSR 
	ALLOCATE(iro(nm1)) !row nr CSR 
	ALLOCATE(di(imax*jmax*px)) !diag nr CSR
	ALLOCATE(di2(imax*jmax*px+1)) !diag nr CSR
	ALLOCATE(beg(imax*jmax*px+1)) !begin new line CSR
	ALLOCATE(LUB(0:nm1)) !LU of LHS
	ALLOCATE(LUBs(0:nm1,kmax/px)) !save LU of LHS
	ALLOCATE(LHS2(nm1)) !LHS
	ALLOCATE(lhs(nm1,kmax/px)) !LHS
	ALLOCATE(wf(0:i1,0:j1,0:k1)) !Wiggle factor blend --> 1=apply AV6 dissipation  0=no AV6 dissipation
	wf = 1.

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
	Ub3old=0.
	Vb3old=0.
	Wb3old=0.
	Ub3new=0.
	Vb3new=0.
	Wb3new=0.	
	
	tmax_inUpuntTSHD=0
	tmax_inVpuntTSHD=0
	tmax_inPpuntTSHD=0

	nu_mol=ekm_mol/rho_b
  
  
	IF (tstart_rms<t_end) then
		ALLOCATE (Uavg(1:imax,1:jmax,1:kmax),Vavg(1:imax,1:jmax,1:kmax),
     1                Wavg(1:imax,1:jmax,1:kmax),Cavg(nfrac,1:imax,1:jmax,1:kmax),
     1                Ravg(1:imax,1:jmax,1:kmax),
     1		      Pavg(1:imax,1:jmax,1:kmax),muavg(1:imax,1:jmax,1:kmax),						
     1		      sigU2(1:imax,1:jmax,1:kmax),sigV2(1:imax,1:jmax,1:kmax),
     1                sigW2(1:imax,1:jmax,1:kmax),sigC2(nfrac,1:imax,1:jmax,1:kmax),
     1                sigR2(1:imax,1:jmax,1:kmax),
     1		      sigUV(1:imax,1:jmax,1:kmax),sigVW(1:imax,1:jmax,1:kmax),
     1                sigUW(1:imax,1:jmax,1:kmax) ,
     1		      sigUC(nfrac,1:imax,1:jmax,1:kmax),sigVC(nfrac,1:imax,1:jmax,1:kmax),
     1                sigWC(nfrac,1:imax,1:jmax,1:kmax),Cmax(nfrac,1:imax,1:jmax,1:kmax),
     1                Cmin(nfrac,1:imax,1:jmax,1:kmax),
     1            Umax(1:imax,1:jmax,1:kmax),Vmax(1:imax,1:jmax,1:kmax),Wmax(1:imax,1:jmax,1:kmax),
     1            Uhormax(1:imax,1:jmax,1:kmax),U3dmax(1:imax,1:jmax,1:kmax),
     1            sig_tau_flow2(1:imax,1:jmax),tau_flow_avg(1:imax,1:jmax),
     1            tau_sl_avg(1:imax,1:jmax),tau_bl_avg(1:imax,1:jmax),tau_mud_avg(1:imax,1:jmax),
     1            sig_tau_sl2(1:imax,1:jmax),sig_tau_bl2(1:imax,1:jmax),sig_tau_mud2(1:imax,1:jmax),	 
     1	          tau_frac_avg(1:nfrac,1:imax,1:jmax),sig_tau_frac2(1:nfrac,1:imax,1:jmax))
!     1            Umin(1:imax,1:jmax,1:kmax),Vmin(1:imax,1:jmax,1:kmax),Wmin(1:imax,1:jmax,1:kmax))
 !    1 			,fUavg(1:imax,1:jmax,1:kmax),sigfU2(1:imax,1:jmax,1:kmax))

		stat_count=0
		stat_time_count=0.
		sigU2=0.
		sigV2=0.
		sigW2=0.
		sigC2=0.
		sigR2=0.
		sigUV=0.
		sigVW=0.
		sigUW=0.
		sigUC=0.	
		sigVC=0.
		sigWC=0.
		sig_tau_flow2=0.

		Uavg=0.
		Vavg=0.
		Wavg=0.
		Cavg=0.
		Cmax=0.
		Cmin=1.e18
		Ravg=0.
		Pavg=0.
		muavg=0.
		Umax=0.
		Vmax=0.
		Wmax=0.
		Uhormax=0.
		U3dmax=0.
		tau_flow_avg=0.
		
		tau_sl_avg=0.
		tau_bl_avg=0.
		tau_mud_avg=0.
		sig_tau_sl2=0.
		sig_tau_bl2=0.
		sig_tau_mud2=0.
		tau_frac_avg=0.
		sig_tau_frac2=0.		
		
		
!		sigfU2=0.
!		fUavg=0.
	endif  

	END SUBROUTINE allocate_global_vars

	SUBROUTINE readtseries(seriesfile,tseries,series)
!	SUBROUTINE readtseries(seriesfile,tseries,series)
!	reads a time series from seriesfile
!	Lynyrd de Wit, October 2010

	CHARACTER*256 seriesfile
	REAL tseries(1:10000),series(1:10000)
	REAL begintime,timestep,endtime
	INTEGER i,k,p,ios

	NAMELIST /timeseries/begintime,timestep,endtime,series
	series=-9999.
	OPEN(10,FILE=seriesfile,IOSTAT=ios,ACTION='read')

		write(*,*) 'File :', seriesfile,'reading'

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

	
		SUBROUTINE readtseries2(seriesfile,tseries,series)
!	SUBROUTINE readtseries(seriesfile,tseries,series)
!	reads a time series from seriesfile
!	Lynyrd de Wit, October 2010

	CHARACTER*256 seriesfile
	REAL tseries(1:10000),series(1:10000)
	REAL begintime,timestep,endtime
	INTEGER i,k,p,ios

	NAMELIST /timeseries/begintime,timestep,endtime,series
	series=-9999.
	OPEN(10,FILE=seriesfile,IOSTAT=ios,ACTION='read')

		write(*,*) 'File :', seriesfile,'reading'

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

	END SUBROUTINE readtseries2

	
		SUBROUTINE readtseries3(seriesfile,tseries,series,nfrac)
!	SUBROUTINE readtseries(seriesfile,tseries,series)
!	reads a time series from seriesfile
!	Lynyrd de Wit, October 2010

	CHARACTER*256 seriesfile
	REAL tseries(1:10000),series(1:nfrac,1:10000)
	REAL begintime,timestep,endtime
	INTEGER i,k,p,ios,nfrac

	NAMELIST /timeseries/begintime,timestep,endtime,series
	series=-9999.
	OPEN(10,FILE=seriesfile,IOSTAT=ios,ACTION='read')

		write(*,*) 'File :', seriesfile,'reading'

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
	DO WHILE (series(1,k).NE.-9999)
		tseries(k)=begintime+REAL(k)*timestep-timestep
		k=k+1
	END DO


	END SUBROUTINE readtseries3
	
	END MODULE nlist
