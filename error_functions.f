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


	MODULE error_functions
	! MODULE error_functions
	! This module contains the suboutine "writeerror" and "finelhelp" 
	! Gerard Dam, 4 March 2007
	!----------------------------------------------------------
	! This module contains the following subroutines/functions:
	!----------------------------------------------------------
	! SUBROUTINE writeerror(error): writes an error on screen
	! SUBROUTINE finelhelp: writes some help on screen 
	!----------------------------------------------------------

	CONTAINS


	subroutine writeerror(error)
	implicit none

	integer error

		write(*,*) ' '
		write(*,*) '##################################################'
		write(*,*) '#                                                #'
		write(*,*) '# ERROR DETECTED:                                #'
		write(*,*) '#                                                #'

	SELECT CASE (error)
	CASE(001)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# px is not defined                              #'			
	CASE(002)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# imax is not defined                            #'			
		write(*,*) '#                                                #'
	CASE(003)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# jmax is not defined                            #'			
		write(*,*) '#                                                #'
	CASE(004)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# kmax is not defined                            #'			
		write(*,*) '#                                                #'
	CASE(005)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# imax_grid is not defined                       #'			
		write(*,*) '#                                                #'
	CASE(006)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# dr_grid is not defined                         #'			
		write(*,*) '#                                                #'
	CASE(007)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# Rmin is not defined                            #'			
		write(*,*) '#                                                #'
	CASE(008)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# schuif_x is not defined                        #'			
		write(*,*) '#                                                #'
	CASE(009)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# dy or dy_grid is not defined                   #'			
		write(*,*) '#                                                #'
	CASE(010)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# depth is not defined                           #'			
		write(*,*) '#                                                #'
	CASE(011)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# jmax is not a multiple of px                   #'			
		write(*,*) '# adjust jmax and/or px                          #'
	CASE(012)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# kmax is not a multiple of px                   #'			
		write(*,*) '# adjust kmax and/or px                          #'
	CASE(013)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# sym_grid_y is not equal to 0,-1,1              #'			
	CASE(030)
		write(*,*) '# namelist &times                                #'			
		write(*,*) '# t_end is not defined                           #'			
	CASE(031)
		write(*,*) '# namelist &times                                #'			
		write(*,*) '# t0_output,dt_output,te_output must be >0       #'			
	CASE(032)
		write(*,*) '# namelist &times                                #'			
		write(*,*) '# tstart_rms is not defined                      #'			
	CASE(033)
		write(*,*) '# namelist &times                                #'			
		write(*,*) '# dt_max is not defined                          #'			
	CASE(034)
		write(*,*) '# namelist &times                                #'			
		write(*,*) '# time_int is not defined                        #'	
		write(*,*) '# choose ''EE1'',''AB2'',''AB3'',''ABv'',or ''RK3''            #'	
		write(*,*) '# for Euler_expl1, AdamsBashforth2,              #'	
		write(*,*) '# for AdamsBashforth3,AdamsBashforth variable dt #'	
		write(*,*) '# (based on AB3), RungeKutta3 time integration   #'	

	CASE(035)
		write(*,*) '# namelist &times                                #'			
		write(*,*) '# CFL is not defined                             #'	
	CASE(036)
		write(*,*) '# namelist &times                                #'			
		write(*,*) '# te_rms must be >tstart_rms                     #'	
	CASE(039)
		write(*,*) '# namelist &ambient                                        #'			
		write(*,*) '# timeseries file mentioned above should be long enough    #'			
	CASE(040)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# U_b is not defined                             #'	
	CASE(041)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# V_b is not defined                             #'	
	CASE(042)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# W_b is not defined                             #'	
	CASE(043)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# rho_b is not defined                           #'	
	CASE(044)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# SEM is not defined                             #'	
	CASE(045)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# nmax2 is not defined                           #'	
	CASE(046)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# nmax1 is not defined                           #'	
	CASE(047)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# lm_min is not defined                          #'	
	CASE(048)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# slip_bot is not defined                        #'	
	CASE(049)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# kn is not defined                              #'	
	CASE(050)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# interaction_bed is not defined                 #'	
	CASE(051)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input bc file does not exist                   #'	
	CASE(052)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input periodicx not defined (0 or 1)           #'	
	CASE(053)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input periodicy not defined (0,1 or 2)         #'	
	CASE(054)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input dpdx not defined but periodicx=1         #'	
	CASE(055)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input Hs should be 0<Hs<depth                  #'	
	CASE(056)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input Tp should be Tp>0                        #'	
	CASE(057)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input wavedirection nx_w,ny_w are incorrect    #'	
		write(*,*) '# sqrt(nx_w^2+ny_w^2) should be 1                #'	
	CASE(058)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# obst%x(1:4) and obst%y(1:4) should be defined  #'	
	CASE(059)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# obst%height should be positive                 #'
	CASE(060)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# W_j is not defined                             #'	
	CASE(061)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# Awjet is not defined                           #'	
	CASE(062)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# Aujet is not defined                           #'	
	CASE(063)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# Avjet is not defined                           #'	
	CASE(064)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# kjet is not defined                            #'	
	CASE(066)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# radius_j is not defined                        #'	
	CASE(067)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# Sc is not defined                              #'
	CASE(068)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# slipvel is not defined, should be 0,1,2        #'
	CASE(069)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# Strouhal is not defined                        #'
	CASE(070)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# azi_n is not defined                           #'
	CASE(071)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# ws is not correct for each fraction            #'		
	CASE(072)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# rho is not correct for each fraction           #'		
	CASE(073)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# c is not correct for each fraction             #'		
	CASE(074)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# c (volume fraction) is larger than 1, physically impossible      #'		
	CASE(075)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# dpart is not correct for each fraction         #'	
	CASE(076)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# tau_e [N/m2] is not correct for each fraction  #'	
	CASE(077)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# tau_d [N/m2] is not correct for each fraction  #'	
	CASE(078)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# M [kg/sm2] is not correct for each fraction    #'	
	CASE(079)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# dfloc is negative or less than dpart           #'	
	CASE(080)
		write(*,*) '# namelist &LESmodel                             #'			
		write(*,*) '# Cs is not defined                              #'	
	CASE(081)
		write(*,*) '# namelist &LESmodel                             #'			
		write(*,*) '# Lmix_type is not defined                       #'	
	CASE(082)
		write(*,*) '# namelist &LESmodel                             #'			
		write(*,*) '# sgs_model is not defined                       #'	
		write(*,*) '# choose ''SSmag'',or ''DSmag'',or ''FSmag'',or''SWALE'',or''Sigma'',or''MixLe'',or''ReaKE''   #'	
		write(*,*) '# for Standard, Filtered Smagorinsky, WALE,      #'
		write(*,*) '# Sigma, Mixing Length RANS model, Realizible K-Eps model  #'
	CASE(083)
		write(*,*) '# namelist &LESmodel                             #'			
		write(*,*) '# nr_HPfilter is not defined                     #'	
	CASE(084)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# kn_sed [m] is not correct for each fraction    #'	
	CASE(085)
		write(*,*) '# namelist &LESmodel                             #'			
		write(*,*) '# damping_drho_dr is not defined                 #'	
		write(*,*) '# choose ''none'',or ''MuAn''                        #'	
		write(*,*) '# for Munk-Anderson damping of eddy viscosity at #'
		write(*,*) '# stably stratification                          #'
	CASE(086)
		write(*,*) '# namelist &LESmodel                             #'			
		write(*,*) '# damping_a1 or a2 [-] is not defined            #'	
	CASE(087)
		write(*,*) '# namelist &LESmodel                             #'			
		write(*,*) '# damping_b1 or b2 [-] is not defined            #'	
	CASE(088)
		write(*,*) '# namelist &LESmodel                             #'			
		write(*,*) '# extra_mix_visc is not defined                  #'	
		write(*,*) '# choose ''none'',or ''Krie''                        #'	
	CASE(089)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# Cs_relax must be between 0 and 1               #'			
	CASE(090)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# kappa is not defined                           #'	
	CASE(091)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# gx is not defined                              #'	
	CASE(092)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# gy is not defined                              #'	
	CASE(093)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# gz is not defined                              #'	
	CASE(094)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# ekm_mol is not defined                         #'
	CASE(095)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# pickup_formula is not defined                  #'
		write(*,*) '# or pickup_formula_swe is not defined           #'
		write(*,*) '# choose ''vanrijn1984'',or ''nielsen1992'',or ''okayasu2010'' #'	
		write(*,*) '# or ''vanrijn2019'', or ''VR2019_Cbed'' #'	
		write(*,*) '# or ''VR1984_Cbed'' #'
	CASE(096)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# kn_d50_multiplier and kn_d50_multiplier_bl must be >0.                  #'		
	CASE(097)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# avalanche_slope must be >=0.                   #'			
	CASE(098)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# calibfac_sand_pickup and calibfac_sand_bedload must be >=0.              #'	
	CASE(099)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# calibfac_Shields_cr and calibfac_Shields_cr_bl must be >=0.               #'			
	CASE(100)
		write(*,*) '# Namelist error occured                         #'			
		write(*,*) '# Namelist involved: HISTORIES                   #'
	CASE(101)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# morfac must be >=0.                            #'	
	CASE(102)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# morfac must be >=1.                            #'			
	CASE(103)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# vwal must be >=0.                              #'			
	CASE(104)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# nl = porosity loose sand after dilantancy      #'			
		write(*,*) '# nl must be >0 and must be >n0                  #'			
		write(*,*) '# n0 = 1-cfixedbed = porosity undisturbed bed    #'		
	CASE(105)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# permeability_kl must be >=0.                   #'		
		write(*,*) '# relation of Kozeny-Carmen could be used for kl #'		
		write(*,*) '# kl = gD15^2/160nu*nl^3/(1-nl)^2                #'		
	CASE(106)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# k_layer_pickup must be >=1 (default 1)         #'		
	CASE(107)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# bl_relax & sl_relax must be between 0 and 1    #'		
	CASE(108)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# fcor must be >=0.                              #'	
	CASE(109)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# bedload_formula is defined incorrectly         #'
		write(*,*) '# choose ''vanrijn2007'',or ''vanrijn2003'',or ''MeyPeMu1947'' #'	
		write(*,*) '# or comment the line and do not define a bedload_formula #'	
	CASE(110)
		write(*,*) '# i history > imax                               #'			
		write(*,*) '#                                                #'
	CASE(111)
		write(*,*) '# i history < 1                                  #'			
		write(*,*) '#                                                #'
	CASE(120)
		write(*,*) '# j history > jmax*px                            #'			
		write(*,*) '#                                                #'
	CASE(121)
		write(*,*) '# j history < 1                                  #'			
		write(*,*) '#                                                #'
	CASE(130)
		write(*,*) '# k history > kmax                               #'			
		write(*,*) '#                                                #'
	CASE(131)
		write(*,*) '# k history < 1                                  #'			
		write(*,*) '#                                                #'
	CASE(141)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# av_slope_z must be increasing and start with 0 #'		
		write(*,*) '# and end with depth                             #'		
	CASE(142)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# avalanche_slope series must be >0 and as long  #'		
		write(*,*) '# av_slope_z                                     #'		
	CASE(143)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# bedslope_effect must be 0,1,2,3,4,5,6,7        #'		
		write(*,*) '# 0 (default) no bedslope effect                 #'		
		write(*,*) '# 1 = adjust Shields_cr following Roulund,       #'		
		write(*,*) '# Fredsoe etal. 2004 for suspension and bedload  #'			
		write(*,*) '# 2 = like 1 but only for bedload                #'		
		write(*,*) '# 3 = bed-slope influence on bedload D3D style   #'		
		write(*,*) '# Bagnold (1966) for longitudinal slope and Ikeda (1982, 1988) for transverse slope   #'	
		write(*,*) '# 4 = like 1 but only for Shields in numerator   #'
		write(*,*) '# 5 = like 2 but only for Shields in numerator   #'
		write(*,*) '# 6 = combi 3 for bedload and 1 for susload but only for Shields in numerator   #'
	CASE(144)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# phi_sediment must be >0. and in degrees        #'		
		write(*,*) '# default phi_sediment=30.                       #'		
	CASE(145)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# wallmodel_tau_sed must be 1,3,4,5,8,9 or 11    #'		
		write(*,*) '# default wallmodel_tau_sed = 1                  #'		
	CASE(146)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# power_VR2019 must be >0                        #'		
	CASE(171)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# type must be 1 (default),2,3 or -1             #'	
	CASE(172)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# fractions must be ordered with increasing dpart#'	
	CASE(173)
		write(*,*) '# namelist &plume -->bedplume                    #'			
		write(*,*) '# move_dx_series,move_dy_series,move_dz_series   #'			
		write(*,*) '# move_nx_series,move_ny_series,                 #'
		write(*,*) '# move_zbed_criterium not same length            #'
		write(*,*) '# make all series same length                    #'			
		write(*,*) '# EXCEPTION: move_zbed_criterium,move_nx_series, #'	
		write(*,*) '# move_ny_series may be length 1 or length series #'	
	CASE(174)
		write(*,*) '# namelist &plume -->bedplume                    #'			
		write(*,*) '# move_filename_series  must have length         #'			
		write(*,*) '# of move_dx_series, or length zero              #'
		write(*,*) '# when length zero then always a outputfile is   #'
		write(*,*) '# generated when bedplume is moved               #'
	CASE(175)
		write(*,*) '# namelist &fractions_in_plume                   #'			
		write(*,*) '# CD must be >0 or 0. for not using it           #'	
		
	CASE(200)
		write(*,*) '#                                                #'			
		write(*,*) '# input file not found                           #'
	CASE(260)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# U_j2 is not defined                            #'	
	CASE(261)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# Awjet2 is not defined                          #'	
	CASE(262)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# Aujet2 is not defined                          #'	
	CASE(263)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# Avjet2 is not defined                          #'	
	CASE(264)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# Strouhal2 is not defined                       #'
	CASE(265)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# azi_n2 is not defined                          #'
	CASE(266)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# radius_j2 is not defined                       #'	
	CASE(267)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# zjet2 is not defined                           #'
	CASE(268)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# zjet2 should be <depth                         #'
	CASE(269)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%x(1:4) and bedplume%y(1:4) should be defined  #'	
	CASE(270)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%height should be positive             #'
	CASE(271)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%c(1:nfrac), bedplume%sedflux(1:nfrac) should be positive         #'
	CASE(272)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%forever should be 0 or 1              #'
	CASE(273)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%t0 and t_end should be >0.            #'
	CASE(274)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%zbottom should be >0. and < height    #'
	CASE(275)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# radius_inner_j must be >0 and < radius_j       #'	
	CASE(276)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# W_j_powerlaw must be >0 (default 7 for 1/7 powerlaw)  #'	
	CASE(277)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%c(1:nfrac), bedplume%sedflux(1:nfrac) cannot both be defined   #'	
	CASE(278)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%changesedsuction should be >0. and < 1. #'		
	CASE(279)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# timeseries file should be long enough          #'		
	CASE(280)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# hindered_settling should be 1,2,3              #'		
	CASE(281)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume covers >100000 cells in x,y plane     #'		
	CASE(282)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# hindered_settling_c should be 0,1              #'		
	CASE(283)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%kn_flow should be >0.                 #'				
	CASE(301)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# please prescribe U_TSHD                        #'	
	CASE(302)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# Lfront<0 or Lfront>LOA                         #'	
	CASE(303)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# Breadth<0                                      #'	
	CASE(304)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# Draught<0                                      #'	
	CASE(305)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# Lback<0 or Lback>LOA                           #'	
	CASE(306)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# Hback<0                                        #'	
	CASE(307)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# xfront>0 (overflow not in TSHD)                #'	
	CASE(308)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# yfront>0.5*Breadth (overflow not in TSHD)      #'	
	CASE(310)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# nprop is not 1 or 2                            #'	
	CASE(311)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# Dprop is not defined                           #'	
	CASE(312)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# xprop is not defined                           #'	
	CASE(313)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# yprop is not defined                           #'
	CASE(314)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# zprop is not between 0 and depth               #'		
	CASE(315)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# Pprop is not defined                           #'	
	CASE(316)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# kn_TSHD is not defined                         #'	
	CASE(317)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# rot_prop is not defined                        #'	
	CASE(318)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# rudder is not defined                          #'	
	CASE(319)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# draghead is not defined                        #'	
		write(*,*) '# choose ''star'',''port'',''both'', or ''none''                #'	
	CASE(320)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# Dsp should be positive                         #'	
	CASE(321)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# softnose should be 0 or 1                      #'	
	CASE(322)
		write(*,*) '# namelist &ship                                 #'			
		write(*,*) '# Hfront should be >0.                           #'	
	CASE(401)
		write(*,*) '# namelist &num_scheme                           #'			
		write(*,*) '# convection is not defined                      #'	
		write(*,*) '# choose ''CDS2'',''CDS6'' for 6th order central #'
		write(*,*) '# combined with 6th diss (numdiff=1/60 gives UPW5) #'	
		write(*,*) '# , or ''COM4'' for compact cds4, or ''CDS4'' for #'
		write(*,*) '# cds4 with 6th diss in which only advected      #'
		write(*,*) '# velocity is discretised with a variant of cds4 #'
		write(*,*) '# , or ''HYB6'' for cds2 with 6th diss.          #'
		write(*,*) '# , or ''HYB4'' for cds2 with 4th diss.          #'
		write(*,*) '# , or ''C4A6'' for cds4 with 6th diss. in which #'
		write(*,*) '# complete stencil (A-Domis) is cds4 for advected and advective vel. #'
		write(*,*) '# ,or ''C2Bl'' for cds2 with blend to upw1 for TVD-r<0; slope blend = numdiff #'		
		write(*,*) '# ,or ''C4Bl'' for cds4 with blend to upw1 for TVD-r<0; slope blend = numdiff #'
		write(*,*) '# (cds4 variant in which only advected velocity is a variant of cds4 #'

	CASE(402)
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# numdiff is not between 0 and 1                 #'	
	CASE(403)
		write(*,*) '# namelist &num_scheme                           #'			
		write(*,*) '# diffusion is not defined                       #'	
		write(*,*) '# choose ''CDS2'',or ''COM4''                        #'	
		write(*,*) '# for 2nd order central scheme (non compact), or 4th order compact   #'
	CASE(404)
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# comp_filter_a is not between 0 and 0.5         #'	
	CASE(405)
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# comp_filter_n is not positive                  #'	
	CASE(406)
		write(*,*) '# namelist &num_scheme                                        #'		
		write(*,*) '# CNdiffz must be 0 (explicit diffusion in z dir)             #'
		write(*,*) '# or 1 (Crank Nicolson semi-implicit diffusion in z dir)      #'
		write(*,*) '# or 2 (Euler backward implicit diffusion in z dir)           #'
		write(*,*) '# or 11 (Crank Nicolson semi-implicit diffusion in all dirs)  #'
		write(*,*) '# or 12 (Euler backward implicit diffusion in all dirs)       #'		
		write(*,*) '# or 31 fully 3D Crank Nicolson diffusion in all dirs         #'
	CASE(407)
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# npresIBM must be >0                            #'
	CASE(408)
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# wiggle_detector should be 0 or 1 or 2 or >20   #'			
	CASE(409)
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# k_ust_tau should be between 1 and kmax         #'	
		write(*,*) '# and k_ust_tau_flow should be between 1 and kmax#'	
	CASE(410)
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# k_pzero should be between 1 and kmax           #'	
	CASE(411)
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# CNdiff_dtfactor should be >=1.                 #'	
	CASE(412)
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# k_ust_tau_sed_range(1) and k_ust_tau_sed_range(2) should be between 1 and kmax #'	
		write(*,*) '# and k_ust_tau_sed_range(1)<k_ust_tau_sed_range(2)                              #'			
	CASE(501)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# imax_grid or jmax_grid is not increasing monotonously       #'	
	CASE(502)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# not enough or negative dr_grid or dy_grid are defined     #'	
	CASE(503)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# last imax_grid is not equal to imax, or        #'	
		write(*,*) '# last jmax_grid is not equal to jmax or half jmax with sym_grid_y=1           #'	
	CASE(504)
		write(*,*) '# namelist &LESmodel                             #'			
		write(*,*) '# Lmix_type is not 1 or 2                        #'	
	CASE(601)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# bc_obst_h should be >0 and <depth              #'
	CASE(602)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# U_b3 should be defined                         #'
	CASE(603)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# V_b3 should be defined                         #'
	CASE(604)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# surf_layer should be >0 and < depth            #'
	CASE(605)
		write(*,*) '# namelist &ambient                              #'
		write(*,*) '# wallup should be 0,1 or 2 (default 0)          #'
	CASE(606)
		write(*,*) '# file does not exist             #'	
	CASE(607)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# cfixedbed must be >=0 and <=1 (default 0.6)    #'			
	CASE(608)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# c_bed must be >0 and <1                        #'			
	CASE(609)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# sum of c_bed must be >=cfixedbed and <=1       #'	
	CASE(610)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# kn_mp monopile roughness must be defined       #'			
	CASE(611)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input U_w must be defined if Hs>0              #'			
		write(*,*) '# please note that unlike U_b, U_w is defined in earth fixed x,y coordinate system !    #'		
	CASE(612)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input V_w must be defined if Hs>0              #'					
		write(*,*) '# please note that unlike V_b, V_w is defined in earth fixed x,y coordinate system !    #'	
	CASE(613)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# istart_morf1, istart_morf2, and/or i_periodicx must be >0 and <imax  #'		
		write(*,*) '# and it should be: istart_morf2(1) > istart_morf1(1)   #'		
		write(*,*) '# and if also no-bed-update zone at right end of domain #'
		write(*,*) '# then it should be: istart_morf2(2) < istart_morf1(2)  #'		
	CASE(614)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# kn_flow_file contains NaN or negative value    #'			
	CASE(620)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# nmax3 is not defined                           #'			
	CASE(621)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# lm_min3 is not defined                         #'	
	CASE(622)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# cbc_perx_j(1) should be <= cbc_perx_j(2)       #'			
		write(*,*) '# cbc_perx_j(1) should be >0                     #'			
		write(*,*) '# cbc_perx_j(2) should be <=jmax/px              #'			
	CASE(623)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# cbc_relax and TBLE_grad_relax and TBLEsl_grad_relax and TBLEbl_grad_relax should be between 0-1     #'	
	CASE(624)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# taulayerTBLE should be 1 (default) or 50       #'			
		write(*,*) '# 1 means using tau_bottom from TBLE as bc CFD   #'			
		write(*,*) '# 50 means using tau_top from TBLE as bc CFD     #'		
	CASE(625)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# obstfile_erodepo should be 1 or 2              #'	
	CASE(626)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# TBLEsed_grad_relax is no longer used           #'			
		write(*,*) '# use TBLEsl_grad_relax and TBLEbl_grad_relax instead  #'	
	CASE(627)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# dpdx_ref_j(1) should be <= dpdx_ref_j(2)       #'			
		write(*,*) '# dpdx_ref_j(1) should be >0                     #'			
		write(*,*) '# dpdx_ref_j(2) should be <=jmax/px              #'			
   	CASE(701)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# restart_dir does not give correct files        #'	
	CASE(702)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# restart_dir does not contain correct size data #'	
	CASE(703)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# avfile does not exist                          #'	
	CASE(704)
		write(*,*) '# namelist &constants                            #'			
		write(*,*) '# avfile does not contain correct data           #'	
	CASE(705)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# obstfile does not contain correct size data    #'	
	CASE(706)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# obstfile does not contain correct length       #'			
		write(*,*) '# ostacle_starttimes vector                      #'			
		write(*,*) '# (should be equal to number of obstacle files)  #'
	CASE(801)															!new error functions
		write(*,*) '# namelist &rheology                             #'
		write(*,*) '# Non_Newtonian should be 0,1 or 2 (default 0)   #'
	CASE(802)
		write(*,*) '# namelist &rheology	                         #'			
		write(*,*) '# Rheological_model is not defined               #'
		write(*,*) '# choose ''SIMPLE'' for simple Bingham#          #'
		write(*,*) '# choose ''JACOBS'' for Jacobs and van Kesteren  #'
		write(*,*) '# choose ''WINTER'' for Winterwerp and Kranenburg#'
		write(*,*) '# choose ''THOMAS'' for Thomas                   #'
	CASE(803)
		write(*,*) '# namelist &rheology                             #'
		write(*,*) '# Apvisc_shear_relax should be between 0 or 1    #'
		write(*,*) '# (default 1 no relaxation)                      #'		
!	CASE(804)
!		write(*,*) '# namelist &rheology                             #'
!		write(*,*) '# thixotropy should be 0 or 1 (default 0)        #'		
	CASE(1001)
		write(*,*) '# namelist &timeseries                           #'			
		write(*,*) '# file cannot be found                           #'	


	CASE(1002)
		write(*,*) '# namelist &timeseries                           #'			
		write(*,*) '# error in reading timeseries                    #'	
		write(*,*) '# should be a NAMELIST /timeseries/ begintime,timestep,endtime,series #'	
		write(*,*) '# ended with a "/"                               #'	
	CASE(1051)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# obst%ero or obst%depo should be between 0-1    #'
	CASE(10000)
		write(*,*) '# timestep becomes too small                     #'			
		write(*,*) '# restart with smaller dt                        #'	

	CASE DEFAULT
		write(*,*) '# error number: ', error
		write(*,*) '# Unexpected error call                          #'	  
		write(*,*) '# check error_functions.for                      #'	



	END SELECT
		write(*,*) '#                                                #'
		write(*,*) '##################################################'
		write(*,*)
		write(*,*) '      ERROR(S) FOUND, DFLOW3D TERMINATED          '

		STOP


	END SUBROUTINE writeerror	

		

	SUBROUTINE finelhelp
	! SUBROUTINE finelhelp
	! This subroutine shows some help information on screen of how
	! to use FINEL2d. 
	! This can be seen by using the command: <FINEL2D_PAR.exe help>
	! Afterwards the program is stopped.
	! Gerard Dam, 4 March 2007
	
	! UPDATE NECCESARRY!!!!
		write(*,*) '##################################################'	
		write(*,*) '#                                                #'	
		write(*,*) '#          FINEL 2D, HELP FUNCTION               #'	
		write(*,*) '#                                                #'	
		write(*,*) '##################################################'	
		write(*,*) '    '	
		write(*,*) ' ================================================'
		write(*,*) ' The following abbreviations are used:'
		write(*,*) ' [M]=mandatory, [O]=optional, DEF=default value'
		write(*,*) ' val=real value, ival=integer value'
		write(*,*) ' str=string variable, iarray=array of integers'
		write(*,*) ' ================================================'
		write(*,*) ' Inputfile should look like this:'	
		write(*,*)
	    write(*,*) '&modtype'
		write(*,*) '[M] modeltype=[1,2,3]; 1=flow,2=flow+silt,3=flow+sand, DEF = 1'
		write(*,*) '/'
		write(*,*) ' '
		write(*,*) '&general'
		write(*,*) '[O] restartfname=[str]; name of the restartfile'
		write(*,*) '[M] meshfname=[str]; name of sepran mesh file'
		write(*,*) '[O] depthfname=[str]; name of depth file [m](pos. downwards)'
		write(*,*) '[O] globaldepth=[val]; global depth [m] (pos.downwards)'
		write(*,*) '[M] roughnesstype=[nikkuradse,chezy,mannning]; type of roughness' 
          write(*,*) '         DEF=nikkuradse'
		write(*,*) '[O] roughnessfname=[str]; name of roughness file [unit of roughnesstype]'
		write(*,*) '[O] globalroughness=[val]; global roughness [unit of roughnesstype]'
		write(*,*) '/'
		write(*,*) ' '

     		write(*,*) '&times'
		write(*,*) '[M] tbeg=[ival]; starttime of calculation[s], DEF=0'
		write(*,*) '[M] tend=[ival]; endtime of calculation [s]'
		write(*,*) '[O] tmap=[ival1 ival2 ival3]; begin,timestep,endtime of writing maps [s]'
		write(*,*) '[O] this=[ival1 ival2 ival3]; begin,timestep,endtime of writing histories [s]'
		write(*,*) '/'
		write(*,*) ' '

		write(*,*) '&flowconstants'
		write(*,*) '[O] gravity=[val]; gravity [m/s2],DEF=9.81'
		write(*,*) '[O] windsp_x=[val]; windspeed in x-dir [m/s],DEF=0'
		write(*,*) '[O] windsp_y=[val]; windspeed in y-dir [m/s],DEF=0'
		write(*,*) '[O] windsp_drag=[val]; winddrag coef., DEF=0'
		write(*,*) '[O] latitude=[val]; Latitude of modelarea [degree],DEF=52 (The Netherlands)'
	    write(*,*) '[O] drycrit=[val]; Critical depth for drying [m],DEF=0.03'
		write(*,*) '/'
		write(*,*) ' '


		write(*,*) '&siltconstants'
		write(*,*) '[O] ws=[val]; fall velocity of silt [m/s],DEF=0.003'
		write(*,*) '[O] M=[val]; ,DEF=0.0003'
		write(*,*) '[O] tau_s=[val]; critical sedimentation shear stress [N/m2], DEF=2.0'
		write(*,*) '[O] tau_e=[val]; critival erosion shear stress [N/m2], DEF=0.8'
		write(*,*) '[O] rho=[val]; dry bulk density [kg/m3], DEF=250'
		write(*,*) '[O] pointsourcefname=[str]; name of pointsource file'
		write(*,*) '/'
		write(*,*) ' '

		write(*,*) '&initcond'
		write(*,*) '[O] hinit=[val]; initial waterdepth [m], DEF=0'
		write(*,*) '[O] uinit=[val]; initial u-velocity [m/s], DEF=0'
		write(*,*) '[O] vinit=[val]; initial v-velocity [m/s], DEF=0'
		write(*,*) '[O] cinit=[val]; initial silt-concentration [mg/l], DEF=0'
		write(*,*) '/'
		write(*,*) ' '
		
		write(*,*) '&histories'
		write(*,*) '/'
		write(*,*) ' '
		write(*,*) '&boundcond'
		write(*,*) '/'
		write(*,*) ' '

		STOP
	END SUBROUTINE finelhelp

	END MODULE error_functions
