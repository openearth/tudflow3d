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
		write(*,*) '# dy is not defined                              #'			
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
		write(*,*) '# slipvel is not defined                         #'
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
		write(*,*) '# choose ''SSmag'',or ''FSmag'',or''SWALE'',or''Sigma'',or''MixLe''  #'	
		write(*,*) '# for Standard, Filtered Smagorinsky, WALE,      #'
		write(*,*) '# Sigma, Mixing Length model                     #'
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
	CASE(100)
		write(*,*) '# Namelist error occured                         #'			
		write(*,*) '# Namelist involved: HISTORIES                   #'
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
		write(*,*) '# bedplume%c(1:nfrac) should be positive         #'
	CASE(272)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%forever should be 0 or 1              #'
	CASE(273)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%t0 should be >0. (default 0.s)        #'
	CASE(274)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# bedplume%zbottom should be >0. and < height    #'
	CASE(275)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# radius_inner_j must be >0 and < radius_j       #'	
	CASE(276)
		write(*,*) '# namelist &plume                                #'			
		write(*,*) '# W_j_powerlaw must be >0 (default 7 for 1/7 powerlaw)  #'	
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
		write(*,*) '# choose ''CDS2'',''CDS6'', or ''COM4'',or              #'	
		write(*,*) '# or ''HYB4'', or ''HYB6'' , or ''C4A6''                 #'
		write(*,*) '# for 2nd order or 6th order central scheme      #'	
		write(*,*) '# or 4th order compact or cds2 with 4th diss.    #'
		write(*,*) '# or cds2 with 6th diss.    #'
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
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# CNdiffz must be 0 (explicit diffusion in z dir)#'
		write(*,*) '# or 1 (Crank Nicolson implicit diffusion in z dir)      #'
	CASE(407)
		write(*,*) '# namelist &num_scheme                           #'		
		write(*,*) '# npresIBM must be >0                            #'
	CASE(501)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# imax_grid is not increasing monotonously       #'	
	CASE(502)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# not enough or negative dr_grid are defined     #'	
	CASE(503)
		write(*,*) '# namelist &simulation                           #'			
		write(*,*) '# last imax_grid is not equal to imax            #'	
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
                write(*,*) '# wallup should be 0 or 1 (default 0)            #'
	CASE(606)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input bedlevel file does not exist             #'	
	CASE(611)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input U_w must be defined if Hs>0              #'			
		write(*,*) '# please note that unlike U_b, U_w is defined in earth fixed x,y coordinate system !    #'		
	CASE(612)
		write(*,*) '# namelist &ambient                              #'			
		write(*,*) '# input V_w must be defined if Hs>0              #'					
		write(*,*) '# please note that unlike V_b, V_w is defined in earth fixed x,y coordinate system !    #'		
	CASE(1001)
		write(*,*) '# namelist &timeseries                           #'			
		write(*,*) '# file cannot be found                           #'	


	CASE(1002)
		write(*,*) '# namelist &timeseries                           #'			
		write(*,*) '# error in reading timeseries                    #'	
		write(*,*) '# should be a NAMELIST /timeseries/ begintime,timestep,endtime,series #'	
		write(*,*) '# ended with a "/"                               #'	

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
