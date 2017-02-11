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

      MODULE sediment
      USE nlist	

      IMPLICIT NONE
      

       CONTAINS

       subroutine init_sediment

       implicit none

 	INTEGER n
	REAL Re_p,atm_pres,rhoair_z,zz,ddzz

	DO n=1,nfrac
	  Re_p=ABS(frac(n)%ws)*frac(n)%dfloc/(ekm_mol/rho_b)
!	hindered settling function acc. to Rowe (1987), 
!	is the smoothed representation of the original Richardson and Zaki (1954) relations:
	  frac(n)%n=(4.7+0.41*Re_p**0.75)/(1.+0.175*Re_p**0.75) 
	ENDDO

	rhocorr_air_z=1.
	atm_pres=101325. !standard atmospheric pressure in Pa 
	DO n=1,nair !correction for comppressible air
		DO k=0,k1
			zz=k*dz-0.5*dz
			ddzz=zz-(depth-frac(nfrac_air(n))%zair_ref_belowsurf) !positive upward and defined wrt zair_ref_belowsurf
			rhoair_z = ((frac(nfrac_air(n))%zair_ref_belowsurf-ddzz)*rho_b*ABS(gz)+atm_pres) / 
     &			       ((frac(nfrac_air(n))%zair_ref_belowsurf)*rho_b*ABS(gz)+atm_pres) * frac(nfrac_air(n))%rho
			rhocorr_air_z(nfrac_air(n),k) = frac(nfrac_air(n))%rho / rhoair_z
		ENDDO
	ENDDO
	
	
      END SUBROUTINE init_sediment

      SUBROUTINE slipvelocity(csed,Wcfd,wsed,rr,k1b,k1e,Wfluid,dt,dz)

	implicit none

	INTEGER n,kp ! local variables!
	REAL    wsed(nfrac,0:i1,0:j1,0:k1)
	REAL    csed(nfrac,0:i1,0:j1,0:k1),csed2(nfrac,0:i1,0:j1,0:k1)
	REAL    Wcfd(0:i1,0:j1,0:k1),Wfluid(0:i1,0:j1,0:k1)
	REAL	rr(0:i1,0:j1,0:k1),rr2(0:i1,0:j1,0:k1)
	REAL	ctot,sum_c_ws,ws(nfrac),ccc(nfrac)
	INTEGER k1b,k1e,kpp,km,iter,kplus
	REAL dt_dzi,noemer,rrpos,rrneg,limiter,dt,dz,ws_basis(nfrac)

	dt_dzi=dt/dz 

	IF (k1b.eq.1.and.k1e.eq.(k1-1)) THEN ! ws in flow domain is determined
	  DO n=1,nfrac
	    ws_basis(n)=frac(n)%ws
	  ENDDO
	ELSE ! ws at bed for ero-depo is determined
	  DO n=1,nfrac
	    ws_basis(n)=frac(n)%ws_dep
	  ENDDO
	ENDIF

	csed2=csed
	rr2=rr 
	if (nobst>0.or.bedlevelfile.ne.''.or.interaction_bed.eq.4.or.interaction_bed.eq.5) then
	 DO i=0,i1
	  DO j=0,j1
		kplus = MIN(kbed(i,j)+1,k1)
		rr2(i,j,kbed(i,j))=rr2(i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
		DO n=1,nfrac
			csed2(n,i,j,kbed(i,j))=csed2(n,i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
	    ENDDO
	  ENDDO
	 ENDDO
	endif


	IF (slipvel.eq.2) THEN ! no hindered settling
	 DO i=0,i1
	  DO j=0,j1
	    DO k=k1b,k1e ! at k=0 wsed should be zero at kmax wsed is calculated
	      ctot=0.
	      sum_c_ws=0.
	      kp = MIN(k+1,k1)
	      DO n=1,nfrac
	        ccc(n) = 0.5*(csed2(n,i,j,k) + csed2(n,i,j,kp))	
	        ws(n)=ws_basis(n)
		! ws is positive downward, wsed is positive upward 
  		sum_c_ws=sum_c_ws+ws(n)*ccc(n)*frac(n)%rho/(0.5*(rr2(i,j,k)+rr2(i,j,kp)))
!		sum_c_ws=sum_c_ws+ws(n)*0.5*(csed(n,i,j,k)*frac(n)%rho/rr(i,j,k)+csed(n,i,j,kp)*frac(n)%rho/rr(i,j,kp))
		! According to drift velocity literature the drift flux is calculated using the mass-fraction in stead of volume-fraction,
		! therefore an extra frac(n)%rho/rho_mix is included
	      ENDDO
	      DO n=1,nfrac
		wsed(n,i,j,k)=Wcfd(i,j,k)+sum_c_ws-ws(n) ! wsed is positive upward
	      ENDDO
	      Wfluid(i,j,k)=sum_c_ws !Wcfd(i,j,k) left out 
	    ENDDO
	  ENDDO
	 ENDDO
	ELSE ! hindered settling
	 DO i=0,i1
	  DO j=0,j1
	    DO k=k1b,k1e ! at k=0 wsed should be zero at kmax wsed is calculated
	      ctot=0.
	      sum_c_ws=0.
	      kp = MIN(k+1,k1)
	      kpp = MIN(k+1,k1-1)
	      km = MAX(k-1,0)
	      DO n=1,nfrac
	        ccc(n) = 0.5*(csed2(n,i,j,k) + csed2(n,i,j,kp))
		ctot=ccc(n)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,kp))+ctot
		! for mud Cfloc must be used to calculate Ctot (Cfloc=Ctot*dfloc/dpart)
		! ccc(n) is used in drift flux correction sum_c_ws which is calculated with mass concentration based on Ctot (not on Cfloc)
	      ENDDO
	      ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
	      DO n=1,nfrac
!		IF (frac(n)%dfloc>frac(n)%dpart) THEN
!			ws(n)=0.002*MIN(5.,ccc(n)*frac(n)%rho) !Whitehouse simple formula for flocculation
!			ws(n)=MAX(0.001,ws(n))
!			ws(n)=frac(n)%ws_dep+MAX(Rp(i)*cos_u(j)-schuif_x,0.)/300.*(frac(n)%ws-frac(n)%ws_dep) ! try time/x-dist dependent settling velocity
!			ws(n)=ws(n)*(1.-ctot)**(frac(n)%n-1.)
!		ELSE
	        ws(n)=ws_basis(n)*(1.-ctot)**(frac(n)%n-1.)
!		ENDIF
		! ws is positive downward, wsed is positive upward 
		! Ri_Za is defined with volume concentration ctot
!		sum_c_ws=sum_c_ws+ws(n)*0.5*(csed(n,i,j,k)*frac(n)%rho/rr(i,j,k)+csed(n,i,j,kp)*frac(n)%rho/rr(i,j,kp))
  		sum_c_ws=sum_c_ws+ws(n)*ccc(n)*frac(n)%rho/(0.5*(rr2(i,j,k)+rr2(i,j,kp)))
		! According to drift velocity literature the drift flux is calculated using the mass-fraction in stead of volume-fraction,
		! therefore an extra frac(n)%rho/rho_mix is included
	      ENDDO
	      DO n=1,nfrac
		wsed(n,i,j,k)=Wcfd(i,j,k)+sum_c_ws-ws(n) ! wsed is positive upward
	      ENDDO


!	      DO iter=1,1
!	        ! TVD interpolation csed for correct return current in drift flux:		     
!		! calculate all again with new TVD interpolation of csed and use first calculation as first guess of direction Wsed
!	      ctot=0.
!	      sum_c_ws=0.
!	      DO n=1,nfrac 
!	        noemer = csed(n,i,j,kp)-csed(n,i,j,k)
!		noemer = MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!		rrpos = (csed(n,i,j,k)-csed(n,i,j,km))/noemer
!		rrneg = (csed(n,i,j,kpp+1)-csed(n,i,j,kp))/noemer
!		if (Wsed(n,i,j,k)>0.) then
!		  ccc(n)=csed(n,i,j,k )+0.5*limiter(rrpos)*(csed(n,i,j,kp)-csed(n,i,j,k ))*(1.-dt_dzi*ABS(Wsed(n,i,j,k)))
!		else
!		  ccc(n)=csed(n,i,j,kp)+0.5*limiter(rrneg)*(csed(n,i,j,k )-csed(n,i,j,kp))*(1.-dt_dzi*ABS(Wsed(n,i,j,k)))
!		endif
!		ctot=ccc(n)*frac(n)%dfloc/frac(n)%dpart+ctot      
!	      ENDDO
!	      ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
!	      DO n=1,nfrac
!	        ws(n)=frac(n)%ws*(1.-ctot)**(frac(n)%n-1.)
!		sum_c_ws=sum_c_ws+ws(n)*ccc(n)*frac(n)%rho/(0.5*(rr(i,j,k)+rr(i,j,kp)))
!	      ENDDO
!	      DO n=1,nfrac
!		wsed(n,i,j,k)=Wcfd(i,j,k)+sum_c_ws-ws(n) ! wsed is positive upward
!	      ENDDO
!	      ENDDO
	      Wfluid(i,j,k)=sum_c_ws !Wcfd(i,j,k) left out 
	    ENDDO
	  ENDDO
	 ENDDO
	ENDIF

	IF (k1b.eq.1.and.k1e.eq.(k1-1)) THEN
	   DO i=0,i1
	     DO j=0,j1
	       DO n=1,nfrac
	    	wsed(n,i,j,kbed(i,j))=Wcfd(i,j,kbed(i,j))  ! prevent sediment to flow through the bed or bed-obstacles
		Wfluid(i,j,kbed(i,j))=0.
		wsed(n,i,j,k1)=Wcfd(i,j,k1) ! prevent sediment to flow out of free surface
		wsed(n,i,j,k1-1)=-wsed(n,i,j,k1-1)*MIN(0.,frac(n)%ws/(ABS(frac(n)%ws)+1.e-12))  !limit wsed at zero for fractions with downward settling velocity--> no transport of sediment in/out domain at top, but air may escape
	   	Wfluid(i,j,k1)=0.
		Wfluid(i,j,k1-1)=0.
	       ENDDO
	     ENDDO
	   ENDDO
	ENDIF

      END SUBROUTINE slipvelocity

       subroutine erosion_deposition(ccnew,cbotnew,ucfd,vcfd,wcfd,rcfd,ccfd,cbotcfd,ddt,dz)

       implicit none

	real     wsed(nfrac,0:i1,0:j1,0:k1),erosion,deposition,uu,vv,absU,z0,ust,yplus,tau !local variables
	integer n,tel,kplus,n1
	REAL     ccnew(nfrac,0:i1,0:j1,0:k1),cbotnew(nfrac,0:i1,0:j1)  ! output
	REAL     ccfd(nfrac,0:i1,0:j1,0:k1),cbotcfd(nfrac,0:i1,0:j1)  ! input
	REAL	 ddt,dz ! input
	REAL     ucfd(0:i1,0:j1,0:k1),vcfd(0:i1,0:j1,0:k1),wcfd(0:i1,0:j1,0:k1),rcfd(0:i1,0:j1,0:k1) ! input
	REAL     wfluid(0:i1,0:j1,0:k1) !dummy
	REAL     cbottot,kn_sed_avg,Mr_avg,tau_e_avg,ctot_firstcel,cbotnewtot
	REAL PSD_bot_sand_massfrac(nfr_sand),PSD_bed_sand_massfrac(nfr_sand),PSD_sand(0:nfr_sand),factor,d50
	REAL cbottot_sand,mbottot_sand,cbedtot,mbedtot_sand,diameter_sand_PSD(0:nfr_sand)
	REAL cbedtot_sand,rho_botsand,rho_bedsand,rho_sand,delta,Dstar,Shields_cr,ustc2,TT,phip,erosion_sand
	REAL erosionf(nfrac),depositionf(nfrac),erosion_avg(nfrac)
	erosion=0.
	deposition=0.


			
	
	IF (nobst>0.or.bedlevelfile.ne.''.or.interaction_bed.eq.4.or.interaction_bed.eq.5) THEN
		call slipvelocity(ccfd,wcfd,wsed,rcfd,0,k1-1,wfluid,ddt,dz) 
		!determine wsed on top of obstacles (actually everywhere in the domain) --> Wsed is not zero on kbed(i,j) in this manner!
	ELSE
		call slipvelocity(ccfd,wcfd,wsed,rcfd,0,0,wfluid,ddt,dz)
	ENDIF

	!write(*,*),'kbed(1,1)=',kbed(1,1)
!	DO n=1,nfrac
!	   Wsed(n,:,:,0)=Wsed(n,:,:,0)+0.5*(Wnew(:,:,1)-Wnew(:,:,0)) !correct for average Wnew at k=1/2 (location of C(:,:,1))
!	ENDDO
!	write(*,*),'Wnew(10,10,0),3xCnew,3xwsed',rank,Wnew(10,10,0),cnew(1:3,10,10,0),wsed(1:3,10,10,0)
	
	IF (interaction_bed.le.3) THEN ! erosion sedimentation without bed update and for each sediment fraction independently
		DO i=1,imax
		  DO j=1,jmax
			!! First Ubot_TSHD and Vbot_TSHD is subtracted to determine tau 
			!! only over ambient velocities not over U_TSHD
			kplus = MIN(kbed(i,j)+1,k1)
			uu=0.5*(ucfd(i,j,kplus)+ucfd(i-1,j,kplus))-Ubot_TSHD(j)
			vv=0.5*(vcfd(i,j,kplus)+vcfd(i,j-1,kplus))-Vbot_TSHD(j)
			absU=sqrt((uu)**2+(vv)**2)
			DO n=1,nfrac
				ust=0.1*absU
				do tel=1,10 ! 10 iter is more than enough
					z0=frac(n)%kn_sed/30.+0.11*nu_mol/MAX(ust,1.e-9) 
					! for tau shear on sediment don't use kn (which is result of bed ripples), use frac(n)%kn_sed 
					! it is adviced to use kn_sed=dfloc
					ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
				enddo
				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
						if (yplus<30.) then
						  do tel=1,10 ! 10 iter is more than enough
								yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
								ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU
						  enddo
						endif
						if (yplus<5.) then !viscous sublayer uplus=yplus
								ust=sqrt(absU*nu_mol/(0.5*dz))
						endif
		!               if (yplus<11.225) then  !viscous sublayer uplus=yplus
		!                       ust=sqrt(absU*nu_mol/(0.5*dz))
		!               endif
				kplus = MIN(kbed(i,j)+1,k1)
				tau=rcfd(i,j,kplus)*ust*ust  

					 ! called before update in correc (in last k3 step of RK3) so Cnewbot and Cnew are used to determine interaction with bed, 
					 ! but effect is added to dcdtbot and dcdt to make superposition on other terms already included in dcdt
					 ! dcdtbot contains sediment volume concentration in bed [-] with ghost bed-cells of dz deep
					 ! dcdt contains sediment volume concentration in water-cells [-] 
				 IF (interaction_bed.ge.2) THEN
					  erosion = frac(n)%M*MAX(0.,(tau/frac(n)%tau_e-1.))*ddt/frac(n)%rho ! m3/m2
				 ENDIF
				 IF (interaction_bed.eq.3) THEN
					   erosion = MIN(erosion,cbotcfd(n,i,j)*dz) ! m3/m2, not more material can be eroded than there was in the bed
				 ENDIF
				 kplus = MIN(kbed(i,j)+1,k1)
				 deposition = MAX(0.,(1.-tau/frac(n)%tau_d))*ccfd(n,i,j,kplus)*MIN(0.,Wsed(n,i,j,kbed(i,j)))*ddt ! m
				 ccnew(n,i,j,kplus)=ccnew(n,i,j,kplus)+(erosion+deposition)/(dz) ! vol conc. [-]
				 cbotnew(n,i,j)=cbotnew(n,i,j)-(erosion+deposition)/(dz) ! vol conc. [-]
			ENDDO
		  ENDDO
		ENDDO
	ELSEIF(interaction_bed.eq.4.or.interaction_bed.eq.5) THEN  ! including bedupdate; erosion based on avg sediment properties top layer
		DO i=1,imax
		  DO j=1,jmax 
			erosionf=0.
			depositionf=0.
			cbotnewtot=0.	
			ctot_firstcel=0.
			!! First Ubot_TSHD and Vbot_TSHD is subtracted to determine tau 
			!! only over ambient velocities not over U_TSHD
			kplus = MIN(kbed(i,j)+1,k1)
			uu=0.5*(ucfd(i,j,kplus)+ucfd(i-1,j,kplus))-Ubot_TSHD(j)
			vv=0.5*(vcfd(i,j,kplus)+vcfd(i,j-1,kplus))-Vbot_TSHD(j)
			absU=sqrt((uu)**2+(vv)**2)			
			cbottot=0.
			cbedtot=0.
			IF (nfr_silt>0) THEN			
				!! 1 determine erosion/sedimentation of mixture of all silt fractions
				DO n1=1,nfr_silt
					n=nfrac_silt(n1)
					cbottot=cbottot+cbotcfd(n,i,j)
					cbedtot=cbedtot+Clivebed(n,i,j,kbed(i,j))
				ENDDO
				
				kn_sed_avg=0.
				Mr_avg=0.
				tau_e_avg=0.
				IF (cbottot>0.) THEN
					DO n1=1,nfr_silt
						n=nfrac_silt(n1)
						kn_sed_avg=kn_sed_avg+(cbotcfd(n,i,j)/cbottot)*frac(n)%kn_sed
						Mr_avg=Mr_avg+(cbotcfd(n,i,j)/cbottot)*frac(n)%M/frac(n)%rho
						tau_e_avg=tau_e_avg+(cbotcfd(n,i,j)/cbottot)*frac(n)%tau_e
					ENDDO
				ELSEIF (cbedtot>0.) THEN ! determine avg sediment characteristics from bed:
					DO n1=1,nfr_silt
						n=nfrac_silt(n1)
						kn_sed_avg=kn_sed_avg+(Clivebed(n,i,j,kbed(i,j))/cbedtot)*frac(n)%kn_sed
						Mr_avg=Mr_avg+(Clivebed(n,i,j,kbed(i,j))/cbedtot)*frac(n)%M/frac(n)%rho
						tau_e_avg=tau_e_avg+(Clivebed(n,i,j,kbed(i,j))/cbedtot)*frac(n)%tau_e
					ENDDO
				ELSE !arbitrarily choose silt frac 1, no effect as there is no erosion possible
					kn_sed_avg=frac(nfrac_silt(1))%kn_sed
					Mr_avg=frac(nfrac_silt(1))%M
					tau_e_avg=frac(nfrac_silt(1))%tau_e
				ENDIF
				
				ust=0.1*absU
				do tel=1,10 ! 10 iter is more than enough
					z0=kn_sed_avg/30.+0.11*nu_mol/MAX(ust,1.e-9) 
					! for tau shear on sediment don't use kn (which is result of bed ripples), use frac(n)%kn_sed 
					! it is adviced to use kn_sed=dfloc
					ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
				enddo
				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
				if (yplus<30.) then
				  do tel=1,10 ! 10 iter is more than enough
						yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
						ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU
				  enddo
				endif
				if (yplus<5.) then !viscous sublayer uplus=yplus
						ust=sqrt(absU*nu_mol/(0.5*dz))
				endif
				kplus = MIN(kbed(i,j)+1,k1)
				tau=rcfd(i,j,kplus)*ust*ust  
				DO n1=1,nfr_silt
					n=nfrac_silt(n1)
					erosion_avg(n) = Mr_avg*MAX(0.,(tau/tau_e_avg-1.))*ddt ! m3/m2	 erosion_avg is filled for silt fractions only with silt erosion			
					IF (cbottot>0.) THEN
						erosionf(n) = erosion_avg(n) * (cbotcfd(n,i,j)/cbottot) !erosion per fraction
					ELSEIF (cbedtot>0.) THEN
						erosionf(n) = erosion_avg(n) * (Clivebed(n,i,j,kbed(i,j))/cbedtot) !erosion per fraction
					ELSE
						erosionf(n) = 0.
					ENDIF
					erosionf(n) = MIN(erosionf(n),(cbotcfd(n,i,j)+Clivebed(n,i,j,kbed(i,j)))*dz) ! m3/m2, not more material can be eroded than there was in top layer cbotcfd+c in top cel bed
					kplus = MIN(kbed(i,j)+1,k1)
					ust=0.1*absU !re-calculate tau with kn_sed for deposition as it is not dependent on avg dpart in mixture
					do tel=1,10 ! 10 iter is more than enough
						z0=frac(n)%kn_sed/30.+0.11*nu_mol/MAX(ust,1.e-9) 
						! for tau shear on sediment don't use kn (which is result of bed ripples), use frac(n)%kn_sed 
						! it is adviced to use kn_sed=dfloc
						ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
					enddo
					yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
					if (yplus<30.) then
					  do tel=1,10 ! 10 iter is more than enough
							yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
							ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU
					  enddo
					endif
					if (yplus<5.) then !viscous sublayer uplus=yplus
							ust=sqrt(absU*nu_mol/(0.5*dz))
					endif
					tau=rcfd(i,j,kplus)*ust*ust  !for deposition apply tau belonging to own frac(n)%kn_sed
	 
					 depositionf(n) = MAX(0.,(1.-tau/frac(n)%tau_d))*ccfd(n,i,j,kplus)*MIN(0.,Wsed(n,i,j,kbed(i,j)))*ddt ! m --> dep is negative due to negative wsed
					 ccnew(n,i,j,kplus)=ccnew(n,i,j,kplus)+(erosionf(n)+depositionf(n))/(dz) ! vol conc. [-]
					 cbotnew(n,i,j)=cbotnew(n,i,j)-(erosionf(n)+depositionf(n))/(dz) ! vol conc. [-]
					 cbotnewtot=cbotnewtot+cbotnew(n,i,j)
					 ctot_firstcel=ccnew(n,i,j,kplus)+ctot_firstcel
				ENDDO
			ENDIF
			!! 2 determine erosion/sedimentation of mixture of all sand fractions
			IF (nfr_sand>0) THEN
				cbottot_sand=0.
				cbedtot_sand=0.
				mbottot_sand=0.
				mbedtot_sand=0.
				rho_botsand=0.
				rho_bedsand=0.
				PSD_sand=0.
				DO n1=1,nfr_sand
					n=nfrac_sand(n1)
					IF (n1<nfr_sand) THEN
						diameter_sand_PSD(n1)=0.5*(frac(nfrac_sand(n1))%dpart+frac(nfrac_sand(n1+1))%dpart)
					ELSEIF (n1>1) THEN
						diameter_sand_PSD(n1)=frac(nfrac_sand(n1))%dpart+( frac(nfrac_sand(n1))%dpart-diameter_sand_PSD(n1-1) )
					ELSE
						diameter_sand_PSD(n1)=frac(n)%dpart
					ENDIF
					IF (nfr_sand>1) THEN
						diameter_sand_PSD(0)=frac(nfrac_sand(1))%dpart+( frac(nfrac_sand(1))%dpart-diameter_sand_PSD(1) )
					ELSE
						diameter_sand_PSD(0)=0.
					ENDIF					
					cbottot_sand=cbottot_sand+cbotcfd(n,i,j)
					cbedtot_sand=cbedtot_sand+Clivebed(n,i,j,kbed(i,j))
					mbottot_sand=mbottot_sand+cbotcfd(n,i,j)*frac(n)%rho
					mbedtot_sand=mbedtot_sand+Clivebed(n,i,j,kbed(i,j))*frac(n)%rho
					PSD_bot_sand_massfrac(n1)=mbottot_sand
					PSD_bed_sand_massfrac(n1)=mbedtot_sand
					rho_botsand=rho_botsand+cbotcfd(n,i,j)*frac(n)%rho
					rho_bedsand=rho_bedsand+Clivebed(n,i,j,kbed(i,j))*frac(n)%rho
				ENDDO
				IF (mbottot_sand>0.) THEN
					PSD_sand(1:nfr_sand)=PSD_bot_sand_massfrac/mbottot_sand
					rho_sand=rho_botsand/cbottot_sand
				ELSEIF (mbedtot_sand>0.) THEN
					PSD_sand(1:nfr_sand)=PSD_bed_sand_massfrac/mbedtot_sand
					rho_sand=rho_bedsand/cbedtot_sand
				ELSE
					PSD_sand=0. ! there is no sand, hence there will be no erosion either therefore choice in PSD and resulting d50 is not important
					rho_sand=frac(nfrac_sand(1))%rho
				ENDIF
				IF (nfr_sand>1) THEN
				! determine d50 based on different input fractions
				! PSD and d50 is calculated under the assumption that each concentration is valid in the range of dpart  AVG(d_n-1;d_n)..AVG(d_n+1;d_n)
				! for the first and last fraction the assumption is that the provided dpart is exactly the middle of the particle size range of this bin
				! Example A: dpart=15,100,200mu and conc=0.06,0.3,0.24 gives PSD=0,0.1,0.6,1 for particle size class edges -27.5,57.5,150,250mu  --> d50=131.5 mu
				! Please note that the first negative number is only for correct interpolation, never will it result in a negative d50.
				! Example B: dpart=dpart=15,100,200mu and conc=0.55,0.05,0.0 gives PSD=0,0.9167,1,1 for particle size class edges -27.5,57.5,150,250mu  --> d50=18.9 mu
				! Example C: dpart=dpart=15,100,200mu and conc=0.05,0.0,0.55 gives PSD=0,0.0833,0.0833,1 for particle size class edges -27.5,57.5,150,250mu  --> d50=195 mu
				! Example D: dpart=dpart=15,100,200mu and conc=0.00,0.0,0.6 gives PSD=0,0,0,1 for particle size class edges -27.5,57.5,150,250mu  --> d50=200 mu
				! Example E: dpart=dpart=15,100,200mu and conc=0.60,0.0,0.0 gives PSD=1,0,0,0 for particle size class edges -27.5,57.5,150,250mu  --> d50=15 mu
				! Example F: dpart=dpart=15,100,200mu and conc=0.00,0.60,0 gives PSD=0,0,1,1 for particle size class edges -27.5,57.5,150,250mu  --> d50=103.75 mu (this is centre of middle PSD, which complies with the median d50 definition however it is not exactly 100 mu because of choice in not perfectly symmetrical bins)
				! Example G: dpart=dpart=15,100,185mu and conc=0.00,0.60,0 gives PSD=0,0,1,1 for particle size class edges -27.5,57.5,142.5,227.5mu  --> d50=100 mu (this new input dpart is symmetric!)

					DO 	n1=1,nfr_sand
						IF (PSD_sand(n1).ge.0.5) THEN
!							IF (n1.eq.1) THEN
!								d50=frac(nfrac_sand(n1))%dpart
!								EXIT
!							ELSE
								factor=(PSD_sand(n1)-0.5)/(PSD_sand(n1)-PSD_sand(n1-1))
								!d50=factor*frac(nfrac_sand(n1-1))%dpart+(1.-factor)*frac(nfrac_sand(n1))%dpart
								d50=factor*diameter_sand_PSD(n1-1)+(1.-factor)*diameter_sand_PSD(n1)
								EXIT
!							ENDIF
						ELSE ! concentration is zero, no d50 found, assign d50 as largest dpart, doesn't matter there can't be erosion anyway:
							d50=frac(nfrac_sand(nfr_sand))%dpart
						ENDIF
					ENDDO
				ELSE
					d50=frac(nfrac_sand(1))%dpart
				ENDIF
				kn_sed_avg=2.*d50 !  kn=2*d50 is mentioned in VanRijn1984 paper, the pickup function which is applied here, however elsewhere vanRijn mentions larger kn_sed like 6*d50...)
				ust=0.1*absU
				do tel=1,10 ! 10 iter is more than enough
					z0=kn_sed_avg/30.+0.11*nu_mol/MAX(ust,1.e-9) 
					ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
				enddo
				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
				if (yplus<30.) then
				  do tel=1,10 ! 10 iter is more than enough
						yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
						ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU
				  enddo
				endif
				if (yplus<5.) then !viscous sublayer uplus=yplus
						ust=sqrt(absU*nu_mol/(0.5*dz))
				endif
				kplus = MIN(kbed(i,j)+1,k1)
				delta = (rho_sand-rho_b)/rho_b !chosen to apply rho_b and not local rcfd because then delta can become very low with very sediment saturated flow
				Dstar = d50 * ((delta*ABS(gz))/nu_mol**2)**(0.333333333333333)
				Shields_cr = 0.3/(1.+1.2*Dstar)+0.055*(1.-exp(-0.02*Dstar))  !Soulsby and Whitehouse 1997 curve through original Shields for threshold of motion sediment particles, Soulsy book Eq. SC(77)
				ustc2 = Shields_cr * ABS(gz)*delta*d50
				
				TT = (ust*ust-ustc2)/(MAX(ustc2,1.e-12))
				TT = MAX(TT,0.) !TT must be positive
				phip = 0.00033*Dstar**0.3*TT**1.5   ! VanRijn1984 pickup function
				
				DO n1=1,nfr_sand
					n=nfrac_sand(n1)			
					erosion_avg(n) = phip * (delta*ABS(gz)*d50)**0.5*ddt  !*rho_sand/rho_sand ! erosion flux in kg/m2/(kg/m3)= m3/m2=m
					IF (cbottot_sand>0.) THEN
						erosionf(n) = erosion_avg(n) * (cbotcfd(n,i,j)/cbottot_sand) !erosion per fraction
					ELSEIF (cbedtot_sand>0.) THEN
						erosionf(n) = erosion_avg(n) * (Clivebed(n,i,j,kbed(i,j))/cbedtot_sand) !erosion per fraction
					ELSE
						erosionf(n) = 0.
					ENDIF
					erosionf(n) = MIN(erosionf(n),(cbotcfd(n,i,j)+Clivebed(n,i,j,kbed(i,j)))*dz) ! m3/m2, not more material can be eroded than there was in top layer cbotcfd+c in top cel bed
					depositionf(n) = ccfd(n,i,j,kplus)*MIN(0.,Wsed(n,i,j,kbed(i,j)))*ddt ! m --> dep is negative due to negative wsed
					ccnew(n,i,j,kplus)=ccnew(n,i,j,kplus)+(erosionf(n)+depositionf(n))/(dz) ! vol conc. [-]
					cbotnew(n,i,j)=cbotnew(n,i,j)-(erosionf(n)+depositionf(n))/(dz) ! vol conc. [-]
					cbotnewtot=cbotnewtot+cbotnew(n,i,j)
					ctot_firstcel=ccnew(n,i,j,kplus)+ctot_firstcel
				ENDDO
			ENDIF				
!update bedlevel for combined erosion deposition silt plus sand fractions: 
			IF (interaction_bed.eq.4) THEN
				IF (cbotnewtot.lt.0.and.(kbed(i,j)-1).ge.0.and.SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed) THEN 
				!add half cell sediment on cbot account for further erosion without lowering 1 dz yet (because otherwise flipflop between ero-1dz and depo+1dz)
					kplus = MIN(kbed(i,j)+1,k1)
					DO n=1,nfrac ! also cbotnew(n,i,j) is le 0:
						cbotnew(n,i,j)=cbotnew(n,i,j)+0.5*Clivebed(n,i,j,kbed(i,j)) 
						Clivebed(n,i,j,kbed(i,j))=0.5*Clivebed(n,i,j,kbed(i,j)) 
					ENDDO			
				ELSEIF (cbotnewtot.lt.0.and.(kbed(i,j)-1).ge.0) THEN !erosion of 1 layer dz:
					kplus = MIN(kbed(i,j)+1,k1)
					drdt(i,j,kbed(i,j))=rho_b
					rnew(i,j,kbed(i,j))=rho_b
					rold(i,j,kbed(i,j))=rho_b					
					DO n=1,nfrac ! also cbotnew(n,i,j) is le 0:
						cbotnew(n,i,j)=cbotnew(n,i,j)+Clivebed(n,i,j,kbed(i,j)) !top layer is now previous top cel bed minus (erosion>top-layer)
						ccnew(n,i,j,kbed(i,j))=(erosionf(n)+depositionf(n))/(dz) !assign all erosion+depo as new sediment concentration new bottom layer fluid
						ccnew(n,i,j,kplus)=ccnew(n,i,j,kplus)-(erosionf(n)+depositionf(n))/(dz) !remove erosion+depo from previous bottom layer fluid
						Clivebed(n,i,j,kbed(i,j))=0. ! not bed anymore but fluid
						drdt(i,j,kbed(i,j)) = drdt(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						rnew(i,j,kbed(i,j)) = rnew(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						rold(i,j,kbed(i,j)) = rold(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
					ENDDO
					!kbed(i,j)=MAX(kbed(i,j)-1,0)  !update bed level at end		
					kbed(i,j)=kbed(i,j)-1
				ELSEIF (ctot_firstcel.ge.cfixedbed.and.kbed(i,j)+1.le.kmax 
     &				.and.(SUM(erosionf(1:nfrac))+SUM(depositionf(1:nfrac))).lt.0.) THEN ! sedimentation of 1 layer dz each time because ctot in fluid already above threshold of bed, only if erosion is less than deposition::
					kbed(i,j)=kbed(i,j)+1
					DO n=1,nfrac 
						Clivebed(n,i,j,kbed(i,j))=ccnew(n,i,j,kbed(i,j))+
     &						(cfixedbed-ctot_firstcel)*ccnew(n,i,j,kbed(i,j))/MAX(ctot_firstcel,1.e-12)  ! reduce ccnew to arrive at ctot = cfixedbed
						cbotnew(n,i,j)=cbotnew(n,i,j)-(cfixedbed-ctot_firstcel)*ccnew(n,i,j,kbed(i,j))/MAX(ctot_firstcel,1.e-12)
						ccnew(n,i,j,kbed(i,j))=Clivebed(n,i,j,kbed(i,j)) ! to have correct density in this cell, this doesn't give incorrect sediment balance as soon as a bed cel becomes fluid again then ccnew is restarted with (ero+depo) and this old concentration is forgotten
					ENDDO
				ELSEIF ((cbotnewtot+ctot_firstcel).ge.cfixedbed.and.kbed(i,j)+1.le.kmax.and.cbotnewtot.gt.1.e-12
     &				.and.(SUM(erosionf(1:nfrac))+SUM(depositionf(1:nfrac))).lt.0.) THEN ! sedimentation of 1 layer dz each time, only if erosion is less than deposition:
					kbed(i,j)=kbed(i,j)+1
					!kbed(i,j)=MIN(kbed(i,j)+1,kmax) !update bed level at start sedimentation 
					DO n=1,nfrac 
						Clivebed(n,i,j,kbed(i,j))=ccnew(n,i,j,kbed(i,j))+
     &						(cfixedbed-ctot_firstcel)*cbotnew(n,i,j)/MAX(cbotnewtot,1.e-12)  ! apply sedimentation ratio between fractions new sediment concentration of cells within bed
						cbotnew(n,i,j)=cbotnew(n,i,j)-(cfixedbed-ctot_firstcel)*cbotnew(n,i,j)/MAX(cbotnewtot,1.e-12)
						ccnew(n,i,j,kbed(i,j))=Clivebed(n,i,j,kbed(i,j)) ! to have correct density in this cell, this doesn't give incorrect sediment balance as soon as a bed cel becomes fluid again then ccnew is restarted with (ero+depo) and this old concentration is forgotten
					ENDDO
				ENDIF
			ENDIF
		  ENDDO
		ENDDO
	ENDIF		
	
	IF (interaction_bed.eq.4) THEN
		call bound_cbot_integer(kbed) ! apply correct boundary conditions for updated kbed
	ENDIF
	

	
      END SUBROUTINE erosion_deposition


       subroutine erosion_deposition_old

       implicit none

	real     wsed(nfrac,0:i1,0:j1,0:k1),erosion,deposition,uu,vv,absU,z0,ust,yplus,tau 
	integer n,tel
	REAL     wdrift(0:i1,0:j1,0:k1),rdrift(0:i1,0:j1,0:k1),wfluid(0:i1,0:j1,0:k1)

	erosion=0.
	deposition=0.
	call slipvelocity(cnew,Wnew,wsed,rnew,0,0,wfluid,0.,1.) !made dt=1 and dz=1
!	DO n=1,nfrac
!	   Wsed(n,:,:,0)=Wsed(n,:,:,0)+0.5*(Wnew(:,:,1)-Wnew(:,:,0)) !correct for average Wnew at k=1/2 (location of C(:,:,1))
!	ENDDO
!	write(*,*),'Wnew(10,10,0),3xCnew,3xwsed',rank,Wnew(10,10,0),cnew(1:3,10,10,0),wsed(1:3,10,10,0)

	DO i=1,imax
	  DO j=1,jmax
	    !! First Ubot_TSHD and Vbot_TSHD is subtracted to determine tau 
	    !! only over ambient velocities not over U_TSHD
  	    uu=Unew(i,j,1)-Ubot_TSHD(j)
	    vv=Vnew(i,j,1)-Vbot_TSHD(j)
	    absU=sqrt((uu)**2+(vv)**2)
	    DO n=1,nfrac
		ust=0.1*absU
		do tel=1,10 ! 10 iter is more than enough
			z0=frac(n)%kn_sed/30.+0.11*nu_mol/MAX(ust,1.e-9) 
			! for tau shear on sediment don't use kn (which is result of bed ripples), use frac(n)%kn_sed 
			! it is adviced to use kn_sed=dfloc
			ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
		enddo
		yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
                if (yplus<30.) then
                  do tel=1,10 ! 10 iter is more than enough
                        yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
                        ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU
                  enddo
                endif
                if (yplus<5.) then !viscous sublayer uplus=yplus
                        ust=sqrt(absU*nu_mol/(0.5*dz))
                endif
!               if (yplus<11.225) then  !viscous sublayer uplus=yplus
!                       ust=sqrt(absU*nu_mol/(0.5*dz))
!               endif
		tau=rnew(i,j,1)*ust*ust  

	      	 ! called before update in correc (in last k3 step of RK3) so Cnewbot and Cnew are used to determine interaction with bed, 
	      	 ! but effect is added to dcdtbot and dcdt to make superposition on other terms already included in dcdt
	      	 ! dcdtbot contains sediment volume concentration in bed [-] with ghost bed-cells of dz deep
	     	 ! dcdt contains sediment volume concentration in water-cells [-] 
		 IF (interaction_bed.ge.2) THEN
	      	  erosion = frac(n)%M*MAX(0.,(tau/frac(n)%tau_e-1.))*dt ! kg/m2
		 ENDIF
		 IF (interaction_bed.eq.3) THEN
	      	   erosion = MIN(erosion,Cnewbot(n,i,j)*frac(n)%rho*dz) ! kg/m2, not more material can be eroded than there was in the bed
		 ENDIF
	     	 deposition = MAX(0.,(1.-tau/frac(n)%tau_d))*Cnew(n,i,j,1)*frac(n)%rho*MIN(0.,Wsed(n,i,j,0))*dt ! kg/m2
	      	 dCdt(n,i,j,1)=dCdt(n,i,j,1)+(erosion+deposition)/(frac(n)%rho*dz) ! vol conc. [-]
	      	 dCdtbot(n,i,j)=dCdtbot(n,i,j)-(erosion+deposition)/(frac(n)%rho*dz) ! vol conc. [-]

!		 if (i.eq.10.and.j.eq.10) then
!	write(*,*),'rank,n,i,j,absU,tau,Cnew,Wsed,ero,dep',rank,n,i,j,absU,tau,Cnew(n,i,j,1),Wsed(n,i,j,0),
!     & frac(n)%M*MAX(0.,(tau/frac(n)%tau_e-1.)),deposition
!		 endif
!		 if (deposition>0.) then
!			write(*,*),'n,i,j,absU,tau,Cnew,Wsed,deposition',n,i,j,absU,tau,Cnew(n,i,j,1),Wsed(n,i,j,0),deposition
!		 endif
			
	    ENDDO
	  ENDDO
	ENDDO

	
      END SUBROUTINE erosion_deposition_old

       subroutine air_bubbles_free_surface

       implicit none

	integer n,t
	REAL wsed(nfrac,0:i1,0:j1,0:k1),Wfluid(0:i1,0:j1,0:k1)

	call slipvelocity(cnew,Wnew,wsed,rnew,kmax,kmax,wfluid,0.,1.) !made dt=0,dz=1
	DO i=1,imax
	  DO j=1,jmax
		DO n=1,nfrac
		  dCdt(n,i,j,kmax)=dCdt(n,i,j,kmax)+Cnew(n,i,j,kmax)*dt*Wsed(n,i,j,kmax)*MIN(0.,frac(n)%ws/(ABS(frac(n)%ws)+1.e-12))/dz 
		  ! air bubbles dissappear at free surface [-]
		ENDDO
	  ENDDO
	ENDDO
       !! Set boundary conditions jet in:
      	DO t=1,tmax_inPpunt
	  i=i_inPpunt(t)
	  j=j_inPpunt(t)
	  DO n=1,nfrac
		dCdt(n,i,j,kmax)=dCdt(n,i,j,kmax)-Cnew(n,i,j,kmax)*dt*Wsed(n,i,j,kmax)*MIN(0.,frac(n)%ws/(ABS(frac(n)%ws)+1.e-12))/dz ! undo for jet
      	  ENDDO
	ENDDO

	end subroutine air_bubbles_free_surface

      END MODULE sediment
