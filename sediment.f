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
	if (nobst>0) then
	 DO i=0,i1
	  DO j=0,j1
	    DO n=1,nfrac
		kplus = MIN(kbed(i,j)+1,k1)
	  	csed2(n,i,j,kbed(i,j))=csed2(n,i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
		rr2(i,j,kbed(i,j))=rr2(i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
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
		wsed(n,i,j,k1-1)=-wsed(n,i,j,k1-1)*MIN(0.,frac(n)%ws/(ABS(frac(n)%ws)+1.e-12))  !limit wsed at zero for fractions with downward settling velocity
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
	integer n,tel,kplus
	REAL     ccnew(nfrac,0:i1,0:j1,0:k1),cbotnew(nfrac,0:i1,0:j1)  ! output
	REAL     ccfd(nfrac,0:i1,0:j1,0:k1),cbotcfd(nfrac,0:i1,0:j1)  ! input
	REAL	 ddt,dz ! input
	REAL     ucfd(0:i1,0:j1,0:k1),vcfd(0:i1,0:j1,0:k1),wcfd(0:i1,0:j1,0:k1),rcfd(0:i1,0:j1,0:k1) ! input
	REAL     wfluid(0:i1,0:j1,0:k1) !dummy

	erosion=0.
	deposition=0.
	
	IF (nobst>0) THEN
		call slipvelocity(ccfd,wcfd,wsed,rcfd,0,k1-1,wfluid,ddt,dz) 
		!determine wsed on top of obstacles (actually everywhere in the domain)
	ELSE
		call slipvelocity(ccfd,wcfd,wsed,rcfd,0,0,wfluid,ddt,dz)
	ENDIF

	
!	DO n=1,nfrac
!	   Wsed(n,:,:,0)=Wsed(n,:,:,0)+0.5*(Wnew(:,:,1)-Wnew(:,:,0)) !correct for average Wnew at k=1/2 (location of C(:,:,1))
!	ENDDO
!	write(*,*),'Wnew(10,10,0),3xCnew,3xwsed',rank,Wnew(10,10,0),cnew(1:3,10,10,0),wsed(1:3,10,10,0)

	DO i=1,imax
	  DO j=1,jmax
	    !! First Ubot_TSHD and Vbot_TSHD is subtracted to determine tau 
	    !! only over ambient velocities not over U_TSHD
		kplus = MIN(kbed(i,j)+1,k1)
  	    uu=ucfd(i,j,kplus)-Ubot_TSHD(j)
	    vv=vcfd(i,j,kplus)-Vbot_TSHD(j)
	    absU=sqrt((uu)**2+(vv)**2)
	    DO n=1,nfrac
		ust=0.1*absU
		do tel=1,10 ! 10 iter is more than enough
			z0=frac(n)%kn_sed/30.+0.11*nu_mol/MAX(ust,1.e-9) 
			! for tau shear on sediment don't use kn (which is result of bed ripples), use frac(n)%kn_sed 
			! it is adviced to use kn_sed=dfloc
			ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
		enddo
		yplus=0.5*dz*ust/nu_mol
                if (yplus<30.) then
                  do tel=1,10 ! 10 iter is more than enough
                        yplus=0.5*dz*ust/nu_mol
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
		yplus=0.5*dz*ust/nu_mol
                if (yplus<30.) then
                  do tel=1,10 ! 10 iter is more than enough
                        yplus=0.5*dz*ust/nu_mol
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
