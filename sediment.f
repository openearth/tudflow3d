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
		USE netcdf
       implicit none

 	INTEGER n,n2
	REAL Re_p,atm_pres,rhoair_z,zz,ddzz,z_rks(1:kmax),interpseries,gvector
	integer :: ncid, rhVarId, status2, ndims, xtype,natts,status
	integer, dimension(nf90_max_var_dims) :: dimids
	integer size1,size2,size3,size4
	

	DO n=1,nfrac
	  Re_p=ABS(frac(n)%ws)*frac(n)%dfloc/(ekm_mol/rho_b)
	  
	  ! hindered_settling = 1	!Hindered settling formula [-] 1=Rowe (1987) (smooth Ri-Za org) (default); 2=Garside (1977); 3=Di Felice (1999)
	  IF (hindered_settling.eq.1) THEN
!	hindered settling function acc. to Rowe (1987), 
!	is the smoothed representation of the original Richardson and Zaki (1954) relations:
!	determined for 0.04<c<0.55 and 0.001<Rep<3e4
	  frac(n)%n=(4.7+0.41*Re_p**0.75)/(1.+0.175*Re_p**0.75) -1. 
	  ELSEIF (hindered_settling.eq.2) THEN
!	hindered settling function acc. to Garside (1977), 
!	higher value of n than Rowe (1987)
!	determined for 0.04<c<0.55 and 0.2<Rep<1e3
	  frac(n)%n=(5.1+0.27*Re_p**0.9)/(1.+0.1*Re_p**0.9)-1.
	  ELSEIF (hindered_settling.eq.3) THEN
!	hindered settling function acc. to Di Felice (1999), 
!	much higher value of n than Rowe (1987) and Garside (1977)
!	determined for 0<c<0.05 and 0.01<Rep<1e3
	  frac(n)%n=(6.5+0.3*Re_p**0.74)/(1.+0.1*Re_p**0.74)-1.
	  ENDIF
	  IF (rank.eq.0) THEN
	  write(*,*),'fraction #,ws,Re_p,n-factor hindered settling:',n,frac(n)%ws,Re_p,frac(n)%n+1.
	  ENDIF
	ENDDO

	rhocorr_air_z=1.
	atm_pres=101325. !standard atmospheric pressure in Pa 
	gvector=sqrt(gx**2+gy**2+gz**2)
	DO n=1,nair !correction for comppressible air
		DO k=0,k1
			zz=k*dz-0.5*dz
			ddzz=zz-(depth-frac(nfrac_air(n))%zair_ref_belowsurf) !positive upward and defined wrt zair_ref_belowsurf
			rhoair_z = ((frac(nfrac_air(n))%zair_ref_belowsurf-ddzz)*rho_b*gvector+atm_pres) / 
     &			       ((frac(nfrac_air(n))%zair_ref_belowsurf)*rho_b*gvector+atm_pres) * frac(nfrac_air(n))%rho
			rhocorr_air_z(nfrac_air(n),k) = frac(nfrac_air(n))%rho / rhoair_z
		ENDDO
	ENDDO
	
	DO k=1,kmax
		z_rks(k)=k*dz-0.5*dz
	ENDDO
	IF (av_slope_z(1)<0.) THEN
	  av_slope(1:imax,1:jmax,0:k1)=avalanche_slope(1)
	ELSE
		IF (av_slope_z(1)>0.) CALL writeerror(141)
		n=1
		DO WHILE (av_slope_z(n+1)>0.)
		  n=n+1
		  IF (av_slope_z(n)-av_slope_z(n-1)<0.) CALL writeerror(141)
		END DO
		IF (av_slope_z(n)<depth) CALL writeerror(141)
		IF (avalanche_slope(n)<0.) CALL writeerror(142)
		n2=0
		DO k=1,kmax
			av_slope(1:imax,1:jmax,k)=interpseries(av_slope_z(1:n),avalanche_slope(1:n),n2,z_rks(k))
		ENDDO
		av_slope(1:imax,1:jmax,0)=av_slope(1:imax,1:jmax,1)
		av_slope(1:imax,1:jmax,k1)=av_slope(1:imax,1:jmax,kmax)
	ENDIF
	
	if (avfile.ne.'') then
		status2 = nf90_open(avfile, nf90_NoWrite, ncid) 
		IF (status2/= nf90_noerr) THEN
			write(*,*),'file =',avfile
			CALL writeerror(703)
		ENDIF
		status = nf90_inq_varid(ncid, "av_slope",rhVarid)
		if (status.eq.nf90_NoErr) then
			call check( nf90_inquire_variable(ncid, rhVarid, dimids = dimIDs))
			call check( nf90_inquire_dimension(ncid, dimIDs(1), len = size1))
			call check( nf90_inquire_dimension(ncid, dimIDs(2), len = size2))
			call check( nf90_inquire_dimension(ncid, dimIDs(3), len = size3))
			IF(size1.ne.imax.or.size2.ne.jmax*px.or.size3.ne.k1+1) CALL writeerror(704)
			call check( nf90_get_var(ncid,rhVarid,av_slope(1:imax,1:jmax,0:k1),
     &                       start=(/1,rank*jmax+1,1/),count=(/imax,jmax,k1+1/)) )
		else
			write(*,*),'file =',avfile,' variable "av_slope" not found '
			CALL writeerror(704)
		endif
		call check( nf90_close(ncid) )
	endif
	

			
	!IF (rank.eq.0) write(*,*),av_slope,z_rks,av_slope_z(1:n),avalanche_slope(1:n)
	
	
      END SUBROUTINE init_sediment

      SUBROUTINE slipvelocity(csed,Wcfd,wsed,rr,k1b,k1e,sumWkm,dt,dz)

	implicit none

	INTEGER n,kp ! local variables!
	REAL    wsed(nfrac,0:i1,0:j1,0:k1)
	REAL    csed(nfrac,0:i1,0:j1,0:k1),csed2(nfrac,0:i1,0:j1,0:k1)
	REAL    Wcfd(0:i1,0:j1,0:k1),sumWkm(0:i1,0:j1,0:k1)
	REAL	rr(0:i1,0:j1,0:k1),rr2(0:i1,0:j1,0:k1)
	REAL	ctot,sum_c_ws,ws(nfrac),ccc(nfrac)
	INTEGER k1b,k1e,kpp,km,iter,kplus,kplus2
	REAL dt_dzi,noemer,rrpos,rrneg,limiter,dt,dz,ws_basis(nfrac),W_km_sum

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

	csed2=cW !csed
	IF (applyVOF.eq.1) THEN
		call state_edges(cW,rhW)
	ENDIF 
	rr2=rhW
	IF (applyVOF.eq.1) THEN
		rhW=rho_b 
	ENDIF 
	if (nobst>0.or.bedlevelfile.ne.''.or.interaction_bed.ge.4) then ! default C(k=0)=C(k=1); therefore only for cases when bed is not necessarily at k=0 this fix is needed:
	 DO i=0,i1
	  DO j=0,j1
		kplus = MIN(kbed(i,j)+1,k1)
		kplus2 = MIN(kbed(i,j)+2,k1)
		!rr2(i,j,kbed(i,j))=rr2(i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
		rr2(i,j,kbed(i,j))=rr(i,j,kplus) ! make rhW(kbed) equal to rho fluid first cell above to get correct drift flux settling
		!rr2(i,j,kbed(i,j))=1.5*rr(i,j,kplus)-0.5*rr(i,j,kplus2)
		DO n=1,nfrac
			!csed2(n,i,j,kbed(i,j))=csed2(n,i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
			csed2(n,i,j,kbed(i,j))=csed(n,i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
			!csed2(n,i,j,kbed(i,j))=1.5*csed(n,i,j,kplus)-0.5*csed(n,i,j,kplus2)
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
	        ccc(n) = csed2(n,i,j,k) !0.5*(csed2(n,i,j,k) + csed2(n,i,j,kp))	
			ctot=ccc(n)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,kp))+ctot
	        ws(n)=ws_basis(n)
		! ws is positive downward, wsed is positive upward 
			sum_c_ws=sum_c_ws+ws(n)*ccc(n)*frac(n)%rho/rr2(i,j,k) !(0.5*(rr2(i,j,k)+rr2(i,j,kp)))
!!		! According to drift velocity literature the drift flux is calculated using the mass-fraction in stead of volume-fraction,
!!		! because sum_c_ws is slipvelocity relative to mixture velocity not relative to volumetric flux!
		! therefore an extra frac(n)%rho/rho_mix is included
	      ENDDO
		  ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
	      DO n=1,nfrac
			wsed(n,i,j,k)=Wcfd(i,j,k)+sum_c_ws-ws(n) ! wsed is positive upward
	      ENDDO
		  W_km_sum=0.
		  do n=1,nfrac
			W_km_sum=W_km_sum+ccc(n)*frac(n)%rho*(wsed(n,i,j,k)-Wcfd(i,j,k)) 
		  enddo 
		  W_km_sum=W_km_sum+(1.-ctot)*rho_b*sum_c_ws !sum_c_ws=fluid return velocity --> slipvelocity of fluid relative to mixture-velocity
		  sumWkm(i,j,k)=W_km_sum !NEW 2-10-2018: contains sum of all fractions and fluid phase drift velocity for correct determination driftflux-force
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
	      DO n=1,nfrac
	        ccc(n) = csed2(n,i,j,k) !0.5*(csed2(n,i,j,k) + csed2(n,i,j,kp))
			ctot=ccc(n)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,kp))+ctot
		! for mud Cfloc must be used to calculate Ctot (Cfloc=Ctot*dfloc/dpart)
		! ccc(n) is used in drift flux correction sum_c_ws which is calculated with mass concentration based on Ctot (not on Cfloc)
	      ENDDO
	      ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
	      DO n=1,nfrac
	        ws(n)=ws_basis(n)*(1.-ctot)**(frac(n)%n) !frac(n)%n lowered with one already
		! ws is positive downward, wsed is positive upward 
		! Ri_Za is defined with volume concentration ctot
  		sum_c_ws=sum_c_ws+ws(n)*ccc(n)*frac(n)%rho/rr2(i,j,k) !(0.5*(rr2(i,j,k)+rr2(i,j,kp)))
!!		! According to drift velocity literature the drift flux is calculated using the mass-fraction in stead of volume-fraction,
!!		! because sum_c_ws is slipvelocity relative to mixture velocity not relative to volumetric flux!
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
!	        ws(n)=frac(n)%ws*(1.-ctot)**(frac(n)%n) !frac(n)%n lowered with one already
!		sum_c_ws=sum_c_ws+ws(n)*ccc(n)*frac(n)%rho/(0.5*(rr(i,j,k)+rr(i,j,kp)))
!	      ENDDO
!	      DO n=1,nfrac
!		wsed(n,i,j,k)=Wcfd(i,j,k)+sum_c_ws-ws(n) ! wsed is positive upward
!	      ENDDO
!	      ENDDO
			 W_km_sum=0.
		  	 do n=1,nfrac
				W_km_sum=W_km_sum+ccc(n)*frac(n)%rho*(wsed(n,i,j,k)-Wcfd(i,j,k)) 
			 enddo 
			 W_km_sum=W_km_sum+(1.-ctot)*rho_b*sum_c_ws !sum_c_ws=fluid return velocity --> slipvelocity of fluid relative to mixture-velocity
		     sumWkm(i,j,k)=W_km_sum !NEW 2-10-2018: contains sum of all fractions and fluid phase drift velocity for correct determination driftflux-force
	    ENDDO
	  ENDDO
	 ENDDO
	ENDIF

	IF (k1b.eq.1.and.k1e.eq.(k1-1)) THEN
	   DO i=0,i1
	     DO j=0,j1
	       DO n=1,nfrac
	    	wsed(n,i,j,kbed(i,j))=Wcfd(i,j,kbed(i,j))  ! prevent sediment to flow through the bed or bed-obstacles
		wsed(n,i,j,k1)=Wcfd(i,j,k1) ! prevent sediment to flow out of free surface
		wsed(n,i,j,k1-1)=-wsed(n,i,j,k1-1)*MIN(0.,frac(n)%ws/(ABS(frac(n)%ws)+1.e-12))  !limit wsed at zero for fractions with downward settling velocity--> no transport of sediment in/out domain at top, but air may escape
	       ENDDO
	     ENDDO
	   ENDDO
	   Wsed(:,:,:,0)=0.
	ENDIF

      END SUBROUTINE slipvelocity
	  
      SUBROUTINE slipvelocity_bed(csed,Wcfd,wsed,rr,sumWkm,dt,dz)

	implicit none

	INTEGER n,kp ! local variables!
	REAL    wsed(nfrac,0:i1,0:j1,0:k1)
	REAL    csed(nfrac,0:i1,0:j1,0:k1),csed2(nfrac,0:i1,0:j1,0:k1)
	REAL    Wcfd(0:i1,0:j1,0:k1),sumWkm(0:i1,0:j1,0:k1)
	REAL	rr(0:i1,0:j1,0:k1),rr2(0:i1,0:j1,0:k1)
	REAL	ctot,sum_c_ws,ws(nfrac),ccc(nfrac)
	INTEGER kpp,km,iter,kplus,kplus2
	REAL dt_dzi,noemer,rrpos,rrneg,limiter,dt,dz,ws_basis(nfrac),W_km_sum

	dt_dzi=dt/dz 

	  DO n=1,nfrac
	    ws_basis(n)=frac(n)%ws_dep
	  ENDDO

	csed2=cW !csed
	IF (applyVOF.eq.1) THEN
		call state_edges(cW,rhW)
	ENDIF 
	rr2=rhW
	IF (applyVOF.eq.1) THEN
		rhW=rho_b 
	ENDIF 
	if (nobst>0.or.bedlevelfile.ne.''.or.interaction_bed.ge.4) then ! default C(k=0)=C(k=1); therefore only for cases when bed is not necessarily at k=0 this fix is needed:
	 DO i=0,i1
	  DO j=0,j1
		kplus = MIN(kbed(i,j)+1,k1)
		kplus2 = MIN(kbed(i,j)+2,k1)
		!rr2(i,j,kbed(i,j))=rr2(i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
		rr2(i,j,kbed(i,j))=rr(i,j,kplus) ! make rhW(kbed) equal to rho fluid first cell above to get correct drift flux settling
		!rr2(i,j,kbed(i,j))=1.5*rr(i,j,kplus)-0.5*rr(i,j,kplus2)
		DO n=1,nfrac
			!csed2(n,i,j,kbed(i,j))=csed2(n,i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
			csed2(n,i,j,kbed(i,j))=csed(n,i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
			!csed2(n,i,j,kbed(i,j))=1.5*csed(n,i,j,kplus)-0.5*csed(n,i,j,kplus2)
	    ENDDO		
	  ENDDO
	 ENDDO
	endif


	IF (slipvel.eq.2) THEN ! no hindered settling
	 DO i=0,i1
	  DO j=0,j1
	    k=kbed(i,j)
	    !DO k=k1b,k1e ! at k=0 wsed should be zero at kmax wsed is calculated
	      ctot=0.
	      sum_c_ws=0.
	      kp = MIN(k+1,k1)
	      DO n=1,nfrac
	        ccc(n) = csed2(n,i,j,k) !0.5*(csed2(n,i,j,k) + csed2(n,i,j,kp))	
			ctot=ccc(n)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,kp))+ctot
	        ws(n)=ws_basis(n)
		! ws is positive downward, wsed is positive upward 
			sum_c_ws=sum_c_ws+ws(n)*ccc(n)*frac(n)%rho/rr2(i,j,k) !(0.5*(rr2(i,j,k)+rr2(i,j,kp)))
!!		! According to drift velocity literature the drift flux is calculated using the mass-fraction in stead of volume-fraction,
!!		! because sum_c_ws is slipvelocity relative to mixture velocity not relative to volumetric flux!
		! therefore an extra frac(n)%rho/rho_mix is included
	      ENDDO
		  ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
	      DO n=1,nfrac
			wsed(n,i,j,k)=Wcfd(i,j,k)+sum_c_ws-ws(n) ! wsed is positive upward
	      ENDDO
		  W_km_sum=0.
		  do n=1,nfrac
			W_km_sum=W_km_sum+ccc(n)*frac(n)%rho*(wsed(n,i,j,k)-Wcfd(i,j,k)) 
		  enddo 
		  W_km_sum=W_km_sum+(1.-ctot)*rho_b*sum_c_ws !sum_c_ws=fluid return velocity --> slipvelocity of fluid relative to mixture-velocity
		  sumWkm(i,j,k)=W_km_sum !NEW 2-10-2018: contains sum of all fractions and fluid phase drift velocity for correct determination driftflux-force
	    !ENDDO
	  ENDDO
	 ENDDO
	ELSE ! hindered settling
	 DO i=0,i1
	  DO j=0,j1
	    !DO k=k1b,k1e ! at k=0 wsed should be zero at kmax wsed is calculated
		k=kbed(i,j)
	      ctot=0.
	      sum_c_ws=0.
	      kp = MIN(k+1,k1)
	      DO n=1,nfrac
	        ccc(n) = csed2(n,i,j,k) !0.5*(csed2(n,i,j,k) + csed2(n,i,j,kp))
			ctot=ccc(n)*frac(n)%dfloc/frac(n)%dpart*0.5*(rhocorr_air_z(n,k)+rhocorr_air_z(n,kp))+ctot
		! for mud Cfloc must be used to calculate Ctot (Cfloc=Ctot*dfloc/dpart)
		! ccc(n) is used in drift flux correction sum_c_ws which is calculated with mass concentration based on Ctot (not on Cfloc)
	      ENDDO
	      ctot=MIN(ctot,1.) ! limit on 1, see also Winterwerp 1999 p.46, because Cfloc can be >1
	      DO n=1,nfrac
	        ws(n)=ws_basis(n)*(1.-ctot)**(frac(n)%n) !frac(n)%n lowered with one already
		! ws is positive downward, wsed is positive upward 
		! Ri_Za is defined with volume concentration ctot
  		sum_c_ws=sum_c_ws+ws(n)*ccc(n)*frac(n)%rho/rr2(i,j,k) !(0.5*(rr2(i,j,k)+rr2(i,j,kp)))
!!		! According to drift velocity literature the drift flux is calculated using the mass-fraction in stead of volume-fraction,
!!		! because sum_c_ws is slipvelocity relative to mixture velocity not relative to volumetric flux!
		! therefore an extra frac(n)%rho/rho_mix is included
	      ENDDO
	      DO n=1,nfrac
		wsed(n,i,j,k)=Wcfd(i,j,k)+sum_c_ws-ws(n) ! wsed is positive upward 
	      ENDDO
			 W_km_sum=0.
		  	 do n=1,nfrac
				W_km_sum=W_km_sum+ccc(n)*frac(n)%rho*(wsed(n,i,j,k)-Wcfd(i,j,k)) 
			 enddo 
			 W_km_sum=W_km_sum+(1.-ctot)*rho_b*sum_c_ws !sum_c_ws=fluid return velocity --> slipvelocity of fluid relative to mixture-velocity
		     sumWkm(i,j,k)=W_km_sum !NEW 2-10-2018: contains sum of all fractions and fluid phase drift velocity for correct determination driftflux-force
	    !ENDDO
	  ENDDO
	 ENDDO
	ENDIF


      END SUBROUTINE slipvelocity_bed

	  

       subroutine erosion_deposition(ccnew,cbotnew,ucfd,vcfd,wcfd,rcfd,ccfd,cbotcfd,ddt,dz)

       implicit none

	   include 'mpif.h'
      integer ierr
	real     wsed(nfrac,0:i1,0:j1,0:k1),erosion,deposition,uu,vv,absU,z0,ust,yplus,tau !local variables
	integer n,tel,kplus,n1,k2,kk
	REAL     ccnew(nfrac,0:i1,0:j1,0:k1),cbotnew(nfrac,0:i1,0:j1)  ! output
	REAL     ccfd(nfrac,0:i1,0:j1,0:k1),cbotcfd(nfrac,0:i1,0:j1)  ! input
	REAL	 ddt,dz ! input
	REAL     ucfd(0:i1,0:j1,0:k1),vcfd(0:i1,0:j1,0:k1),wcfd(0:i1,0:j1,0:k1),rcfd(0:i1,0:j1,0:k1) ! input
	REAL     sumWkm(0:i1,0:j1,0:k1) !dummy
	REAL     cbottot,kn_sed_avg,Mr_avg,tau_e_avg,ctot_firstcel,cbotnewtot,cbotnewtot_pos
	REAL PSD_bot_sand_massfrac(nfr_sand),PSD_bed_sand_massfrac(nfr_sand),PSD_sand(0:nfr_sand),factor,d50
	REAL cbottot_sand,mbottot_sand,cbedtot,mbedtot_sand,diameter_sand_PSD(0:nfr_sand)
	REAL cbedtot_sand,rho_botsand,rho_bedsand,rho_sand,delta,Dstar,Shields_cr,ustc2,TT,phipp,erosion_sand
	REAL ws_botsand,ws_bedsand,ws_sand,Re_p,CD,Shields_eff,dubdt,Fd,Fi,Fl,W
	REAL erosionf(nfrac),depositionf(nfrac),erosion_avg(nfrac)
	REAL zb_all(0:i1,0:j1),maxbedslope(0:i1,0:j1),sl1,sl2,sl3,sl4,sl5,sl6,sl7,sl8
	REAL d_cbotnew(nfrac,0:i1,0:j1),dbed,dl,dbed_allowed,dbed_adjust,dz_botlayer,c_adjust,c_adjustA,c_adjustB
	REAL*8 cbf(nfrac,0:i1),cbb(nfrac,0:i1),zbf(0:i1),zbb(0:i1),reduced_sed
	INTEGER itrgt,jtrgt,nav,n_av,kplus2,kpp
	REAL ws_botsand2,rho_botsand2,mbottot_sand2,PSD_bot_sand_massfrac2(nfr_sand),have_avalanched,have_avalanched_tmp,cctot
	REAL ccfdtot_firstcel,wsedbed,distance_to_bed,zb_W,gvector
	REAL pickup_random(1:imax,1:jmax),vs,ve,Uhor(1:imax,1:jmax,1:kmax),cbed,uu2,vv2,bed_slope,facx,facy,bs_geo
	REAL qb,MMM,MME,flux,ucr,qbf(1:nfrac),uuRrel,uuLrel,vvRrel,vvLrel,Shields,absUbl,ve_check
      integer clock,nnn,k_maxU
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed	
	  INTEGER kbedp(0:i1,0:j1),ibeg,iend 
	
	erosion=0.
	deposition=0.

	
!	IF (nobst>0.or.bedlevelfile.ne.''.or.interaction_bed.eq.4.or.interaction_bed.eq.6) THEN
!		call slipvelocity(ccfd,wcfd,wsed,rcfd,0,k1-1,sumWkm,ddt,dz) 
!		!determine wsed on top of obstacles (actually everywhere in the domain) --> Wsed is not zero on kbed(i,j) in this manner!
!	ELSE
!		call slipvelocity(ccfd,wcfd,wsed,rcfd,0,0,sumWkm,ddt,dz)
!	ENDIF
	!call slipvelocity_bed(ccfd,wcfd,wsed,rcfd,sumWkm,ddt,dz) 
	call slipvelocity_bed(ccfd,0.*wcfd,wsed,rcfd,sumWkm,ddt,dz) 
	
	DO i=0,i1
	  DO j=0,j1
		kbedp(i,j)=MIN(kbed(i,j)+2,k1) 
	  ENDDO 
	ENDDO 
	IF (interaction_bed.ge.4) THEN
	 DO i=1,imax
	  DO j=1,jmax
		zbed(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(cbotcfd(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
	  ENDDO 
	 ENDDO 
	 call bound_cbot(zbed)
	ENDIF 
				  
	d_cbotnew = 0.	  
	IF (interaction_bed.le.3.and.time_n.ge.tstart_morf) THEN ! erosion sedimentation without bed update and for each sediment fraction independently
		DO i=1,imax
		  DO j=1,jmax
			!! First Ubot_TSHD and Vbot_TSHD is subtracted to determine tau 
			!! only over ambient velocities not over U_TSHD
			IF (IBMorder.eq.2) THEN
				zb_W=zbed(i,j)
				kpp=MIN(CEILING(zb_W/dz+0.5)+1,k1)		!kpp is between 1*dz-2*dz distance from bed 	
				!kpp=MIN(CEILING(zb_W/dz+0.5),k1)		!kpp is between 0-dz distance from bed 	
				distance_to_bed=(REAL(kpp)-0.5)*dz-zb_W
				!as start use first cell (0.5-1)*dz distance from bed 
!				IF (distance_to_bed<0.5*dz) THEN !first cell too close to bed, therefore use second cell (1-1.5)*dz distance from bed 
!				  kpp=MIN(CEILING(zb_W/dz+0.5)+1,k1)		!kpp is in principle between 1*dz-2*dz distance from bed, but due to this if-statement only 1-1.5 from bed
!				  distance_to_bed=(REAL(kpp)-0.5)*dz-zb_W
!				ENDIF 
			ELSE 
				kpp = MIN(kbedp(i,j),k1)                !kpp is 1.5*dz from 0-order ibm bed
				distance_to_bed=1.5*dz				
			ENDIF
			IF (kbedp(i,j)<kbedp(i+1,j).and.kbedp(i,j)<kbedp(i-1,j)) THEN ! pit
				uu=0.5*(ucfd(i,j,kpp)+ucfd(i-1,j,kpp))-Ubot_TSHD(j)
			ELSE 
				uu=0.5*(ucfd(i,j,MAX(kbedp(i,j),kbedp(i+1,j),kpp))+ucfd(i-1,j,MAX(kbedp(i,j),kbedp(i-1,j),kpp)))-Ubot_TSHD(j) !choice to not alter distance_to_bed at slopes
			ENDIF
			IF (kbedp(i,j)<kbedp(i,j+1).and.kbedp(i,j)<kbedp(i,j-1)) THEN ! pit
				vv=0.5*(vcfd(i,j,kpp)+vcfd(i,j-1,kpp))-Vbot_TSHD(j)
			ELSE 
				vv=0.5*(vcfd(i,j,MAX(kbedp(i,j),kbedp(i,j+1),kpp))+vcfd(i,j-1,MAX(kbedp(i,j),kbedp(i,j-1),kpp)))-Vbot_TSHD(j) !choice to not alter distance_to_bed at slopes
			ENDIF
			IF (pickup_bedslope_geo.eq.1) THEN 
				bed_slope = atan((zbed(i+1,j)-zbed(i-1,j))/(Rp(i+1)-Rp(i-1)))
				uu2 = uu*cos(bed_slope)+wcfd(i,j,kpp)*sin(bed_slope)
				facx = 1./cos(bed_slope)
				bed_slope = atan((zbed(i,j+1)-zbed(i,j-1))/(Rp(i)*(phip(j+1)-phip(j-1))))
				vv2 = vv*cos(bed_slope)+wcfd(i,j,kpp)*sin(bed_slope)
				facy = 1./cos(bed_slope)
				bs_geo = facx*facy ! increase in dx and dy (area) over which pickup and deposition take place
				absU = sqrt((uu2)**2+(vv2)**2)	
			ELSE 
				absU=sqrt((uu)**2+(vv)**2)	
				bs_geo = 1.
			ENDIF 
			ust=0.1*absU
			do tel=1,10 ! 10 iter is more than enough
				z0=MAX(kn/30.+0.11*nu_mol/MAX(ust,1.e-9),1e-9) 
				ust=absU/MAX(1./kappa*log(MAX(distance_to_bed/z0,1.001)),2.) !ust maximal 0.5*absU
			enddo
			distance_to_bed=z_tau_sed
			absU=ust/kappa*log(distance_to_bed/z0) ! replace absU with velocity that is valid at z_tau_sed from bed (user input to make result less dependent of grid resolution)

			DO n=1,nfrac
				ust=0.1*absU
				do tel=1,10 ! 10 iter is more than enough
					z0=frac(n)%kn_sed/30.+0.11*nu_mol/MAX(ust,1.e-9) 
					! for tau shear on sediment don't use kn (which is result of bed ripples), use frac(n)%kn_sed; it is adviced to use kn_sed=dfloc
					ust=absU/MAX(1./kappa*log(distance_to_bed/z0),2.) !ust maximal 0.5*absU
				enddo

				kplus = MIN(kbed(i,j)+1,k1)
				kplus2 = MIN(kbed(i,j)+2,k1)
				tau=rho_b*ust*ust  

					 ! called before update in correc (in last k3 step of RK3) so Cnewbot and Cnew are used to determine interaction with bed, 
					 ! but effect is added to dcdtbot and dcdt to make superposition on other terms already included in dcdt
					 ! dcdtbot contains sediment volume concentration in bed [-] with ghost bed-cells of dz deep
					 ! dcdt contains sediment volume concentration in water-cells [-] 
				 IF (interaction_bed.ge.2) THEN
					  erosion = frac(n)%M*MAX(0.,(tau/frac(n)%tau_e-1.))*ddt/frac(n)%rho*bs_geo ! m3/m2
				 ENDIF
!				 erosion = erosion -MIN(0.,ddt*0.5*(Diffcof(i,j,kplus)+Diffcof(i,j,kplus2))*(ccnew(n,i,j,kplus2)-ccnew(n,i,j,kplus))/dz)  !add upward diffusion to erosion flux
				 IF (interaction_bed.eq.3) THEN
					   erosion = MIN(erosion,cbotnew(n,i,j)*dz) ! m3/m2, not more material can be eroded than there was in the bed
				 ENDIF
				 IF (depo_implicit.eq.1) THEN  !determine deposition as sink implicit
					 kplus = MIN(kbed(i,j)+1,k1) 
					 ccnew(n,i,j,kplus)=(ccnew(n,i,j,kplus)+erosion/dz)/
     &					 (1.-MAX(0.,(1.-tau/frac(n)%tau_d))*MIN(0.,Wsed(n,i,j,kbed(i,j)))*ddt/dz) ! vol conc. [-]
					 deposition = MAX(0.,(1.-tau/frac(n)%tau_d))*ccnew(n,i,j,kplus)*MIN(0.,Wsed(n,i,j,kbed(i,j)))*ddt ! m !ccfd
					 cbotnew(n,i,j)=cbotnew(n,i,j)-(erosion+deposition)/(dz) ! vol conc. [-]				 
				 ELSE
					 kplus = MIN(kbed(i,j)+1,k1)
					 deposition = MAX(0.,(1.-tau/frac(n)%tau_d))*ccnew(n,i,j,kplus)*MIN(0.,Wsed(n,i,j,kbed(i,j)))*ddt ! m !ccfd
					 ccnew(n,i,j,kplus)=ccnew(n,i,j,kplus)+(erosion+deposition)/(dz) ! vol conc. [-]
					 cbotnew(n,i,j)=cbotnew(n,i,j)-(erosion+deposition)/(dz) ! vol conc. [-]
			     ENDIF
			ENDDO
		  ENDDO
		ENDDO
	ELSEIF((interaction_bed.ge.4).and.time_n.ge.tstart_morf) THEN  ! including bedupdate; erosion based on avg sediment properties top layer
		IF (pickup_fluctuations.eq.1) THEN
			!1 add white noise to pickup
		  CALL SYSTEM_CLOCK(COUNT=clock)
		  CALL RANDOM_SEED(size = nnn)
		  ALLOCATE(seed(nnn))
		  CALL SYSTEM_CLOCK(COUNT=clock)
		  seed = clock + 37 * (/ (i - 1, i = 1, nnn) /)
		  CALL RANDOM_SEED(PUT = seed)			
		  call random_number(pickup_random) ! uniform distribution 0,1
 		  pickup_random=1.+2.*(pickup_random-0.5)*pickup_fluctuations_ampl
		ENDIF	
		IF (cbed_method.eq.2) THEN
			Uhor(1:imax,1:jmax,1:kmax)=sqrt((0.5*(ucfd(0:imax-1,1:jmax,1:kmax)+ucfd(1:imax,1:jmax,1:kmax)))**2 + 
     &    				(0.5*(vcfd(1:imax,0:jmax-1,1:kmax)+vcfd(1:imax,1:jmax,1:kmax)))**2)	
		ENDIF 
		
		IF (periodicx.eq.0.and.monopile<0) THEN
		  ibeg=2
		  iend=imax-1
		ELSEIF (periodicx.eq.0.and.monopile>0) THEN
		  ibeg=1
		  iend=imax-1 
		ELSE 
		  ibeg=1
		  iend=imax 
		ENDIF 
		DO i=ibeg,iend
		  DO j=1,jmax 
			erosionf=0.
			depositionf=0.
			cbotnewtot=0.
			cbotnewtot_pos=0.
			ctot_firstcel=0.
			!! First Ubot_TSHD and Vbot_TSHD is subtracted to determine tau 
			!! only over ambient velocities not over U_TSHD
			IF (IBMorder.eq.2) THEN
				zb_W=zbed(i,j)
				kpp=MIN(CEILING(zb_W/dz+0.5)+1,k1)		!kpp is between 1*dz-2*dz distance from bed 	
				!kpp=MIN(CEILING(zb_W/dz+0.5),k1)		!kpp is between 0-dz distance from bed 	
				distance_to_bed=(REAL(kpp)-0.5)*dz-zb_W
				!as start use first cell (0.5-1)*dz distance from bed 
!				IF (distance_to_bed<0.5*dz) THEN !first cell too close to bed, therefore use second cell (1-1.5)*dz distance from bed 
!				  kpp=MIN(CEILING(zb_W/dz+0.5)+1,k1)		!kpp is in principle between 1*dz-2*dz distance from bed, but due to this if-statement only 1-1.5 from bed
!				  distance_to_bed=(REAL(kpp)-0.5)*dz-zb_W
!				ENDIF 
			ELSE 
				kpp = MIN(kbedp(i,j),k1)                !kpp is 1.5*dz from 0-order ibm bed
				distance_to_bed=1.5*dz				
			ENDIF
			IF (kbedp(i,j)<kbedp(i+1,j).and.kbedp(i,j)<kbedp(i-1,j)) THEN ! pit
				uu=0.5*(ucfd(i,j,kpp)+ucfd(i-1,j,kpp))-Ubot_TSHD(j)
			ELSE 
				uu=0.5*(ucfd(i,j,MAX(kbedp(i,j),kbedp(i+1,j),kpp))+ucfd(i-1,j,MAX(kbedp(i,j),kbedp(i-1,j),kpp)))-Ubot_TSHD(j) !choice to not alter distance_to_bed at slopes
			ENDIF
			IF (kbedp(i,j)<kbedp(i,j+1).and.kbedp(i,j)<kbedp(i,j-1)) THEN ! pit
				vv=0.5*(vcfd(i,j,kpp)+vcfd(i,j-1,kpp))-Vbot_TSHD(j)
			ELSE 
				vv=0.5*(vcfd(i,j,MAX(kbedp(i,j),kbedp(i,j+1),kpp))+vcfd(i,j-1,MAX(kbedp(i,j),kbedp(i,j-1),kpp)))-Vbot_TSHD(j) !choice to not alter distance_to_bed at slopes
			ENDIF
			uuRrel = ucfd(i  ,j,MAX(kbedp(i,j),kbedp(i+1,j),kpp))-Ubot_TSHD(j)  !no correction for pit because uuR from cell i always needs to be same as uuL from i+1 in bedload otherwise interuption and pit may never fill up
			uuLrel = ucfd(i-1,j,MAX(kbedp(i,j),kbedp(i-1,j),kpp))-Ubot_TSHD(j)			
			vvRrel = vcfd(i,j  ,MAX(kbedp(i,j),kbedp(i,j+1),kpp))-Vbot_TSHD(j)
			vvLrel = vcfd(i,j-1,MAX(kbedp(i,j),kbedp(i,j-1),kpp))-Vbot_TSHD(j)			
			uuR_relax(i,j)  = bl_relax*uuRrel+(1.-bl_relax)*uuR_relax(i,j)  	!needed for bedload-fluxes
			uuL_relax(i,j)  = bl_relax*uuLrel+(1.-bl_relax)*uuL_relax(i,j)
			vvR_relax(i,j)  = bl_relax*vvRrel+(1.-bl_relax)*vvR_relax(i,j)
			vvL_relax(i,j)  = bl_relax*vvLrel+(1.-bl_relax)*vvL_relax(i,j)
			absUbl = MAX(sqrt((0.5*(uuR_relax(i,j)+uuL_relax(i,j)))**2+(0.5*(vvR_relax(i,j)+vvL_relax(i,j)))**2),1.e-6)
			uuRrel = uuR_relax(i,j)/absUbl
			uuLrel = uuL_relax(i,j)/absUbl
			vvRrel = vvR_relax(i,j)/absUbl
			vvLrel = vvL_relax(i,j)/absUbl

			IF (pickup_bedslope_geo.eq.1) THEN 
				bed_slope = atan((zbed(i+1,j)-zbed(i-1,j))/(Rp(i+1)-Rp(i-1)))
				uu2 = uu*cos(bed_slope)+wcfd(i,j,kpp)*sin(bed_slope)
				facx = 1./cos(bed_slope)
				bed_slope = atan((zbed(i,j+1)-zbed(i,j-1))/(Rp(i)*(phip(j+1)-phip(j-1))))
				vv2 = vv*cos(bed_slope)+wcfd(i,j,kpp)*sin(bed_slope)
				facy = 1./cos(bed_slope)
				bs_geo = facx*facy ! increase in dx and dy (area) over which pickup and deposition take place
				absU = sqrt((uu2)**2+(vv2)**2)	
			ELSE 
				absU=sqrt((uu)**2+(vv)**2)	
				bs_geo = 1.
			ENDIF 
			ust=0.1*absU
			do tel=1,10 ! 10 iter is more than enough
				z0=MAX(kn/30.+0.11*nu_mol/MAX(ust,1.e-9),1e-9) 
				ust=absU/MAX(1./kappa*log(MAX(distance_to_bed/z0,1.001)),2.) !ust maximal 0.5*absU
			enddo
			distance_to_bed=z_tau_sed
			absU=ust/kappa*log(distance_to_bed/z0) ! replace absU with velocity that is valid at z_tau_sed from bed (user input to make result less dependent of grid resolution)			
			cbottot=0.
			cbedtot=0.
			IF (nfr_silt>0) THEN			
				!! 1 determine erosion/sedimentation of mixture of all silt fractions
				DO n1=1,nfr_silt
					n=nfrac_silt(n1)
					cbottot=cbottot+MAX(cbotnew(n,i,j),0.)
					cbedtot=cbedtot+Clivebed(n,i,j,kbed(i,j))
				ENDDO
				
				kn_sed_avg=0.
				Mr_avg=0.
				tau_e_avg=0.
				IF ((interaction_bed.ge.6.and.kbed(i,j).eq.0).or.interaction_bed.eq.7) THEN !unlimited erosion in case kbed.eq.0
					DO n1=1,nfr_silt
						n=nfrac_silt(n1)
						kn_sed_avg=kn_sed_avg+(c_bed(n)/cfixedbed)*frac(n)%kn_sed
						Mr_avg=Mr_avg+(c_bed(n)/cfixedbed)*frac(n)%M/frac(n)%rho
						tau_e_avg=tau_e_avg+(c_bed(n)/cfixedbed)*frac(n)%tau_e
					ENDDO					
				ELSEIF (cbottot>0.) THEN
					DO n1=1,nfr_silt
						n=nfrac_silt(n1)
						kn_sed_avg=kn_sed_avg+(MAX(cbotnew(n,i,j),0.)/cbottot)*frac(n)%kn_sed
						Mr_avg=Mr_avg+(MAX(cbotnew(n,i,j),0.)/cbottot)*frac(n)%M/frac(n)%rho
						tau_e_avg=tau_e_avg+(MAX(cbotnew(n,i,j),0.)/cbottot)*frac(n)%tau_e
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
					Mr_avg=frac(nfrac_silt(1))%M/frac(n)%rho
					tau_e_avg=frac(nfrac_silt(1))%tau_e
				ENDIF
				
				ust=0.1*absU
				do tel=1,10 ! 10 iter is more than enough
					z0=kn_sed_avg/30.+0.11*nu_mol/MAX(ust,1.e-9) 
					! for tau shear on sediment don't use kn (which is result of bed ripples), use frac(n)%kn_sed; it is adviced to use kn_sed=dfloc
					ust=absU/MAX(1./kappa*log(distance_to_bed/z0),2.) !ust maximal 0.5*absU
				enddo
				kplus = MIN(kbed(i,j)+1,k1)
				kplus2 = MIN(kbed(i,j)+2,k1)
				tau=rho_b*ust*ust  
				DO n1=1,nfr_silt
					n=nfrac_silt(n1)
					erosion_avg(n) = Mr_avg*MAX(0.,(tau/tau_e_avg-1.))*ddt*bednotfixed(i,j,kbed(i,j))*morfac*bs_geo ! m3/m2	 erosion_avg is filled for silt fractions only with silt erosion
					IF ((interaction_bed.ge.6.and.kbed(i,j).eq.0).or.interaction_bed.eq.7) THEN !unlimited erosion in case kbed.eq.0
						erosionf(n) = erosion_avg(n) * (c_bed(n)/cfixedbed) !erosion per fraction
!						erosionf(n) = erosionf(n) 
!     &					-MIN(0.,ddt*0.5*(Diffcof(i,j,kplus)+Diffcof(i,j,kplus2))*(ccnew(n,i,j,kplus2)-ccnew(n,i,j,kplus))/dz)  !add upward diffusion to erosion flux
					ELSEIF (cbottot>0.) THEN
						erosionf(n) = erosion_avg(n) * (cbotnew(n,i,j)/cbottot) !erosion per fraction
!						erosionf(n) = erosionf(n) 
!     &					-MIN(0.,ddt*0.5*(Diffcof(i,j,kplus)+Diffcof(i,j,kplus2))*(ccnew(n,i,j,kplus2)-ccnew(n,i,j,kplus))/dz)  !add upward diffusion to erosion flux
						erosionf(n) = MIN(erosionf(n),(cbotnew(n,i,j)+SUM(Clivebed(n,i,j,0:kbed(i,j))))*dz/morfac2) ! m3/m2, not more material can be eroded as available 
						erosionf(n) = MAX(erosionf(n),0.)
					ELSEIF (cbedtot>0.) THEN
						erosionf(n) = erosion_avg(n) * (Clivebed(n,i,j,kbed(i,j))/cbedtot) !erosion per fraction
!						erosionf(n) = erosionf(n) 
!     &					-MIN(0.,ddt*0.5*(Diffcof(i,j,kplus)+Diffcof(i,j,kplus2))*(ccnew(n,i,j,kplus2)-ccnew(n,i,j,kplus))/dz)  !add upward diffusion to erosion flux
						erosionf(n) = MIN(erosionf(n),(cbotnew(n,i,j)+SUM(Clivebed(n,i,j,0:kbed(i,j))))*dz/morfac2) ! m3/m2, not more material can be eroded as available 
						erosionf(n) = MAX(erosionf(n),0.)
					ELSE
						erosionf(n) = 0.
					ENDIF
					
					kplus = MIN(kbed(i,j)+1,k1)
					ust=0.1*absU !re-calculate tau with kn_sed for deposition as it is not dependent on avg dpart in mixture
					do tel=1,10 ! 10 iter is more than enough
						z0=frac(n)%kn_sed/30.+0.11*nu_mol/MAX(ust,1.e-9) 
						! for tau shear on sediment don't use kn (which is result of bed ripples), use frac(n)%kn_sed; it is adviced to use kn_sed=dfloc
						ust=absU/MAX(1./kappa*log(distance_to_bed/z0),2.) !ust maximal 0.5*absU
					enddo
					tau=rho_b*ust*ust  !for deposition apply tau belonging to own frac(n)%kn_sed
					IF (depo_implicit.eq.1) THEN  !determine deposition as sink implicit
					 ccnew(n,i,j,kplus)=(ccnew(n,i,j,kplus)+erosionf(n)/dz)/ ! vol conc. [-]
     &				(1.-MAX(0.,(1.-tau/frac(n)%tau_d))*MIN(0.,Wsed(n,i,j,kbed(i,j)))*ddt/dz*bednotfixed_depo(i,j,kbed(i,j))*morfac)
					 depositionf(n) = MAX(0.,(1.-tau/frac(n)%tau_d))*ccnew(n,i,j,kplus)*MIN(0.,Wsed(n,i,j,kbed(i,j)))*ddt !ccfd
     & *bednotfixed_depo(i,j,kbed(i,j))*morfac				 ! m --> dep is negative due to negative wsed					 
					 cbotnew(n,i,j)=cbotnew(n,i,j)-b_update*morfac2*(erosionf(n)+depositionf(n))/(dz) ! vol conc. [-]  !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
					 cbotnewtot=cbotnewtot+cbotnew(n,i,j)
					 cbotnewtot_pos=cbotnewtot_pos+MAX(cbotnew(n,i,j),0.)
					 ctot_firstcel=ccnew(n,i,j,kplus)+ctot_firstcel					
					ELSE
					 depositionf(n) = MAX(0.,(1.-tau/frac(n)%tau_d))*ccnew(n,i,j,kplus)*MIN(0.,Wsed(n,i,j,kbed(i,j)))*ddt !ccfd
     & *bednotfixed_depo(i,j,kbed(i,j))*morfac			 ! m --> dep is negative due to negative wsed
					 ccnew(n,i,j,kplus)=ccnew(n,i,j,kplus)+(erosionf(n)+depositionf(n))/(dz) ! vol conc. [-]
					 cbotnew(n,i,j)=cbotnew(n,i,j)-b_update*morfac2*(erosionf(n)+depositionf(n))/(dz) ! vol conc. [-]  !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
					 cbotnewtot=cbotnewtot+cbotnew(n,i,j)
					 cbotnewtot_pos=cbotnewtot_pos+MAX(cbotnew(n,i,j),0.)
					 ctot_firstcel=ccnew(n,i,j,kplus)+ctot_firstcel
					ENDIF
				ENDDO
			ENDIF
			!! 2 determine erosion/sedimentation of mixture of all sand fractions
			IF (nfr_sand>0) THEN
				cbottot_sand=0.
				cbedtot_sand=0.
				mbottot_sand=0.
				mbottot_sand2=0.
				mbedtot_sand=0.
				rho_botsand=0.
				rho_botsand2=0.
				rho_bedsand=0.
				ws_botsand=0.
				ws_botsand2=0.
				ws_bedsand=0.
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
					cbottot_sand=cbottot_sand+MAX(cbotnew(n,i,j),0.)
					cbedtot_sand=cbedtot_sand+Clivebed(n,i,j,kbed(i,j))
					mbottot_sand=mbottot_sand+MAX(cbotnew(n,i,j),0.)*frac(n)%rho
					mbottot_sand2=mbottot_sand2+c_bed(n)*frac(n)%rho
					mbedtot_sand=mbedtot_sand+Clivebed(n,i,j,kbed(i,j))*frac(n)%rho
					PSD_bot_sand_massfrac(n1)=mbottot_sand
					PSD_bot_sand_massfrac2(n1)=mbottot_sand2
					PSD_bed_sand_massfrac(n1)=mbedtot_sand
					rho_botsand=rho_botsand+MAX(cbotnew(n,i,j),0.)*frac(n)%rho
					rho_botsand2=rho_botsand2+c_bed(n)*frac(n)%rho
					rho_bedsand=rho_bedsand+Clivebed(n,i,j,kbed(i,j))*frac(n)%rho
					ws_botsand=ws_botsand+MAX(cbotnew(n,i,j),0.)*frac(n)%ws
					ws_botsand2=ws_botsand2+c_bed(n)*frac(n)%ws
					ws_bedsand=ws_bedsand+Clivebed(n,i,j,kbed(i,j))*frac(n)%ws					
				ENDDO
				IF ((interaction_bed.ge.6.and.kbed(i,j).eq.0).or.interaction_bed.eq.7) THEN !unlimited erosion in case kbed.eq.0
					PSD_sand(1:nfr_sand)=PSD_bot_sand_massfrac2/mbottot_sand2
					rho_sand=rho_botsand2/cfixedbed
					ws_sand=ws_botsand2/cfixedbed				
				ELSEIF (mbottot_sand>0.) THEN
					PSD_sand(1:nfr_sand)=PSD_bot_sand_massfrac/mbottot_sand
					rho_sand=rho_botsand/cbottot_sand
					ws_sand=ws_botsand/cbottot_sand
				ELSEIF (mbedtot_sand>0.) THEN
					PSD_sand(1:nfr_sand)=PSD_bed_sand_massfrac/mbedtot_sand
					rho_sand=rho_bedsand/cbedtot_sand
					ws_sand=ws_bedsand/cbedtot_sand
				ELSE
					PSD_sand=0. ! there is no sand, hence there will be no erosion either therefore choice in PSD and resulting d50 is not important
					rho_sand=frac(nfrac_sand(1))%rho
					ws_sand=frac(nfrac_sand(1))%ws
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
				
				kplus = MIN(kbed(i,j)+1,k1)
				kplus2 = MIN(kbed(i,j)+2,k1)
				delta = (rho_sand-rho_b)/rho_b !(rho_sand-rcfd(i,j,kplus))/rcfd(i,j,kplus) !switched to using rho_b instead of rcfd 13-2-2019
				delta = MAX(delta,0.1) ! rho_sand must be > 2*rho_fluid
				gvector=sqrt(gx**2+gy**2+gz**2)
				Dstar = d50 * ((delta*gvector)/nu_mol**2)**(0.333333333333333)
				Shields_cr = 0.3/(1.+1.2*Dstar)+0.055*(1.-exp(-0.02*Dstar))  !Soulsby and Whitehouse 1997 curve through original Shields for threshold of motion sediment particles, Soulsy book Eq. SC(77)	
				Shields_cr = Shields_cr*calibfac_Shields_cr
				kn_sed_avg=kn_d50_multiplier*d50 !  kn=2*d50 is mentioned in VanRijn1984 paper, the pickup function which is applied here, however elsewhere vanRijn mentions larger kn_sed like 6*d50...)
				ust=0.1*absU
				do tel=1,10 ! 10 iter is more than enough
					z0=kn_sed_avg/30.+0.11*nu_mol/MAX(ust,1.e-9) 
					ust=absU/MAX(1./kappa*log(distance_to_bed/z0),2.) !ust maximal 0.5*absU
				enddo
				
				
				IF (cbed_method.eq.2) THEN
					k_maxU = MAXLOC(Uhor(i,j,1:kmax),DIM=1)
					IF (k_maxU.le.kbed(i,j)) THEN 
						cbed = MIN(cfixedbed,SUM(ccfd(1:nfrac,i,j,kplus)))
					ELSE 
						cbed = 0.
						DO k=kplus,k_maxU
							cbed = cbed + SUM(ccfd(1:nfrac,i,j,k))
						ENDDO 					
						cbed = cbed /MAX(DBLE(k_maxU-kbed(i,j)),1.)
						cbed = MIN(cfixedbed,cbed)	
					ENDIF 
				ELSE 
					cbed = MIN(cfixedbed,SUM(ccfd(1:nfrac,i,j,kplus)))
				ENDIF 
				IF (pickup_formula.eq.'vanrijn1984') THEN
					ustc2 = Shields_cr * gvector*delta*d50
					TT = (ust*ust-ustc2)/(MAX(ustc2,1.e-12))
					TT = MAX(TT,0.) !TT must be positive
					phipp = calibfac_sand_pickup*0.00033*Dstar**0.3*TT**1.5   ! general pickup function					
				ELSEIF (pickup_formula.eq.'nielsen1992') THEN
					Re_p=ABS(ws_sand)*d50/nu_mol
					CD = MAX(0.4,24./Re_p*(1.+0.15*Re_p**0.687)) ! Okayasu et al 2010 Eq.4 combined with van Rhee p28 eq3.5 0.4 limit for Re>2000 	
					!ust = kappa * absU / (log(15.05*dz/d50)) ! in Okayasu 2010 hydraulic rough flow is assumed with kn_sed=d50, in TUDflow3d relations valid for both smooth and rough flow are used and kn_sed is user defined 
					Shields_eff = 0.75 * CD * ust**2/(delta*gvector*d50)
					TT = (Shields_eff - Shields_cr)/Shields_cr
					TT = MAX(TT,0.) !TT must be positive
					phipp = calibfac_sand_pickup*0.00033*Dstar**0.3*TT**1.5   ! general pickup function						
				ELSEIF (pickup_formula.eq.'okayasu2010') THEN
					Re_p=ABS(ws_sand)*d50/nu_mol
					CD = MAX(0.4,24./Re_p*(1.+0.15*Re_p**0.687)) ! Okayasu et al 2010 Eq.4 combined with van Rhee p28 eq3.5 0.4 limit for Re>2000 	
					!ust = kappa * absU / (log(15.05*dz/d50)) ! in Okayasu 2010 hydraulic rough flow is assumed with kn_sed=d50, in TUDflow3d relations valid for both smooth and rough flow are used and kn_sed is user defined 
					Fd = 0.125*CD*rcfd(i,j,kplus)*pi*d50**2*ust**2
					dubdt = (absU-ubot(i,j))/ddt
					ubot(i,j)=absU
					Fi = 0.25*rcfd(i,j,kplus)*pi*d50**3*dubdt ! Ci=1.5 --> 1.5/6=0.25
					W = 0.166667*gvector*(rho_sand-rcfd(i,j,kplus))*pi*d50**3
					Fl = 0.025*rcfd(i,j,kplus)*pi*d50**2*ust**2 ! Cl=0.2 --> 0.2*0.5/4=0.025
					Shields_eff = (Fd+Fi)/MAX(W-Fl,1.e-12)
					TT = (Shields_eff - Shields_cr)/Shields_cr	
					TT = MAX(TT,0.) !TT must be positive
					phipp = calibfac_sand_pickup*0.00033*Dstar**0.3*TT**1.5   ! general pickup function	
				ELSEIF (pickup_formula.eq.'vanrijn2019') THEN
					ustc2 = Shields_cr * gvector*delta*d50
					TT = (ust*ust-ustc2)/(MAX(ustc2,1.e-12))
					TT = MAX(TT,0.) !TT must be positive
					phipp = calibfac_sand_pickup*0.00033*Dstar**0.3*TT**1.5   ! general pickup function	
					Shields_eff = ust**2/(delta*gvector*d50)
					phipp = phipp/(MAX(Shields_eff,1.)) ! correction: reduced pickup for high speed erosion
				ELSEIF (pickup_formula.eq.'VR2019_Cbed') THEN
					ustc2 = Shields_cr * gvector*delta*d50
					TT = (ust*ust-ustc2)/(MAX(ustc2,1.e-12))
					TT = MAX(TT,0.) !TT must be positive
					phipp = calibfac_sand_pickup*0.00033*Dstar**0.3*TT**1.5   ! general pickup function	
					Shields_eff = ust**2/(delta*gvector*d50)
					phipp = phipp/(MAX(Shields_eff,1.))*(cfixedbed-cbed)/(cfixedbed) ! correction: reduced pickup for high speed erosion and reduction for cbed according to VanRhee and Talmon 2010
				ELSEIF (pickup_formula.eq.'VR1984_Cbed') THEN
					ustc2 = Shields_cr * gvector*delta*d50
					TT = (ust*ust-ustc2)/(MAX(ustc2,1.e-12))
					TT = MAX(TT,0.) !TT must be positive
					phipp = calibfac_sand_pickup*0.00033*Dstar**0.3*TT**1.5   ! general pickup function	
					phipp = phipp*(cfixedbed-cbed)/(cfixedbed) ! correction: reduced pickup for for cbed according to VanRhee and Talmon 2010					
				ELSE
					TT=0.
					phipp = 0.  ! general pickup function					
				ENDIF
				IF (pickup_correction.eq.'MastBergenvdBerg2003') THEN 
					vs = MAX(sqrt(gvector*delta*d50),1.e-12)
					!vwal = (-cfixedbed*delta*sin(phi-alpha)/sin(alpha))/(delta_nsed/permeability_kt) !vwal is user input, here Eq.6 from MastBergenvdBerg2003 is mentioned to know how to calculate it (in Eq. 6 the minus sign was forgotten)
					ve = 0.5*vwal+vs*sqrt((0.5*vwal/vs)**2+fcor*phipp*delta/delta_nsed*(permeability_kl/vs))
					phipp = ve*cfixedbed/vs 
					
					ve_check = phipp*(delta*gvector*d50)**0.5*morfac*bs_geo +  ! this is erosion (positive value)
     & 					SUM(ccnew(1:nfrac,i,j,kplus))*(MIN(0.,Wsed(n,i,j,kbed(i,j))))*morfac ! this is depo (neg value) in m/s 
					ve_check = ve_check/cfixedbed  ! correction needed for pore volume 
					phipp = phipp + MAX(vwal*bs_geo-ve_check,0.)/((delta*gvector*d50)**0.5*morfac*bs_geo)*cfixedbed
					IF (wbed_correction.eq.1) Wbed(i,j)=ve_check+MAX(vwal*bs_geo-ve_check,0.)
				ENDIF
				IF (pickup_fluctuations.eq.1) THEN
					!1 add white noise to pickup
					phipp = phipp*pickup_random(i,j)					
				ENDIF
					
				IF (reduction_sedimentation_shields>0) THEN ! PhD thesis vRhee p146 eq 7.74
					reduced_sed = 1.-(ust*ust/(delta*gvector*d50))/reduction_sedimentation_shields
					reduced_sed = MAX(reduced_sed,0.)
				ELSE
					reduced_sed  = 1.
				ENDIF

				ccfdtot_firstcel=0.
				wsedbed=0.
				DO n1=1,nfr_sand
					n=nfrac_sand(n1)
					ccfdtot_firstcel=ccfd(n,i,j,kplus)*DBLE(depo_cbed_option)+ccfdtot_firstcel	 !depo_cbed_option=0 means correction (1-ccfdtot_firstcel) not needed
					wsedbed=wsedbed-ccnew(n,i,j,kplus)*MIN(0.,reduced_sed*Wsed(n,i,j,kbed(i,j)))*DBLE(depo_cbed_option) !wsed upward
				ENDDO
				wsedbed=wsedbed/MAX(1.-cfixedbed-ccfdtot_firstcel,1.e-6)
				wsedbed=MIN(wsedbed,0.1*dz/ddt) !to keep concentration positive: never may wsed be more than 0.1 grid cell in one time step
				IF (bedload_formula.eq.'vanrijn2007'.and.time_n.ge.tstart_morf2) THEN 
					kn_sed_avg=kn_d50_multiplier_bl*d50 
					ust=0.1*absUbl
					do tel=1,10 ! 10 iter is more than enough
						z0=kn_sed_avg/30.+0.11*nu_mol/MAX(ust,1.e-9) 
						ust=absUbl/MAX(1./kappa*log(distance_to_bed/z0),2.) !ust maximal 0.5*absU
					enddo	
					ustc2 = Shields_cr * gvector*delta*d50
					MME = MAX(0.,(ust*ust-ustc2)/ustc2)
					qb = calibfac_sand_bedload*0.5*rho_sand*d50*Dstar**(-0.3)*ust*MME ! [kg/m/s]	
!		IF (ABS(SUM(cbotnew(1:nfrac,i,j)))>0.5*cfixedbed) write(*,'(a,i6.1,i6.1,i6.1,f11.6,f11.6,f11.6,f11.6,f11.6,f11.6,f11.6)'),
!     & 'rank,i,j,qb,absU,uuR,uuL,vvR,vvL,cbotnew'
!     & ,rank,i,j,qb,absUbl,uuR_relax(i,j),uuL_relax(i,j),vvR_relax(i,j),vvL_relax(i,j),SUM(cbotnew(1:nfrac,i,j))
		IF ((rank.eq.1.or.rank.eq.0).and.i>96.and.i<100.and.j>1.and.j<8) write
     & (*,'(a,i6.1,i6.1,i6.1,i6.1,f11.6,f11.6,f11.6,f11.6,f11.6,f11.6,f11.6)'),'istep,rank,i,j,qb,absU,uuR,uuL,vvR,vvL,cbotnew'
     & ,istep,rank,i,j,qb,absUbl,uuR_relax(i,j),uuL_relax(i,j),vvR_relax(i,j),vvL_relax(i,j),SUM(cbotnew(1:nfrac,i,j))	 
					qb = qb*bednotfixed(i,j,kbed(i,j))*morfac*morfac2 !kg/m/s				
				ELSEIF (bedload_formula.eq.'vanrijn2003'.and.time_n.ge.tstart_morf2) THEN !as taken from D3D manual
					kn_sed_avg=kn_d50_multiplier_bl*d50 
					ust=0.1*absUbl
					do tel=1,10 ! 10 iter is more than enough
						z0=kn_sed_avg/30.+0.11*nu_mol/MAX(ust,1.e-9) 
						ust=absUbl/MAX(1./kappa*log(distance_to_bed/z0),2.) !ust maximal 0.5*absU
					enddo	
					ucr = sqrt(Shields_cr * gvector*delta*d50)*log(distance_to_bed/z0)/kappa
					MMM = absUbl**2/(delta*gvector*d50)
					MME = (MAX(absUbl-ucr,0.))**2/(delta*gvector*d50)
					qb = calibfac_sand_bedload*0.006*rho_sand*ws_sand*d50*MMM**0.5*MME**0.7 ! [kg/m/s] 
					qb = qb*bednotfixed(i,j,kbed(i,j))*morfac*morfac2 !kg/m/s				
				ELSEIF (bedload_formula.eq.'MeyPeMu1947'.and.time_n.ge.tstart_morf2) THEN !as taken from D3D manual
					kn_sed_avg=kn_d50_multiplier_bl*d50 
					ust=0.1*absUbl
					do tel=1,10 ! 10 iter is more than enough
						z0=kn_sed_avg/30.+0.11*nu_mol/MAX(ust,1.e-9) 
						ust=absUbl/MAX(1./kappa*log(distance_to_bed/z0),2.) !ust maximal 0.5*absU
					enddo	
					MME = (MAX(Shields-Shields_cr,0.))**1.5
					qb = calibfac_sand_bedload*8.*rho_sand*d50*sqrt(gvector*delta*d50)*MME ! [kg/m/s] 
					qb = qb*bednotfixed(i,j,kbed(i,j))*morfac*morfac2 !kg/m/s
				ENDIF 
				
				DO n1=1,nfr_sand
					n=nfrac_sand(n1)			
					erosion_avg(n) = phipp * (delta*gvector*d50)**0.5*ddt*bednotfixed(i,j,kbed(i,j))*morfac*bs_geo  !*rho_sand/rho_sand ! erosion flux in kg/m2/(kg/m3)= m3/m2=m
					IF ((interaction_bed.ge.6.and.kbed(i,j).eq.0).or.interaction_bed.eq.7) THEN !unlimited erosion in case kbed.eq.0
						erosionf(n) = erosion_avg(n) * (c_bed(n)/cfixedbed) !erosion per fraction
						erosionf(n) = erosionf(n) + ccnew(n,i,j,kplus)*MAX(0.,reduced_sed*Wsed(n,i,j,kbed(i,j)))*ddt !when wsed upward (e.g. fine fractions) then add negative deposition to erosion
!						erosionf(n) = erosionf(n) 
!     &					-MIN(0.,ddt*0.5*(Diffcof(i,j,kplus)+Diffcof(i,j,kplus2))*(ccnew(n,i,j,kplus2)-ccnew(n,i,j,kplus))/dz)  !add upward diffusion to erosion flux
						qbf(n) = qb * (c_bed(n)/cfixedbed)
					ELSEIF (cbottot_sand>0.) THEN
						erosionf(n) = erosion_avg(n) * (cbotnew(n,i,j)/cbottot_sand) !erosion per fraction
						erosionf(n) = erosionf(n) + ccnew(n,i,j,kplus)*MAX(0.,reduced_sed*Wsed(n,i,j,kbed(i,j)))*ddt !when wsed upward (e.g. fine fractions) then add negative deposition to erosion
!						erosionf(n) = erosionf(n) 
!     &					-MIN(0.,ddt*0.5*(Diffcof(i,j,kplus)+Diffcof(i,j,kplus2))*(ccnew(n,i,j,kplus2)-ccnew(n,i,j,kplus))/dz)  !add upward diffusion to erosion flux
						erosionf(n) = MIN(erosionf(n),(cbotnew(n,i,j)+SUM(Clivebed(n,i,j,0:kbed(i,j))))*dz/morfac2) ! m3/m2, not more material can be eroded as available 
						erosionf(n) = MAX(erosionf(n),0.)
						qbf(n) = qb * (cbotnew(n,i,j)/cbottot_sand)
					ELSEIF (cbedtot_sand>0.0) THEN
						erosionf(n) = erosion_avg(n) * (Clivebed(n,i,j,kbed(i,j))/cbedtot_sand) !erosion per fraction
						erosionf(n) = erosionf(n) + ccnew(n,i,j,kplus)*MAX(0.,reduced_sed*Wsed(n,i,j,kbed(i,j)))*ddt !when wsed upward (e.g. fine fractions) then add negative deposition to erosion
!						erosionf(n) = erosionf(n) 
!     &					-MIN(0.,ddt*0.5*(Diffcof(i,j,kplus)+Diffcof(i,j,kplus2))*(ccnew(n,i,j,kplus2)-ccnew(n,i,j,kplus))/dz)  !add upward diffusion to erosion flux
						erosionf(n) = MIN(erosionf(n),(cbotnew(n,i,j)+SUM(Clivebed(n,i,j,0:kbed(i,j))))*dz/morfac2) ! m3/m2, not more material can be eroded as available 
						erosionf(n) = MAX(erosionf(n),0.)
						qbf(n) = qb * (Clivebed(n,i,j,kbed(i,j))/cbedtot_sand)
					ELSE
						erosionf(n) = 0.
						qbf(n) = qb * (c_bed(n)/cfixedbed) !don't make qb zero this is not robust 
					ENDIF
					qbU(n,i,j) = qbf(n)*uuRrel
					qbV(n,i,j) = qbf(n)*vvRrel	
					cctot=0.
					DO k=kplus,kbed(i,j)+k_layer_pickup
						cctot=cctot+MAX(ccfd(n,i,j,k)+1.e-12,0.)
					ENDDO
					DO k=kplus+1,kplus+k_layer_pickup-1 
						ccnew(n,i,j,k) = ccnew(n,i,j,k) + erosionf(n)*MAX(ccfd(n,i,j,k)+1.e-12,0.)/cctot/dz !pickup is spread over multiple k-layers 
					ENDDO 
					IF (depo_implicit.eq.1) THEN  !determine deposition as sink implicit
					ccnew(n,i,j,kplus)=(ccnew(n,i,j,kplus)+erosionf(n)*MAX(ccfd(n,i,j,kplus)+1.e-12,0.)/cctot/dz)/ ! vol conc. [-]
     &      		(1.-(MIN(0.,reduced_sed*Wsed(n,i,j,kbed(i,j)))-wsedbed)*ddt/dz*bednotfixed_depo(i,j,kbed(i,j))*morfac)					
					depositionf(n) = ccnew(n,i,j,kplus)*(MIN(0.,reduced_sed*Wsed(n,i,j,kbed(i,j)))-wsedbed)*ddt !ccfd
     &     *bednotfixed_depo(i,j,kbed(i,j))*morfac! m --> dep is negative due to negative wsed					
					cbotnew(n,i,j)=cbotnew(n,i,j)-b_update*morfac2*(erosionf(n)+depositionf(n))/(dz) ! vol conc. [-] !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
					cbotnewtot=cbotnewtot+cbotnew(n,i,j)
					cbotnewtot_pos=cbotnewtot_pos+MAX(cbotnew(n,i,j),0.)
					ctot_firstcel=ccnew(n,i,j,kplus)+ctot_firstcel
					ELSE
					depositionf(n) = ccnew(n,i,j,kplus)*(MIN(0.,reduced_sed*Wsed(n,i,j,kbed(i,j)))-wsedbed)*ddt !ccfd
     &     *bednotfixed_depo(i,j,kbed(i,j))*morfac ! m --> dep is negative due to negative wsed
					ccnew(n,i,j,kplus)=ccnew(n,i,j,kplus)+(erosionf(n)*MAX(ccfd(n,i,j,kplus)+1.e-12,0.)/cctot+depositionf(n))/(dz) ! vol conc. [-]
					cbotnew(n,i,j)=cbotnew(n,i,j)-b_update*morfac2*(erosionf(n)+depositionf(n))/(dz) ! vol conc. [-] !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
					cbotnewtot=cbotnewtot+cbotnew(n,i,j)
					cbotnewtot_pos=cbotnewtot_pos+MAX(cbotnew(n,i,j),0.)
					ctot_firstcel=ccnew(n,i,j,kplus)+ctot_firstcel
					ENDIF
					IF ((bedload_formula.ne.'nonenon0000').and.time_n.ge.tstart_morf2) THEN
						flux = uuRrel*qbf(n)*Ru(i)*dphi2(j)*ddt ![kg] change in mass over this edge
				!write(*,*),'i,j,kbed(i,j),flux,uuRrel,qbf(n),Ru(i)*dphi2(j)',i,j,kbed(i,j),flux,uuRrel,qbf(n),Ru(i)*dphi2(j)
						IF (flux>0.) THEN
							d_cbotnew(n,i,j) = d_cbotnew(n,i,j) - flux/(rho_sand*dr(i)*Rp(i)*dphi2(j)*dz)
							d_cbotnew(n,i+1,j) = d_cbotnew(n,i+1,j) + flux/(rho_sand*dr(i+1)*Rp(i+1)*dphi2(j)*dz)
						ENDIF 
						flux = uuLrel*qbf(n)*Ru(i-1)*dphi2(j)*ddt ![kg] change in mass over this edge
				!write(*,*),'i,j,kbed(i,j),flux,uuLrel,qbf(n),Ru(i-1)*dphi2(j)',i,j,kbed(i,j),flux,uuLrel,qbf(n),Ru(i-1)*dphi2(j)
						IF (flux<0.) THEN
							d_cbotnew(n,i,j) = d_cbotnew(n,i,j) + flux/(rho_sand*dr(i)*Rp(i)*dphi2(j)*dz)
							d_cbotnew(n,i-1,j) = d_cbotnew(n,i-1,j) - flux/(rho_sand*dr(i-1)*Rp(i-1)*dphi2(j)*dz)
						ENDIF 						
						flux = vvRrel*qbf(n)*dr(i)*ddt ![kg] change in mass over this edge
				!write(*,*),'i,j,kbed(i,j),flux,vvRrel,qbf(n),dr(i)',i,j,kbed(i,j),flux,vvRrel,qbf(n),dr(i)
						IF (flux>0.) THEN 
							d_cbotnew(n,i,j) = d_cbotnew(n,i,j) - flux/(rho_sand*dr(i)*Rp(i)*dphi2(j)*dz)
							d_cbotnew(n,i,j+1) = d_cbotnew(n,i,j+1) + flux/(rho_sand*dr(i)*Rp(i)*dphi2(j+1)*dz)
						ENDIF 
						flux = vvLrel*qbf(n)*dr(i)*ddt ![kg] change in mass over this edge
				!write(*,*),'i,j,kbed(i,j),flux,vvLrel,qbf(n),dr(i)',i,j,kbed(i,j),flux,vvLrel,qbf(n),dr(i)
						IF (flux<0.) THEN 
							d_cbotnew(n,i,j) = d_cbotnew(n,i,j) + flux/(rho_sand*dr(i)*Rp(i)*dphi2(j)*dz)
							d_cbotnew(n,i,j-1) = d_cbotnew(n,i,j-1) - flux/(rho_sand*dr(i)*Rp(i)*dphi2(j-1)*dz)
						ENDIF						
					ENDIF 
				ENDDO
			ENDIF	
	
!update bedlevel for combined erosion deposition silt plus sand fractions: 
!combination silt and sand gives that in eroding cases sand will erode faster than silt from cbotnew and once all sand is eroded 
!it can take a while before all silt has eroded and next bed-cell is reached with sand again 
!--> this is a rather strong sheltering effect of even a small amount of silt --> improvement is possible
! individual fraction in cbotnew can become <0 but cbotnewtot can still be >0 this may not enter Clivebed and is prevented below
! negative cbotnew still cannot enter Clivebed, but bedupdate (kbed+1) is allowed when total cbotnew(fracs>0) is enough even when individual fractions are negative in cbotnew
			IF ((interaction_bed.eq.4.or.interaction_bed.eq.6).and.time_n.ge.tstart_morf2) THEN
				IF (cbotnewtot.lt.0.and.(kbed(i,j)-1).ge.0.and.SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed-1.e-9) THEN 
				!add half cell sediment on cbot account for further erosion without lowering 1 dz yet (because otherwise flipflop between ero-1dz and depo+1dz)
				! this is at 1*dz			
					DO n=1,nfrac ! also cbotnew(n,i,j) is le 0:
						cbotnew(n,i,j)=cbotnew(n,i,j)+0.5*Clivebed(n,i,j,kbed(i,j)) 
						Clivebed(n,i,j,kbed(i,j))=0.5*Clivebed(n,i,j,kbed(i,j)) 
					ENDDO	
				ELSEIF (cbotnewtot.lt.0.and.(kbed(i,j)-1).ge.0) THEN !erosion of 1 layer dz:
				! this is at 0.5*dz
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
					kbedt(i,j)=kbed(i,j)	 
				ELSEIF ((cbotnewtot_pos).ge.0.5*cfixedbed.and.kbed(i,j)+1.le.kmax.and.
     &             (SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed-1.e-9.or.kbed(i,j).eq.0)) THEN
					kbed(i,j)=kbed(i,j)+1
					kbedt(i,j)=kbed(i,j)
					drdt(i,j,kbed(i,j))=rho_b
					rnew(i,j,kbed(i,j))=rho_b
					rold(i,j,kbed(i,j))=rho_b
					drdt(i,j,kbed(i,j)+1)=rho_b
					rnew(i,j,kbed(i,j)+1)=rho_b
					rold(i,j,kbed(i,j)+1)=rho_b					 
					c_adjustA = 0.5*cfixedbed/cbotnewtot_pos
					! fill Clivebed to 0.5 cfixedbed from cbotnew and place all ccnew from 1st fluid cell into cbotnew as buffer against flip-flop erosive state
					
					DO n=1,nfrac 
						Clivebed(n,i,j,kbed(i,j))=c_adjustA*MAX(cbotnew(n,i,j),0.)  
						cbotnew(n,i,j)=cbotnew(n,i,j)-c_adjustA*MAX(cbotnew(n,i,j),0.)+ccnew(n,i,j,kbed(i,j))
						IF (morfac2.gt.1.0000001) THEN
						 cctot=0.
						 DO k=kbed(i,j)+1,kmax 
							cctot=cctot+MAX(ccfd(n,i,j,k),0.)
						 ENDDO						
						 IF (cctot<1e-9) THEN
						  ccnew(n,i,j,kbed(i,j)+1)=ccnew(n,i,j,kbed(i,j)+1)+(morfac2-1.)/morfac2*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
						  drdt(i,j,kbed(i,j)+1) = drdt(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						  rnew(i,j,kbed(i,j)+1) = rnew(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						  rold(i,j,kbed(i,j)+1) = rold(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density												
						 ELSE
						  DO k=kbed(i,j)+1,kmax
						   drdt(i,j,k)=rho_b
						   rnew(i,j,k)=rho_b
						   rold(i,j,k)=rho_b
						  ENDDO
						  DO k=kbed(i,j)+1,kmax !redistribute morfac2 buried sediment over water column above 
						   ccnew(n,i,j,k)=ccnew(n,i,j,k)+(morfac2-1.)/morfac2*MAX(ccfd(n,i,j,k),0.)/cctot*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
						   drdt(i,j,k) = drdt(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						   rnew(i,j,k) = rnew(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						   rold(i,j,k) = rold(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density	
						  ENDDO 
						 ENDIF
						ENDIF
						ccnew(n,i,j,kbed(i,j))=0. 
						drdt(i,j,kbed(i,j)) = drdt(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						rnew(i,j,kbed(i,j)) = rnew(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						rold(i,j,kbed(i,j)) = rold(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density							
					ENDDO	
				ELSEIF (ctot_firstcel.ge.cfixedbed.and.kbed(i,j)+1.le.kmax.and.
     &             (SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed-1.e-9.or.kbed(i,j).eq.0)) THEN
					kbed(i,j)=kbed(i,j)+1
					kbedt(i,j)=kbed(i,j)
					drdt(i,j,kbed(i,j))=rho_b
					rnew(i,j,kbed(i,j))=rho_b
					rold(i,j,kbed(i,j))=rho_b
					drdt(i,j,kbed(i,j)+1)=rho_b
					rnew(i,j,kbed(i,j)+1)=rho_b
					rold(i,j,kbed(i,j)+1)=rho_b					 
					c_adjustA = (0.5*cfixedbed-cbotnewtot_pos)/ctot_firstcel
					! cbotnewtot_pos is smaller as 0.5*cfixedbed; otherwise captured in ELSEIF above 
					! fill Clivebed to 0.5 cfixedbed, first from cbotnew and top off to 0.5*cfixedbed from ccnew, move remainder of ccnew into cbotnew as buffer against flip-flop erosive state
					
					DO n=1,nfrac 
						Clivebed(n,i,j,kbed(i,j))=MAX(cbotnew(n,i,j),0.)+c_adjustA*ccnew(n,i,j,kbed(i,j)) 
						cbotnew(n,i,j)=(1.-c_adjustA)*ccnew(n,i,j,kbed(i,j))
						IF (morfac2.gt.1.0000001) THEN
						 cctot=0.
						 DO k=kbed(i,j)+1,kmax 
							cctot=cctot+MAX(ccfd(n,i,j,k),0.)
						 ENDDO						
						 IF (cctot<1e-9) THEN
						  ccnew(n,i,j,kbed(i,j)+1)=ccnew(n,i,j,kbed(i,j)+1)+(morfac2-1.)/morfac2*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
						  drdt(i,j,kbed(i,j)+1) = drdt(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						  rnew(i,j,kbed(i,j)+1) = rnew(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						  rold(i,j,kbed(i,j)+1) = rold(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density												
						 ELSE
						  DO k=kbed(i,j)+1,kmax
						   drdt(i,j,k)=rho_b
						   rnew(i,j,k)=rho_b
						   rold(i,j,k)=rho_b
						  ENDDO
						  DO k=kbed(i,j)+1,kmax !redistribute morfac2 buried sediment over water column above 
						   ccnew(n,i,j,k)=ccnew(n,i,j,k)+(morfac2-1.)/morfac2*MAX(ccfd(n,i,j,k),0.)/cctot*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
						   drdt(i,j,k) = drdt(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						   rnew(i,j,k) = rnew(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						   rold(i,j,k) = rold(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density	
						  ENDDO 
						 ENDIF
						ENDIF
						ccnew(n,i,j,kbed(i,j))=0. 
						drdt(i,j,kbed(i,j)) = drdt(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						rnew(i,j,kbed(i,j)) = rnew(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
						rold(i,j,kbed(i,j)) = rold(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density							
					ENDDO	
				ELSEIF (cbotnewtot_pos.ge.0.5*cfixedbed.and.kbed(i,j)+1.le.kmax) THEN
				! Clivebed not completely full; therefore no bedupdate kbed+1, but only redistribution from cbotnew tot Clivebed 
				! delay in this elseif is no longer needed as the potential flip-flop is with the first IF which is only a bookkeeping flip-flop; not changing kbed up and down unneeded			
					c_adjustA = (cfixedbed-SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cbotnewtot_pos
					! top off Clivebed to cfixedbed from cbotnew and leave remainder into cbotnew 
					DO n=1,nfrac 
						Clivebed(n,i,j,kbed(i,j))=Clivebed(n,i,j,kbed(i,j))+c_adjustA*MAX(cbotnew(n,i,j),0.)  
						cbotnew(n,i,j)=cbotnew(n,i,j)-c_adjustA*MAX(cbotnew(n,i,j),0.)						
					ENDDO																			
				ENDIF

			ENDIF
			
			
!!!			IF (interaction_bed.eq.4.or.interaction_bed.eq.6) THEN
!!!				IF (cbotnewtot.lt.0.and.(kbed(i,j)-1).ge.0.and.SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed) THEN 
!!!				!add half cell sediment on cbot account for further erosion without lowering 1 dz yet (because otherwise flipflop between ero-1dz and depo+1dz)
!!!					DO n=1,nfrac ! also cbotnew(n,i,j) is le 0:
!!!						cbotnew(n,i,j)=cbotnew(n,i,j)+0.5*Clivebed(n,i,j,kbed(i,j)) 
!!!						Clivebed(n,i,j,kbed(i,j))=0.5*Clivebed(n,i,j,kbed(i,j)) 
!!!					ENDDO	
!!!				ELSEIF (cbotnewtot.lt.0.and.(kbed(i,j)-1).ge.0) THEN !erosion of 1 layer dz:
!!!					kplus = MIN(kbed(i,j)+1,k1)
!!!					drdt(i,j,kbed(i,j))=rho_b
!!!					rnew(i,j,kbed(i,j))=rho_b
!!!					rold(i,j,kbed(i,j))=rho_b					
!!!					DO n=1,nfrac ! also cbotnew(n,i,j) is le 0:
!!!						cbotnew(n,i,j)=cbotnew(n,i,j)+Clivebed(n,i,j,kbed(i,j)) !top layer is now previous top cel bed minus (erosion>top-layer)
!!!						ccnew(n,i,j,kbed(i,j))=(erosionf(n)+depositionf(n))/(dz) !assign all erosion+depo as new sediment concentration new bottom layer fluid
!!!						ccnew(n,i,j,kplus)=ccnew(n,i,j,kplus)-(erosionf(n)+depositionf(n))/(dz) !remove erosion+depo from previous bottom layer fluid
!!!						Clivebed(n,i,j,kbed(i,j))=0. ! not bed anymore but fluid
!!!						drdt(i,j,kbed(i,j)) = drdt(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						rnew(i,j,kbed(i,j)) = rnew(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						rold(i,j,kbed(i,j)) = rold(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!					ENDDO
!!!					!kbed(i,j)=MAX(kbed(i,j)-1,0)  !update bed level at end		
!!!					kbed(i,j)=kbed(i,j)-1
!!!					kbedt(i,j)=kbed(i,j)
!!!				!! 31-8-2018 switched top 2 lines ELSEIF on instead of bottom 2 lines because ctot_firstcel can become >>cbed in TSHD placement sim; however previous sims diffuser depostion were done with bottom 2 lines!!
!!!!!				ELSEIF ((cbotnewtot_pos+ctot_firstcel).ge.cfixedbed.and.kbed(i,j)+1.le.kmax.and.cbotnewtot_pos.gt.1.e-12.and.
!!!!!     &             (SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed-1.e-12.or.kbed(i,j).eq.0).and.depo_cbed_option.eq.0) THEN
!!!!! test 8-may-2019 with this elseif instead of 2 lines above (=better!; no large zone with near-zero concentration first cell fluid):	 
!!!				ELSEIF ((cbotnewtot_pos).ge.cfixedbed.and.kbed(i,j)+1.le.kmax.and.
!!!     &             (SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed-1.e-12.or.kbed(i,j).eq.0).and.depo_cbed_option.eq.0) THEN
!!!	 
!!!!     &               .and.MINVAL(cbotnew(1:nfrac,i,j)).ge.0) THEN !if kbed=0 then sedimentation can happen even if Clivebed empty
!!!!				ELSEIF (cbotnewtot.ge.cfixedbed.and.kbed(i,j)+1.le.kmax.and.cbotnewtot.gt.1.e-12.and.
!!!!     &             (SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed-1.e-12.or.kbed(i,j).eq.0)) THEN !if kbed=0 then sedimentation can happen even if Clivebed empty	 
!!!					kbed(i,j)=kbed(i,j)+1
!!!					kbedt(i,j)=kbed(i,j)
!!!					drdt(i,j,kbed(i,j))=rho_b
!!!					rnew(i,j,kbed(i,j))=rho_b
!!!					rold(i,j,kbed(i,j))=rho_b
!!!					drdt(i,j,kbed(i,j)+1)=rho_b
!!!					rnew(i,j,kbed(i,j)+1)=rho_b
!!!					rold(i,j,kbed(i,j)+1)=rho_b					
!!!					!kbed(i,j)=MIN(kbed(i,j)+1,kmax) !update bed level at start sedimentation 
!!!					c_adjustA = MAX(cfixedbed-ctot_firstcel,0.)/MAX(cbotnewtot_pos,1.e-12)    !first fluid cel not yet filled up --> c_adjustA>0 & c_adjustB=0--> fluid cel transformed into bed and sediment from cbotnew moved to bed
!!!					c_adjustB = MIN(cfixedbed-ctot_firstcel,0.)/MAX(ctot_firstcel,1.e-12) !first fluid cel already more than filled up --> c_adjustB<0&c_adjustA=0 --> fluid cel transformed into bed and excess sediment to cbotnew
!!!!					DO n=1,nfrac 
!!!!						Clivebed(n,i,j,kbed(i,j))=ccnew(n,i,j,kbed(i,j))+
!!!!     &						c_adjustA*cbotnew(n,i,j)+c_adjustB*ccnew(n,i,j,kbed(i,j))  ! apply sedimentation ratio between fractions new sediment concentration of cells within bed
!!!!						cbotnew(n,i,j)=cbotnew(n,i,j)-c_adjustA*cbotnew(n,i,j)-c_adjustB*ccnew(n,i,j,kbed(i,j))
!!!!!     &						+(morfac2-1.)*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
!!!!						ccnew(n,i,j,kbed(i,j)+1)=ccnew(n,i,j,kbed(i,j)+1)+(morfac2-1.)/morfac2*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
!!!!						drdt(i,j,kbed(i,j)+1) = drdt(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!!						rnew(i,j,kbed(i,j)+1) = rnew(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!!						rold(i,j,kbed(i,j)+1) = rold(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density						
!!!!						ccnew(n,i,j,kbed(i,j))=0. 
!!!!						drdt(i,j,kbed(i,j)) = drdt(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!!						rnew(i,j,kbed(i,j)) = rnew(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!!						rold(i,j,kbed(i,j)) = rold(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density							
!!!!					ENDDO
!!!!					cctot=0.
!!!!					DO k=kbed(i,j)+1,kmax 
!!!!						cctot=cctot+SUM(ccfd(1:nfrac,i,j,k))
!!!!					ENDDO
!!!					DO n=1,nfrac 
!!!						Clivebed(n,i,j,kbed(i,j))=ccnew(n,i,j,kbed(i,j))+
!!!     &						c_adjustA*MAX(cbotnew(n,i,j),0.)+c_adjustB*ccnew(n,i,j,kbed(i,j))  ! apply sedimentation ratio between fractions new sediment concentration of cells within bed
!!!						cbotnew(n,i,j)=cbotnew(n,i,j)-c_adjustA*MAX(cbotnew(n,i,j),0.)-c_adjustB*ccnew(n,i,j,kbed(i,j))
!!!!     &						+(morfac2-1.)*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
!!!						IF (morfac2.gt.1.0000001) THEN
!!!						 cctot=0.
!!!						 DO k=kbed(i,j)+1,kmax 
!!!							cctot=cctot+MAX(ccfd(n,i,j,k),0.)
!!!						 ENDDO						
!!!						 IF (cctot<1e-12) THEN
!!!						  ccnew(n,i,j,kbed(i,j)+1)=ccnew(n,i,j,kbed(i,j)+1)+(morfac2-1.)/morfac2*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
!!!						  drdt(i,j,kbed(i,j)+1) = drdt(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						  rnew(i,j,kbed(i,j)+1) = rnew(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						  rold(i,j,kbed(i,j)+1) = rold(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density												
!!!						 ELSE
!!!						  DO k=kbed(i,j)+1,kmax
!!!						   drdt(i,j,k)=rho_b
!!!						   rnew(i,j,k)=rho_b
!!!						   rold(i,j,k)=rho_b
!!!						  ENDDO
!!!						  DO k=kbed(i,j)+1,kmax !redistribute morfac2 buried sediment over water column above 
!!!						   ccnew(n,i,j,k)=ccnew(n,i,j,k)+(morfac2-1.)/morfac2*MAX(ccfd(n,i,j,k),0.)/cctot*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
!!!						   drdt(i,j,k) = drdt(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						   rnew(i,j,k) = rnew(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						   rold(i,j,k) = rold(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density	
!!!						  ENDDO 
!!!						 ENDIF
!!!						ENDIF
!!!						ccnew(n,i,j,kbed(i,j))=0. 
!!!						drdt(i,j,kbed(i,j)) = drdt(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						rnew(i,j,kbed(i,j)) = rnew(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						rold(i,j,kbed(i,j)) = rold(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density							
!!!					ENDDO	
!!!				ELSEIF ((ctot_firstcel).ge.cfixedbed.and.kbed(i,j)+1.le.kmax.and.
!!!     &             (SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed-1.e-12.or.kbed(i,j).eq.0).and.depo_cbed_option.eq.0) THEN
!!!	 
!!!					kbed(i,j)=kbed(i,j)+1
!!!					kbedt(i,j)=kbed(i,j)
!!!					drdt(i,j,kbed(i,j))=rho_b
!!!					rnew(i,j,kbed(i,j))=rho_b
!!!					rold(i,j,kbed(i,j))=rho_b
!!!					drdt(i,j,kbed(i,j)+1)=rho_b
!!!					rnew(i,j,kbed(i,j)+1)=rho_b
!!!					rold(i,j,kbed(i,j)+1)=rho_b					
!!!					!kbed(i,j)=MIN(kbed(i,j)+1,kmax) !update bed level at start sedimentation 
!!!					c_adjustA = MAX(cfixedbed-ctot_firstcel,0.)/MAX(cbotnewtot_pos,1.e-12)    !first fluid cel not yet filled up --> c_adjustA>0 & c_adjustB=0--> fluid cel transformed into bed and sediment from cbotnew moved to bed
!!!					c_adjustB = MIN(cfixedbed-ctot_firstcel,0.)/MAX(ctot_firstcel,1.e-12) !first fluid cel already more than filled up --> c_adjustB<0&c_adjustA=0 --> fluid cel transformed into bed and excess sediment to cbotnew
!!!					DO n=1,nfrac 
!!!						Clivebed(n,i,j,kbed(i,j))=ccnew(n,i,j,kbed(i,j))+
!!!     &						c_adjustA*MAX(cbotnew(n,i,j),0.)+c_adjustB*ccnew(n,i,j,kbed(i,j))  ! apply sedimentation ratio between fractions new sediment concentration of cells within bed
!!!						cbotnew(n,i,j)=cbotnew(n,i,j)-c_adjustA*MAX(cbotnew(n,i,j),0.)-c_adjustB*ccnew(n,i,j,kbed(i,j))
!!!!     &						+(morfac2-1.)*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
!!!
!!!						IF (morfac2.gt.1.0000001) THEN
!!!						 cctot=0.
!!!						 DO k=kbed(i,j)+1,kmax 
!!!							cctot=cctot+MAX(ccfd(n,i,j,k),0.)
!!!						 ENDDO						
!!!						 IF (cctot<1e-12) THEN
!!!						  ccnew(n,i,j,kbed(i,j)+1)=ccnew(n,i,j,kbed(i,j)+1)+(morfac2-1.)/morfac2*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
!!!						  drdt(i,j,kbed(i,j)+1) = drdt(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						  rnew(i,j,kbed(i,j)+1) = rnew(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						  rold(i,j,kbed(i,j)+1) = rold(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density												
!!!						 ELSE
!!!						  DO k=kbed(i,j)+1,kmax
!!!						   drdt(i,j,k)=rho_b
!!!						   rnew(i,j,k)=rho_b
!!!						   rold(i,j,k)=rho_b
!!!						  ENDDO
!!!						  DO k=kbed(i,j)+1,kmax !redistribute morfac2 buried sediment over water column above 
!!!						   ccnew(n,i,j,k)=ccnew(n,i,j,k)+(morfac2-1.)/morfac2*MAX(ccfd(n,i,j,k),0.)/cctot*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
!!!						   drdt(i,j,k) = drdt(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						   rnew(i,j,k) = rnew(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						   rold(i,j,k) = rold(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density	
!!!						  ENDDO 
!!!						 ENDIF
!!!						ENDIF
!!!						ccnew(n,i,j,kbed(i,j))=0. 
!!!						drdt(i,j,kbed(i,j)) = drdt(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						rnew(i,j,kbed(i,j)) = rnew(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						rold(i,j,kbed(i,j)) = rold(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density							
!!!					ENDDO					
!!!				!! 24-2-2020: added option to increase kbed already when 0.5*dz is filled instead of 1*dz (erosive and depositing state in this code are more similar then+advantage 2nd order ibm)
!!!!				ELSEIF ((cbotnewtot_pos).ge.0.5*cfixedbed.and.kbed(i,j)+1.le.kmax.and.
!!!!     &             (SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed-1.e-12.or.kbed(i,j).eq.0).and.depo_cbed_option.eq.0) THEN
!!!
!!!	 
!!!				!! 27-2-2019 test update bed without burying ctot_firstcel because that makes total sediment in fluid smaller abruptly
!!!				ELSEIF ((cbotnewtot_pos).ge.cfixedbed.and.kbed(i,j)+1.le.kmax.and.
!!!     &             (SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed-1.e-12.or.kbed(i,j).eq.0).and.depo_cbed_option.eq.1) THEN
!!!
!!!					kbed(i,j)=kbed(i,j)+1
!!!					kbedt(i,j)=kbed(i,j)
!!!					drdt(i,j,kbed(i,j))=rho_b
!!!					rnew(i,j,kbed(i,j))=rho_b
!!!					rold(i,j,kbed(i,j))=rho_b
!!!					drdt(i,j,kbed(i,j)+1)=rho_b
!!!					rnew(i,j,kbed(i,j)+1)=rho_b
!!!					rold(i,j,kbed(i,j)+1)=rho_b					
!!!					!kbed(i,j)=MIN(kbed(i,j)+1,kmax) !update bed level at start sedimentation 
!!!					!!!c_adjustA = MAX(cfixedbed-ctot_firstcel,0.)/MAX(cbotnewtot_pos,1.e-12)    !first fluid cel not yet filled up --> c_adjustA>0 & c_adjustB=0--> fluid cel transformed into bed and sediment from cbotnew moved to bed
!!!					!!!c_adjustB = MIN(cfixedbed-ctot_firstcel,0.)/MAX(ctot_firstcel,1.e-12) !first fluid cel already more than filled up --> c_adjustB<0&c_adjustA=0 --> fluid cel transformed into bed and excess sediment to cbotnew
!!!
!!!					DO n=1,nfrac 
!!!						Clivebed(n,i,j,kbed(i,j))=cfixedbed/cbotnewtot_pos*MAX(cbotnew(n,i,j),0.) 
!!!						cbotnew(n,i,j)=cbotnew(n,i,j)-cfixedbed/cbotnewtot_pos*MAX(cbotnew(n,i,j),0.) 
!!!						cctot=0.
!!!						DO k=kbed(i,j)+1,kmax 
!!!							cctot=cctot+MAX(ccfd(n,i,j,k),0.)
!!!						ENDDO						
!!!						IF (cctot<1e-12) THEN !no sediment-Rouse profile available, therefore place ccnew of cell converted into bed in first cell above:
!!!						 ccnew(n,i,j,kbed(i,j)+1)=ccnew(n,i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
!!!						 drdt(i,j,kbed(i,j)+1) = drdt(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						 rnew(i,j,kbed(i,j)+1) = rnew(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						 rold(i,j,kbed(i,j)+1) = rold(i,j,kbed(i,j)+1)+ccnew(n,i,j,kbed(i,j)+1)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density												
!!!						ELSE
!!!						 DO k=kbed(i,j)+1,kmax
!!!						  drdt(i,j,k)=rho_b
!!!						  rnew(i,j,k)=rho_b
!!!						  rold(i,j,k)=rho_b
!!!						 ENDDO
!!!						 DO k=kbed(i,j)+1,kmax !redistribute ccnew of cell converted into bed over full water column:
!!!						  ccnew(n,i,j,k)=ccnew(n,i,j,k)+MAX(ccfd(n,i,j,k),0.)/cctot*ccnew(n,i,j,kbed(i,j)) !morfac2 makes bed changes faster but leaves c-fluid same: every m3 sediment in fluid corresponds to morfac2 m3 in bed! 
!!!						  drdt(i,j,k) = drdt(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						  rnew(i,j,k) = rnew(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						  rold(i,j,k) = rold(i,j,k)+ccnew(n,i,j,k)*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density	
!!!						 ENDDO 
!!!						ENDIF
!!!						ccnew(n,i,j,kbed(i,j))=0. 
!!!						drdt(i,j,kbed(i,j)) = drdt(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						rnew(i,j,kbed(i,j)) = rnew(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
!!!						rold(i,j,kbed(i,j)) = rold(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density							
!!!					ENDDO										
!!!				ELSEIF (cbotnewtot_pos.gt.1.e-12.and.SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).lt.cfixedbed.and.kbed(i,j).gt.0) THEN
!!!!     &                 .and.MINVAL(cbotnew(1:nfrac,i,j)).ge.0) THEN !only allowed if kbed>0, because Clivebed(:,:,0)=0
!!!					! add sediment to Clivebed to bring it to cfixedbed again 
!!!					c_adjust=MIN((cfixedbed-SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cbotnewtot_pos,1.)				
!!!					DO n=1,nfrac ! also cbotnew(n,i,j) is le 0:
!!!						Clivebed(n,i,j,kbed(i,j))=Clivebed(n,i,j,kbed(i,j))+c_adjust*MAX(cbotnew(n,i,j),0.)
!!!						cbotnew(n,i,j)=cbotnew(n,i,j)-c_adjust*MAX(cbotnew(n,i,j),0.)
!!!					ENDDO				
!!!				ENDIF
!!!			ENDIF			
		  ENDDO
		ENDDO
		! kbed is updated, hence also zbed needs to be updated:
		DO i=0,i1
			DO j=0,j1
				zbed(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(cbotnew(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
			ENDDO
		ENDDO		
		IF (wbed_correction.eq.1) call bound_cbot(Wbed)
	ENDIF		
	
	IF ((interaction_bed.eq.4.or.interaction_bed.eq.6).and.time_n.ge.tstart_morf.and.time_n.ge.tstart_morf2) THEN
		IF (bedload_formula.ne.'nonenon0000') THEN 
!! mpi transfer sum d_cbotnew over edges:
			call shiftf_lreverse(d_cbotnew,cbf) 
			call shiftb_lreverse(d_cbotnew,cbb) 

			if (periodicy.eq.0.or.periodicy.eq.2) then
				if (rank.eq.0) then ! boundaries in j-direction
				  do n=1,nfrac
				   do i=1,imax
					   d_cbotnew(n,i,1) = 0. !bedload cannot fill or empty first grid cell: d.dn flux = 0.
					   d_cbotnew(n,i,jmax) = d_cbotnew(n,i,jmax) + cbb(n,i) 
				   enddo
				  enddo
				elseif (rank.eq.px-1) then
				  do n=1,nfrac
				   do i=1,imax
					   d_cbotnew(n,i,1) = d_cbotnew(n,i,1) + cbf(n,i)
					   d_cbotnew(n,i,jmax) = 0. !bedload cannot fill or empty first grid cell: d.dn flux = 0.
				   enddo
				  enddo		   
				else
				  do n=1,nfrac
				   do i=1,imax
					   d_cbotnew(n,i,1) = d_cbotnew(n,i,1) + cbf(n,i)
					   d_cbotnew(n,i,jmax) = d_cbotnew(n,i,jmax) + cbb(n,i) 
				   enddo
				  enddo
				endif
			else
			  do n=1,nfrac
			   do i=1,imax
				   d_cbotnew(n,i,1) = d_cbotnew(n,i,1) + cbf(n,i)
				   d_cbotnew(n,i,jmax) = d_cbotnew(n,i,jmax) + cbb(n,i) 
			   enddo
			  enddo
			endif
			  do n=1,nfrac
			   do j=1,jmax
				   d_cbotnew(n,1,j) = 0. !bedload cannot fill or empty first grid cell: d.dn flux = 0.
				   d_cbotnew(n,imax,j) = 0.
			   enddo
			  enddo
	
!! add c_cbotnew with original cbotnew:
			DO i=1,imax
				DO j=1,jmax
					DO n=1,nfrac
						cbotnew(n,i,j)=cbotnew(n,i,j)+d_cbotnew(n,i,j)
					ENDDO
				ENDDO
			ENDDO	
		ENDIF 
		IF (U_TSHD>0.) THEN !reset "inflow" bedlevel
			kbed(1,0:j1)=kbedin(0:j1)
			kbedt(1,0:j1)=kbedin(0:j1)
			kbed(0,0:j1)=kbedin(0:j1)
			kbedt(0,0:j1)=kbedin(0:j1)			
			do j=0,j1 
				do i=0,1
					if (kbed(i,j).eq.0) then !if "inflow" bed is zero than force Clivebed and cbotnew to zero
							Clivebed(n,i,j,0)=0.
							cbotnew(n,i,j)=0. !make buffer layer empty 
					endif
					do k=1,kbed(i,j) ! assign initial bed concentrations; k=0 remains empty
						do n=1,nfrac
							Clivebed(n,i,j,k)=c_bed(n)
							cbotnew(n,i,j)=0. !make buffer layer empty 
						enddo
					enddo
					do k=kbed(i,j)+1,kmax ! remove all bed-sediment above new bed
						do n=1,nfrac
							Clivebed(n,i,j,k)=0.
						enddo
					enddo					
				enddo
			enddo		  			
		ENDIF	
		call bound_cbot_integer(kbed) ! apply correct boundary conditions for updated kbed	
		DO i=0,i1
			DO j=0,j1
				zbed(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(cbotnew(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
			ENDDO
		ENDDO
		
		nav=0
		have_avalanched=1  ! default do avalanche, during avalanche procedure this switch can be turned into 0 to stop avalanching
		IF (avalanche_until_done.eq.1) THEN
			n_av=MAX(NINT(morfac2),NINT(morfac),kmax*1000)
		ELSE 
			n_av=MAX(NINT(morfac2),NINT(morfac))
		ENDIF		
		IF (MAXVAL(av_slope).gt.0) THEN
		 DO tel=1,n_av ! normally avalanche 1 time every timestep; with morfac more times avalanche every timestep
		  IF (have_avalanched>0.1) then !.true.) then !have_avalanched>0.1) THEN
			nav=nav+1
			d_cbotnew = 0.
			DO i=1,imax
				DO j=1,jmax
					zb_all(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+ (SUM(cbotnew(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
				ENDDO
			ENDDO
		! mpi boundary for zb_all:
			call shiftf_l(zb_all,zbf) 
			call shiftb_l(zb_all,zbb) 

			if (periodicy.eq.0.or.periodicy.eq.2) then
				if (rank.eq.0) then ! boundaries in j-direction
				   do i=1,imax
				   zb_all(i,0) = zb_all(i,1) 
				   zb_all(i,j1) =zbb(i) 
				   enddo
				elseif (rank.eq.px-1) then
				   do i=1,imax
				   zb_all(i,0) = zbf(i)
				   zb_all(i,j1) =zb_all(i,jmax)  
				   enddo
				else
				   do i=1,imax
				   zb_all(i,0) = zbf(i)
				   zb_all(i,j1) =zbb(i) 
				   enddo
				endif
			else
			   do i=1,imax
				   zb_all(i,0) = zbf(i)
				   zb_all(i,j1) =zbb(i) 
			   enddo
			endif
			 ! boundaries in i-direction
			if (periodicx.eq.0.or.periodicx.eq.2) then
				 do j=0,j1
				   zb_all(0,j)    =    zb_all(1,j)
				   zb_all(i1,j)   =    zb_all(imax,j)
				 enddo   
			else 
				 do j=0,j1
				   zb_all(0,j)    =    zb_all(imax,j)
				   zb_all(i1,j)   =    zb_all(1,j)
				 enddo   
			endif			
			
			IF (avalanche_until_done.eq.1) THEN
				have_avalanched=0.
			ELSE 
				have_avalanched=1.
			ENDIF
			DO i=1,imax
				DO j=1,jmax
					sl1=(Rp(i)-Rp(i-1))/MAX(zb_all(i,j)-zb_all(i-1,j),1.e-18)
					sl2=(Rp(i+1)-Rp(i))/MAX(zb_all(i,j)-zb_all(i+1,j),1.e-18)
					sl3=(Rp(i)*phip(j)-Rp(i)*phip(j-1))/MAX(zb_all(i,j)-zb_all(i,j-1),1.e-18)
					sl4=(Rp(i)*phip(j+1)-Rp(i)*phip(j))/MAX(zb_all(i,j)-zb_all(i,j+1),1.e-18)					
					sl5=SQRT((Rp(i)*phip(j+1)-Rp(i)*phip(j))**2+(Rp(i)-Rp(i-1))**2)/MAX(zb_all(i,j)-zb_all(i-1,j+1),1.e-18)
					sl6=SQRT((Rp(i)*phip(j+1)-Rp(i)*phip(j))**2+(Rp(i)-Rp(i+1))**2)/MAX(zb_all(i,j)-zb_all(i+1,j+1),1.e-18)
					sl7=SQRT((Rp(i)*phip(j-1)-Rp(i)*phip(j))**2+(Rp(i)-Rp(i+1))**2)/MAX(zb_all(i,j)-zb_all(i+1,j-1),1.e-18)
					sl8=SQRT((Rp(i)*phip(j-1)-Rp(i)*phip(j))**2+(Rp(i)-Rp(i-1))**2)/MAX(zb_all(i,j)-zb_all(i-1,j-1),1.e-18)
					maxbedslope(i,j)=MIN(sl1,sl2,sl3,sl4,sl5,sl6,sl7,sl8)
					IF (maxbedslope(i,j).lt.0.9999*av_slope(i,j,kbed(i,j))*bednotfixed(i,j,kbed(i,j)).and.kbed(i,j).ge.1) THEN
					! avalanche...
						have_avalanched=have_avalanched+1.
						IF (sl1.le.maxbedslope(i,j)) THEN
							itrgt=i-1
							jtrgt=j
							dbed = zb_all(i,j)-zb_all(itrgt,jtrgt)
							dl = Rp(i)-Rp(i-1)
						ELSEIF (sl2.le.maxbedslope(i,j)) THEN
							itrgt=i+1
							jtrgt=j
							dbed = zb_all(i,j)-zb_all(itrgt,jtrgt)
							dl = Rp(i+1)-Rp(i)							
						ELSEIF (sl3.le.maxbedslope(i,j)) THEN
							itrgt=i
							jtrgt=j-1
							dbed = zb_all(i,j)-zb_all(itrgt,jtrgt)
							dl = Rp(i)*phip(j)-Rp(i)*phip(j-1)							
						ELSEIF (sl4.le.maxbedslope(i,j)) THEN
							itrgt=i
							jtrgt=j+1
							dbed = zb_all(i,j)-zb_all(itrgt,jtrgt)
							dl = Rp(i)*phip(j+1)-Rp(i)*phip(j)
						ELSEIF (sl5.le.maxbedslope(i,j)) THEN
							itrgt=i-1
							jtrgt=j+1
							dbed = zb_all(i,j)-zb_all(itrgt,jtrgt)
							dl = SQRT((Rp(i)*phip(j)-Rp(i)*phip(jtrgt))**2+(Rp(i)-Rp(itrgt))**2)	
						ELSEIF (sl6.le.maxbedslope(i,j)) THEN
							itrgt=i+1
							jtrgt=j+1
							dbed = zb_all(i,j)-zb_all(itrgt,jtrgt)
							dl = SQRT((Rp(i)*phip(j)-Rp(i)*phip(jtrgt))**2+(Rp(i)-Rp(itrgt))**2)		
						ELSEIF (sl7.le.maxbedslope(i,j)) THEN
							itrgt=i+1
							jtrgt=j-1
							dbed = zb_all(i,j)-zb_all(itrgt,jtrgt)
							dl = SQRT((Rp(i)*phip(j)-Rp(i)*phip(jtrgt))**2+(Rp(i)-Rp(itrgt))**2)
						ELSEIF (sl8.le.maxbedslope(i,j)) THEN
							itrgt=i-1
							jtrgt=j-1
							dbed = zb_all(i,j)-zb_all(itrgt,jtrgt)
							dl = SQRT((Rp(i)*phip(j)-Rp(i)*phip(jtrgt))**2+(Rp(i)-Rp(itrgt))**2)								
						ELSE
							write(*,*),'Warning avalanche cell not found,i,j:',i,j
							CYCLE 
						ENDIF
						dbed_allowed = dl/MAX(av_slope(i,j,kbed(i,j))*bednotfixed(i,j,kbed(i,j)),1.e-18)
						dbed_adjust = vol_Vp(itrgt,jtrgt)/(vol_Vp(i,j)+vol_Vp(itrgt,jtrgt))*(dbed - dbed_allowed)
						dz_botlayer =SUM(cbotnew(1:nfrac,i,j))/cfixedbed*dz
      ! write(*,*),'rank,i,j:',rank,i,j,dbed,dbed_allowed,dbed_adjust,dz_botlayer,maxbedslope(i,j),av_slope(kbed(i,j))	
						IF (dbed_adjust.le.dz_botlayer) THEN ! only cbotnew adjusted
							!DO n1=1,nfr_sand !avalanche sand fractions only, not silt
								!n=nfrac_sand(n1)
							DO n=1,nfrac !avalanche all fractions, sand and silt
								c_adjust = dbed_adjust/MAX(dz_botlayer,1.e-18)*cbotnew(n,i,j)
								d_cbotnew(n,i,j)=d_cbotnew(n,i,j) - c_adjust
								d_cbotnew(n,itrgt,jtrgt)=d_cbotnew(n,itrgt,jtrgt)+c_adjust*vol_Vp(i,j)/vol_Vp(itrgt,jtrgt) ! dump all avalanche in cbotnew, next timestep it can be added to fixed bed in routine above
							ENDDO
							
						ELSE !avalanche full cbotnew
							kplus = MIN(kbed(i,j)+1,k1)
							drdt(i,j,kbed(i,j))=rho_b
							rnew(i,j,kbed(i,j))=rho_b
							rold(i,j,kbed(i,j))=rho_b							
							!DO n1=1,nfr_sand !avalanche sand fractions only, not silt
							!	n=nfrac_sand(n1)
							DO n=1,nfrac !avalanche all fractions, sand and silt
								c_adjust = cbotnew(n,i,j) !avalanche complete cbotnew 
								d_cbotnew(n,i,j)=d_cbotnew(n,i,j) - c_adjust
								d_cbotnew(n,itrgt,jtrgt)=d_cbotnew(n,itrgt,jtrgt) + c_adjust*vol_Vp(i,j)/vol_Vp(itrgt,jtrgt)
								c_adjust = MIN(dz,dbed_adjust-dz_botlayer)/dz*Clivebed(n,i,j,kbed(i,j)) !never more than one dz layer is eroded by avalanche
								d_cbotnew(n,i,j) = d_cbotnew(n,i,j) + Clivebed(n,i,j,kbed(i,j)) - c_adjust ! erosion one layer dz, after that avalanche
								Clivebed(n,i,j,kbed(i,j))=0. ! not bed anymore but fluid --> old Clivebed is added tot cbotnew [sediment budget OK]
								d_cbotnew(n,itrgt,jtrgt)=d_cbotnew(n,itrgt,jtrgt) + c_adjust*vol_Vp(i,j)/vol_Vp(itrgt,jtrgt) ! dump all avalanche in cbotnew, next timestep it can be added to fixed bed in routine above	
								ccnew(n,i,j,kbed(i,j))= 0. !start with fluid cell without sediment concentration
								drdt(i,j,kbed(i,j)) = drdt(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
								rnew(i,j,kbed(i,j)) = rnew(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density
								rold(i,j,kbed(i,j)) = rold(i,j,kbed(i,j))+ccnew(n,i,j,kbed(i,j))*(frac(n)%rho-rho_b) ! prevent large source in pres-corr by sudden increase in density					
							ENDDO	
							kbed(i,j)=kbed(i,j)-1  !update bed level at end	
							kbedt(i,j)=kbed(i,j)
							have_avalanched=have_avalanched+1.	
						ENDIF
					ENDIF
				ENDDO
			ENDDO
			!write(*,*),'rank,have_avalanched A:',rank,have_avalanched
			have_avalanched_tmp=have_avalanched
			call mpi_allreduce(have_avalanched_tmp,have_avalanched,1,mpi_double_precision,mpi_max,mpi_comm_world,ierr)
			
			call bound_cbot_integer(kbed) ! apply correct boundary conditions for updated kbed
			DO i=0,i1
				DO j=0,j1
					zbed(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(cbotnew(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
				ENDDO
			ENDDO			
			
!! mpi transfer sum d_cbotnew over edges:
			call shiftf_lreverse(d_cbotnew,cbf) 
			call shiftb_lreverse(d_cbotnew,cbb) 

			if (periodicy.eq.0.or.periodicy.eq.2) then
				if (rank.eq.0) then ! boundaries in j-direction
				  do n=1,nfrac
				   do i=1,imax
					   d_cbotnew(n,i,1) = d_cbotnew(n,i,1) - d_cbotnew(n,i,0) ! undo avalanche to j=0
					   d_cbotnew(n,i,jmax) = d_cbotnew(n,i,jmax) + cbb(n,i) 
				   enddo
				  enddo
				elseif (rank.eq.px-1) then
				  do n=1,nfrac
				   do i=1,imax
					   d_cbotnew(n,i,1) = d_cbotnew(n,i,1) + cbf(n,i)
					   d_cbotnew(n,i,jmax) = d_cbotnew(n,i,jmax) - d_cbotnew(n,i,j1) ! undo avalanche to j1
				   enddo
				  enddo		   
				else
				  do n=1,nfrac
				   do i=1,imax
					   d_cbotnew(n,i,1) = d_cbotnew(n,i,1) + cbf(n,i)
					   d_cbotnew(n,i,jmax) = d_cbotnew(n,i,jmax) + cbb(n,i) 
				   enddo
				  enddo
				endif
			else
			  do n=1,nfrac
			   do i=1,imax
				   d_cbotnew(n,i,1) = d_cbotnew(n,i,1) + cbf(n,i)
				   d_cbotnew(n,i,jmax) = d_cbotnew(n,i,jmax) + cbb(n,i) 
			   enddo
			  enddo
			endif
	
!! add c_cbotnew with original cbotnew:
			DO i=1,imax
				DO j=1,jmax
					DO n=1,nfrac
						cbotnew(n,i,j)=cbotnew(n,i,j)+d_cbotnew(n,i,j)
					ENDDO
				ENDDO
			ENDDO
		  ENDIF  ! end have_avalanched 
		 ENDDO !# avalanche steps
		 !write(*,*),'rank,have_avalanched:',rank,have_avalanched_tmp
		 IF (rank.eq.0.and.MOD(istep,10).eq.0) THEN
			write(*,*),'istep,avalanche steps:',istep,nav
		 ENDIF		 
		ENDIF
			! update ibm2 kbed22 which is subtly different from kbed:
		if (IBMorder.eq.2) then
			DO i=1,imax
				DO j=1,jmax			
					zb_W=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(cbotnew(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
					kbed22(i,j)=FLOOR(zb_W/dz)
				ENDDO 
			ENDDO 
			call bound_cbot_integer(kbed22)			
		else 
			kbed22=kbed
		endif 		
	ENDIF

	
      END SUBROUTINE erosion_deposition

       subroutine advec_update_Clivebed(ccnew,cbotnew,ddtt)

       implicit none
	   
      include 'mpif.h'
      integer ierr
	integer n,ib,ie,jb,je,kb,ke,k2,kk,ketmp
	REAL ccnew(1:nfrac,0:i1,0:j1,0:k1),cbotnew(1:nfrac,0:i1,0:j1),Cadvec(0:i1,0:j1,0:k1),ddtt
	REAL c_adjust,cbotnewtot,ctot_firstcel,c_adjustA,c_adjustB

	 IF ((interaction_bed.eq.4.or.interaction_bed.eq.6).and.ABS(U_TSHD).gt.1e-6) THEN ! bring Clivebed to lowest possible gridcell to adjust for advected with U_TSHD	
	ib=1
	ie=imax
	jb=1
	je=jmax
	kb=1
	ke=MAXVAL(kbed(ib:ie,jb:je)) !kmax
	ketmp=ke
	call mpi_allreduce(ketmp,ke,1,mpi_integer,mpi_max,mpi_comm_world,ierr)
	DO n=1,nfrac
		Cadvec=0.
		call adveccbot3d_TVD(Cadvec(:,:,:),Clivebed(n,:,:,:),Ubot_TSHD,Vbot_TSHD,Ru,Rp,dr,phiv,phipt,dz,
     +            	i1,j1,k1,ib,ie,jb,je,kb,ke,ddtt,rank,px,periodicx,periodicy)
		 DO i=ib,ie
		   DO j=jb,je
			 DO k=kb,ke 
				Clivebed(n,i,j,k)= Clivebed(n,i,j,k) + ddtt*Cadvec(i,j,k) ! time update EE1 (stable and conserving TVD advec scheme)
			 ENDDO
		  ENDDO
		ENDDO
	ENDDO

		 DO i=1,imax
		   DO j=1,jmax
			 DO k=1,ke !kmax 
				DO k2=k+1,ke !kmax
				  c_adjust= MAX(MIN(cfixedbed-SUM(Clivebed(1:nfrac,i,j,k)),SUM(Clivebed(1:nfrac,i,j,k2))),0.)
     &					/(MAX(SUM(Clivebed(1:nfrac,i,j,k2)),1e-18))
				  DO n=1,nfrac
					Clivebed(n,i,j,k)=Clivebed(n,i,j,k)+c_adjust*Clivebed(n,i,j,k2)      !top op Clivebed(k) up to cfixedbed from all above cells
					Clivebed(n,i,j,k2)=Clivebed(n,i,j,k2)-c_adjust*Clivebed(n,i,j,k2)
				  ENDDO
				ENDDO
			 ENDDO
		 
			 IF (SUM(Clivebed(1:nfrac,i,j,1)).ge.cfixedbed) THEN
			  DO k=1,ke !kmax
			   IF (SUM(Clivebed(1:nfrac,i,j,k)).ge.cfixedbed.and.SUM(Clivebed(1:nfrac,i,j,k+1)).lt.cfixedbed) THEN
!			   IF (SUM(Clivebed(1:nfrac,i,j,k))+SUM(cbotnew(1:nfrac,i,j)).ge.cfixedbed.and.
!     &			   SUM(Clivebed(1:nfrac,i,j,k+1))+SUM(cbotnew(1:nfrac,i,j)).lt.cfixedbed) THEN 
			   
			      IF (k>kbed(i,j)) THEN
				    DO n=1,nfrac
					  cbotnew(n,i,j)=cbotnew(n,i,j)+SUM(ccnew(n,i,j,kbed(i,j):k))  ! add all suspended sediment of fluidcells now covered inside bed to cbotnew
					  ccnew(n,i,j,kbed(i,j):k)=0.
					ENDDO
					drdt(i,j,kbed(i,j):k) = rho_b  ! prevent large source in pres-corr by sudden increase in density
					rnew(i,j,kbed(i,j):k) = rho_b  ! prevent large source in pres-corr by sudden increase in density
					rold(i,j,kbed(i,j):k) = rho_b  ! prevent large source in pres-corr by sudden increase in density						
				  ENDIF
				  DO n=1,nfrac
				    cbotnew(n,i,j)=cbotnew(n,i,j)+SUM(Clivebed(n,i,j,k+1:kmax))! add all partially filled bed-cells above the bed to cbotnew
					Clivebed(n,i,j,k+1:kmax)=0.
				  ENDDO
			      kbed(i,j)=k ! adjust kbed
				  kbedt(i,j)=k
				  EXIT
			   ENDIF
			  ENDDO
			 ELSE
			  DO n=1,nfrac
				  cbotnew(n,i,j)=cbotnew(n,i,j)+SUM(Clivebed(n,i,j,1:kmax))! add all partially filled bed-cells above the bed to cbotnew
				  Clivebed(n,i,j,1:kmax)=0.					  
			  ENDDO
			  kbed(i,j)=0 
			  kbedt(i,j)=0
			 ENDIF
 			 DO k=1,kbed(i,j)
			   c_adjust= MIN(cfixedbed-SUM(Clivebed(1:nfrac,i,j,k)),0.) !negative in case Clivebed contains too much sediment otherwise zero
     &					/(SUM(Clivebed(1:nfrac,i,j,k)))			 
			   DO n=1,nfrac
					cbotnew(n,i,j)=cbotnew(n,i,j)-c_adjust*Clivebed(n,i,j,k) 
					Clivebed(n,i,j,k)=Clivebed(n,i,j,k)+c_adjust*Clivebed(n,i,j,k)       
			   ENDDO
			 ENDDO 

!			 IF (SUM(cbotnew(1:nfrac,i,j)).ge.cfixedbed.and.kbed(i,j)+1.le.kmax.and.
!     &			 (SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.cfixedbed.or.kbed(i,j).eq.0)) THEN
			  DO WHILE (SUM(cbotnew(1:nfrac,i,j)).ge.(cfixedbed-1.e-12).and.kbed(i,j)+1.le.kmax.and.
     &			 (SUM(Clivebed(1:nfrac,i,j,kbed(i,j))).ge.(cfixedbed-1.e-12).or.kbed(i,j).eq.0))
				kbed(i,j)=kbed(i,j)+1
				kbedt(i,j)=kbed(i,j)
					drdt(i,j,kbed(i,j))=rho_b
					rnew(i,j,kbed(i,j))=rho_b
					rold(i,j,kbed(i,j))=rho_b	
					cbotnewtot=SUM(cbotnew(1:nfrac,i,j))
					ctot_firstcel=SUM(ccnew(1:nfrac,i,j,kbed(i,j)))
					c_adjustA = MAX(cfixedbed-ctot_firstcel,0.)/MAX(cbotnewtot,1.e-12)
					c_adjustB = MIN(cfixedbed-ctot_firstcel,0.)/MAX(ctot_firstcel,1.e-12)					
					DO n=1,nfrac 
						Clivebed(n,i,j,kbed(i,j))=ccnew(n,i,j,kbed(i,j))+c_adjustA*cbotnew(n,i,j)+c_adjustB*ccnew(n,i,j,kbed(i,j))
						cbotnew(n,i,j)=cbotnew(n,i,j)-c_adjustA*cbotnew(n,i,j)-c_adjustB*ccnew(n,i,j,kbed(i,j))
						ccnew(n,i,j,kbed(i,j))=0.
					ENDDO
!				write(*,*),'WHILE rank,i,j,kbed',rank,i,j,kbed(i,j),SUM(Clivebed(1:nfrac,i,j,kbed(i,j))),SUM(cbotnew(1:nfrac,i,j))	
			  ENDDO
			 !ENDIF

			 !IF (kbed(i,j)>1) THEN
			 IF (SUM(cbotnew(1:nfrac,i,j))>cfixedbed) THEN
				  write(*,*),'rank,i,j,kbed',rank,i,j,kbed(i,j),SUM(Clivebed(1:nfrac,i,j,kbed(i,j))),SUM(cbotnew(1:nfrac,i,j))
				  DO kk=1,kbed(i,j)+1
				   write(*,*)'Cb',kk,SUM(Clivebed(1:nfrac,i,j,kk))
				  ENDDO
			 ENDIF
		  ENDDO
		ENDDO
		call bound_cbot_integer(kbed) ! apply correct boundary conditions for updated kbed	
		DO i=0,i1
			DO j=0,j1
				zbed(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(cbotnew(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
			ENDDO
		ENDDO		
	 ENDIF	
	 


	end subroutine advec_update_Clivebed
	
       subroutine air_bubbles_free_surface

       implicit none

	integer n,t
	REAL wsed(nfrac,0:i1,0:j1,0:k1),sumWkm(0:i1,0:j1,0:k1)

	call slipvelocity(cnew,Wnew,wsed,rnew,kmax,kmax,sumWkm,0.,1.) !made dt=0,dz=1
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
