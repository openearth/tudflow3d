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


      subroutine mkgrid

      USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      real  delta(imax),Yplus,X,Rmax,dx,xx,yy,theta_V,theta_U,phi,maxh_obst
      integer ngridsteps,begin,n,ii
	  real dphi_shft,schuifphi,phiv2t(0:jmax*px+1),factor
      real cbf(0:i1)
      real cbb(0:i1)	  
c******************************************************************
      pi   = 4.0 * atan(1.0)
	
	!dphi = atan(dy/schuif_x)
	IF (LOA>0.) THEN ! ship:
	      dz   = depth     / kmax
	ELSE ! flat plate:
	      dz   = depth     / (kmax-kjet)
	ENDIF
	
	IF (z_tau_sed<0.) THEN 
	  z_tau_sed = 0.5*dz
	ENDIF 
        kbed_bc=0 ! !FLOOR(bc_obst_h/dz) ! kbed_bc no longer needed now kbed(i,j) is used
		
	ksurf_bc=kmax-FLOOR(surf_layer/dz) ! if surf_layer is zero than ksurf_bc=kmax	

	maxh_obst=0.
	DO n=1,nobst
	   maxh_obst=MAX(maxh_obst,ob(n)%height)
	ENDDO
	kmaxTSHD_ind=kmax+2 !FLOOR(Draught/dz)+CEILING(maxh_obst/dz)+3 
!	kmaxTSHD_ind=MIN(kmax+2,kmaxTSHD_ind) !used to conservatively allocate dummy_ind TSHD hull
!	kmaxTSHD_ind=2*kmaxTSHD_ind

        ii=(i1+1)*(j1+1)*(kmaxTSHD_ind)*2 ! per partition it is possible that obstacles overlap
		
	IF (LOA>0.or.nobst>0) THEN
        ALLOCATE(i_inPpuntTSHD(ii))
        ALLOCATE(j_inPpuntTSHD(ii))
        ALLOCATE(k_inPpuntTSHD(ii))
        ALLOCATE(i_inUpuntTSHD(ii))
        ALLOCATE(j_inUpuntTSHD(ii))
        ALLOCATE(k_inUpuntTSHD(ii))
        ALLOCATE(i_inVpuntTSHD(ii))
        ALLOCATE(j_inVpuntTSHD(ii))
        ALLOCATE(k_inVpuntTSHD(ii))
        ALLOCATE(i_inWpuntTSHD(ii))
        ALLOCATE(j_inWpuntTSHD(ii))
        ALLOCATE(k_inWpuntTSHD(ii))
        !ALLOCATE(facIBMu(ii))
		!ALLOCATE(facIBMv(ii))
		!ALLOCATE(facIBMw(ii))
	ENDIF
	
	Ru(0)=Rmin

	ngridsteps=0
	DO WHILE (imax_grid(ngridsteps+1).NE.0)
	  ngridsteps=ngridsteps+1
	END DO

	!! Checks on input imax_grid:
	DO n=2,ngridsteps
	  IF(imax_grid(n)<imax_grid(n-1)) CALL writeerror(501)
	ENDDO
	DO n=1,ngridsteps
	  IF(dr_grid(n)<0.) CALL writeerror(502)
	  IF(lim_r_grid(n)<0) lim_r_grid(n)=dr_grid(n)
	  IF(fac_r_grid(n)<0) fac_r_grid(n)=1.	  
	ENDDO
	IF(imax_grid(ngridsteps).ne.imax) CALL writeerror(503)

	begin=1
	DO n=1,ngridsteps
	  dx=dr_grid(n)
	  dx=dx/fac_r_grid(n) !to make first grid cell desired dx
	  DO i = begin , imax_grid(n) ! build up grid with user defined dr_grid
	    dx=dx*fac_r_grid(n)
        if (fac_r_grid(n)<1.) then
            dx=max(dx,lim_r_grid(n))
        elseif (fac_r_grid(n)>1.) then
            dx=min(dx,lim_r_grid(n))
        endif		
		Ru(i)= Ru(i-1) + dx
	  ENDDO	
	  begin = imax_grid(n)+1
	ENDDO
      do i = 1 , imax
         Rp(i) = ( Ru(i) + Ru(i-1) ) / 2.0
         dr(i) = ( Ru(i) - Ru(i-1) )
      enddo
      dr(i1) = dr(imax)
      Ru(i1) = Ru(imax) + dr(i1)
      Rp(i1) = Ru(imax) + dr(i1) / 2.0
      dr(0)  = dr(1)
      Rp(0)  = Ru(0) - dr(0) / 2.0
	  
	  
	phivt(0)=0.
	ngridsteps=0
	DO WHILE (jmax_grid(ngridsteps+1).NE.0)
	  ngridsteps=ngridsteps+1
	END DO
	IF (ngridsteps.eq.1.and.jmax_grid(1).eq.1) THEN
	  jmax_grid(1)=jmax*px
	  dy_grid(1)=dy
	ENDIF
	!! Checks on input imax_grid:
	DO n=2,ngridsteps
	  IF(jmax_grid(n)<jmax_grid(n-1)) CALL writeerror(501)
	ENDDO
	DO n=1,ngridsteps
	  IF(dy_grid(n)<0.) CALL writeerror(502)
	  IF(lim_y_grid(n)<0) lim_y_grid(n)=dy_grid(n)
	  IF(fac_y_grid(n)<0) fac_y_grid(n)=1.	  
	ENDDO

	IF (sym_grid_y.eq.1) THEN
	  IF(jmax*px-2*jmax_grid(ngridsteps).ne.0.and.jmax*px-2*jmax_grid(ngridsteps).ne.1) CALL writeerror(503)
	ELSE
	  IF(jmax_grid(ngridsteps).ne.jmax*px) CALL writeerror(503)
	ENDIF
	IF (ngridsteps.eq.1.and.fac_y_grid(1).eq.1.) THEN
	   poissolver=0 !0=FFT pressure-poisson solver
	ELSE
	   poissolver=3 !3=PARDISO pressure-poisson solver
	ENDIF
	
	write(*,*) 'poissolver=', poissolver
	
	begin=1
	DO n=1,ngridsteps
	  dx=dy_grid(n)
	  dx=dx/fac_y_grid(n) !to make first grid cell desired dx
	  DO j = begin , jmax_grid(n) ! build up grid with user defined dr_grid
	    dx=dx*fac_y_grid(n)
        if (fac_y_grid(n)<1.) then
            dx=max(dx,lim_y_grid(n))
        elseif (fac_y_grid(n)>1.) then
            dx=min(dx,lim_y_grid(n))
        endif		
		phivt(j)= phivt(j-1) + atan(dx/schuif_x)
	  ENDDO	
	  begin = jmax_grid(n)+1
	ENDDO
	
	n=ngridsteps
	IF (sym_grid_y.eq.1) THEN
		IF (MOD(jmax*px,2).eq.0) THEN !even
			DO i=0,jmax_grid(n)-1
				phiv2t(i)=-phivt(jmax_grid(n)-i)
				phiv2t(i+jmax_grid(n)+1)=phivt(i+1)
			ENDDO
			phiv2t(jmax_grid(n))=phivt(0) !center
		ELSE !odd
			dphi_shft=0.5*(phivt(1)-phivt(0))
			DO i=0,jmax_grid(n)-1
				phiv2t(i)=-phivt(jmax_grid(n)-i)-dphi_shft;
				phiv2t(i+jmax_grid(n)+2)=phivt(i+1)+dphi_shft;
			ENDDO
			phiv2t(jmax_grid(n))=phivt(0)-dphi_shft    !center
			phiv2t(jmax_grid(n)+1)=phivt(0)+dphi_shft  !center
		ENDIF
		phivt=phiv2t
	ENDIF

	IF (monopile>0) THEN 
		factor = 2.*pi/(phivt(jmax*px)-phivt(0))
		phivt = phivt*factor
	ENDIF 
	
	DO j = 1 , jmax*px
		phipt((j)) = ( phivt(j) + phivt(j-1) ) / 2.0
		dphi2t((j)) = ( phivt(j) - phivt(j-1) )
	ENDDO
	dphi2t(jmax*px+1) = dphi2t(jmax*px)
	phivt(jmax*px+1) = phivt(jmax*px) + dphi2t(jmax*px+1)
	phipt(jmax*px+1) = phivt(jmax*px) + dphi2t(jmax*px+1) / 2.0
	dphi2t(0)  = dphi2t(1)
	phipt(0)  = phivt(0) - dphi2t(0) / 2.0

	IF (monopile>0) THEN 
		schuifphi=(phipt(jmax*px+1)+phipt(0))/2.+pi
	ELSE 
		schuifphi=(phipt(jmax*px+1)+phipt(0))/2.
	ENDIF 
	phivt=phivt-schuifphi
	phipt=phipt-schuifphi
	phiv(0:j1)=phivt(0+rank*jmax:j1+rank*jmax)
	phip(0:j1)=phipt(0+rank*jmax:j1+rank*jmax)
	dphi2(0:j1)=dphi2t(0+rank*jmax:j1+rank*jmax)
	
      IF (Lmix_type.eq.1) THEN
	do i=1,imax
	 do j=1,jmax
	  Lmix2(i,j) = (Rp(i)*dphi2(j)*dr(i)*dz)**(2./3.)
	  Lmix2hat(i,j) = (Rp(i)*(dphi2(j-1)+dphi2(j)+dphi2(j+1))*(dr(i-1)+dr(i)+dr(i+1))*3.*dz)**(2./3.)
	 enddo
	enddo
      ELSEIF (Lmix_type.eq.2) THEN
	do i=1,imax
	 do j=1,jmax	
           Lmix2(i,j) = ( ( (Rp(i)*dphi2(j))**2 +(dr(i))**2 + dz**2  ) / 3.0 )
		   Lmix2hat(i,j) = ( ( (Rp(i)*(dphi2(j-1)+dphi2(j)+dphi2(j+1)))**2 +(dr(i-1)+dr(i)+dr(i+1))**2 + (3.*dz)**2  ) / 3.0 )
	 enddo
	enddo
      ELSE
	CALL writeerror(504)
      ENDIF



      do i=0,i1
	do j=0,j1
	  theta_U = phip(j) !( j +(rank*jmax)-0.5*jmax*px-0.5)*dphi
	  theta_V = phiv(j) !( j +(rank*jmax)-0.5*jmax*px)*dphi
	  cos_u(j) = cos(theta_U)
	  cos_v(j) = cos(theta_V)
	  sin_u(j) = sin(theta_U)	
	  sin_v(j) = sin(theta_V)
	  xx=Ru(i)*cos(theta_U)-schuif_x
	  yy=Ru(i)*sin(theta_U)
	  azi_angle_u(i,j) = atan2(yy,xx)
	  xx=Rp(i)*cos(theta_U)-schuif_x
	  yy=Rp(i)*sin(theta_U)
	  azi_angle_p(i,j) = atan2(yy,xx)
	  xx=Rp(i)*cos(theta_V)-schuif_x
	  yy=Rp(i)*sin(theta_V)
	  azi_angle_v(i,j) = atan2(yy,xx)
	enddo
	enddo
      
	do j=0,jmax*px+1
	  theta_U = phipt(j) 
	  theta_V = phivt(j) 
	  cos_ut(j) = cos(theta_U)
	  cos_vt(j) = cos(theta_V)
	  sin_ut(j) = sin(theta_U)	
	  sin_vt(j) = sin(theta_V)
	enddo	  
	!! wave influence parameters:
	Lw=gz*Tp**2/(2.*pi) ! start with deep water wave length
	do k=1,20
	  Lw=gz*Tp**2/(2.*pi)*tanh(2.*pi*(depth-bc_obst_h)/Lw)
	enddo
	kabs_w=2.*pi/Lw
	om_w=2.*pi/Tp
	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
	kx_w=kabs_w*nx_w*cos(phi)+kabs_w*ny_w*sin(phi)   !nx_w pos in TSHD dir, ny_w pos waves from starboard (is negative y-dir!)
	ky_w=-kabs_w*ny_w*cos(phi)+kabs_w*nx_w*sin(phi) 
	do i=1,imax
	  do j=1,jmax*px
	    vol_V(i,j)=(Ru(i)-Ru(i-1))*(phivt(j)-phivt(j-1))*Rp(i)*dz 
	  enddo
	enddo
	do i=1,imax
	  do j=1,jmax
	    vol_Vp(i,j)=(Ru(i)-Ru(i-1))*(phiv(j)-phiv(j-1))*Rp(i)*dz
	  enddo
	enddo	

	call shiftf_l(vol_Vp,cbf) 
	call shiftb_l(vol_Vp,cbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
		if (rank.eq.0) then ! boundaries in j-direction
		   do i=1,imax
		   vol_Vp(i,0) = vol_Vp(i,1) 
		   vol_Vp(i,j1) =cbb(i) 
		   enddo
		elseif (rank.eq.px-1) then
		   do i=1,imax
		   vol_Vp(i,0) = cbf(i)
		   vol_Vp(i,j1) =vol_Vp(i,jmax)  
		   enddo
		else
		   do i=1,imax
		   vol_Vp(i,0) = cbf(i)
		   vol_Vp(i,j1) =cbb(i) 
		   enddo
		endif
	else
	   do i=1,imax
		   vol_Vp(i,0) = cbf(i)
		   vol_Vp(i,j1) =cbb(i) 
	   enddo
	endif
	 ! boundaries in i-direction
	if (periodicx.eq.0.or.periodicx.eq.2) then
         do j=0,j1
		   vol_Vp(0,j)    =    vol_Vp(1,j)
		   vol_Vp(i1,j)   =    vol_Vp(imax,j)
         enddo   
	else 
         do j=0,j1
		   vol_Vp(0,j)    =    vol_Vp(imax,j)
		   vol_Vp(i1,j)   =    vol_Vp(1,j)
         enddo   
	endif

	
	
      end


      subroutine fkdat
      USE nlist
	  USE netcdf
	  USE sediment
      implicit none
	  
	  include 'mpif.h'

	integer n,n2,ii,inout,load_var,load_var2
	real ampli,phi,uu,vv,xTSHD(4),yTSHD(4),z,xx,yy
	real U2,V2,z0_U,ust_U_b,z0_V,ust_V_b,interpseries
	  	integer :: ncid, rhVarId, status2, ndims, xtype,natts,status
		integer, dimension(nf90_max_var_dims) :: dimids
		REAL dummy_var(1:imax,1:px*jmax,1:kmax),dummy_var2(1:imax,1:px*jmax,1:kmax) 
		CHARACTER(len=256) :: command,dummy,restart_file(1000)
		INTEGER ressystem, io, nfound,jpx,size1,size2,size3,size4,load_Cbed
		real*8 fc_local(1:imax,1:jmax,1:kmax),fc_local_vec(imax*jmax*kmax),fc_global_vec(imax*jmax*px*kmax)
		integer tel,ios,j2,r,status4(MPI_STATUS_SIZE),ierr,obstacle(imax,jmax,kmax),t

!       include 'param.txt'
!       include 'common.txt'
           Unew =0.
           Vnew =0.
           Wnew =0.
           Uold =0.
           Vold =0.
           Wold =0.
	   dUdt =0.
	   dVdt =0.
	   dWdt =0.
        rold =rho_b
        rnew =rho_b
	drdt=rho_b
		rhoU=0.
		rhoV=0.
		rhoW=0.


	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif

	if (periodicx.eq.0.or.periodicx.eq.2) then
		cnew =0.
		cold =0.
		dcdt=0.
		IF (time_int.eq.'AB2'.or.time_int.eq.'AB3'.or.time_int.eq.'ABv') THEN
			cc=0.
			ccold=0.
			ccbot=0.
			ccoldbot=0.
		ENDIF
!		IF (surf_layer>0.) THEN
!		do i=0,i1
!	         do j=0,j1
!		  do k=0,k1
!			if (k.gt.ksurf_bc) then
!			  uu=-U_b3*cos(phi)+V_b3*sin(phi)+U_TSHD*cos(phi)
!			  vv=-V_b3*cos(phi)-U_b3*sin(phi)+U_TSHD*sin(phi)
!			  Uold(i,j,k)=uu*cos_u(j)+vv*sin_u(j)
!			  Vold(i,j,k)=vv*cos_v(j)-uu*sin_v(j)
!			endif
!		    enddo
!		   enddo
!		  enddo
!		  Unew=Uold
!		  dUdt=Uold*rold
!		  Vnew=Vold
!		  dVdt=Vold*rold
!		ENDIF

	else !periodic in x direction:
		!if (LOA>0.) then ! apply U_b as bc also for periodic sims, so to keep random fluctuations alive
		!	Unew=-U_init !new steering dpdx dynamically uses U_b without +/- switch for LOA>0. 
		!else
			Unew=U_init
		!endif
		Vnew=V_init
		Wnew=W_b !Wnew*sqrt(dpdx*(depth-bc_obst_h)/rho_b)*4. !scale with 20u_tau
		Uold=Unew
		Vold=Vnew
		Wold=Wnew
		dUdt=Unew*rnew
		dVdt=Vnew*rnew
		dWdt=Wnew*rnew
	  	do n=1,nfrac
			cnew(n,:,:,:) = frac(n)%c
			cold(n,:,:,:) = frac(n)%c
			dcdt(n,:,:,:) = frac(n)%c
			IF (time_int.eq.'AB2'.or.time_int.eq.'AB3'.or.time_int.eq.'ABv') THEN
				cc(n,:,:,:) = frac(n)%c
				ccold(n,:,:,:) = frac(n)%c
				ccbot=0.
				ccoldbot=0.
			ENDIF
	  	enddo
	endif

	IF (nbedplume>0) THEN

	  ! start with Log-profile initial condition when bedplume:
        IF (.false.) THEN ! next bit is not needed 3-11-2014
	  DO k=0,k1
	    ust_U_b=MAX(ABS(U_b),1.e-6)
	    ust_V_b=MAX(ABS(V_b),1.e-6)
	    IF (slip_bot.eq.1) THEN
	      do ii=1,10
	        z0_U=0.11*nu_mol/MAX(ABS(ust_U_b),1.e-9)+kn/30
	        ust_U_b=ABS(U_b)*kappa/(log((depth)/z0_U)-1);
	        z0_V=0.11*nu_mol/MAX(ABS(ust_V_b),1.e-9)+kn/30
	        ust_V_b=ABS(V_b)*kappa/(log((depth)/z0_V)-1);
	      enddo
	      z=k*dz-0.5*dz
	      IF (LOA>0.) THEN
	        U2=-ust_U_b/kappa*log(z/z0_U)*signU_b*cos(phi)+ust_V_b/kappa*log(z/z0_V)*signV_b*sin(phi)+U_TSHD*cos(phi)
	        V2=-ust_V_b/kappa*log(z/z0_V)*signV_b*cos(phi)-ust_U_b/kappa*log(z/z0_U)*signU_b*sin(phi)+U_TSHD*sin(phi)
	      ELSE
	        U2=ust_U_b/kappa*log(z/z0_U)*signU_b
	        V2=ust_V_b/kappa*log(z/z0_V)*signV_b
	      ENDIF
	      DO j=0,j1 
		  DO i=0,i1		
		  	Unew(i,j,k)=U2*cos_u(j)+V2*sin_u(j)
		  	Vnew(i,j,k)=-U2*sin_v(j)+V2*cos_v(j) 
		  	Uold(i,j,k)=U2*cos_u(j)+V2*sin_u(j)
		  	Vold(i,j,k)=-U2*sin_v(j)+V2*cos_v(j) 
		  	dUdt(i,j,k)=(U2*cos_u(j)+V2*sin_u(j) ) *rnew(i,j,k)
		  	dVdt(i,j,k)=(-U2*sin_v(j)+V2*cos_v(j)) *rnew(i,j,k)
			Wnew(i,j,k)=W_b
			Wold(i,j,k)=W_b
			dWdt(i,j,k)=W_b*rnew(i,j,k)
		  ENDDO
		ENDDO
	    ELSE
		IF (LOA>0.) THEN
	    	  U2=-U_b*cos(phi)+V_b*sin(phi)+U_TSHD*cos(phi)
	    	  V2=-V_b*cos(phi)-U_b*sin(phi)+U_TSHD*sin(phi)
		ELSE
	    	  U2=U_b
	    	  V2=V_b
		ENDIF
		DO j=0,j1 
		  DO i=0,i1		
		  	Unew(i,j,k)=U2*cos_u(j)+V2*sin_u(j)
		  	Vnew(i,j,k)=-U2*sin_v(j)+V2*cos_v(j) 
		  	Uold(i,j,k)=U2*cos_u(j)+V2*sin_u(j)
		  	Vold(i,j,k)=-U2*sin_v(j)+V2*cos_v(j) 
		  	dUdt(i,j,k)=(U2*cos_u(j)+V2*sin_u(j) )*rnew(i,j,k)
		  	dVdt(i,j,k)=(-U2*sin_v(j)+V2*cos_v(j))*rnew(i,j,k)
			Wnew(i,j,k)=W_b
			Wold(i,j,k)=W_b
			dWdt(i,j,k)=W_b*rnew(i,j,k)
		  ENDDO
		ENDDO
	    ENDIF
	  ENDDO
	ENDIF ! this bit is not needed 3-11-2014	  


	DO n=1,nbedplume
      !! Search for P,V:
      do k=0,k1
       do i=0,i1  
         do j=j1,0,-1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  IF (k.le.FLOOR(bp(n)%height/dz).and.k.ge.CEILING(bp(n)%zbottom/dz)) THEN ! obstacle: 
		xTSHD(1:4)=bp(n)%x*cos(phi)-bp(n)%y*sin(phi)
		yTSHD(1:4)=bp(n)%x*sin(phi)+bp(n)%y*cos(phi)
		if (bp(n)%radius.gt.0.) then
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
	   IF (bp(n)%t0.eq.0.) THEN
	    DO n2=1,nfrac
		Cold(n2,i,j,k)=bp(n)%c(n2)
		Cnew(n2,i,j,k)=bp(n)%c(n2)
		dCdt(n2,i,j,k)=bp(n)%c(n2)
		! rho is calculated in state called after fkdat
	    ENDDO
	   ENDIF
	    IF (bp(n)%u.ne.-99999.and.bp(n).t0.eq.0.) THEN
	      Uold(i,j,k)=bp(n)%u*cos_u(j)+bp(n)%v*sin_u(j)
		  Uold(MAX(i-1,0),j,k)=bp(n)%u*cos_u(j)+bp(n)%v*sin_u(j)
	      Vold(i,j,k)=-bp(n)%u*sin_v(j)+bp(n)%v*cos_v(j) 
		  Vold(i,MAX(j-1,0),k)=-bp(n)%u*sin_v(MAX(j-1,0))+bp(n)%v*cos_v(MAX(j-1,0)) 
	      Wold(i,j,k)=bp(n)%w
		  Wold(i,j,MAX(k-1,0))=bp(n)%w
	      Unew(i,j,k)=Uold(i,j,k)
	      Vnew(i,j,k)=Vold(i,j,k)
	      Wnew(i,j,k)=bp(n)%w
	      dUdt(i,j,k)=Uold(i,j,k)*rnew(i,j,k)
		  dUdt(MAX(i-1,0),j,k)=Uold(MAX(i-1,0),j,k)*rnew(MAX(i-1,0),j,k)
	      dVdt(i,j,k)=Vold(i,j,k)*rnew(i,j,k)
		  dVdt(i,MAX(j-1,0),k)=Vold(i,MAX(j-1,0),k)*rnew(i,MAX(j-1,0),k)
	      dWdt(i,j,k)=bp(n)%w*rnew(i,j,k)
		  dWdt(i,j,MAX(k-1,0))=bp(n)%w*rnew(i,j,MAX(k-1,0))
	    ENDIF
	   endif
	  enddo
	 enddo
	enddo
	ENDDO ! bedplume loop
	
	
	obstacle=0
	  do t=1,tmax_inPpuntTSHD
 	    i=i_inPpuntTSHD(t)
 	    j=j_inPpuntTSHD(t)		
 	    k=k_inPpuntTSHD(t)		
		IF (i.ge.1.and.i.le.imax.and.j.ge.1.and.j.le.jmax.and.k.ge.1.and.k.le.kmax) THEN
			obstacle(i,j,k)=1
		ENDIF
	  enddo
      
		do j=1,jmax
			do i=1,imax
				do k=1,kbed(i,j)
					obstacle(i,j,k)=1
				enddo
			enddo
		enddo	

		do i=1,imax
			do j=1,jmax
				do k=1,kmax	
					fc_local(i,j,k)=1.-DBLE(obstacle(i,j,k))
					tel=k+(j-1)*kmax+(i-1)*kmax*jmax
					fc_local_vec(tel)=fc_local(i,j,k)
				enddo 
			enddo
		enddo
					
		call MPI_Allgather(fc_local_vec,imax*jmax*kmax,MPI_REAL8,fc_global_vec,imax*jmax*kmax,MPI_REAL8,MPI_COMM_WORLD,status4,ierr)
		do r=1,px
			do i=1,imax
				do j=1,jmax
					do k=1,kmax	
						tel=k+(j-1)*kmax+(i-1)*kmax*jmax+(r-1)*(imax*jmax*kmax) 
						j2=j+(r-1)*jmax
						fc_global(i,j2,k)=fc_global_vec(tel)
					enddo 
				enddo
			enddo	
		enddo
		if (periodicy.eq.1) then
		  fc_global(1:imax,0,1:kmax)=fc_global(1:imax,jmax*px,1:kmax)
		  fc_global(1:imax,jmax*px+1,1:kmax)=fc_global(1:imax,1,1:kmax)
		else
		  fc_global(1:imax,0,1:kmax)=fc_global(1:imax,1,1:kmax)
		  fc_global(1:imax,jmax*px+1,1:kmax)=fc_global(1:imax,jmax*px,1:kmax)		
		endif
		if (periodicx.eq.1) then
		  fc_global(0,0:jmax*px+1,1:kmax)=fc_global(imax,0:jmax*px+1,1:kmax)
		  fc_global(i1,0:jmax*px+1,1:kmax)=fc_global(1,0:jmax*px+1,1:kmax)		
		else
		  fc_global(0,0:jmax*px+1,1:kmax)=fc_global(1,0:jmax*px+1,1:kmax)
		  fc_global(i1,0:jmax*px+1,1:kmax)=fc_global(imax,0:jmax*px+1,1:kmax)
		endif
		
	DO n=1,nbedplume
	  bp(n)%volncells=0.
      do k=1,kmax
       do i=1,imax 
         do j=1,jmax*px       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
	  xx=Rp(i)*cos_ut(j)-schuif_x !global xx over different processors
	  yy=Rp(i)*sin_ut(j)          !global yy over different processors
	  IF (k.le.FLOOR(bp(n)%height/dz).and.k.ge.CEILING(bp(n)%zbottom/dz)) THEN ! obstacle: 
		xTSHD(1:4)=bp(n)%x*cos(phi)-bp(n)%y*sin(phi)
		yTSHD(1:4)=bp(n)%x*sin(phi)+bp(n)%y*cos(phi)
		if (bp(n)%radius.gt.0.) then
		  inout=0
		  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n)%radius)**2) THEN
			inout=1
		  ENDIF
		else 
		  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		endif 		
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
	   bp(n)%volncells=bp(n)%volncells+vol_V(i,j)*fc_global(i,j,k)
	   endif
	  enddo
	 enddo
	enddo	  
	 if (bp(n)%volncells.le.0.) then
	   write(*,*),'WARNING, no cells found for bedplume number: ',n,bp(n)%volncells,rank
	   write(*,*),'In case a sedflux or Q has been defined the code will crash because of division by a zero volncells'
	 endif
	ENDDO ! bedplume loop	

	ENDIF

	if (restart_dir.ne.'') then
		WRITE(*,*) 'Listing restart files with system call'
		WRITE(command,'(a,a,a,i4.4,a)')'ls ',TRIM(restart_dir),' > restart_temp',INT(rank),'.txt'
		CALL SYSTEM(command) 
		
		WRITE(command,'(a,i4.4,a)')'restart_temp',INT(rank),'.txt'
		OPEN(unit=25,file=command) 
		
		nfound = 0
		DO
			READ(25,'(a)',iostat=io) dummy
			IF (io/=0) EXIT
			nfound = nfound + 1
			restart_file(nfound)=TRIM(dummy)
		ENDDO
		IF (rank.eq.0) THEN
			WRITE(*,'(A,I0)') ' Number of restart files found: ',nfound	
			DO n=1,nfound 
			  write(*,*),restart_file(n)
			ENDDO
		ENDIF
		CLOSE(25, STATUS = 'DELETE')
		!CLOSE(25, STATUS = 'KEEP')
		
		jpx = (jmax*px/nfound)

		IF (nfound.gt.0) THEN
			n2=1 
			status2 = nf90_open(restart_file(n2), nf90_NoWrite, ncid) 
			IF (status2/= nf90_noerr) THEN
				write(*,*),'initconditionsfile =',restart_file(n2)
				CALL writeerror(701)
			ENDIF
			status = nf90_inq_varid(ncid, "time",rhVarid)
			if (status.eq.nf90_NoErr) then
				call check( nf90_get_var(ncid,rhVarid,time_n))
				IF (rank.eq.0) write(*,*),'restart time =',time_n
				trestart=time_n
			endif
			call check( nf90_close(ncid) )
		ENDIF
		load_var=0
		DO n2=1,nfound
			status2 = nf90_open(restart_file(n2), nf90_NoWrite, ncid) 
			IF (status2/= nf90_noerr) THEN
				write(*,*),'initconditionsfile =',restart_file(n2)
				CALL writeerror(701)
			ENDIF
			IF (rank.eq.0) write(*,*),'reading initconditionsfile =',restart_file(n2)
			status = nf90_inq_varid(ncid, "U",rhVarid)
			if (status.eq.nf90_NoErr) then
				call check( nf90_inquire_variable(ncid, rhVarid, dimids = dimIDs))
				call check( nf90_inquire_dimension(ncid, dimIDs(1), len = size1))
				call check( nf90_inquire_dimension(ncid, dimIDs(2), len = size2))
				call check( nf90_inquire_dimension(ncid, dimIDs(3), len = size3))
				IF(size1.ne.imax.or.size2*nfound.ne.jmax*px.or.size3.ne.kmax) CALL writeerror(702)
				call check( nf90_get_var(ncid,rhVarid,dummy_var(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1:kmax),
     &                       start=(/1,1,1/),count=(/imax,jpx,kmax/)) )
				load_var=1
			else
				write(*,*),'initconditionsfile =',restart_file(n2),' variable "U" not found and not used as initial condition'
				load_var=0
				!dummy_var(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1:kmax)=0.
			endif
			call check( nf90_close(ncid) )
		ENDDO
		load_var2=0
		DO n2=1,nfound
			status2 = nf90_open(restart_file(n2), nf90_NoWrite, ncid) 
			IF (status2/= nf90_noerr) THEN
				write(*,*),'initconditionsfile =',restart_file(n2)
				CALL writeerror(701)
			ENDIF
			status = nf90_inq_varid(ncid, "V",rhVarid)
			if (status.eq.nf90_NoErr) then
				call check( nf90_inquire_variable(ncid, rhVarid, dimids = dimIDs))
				call check( nf90_inquire_dimension(ncid, dimIDs(1), len = size1))
				call check( nf90_inquire_dimension(ncid, dimIDs(2), len = size2))
				call check( nf90_inquire_dimension(ncid, dimIDs(3), len = size3))
				IF(size1.ne.imax.or.size2*nfound.ne.jmax*px.or.size3.ne.kmax) CALL writeerror(702)
				call check( nf90_get_var(ncid,rhVarid,dummy_var2(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1:kmax),
     &                       start=(/1,1,1/),count=(/imax,jpx,kmax/)) )
				load_var2=1
			else
				write(*,*),'initconditionsfile =',restart_file(n2),' variable "V" not found and not used as initial condition'
				load_var2=0
				!dummy_var2(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1:kmax)=0.
			endif
			call check( nf90_close(ncid) )
		ENDDO
		IF (load_var.eq.1.and.load_var2.eq.1) THEN
		 do i=1,imax
		  do j=1,jmax
		    do k=1,kmax 
			  Unew(i,j,k)=dummy_var(i,j+rank*jmax,k)*cos_u(j)+dummy_var2(i,j+rank*jmax,k)*sin_u(j)
			  Vnew(i,j,k)=dummy_var2(i,j+rank*jmax,k)*cos_v(j)-dummy_var(i,j+rank*jmax,k)*sin_v(j)
			enddo 
		  enddo
		 enddo
		ENDIF
		load_var=0
		DO n2=1,nfound
			status2 = nf90_open(restart_file(n2), nf90_NoWrite, ncid) 
			IF (status2/= nf90_noerr) THEN
				write(*,*),'initconditionsfile =',restart_file(n2)
				CALL writeerror(701)
			ENDIF
			status = nf90_inq_varid(ncid, "W",rhVarid)
			if (status.eq.nf90_NoErr) then
				call check( nf90_inquire_variable(ncid, rhVarid, dimids = dimIDs))
				call check( nf90_inquire_dimension(ncid, dimIDs(1), len = size1))
				call check( nf90_inquire_dimension(ncid, dimIDs(2), len = size2))
				call check( nf90_inquire_dimension(ncid, dimIDs(3), len = size3))
				IF(size1.ne.imax.or.size2*nfound.ne.jmax*px.or.size3.ne.kmax) CALL writeerror(702)
				call check( nf90_get_var(ncid,rhVarid,dummy_var(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1:kmax),
     &                       start=(/1,1,1/),count=(/imax,jpx,kmax/)) )
				load_var=1
			else
				write(*,*),'initconditionsfile =',restart_file(n2),' variable "W" not found and not used as initial condition'
				load_var=0
				!dummy_var(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1:kmax)=0.
			endif
			call check( nf90_close(ncid) )
		ENDDO
		IF (load_var.eq.1) THEN
		 do i=1,imax
		  do j=1,jmax
		    do k=1,kmax 
			  Wnew(i,j,k)=dummy_var(i,j+rank*jmax,k)
			enddo 
		  enddo
		 enddo
		ENDIF
		load_var=0
		DO n=1,nfrac
			DO n2=1,nfound
				status2 = nf90_open(restart_file(n2), nf90_NoWrite, ncid) 
				IF (status2/= nf90_noerr) THEN
					write(*,*),'initconditionsfile =',restart_file(n2)
					CALL writeerror(701)
				ENDIF
				status = nf90_inq_varid(ncid, "C",rhVarid)
				if (status.eq.nf90_NoErr) then
					call check( nf90_inquire_variable(ncid, rhVarid, dimids = dimIDs))
					call check( nf90_inquire_dimension(ncid, dimIDs(1), len = size1))
					call check( nf90_inquire_dimension(ncid, dimIDs(2), len = size2))
					call check( nf90_inquire_dimension(ncid, dimIDs(3), len = size3))
					call check( nf90_inquire_dimension(ncid, dimIDs(4), len = size4))
					IF(size1.ne.nfrac.or.size2.ne.imax.or.size3*nfound.ne.jmax*px.or.size4.ne.kmax) CALL writeerror(702)
					call check( nf90_get_var(ncid,rhVarid,dummy_var(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1:kmax),
     &                       start=(/n,1,1,1/),count=(/1,imax,jpx,kmax/)) )
					load_var=1
				else
					write(*,*),'initconditionsfile =',restart_file(n2),' variable "C" not found and not used as initial condition'
					load_var=0
					!dummy_var(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1:kmax)=0.
				endif
				call check( nf90_close(ncid) )
			ENDDO
			IF (load_var.eq.1) THEN
			 do i=1,imax
			  do j=1,jmax
				do k=1,kmax 
				  Cnew(n,i,j,k)=dummy_var(i,j+rank*jmax,k)
				enddo 
			  enddo
			 enddo
			ENDIF
		ENDDO
		load_Cbed=0
		load_var=0
		IF (interaction_bed.ge.4) THEN
		 DO n=1,nfrac
			DO n2=1,nfound
				status2 = nf90_open(restart_file(n2), nf90_NoWrite, ncid) 
				IF (status2/= nf90_noerr) THEN
					write(*,*),'initconditionsfile =',restart_file(n2)
					CALL writeerror(701)
				ENDIF
				status = nf90_inq_varid(ncid, "Cbed",rhVarid)
				if (status.eq.nf90_NoErr) then
					load_Cbed=1
					call check( nf90_inquire_variable(ncid, rhVarid, dimids = dimIDs))
					call check( nf90_inquire_dimension(ncid, dimIDs(1), len = size1))
					call check( nf90_inquire_dimension(ncid, dimIDs(2), len = size2))
					call check( nf90_inquire_dimension(ncid, dimIDs(3), len = size3))
					call check( nf90_inquire_dimension(ncid, dimIDs(4), len = size4))
					IF(size1.ne.nfrac.or.size2.ne.imax.or.size3*nfound.ne.jmax*px.or.size4.ne.kmax) CALL writeerror(702)
					call check( nf90_get_var(ncid,rhVarid,dummy_var(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1:kmax),
     &                       start=(/n,1,1,1/),count=(/1,imax,jpx,kmax/)) )
					load_var=1
				else
					write(*,*),'initconditionsfile =',restart_file(n2),' variable "Cbed" not found and not used as initial condition'
					load_var=0
					!dummy_var(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1:kmax)=0.
				endif
				call check( nf90_close(ncid) )
			ENDDO
			IF (load_var.eq.1) THEN
			 do i=1,imax
			  do j=1,jmax
				do k=1,kmax 
				  Clivebed(n,i,j,k)=dummy_var(i,j+rank*jmax,k)
				enddo 
			  enddo
			 enddo
			ENDIF
		 ENDDO
		ENDIF
		load_var=0
		DO n=1,nfrac
			DO n2=1,nfound
				status2 = nf90_open(restart_file(n2), nf90_NoWrite, ncid) 
				IF (status2/= nf90_noerr) THEN
					write(*,*),'initconditionsfile =',restart_file(n2)
					CALL writeerror(701)
				ENDIF
				status = nf90_inq_varid(ncid, "mass_bed",rhVarid)
				if (status.eq.nf90_NoErr) then
					call check( nf90_inquire_variable(ncid, rhVarid, dimids = dimIDs))
					call check( nf90_inquire_dimension(ncid, dimIDs(1), len = size1))
					call check( nf90_inquire_dimension(ncid, dimIDs(2), len = size2))
					call check( nf90_inquire_dimension(ncid, dimIDs(3), len = size3))
					IF(size1.ne.nfrac.or.size2.ne.imax.or.size3*nfound.ne.jmax*px) CALL writeerror(702)
					call check( nf90_get_var(ncid,rhVarid,dummy_var(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1),
     &                       start=(/n,1,1/),count=(/1,imax,jpx/)) )
					load_var=1
				else
					write(*,*),'initconditionsfile =',restart_file(n2),' variable "mass_bed" not found and not used as initial condition'
					load_var=0
					!dummy_var(1:imax,(n2-1)*jpx+1:(n2-1)*jpx+jpx,1)=0.
				endif
				call check( nf90_close(ncid) )
			ENDDO
			IF (load_var.eq.1) THEN
			 do i=1,imax
			  do j=1,jmax
				  Cnewbot(n,i,j)=dummy_var(i,j+rank*jmax,1)/dz/frac(n)%rho 
			  enddo
			 enddo
			ENDIF
		ENDDO		
		
		Uold=Unew
		Vold=Vnew 
		Wold=Wnew 
		Cold=Cnew
		
		IF (interaction_bed.ge.4.and.load_Cbed.eq.1) THEN
			kbed=0
			kbedt=0
			do i=1,imax
			  do j=1,jmax
				k=1
				do WHILE (SUM(Clivebed(1:nfrac,i,j,k)).gt.0.)
				  k=k+1
				  kbed(i,j)=k
				  kbedt(i,j)=k
				enddo
			  enddo
			enddo
		  call bound_cbot_integer(kbed) 
		  call bound_cbot_integer(kbedt)	

	!		call state(cnew,rnew)
	!		CALL erosion_deposition(Cnew,cnewbot,Unew,Vnew,Wnew,Rnew,Cnew,Cnewbot,dt,dz)

		  do n=1,nfrac
			 IF (interaction_bed.ge.4) THEN
			   call bound_c(Clivebed(n,:,:,:),0.,n,0.)
			 ENDIF
			 call bound_cbot(Cnewbot(n,:,:))
		  enddo	
	  ENDIF

	
	ENDIF
	
	END 

      SUBROUTINE init_transpose
      USE nlist
!         INCLUDE 'param.txt'
!         parameter (nr=imax,mt=jmax,nx=kmax,mx=kmax/px,NT=jmax*px)
! 	integer Xii(Nx),Xkk(Nx,Mt)
! 	common /XPOSE/ Xii,Xkk
! 	integer Tkk(Nt),Tii(Nt,Mx)
! 	common /TPOSE/ Tkk,Tii

! 	integer Xii(kmax),Xkk(kmax,jmax)
! 	common /XPOSE/ Xii,Xkk
! 	integer Tkk(jmax*px),Tii(jmax*px,kmax/px)
! 	common /TPOSE/ Tkk,Tii
	integer nr,mt,nx,mx,nt

C ...  Locals
	  
      
	nr=imax
	mt=jmax
	nx=kmax
	mx=kmax/px
	NT=jmax*px

	do k = 1,Mt
	  do i = 1,Nx
	    Xii(i) = MOD(i-1,Mx)+1
	    Xkk(i,k) = INT((i-1)/Mx)*Mt+k
	  end do
	end do
	
	do i = 1,Mx
	  do k = 1,Nt
	    Tkk(k) = MOD(k-1,Mt) + 1
	    Tii(k,i) = INT((k-1)/Mt)*Mx + i
	  end do
	end do
	  
      RETURN
      END

	
	subroutine read_bc_from_coarse_sim(filenm)
      	USE netcdf
      	USE nlist	
	implicit none

	CHARACTER*256 filenm,varname	
	integer :: ncid, rhVarId, status, ncid2, rhVarId2, status2, ndims, xtype,natts,nfrc,n
	integer, dimension(nf90_max_var_dims) :: dimids
	real cbc_lateral_dummy(nfrac,1:imax,1:kmax)
	real cbc_front_dummy(nfrac,0:jmax+1,1:kmax)

	nfrc=-1
	IF (rank.eq.0) THEN
	       	status = nf90_open(filenm, nf90_NoWrite, ncid) 
		IF (status/= nf90_noerr) THEN
			write(*,*),'bcfile =',filenm
			CALL writeerror(51)
		ENDIF
		call check( nf90_inq_varid(ncid, "ubc1a",rhVarid) )
		call check( nf90_get_var(ncid,rhVarid,Ubcoarse1(1:imax,1:kmax),start=(/1,1/),count=(/imax,kmax/)) )
		call check( nf90_inq_varid(ncid, "vbc1a",rhVarid) )
		call check( nf90_get_var(ncid,rhVarid,Vbcoarse1(1:imax,1:kmax),start=(/1,1/),count=(/imax,kmax/)) )
		call check( nf90_inq_varid(ncid, "wbc1a",rhVarid) )
		call check( nf90_get_var(ncid,rhVarid,Wbcoarse1(1:imax,1:kmax),start=(/1,1/),count=(/imax,kmax/)) )
		status = nf90_inq_varid(ncid, "cbc1a",rhVarid)
		IF (status== nf90_noerr) THEN
		  call check( nf90_inquire_variable(ncid, RhVarId, varname, xtype, ndims, dimids, natts))
		  status = nf90_inquire_dimension(ncid, dimids(1), len=nfrc)
		  IF (nfrc>nfrac) THEN
		     write(*,*) 'WARNING, number of fractions in this file is larger than nfrac'
		     write(*,*) 'extra fractions in bcfilename are ignored!'
		     write(*,*) 'bcfilename : ',filenm
		     write(*,*) 'nfrac in bcfile : ',nfrc
		     write(*,*) 'nfrac in simulation : ',nfrac
	    	  ENDIF
		  IF (nfrc>0) THEN
		    call check( nf90_inq_varid(ncid, "cbc1a",rhVarid) )
		    call check( nf90_get_var(ncid,rhVarid,cbc_lateral_dummy(1:nfrc,:,:),start=(/1,1,1/),count=(/nfrc,imax,kmax/)) )
		    DO n=1,nfrc
			Cbcoarse1(n,:,:)=cbc_lateral_dummy(n,:,:)
		    ENDDO
			write(*,*),'Cbc1 lateral boundary read, nfrac=',nfrc
		  ENDIF		    
		ENDIF
		call check( nf90_close(ncid) )

!	write(*,*),'rank,Ubc1(1,1:3),Ubce1(1:3,1)',rank,Ubcoarse1(1,1:3),Ubcoarse1(1:3,1)
!     	write(*,*),'rank,Ubc1(jmax,kmax-2:kmax),Ubc1(jmax-2:jmax,kmax)',rank,Ubcoarse1(jmax,kmax-2:kmax),Ubcoarse1(jmax-2:jmax,kmax)
	ENDIF
	IF (rank.eq.px-1) THEN
	       	status = nf90_open(filenm, nf90_NoWrite, ncid) 
		IF (status/= nf90_noerr) THEN
			write(*,*),'bcfile =',filenm
			CALL writeerror(51)
		ENDIF
		call check( nf90_inq_varid(ncid, "ubc1b",rhVarid) )
		call check( nf90_get_var(ncid,rhVarid,Ubcoarse1(1:imax,1:kmax),start=(/1,1/),count=(/imax,kmax/)) )
		call check( nf90_inq_varid(ncid, "vbc1b",rhVarid) )
		call check( nf90_get_var(ncid,rhVarid,Vbcoarse1(1:imax,1:kmax),start=(/1,1/),count=(/imax,kmax/)) )
		call check( nf90_inq_varid(ncid, "wbc1b",rhVarid) )
		call check( nf90_get_var(ncid,rhVarid,Wbcoarse1(1:imax,1:kmax),start=(/1,1/),count=(/imax,kmax/)) )
		status = nf90_inq_varid(ncid, "cbc1b",rhVarid)
		IF (status== nf90_noerr) THEN
		  call check( nf90_inquire_variable(ncid, RhVarId, varname, xtype, ndims, dimids, natts))
		  status = nf90_inquire_dimension(ncid, dimids(1), len=nfrc)
		  IF (nfrc>nfrac) THEN
		     write(*,*) 'WARNING, number of fractions in this file is larger than nfrac'
		     write(*,*) 'extra fractions in bcfilename are ignored!'
		     write(*,*) 'bcfilename : ',filenm
		     write(*,*) 'nfrac in bcfile : ',nfrc
		     write(*,*) 'nfrac in simulation : ',nfrac
	    	  ENDIF
		  IF (nfrc>0) THEN
		    call check( nf90_inq_varid(ncid, "cbc1b",rhVarid) )
		    call check( nf90_get_var(ncid,rhVarid,cbc_lateral_dummy(1:nfrc,:,:),start=(/1,1,1/),count=(/nfrc,imax,kmax/)) )
		    DO n=1,nfrc
			Cbcoarse1(n,:,:)=cbc_lateral_dummy(n,:,:)
		    ENDDO
			write(*,*),'Cbc1 lateral boundary read, nfrac=',nfrc
		  ENDIF		    
		ENDIF
		call check( nf90_close(ncid) )

!	write(*,*),'rank,Ubc1(1,1:3),Ubce1(1:3,1)',rank,Ubcoarse1(1,1:3),Ubcoarse1(1:3,1)
!     	write(*,*),'rank,Ubc1(jmax,kmax-2:kmax),Ubc1(jmax-2:jmax,kmax)',rank,Ubcoarse1(jmax,kmax-2:kmax),Ubcoarse1(jmax-2:jmax,kmax)

	ENDIF

	! read ubcoarse2 (0:j1,1:kmax) -> in matlab index 0 does not exist, so netcdffile ubc2(1:j1+1,1:kmax)
	! ubc2 consists inflow for all px domains!
       	status2 = nf90_open(filenm, nf90_NoWrite, ncid2) 
	IF (status2/= nf90_noerr) THEN
		write(*,*),'bcfile =',filenm
		CALL writeerror(51)
	ENDIF

	call check( nf90_inq_varid(ncid2, "ubc2",rhVarid2) )
	call check( nf90_get_var(ncid2,rhVarid2,Ubcoarse2(0:j1,1:kmax),start=(/rank*jmax+1,1/),count=(/j1+1,kmax/)) )
	call check( nf90_inq_varid(ncid2, "vbc2",rhVarid2) )
	call check( nf90_get_var(ncid2,rhVarid2,Vbcoarse2(0:j1,1:kmax),start=(/rank*jmax+1,1/),count=(/j1+1,kmax/)) )
	call check( nf90_inq_varid(ncid2, "wbc2",rhVarid2) )
	call check( nf90_get_var(ncid2,rhVarid2,Wbcoarse2(0:j1,1:kmax),start=(/rank*jmax+1,1/),count=(/j1+1,kmax/)) )
		status = nf90_inq_varid(ncid2, "cbc2",rhVarid)
		IF (status== nf90_noerr) THEN
		  call check( nf90_inquire_variable(ncid2, RhVarId, varname, xtype, ndims, dimids, natts))
		  status = nf90_inquire_dimension(ncid2, dimids(1), len=nfrc)
		  IF (nfrc>nfrac) THEN
		     write(*,*) 'WARNING, number of fractions in this file is larger than nfrac'
		     write(*,*) 'extra fractions in bcfilename are ignored!'
		     write(*,*) 'bcfilename : ',filenm
		     write(*,*) 'nfrac in bcfile : ',nfrc
		     write(*,*) 'nfrac in simulation : ',nfrac
	    	  ENDIF
		  IF (nfrc>0) THEN
		    call check( nf90_inq_varid(ncid2, "cbc2",rhVarid) )
		    call check( nf90_get_var(ncid2,rhVarid,cbc_front_dummy(1:nfrc,:,:),start=(/1,rank*jmax+1,1/),count=(/nfrc,j1+1,kmax/)) )
		    DO n=1,nfrc
			Cbcoarse2(n,:,:)=cbc_front_dummy(n,:,:)
		    ENDDO
			write(*,*),'Cbc2 front inflow boundary read, nfrac=',nfrc
		  ENDIF		    
		ENDIF

	call check( nf90_close(ncid2) )


!	write(*,*),'rank,Ubcoarse2(0,1:3),Ubcoarse2(0:2,1)',rank,Ubcoarse2(0,1:3),Ubcoarse2(0:2,1)
!	write(*,*),'rank,Vbcoarse2(0,1:3),Vbcoarse2(0:2,1)',rank,Vbcoarse2(0,1:3),Vbcoarse2(0:2,1)
!	write(*,*),'rank,Wbcoarse2(0,1:3),Wbcoarse2(0:2,1)',rank,Wbcoarse2(0,1:3),Wbcoarse2(0:2,1)
!	write(*,*),'rank,Ubcoarse2(j1,kmax-2:kmax),Ubcoarse2(j1-2:j1,kmax)',rank,Ubcoarse2(j1,kmax-2:kmax),Ubcoarse2(j1-2:j1,kmax)
!	write(*,*),'rank,Vbcoarse2(j1,kmax-2:kmax),Vbcoarse2(j1-2:j1,kmax)',rank,Vbcoarse2(j1,kmax-2:kmax),Vbcoarse2(j1-2:j1,kmax)
!	write(*,*),'rank,Wbcoarse2(j1,kmax-2:kmax),Wbcoarse2(j1-2:j1,kmax)',rank,Wbcoarse2(j1,kmax-2:kmax),Wbcoarse2(j1-2:j1,kmax)


	end subroutine read_bc_from_coarse_sim





	subroutine determine_indices_jet_in

      USE nlist
	implicit none
! 	include 'param.txt'
! 	include 'common.txt'

	REAL xx,yy,r_orifice2,zz
	INTEGER tel1,tel2,tel3,t,inout
	INTEGER*2 j_inPpunt_dummy((i1+1)*(j1+1)*px),j_inVpunt_dummy((i1+1)*(j1+1)*px),j_inPpuntrand_dummy((i1+1)*(j1+1)*px)
	INTEGER*2 i_inPpunt_dummy((i1+1)*(j1+1)*px),i_inVpunt_dummy((i1+1)*(j1+1)*px),i_inPpuntrand_dummy((i1+1)*(j1+1)*px)
	INTEGER*2 j_inPpunt2_dummy((k1+1)*(j1+1)*px),j_inVpunt2_dummy((k1+1)*(j1+1)*px)
	INTEGER*2 k_inPpunt2_dummy((k1+1)*(j1+1)*px),k_inVpunt2_dummy((k1+1)*(j1+1)*px)
	LOGICAL in
	REAL rr,jetcorr

      
      r_orifice2=radius_j**2
      tel1=0
      tel2=0
      tel3=0
	tmax_inUpunt=0
	tmax_inVpunt=0
	tmax_inPpunt=0
	tmax_inPpuntrand=0
	
	
	if (radius_j>0.) then
		Aplume=0.
      do i=0,i1  
	in=.false.
        do j=jmax*px,1,-1 ! first search in all partitions
	  xx=Rp(i)*cos_ut(j)-schuif_x
	  yy=Rp(i)*sin_ut(j)
	  if ((xx*xx+yy*yy).le.r_orifice2.and.(xx*xx+yy*yy).gt.radius_inner_j**2 ) then
	   CALL PNPOLY (xx,yy, xj(1:4), yj(1:4), 4, inout ) ! only grid inside xj,yj box is jet (xj,yj default is complete domain)
	   if (inout.eq.1) then
	    tel1=tel1+1
	    i_inPpunt_dummy(tel1)=i ! Ppunt is U
	    j_inPpunt_dummy(tel1)=j
		
		if (Q_j>0.or.plumeQtseriesfile.ne.'') THEN
			jetcorr=pi/(2.*pi*(1/(1./W_j_powerlaw+1)-1./(1./W_j_powerlaw+2.))) !=1.22449 for W_j_powerlaw=7
			rr=1.-sqrt((xx**2+yy**2)/MAX(radius_j**2,1.e-12))
			rr=MAX(rr,0.)
			Aplume=Aplume+vol_V(i,j)/dz*jetcorr*rr**(1./W_j_powerlaw) !m2 inflow area weighted with inflow velocity profile --> W_j==Q_j/Aplume
		endif

	    if (in) then ! rand Vpunt niet binnen jet buisje, dus alleen zodra al een ppunt gevonden is dan vpunt ook vinden
	    	tel2=tel2+1
		i_inVpunt_dummy(tel2)=i ! Vpunt is V
		j_inVpunt_dummy(tel2)=j
	    else !eerste Ppunt is rand
		tel3=tel3+1
		i_inPpuntrand_dummy(tel3)=i
		j_inPpuntrand_dummy(tel3)=j
	    endif
	    in=.true.
	   endif
	  endif
	enddo

	!laatste Ppunt is ook rand punt:
	if (tel1>0) then
	  tel3=tel3+1
	  i_inPpuntrand_dummy(tel3)=i_inPpunt_dummy(tel1)
          j_inPpuntrand_dummy(tel3)=j_inPpunt_dummy(tel1)
	endif
      enddo
      tmax_inPpunt=tel1
      tmax_inVpunt=tel2

	IF (Q_j>0.or.plumeQtseriesfile.ne.'') THEN
		IF (Aplume.le.0.) THEN
			write(*,*),'WARNING no plume inflow cells are found, but Q_j or plumeQtseriesfile is defined'
		ENDIF
		IF (Q_j>0.and.plumeQtseriesfile.eq.'') THEN
			W_j=-Q_j/Aplume
		ELSE
			W_j=0.			
		ENDIF
	ENDIF
	
	  
      tel1=0
      do t=1,tmax_inPpunt ! now search for all local j-indices between 0 and j1:
	if ((j_inPpunt_dummy(t)-rank*jmax).ge.0.and.(j_inPpunt_dummy(t)-rank*jmax).le.j1) then
		tel1=tel1+1
		j_inPpunt(tel1)=j_inPpunt_dummy(t)-rank*jmax
		i_inPpunt(tel1)=i_inPpunt_dummy(t)
	endif
      enddo
      tmax_inPpunt=tel1
      tel1=0
      do t=1,tmax_inVpunt ! now search for all local j-indices between 0 and j1:
	if ((j_inVpunt_dummy(t)-rank*jmax).ge.0.and.(j_inVpunt_dummy(t)-rank*jmax).le.j1) then
		tel1=tel1+1
		j_inVpunt(tel1)=j_inVpunt_dummy(t)-rank*jmax
		i_inVpunt(tel1)=i_inVpunt_dummy(t)
	endif
      enddo
      tmax_inVpunt=tel1
      tel1=0
      do t=1,tel3 ! now search for all local j-indices between 0 and j1:
	if ((j_inPpuntrand_dummy(t)-rank*jmax).ge.0.and.(j_inPpuntrand_dummy(t)-rank*jmax).le.j1) then
		tel1=tel1+1
		j_inPpuntrand(tel1)=j_inPpuntrand_dummy(t)-rank*jmax
		i_inPpuntrand(tel1)=i_inPpuntrand_dummy(t)
	endif
      enddo
      tel3=tel1
	 

!      do i=0,i1  
!	in=.false.
!        do j=j1,0,-1
!	  xx=Rp(i)*cos_u(j)-schuif_x
!	  yy=Rp(i)*sin_u(j)
!	  if ((xx*xx+yy*yy).le.r_orifice2) then
!	    tel1=tel1+1
!	    i_inPpunt(tel1)=i ! Ppunt is U
!	    j_inPpunt(tel1)=j
!	    
!	    if (in.or.j.eq.j1) then ! rand Vpunt niet binnen jet buisje, dus alleen zodra al een ppunt gevonden is dan vpunt ook vinden
!	    	tel2=tel2+1
!		i_inVpunt(tel2)=i ! Vpunt is V
!		j_inVpunt(tel2)=j
!	    else !eerste Ppunt is rand
!		tel3=tel3+1
!		i_inPpuntrand(tel3)=i
!		j_inPpuntrand(tel3)=j
!	    endif
!	    in=.true.
!	  endif
!	enddo
!	!laatste Ppunt is ook rand punt:
!	if (tel1>0) then
!	  tel3=tel3+1
!	  i_inPpuntrand(tel3)=i_inPpunt(tel1)
!          j_inPpuntrand(tel3)=j_inPpunt(tel1)
!	endif
!!!	laatste Vpunt ook overslaan:
!!	if (in.and.j_inVpunt(tel2).gt.0) then !always if a jet-point is found an extra Vpunt must included
!!	  tel2=tel2+1
!!	  i_inVpunt(tel2)=i_inVpunt(tel2-1)
!!	  j_inVpunt(tel2)=j_inVpunt(tel2-1)-1
!!	endif
!      enddo
!      tmax_inPpunt=tel1
!      tmax_inVpunt=tel2

      !! Same, but now shifted i,j index to find all Upunt
      tel1=0
      tel2=0
      do j=0,j1
	in=.false.
        do i=i1,0,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  if ((xx*xx+yy*yy).le.r_orifice2.and.(xx*xx+yy*yy).gt.radius_inner_j**2) then
	   CALL PNPOLY (xx,yy, xj(1:4), yj(1:4), 4, inout )  ! only grid inside xj,yj box is jet (xj,yj default is complete domain)
	   if (inout.eq.1) then
	    if (in) then ! rand Upunt niet binnen jet buisje, dus eerste overslaan
	    	    tel2=tel2+1
		    i_inUpunt(tel2)=i  ! Upunt is W
		    j_inUpunt(tel2)=j
	    else ! eerste Upunt is rand Ppunt:
                tel3=tel3+1
                i_inPpuntrand(tel3)=i
                j_inPpuntrand(tel3)=j
	    endif
	    in=.true.
	   endif
	  endif
	enddo
        !laatste Upunt is ook rand punt:
	if (tel2>0) then
          tel3=tel3+1
          i_inPpuntrand(tel3)=i_inUpunt(tel2)
          j_inPpuntrand(tel3)=j_inUpunt(tel2)
	endif
!!	laatste Upunt ook overslaan:
!	if (in.and.i_inUpunt(tel2).gt.0) then !always if a jet-point is found an extra Upunt must included
!	  tel2=tel2+1
!	  i_inUpunt(tel2)=i_inUpunt(tel2-1)-1
!	  j_inUpunt(tel2)=j_inUpunt(tel2-1)
!	endif
      enddo
      tmax_inUpunt=tel2
      tmax_inPpuntrand=tel3
	endif
	!! Er zitten nu wel dubbele cellen in Ppuntrand, maar dat is is voor snelheid van code niet erg

      !! jet front side:
      r_orifice2=radius_j2**2
      tel1=0
      tel2=0
      tel3=0
	tmax_inWpunt2=0
	tmax_inVpunt2=0
	tmax_inPpunt2=0

	if (radius_j2>0.) then
      do k=0,k1   !!! for this one simulation i is k  !!!
	in=.false.
        do j=jmax*px+1,0,-1
	  zz=k*dz-0.5*dz-zjet2 !Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(0)*sin_ut(j)
	  if ((zz*zz+yy*yy).le.r_orifice2) then
	    tel1=tel1+1
	    k_inPpunt2_dummy(tel1)=k ! Ppunt is U
	    j_inPpunt2_dummy(tel1)=j

	    if (in) then ! rand Vpunt niet binnen jet buisje, dus alleen zodra al een ppunt gevonden is dan vpunt ook vinden
	    	tel2=tel2+1
		k_inVpunt2_dummy(tel2)=k ! Vpunt is V
		j_inVpunt2_dummy(tel2)=j
	    else !eerste Ppunt is rand
!		tel3=tel3+1
!		i_inPpuntrand(tel3)=i
!		j_inPpuntrand(tel3)=j
	    endif
	    in=.true.
	  endif
	enddo
      enddo
      tmax_inPpunt2=tel1
      tmax_inVpunt2=tel2
      tel1=0
      do t=1,tmax_inPpunt2 ! now search for all local j-indices between 0 and j1:
	if ((j_inPpunt2_dummy(t)-rank*jmax).ge.0.and.(j_inPpunt2_dummy(t)-rank*jmax).le.j1) then
		tel1=tel1+1
		j_inPpunt2(tel1)=j_inPpunt2_dummy(t)-rank*jmax
		k_inPpunt2(tel1)=k_inPpunt2_dummy(t)
	endif
      enddo
      tmax_inPpunt2=tel1
      tel1=0
      do t=1,tmax_inVpunt2 ! now search for all local j-indices between 0 and j1:
	if ((j_inVpunt2_dummy(t)-rank*jmax).ge.0.and.(j_inVpunt2_dummy(t)-rank*jmax).le.j1) then
		tel1=tel1+1
		j_inVpunt2(tel1)=j_inVpunt2_dummy(t)-rank*jmax
		k_inVpunt2(tel1)=k_inVpunt2_dummy(t)
	endif
      enddo
      tmax_inVpunt2=tel1


      !! Same, but now shifted i,j index to find all Upunt
      tel1=0
      tel2=0
      do j=0,j1
	in=.false.
        do k=k1,0,-1 !!! for this one simulation i is k  !!!
	  zz=k*dz-0.5*dz-zjet2
	  yy=Rp(0)*sin_u(j)
	  if ((zz*zz+yy*yy).le.r_orifice2) then
	    if (in) then ! rand Upunt niet binnen jet buisje, dus eerste overslaan
	    	    tel2=tel2+1
		    k_inWpunt2(tel2)=k 
		    j_inWpunt2(tel2)=j
	    else ! eerste Wpunt is rand Ppunt:
!                tel3=tel3+1
!                i_inPpuntrand(tel3)=i
!                j_inPpuntrand(tel3)=j
	    endif
	    in=.true.
	  endif
	enddo
        !laatste Upunt is ook rand punt:
!        tel3=tel3+1
!        i_inPpuntrand(tel3)=i_inUpunt(tel2)
!        j_inPpuntrand(tel3)=j_inUpunt(tel2)
      enddo
      tmax_inWpunt2=tel2
!      tmax_inPpuntrand=tel3
	endif
      end


	subroutine init_propeller

      USE nlist
	implicit none
      include 'mpif.h' 

	REAL xx,yy,zz,dxi,dxip
	REAL xprop2(2),yprop2(2),xprop3(2),yprop3(2),phi,uprop0
	INTEGER i_min(2),j_min(2),n,tel1,t,nprop2,found(2)
	REAL Ua,YYY,Ax,Aphi,fbx,fbphi,fbx2,fby,fby2,fbz,km
	REAL sum_profile,tel,dzz,dyy,dxi_min !,dyprop2_min(2)
!	REAL*4 Ppropx_dummy(0:i1,0:px*jmax+1,0:k1)
!	REAL*4 Ppropy_dummy(0:i1,0:px*jmax+1,0:k1)
!	REAL*4 Ppropz_dummy(0:i1,0:px*jmax+1,0:k1)
	INTEGER i_inVpunt_rudder_dummy(i1*px*k1)
	INTEGER j_inVpunt_rudder_dummy(i1*px*k1)
	INTEGER k_inVpunt_rudder_dummy(i1*px*k1)
        integer ileng,ierr,itag,status(MPI_STATUS_SIZE)

	
	
	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
	if (nprop.eq.1) then
	  yprop=0.
	endif

	xprop3(1)=xprop+xfront
	xprop3(2)=xprop+xfront
	yprop3(1)=yprop+yfront
	yprop3(2)=-yprop+yfront
	xprop2(1)=xprop3(1)*cos(phi)-yprop3(1)*sin(phi)
	yprop2(1)=xprop3(1)*sin(phi)+yprop3(1)*cos(phi)
	xprop2(2)=xprop3(2)*cos(phi)-yprop3(2)*sin(phi)
	yprop2(2)=xprop3(2)*sin(phi)+yprop3(2)*cos(phi)

	Ppropx=0.
	Ppropy=0.
	Ppropz=0.
	tmax_inVpunt_rudder=0
	sum_profile=0.
	tel=0.
	IF (rank.eq.0) THEN !only search for ship grid cells on first proc to spare memory!
	ALLOCATE(Ppropx_dummy(0:i1,0:px*jmax+1,0:k1))
	ALLOCATE(Ppropy_dummy(0:i1,0:px*jmax+1,0:k1))
	ALLOCATE(Ppropz_dummy(0:i1,0:px*jmax+1,0:k1))
	Ppropx_dummy=0.
	Ppropy_dummy=0.
	Ppropz_dummy=0.


	IF (LOA>0.) THEN

	
	uprop0=1.15*(Pprop/MAX(REAL(nprop),1.)/(MAX(rho_b*Dprop*Dprop,1.e-12)))**(1./3.) ! factor 0.7 left out because full Dprop is forced by propeller
	Ua=U_TSHD-U_b	! ambient velocity (positive flowing out of propeller) 

	tel1=0

		nprop2=0
		found=0
	    do n=1,nprop
                dxi_min=100000.*dz
		do i=1,imax
		  do j=1,jmax*px
 		    xx=Ru(i)*cos_ut(j)-schuif_x
		    yy=Ru(i)*sin_ut(j)
		    dxi=sqrt((xx-xprop2(n))**2+(yy-yprop2(n))**2)
		    IF (dxi<dxi_min.and.dxi<3.*sqrt(dr(i)**2+(Rp(i)*dphi2t(j))**2)) THEN
			  IF (found(n).eq.0) THEN
				nprop2=nprop2+1
				found(n)=1
			  ENDIF
		      dxi_min=dxi
			  i_min(nprop2)=i
			  j_min(nprop2)=j			  
		    ENDIF
		  enddo
		enddo
	    enddo

	!! search for indices nearest to y-hart of prop:

	      do n=1,nprop2
		!dyprop2_min(n)=Ru(i_min(n)+2)*dphi !minimum distance to find hart prop is dy
		do j=1,jmax*px
		  yy=Ru(i_min(n))*sin_ut(j)
		  dyy=ABS(yy-yprop2(n))
		  !dyprop2_min(n)=MIN(dyprop2_min(n),dyy)
	        enddo
	      enddo
		IF (cutter.eq.'yes') THEN
	      do n=1,nprop2
		do i=1,imax
		do j=1,jmax*px
			do k=1,kmax
			  xx=Ru(i)*cos_ut(j)-schuif_x
			  yy=Rp(i)*sin_ut(j)
			  zz=k*dz-0.5*dz-(depth-zprop)
			  dzz=zz
			  dyy=yy-yprop2(n)
			  if ((dzz*dzz+dyy*dyy).le.(Dprop/2.)**2.and.xx.ge.xprop2(n)-Dprop/2.and.xx.le.xprop2(n)+Dprop/2.) then
			    YYY    = sqrt(dzz*dzz+dyy*dyy)/(Dprop/2.)
!			    Ax    = Pprop/(pi/4.*Dprop**2*uprop0*vol_V(i,j)/(Ru(i)*(phivt(j)-phivt(j-1))*dz))
!     &                              /REAL(nprop)
			    Ax    = Pprop/(pi/4.*Dprop**2*uprop0*vol_V(i,j)/(Ru(i)*(phivt(j)-phivt(j-1))*dz))
     &                              /REAL(nprop)*dr(i)/Dprop !last term to make sure to divide Pprop over complete dx of cutter wheel, otherwise too much power is added to flow
				  fbx2 = 0. !only rotation
			      !fbx2   = 3./4.*Ax
			    if (n.eq.1) then
				    fbphi = 1./4.*Ax*YYY*1.5   * (-rot_prop)/ABS(rot_prop)
			    else
				    fbphi =-1./4.*Ax*YYY*1.5   * (-rot_prop)/ABS(rot_prop)
			    endif
					!! right rotation positive 
					!! correction 1.5 is for integral linear function YYY on circilar area
					!! fbphi=1/4 of total thrust, results in about Utangent=50%Uprop 
					!! (fbx is reduced to keep total thrust/power similar)

		 	    fby2   = -fbphi * sin(atan2(dzz,dyy))
			    if (rudder.eq.1) then 
!				! add downward force at one side of rudder and upward at other side: 
!				! force is zero at z=D/2 and max at z=-D/2 (right part, vice versa left part)
				! test show this does not work as a rudder...
!				fbz    =  fbphi * cos(atan2(dzz,dyy)) 
!     &					  + 0.4*fbx2 * SIGN(1.,dyy)*(rot_prop)/ABS(rot_prop)*(SIGN(1.,-dyy)*dzz+Dprop/2.)/Dprop
!				! add downward force bottom and upward on top: 
				fbz    =  fbphi * cos(atan2(dzz,dyy)) 
     &					  - 0.4*fbx2 * dzz/(0.5*Dprop)
			    else
				fbz    =  fbphi * cos(atan2(dzz,dyy)) 
			    endif

			    fbx   =  fbx2 * cos(phi) - fby2 * sin(phi)
			    fby   =  fbx2 * sin(phi) + fby2 * cos(phi)
			
			   Ppropx_dummy(i,j  ,k)=cos_ut(j)*fbx+sin_ut(j)*fby	
			   Ppropy_dummy(i,j  ,k)=Ppropy_dummy(i,j  ,k)
     & -0.5*sin_vt(j)*fbx+0.5*cos_vt(j)*fby
			   Ppropy_dummy(i,j-1,k)=Ppropy_dummy(i,j-1,k)
     & -0.5*sin_vt(j-1)*fbx+0.5*cos_vt(j-1)*fby
			   Ppropz_dummy(i,j,k  )=Ppropz_dummy(i,j,k)+0.5*fbz
			   Ppropz_dummy(i,j,k-1)=Ppropz_dummy(i,j,k-1)+0.5*fbz
			  endif			   
			enddo 
		enddo
	      enddo
		  enddo
		ELSE
	      do n=1,nprop2
		do j=1,jmax*px
			do k=1,kmax
			  yy=Rp(i_min(n))*sin_ut(j)
			  zz=k*dz-0.5*dz-(depth-zprop)
			  dzz=zz
			  dyy=yy-yprop2(n)
			  if ((dzz*dzz+dyy*dyy).le.(Dprop/2.)**2) then
			    YYY    = sqrt(dzz*dzz+dyy*dyy)/(Dprop/2.)
			    Ax    = Pprop/(pi/4.*Dprop**2*uprop0*vol_V(i_min(n),j)/(Ru(i_min(n))*(phivt(j)-phivt(j-1))*dz))
     &                              /REAL(nprop)!*(uprop0+Ua)/uprop0
			    fbx2   = 3./4.*Ax
				
			    if (n.eq.1) then
				    fbphi = 1./4.*Ax*YYY*1.5   * (-rot_prop)/ABS(rot_prop)
			    else
				    fbphi =-1./4.*Ax*YYY*1.5   * (-rot_prop)/ABS(rot_prop)
			    endif
					!! right rotation positive 
					!! correction 1.5 is for integral linear function YYY on circilar area
					!! fbphi=1/4 of total thrust, results in about Utangent=50%Uprop 
					!! (fbx is reduced to keep total thrust/power similar)

		 	    fby2   = -fbphi * sin(atan2(dzz,dyy))
			    if (rudder.eq.1) then 
!				! add downward force at one side of rudder and upward at other side: 
!				! force is zero at z=D/2 and max at z=-D/2 (right part, vice versa left part)
				! test show this does not work as a rudder...
!				fbz    =  fbphi * cos(atan2(dzz,dyy)) 
!     &					  + 0.4*fbx2 * SIGN(1.,dyy)*(rot_prop)/ABS(rot_prop)*(SIGN(1.,-dyy)*dzz+Dprop/2.)/Dprop
!				! add downward force bottom and upward on top: 
				fbz    =  fbphi * cos(atan2(dzz,dyy)) 
     &					  - 0.4*fbx2 * dzz/(0.5*Dprop)
			    else
				fbz    =  fbphi * cos(atan2(dzz,dyy)) 
			    endif

			    fbx   =  fbx2 * cos(phi) - fby2 * sin(phi)
			    fby   =  fbx2 * sin(phi) + fby2 * cos(phi)
			
			   Ppropx_dummy(i_min(n),j  ,k)=cos_ut(j)*fbx+sin_ut(j)*fby	
			   Ppropy_dummy(i_min(n),j  ,k)=Ppropy_dummy(i_min(n),j  ,k)
     & -0.5*sin_vt(j)*fbx+0.5*cos_vt(j)*fby
			   Ppropy_dummy(i_min(n),j-1,k)=Ppropy_dummy(i_min(n),j-1,k)
     & -0.5*sin_vt(j-1)*fbx+0.5*cos_vt(j-1)*fby
			   Ppropz_dummy(i_min(n),j,k  )=Ppropz_dummy(i_min(n),j,k)+0.5*fbz
			   Ppropz_dummy(i_min(n),j,k-1)=Ppropz_dummy(i_min(n),j,k-1)+0.5*fbz

!			  if (rudder.eq.1.and.(i_min(n)+7)<i1) then
!			    if (ABS(dyy).le.(1.2*dyprop2_min(n))) then !! Middle of prop is found for rudder
!			      if (ABS(dzz).le.(Dprop/2.)) then
!				tel1=tel1+1
!				i_inVpunt_rudder_dummy(tel1)=i_min(n)+6
!				j_inVpunt_rudder_dummy(tel1)=j
!				k_inVpunt_rudder_dummy(tel1)=k
!				tel1=tel1+1 ! rudder is 2dx long
!				i_inVpunt_rudder_dummy(tel1)=i_min(n)+7
!				j_inVpunt_rudder_dummy(tel1)=j
!				k_inVpunt_rudder_dummy(tel1)=k
!			      endif
!			    endif
			  endif
			enddo
		enddo
	      enddo		
		ENDIF
		  
		  
	ENDIF !IF LOA>0
	ENDIF !rank=0


	IF (LOA>0.) THEN
	call mpi_bcast(tel1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call mpi_bcast(i_inVpunt_rudder_dummy,i1*px*k1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call mpi_bcast(j_inVpunt_rudder_dummy,i1*px*k1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call mpi_bcast(k_inVpunt_rudder_dummy,i1*px*k1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	
	IF (rank.eq.0) THEN
	  do i=1,px-1
      	   call mpi_send(Ppropx_dummy(:,0+i*jmax:j1+i*jmax,:),(i1+1)*(j1+1)*(k1+1),MPI_REAL8,i,i+100,MPI_COMM_WORLD,status,ierr)
      	   call mpi_send(Ppropy_dummy(:,0+i*jmax:j1+i*jmax,:),(i1+1)*(j1+1)*(k1+1),MPI_REAL8,i,i+101,MPI_COMM_WORLD,status,ierr)
      	   call mpi_send(Ppropz_dummy(:,0+i*jmax:j1+i*jmax,:),(i1+1)*(j1+1)*(k1+1),MPI_REAL8,i,i+102,MPI_COMM_WORLD,status,ierr)
	  enddo
	  Ppropx(:,:,:)=Ppropx_dummy(:,0+rank*jmax:j1+rank*jmax,:)
	  Ppropy(:,:,:)=Ppropy_dummy(:,0+rank*jmax:j1+rank*jmax,:) 
	  Ppropz(:,:,:)=Ppropz_dummy(:,0+rank*jmax:j1+rank*jmax,:)
	else
		call mpi_recv(Ppropx(:,:,:),(i1+1)*(j1+1)*(k1+1),MPI_REAL8,0,rank+100,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(Ppropy(:,:,:),(i1+1)*(j1+1)*(k1+1),MPI_REAL8,0,rank+101,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(Ppropz(:,:,:),(i1+1)*(j1+1)*(k1+1),MPI_REAL8,0,rank+102,MPI_COMM_WORLD,status,ierr)
	endif


	IF (rank.eq.0) THEN
	  DEALLOCATE(Ppropx_dummy)
	  DEALLOCATE(Ppropy_dummy)
	  DEALLOCATE(Ppropz_dummy)
	ENDIF		

	!! all ranks:
          tmax_inVpunt_rudder=tel1	
	      do n=1,nprop
		      tel1=0
		      do t=1,tmax_inVpunt_rudder ! now search for all local j-indices between 0 and j1:
			if ((j_inVpunt_rudder_dummy(t)-rank*jmax).ge.0.and.(j_inVpunt_rudder_dummy(t)-rank*jmax).le.j1) then
				tel1=tel1+1
				j_inVpunt_rudder(tel1)=j_inVpunt_rudder_dummy(t)-rank*jmax
				i_inVpunt_rudder(tel1)=i_inVpunt_rudder_dummy(t)
				k_inVpunt_rudder(tel1)=k_inVpunt_rudder_dummy(t)
			endif
		      enddo
	      enddo
	      tmax_inVpunt_rudder=tel1
	ENDIF !IF LOA>0

	! in case periodic x flow then Ppropx is used as driving force:
	if (periodicx.eq.1) then
		Ppropx=dpdx
	endif
	! in case periodic y flow then Ppropy is used as driving force:
	if (periodicy.eq.1) then
		Ppropy=dpdy
	endif

	end

	subroutine init_propeller_old

      USE nlist
	implicit none

	REAL xx,yy,zz,dxi,dxip
	REAL xprop2(2),yprop2(2),phi,uprop0
	INTEGER i_min(2),j_min(2),n,tel1,t
	REAL Ua,YYY,Ax,Aphi,fbx,fbphi,fbx2,fby,fby2,fbz,km
	REAL sum_profile,tel,dzz,dyy,dxi_min !,dyprop2_min(2)
!	REAL*4 Ppropx_dummy(0:i1,0:px*jmax+1,0:k1)
!	REAL*4 Ppropy_dummy(0:i1,0:px*jmax+1,0:k1)
!	REAL*4 Ppropz_dummy(0:i1,0:px*jmax+1,0:k1)
	INTEGER*2 i_inVpunt_rudder_dummy(i1*px*k1)
	INTEGER*2 j_inVpunt_rudder_dummy(i1*px*k1)
	INTEGER*2 k_inVpunt_rudder_dummy(i1*px*k1)

	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
	if (nprop.eq.1) then
	  yprop=0.
	endif

	ALLOCATE(Ppropx_dummy(0:i1,0:px*jmax+1,0:k1))
	ALLOCATE(Ppropy_dummy(0:i1,0:px*jmax+1,0:k1))
	ALLOCATE(Ppropz_dummy(0:i1,0:px*jmax+1,0:k1))


	xprop2(1)=xprop+xfront
	xprop2(2)=xprop+xfront
	yprop2(1)=yprop+yfront
	yprop2(2)=-yprop+yfront
	xprop2(1)=xprop2(1)*cos(phi)-yprop2(1)*sin(phi)
	yprop2(1)=xprop2(1)*sin(phi)+yprop2(1)*cos(phi)
	xprop2(2)=xprop2(2)*cos(phi)-yprop2(2)*sin(phi)
	yprop2(2)=xprop2(2)*sin(phi)+yprop2(2)*cos(phi)

	Ppropx=0.
	Ppropy=0.
	Ppropz=0.
	Ppropx_dummy=0.
	Ppropy_dummy=0.
	Ppropz_dummy=0.
	sum_profile=0.
	tel=0.
	tmax_inVpunt_rudder=0

	IF (LOA>0.) THEN

	
	uprop0=1.15*(Pprop/REAL(nprop)/(rho_b*Dprop*Dprop))**(1./3.) ! factor 0.7 left out because full Dprop is forced by propeller
	Ua=U_TSHD-U_b	! ambient velocity (positive flowing out of propeller) 

	tel1=0

	    do n=1,nprop
                dxi_min=100000.*dz
		do i=1,imax
		  do j=1,j1
 		    xx=Ru(i)*cos_u(j)-schuif_x
		    yy=Ru(i)*sin_u(j)
		    dxi=sqrt((xx-xprop2(n))**2+(yy-yprop2(n))**2)
		    IF (dxi<dxi_min) THEN
		        dxi_min=dxi
			i_min(n)=i
			j_min(n)=j
		    ENDIF
		  enddo
		enddo
	    enddo

	!! search for indices nearest to y-hart of prop:

	      do n=1,nprop
		!dyprop2_min(n)=Ru(i_min(n)+2)*dphi !minimum distance to find hart prop is dy
		do j=1,jmax*px
		  yy=Ru(i_min(n))*sin_ut(j)
		  dyy=ABS(yy-yprop2(n))
		  !dyprop2_min(n)=MIN(dyprop2_min(n),dyy)
	        enddo
	      enddo

	      do n=1,nprop
		do j=1,jmax*px
			do k=1,kmax
			  yy=Rp(i_min(n))*sin_ut(j)
			  zz=k*dz-0.5*dz-(depth-zprop)
			  dzz=zz
			  dyy=yy-yprop2(n)
			  if ((dzz*dzz+dyy*dyy).le.(Dprop/2.)**2) then
			    YYY    = sqrt(dzz*dzz+dyy*dyy)/(Dprop/2.)
			    Ax    = Pprop/(pi/4.*Dprop**2*uprop0*vol_V(i_min(n),j)/(Ru(i_min(n))*(phivt(j)-phivt(j-1))*dz))
     &                              /REAL(nprop)!*(uprop0+Ua)/uprop0
			    fbx2   = 3./4.*Ax 
			    if (n.eq.1) then
				    fbphi = 1./4.*Ax*YYY*1.5   * (-rot_prop)/ABS(rot_prop)
			    else
				    fbphi =-1./4.*Ax*YYY*1.5   * (-rot_prop)/ABS(rot_prop)
			    endif
					!! right rotation positive 
					!! correction 1.5 is for integral linear function YYY on circilar area
					!! fbphi=1/4 of total thrust, results in about Utangent=50%Uprop 
					!! (fbx is reduced to keep total thrust/power similar)

		 	    fby2   = -fbphi * sin(atan2(dzz,dyy))
			    if (rudder.eq.1) then 
!				! add downward force bottom and upward on top: 
				fbz    =  fbphi * cos(atan2(dzz,dyy)) 
     &					  - 0.4*fbx2 * dzz/(0.5*Dprop)
			    else
				fbz    =  fbphi * cos(atan2(dzz,dyy)) 
			    endif

			    fbx   =  fbx2 * cos(phi) - fby2 * sin(phi)
			    fby   =  fbx2 * sin(phi) + fby2 * cos(phi)
			
			   Ppropx_dummy(i_min(n),j  ,k)=cos_ut(j)*fbx+sin_ut(j)*fby	
			   Ppropy_dummy(i_min(n),j  ,k)=Ppropy_dummy(i_min(n),j  ,k)
     & -0.5*sin_vt(j)*fbx+0.5*cos_vt(j)*fby
			   Ppropy_dummy(i_min(n),j-1,k)=Ppropy_dummy(i_min(n),j-1,k)
     & -0.5*sin_vt(j-1)*fbx+0.5*cos_vt(j-1)*fby
			   Ppropz_dummy(i_min(n),j,k  )=Ppropz_dummy(i_min(n),j,k)+0.5*fbz
			   Ppropz_dummy(i_min(n),j,k-1)=Ppropz_dummy(i_min(n),j,k-1)+0.5*fbz
			  endif
!			  if (rudder.eq.1.and.(i_min(n)+7)<i1) then
!			    if (ABS(dyy).le.(1.2*dyprop2_min(n))) then !! Middle of prop is found for rudder
!			      if (ABS(dzz).le.(Dprop/2.)) then
!				tel1=tel1+1
!				i_inVpunt_rudder_dummy(tel1)=i_min(n)+6
!				j_inVpunt_rudder_dummy(tel1)=j
!				k_inVpunt_rudder_dummy(tel1)=k
!				tel1=tel1+1 ! rudder is 2dx long
!				i_inVpunt_rudder_dummy(tel1)=i_min(n)+7
!				j_inVpunt_rudder_dummy(tel1)=j
!				k_inVpunt_rudder_dummy(tel1)=k
!			      endif
!			    endif
!			  endif
			enddo
		enddo
	      enddo
          tmax_inVpunt_rudder=tel1	
	      do n=1,nprop
		      tel1=0
		      do t=1,tmax_inVpunt_rudder ! now search for all local j-indices between 0 and j1:
			if ((j_inVpunt_rudder_dummy(t)-rank*jmax).ge.0.and.(j_inVpunt_rudder_dummy(t)-rank*jmax).le.j1) then
				tel1=tel1+1
				j_inVpunt_rudder(tel1)=j_inVpunt_rudder_dummy(t)-rank*jmax
				i_inVpunt_rudder(tel1)=i_inVpunt_rudder_dummy(t)
				k_inVpunt_rudder(tel1)=k_inVpunt_rudder_dummy(t)
			endif
		      enddo
	      enddo
	      tmax_inVpunt_rudder=tel1
	      Ppropx(:,:,:)=DBLE(Ppropx_dummy(:,0+rank*jmax:j1+rank*jmax,:))
	      Ppropy(:,:,:)=DBLE(Ppropy_dummy(:,0+rank*jmax:j1+rank*jmax,:))
	      Ppropz(:,:,:)=DBLE(Ppropz_dummy(:,0+rank*jmax:j1+rank*jmax,:))
	ENDIF !IF LOA>0

	  DEALLOCATE(Ppropx_dummy)
	  DEALLOCATE(Ppropy_dummy)
	  DEALLOCATE(Ppropz_dummy)

	! in case periodic x flow then Ppropx is used as driving force:
	if (periodicx.eq.1) then
		Ppropx=dpdx
	endif
	! in case periodic y flow then Ppropy is used as driving force:
	if (periodicy.eq.1) then
		Ppropy=dpdy
	endif

	end


	subroutine determine_indices_ship_in

      	USE netcdf
      USE nlist
	implicit none
! 	include 'param.txt'
! 	include 'common.txt'

	integer :: ncid, rhVarId, status2, ndims, xtype,natts
	integer, dimension(nf90_max_var_dims) :: dimids

	REAL xx,yy,zz,phi
	INTEGER tel1,tel2,tel3,inout,kmaxTSHD,t,n
	LOGICAL in

	REAL xnoseS(20),ynoseS(20),xnoseP(20),ynoseP(20)
	REAL xTSHD(1:42),yTSHD(1:42),xTSHD2(1:42),yTSHD2(1:42)
	REAL Utshd,Vtshd
	REAL ydh,ydh_rot2,xdh2,xdh_rot2
	INTEGER*2 j_inPpuntTSHD_dummy((i1+1)*(j1+1)*px*kmaxTSHD_ind),j_inVpuntTSHD_dummy((i1+1)*(j1+1)*px*kmaxTSHD_ind)
	INTEGER*2 j_inVpunt_tauTSHD_dummy((i1+1)*(j1+1)*px)
	INTEGER*2 i_inPpuntTSHD_dummy((i1+1)*(j1+1)*px*kmaxTSHD_ind),i_inVpuntTSHD_dummy((i1+1)*(j1+1)*px*kmaxTSHD_ind)
	INTEGER*2 i_inVpunt_tauTSHD_dummy((i1+1)*(j1+1)*px)
	INTEGER*2 k_inPpuntTSHD_dummy((i1+1)*(j1+1)*px*kmaxTSHD_ind),k_inVpuntTSHD_dummy((i1+1)*(j1+1)*px*kmaxTSHD_ind)
	INTEGER*2 k_inVpunt_tauTSHD_dummy((i1+1)*(j1+1)*px),kbed3(0:i1,0:j1)
	REAL adj_ment,zbed3(0:i1,0:j1),zbed4(0:i1+1,0:j1+1),cbb(0:i1),zbedcell,zbotcell,zb_W

	tmax_inUpuntTSHD=0
	tmax_inVpuntTSHD=0
	tmax_inWpuntTSHD=0
	tmax_inUpunt_tauTSHD=0
	tmax_inVpunt_tauTSHD=0
	tmax_inWpunt_suction=0


	IF (LOA>0.) THEN

	! build shape TSHD hull:
	!! build starboard nose shape by hyperbool
	xnoseS(1)=0.2 
	do i=2,20
	  xnoseS(i)=xnoseS(i-1)+1./20. !Lfront/20.
	enddo
	ynoseS=-1./xnoseS
	ynoseS=ynoseS-ynoseS(1)
	ynoseS=ynoseS*Breadth*0.5/(ynoseS(20)-ynoseS(1))

	xnoseS=xnoseS-0.2 
	xnoseS=xnoseS*Lfront
	!! build portside nose shape by hyperbool
	xnoseP(1)=1+0.15 ! plus 0.15 to make xnoseP and xnoseS exactly similar!!!
	do i=2,20
	  xnoseP(i)=xnoseP(i-1)-1./20. !Lfront/20.
	enddo
	ynoseP=1./xnoseP
	ynoseP=ynoseP-ynoseP(20)
	ynoseP=-ynoseP*Breadth*0.5/(ynoseP(1)-ynoseP(20)) 

	xnoseP=xnoseP-0.2 
	xnoseP=xnoseP*Lfront !xnoseP(20)

	xTSHD(1:20)=xnoseS
	xTSHD(21)=LOA
	xTSHD(22)=LOA
	xTSHD(23:42)=xnoseP
	yTSHD(1:20)=ynoseS
	yTSHD(21)=Breadth*0.5
	yTSHD(22)=-Breadth*0.5
	yTSHD(23:42)=ynoseP

	xTSHD=xTSHD+xfront
	yTSHD=yTSHD+yfront
	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
	xTSHD2=xTSHD
	yTSHD2=yTSHD
	xTSHD=xTSHD2*cos(phi)-yTSHD2*sin(phi)
	yTSHD=xTSHD2*sin(phi)+yTSHD2*cos(phi)

!	write(*,*),'xTSHD:'
!	write(*,*),xTSHD
!	write(*,*),'yTSHD:'
!	write(*,*),yTSHD

	kmaxTSHD=FLOOR(Draught/dz)
	IF (plume_z_outflow_belowsurf>0.) THEN
		kjet=CEILING(plume_z_outflow_belowsurf/dz)+1
	ELSE
		kjet=kmaxTSHD+1  
	ENDIF

!      !! Search for P,V:
!      tel1=0
!      tel2=0
!      tel3=0
!      do k=kmax-kmaxTSHD,k1   
!       do i=0,i1  
!	in=.false.
!        do j=j1,0,-1                   ! search in j dir, search on rho-loc, when rho-loc is in then also V left right are in (made zero) 
!	  xx=Rp(i)*cos_u(j)-schuif_x
!	  yy=Rp(i)*sin_u(j)
! 	  zz=MIN(k*dz,depth) !k1 would lead to zz>depth
!	  if (Hback>0.) then ! make shape shorter for slope at back:
!	    xTSHD2(21:22)=LOA-(1.+MAX((depth-DRAUGHT-zz)/Hback,-1.))*Lback+xfront
!	    xTSHD(21:22)=xTSHD2(21:22)*cos(phi)-yTSHD2(21:22)*sin(phi)
!	    yTSHD(21:22)=xTSHD2(21:22)*sin(phi)+yTSHD2(21:22)*cos(phi)
!	  else ! leave shape same length:
!	    xTSHD2(21:22)=LOA+xfront
!	    xTSHD(21:22)=xTSHD2(21:22)*cos(phi)-yTSHD2(21:22)*sin(phi)
!	    yTSHD(21:22)=xTSHD2(21:22)*sin(phi)+yTSHD2(21:22)*cos(phi)
!	  endif
!	  CALL PNPOLY (xx,yy, xTSHD, yTSHD, 42, inout ) 
!	  if (inout.eq.1) then
!	      tel1=tel1+1
!	      i_inPpuntTSHD(tel1)=i ! Ppunt
!	      j_inPpuntTSHD(tel1)=j
!  	      k_inPpuntTSHD(tel1)=k
!      	      tel2=tel2+1
!	      i_inVpuntTSHD(tel2)=i ! Vpunt
!	      j_inVpuntTSHD(tel2)=j
!	      k_inVpuntTSHD(tel2)=k
!	      if (k.eq.kmax-kmaxTSHD.and.kn_TSHD>0.and.j>0.and.j<j1.and.i>0.and.i<i1) then
!		tel3=tel3+1
!	        i_inVpunt_tauTSHD(tel3)=i ! Wpunt
!	        j_inVpunt_tauTSHD(tel3)=j
!		k_inVpunt_tauTSHD(tel3)=kmax-kmaxTSHD-1
!	      endif
!	    in=.true.
!	  endif
!	enddo

!!	add one Vpunt extra because staggered:
!	if (in.and.j_inVpuntTSHD(tel2)>0) then !always if a jet-point is found an extra Vpunt must be included, except when at border between processors
!	  tel2=tel2+1
!	  i_inVpuntTSHD(tel2)=i_inVpuntTSHD(tel2-1)
!	  j_inVpuntTSHD(tel2)=j_inVpuntTSHD(tel2-1)-1
!	  k_inVpuntTSHD(tel2)=k_inVpuntTSHD(tel2-1)
!	endif
!       enddo
!      enddo
!      tmax_inPpuntTSHD=tel1
!      tmax_inVpuntTSHD=tel2
!      tmax_inVpunt_tauTSHD=tel3
	
!      !! Search for P,V on all partitions:
      tel1=0
      tel2=0
      tel3=0
      do k=kmax-kmaxTSHD,k1   
       do i=0,i1  
	in=.false.
        do j=jmax*px+1,0,-1       ! global search in j dir, search on rho-loc, when rho-loc is in then also V left right are in (made zero) 
	  xx=Rp(i)*cos_ut(j)-schuif_x
	  yy=Rp(i)*sin_ut(j)
 	  zz=MIN(k*dz,depth) !k1 would lead to zz>depth
	  if (softnose.eq.1) then !make nose shorter gradually
		adj_ment=(1.+MAX(MIN(depth-DRAUGHT-zz,0.)/(Hfront),-1.)) !grows linear from 0. at Hfront to 1 at keel TSHD 
		adj_ment=1.-0.25*adj_ment !1 at Hfront and 0.75 at keel TSHD
		xTSHD2(1:20)=xnoseS*adj_ment+Lfront*(1.-adj_ment)+xfront ! make nose shorter
		xTSHD2(23:42)=xnoseP*adj_ment+Lfront*(1.-adj_ment)+xfront ! make nose shorter
		xTSHD=xTSHD2*cos(phi)-yTSHD2*sin(phi)
		yTSHD=xTSHD2*sin(phi)+yTSHD2*cos(phi)
	  endif

	  if (Hback>0.) then ! make shape shorter for slope at back:
	    xTSHD2(21:22)=LOA-(1.+MAX((depth-DRAUGHT-zz)/Hback,-1.))*Lback+xfront
	    xTSHD(21:22)=xTSHD2(21:22)*cos(phi)-yTSHD2(21:22)*sin(phi)
	    yTSHD(21:22)=xTSHD2(21:22)*sin(phi)+yTSHD2(21:22)*cos(phi)
	  else ! leave shape same length:
	    xTSHD2(21:22)=LOA+xfront
	    xTSHD(21:22)=xTSHD2(21:22)*cos(phi)-yTSHD2(21:22)*sin(phi)
	    yTSHD(21:22)=xTSHD2(21:22)*sin(phi)+yTSHD2(21:22)*cos(phi)
	  endif
	  CALL PNPOLY (xx,yy, xTSHD, yTSHD, 42, inout ) 
	  if (inout.eq.1) then
	      tel1=tel1+1
	      i_inPpuntTSHD_dummy(tel1)=i ! Ppunt
	      j_inPpuntTSHD_dummy(tel1)=j
  	      k_inPpuntTSHD_dummy(tel1)=k
      	      tel2=tel2+1
	      i_inVpuntTSHD_dummy(tel2)=i ! Vpunt
	      j_inVpuntTSHD_dummy(tel2)=j
	      k_inVpuntTSHD_dummy(tel2)=k
	      if (k.eq.kmax-kmaxTSHD.and.(kmax-kmaxTSHD)>0.and.kn_TSHD.ge.0.0) then
		tel3=tel3+1
	        i_inVpunt_tauTSHD_dummy(tel3)=i ! Wpunt
	        j_inVpunt_tauTSHD_dummy(tel3)=j
		k_inVpunt_tauTSHD_dummy(tel3)=kmax-kmaxTSHD-1
	      endif
	    in=.true.
	  endif
	enddo

!	add one Vpunt extra because staggered:
	if (tel2>0) then
		if (in.and.j_inVpuntTSHD_dummy(tel2).gt.0) then !always if a jet-point is found an extra Vpunt must be included
		  tel2=tel2+1
		  i_inVpuntTSHD_dummy(tel2)=i_inVpuntTSHD_dummy(tel2-1)
		  j_inVpuntTSHD_dummy(tel2)=j_inVpuntTSHD_dummy(tel2-1)-1
		  k_inVpuntTSHD_dummy(tel2)=k_inVpuntTSHD_dummy(tel2-1)
		endif
	endif
       enddo
      enddo
      tmax_inPpuntTSHD=tel1
      tmax_inVpuntTSHD=tel2
      tmax_inVpunt_tauTSHD=tel3
      tel1=0
      do t=1,tmax_inPpuntTSHD ! now search for all local j-indices between 0 and j1:
	if ((j_inPpuntTSHD_dummy(t)-rank*jmax).ge.0.and.(j_inPpuntTSHD_dummy(t)-rank*jmax).le.j1) then
		tel1=tel1+1
		j_inPpuntTSHD(tel1)=j_inPpuntTSHD_dummy(t)-rank*jmax
		i_inPpuntTSHD(tel1)=i_inPpuntTSHD_dummy(t)
		k_inPpuntTSHD(tel1)=k_inPpuntTSHD_dummy(t)
	endif
      enddo
      tmax_inPpuntTSHD=tel1
      tel1=0
      do t=1,tmax_inVpuntTSHD ! now search for all local j-indices between 0 and j1:
	if ((j_inVpuntTSHD_dummy(t)-rank*jmax).ge.0.and.(j_inVpuntTSHD_dummy(t)-rank*jmax).le.j1) then
		tel1=tel1+1
		j_inVpuntTSHD(tel1)=j_inVpuntTSHD_dummy(t)-rank*jmax
		i_inVpuntTSHD(tel1)=i_inVpuntTSHD_dummy(t)
		k_inVpuntTSHD(tel1)=k_inVpuntTSHD_dummy(t)
	endif
      enddo
      tmax_inVpuntTSHD=tel1
      tel1=0
      do t=1,tmax_inVpunt_tauTSHD ! now search for all local j-indices between 0 and j1:
	if ((j_inVpunt_tauTSHD_dummy(t)-rank*jmax).ge.0.and.(j_inVpunt_tauTSHD_dummy(t)-rank*jmax).le.j1) then
		tel1=tel1+1
		j_inVpunt_tauTSHD(tel1)=j_inVpunt_tauTSHD_dummy(t)-rank*jmax
		i_inVpunt_tauTSHD(tel1)=i_inVpunt_tauTSHD_dummy(t)
		k_inVpunt_tauTSHD(tel1)=k_inVpunt_tauTSHD_dummy(t)
	endif
      enddo
      tmax_inVpunt_tauTSHD=tel1

!	write(*,*) 'xTSHD',xTSHD
!	write(*,*) 'yTSHD',yTSHD

      !! Search for W:
      tel1=0
      do j=0,j1 
       do i=0,i1  
	in=.false.
        do k=k1,kmax-kmaxTSHD,-1   ! search on rho-loc, when rho-loc is in then also W down are in (made zero) 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
 	  zz=MIN(k*dz,depth) !k1 would lead to zz>depth
	  if (softnose.eq.1) then !make nose shorter gradually
		adj_ment=(1.+MAX(MIN(depth-DRAUGHT-zz,0.)/(Hfront),-1.)) !grows linear from 0. at Hfront to 1 at keel TSHD 
		adj_ment=1.-0.25*adj_ment !1 at Hfront and 0.75 at keel TSHD
		xTSHD2(1:20)=xnoseS*adj_ment+Lfront*(1.-adj_ment)+xfront ! make nose shorter
		xTSHD2(23:42)=xnoseP*adj_ment+Lfront*(1.-adj_ment)+xfront ! make nose shorter
		xTSHD=xTSHD2*cos(phi)-yTSHD2*sin(phi)
		yTSHD=xTSHD2*sin(phi)+yTSHD2*cos(phi)
	  endif
	  if (Hback>0.) then ! make shape shorter for slope at back:
	    xTSHD2(21:22)=LOA-(1.+MAX((depth-DRAUGHT-zz)/Hback,-1.))*Lback+xfront
	    xTSHD(21:22)=xTSHD2(21:22)*cos(phi)-yTSHD2(21:22)*sin(phi)
	    yTSHD(21:22)=xTSHD2(21:22)*sin(phi)+yTSHD2(21:22)*cos(phi)
	  else ! leave shape same length:
	    xTSHD2(21:22)=LOA+xfront
	    xTSHD(21:22)=xTSHD2(21:22)*cos(phi)-yTSHD2(21:22)*sin(phi)
	    yTSHD(21:22)=xTSHD2(21:22)*sin(phi)+yTSHD2(21:22)*cos(phi)
	  endif
	  CALL PNPOLY (xx,yy, xTSHD, yTSHD, 42, inout ) 
	  if (inout.eq.1) then
	      tel1=tel1+1
	      i_inWpuntTSHD(tel1)=i ! Wpunt
	      j_inWpuntTSHD(tel1)=j
  	      k_inWpuntTSHD(tel1)=k
	      in=.true.
	  endif
	enddo
!	add one Wpunt extra because staggered:
	if (tel1>0) then
		if (in.and.k_inWpuntTSHD(tel1)>0) then !always if a TSHD-point is found an extra Wpunt must included
		  tel1=tel1+1
		  i_inWpuntTSHD(tel1)=i_inWpuntTSHD(tel1-1)
		  j_inWpuntTSHD(tel1)=j_inWpuntTSHD(tel1-1)
		  k_inWpuntTSHD(tel1)=k_inWpuntTSHD(tel1-1)-1
		endif
	endif
       enddo
      enddo
      tmax_inWpuntTSHD=tel1

      !! Search for U: (search in i-dir)
      tel2=0
      tel3=0
      do k=kmax-kmaxTSHD,k1   
       do j=0,j1
	in=.false.
        do i=i1,0,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
 	  zz=MIN(k*dz,depth) !k1 would lead to zz>depth
	  if (softnose.eq.1) then !make nose shorter gradually
		adj_ment=(1.+MAX(MIN(depth-DRAUGHT-zz,0.)/(Hfront),-1.)) !grows linear from 0. at Hfront to 1 at keel TSHD 
		adj_ment=1.-0.25*adj_ment !1 at Hfront and 0.75 at keel TSHD
		xTSHD2(1:20)=xnoseS*adj_ment+Lfront*(1.-adj_ment)+xfront ! make nose shorter
		xTSHD2(23:42)=xnoseP*adj_ment+Lfront*(1.-adj_ment)+xfront ! make nose shorter
		xTSHD=xTSHD2*cos(phi)-yTSHD2*sin(phi)
		yTSHD=xTSHD2*sin(phi)+yTSHD2*cos(phi)
	  endif
	  if (Hback>0.) then ! make shape shorter for slope at back:
	    xTSHD2(21:22)=LOA-(1.+MAX((depth-DRAUGHT-zz)/Hback,-1.))*Lback+xfront
	    xTSHD(21:22)=xTSHD2(21:22)*cos(phi)-yTSHD2(21:22)*sin(phi)
	    yTSHD(21:22)=xTSHD2(21:22)*sin(phi)+yTSHD2(21:22)*cos(phi)
	  else ! leave shape same length:
	    xTSHD2(21:22)=LOA+xfront
	    xTSHD(21:22)=xTSHD2(21:22)*cos(phi)-yTSHD2(21:22)*sin(phi)
	    yTSHD(21:22)=xTSHD2(21:22)*sin(phi)+yTSHD2(21:22)*cos(phi)
	  endif
	  CALL PNPOLY (xx,yy, xTSHD, yTSHD, 42, inout ) 
	  if (inout.eq.1) then
	    	    tel2=tel2+1
		    i_inUpuntTSHD(tel2)=i  ! Upunt is U
		    j_inUpuntTSHD(tel2)=j
		    k_inUpuntTSHD(tel2)=k
	      if (k.eq.kmax-kmaxTSHD.and.(kmax-kmaxTSHD)>0.and.kn_TSHD.ge.0.0) then
		tel3=tel3+1
	        i_inUpunt_tauTSHD(tel3)=i ! Upunt-tau
	        j_inUpunt_tauTSHD(tel3)=j
		k_inUpunt_tauTSHD(tel3)=kmax-kmaxTSHD-1
	      endif
	    in=.true.
	  endif
	enddo

!	add one Upunt extra because staggered:
	if (tel2>0) then
		if (in.and.i_inUpuntTSHD(tel2)>0) then !always if a TSHD-point is found an extra Upunt must included
		  tel2=tel2+1
		  i_inUpuntTSHD(tel2)=i_inUpuntTSHD(tel2-1)-1
		  j_inUpuntTSHD(tel2)=j_inUpuntTSHD(tel2-1)
		  k_inUpuntTSHD(tel2)=k_inUpuntTSHD(tel2-1)
		endif
	endif
       enddo
      enddo
      tmax_inUpuntTSHD=tel2
      tmax_inUpunt_tauTSHD=tel3

	!! Determine U_TSHD split up in u,v along r,phi axes to calculate correct tau_wall
	!! Tau_wall only works for U_b,V_b
	Utshd=U_TSHD*cos(phi) ! in x,y coordinates
	Vtshd=U_TSHD*sin(phi) ! in x,y coordinates
	do j=0,j1
	  Ubot_TSHD(j)    =    Utshd*cos_u(j)+Vtshd*sin_u(j)   ! in r,phi coordinates
	  Vbot_TSHD(j)    =   -Utshd*sin_v(j)+Vtshd*cos_v(j)   ! in r,phi coordinates
	enddo

	xdh = xfront + Lfront + xdh !draghead starts at end of nose
	IF (draghead.eq.'star'.or.draghead.eq.'both') THEN
	! search for starboard draghead:

	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif

	ydh=0.5*Breadth+1.5*Dsp

      !! Search for P,V:

      tel1=0 !tmax_inPpuntTSHD
      tel2=0 !tmax_inVpuntTSHD
      do k=0,k1 
       do i=0,i1  
	in=.false.
        do j=jmax*px+1,0,-1         ! search in j dir, search on rho-loc, when rho-loc is in then also V left right are in (made zero) 
	  xx=Rp(i)*cos_ut(j)-schuif_x
	  yy=Rp(i)*sin_ut(j)
	  IF (k.le.FLOOR(Dsp/dz)) THEN ! draghead:
	xTSHD2(1)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(2)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(3)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(4)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	yTSHD2(1)=0.5*Breadth
	yTSHD2(2)=0.5*Breadth
	yTSHD2(3)=0.5*Breadth+3.*Dsp
	yTSHD2(4)=0.5*Breadth+3.*Dsp
	xTSHD=xTSHD2*cos(phi)-yTSHD2*sin(phi)
	yTSHD=xTSHD2*sin(phi)+yTSHD2*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE ! suction pipe:
	 	inout=0
				xdh2=xdh-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
		xdh_rot2=xdh2*cos(phi)-ydh*sin(phi)
		ydh_rot2=xdh2*sin(phi)+ydh*cos(phi)
		IF((xx-xdh_rot2)**2+(yy-ydh_rot2)**2.le.(0.5*Dsp)**2) THEN
			inout=1
		ENDIF
	  ENDIF
	  if (inout.eq.1) then
	      tel1=tel1+1
	      i_inPpuntTSHD_dummy(tel1)=i ! Ppunt
	      j_inPpuntTSHD_dummy(tel1)=j
  	      k_inPpuntTSHD_dummy(tel1)=k
      	      tel2=tel2+1
	      i_inVpuntTSHD_dummy(tel2)=i ! Vpunt
	      j_inVpuntTSHD_dummy(tel2)=j
	      k_inVpuntTSHD_dummy(tel2)=k
	    in=.true.
	  endif
	enddo
!	add one Vpunt extra because staggered:
	if (tel2>0) then
		if (in.and.j_inVpuntTSHD_dummy(tel2)>0) then !always if a jet-point is found an extra Vpunt must be included
		  tel2=tel2+1
		  i_inVpuntTSHD_dummy(tel2)=i_inVpuntTSHD_dummy(tel2-1)
		  j_inVpuntTSHD_dummy(tel2)=j_inVpuntTSHD_dummy(tel2-1)-1
		  k_inVpuntTSHD_dummy(tel2)=k_inVpuntTSHD_dummy(tel2-1)
		endif
	endif
       enddo
      enddo
!      tmax_inPpuntTSHD=tel1
!      tmax_inVpuntTSHD=tel2

      tel3=tmax_inPpuntTSHD
      do t=1,tel1 ! now search for all local j-indices between 0 and j1:
	if ((j_inPpuntTSHD_dummy(t)-rank*jmax).ge.0.and.(j_inPpuntTSHD_dummy(t)-rank*jmax).le.j1) then
		tel3=tel3+1
		j_inPpuntTSHD(tel3)=j_inPpuntTSHD_dummy(t)-rank*jmax
		i_inPpuntTSHD(tel3)=i_inPpuntTSHD_dummy(t)
		k_inPpuntTSHD(tel3)=k_inPpuntTSHD_dummy(t)
	endif
      enddo
      tmax_inPpuntTSHD=tel3	
      tel3=tmax_inVpuntTSHD
      do t=1,tel2 ! now search for all local j-indices between 0 and j1:
	if ((j_inVpuntTSHD_dummy(t)-rank*jmax).ge.0.and.(j_inVpuntTSHD_dummy(t)-rank*jmax).le.j1) then
		tel3=tel3+1
		j_inVpuntTSHD(tel3)=j_inVpuntTSHD_dummy(t)-rank*jmax
		i_inVpuntTSHD(tel3)=i_inVpuntTSHD_dummy(t)
		k_inVpuntTSHD(tel3)=k_inVpuntTSHD_dummy(t)
	endif
      enddo
      tmax_inVpuntTSHD=tel3


      !! Search for W:
      tel1=tmax_inWpuntTSHD
      do j=0,j1 
       do i=0,i1  
	in=.false.
        do k=0,k1 ! one extra in vertical direction for W  
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  IF (k.le.FLOOR(Dsp/dz)) THEN ! draghead:
	xTSHD2(1)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(2)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(3)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(4)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	yTSHD2(1)=0.5*Breadth
	yTSHD2(2)=0.5*Breadth
	yTSHD2(3)=0.5*Breadth+3.*Dsp
	yTSHD2(4)=0.5*Breadth+3.*Dsp
	xTSHD=xTSHD2*cos(phi)-yTSHD2*sin(phi)
	yTSHD=xTSHD2*sin(phi)+yTSHD2*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE ! suction pipe:
	 	inout=0
				xdh2=xdh-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
		xdh_rot2=xdh2*cos(phi)-ydh*sin(phi)
		ydh_rot2=xdh2*sin(phi)+ydh*cos(phi)
		IF((xx-xdh_rot2)**2+(yy-ydh_rot2)**2.le.(0.5*Dsp)**2) THEN
			inout=1
		ENDIF
	  ENDIF
	  if (inout.eq.1) then
	      tel1=tel1+1
	      i_inWpuntTSHD(tel1)=i ! Wpunt
	      j_inWpuntTSHD(tel1)=j
  	      k_inWpuntTSHD(tel1)=k
	      in=.true.
	  endif
	enddo
       enddo
      enddo
      tmax_inWpuntTSHD=tel1

      !! Search for U: (search in i-dir)
      tel2=tmax_inUpuntTSHD
      tel3=tmax_inWpunt_suction
      do k=0,k1
       do j=0,j1
	in=.false.
        do i=i1,0,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  IF (k.le.FLOOR(Dsp/dz)) THEN ! draghead:
	xTSHD2(1)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(2)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(3)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(4)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	yTSHD2(1)=0.5*Breadth
	yTSHD2(2)=0.5*Breadth
	yTSHD2(3)=0.5*Breadth+3.*Dsp
	yTSHD2(4)=0.5*Breadth+3.*Dsp
	xTSHD=xTSHD2*cos(phi)-yTSHD2*sin(phi)
	yTSHD=xTSHD2*sin(phi)+yTSHD2*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE ! suction pipe:
	 	inout=0
				xdh2=xdh-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
		xdh_rot2=xdh2*cos(phi)-ydh*sin(phi)
		ydh_rot2=xdh2*sin(phi)+ydh*cos(phi)
		IF((xx-xdh_rot2)**2+(yy-ydh_rot2)**2.le.(0.5*Dsp)**2) THEN
			inout=1
		ENDIF
	  ENDIF
	  if (inout.eq.1) then
	    	tel2=tel2+1
		i_inUpuntTSHD(tel2)=i  ! Upunt is U
		j_inUpuntTSHD(tel2)=j
		k_inUpuntTSHD(tel2)=k
		if (k.eq.0.and..NOT.in.and.i_inUpuntTSHD(tel2)<i1) then ! one row behind draghead should get suction velocity:
		  tel3=tel3+1
		  i_inWpunt_suction(tel3)=i_inUpuntTSHD(tel2)+1
		  j_inWpunt_suction(tel3)=j_inUpuntTSHD(tel2)
		  k_inWpunt_suction(tel3)=k_inUpuntTSHD(tel2)		
		endif
	        in=.true.


	  endif
	enddo
!	add one Upunt extra because staggered:
	if (tel2>0) then
		if (in.and.i_inUpuntTSHD(tel2)>0) then !always if a TSHD-point is found an extra Upunt must included
		  tel2=tel2+1
		  i_inUpuntTSHD(tel2)=i_inUpuntTSHD(tel2-1)-1
		  j_inUpuntTSHD(tel2)=j_inUpuntTSHD(tel2-1)
		  k_inUpuntTSHD(tel2)=k_inUpuntTSHD(tel2-1)
		endif
	endif
       enddo
      enddo
      tmax_inUpuntTSHD=tel2
      tmax_inWpunt_suction=tel3

	ENDIF !endif starboard draghead

	IF (draghead.eq.'port'.or.draghead.eq.'both') THEN
	! search for portside draghead:

	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif

	ydh=-0.5*Breadth-1.5*Dsp

      !! Search for P,V:
      tel1=0 !tmax_inPpuntTSHD
      tel2=0 !tmax_inVpuntTSHD
      do k=0,k1
       do i=0,i1  
	in=.false.
        do j=jmax*px+1,0,-1       ! global search in j dir, search on rho-loc, when rho-loc is in then also V left right are in (made zero) 
	  xx=Rp(i)*cos_ut(j)-schuif_x
	  yy=Rp(i)*sin_ut(j)
	  IF (k.le.FLOOR(Dsp/dz)) THEN ! draghead:
	xTSHD2(1)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(2)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(3)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(4)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	yTSHD2(1)=-0.5*Breadth-3.*Dsp
	yTSHD2(2)=-0.5*Breadth-3.*Dsp
	yTSHD2(3)=-0.5*Breadth
	yTSHD2(4)=-0.5*Breadth
	xTSHD=xTSHD2*cos(phi)-yTSHD2*sin(phi)
	yTSHD=xTSHD2*sin(phi)+yTSHD2*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE ! suction pipe:
	 	inout=0
		xdh2=xdh-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
		xdh_rot2=xdh2*cos(phi)-ydh*sin(phi)
		ydh_rot2=xdh2*sin(phi)+ydh*cos(phi)
		IF((xx-xdh_rot2)**2+(yy-ydh_rot2)**2.le.(0.5*Dsp)**2) THEN
			inout=1
		ENDIF
	  ENDIF
	  if (inout.eq.1) then
	      tel1=tel1+1
	      i_inPpuntTSHD_dummy(tel1)=i ! Ppunt
	      j_inPpuntTSHD_dummy(tel1)=j
  	      k_inPpuntTSHD_dummy(tel1)=k
      	      tel2=tel2+1
	      i_inVpuntTSHD_dummy(tel2)=i ! Vpunt
	      j_inVpuntTSHD_dummy(tel2)=j
	      k_inVpuntTSHD_dummy(tel2)=k
	    in=.true.
	  endif
	enddo

!	add one Vpunt extra because staggered:
	if (tel2>0) then
		if (in.and.j_inVpuntTSHD_dummy(tel2)>0) then !always if a jet-point is found an extra Vpunt must be included
		  tel2=tel2+1
		  i_inVpuntTSHD_dummy(tel2)=i_inVpuntTSHD_dummy(tel2-1)
		  j_inVpuntTSHD_dummy(tel2)=j_inVpuntTSHD_dummy(tel2-1)-1
		  k_inVpuntTSHD_dummy(tel2)=k_inVpuntTSHD_dummy(tel2-1)
		endif
	endif
       enddo
      enddo
!      tmax_inPpuntTSHD=tel1
!      tmax_inVpuntTSHD=tel2

      tel3=tmax_inPpuntTSHD
      do t=1,tel1 ! now search for all local j-indices between 0 and j1:
	if ((j_inPpuntTSHD_dummy(t)-rank*jmax).ge.0.and.(j_inPpuntTSHD_dummy(t)-rank*jmax).le.j1) then
		tel3=tel3+1
		j_inPpuntTSHD(tel3)=j_inPpuntTSHD_dummy(t)-rank*jmax
		i_inPpuntTSHD(tel3)=i_inPpuntTSHD_dummy(t)
		k_inPpuntTSHD(tel3)=k_inPpuntTSHD_dummy(t)
	endif
      enddo
      tmax_inPpuntTSHD=tel3	
      tel3=tmax_inVpuntTSHD
      do t=1,tel2 ! now search for all local j-indices between 0 and j1:
	if ((j_inVpuntTSHD_dummy(t)-rank*jmax).ge.0.and.(j_inVpuntTSHD_dummy(t)-rank*jmax).le.j1) then
		tel3=tel3+1
		j_inVpuntTSHD(tel3)=j_inVpuntTSHD_dummy(t)-rank*jmax
		i_inVpuntTSHD(tel3)=i_inVpuntTSHD_dummy(t)
		k_inVpuntTSHD(tel3)=k_inVpuntTSHD_dummy(t)
	endif
      enddo
      tmax_inVpuntTSHD=tel3


      !! Search for W:
      tel1=tmax_inWpuntTSHD
      do j=0,j1 
       do i=0,i1  
	in=.false.
        do k=0,k1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  IF (k.le.FLOOR(Dsp/dz)) THEN ! draghead:
	xTSHD2(1)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(2)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(3)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(4)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	yTSHD2(1)=-0.5*Breadth-3.*Dsp
	yTSHD2(2)=-0.5*Breadth-3.*Dsp
	yTSHD2(3)=-0.5*Breadth
	yTSHD2(4)=-0.5*Breadth
	xTSHD=xTSHD2*cos(phi)-yTSHD2*sin(phi)
	yTSHD=xTSHD2*sin(phi)+yTSHD2*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE ! suction pipe:
	 	inout=0
				xdh2=xdh-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
		xdh_rot2=xdh2*cos(phi)-ydh*sin(phi)
		ydh_rot2=xdh2*sin(phi)+ydh*cos(phi)
		IF((xx-xdh_rot2)**2+(yy-ydh_rot2)**2.le.(0.5*Dsp)**2) THEN
			inout=1
		ENDIF
	  ENDIF
	  if (inout.eq.1) then
	      tel1=tel1+1
	      i_inWpuntTSHD(tel1)=i ! Wpunt
	      j_inWpuntTSHD(tel1)=j
  	      k_inWpuntTSHD(tel1)=k
	      in=.true.
	  endif
	enddo
       enddo
      enddo
      tmax_inWpuntTSHD=tel1

      !! Search for U: (search in i-dir)
      tel2=tmax_inUpuntTSHD
      tel3=tmax_inWpunt_suction
      do k=0,k1
       do j=0,j1
	in=.false.
        do i=i1,0,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  IF (k.le.FLOOR(Dsp/dz)) THEN ! draghead:
	xTSHD2(1)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(2)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(3)=xdh+0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	xTSHD2(4)=xdh-0.5*Dsp-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
	yTSHD2(1)=-0.5*Breadth-3.*Dsp
	yTSHD2(2)=-0.5*Breadth-3.*Dsp
	yTSHD2(3)=-0.5*Breadth
	yTSHD2(4)=-0.5*Breadth
	xTSHD=xTSHD2*cos(phi)-yTSHD2*sin(phi)
	yTSHD=xTSHD2*sin(phi)+yTSHD2*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE ! suction pipe:
	 	inout=0
				xdh2=xdh-REAL(k)/REAL(k1)*(xdh-Lfront-xfront)
		xdh_rot2=xdh2*cos(phi)-ydh*sin(phi)
		ydh_rot2=xdh2*sin(phi)+ydh*cos(phi)
		IF((xx-xdh_rot2)**2+(yy-ydh_rot2)**2.le.(0.5*Dsp)**2) THEN
			inout=1
		ENDIF
	  ENDIF
	  if (inout.eq.1) then
	    	tel2=tel2+1
		i_inUpuntTSHD(tel2)=i  ! Upunt is U
		j_inUpuntTSHD(tel2)=j
		k_inUpuntTSHD(tel2)=k
		if (k.eq.0.and..NOT.in.and.i_inUpuntTSHD(tel2)<i1) then ! one row behind draghead should get suction velocity:
		  tel3=tel3+1
		  i_inWpunt_suction(tel3)=i_inUpuntTSHD(tel2)+1
		  j_inWpunt_suction(tel3)=j_inUpuntTSHD(tel2)
		  k_inWpunt_suction(tel3)=k_inUpuntTSHD(tel2)		
		endif
                in=.true.
	  endif
	enddo
!	add one Upunt extra because staggered:
	if (tel2>0) then
		if (in.and.i_inUpuntTSHD(tel2)>0) then !always if a TSHD-point is found an extra Upunt must included
		  tel2=tel2+1
		  i_inUpuntTSHD(tel2)=i_inUpuntTSHD(tel2-1)-1
		  j_inUpuntTSHD(tel2)=j_inUpuntTSHD(tel2-1)
		  k_inUpuntTSHD(tel2)=k_inUpuntTSHD(tel2-1)
		endif
	endif
       enddo
      enddo
      tmax_inUpuntTSHD=tel2
      tmax_inWpunt_suction=tel3 

!	do i=1,tmax_inWpunt_suction
!	  if(j_inWpunt_suction(i)>0.and.j_inWpunt_suction(i)<j1) then
!	    write(*,*),'rank,i,j,k',rank,i_inWpunt_suction(i),j_inWpunt_suction(i)+rank*jmax,k_inWpunt_suction(i)
!	  endif
!	enddo


      ENDIF !endif portside draghead


      ENDIF !endif LOA>0



	!! Search for obstacles near bed:  
	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
	kbed=0
	zbed=0.

	DO n=1,nobst
	!! search for zbed and kbed on each proc (j=0,j1): (used for determination of obstacle)
	      do j=0,j1 
	       do i=0,i1  
		in=.false.
		do k=0,k1 ! one extra in vertical direction for W  
		  xx=Rp(i)*cos_u(j)-schuif_x
		  yy=Rp(i)*sin_u(j)
		  IF (k.le.FLOOR(ob(n)%height/dz).and.ob(n)%zbottom.le.0.) THEN ! obstacle part of bed:
			xTSHD(1:4)=ob(n)%x*cos(phi)-ob(n)%y*sin(phi)
			yTSHD(1:4)=ob(n)%x*sin(phi)+ob(n)%y*cos(phi)
			CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		  ELSE 
		 	inout=0
		  ENDIF
		  if (inout.eq.1) then
		      kbed(i,j)=MAX(kbed(i,j),FLOOR(ob(n)%height/dz)) !zero without obstacle, otherwise max of all obstacles at i,j
			  kbed(i,j)=MIN(kbed(i,j),kmax)
		      zbed(i,j)=MAX(zbed(i,j),ob(n)%height) !zero without obstacle, otherwise max of all obstacles at i,j  
			  IF (interaction_bed.ge.4) THEN
			    bednotfixed(i,j,k)=ob(n)%ero ! obstacle in bed cannot be avalanched or eroded if ob(n)%ero=0. (default) user can choose to have erosion: ob(n)%ero=1.
			    bednotfixed_depo(i,j,k)=ob(n)%depo ! obstacle in bed cannot have deposition if ob(n)%depo=0. (default)
			  ENDIF
		  endif
		enddo
	       enddo
	      enddo
	ENDDO

	kbed2=0
	DO n=1,nobst
	!! search for zbed and kbed for dummy array (j=0:px*jmax+1):
	      do j=0,jmax*px+1 
	       do i=0,i1  
		in=.false.
		do k=0,k1 ! one extra in vertical direction for W  
		  xx=Rp(i)*cos_ut(j)-schuif_x
		  yy=Rp(i)*sin_ut(j)
		  IF (k.le.FLOOR(ob(n)%height/dz).and.ob(n)%zbottom.le.0.) THEN ! obstacle part of bed:
			xTSHD(1:4)=ob(n)%x*cos(phi)-ob(n)%y*sin(phi)
			yTSHD(1:4)=ob(n)%x*sin(phi)+ob(n)%y*cos(phi)
			CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
		  ELSE 
		 	inout=0
		  ENDIF
		  if (inout.eq.1) then
		      kbed2(i,j)=MAX(kbed2(i,j),FLOOR(ob(n)%height/dz)) !zero without obstacle, otherwise max of all obstacles at i,j
		  endif
		enddo
	       enddo
	      enddo
	ENDDO

	DO n=1,nobst
      !! Search for P,V:
      tel1=0 !tmax_inPpuntTSHD
      tel2=0 !tmax_inVpuntTSHD
      do k=0,k1
       do i=0,i1  
	in=.false.
        do j=jmax*px+1,0,-1       ! global search in j dir, search on rho-loc, when rho-loc is in then also V left right are in (made zero) 
	  xx=Rp(i)*cos_ut(j)-schuif_x
	  yy=Rp(i)*sin_ut(j)

!	  IF (k.le.FLOOR(ob(n)%height/dz)) THEN ! obstacle:
	  !IF ((k.ge.CEILING(ob(n)%zbottom/dz)).and.(k.le.FLOOR(ob(n)%height/dz)).and.(FLOOR(ob(n)%height/dz).eq.kbed2(i,j))) THEN ! obstacle:
	  IF ((k.ge.CEILING(ob(n)%zbottom/dz)).and.(k.le.FLOOR(ob(n)%height/dz)).and.ob(n)%zbottom.gt.0.) THEN ! obstacle: adjustment 10-1-2017 because some obstacles are missed with check on kbed
		xTSHD(1:4)=ob(n)%x*cos(phi)-ob(n)%y*sin(phi)
		yTSHD(1:4)=ob(n)%x*sin(phi)+ob(n)%y*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
	      tel1=tel1+1
	      i_inPpuntTSHD_dummy(tel1)=i ! Ppunt
	      j_inPpuntTSHD_dummy(tel1)=j
  	      k_inPpuntTSHD_dummy(tel1)=k
      	      tel2=tel2+1
	      i_inVpuntTSHD_dummy(tel2)=i ! Vpunt
	      j_inVpuntTSHD_dummy(tel2)=j
	      k_inVpuntTSHD_dummy(tel2)=k
	    in=.true.
	  endif
	enddo

!	add one Vpunt extra because staggered:
	if (tel2>0) then
		if (in.and.j_inVpuntTSHD_dummy(tel2)>0) then !always if a jet-point is found an extra Vpunt must be included
		  tel2=tel2+1
		  i_inVpuntTSHD_dummy(tel2)=i_inVpuntTSHD_dummy(tel2-1)
		  j_inVpuntTSHD_dummy(tel2)=j_inVpuntTSHD_dummy(tel2-1)-1
		  k_inVpuntTSHD_dummy(tel2)=k_inVpuntTSHD_dummy(tel2-1)
		endif
	endif
       enddo
      enddo
!      tmax_inPpuntTSHD=tel1
!      tmax_inVpuntTSHD=tel2

	
      tel3=tmax_inPpuntTSHD
      do t=1,tel1 ! now search for all local j-indices between 0 and j1:
	if ((j_inPpuntTSHD_dummy(t)-rank*jmax).ge.0.and.(j_inPpuntTSHD_dummy(t)-rank*jmax).le.j1) then
		tel3=tel3+1
		j_inPpuntTSHD(tel3)=j_inPpuntTSHD_dummy(t)-rank*jmax
		i_inPpuntTSHD(tel3)=i_inPpuntTSHD_dummy(t)
		k_inPpuntTSHD(tel3)=k_inPpuntTSHD_dummy(t)
	endif
      enddo
      tmax_inPpuntTSHD=tel3	
      tel3=tmax_inVpuntTSHD
      do t=1,tel2 ! now search for all local j-indices between 0 and j1:
	if ((j_inVpuntTSHD_dummy(t)-rank*jmax).ge.0.and.(j_inVpuntTSHD_dummy(t)-rank*jmax).le.j1) then
		tel3=tel3+1
		j_inVpuntTSHD(tel3)=j_inVpuntTSHD_dummy(t)-rank*jmax
		i_inVpuntTSHD(tel3)=i_inVpuntTSHD_dummy(t)
		k_inVpuntTSHD(tel3)=k_inVpuntTSHD_dummy(t)
	endif
      enddo
      tmax_inVpuntTSHD=tel3

      !! Search for W:
      tel1=tmax_inWpuntTSHD
      do j=0,j1 
       do i=0,i1  
	in=.false.
        do k=0,k1 ! one extra in vertical direction for W  
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(ob(n)%height/dz)) THEN ! obstacle:
	  !IF ((k.ge.CEILING(ob(n)%zbottom/dz)).and.(k.le.FLOOR(ob(n)%height/dz)).and.(FLOOR(ob(n)%height/dz).eq.kbed(i,j))) THEN ! obstacle:
	  IF ((k.ge.CEILING(ob(n)%zbottom/dz)).and.(k.le.FLOOR(ob(n)%height/dz)).and.ob(n)%zbottom.gt.0.) THEN ! obstacle: adjustment 10-1-2017 because some obstacles are missed with check on kbed
		xTSHD(1:4)=ob(n)%x*cos(phi)-ob(n)%y*sin(phi)
		yTSHD(1:4)=ob(n)%x*sin(phi)+ob(n)%y*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
	      tel1=tel1+1
	      i_inWpuntTSHD(tel1)=i ! Wpunt
	      j_inWpuntTSHD(tel1)=j
  	      k_inWpuntTSHD(tel1)=k
	      in=.true.
		  IF (interaction_bed.ge.4) THEN
			  bednotfixed(i,j,k)=ob(n)%ero ! obstacle in bed cannot be avalanched or eroded if ob(n)%ero=0. (default) user can choose to have erosion: ob(n)%ero=1.
			  bednotfixed_depo(i,j,k)=ob(n)%depo ! obstacle in bed cannot have deposition if ob(n)%depo=0. (default)		  
		  ENDIF
		  
!	      kbed(i,j)=MAX(kbed(i,j),FLOOR(ob(n)%height/dz)) !zero without obstacle, otherwise max of all obstacles at i,j
!	      zbed(i,j)=MAX(zbed(i,j),ob(n)%height) !zero without obstacle, otherwise max of all obstacles at i,j  
	  endif
	enddo
       enddo
      enddo
      tmax_inWpuntTSHD=tel1

      !! Search for U: (search in i-dir)
      tel2=tmax_inUpuntTSHD
      do k=0,k1
       do j=0,j1
	in=.false.
        do i=i1,0,-1 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
!	  IF (k.le.FLOOR(ob(n)%height/dz)) THEN ! obstacle:
!	  IF ((k.ge.CEILING(ob(n)%zbottom/dz)).and.(k.le.FLOOR(ob(n)%height/dz)).and.(FLOOR(ob(n)%height/dz).eq.kbed(i,j))) THEN ! obstacle:
	  IF ((k.ge.CEILING(ob(n)%zbottom/dz)).and.(k.le.FLOOR(ob(n)%height/dz)).and.ob(n)%zbottom.gt.0.) THEN ! obstacle: adjustment 10-1-2017 because some obstacles are missed with check on kbed
		xTSHD(1:4)=ob(n)%x*cos(phi)-ob(n)%y*sin(phi)
		yTSHD(1:4)=ob(n)%x*sin(phi)+ob(n)%y*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
	    	tel2=tel2+1
		i_inUpuntTSHD(tel2)=i  ! Upunt is U
		j_inUpuntTSHD(tel2)=j
		k_inUpuntTSHD(tel2)=k
                in=.true.
	  endif
	enddo
!	add one Upunt extra because staggered:
	if (tel2>0) then
		if (in.and.i_inUpuntTSHD(tel2)>0) then !always if a TSHD-point is found an extra Upunt must included
		  tel2=tel2+1
		  i_inUpuntTSHD(tel2)=i_inUpuntTSHD(tel2-1)-1
		  j_inUpuntTSHD(tel2)=j_inUpuntTSHD(tel2-1)
		  k_inUpuntTSHD(tel2)=k_inUpuntTSHD(tel2-1)
		endif
	endif
       enddo
      enddo
      tmax_inUpuntTSHD=tel2
	ENDDO ! obstacles loop

	
	!! waarom nog een keer kbed bepalen??
!	kbed=0
!	zbed=0.
!	DO n=1,nobst
!	!! search for zbed and kbed on each proc (j=0,j1): (used for deposition and bc in solver [adjusted for ob(n)%zbottom])
!	      do j=0,j1 
!	       do i=0,i1  
!		in=.false.
!		do k=0,k1 
!		  xx=Rp(i)*cos_u(j)-schuif_x
!		  yy=Rp(i)*sin_u(j)
!		  IF (k.le.FLOOR(ob(n)%height/dz)) THEN ! obstacle:
!			xTSHD(1:4)=ob(n)%x*cos(phi)-ob(n)%y*sin(phi)
!			yTSHD(1:4)=ob(n)%x*sin(phi)+ob(n)%y*cos(phi)
!			CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
!		  ELSE 
!		 	inout=0
!		  ENDIF
!		  if (inout.eq.1.and.ob(n)%zbottom.le.0.) then ! only when zbottom.le.0, then kbed and zbed should deviate from zero, zbottom default<0
!		      kbed(i,j)=MAX(kbed(i,j),FLOOR(ob(n)%height/dz)) !zero without obstacle, otherwise max of all obstacles at i,j
!		      zbed(i,j)=MAX(zbed(i,j),ob(n)%height) !zero without obstacle, otherwise max of all obstacles at i,j  
!		  endif
!		enddo
!	       enddo
!	      enddo
!	ENDDO

	
	if (.true.) then ! true = 0th order IBM staircase manner, false = 1st order IBM scaling with volume grid cell inside object
	if (bedlevelfile.ne.'') then

       	status2 = nf90_open(bedlevelfile, nf90_NoWrite, ncid) 
	IF (status2/= nf90_noerr) THEN
		write(*,*),'bedlevelfile =',bedlevelfile
		CALL writeerror(606)
	ENDIF
	call check( nf90_inq_varid(ncid, "zbed",rhVarid) )
	call check( nf90_get_var(ncid,rhVarid,zbed3(1:imax,0:j1),start=(/1,rank*jmax+1/),count=(/imax,jmax+2/)) )
	call check( nf90_close(ncid) )
	
	! bedlevel file obstacles have no i=0 or i=i1 in zbed3, but do have j=0 and j=j1
	kbed3=0
	zbed3(0,0:j1)=zbed3(1,0:j1)
	zbed3(i1,0:j1)=zbed3(imax,0:j1)
	!! search for zbed and kbed on each proc (j=0,j1): (used for deposition and bc in solver [adjusted for ob(n)%zbottom])
	      do j=0,j1 
	       do i=0,i1 !imax  !including i1 strangely gives crash (13/4/15) !1,imax !0,i1
				kbed3(i,j)=FLOOR(zbed3(i,j)/dz+0.5) !added+0.5 10-3-2020 because now top kbed level is within +0.5 and -0.5dz from real bed (cnewbot can be positive and negative) 
				kbed3(i,j)=MAX(0,kbed3(i,j))
				kbed3(i,j)=MIN(kmax,kbed3(i,j))
				zbed(i,j)=MAX(zbed(i,j),zbed3(i,j)) !zero without obstacle, otherwise max of all obstacles at i,j  
				
		enddo
	       enddo
		
		
		IF(.false.) THEN ! 11-2-2017 adjusted; bedlevelfile is now immersed boundary via kbed not via TSHD arrays
      !! Search for P,V:
      tel1=tmax_inPpuntTSHD
      tel2=tmax_inVpuntTSHD
      do k=0,k1
       do i=0,i1  
	in=.false.
        do j=j1,0,-1       
	  IF (k.le.kbed3(i,j).and.kbed3(i,j).gt.kbed(i,j)) THEN ! new obstacle:
		inout=1
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
	      tel1=tel1+1
	      i_inPpuntTSHD(tel1)=i ! Ppunt
	      j_inPpuntTSHD(tel1)=j
  	      k_inPpuntTSHD(tel1)=k
      	      tel2=tel2+1
	      i_inVpuntTSHD(tel2)=i ! Vpunt
	      j_inVpuntTSHD(tel2)=j
	      k_inVpuntTSHD(tel2)=k
	    in=.true.
	  endif
	enddo
!	add one Vpunt extra because staggered:
	if (tel2>0) then
		if (in.and.j_inVpuntTSHD_dummy(tel2)>0) then !always if a jet-point is found an extra Vpunt must be included
		  tel2=tel2+1
		  i_inVpuntTSHD(tel2)=i_inVpuntTSHD(tel2-1)
		  j_inVpuntTSHD(tel2)=j_inVpuntTSHD(tel2-1)-1
		  k_inVpuntTSHD(tel2)=k_inVpuntTSHD(tel2-1)
		endif
	endif
       enddo
      enddo
      tmax_inPpuntTSHD=tel1
      tmax_inVpuntTSHD=tel2

      !! Search for W:
      tel1=tmax_inWpuntTSHD
      do j=0,j1 
       do i=0,i1  
	in=.false.
        do k=0,k1 ! one extra in vertical direction for W  
	  IF (k.le.kbed3(i,j).and.kbed3(i,j).gt.kbed(i,j)) THEN ! new obstacle:
		inout=1
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
	      tel1=tel1+1
	      i_inWpuntTSHD(tel1)=i ! Wpunt
	      j_inWpuntTSHD(tel1)=j
  	      k_inWpuntTSHD(tel1)=k
	      in=.true.
!	      kbed(i,j)=MAX(kbed(i,j),FLOOR(ob(n)%height/dz)) !zero without obstacle, otherwise max of all obstacles at i,j
!	      zbed(i,j)=MAX(zbed(i,j),ob(n)%height) !zero without obstacle, otherwise max of all obstacles at i,j  
	  endif
	enddo
       enddo
      enddo
      tmax_inWpuntTSHD=tel1

      !! Search for U: (search in i-dir)
      tel2=tmax_inUpuntTSHD
      do k=0,k1
       do j=0,j1
	in=.false.
        do i=i1,0,-1 
	  IF (k.le.kbed3(i,j).and.kbed3(i,j).gt.kbed(i,j)) THEN ! new obstacle:
		inout=1
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
	    	tel2=tel2+1
		i_inUpuntTSHD(tel2)=i  ! Upunt is U
		j_inUpuntTSHD(tel2)=j
		k_inUpuntTSHD(tel2)=k
                in=.true.
	  endif
	enddo
!	add one Upunt extra because staggered:
	if (tel2>0) then
		if (in.and.i_inUpuntTSHD(tel2)>0) then !always if a TSHD-point is found an extra Upunt must included
		  tel2=tel2+1
		  i_inUpuntTSHD(tel2)=i_inUpuntTSHD(tel2-1)-1
		  j_inUpuntTSHD(tel2)=j_inUpuntTSHD(tel2-1)
		  k_inUpuntTSHD(tel2)=k_inUpuntTSHD(tel2-1)
		endif
	endif
       enddo
      enddo
      tmax_inUpuntTSHD=tel2
	ENDIF ! endif (.false.) ! 11-2-2017 adjusted; bedlevelfile is now immersed boundary via kbed not via TSHD arrays
		!! update kbed on each proc (j=0,j1) with kbed3 (used for deposition and bc in solver [adjusted for ob(n)%zbottom])
		do j=0,j1 
			do i=0,i1 !imax !including i1 strangely gives crash (13/4/15) !1,imax !0,i1
				kbed(i,j)=MAX(kbed(i,j),kbed3(i,j)) !zero without obstacle, otherwise max of all obstacles at i,j  
			enddo
		enddo
	endif
	else ! 1st order IBM with facIBM scaling with percentage volume of grid cell captured in immersed object
!	if (bedlevelfile.ne.'') then
!
!       	status2 = nf90_open(bedlevelfile, nf90_NoWrite, ncid) 
!	IF (status2/= nf90_noerr) THEN
!		write(*,*),'bedlevelfile =',bedlevelfile
!		CALL writeerror(606)
!	ENDIF
!	call check( nf90_inq_varid(ncid, "zbed",rhVarid) )
!	call check( nf90_get_var(ncid,rhVarid,zbed3(1:imax,0:j1),start=(/1,rank*jmax+1/),count=(/imax,jmax+2/)) )
!	call check( nf90_close(ncid) )
!	
!	!write(*,*),'rank,zbed3(1,0:j1)',rank,zbed3(1,0:j1)
!	
!	! bedlevel file obstacles have no i=0 or i=i1 in zbed3, but do have j=0 and j=j1
!	kbed3=0
!	!! search for zbed and kbed on each proc (j=0,j1): (used for deposition and bc in solver [adjusted for ob(n)%zbottom])
!	      do j=0,j1 
!	       do i=0,imax  !including i1 strangely gives crash (13/4/15) !1,imax !0,i1
!				kbed3(i,j)=FLOOR(zbed3(i,j)/dz)
!				kbed3(i,j)=MAX(0,kbed3(i,j))
!				kbed3(i,j)=MIN(kmax,kbed3(i,j))
!				zbed(i,j)=MAX(zbed(i,j),zbed3(i,j)) !zero without obstacle, otherwise max of all obstacles at i,j  
!		enddo
!	       enddo
!		facIBMu=0. !initialise all facIBM zero because for all ibm cells velocity must be zero
!		facIBMv=0.
!		facIBMw=0.
!!c get stuff from other CPU's
!
!	zbed4(0:i1,0:j1)=zbed(0:i1,0:j1)	  
!	!call shiftf_l2(zbed4,cbf) 
!	call shiftb_l2(zbed4(0:i1,0:j1),cbb) 
!	if (periodicy.eq.0.or.periodicy.eq.2) then
!		if (rank.eq.0) then ! boundaries in j-direction
!		   do i=1,imax
!		   !zbed4(i,-1) = zbed4(i,1) !cbf(i,k)
!		   zbed4(i,j1+1) =cbb(i) 
!		   enddo
!		elseif (rank.eq.px-1) then
!		   do i=1,imax
!		   !zbed4(i,-1) = cbf(i)
!		   zbed4(i,j1+1) =zbed4(i,jmax) !cbb(i,k) 
!		   enddo
!		else
!		   do i=1,imax
!		   !zbed4(i,-1) = cbf(i)
!		   zbed4(i,j1+1) =cbb(i) 
!		   enddo
!		endif
!	else
!	   do i=1,imax
!		   !zbed4(i,-1) = cbf(i)
!		   zbed4(i,j1+1) =cbb(i) 
!	   enddo
!	endif
!	if (periodicx.eq.1) then
!		zbed4(0,0:j1)=zbed4(imax,0:j1)
!		zbed4(i1,0:j1)=zbed4(1,0:j1)
!		zbed4(i1+1,0:j1)=zbed4(2,0:j1)
!	else
!		zbed4(0,0:j1)=zbed4(1,0:j1)
!		zbed4(i1,0:j1)=zbed4(imax,0:j1)
!		zbed4(i1+1,0:j1)=zbed4(i1,0:j1)
!	endif
!	!write(*,*),'rank,zbed4(0:i1+1,0:j1)',rank,zbed4(0:i1+1,0:j1)
!	
!		
!      !! Search for P:
!      tel1=tmax_inPpuntTSHD
!      do k=0,k1
!       do i=0,i1  
!	in=.false.
!        do j=j1,0,-1       
!	  IF (k.le.kbed3(i,j).and.kbed3(i,j).gt.kbed(i,j)) THEN ! new obstacle:
!		inout=1
!	  ELSE 
!	 	inout=0
!	  ENDIF
!	  if (inout.eq.1) then
!	      tel1=tel1+1
!	      i_inPpuntTSHD(tel1)=i ! Ppunt
!	      j_inPpuntTSHD(tel1)=j
!  	      k_inPpuntTSHD(tel1)=k
!	    in=.true.
!	  endif
!	enddo
!       enddo
!      enddo	
!	
!      !! Search for V:
!      tel2=tmax_inVpuntTSHD
!      do k=0,k1
!       do i=0,i1  
!	in=.false.
!        do j=j1,0,-1       
!	  IF (k-3.le.kbed3(i,j).and.kbed3(i,j).gt.kbed(i,j)) THEN ! new obstacle:
!		inout=1
!	  ELSE 
!	 	inout=0
!	  ENDIF
!	  if (inout.eq.1) then
!      	      tel2=tel2+1
!	      i_inVpuntTSHD(tel2)=i ! Vpunt
!	      j_inVpuntTSHD(tel2)=j
!	      k_inVpuntTSHD(tel2)=k
!		  zbotcell=(k-1)*dz+0.5*dz
!		  zbedcell=0.5*(zbed4(i,j)+zbed4(i,j+1))
!		  facIBMv(tel2)=max(min((zbotcell-zbedcell)/dz,1.),0.)
!	    in=.true.
!	  endif
!	enddo
!
!	
!!	add one Vpunt extra because staggered: !not needed because each individual point is being checked with 1st or 2nd order ibm
!!	if (in.and.j_inVpuntTSHD_dummy(tel2)>0) then !always if a jet-point is found an extra Vpunt must be included
!!	  tel2=tel2+1
!!	  i_inVpuntTSHD(tel2)=i_inVpuntTSHD(tel2-1)
!!	  j_inVpuntTSHD(tel2)=j_inVpuntTSHD(tel2-1)-1
!!	  k_inVpuntTSHD(tel2)=k_inVpuntTSHD(tel2-1)
!!		  zbotcell=(k-1)*dz+0.5*dz
!!		  zbedcell=0.5*(zbed4(i,j-1)+zbed4(i,j))
!!		  facIBMv(tel2)=max(min((zbotcell-zbedcell)/dz,1.),0.)
!!	endif
!       enddo
!      enddo
!      tmax_inPpuntTSHD=tel1
!      tmax_inVpuntTSHD=tel2
!
!      !! Search for W:
!      tel1=tmax_inWpuntTSHD
!      do j=0,j1 
!       do i=0,i1  
!	in=.false.
!        do k=0,k1 ! one extra in vertical direction for W  
!	  IF (k-3.le.kbed3(i,j).and.kbed3(i,j).gt.kbed(i,j)) THEN ! new obstacle:
!		inout=1
!	  ELSE 
!	 	inout=0
!	  ENDIF
!	  if (inout.eq.1) then
!	      tel1=tel1+1
!	      i_inWpuntTSHD(tel1)=i ! Wpunt
!	      j_inWpuntTSHD(tel1)=j
!  	      k_inWpuntTSHD(tel1)=k
!		  zbotcell=k*dz
!		  zbedcell=zbed4(i,j)
!		  facIBMw(tel1)=max(min((zbotcell-zbedcell)/dz,1.),0.)
!		  
!		!  write(*,*),'i,j,k,zbed4(i,j),zbotcell,facIBMw(tel1)',i,j,k,zbed4(i,j),zbotcell,facIBMw(tel1)
!		  
!	      in=.true.
!!	      kbed(i,j)=MAX(kbed(i,j),FLOOR(ob(n)%height/dz)) !zero without obstacle, otherwise max of all obstacles at i,j
!!	      zbed(i,j)=MAX(zbed(i,j),ob(n)%height) !zero without obstacle, otherwise max of all obstacles at i,j  
!	  endif
!	enddo
!       enddo
!      enddo
!      tmax_inWpuntTSHD=tel1
!
!      !! Search for U: (search in i-dir)
!      tel2=tmax_inUpuntTSHD
!      do k=0,k1
!       do j=0,j1
!	in=.false.
!        do i=i1,0,-1 
!	  IF (k-3.le.kbed3(i,j).and.kbed3(i,j).gt.kbed(i,j)) THEN ! new obstacle:
!		inout=1
!	  ELSE 
!	 	inout=0
!		
!	  ENDIF
!	  if (inout.eq.1) then
!	    	tel2=tel2+1
!		i_inUpuntTSHD(tel2)=i  ! Upunt is U
!		j_inUpuntTSHD(tel2)=j
!		k_inUpuntTSHD(tel2)=k
!		  zbotcell=(k-1)*dz+0.5*dz
!		  zbedcell=0.5*(zbed4(i,j)+zbed4(i+1,j))
!		  facIBMu(tel2)=max(min((zbotcell-zbedcell)/dz,1.),0.)		
!                in=.true.
!			!	write(*,*),'i,j,k,zbedcell,zbotcell,facIBMu(tel2)',i,j,k,zbedcell,zbotcell,facIBMu(tel2)
!	  endif
!	enddo
!!	add one Upunt extra because staggered: ! not needed with 1st or 2nd order ibm, then each individual point is checked
!!	if (in.and.i_inUpuntTSHD(tel2)>0) then !always if a TSHD-point is found an extra Upunt must included
!!	  tel2=tel2+1
!!	  i_inUpuntTSHD(tel2)=i_inUpuntTSHD(tel2-1)-1
!!	  j_inUpuntTSHD(tel2)=j_inUpuntTSHD(tel2-1)
!!	  k_inUpuntTSHD(tel2)=k_inUpuntTSHD(tel2-1)
!!		  zbotcell=(k-1)*dz+0.5*dz
!!		  zbedcell=0.5*(zbed4(i_inUpuntTSHD(tel2), j_inUpuntTSHD(tel2))+zbed4( i_inUpuntTSHD(tel2)+1, j_inUpuntTSHD(tel2)))
!!		  facIBMu(tel2)=max(min((zbotcell-zbedcell)/dz,1.),0.)			  
!!		  write(*,*),'i,j,k,zbedcell,zbotcell,facIBMu(tel2)', i_inUpuntTSHD(tel2), j_inUpuntTSHD(tel2),k,zbedcell,zbotcell,facIBMu(tel2)
!!	endif
!       enddo
!      enddo
!      tmax_inUpuntTSHD=tel2
!	  
!	
!		!! update kbed on each proc (j=0,j1) with kbed3 (used for deposition and bc in solver [adjusted for ob(n)%zbottom])
!	      do j=0,j1 
!	       do i=0,imax !including i1 strangely gives crash (13/4/15) !1,imax !0,i1
!	          kbed(i,j)=MAX(kbed(i,j),kbed3(i,j)) !zero without obstacle, otherwise max of all obstacles at i,j  
!		enddo
!	       enddo
!	endif
	endif !endif 1st order ibm	
	
	      do j=0,j1 
	       do i=0,imax 
	          kbed(i,j)=MIN(kbed(i,j),kmax) ! make sure kbed never is larger than kmax
	       enddo
	      enddo


		IF (interaction_bed.eq.4.or.interaction_bed.eq.6) THEN   
			do j=0,j1 
				do i=0,i1 !imax !including i1 strangely gives crash (13/4/15) !1,imax !0,i1
					do k=1,kbed(i,j) ! assign initial bed concentrations; k=0 remains empty
						do n=1,nfrac
							Clivebed(n,i,j,k)=c_bed(n)*bednotfixed(i,j,k) ! in case obstacle no erosion (bednotfixed=0) then no cbed
							!Cold(n,i,j,k)=c_bed(n)
							!Cnew(n,i,j,k)=c_bed(n)
							!dCdt(n,i,j,k)=c_bed(n)	
						enddo
					enddo
					do n=1,nfrac
					  ! place remainder which does not fit exactly in 1dz into cnewbot; this can be positive or negative as cnewbot varies between +/-0.5*cfixedbed 
					  Cnewbot(n,i,j)=(zbed(i,j)-DBLE(kbed(i,j))*dz)/dz*c_bed(n)*bednotfixed(i,j,kbed(i,j))
					enddo 
				enddo
			enddo
			Coldbot = Cnewbot 
			dCdtbot = Cnewbot
		ENDIF	

		IF (U_TSHD>0.) THEN
			kbedin(0:j1)=kbed(1,0:j1)
		ENDIF
		
		  IF (IBMorder.ne.2) THEN 
		    kbed22=kbed 
		  ELSE 
		  	do j=1,jmax 
				do i=1,imax 
					IF (interaction_bed.ge.4) THEN 
					  zb_W=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
					ELSE 
					  zb_W=zbed(i,j)
					ENDIF
					kbed22(i,j)=FLOOR(zb_W/dz) 
				ENDDO
			ENDDO 
			call bound_cbot_integer(kbed22) 
		  ENDIF 
      end

	subroutine bedroughness_init

        USE nlist
	implicit none

	real phi,z,ust_U_b,ust_V_b,z0_U,z0_V,correction
	integer ii


	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif

	!! velocity U_b is determined for (depth-bc_obst_h): 
	!! when local depth differs the Ubc is adjusted such that total flux remains equal !13-4-15 correction switched off !

	!! fill Ubc1,Vbc1 (lateral boundaries at j=0 or j=jmax):
	    IF (rank.eq.0) THEN
		j=0
	    ELSE ! also for rank between 0 and px-1 now Ubc1 is filled, but this is not used anyways
		j=jmax 
	    ENDIF	
	    do k=1,kmax
		 do i=1,imax
		  if (zbed(i,j).lt.depth) then
	  	    correction=1. !13-4-15 correction switched off !(depth-bc_obst_h)/(depth-zbed(i,j))
	          else !if obstacle persists through full watercolumn then U_b should be zero at this location
		    correction=0.
		  endif
		  if (LOA>0.) then !ship:
		   ust_U_b=MAX(ABS(U_b),1.e-6)
		   ust_V_b=MAX(ABS(V_b),1.e-6)
		   if (slip_bot.eq.1) then
			    do ii=1,10
			      z0_U=0.11*nu_mol/MAX(ust_U_b,1.e-9)+kn/30
			      ust_U_b=correction*ABS(U_b)*kappa/(log(MAX(depth-zbed(i,j),1.e-12)/z0_U)-1);
			      z0_V=0.11*nu_mol/MAX(ust_V_b,1.e-9)+kn/30
			      ust_V_b=correction*ABS(V_b)*kappa/(log(MAX(depth-zbed(i,j),1.e-12)/z0_V)-1);
			    enddo
		   endif
		    if (k.gt.kbed(i,j)) then
		      if (slip_bot.eq.1) then
				if (wallup.eq.1) then
					  z=depth-((k-kbed(i,j))*dz-0.5*dz)
				else
				  z=(k-kbed(i,j))*dz-0.5*dz
				endif
		    	Ubc1(i,k)=-ust_U_b/kappa*log(z/z0_U)*signU_b*cos(phi)+ust_V_b/kappa*log(z/z0_V)*signV_b*sin(phi)
     &                   +U_TSHD*cos(phi)*correction
		    	Vbc1(i,k)=-ust_V_b/kappa*log(z/z0_V)*signV_b*cos(phi)-ust_U_b/kappa*log(z/z0_U)*signU_b*sin(phi)
     &                   +U_TSHD*sin(phi)*correction
	 	      else
				if (k.gt.ksurf_bc) then
				  Ubc1(i,k)=(-U_b3*cos(phi)+V_b3*sin(phi)+U_TSHD*cos(phi))*correction
				  Vbc1(i,k)=(-V_b3*cos(phi)-U_b3*sin(phi)+U_TSHD*sin(phi))*correction
				else
				  Ubc1(i,k)=sqrt((U_TSHD-U_b)**2+V_b**2)*correction
				  Vbc1(i,k)=0.
				endif
		      endif
		    else
			  Ubc1(i,k)=Ubot_TSHD(j) !0.
			  Vbc1(i,k)=Vbot_TSHD(j) !0.
		    endif
		  else ! no ship, then plate:
		    ust_U_b=MAX(ABS(U_b),1.e-6)
		    ust_V_b=MAX(ABS(V_b),1.e-6)
		    if (slip_bot.eq.1) then
			    do ii=1,10
			      z0_U=0.11*nu_mol/MAX(ABS(ust_U_b),1.e-9)+kn/30
			      ust_U_b=correction*U_b*kappa/(log(MAX(depth-zbed(i,j),1.e-12)/z0_U)-1);
			      z0_V=0.11*nu_mol/MAX(ABS(ust_V_b),1.e-9)+kn/30
			      ust_V_b=correction*V_b*kappa/(log(MAX(depth-zbed(i,j),1.e-12)/z0_V)-1);
			    enddo
		    endif
		      if (k.le.(kmax-kjet).and.k.gt.kbed(i,j)) then
			    if (slip_bot.eq.1) then
				  if (wallup.eq.1) then
					z=depth-((k-kbed(i,j))*dz-0.5*dz)
				  else
					z=(k-kbed(i,j))*dz-0.5*dz
				  endif
			      Ubc1(i,k)=ust_U_b/kappa*log(z/z0_U) 
			      Vbc1(i,k)=ust_V_b/kappa*log(z/z0_V)
		 	    else
				  Ubc1(i,k)=U_b*correction 
				  Vbc1(i,k)=0.
			    endif
		      else
			    Ubc1(i,k)=0.
			    Vbc1(i,k)=0.
		      endif
		  endif
		 enddo
	    enddo 

	!! fill Ubc2,Vbc2 (front boundary at i=0):
		if (monopile>0) then !inflow Ubc2,Vbc2 defined at imax instead at i=0
		  i=imax 
		else 
	      i=0
		endif 
	    do k=1,kmax
		 do j=0,j1
		  if (zbed(i,j).lt.depth) then
	  	    correction=1. !13-4-15 correction switched off !(depth-bc_obst_h)/(depth-zbed(i,j))
	          else !if obstacle persists through full watercolumn then U_b should be zero at this location
		    correction=0.
		  endif
		  if (LOA>0.) then !ship:
		   ust_U_b=MAX(ABS(U_b),1.e-6)
		   ust_V_b=MAX(ABS(V_b),1.e-6)
		   if (slip_bot.eq.1) then
		    do ii=1,10
		      z0_U=0.11*nu_mol/MAX(ust_U_b,1.e-9)+kn/30
		      ust_U_b=correction*ABS(U_b)*kappa/(log(MAX(depth-zbed(i,j),1.e-12)/z0_U)-1);
		      z0_V=0.11*nu_mol/MAX(ust_V_b,1.e-9)+kn/30
		      ust_V_b=correction*ABS(V_b)*kappa/(log(MAX(depth-zbed(i,j),1.e-12)/z0_V)-1);
		    enddo
		   endif
		    if (k.gt.kbed(i,j)) then
		     if (slip_bot.eq.1) then
				if (wallup.eq.1) then
				  z=depth-((k-kbed(i,j))*dz-0.5*dz)
				else
				  z=(k-kbed(i,j))*dz-0.5*dz
				endif
		    	Ubc2(j,k)=-ust_U_b/kappa*log(z/z0_U)*signU_b*cos(phi)+ust_V_b/kappa*log(z/z0_V)*signV_b*sin(phi)
     &                   +U_TSHD*cos(phi)*correction
		    	Vbc2(j,k)=-ust_V_b/kappa*log(z/z0_V)*signV_b*cos(phi)-ust_U_b/kappa*log(z/z0_U)*signU_b*sin(phi)
     &                   +U_TSHD*sin(phi)*correction
	 	     else
				if (k.gt.ksurf_bc) then
				  Ubc2(j,k)=(-U_b3*cos(phi)+V_b3*sin(phi)+U_TSHD*cos(phi))*correction
				  Vbc2(j,k)=(-V_b3*cos(phi)-U_b3*sin(phi)+U_TSHD*sin(phi))*correction
				else
				  Ubc2(j,k)=sqrt((U_TSHD-U_b)**2+V_b**2)*correction
				  Vbc2(j,k)=0.
				endif
		     endif
		    else
			  Ubc2(j,k)=Ubot_TSHD(j) !0.
			  Vbc2(j,k)=Vbot_TSHD(j) !0.
		    endif
		  else ! no ship, then plate:
		   ust_U_b=MAX(ABS(U_b),1.e-6)
		   ust_V_b=MAX(ABS(V_b),1.e-6)
		   if (slip_bot.eq.1) then
		    do ii=1,10
		      z0_U=0.11*nu_mol/MAX(ABS(ust_U_b),1.e-9)+kn/30
		      ust_U_b=correction*U_b*kappa/(log(MAX(depth-zbed(i,j),1.e-12)/z0_U)-1);
		      z0_V=0.11*nu_mol/MAX(ABS(ust_V_b),1.e-9)+kn/30
		      ust_V_b=correction*V_b*kappa/(log(MAX(depth-zbed(i,j),1.e-12)/z0_V)-1);
		    enddo
		   endif
		      if (k.le.(kmax-kjet).and.k.gt.kbed(i,j)) then
			    if (slip_bot.eq.1) then
                        	if (wallup.eq.1) then
                          	  z=depth-((k-kbed(i,j))*dz-0.5*dz)
                        	else
                          	  z=(k-kbed(i,j))*dz-0.5*dz
                                endif
			    	Ubc2(j,k)=ust_U_b/kappa*log(z/z0_U) 
			    	Vbc2(j,k)=ust_V_b/kappa*log(z/z0_V)
		 	    else
				Ubc2(j,k)=U_b*correction 
				Vbc2(j,k)=0.
			    endif
		      else
			Ubc2(j,k)=0.
			Vbc2(j,k)=0.
		      endif
		  endif
		 enddo
	    enddo

	
	
	call bound_cbot_integer(kbed)
	! called as last therefore now kbedt (used to apply tau wall) can be defined:
	IF (wallup.eq.1) THEN
	  kbedt=kmax-1 ! apply tau wall at upper boundary kmax
	ELSE
	  kbedt=kbed ! apply tau wall at lower boundary kbed
	ENDIF
	

	end 



C .................................................................. 
C 
C SUBROUTINE PNPOLY 
C 
C PURPOSE 
C TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON 
C 
C USAGE 
C CALL PNPOLY (PX, PY, XX, YY, N, INOUT ) 
C 
C DESCRIPTION OF THE PARAMETERS 
C PX - X-COORDINATE OF POINT IN QUESTION. 
C PY - Y-COORDINATE OF POINT IN QUESTION. 
C XX - N LONG VECTOR CONTAINING X-COORDINATES OF 
C VERTICES OF POLYGON. 
C YY - N LONG VECTOR CONTAING Y-COORDINATES OF 
C VERTICES OF POLYGON. 
C N - NUMBER OF VERTICES IN THE POLYGON. 
C INOUT - THE SIGNAL RETURNED: 
C -1 IF THE POINT IS OUTSIDE OF THE POLYGON, 
C 0 IF THE POINT IS ON AN EDGE OR AT A VERTEX, 
C 1 IF THE POINT IS INSIDE OF THE POLYGON. 
C 
C REMARKS 
C THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE. 
C THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY 
C OPTIONALLY BE INCREASED BY 1. 
C THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING 
C OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX 
C OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING 
C N, THESE FIRST VERTICES MUST BE COUNTED TWICE. 
C INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED. 
C THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM 
C WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70. 
C 
C SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED 
C NONE 
C 
C METHOD 
C A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT 
C CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE 
C POINT IS INSIDE OF THE POLYGON. 
C 
C .................................................................. 
C 
	SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT) 
	implicit none
	REAL X(200),Y(200),XX(N),YY(N),PX,PY
	LOGICAL MX,MY,NX,NY 
	INTEGER O,INOUT,I,N,J,MAXDIM
C OUTPUT UNIT FOR PRINTED MESSAGES 
	DATA O/6/ 
	MAXDIM=200 

	IF(N.LE.MAXDIM)GO TO 6 
	WRITE(O,7) 
    7 FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY. 
     & RESULTS INVALID') 
	RETURN 
    6 DO 1 I=1,N 
	X(I)=XX(I)-PX 
    1 Y(I)=YY(I)-PY 
	INOUT=-1 
	DO 2 I=1,N 
	  J=1+MOD(I,N) 
	  MX=X(I).GE.0.0 
	  NX=X(J).GE.0.0 
	  MY=Y(I).GE.0.0 
	  NY=Y(J).GE.0.0 
	  IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2 
	    IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3 
	      INOUT=-INOUT 
	      GO TO 2 
    3 IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5 
    4 INOUT=0 
	RETURN 
    5 INOUT=-INOUT 
    2 CONTINUE 
	RETURN 
	END 




