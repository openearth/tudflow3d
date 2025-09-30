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

	subroutine Kinetic_equation(ib,ie,jb,je,kb,ke)
	
	USE nlist
	implicit none
	!include 'mpif.h'
	
	integer ib,ie,jb,je,kb,ke
	integer*8 kbed_dummy(0:i1,0:j1)
	real lambda_advec(0:i1,0:j1,0:k1)
	
	lambda_advec(:,:,:)=0.
	kbed_dummy(:,:)=0

	call advecc_VLE(lambda_advec(:,:,:),lambda_old(:,:,:),Unew,Vnew,Wnew,rnew,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy,'volufrac',kbed_dummy) !advection of lambda

	do i=1,imax
		do j=1,jmax
			do k=1,kmax
				lambda_new(i,j,k)=(lambda_old(i,j,k)+dt*(Kin_eq_a*Kin_eq_lambda_0+(lambda_advec(i,j,k))))/
     &							  (1.+dt*(Kin_eq_a+Kin_eq_b*strain(i,j,k)))
			enddo
		enddo
	enddo
	call bound_3D(lambda_new) !Neumann B.C.
		
	end subroutine Kinetic_equation
	
	subroutine Houska_Papanastasiou(ib,ie,jb,je,kb,ke)
	
	USE nlist
	implicit none
	!include 'mpif.h'
	
	integer ib,ie,jb,je,kb,ke

	call strain_magnitude(Unew,Vnew,Wnew)
	call Kinetic_equation(ib,ie,jb,je,kb,ke)
	
	do i=1,imax						!start loop, in r-direction
		do j=1,jmax						!start loop, in phi-direction
			do k=1,kmax					!start loop, in z-direction
				
				muA(i,j,k)= ((HOUSKA_tauy_inf+lambda_old(i,j,k)*(HOUSKA_tauy_0-HOUSKA_tauy_inf))/(1.e-12+strain(i,j,k)))*
     &						(1.-exp(-PAPANASTASIOUS_m*(1.e-12+strain(i,j,k))))+
     &						(HOUSKA_eta_inf+(lambda_old(i,j,k)*(HOUSKA_eta_0-HOUSKA_eta_inf)))*(strain(i,j,k)**(HOUSKA_n-1.))
				ekm(i,j,k)= ekm(i,j,k) + muA(i,j,k)
			enddo
		enddo
	enddo
	
	lambda_old(:,:,:)=lambda_new(:,:,:)

	end subroutine Houska_Papanastasiou

	subroutine strain_magnitude(Uvel,Vvel,Wvel)
			
	USE nlist
	implicit none
	!include 'mpif.h'
	
	real dzi
	real dRpp_i,dRp_i
	real Uvel(0:i1,0:j1,0:k1),Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1)
	real Uvel2(0:i1,0:j1,0:k1),Vvel2(0:i1,0:j1,0:k1),Wvel2(0:i1,0:j1,0:k1)	 
	real dudx,dudy,dudz,strain2(0:i1,0:j1,0:k1)
	real dvdx,dvdy,dvdz
	real dwdx,dwdy,dwdz
	real S11,S12,S13,S22,S23,S33,SijSij,shear,divergentie
	integer n,im,ip,jm,jp,km,kp

	SijSij=0.

	IF (rheo_shear_method.eq.2) THEN !determine shear on rho*U and divide by rho at the end 
	  Uvel = Uvel * rhU !before bound_rhoU so rhU belongs to Unew 
	  Vvel = Vvel * rhV !before bound_rhoU so rhU belongs to Unew 
	  Wvel = Wvel * rhW !before bound_rhoU so rhU belongs to Unew 
	  dzi=1./dz
	  do i=1,imax
		dRpp_i=1./(Rp(i+1)-Rp(i))
		dRp_i=1./(Rp(i)-Rp(i-1))
		do j=1,jmax
			do k=1,kmax
				
				shear = 2.0*(
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)) +	
     1			((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi) )
     
     
				shear = shear + 0.25*(
	1			((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))/(Ru(i)*(phip(j+1)-phip(j))) +
     2			( Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))/(Ru(i)*(phip(j)-phip(j-1))) +
     2			( Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))/(Ru(i-1)*(phip(j+1)-phip(j))) +
     2			( Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2 +
	1			((Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))/(Ru(i-1)*(phip(j)-phip(j-1))) +
     2			( Vvel(i  ,j-1,k  )/Rp(i)-Vvel(i-1,j-1,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2
     e				)

				shear = shear + 0.25*(
	1			((Uvel(i  ,j  ,k+1)-Uvel(i  ,j  ,k  ))*dzi +
     2			( Wvel(i+1,j  ,k  )-Wvel(i  ,j  ,k  ))*dRpp_i)**2 +
	1			((Uvel(i  ,j  ,k  )-Uvel(i  ,j  ,k-1))*dzi +
     2			( Wvel(i+1,j  ,k-1)-Wvel(i  ,j  ,k-1))*dRpp_i)**2 +
	1			((Uvel(i-1,j  ,k+1)-Uvel(i-1,j  ,k  ))*dzi +
     2			( Wvel(i  ,j  ,k  )-Wvel(i-1,j  ,k  ))*dRp_i)**2 +
	1			((Uvel(i-1,j  ,k  )-Uvel(i-1,j  ,k-1))*dzi +
     2			( Wvel(i  ,j  ,k-1)-Wvel(i-1,j  ,k-1))*dRp_i)**2 
     e				)
	 
				shear = shear + 0.25*(
	1			((Vvel(i  ,j  ,k+1)-Vvel(i  ,j  ,k  ))*dzi +
     2			( Wvel(i  ,j+1,k  )-Wvel(i  ,j  ,k  ))/(Rp(i)*(phip(j+1)-phip(j))) )**2 +
	1			((Vvel(i  ,j  ,k  )-Vvel(i  ,j  ,k-1))*dzi +
     2			( Wvel(i  ,j+1,k-1)-Wvel(i  ,j  ,k-1))/(Rp(i)*(phip(j+1)-phip(j))) )**2 +
	1			((Vvel(i  ,j-1,k+1)-Vvel(i  ,j-1,k  ))*dzi +
     2			( Wvel(i  ,j  ,k  )-Wvel(i  ,j-1,k  ))/(Rp(i)*(phip(j)-phip(j-1))) )**2 +
	1			((Vvel(i  ,j-1,k  )-Vvel(i  ,j-1,k-1))*dzi +
     2			( Wvel(i  ,j  ,k-1)-Wvel(i  ,j-1,k-1))/(Rp(i)*(phip(j)-phip(j-1))) )**2 
     e				)

!! 
      divergentie= 2./3.*(( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) ) / ( dz ))   
       
		shear = shear 
     +             +1.5*divergentie*divergentie
     +             -divergentie*2.*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))
     +             -divergentie*2.*((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1)))+0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))
     +             -divergentie*2.*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)

				strain(i,j,k) = Apvisc_shear_relax*sqrt(shear)/rnew(i,j,k) + (1.-Apvisc_shear_relax)*strain(i,j,k)
			enddo
		enddo
	  enddo
	  Uvel = Uvel / rhU !before bound_rhoU so rhU belongs to Unew 
	  Vvel = Vvel / rhV !before bound_rhoU so rhU belongs to Unew 
	  Wvel = Wvel / rhW !before bound_rhoU so rhU belongs to Unew 	  
	ELSEIF (rheo_shear_method.eq.3) THEN ! determine shear with horizontal velocity components only 
	  dzi=1./dz
	  do i=1,imax
		dRpp_i=1./(Rp(i+1)-Rp(i))
		dRp_i=1./(Rp(i)-Rp(i-1))
		do j=1,jmax
			do k=1,kmax
				
				shear = 2.0*(
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)))
     
     
				shear = shear + 0.25*(
	1			((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))/(Ru(i)*(phip(j+1)-phip(j))) +
     2			( Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))/(Ru(i)*(phip(j)-phip(j-1))) +
     2			( Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))/(Ru(i-1)*(phip(j+1)-phip(j))) +
     2			( Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2 +
	1			((Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))/(Ru(i-1)*(phip(j)-phip(j-1))) +
     2			( Vvel(i  ,j-1,k  )/Rp(i)-Vvel(i-1,j-1,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2
     e				)

				strain(i,j,k) = Apvisc_shear_relax*sqrt(shear) + (1.-Apvisc_shear_relax)*strain(i,j,k)
			enddo
		enddo
	  enddo
	ELSEIF (rheo_shear_method.eq.4) THEN ! determine shear with horizontal velocity components only and including rho*U
	  Uvel = Uvel * rhU !before bound_rhoU so rhU belongs to Unew 
	  Vvel = Vvel * rhV !before bound_rhoU so rhU belongs to Unew 
	 ! Wvel = Wvel * rhW !before bound_rhoU so rhU belongs to Unew 
	  dzi=1./dz
	  do i=1,imax
		dRpp_i=1./(Rp(i+1)-Rp(i))
		dRp_i=1./(Rp(i)-Rp(i-1))
		do j=1,jmax
			do k=1,kmax
				
				shear = 2.0*(
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)))
     
     
				shear = shear + 0.25*(
	1			((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))/(Ru(i)*(phip(j+1)-phip(j))) +
     2			( Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))/(Ru(i)*(phip(j)-phip(j-1))) +
     2			( Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))/(Ru(i-1)*(phip(j+1)-phip(j))) +
     2			( Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2 +
	1			((Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))/(Ru(i-1)*(phip(j)-phip(j-1))) +
     2			( Vvel(i  ,j-1,k  )/Rp(i)-Vvel(i-1,j-1,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2
     e				)

				strain(i,j,k) = Apvisc_shear_relax*sqrt(shear)/rnew(i,j,k) + (1.-Apvisc_shear_relax)*strain(i,j,k)
			enddo
		enddo
	  enddo
	  Uvel = Uvel / rhU !before bound_rhoU so rhU belongs to Unew 
	  Vvel = Vvel / rhV !before bound_rhoU so rhU belongs to Unew 
!	  Wvel = Wvel / rhW !before bound_rhoU so rhU belongs to Unew 
	ELSEIF (rheo_shear_method.eq.5) THEN !determine shear on only dudx dvdy dwdz 
	  dzi=1./dz
	  do i=1,imax
		dRpp_i=1./(Rp(i+1)-Rp(i))
		dRp_i=1./(Rp(i)-Rp(i-1))
		do j=1,jmax
			do k=1,kmax
				
				shear = 2.0*(
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)) +	
     1			((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi) )
				strain(i,j,k) = Apvisc_shear_relax*sqrt(shear) + (1.-Apvisc_shear_relax)*strain(i,j,k)
			enddo
		enddo
	  enddo
	ELSEIF (rheo_shear_method.eq.6) THEN !determine shear on rho*U and divide by rho at the end + only dudx dvdy dwdz 
	  Uvel = Uvel * rhU !before bound_rhoU so rhU belongs to Unew 
	  Vvel = Vvel * rhV !before bound_rhoU so rhU belongs to Unew 
	  Wvel = Wvel * rhW !before bound_rhoU so rhU belongs to Unew 
	  dzi=1./dz
	  do i=1,imax
		dRpp_i=1./(Rp(i+1)-Rp(i))
		dRp_i=1./(Rp(i)-Rp(i-1))
		do j=1,jmax
			do k=1,kmax
				
				shear = 2.0*(
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)) +	
     1			((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi) )
				strain(i,j,k) = Apvisc_shear_relax*sqrt(shear)/rnew(i,j,k) + (1.-Apvisc_shear_relax)*strain(i,j,k)
			enddo
		enddo
	  enddo
	  Uvel = Uvel / rhU !before bound_rhoU so rhU belongs to Unew 
	  Vvel = Vvel / rhV !before bound_rhoU so rhU belongs to Unew 
	  Wvel = Wvel / rhW !before bound_rhoU so rhU belongs to Unew 
	ELSEIF (rheo_shear_method.eq.7) THEN !determine shear on only dudx dvdy
	  dzi=1./dz
	  do i=1,imax
		dRpp_i=1./(Rp(i+1)-Rp(i))
		dRp_i=1./(Rp(i)-Rp(i-1))
		do j=1,jmax
			do k=1,kmax
				
				shear = 2.0*(
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))) 	

				strain(i,j,k) = Apvisc_shear_relax*sqrt(shear) + (1.-Apvisc_shear_relax)*strain(i,j,k)
			enddo
		enddo
	  enddo
	ELSEIF (rheo_shear_method.eq.8) THEN !determine shear on rho*U and divide by rho at the end + only dudx dvdy dwdz 
	  Uvel = Uvel * rhU !before bound_rhoU so rhU belongs to Unew 
	  Vvel = Vvel * rhV !before bound_rhoU so rhU belongs to Unew 
	  dzi=1./dz
	  do i=1,imax
		dRpp_i=1./(Rp(i+1)-Rp(i))
		dRp_i=1./(Rp(i)-Rp(i-1))
		do j=1,jmax
			do k=1,kmax
				
				shear = 2.0*(
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)))
				strain(i,j,k) = Apvisc_shear_relax*sqrt(shear)/rnew(i,j,k) + (1.-Apvisc_shear_relax)*strain(i,j,k)
			enddo
		enddo
	  enddo
	  Uvel = Uvel / rhU !before bound_rhoU so rhU belongs to Unew 
	  Vvel = Vvel / rhV !before bound_rhoU so rhU belongs to Unew 
	ELSEIF (rheo_shear_method.eq.12) THEN !determine shear on Ctot*U and divide by Ctot at the end 
	  
	  Uvel2 = Uvel * SUM(cU(1:nfrac,:,:,:),DIM=1) !before bound_rhoU so cU belongs to Unew 
	  Vvel2 = Vvel * SUM(cV(1:nfrac,:,:,:),DIM=1) !before bound_rhoU so cU belongs to Unew 
	  Wvel2 = Wvel * SUM(cW(1:nfrac,:,:,:),DIM=1) !before bound_rhoU so cU belongs to Unew 
	  dzi=1./dz
	  do i=1,imax
		dRpp_i=1./(Rp(i+1)-Rp(i))
		dRp_i=1./(Rp(i)-Rp(i-1))
		do j=1,jmax
			do k=1,kmax
				
				shear = 2.0*(
     1			((Uvel2(i,j,k)-Uvel2(i-1,j,k))/dr(i))*((Uvel2(i,j,k)-Uvel2(i-1,j,k))/dr(i)) + 
     1			((Vvel2(i,j,k)-Vvel2(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel2(i,j,k)+Uvel2(i-1,j,k))/Rp(i))*
     1			((Vvel2(i,j,k)-Vvel2(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel2(i,j,k)+Uvel2(i-1,j,k))/Rp(i)) +	
     1			((Wvel2(i,j,k)-Wvel2(i,j,k-1))*dzi)*((Wvel2(i,j,k)-Wvel2(i,j,k-1))*dzi) )
     
     
				shear = shear + 0.25*(
	1			((Uvel2(i  ,j+1,k  )-Uvel2(i  ,j  ,k  ))/(Ru(i)*(phip(j+1)-phip(j))) +
     2			( Vvel2(i+1,j  ,k  )/Rp(i+1)-Vvel2(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel2(i  ,j  ,k  )-Uvel2(i  ,j-1,k  ))/(Ru(i)*(phip(j)-phip(j-1))) +
     2			( Vvel2(i+1,j-1,k  )/Rp(i+1)-Vvel2(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel2(i-1,j+1,k  )-Uvel2(i-1,j  ,k  ))/(Ru(i-1)*(phip(j+1)-phip(j))) +
     2			( Vvel2(i  ,j  ,k  )/Rp(i)-Vvel2(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2 +
	1			((Uvel2(i-1,j  ,k  )-Uvel2(i-1,j-1,k  ))/(Ru(i-1)*(phip(j)-phip(j-1))) +
     2			( Vvel2(i  ,j-1,k  )/Rp(i)-Vvel2(i-1,j-1,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2
     e				)

				shear = shear + 0.25*(
	1			((Uvel2(i  ,j  ,k+1)-Uvel2(i  ,j  ,k  ))*dzi +
     2			( Wvel2(i+1,j  ,k  )-Wvel2(i  ,j  ,k  ))*dRpp_i)**2 +
	1			((Uvel2(i  ,j  ,k  )-Uvel2(i  ,j  ,k-1))*dzi +
     2			( Wvel2(i+1,j  ,k-1)-Wvel2(i  ,j  ,k-1))*dRpp_i)**2 +
	1			((Uvel2(i-1,j  ,k+1)-Uvel2(i-1,j  ,k  ))*dzi +
     2			( Wvel2(i  ,j  ,k  )-Wvel2(i-1,j  ,k  ))*dRp_i)**2 +
	1			((Uvel2(i-1,j  ,k  )-Uvel2(i-1,j  ,k-1))*dzi +
     2			( Wvel2(i  ,j  ,k-1)-Wvel2(i-1,j  ,k-1))*dRp_i)**2 
     e				)
	 
				shear = shear + 0.25*(
	1			((Vvel2(i  ,j  ,k+1)-Vvel2(i  ,j  ,k  ))*dzi +
     2			( Wvel2(i  ,j+1,k  )-Wvel2(i  ,j  ,k  ))/(Rp(i)*(phip(j+1)-phip(j))) )**2 +
	1			((Vvel2(i  ,j  ,k  )-Vvel2(i  ,j  ,k-1))*dzi +
     2			( Wvel2(i  ,j+1,k-1)-Wvel2(i  ,j  ,k-1))/(Rp(i)*(phip(j+1)-phip(j))) )**2 +
	1			((Vvel2(i  ,j-1,k+1)-Vvel2(i  ,j-1,k  ))*dzi +
     2			( Wvel2(i  ,j  ,k  )-Wvel2(i  ,j-1,k  ))/(Rp(i)*(phip(j)-phip(j-1))) )**2 +
	1			((Vvel2(i  ,j-1,k  )-Vvel2(i  ,j-1,k-1))*dzi +
     2			( Wvel2(i  ,j  ,k-1)-Wvel2(i  ,j-1,k-1))/(Rp(i)*(phip(j)-phip(j-1))) )**2 
     e				)

!! 
      divergentie= 2./3.*(( Ru(i)*Uvel2(i,j,k) - Ru(i-1)*Uvel2(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel2(i,j,k) -         Vvel2(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       Wvel2(i,j,k) -         Wvel2(i,j,k-1) ) / ( dz ))   
       
		shear = shear 
     +       +1.5*divergentie*divergentie
     +       -divergentie*2.*((Uvel2(i,j,k)-Uvel2(i-1,j,k))/dr(i))
     +       -divergentie*2.*((Vvel2(i,j,k)-Vvel2(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1)))+0.5*(Uvel2(i,j,k)+Uvel2(i-1,j,k))/Rp(i))
     +       -divergentie*2.*((Wvel2(i,j,k)-Wvel2(i,j,k-1))*dzi)

	    strain(i,j,k) = Apvisc_shear_relax*sqrt(shear)/(SUM(cnew(1:nfrac,i,j,k))+1.e-12) + (1.-Apvisc_shear_relax)*strain(i,j,k)
			enddo
		enddo
	  enddo
	ELSE 
	  dzi=1./dz
	  do i=1,imax
		dRpp_i=1./(Rp(i+1)-Rp(i))
		dRp_i=1./(Rp(i)-Rp(i-1))
		do j=1,jmax
			do k=1,kmax
!				dudx = (Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i) 																	!du/dx
!				dvdy = (Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)		!dv/dy
!				dwdz = (Wvel(i,j,k)-Wvel(i,j,k-1))*dzi																		!dw/dz
!
!				dudy = ((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))/(Ru(i)*(phip(j+1)-phip(j))) + 
!     &		        	(Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))/(Ru(i)*(phip(j)-phip(j-1))) +
!     &		        	(Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))/(Ru(i-1)*(phip(j+1)-phip(j))) +
!     &		        	(Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))/(Ru(i-1)*(phip(j)-phip(j-1))) ) * 0.25
!				dudz = ((Uvel(i  ,j  ,k+1)-Uvel(i  ,j  ,k  ))*dzi +
!     &		        	(Uvel(i  ,j  ,k  )-Uvel(i  ,j  ,k-1))*dzi +
!     &		        	(Uvel(i-1,j  ,k+1)-Uvel(i-1,j  ,k  ))*dzi +
!     &		        	(Uvel(i-1,j  ,k  )-Uvel(i-1,j  ,k-1))*dzi ) * 0.25
!				dvdx = ((Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) +
!     &		       		(Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) + 
!     &		       		(Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1)  +
!     &		       		(Vvel(i  ,j-1,k  )/Rp(i)-Vvel(i-1,j-1,k  )/Rp(i-1))*dRp_i*Ru(i-1) ) * 0.25
!				dvdz = ((Vvel(i  ,j  ,k+1)-Vvel(i  ,j  ,k  ))*dzi +
!     &		        	(Vvel(i  ,j  ,k  )-Vvel(i  ,j  ,k-1))*dzi +
!     &		        	(Vvel(i  ,j-1,k+1)-Vvel(i  ,j-1,k  ))*dzi +
!     &		        	(Vvel(i  ,j-1,k  )-Vvel(i  ,j-1,k-1))*dzi ) * 0.25
!				dwdx = ((Wvel(i+1,j  ,k  )-Wvel(i  ,j  ,k  ))*dRpp_i + 
!     &		       		(Wvel(i+1,j  ,k-1)-Wvel(i  ,j  ,k-1))*dRpp_i +
!     &		       		(Wvel(i  ,j  ,k  )-Wvel(i-1,j  ,k  ))*dRp_i +
!     &		       		(Wvel(i  ,j  ,k-1)-Wvel(i-1,j  ,k-1))*dRp_i ) * 0.25
!				dwdy = ((Wvel(i  ,j+1,k  )-Wvel(i  ,j  ,k  ))/(Rp(i)*(phip(j+1)-phip(j))) +
!     &		    		(Wvel(i  ,j+1,k-1)-Wvel(i  ,j  ,k-1))/(Rp(i)*(phip(j+1)-phip(j))) +
!     &		       		(Wvel(i  ,j  ,k  )-Wvel(i  ,j-1,k  ))/(Rp(i)*(phip(j)-phip(j-1))) +
!     &		        	(Wvel(i  ,j  ,k-1)-Wvel(i  ,j-1,k-1))/(Rp(i)*(phip(j)-phip(j-1))) ) * 0.25
!
!				S11 =      dudx			!strainrate tensor elements
!				S22 =      dvdy
!				S33 =      dwdz
!				S12 = 0.5*(dvdx+dudy)
!				S13 = 0.5*(dudz+dwdx)
!				S23 = 0.5*(dvdz+dwdy)
!				SijSij = S11*S11 + S22*S22 + S33*S33 + 2.*S12*S12 + 2.*S13*S13 + 2.*S23*S23
!				strain(i,j,k)=sqrt(2*SijSij)
				
				shear = 2.0*(
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)) +	
     1			((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi) )
     
     
				shear = shear + 0.25*(
	1			((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))/(Ru(i)*(phip(j+1)-phip(j))) +
     2			( Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))/(Ru(i)*(phip(j)-phip(j-1))) +
     2			( Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))/(Ru(i-1)*(phip(j+1)-phip(j))) +
     2			( Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2 +
	1			((Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))/(Ru(i-1)*(phip(j)-phip(j-1))) +
     2			( Vvel(i  ,j-1,k  )/Rp(i)-Vvel(i-1,j-1,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2
     e				)

				shear = shear + 0.25*(
	1			((Uvel(i  ,j  ,k+1)-Uvel(i  ,j  ,k  ))*dzi +
     2			( Wvel(i+1,j  ,k  )-Wvel(i  ,j  ,k  ))*dRpp_i)**2 +
	1			((Uvel(i  ,j  ,k  )-Uvel(i  ,j  ,k-1))*dzi +
     2			( Wvel(i+1,j  ,k-1)-Wvel(i  ,j  ,k-1))*dRpp_i)**2 +
	1			((Uvel(i-1,j  ,k+1)-Uvel(i-1,j  ,k  ))*dzi +
     2			( Wvel(i  ,j  ,k  )-Wvel(i-1,j  ,k  ))*dRp_i)**2 +
	1			((Uvel(i-1,j  ,k  )-Uvel(i-1,j  ,k-1))*dzi +
     2			( Wvel(i  ,j  ,k-1)-Wvel(i-1,j  ,k-1))*dRp_i)**2 
     e				)
	 
				shear = shear + 0.25*(
	1			((Vvel(i  ,j  ,k+1)-Vvel(i  ,j  ,k  ))*dzi +
     2			( Wvel(i  ,j+1,k  )-Wvel(i  ,j  ,k  ))/(Rp(i)*(phip(j+1)-phip(j))) )**2 +
	1			((Vvel(i  ,j  ,k  )-Vvel(i  ,j  ,k-1))*dzi +
     2			( Wvel(i  ,j+1,k-1)-Wvel(i  ,j  ,k-1))/(Rp(i)*(phip(j+1)-phip(j))) )**2 +
	1			((Vvel(i  ,j-1,k+1)-Vvel(i  ,j-1,k  ))*dzi +
     2			( Wvel(i  ,j  ,k  )-Wvel(i  ,j-1,k  ))/(Rp(i)*(phip(j)-phip(j-1))) )**2 +
	1			((Vvel(i  ,j-1,k  )-Vvel(i  ,j-1,k-1))*dzi +
     2			( Wvel(i  ,j  ,k-1)-Wvel(i  ,j-1,k-1))/(Rp(i)*(phip(j)-phip(j-1))) )**2 
     e				)

!! 
      divergentie= 2./3.*(( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) ) / ( dz ))   
       
		shear = shear 
     +             +1.5*divergentie*divergentie
     +             -divergentie*2.*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))
     +             -divergentie*2.*((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1)))+0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))
     +             -divergentie*2.*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)

				strain(i,j,k) = Apvisc_shear_relax*sqrt(shear) + (1.-Apvisc_shear_relax)*strain(i,j,k) 
			enddo
		enddo
	  enddo	
	ENDIF 
	
		IF (rheo_shear_method.eq.21) THEN ! re-use rheo_shear_method also for other possibility in determination shear
		  strain2=strain
		  call bound_3D(strain2) 
		  do i=1,imax
			do j=1,jmax
				do k=1,kmax	
					im=i-1 
					ip=i+1 
					jm=j-1 
					jp=j+1 
					km=k-1 
					kp=k+1 
					strain(i,j,k) = MIN(strain2(i,j,k),strain2(im,j,k),strain2(ip,j,k),strain2(i,jm,k),strain2(i,jp,k),
     &				strain2(i,j,km),strain2(i,j,kp))
				enddo
			enddo
		  enddo	
		ELSEIF (rheo_shear_method.eq.31) THEN
		  do i=1,imax
			do j=1,jmax
				do k=1,kmax	
				   IF (SUM(cnew(1:nfrac,i,j,k))<1.e-6) THEN
				     strain(i,j,k)=0.  !forces strain to 0 for cells outside patch of rheology material
				   ENDIF 
				enddo
			enddo
		  enddo	
		ENDIF 
	
!!	IF (Apvisc_interp.eq.3.or.Apvisc_interp.eq.4) THEN ! apply high muA at edges patch for staggered arrangement:
!!		strain2=strain
!!		call bound_3D(strain2) 
!!		do i=1,imax
!!			do j=1,jmax
!!				do k=1,kmax	
!!					im=i-1 
!!					ip=i+1 
!!					jm=j-1 
!!					jp=j+1 
!!					km=k-1 
!!					kp=k+1 
!!					strain(i,j,k) = MIN(
!!     &					strain2(i,j,k),strain2(i,jm,k),strain2(i,jp,k),strain2(i,j,km),strain2(i,j,kp),
!!     &					strain2(i,jm,km),strain2(i,jp,km),strain2(i,jm,kp),strain2(i,jp,kp),
!!     &					strain2(im,j,k),strain2(im,jm,k),strain2(im,jp,k),strain2(im,j,km),strain2(im,j,kp),
!!     &					strain2(im,jm,km),strain2(im,jp,km),strain2(im,jm,kp),strain2(im,jp,kp),
!!     &					strain2(ip,j,k),strain2(ip,jm,k),strain2(ip,jp,k),strain2(ip,j,km),strain2(ip,j,kp),
!!     &					strain2(ip,jm,km),strain2(ip,jp,km),strain2(ip,jm,kp),strain2(ip,jp,kp))
!!				enddo
!!			enddo
!!		enddo
!!	ENDIF 
	
	end subroutine strain_magnitude
	
	subroutine stress_magnitude
	USE nlist
	implicit none
	!include 'mpif.h'
	
	do i=1,imax
		do j=1,jmax
			do k=1,kmax
				stress(i,j,k)=muA(i,j,k)*strain(i,j,k)
			enddo
		enddo
	enddo
	end subroutine stress_magnitude

	subroutine Simple_Bingham
	
	USE nlist
	implicit none
	!include 'mpif.h'
	
	integer n
	
	IF (MAXVAL(SIMPLE_climit)>1.e-12) THEN 
		muB=0.
		tauY=0.
		do n=1,nfrac
			do i=1,imax
				do j=1,jmax
					do k=1,kmax
						IF (Cnew(n,i,j,k)>SIMPLE_climit(n)) THEN
							muB(i,j,k) = SIMPLE_muB
							tauY(i,j,k) = SIMPLE_tauy
						ENDIF
					enddo
				enddo
			enddo
		enddo
	ELSE
		muB(1:imax,1:jmax,1:kmax) = SIMPLE_muB
		tauY(1:imax,1:jmax,1:kmax) = SIMPLE_tauy
	ENDIF 

	end subroutine Simple_Bingham


	subroutine Rheo_Jacobs_and_vKesteren
	
	USE nlist
	implicit none
	!include 'mpif.h'
	
	real phi_clay(0:i1,0:j1,0:k1),phi_silt(0:i1,0:j1,0:k1),phi_sand(0:i1,0:j1,0:k1)
	real phi_solids(0:i1,0:j1,0:k1),phi_fines(0:i1,0:j1,0:k1)
	real sol_frac(0:i1,0:j1,0:k1),W_rel(0:i1,0:j1,0:k1)
	real lambda_b(0:i1,0:j1,0:k1),sol_eff(0:i1,0:j1,0:k1)
	integer n

	real rho_sed			!temporary
	
	rho_sed=2650.
	phi_clay(:,:,:)=0.
	phi_silt(:,:,:)=0.
	phi_sand(:,:,:)=0.
	phi_solids(:,:,:)=0.
	phi_fines(:,:,:)=0.
	
	do i=1,imax 
		do j=1,jmax
			do k=1,kmax
				do n=1,nfrac
					if (frac(n)%dpart.le.2) then								!check for clay fractions
					phi_clay(i,j,k)= phi_clay(i,j,k) + cnew(n,i,j,k)
					elseif (frac(n)%dpart.gt.2.and.frac(n)%dpart.le.63) then 	!check for silt fractions
					phi_silt(i,j,k)= phi_silt(i,j,k) + cnew(n,i,j,k)
					else
					phi_sand(i,j,k)= phi_sand(i,j,k) + cnew(n,i,j,k)
					endif
				enddo
				phi_fines(i,j,k)=phi_clay(i,j,k)+phi_silt(i,j,k)
				phi_solids(i,j,k)=phi_fines(i,j,k)+phi_sand(i,j,k)
				lambda_b(i,j,k)=1./(((BAGNOLD_phi_max/(1.e-12+phi_sand(i,j,k)))**(1./3.))-1.)				!linear sand concentration Bagnold 1956
				sol_eff(i,j,k)=exp(BAGNOLD_beta*lambda_b(i,j,k))											!exponential term expressing influens of sand and silt
				W_rel(i,j,k)=(rho_b/(JACOBS_Aclay*rho_sed))*((1.-phi_solids(i,j,k))/(1.e-12+phi_fines(i,j,k)))		!relative water content
				tauY(i,j,k)=MIN(MAX_tauy,sol_eff(i,j,k)*JACOBS_Ky*(W_rel(i,j,k)**JACOBS_By))
				muB(i,j,k)=MIN(MAX_mu,sol_eff(i,j,k)*(JACOBS_muw+JACOBS_Kmu*(W_rel(i,j,k)**JACOBS_Bmu)))
			enddo
		enddo
	enddo

	end subroutine Rheo_Jacobs_and_vKesteren

	subroutine Rheo_Winterwerp_and_Kranenburg
	
	USE nlist
	implicit none
	!include 'mpif.h'

	real phi_clay(0:i1,0:j1,0:k1),phi_silt(0:i1,0:j1,0:k1),phi_sand(0:i1,0:j1,0:k1)
	real phi_solids(0:i1,0:j1,0:k1),phi_fines(0:i1,0:j1,0:k1)
	real lambda_b(0:i1,0:j1,0:k1),sol_eff(0:i1,0:j1,0:k1),sol_frac(0:i1,0:j1,0:k1)
	real shear_thin(0:i1,0:j1,0:k1)
	integer n
	
	phi_clay(:,:,:)=0.
	phi_silt(:,:,:)=0.
	phi_sand(:,:,:)=0.
	phi_solids(:,:,:)=0.
	phi_fines(:,:,:)=0.
	
	call strain_magnitude(Unew,Vnew,Wnew)
	do i=1,imax
		do j=1,jmax
			do k=1,kmax
				do n=1,nfrac
					if (frac(n)%dpart.le.2) then								!check for clay fractions
					phi_clay(i,j,k)= phi_clay(i,j,k) + cnew(n,i,j,k)
					elseif (frac(n)%dpart.gt.2.and.frac(n)%dpart.le.63) then 	!check for silt fractions
					phi_silt(i,j,k)= phi_silt(i,j,k) + cnew(n,i,j,k)
					else
					phi_sand(i,j,k)= phi_sand(i,j,k) + cnew(n,i,j,k)
					endif
				enddo
				phi_fines(i,j,k)=phi_clay(i,j,k)+phi_silt(i,j,k)
				phi_solids(i,j,k)=phi_fines(i,j,k)+phi_sand(i,j,k)
				lambda_b(i,j,k)=1./(((BAGNOLD_phi_max/(1.e-12+phi_sand(i,j,k)))**(1./3.))-1.)				!linear sand concentration Bagnold 1956
				sol_eff(i,j,k)=exp(BAGNOLD_beta*lambda_b(i,j,k))									!exponential term expressing influens of sand and silt
				sol_frac(i,j,k)=phi_fines(i,j,k)/(1.-phi_sand(i,j,k))
				tauY(i,j,k)=MIN(MAX_tauy,sol_eff(i,j,k)*WINTER_Ay*sol_frac(i,j,k)**(2./(3.-WINTER_nf)))
				shear_thin(i,j,k)=((1./(1.e-12+strain(i,j,k)))**(((WINTER_af+1.)*(3.-WINTER_nf))/3.))
				muB(i,j,k)= sol_eff(i,j,k)*(WINTER_muw+WINTER_Amu*(sol_frac(i,j,k)**(((2.*(WINTER_af+1.))/3.)))*shear_thin(i,j,k))
				muB(i,j,k)=MIN(MAX_mu,muB(i,j,k))
			enddo
		enddo
	enddo
	end subroutine Rheo_Winterwerp_and_Kranenburg

	subroutine Rheo_Thomas

	USE nlist
	implicit none
	!include 'mpif.h'
	
	real phi_sand(0:i1,0:j1,0:k1),phi_fines(0:i1,0:j1,0:k1),phi_solids(0:i1,0:j1,0:k1)
	real sol_eff_y(0:i1,0:j1,0:k1),sol_eff_mu(0:i1,0:j1,0:k1)
	integer n,n1
	
	phi_solids(:,:,:)=0.
	phi_fines(:,:,:)=0.
	phi_sand(:,:,:)=0.

	do i=1,imax 
		do j=1,jmax
			do k=1,kmax
				do n1=1,nfr_silt 
					n=nfrac_silt(n1)				
					phi_fines(i,j,k)= phi_fines(i,j,k)+cnew(n,i,j,k)
				enddo 
				do n1=1,nfr_sand 
					n=nfrac_sand(n1)				
					phi_sand(i,j,k)= phi_sand(i,j,k)+cnew(n,i,j,k)
				enddo 
!				do n=1,nfrac
!					if (frac(n)%type.eq.1) then								!fines < 45 micron
!					elseif (frac(n)%type.eq.2) then 							!To include all fractions we let sand > 44 micron (instead of 63 micron)
!					phi_sand(i,j,k)= phi_sand(i,j,k)+cnew(n,i,j,k)
!					endif
!				enddo
				phi_solids(i,j,k)=MIN(THOMAS_phi_sand_max,phi_fines(i,j,k)+phi_sand(i,j,k))
				phi_sand(i,j,k)=MIN(THOMAS_phi_sand_max,phi_sand(i,j,k))
				phi_fines(i,j,k)=MIN(THOMAS_phi_sand_max,phi_fines(i,j,k))
				sol_eff_y(i,j,k)=(1.-(phi_sand(i,j,k)/(THOMAS_ky*THOMAS_phi_sand_max)))**-2.5		!solids effect
				sol_eff_mu(i,j,k)=(1.-(phi_sand(i,j,k)/(THOMAS_kmu*THOMAS_phi_sand_max)))**-2.5
				tauY(i,j,k)=MIN(MAX_tauy,THOMAS_Cy*((phi_fines(i,j,k)/(1.-phi_sand(i,j,k)))**THOMAS_Py)*sol_eff_y(i,j,k))		!depends on definition of fines (in-/excluding silt)
				muB(i,j,k)=MIN(MAX_mu,THOMAS_Cmu*exp(THOMAS_Pmu*(phi_fines(i,j,k)/(1.-phi_solids(i,j,k))))*sol_eff_mu(i,j,k))
			enddo
		enddo
	enddo
	end subroutine Rheo_Thomas

	subroutine Bingham_Papanastasiou							!tauY and muB, yield stress and bingham viscosity (rr = density)
	USE nlist
	implicit none
	!include 'mpif.h'

	real ebb(0:i1,0:k1)
	real ebf(0:i1,0:k1)
	real muA2(0:i1,0:j1,0:k1)
	real tau_app(0:i1,0:j1,0:k1),tau_app_u(0:i1,0:j1,0:k1),tau_app_v(0:i1,0:j1,0:k1)
	real dppdx,dppdy,dy2,dx2,pppp(0:i1,0:j1,0:k1)
	real driving_shearstressz,driving_shearstressy,driving_shearstressx,driving_shearstress,dpdz_hydr
	integer im,ip,jm,jp,km,kp,n

	if (Rheological_model.ne.'WINTER') then
		call strain_magnitude(Unew,Vnew,Wnew)				!determine the magnite of the strain rate
	endif


	tau_app=0.
	if (Apvisc_force_eq.eq.1) then 
	!! try to split off horizontal dpdx,dpdy part of BYS to arrive at horizontal equilibrium numerically more stable
	    pppp(1:imax,1:jmax,1:kmax)=Pold+p !now full pold not just dp
		call bound_3D(pppp) 
		call bound_3D(tauY) 
		do i=1,imax
			do j=1,jmax
				do k=1,kmax	
					im=i-1 
					ip=i+1 
					jm=j-1 
					jp=j+1 
					!! at U-loc:
					dppdx = (pppp(ip,j,k)-pppp(i,j,k))/(Rp(ip)-Rp(i))
					dppdy = (0.5*(pppp(i,jp,k)+pppp(ip,jp,k))-0.5*(pppp(i,jm,k)+pppp(ip,jm,k)))/( Ru(i) * (phip(jp)-phip(jm)) )
					tau_app_u(i,j,k) = MIN(tauY(i,j,k),tauY(ip,j,k),dz*sqrt(dppdx**2+dppdy**2))
					Ppropx(i,j,k) = Ppropx(i,j,k)+tau_app_u(i,j,k)/dz*dppdx/sqrt((dppdx)**2+(dppdy)**2+1.e-18) !sign of dppdx automatically gives right addition or subtraction to slow down driving force
					!Ppropx(i,j,k) = Ppropx(i,j,k)+dppdx !check -> should give near-zero velocity -> yes it does
					!! at V-loc:
					dppdy = (pppp(i,jp,k)-pppp(i,j,k))/( Rp(i) * (phip(jp)-phip(j)) ) !
					dppdx = (0.5*(pppp(ip,j,k)+pppp(ip,jp,k))-0.5*(pppp(im,j,k)+pppp(im,jp,k)))/(Rp(ip)-Rp(im))
					tau_app_v(i,j,k) = MIN(tauY(i,j,k),tauY(i,jp,k),dz*sqrt(dppdx**2+dppdy**2))
					Ppropy(i,j,k) = Ppropy(i,j,k)+tau_app_v(i,j,k)/dz*dppdy/sqrt((dppdx)**2+(dppdy)**2+1.e-18) !sign of dppdy automatically gives right addition or subtraction to slow down driving force	
					!Ppropy(i,j,k) = Ppropy(i,j,k)+dppdy !check -> should give near-zero velocity -> yes it does
				enddo
			enddo
		enddo		

		call bound_3D(tau_app_u)
		call bound_3D(tau_app_v)
		do i=1,imax
			do j=1,jmax
				do k=1,kmax	
					im=i-1 
					jm=j-1 
					tau_app(i,j,k) = 0.25*(tau_app_u(i,j,k)+tau_app_u(im,j,k)+tau_app_v(i,j,k)+tau_app_v(i,jm,k))
				enddo
			enddo
		enddo
	endif 
				
				
				
	
	
	IF (Apvisc_interp.eq.3.or.Apvisc_interp.eq.4) THEN ! apply high muA at edges patch for staggered arrangement:
	  DO n=1,1 !10
		muA2=tauY
		call bound_3D(muA2) 
		do i=1,imax
			do j=1,jmax
				do k=1,kmax	
					im=i-1 
					ip=i+1 
					jm=j-1 
					jp=j+1 
					km=k-1 
					kp=k+1 
					tauY(i,j,k) = MAX(
     &					muA2(i,j,k),muA2(i,jm,k),muA2(i,jp,k),muA2(i,j,km),muA2(i,j,kp),
     &					muA2(i,jm,km),muA2(i,jp,km),muA2(i,jm,kp),muA2(i,jp,kp),
     &					muA2(im,j,k),muA2(im,jm,k),muA2(im,jp,k),muA2(im,j,km),muA2(im,j,kp),
     &					muA2(im,jm,km),muA2(im,jp,km),muA2(im,jm,kp),muA2(im,jp,kp),
     &					muA2(ip,j,k),muA2(ip,jm,k),muA2(ip,jp,k),muA2(ip,j,km),muA2(ip,j,kp),
     &					muA2(ip,jm,km),muA2(ip,jp,km),muA2(ip,jm,kp),muA2(ip,jp,kp))
				enddo
			enddo
		enddo
	  ENDDO 
	ENDIF 
	
	
	if (Apvisc_force_eq.eq.2) then 
	    pppp(1:imax,1:jmax,1:kmax)=Pold+p !now full pold not just dp
		call bound_3D(pppp) 	
		do i=1,imax						!start loop, in r-direction
			do j=1,jmax						!start loop, in phi-direction
				do k=1,kmax					!start loop, in z-direction
				  ! unyielded when: delta_Pz*dx*dy < BYS*dz*dy  or delta_Pz*dx*dy < BYS*dz*dx & 
				  !                 delta_Py*dz*dx < BYS*dy*dx  or delta_Py*dz*dx < BYS*dy*dz  
				  !	                delta_Px*dz*dy < BYS*dx*dy  or delta_Px*dz*dy < BYS*dx*dz & 		
				  ! assume very wide cells in laterally uniform situation: dy>>> then check #2 and #6 give yield
				  ! or assume very wide cells in x-dir for uniform situation in x-dir with dx>>> then check #1 and #4 give yield
				  ! or assume very high cells dz>>> then check #3 and #5 give yield 
				  ! each force-dir needs to be taken up by either one of the edges where BYS can counteract the pressure gradient, doesn't need to be able to counteract on both edges for each dir 
				  dpdz_hydr = (rho_b-rnew(i,j,k))*gz !negative when rnew>rho_b
				  driving_shearstressz = ABS((pppp(i,j,k+1)-pppp(i,j,k-1))/(2.*dz)-dpdz_hydr)*dr(i)
				  driving_shearstressz = MIN(driving_shearstressz,ABS((pppp(i,j,k+1)-pppp(i,j,k-1))/(2.*dz)-dpdz_hydr)*Rp(i)*dphi2(j))
				  
				  driving_shearstressy = ABS(pppp(i,j+1,k)-pppp(i,j-1,k))*dr(i)/(2.*Rp(i)*dphi2(j))
				  driving_shearstressy = MIN(driving_shearstressy,ABS(pppp(i,j+1,k)-pppp(i,j-1,k))*dz/(2.*Rp(i)*dphi2(j)))
				  
				  driving_shearstressx = ABS(pppp(i+1,j,k)-pppp(i-1,j,k))*dz/(2.*dr(i))
				  driving_shearstressx = MIN(driving_shearstressx,ABS(pppp(i+1,j,k)-pppp(i-1,j,k))*Rp(i)*dphi2(j)/(2.*dr(i)))
				  driving_shearstress = MAX(driving_shearstressx,driving_shearstressy,driving_shearstressz) 
				  if (driving_shearstress<tauY(i,j,k)) then !unyielded zone with maximum appararent viscosity 
				    muA(i,j,k)= muB(i,j,k)+(tauY(i,j,k)-tau_app(i,j,k))*PAPANASTASIOUS_m
				  elseif (strain(i,j,k)>1.e-8) then 
					muA(i,j,k)= muB(i,j,k)+MAX(0.,tauY(i,j,k)-tau_app(i,j,k))/(max(1.e-12,strain(i,j,k)-shear0limit))*
     &				(1.-exp(-PAPANASTASIOUS_m*(max(1.e-12,strain(i,j,k)-shear0limit))))
				  else 
					muA(i,j,k)= muB(i,j,k)+(tauY(i,j,k)-tau_app(i,j,k))*PAPANASTASIOUS_m
				  endif 
				enddo
			enddo
		enddo
	else 
		do i=1,imax						!start loop, in r-direction
			do j=1,jmax						!start loop, in phi-direction
				do k=1,kmax					!start loop, in z-direction
				  if (strain(i,j,k)>1.e-8) then 
					muA(i,j,k)= muB(i,j,k)+MAX(0.,tauY(i,j,k)-tau_app(i,j,k))/(max(1.e-12,strain(i,j,k)-shear0limit))*
     &				(1.-exp(-PAPANASTASIOUS_m*(max(1.e-12,strain(i,j,k)-shear0limit))))
				  else 
					muA(i,j,k)= muB(i,j,k)+(tauY(i,j,k)-tau_app(i,j,k))*PAPANASTASIOUS_m
				  endif 
				enddo
			enddo
		enddo
	endif 
!!	IF (Apvisc_interp.eq.3.or.Apvisc_interp.eq.4) THEN ! apply high muA at edges patch for staggered arrangement:
!!	  DO n=1,5
!!		muA2=muA
!!		call bound_3D(muA2) 
!!		do i=1,imax
!!			do j=1,jmax
!!				do k=1,kmax	
!!					im=i-1 
!!					ip=i+1 
!!					jm=j-1 
!!					jp=j+1 
!!					km=k-1 
!!					kp=k+1 
!!					muA(i,j,k) = MAX(
!!     &					muA2(i,j,k),muA2(i,jm,k),muA2(i,jp,k),muA2(i,j,km),muA2(i,j,kp),
!!     &					muA2(i,jm,km),muA2(i,jp,km),muA2(i,jm,kp),muA2(i,jp,kp),
!!     &					muA2(im,j,k),muA2(im,jm,k),muA2(im,jp,k),muA2(im,j,km),muA2(im,j,kp),
!!     &					muA2(im,jm,km),muA2(im,jp,km),muA2(im,jm,kp),muA2(im,jp,kp),
!!     &					muA2(ip,j,k),muA2(ip,jm,k),muA2(ip,jp,k),muA2(ip,j,km),muA2(ip,j,kp),
!!     &					muA2(ip,jm,km),muA2(ip,jp,km),muA2(ip,jm,kp),muA2(ip,jp,kp))
!!				enddo
!!			enddo
!!		enddo
!!	  ENDDO 
!!	ENDIF 
	ekm= ekm + muA - ekm_mol !in muB (part of muA) ekm_mol is included  
	call stress_magnitude !calculate the stress magnitude
	call bound_3D(ekm)
	call bound_3D(muA) 
	
	END subroutine Bingham_Papanastasiou
	
	
	subroutine Bingham_Fluent_manner							!tauY and muB, yield stress and bingham viscosity (rr = density)
	USE nlist
	implicit none
	!include 'mpif.h'

	real ebb(0:i1,0:k1)
	real ebf(0:i1,0:k1)
	real muA2(0:i1,0:j1,0:k1)
	real tau_app(0:i1,0:j1,0:k1),tau_app_u(0:i1,0:j1,0:k1),tau_app_v(0:i1,0:j1,0:k1)
	real dppdx,dppdy,dy2,dx2,pppp(0:i1,0:j1,0:k1)
	real driving_shearstressz,driving_shearstressy,driving_shearstressx,driving_shearstress,dpdz_hydr
	integer im,ip,jm,jp,km,kp,n

	if (Rheological_model.ne.'WINTER') then
		call strain_magnitude(Unew,Vnew,Wnew)				!determine the magnite of the strain rate
	endif


	tau_app=0.
	if (Apvisc_force_eq.eq.1) then 
	!! try to split off horizontal dpdx,dpdy part of BYS to arrive at horizontal equilibrium numerically more stable
	    pppp(1:imax,1:jmax,1:kmax)=Pold+p !now full pold not just dp
		call bound_3D(pppp) 
		call bound_3D(tauY) 
		do i=1,imax
			do j=1,jmax
				do k=1,kmax	
					im=i-1 
					ip=i+1 
					jm=j-1 
					jp=j+1 
					!! at U-loc:
					dppdx = (pppp(ip,j,k)-pppp(i,j,k))/(Rp(ip)-Rp(i))
					dppdy = (0.5*(pppp(i,jp,k)+pppp(ip,jp,k))-0.5*(pppp(i,jm,k)+pppp(ip,jm,k)))/( Ru(i) * (phip(jp)-phip(jm)) )
					tau_app_u(i,j,k) = MIN(tauY(i,j,k),tauY(ip,j,k),dz*sqrt(dppdx**2+dppdy**2))
					Ppropx(i,j,k) = Ppropx(i,j,k)+tau_app_u(i,j,k)/dz*dppdx/sqrt((dppdx)**2+(dppdy)**2+1.e-18) !sign of dppdx automatically gives right addition or subtraction to slow down driving force
					!Ppropx(i,j,k) = Ppropx(i,j,k)+dppdx !check -> should give near-zero velocity -> yes it does
					!! at V-loc:
					dppdy = (pppp(i,jp,k)-pppp(i,j,k))/( Rp(i) * (phip(jp)-phip(j)) ) !
					dppdx = (0.5*(pppp(ip,j,k)+pppp(ip,jp,k))-0.5*(pppp(im,j,k)+pppp(im,jp,k)))/(Rp(ip)-Rp(im))
					tau_app_v(i,j,k) = MIN(tauY(i,j,k),tauY(i,jp,k),dz*sqrt(dppdx**2+dppdy**2))
					Ppropy(i,j,k) = Ppropy(i,j,k)+tau_app_v(i,j,k)/dz*dppdy/sqrt((dppdx)**2+(dppdy)**2+1.e-18) !sign of dppdy automatically gives right addition or subtraction to slow down driving force	
					!Ppropy(i,j,k) = Ppropy(i,j,k)+dppdy !check -> should give near-zero velocity -> yes it does
				enddo
			enddo
		enddo		

		call bound_3D(tau_app_u)
		call bound_3D(tau_app_v)
		do i=1,imax
			do j=1,jmax
				do k=1,kmax	
					im=i-1 
					jm=j-1 
					tau_app(i,j,k) = 0.25*(tau_app_u(i,j,k)+tau_app_u(im,j,k)+tau_app_v(i,j,k)+tau_app_v(i,jm,k))
				enddo
			enddo
		enddo
	endif 
					
	
	IF (Apvisc_interp.eq.3.or.Apvisc_interp.eq.4) THEN ! apply high muA at edges patch for staggered arrangement:
	  DO n=1,1 !10
		muA2=tauY
		call bound_3D(muA2) 
		do i=1,imax
			do j=1,jmax
				do k=1,kmax	
					im=i-1 
					ip=i+1 
					jm=j-1 
					jp=j+1 
					km=k-1 
					kp=k+1 
					tauY(i,j,k) = MAX(
     &					muA2(i,j,k),muA2(i,jm,k),muA2(i,jp,k),muA2(i,j,km),muA2(i,j,kp),
     &					muA2(i,jm,km),muA2(i,jp,km),muA2(i,jm,kp),muA2(i,jp,kp),
     &					muA2(im,j,k),muA2(im,jm,k),muA2(im,jp,k),muA2(im,j,km),muA2(im,j,kp),
     &					muA2(im,jm,km),muA2(im,jp,km),muA2(im,jm,kp),muA2(im,jp,kp),
     &					muA2(ip,j,k),muA2(ip,jm,k),muA2(ip,jp,k),muA2(ip,j,km),muA2(ip,j,kp),
     &					muA2(ip,jm,km),muA2(ip,jp,km),muA2(ip,jm,kp),muA2(ip,jp,kp))
				enddo
			enddo
		enddo
	  ENDDO 
	ENDIF 
	
	
	if (Apvisc_force_eq.eq.2) then 
	    pppp(1:imax,1:jmax,1:kmax)=Pold+p !now full pold not just dp
		call bound_3D(pppp) 	
		do i=1,imax						!start loop, in r-direction
			do j=1,jmax						!start loop, in phi-direction
				do k=1,kmax					!start loop, in z-direction
				  ! unyielded when: delta_Pz*dx*dy < BYS*dz*dy  or delta_Pz*dx*dy < BYS*dz*dx & 
				  !                 delta_Py*dz*dx < BYS*dy*dx  or delta_Py*dz*dx < BYS*dy*dz  
				  !	                delta_Px*dz*dy < BYS*dx*dy  or delta_Px*dz*dy < BYS*dx*dz & 		
				  ! assume very wide cells in laterally uniform situation: dy>>> then check #2 and #6 give yield
				  ! or assume very wide cells in x-dir for uniform situation in x-dir with dx>>> then check #1 and #4 give yield
				  ! or assume very high cells dz>>> then check #3 and #5 give yield 
				  ! each force-dir needs to be taken up by either one of the edges where BYS can counteract the pressure gradient, doesn't need to be able to counteract on both edges for each dir 
				  dpdz_hydr = (rho_b-rnew(i,j,k))*gz !negative when rnew>rho_b
				  driving_shearstressz = ABS((pppp(i,j,k+1)-pppp(i,j,k-1))/(2.*dz)-dpdz_hydr)*dr(i)
				  driving_shearstressz = MIN(driving_shearstressz,ABS((pppp(i,j,k+1)-pppp(i,j,k-1))/(2.*dz)-dpdz_hydr)*Rp(i)*dphi2(j))
				  
				  driving_shearstressy = ABS(pppp(i,j+1,k)-pppp(i,j-1,k))*dr(i)/(2.*Rp(i)*dphi2(j))
				  driving_shearstressy = MIN(driving_shearstressy,ABS(pppp(i,j+1,k)-pppp(i,j-1,k))*dz/(2.*Rp(i)*dphi2(j)))
				  
				  driving_shearstressx = ABS(pppp(i+1,j,k)-pppp(i-1,j,k))*dz/(2.*dr(i))
				  driving_shearstressx = MIN(driving_shearstressx,ABS(pppp(i+1,j,k)-pppp(i-1,j,k))*Rp(i)*dphi2(j)/(2.*dr(i)))
				  driving_shearstress = MAX(driving_shearstressx,driving_shearstressy,driving_shearstressz) 
				  if (driving_shearstress<tauY(i,j,k)) then !unyielded zone with maximum appararent viscosity 
					muA(i,j,k)= muB(i,j,k)+2.*MAX(0.,tauY(i,j,k)-tau_app(i,j,k))/shear0limit					
				  elseif (strain(i,j,k)>shear0limit) then 
				    muA(i,j,k)= muB(i,j,k)+MAX(0.,tauY(i,j,k)-tau_app(i,j,k))/strain(i,j,k)
				  else 
					muA(i,j,k)= muB(i,j,k)+(2.-strain(i,j,k)/shear0limit)*MAX(0.,tauY(i,j,k)-tau_app(i,j,k))/shear0limit
				  endif 
				enddo
			enddo
		enddo
	else 
		do i=1,imax						!start loop, in r-direction
			do j=1,jmax						!start loop, in phi-direction
				do k=1,kmax					!start loop, in z-direction
				  if (strain(i,j,k)>shear0limit) then 
				    muA(i,j,k)= muB(i,j,k)+MAX(0.,tauY(i,j,k)-tau_app(i,j,k))/strain(i,j,k)
				  else 
					muA(i,j,k)= muB(i,j,k)+(2.-strain(i,j,k)/shear0limit)*MAX(0.,tauY(i,j,k)-tau_app(i,j,k))/shear0limit
				  endif 
				enddo
			enddo
		enddo
	endif 

	ekm= ekm + muA - ekm_mol !in muB (part of muA) ekm_mol is included  
	call stress_magnitude !calculate the stress magnitude
	call bound_3D(ekm)
	call bound_3D(muA) 
	
	END subroutine Bingham_Fluent_manner	
	

