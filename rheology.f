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
	integer kbed_dummy(0:i1,0:j1)
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
	call bound_p(lambda_new) !Neumann B.C.
		
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
	
	real dzi,n
	real dRpp_i,dRp_i
	real Uvel(0:i1,0:j1,0:k1),Vvel(0:i1,0:j1,0:k1),
     +   Wvel(0:i1,0:j1,0:k1),rr(0:i1,0:j1,0:k1)
	real dudx,dudy,dudz
	real dvdx,dvdy,dvdz
	real dwdx,dwdy,dwdz
	real S11,S12,S13,S22,S23,S33,SijSij

	SijSij=0.

	dzi=1./dz
	do i=1,imax
		dRpp_i=1./(Rp(i+1)-Rp(i))
		dRp_i=1./(Rp(i)-Rp(i-1))	  
		do j=1,jmax
			do k=1,kmax
				dudx = (Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i) 																	!du/dx
				dvdy = (Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*(phiv(j)-phiv(j-1))) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)		!dv/dy
				dwdz = (Wvel(i,j,k)-Wvel(i,j,k-1))*dzi																		!dw/dz

				dudy = ((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))/(Ru(i)*(phip(j+1)-phip(j))) + 
     &		        	(Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))/(Ru(i)*(phip(j)-phip(j-1))) +
     &		        	(Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))/(Ru(i-1)*(phip(j+1)-phip(j))) +
     &		        	(Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))/(Ru(i-1)*(phip(j)-phip(j-1))) ) * 0.25
				dudz = ((Uvel(i  ,j  ,k+1)-Uvel(i  ,j  ,k  ))*dzi +
     &		        	(Uvel(i  ,j  ,k  )-Uvel(i  ,j  ,k-1))*dzi +
     &		        	(Uvel(i-1,j  ,k+1)-Uvel(i-1,j  ,k  ))*dzi +
     &		        	(Uvel(i-1,j  ,k  )-Uvel(i-1,j  ,k-1))*dzi ) * 0.25
				dvdx = ((Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) +
     &		       		(Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) + 
     &		       		(Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1)  +
     &		       		(Vvel(i  ,j-1,k  )/Rp(i)-Vvel(i-1,j-1,k  )/Rp(i-1))*dRp_i*Ru(i-1) ) * 0.25
				dvdz = ((Vvel(i  ,j  ,k+1)-Vvel(i  ,j  ,k  ))*dzi +
     &		        	(Vvel(i  ,j  ,k  )-Vvel(i  ,j  ,k-1))*dzi +
     &		        	(Vvel(i  ,j-1,k+1)-Vvel(i  ,j-1,k  ))*dzi +
     &		        	(Vvel(i  ,j-1,k  )-Vvel(i  ,j-1,k-1))*dzi ) * 0.25
				dwdx = ((Wvel(i+1,j  ,k  )-Wvel(i  ,j  ,k  ))*dRpp_i + 
     &		       		(Wvel(i+1,j  ,k-1)-Wvel(i  ,j  ,k-1))*dRpp_i +
     &		       		(Wvel(i  ,j  ,k  )-Wvel(i-1,j  ,k  ))*dRp_i +
     &		       		(Wvel(i  ,j  ,k-1)-Wvel(i-1,j  ,k-1))*dRp_i ) * 0.25
				dwdy = ((Wvel(i  ,j+1,k  )-Wvel(i  ,j  ,k  ))/(Rp(i)*(phip(j+1)-phip(j))) +
     &		    		(Wvel(i  ,j+1,k-1)-Wvel(i  ,j  ,k-1))/(Rp(i)*(phip(j+1)-phip(j))) +
     &		       		(Wvel(i  ,j  ,k  )-Wvel(i  ,j-1,k  ))/(Rp(i)*(phip(j)-phip(j-1))) +
     &		        	(Wvel(i  ,j  ,k-1)-Wvel(i  ,j-1,k-1))/(Rp(i)*(phip(j)-phip(j-1))) ) * 0.25

				S11 =      dudx			!strainrate tensor elements
				S22 =      dvdy
				S33 =      dwdz
				S12 = 0.5*(dvdx+dudy)
				S13 = 0.5*(dudz+dwdx)
				S23 = 0.5*(dvdz+dwdy)
				SijSij = S11*S11 + S22*S22 + S33*S33 + 2.*S12*S12 + 2.*S13*S13 + 2.*S23*S23
				strain(i,j,k)=sqrt(2*SijSij)
			enddo
		enddo
	enddo
	
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

	do i=1,imax
		do j=1,jmax
			do k=1,kmax
				muB(i,j,k) = SIMPLE_muB
				tauY(i,j,k) = SIMPLE_tauy
			enddo
		enddo
	enddo

	end subroutine Simple_Bingham


	subroutine Rheo_Jacobs_and_vKesteren
	
	USE nlist
	implicit none
	!include 'mpif.h'
	
	real phi_clay(0:i1,0:j1,0:k1),phi_silt(0:i1,0:j1,0:k1),phi_sand(0:i1,0:j1,0:k1)
	real phi_solids(0:i1,0:j1,0:k1),phi_fines(0:i1,0:j1,0:k1)
	real sol_frac(0:i1,0:j1,0:k1),W_rel(0:i1,0:j1,0:k1)
	real lambda_b(0:i1,0:j1,0:k1),sol_eff(0:i1,0:j1,0:k1)
	real n

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
				tauY(i,j,k)=sol_eff(i,j,k)*JACOBS_Ky*(W_rel(i,j,k)**JACOBS_By)
				muB(i,j,k)=sol_eff(i,j,k)*(JACOBS_muw+JACOBS_Kmu*(W_rel(i,j,k)**JACOBS_Bmu))
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
	real n
	
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
				tauY(i,j,k)=sol_eff(i,j,k)*WINTER_Ay*sol_frac(i,j,k)**(2./(3.-WINTER_nf))
				shear_thin(i,j,k)=((1./(1.e-12+strain(i,j,k)))**(((WINTER_af+1.)*(3.-WINTER_nf))/3.))
				muB(i,j,k)= sol_eff(i,j,k)*(WINTER_muw+WINTER_Amu*(sol_frac(i,j,k)**(((2.*(WINTER_af+1.))/3.)))*shear_thin(i,j,k))
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
	real n
	
	phi_solids(:,:,:)=0.
	phi_fines(:,:,:)=0.
	phi_sand(:,:,:)=0.

	do i=1,imax 
		do j=1,jmax
			do k=1,kmax
				do n=1,nfrac
					if (frac(n)%dpart.le.44) then								!fines < 45 micron
					phi_fines(i,j,k)= phi_fines(i,j,k)+cnew(n,i,j,k)
					elseif (frac(n)%dpart.gt.44) then 							!To include all fractions we let sand > 44 micron (instead of 63 micron)
					phi_sand(i,j,k)= phi_sand(i,j,k)+cnew(n,i,j,k)
					endif
				enddo
				phi_solids(i,j,k)=phi_fines(i,j,k)+phi_sand(i,j,k)
				sol_eff_y(i,j,k)=(1.-(phi_sand(i,j,k)/(THOMAS_ky*THOMAS_phi_sand_max)))**-2.5		!solids effect
				sol_eff_mu(i,j,k)=(1.-(phi_sand(i,j,k)/(THOMAS_kmu*THOMAS_phi_sand_max)))**-2.5
				tauY(i,j,k)=THOMAS_Cy*((phi_fines(i,j,k)/(1.-phi_sand(i,j,k)))**THOMAS_Py)*sol_eff_y(i,j,k)		!depends on definition of fines (in-/excluding silt)
				muB(i,j,k)=THOMAS_Cmu*exp(THOMAS_Pmu*(phi_fines(i,j,k)/(1.-phi_solids(i,j,k))))*sol_eff_mu(i,j,k)
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

	if (Rheological_model.ne.'WINTER') then
		call strain_magnitude(Unew,Vnew,Wnew)				!determine the magnite of the strain rate
	endif


	do i=1,imax						!start loop, in r-direction
		do j=1,jmax						!start loop, in phi-direction
			do k=1,kmax					!start loop, in z-direction
				muA(i,j,k)= muB(i,j,k)+(tauY(i,j,k)/(1.e-12+strain(i,j,k)))*(1.-exp(-PAPANASTASIOUS_m*(1.e-12+strain(i,j,k))))
				ekm(i,j,k)= ekm(i,j,k) + muA(i,j,k)
			enddo
		enddo
	enddo
	call stress_magnitude !calculate the stress magnitude

!!  Boundary conditions Neumann
	call shiftf(ekm,ebf) 
	call shiftb(ekm,ebb) 
	if (periodicy.eq.0.or.periodicy.eq.2) then
		if (rank.eq.0) then ! boundaries in j-direction
			do k=1,kmax
				do i=1,imax
				ekm(i,0,k) = ekm(i,1,k) 
				ekm(i,j1,k)= ebb(i,k) 
				enddo
			enddo
		elseif (rank.eq.px-1) then
			do k=1,kmax
				do i=1,imax
				ekm(i,0,k) = ebf(i,k)
				ekm(i,j1,k)= ekm(i,jmax,k) 
				enddo
			enddo
		else
			do k=1,kmax
				do i=1,imax
				ekm(i,0,k) = ebf(i,k)
				ekm(i,j1,k) =ebb(i,k) 
				enddo
			enddo
		endif
	else
		do k=1,kmax
			do i=1,imax
				ekm(i,0,k) = ebf(i,k)
				ekm(i,j1,k) =ebb(i,k) 
			enddo
		enddo
	endif
	if (periodicx.eq.0.or.periodicx.eq.2) then
		do k=1,kmax ! boundaries in i-direction
			do j=0,j1
				ekm(0,j,k) = ekm(1,j,k)
				ekm(i1,j,k) = ekm(imax,j,k)
			enddo
		enddo
	else
		do k=1,kmax ! boundaries in i-direction
			do j=0,j1
				ekm(0,j,k) = ekm(imax,j,k)
				ekm(i1,j,k) = ekm(1,j,k)
			enddo
		enddo
	endif
	do i=0,i1 ! boundaries in k-direction
		do j=0,j1
			ekm(i,j,0) = ekm(i,j,1)
			ekm(i,j,k1) = ekm(i,j,kmax)
		enddo
	enddo
	END subroutine Bingham_Papanastasiou
