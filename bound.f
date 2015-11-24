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


      subroutine bound(Ubound,Vbound,Wbound,rho,botstress,tt,Ub1in,Vb1in,Wb1in,Ub2in,Vb2in,Wb2in)
     
      USE nlist

      implicit none
c
!       include 'param.txt'
!       include 'common.txt'

c
      integer jtmp,t
c
      real  Ubound(0:i1,0:j1,0:k1),Vbound(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +      Wbound(0:i1,0:j1,0:k1)
      real ubb(0:i1,0:k1),val,theta,Ubc,Vbc,Wbc,xx,yy,r_orifice2,Wjet,theta_U,theta_V
      real vbb(0:i1,0:k1)
      real wbb(0:i1,0:k1)
      real ubf(0:i1,0:k1)
      real vbf(0:i1,0:k1)
      real wbf(0:i1,0:k1)
      real cbf(0:i1,0:k1)
      real cbb(0:i1,0:k1)
	real ust_U_b,ust_V_b,Chezy,fluc,f,tt,z0_U,z0_V,phi
	integer botstress,n
      real Ub1(0:i1,0:k1),Vb1(0:i1,0:k1),Wb1(0:i1,0:k1),Ub2(0:j1,0:k1),Vb2(0:j1+1,0:k1),Wb2(0:j1,0:k1)
      real,intent(in) :: Ub1in(0:i1,0:k1),Vb1in(0:i1,0:k1),Wb1in(0:i1,0:k1),Ub2in(0:j1,0:k1),Vb2in(0:j1+1,0:k1),Wb2in(0:j1,0:k1)
	real x,y,z!,u(1:kmax-2),v(1:kmax-2),w(1:kmax-2),zz(1:kmax-2)
	real Wbc1,Wbc2!,Ubc1,Vbc1
	real Ujetbc,Vjetbc,Wjetbc,Ujetbc2,Vjetbc2
	real rr,interpseries
      real Ujet2,zzz,z2,val2,jetcorr

! 	pi   = 4.0 * atan(1.0)
c
c
c*************************************************************
c
c     Subroutine bound sets the boundary conditions for all variables,
c     except for the diffusion coefficients. These are set in submod.
c     The common boundary conditions for the pressure are set in mkgrid.
c
c*************************************************************
c
c*************************************************************
c	Ubc,Vbc,Wbc are boundary velocities in carthesian
c	x,y,z coordinate system, not in r,theta,z like this code

	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif

c 	influence of waves on lateral boundaries:
	IF(Hs>0.) THEN
	 DO k=0,k1
	  DO i=0,i1
	    IF (rank.eq.0) THEN
		j=0
	    ELSE ! also for rank between 0 and px-1 now Ub1 is filled, but this is not used anyways
		j=jmax 
	    ENDIF
	    z=(k-0.5)*dz-zbed(i,j)
	    z=MIN(z,depth-zbed(i,j))
	    z=MAX(z,0.)
	    z2=k*dz-zbed(i,j)
	    z2=MIN(z2,depth-zbed(i,j))
	    z2=MAX(z2,0.)
	    val=z/(ABS(z)+1e-12)*cosh(kabs_w*z)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.
	    val2=sinh(kabs_w*z2)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.	
	    ! force val to zero when under obstacle height, val2 is automatically
	    x=Ru(i)*cos_u(j)-schuif_x
	    y=Ru(i)*sin_u(j)
	    Ub1(i,k)=Ub1in(i,k)+kx_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_v(j)-schuif_x
	    y=Rp(i)*sin_v(j)
	    Vb1(i,k)=Vb1in(i,k)+ky_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_u(j)-schuif_x
	    y=Rp(i)*sin_u(j)
	    Wb1(i,k)=Wb1in(i,k)+val2*sin(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	  ENDDO
	  DO j=0,j1
	    i=0
	    z=(k-0.5)*dz-zbed(i,j)
	    z=MIN(z,depth-zbed(i,j))
	    z=MAX(z,0.)
	    z2=k*dz-zbed(i,j)
	    z2=MIN(z2,depth-zbed(i,j))
	    z2=MAX(z2,0.)
		val=z/(ABS(z)+1e-12)*cosh(kabs_w*z)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.
		val2=sinh(kabs_w*z2)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.	
	    ! force val to zero when under obstacle height, val2 is automatically
	    x=Ru(i)*cos_u(j)-schuif_x
	    y=Ru(i)*sin_u(j)
	    Ub2(j,k)=Ub2in(j,k)+kx_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_v(j)-schuif_x
	    y=Rp(i)*sin_v(j)
	    Vb2(j,k)=Vb2in(j,k)+ky_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_u(j)-schuif_x
	    y=Rp(i)*sin_u(j)
	    Wb2(j,k)=Wb2in(j,k)+val2*sin(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	  ENDDO	
	 ENDDO
	 Vb2(j1+1,:)=Vb2(j1,:)
	ELSE
	  Ub1=Ub1in
	  Vb1=Vb1in
	  Wb1=Wb1in
	  Ub2=Ub2in
	  Vb2=Vb2in
	  Wb2=Wb2in
	ENDIF
	


c get stuff from other CPU's
	  call shiftf(Ubound,ubf)
	  call shiftf(Vbound,vbf) 
	  call shiftf(Wbound,wbf) 
	  call shiftb(Ubound,ubb) 
	  call shiftb(Vbound,vbb) 
	  call shiftb(Wbound,wbb) 
 
	Wbc1=W_b
	if (periodicy.eq.0) then
	  if (rank.eq.0) then ! boundaries in j-direction
!	    write(*,*),'rank (moet 0 zijn)',rank
	    j=0
	    do k=1,kmax
		do i=1,imax
		    if (bcfile.ne.'') then
			Ubc1(i,k)=Ubcoarse1(i,k)
			Vbc1(i,k)=Vbcoarse1(i,k)
			Wbc1=Wbcoarse1(i,k)
		    endif
		    Ubc=0.5*(Ub1(i,k)+Ub1(i+1,k))+Ubc1(i,k)
		    Vbc=Vb1(i,k)+Vbc1(i,k)
		    Wbc=0.5*(Wb1(i,k)+Wb1(i,k+1))+Wbc1		    
 		   Ubound(i,0,k) = Ubc*cos_u(0)+Vbc*sin_u(0) !Ubf(i,k)
!		   Ubound(i,0,k) = 2.*(Ubc*cos_u(0)+Vbc*sin_u(0))-Ubound(i,1,k) !Ubf(i,k)
		   Vbound(i,0,k) = -Ubc*sin_v(0)+Vbc*cos_v(0) !Vbf(i,k)
!		   Vbound(i,0,k) = 2.*(-Ubc*sin_v(0)+Vbc*cos_v(0))-Vbound(i,1,k) !Vbf(i,k)

		   Wbound(i,0,k) = Wbc !2.*Wbc-Wbound(i,1,k)!Wbc !Wbf(i,k)
		   Ubound(i,j1,k) =Ubb(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	elseif (rank.eq.px-1) then
!		write(*,*),'rank (moet px-1 zijn)',rank
		j=jmax
		do k=1,kmax
		   do i=1,imax
		    if (bcfile.ne.'') then
			Ubc1(i,k)=Ubcoarse1(i,k)
			Vbc1(i,k)=Vbcoarse1(i,k)
			Wbc1=Wbcoarse1(i,k)
		    endif
		    Ubc=0.5*(Ub1(i,k)+Ub1(i+1,k))+Ubc1(i,k)
		    Vbc=Vb1(i,k)+Vbc1(i,k)
		    Wbc=0.5*(Wb1(i,k)+Wb1(i,k+1))+Wbc1	
		   Ubound(i,0,k) = Ubf(i,k)
		   Vbound(i,0,k) = Vbf(i,k)
		   Wbound(i,0,k) = Wbf(i,k)

 		   Ubound(i,j1,k) =Ubc*cos_u(j1)+Vbc*sin_u(j1) !Ubb(i,k)
! 		   Ubound(i,j1,k) =2.*(Ubc*cos_u(j1)+Vbc*sin_u(j1))-Ubound(i,jmax,k) !Ubb(i,k)
		   Vbound(i,jmax,k) =-Ubc*sin_v(jmax)+Vbc*cos_v(jmax) !Vbb(i,k)
!		   Vbound(i,jmax,k) =2.*(-Ubc*sin_v(jmax)+Vbc*cos_v(jmax))-Vbound(i,jmax-1,k) !Vbb(i,k)
		   Wbound(i,j1,k) = Wbc !2.*Wbc-Wbound(i,jmax,k) !Wbb(i,k) 
		   enddo
		enddo	
	  else
		do k=1,kmax
		   do i=1,imax
		     Ubound(i,0,k) = Ubf(i,k)
		     Vbound(i,0,k) = Vbf(i,k)
		     Wbound(i,0,k) = Wbf(i,k)
		     Ubound(i,j1,k) =Ubb(i,k)
		     Vbound(i,j1,k) =Vbb(i,k)
		     Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	  endif
	elseif (periodicy.eq.1) then !periodic in y:
	  do k=1,kmax
	   do i=1,imax
	     Ubound(i,0,k) = Ubf(i,k)
	     Vbound(i,0,k) = Vbf(i,k)
	     Wbound(i,0,k) = Wbf(i,k)
	     Ubound(i,j1,k) =Ubb(i,k)
	     Vbound(i,j1,k) =Vbb(i,k)
	     Wbound(i,j1,k) =Wbb(i,k) 
	   enddo
	  enddo
	elseif (periodicy.eq.2) then !free slip lateral boundaries
	  if (rank.eq.0) then ! boundaries in j-direction
	    j=0
	    do k=1,kmax
		do i=1,imax
 		   Ubound(i,0,k) = Ubound(i,1,k) 
		   Vbound(i,0,k) = 0. 
		   Wbound(i,0,k) = Wbound(i,1,k)
		   Ubound(i,j1,k) =Ubb(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	elseif (rank.eq.px-1) then
!		write(*,*),'rank (moet px-1 zijn)',rank
		j=jmax
		do k=1,kmax
		   do i=1,imax
		   Ubound(i,0,k) = Ubf(i,k)
		   Vbound(i,0,k) = Vbf(i,k)
		   Wbound(i,0,k) = Wbf(i,k)
 		   Ubound(i,j1,k) =Ubound(i,jmax,k)
		   Vbound(i,jmax,k) = 0. 
		   Wbound(i,j1,k) = Wbound(i,jmax,k) 
		   enddo
		enddo	
	  else
		do k=1,kmax
		   do i=1,imax
		     Ubound(i,0,k) = Ubf(i,k)
		     Vbound(i,0,k) = Vbf(i,k)
		     Wbound(i,0,k) = Wbf(i,k)
		     Ubound(i,j1,k) =Ubb(i,k)
		     Vbound(i,j1,k) =Vbb(i,k)
		     Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	  endif
	endif

!	Ubc1=0.
!	Vbc1=0.
	Wbc1=W_b

	if (periodicx.eq.0) then
          do k=1,kmax ! boundaries in i-direction
           do j=0,j1
		    if (bcfile.ne.'') then
			Ubc2(j,k)=Ubcoarse2(j,k)
			Vbc2(j,k)=Vbcoarse2(j,k)
			Wbc2=Wbcoarse2(j,k)
		    endif
		    Ubc=Ub2(j,k)+Ubc2(j,k)
		    Vbc=0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k)
		    Wbc=0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2	
		   Ubound(0,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
		   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
		   Vbound(0,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
		   Vbound(i1,j,k)   =    Vbound(imax,j,k)
		   Wbound(0,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
		   Wbound(i1,j,k)   =    Wbound(imax,j,k)
           enddo   
        enddo
       elseif (periodicx.eq.2) then ! no outflow in x direction:
         do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    0.
		   Ubound(imax,j,k) =    0.
		   Ubound(i1,j,k)   =    0.
		   Vbound(0,j,k)    =    -Vbound(1,j,k)
		   Vbound(i1,j,k)   =    -Vbound(imax,j,k)
		   Wbound(0,j,k)    =    -Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
		   Wbound(i1,j,k)   =    2.*W_ox - Wbound(imax,j,k)
           enddo   
         enddo
       else ! periodic in x:
         do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    Ubound(imax,j,k)
		   Ubound(i1,j,k)   =    Ubound(1,j,k)
		   Vbound(0,j,k)    =    Vbound(imax,j,k)
		   Vbound(i1,j,k)   =    Vbound(1,j,k)
		   Wbound(0,j,k)    =    Wbound(imax,j,k)
		   Wbound(i1,j,k)   =    Wbound(1,j,k)
           enddo   
         enddo
       endif

       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         Ubound(i,j,k1)   = Ubound(i,j,kmax)
         Vbound(i,j,k1)   = Vbound(i,j,kmax)
         Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1)   = 0.
         Ubound(i,j,0)    = Ubound(i,j,1) 
         Vbound(i,j,0)    = Vbound(i,j,1)
         Wbound(i,j,0)    = 0.

         xx=Rp(i)*cos_u(j)-schuif_x
	 yy=Rp(i)*sin_u(j)
	Wbound(i,j,kmax)=(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo

      !! Set boundary conditions vertical jet in:
       	if (plumetseriesfile.eq.'') then
       		Wjet=W_j
       		f=Strouhal*ABS(W_j)/(radius_j*2.) !Strouholt number is 0.3
	else
       		Wjet=interpseries(plumetseries,plumeUseries,plumeseriesloc,tt)
       		f=Strouhal*ABS(W_j)/(radius_j*2.) !Strouholt number is 0.3
	endif
	jetcorr=pi/(2.*pi*(1/(1./W_j_powerlaw+1)-1./(1./W_j_powerlaw+2.))) !=1.22449 for W_j_powerlaw=7
       do t=1,tmax_inPpunt
 	i=i_inPpunt(t)
 	j=j_inPpunt(t)
 	fluc=0.
 	do n=1,azi_n
 	  fluc=fluc+sin(2.*pi*f*tt/n+azi_angle_p(i,j))
 	enddo
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  rr=1.-sqrt((xx**2+yy**2)/radius_j**2)
	  rr=MAX(rr,0.)
 	Wbound(i,j,kmax)=jetcorr*Wjet*rr**(1./W_j_powerlaw)+Awjet/REAL(azi_n)*Wjet*fluc 
! 	Wbound(i,j,k1)=jetcorr*Wjet*rr**(1./W_j_powerlaw)+Awjet/REAL(azi_n)*Wjet*fluc 
       enddo
       do t=1,tmax_inUpunt
 	i=i_inUpunt(t)
 	j=j_inUpunt(t)
 	fluc=0.
 	do n=1,azi_n
 	  fluc=fluc+sin(2.*pi*f*tt/n+azi_angle_u(i,j))
 	enddo
	Ujetbc=Aujet/azi_n*Wjet*fluc
	Vjetbc=Avjet/azi_n*Wjet*fluc
	Ujetbc2=Ujetbc*cos(azi_angle_u(i,j))-Vjetbc*sin(azi_angle_u(i,j))
	Vjetbc2=Ujetbc*sin(azi_angle_u(i,j))+Vjetbc*cos(azi_angle_u(i,j))

  	Ubound(i,j,k1)=2.*(Ujetbc2*cos_u(j)+Vjetbc2*sin_u(j))-Ubound(i,j,kmax) !(Aujet/6.*Wjet*fluc)	
       enddo
       do t=1,tmax_inVpunt
 	i=i_inVpunt(t)
 	j=j_inVpunt(t)
 	fluc=0.
 	do n=1,azi_n
 	  fluc=fluc+sin(2.*pi*f*tt/n+azi_angle_v(i,j))
 	enddo
	Ujetbc=Aujet/azi_n*Wjet*fluc
	Vjetbc=Avjet/azi_n*Wjet*fluc
	Ujetbc2=Ujetbc*cos(azi_angle_v(i,j))-Vjetbc*sin(azi_angle_v(i,j))
	Vjetbc2=Ujetbc*sin(azi_angle_v(i,j))+Vjetbc*cos(azi_angle_v(i,j))

  	Vbound(i,j,k1)=2.*(Ujetbc2*sin_v(j)+Vjetbc2*cos_v(j))-Vbound(i,j,kmax) !(Aujet/6.*Wjet*fluc)
       enddo

      !! Set boundary conditions horizontal jet2 in:
       	if (plumetseriesfile2.eq.'') then
       		Ujet2=U_j2
       		f=Strouhal*ABS(U_j2)/(radius_j2*2.) !Strouholt number is 0.3
	else
       		Ujet2=interpseries(plumetseries2,plumeUseries2,plumeseriesloc2,tt)
       		f=Strouhal2*ABS(U_j2)/(radius_j2*2.) !Strouholt number is 0.3
	endif
       do t=1,tmax_inPpunt2
 	  k=k_inPpunt2(t)
 	  j=j_inPpunt2(t)
	  zzz=k*dz-0.5*dz-zjet2
	  yy=Rp(0)*sin_u(j)
 	  fluc=0.
 	  do n=1,azi_n2
 	    fluc=fluc+sin(2.*pi*f*tt/n+atan2(zzz,yy))
 	  enddo
	    rr=1.-sqrt((zzz**2+yy**2)/radius_j2**2)
	    rr=MAX(rr,0.)
	    Ujetbc=jetcorr*Ujet2*rr**(1./W_j_powerlaw)+Aujet2/azi_n2*Ujet2*fluc 
	    Vjetbc=Avjet/azi_n*Ujet2*fluc
	    Ubound(0,j,k)=(Ujetbc*cos_u(j)+Vjetbc*sin_u(j))
       enddo
       do t=1,tmax_inWpunt2
 	  k=k_inWpunt2(t)
 	  j=j_inWpunt2(t)
	  zzz=k*dz-zjet2
	  yy=Rp(0)*sin_u(j)
 	  fluc=0.
 	  do n=1,azi_n
 	    fluc=fluc+sin(2.*pi*f*tt/n+atan2(zzz,yy))
 	  enddo
	  Wjetbc=Awjet/azi_n*Ujet2*fluc
	  Vjetbc=Avjet/azi_n*Ujet2*fluc
	  Wbound(0,j,k)=2.*(cos(atan2(zzz,yy))*Wjetbc+sin(atan2(zzz,yy))*Vjetbc)-Wbound(1,j,k)
       enddo
       do t=1,tmax_inVpunt2
 	  k=k_inVpunt2(t)
 	  j=j_inVpunt2(t)
	  zzz=k*dz-0.5*dz-zjet2
	  yy=Rp(0)*sin_v(j)
	  rr=1.-sqrt((zzz**2+yy**2)/radius_j**2)
	  rr=MAX(rr,0.)
 	  fluc=0.
 	  do n=1,azi_n
 	    fluc=fluc+sin(2.*pi*f*tt/n+atan2(zzz,yy))
 	  enddo
	  Ujetbc=Ujet2*rr**(1./W_j_powerlaw)+Aujet2/azi_n2*Ujet2*fluc
	  Vjetbc=Avjet/azi_n*Ujet2*fluc
	  Wjetbc=Awjet/azi_n*Ujet2*fluc
	  Vbound(0,j,k)=2.*((-Ujetbc*sin_v(j)+Vjetbc*cos_v(j))*cos(atan2(zzz,yy))-sin(atan2(zzz,yy))*Wjetbc)-Vbound(1,j,k)
       enddo

	
       do t=1,tmax_inWpunt_suction ! add suction in front of draghead
 	  k=k_inWpunt_suction(t)
 	  j=j_inWpunt_suction(t)
 	  i=i_inWpunt_suction(t)
	  Wjetbc=pi*radius_j**2*W_j*perc_dh_suction/(dr(i)*3.*Dsp) !(m3/s)/m2=m/s per suction cell (perc_dh_suction is corrected for number of dragheads)
	  Wbound(i,j,k)=Wjetbc
       enddo

      end

      subroutine bound_incljet(Ubound,Vbound,Wbound,rho,botstress,tt,Ub1in,Vb1in,Wb1in,Ub2in,Vb2in,Wb2in)
     
      USE nlist

      implicit none
c
!       include 'param.txt'
!       include 'common.txt'

c
      integer jtmp,t,inout
	real xTSHD(1:4),yTSHD(1:4)
c
      real  Ubound(0:i1,0:j1,0:k1),Vbound(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +      Wbound(0:i1,0:j1,0:k1)
      real ubb(0:i1,0:k1),val,theta,Ubc,Vbc,Wbc,xx,yy,r_orifice2,Wjet,theta_U,theta_V
      real vbb(0:i1,0:k1)
      real wbb(0:i1,0:k1)
      real ubf(0:i1,0:k1)
      real vbf(0:i1,0:k1)
      real wbf(0:i1,0:k1)
      real cbf(0:i1,0:k1)
      real cbb(0:i1,0:k1)
	real ust_U_b,ust_V_b,Chezy,fluc,f,tt,z0_U,z0_V,phi
	integer botstress,n,jp
      real Ub1(0:i1,0:k1),Vb1(0:i1,0:k1),Wb1(0:i1,0:k1),Ub2(0:j1,0:k1),Vb2(0:j1+1,0:k1),Wb2(0:j1,0:k1)
      real,intent(in):: Ub1in(0:i1,0:k1),Vb1in(0:i1,0:k1),Wb1in(0:i1,0:k1),Ub2in(0:j1,0:k1),Vb2in(0:j1+1,0:k1),Wb2in(0:j1,0:k1)
	real x,y,z
	real Wbc1,Wbc2
	real Ujetbc,Vjetbc,Wjetbc,Ujetbc2,Vjetbc2
	real rr,interpseries
	real uu1,uu2,vv1,vv2
      real zzz,Ujet2,z2,val2,jetcorr

      real  Ubound2(0:i1,0:j1,0:k1),Vbound2(0:i1,0:j1,0:k1),Wbound2(0:i1,0:j1,0:k1) 

! 	pi   = 4.0 * atan(1.0)
c
c
c*************************************************************
c
c     Subroutine bound sets the boundary conditions for all variables,
c     except for the diffusion coefficients. These are set in submod.
c     The common boundary conditions for the pressure are set in mkgrid.
c
c*************************************************************
c
c*************************************************************
c	Ubc,Vbc,Wbc are boundary velocities in carthesian
c	x,y,z coordinate system, not in r,theta,z like this code

	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
	
c 	influence of waves on lateral boundaries:
	IF(Hs>0.) THEN
	 DO k=0,k1
	  DO i=0,i1
	    IF (rank.eq.0) THEN
		j=0
	    ELSE ! also for rank between 0 and px-1 now Ub1 is filled, but this is not used anyways
		j=jmax 
	    ENDIF
	    z=(k-0.5)*dz-zbed(i,j)
	    z=MIN(z,depth-zbed(i,j))
	    z=MAX(z,0.)
	    z2=k*dz-zbed(i,j)
	    z2=MIN(z2,depth-zbed(i,j))
	    z2=MAX(z2,0.)
		val=z/(ABS(z)+1e-12)*cosh(kabs_w*z)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.
		val2=sinh(kabs_w*z2)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.			
	    ! force val to zero when under obstacle height, val2 is automatically
	    x=Ru(i)*cos_u(j)-schuif_x
	    y=Ru(i)*sin_u(j)
	    Ub1(i,k)=Ub1in(i,k)+kx_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_v(j)-schuif_x
	    y=Rp(i)*sin_v(j)
	    Vb1(i,k)=Vb1in(i,k)+ky_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_u(j)-schuif_x
	    y=Rp(i)*sin_u(j)
	    Wb1(i,k)=Wb1in(i,k)+val2*sin(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	  ENDDO
	  DO j=0,j1
	    i=0
	    z=(k-0.5)*dz-zbed(i,j)
	    z=MIN(z,depth-zbed(i,j))
	    z=MAX(z,0.)
	    z2=k*dz-zbed(i,j)
	    z2=MIN(z2,depth-zbed(i,j))
	    z2=MAX(z2,0.)
		val=z/(ABS(z)+1e-12)*cosh(kabs_w*z)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.
	    val2=sinh(kabs_w*z2)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.		
	    ! force val to zero when under obstacle height, val2 is automatically
	    x=Ru(i)*cos_u(j)-schuif_x
	    y=Ru(i)*sin_u(j)
	    Ub2(j,k)=Ub2in(j,k)+kx_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_v(j)-schuif_x
	    y=Rp(i)*sin_v(j)
	    Vb2(j,k)=Vb2in(j,k)+ky_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_u(j)-schuif_x
	    y=Rp(i)*sin_u(j)
	    Wb2(j,k)=Wb2in(j,k)+val2*sin(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	  ENDDO	
	 ENDDO
	 Vb2(j1+1,:)=Vb2(j1,:)
	ELSE
	  Ub1=Ub1in
	  Vb1=Vb1in
	  Wb1=Wb1in
	  Ub2=Ub2in
	  Vb2=Vb2in
	  Wb2=Wb2in
	ENDIF

c get stuff from other CPU's
	  call shiftf(Ubound,ubf)
	  call shiftf(Vbound,vbf) 
	  call shiftf(Wbound,wbf) 
	  call shiftb(Ubound,ubb) 
	  call shiftb(Vbound,vbb) 
	  call shiftb(Wbound,wbb) 

	Wbc1=W_b
	if (periodicy.eq.0) then
	  if (rank.eq.0) then ! boundaries in j-direction
!		write(*,*),'rank (moet 0 zijn)',rank
		j=0
		do k=1,kmax
		   do i=1,imax
		    if (bcfile.ne.'') then
			Ubc1(i,k)=Ubcoarse1(i,k)
			Vbc1(i,k)=Vbcoarse1(i,k)
			Wbc1=Wbcoarse1(i,k)
		    endif
		    Ubc=0.5*(Ub1(i,k)+Ub1(i+1,k))+Ubc1(i,k)
		    Vbc=Vb1(i,k)+Vbc1(i,k)
		    Wbc=0.5*(Wb1(i,k)+Wb1(i,k+1))+Wbc1		    
 		   Ubound(i,0,k) = Ubc*cos_u(0)+Vbc*sin_u(0) !Ubf(i,k)
!		   Ubound(i,0,k) = 2.*(Ubc*cos_u(0)+Vbc*sin_u(0))-Ubound(i,1,k) !Ubf(i,k)
		   Vbound(i,0,k) = -Ubc*sin_v(0)+Vbc*cos_v(0) !Vbf(i,k)

!		   Vbound(i,0,k) = 2.*(-Ubc*sin_v(0)+Vbc*cos_v(0))-Vbound(i,1,k) !Vbf(i,k)
		   Wbound(i,0,k) = Wbc !2.*Wbc-Wbound(i,1,k)!Wbc !Wbf(i,k)
		   Ubound(i,j1,k) =Ubb(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
!		write(*,*),'rank (moet px-1 zijn)',rank
		j=jmax
		do k=1,kmax
		   do i=1,imax
		    if (bcfile.ne.'') then
			Ubc1(i,k)=Ubcoarse1(i,k)
			Vbc1(i,k)=Vbcoarse1(i,k)
			Wbc1=Wbcoarse1(i,k)
		    endif
		    Ubc=0.5*(Ub1(i,k)+Ub1(i+1,k))+Ubc1(i,k)
		    Vbc=Vb1(i,k)+Vbc1(i,k)
		    Wbc=0.5*(Wb1(i,k)+Wb1(i,k+1))+Wbc1	

		   Ubound(i,0,k) = Ubf(i,k)
		   Vbound(i,0,k) = Vbf(i,k)
		   Wbound(i,0,k) = Wbf(i,k)
 		   Ubound(i,j1,k) =Ubc*cos_u(j1)+Vbc*sin_u(j1) !Ubb(i,k)
! 		   Ubound(i,j1,k) =2.*(Ubc*cos_u(j1)+Vbc*sin_u(j1))-Ubound(i,jmax,k) !Ubb(i,k)
		   Vbound(i,jmax,k) =-Ubc*sin_v(jmax)+Vbc*cos_v(jmax) !Vbb(i,k)

!		   Vbound(i,jmax,k) =2.*(-Ubc*sin_v(jmax)+Vbc*cos_v(jmax))-Vbound(i,jmax-1,k) !Vbb(i,k)

		   Wbound(i,j1,k) = Wbc !2.*Wbc-Wbound(i,jmax,k) !Wbb(i,k) 

		   enddo
		enddo	
	  else
!		write(*,*),'rank (moet geen 0 of px zijn)',rank
		do k=1,kmax
		   do i=1,imax
		     Ubound(i,0,k) = Ubf(i,k)
		     Vbound(i,0,k) = Vbf(i,k)
		     Wbound(i,0,k) = Wbf(i,k)
		     Ubound(i,j1,k) =Ubb(i,k)
		     Vbound(i,j1,k) =Vbb(i,k)
		     Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	  endif
	elseif (periodicy.eq.1) then !periodic in y:
	  do k=1,kmax
	   do i=1,imax
	     Ubound(i,0,k) = Ubf(i,k)
	     Vbound(i,0,k) = Vbf(i,k)
	     Wbound(i,0,k) = Wbf(i,k)
	     Ubound(i,j1,k) =Ubb(i,k)
	     Vbound(i,j1,k) =Vbb(i,k)
	     Wbound(i,j1,k) =Wbb(i,k) 
	   enddo
	  enddo
	elseif (periodicy.eq.2) then !free slip lateral boundaries
	  if (rank.eq.0) then ! boundaries in j-direction
	    j=0
	    do k=1,kmax
		do i=1,imax
 		   Ubound(i,0,k) = Ubound(i,1,k) 
		   Vbound(i,0,k) = 0. 
		   Wbound(i,0,k) = Wbound(i,1,k)
		   Ubound(i,j1,k) =Ubb(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	elseif (rank.eq.px-1) then
!		write(*,*),'rank (moet px-1 zijn)',rank
		j=jmax
		do k=1,kmax
		   do i=1,imax
		   Ubound(i,0,k) = Ubf(i,k)
		   Vbound(i,0,k) = Vbf(i,k)
		   Wbound(i,0,k) = Wbf(i,k)
 		   Ubound(i,j1,k) =Ubound(i,jmax,k)
		   Vbound(i,jmax,k) = 0. 
		   Wbound(i,j1,k) = Wbound(i,jmax,k) 
		   enddo
		enddo	
	  else
		do k=1,kmax
		   do i=1,imax
		     Ubound(i,0,k) = Ubf(i,k)
		     Vbound(i,0,k) = Vbf(i,k)
		     Wbound(i,0,k) = Wbf(i,k)
		     Ubound(i,j1,k) =Ubb(i,k)
		     Vbound(i,j1,k) =Vbb(i,k)
		     Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	  endif
	endif

	Wbc1=W_b
	if (periodicx.eq.0) then
	      do k=1,kmax ! boundaries in i-direction
         	do j=0,j1
		    if (bcfile.ne.'') then
			Ubc2(j,k)=Ubcoarse2(j,k)
			Vbc2(j,k)=Vbcoarse2(j,k)
			Wbc2=Wbcoarse2(j,k)
		    endif
		    Ubc=Ub2(j,k)+Ubc2(j,k)
		    Vbc=0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k)
		    Wbc=0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2	
		   Ubound(0,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
		   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
		   Vbound(0,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
		   Vbound(i1,j,k)   =    Vbound(imax,j,k)
		   Wbound(0,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
		   Wbound(i1,j,k)   =    Wbound(imax,j,k)
         	enddo   
      	    enddo
       elseif (periodicx.eq.2) then ! no outflow in x direction:
         do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    0.
		   Ubound(imax,j,k) =    0.
		   Ubound(i1,j,k)   =    0.
		   Vbound(0,j,k)    =    -Vbound(1,j,k)
		   Vbound(i1,j,k)   =    -Vbound(imax,j,k)
		   Wbound(0,j,k)    =    -Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
		   Wbound(i1,j,k)   =    2.*W_ox - Wbound(imax,j,k)
           enddo   
         enddo
       else ! periodic in x:
         do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    Ubound(imax,j,k)
		   Ubound(i1,j,k)   =    Ubound(1,j,k)
		   Vbound(0,j,k)    =    Vbound(imax,j,k)
		   Vbound(i1,j,k)   =    Vbound(1,j,k)
		   Wbound(0,j,k)    =    Wbound(imax,j,k)
		   Wbound(i1,j,k)   =    Wbound(1,j,k)
           enddo   
         enddo
       endif

       !! Get UVW in pipe before wall stresses are taken into account, later these UVW2 are put back in pipe
	Ubound2=Ubound
	Vbound2=Vbound
	Wbound2=Wbound
	
!	if (rank.eq.INT(px/2)) then ! px must be odd, middle piece contains jet in
!	write(*,*),'rank (moet 1 zijn)',rank
       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         Ubound(i,j,k1)   = Ubound(i,j,kmax)
         Vbound(i,j,k1)   = Vbound(i,j,kmax)
         Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1) = 0.
         Ubound(i,j,0)    = Ubound(i,j,1) 
         Vbound(i,j,0)    = Vbound(i,j,1)
         Wbound(i,j,0)    = 0.

         xx=Rp(i)*cos_u(j)-schuif_x
	 yy=Rp(i)*sin_u(j)
	Wbound(i,j,kmax)=(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo

	DO n=1,nbedplume
	IF (bp(n)%u.ne.-99999.) THEN ! apply bedplume velocity boundary condition:
	IF ((bp(n)%forever.eq.1.and.time_np.gt.bp(n)%t0)
     &     .or.(bp(n)%forever.eq.0.and.time_n.lt.bp(n)%t0.and.time_np.gt.bp(n)%t0)) THEN
	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
      !! Search for P,V:
      do k=k1,0,-1 !from top to bottom
       do i=0,i1  
         do j=j1,0,-1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  IF (k.le.FLOOR(bp(n)%height/dz).and.k.ge.CEILING(bp(n)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n)%x*cos(phi)-bp(n)%y*sin(phi)
		yTSHD(1:4)=bp(n)%x*sin(phi)+bp(n)%y*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
		Ubound(i,j,k)=bp(n)%u*cos_u(j)+bp(n)%v*sin_u(j)
		Vbound(i,j,k)=-bp(n)%u*sin_v(j)+bp(n)%v*cos_v(j)
		Wbound(i,j,k)=bp(n)%w
		Wbound(i,j,MAX(k-1,0))=bp(n)%w
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDIF
	ENDDO ! bedplume loop
	
	 
	IF (LOA>0.) THEN ! ship:
	  do t=1,tmax_inWpuntTSHD
 	    i=i_inWpuntTSHD(t)
 	    j=j_inWpuntTSHD(t)		
 	    k=k_inWpuntTSHD(t)		
	    Wbound(i,j,k)=0.
	  enddo
	  do t=1,tmax_inUpuntTSHD
 	    i=i_inUpuntTSHD(t)
 	    j=j_inUpuntTSHD(t)		
 	    k=k_inUpuntTSHD(t)		
	    Ubound(i,j,k)=0.
	  enddo
	  do t=1,tmax_inVpuntTSHD
 	    i=i_inVpuntTSHD(t)
 	    j=j_inVpuntTSHD(t)		
 	    k=k_inVpuntTSHD(t)		
	    Vbound(i,j,k)=0.
	  enddo
	  do t=1,tmax_inVpunt_rudder ! apply rudder 
 	    i=i_inVpunt_rudder(t)
 	    j=j_inVpunt_rudder(t)		
 	    k=k_inVpunt_rudder(t)
	    jp=MIN(j+1,j1) ! if j is j1, then jp=j1 --> uu2 is incorrect, but this is no problem as Ubound(i,j1,k) is not used	
	    uu1 = -sin_v(j)*Vbound(i,j,k)+(Ubound(i,j,k))*cos_v(j) ! u wanted along rudder (v wanted = zero)
	    uu2 = -sin_v(j)*Vbound(i,j,k)+(Ubound(i,jp,k))*cos_v(j) ! u wanted along rudder (v wanted = zero)
	    Vbound(i,j,k)= -sin_v(j)*0.5*(uu1+uu2)
	    Ubound(i,j,k)= cos_u(j)*uu1
	    Ubound(i,jp,k)= cos_u(jp)*uu2
	  enddo
	  do i=0,i1 ! prescribe UTSHD in obstacle:
	    do j=0,j1
	      do k=1,kbed(i,j)
	        Ubound(i,j,k)=Ubot_TSHD(j)
	        Vbound(i,j,k)=Vbot_TSHD(j)
	      enddo
	    enddo
	  enddo
	ELSE 
          !! Set boundary conditions plate:
	  if (kjet>0) then
	  Ubound(0:i1,0:j1,kmax-kjet+1:k1)=0. ! top layer (kjet) is NO slip boundary, let develop freely no slip is obtained through wall_fun
	  Vbound(0:i1,0:j1,kmax-kjet+1:k1)=0. ! top layer (kjet) is NO slip boundary, let develop freely no slip is obtained through wall_fun
	  Wbound(0:i1,0:j1,kmax-kjet:kmax)=0.
	  endif
	ENDIF
	! set boundary conditions vertical jet:
       	if (plumetseriesfile.eq.'') then
       		Wjet=W_j
       		f=Strouhal*ABS(W_j)/(radius_j*2.) !Strouholt number is 0.3
	else
       		Wjet=interpseries(plumetseries,plumeUseries,plumeseriesloc,tt)
       		f=Strouhal*ABS(W_j)/(radius_j*2.) !Strouholt number is 0.3
	endif
	jetcorr=pi/(2.*pi*(1/(1./W_j_powerlaw+1)-1./(1./W_j_powerlaw+2.))) !=1.22449 for W_j_powerlaw=7
       do t=1,tmax_inPpunt
 	i=i_inPpunt(t)
 	j=j_inPpunt(t)
 	fluc=0.
 	do n=1,azi_n
 	  fluc=fluc+sin(2.*pi*f*tt/n+azi_angle_p(i,j))
 	enddo
	xx=Rp(i)*cos_u(j)-schuif_x
	yy=Rp(i)*sin_u(j)
	rr=1.-sqrt((xx**2+yy**2)/radius_j**2)
	rr=MAX(rr,0.)
	if (outflow_overflow_down.eq.1) then
	 	do k=kmax-kjet,kmax-kjet
	 	  Wbound(i,j,k)=Wbound2(i,j,k)	
	 	enddo
	 	do k=kmax-kjet+1,kmax
	 	  Wbound(i,j,k)=(jetcorr*rr**(1./W_j_powerlaw)+Awjet/REAL(azi_n)*fluc)*Wjet
	  	enddo
	else
	 	do k=kmax-kjet,kmax
	 	  Wbound(i,j,k)=Wbound2(i,j,k) !Wjet+Awjet/REAL(azi_n)*Wjet*fluc
	 	enddo
	 	Wbound(i,j,kmax)=(jetcorr*rr**(1./W_j_powerlaw)+Awjet/REAL(azi_n)*fluc)*Wjet
	! 	Wbound(i,j,k1)=jetcorr*Wjet*rr**(1./W_j_powerlaw)+Awjet/REAL(azi_n)*Wjet*fluc 
	endif
       enddo
       do t=1,tmax_inUpunt
 	i=i_inUpunt(t)
 	j=j_inUpunt(t)
 	fluc=0.
 	do n=1,azi_n
 	  fluc=fluc+sin(2.*pi*f*tt/n+azi_angle_u(i,j))
 	enddo
	Ujetbc=Aujet/azi_n*Wjet*fluc
	Vjetbc=Avjet/azi_n*Wjet*fluc
	Ujetbc2=Ujetbc*cos(azi_angle_u(i,j))-Vjetbc*sin(azi_angle_u(i,j))
	Vjetbc2=Ujetbc*sin(azi_angle_u(i,j))+Vjetbc*cos(azi_angle_u(i,j))
	if (outflow_overflow_down.eq.1) then
		do k=kmax-kjet,kmax-kjet+1
		  Ubound(i,j,k)=Ubound2(i,j,k)
		enddo	
		do k=kmax-kjet+2,k1
		  Ubound(i,j,k)=0.
		enddo	
	else
		do k=kmax-kjet,k1
		  Ubound(i,j,k)=Ubound2(i,j,k)
		enddo	
	  	Ubound(i,j,k1)=2.*(Ujetbc2*cos_u(j)+Vjetbc2*sin_u(j))-Ubound(i,j,kmax) !(Aujet/6.*Wjet*fluc)	
	endif
	
       enddo
       do t=1,tmax_inVpunt
 	i=i_inVpunt(t)
 	j=j_inVpunt(t)
 	fluc=0.
 	do n=1,azi_n
 	  fluc=fluc+sin(2.*pi*f*tt/n+azi_angle_v(i,j))
 	enddo
	Ujetbc=Aujet/azi_n*Wjet*fluc
	Vjetbc=Avjet/azi_n*Wjet*fluc
	Ujetbc2=Ujetbc*cos(azi_angle_v(i,j))-Vjetbc*sin(azi_angle_v(i,j))
	Vjetbc2=Ujetbc*sin(azi_angle_v(i,j))+Vjetbc*cos(azi_angle_v(i,j))
	if (outflow_overflow_down.eq.1) then
		do k=kmax-kjet,kmax-kjet+1
		  Vbound(i,j,k)=Vbound2(i,j,k)
		enddo
		do k=kmax-kjet+2,k1
		  Vbound(i,j,k)=0.
		enddo
	else
		do k=kmax-kjet,k1
		  Vbound(i,j,k)=Vbound2(i,j,k)
		enddo
	  	Vbound(i,j,k1)=2.*(Ujetbc2*sin_v(j)+Vjetbc2*cos_v(j))-Vbound(i,j,kmax) !(Aujet/6.*Wjet*fluc)
	endif

       enddo
      !! Set boundary conditions horizontal jet2 in:
       	if (plumetseriesfile2.eq.'') then
       		Ujet2=U_j2
       		f=Strouhal*ABS(U_j2)/(radius_j2*2.) !Strouholt number is 0.3
	else
       		Ujet2=interpseries(plumetseries2,plumeUseries2,plumeseriesloc2,tt)
       		f=Strouhal2*ABS(U_j2)/(radius_j2*2.) !Strouholt number is 0.3
	endif
       do t=1,tmax_inPpunt2
 	  k=k_inPpunt2(t)
 	  j=j_inPpunt2(t)
	  zzz=k*dz-0.5*dz-zjet2
	  yy=Rp(0)*sin_u(j)
 	  fluc=0.
 	  do n=1,azi_n2
 	    fluc=fluc+sin(2.*pi*f*tt/n+atan2(zzz,yy))
 	  enddo
	    rr=1.-sqrt((zzz**2+yy**2)/radius_j2**2)
	    rr=MAX(rr,0.)
	    Ujetbc=jetcorr*Ujet2*rr**(1./W_j_powerlaw)+Aujet2/azi_n2*Ujet2*fluc 
	    Vjetbc=Avjet/azi_n*Ujet2*fluc
	    Ubound(0,j,k)=(Ujetbc*cos_u(j)+Vjetbc*sin_u(j))
       enddo
       do t=1,tmax_inWpunt2
 	  k=k_inWpunt2(t)
 	  j=j_inWpunt2(t)
	  zzz=k*dz-zjet2
	  yy=Rp(0)*sin_u(j)
 	  fluc=0.
 	  do n=1,azi_n
 	    fluc=fluc+sin(2.*pi*f*tt/n+atan2(zzz,yy))
 	  enddo
	  Wjetbc=Awjet/azi_n*Ujet2*fluc
	  Vjetbc=Avjet/azi_n*Ujet2*fluc
	  Wbound(0,j,k)=2.*(cos(atan2(zzz,yy))*Wjetbc+sin(atan2(zzz,yy))*Vjetbc)-Wbound(1,j,k)
       enddo
       do t=1,tmax_inVpunt2
 	  k=k_inVpunt2(t)
 	  j=j_inVpunt2(t)
	  zzz=k*dz-0.5*dz-zjet2
	  yy=Rp(0)*sin_v(j)
	  rr=1.-sqrt((zzz**2+yy**2)/radius_j**2)
	  rr=MAX(rr,0.)
 	  fluc=0.
 	  do n=1,azi_n
 	    fluc=fluc+sin(2.*pi*f*tt/n+atan2(zzz,yy))
 	  enddo
	  Ujetbc=Ujet2*rr**(1./W_j_powerlaw)+Aujet2/azi_n2*Ujet2*fluc
	  Vjetbc=Avjet/azi_n*Ujet2*fluc
	  Wjetbc=Awjet/azi_n*Ujet2*fluc
	  Vbound(0,j,k)=2.*((-Ujetbc*sin_v(j)+Vjetbc*cos_v(j))*cos(atan2(zzz,yy))-sin(atan2(zzz,yy))*Wjetbc)-Vbound(1,j,k)
       enddo

       do t=1,tmax_inWpunt_suction ! add suction in front of draghead
 	  k=k_inWpunt_suction(t)
 	  j=j_inWpunt_suction(t)
 	  i=i_inWpunt_suction(t)
	  Wjetbc=pi*radius_j**2*W_j*perc_dh_suction/(dr(i)*3.*Dsp) !(m3/s)/m2=m/s per suction cell (perc_dh_suction is corrected for number of dragheads)
	  Wbound(i,j,k)=Wjetbc
       enddo

	! apply boundary conditions for Vbound again to fix small inconsistency in Vbound(:,j1,:) made in rudder
	call shiftf(Vbound,vbf) 
	call shiftb(Vbound,vbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then ! boundaries in j-direction

	  elseif (rank.eq.px-1) then
	
	  else
		do k=0,k1
		   do i=1,imax
		   Vbound(i,0,k) = Vbf(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   enddo
		enddo
	  endif
	else !periodic in y:
	  do k=0,k1
	   do i=1,imax
	     Vbound(i,0,k) = Vbf(i,k)
	     Vbound(i,j1,k) =Vbb(i,k)
	   enddo
	  enddo
	endif


      end

      subroutine bound_rhoU(Ubound,Vbound,Wbound,rho,botstress,tt,Ub1in,Vb1in,Wb1in,Ub2in,Vb2in,Wb2in)

      USE nlist

      implicit none
c
!       include 'param.txt'
!       include 'common.txt'

c
      integer jtmp,botstress,n,t,jp,inout
	real xTSHD(1:4),yTSHD(1:4)
c
      real  Ubound(0:i1,0:j1,0:k1),Vbound(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +      Wbound(0:i1,0:j1,0:k1)
      real ubb(0:i1,0:k1),val,theta,Ubc,Vbc,Wbc,xx,yy,Ujet,Vjet,Wjet,theta_U,theta_V
      real vbb(0:i1,0:k1)
      real wbb(0:i1,0:k1)
      real ubf(0:i1,0:k1)
      real vbf(0:i1,0:k1)
      real wbf(0:i1,0:k1)
      real cbf(0:i1,0:k1)
      real cbb(0:i1,0:k1)
	real ust_U_b,ust_V_b,Chezy,fluc,f,tt,z0_U,z0_V,phi,z,x,y
      real Ub1(0:i1,0:k1),Vb1(0:i1,0:k1),Wb1(0:i1,0:k1),Ub2(0:j1,0:k1),Vb2(0:j1+1,0:k1),Wb2(0:j1,0:k1)
      real,intent(in):: Ub1in(0:i1,0:k1),Vb1in(0:i1,0:k1),Wb1in(0:i1,0:k1),Ub2in(0:j1,0:k1),Vb2in(0:j1+1,0:k1),Wb2in(0:j1,0:k1)
	real Wbc1,Wbc2
	real Ujetbc,Vjetbc,Wjetbc,Ujetbc2,Vjetbc2
	real rr,interpseries
	real uu1,uu2,vv1,vv2,test
      real zzz,Ujet2,z2,val2,jetcorr

!      real  Ubound2(0:i1,0:j1,0:k1),Vbound2(0:i1,0:j1,0:k1),Wbound2(0:i1,0:j1,0:k1)
      real  Ubound2(-1:i1+1,-1:j1+1,0:k1),Vbound2(-1:i1+1,-1:j1+1,0:k1),Wbound2(0:i1,0:j1,0:k1)
c
c
c*************************************************************
c
c     Subroutine bound sets the boundary conditions for all variables,
c     except for the diffusion coefficients. These are set in submod.
c     The common boundary conditions for the pressure are set in mkgrid.
c
c*************************************************************
c
c*************************************************************
c	Ubc,Vbc,Wbc are boundary velocities in carthesian
c	x,y,z coordinate system, not in r,theta,z like this code

	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif

c 	influence of waves on lateral boundaries:
	IF(Hs>0.) THEN
	 DO k=0,k1
	  DO i=0,i1
	    IF (rank.eq.0) THEN
		j=0
	    ELSE ! also for rank between 0 and px-1 now Ub1 is filled, but this is not used anyways
		j=jmax 
	    ENDIF
	    z=(k-0.5)*dz-zbed(i,j)
	    z=MIN(z,depth-zbed(i,j))
	    z=MAX(z,0.)
	    z2=k*dz-zbed(i,j)
	    z2=MIN(z2,depth-zbed(i,j))
	    z2=MAX(z2,0.)
		val=z/(ABS(z)+1e-12)*cosh(kabs_w*z)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.
	    val2=sinh(kabs_w*z2)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.	
	    ! force val to zero when under obstacle height, val2 is automatically
	    x=Ru(i)*cos_u(j)-schuif_x
	    y=Ru(i)*sin_u(j)
	    Ub1(i,k)=Ub1in(i,k)+kx_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_v(j)-schuif_x
	    y=Rp(i)*sin_v(j)
	    Vb1(i,k)=Vb1in(i,k)+ky_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_u(j)-schuif_x
	    y=Rp(i)*sin_u(j)
	    Wb1(i,k)=Wb1in(i,k)+val2*sin(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	  ENDDO
	  DO j=0,j1
	    i=0
	    z=(k-0.5)*dz-zbed(i,j)
	    z=MIN(z,depth-zbed(i,j))
	    z=MAX(z,0.)
	    z2=k*dz-zbed(i,j)
	    z2=MIN(z2,depth-zbed(i,j))
	    z2=MAX(z2,0.)
		val=z/(ABS(z)+1e-12)*cosh(kabs_w*z)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.
	    val2=sinh(kabs_w*z2)/sinh(kabs_w*(depth-zbed(i,j)))*(om_w-kx_w*U_w-ky_w*V_w)*Hs/2.	
	    ! force val to zero when under obstacle height, val2 is automatically
	    x=Ru(i)*cos_u(j)-schuif_x
	    y=Ru(i)*sin_u(j)
	    Ub2(j,k)=Ub2in(j,k)+kx_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_v(j)-schuif_x
	    y=Rp(i)*sin_v(j)
	    Vb2(j,k)=Vb2in(j,k)+ky_w/kabs_w*val*cos(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	    x=Rp(i)*cos_u(j)-schuif_x
	    y=Rp(i)*sin_u(j)
	    Wb2(j,k)=Wb2in(j,k)+val2*sin(kx_w*(x-U_TSHD*cos(phi)*tt)+ky_w*(y-U_TSHD*sin(phi)*tt)-om_w*tt)
	  ENDDO	
	 ENDDO
	 Vb2(j1+1,:)=Vb2(j1,:)
	ELSE
	  Ub1=Ub1in
	  Vb1=Vb1in
	  Wb1=Wb1in
	  Ub2=Ub2in
	  Vb2=Vb2in
	  Wb2=Wb2in
	ENDIF

c get stuff from other CPU's
	call shiftf(Ubound,ubf)
	call shiftf(Vbound,vbf) 
	call shiftf(Wbound,wbf) 
	call shiftb(Ubound,ubb) 
	call shiftb(Vbound,vbb) 
	call shiftb(Wbound,wbb) 

	Wbc1=W_b

	if (periodicy.eq.0) then
	  if (rank.eq.0) then ! boundaries in j-direction
		do k=1,kmax
		   do i=1,imax
		    if (bcfile.ne.'') then
			Ubc1(i,k)=Ubcoarse1(i,k)
			Vbc1(i,k)=Vbcoarse1(i,k)
			Wbc1=Wbcoarse1(i,k)
		    endif
		    Ubc=0.5*(rho(i,0,k)+rho(i,1,k))*(0.5*(Ub1(i,k)+Ub1(i+1,k))+Ubc1(i,k))
		    Vbc=0.5*(rho(i,0,k)+rho(i,1,k))*(Vb1(i,k)+Vbc1(i,k))
		    Wbc=0.5*(rho(i,0,k)+rho(i,1,k))*(0.5*(Wb1(i,k)+Wb1(i,k+1))+Wbc1)

 		   Ubound(i,0,k) = Ubc*cos_u(0)+Vbc*sin_u(0) !Ubf(i,k)
!		   Ubound(i,0,k) = 2.*(Ubc*cos_u(0)+Vbc*sin_u(0))-Ubound(i,1,k) !Ubf(i,k)
		   Vbound(i,0,k) = -Ubc*sin_v(0)+Vbc*cos_v(0) !Vbf(i,k)

!		   Vbound(i,0,k) = 2.*(-Ubc*sin_v(0)+Vbc*cos_v(0))-Vbound(i,1,k) !Vbf(i,k)
		   Wbound(i,0,k) = Wbc !2.*Wbc-Wbound(i,1,k)!Wbc !Wbf(i,k)
		   
		   Ubound(i,j1,k) =Ubb(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,kmax
		   do i=1,imax
		    if (bcfile.ne.'') then
			Ubc1(i,k)=Ubcoarse1(i,k)
			Vbc1(i,k)=Vbcoarse1(i,k)
			Wbc1=Wbcoarse1(i,k)
		    endif
		    Ubc=0.5*(rho(i,jmax,k)+rho(i,j1,k))*(0.5*(Ub1(i,k)+Ub1(i+1,k))+Ubc1(i,k))
		    Vbc=0.5*(rho(i,jmax,k)+rho(i,j1,k))*(Vb1(i,k)+Vbc1(i,k))
		    Wbc=0.5*(rho(i,jmax,k)+rho(i,j1,k))*(0.5*(Wb1(i,k)+Wb1(i,k+1))+Wbc1)

!		    Ubc=rho_b*(0.5*(Ub1(i,k)+Ub1(i+1,k))+Ubc1)
!		    Vbc=rho_b*(Vb1(i,k)+Vbc1)
!		    Wbc=rho_b*(0.5*(Wb1(i,k)+Wb1(i,k+1))+Wbc1)
		   Ubound(i,0,k) = Ubf(i,k)
		   Vbound(i,0,k) = Vbf(i,k)
		   Wbound(i,0,k) = Wbf(i,k)

 		   Ubound(i,j1,k) =Ubc*cos_u(j1)+Vbc*sin_u(j1) !Ubb(i,k)
! 		   Ubound(i,j1,k) =2.*(Ubc*cos_u(j1)+Vbc*sin_u(j1))-Ubound(i,jmax,k) !Ubb(i,k)
		   Vbound(i,jmax,k) =-Ubc*sin_v(jmax)+Vbc*cos_v(jmax) !Vbb(i,k)

!		   Vbound(i,jmax,k) =2.*(-Ubc*sin_v(jmax)+Vbc*cos_v(jmax))-Vbound(i,jmax-1,k) !Vbb(i,k)
		   Wbound(i,j1,k) =Wbc !2.*Wbc-Wbound(i,jmax,k) !Wbb(i,k) 

		   enddo
		enddo	
	  else
		do k=1,kmax
		   do i=1,imax
		   Ubound(i,0,k) = Ubf(i,k)
		   Vbound(i,0,k) = Vbf(i,k)
		   Wbound(i,0,k) = Wbf(i,k)
		   Ubound(i,j1,k) =Ubb(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	  endif
	elseif (periodicy.eq.1) then !periodic in y:
	  do k=1,kmax
	   do i=1,imax
	     Ubound(i,0,k) = Ubf(i,k)
	     Vbound(i,0,k) = Vbf(i,k)
	     Wbound(i,0,k) = Wbf(i,k)
	     Ubound(i,j1,k) =Ubb(i,k)
	     Vbound(i,j1,k) =Vbb(i,k)
	     Wbound(i,j1,k) =Wbb(i,k) 
	   enddo
	  enddo
	elseif (periodicy.eq.2) then !free slip lateral boundaries
	  if (rank.eq.0) then ! boundaries in j-direction
	    j=0
	    do k=1,kmax
		do i=1,imax
 		   Ubound(i,0,k) = Ubound(i,1,k) 
		   Vbound(i,0,k) = 0. 
		   Wbound(i,0,k) = Wbound(i,1,k)
		   Ubound(i,j1,k) =Ubb(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	elseif (rank.eq.px-1) then
!		write(*,*),'rank (moet px-1 zijn)',rank
		j=jmax
		do k=1,kmax
		   do i=1,imax
		   Ubound(i,0,k) = Ubf(i,k)
		   Vbound(i,0,k) = Vbf(i,k)
		   Wbound(i,0,k) = Wbf(i,k)
 		   Ubound(i,j1,k) =Ubound(i,jmax,k)
		   Vbound(i,jmax,k) = 0. 
		   Wbound(i,j1,k) = Wbound(i,jmax,k) 
		   enddo
		enddo	
	  else
		do k=1,kmax
		   do i=1,imax
		     Ubound(i,0,k) = Ubf(i,k)
		     Vbound(i,0,k) = Vbf(i,k)
		     Wbound(i,0,k) = Wbf(i,k)
		     Ubound(i,j1,k) =Ubb(i,k)
		     Vbound(i,j1,k) =Vbb(i,k)
		     Wbound(i,j1,k) =Wbb(i,k) 
		   enddo
		enddo
	  endif
	endif

	Wbc1=W_b
	if (periodicx.eq.0) then	
      	   do k=1,kmax ! boundaries in i-direction
	   	do j=0,j1
		    if (bcfile.ne.'') then
			Ubc2(j,k)=Ubcoarse2(j,k)
			Vbc2(j,k)=Vbcoarse2(j,k)
			Wbc2=Wbcoarse2(j,k)
		    endif
		    Ubc=0.5*(rho(0,j,k)+rho(1,j,k))*(Ub2(j,k)+Ubc2(j,k))
		    Vbc=0.5*(rho(0,j,k)+rho(1,j,k))*(0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k))
		    Wbc=0.5*(rho(0,j,k)+rho(1,j,k))*(0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2)	
		   Ubound(0,j,k)    =    (Ubc*cos_u(j)+Vbc*sin_u(j))
		   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
		   Vbound(0,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j)) !-Vbound(1,j,k)
		   Vbound(i1,j,k)   =    Vbound(imax,j,k)
		   Wbound(0,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
		   Wbound(i1,j,k)   =    Wbound(imax,j,k)
         	enddo   
      	    enddo
       elseif (periodicx.eq.2) then ! no outflow in x direction:
         do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    0.
		   Ubound(imax,j,k) =    0.
		   Ubound(i1,j,k)   =    0.
		   Vbound(0,j,k)    =    -Vbound(1,j,k)
		   Vbound(i1,j,k)   =    -Vbound(imax,j,k)
		   Wbound(0,j,k)    =    -Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
		   Wbound(i1,j,k)   =    2.*W_ox*rho_b - Wbound(imax,j,k)
           enddo   
         enddo
       else ! periodic in x:
         do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    Ubound(imax,j,k)
		   Ubound(i1,j,k)   =    Ubound(1,j,k)
		   Vbound(0,j,k)    =    Vbound(imax,j,k)
		   Vbound(i1,j,k)   =    Vbound(1,j,k)
		   Wbound(0,j,k)    =    Wbound(imax,j,k)
		   Wbound(i1,j,k)   =    Wbound(1,j,k)
           enddo   
         enddo
       endif

       !! Get UVW in pipe before wall stresses are taken into account, later these UVW2 are put back in pipe
	Ubound2(0:i1,0:j1,0:k1)=Ubound
	Vbound2(0:i1,0:j1,0:k1)=Vbound
	Wbound2(0:i1,0:j1,0:k1)=Wbound
	! fill Ubound2,Vbound2 with one extra row positive and negative in i,j direction for shear stress 
	if (periodicx.eq.0) then
		Ubound2(-1,0:j1,0:k1)=Ubound(0,0:j1,0:k1)
		Ubound2(i1+1,0:j1,0:k1)=Ubound(imax,0:j1,0:k1)
		Vbound2(-1,0:j1,0:k1)=Vbound(0,0:j1,0:k1)
		Vbound2(i1+1,0:j1,0:k1)=Vbound(imax,0:j1,0:k1)
	else 
		Ubound2(-1,0:j1,0:k1)=Ubound(imax-1,0:j1,0:k1)
		Ubound2(i1+1,0:j1,0:k1)=Ubound(1,0:j1,0:k1)
		Vbound2(-1,0:j1,0:k1)=Vbound(imax-1,0:j1,0:k1)
		Vbound2(i1+1,0:j1,0:k1)=Vbound(1,0:j1,0:k1)
	endif

c get stuff from other CPU's
	  call shiftf2(Ubound,ubf)
	  call shiftb2(Ubound,ubb) 
	  call shiftf2(Vbound,vbf)
	  call shiftb2(Vbound,vbb)
	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1
		   Ubound2(i,-1,k) = Ubound(i,0,k)
		   Ubound2(i,j1+1,k) =Ubb(i,k)
		   Vbound2(i,-1,k) = Vbound(i,0,k)
		   Vbound2(i,j1+1,k) =Vbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1
		   Ubound2(i,-1,k) = Ubf(i,k)
		   Ubound2(i,j1+1,k) =Ubound(i,j1,k)
		   Vbound2(i,-1,k) = Vbf(i,k)
		   Vbound2(i,j1+1,k) =Vbound(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1
		   Ubound2(i,-1,k) = Ubf(i,k)
		   Ubound2(i,j1+1,k) =Ubb(i,k)
		   Vbound2(i,-1,k) = Vbf(i,k)
		   Vbound2(i,j1+1,k) =Vbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=0,k1
		   do i=0,i1
		   Ubound2(i,-1,k) = Ubf(i,k)
		   Ubound2(i,j1+1,k) =Ubb(i,k)
		   Vbound2(i,-1,k) = Vbf(i,k)
		   Vbound2(i,j1+1,k) =Vbb(i,k)
		   enddo
		enddo
	endif


       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         if (botstress.ge.1) then  !((botstress.eq.1.and.LOA>0.).or.(botstress.eq.1.and.periodicx.eq.1)) then ! with ship or periodic apply shear stress at seabed:
			!! First Ubot_TSHD and Vbot_TSHD is subtracted to determine tau 
			!! only over ambient velocities not over U_TSHD
			!! Then Ubound,Vbound are reduced by tau, 
			!! now Ubot_TSHD and Vbot_TSHD are added again to get correct total velocity
			uu1=Ubound(i,j,kbedt(i,j)+1)-Ubot_TSHD(j)*rho(i,j,kbedt(i,j)+1)
			vv1=0.25*(Vbound2(i,j,kbedt(i,j)+1)+Vbound2(i+1,j,kbedt(i,j)+1)+Vbound2(i,j-1,kbedt(i,j)+1)
     &                      +Vbound2(i+1,j-1,kbedt(i,j)+1))-Vbot_TSHD(j)*rho(i,j,kbedt(i,j)+1)
			uu2=0.25*(Ubound2(i,j,kbedt(i,j)+1)+Ubound2(i,j+1,kbedt(i,j)+1)+Ubound2(i-1,j,kbedt(i,j)+1)
     &                      +Ubound2(i-1,j+1,kbedt(i,j)+1))-Ubot_TSHD(j)*rho(i,j,kbedt(i,j)+1)
			vv2=Vbound(i,j,kbedt(i,j)+1)-Vbot_TSHD(j)*rho(i,j,kbedt(i,j)+1)
			call wall_fun_rho(uu1,vv1
     &  ,rho(i,j,1),dz,dt,kn,kappa,(depth-bc_obst_h),U_b,nu_mol)
			call wall_fun_rho(vv2,uu2
     &  ,rho(i,j,1),dz,dt,kn,kappa,(depth-bc_obst_h),U_b,nu_mol)
			Ubound(i,j,kbedt(i,j)+1)=uu1+Ubot_TSHD(j)*rho(i,j,kbedt(i,j)+1)
			Vbound(i,j,kbedt(i,j)+1)=vv2+Vbot_TSHD(j)*rho(i,j,kbedt(i,j)+1)
	   elseif (botstress.ge.1.and.kjet>0.and.LOA<0.) then ! with flat plate shear stress must be applied:
!	    call wall_fun_rho(Vbound(i,j,kmax-kjet-1),Ubound(i,j,kmax-kjet-1),rho(i,j,kmax-kjet-1),dz,dt,kn,kappa,(depth-bc_obst_h),U_b,nu_mol)
!	    call wall_fun_rho(Ubound(i,j,kmax-kjet-1),Vbound(i,j,kmax-kjet-1),rho(i,j,kmax-kjet-1),dz,dt,kn,kappa,(depth-bc_obst_h),U_b,nu_mol)
		uu2=0.25*(Ubound2(i,j,kmax-kjet-1)+Ubound2(i,j+1,kmax-kjet-1)+Ubound2(i-1,j,kmax-kjet-1)+Ubound2(i-1,j+1,kmax-kjet-1))
		vv2=0.25*(Vbound2(i,j,kmax-kjet-1)+Vbound2(i+1,j,kmax-kjet-1)+Vbound2(i,j-1,kmax-kjet-1)+Vbound2(i+1,j-1,kmax-kjet-1))

	    call wall_fun_rho(Vbound(i,j,kmax-kjet-1),uu2,rho(i,j,kmax-kjet-1),dz,dt,kn,kappa,(depth-bc_obst_h),U_b,nu_mol)
	    call wall_fun_rho(Ubound(i,j,kmax-kjet-1),vv2,rho(i,j,kmax-kjet-1),dz,dt,kn,kappa,(depth-bc_obst_h),U_b,nu_mol)

	   endif
	
         Ubound(i,j,k1)   = Ubound(i,j,kmax)
         Vbound(i,j,k1)   = Vbound(i,j,kmax)
         Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1)   = 0.
         Ubound(i,j,0)    = Ubound(i,j,1) 
         Vbound(i,j,0)    = Vbound(i,j,1)
         Wbound(i,j,0)    = 0.

         xx=Rp(i)*cos_u(j)-schuif_x
	 yy=Rp(i)*sin_u(j)
	! apply vertical boundary condition belonging to waves:
	Wbound(i,j,kmax)=rho_b*(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)

        enddo
       enddo

	DO n=1,nbedplume
	IF (bp(n)%u.ne.-99999.) THEN ! apply bedplume velocity boundary condition:
	IF ((bp(n)%forever.eq.1.and.time_np.gt.bp(n)%t0)
     &     .or.(bp(n)%forever.eq.0.and.time_n.lt.bp(n)%t0.and.time_np.gt.bp(n)%t0)) THEN
	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
      !! Search for P,V:
      do k=k1,0,-1 !from top to bottom
       do i=0,i1  
         do j=j1,0,-1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  IF (k.le.FLOOR(bp(n)%height/dz).and.k.ge.CEILING(bp(n)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n)%x*cos(phi)-bp(n)%y*sin(phi)
		yTSHD(1:4)=bp(n)%x*sin(phi)+bp(n)%y*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
		Ubound(i,j,k)=(bp(n)%u*cos_u(j)+bp(n)%v*sin_u(j))*rho(i,j,k)
		Vbound(i,j,k)=(-bp(n)%u*sin_v(j)+bp(n)%v*cos_v(j))*rho(i,j,k)
		Wbound(i,j,k)=bp(n)%w*rho(i,j,k)
		Wbound(i,j,MAX(k-1,0))=bp(n)%w*rho(i,j,k)
	   endif
	  enddo
	 enddo
	enddo
	ENDIF
	ENDIF
	ENDDO ! bedplume loop


	IF (LOA>0.) THEN ! ship:
	  do t=1,tmax_inWpuntTSHD
 	    i=i_inWpuntTSHD(t)
 	    j=j_inWpuntTSHD(t)		
 	    k=k_inWpuntTSHD(t)		
	    Wbound(i,j,k)=0.
	  enddo
	  do t=1,tmax_inUpuntTSHD
 	    i=i_inUpuntTSHD(t)
 	    j=j_inUpuntTSHD(t)		
 	    k=k_inUpuntTSHD(t)		
	    Ubound(i,j,k)=0.
	  enddo
	  do t=1,tmax_inVpuntTSHD
 	    i=i_inVpuntTSHD(t)
 	    j=j_inVpuntTSHD(t)		
 	    k=k_inVpuntTSHD(t)		
	    Vbound(i,j,k)=0.
	  enddo
	  do i=0,i1 ! prescribe UTSHD in obstacle:
	    do j=0,j1
	      do k=1,kbed(i,j)
	        Ubound(i,j,k)=Ubot_TSHD(j)*rho(i,j,k)
	        Vbound(i,j,k)=Vbot_TSHD(j)*rho(i,j,k)
	      enddo
	    enddo
	  enddo
	  !if (botstress.ge.1) then !if kn_TSHD>0 then tauTSHD independent of botstress
	  if (kn_TSHD.ge.0.0) then !apply tauTSHD independent of botstress
	   do t=1,tmax_inVpunt_tauTSHD
 	    i=i_inVpunt_tauTSHD(t)
 	    j=j_inVpunt_tauTSHD(t)		
 	    k=k_inVpunt_tauTSHD(t)
	    uu2=0.25*(Ubound2(i,j,k)+Ubound2(i,j+1,k)+Ubound2(i-1,j,k)+Ubound2(i-1,j+1,k))		
!	    call wall_fun_rho(Vbound(i,j,k),Ubound2(i,j,k),rho(i,j,k),dz,dt,kn_TSHD,kappa,(depth-bc_obst_h),U_b,nu_mol)
	    call wall_fun_rho(Vbound(i,j,k),uu2,rho(i,j,k),dz,dt,kn_TSHD,kappa,(depth-bc_obst_h),U_b,nu_mol)
	   enddo
	  endif
	  !if (botstress.ge.1) then !if kn_TSHD>0 then tauTSHD independent of botstress
          if (kn_TSHD.ge.0.0) then !apply tauTSHD independent of botstress
	   do t=1,tmax_inUpunt_tauTSHD
 	    i=i_inUpunt_tauTSHD(t)
 	    j=j_inUpunt_tauTSHD(t)		
 	    k=k_inUpunt_tauTSHD(t)		
	    vv2=0.25*(Vbound2(i,j,k)+Vbound2(i+1,j,k)+Vbound2(i,j-1,k)+Vbound2(i+1,j-1,k))
!	    call wall_fun_rho(Ubound(i,j,k),Vbound2(i,j,k),rho(i,j,k),dz,dt,kn_TSHD,kappa,(depth-bc_obst_h),U_b,nu_mol)
	    call wall_fun_rho(Ubound(i,j,k),vv2,rho(i,j,k),dz,dt,kn_TSHD,kappa,(depth-bc_obst_h),U_b,nu_mol)
	   enddo
	  endif
	  do t=1,tmax_inVpunt_rudder ! apply rudder 
 	    i=i_inVpunt_rudder(t)
 	    j=j_inVpunt_rudder(t)		
 	    k=k_inVpunt_rudder(t)
	    jp=MIN(j+1,j1) ! if j is j1, then jp=j1 --> uu2 is incorrect, but this is no problem as Ubound(i,j1,k) is not used	
	    uu1 = -sin_v(j)*Vbound(i,j,k)+(Ubound(i,j,k))*cos_v(j) ! u wanted along rudder (v wanted = zero)
	    uu2 = -sin_v(j)*Vbound(i,j,k)+(Ubound(i,jp,k))*cos_v(j) ! u wanted along rudder (v wanted = zero)
	    Vbound(i,j,k)= -sin_v(j)*0.5*(uu1+uu2)
	    Ubound(i,j,k)= cos_u(j)*uu1
	    Ubound(i,jp,k)= cos_u(jp)*uu2
	  enddo
	ELSE 
          !! Set boundary conditions plate:
	  if (kjet>0) then
	  Ubound(0:i1,0:j1,kmax-kjet+1:k1)=0. ! top layer (kjet) is NO slip boundary, let develop freely no slip is obtained through wall_fun
	  Vbound(0:i1,0:j1,kmax-kjet+1:k1)=0. ! top layer (kjet) is NO slip boundary, let develop freely no slip is obtained through wall_fun
	  Wbound(0:i1,0:j1,kmax-kjet:kmax)=0.
	  endif
	ENDIF

	! set boundary conditions vertical jet:
       	if (plumetseriesfile.eq.'') then
       		Wjet=W_j
       		f=Strouhal*ABS(W_j)/(radius_j*2.) !Strouholt number is 0.3
	else
       		Wjet=interpseries(plumetseries,plumeUseries,plumeseriesloc,tt)
       		f=Strouhal*ABS(W_j)/(radius_j*2.) !Strouholt number is 0.3
	endif
	jetcorr=pi/(2.*pi*(1/(1./W_j_powerlaw+1)-1./(1./W_j_powerlaw+2.))) !=1.22449 for W_j_powerlaw=7
       do t=1,tmax_inPpunt
 	i=i_inPpunt(t)
 	j=j_inPpunt(t)
 	fluc=0.
 	do n=1,azi_n
 	  fluc=fluc+sin(2.*pi*f*tt/n+azi_angle_p(i,j))
 	enddo
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  rr=1.-sqrt((xx**2+yy**2)/radius_j**2)
	  rr=MAX(rr,0.)
	if (outflow_overflow_down.eq.1) then
 		do k=kmax-kjet,kmax-kjet
 		  Wbound(i,j,k)=Wbound2(i,j,k) !Wjet*0.5*(rho(i,j,k)+rho(i,j,k+1))+Awjet/6.*Wjet*0.5*(rho(i,j,k)+rho(i,j,k+1))*fluc
 		enddo
 		do k=kmax-kjet+1,kmax
 		  Wbound(i,j,k)=(jetcorr*rr**(1./W_j_powerlaw)+fluc*Awjet/REAL(azi_n))*Wjet*0.5*(rho(i,j,k)+rho(i,j,k+1))
 		enddo
	else
 		do k=kmax-kjet,kmax
 		  Wbound(i,j,k)=Wbound2(i,j,k) !Wjet*0.5*(rho(i,j,k)+rho(i,j,k+1))+Awjet/6.*Wjet*0.5*(rho(i,j,k)+rho(i,j,k+1))*fluc
 		enddo
	 	Wbound(i,j,kmax)=(jetcorr*rr**(1./W_j_powerlaw)+fluc*Awjet/REAL(azi_n))*Wjet*0.5*(rho(i,j,kmax)+rho(i,j,k1))
	endif
       enddo

       do t=1,tmax_inUpunt
 	i=i_inUpunt(t)
 	j=j_inUpunt(t)
 	fluc=0.
 	do n=1,azi_n
 	  fluc=fluc+sin(2.*pi*f*tt/n+azi_angle_u(i,j))
 	enddo
	Ujetbc=Aujet/azi_n*Wjet*fluc
	Vjetbc=Avjet/azi_n*Wjet*fluc
	Ujetbc2=Ujetbc*cos(azi_angle_u(i,j))-Vjetbc*sin(azi_angle_u(i,j))
	Vjetbc2=Ujetbc*sin(azi_angle_u(i,j))+Vjetbc*cos(azi_angle_u(i,j))
	if (outflow_overflow_down.eq.1) then
		do k=kmax-kjet,kmax-kjet+1
		  Ubound(i,j,k)=Ubound2(i,j,k)
		enddo
		do k=kmax-kjet+2,k1
		  Ubound(i,j,k)=0.
		enddo	
	else
		do k=kmax-kjet,k1
		  Ubound(i,j,k)=Ubound2(i,j,k)
		enddo	
	  	Ubound(i,j,k1)=(rho(i,j,k1)+rho(i+1,j,k1))*(Ujetbc2*cos_u(j)+Vjetbc2*sin_u(j))-Ubound(i,j,kmax)
	endif
       enddo
       do t=1,tmax_inVpunt
 	i=i_inVpunt(t)
 	j=j_inVpunt(t)
 	fluc=0.
 	do n=1,azi_n
 	  fluc=fluc+sin(2.*pi*f*tt/n+azi_angle_v(i,j))
 	enddo
	Ujetbc=Aujet/azi_n*Wjet*fluc
	Vjetbc=Avjet/azi_n*Wjet*fluc
	Ujetbc2=Ujetbc*cos(azi_angle_v(i,j))-Vjetbc*sin(azi_angle_v(i,j))
	Vjetbc2=Ujetbc*sin(azi_angle_v(i,j))+Vjetbc*cos(azi_angle_v(i,j))
	jp=min(j+1,j1)	!this error is no issue as Vbound(j1) is not used
	if (outflow_overflow_down.eq.1) then
		do k=kmax-kjet,kmax-kjet+1
		  Vbound(i,j,k)=Vbound2(i,j,k)
		enddo
		do k=kmax-kjet+2,k1
		  Vbound(i,j,k)=0.
		enddo
	else
		do k=kmax-kjet,k1
		  Vbound(i,j,k)=Vbound2(i,j,k)
		enddo
	  	Vbound(i,j,k1)=(rho(i,j,k1)+rho(i,jp,k1))*(Ujetbc2*sin_v(j)+Vjetbc2*cos_v(j))-Vbound(i,j,kmax) !(Aujet/6.*Wjet*fluc)
	endif
       enddo

      !! Set boundary conditions horizontal jet2 in:
       	if (plumetseriesfile2.eq.'') then
       		Ujet2=U_j2
       		f=Strouhal*ABS(U_j2)/(radius_j2*2.) !Strouholt number is 0.3
	else
       		Ujet2=interpseries(plumetseries2,plumeUseries2,plumeseriesloc2,tt)
       		f=Strouhal2*ABS(U_j2)/(radius_j2*2.) !Strouholt number is 0.3
	endif
       do t=1,tmax_inPpunt2
 	  k=k_inPpunt2(t)
 	  j=j_inPpunt2(t)
	  zzz=k*dz-0.5*dz-zjet2
	  yy=Rp(0)*sin_u(j)
 	  fluc=0.
 	  do n=1,azi_n2
 	    fluc=fluc+sin(2.*pi*f*tt/n+atan2(zzz,yy))
 	  enddo
	    rr=1.-sqrt((zzz**2+yy**2)/radius_j2**2)
	    rr=MAX(rr,0.)
	    Ujetbc=jetcorr*Ujet2*rr**(1./W_j_powerlaw)+Aujet2/azi_n2*Ujet2*fluc 
	    Vjetbc=Avjet/azi_n*Ujet2*fluc
	    Ubound(0,j,k)=(Ujetbc*cos_u(j)+Vjetbc*sin_u(j))*0.5*(rho(0,j,k)+rho(1,j,k))
       enddo
       do t=1,tmax_inWpunt2
 	  k=k_inWpunt2(t)
 	  j=j_inWpunt2(t)
	  zzz=k*dz-zjet2
	  yy=Rp(0)*sin_u(j)
 	  fluc=0.
 	  do n=1,azi_n
 	    fluc=fluc+sin(2.*pi*f*tt/n+atan2(zzz,yy))
 	  enddo
	  Wjetbc=Awjet/azi_n*Ujet2*fluc
	  Vjetbc=Avjet/azi_n*Ujet2*fluc
	  Wbound(0,j,k)=(rho(0,j,k)+rho(0,j,k+1))*(cos(atan2(zzz,yy))*Wjetbc+sin(atan2(zzz,yy))*Vjetbc)-Wbound(1,j,k)
       enddo
       do t=1,tmax_inVpunt2
 	  k=k_inVpunt2(t)
 	  j=j_inVpunt2(t)
	  zzz=k*dz-0.5*dz-zjet2
	  yy=Rp(0)*sin_v(j)
	  rr=1.-sqrt((zzz**2+yy**2)/radius_j**2)
	  rr=MAX(rr,0.)
 	  fluc=0.
 	  do n=1,azi_n
 	    fluc=fluc+sin(2.*pi*f*tt/n+atan2(zzz,yy))
 	  enddo
	  Ujetbc=Ujet2*rr**(1./W_j_powerlaw)+Aujet2/azi_n2*Ujet2*fluc
	  Vjetbc=Avjet/azi_n*Ujet2*fluc
	  Wjetbc=Awjet/azi_n*Ujet2*fluc
	  jp=min(j+1,j1) !this error is no issue as Vbound(j1) is not used
	  Vbound(0,j,k)=(rho(0,j,k)+rho(0,jp,k))*
     &                  ((-Ujetbc*sin_v(j)+Vjetbc*cos_v(j))*cos(atan2(zzz,yy))-sin(atan2(zzz,yy))*Wjetbc)-Vbound(1,j,k)
       enddo

       do t=1,tmax_inWpunt_suction ! add suction in front of draghead
 	  k=k_inWpunt_suction(t)
 	  j=j_inWpunt_suction(t)
 	  i=i_inWpunt_suction(t)
	  Wjetbc=pi*radius_j**2*W_j*perc_dh_suction/(dr(i)*3.*Dsp) !(m3/s)/m2=m/s per suction cell (perc_dh_suction is corrected for number of dragheads)
	  Wbound(i,j,k)=Wjetbc*0.5*(rho(i,j,k)+rho(i,j,k+1))
       enddo


	! apply boundary conditions for Vbound again to fix small inconsistency in Vbound(:,j1,:) made in jet, jet2 and rudder
	call shiftf(Vbound,vbf) 
	call shiftb(Vbound,vbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then ! boundaries in j-direction

	  elseif (rank.eq.px-1) then
	
	  else
		do k=0,k1
		   do i=1,imax
		   Vbound(i,0,k) = Vbf(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   enddo
		enddo
	  endif
	else !periodic in y:
	  do k=0,k1
	   do i=1,imax
	     Vbound(i,0,k) = Vbf(i,k)
	     Vbound(i,j1,k) =Vbb(i,k)
	   enddo
	  enddo
	endif




      end

      subroutine make_UtransC_zeroIBM(Ubound,Vbound,Wbound)

      USE nlist

      implicit none

	integer t
	real Ubound(0:i1,0:j1,0:k1),Vbound(0:i1,0:j1,0:k1),Wbound(0:i1,0:j1,0:k1)
	real Ubound2(0:i1,0:j1,0:k1),Vbound2(0:i1,0:j1,0:k1),Wbound2(0:i1,0:j1,0:k1)


	Ubound2=Ubound
	Vbound2=Vbound
	Wbound2=Wbound
	  do t=1,tmax_inWpuntTSHD
 	    i=i_inWpuntTSHD(t)
 	    j=j_inWpuntTSHD(t)		
 	    k=k_inWpuntTSHD(t)		
	    Wbound(i,j,k)=0.
	  enddo
	  do t=1,tmax_inUpuntTSHD
 	    i=i_inUpuntTSHD(t)
 	    j=j_inUpuntTSHD(t)		
 	    k=k_inUpuntTSHD(t)		
	    Ubound(i,j,k)=0.
	  enddo
	  do t=1,tmax_inVpuntTSHD
 	    i=i_inVpuntTSHD(t)
 	    j=j_inVpuntTSHD(t)		
 	    k=k_inVpuntTSHD(t)		
	    Vbound(i,j,k)=0.
	  enddo

        do t=1,tmax_inPpunt ! make vel jet not zero
 	 i=i_inPpunt(t)
 	 j=j_inPpunt(t)
 	 do k=kmax-kjet,kmax
 	   Wbound(i,j,k)=Wbound2(i,j,k) 
 	 enddo
	enddo
        do t=1,tmax_inVpunt ! make vel jet not zero
 	 i=i_inVpunt(t)
 	 j=j_inVpunt(t)
 	 do k=kmax-kjet,k1
 	   Vbound(i,j,k)=Vbound2(i,j,k) 
 	 enddo
	enddo
        do t=1,tmax_inUpunt ! make vel jet not zero
 	 i=i_inUpunt(t)
 	 j=j_inUpunt(t)
 	 do k=kmax-kjet,k1
 	   Ubound(i,j,k)=Ubound2(i,j,k) 
 	 enddo
	enddo


	end


      subroutine bound_c(Cbound,cjet,n)
      
      USE nlist

      implicit none
c
!       include 'param.txt'
!       include 'common.txt'

c
      integer jtmp,t,kcheck,n,n2,inout
c
      real Cbound(0:i1,0:j1,0:k1),cjet,Cbound2(0:i1,0:j1,0:k1)
      real ubb(0:i1,0:k1),val,theta,rbc,xx,yy,r_orifice2,rjet,theta_U,theta_V
      real cbf(0:i1,0:k1)
      real cbb(0:i1,0:k1)
	real xTSHD(4),yTSHD(4),phi
c
c
c*************************************************************
c
c     Subroutine bound sets the boundary conditions for Cbound
c     except for the diffusion coefficients. These are set in submod.
c     The common boundary conditions for the pressure are set in mkgrid.
c
c     Set boundary conditions for j=0 and j=j1. Because of the
c     fact that the index j denotes the tangential direction,
c     we have to set the values at j=0 equal to the values at
c     j=jmax and the values at j=j1 equal to the values at j=1.
c
c*************************************************************
c
c*************************************************************
c get stuff from other CPU's

	
	call shiftf(Cbound,cbf) 
	call shiftb(Cbound,cbb) 

!	write(*,*),'rank,px,INT(px/2)',rank,px,INT(px/2)

	! Cbcoarse1 and Cbcoarse2 are zero when no bcfile is used, or when nfrac is zero in bcfile;
	! in that case simply zero lateral and inflow bc are applied to C

	
	if (periodicy.eq.0) then
		if (rank.eq.0) then ! boundaries in j-direction
			do k=1,kmax
			   do i=1,imax
			   Cbound(i,0,k) = Cbcoarse1(n,i,k) !Cbound(i,1,k) !cbf(i,k)
			   Cbound(i,j1,k) =cbb(i,k) 
			   enddo
			enddo
		elseif (rank.eq.px-1) then
			do k=1,kmax
			   do i=1,imax
			   Cbound(i,0,k) = cbf(i,k)
			   Cbound(i,j1,k) = Cbcoarse1(n,i,k) !Cbound(i,jmax,k) !cbb(i,k) 
			   enddo
			enddo	
		else
			do k=1,kmax
			   do i=1,imax
			   Cbound(i,0,k) = cbf(i,k)
			   Cbound(i,j1,k) =cbb(i,k) 
			   enddo
			enddo
		endif
	elseif (periodicy.eq.1) then ! periodic in y:
		do k=1,kmax
		   do i=1,imax
			   Cbound(i,0,k) = cbf(i,k)
			   Cbound(i,j1,k) =cbb(i,k) 
		   enddo
		enddo
	elseif (periodicy.eq.2) then ! free slip in y:
		if (rank.eq.0) then ! boundaries in j-direction
			do k=1,kmax
			   do i=1,imax
			   Cbound(i,0,k) = Cbound(i,1,k) !cbf(i,k)
			   Cbound(i,j1,k) =cbb(i,k) 
			   enddo
			enddo
		elseif (rank.eq.px-1) then
			do k=1,kmax
			   do i=1,imax
			   Cbound(i,0,k) = cbf(i,k)
			   Cbound(i,j1,k) =Cbound(i,jmax,k) !cbb(i,k) 
			   enddo
			enddo	
		else
			do k=1,kmax
			   do i=1,imax
			   Cbound(i,0,k) = cbf(i,k)
			   Cbound(i,j1,k) =cbb(i,k) 
			   enddo
			enddo
		endif
	endif

	if (periodicx.eq.0) then
	      do k=1,kmax ! boundaries in i-direction
		 do j=0,j1
			   Cbound(0,j,k)    =    Cbcoarse2(n,j,k) !Cbound(1,j,k)
			   Cbound(i1,j,k)   =    Cbound(imax,j,k)
		 enddo   
	      enddo
       elseif (periodicx.eq.2) then ! no outflow in x direction:
	      do k=1,kmax ! boundaries in i-direction
		 do j=0,j1
			   Cbound(0,j,k)    =    0. 
			   Cbound(i1,j,k)   =    0. 
		 enddo   
	      enddo
	else
	      do k=1,kmax ! boundaries in i-direction
		 do j=0,j1
			   Cbound(0,j,k)    =    Cbound(imax,j,k)
			   Cbound(i1,j,k)   =    Cbound(1,j,k)
		 enddo   
	      enddo
	endif
	
      do j=0,j1 ! boundaries in k-direction
         do i=0,i1
         Cbound(i,j,k1)   = Cbound(i,j,kmax)
         Cbound(i,j,0)    = Cbound(i,j,1)
         enddo
       enddo
!	endif

       !! Set boundary conditions jet in:
	IF (outflow_overflow_down.eq.1) THEN  
      do t=1,tmax_inPpunt
		i=i_inPpunt(t)
		j=j_inPpunt(t)
		do k=k1-kjet+1,k1
	      Cbound(i,j,k)=cjet !0.
		enddo
	  enddo
	ELSE
      do t=1,tmax_inPpunt
		i=i_inPpunt(t)
		j=j_inPpunt(t)	
	    Cbound(i,j,k1)=cjet !Cbound(i,j,k1-kjet)=cjet	
	  enddo
	ENDIF
      
       !! Set boundary conditions jet2 in:
      do t=1,tmax_inPpunt2
	k=k_inPpunt2(t)
	j=j_inPpunt2(t)
!	do k=1,kjet !k1-kjet+1,k1
!	  Cbound(i,j,k)=cjet !0.
!	enddo
	Cbound(0,j,k)=cjet !Cbound(i,j,k1-kjet)=cjet	
      enddo

	Cbound2=Cbound

	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0)
     &     .or.(bp(n2)%forever.eq.0.and.time_n.lt.bp(n2)%t0.and.time_np.gt.bp(n2)%t0)) THEN
	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
      !! Search for P,V:
      do k=0,k1
       do i=0,i1  
         do j=j1,0,-1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
	  xx=Rp(i)*cos_u(j)-schuif_x
	  yy=Rp(i)*sin_u(j)
	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
	  ELSE 
	 	inout=0
	  ENDIF
	  if (inout.eq.1) then
		Cbound(i,j,k)=bp(n2)%c(n)
		! rho is calculated in state called after fkdat
	   endif
	  enddo
	 enddo
	enddo
	  ! remove sediment from obstacles/TSHD after placement of bedplume:
	  do t=1,tmax_inPpuntTSHD ! when no TSHD then this loop is skipped
 	    k=k_inPpuntTSHD(t)		
	    i=i_inPpuntTSHD(t)
            j=j_inPpuntTSHD(t)
            Cbound(i,j,k)=0.  ! remove sediment from hull
	  enddo
	ENDIF
	ENDDO ! bedplume loop

	kcheck=kmax-(kjet-1) 
	  do t=1,tmax_inPpuntTSHD ! when no TSHD then this loop is skipped
 	    k=k_inPpuntTSHD(t)		
	    i=i_inPpuntTSHD(t)
            j=j_inPpuntTSHD(t)
	      if (k.eq.kcheck) then 
	        Cbound(i,j,k-1)=Cbound(i,j,k-1)+Cbound(i,j,k) ! move sediment from first (lowest) cell in TSHD-hull to first cell in fluid
		!! Idea is to undo the effect of diffusion and air bubble rising up into hull
	      endif
	      if (k.gt.0.and.k.eq.kbed(i,j)) then
	        Cbound(i,j,kbed(i,j)+1)=Cbound(i,j,kbed(i,j)+1)+Cbound(i,j,k) ! move sediment from heightest cell in obstacle to first cell in fluid
		!! Idea is to undo the effect of diffusion into the obstacle
	      endif
            Cbound(i,j,k)=0.  ! remove sediment from hull, not only in lowest line of cells in hull
	  enddo
      do t=1,tmax_inPpunt
	i=i_inPpunt(t)
	j=j_inPpunt(t)
	do k=k1-kjet-2,k1
	  Cbound(i,j,k)=Cbound2(i,j,k)
	enddo
      enddo

      end

	
	subroutine bound_cbot(Cbound)
      
      USE nlist

      implicit none
c
!       include 'param.txt'
!       include 'common.txt'

c
      integer jtmp,t,kcheck
c
      real Cbound(0:i1,0:j1)
      real cbf(0:i1)
      real cbb(0:i1)

c
c
c*************************************************************
c
c     Subroutine bound sets the boundary conditions for Cbound
c     except for the diffusion coefficients. These are set in submod.
c     The common boundary conditions for the pressure are set in mkgrid.
c
c     Set boundary conditions for j=0 and j=j1. Because of the
c     fact that the index j denotes the tangential direction,
c     we have to set the values at j=0 equal to the values at
c     j=jmax and the values at j=j1 equal to the values at j=1.
c
c*************************************************************
c
c*************************************************************
c get stuff from other CPU's
	  
	call shiftf_l(Cbound,cbf) 
	call shiftb_l(Cbound,cbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
		if (rank.eq.0) then ! boundaries in j-direction
		   do i=1,imax
		   Cbound(i,0) = Cbound(i,1) !cbf(i,k)
		   Cbound(i,j1) =cbb(i) 
		   enddo
		elseif (rank.eq.px-1) then
		   do i=1,imax
		   Cbound(i,0) = cbf(i)
		   Cbound(i,j1) =Cbound(i,jmax) !cbb(i,k) 
		   enddo
		else
		   do i=1,imax
		   Cbound(i,0) = cbf(i)
		   Cbound(i,j1) =cbb(i) 
		   enddo
		endif
	else
	   do i=1,imax
		   Cbound(i,0) = cbf(i)
		   Cbound(i,j1) =cbb(i) 
	   enddo
	endif
	
	 ! boundaries in i-direction
	if (periodicx.eq.0) then
         do j=0,j1
		   Cbound(0,j)    =    0. !Cbound(1,j)
		   Cbound(i1,j)   =    Cbound(imax,j)
         enddo   
	else 
         do j=0,j1
		   Cbound(0,j)    =    Cbound(imax,j)
		   Cbound(i1,j)   =    Cbound(1,j)
         enddo   
	endif

		end

	subroutine wall_fun(uu,vv,rr,dz,dt,kn,kappa,depth,U_b,nu_mol)
		
	implicit none

c*************************************************************
c	Wall function
c	Determine tau-wall with sqrt(uu**2,vv**2) and adapt uu	
c	Subtract tau_wall/rho*dt*dx*dy/(dx*dy*dz) from uu (m/s)
c*************************************************************
	real uu,vv,rr,dz,dt,absU,ust,z0,kn,kappa,yplus,tau,depth,U_b,nu_mol
	integer tel

	absU=sqrt((uu)**2+(vv)**2)
	ust=0.1*absU
	if (kn>0.) then !walls with rougness (log-law used which can be hydr smooth or hydr rough):
		do tel=1,10 ! 10 iter is more than enough
			z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
			ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
		enddo
!		yplus=0.5*dz*ust/nu_mol
!	   	if (yplus<30.) then
!		  do tel=1,10 ! 10 iter is more than enough
!			yplus=0.5*dz*ust/nu_mol
!			ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
!		  enddo	
!		endif
!		if (yplus<5.) then !viscous sublayer uplus=yplus
!			ust=sqrt(absU*nu_mol/(0.5*dz))
!		endif
!		if (yplus<11.225) then	!viscous sublayer uplus=yplus
!			ust=sqrt(absU*nu_mol/(0.5*dz))
!		endif
	else
		do tel=1,10 ! 10 iter is more than enough
			yplus=0.5*dz*ust/nu_mol
			ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
		enddo
	   	if (yplus<30.) then
		  do tel=1,10 ! 10 iter is more than enough
			yplus=0.5*dz*ust/nu_mol
			ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
		  enddo	
		endif
		if (yplus<5.) then !viscous sublayer uplus=yplus
			ust=sqrt(absU*nu_mol/(0.5*dz))
		endif

	endif
! 	write(*,*) 'ust wall function',ust

	tau=ust*ust  ! omit rho otherwise first *rho then /rho to get dimensions right
	uu = uu - tau*dt/dz*uu/MAX(absU,1.e-6) 	!! only uu is adjusted
	end

	subroutine wall_fun_rho(uu,vv,rr,dz,dt,kn,kappa,depth,U_b,nu_mol)
		
	implicit none

c*************************************************************
c	Wall function
c	Determine tau-wall with sqrt(uu**2,vv**2) and adapt uu	
c	Subtract tau_wall*dt*dx*dy/(dx*dy*dz) from uu (kg/m3*m/s)
c*************************************************************
	real uu,vv,rr,dz,dt,absU,ust,z0,kn,kappa,yplus,tau,depth,U_b,nu_mol
	integer tel

	absU=sqrt((uu/rr)**2+(vv/rr)**2)
	ust=0.1*absU
	if (kn>0.) then !walls with rougness (log-law used which can be hydr smooth or hydr rough):
		do tel=1,10 ! 10 iter is more than enough
			z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
			ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
		enddo
!		yplus=0.5*dz*ust/nu_mol
!	   	if (yplus<30.) then
!		  do tel=1,10 ! 10 iter is more than enough
!			yplus=0.5*dz*ust/nu_mol
!			ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
!		  enddo	
!		endif
!		if (yplus<5.) then !viscous sublayer uplus=yplus
!			ust=sqrt(absU*nu_mol/(0.5*dz))
!		endif
!		if (yplus<11.225) then	!viscous sublayer uplus=yplus
!			ust=sqrt(absU*nu_mol/(0.5*dz))
!		endif
	else
		do tel=1,10 ! 10 iter is more than enough
			yplus=0.5*dz*ust/nu_mol
			ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
		enddo
	   	if (yplus<30.) then
		  do tel=1,10 ! 10 iter is more than enough
			yplus=0.5*dz*ust/nu_mol
			ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
		  enddo	
		endif
		if (yplus<5.) then !viscous sublayer uplus=yplus
			ust=sqrt(absU*nu_mol/(0.5*dz))
		endif
	endif

	tau=rr*ust*ust
	uu = uu - tau*dt/dz*uu/rr/MAX(absU,1.e-6) 	!! only uu is adjusted
	end


	REAL function interpseries(tseries,series,sloc,tt)

	REAL tseries(1:10000),series(1:10000),tt
	INTEGER sloc

	IF (tt>tseries(sloc+1)) THEN
		sloc=sloc+1
	ENDIF

!	write(*,*)'sloc,tseries(sloc),tseries(sloc+1)',sloc,tseries(sloc),tseries(sloc+1)
!	write(*,*)'series(sloc),series(sloc+1),interpval',series(sloc),series(sloc+1),
!     & (tt-tseries(sloc))/(tseries(sloc+1)-tseries(sloc))*(series(sloc+1)-series(sloc))+series(sloc)


	interpseries=(tt-tseries(sloc))/(tseries(sloc+1)-tseries(sloc))*(series(sloc+1)-series(sloc))+series(sloc)	



	end function interpseries	
	  
	subroutine shiftb(UT,UP)

      USE nlist


      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf,ierr
!       parameter (ileng= (k1+1)*(i1+1) )
      include 'mpif.h'
      integer itag,status(MPI_STATUS_SIZE),l
      real*8 ut(0:i1,0:j1,0:k1)
      real*8 up(0:i1,0:k1),UTMP(0:I1,0:K1)
	!real up(0:i1,0:k1),UTMP(0:I1,0:K1)
      do i=0,i1
	 do k=0,k1
	  utmp(i,k) =UT(i,1,k)
          enddo
      enddo
      itag = 10
      ileng = (k1+1)*(i1+1)
      rankf=rank+1
      rankb=rank-1
	if (periodicy.eq.0.or.periodicy.eq.2) then
           if(rank.eq.px-1)rankf=MPI_PROC_NULL
           if(rank.eq.   0)rankb=MPI_PROC_NULL
	else 
      	   if(rank.eq.px-1)rankf=0 ! MPI_PROC_NULL
	   if(rank.eq.   0)rankb=px-1 !MPI_PROC_NULL
	endif
      call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankb,itag,
     $                  up ,ileng,MPI_REAL8,rankf,itag, MPI_COMM_WORLD,status,ierr)
!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankb,itag,
!     $                  up ,ileng,MPI_REAL,rankf,itag, MPI_COMM_WORLD,status,ierr)
c      if (rank.eq.   0) then
c         call MPI_SEND(utmp  ,ileng,MPI_REAL8,px-1,itag,MPI_COMM_WORLD,ierr)
c        endif
c       if (rank.eq.px-1) then
c           call MPI_RECV(up,ileng,MPI_REAL8,0 ,itag,MPI_COMM_WORLD,status,ierr)
c          endif

      end

	subroutine shiftf(UT,UP)
      USE nlist

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf
!       parameter (ileng= (k1+1)*(i1+1))
      include 'mpif.h'
      real*8 UT(0:i1,0:j1,0:k1),UP(0:i1,0:k1),UTMP(0:i1,0:k1)
	!real UT(0:i1,0:j1,0:k1),UP(0:i1,0:k1),UTMP(0:i1,0:k1)
      integer  itag,status(MPI_STATUS_SIZE),l,ierr
      itag = 11
      ileng = (k1+1)*(i1+1)
      do k=0,k1
        do i=0,i1
	  UTMP(i,k) =UT(i,jmax,k)
          enddo
      enddo
      rankf=rank+1
      rankb=rank-1
	if (periodicy.eq.0.or.periodicy.eq.2) then
           if(rank.eq.px-1)rankf=MPI_PROC_NULL
           if(rank.eq.   0)rankb=MPI_PROC_NULL
	else 
      	   if(rank.eq.px-1)rankf=0 ! MPI_PROC_NULL
	   if(rank.eq.   0)rankb=px-1 !MPI_PROC_NULL
	endif

      call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankf,itag,
     $                  up   ,ileng,MPI_REAL8,rankb,itag, MPI_COMM_WORLD,status,ierr)

!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankf,itag,
!     $                  up   ,ileng,MPI_REAL,rankb,itag, MPI_COMM_WORLD,status,ierr)

c      if (rank.eq.px-1) then
c       call MPI_SEND(UTMP,ileng,MPI_REAL8,0,itag,MPI_COMM_WORLD,ierr)
c       else
c       if (rank.eq.0) then
c       call MPI_RECV(UP ,ileng,MPI_REAL8,px-1,itag,MPI_COMM_WORLD,status,ierr)
c       endif
c      endif
      end

	subroutine shiftb_l(UT1,UP1)

      USE nlist


      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf,ierr
!       parameter (ileng= (k1+1)*(i1+1) )
      include 'mpif.h'
      integer itag1,status(MPI_STATUS_SIZE),l
      real*8 ut1(0:i1,0:j1)
      real*8 up1(0:i1),UTMP1(0:I1)
	!real up(0:i1,0:k1),UTMP(0:I1,0:K1)
      do i=0,i1
	  utmp1(i) =UT1(i,1)
      enddo
      itag1 = 12
      ileng = i1+1
      rankf=rank+1
      rankb=rank-1
	if (periodicy.eq.0.or.periodicy.eq.2) then
           if(rank.eq.px-1)rankf=MPI_PROC_NULL
           if(rank.eq.   0)rankb=MPI_PROC_NULL
	else 
      	   if(rank.eq.px-1)rankf=0 ! MPI_PROC_NULL
	   if(rank.eq.   0)rankb=px-1 !MPI_PROC_NULL
	endif

      call mpi_sendrecv(utmp1 ,ileng,MPI_REAL8,rankb,itag1,
     $                  up1 ,ileng,MPI_REAL8,rankf,itag1, MPI_COMM_WORLD,status,ierr)

!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankb,itag,
!     $                  up ,ileng,MPI_REAL,rankf,itag, MPI_COMM_WORLD,status,ierr)
c      if (rank.eq.   0) then
c         call MPI_SEND(utmp  ,ileng,MPI_REAL8,px-1,itag,MPI_COMM_WORLD,ierr)
c        endif
c       if (rank.eq.px-1) then
c           call MPI_RECV(up,ileng,MPI_REAL8,0 ,itag,MPI_COMM_WORLD,status,ierr)
c          endif

      end

	subroutine shiftf_l(UT,UP)
      USE nlist

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf
!       parameter (ileng= (k1+1)*(i1+1))
      include 'mpif.h'
      real*8 UT(0:i1,0:j1),UP(0:i1),UTMP(0:i1)
	!real UT(0:i1,0:j1,0:k1),UP(0:i1,0:k1),UTMP(0:i1,0:k1)
      integer  itag2,status(MPI_STATUS_SIZE),l,ierr
      itag2 = 13
      ileng = i1+1
        do i=0,i1
	  UTMP(i) =UT(i,jmax)
          enddo
      rankf=rank+1
      rankb=rank-1
	if (periodicy.eq.0.or.periodicy.eq.2) then
           if(rank.eq.px-1)rankf=MPI_PROC_NULL
           if(rank.eq.   0)rankb=MPI_PROC_NULL
	else 
      	   if(rank.eq.px-1)rankf=0 ! MPI_PROC_NULL
	   if(rank.eq.   0)rankb=px-1 !MPI_PROC_NULL
	endif

      call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankf,itag2,
     $                  up   ,ileng,MPI_REAL8,rankb,itag2, MPI_COMM_WORLD,status,ierr)
!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankf,itag,
!     $                  up   ,ileng,MPI_REAL,rankb,itag, MPI_COMM_WORLD,status,ierr)

c      if (rank.eq.px-1) then
c       call MPI_SEND(UTMP,ileng,MPI_REAL8,0,itag,MPI_COMM_WORLD,ierr)
c       else
c       if (rank.eq.0) then
c       call MPI_RECV(UP ,ileng,MPI_REAL8,px-1,itag,MPI_COMM_WORLD,status,ierr)
c       endif
c      endif


      end

	subroutine shiftb_l2(UT1,UP1)

      USE nlist


      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf,ierr
!       parameter (ileng= (k1+1)*(i1+1) )
      include 'mpif.h'
      integer itag1,status(MPI_STATUS_SIZE),l
      real*8 ut1(0:i1,0:j1)
      real*8 up1(0:i1),UTMP1(0:I1)
	!real up(0:i1,0:k1),UTMP(0:I1,0:K1)
      do i=0,i1
	  utmp1(i) =UT1(i,2)
      enddo
      itag1 = 12
      ileng = i1+1
      rankf=rank+1
      rankb=rank-1
	if (periodicy.eq.0.or.periodicy.eq.2) then
           if(rank.eq.px-1)rankf=MPI_PROC_NULL
           if(rank.eq.   0)rankb=MPI_PROC_NULL
	else 
      	   if(rank.eq.px-1)rankf=0 ! MPI_PROC_NULL
	   if(rank.eq.   0)rankb=px-1 !MPI_PROC_NULL
	endif

      call mpi_sendrecv(utmp1 ,ileng,MPI_REAL8,rankb,itag1,
     $                  up1 ,ileng,MPI_REAL8,rankf,itag1, MPI_COMM_WORLD,status,ierr)

!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankb,itag,
!     $                  up ,ileng,MPI_REAL,rankf,itag, MPI_COMM_WORLD,status,ierr)
c      if (rank.eq.   0) then
c         call MPI_SEND(utmp  ,ileng,MPI_REAL8,px-1,itag,MPI_COMM_WORLD,ierr)
c        endif
c       if (rank.eq.px-1) then
c           call MPI_RECV(up,ileng,MPI_REAL8,0 ,itag,MPI_COMM_WORLD,status,ierr)
c          endif

      end

	subroutine shiftf_l2(UT,UP)
      USE nlist

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf
!       parameter (ileng= (k1+1)*(i1+1))
      include 'mpif.h'
      real*8 UT(0:i1,0:j1),UP(0:i1),UTMP(0:i1)
	!real UT(0:i1,0:j1,0:k1),UP(0:i1,0:k1),UTMP(0:i1,0:k1)
      integer  itag2,status(MPI_STATUS_SIZE),l,ierr
      itag2 = 13
      ileng = i1+1
        do i=0,i1
	  UTMP(i) =UT(i,jmax-1)
          enddo
      rankf=rank+1
      rankb=rank-1
	if (periodicy.eq.0.or.periodicy.eq.2) then
           if(rank.eq.px-1)rankf=MPI_PROC_NULL
           if(rank.eq.   0)rankb=MPI_PROC_NULL
	else 
      	   if(rank.eq.px-1)rankf=0 ! MPI_PROC_NULL
	   if(rank.eq.   0)rankb=px-1 !MPI_PROC_NULL
	endif

      call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankf,itag2,
     $                  up   ,ileng,MPI_REAL8,rankb,itag2, MPI_COMM_WORLD,status,ierr)
!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankf,itag,
!     $                  up   ,ileng,MPI_REAL,rankb,itag, MPI_COMM_WORLD,status,ierr)

c      if (rank.eq.px-1) then
c       call MPI_SEND(UTMP,ileng,MPI_REAL8,0,itag,MPI_COMM_WORLD,ierr)
c       else
c       if (rank.eq.0) then
c       call MPI_RECV(UP ,ileng,MPI_REAL8,px-1,itag,MPI_COMM_WORLD,status,ierr)
c       endif
c      endif
	end 

