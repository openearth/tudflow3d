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


      subroutine bound(Ubound,Vbound,Wbound,rho,botstress,tt,Ub1in,Vb1in,Wb1in,Ub2in,Vb2in,Wb2in,Ub3in,Vb3in,Wb3in)
     
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
	integer botstress,n,jbeg,jend
      real Ub1(0:i1,0:k1),Vb1(0:i1,0:k1),Wb1(0:i1,0:k1),Ub2(0:j1,0:k1),Vb2(0:j1+1,0:k1),Wb2(0:j1,0:k1)
      real,intent(in) :: Ub1in(0:i1,0:k1),Vb1in(0:i1,0:k1),Wb1in(0:i1,0:k1),Ub2in(0:j1,0:k1),Vb2in(0:j1+1,0:k1),Wb2in(0:j1,0:k1)
	  real,intent(in) :: Ub3in(0:i1,0:j1),Vb3in(0:i1,0:j1),Wb3in(0:i1,0:j1)
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
	! 10-3-2020 noted that perhaps exact location wave inflow bc is not consistent with staggered layout near boundary, have to research later [Ub1in,Vb1in,Wb1in are defined at Uloc; Ub2in,Vb2in,Wb2in are defined at Vloc; but Ub1,Vb1,Wb1,Ub2,Vb2,Wb2 defined at u,v,w loc respectively and they are added without correcting
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
	    IF (periodicx.eq.3.or.periodicx.eq.4) THEN
			i=i1
		ELSE 
			i=0
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
	


!c get stuff from other CPU's
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
		   Vbound(i,j1,k) = Vbound(i,jmax,k)
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
		   Vbound(i,j1,k) = 0.
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

	if (periodicx.eq.0.and.monopile.eq.-1) then
		if (Uoutflow.eq.2) then ! convective outflow boundary
          do k=1,kmax ! boundaries in i-direction
           do j=0,j1
		    Wbc2=0.
		    if (bcfile.ne.'') then
			Ubc2(j,k)=Ubcoarse2(j,k)
			Vbc2(j,k)=Vbcoarse2(j,k)
			Wbc2=Wbcoarse2(j,k)
		    endif
		    Ubc=Ub2(j,k)+Ubc2(j,k)
		    Vbc=0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k)
		    Wbc=0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2	
		   Ubound(0,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
		   Vbound(0,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
		   Wbound(0,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
		   Ubc = MAX(Ubc2(j,k)*cos_u(j)+Vbc2(j,k)*sin_u(j),0.)  !Unormal outflow boundary at U-loc; should be >0 otherwise no outflow and divide by zero 
!		   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k)-dr(imax)/(Ubc*dt)*(Unew(imax-1,j,k)-Uold(imax-1,j,k)),0.)
!		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
!		   Vbound(i1,j,k)   =    Vbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Vnew(imax,j,k)-Vold(imax,j,k))		   
!		   Wbound(i1,j,k)   =    Wbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Wnew(imax,j,k)-Wold(imax,j,k))
		   
		   ! solving Ui+1^n+1-Ui+1^n)/dt+Unormal(Ui+1^n+1-Ui^n+1)/dx=0. !implicit treatment of dUdn term
		   Ubound(imax,j,k) = (Unew(imax,j,k)+Ubc*dt/dr(imax)*Ubound(imax-1,j,k))/(1.+Ubc*dt/dr(imax)) 
		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
		   Vbound(i1,j,k) = (Vnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Vbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
		   Wbound(i1,j,k) = (Wnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Wbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 

           enddo   
          enddo		
		else !Neumann outflow boundary 
          do k=1,kmax ! boundaries in i-direction
           do j=0,j1
		    Wbc2=0.
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
		endif 
       elseif (periodicx.eq.0.and.monopile>0) then ! closed boundary at i=0 (flow past circular cylinder); mixed inflow/outflow bc at imax; monopile=1=partial slip with tau-wall monopile2=free slip
		if (monopile.eq.4) then !no slip monopile wall 
         do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    0.
		   Vbound(0,j,k)    =    -Vbound(1,j,k)
		   Wbound(0,j,k)    =    -Wbound(1,j,k)
           enddo   
         enddo	
		 else 
		  do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    0.
		   Vbound(0,j,k)    =    Vbound(1,j,k)
		   Wbound(0,j,k)    =    Wbound(1,j,k)
           enddo   
          enddo	
		 endif 
		 if (rank.le.CEILING(REAL(px)/4.)-1) then ! first 1/4 of domain is outflow bc 
			jend = CEILING(REAL(jmax*px)*0.25) !end of first 1/4 global j-index
			jend = MIN(jend-rank*jmax,j1) !end of first 1/4 local j-index 
			if (Uoutflow.eq.2) then ! convective outflow boundary
			 do j=0,jend
			   do k=1,kmax 
				   Ubc = MAX(Ubc2(j,k)*cos_u(j)+Vbc2(j,k)*sin_u(j),0.)  !Unormal outflow boundary at U-loc; should be >0 otherwise no outflow and divide by zero 
				    !if (k.eq.1) write(*,*),'rank,j,Ubc',rank,j,Ubc
!				   Wbound(i1,j,k)   =    Wbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Wnew(imax,j,k)-Wold(imax,j,k))
!				   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k)-dr(imax)/(Ubc*dt)*(Unew(imax-1,j,k)-Uold(imax-1,j,k)),0.)
!				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
!				   Vbound(i1,j,k)   =    Vbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Vnew(imax,j,k)-Vold(imax,j,k))		   
				   
				   ! solving Ui+1^n+1-Ui+1^n)/dt+Unormal(Ui+1^n+1-Ui^n+1)/dx=0. !implicit treatment of dUdn term
				   Ubound(imax,j,k) = (Unew(imax,j,k)+Ubc*dt/dr(imax)*Ubound(imax-1,j,k))/(1.+Ubc*dt/dr(imax)) 
				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
				   Vbound(i1,j,k) = (Vnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Vbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
				   Wbound(i1,j,k) = (Wnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Wbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 				   
			   enddo   
			 enddo	
			else !Neumann outflow 
			 do j=0,jend
			   do k=1,kmax 
			     Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
			     Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
			     Vbound(i1,j,k)    =    Vbound(imax,j,k)
			     Wbound(i1,j,k)    =    Wbound(imax,j,k)
			   enddo   
			 enddo			
			endif
			 do j=jend+1,j1  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
					Wbc2=0.
					if (bcfile.ne.'') then
						Ubc2(j,k)=Ubcoarse2(j,k) 
						Vbc2(j,k)=Vbcoarse2(j,k)
						Wbc2=Wbcoarse2(j,k)
					endif
					Ubc=Ub2(j,k)+Ubc2(j,k)
					Vbc=0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k)
					Wbc=0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2	
				    Ubound(imax,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Ubound(i1,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Vbound(i1,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
				    Wbound(i1,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
			   enddo   
			 enddo
		 elseif (rank.ge.px-CEILING(REAL(px)/4.)-1) then ! last 1/4 of domain is outflow bc 
			jbeg = jmax*px - CEILING(REAL(jmax*px)*0.25) ! begin of last 1/4 global j-index
			jbeg = MAX(jbeg-rank*jmax,0) !begin of last 1/4 local j-index
			if (Uoutflow.eq.2) then ! convective outflow boundary
			 do j=jbeg,j1 ! outflow 
			   do k=1,kmax 
				   Ubc = MAX(Ubc2(j,k)*cos_u(j)+Vbc2(j,k)*sin_u(j),0.)  !Unormal outflow boundary at U-loc; should be >0 otherwise no outflow and divide by zero 
				   				   ! if (k.eq.1) write(*,*),'rank,j,Ubc',rank,j,Ubc
!				   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k)-dr(imax)/(Ubc*dt)*(Unew(imax-1,j,k)-Uold(imax-1,j,k)),0.)
!				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
!				   Vbound(i1,j,k)   =    Vbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Vnew(imax,j,k)-Vold(imax,j,k))		   
!				   Wbound(i1,j,k)   =    Wbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Wnew(imax,j,k)-Wold(imax,j,k))
				   
				   ! solving Ui+1^n+1-Ui+1^n)/dt+Unormal(Ui+1^n+1-Ui^n+1)/dx=0. !implicit treatment of dUdn term
				   Ubound(imax,j,k) = (Unew(imax,j,k)+Ubc*dt/dr(imax)*Ubound(imax-1,j,k))/(1.+Ubc*dt/dr(imax)) 
				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
				   Vbound(i1,j,k) = (Vnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Vbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
				   Wbound(i1,j,k) = (Wnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Wbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 				   
			   enddo   
			 enddo	
			else !Neumann outflow condition 
			 do j=jbeg,j1 ! outflow 
			   do k=1,kmax 
			     Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
			     Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
			     Vbound(i1,j,k)    =    Vbound(imax,j,k)
			     Wbound(i1,j,k)    =    Wbound(imax,j,k)
			   enddo   
			 enddo				
			endif
			 do j=0,MIN(j1,jbeg-1)  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
					Wbc2=0.
					if (bcfile.ne.'') then
						Ubc2(j,k)=Ubcoarse2(j,k) 
						Vbc2(j,k)=Vbcoarse2(j,k)
						Wbc2=Wbcoarse2(j,k)
					endif
					Ubc=Ub2(j,k)+Ubc2(j,k)
					Vbc=0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k)
					Wbc=0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2	
				    Ubound(imax,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Ubound(i1,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Vbound(i1,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
				    Wbound(i1,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
			   enddo   
			 enddo			 
		 else !middle 2/4 is inflow bc 
			 do j=0,j1  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
					Wbc2=0.
					if (bcfile.ne.'') then
						Ubc2(j,k)=Ubcoarse2(j,k) 
						Vbc2(j,k)=Vbcoarse2(j,k)
						Wbc2=Wbcoarse2(j,k)
					endif
					Ubc=Ub2(j,k)+Ubc2(j,k)
					Vbc=0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k)
					Wbc=0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2	
				    Ubound(imax,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Ubound(i1,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Vbound(i1,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
				    Wbound(i1,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
			   enddo   
			 enddo	
		 endif 		
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

	  if (botstress.eq.-2) then
       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         Ubound(i,j,k1)   = -Ubound(i,j,kmax)
         Vbound(i,j,k1)   = -Vbound(i,j,kmax)
         !Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1)   = 0.
         Ubound(i,j,0)    = -Ubound(i,j,1) 
         Vbound(i,j,0)    = -Vbound(i,j,1)
         Wbound(i,j,0)    = Wbed(i,j) !0.
         xx=Rp(i)*cos_u(j)-schuif_x
	     yy=Rp(i)*sin_u(j)
	Wbound(i,j,kmax)=(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo	
	  elseif (botstress.eq.-1) then
       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         Ubound(i,j,k1)   = Ubound(i,j,kmax)
         Vbound(i,j,k1)   = Vbound(i,j,kmax)
         !Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1) = 0.
         Ubound(i,j,0)    = -Ubound(i,j,1) 
         Vbound(i,j,0)    = -Vbound(i,j,1)
         Wbound(i,j,0)    = Wbed(i,j) !0.

         xx=Rp(i)*cos_u(j)-schuif_x
	     yy=Rp(i)*sin_u(j)
	     Wbound(i,j,kmax)=(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo	   
	  else
       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         Ubound(i,j,k1)   = Ubound(i,j,kmax)
         Vbound(i,j,k1)   = Vbound(i,j,kmax)
         !Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1)   = 0.
         Ubound(i,j,0)    = Ubound(i,j,1) 
         Vbound(i,j,0)    = Vbound(i,j,1)
         Wbound(i,j,0)    = Wbed(i,j) !0.

         xx=Rp(i)*cos_u(j)-schuif_x
	     yy=Rp(i)*sin_u(j)
	Wbound(i,j,kmax)=(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo
	  endif

      !! Set boundary conditions vertical jet in:
       	if (plumetseriesfile.eq.'') then
       		Wjet=W_j
       		f=Strouhal*ABS(W_j)/MAX(radius_j*2.,1.e-12) !Strouholt number is 0.3
	else
       		Wjet=interpseries(plumetseries,plumeUseries,plumeseriesloc,tt)
       		f=Strouhal*ABS(W_j)/MAX(radius_j*2.,1.e-12) !Strouholt number is 0.3
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
	  rr=1.-sqrt((xx**2+yy**2)/MAX(radius_j**2,1.e-12))
	  rr=MAX(rr,0.)
	  if (outflow_overflow_down.eq.1) then
		Wbound(i,j,kmax)=jetcorr*Wjet*rr**(1./W_j_powerlaw)+Awjet/REAL(azi_n)*Wjet*fluc !no turb SEM fluc
	  else
		Wbound(i,j,kmax)=jetcorr*Wjet*rr**(1./W_j_powerlaw)+Awjet/REAL(azi_n)*Wjet*fluc+Wb3in(i,j) !including turb SEM fluc
	  endif
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
	if (outflow_overflow_down.eq.1) then
		Ujetbc2=Ujetbc*cos(azi_angle_u(i,j))-Vjetbc*sin(azi_angle_u(i,j)) !no turb SEM fluc
		Vjetbc2=Ujetbc*sin(azi_angle_u(i,j))+Vjetbc*cos(azi_angle_u(i,j))
	else
		Ujetbc2=Ujetbc*cos(azi_angle_u(i,j))-Vjetbc*sin(azi_angle_u(i,j))+Ub3in(i,j) !including turb SEM fluc
		Vjetbc2=Ujetbc*sin(azi_angle_u(i,j))+Vjetbc*cos(azi_angle_u(i,j))+Vb3in(i,j)
	endif
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
	if (outflow_overflow_down.eq.1) then
		Ujetbc2=Ujetbc*cos(azi_angle_v(i,j))-Vjetbc*sin(azi_angle_v(i,j)) !no turb SEM fluc
		Vjetbc2=Ujetbc*sin(azi_angle_v(i,j))+Vjetbc*cos(azi_angle_v(i,j))
	else
		Ujetbc2=Ujetbc*cos(azi_angle_v(i,j))-Vjetbc*sin(azi_angle_v(i,j))+Ub3in(i,j) !including turb SEM fluc
		Vjetbc2=Ujetbc*sin(azi_angle_v(i,j))+Vjetbc*cos(azi_angle_v(i,j))+Vb3in(i,j) !including turb SEM fluc
	endif
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
	  rr=1.-sqrt((zzz**2+yy**2)/MAX(radius_j**2,1.e-12))
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
	  Wjetbc=pi*MAX(radius_j**2,1.e-12)*W_j*perc_dh_suction/(dr(i)*3.*Dsp) !(m3/s)/m2=m/s per suction cell (perc_dh_suction is corrected for number of dragheads)
	  Wbound(i,j,k)=Wjetbc
       enddo

		if (i_periodicx>0) then 
         do k=0,k1 
           do j=0,j1
		   Ubound(0,j,k)    =    Ubound(i_periodicx,j,k)
		   Vbound(0,j,k)    =    Vbound(i_periodicx,j,k)
		   Wbound(0,j,k)    =    Wbound(i_periodicx,j,k)
           enddo   
         enddo	
		endif 

      end

      subroutine bound_incljet(Ubound,Vbound,Wbound,rho,botstress,tt,Ub1in,Vb1in,Wb1in,Ub2in,Vb2in,Wb2in,Ub3in,Vb3in,Wb3in)
     
      USE nlist

      implicit none
c
!       include 'param.txt'
!       include 'common.txt'

c
      integer jtmp,t,inout,im,jm
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
	  real,intent(in) :: Ub3in(0:i1,0:j1),Vb3in(0:i1,0:j1),Wb3in(0:i1,0:j1)
	real x,y,z
	real Wbc1,Wbc2
	real Ujetbc,Vjetbc,Wjetbc,Ujetbc2,Vjetbc2
	real rr,interpseries
	real uu1,uu2,vv1,vv2
      real zzz,Ujet2,z2,val2,jetcorr
	  integer kp,kpp,kb,tel,jbeg,jend 
	  real zb_W,zb_U,zb_V,vel_ibm2,distance_to_bed_kp,distance_to_bed_kpp,yplus,absU,z0,ust

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
	    IF (periodicx.eq.3.or.periodicx.eq.4) THEN
			i=i1
		ELSE 
			i=0
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
	
!c get stuff from other CPU's
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
		   Vbound(i,j1,k) = Vbound(i,jmax,k)
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
		   Vbound(i,j1,k) = Vbound(i,jmax,k)
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
	if (periodicx.eq.0.and.monopile.eq.-1) then
		if (Uoutflow.eq.2) then ! convective outflow boundary
	      do k=1,kmax ! boundaries in i-direction
         	do j=0,j1
			Wbc2=0.
		    if (bcfile.ne.'') then
			Ubc2(j,k)=Ubcoarse2(j,k)
			Vbc2(j,k)=Vbcoarse2(j,k)
			Wbc2=Wbcoarse2(j,k)
		    endif
		    Ubc=Ub2(j,k)+Ubc2(j,k)
		    Vbc=0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k)
		    Wbc=0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2	
		   Ubound(0,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
		   Vbound(0,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
		   Wbound(0,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
		   Ubc = MAX(Ubc2(j,k)*cos_u(j)+Vbc2(j,k)*sin_u(j),0.)  !Unormal outflow boundary at U-loc; should be >0 otherwise no outflow and divide by zero 
!		   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k)-dr(imax)/(Ubc*dt)*(Unew(imax-1,j,k)-Uold(imax-1,j,k)),0.)
!		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
!		   Vbound(i1,j,k)   =    Vbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Vnew(imax,j,k)-Vold(imax,j,k))		   
!		   Wbound(i1,j,k)   =    Wbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Wnew(imax,j,k)-Wold(imax,j,k))
		   
		   ! solving Ui+1^n+1-Ui+1^n)/dt+Unormal(Ui+1^n+1-Ui^n+1)/dx=0. !implicit treatment of dUdn term
		   Ubound(imax,j,k) = (Unew(imax,j,k)+Ubc*dt/dr(imax)*Ubound(imax-1,j,k))/(1.+Ubc*dt/dr(imax)) 
		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
		   Vbound(i1,j,k) = (Vnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Vbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
		   Wbound(i1,j,k) = (Wnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Wbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 		   
         	enddo   
      	    enddo
		else !Neumann outflow boundary 
	      do k=1,kmax ! boundaries in i-direction
         	do j=0,j1
			Wbc2=0.
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
		endif 
       elseif (periodicx.eq.0.and.monopile>0) then ! closed boundary at i=0 (flow past circular cylinder); mixed inflow/outflow bc at imax; 3 = partial slip with tau-wall/4=free slip
		if (monopile.eq.4) then !no slip monopile wall 
         do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    0.
		   Vbound(0,j,k)    =    -Vbound(1,j,k)
		   Wbound(0,j,k)    =    -Wbound(1,j,k)
           enddo   
         enddo	
		 else 
		  do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    0.
		   Vbound(0,j,k)    =    Vbound(1,j,k)
		   Wbound(0,j,k)    =    Wbound(1,j,k)
           enddo   
          enddo	
		 endif 		 
		 if (rank.le.CEILING(REAL(px)/4.)-1) then ! first 1/4 of domain is outflow bc 
			jend = CEILING(REAL(jmax*px)*0.25) !end of first 1/4 global j-index
			jend = MIN(jend-rank*jmax,j1) !end of first 1/4 local j-index 
			if (Uoutflow.eq.2) then ! convective outflow boundary
			 do j=0,jend
			   do k=1,kmax 
				   Ubc = MAX(Ubc2(j,k)*cos_u(j)+Vbc2(j,k)*sin_u(j),0.)  !Unormal outflow boundary at U-loc; should be >0 otherwise no outflow and divide by zero 
!				   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k)-dr(imax)/(Ubc*dt)*(Unew(imax-1,j,k)-Uold(imax-1,j,k)),0.)
!				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
!				   Vbound(i1,j,k)   =    Vbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Vnew(imax,j,k)-Vold(imax,j,k))		   
!				   Wbound(i1,j,k)   =    Wbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Wnew(imax,j,k)-Wold(imax,j,k))
				   
				   ! solving Ui+1^n+1-Ui+1^n)/dt+Unormal(Ui+1^n+1-Ui^n+1)/dx=0. !implicit treatment of dUdn term
				   Ubound(imax,j,k) = (Unew(imax,j,k)+Ubc*dt/dr(imax)*Ubound(imax-1,j,k))/(1.+Ubc*dt/dr(imax)) 
				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
				   Vbound(i1,j,k) = (Vnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Vbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
				   Wbound(i1,j,k) = (Wnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Wbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 				   
			   enddo   
			 enddo		
			else !Neumann outflow boundary 
			 do j=0,jend
			   do k=1,kmax 
			     Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
			     Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
			     Vbound(i1,j,k)    =    Vbound(imax,j,k)
			     Wbound(i1,j,k)    =    Wbound(imax,j,k)
			   enddo   
			 enddo			
			endif 
			 do j=jend+1,j1  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
					Wbc2=0.
					if (bcfile.ne.'') then
						Ubc2(j,k)=Ubcoarse2(j,k) 
						Vbc2(j,k)=Vbcoarse2(j,k)
						Wbc2=Wbcoarse2(j,k)
					endif
					Ubc=Ub2(j,k)+Ubc2(j,k)
					Vbc=0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k)
					Wbc=0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2	
				    Ubound(imax,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Ubound(i1,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Vbound(i1,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
				    Wbound(i1,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
			   enddo   
			 enddo
		 elseif (rank.ge.px-CEILING(REAL(px)/4.)-1) then ! last 1/4 of domain is outflow bc 
			jbeg = jmax*px - CEILING(REAL(jmax*px)*0.25) ! begin of last 1/4 global j-index
			jbeg = MAX(jbeg-rank*jmax,0) !begin of last 1/4 local j-index
			if (Uoutflow.eq.2) then ! convective outflow boundary
			 do j=jbeg,j1 ! outflow 
			   do k=1,kmax 
				   Ubc = MAX(Ubc2(j,k)*cos_u(j)+Vbc2(j,k)*sin_u(j),0.)  !Unormal outflow boundary at U-loc; should be >0 otherwise no outflow and divide by zero 
!				   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k)-dr(imax)/(Ubc*dt)*(Unew(imax-1,j,k)-Uold(imax-1,j,k)),0.)
!				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
!				   Vbound(i1,j,k)   =    Vbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Vnew(imax,j,k)-Vold(imax,j,k))		   
!				   Wbound(i1,j,k)   =    Wbound(imax,j,k)-(Rp(i1)-Rp(imax))/(Ubc*dt)*(Wnew(imax,j,k)-Wold(imax,j,k))
				   
				   ! solving Ui+1^n+1-Ui+1^n)/dt+Unormal(Ui+1^n+1-Ui^n+1)/dx=0. !implicit treatment of dUdn term
				   Ubound(imax,j,k) = (Unew(imax,j,k)+Ubc*dt/dr(imax)*Ubound(imax-1,j,k))/(1.+Ubc*dt/dr(imax)) 
				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
				   Vbound(i1,j,k) = (Vnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Vbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
				   Wbound(i1,j,k) = (Wnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Wbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 				   
			   enddo   
			 enddo
			else !Neumann outflow boundary 
			 do j=jbeg,j1 ! outflow 
			   do k=1,kmax 
			     Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
			     Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
			     Vbound(i1,j,k)    =    Vbound(imax,j,k)
			     Wbound(i1,j,k)    =    Wbound(imax,j,k)
			   enddo   
			 enddo			
			endif 
			 do j=0,MIN(j1,jbeg-1)  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
					Wbc2=0.
					if (bcfile.ne.'') then
						Ubc2(j,k)=Ubcoarse2(j,k) 
						Vbc2(j,k)=Vbcoarse2(j,k)
						Wbc2=Wbcoarse2(j,k)
					endif
					Ubc=Ub2(j,k)+Ubc2(j,k)
					Vbc=0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k)
					Wbc=0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2	
				    Ubound(imax,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Ubound(i1,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Vbound(i1,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
				    Wbound(i1,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
			   enddo   
			 enddo			 
		 else !middle 2/4 is inflow bc 
			 do j=0,j1  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
					Wbc2=0.
					if (bcfile.ne.'') then
						Ubc2(j,k)=Ubcoarse2(j,k) 
						Vbc2(j,k)=Vbcoarse2(j,k)
						Wbc2=Wbcoarse2(j,k)
					endif
					Ubc=Ub2(j,k)+Ubc2(j,k)
					Vbc=0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k)
					Wbc=0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2	
				    Ubound(imax,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Ubound(i1,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Vbound(i1,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
				    Wbound(i1,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
			   enddo   
			 enddo	
		 endif 		
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
	
	  if (botstress.eq.-2) then
       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         Ubound(i,j,k1)   = -Ubound(i,j,kmax)
         Vbound(i,j,k1)   = -Vbound(i,j,kmax)
         !Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1) = 0.
         Ubound(i,j,0)    = -Ubound(i,j,1) 
         Vbound(i,j,0)    = -Vbound(i,j,1)
         Wbound(i,j,0)    = Wbed(i,j) !0.

         xx=Rp(i)*cos_u(j)-schuif_x
	     yy=Rp(i)*sin_u(j)
	     Wbound(i,j,kmax)=(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo	
	  elseif (botstress.eq.-1) then
       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         Ubound(i,j,k1)   = Ubound(i,j,kmax)
         Vbound(i,j,k1)   = Vbound(i,j,kmax)
         !Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1) = 0.
         Ubound(i,j,0)    = -Ubound(i,j,1) 
         Vbound(i,j,0)    = -Vbound(i,j,1)
         Wbound(i,j,0)    = Wbed(i,j) !0.

         xx=Rp(i)*cos_u(j)-schuif_x
	     yy=Rp(i)*sin_u(j)
	     Wbound(i,j,kmax)=(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo	   
	  else
       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         Ubound(i,j,k1)   = Ubound(i,j,kmax)
         Vbound(i,j,k1)   = Vbound(i,j,kmax)
         !Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1) = 0.
         Ubound(i,j,0)    = Ubound(i,j,1) 
         Vbound(i,j,0)    = Vbound(i,j,1)
         Wbound(i,j,0)    = Wbed(i,j) !0.

         xx=Rp(i)*cos_u(j)-schuif_x
	     yy=Rp(i)*sin_u(j)
	     Wbound(i,j,kmax)=(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo
	  endif

!! there used to be bedplume loop here to assign bp()u,v,w directly, but now with determination of force in bedplume in bound_rhoU this is no longer needed
	
	 
	IF (LOA>0.) THEN ! ship:
	  do t=1,tmax_inWpuntTSHD
 	    i=i_inWpuntTSHD(t)
 	    j=j_inWpuntTSHD(t)		
 	    k=k_inWpuntTSHD(t)		
		Wbound(i,j,k)=0.
	    !Wbound(i,j,k)=facIBMw(t)*Wbound(i,j,k) !0.
		!if (facIBMw(t)<1.) then
		!  Wbound(i,j,k)=facIBMw(t)/(1.+facIBMw(t))*Wbound(i,j,k+1)
		!endif		
	  enddo
	  do t=1,tmax_inUpuntTSHD
 	    i=i_inUpuntTSHD(t)
 	    j=j_inUpuntTSHD(t)		
 	    k=k_inUpuntTSHD(t)		
		Ubound(i,j,k)=0.
	    !Ubound(i,j,k)=facIBMu(t)*Ubound(i,j,k) !0.
		!if (facIBMu(t)<1.) then
		!  ubound(i,j,k)=facIBMu(t)/(1.+facIBMu(t))*ubound(i,j,k+1)
		!endif		
	  enddo
	  do t=1,tmax_inVpuntTSHD
 	    i=i_inVpuntTSHD(t)
 	    j=j_inVpuntTSHD(t)		
 	    k=k_inVpuntTSHD(t)		
		Vbound(i,j,k)=0.
	    !Vbound(i,j,k)=facIBMv(t)*Vbound(i,j,k) !0.
		!if (facIBMv(t)<1.) then
		!  vbound(i,j,k)=facIBMv(t)/(1.+facIBMv(t))*vbound(i,j,k+1)
		!endif		
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
	  IF (IBMorder.ne.2) THEN
	  do i=0,i1 ! prescribe UTSHD in obstacle:
	    do j=0,j1
	      do k=1,kbed(i,j) 
	        Ubound(i,j,k)=Ubot_TSHD(j)
	        Vbound(i,j,k)=Vbot_TSHD(j)
	      enddo
	    enddo
	  enddo
	  ENDIF 
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
       		f=Strouhal*ABS(W_j)/MAX(radius_j*2.,1.e-12) !Strouholt number is 0.3
	else
       		Wjet=interpseries(plumetseries,plumeUseries,plumeseriesloc,tt)
       		f=Strouhal*ABS(W_j)/MAX(radius_j*2.,1.e-12) !Strouholt number is 0.3
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
	rr=1.-sqrt((xx**2+yy**2)/MAX(radius_j**2,1.e-12))
	rr=MAX(rr,0.)
	if (outflow_overflow_down.eq.1) then
	 	do k=kmax-kjet,kmax-kjet
	 	  Wbound(i,j,k)=Wbound2(i,j,k)	
	 	enddo
	 	do k=kmax-kjet+1,kmax
	 	  Wbound(i,j,k)=(jetcorr*rr**(1./W_j_powerlaw)+Awjet/REAL(azi_n)*fluc)*Wjet+Wb3in(i,j)
	  	enddo
	else
	 	do k=kmax-kjet,kmax
	 	  Wbound(i,j,k)=Wbound2(i,j,k) !Wjet+Awjet/REAL(azi_n)*Wjet*fluc
	 	enddo
	 	Wbound(i,j,kmax)=(jetcorr*rr**(1./W_j_powerlaw)+Awjet/REAL(azi_n)*fluc)*Wjet+Wb3in(i,j)
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
	Ujetbc2=Ujetbc*cos(azi_angle_u(i,j))-Vjetbc*sin(azi_angle_u(i,j))+Ub3in(i,j)
	Vjetbc2=Ujetbc*sin(azi_angle_u(i,j))+Vjetbc*cos(azi_angle_u(i,j))+Vb3in(i,j)
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
	Ujetbc2=Ujetbc*cos(azi_angle_v(i,j))-Vjetbc*sin(azi_angle_v(i,j))+Ub3in(i,j)
	Vjetbc2=Ujetbc*sin(azi_angle_v(i,j))+Vjetbc*cos(azi_angle_v(i,j))+Vb3in(i,j)
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
	  rr=1.-sqrt((zzz**2+yy**2)/MAX(radius_j**2,1.e-12))
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
	  Wjetbc=pi*MAX(radius_j**2,1.e-12)*W_j*perc_dh_suction/(dr(i)*3.*Dsp) !(m3/s)/m2=m/s per suction cell (perc_dh_suction is corrected for number of dragheads)
	  Wbound(i,j,k)=Wjetbc
       enddo
	IF (interaction_bed.ge.4) THEN
	 DO i=1,imax
	  DO j=1,jmax
		zbed(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
	  ENDDO 
	 ENDDO 
	 call bound_cbot(zbed)
	ENDIF 	
	IF ((interaction_bed.ge.4.or.bedlevelfile.ne.''.or.nobst>0).and.IBMorder.eq.2) THEN ! order-2 IBM before information is exchanged between partitions (hence only j=1-jmax); order-0 IBM is done at and of this subroutine 
		DO i=1,imax
			DO j=1,jmax
				zb_W=zbed(i,j)
				kb=FLOOR(zb_W/dz)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
				kp=MIN(CEILING(zb_W/dz),kmax)			!location velocity which must be adjusted 2nd order IBM
				kpp=MIN(CEILING(zb_W/dz)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity
				distance_to_bed_kp=REAL(kp)*dz-zb_W
				distance_to_bed_kpp=REAL(kpp)*dz-zb_W
				vel_ibm2=distance_to_bed_kp/distance_to_bed_kpp*Wbound(i,j,kpp) !linear interpolation with zero velocity at bed	
				!vel_ibm2=distance_to_bed_kp/distance_to_bed_kpp*Wbound(i,j,kp) !apply 'force' to W(kp) to make zero when bed approaches location kp				
				  DO k=1,kb
					Wbound(i,j,k)=0.
				  ENDDO						
				Wbound(i,j,kp)=vel_ibm2 			

				zb_U=0.5*(zbed(i,j)+zbed(i+1,j))
				kb=FLOOR(zb_U/dz+0.5)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
				DO k=1,kb
					Ubound(i,j,k)=0.
				ENDDO
				kp=MIN(CEILING(zb_U/dz+0.5),kmax)			!location velocity which must be adjusted 2nd order IBM -->(0-1)*dz distance from bed
				distance_to_bed_kp=(REAL(kp)-0.5)*dz-zb_U				  
				IF (dUVdn_IBMbed.eq.0.and.distance_to_bed_kp>0.1*dz) THEN
					Ubound(i,j,kb) = Ubound(i,j,kp)
				ENDIF				  
				IF (slip_bot.ge.1) THEN	! tests showed it is best to apply tau shear stress and not prescribe U,V velocity straightaway including influence tau (latter gives too large near bed velocities)							  
					! tau is only applied in bound_rhoU				
				ELSE ! without partial slip bc, simply prescribe 2nd order accurate velocity:		
					kp=MIN(CEILING(zb_U/dz+0.5),kmax)			!location velocity which must be adjusted 2nd order IBM
					kpp=MIN(CEILING(zb_U/dz+0.5)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity
					distance_to_bed_kp=(REAL(kp)-0.5)*dz-zb_U
					distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_U	
					Ubound(i,j,kp)=distance_to_bed_kp/distance_to_bed_kpp*Ubound(i,j,kpp) !linear interpolation with zero velocity at bed			
				ENDIF 

				zb_V=0.5*(zbed(i,j)+zbed(i,j+1))
				kb=FLOOR(zb_V/dz+0.5)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
				DO k=1,kb
					Vbound(i,j,k)=0.
				ENDDO
				kp=MIN(CEILING(zb_V/dz+0.5),kmax)			!location velocity which must be adjusted 2nd order IBM
				distance_to_bed_kp=(REAL(kp)-0.5)*dz-zb_V
				IF (dUVdn_IBMbed.eq.0.and.distance_to_bed_kp>0.1*dz) THEN
					Vbound(i,j,kb) = Vbound(i,j,kp)
				ENDIF					  
				IF (slip_bot.ge.1) THEN				
					! tau is only applied in bound_rhoU
				ELSE ! without partial slip bc, simply prescribe 2nd order accurate velocity:				  
					kp=MIN(CEILING(zb_V/dz+0.5),kmax)			!location velocity which must be adjusted 2nd order IBM
					kpp=MIN(CEILING(zb_V/dz+0.5)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity
					distance_to_bed_kp=(REAL(kp)-0.5)*dz-zb_V
					distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_V						
					Vbound(i,j,kp)=distance_to_bed_kp/distance_to_bed_kpp*Vbound(i,j,kpp) !linear interpolation with zero velocity at bed
				ENDIF
			ENDDO
		ENDDO
	ENDIF	
	
	IF (IBMorder.ne.2.or.dUVdn_IBMbed.eq.-2) THEN !19-2-2021 apply for IBMorder.eq.2 to prevent SSC issues near bed when dUVdn_IBMbed=-2
		IF (interaction_bed.ge.4.or.bedlevelfile.ne.''.or.nobst>0) THEN ! dynamic bed update, IBM bed re-defined every timestep: !order-0 IBM: make vel zero for all cells within bed
			DO i=0,i1 !1,imax
				DO j=0,j1 !1,jmax
					DO k=1,kbed(i,j)
						Ubound(i,j,k)=Ubot_TSHD(j) !0.
						Vbound(i,j,k)=Vbot_TSHD(j) !0.
						im=MAX(i-1,0)
						jm=MAX(j-1,0)
						Ubound(im,j,k)=Ubot_TSHD(j) !0.
						Vbound(i,jm,k)=Vbot_TSHD(jm) !0.
						Wbound(i,j,k)=0.
					ENDDO 
				ENDDO
			ENDDO
			IF (dUVdn_IBMbed.eq.0) THEN 
				DO i=1,imax
					DO j=1,jmax
						k=kbed(i,j) 
						IF (kbed(i+1,j).eq.k) Ubound(i,j,k)=Ubound(i,j,MIN(k+1,k1))
						IF (kbed(i,j+1).eq.k) Vbound(i,j,k)=Vbound(i,j,MIN(k+1,k1))
						im=MAX(i-1,0)
						jm=MAX(j-1,0)
						IF (kbed(i-1,j).eq.k) Ubound(im,j,k)=Ubound(im,j,MIN(k+1,k1))
						IF (kbed(i,j-1).eq.k) Vbound(i,jm,k)=Vbound(i,jm,MIN(k+1,k1))				
					ENDDO
				ENDDO	 
			ENDIF 
		ENDIF		
	ENDIF 
	
	! apply j-boundary conditions for U,V,W again for application tau and ibm2 on j=1:jmax and not 0:j1 + 	
	! to fix small inconsistency in Vbound(:,j1,:) made in jet, jet2 and rudder
	call shiftf(Vbound,vbf) 
	call shiftb(Vbound,vbb) 
	call shiftf(Ubound,ubf) 
	call shiftb(Ubound,ubb) 
	call shiftf(Wbound,wbf) 
	call shiftb(Wbound,wbb) 	
	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then ! boundaries in j-direction

	  elseif (rank.eq.px-1) then
	
	  else
		do k=0,k1
		   do i=1,imax
		   Vbound(i,0,k) = Vbf(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   Ubound(i,0,k) = Ubf(i,k)
		   Ubound(i,j1,k) =Ubb(i,k)	
		   Wbound(i,0,k) = Wbf(i,k)
		   Wbound(i,j1,k) =Wbb(i,k)			   
		   enddo
		enddo
	  endif
	else !periodic in y:
	  do k=0,k1
	   do i=1,imax
	     Vbound(i,0,k) = Vbf(i,k)
	     Vbound(i,j1,k) =Vbb(i,k)
	     Ubound(i,0,k) = Ubf(i,k)
	     Ubound(i,j1,k) =Ubb(i,k)	
	     Wbound(i,0,k) = Wbf(i,k)
	     Wbound(i,j1,k) =Wbb(i,k)			 
	   enddo
	  enddo
	endif
 
		if (i_periodicx>0) then 
         do k=0,k1 
           do j=0,j1
		   Ubound(0,j,k)    =    Ubound(i_periodicx,j,k)
		   Vbound(0,j,k)    =    Vbound(i_periodicx,j,k)
		   Wbound(0,j,k)    =    Wbound(i_periodicx,j,k)
           enddo   
         enddo	
		endif 
	

      end

      subroutine bound_rhoU(Ubound,Vbound,Wbound,rho,botstress,mpstress,tt,Ub1in,Vb1in,Wb1in,Ub2in,Vb2in,Wb2in,Ub3in,Vb3in,Wb3in)

      USE nlist

      implicit none
c
!       include 'param.txt'
!       include 'common.txt'

c
      integer jtmp,botstress,n,t,jp,inout,im,jm,kplus,kplus2,tel,mpstress
	real xTSHD(1:4),yTSHD(1:4)
c
      real  Ubound(0:i1,0:j1,0:k1),Vbound(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +      Wbound(0:i1,0:j1,0:k1)
      real ubb(0:i1,0:k1),val,theta,Ubc,Vbc,Wbc,xx,yy,Ujet,Vjet,Wjet,theta_U,theta_V
      real vbb(0:i1,0:k1)
      real wbb(0:i1,0:k1),rbb(0:i1,0:k1)
      real ubf(0:i1,0:k1)
      real vbf(0:i1,0:k1)
      real wbf(0:i1,0:k1),rbf(0:i1,0:k1)
      real cbf(0:i1,0:k1)
      real cbb(0:i1,0:k1)
	real ust_U_b,ust_V_b,Chezy,fluc,f,tt,z0_U,z0_V,phi,z,x,y
      real Ub1(0:i1,0:k1),Vb1(0:i1,0:k1),Wb1(0:i1,0:k1),Ub2(0:j1,0:k1),Vb2(0:j1+1,0:k1),Wb2(0:j1,0:k1)
      real,intent(in):: Ub1in(0:i1,0:k1),Vb1in(0:i1,0:k1),Wb1in(0:i1,0:k1),Ub2in(0:j1,0:k1),Vb2in(0:j1+1,0:k1),Wb2in(0:j1,0:k1)
	  real,intent(in) :: Ub3in(0:i1,0:j1),Vb3in(0:i1,0:j1),Wb3in(0:i1,0:j1)
	real Wbc1,Wbc2
	real Ujetbc,Vjetbc,Wjetbc,Ujetbc2,Vjetbc2
	real rr,interpseries
	real uu1,uu2,vv1,vv2,test
      real zzz,Ujet2,z2,val2,jetcorr
	  real duu,du,Uav,Tadapt,fbx2,fby2,fbz2,fbx,fby,ppp(0:i1,0:j1,0:k1),scal

!      real  Ubound2(0:i1,0:j1,0:k1),Vbound2(0:i1,0:j1,0:k1),Wbound2(0:i1,0:j1,0:k1)
      real  Ubound2(-1:i1+1,-1:j1+1,0:k1),Vbound2(-1:i1+1,-1:j1+1,0:k1),Wbound2(-1:i1+1,-1:j1+1,0:k1),rho2(-1:i1+1,-1:j1+1,0:k1),rr1,rr2
	real Propx_dummy(0:i1,0:px*jmax+1,1:kmax)
	real Propy_dummy(0:i1,0:px*jmax+1,1:kmax)
	real Propz_dummy(0:i1,0:px*jmax+1,1:kmax)	  
	  integer kp,kpp,kb,jbeg,jend
	  real zb_W,zb_U,zb_V,vel_ibm2,distance_to_bed_kp,distance_to_bed_kpp,yplus,absU,z0,ust,rr2,ww1,ww2
	  real dpdz1,dpdy2,tauw1,tauv1,tauw2,tauv2,ust1,ust2,tauu1,tauu2,dpdx11,dpdy22,tau

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

	
	 if (split_rho_cont.eq.'VL2'.or.split_rho_cont.eq.'SB2') then
	    if (split_rho_cont.eq.'VL2') THEN
		 do n=1,nfrac
	      call c_edges_VL2_nocfl(cU(n,:,:,:),cV(n,:,:,:),cW(n,:,:,:),dcdt(n,:,:,:),Ubound,Vbound,Wbound,rho,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,1,imax,1,jmax,1,kmax,dt,rank,px,periodicx,periodicy) 
		 enddo
		endif
	    if (split_rho_cont.eq.'SB2') THEN
		 do n=1,nfrac
	      call c_edges_SB2_nocfl(cU(n,:,:,:),cV(n,:,:,:),cW(n,:,:,:),dcdt(n,:,:,:),Ubound,Vbound,Wbound,rho,Ru,Rp,dr,phiv,phipt,dz,
     +            i1,j1,k1,1,imax,1,jmax,1,kmax,dt,rank,px,periodicx,periodicy) 
		 enddo
		endif	
	 ! bound_rhoU always called for n+1 timestep which always corresponds to dcdt (but this is not passed explicitly to this subroutine)
	 ! this is the first moment in a timestep that c_edges_TVD is called in the code		
		call state_edges(cU,rhU)
		call state_edges(cV,rhV)
		call state_edges(cW,rhW)
		if (nobst>0.or.bedlevelfile.ne.''.or.interaction_bed.ge.4) then ! default C(k=0)=C(k=1); therefore only for cases when bed is not necessarily at k=0 this fix is needed:
		 DO i=0,i1
		  DO j=0,j1
			kplus = MIN(kbed(i,j)+1,k1)
			kplus2 = MIN(kbed(i,j)+2,k1)
			rhW(i,j,kbed(i,j))=rho(i,j,kplus) ! make rhW(kbed) equal to rho fluid first cell above to get correct drift flux settling
			!rhW(i,j,kbed(i,j))=1.5*rho(i,j,kplus)-0.5*rho(i,j,kplus2)
			DO n=1,nfrac
				cW(n,i,j,kbed(i,j))=dcdt(n,i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
				!cW(n,i,j,kbed(i,j))=1.5*dcdt(n,i,j,kplus)-0.5*dcdt(n,i,j,kplus2)
			ENDDO
		  ENDDO
		 ENDDO
		endif		
	 else 
		cU(1:nfrac,0:imax,0:j1,0:k1)=0.5*(dcdt(1:nfrac,0:imax,0:j1,0:k1)+dcdt(1:nfrac,1:imax+1,0:j1,0:k1)) !er stond hier cnew ipv dcdt, maar rest gebruikt rho==drdt <--> dcdt
		cV(1:nfrac,0:i1,0:jmax,0:k1)=0.5*(dcdt(1:nfrac,0:i1,0:jmax,0:k1)+dcdt(1:nfrac,0:i1,1:jmax+1,0:k1)) !lijkt typo geweest te zijn (hierboven 'VL2' en 'SB2' wordt ook dcdt gebruikt)
		IF (hindered_settling_c.eq.1) THEN 
			cW(1:nfrac,0:i1,0:j1,0:kmax)=MAX(dcdt(1:nfrac,0:i1,0:j1,0:kmax),dcdt(1:nfrac,0:i1,0:j1,1:kmax+1))	 
		ELSE 
			cW(1:nfrac,0:i1,0:j1,0:kmax)=0.5*(dcdt(1:nfrac,0:i1,0:j1,0:kmax)+dcdt(1:nfrac,0:i1,0:j1,1:kmax+1))	 
		ENDIF 
		rhU(0:imax,0:j1,0:k1)=0.5*(rho(0:imax,0:j1,0:k1)+rho(1:imax+1,0:j1,0:k1))
		rhV(0:i1,0:jmax,0:k1)=0.5*(rho(0:i1,0:jmax,0:k1)+rho(0:i1,1:jmax+1,0:k1))
		rhW(0:i1,0:j1,0:kmax)=0.5*(rho(0:i1,0:j1,0:kmax)+rho(0:i1,0:j1,1:kmax+1))
      
		if (nobst>0.or.bedlevelfile.ne.''.or.interaction_bed.ge.4) then ! default C(k=0)=C(k=1); therefore only for cases when bed is not necessarily at k=0 this fix is needed:
		 DO i=0,i1
		  DO j=0,j1
			kplus = MIN(kbed(i,j)+1,k1)
			kplus2 = MIN(kbed(i,j)+2,k1)
			rhW(i,j,kbed(i,j))=rho(i,j,kplus) ! make rhW(kbed) equal to rho fluid first cell above to get correct drift flux settling
			!rhW(i,j,kbed(i,j))=1.5*rho(i,j,kplus)-0.5*rho(i,j,kplus2)
			DO n=1,nfrac
				cW(n,i,j,kbed(i,j))=dcdt(n,i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
				!cW(n,i,j,kbed(i,j))=1.5*dcdt(n,i,j,kplus)-0.5*dcdt(n,i,j,kplus2)
			ENDDO
		  ENDDO
		 ENDDO
		endif
	 IF (applyVOF.eq.1) THEN !with apply_VOF=1 density on edges should not be real density but constant 
		rhU=rho_b
		rhV=rho_b
		rhW=rho_b
	 endif
	 
	 
	 if (periodicx.eq.1) then
		rhU(i1,0:j1,0:k1)=  rhU(1,0:j1,0:k1)
		rhV(i1,0:j1,0:k1)=  rhV(1,0:j1,0:k1)
		rhW(i1,0:j1,0:k1)=  rhW(1,0:j1,0:k1)
		rhU(0,0:j1,0:k1)=  rhU(imax,0:j1,0:k1)
		rhV(0,0:j1,0:k1)=  rhV(imax,0:j1,0:k1)
		rhW(0,0:j1,0:k1)=  rhW(imax,0:j1,0:k1)		
	 else
		rhU(i1,0:j1,0:k1)=  rho(i1,0:j1,0:k1) 
		rhU(imax,0:j1,0:k1)=  	 0.5*(rho(imax,0:j1,0:k1)+rho(imax+1,0:j1,0:k1))
		rhV(i1,0:jmax,0:k1)=  0.5*(rho(i1,0:jmax,0:k1)+rho(i1,1:j1,0:k1))
		rhW(i1,0:j1,0:kmax)=  0.5*(rho(i1,0:j1,0:kmax)+rho(i1,0:j1,1:k1))
		rhU(0,0:j1,0:k1)=  	 0.5*(rho(0,0:j1,0:k1)+rho(1,0:j1,0:k1))
		rhV(0,0:jmax,0:k1)=  0.5*(rho(0,0:jmax,0:k1)+rho(0,1:j1,0:k1))
		rhW(0,0:j1,0:kmax)=  0.5*(rho(0,0:j1,0:kmax)+rho(0,0:j1,1:k1))		
	 endif
	 rhU(0:imax,0:j1,k1)=  0.5*(rho(0:imax,0:j1,k1)+rho(1:i1,0:j1,k1))
	 rhV(0:i1,0:jmax,k1)=  0.5*(rho(0:i1,0:jmax,k1)+rho(0:i1,1:j1,k1))
	 rhW(0:i1,0:j1,k1)=    rho(0:i1,0:j1,k1)
	 rhW(0:i1,0:j1,kmax)=    0.5*(rho(0:i1,0:j1,kmax)+rho(0:i1,0:j1,kmax+1))
	 rhU(0:imax,0:j1,0)=  0.5*(rho(0:imax,0:j1,0)+rho(1:i1,0:j1,0))
	 rhV(0:i1,0:jmax,0)=  0.5*(rho(0:i1,0:jmax,0)+rho(0:i1,1:j1,0))
	 rhW(0:i1,0:j1,0)=    0.5*(rho(0:i1,0:j1,0)+rho(0:i1,0:j1,1))
	
		!! boundary conditions in y-dir as last to prevail (especially for internal partitions this is essential)
!c get stuff from other CPU's
	 call shiftf(rhU,ubf) 
	 call shiftf(rhV,vbf)
	 call shiftf(rhW,wbf)
	 call shiftb(rhU,ubb)
	 call shiftb(rhV,vbb)
	 call shiftb(rhW,wbb)
	 if (periodicy.eq.0.or.periodicy.eq.2) then !for rho it does not matter periodicy is 0 or 2
	  if (rank.eq.0) then ! boundaries in j-direction
		do k=0,k1
		   do i=0,i1
		   rhU(i,0,k) =0.5*(rho(i,0,k)+rho(MIN(i+1,i1),0,k))		   
		   rhV(i,0,k) =0.5*(rho(i,0,k)+rho(i,1,k))
		   rhW(i,0,k) =0.5*(rho(i,0,k)+rho(i,0,MIN(k+1,k1)))
		   rhU(i,j1,k) =Ubb(i,k)
		   rhV(i,j1,k) =Vbb(i,k)
		   rhW(i,j1,k) =Wbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1
		   rhU(i,j1,k) =0.5*(rho(i,j1,k)+rho(MIN(i+1,i1),j1,k))		   
		   rhV(i,j1,k) =rho(i,j1,k)
		   rhV(i,jmax,k) =0.5*(rho(i,jmax,k)+rho(i,jmax+1,k))
		   rhW(i,j1,k) =0.5*(rho(i,j1,k)+rho(i,j1,MIN(k+1,k1)))
		   rhU(i,0,k) =Ubf(i,k)
		   rhV(i,0,k) =Vbf(i,k)
		   rhW(i,0,k) =Wbf(i,k)
		   enddo
		enddo	
	  else
		do k=0,k1
		   do i=0,i1
		   rhU(i,0,k) = Ubf(i,k)
		   rhV(i,0,k) = Vbf(i,k)
		   rhW(i,0,k) = Wbf(i,k)
		   rhU(i,j1,k) =Ubb(i,k)
		   rhV(i,j1,k) =Vbb(i,k)
		   rhW(i,j1,k) =Wbb(i,k)
		   enddo
		enddo
	  endif
	 elseif (periodicy.eq.1) then !periodic in y:
		do k=0,k1
		   do i=0,i1
		   rhU(i,0,k) = Ubf(i,k)
		   rhV(i,0,k) = Vbf(i,k)
		   rhW(i,0,k) = Wbf(i,k)
		   rhU(i,j1,k) =Ubb(i,k)
		   rhV(i,j1,k) =Vbb(i,k)
		   rhW(i,j1,k) =Wbb(i,k)
	   enddo
	  enddo
	 endif
		 DO i=0,i1
		  DO j=0,j1
		   DO k=1,kbed(i,j) !prevent source term at immersed boundary in PPE --> rhW not adapted because dpdn=0 has been enforced in bound_p already
		    rhU(i,j,k)=rho_b 
			rhU(MAX(0,i-1),j,k)=rho_b 
		    rhV(i,j,k)=rho_b 
			rhV(i,MAX(0,j-1),k)=rho_b 
		  ENDDO
		 ENDDO
		ENDDO  
	ENDIF 
	
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
	    IF (periodicx.eq.3.or.periodicx.eq.4) THEN
			i=i1
		ELSE 
			i=0
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

	IF (botstress.eq.3.or.mpstress.eq.3) THEN 
	  IF (pres_in_predictor_step.eq.0) THEN 
		ppp(1:imax,1:jmax,1:kmax)=p 
	  ELSE 
		ppp(1:imax,1:jmax,1:kmax)=pold
	  ENDIF 
	  call bound_p(ppp)
	ENDIF

!	IF ((interaction_bed.ge.4.or.bedlevelfile.ne.''.or.nobst>0).and.IBMorder.eq.2) THEN ! order-2 IBM before information is exchanged between partitions (hence only j=1-jmax); order-0 IBM is done at and of this subroutine 
!		DO i=1,imax
!			DO j=1,jmax
!				IF (interaction_bed.ge.4) THEN
!				  !zb_W=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(cnewbot(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
!				  zb_W=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
!				ELSE 
!				  zb_W=zbed(i,j)
!				ENDIF
!				kb=FLOOR(zb_W/dz)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
!				kp=MIN(CEILING(zb_W/dz),kmax)			!location velocity which must be adjusted 2nd order IBM
!				kpp=MIN(CEILING(zb_W/dz)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity
!				distance_to_bed_kp=REAL(kp)*dz-zb_W
!				distance_to_bed_kpp=REAL(kpp)*dz-zb_W
!				vel_ibm2=distance_to_bed_kp/distance_to_bed_kpp*Wbound(i,j,kpp) !linear interpolation with zero velocity at bed		
!				DO k=1,kb
!				  Wbound(i,j,k)=0.
!				ENDDO					
!				Wbound(i,j,kp)=vel_ibm2 
!				IF (interaction_bed.ge.4) THEN
!!				  zb_U=0.5*(zb_W+REAL(MAX(kbed(i+1,j)-1,0))*dz+(SUM(cnewbot(1:nfrac,i+1,j))+
!!     &				SUM(Clivebed(1:nfrac,i+1,j,kbed(i+1,j))))/cfixedbed*dz)
!				  zb_U=0.5*(zb_W+REAL(MAX(kbed(i+1,j)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i+1,j))+
!     &				SUM(Clivebed(1:nfrac,i+1,j,kbed(i+1,j))))/cfixedbed*dz)	 
!				ELSE 
!				  zb_U=0.5*(zbed(i,j)+zbed(i+1,j))
!				ENDIF
!				kb=FLOOR(zb_U/dz+0.5)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
!				kp=MIN(CEILING(zb_U/dz+0.5),kmax)			!location velocity which must be adjusted 2nd order IBM -->(0-1)*dz distance from bed
!				kpp=MIN(CEILING(zb_U/dz+0.5)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity -->(1-2)*dz distance from bed
!				distance_to_bed_kp=(REAL(kp)-0.5)*dz-zb_U
!				distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_U				
!				IF (distance_to_bed_kp>0.1*dz) THEN !apply tau on first cell (0.1-1)*dz distance from bed based on velocity of that cell
!					kpp=kp
!					distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_U
!				ELSE !apply tau on second cell (1-1.1)*dz distance from bed based on velocity of that cell, do nothing for cell below: (0-0.1)dz from bed 
!					kp=kpp 
!!				  ELSE !apply tau on first cell (0-0.5)*dz distance from bed based on ustar of second cell (kpp) (1-1.5)*dz distance from bed 
!				ENDIF				
!				IF (botstress.eq.1.or.botstress.eq.2) THEN	! tests showed it is best to apply tau shear stress at first cell above ibm bed and not prescribe U,V velocity straightaway including influence tau (latter gives too large near bed velocities)			
!				  absU=sqrt(Ubound(i,j,kpp)**2+(0.25*(Vbound(i,j,kpp)+Vbound(i,j-1,kpp)+Vbound(i+1,j,kpp)+Vbound(i+1,j-1,kpp)))**2)
!				  absU=absU/rhU(i,j,kpp)
!				  ust=0.1*absU
!				  if (kn>0.) then !walls with rougness (log-law used which can be hydr smooth or hydr rough):
!					do tel=1,10 ! 10 iter is more than enough
!						z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
!						ust=absU/MAX(1./kappa*log(MAX(distance_to_bed_kpp/z0,1.001)),2.) !ust maximal 0.5*absU
!					enddo
!					!vel_ibm2=Ubound(i,j,kpp)/MAX(absU,1e-9)/rhU(i,j,kpp)* (ust/kappa*log(MAX(distance_to_bed_kp/z0,1.001)))
!				  else
!					do tel=1,10 ! 10 iter is more than enough
!						yplus=MAX(distance_to_bed_kpp*ust/nu_mol,1e-12)
!						ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
!					enddo
!					!vel_ibm2=Ubound(i,j,kpp)/MAX(absU,1e-9)/rhU(i,j,kpp)* ust*(2.5*log(MAX(distance_to_bed_kp*ust/nu_mol,1.001))+5.5)
!					if (yplus<30.) then
!					  do tel=1,10 ! 10 iter is more than enough
!						yplus=MAX(distance_to_bed_kpp*ust/nu_mol,1e-12)
!						ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
!					  enddo	
!					  !vel_ibm2=Ubound(i,j,kpp)/MAX(absU,1e-9)/rhU(i,j,kpp)* ust*(5*log(MAX(distance_to_bed_kp*ust/nu_mol,1.001))-3.05)
!					endif
!					if (yplus<5.) then !viscous sublayer uplus=yplus
!						ust=sqrt(absU*nu_mol/(distance_to_bed_kpp))
!						!vel_ibm2=Ubound(i,j,kpp)/MAX(absU,1e-9)/rhU(i,j,kpp)* ust*ust*distance_to_bed_kp/nu_mol
!					endif
!				  endif
!				  tau_fl_Unew(i,j) = rhU(i,j,kp)*ust*ust
!				  DO k=1,kb
!					Ubound(i,j,k)=0.
!				  ENDDO
!				  !Ubound(i,j,kp)=vel_ibm2*rhU(i,j,kp)	
!				  absU=sqrt(Ubound(i,j,kp)**2+(0.25*(Vbound(i,j,kp)+Vbound(i,j-1,kp)+Vbound(i+1,j,kp)+Vbound(i+1,j-1,kp)))**2)
!				  absU=absU/rhU(i,j,kp)				
!				  Ubound(i,j,kp) = Ubound(i,j,kp) / (1. + ust*ust*dt/dz/MAX(absU,1.e-9))  ! implicit = more stable	
!				  IF (dUVdn_IBMbed.eq.0.and.distance_to_bed_kp>0.1*dz) THEN
!					Ubound(i,j,kb) = Ubound(i,j,kp)
!				  ENDIF
!				ELSEIF(botstress.eq.3) THEN 
!				  rr1=rhU(i,j,kp) 																				!at U-gridpoint
!				  dpdx11=(ppp(i+1,j,kp)-ppp(i,j,kp))/(Rp(i+1)-Rp(i))/rr1										!at U-gridpoint
!				  tauu1=tau_fl_Uold(i,j)																		!at U-gridpoint
!				  tauv1=0.25*(tau_fl_Vold(i,j)+tau_fl_Vold(i+1,j)+tau_fl_Vold(i,j-1)+tau_fl_Vold(i+1,j-1))		!at U-gridpoint
!				  uu1=Ubound(i,j,kp)
!				  vv1=0.25*(Vbound(i,j,kp)+Vbound(i+1,j,kp)+Vbound(i,j-1,kp)+Vbound(i+1,j-1,kp))
!				  ust1 = ((tauu1/rr1)**2+(tauv1/rr1)**2)**0.25													!at U-gridpoint
!				  call wall_fun_rho_TBL(uu1,vv1,dpdx11,rr1,2.*distance_to_bed_kp,1.,dt,kn,kappa,nu_mol,ust1,tau_fl_Unew(i,j))
!				  Ubound(i,j,kp) = Ubound(i,j,kp) / (1.+tau_fl_Unew(i,j)*dt/dz/MAX(ABS(Ubound(i,j,kp)),1.e-9))  
!				  ! don't use uu1 directly as dz and distance_to_bed_kp are not equal
!				  IF (dUVdn_IBMbed.eq.0.and.distance_to_bed_kp>0.1*dz) THEN
!					Ubound(i,j,kb) = Ubound(i,j,kp)
!				  ENDIF				  
!				ELSE ! without partial slip bc, simply prescribe 2nd order accurate velocity:
!				  DO k=1,kb
!					Ubound(i,j,k)=0.
!				  ENDDO
!				  Ubound(i,j,kp)=distance_to_bed_kp/distance_to_bed_kpp*Ubound(i,j,kpp) !linear interpolation with zero velocity at bed			
!				ENDIF 
!	 			IF (interaction_bed.ge.4) THEN
!!				  zb_V=0.5*(zb_W+REAL(MAX(kbed(i,j+1)-1,0))*dz+(SUM(cnewbot(1:nfrac,i,j+1))+
!!     &				SUM(Clivebed(1:nfrac,i,j+1,kbed(i,j+1))))/cfixedbed*dz)
!				  zb_V=0.5*(zb_W+REAL(MAX(kbed(i,j+1)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i,j+1))+
!     &				SUM(Clivebed(1:nfrac,i,j+1,kbed(i,j+1))))/cfixedbed*dz)	 
!				ELSE 
!				  zb_V=0.5*(zbed(i,j)+zbed(i,j+1))
!				ENDIF	
!				
!				kb=FLOOR(zb_V/dz+0.5)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
!				kp=MIN(CEILING(zb_V/dz+0.5),kmax)			!location velocity which must be adjusted 2nd order IBM
!				kpp=MIN(CEILING(zb_V/dz+0.5)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity
!				distance_to_bed_kp=(REAL(kp)-0.5)*dz-zb_V
!				distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_V
!				IF (distance_to_bed_kp>0.1*dz) THEN !apply tau on first cell (0.1-1)*dz distance from bed based on velocity of that cell
!					kpp=kp
!					distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_V
!				ELSE !apply tau on second cell (1-1.1)*dz distance from bed based on velocity of that cell, do nothing for cell below: (0-0.1)dz from bed 
!					kp=kpp 
!!				  ELSE !apply tau on first cell (0-0.1)*dz distance from bed based on ustar of second cell (kpp) (1-1.1)*dz distance from bed 
!				ENDIF				
!				IF (botstress.eq.1.or.botstress.eq.2) THEN	
!				  absU=sqrt(Vbound(i,j,kpp)**2+(0.25*(Ubound(i,j,kpp)+Ubound(i,j+1,kpp)+Ubound(i-1,j,kpp)+Ubound(i-1,j+1,kpp)))**2)
!				  absU=absU/rhV(i,j,kpp)
!				  ust=0.1*absU
!				  if (kn>0.) then !walls with rougness (log-law used which can be hydr smooth or hydr rough):
!					do tel=1,10 ! 10 iter is more than enough
!						z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
!						ust=absU/MAX(1./kappa*log(MAX(distance_to_bed_kpp/z0,1.001)),2.) !ust maximal 0.5*absU
!					enddo
!					!vel_ibm2=Vbound(i,j,kpp)/MAX(absU,1e-9)/rhV(i,j,kpp)* (ust/kappa*log(MAX(distance_to_bed_kp/z0,1.001)))
!				  else
!					do tel=1,10 ! 10 iter is more than enough
!						yplus=MAX(distance_to_bed_kpp*ust/nu_mol,1e-12)
!						ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
!					enddo
!					!vel_ibm2=Vbound(i,j,kpp)/MAX(absU,1e-9)/rhV(i,j,kpp)* ust*(2.5*log(MAX(distance_to_bed_kp*ust/nu_mol,1.001))+5.5)
!					if (yplus<30.) then
!					  do tel=1,10 ! 10 iter is more than enough
!						yplus=MAX(distance_to_bed_kpp*ust/nu_mol,1e-12)
!						ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
!					  enddo	
!					  !vel_ibm2=Vbound(i,j,kpp)/MAX(absU,1e-9)/rhV(i,j,kpp)* ust*(5*log(MAX(distance_to_bed_kp*ust/nu_mol,1.001))-3.05)
!					endif
!					if (yplus<5.) then !viscous sublayer uplus=yplus
!						ust=sqrt(absU*nu_mol/(distance_to_bed_kpp))
!						!vel_ibm2=Vbound(i,j,kpp)/MAX(absU,1e-9)/rhV(i,j,kpp)* ust*ust*distance_to_bed_kp/nu_mol
!					endif
!				  endif
!				  tau_fl_Vnew(i,j) = rhV(i,j,kp)*ust*ust  
!				  DO k=1,kb
!					Vbound(i,j,k)=0.
!				  ENDDO
!				  !Vbound(i,j,kp)=vel_ibm2*rhV(i,j,kp)
!				  absU=sqrt(Vbound(i,j,kp)**2+(0.25*(Ubound(i,j,kp)+Ubound(i,j+1,kp)+Ubound(i-1,j,kp)+Ubound(i-1,j+1,kp)))**2)
!				  absU=absU/rhV(i,j,kp)				  
!				  Vbound(i,j,kp) = Vbound(i,j,kp) / (1. + ust*ust*dt/dz/MAX(absU,1.e-9))  ! implicit = more stable	
!				  IF (dUVdn_IBMbed.eq.0.and.distance_to_bed_kp>0.1*dz) THEN
!					Vbound(i,j,kb) = Vbound(i,j,kp)
!				  ENDIF	
!				ELSEIF(botstress.eq.3) THEN 
!				  rr2=rhV(i,j,kp) 																				!at V-gridpoint
!				  dpdy22=(ppp(i,j+1,kp)-ppp(i,j,kp))/(Rp(i)*(phip(j+1)-phip(j)))/rr2							!at V-gridpoint	
!				  tauu2=0.25*(tau_fl_Uold(i,j)+tau_fl_Uold(i,j+1)+tau_fl_Uold(i-1,j)+tau_fl_Uold(i-1,j+1)) 		!at V-gridpoint
!				  tauv2=tau_fl_Vold(i,j) 																		!at V-gridpoint
!				  uu2=0.25*(Ubound(i,j,kp)+Ubound(i,j+1,kp)+Ubound(i-1,j,kp)+Ubound(i-1,j+1,kp))
!				  vv2=Vbound(i,j,kp)
!				  ust2 = ((tauu2/rr2)**2+(tauv2/rr2)**2)**0.25 								!at V-gridpoint	
!				  call wall_fun_rho_TBL(vv2,uu2,dpdy22,rr2,2.*distance_to_bed_kp,1.,dt,kn,kappa,nu_mol,ust2,tau_fl_Vnew(i,j))
!				  Vbound(i,j,kp) = Vbound(i,j,kp) / (1.+tau_fl_Vnew(i,j)*dt/dz/MAX(ABS(Vbound(i,j,kp)),1.e-9))  
!				  ! don't use vv2 directly as dz and distance_to_bed_kp are not equal
!				  IF (dUVdn_IBMbed.eq.0.and.distance_to_bed_kp>0.1*dz) THEN
!					Vbound(i,j,kb) = Vbound(i,j,kp)
!				  ENDIF				
!				ELSE ! without partial slip bc, simply prescribe 2nd order accurate velocity:
!				  DO k=1,kb
!					Vbound(i,j,k)=0.
!				  ENDDO
!				  Vbound(i,j,kp)=distance_to_bed_kp/distance_to_bed_kpp*Vbound(i,j,kpp) !linear interpolation with zero velocity at bed
!				ENDIF
!			ENDDO
!		ENDDO
!		call bound_cbot(tau_fl_Unew)
!		call bound_cbot(tau_fl_Vnew)
!		tau_fl_Uold = tau_fl_Unew 
!		tau_fl_Vold = tau_fl_Vnew		
!	ENDIF

!c get stuff from other CPU's
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
		    Ubc=rhV(i,0,k)*(0.5*(Ub1(i,k)+Ub1(i+1,k))+Ubc1(i,k))    !0.5*(rho(i,0,k)+rho(i,1,k))*
		    Vbc=rhV(i,0,k)*(Vb1(i,k)+Vbc1(i,k))                     !0.5*(rho(i,0,k)+rho(i,1,k))*
		    Wbc=rhV(i,0,k)*(0.5*(Wb1(i,k)+Wb1(i,k+1))+Wbc1)         !0.5*(rho(i,0,k)+rho(i,1,k))*

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
		    Ubc=rhV(i,jmax,k)*(0.5*(Ub1(i,k)+Ub1(i+1,k))+Ubc1(i,k))                 !0.5*(rho(i,jmax,k)+rho(i,j1,k))
		    Vbc=rhV(i,jmax,k)*(Vb1(i,k)+Vbc1(i,k))                                  !0.5*(rho(i,jmax,k)+rho(i,j1,k))
		    Wbc=rhV(i,jmax,k)*(0.5*(Wb1(i,k)+Wb1(i,k+1))+Wbc1)                      !0.5*(rho(i,jmax,k)+rho(i,j1,k))

!		    Ubc=rho_b*(0.5*(Ub1(i,k)+Ub1(i+1,k))+Ubc1)
!		    Vbc=rho_b*(Vb1(i,k)+Vbc1)
!		    Wbc=rho_b*(0.5*(Wb1(i,k)+Wb1(i,k+1))+Wbc1)
		   Ubound(i,0,k) = Ubf(i,k)
		   Vbound(i,0,k) = Vbf(i,k)
		   Wbound(i,0,k) = Wbf(i,k)

 		   Ubound(i,j1,k) =Ubc*cos_u(j1)+Vbc*sin_u(j1) !Ubb(i,k)
! 		   Ubound(i,j1,k) =2.*(Ubc*cos_u(j1)+Vbc*sin_u(j1))-Ubound(i,jmax,k) !Ubb(i,k)
		   Vbound(i,jmax,k) =-Ubc*sin_v(jmax)+Vbc*cos_v(jmax) !Vbb(i,k)
		   Vbound(i,j1,k) = Vbound(i,jmax,k)
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
			if (botstress.ge.1.and.kn_sidewalls>0.) then  !bed shear stress at lateral walls
				j=1
				rr1=0.5*(rho(i,j,k)+rho(i,j,k+1)) !at W-gridpoint
				rr2=0.5*(rho(i,j,k)+rho(i+1,j,k)) !at U-gridpoint
				ww1=Wbound(i,j,k)
				uu1=0.25*(Ubound(i,j,k)+Ubound(i-1,j,k)+Ubound(i,j,k+1)+Ubound(i-1,j,k+1))	!at W-gridpoint
				ww2=0.25*(Wbound(i,j,k)+Wbound(i,j,k-1)+Wbound(i+1,j,k)+Wbound(i+1,j,k-1)) !at U-gridpoint
				uu2=Ubound(i,j,k)
				call wall_fun_rho(ww1,uu1,rr1,Rp(i)*dphi2(j),1.,dt,kn_sidewalls,kappa,nu_mol,tau)
				call wall_fun_rho(uu2,ww2,rr2,Rp(i)*dphi2(j),1.,dt,kn_sidewalls,kappa,nu_mol,tau)
				Wbound(i,j,k)=ww1
				Ubound(i,j,k)=uu2
				j=0
			endif 		
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
			if (botstress.ge.1.and.kn_sidewalls>0.) then  !bed shear stress at lateral walls
				j=jmax
				rr1=0.5*(rho(i,j,k)+rho(i,j,k+1)) !at W-gridpoint
				rr2=0.5*(rho(i,j,k)+rho(i+1,j,k)) !at U-gridpoint
				ww1=Wbound(i,j,k)
				uu1=0.25*(Ubound(i,j,k)+Ubound(i-1,j,k)+Ubound(i,j,k+1)+Ubound(i-1,j,k+1))	!at W-gridpoint
				ww2=0.25*(Wbound(i,j,k)+Wbound(i,j,k-1)+Wbound(i+1,j,k)+Wbound(i+1,j,k-1)) !at U-gridpoint
				uu2=Ubound(i,j,k)
				call wall_fun_rho(ww1,uu1,rr1,Rp(i)*dphi2(j),1.,dt,kn_sidewalls,kappa,nu_mol,tau)
				call wall_fun_rho(uu2,ww2,rr2,Rp(i)*dphi2(j),1.,dt,kn_sidewalls,kappa,nu_mol,tau)
				Wbound(i,j,k)=ww1
				Ubound(i,j,k)=uu2
			endif 		   
		   Ubound(i,0,k) = Ubf(i,k)
		   Vbound(i,0,k) = Vbf(i,k)
		   Wbound(i,0,k) = Wbf(i,k)
 		   Ubound(i,j1,k) =Ubound(i,jmax,k)
		   Vbound(i,jmax,k) = 0. 
		   Vbound(i,j1,k) = Vbound(i,jmax,k)
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


       !! Get UVW in pipe before wall stresses are taken into account, later these UVW2 are put back in pipe
	Ubound2(0:i1,0:j1,0:k1)=Ubound
	Vbound2(0:i1,0:j1,0:k1)=Vbound
	Wbound2(0:i1,0:j1,0:k1)=Wbound
	rho2(0:i1,0:j1,0:k1)=rho
	! fill Ubound2,Vbound2 with one extra row positive and negative in i,j direction for shear stress 
	if (periodicx.eq.0.or.periodicx.eq.2) then
		Ubound2(-1,0:j1,0:k1)=Ubound(0,0:j1,0:k1)
		Ubound2(i1+1,0:j1,0:k1)=Ubound(imax,0:j1,0:k1)
		Vbound2(-1,0:j1,0:k1)=Vbound(0,0:j1,0:k1)
		Vbound2(i1+1,0:j1,0:k1)=Vbound(imax,0:j1,0:k1)
		Wbound2(-1,0:j1,0:k1)=Wbound(0,0:j1,0:k1)
		Wbound2(i1+1,0:j1,0:k1)=Wbound(imax,0:j1,0:k1)		
		rho2(-1,0:j1,0:k1)=rho(0,0:j1,0:k1)		
		rho2(i1+1,0:j1,0:k1)=rho(imax,0:j1,0:k1)
	else 
		Ubound2(-1,0:j1,0:k1)=Ubound(imax-1,0:j1,0:k1)
		Ubound2(i1+1,0:j1,0:k1)=Ubound(2,0:j1,0:k1)
		Vbound2(-1,0:j1,0:k1)=Vbound(imax-1,0:j1,0:k1)
		Vbound2(i1+1,0:j1,0:k1)=Vbound(2,0:j1,0:k1)
		Wbound2(-1,0:j1,0:k1)=Wbound(imax-1,0:j1,0:k1)
		Wbound2(i1+1,0:j1,0:k1)=Wbound(2,0:j1,0:k1)		
		rho2(-1,0:j1,0:k1)=rho(imax-1,0:j1,0:k1)
		rho2(i1+1,0:j1,0:k1)=rho(2,0:j1,0:k1)		
	endif

!c get stuff from other CPU's
	  call shiftf2(Ubound,ubf)
	  call shiftb2(Ubound,ubb) 
	  call shiftf2(Vbound,vbf)
	  call shiftb2(Vbound,vbb)
	  call shiftf2(Wbound,wbf)
	  call shiftb2(Wbound,wbb) 	  
	  call shiftf2(rho,rbf)
	  call shiftb2(rho,rbb)	  
	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1
		   Ubound2(i,-1,k) = Ubound(i,0,k)
		   Ubound2(i,j1+1,k) =Ubb(i,k)
		   Vbound2(i,-1,k) = Vbound(i,0,k)
		   Vbound2(i,j1+1,k) =Vbb(i,k)
		   Wbound2(i,-1,k) = Wbound(i,0,k)
		   Wbound2(i,j1+1,k) =Wbb(i,k)		   
		   rho2(i,-1,k) = rho(i,0,k)
		   rho2(i,j1+1,k) =rbb(i,k)		   
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1
		   Ubound2(i,-1,k) = Ubf(i,k)
		   Ubound2(i,j1+1,k) =Ubound(i,j1,k)
		   Vbound2(i,-1,k) = Vbf(i,k)
		   Vbound2(i,j1+1,k) =Vbound(i,j1,k)
		   Wbound2(i,-1,k) = Wbf(i,k)
		   Wbound2(i,j1+1,k) =Wbound(i,j1,k)		   
		   rho2(i,-1,k) = rbf(i,k)
		   rho2(i,j1+1,k) =rho(i,j1,k)		   
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1
		   Ubound2(i,-1,k) = Ubf(i,k)
		   Ubound2(i,j1+1,k) =Ubb(i,k)
		   Vbound2(i,-1,k) = Vbf(i,k)
		   Vbound2(i,j1+1,k) =Vbb(i,k)
		   Wbound2(i,-1,k) = Wbf(i,k)
		   Wbound2(i,j1+1,k) =Wbb(i,k)		   
		   rho2(i,-1,k) = rbf(i,k)
		   rho2(i,j1+1,k) =rbb(i,k)		   
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
		   Wbound2(i,-1,k) = Wbf(i,k)
		   Wbound2(i,j1+1,k) =Wbb(i,k)		   
		   rho2(i,-1,k) = rbf(i,k)
		   rho2(i,j1+1,k) =rbb(i,k)		   
		   enddo
		enddo
	endif


	! apply wall stress at monopile wall r=0 before information is exchanged between partitions (hence only j=1-jmax)
	IF (monopile>0.and.mpstress.eq.1) then ! partial slip tau based on log-wall 
         do k=1,kmax 
           do j=1,jmax
			rr1=0.5*(rho(1,j,k)+rho(1,j,k+1)) !at W-gridpoint
			rr2=0.5*(rho(1,j,k)+rho(1,j+1,k)) !at V-gridpoint
			ww1=Wbound(1,j,k)
			vv1=0.25*(Vbound(1,j,k)+Vbound(1,j-1,k)+Vbound(1,j,k+1)+Vbound(1,j-1,k+1))	!at W-gridpoint
			ww2=0.25*(Wbound(1,j,k)+Wbound(1,j,k-1)+Wbound(1,j+1,k)+Wbound(1,j+1,k-1)) !at V-gridpoint
			vv2=Vbound(1,j,k)
			call wall_fun_rho(ww1,vv1,rr1,dr(1),1.,dt,kn_mp,kappa,nu_mol,tau)
			call wall_fun_rho(vv2,ww2,rr2,dr(1),1.,dt,kn_mp,kappa,nu_mol,tau)
			Wbound(1,j,k)=ww1
			Vbound(1,j,k)=vv2
		  enddo
		 enddo 
	ELSEIF (monopile>0.and.mpstress.eq.3) then ! partial slip tau based on simplified TBL wall-law taking pressure gradient into account
         do k=1,kmax 
           do j=1,jmax
			rr1=0.5*(rho(1,j,k)+rho(1,j,k+1)) 											!at W-gridpoint
			rr2=0.5*(rho(1,j,k)+rho(1,j+1,k)) 											!at V-gridpoint 
			dpdz1=(ppp(1,j,k+1)-ppp(1,j,k))/dz/rr1										!at W-gridpoint
			dpdy2=(ppp(1,j+1,k)-ppp(1,j,k))/(Rp(1)*(phip(j+1)-phip(j)))/rr2				!at V-gridpoint
			tauw1=tau2Wold(j,k)
			tauv1=0.25*(tau2Vold(j,k)+tau2Vold(j-1,k)+tau2Vold(j,k+1)+tau2Vold(j-1,k+1))	!at W-gridpoint
			tauw2=0.25*(tau2Wold(j,k)+tau2Wold(j,k-1)+tau2Wold(j+1,k)+tau2Wold(j+1,k-1)) 	!at V-gridpoint
			tauv2=tau2Vold(j,k)
			ust1 = ((tauw1/rr1)**2+(tauv1/rr1)**2)**0.25								!at W-gridpoint
			ust2 = ((tauw2/rr2)**2+(tauv2/rr2)**2)**0.25 								!at V-gridpoint			
			
			ww1=Wbound(1,j,k)
			vv1=0.25*(Vbound(1,j,k)+Vbound(1,j-1,k)+Vbound(1,j,k+1)+Vbound(1,j-1,k+1))	!at W-gridpoint
			ww2=0.25*(Wbound(1,j,k)+Wbound(1,j,k-1)+Wbound(1,j+1,k)+Wbound(1,j+1,k-1)) !at V-gridpoint
			vv2=Vbound(1,j,k)
			call wall_fun_rho_TBL(ww1,vv1,dpdz1,rr1,dr(1),Ru(0)/Rp(1),dt,kn_mp,kappa,nu_mol,ust1,tau2Wnew(j,k))
			call wall_fun_rho_TBL(vv2,ww2,dpdy2,rr2,dr(1),Ru(0)/Rp(1),dt,kn_mp,kappa,nu_mol,ust2,tau2Vnew(j,k))
			Wbound(1,j,k)=ww1
			Vbound(1,j,k)=vv2
		  enddo
		 enddo 
		call bound_2Djk(tau2Wnew)
		call bound_2Djk(tau2Vnew)
		tau2Wold = tau2Wnew 
		tau2Vold = tau2Vnew
		
	ENDIF

	Wbc1=W_b
	if (periodicx.eq.0.and.monopile.eq.-1) then	
		if (Uoutflow.eq.2) then ! convective outflow boundary
      	   do k=1,kmax ! boundaries in i-direction
	   	do j=0,j1
			Wbc2=0.
		    if (bcfile.ne.'') then
			Ubc2(j,k)=Ubcoarse2(j,k)
			Vbc2(j,k)=Vbcoarse2(j,k)
			Wbc2=Wbcoarse2(j,k)
		    endif
		    Ubc=rhU(0,j,k)*(Ub2(j,k)+Ubc2(j,k))                                     !0.5*(rho(0,j,k)+rho(1,j,k))*
		    Vbc=rhU(0,j,k)*(0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k))                    !0.5*(rho(0,j,k)+rho(1,j,k))*
		    Wbc=rhU(0,j,k)*(0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2)	                     !0.5*(rho(0,j,k)+rho(1,j,k))*
		   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
		   Vbound(i1,j,k)   =    Vbound(imax,j,k)
		   Wbound(i1,j,k)   =    Wbound(imax,j,k)
		   Ubc = MAX(Ubc2(j,k)*cos_u(j)+Vbc2(j,k)*sin_u(j),0.)  !Unormal outflow boundary at U-loc; should be >0 otherwise no outflow and divide by zero 
!		   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k)-rhU(imax,j,k)*dr(imax)/(Ubc*dt)*(Unew(imax-1,j,k)-Uold(imax-1,j,k)),0.)
!		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
!		   Vbound(i1,j,k)   =    Vbound(imax,j,k)-rhV(imax,j,k)*(Rp(i1)-Rp(imax))/(Ubc*dt)*(Vnew(imax,j,k)-Vold(imax,j,k))		   
!		   Wbound(i1,j,k)   =    Wbound(imax,j,k)-rhW(imax,j,k)*(Rp(i1)-Rp(imax))/(Ubc*dt)*(Wnew(imax,j,k)-Wold(imax,j,k))		

		   ! solving Ui+1^n+1-Ui+1^n)/dt+Unormal(Ui+1^n+1-Ui^n+1)/dx=0. !implicit treatment of dUdn term
		   Ubound(imax,j,k) = (rhU(imax,j,k)*Unew(imax,j,k)+Ubc*dt/dr(imax)*Ubound(imax-1,j,k))/(1.+Ubc*dt/dr(imax)) 
		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
		   Vbound(i1,j,k) = (rhV(i1,j,k)*Vnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Vbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
		   Wbound(i1,j,k) = (rhW(i1,j,k)*Wnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Wbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
	   
         	enddo   
      	    enddo
		else !Neumann outflow boundary 
      	   do k=1,kmax ! boundaries in i-direction
	   	do j=0,j1
			Wbc2=0.
		    if (bcfile.ne.'') then
			Ubc2(j,k)=Ubcoarse2(j,k)
			Vbc2(j,k)=Vbcoarse2(j,k)
			Wbc2=Wbcoarse2(j,k)
		    endif
		    Ubc=rhU(0,j,k)*(Ub2(j,k)+Ubc2(j,k))                                     !0.5*(rho(0,j,k)+rho(1,j,k))*
		    Vbc=rhU(0,j,k)*(0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k))                    !0.5*(rho(0,j,k)+rho(1,j,k))*
		    Wbc=rhU(0,j,k)*(0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2)	                     !0.5*(rho(0,j,k)+rho(1,j,k))*
		   Ubound(0,j,k)    =    (Ubc*cos_u(j)+Vbc*sin_u(j))
		   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
		   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
		   Vbound(0,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j)) !-Vbound(1,j,k)
		   Vbound(i1,j,k)   =    Vbound(imax,j,k)
		   Wbound(0,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
		   Wbound(i1,j,k)   =    Wbound(imax,j,k)
         	enddo   
      	    enddo		
		endif 
       elseif (periodicx.eq.0.and.monopile>0) then ! closed boundary at i=0 (flow past circular cylinder); mixed inflow/outflow bc at imax; 
		if (monopile.eq.4) then !no slip monopile wall 
         do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    0.
		   Vbound(0,j,k)    =    -Vbound(1,j,k)
		   Wbound(0,j,k)    =    -Wbound(1,j,k)
           enddo   
         enddo	
		 else 
		  do k=1,kmax 
           do j=0,j1
		   Ubound(0,j,k)    =    0.
		   Vbound(0,j,k)    =    Vbound(1,j,k)
		   Wbound(0,j,k)    =    Wbound(1,j,k)
           enddo   
          enddo	
		 endif 
		 if (rank.le.CEILING(REAL(px)/4.)-1) then ! first 1/4 of domain is outflow bc 
			jend = CEILING(REAL(jmax*px)*0.25) !end of first 1/4 global j-index
			jend = MIN(jend-rank*jmax,j1) !end of first 1/4 local j-index 
			if (Uoutflow.eq.2) then ! convective outflow boundary
			 do j=0,jend
			   do k=1,kmax 
				   Ubc = MAX(Ubc2(j,k)*cos_u(j)+Vbc2(j,k)*sin_u(j),0.)  !Unormal outflow boundary at U-loc; should be >0 otherwise no outflow and divide by zero 
!				   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k)-rhU(imax,j,k)*dr(imax)/(Ubc*dt)*(Unew(imax-1,j,k)-Uold(imax-1,j,k)),0.)
!				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
!				   Vbound(i1,j,k)   =    Vbound(imax,j,k)-rhV(imax,j,k)*(Rp(i1)-Rp(imax))/(Ubc*dt)*(Vnew(imax,j,k)-Vold(imax,j,k))	   
!				   Wbound(i1,j,k)   =    Wbound(imax,j,k)-rhW(imax,j,k)*(Rp(i1)-Rp(imax))/(Ubc*dt)*(Wnew(imax,j,k)-Wold(imax,j,k))
				   
				   ! solving Ui+1^n+1-Ui+1^n)/dt+Unormal(Ui+1^n+1-Ui^n+1)/dx=0. !implicit treatment of dUdn term
				   Ubound(imax,j,k) = (rhU(imax,j,k)*Unew(imax,j,k)+Ubc*dt/dr(imax)*Ubound(imax-1,j,k))/(1.+Ubc*dt/dr(imax)) 
				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
				   Vbound(i1,j,k) = (rhV(i1,j,k)*Vnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Vbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
				   Wbound(i1,j,k) = (rhW(i1,j,k)*Wnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Wbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
				   
			   enddo   
			 enddo	
			else !Neumann outflow boundary 
			 do j=0,jend
			   do k=1,kmax 
			     Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
			     Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
			     Vbound(i1,j,k)    =    Vbound(imax,j,k)
			     Wbound(i1,j,k)    =    Wbound(imax,j,k)
			   enddo   
			 enddo				
			endif 
			 do j=jend+1,j1  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
					Wbc2=0.
					if (bcfile.ne.'') then
						Ubc2(j,k)=Ubcoarse2(j,k) 
						Vbc2(j,k)=Vbcoarse2(j,k)
						Wbc2=Wbcoarse2(j,k)
					endif
					Ubc=rhU(imax,j,k)*(Ub2(j,k)+Ubc2(j,k))
					Vbc=rhU(imax,j,k)*(0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k))
					Wbc=rhU(imax,j,k)*(0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2)	
				    Ubound(imax,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Ubound(i1,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Vbound(i1,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
				    Wbound(i1,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
			   enddo   
			 enddo
		 elseif (rank.ge.px-CEILING(REAL(px)/4.)-1) then ! last 1/4 of domain is outflow bc 
			jbeg = jmax*px - CEILING(REAL(jmax*px)*0.25) ! begin of last 1/4 global j-index
			jbeg = MAX(jbeg-rank*jmax,0) !begin of last 1/4 local j-index
			if (Uoutflow.eq.2) then ! convective outflow boundary
			 do j=jbeg,j1 ! outflow 
			   do k=1,kmax 
				   Ubc = MAX(Ubc2(j,k)*cos_u(j)+Vbc2(j,k)*sin_u(j),0.)  !Unormal outflow boundary at U-loc; should be >0 otherwise no outflow and divide by zero 
!				   Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k)-rhU(imax,j,k)*dr(imax)/(Ubc*dt)*(Unew(imax-1,j,k)-Uold(imax-1,j,k)),0.)
!				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
!				   Vbound(i1,j,k)   =    Vbound(imax,j,k)-rhV(imax,j,k)*(Rp(i1)-Rp(imax))/(Ubc*dt)*(Vnew(imax,j,k)-Vold(imax,j,k))	   
!				   Wbound(i1,j,k)   =    Wbound(imax,j,k)-rhW(imax,j,k)*(Rp(i1)-Rp(imax))/(Ubc*dt)*(Wnew(imax,j,k)-Wold(imax,j,k))
				   
				   ! solving Ui+1^n+1-Ui+1^n)/dt+Unormal(Ui+1^n+1-Ui^n+1)/dx=0. !implicit treatment of dUdn term
				   Ubound(imax,j,k) = (rhU(imax,j,k)*Unew(imax,j,k)+Ubc*dt/dr(imax)*Ubound(imax-1,j,k))/(1.+Ubc*dt/dr(imax)) 
				   Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)	
				   Vbound(i1,j,k) = (rhV(i1,j,k)*Vnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Vbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
				   Wbound(i1,j,k) = (rhW(i1,j,k)*Wnew(i1,j,k)+Ubc*dt/(Rp(i1)-Rp(imax))*Wbound(imax,j,k))/(1.+Ubc*dt/(Rp(i1)-Rp(imax))) 
				   
			   enddo   
			 enddo
			else !Neumann outflow boundary 
			 do j=jbeg,j1 ! outflow 
			   do k=1,kmax 
			     Ubound(imax,j,k) =    MAX(Ubound(imax-1,j,k),0.)
			     Ubound(i1,j,k)   =    MAX(Ubound(imax,j,k),0.)
			     Vbound(i1,j,k)    =    Vbound(imax,j,k)
			     Wbound(i1,j,k)    =    Wbound(imax,j,k)
			   enddo   
			 enddo			
			endif 
			 do j=0,MIN(j1,jbeg-1)  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
					Wbc2=0.
					if (bcfile.ne.'') then
						Ubc2(j,k)=Ubcoarse2(j,k) 
						Vbc2(j,k)=Vbcoarse2(j,k)
						Wbc2=Wbcoarse2(j,k)
					endif
					Ubc=rhU(imax,j,k)*(Ub2(j,k)+Ubc2(j,k))
					Vbc=rhU(imax,j,k)*(0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k))
					Wbc=rhU(imax,j,k)*(0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2)	
				    Ubound(imax,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Ubound(i1,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Vbound(i1,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
				    Wbound(i1,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
			   enddo   
			 enddo			 
		 else !middle 2/4 is inflow bc 
			 do j=0,j1  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
					Wbc2=0.
					if (bcfile.ne.'') then
						Ubc2(j,k)=Ubcoarse2(j,k) 
						Vbc2(j,k)=Vbcoarse2(j,k)
						Wbc2=Wbcoarse2(j,k)
					endif
					Ubc=rhU(imax,j,k)*(Ub2(j,k)+Ubc2(j,k))
					Vbc=rhU(imax,j,k)*(0.5*(Vb2(j,k)+Vb2(j+1,k))+Vbc2(j,k))
					Wbc=rhU(imax,j,k)*(0.5*(Wb2(j,k)+Wb2(j,k+1))+Wbc2)	
				    Ubound(imax,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Ubound(i1,j,k)    =    Ubc*cos_u(j)+Vbc*sin_u(j)
				    Vbound(i1,j,k)    =    (-Ubc*sin_v(j)+Vbc*cos_v(j)) !2.*(-Ubc*sin_v(j)+Vbc*cos_v(j))-Vbound(1,j,k)
				    Wbound(i1,j,k)    =    Wbc !2.*Wbc - Wbound(1,j,k) !Wbound(0,j,k)    =  Wbc !- Wbound(1,j,k)
			   enddo   
			 enddo	
		 endif 			
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

	  if (botstress.eq.-2) then
       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         Ubound(i,j,k1)   = -Ubound(i,j,kmax)
         Vbound(i,j,k1)   = -Vbound(i,j,kmax)
         !Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1)   = 0.
         Ubound(i,j,0)    = -Ubound(i,j,1) 
         Vbound(i,j,0)    = -Vbound(i,j,1)
         Wbound(i,j,0)    = Wbed(i,j)*rhW(i,j,0) !0.
         xx=Rp(i)*cos_u(j)-schuif_x
	     yy=Rp(i)*sin_u(j)
	     Wbound(i,j,kmax)=(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo	
	  elseif (botstress.eq.-1) then
       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         Ubound(i,j,k1)   = Ubound(i,j,kmax)
         Vbound(i,j,k1)   = Vbound(i,j,kmax)
         !Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1) = 0.
         Ubound(i,j,0)    = -Ubound(i,j,1) 
         Vbound(i,j,0)    = -Vbound(i,j,1)
         Wbound(i,j,0)    = Wbed(i,j)*rhW(i,j,0) !0.

         xx=Rp(i)*cos_u(j)-schuif_x
	     yy=Rp(i)*sin_u(j)
	     Wbound(i,j,kmax)=(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo
	  elseif (botstress.eq.3.and.IBMorder.ne.2) then 
		do j=1,jmax ! boundaries in k-direction
		  do i=1,imax 
			rr1=rhU(i,j,kbedt(i,j)+1) 					 												!at U-gridpoint
			rr2=rhV(i,j,kbedt(i,j)+1) 																	!at V-gridpoint
			dpdx11=(ppp(i+1,j,kbedt(i,j)+1)-ppp(i,j,kbedt(i,j)+1))/(Rp(i+1)-Rp(i))/rr1					!at U-gridpoint
			dpdy22=(ppp(i,j+1,kbedt(i,j)+1)-ppp(i,j,kbedt(i,j)+1))/(Rp(i)*(phip(j+1)-phip(j)))/rr2		!at V-gridpoint	
			tauu1=tau_fl_Uold(i,j)																		!at U-gridpoint
			tauv1=0.25*(tau_fl_Vold(i,j)+tau_fl_Vold(i+1,j)+tau_fl_Vold(i,j-1)+tau_fl_Vold(i+1,j-1))	!at U-gridpoint
			tauu2=0.25*(tau_fl_Uold(i,j)+tau_fl_Uold(i,j+1)+tau_fl_Uold(i-1,j)+tau_fl_Uold(i-1,j+1)) 	!at V-gridpoint
			tauv2=tau_fl_Vold(i,j) 																		!at V-gridpoint
			
			uu1=Ubound(i,j,kbedt(i,j)+1)-Ubot_TSHD(j)*rr1
			vv1=0.25*(Vbound2(i,j,kbedt(i,j)+1)+Vbound2(i+1,j,kbedt(i,j)+1)+Vbound2(i,j-1,kbedt(i,j)+1)
     &                      +Vbound2(i+1,j-1,kbedt(i,j)+1))-Vbot_TSHD(j)*rr1
			uu2=0.25*(Ubound2(i,j,kbedt(i,j)+1)+Ubound2(i,j+1,kbedt(i,j)+1)+Ubound2(i-1,j,kbedt(i,j)+1)
     &                      +Ubound2(i-1,j+1,kbedt(i,j)+1))-Ubot_TSHD(j)*rr2
			vv2=Vbound(i,j,kbedt(i,j)+1)-Vbot_TSHD(j)*rr2
			ust1 = ((tauu1/rr1)**2+(tauv1/rr1)**2)**0.25								!at U-gridpoint
			ust2 = ((tauu2/rr2)**2+(tauv2/rr2)**2)**0.25 								!at V-gridpoint	
			IF (slip_bot.eq.3) THEN 
			 call wall_fun_rho_TBL(uu1,vv1,dpdx11,rr1,dz,1.,dt,kn,kappa,nu_mol,ust1,tau_fl_Unew(i,j))
			 call wall_fun_rho_TBL(vv2,uu2,dpdy22,rr2,dz,1.,dt,kn,kappa,nu_mol,ust2,tau_fl_Vnew(i,j))
			ELSEIF (slip_bot.eq.4) THEN 
			 call wall_fun_rho_TBL(uu1,vv1,0.,rr1,dz,1.,dt,kn,kappa,nu_mol,ust1,tau_fl_Unew(i,j))
			 call wall_fun_rho_TBL(vv2,uu2,0.,rr2,dz,1.,dt,kn,kappa,nu_mol,ust2,tau_fl_Vnew(i,j))
			ELSEIF (slip_bot.eq.5) THEN 
			 call wall_fun_rho_VD(uu1,vv1,rr1,dz,1.,dt,kn,kappa,nu_mol,ust1,tau_fl_Unew(i,j))
			 call wall_fun_rho_VD(vv2,uu2,rr2,dz,1.,dt,kn,kappa,nu_mol,ust2,tau_fl_Vnew(i,j))
			ENDIF 
			Ubound(i,j,kbedt(i,j)+1)=uu1+Ubot_TSHD(j)*rr1 !Ubound adjusted 1:imax,1:jmax; but later bc 0,i1,0,0,j1 applied
			Vbound(i,j,kbedt(i,j)+1)=vv2+Vbot_TSHD(j)*rr2 !Vbound adjusted 1:imax,1:jmax; but later bc 0,i1,0,0,j1 applied
		  enddo
		 enddo 
		do j=0,j1 ! boundaries in k-direction
		  do i=0,i1		
         Ubound(i,j,k1)   = Ubound(i,j,kmax)
         Vbound(i,j,k1)   = Vbound(i,j,kmax)
         !Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1)   = 0.
         Ubound(i,j,0)    = Ubound(i,j,1) 
         Vbound(i,j,0)    = Vbound(i,j,1)
         Wbound(i,j,0)    = Wbed(i,j)*rhW(i,j,0) !0.

         xx=Rp(i)*cos_u(j)-schuif_x
	     yy=Rp(i)*sin_u(j)
	! apply vertical boundary condition belonging to waves:
	     Wbound(i,j,kmax)=rho_b*(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo	  
	  else
       do j=0,j1 ! boundaries in k-direction
        do i=0,i1
         if ((botstress.eq.1.or.botstress.eq.2).and.(kjet.eq.0.or.LOA>0.).and.IBMorder.ne.2) then  !bed shear stress at bed below
			!! First Ubot_TSHD and Vbot_TSHD is subtracted to determine tau 
			!! only over ambient velocities not over U_TSHD
			!! Then Ubound,Vbound are reduced by tau, 
			!! now Ubot_TSHD and Vbot_TSHD are added again to get correct total velocity
			rr1=0.5*(rho2(i,j,kbedt(i,j)+1)+rho2(i+1,j,kbedt(i,j)+1)) !at U-gridpoint
			rr2=0.5*(rho2(i,j,kbedt(i,j)+1)+rho2(i,j+1,kbedt(i,j)+1)) !at V-gridpoint
			uu1=Ubound(i,j,kbedt(i,j)+1)-Ubot_TSHD(j)*rr1
			vv1=0.25*(Vbound2(i,j,kbedt(i,j)+1)+Vbound2(i+1,j,kbedt(i,j)+1)+Vbound2(i,j-1,kbedt(i,j)+1)
     &                      +Vbound2(i+1,j-1,kbedt(i,j)+1))-Vbot_TSHD(j)*rr1
			uu2=0.25*(Ubound2(i,j,kbedt(i,j)+1)+Ubound2(i,j+1,kbedt(i,j)+1)+Ubound2(i-1,j,kbedt(i,j)+1)
     &                      +Ubound2(i-1,j+1,kbedt(i,j)+1))-Ubot_TSHD(j)*rr2
			vv2=Vbound(i,j,kbedt(i,j)+1)-Vbot_TSHD(j)*rr2
			call wall_fun_rho(uu1,vv1,rr1,dz,1.,dt,kn,kappa,nu_mol,tau)
			tau_fl_Unew(i,j)=tau
			call wall_fun_rho(vv2,uu2,rr2,dz,1.,dt,kn,kappa,nu_mol,tau)
			tau_fl_Vnew(i,j)=tau
			Ubound(i,j,kbedt(i,j)+1)=uu1+Ubot_TSHD(j)*rr1
			Vbound(i,j,kbedt(i,j)+1)=vv2+Vbot_TSHD(j)*rr2
			
	   elseif ((botstress.eq.1.or.botstress.eq.2).and.kjet>0.and.LOA<0.and.IBMorder.ne.2) then ! with flat plate shear stress must be applied:
!	    call wall_fun_rho(Vbound(i,j,kmax-kjet-1),Ubound(i,j,kmax-kjet-1),rho(i,j,kmax-kjet-1),dz,dt,kn,kappa,(depth-bc_obst_h),U_b,nu_mol)
!	    call wall_fun_rho(Ubound(i,j,kmax-kjet-1),Vbound(i,j,kmax-kjet-1),rho(i,j,kmax-kjet-1),dz,dt,kn,kappa,(depth-bc_obst_h),U_b,nu_mol)
			rr1=0.5*(rho2(i,j,kmax-kjet-1)+rho2(i+1,j,kmax-kjet-1)) !at U-gridpoint
			rr2=0.5*(rho2(i,j,kmax-kjet-1)+rho2(i,j+1,kmax-kjet-1)) !at V-gridpoint
		uu2=0.25*(Ubound2(i,j,kmax-kjet-1)+Ubound2(i,j+1,kmax-kjet-1)+Ubound2(i-1,j,kmax-kjet-1)+Ubound2(i-1,j+1,kmax-kjet-1))
		vv2=0.25*(Vbound2(i,j,kmax-kjet-1)+Vbound2(i+1,j,kmax-kjet-1)+Vbound2(i,j-1,kmax-kjet-1)+Vbound2(i+1,j-1,kmax-kjet-1))

	    call wall_fun_rho(Vbound(i,j,kmax-kjet-1),uu2,rr2,dz,1.,dt,kn,kappa,nu_mol,tau)
		tau_fl_Vnew(i,j)=tau
	    call wall_fun_rho(Ubound(i,j,kmax-kjet-1),vv2,rr1,dz,1.,dt,kn,kappa,nu_mol,tau)
		tau_fl_Unew(i,j)=tau

	   endif
         if (botstress.ge.1.and.wallup.eq.2) then  !bed shear stress at top domain
			rr1=0.5*(rho2(i,j,kmax)+rho2(i+1,j,kmax)) !at U-gridpoint
			rr2=0.5*(rho2(i,j,kmax)+rho2(i,j+1,kmax)) !at V-gridpoint
			uu1=Ubound(i,j,kmax)-Ubot_TSHD(j)*rr1
			vv1=0.25*(Vbound2(i,j,kmax)+Vbound2(i+1,j,kmax)+Vbound2(i,j-1,kmax)
     &                      +Vbound2(i+1,j-1,kmax))-Vbot_TSHD(j)*rr1
			uu2=0.25*(Ubound2(i,j,kmax)+Ubound2(i,j+1,kmax)+Ubound2(i-1,j,kmax)
     &                      +Ubound2(i-1,j+1,kmax))-Ubot_TSHD(j)*rr2
			vv2=Vbound(i,j,kmax)-Vbot_TSHD(j)*rr2
			call wall_fun_rho(uu1,vv1,rr1,dz,1.,dt,kn,kappa,nu_mol,tau)
			call wall_fun_rho(vv2,uu2,rr2,dz,1.,dt,kn,kappa,nu_mol,tau)
			Ubound(i,j,kmax)=uu1+Ubot_TSHD(j)*rr1
			Vbound(i,j,kmax)=vv2+Vbot_TSHD(j)*rr2
		endif 
	
         Ubound(i,j,k1)   = Ubound(i,j,kmax)
         Vbound(i,j,k1)   = Vbound(i,j,kmax)
         !Wbound(i,j,kmax) = 0.
         Wbound(i,j,k1)   = 0.
         Ubound(i,j,0)    = Ubound(i,j,1) 
         Vbound(i,j,0)    = Vbound(i,j,1)
         Wbound(i,j,0)    = Wbed(i,j)*rhW(i,j,0) !0.

         xx=Rp(i)*cos_u(j)-schuif_x
	     yy=Rp(i)*sin_u(j)
	! apply vertical boundary condition belonging to waves:
	     Wbound(i,j,kmax)=rho_b*(om_w+kx_w*U_TSHD*cos(phi)+ky_w*U_TSHD*sin(phi))*Hs/2.
     &    *sin(kx_w*(xx-U_TSHD*cos(phi)*tt)+ky_w*(yy-U_TSHD*sin(phi)*tt)-om_w*tt)
        enddo
       enddo
	  endif
		call bound_cbot(tau_fl_Unew)
		call bound_cbot(tau_fl_Vnew)
		tau_fl_Uold = tau_fl_Unew 
		tau_fl_Vold = tau_fl_Vnew		
!	DO n=1,nbedplume !make all forces zero before new forces are applied for bedplume
!		IF (bp(n)%u.ne.-99999.and.bp(n)%velocity_force.ne.0.) THEN ! apply bedplume velocity boundary condition:
!			IF ((bp(n)%forever.eq.1.and.time_np.gt.bp(n)%t0.and.time_np.lt.bp(n)%t_end)
!     &     .or.(bp(n)%forever.eq.0.and.time_n.lt.bp(n)%t0.and.time_np.gt.bp(n)%t0)) THEN
!				Propx_dummy=0.
!				Propy_dummy=0.
!				Propz_dummy=0.
!			ENDIF
!		ENDIF 
!	ENDDO
	 
	   
	DO n=1,nbedplume
	IF (bp(n)%u.ne.-99999.) THEN ! apply bedplume velocity boundary condition:
	IF ((bp(n)%forever.eq.1.and.time_np.gt.bp(n)%t0.and.time_np.lt.bp(n)%t_end)
     &     .or.(bp(n)%forever.eq.0.and.time_n.lt.bp(n)%t0.and.time_np.gt.bp(n)%t0)) THEN
	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif
	 IF (bp(n)%velocity_force.eq.0.) THEN
      do k=MAX(1,CEILING(bp(n)%zbottom/dz)),MIN(kmax,FLOOR(bp(n)%height/dz))! do k=k1,0,-1 !from top to bottom
	   do tel=1,bp(n)%tmax 
	     i=bp(n)%i(tel) 
		 j=bp(n)%j(tel) 
!!       do i=0,i1  
!!         do j=0,j1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
!!			xx=Rp(i)*cos_u(j)-schuif_x
!!			yy=Rp(i)*sin_u(j)
!!!	  		IF (k.le.FLOOR(bp(n)%height/dz).and.k.ge.CEILING(bp(n)%zbottom/dz)) THEN ! obstacle:
!!			xTSHD(1:4)=bp(n)%x*cos(phi)-bp(n)%y*sin(phi)
!!			yTSHD(1:4)=bp(n)%x*sin(phi)+bp(n)%y*cos(phi)
!!			CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
!!!	  		ELSE 
!!!	 			inout=0
!!!	  		ENDIF
!!			if (inout.eq.1) then
			  Ubound(i,j,k)=(bp(n)%u*cos_u(j)+bp(n)%v*sin_u(j))*rhU(i,j,k) !rho(i,j,k)
			  Ubound(MAX(i-1,0),j,k)=(bp(n)%u*cos_u(j)+bp(n)%v*sin_u(j))*rhU(MAX(i-1,0),j,k) !rho(MAX(i-1,0),j,k)
			  Vbound(i,j,k)=(-bp(n)%u*sin_v(j)+bp(n)%v*cos_v(j))*rhV(i,j,k) !rho(i,j,k)
			  Vbound(i,MAX(j-1,0),k)=(-bp(n)%u*sin_v(MAX(j-1,0))+bp(n)%v*cos_v(MAX(j-1,0)))*rhV(i,MAX(j-1,0),k) !rho(i,MAX(j-1,0),k)
			  Wbound(i,j,k)=bp(n)%w*rhW(i,j,k) !rho(i,j,k)
			  Wbound(i,j,MAX(k-1,0))=bp(n)%w*rhW(i,j,MAX(k-1,0)) !rho(i,j,MAX(k-1,0))
!			endif
!		 enddo
	   enddo
	  enddo
	 ELSE
!      do k=MAX(1,CEILING(bp(n)%zbottom/dz)),MIN(kmax,FLOOR(bp(n)%height/dz))! do k=k1,0,-1 !from top to bottom
!       do i=0,i1  
!         do j=0,jmax*px+1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
!			xx=Rp(i)*cos_ut(j)-schuif_x
!			yy=Rp(i)*sin_ut(j)
!!	  		IF (k.le.FLOOR(bp(n)%height/dz).and.k.ge.CEILING(bp(n)%zbottom/dz)) THEN ! obstacle:
!			xTSHD(1:4)=bp(n)%x*cos(phi)-bp(n)%y*sin(phi)
!			yTSHD(1:4)=bp(n)%x*sin(phi)+bp(n)%y*cos(phi)
!			CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
!!	  		ELSE 
!!	 			inout=0
!!	  		ENDIF
!			if (inout.eq.1) then
!			  fbx2=ABS(bp(n)%u)*bp(n)%u/(bp(n)%x(2)-bp(n)%x(1))
!			  fby2=ABS(bp(n)%v)*bp(n)%v/(bp(n)%y(3)-bp(n)%y(2))
!			  fbz2=ABS(bp(n)%w)*bp(n)%w/(bp(n)%height-bp(n)%zbottom)
!			  fbx   =  fbx2 * cos(phi) - fby2 * sin(phi)
!			  fby   =  fbx2 * sin(phi) + fby2 * cos(phi)		
!			  Propx_dummy(i,j,k) = Propx_dummy(i,j,k)+0.5*(cos_ut(j)*fbx+sin_ut(j)*fby)
!			  Propx_dummy(MAX(i-1,0),j,k) = Propx_dummy(MAX(i-1,0),j,k) + 0.5*(cos_ut(j)*fbx+sin_ut(j)*fby)
!			  Propy_dummy(i,j,k) = Propy_dummy(i,j,k)+0.5*(-sin_vt(j)*fbx+cos_vt(j)*fby)
!			  Propy_dummy(i,MAX(j-1,0),k) = Propy_dummy(i,MAX(j-1,0),k) + 0.5*(-sin_vt(MAX(j-1,0))*fbx+cos_vt(MAX(j-1,0))*fby)
!			  Propz_dummy(i,j,k) = Propz_dummy(i,j,k)+0.5*fbz2
!			  Propz_dummy(i,j,MAX(k-1,0)) = Propz_dummy(i,j,MAX(k-1,0)) + 0.5*fbz2
!			endif
!		 enddo
!	   enddo
!	  enddo
	 ENDIF
	ENDIF
	ENDIF
	ENDDO ! bedplume loop
!	DO n=1,nbedplume !make all forces zero before new forces are applied for bedplume
!		IF (bp(n)%u.ne.-99999.and.bp(n)%velocity_force.ne.0.) THEN ! apply bedplume velocity boundary condition:
!			IF ((bp(n)%forever.eq.1.and.time_np.gt.bp(n)%t0.and.time_np.lt.bp(n)%t_end)
!     &     .or.(bp(n)%forever.eq.0.and.time_n.lt.bp(n)%t0.and.time_np.gt.bp(n)%t0)) THEN
!				Ppropx(0:i1,0:j1,1:kmax)=Propx_dummy(0:i1,rank*jmax:rank*jmax+j1,1:kmax)*rhU(0:i1,0:j1,1:kmax)
!				Ppropy(0:i1,0:j1,1:kmax)=Propy_dummy(0:i1,rank*jmax:rank*jmax+j1,1:kmax)*rhV(0:i1,0:j1,1:kmax)
!				Ppropz(0:i1,0:j1,1:kmax)=Propz_dummy(0:i1,rank*jmax:rank*jmax+j1,1:kmax)*rhW(0:i1,0:j1,1:kmax)	
!			ENDIF
!		ENDIF 
!	ENDDO	



	IF (LOA>0.) THEN ! ship:
	  do t=1,tmax_inWpuntTSHD
 	    i=i_inWpuntTSHD(t)
 	    j=j_inWpuntTSHD(t)		
 	    k=k_inWpuntTSHD(t)		
		Wbound(i,j,k)=0.
	    !Wbound(i,j,k)=facIBMw(t)*Wbound(i,j,k) !0.
		!if (facIBMw(t)<1.) then
		!  Wbound(i,j,k)=facIBMw(t)/(1.+facIBMw(t))*Wbound(i,j,k+1)
		!endif
		
	  enddo
	  do t=1,tmax_inUpuntTSHD
 	    i=i_inUpuntTSHD(t)
 	    j=j_inUpuntTSHD(t)		
 	    k=k_inUpuntTSHD(t)		
		Ubound(i,j,k)=0.
	    !Ubound(i,j,k)=facIBMu(t)*Ubound(i,j,k) !0.
		!if (facIBMu(t)<1.) then
		!  Ubound(i,j,k)=facIBMu(t)/(1.+facIBMu(t))*Ubound(i,j,k+1)
		!endif		
	  enddo
	  do t=1,tmax_inVpuntTSHD
 	    i=i_inVpuntTSHD(t)
 	    j=j_inVpuntTSHD(t)		
 	    k=k_inVpuntTSHD(t)		
		Vbound(i,j,k)=0.
	    !Vbound(i,j,k)=facIBMv(t)*Vbound(i,j,k) !0.
		!if (facIBMv(t)<1.) then
		!  Vbound(i,j,k)=facIBMv(t)/(1.+facIBMv(t))*Vbound(i,j,k+1)
		!endif		
	  enddo
	  IF (IBMorder.ne.2) THEN
	  do i=0,i1 ! prescribe UTSHD in obstacle:
	    do j=0,j1
	      do k=1,kbed(i,j) 
	        Ubound(i,j,k)=Ubot_TSHD(j)*rhU(i,j,k) !*rho(i,j,k)
	        Vbound(i,j,k)=Vbot_TSHD(j)*rhV(i,j,k) !*rho(i,j,k)
	      enddo
	    enddo
	  enddo
	  ENDIF 
	  !if (botstress.ge.1) then !if kn_TSHD>0 then tauTSHD independent of botstress
	  if (kn_TSHD.ge.0.0) then !apply tauTSHD independent of botstress
	   do t=1,tmax_inVpunt_tauTSHD
 	    i=i_inVpunt_tauTSHD(t)
 	    j=j_inVpunt_tauTSHD(t)		
 	    k=k_inVpunt_tauTSHD(t)
	    uu2=0.25*(Ubound2(i,j,k)+Ubound2(i,j+1,k)+Ubound2(i-1,j,k)+Ubound2(i-1,j+1,k))	
		rr2=0.5*(rho2(i,j,k)+rho2(i,j+1,k)) !at V-gridpoint		
		
!	    call wall_fun_rho(Vbound(i,j,k),Ubound2(i,j,k),rho(i,j,k),dz,dt,kn_TSHD,kappa,(depth-bc_obst_h),U_b,nu_mol)
	    call wall_fun_rho(Vbound(i,j,k),uu2,rr2,dz,1.,dt,kn_TSHD,kappa,nu_mol,tau)
	   enddo
	  endif
	  !if (botstress.ge.1) then !if kn_TSHD>0 then tauTSHD independent of botstress
          if (kn_TSHD.ge.0.0) then !apply tauTSHD independent of botstress
	   do t=1,tmax_inUpunt_tauTSHD
 	    i=i_inUpunt_tauTSHD(t)
 	    j=j_inUpunt_tauTSHD(t)		
 	    k=k_inUpunt_tauTSHD(t)		
	    vv2=0.25*(Vbound2(i,j,k)+Vbound2(i+1,j,k)+Vbound2(i,j-1,k)+Vbound2(i+1,j-1,k))
		rr1=0.5*(rho2(i,j,k)+rho2(i+1,j,k)) !at U-gridpoint
!	    call wall_fun_rho(Ubound(i,j,k),Vbound2(i,j,k),rho(i,j,k),dz,dt,kn_TSHD,kappa,(depth-bc_obst_h),U_b,nu_mol)
	    call wall_fun_rho(Ubound(i,j,k),vv2,rr1,dz,1.,dt,kn_TSHD,kappa,nu_mol,tau)
	   enddo
	  endif
	  do t=1,tmax_inVpunt_rudder ! apply rudder 
 	    i=i_inVpunt_rudder(t)
 	    j=j_inVpunt_rudder(t)		
 	    k=k_inVpunt_rudder(t)
	    jp=MIN(j+1,j1) ! if j is j1, then jp=j1 --> uu2 is incorrect, but this is no problem as Ubound(i,j1,k) is not used	
	    uu1 = -sin_v(j)*Vbound(i,j,k)+(Ubound(i,j,k))*cos_v(j) ! u wanted along rudder (v wanted = zero)
	    uu2 = -sin_v(j)*Vbound(i,j,k)+(Ubound(i,jp,k))*cos_v(j) ! u wanted along rudder (v wanted = zero)
	    Vbound(i,j,k)= -sin_v(j)*0.5*(uu1+uu2)*rhV(i,j,k)
	    Ubound(i,j,k)= cos_u(j)*uu1*rhU(i,j,k)
	    Ubound(i,jp,k)= cos_u(jp)*uu2*rhU(i,jp,k)
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
       		f=Strouhal*ABS(W_j)/MAX(radius_j*2.,1.e-12) !Strouholt number is 0.3
	else
       		Wjet=interpseries(plumetseries,plumeUseries,plumeseriesloc,tt)
       		f=Strouhal*ABS(W_j)/MAX(radius_j*2.,1.e-12) !Strouholt number is 0.3
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
	  rr=1.-sqrt((xx**2+yy**2)/MAX(radius_j**2,1.e-12))
	  rr=MAX(rr,0.)
	if (outflow_overflow_down.eq.1) then
 		do k=kmax-kjet,kmax-kjet
 		  Wbound(i,j,k)=Wbound2(i,j,k) !Wjet*0.5*(rho(i,j,k)+rho(i,j,k+1))+Awjet/6.*Wjet*0.5*(rho(i,j,k)+rho(i,j,k+1))*fluc
 		enddo
 		do k=kmax-kjet+1,kmax
 		  Wbound(i,j,k)=((jetcorr*rr**(1./W_j_powerlaw)+fluc*Awjet/REAL(azi_n))*Wjet+Wb3in(i,j))*rhW(i,j,k) !*0.5*(rho(i,j,k)+rho(i,j,k+1))
 		enddo
	else
 		do k=kmax-kjet,kmax
 		  Wbound(i,j,k)=Wbound2(i,j,k) !Wjet*0.5*(rho(i,j,k)+rho(i,j,k+1))+Awjet/6.*Wjet*0.5*(rho(i,j,k)+rho(i,j,k+1))*fluc
 		enddo
	 	Wbound(i,j,kmax)=((jetcorr*rr**(1./W_j_powerlaw)+fluc*Awjet/REAL(azi_n))*Wjet+Wb3in(i,j))*rhW(i,j,kmax) !*0.5*(rho(i,j,kmax)+rho(i,j,k1))
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
	Ujetbc2=Ujetbc*cos(azi_angle_u(i,j))-Vjetbc*sin(azi_angle_u(i,j))+Ub3in(i,j)
	Vjetbc2=Ujetbc*sin(azi_angle_u(i,j))+Vjetbc*cos(azi_angle_u(i,j))+Vb3in(i,j)
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
	  	Ubound(i,j,k1)=2.*rhU(i,j,k1)*(Ujetbc2*cos_u(j)+Vjetbc2*sin_u(j))-Ubound(i,j,kmax)   !(rho(i,j,k1)+rho(i+1,j,k1))*
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
	Ujetbc2=Ujetbc*cos(azi_angle_v(i,j))-Vjetbc*sin(azi_angle_v(i,j))+Ub3in(i,j)
	Vjetbc2=Ujetbc*sin(azi_angle_v(i,j))+Vjetbc*cos(azi_angle_v(i,j))+Vb3in(i,j)
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
	  	Vbound(i,j,k1)=(2.*rhV(i,j,k1)*Ujetbc2*sin_v(j)+Vjetbc2*cos_v(j))-Vbound(i,j,kmax) !(Aujet/6.*Wjet*fluc)   (rho(i,j,k1)+rho(i,jp,k1))*(
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
	    Ubound(0,j,k)=(Ujetbc*cos_u(j)+Vjetbc*sin_u(j))*rhU(0,j,k) !*0.5*(rho(0,j,k)+rho(1,j,k))
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
	  Wbound(0,j,k)=2.*rhW(0,j,k)*(cos(atan2(zzz,yy))*Wjetbc+sin(atan2(zzz,yy))*Vjetbc)-Wbound(1,j,k)   !(rho(0,j,k)+rho(0,j,k+1))*
       enddo
       do t=1,tmax_inVpunt2
 	  k=k_inVpunt2(t)
 	  j=j_inVpunt2(t)
	  zzz=k*dz-0.5*dz-zjet2
	  yy=Rp(0)*sin_v(j)
	  rr=1.-sqrt((zzz**2+yy**2)/MAX(radius_j**2,1.e-12))
	  rr=MAX(rr,0.)
 	  fluc=0.
 	  do n=1,azi_n
 	    fluc=fluc+sin(2.*pi*f*tt/n+atan2(zzz,yy))
 	  enddo
	  Ujetbc=Ujet2*rr**(1./W_j_powerlaw)+Aujet2/azi_n2*Ujet2*fluc
	  Vjetbc=Avjet/azi_n*Ujet2*fluc
	  Wjetbc=Awjet/azi_n*Ujet2*fluc
	  jp=min(j+1,j1) !this error is no issue as Vbound(j1) is not used
	  Vbound(0,j,k)=2.*rhV(0,j,k)*
     &                  ((-Ujetbc*sin_v(j)+Vjetbc*cos_v(j))*cos(atan2(zzz,yy))-sin(atan2(zzz,yy))*Wjetbc)-Vbound(1,j,k)        !(rho(0,j,k)+rho(0,jp,k))*
       enddo

       do t=1,tmax_inWpunt_suction ! add suction in front of draghead
 	  k=k_inWpunt_suction(t)
 	  j=j_inWpunt_suction(t)
 	  i=i_inWpunt_suction(t)
	  Wjetbc=pi*MAX(radius_j**2,1.e-12)*W_j*perc_dh_suction/(dr(i)*3.*Dsp) !(m3/s)/m2=m/s per suction cell (perc_dh_suction is corrected for number of dragheads)
	  Wbound(i,j,k)=Wjetbc*rhW(i,j,k) !0.5*(rho(i,j,k)+rho(i,j,k+1))
       enddo

	IF (interaction_bed.ge.4) THEN
	 DO i=1,imax
	  DO j=1,jmax
		zbed(i,j)=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
	  ENDDO 
	 ENDDO 
	 call bound_cbot(zbed)
	ENDIF 
	IF ((interaction_bed.ge.4.or.bedlevelfile.ne.''.or.nobst>0).and.IBMorder.eq.2) THEN ! order-2 IBM j=1,jmax --> at end bound_rhoU lateral boundaries are exchanged to get j=0 and j=j1 right 
		DO i=1,imax
			DO j=1,jmax
			    zb_W=zbed(i,j)
				kb=FLOOR(zb_W/dz)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
				kp=MIN(CEILING(zb_W/dz),kmax)			!location velocity which must be adjusted 2nd order IBM W=(0-1)dz distance from bed
				kpp=MIN(CEILING(zb_W/dz)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity
				distance_to_bed_kp=REAL(kp)*dz-zb_W
				distance_to_bed_kpp=REAL(kpp)*dz-zb_W
				vel_ibm2=distance_to_bed_kp/distance_to_bed_kpp*Wbound(i,j,kpp) !linear interpolation with zero velocity at bed	
				!vel_ibm2=distance_to_bed_kp/distance_to_bed_kpp*Wbound(i,j,kp) !apply 'force' to W(kp) to make zero when bed approaches location kp
				DO k=1,kb
				  Wbound(i,j,k)=0.
				ENDDO					
				Wbound(i,j,kp)=vel_ibm2 
				
			    zb_U=0.5*(zbed(i,j)+zbed(i+1,j))
				kb=FLOOR(zb_U/dz+0.5)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
				kp=MIN(CEILING(zb_U/dz+0.5),kmax)			!location velocity which must be adjusted 2nd order IBM -->(0-1)*dz distance from bed
				kpp=MIN(CEILING(zb_U/dz+0.5)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity -->(1-2)*dz distance from bed
				distance_to_bed_kp=(REAL(kp)-0.5)*dz-zb_U
				distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_U				
				IF (distance_to_bed_kp>0.1*dz) THEN !apply tau on first cell (0.1-1)*dz distance from bed based on velocity of that cell
					kpp=kp
					distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_U
				ELSE !apply tau on second cell (1-1.1)*dz distance from bed based on velocity of that cell, do nothing for cell below: (0-0.1)dz from bed 
					kp=kpp
!				  ELSE !apply tau on first cell (0-0.5)*dz distance from bed based on ustar of second cell (kpp) (1-1.5)*dz distance from bed 
				ENDIF				
				IF (botstress.eq.1.or.botstress.eq.2) THEN	! tests showed it is best to apply tau shear stress at first cell above ibm bed and not prescribe U,V velocity straightaway including influence tau (latter gives too large near bed velocities)			
				  absU=sqrt(Ubound(i,j,kpp)**2+(0.25*(Vbound(i,j,kpp)+Vbound(i,j-1,kpp)+Vbound(i+1,j,kpp)+Vbound(i+1,j-1,kpp)))**2)
				  absU=absU/rhU(i,j,kpp)
				  ust=0.1*absU
				  if (kn>0.) then !walls with rougness (log-law used which can be hydr smooth or hydr rough):
					do tel=1,10 ! 10 iter is more than enough
						z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
						ust=absU/MAX(1./kappa*log(MAX(distance_to_bed_kpp/z0,1.001)),2.) !ust maximal 0.5*absU
					enddo
					!vel_ibm2=Ubound(i,j,kpp)/MAX(absU,1e-9)/rhU(i,j,kpp)* (ust/kappa*log(MAX(distance_to_bed_kp/z0,1.001)))
				  else
					do tel=1,10 ! 10 iter is more than enough
						yplus=MAX(distance_to_bed_kpp*ust/nu_mol,1e-12)
						ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
					enddo
					!vel_ibm2=Ubound(i,j,kpp)/MAX(absU,1e-9)/rhU(i,j,kpp)* ust*(2.5*log(MAX(distance_to_bed_kp*ust/nu_mol,1.001))+5.5)
					if (yplus<30.) then
					  do tel=1,10 ! 10 iter is more than enough
						yplus=MAX(distance_to_bed_kpp*ust/nu_mol,1e-12)
						ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
					  enddo	
					  !vel_ibm2=Ubound(i,j,kpp)/MAX(absU,1e-9)/rhU(i,j,kpp)* ust*(5*log(MAX(distance_to_bed_kp*ust/nu_mol,1.001))-3.05)
					endif
					if (yplus<5.) then !viscous sublayer uplus=yplus
						ust=sqrt(absU*nu_mol/(distance_to_bed_kpp))
						!vel_ibm2=Ubound(i,j,kpp)/MAX(absU,1e-9)/rhU(i,j,kpp)* ust*ust*distance_to_bed_kp/nu_mol
					endif
				  endif
				  tau_fl_Unew(i,j) = rhU(i,j,kp)*ust*ust
				  DO k=1,kb
					Ubound(i,j,k)=0.
				  ENDDO
				  !Ubound(i,j,kp)=vel_ibm2*rhU(i,j,kp)	
				  absU=sqrt(Ubound(i,j,kp)**2+(0.25*(Vbound(i,j,kp)+Vbound(i,j-1,kp)+Vbound(i+1,j,kp)+Vbound(i+1,j-1,kp)))**2)
				  absU=absU/rhU(i,j,kp)				
				  Ubound(i,j,kp) = Ubound(i,j,kp) / (1. + ust*ust*dt/dz/MAX(absU,1.e-9))  ! implicit = more stable	
				  IF (dUVdn_IBMbed.eq.0.and.distance_to_bed_kp>0.1*dz) THEN
					Ubound(i,j,kb) = Ubound(i,j,kp)
				  ENDIF
				ELSEIF(botstress.eq.3) THEN 
				  rr1=rhU(i,j,kpp) 																				!at U-gridpoint
				  dpdx11=(ppp(i+1,j,kpp)-ppp(i,j,kpp))/(Rp(i+1)-Rp(i))/rr1										!at U-gridpoint
				  tauu1=tau_fl_Uold(i,j)																		!at U-gridpoint
				  tauv1=0.25*(tau_fl_Vold(i,j)+tau_fl_Vold(i+1,j)+tau_fl_Vold(i,j-1)+tau_fl_Vold(i+1,j-1))		!at U-gridpoint
				  uu1=Ubound(i,j,kpp)
				  vv1=0.25*(Vbound(i,j,kpp)+Vbound(i+1,j,kpp)+Vbound(i,j-1,kpp)+Vbound(i+1,j-1,kpp))
				  ust1 = ((tauu1/rr1)**2+(tauv1/rr1)**2)**0.25													!at U-gridpoint
				  scal=2.*distance_to_bed_kpp/dz
				  IF (slip_bot.eq.3) THEN
				   call wall_fun_rho_TBL(uu1,vv1,dpdx11,rr1,2.*distance_to_bed_kpp,scal,dt,kn,kappa,nu_mol,ust1,tau_fl_Unew(i,j))
				  ELSEIF (slip_bot.eq.4) THEN 
				   call wall_fun_rho_TBL(uu1,vv1,0.,rr1,2.*distance_to_bed_kpp,scal,dt,kn,kappa,nu_mol,ust1,tau_fl_Unew(i,j))
				  ELSEIF (slip_bot.eq.5) THEN 
				   call wall_fun_rho_VD(uu1,vv1,rr1,2.*distance_to_bed_kpp,scal,dt,kn,kappa,nu_mol,ust1,tau_fl_Unew(i,j))				   
				  ENDIF
				  DO k=1,kb
					Ubound(i,j,k)=0.
				  ENDDO				  
				  Ubound(i,j,kp) = uu1 !Ubound(i,j,kp) / (1.+tau_fl_Unew(i,j)*dt/dz/MAX(ABS(Ubound(i,j,kp)),1.e-9)) 
				  IF (dUVdn_IBMbed.eq.0.and.distance_to_bed_kp>0.1*dz) THEN
					Ubound(i,j,kb) = Ubound(i,j,kp)
				  ENDIF				  
				ELSE ! without partial slip bc, simply prescribe 2nd order accurate velocity:
				  DO k=1,kb
					Ubound(i,j,k)=0.
				  ENDDO
				  Ubound(i,j,kp)=distance_to_bed_kp/distance_to_bed_kpp*Ubound(i,j,kpp) !linear interpolation with zero velocity at bed			
				ENDIF 
				zb_V=0.5*(zbed(i,j)+zbed(i,j+1))
				
				kb=FLOOR(zb_V/dz+0.5)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
				kp=MIN(CEILING(zb_V/dz+0.5),kmax)			!location velocity which must be adjusted 2nd order IBM
				kpp=MIN(CEILING(zb_V/dz+0.5)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity
				distance_to_bed_kp=(REAL(kp)-0.5)*dz-zb_V
				distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_V
				IF (distance_to_bed_kp>0.1*dz) THEN !apply tau on first cell (0.1-1)*dz distance from bed based on velocity of that cell
					kpp=kp
					distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_V
				ELSE !apply tau on second cell (1-1.1)*dz distance from bed based on velocity of that cell, do nothing for cell below: (0-0.1)dz from bed 
					kp=kpp 
!				  ELSE !apply tau on first cell (0-0.1)*dz distance from bed based on ustar of second cell (kpp) (1-1.1)*dz distance from bed 
				ENDIF				
				IF (botstress.eq.1.or.botstress.eq.2) THEN	
				  absU=sqrt(Vbound(i,j,kpp)**2+(0.25*(Ubound(i,j,kpp)+Ubound(i,j+1,kpp)+Ubound(i-1,j,kpp)+Ubound(i-1,j+1,kpp)))**2)
				  absU=absU/rhV(i,j,kpp)
				  ust=0.1*absU
				  if (kn>0.) then !walls with rougness (log-law used which can be hydr smooth or hydr rough):
					do tel=1,10 ! 10 iter is more than enough
						z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
						ust=absU/MAX(1./kappa*log(MAX(distance_to_bed_kpp/z0,1.001)),2.) !ust maximal 0.5*absU
					enddo
					!vel_ibm2=Vbound(i,j,kpp)/MAX(absU,1e-9)/rhV(i,j,kpp)* (ust/kappa*log(MAX(distance_to_bed_kp/z0,1.001)))
				  else
					do tel=1,10 ! 10 iter is more than enough
						yplus=MAX(distance_to_bed_kpp*ust/nu_mol,1e-12)
						ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
					enddo
					!vel_ibm2=Vbound(i,j,kpp)/MAX(absU,1e-9)/rhV(i,j,kpp)* ust*(2.5*log(MAX(distance_to_bed_kp*ust/nu_mol,1.001))+5.5)
					if (yplus<30.) then
					  do tel=1,10 ! 10 iter is more than enough
						yplus=MAX(distance_to_bed_kpp*ust/nu_mol,1e-12)
						ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
					  enddo	
					  !vel_ibm2=Vbound(i,j,kpp)/MAX(absU,1e-9)/rhV(i,j,kpp)* ust*(5*log(MAX(distance_to_bed_kp*ust/nu_mol,1.001))-3.05)
					endif
					if (yplus<5.) then !viscous sublayer uplus=yplus
						ust=sqrt(absU*nu_mol/(distance_to_bed_kpp))
						!vel_ibm2=Vbound(i,j,kpp)/MAX(absU,1e-9)/rhV(i,j,kpp)* ust*ust*distance_to_bed_kp/nu_mol
					endif
				  endif
				  tau_fl_Vnew(i,j) = rhV(i,j,kp)*ust*ust  
				  DO k=1,kb
					Vbound(i,j,k)=0.
				  ENDDO
				  !Vbound(i,j,kp)=vel_ibm2*rhV(i,j,kp)
				  absU=sqrt(Vbound(i,j,kp)**2+(0.25*(Ubound(i,j,kp)+Ubound(i,j+1,kp)+Ubound(i-1,j,kp)+Ubound(i-1,j+1,kp)))**2)
				  absU=absU/rhV(i,j,kp)				  
				  Vbound(i,j,kp) = Vbound(i,j,kp) / (1. + ust*ust*dt/dz/MAX(absU,1.e-9))  ! implicit = more stable	
				  IF (dUVdn_IBMbed.eq.0.and.distance_to_bed_kp>0.1*dz) THEN
					Vbound(i,j,kb) = Vbound(i,j,kp)
				  ENDIF	
				ELSEIF(botstress.eq.3) THEN 
				  rr2=rhV(i,j,kpp) 																				!at V-gridpoint
				  dpdy22=(ppp(i,j+1,kpp)-ppp(i,j,kpp))/(Rp(i)*(phip(j+1)-phip(j)))/rr2							!at V-gridpoint	
				  tauu2=0.25*(tau_fl_Uold(i,j)+tau_fl_Uold(i,j+1)+tau_fl_Uold(i-1,j)+tau_fl_Uold(i-1,j+1)) 		!at V-gridpoint
				  tauv2=tau_fl_Vold(i,j) 																		!at V-gridpoint
				  uu2=0.25*(Ubound(i,j,kpp)+Ubound(i,j+1,kpp)+Ubound(i-1,j,kpp)+Ubound(i-1,j+1,kpp))
				  vv2=Vbound(i,j,kpp)
				  ust2 = ((tauu2/rr2)**2+(tauv2/rr2)**2)**0.25 								!at V-gridpoint	
				  scal=2.*distance_to_bed_kpp/dz
				  IF (slip_bot.eq.3) THEN
				   call wall_fun_rho_TBL(vv2,uu2,dpdy22,rr2,2.*distance_to_bed_kpp,scal,dt,kn,kappa,nu_mol,ust2,tau_fl_Vnew(i,j))
				  ELSEIF (slip_bot.eq.4) THEN 
				   call wall_fun_rho_TBL(vv2,uu2,0.,rr2,2.*distance_to_bed_kpp,scal,dt,kn,kappa,nu_mol,ust2,tau_fl_Vnew(i,j))
				  ELSEIF (slip_bot.eq.5) THEN 
				   call wall_fun_rho_VD(vv2,uu2,rr2,2.*distance_to_bed_kpp,scal,dt,kn,kappa,nu_mol,ust2,tau_fl_Vnew(i,j))
				  ENDIF 
				  DO k=1,kb
					Vbound(i,j,k)=0.
				  ENDDO				  
				  Vbound(i,j,kp) = vv2 !Vbound(i,j,kp) / (1.+tau_fl_Vnew(i,j)*dt/dz/MAX(ABS(Vbound(i,j,kp)),1.e-9))  
				  IF (dUVdn_IBMbed.eq.0.and.distance_to_bed_kp>0.1*dz) THEN
					Vbound(i,j,kb) = Vbound(i,j,kp)
				  ENDIF				
				ELSE ! without partial slip bc, simply prescribe 2nd order accurate velocity:
				  DO k=1,kb
					Vbound(i,j,k)=0.
				  ENDDO
				  Vbound(i,j,kp)=distance_to_bed_kp/distance_to_bed_kpp*Vbound(i,j,kpp) !linear interpolation with zero velocity at bed
				ENDIF
			ENDDO
		ENDDO
		call bound_cbot(tau_fl_Unew)
		call bound_cbot(tau_fl_Vnew)
		tau_fl_Uold = tau_fl_Unew 
		tau_fl_Vold = tau_fl_Vnew		
	ENDIF

	IF (IBMorder.ne.2.or.dUVdn_IBMbed.eq.-2) THEN !19-2-2021 apply for IBMorder.eq.2 to prevent SSC issues near bed when dUVdn_IBMbed=-2
	 IF (interaction_bed.ge.4.or.bedlevelfile.ne.''.or.nobst>0) THEN ! dynamic bed update, IBM bed re-defined every timestep: !order-0 IBM: make vel zero for all cells within bed
		DO i=0,i1
			DO j=0,j1
				DO k=1,kbed(i,j) 
					Ubound(i,j,k)=Ubot_TSHD(j)*rhU(i,j,k) !*rho(i,j,k) !0.
					Vbound(i,j,k)=Vbot_TSHD(j)*rhV(i,j,k) !*rho(i,j,k) !0.
					im=MAX(i-1,0)
					jm=MAX(j-1,0)
					Ubound(im,j,k)=Ubot_TSHD(j)*rhU(im,j,k) !*rho(im,j,k) !0.
					Vbound(i,jm,k)=Vbot_TSHD(jm)*rhV(i,jm,k) !*rho(i,jm,k) !0.					
					Wbound(i,j,k)=0.
				ENDDO 
			ENDDO
		ENDDO
	  !IF (IBMorder.ne.2) THEN
	  IF (dUVdn_IBMbed.eq.0) THEN 
		DO i=1,imax
			DO j=1,jmax
				k=kbed(i,j) 
				IF (kbed(i+1,j).eq.k) Ubound(i,j,k)=Ubound(i,j,MIN(k+1,k1))
				IF (kbed(i,j+1).eq.k) Vbound(i,j,k)=Vbound(i,j,MIN(k+1,k1))
				im=MAX(i-1,0)
				jm=MAX(j-1,0)
				IF (kbed(i-1,j).eq.k) Ubound(im,j,k)=Ubound(im,j,MIN(k+1,k1))
				IF (kbed(i,j-1).eq.k) Vbound(i,jm,k)=Vbound(i,jm,MIN(k+1,k1))				
			ENDDO
		ENDDO	 
	  ENDIF 
	 ENDIF	
	ENDIF
	

	! apply j-boundary conditions for U,V,W again for application tau and ibm2 on j=1:jmax and not 0:j1 + 	
	! to fix small inconsistency in Vbound(:,j1,:) made in jet, jet2 and rudder
	call shiftf(Vbound,vbf) 
	call shiftb(Vbound,vbb) 
	call shiftf(Ubound,ubf) 
	call shiftb(Ubound,ubb) 
	call shiftf(Wbound,wbf) 
	call shiftb(Wbound,wbb) 	
	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then ! boundaries in j-direction

	  elseif (rank.eq.px-1) then
	
	  else
		do k=0,k1
		   do i=1,imax
		   Vbound(i,0,k) = Vbf(i,k)
		   Vbound(i,j1,k) =Vbb(i,k)
		   Ubound(i,0,k) = Ubf(i,k)
		   Ubound(i,j1,k) =Ubb(i,k)	
		   Wbound(i,0,k) = Wbf(i,k)
		   Wbound(i,j1,k) =Wbb(i,k)			   
		   enddo
		enddo
	  endif
	else !periodic in y:
	  do k=0,k1
	   do i=1,imax
	     Vbound(i,0,k) = Vbf(i,k)
	     Vbound(i,j1,k) =Vbb(i,k)
	     Ubound(i,0,k) = Ubf(i,k)
	     Ubound(i,j1,k) =Ubb(i,k)	
	     Wbound(i,0,k) = Wbf(i,k)
	     Wbound(i,j1,k) =Wbb(i,k)			 
	   enddo
	  enddo
	endif

		if (i_periodicx>0) then 
         do k=0,k1 
           do j=0,j1
		   Ubound(0,j,k)    =    Ubound(i_periodicx,j,k)
		   Vbound(0,j,k)    =    Vbound(i_periodicx,j,k)
		   Wbound(0,j,k)    =    Wbound(i_periodicx,j,k)
           enddo   
         enddo	
		endif 

      end

      subroutine make_UtransC_zeroIBM(Ubound,Vbound,Wbound)

      USE nlist

      implicit none

	integer t,im,jm,kb,kp,kpp
	real Ubound(0:i1,0:j1,0:k1),Vbound(0:i1,0:j1,0:k1),Wbound(0:i1,0:j1,0:k1)
	real Ubound2(0:i1,0:j1,0:k1),Vbound2(0:i1,0:j1,0:k1),Wbound2(0:i1,0:j1,0:k1)
	real zb_W,zb_U,zb_V,vel_ibm2,distance_to_bed_kp,distance_to_bed_kpp,yplus,absU,z0,ust


	Ubound2=Ubound
	Vbound2=Vbound
	Wbound2=Wbound
	  do t=1,tmax_inWpuntTSHD
 	    i=i_inWpuntTSHD(t)
 	    j=j_inWpuntTSHD(t)		
 	    k=k_inWpuntTSHD(t)		
	    Wbound(i,j,k)=0.
		!Wbound(i,j,k)=facIBMw(t)*Wbound(i,j,k) !0.
		!if (facIBMw(t)<1.) then
		!  Wbound(i,j,k)=facIBMw(t)/(1.+facIBMw(t))*Wbound(i,j,k+1)
		!endif			
	  enddo
	  do t=1,tmax_inUpuntTSHD
 	    i=i_inUpuntTSHD(t)
 	    j=j_inUpuntTSHD(t)		
 	    k=k_inUpuntTSHD(t)		
	    Ubound(i,j,k)=0.
	    !Ubound(i,j,k)=facIBMu(t)*Ubound(i,j,k) !0.
		!if (facIBMu(t)<1.) then
		!  Ubound(i,j,k)=facIBMu(t)/(1.+facIBMu(t))*Ubound(i,j,k+1)
		!endif			
	  enddo
	  do t=1,tmax_inVpuntTSHD
 	    i=i_inVpuntTSHD(t)
 	    j=j_inVpuntTSHD(t)		
 	    k=k_inVpuntTSHD(t)		
	    Vbound(i,j,k)=0.
	    !Vbound(i,j,k)=facIBMv(t)*Vbound(i,j,k) !0.
		!if (facIBMv(t)<1.) then
		!  Vbound(i,j,k)=facIBMv(t)/(1.+facIBMv(t))*Vbound(i,j,k+1)
		!endif			
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

	IF (interaction_bed.ge.4.or.bedlevelfile.ne.''.or.nobst>0) THEN ! dynamic bed update, IBM bed re-defined every timestep:
		DO i=0,i1
			DO j=0,j1
				DO k=1,kbed(i,j) 
					Ubound(i,j,k)=0. !Ubot_TSHD(j) ! 0.
					Vbound(i,j,k)=0. !Vbot_TSHD(j) ! 0.
					im=MAX(i-1,0)
					jm=MAX(j-1,0)
					Ubound(im,j,k)=0. !Ubot_TSHD(j) !0.
					Vbound(i,jm,k)=0. !Vbot_TSHD(jm) !0.					
					Wbound(i,j,k)=0.
				ENDDO
			ENDDO
		ENDDO
	ENDIF	 
!	IF ((interaction_bed.ge.4.or.bedlevelfile.ne.''.or.nobst>0).and.IBMorder.eq.2) THEN ! order-2 IBM before information is exchanged between partitions (hence only j=1-jmax); order-0 IBM is done at and of this subroutine 
!		DO i=1,imax
!			DO j=1,jmax
!				IF (interaction_bed.ge.4) THEN
!				  !zb_W=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(cnewbot(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
!				  zb_W=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
!				ELSE 
!				  zb_W=zbed(i,j)
!				ENDIF
!				kb=FLOOR(zb_W/dz)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
!				kp=MIN(CEILING(zb_W/dz),kmax)			!location velocity which must be adjusted 2nd order IBM
!				kpp=MIN(CEILING(zb_W/dz)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity
!				distance_to_bed_kp=REAL(kp)*dz-zb_W
!				distance_to_bed_kpp=REAL(kpp)*dz-zb_W
!				vel_ibm2=distance_to_bed_kp/distance_to_bed_kpp*Wbound(i,j,kpp) !linear interpolation with zero velocity at bed		
!				  DO k=kbed(i,j),kb
!					Wbound(i,j,k)=0.
!				  ENDDO						
!				Wbound(i,j,kp)=vel_ibm2 				
!
!				IF (interaction_bed.ge.4) THEN
!!					  zb_U=0.5*(zb_W+REAL(MAX(kbed(i+1,j)-1,0))*dz+(SUM(cnewbot(1:nfrac,i+1,j))+
!!     &				SUM(Clivebed(1:nfrac,i+1,j,kbed(i+1,j))))/cfixedbed*dz)
!				  zb_U=0.5*(zb_W+REAL(MAX(kbed(i+1,j)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i+1,j))+
!     &			  SUM(Clivebed(1:nfrac,i+1,j,kbed(i+1,j))))/cfixedbed*dz)	 
!				ELSE 
!				  zb_U=0.5*(zbed(i,j)+zbed(i+1,j))
!				ENDIF
!				kb=FLOOR(zb_U/dz+0.5)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
!				  DO k=kbed(i,j),kb
!					Ubound(i,j,k)=0.
!				  ENDDO
!				IF (slip_bot.ge.1) THEN	! tests showed it is best to apply tau shear stress at first cell above ibm bed and not prescribe U,V velocity straightaway including influence tau (latter gives too large near bed velocities)							  
!					! tau is only applied in bound_rhoU				
!				ELSE ! without partial slip bc, simply prescribe 2nd order accurate velocity:		
!					kp=MIN(CEILING(zb_U/dz+0.5),kmax)			!location velocity which must be adjusted 2nd order IBM
!					kpp=MIN(CEILING(zb_U/dz+0.5)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity
!					distance_to_bed_kp=(REAL(kp)-0.5)*dz-zb_U
!					distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_U	
!					Ubound(i,j,kp)=distance_to_bed_kp/distance_to_bed_kpp*Ubound(i,j,kpp) !linear interpolation with zero velocity at bed			
!				ENDIF 
!
!					IF (interaction_bed.ge.4) THEN
!!				  zb_V=0.5*(zb_W+REAL(MAX(kbed(i,j+1)-1,0))*dz+(SUM(cnewbot(1:nfrac,i,j+1))+
!!     &				SUM(Clivebed(1:nfrac,i,j+1,kbed(i,j+1))))/cfixedbed*dz)
!				  zb_V=0.5*(zb_W+REAL(MAX(kbed(i,j+1)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i,j+1))+
!     &				SUM(Clivebed(1:nfrac,i,j+1,kbed(i,j+1))))/cfixedbed*dz)	 
!					ELSE 
!					  zb_V=0.5*(zbed(i,j)+zbed(i,j+1))
!					ENDIF	
!					kb=FLOOR(zb_V/dz+0.5)						!location vel=0 for 2nd order IBM because below 2nd-order zbed
!				  DO k=kbed(i,j),kb
!					Vbound(i,j,k)=0.
!				  ENDDO
!				IF (slip_bot.ge.1) THEN				
!					! tau is only applied in bound_rhoU
!				ELSE ! without partial slip bc, simply prescribe 2nd order accurate velocity:				  
!					kp=MIN(CEILING(zb_V/dz+0.5),kmax)			!location velocity which must be adjusted 2nd order IBM
!					kpp=MIN(CEILING(zb_V/dz+0.5)+1,kmax+1)		!location one cell above to determine 2nd order ibm velocity
!					distance_to_bed_kp=(REAL(kp)-0.5)*dz-zb_V
!					distance_to_bed_kpp=(REAL(kpp)-0.5)*dz-zb_V						
!					Vbound(i,j,kp)=distance_to_bed_kp/distance_to_bed_kpp*Vbound(i,j,kpp) !linear interpolation with zero velocity at bed
!				ENDIF
!			ENDDO
!		ENDDO
!	ENDIF	
	
	end

	

      subroutine bound_c(Cbound,cjet,n,ddtt)
      
      USE nlist

      implicit none
c
!       include 'param.txt'
!       include 'common.txt'

c
      integer jtmp,t,kcheck,n,n2,inout,tel,jbeg,jend
c
      real Cbound(0:i1,0:j1,0:k1),cjet,Cbound2(0:i1,0:j1,0:k1)
      real ubb(0:i1,0:k1),val,theta,rbc,xx,yy,r_orifice2,rjet,theta_U,theta_V
      real cbf(0:i1,0:k1)
      real cbb(0:i1,0:k1)
	real xTSHD(4),yTSHD(4),phi,ddtt,interpseries
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

	! start with placing bedplume concentrations in order to get lateral boundaries j=0 and j=jmax*px+1 correct (no error in diffusion out of domain)

	DO n2=1,nbedplume
	IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end)
     &     .or.(bp(n2)%forever.eq.0.and.time_n.lt.bp(n2)%t0.and.time_np.gt.bp(n2)%t0)) THEN
	! rotation ship for ambient side current
!!	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
!!	  phi=atan2(V_b,1.e-12)
!!	else
!!	  phi=atan2(V_b,(U_TSHD-U_b))
!!	endif
      do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz))! do k=0,k1
	   do tel=1,bp(n2)%tmax 
	     i=bp(n2)%i(tel) 
		 j=bp(n2)%j(tel) 	  
!!       do i=0,i1  
!!         do j=j1,0,-1       
!!	  xx=Rp(i)*cos_u(j)-schuif_x
!!	  yy=Rp(i)*sin_u(j)
!!!	  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle:
!!		xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
!!		yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!		CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
!!!	  ELSE 
!!!	 	inout=0
!!!	  ENDIF
!!	  if (inout.eq.1.and.k>kbed(i,j)) then
      if (k>kbed(i,j)) then
	    if (bp(n2)%c(n)>0.) then
		  Cbound(i,j,k)=bp(n2)%c(n)
		else
		  Cbound(i,j,k)=(Cbound(i,j,k)+bp(n2)%sedflux(n)*ddtt*fc_global(i,j+jmax*rank,k)/bp(n2)%volncells/frac(n)%rho)
     &  /(1.-MIN(0.,bp(n2)%Q*bp(n2)%changesedsuction*fc_global(i,j+jmax*rank,k))*ddtt/bp(n2)%volncells) !when Q negative, remove sediment from cell as well   IMPLICIT 
	 ! IMPLICIT: c^n+1-c^n=-Qout_cel/Vol_cel*dt*c^n+1 --> c^n+1 = c^n/(1+Qout_cel/Vol_cel*dt)
		endif
		! rho is calculated in state called after fkdat
	   endif
!	  enddo
	 enddo
	enddo
	  ! remove sediment from obstacles/TSHD after placement of bedplume:
!	  do t=1,tmax_inPpuntTSHD ! when no TSHD then this loop is skipped
!	    k=k_inPpuntTSHD(t)		
!	    i=i_inPpuntTSHD(t)
!            j=j_inPpuntTSHD(t)
!            Cbound(i,j,k)=0.  ! remove sediment from hull
!	  enddo
	ENDIF
	ENDDO ! bedplume loop
	

!c get stuff from other CPU's

	
	call shiftf(Cbound,cbf) 
	call shiftb(Cbound,cbb) 

!	write(*,*),'rank,px,INT(px/2)',rank,px,INT(px/2)

	! Cbcoarse1 and Cbcoarse2 are zero when no bcfile is used, or when nfrac is zero in bcfile;
	! in that case simply zero lateral and inflow bc are applied to C

	if (bcfile.eq.'') then
		if (rank.eq.0) then
			Cbcoarse1(n,1:imax,1:kmax)=0. !Cbound(1:imax,1,1:kmax)
		elseif (rank.eq.px-1) then	
			Cbcoarse1(n,1:imax,1:kmax)=0. !Cbound(1:imax,jmax,1:kmax)
		endif
		Cbcoarse2(n,0:j1,1:kmax)=0. !Cbound(1,0:j1,1:kmax)
	endif
		
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
	elseif (periodicy.eq.2) then ! free slip in y: no advective flux through lateral walls; hence dcdn prevents loss of sediment via diff over walls
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

	if (periodicx.eq.0.and.monopile.eq.-1) then
	      do k=1,kmax ! boundaries in i-direction
		 do j=0,j1
			   Cbound(0,j,k)    =    Cbcoarse2(n,j,k) !Cbound(1,j,k)
			   Cbound(i1,j,k)   =    Cbound(imax,j,k)
		 enddo   
	      enddo
       elseif (periodicx.eq.0.and.monopile>0) then ! closed boundary at i=0 (flow past circular cylinder); mixed inflow/outflow bc at imax; 3 = partial slip with tau-wall/4=free slip
         do k=1,kmax 
           do j=0,j1
		   Cbound(0,j,k)    =    Cbound(1,j,k) !dcdn=0. to have no loss of sediment at wall monopile 
           enddo   
         enddo	
		 if (rank.le.CEILING(REAL(px)/4.)-1) then ! first 1/4 of domain is outflow bc 
			jend = CEILING(REAL(jmax*px)*0.25) !end of first 1/4 global j-index
			jend = MIN(jend-rank*jmax,j1) !end of first 1/4 local j-index 
			 do j=0,jend
			   do k=1,kmax 
			     Cbound(i1,j,k)   =    Cbound(imax,j,k)
			   enddo   
			 enddo			
			 do j=jend+1,j1  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
				    Cbound(i1,j,k)    =    Cbcoarse2(n,j,k)
			   enddo   
			 enddo
		 elseif (rank.ge.px-CEILING(REAL(px)/4.)-1) then ! last 1/4 of domain is outflow bc 
			jbeg = jmax*px - CEILING(REAL(jmax*px)*0.25) ! begin of last 1/4 global j-index
			jbeg = MAX(jbeg-rank*jmax,0) !begin of last 1/4 local j-index
			 do j=jbeg,j1 ! outflow 
			   do k=1,kmax 
			     Cbound(i1,j,k)    =    Cbound(imax,j,k)
			   enddo   
			 enddo	
			 do j=0,MIN(j1,jbeg-1)  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
				    Cbound(i1,j,k)    =    Cbcoarse2(n,j,k)
			   enddo   
			 enddo			 
		 else !middle 2/4 is inflow bc 
			 do j=0,j1  !remainder belongs to second 1/4 of domain is inflow 
			   do k=1,kmax 
				    Cbound(i1,j,k)    =    Cbcoarse2(n,j,k)
			   enddo   
			 enddo	
		 endif 			  
       elseif (periodicx.eq.2) then ! no outflow in x direction:
	      do k=1,kmax ! boundaries in i-direction
		 do j=0,j1
			   Cbound(0,j,k)    =    Cbound(1,j,k) !0. 
			   Cbound(i1,j,k)   =    Cbound(imax,j,k) !0. 
		 enddo   
	      enddo
	else ! periodic x boundaries
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

!			IF (interaction_bed.ge.4) THEN ! dynamic bed update, IBM bed re-defined every timestep, make c inside bed zero (diffcof is made zero in turbulence.f to eliminate incorrect diffusion into bed): ! 3-10-2017:why needed, hence switched off; might give trouble for some simulations (hopper filling?) 
!			DO i=0,i1 
!				DO j=0,j1
!					DO k=0,kbed(i,j)
!						Cbound(i,j,k)=0. 
!					ENDDO
!				ENDDO
!			ENDDO
!		ENDIF	
		!! from 17-2-2017 bedlevelfile is dealt with via kbed, not via TSHD immersed boundary; with utr,vtr,wtr and diffcof zero for all sides of bed-cell concentration in bed should stay exactly zero and tricks like below to remove sediment from bed and add it to lowest fluid cell shouldn't be necessary
		if (.false.) then !25-1-2018 switched off because not needed anymore and potentially dangerous
	kcheck=kmax-(kjet-1) 
	  do t=1,tmax_inPpuntTSHD ! when no TSHD then this loop is skipped
 	    k=k_inPpuntTSHD(t)		
	    i=i_inPpuntTSHD(t)
            j=j_inPpuntTSHD(t)
	      if (k.eq.kcheck.and.frac(n)%ws.lt.0.and.kjet>0.) then !and.interaction_bed.ne.4) then 
	        Cbound(i,j,k-1)=Cbound(i,j,k-1)+Cbound(i,j,k) ! move air from first (lowest) cell in TSHD-hull to first cell in fluid
		!! Idea is to undo the effect of diffusion and air bubble rising up into hull
	      endif
		  !! switched off 5-10-2017:
!	      if (k.gt.0.and.k.eq.kbed(i,j)) then !.and.interaction_bed.ne.4) then ! does not happen for bed with interaction_bed=4, or for bedlevelfile 
!	        Cbound(i,j,kbed(i,j)+1)=Cbound(i,j,kbed(i,j)+1)+Cbound(i,j,k) ! move sediment from heightest cell in obstacle to first cell in fluid
!		!! Idea is to undo the effect of diffusion into the obstacle
!	      endif
		    !IF (interaction_bed.ne.4) then
            Cbound(i,j,k)=0.  ! remove sediment from hull, not only in lowest line of cells in hull
			!endif
	  enddo
	  endif
	  
	  
      do t=1,tmax_inPpunt
	i=i_inPpunt(t)
	j=j_inPpunt(t)
	do k=k1-kjet-2,k1
	  Cbound(i,j,k)=Cbound2(i,j,k)
	enddo
      enddo

		if (i_periodicx>0) then 
         do k=0,k1 
           do j=0,j1
		   Cbound(0,j,k)    =    Cbound(i_periodicx,j,k)
           enddo   
         enddo	
		endif 
		
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
!c get stuff from other CPU's
	  
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
		   Cbound(0,j)    =    Cbound(1,j)
		   Cbound(i1,j)   =    Cbound(imax,j)
         enddo   
	elseif (periodicx.eq.2) then
         do j=0,j1
		   Cbound(0,j)    =    Cbound(1,j)
		   Cbound(i1,j)   =    Cbound(imax,j)
         enddo   	
	else 
         do j=0,j1
		   Cbound(0,j)    =    Cbound(imax,j)
		   Cbound(i1,j)   =    Cbound(1,j)
         enddo   
	endif

		end

	subroutine bound_cbot_integer(Cbound) 
      
      USE nlist

      implicit none
c
      integer jtmp,t,kcheck
c
      integer*8 Cbound(0:i1,0:j1)
      integer*8 cbf(0:i1)
      integer*8 cbb(0:i1)
c
c
c*************************************************************
c
c     Subroutine sets the boundary conditions for Cbound for integer array
c
c
c*************************************************************
!c get stuff from other CPU's
	  
	call shiftf_intl(Cbound,cbf) 
	call shiftb_intl(Cbound,cbb) 

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
	if (periodicx.eq.0.or.periodicx.eq.2) then
         do j=0,j1
		   Cbound(0,j)    =    Cbound(1,j)
		   Cbound(i1,j)   =    Cbound(imax,j)
         enddo   
	else 
         do j=0,j1
		   Cbound(0,j)    =    Cbound(imax,j)
		   Cbound(i1,j)   =    Cbound(1,j)
         enddo   
	endif

		end		
		


	subroutine bound_2Djk(Cbound)
      
      USE nlist

      implicit none
c
!       include 'param.txt'
!       include 'common.txt'

c
      integer jtmp,t,kcheck
c
      real Cbound(0:j1,0:k1)
      real cbf(0:k1)
      real cbb(0:k1)

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
!c get stuff from other CPU's
	  
	call shiftf_k(Cbound,cbf) 
	call shiftb_k(Cbound,cbb) 

	  if (periodicy.eq.0.or.periodicy.eq.2) then
		if (rank.eq.0) then ! boundaries in j-direction
			do k=1,kmax
			   Cbound(0,k) = Cbound(1,k) !cbf(i,k)
			   Cbound(j1,k) =cbb(k) 
			enddo
		elseif (rank.eq.px-1) then
			do k=1,kmax
			   Cbound(0,k) = cbf(k)
			   Cbound(j1,k) = Cbound(jmax,k) !cbb(i,k) 
			enddo	
		else
			do k=1,kmax
			   Cbound(0,k) = cbf(k)
			   Cbound(j1,k) =cbb(k) 
			enddo
		endif
	  elseif (periodicy.eq.1) then ! periodic in y:
		do k=1,kmax
			   Cbound(0,k) = cbf(k)
			   Cbound(j1,k) =cbb(k) 
		enddo
	  endif


      do j=0,j1 ! boundaries in k-direction
         Cbound(j,k1)   = Cbound(j,kmax)
         Cbound(j,0)    = Cbound(j,1)
       enddo

		end
		

      subroutine bound_p(Cbound)
      
      USE nlist

      implicit none

      real Cbound(0:i1,0:j1,0:k1)
      real cbf(0:i1,0:k1)
      real cbb(0:i1,0:k1)
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

	!c get stuff from other CPU's
	  call shiftf(Cbound,cbf) 
	  call shiftb(Cbound,cbb) 

	  if (periodicy.eq.0.or.periodicy.eq.2) then
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
			   Cbound(i,j1,k) = Cbound(i,jmax,k) !cbb(i,k) 
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
	  endif

	  if (periodicx.eq.0.or.periodicx.eq.2) then
	      do k=1,kmax ! boundaries in i-direction
		 do j=0,j1
			   Cbound(0,j,k)    =    Cbound(1,j,k)
			   Cbound(i1,j,k)   =    Cbound(imax,j,k)
		 enddo   
	      enddo
	  else ! periodic x boundaries
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
		
!	 DO i=0,i1
!	  DO j=0,j1
!		Cbound(i,j,kbed(i,j))=Cbound(i,j,MIN(kbed(i,j)+1,k1)) ! make dpdn=0 at bed
!	  ENDDO
!	 ENDDO
	 
      end


		
	subroutine wall_fun(uu,vv,rr,dz,dt,kn,kappa,nu_mol)
		
	implicit none

c*************************************************************
c	Wall function
c	Determine tau-wall with sqrt(uu**2,vv**2) and adapt uu	
c	Subtract tau_wall/rho*dt*dx*dy/(dx*dy*dz) from uu (m/s)
c*************************************************************
	real uu,vv,rr,dz,dt,absU,ust,z0,kn,kappa,yplus,tau,nu_mol
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
			yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
			ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
		enddo
	   	if (yplus<30.) then
		  do tel=1,10 ! 10 iter is more than enough
			yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
			ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
		  enddo	
		endif
		if (yplus<5.) then !viscous sublayer uplus=yplus
			ust=sqrt(absU*nu_mol/(0.5*dz))
		endif

	endif
! 	write(*,*) 'ust wall function',ust

	tau=ust*ust  ! omit rho otherwise first *rho then /rho to get dimensions right
	!uu = uu - tau*dt/dz*uu/MAX(absU,1.e-6) 	!! only uu is adjusted
	uu = uu / (1. + tau*dt/dz/MAX(absU,1.e-12)) 	!! only uu is adjusted !! implicit = more stable
	end

	subroutine wall_fun_rho(uu,vv,rr,dz,scal,dt,kn,kappa,nu_mol,tau)
		
	implicit none

c*************************************************************
c	Wall function
c	Determine tau-wall with sqrt(uu**2,vv**2) and adapt uu	
c	Substract tau_wall*dt*dx*dy/(dx*dy*dz) from uu (kg/m3*m/s)
c*************************************************************
	real uu,vv,rr,dz,dt,absU,ust,z0,kn,kappa,yplus,tau,nu_mol,scal
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
			yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
			ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
		enddo
	   	if (yplus<30.) then
		  do tel=1,10 ! 10 iter is more than enough
			yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
			ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
		  enddo	
		endif
		if (yplus<5.) then !viscous sublayer uplus=yplus
			ust=sqrt(absU*nu_mol/(0.5*dz))
		endif
	endif

	tau=rr*ust*ust
	!uu = uu - tau*dt/dz*uu/rr/MAX(absU,1.e-6) 	!! only uu is adjusted
	uu = uu / (1. + tau*scal*dt/dz/rr/MAX(absU,1.e-9)) 	!! only uu is adjusted ! implicit = more stable
	end

	subroutine wall_fun_rho_TBL(uu,vv,Fi,rr,dz,scal,dt,kn,kappa,nu_mol,ust,tau)
		
	implicit none

c*************************************************************
c	Wall function based on simplified TBL with influence of dpdx 
c   Based on (Wang and Moin 2002) Dynamic wall modeling for large-eddy simulation of complex turbulent flows
c	Determine tau-wall and adapt uu	
c	Substract tau_wall*dt*dx*dy/(dx*dy*dz) from uu (kg/m3*m/s)
c*************************************************************
	real uu,vv,Fi,rr,dz,dt,ust,z0,kn,kappa,yplus,tau,nu_mol,absU
	real int1,int2,zzz,AA,nu_tot,zplus,ddzz(1:50),zz(1:50),scal,fac,kplus
	integer tel,tel0,nlayer
	
	nlayer = 50
	AA = 19. ! Van Driest damping factor (Wang and Moin 2002) and VanBalen use value 19, others often use 26. In my matlab tests 26 gives better velocity profiles, but 19 gives better (closer to log law) tau based on velocity at certain distance.
	ddzz(1) = dz/100. !0.5*dz/50. -->input velocity uu,vv defined at 0.5*dz from wall
	do tel =2,nlayer ! apply 50 grid layers between wall z=0 and first velocity point at 0.5*dz with grow factor 
		ddzz(tel)=ddzz(tel-1)*1.1 
	enddo
	ddzz=0.5*dz/SUM(ddzz)*ddzz 
	zz(1) = 0.5*ddzz(1)
	do tel =2,nlayer 
		zz(tel)=zz(tel-1)+0.5*ddzz(tel-1)+0.5*ddzz(tel) 
	enddo	
	int1=0.
	int2=0. 
!	absU=sqrt((uu/rr)**2+(vv/rr)**2)
!	if (ust<0.0001*absU) then 		
!		ust=0.1*absU
!		if (kn>0.) then !walls with rougness (log-law used which can be hydr smooth or hydr rough):
!			do tel=1,10 ! 10 iter is more than enough
!				z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
!				ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
!			enddo
!		else
!			do tel=1,10 ! 10 iter is more than enough
!				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
!				ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
!			enddo
!			if (yplus<30.) then
!			  do tel=1,10 ! 10 iter is more than enough
!				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
!				ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
!			  enddo	
!			endif
!			if (yplus<5.) then !viscous sublayer uplus=yplus
!				ust=sqrt(absU*nu_mol/(0.5*dz))
!			endif
!		endif	  
!	endif 
	if (kn>0.) then
		z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
	else
		z0=0.
	endif
	! I found a nice relation between 1/kappa*log(y/y0) for rough walls and VanDriest velocity profile with z0 added to z in zplus in combination with Crimaldi 2006 correction	
	kplus = kn*ust/nu_mol ! no limiting for negative kn needed 
	fac = MIN(MAX((60.-kplus)/55.,0.),1.) ! Crimaldi 2006 correction for rough surfaces	
	do tel =1,nlayer ! apply 50 grid layers between wall z=0 and first velocity point at 0.5*dz 
		zplus = (zz(tel)+z0)*ust/nu_mol 
		nu_tot = nu_mol+nu_mol*kappa*zplus*(1.-fac*exp(-zplus/AA))**2 
		int1 = int1+1./nu_tot*ddzz(tel)
		int2 = int2+zzz/nu_tot*ddzz(tel) 
	enddo
	tau = rr/int1*(uu/rr-Fi*int2) 
!	if (tau*uu<0.) then ! they have different signs which would result in tau enhancing the flow instead of slowing it down
!		tau = ABS(rr/int1*(uu/rr))
!	else 
!		tau = ABS(tau)
!	endif 
	tau = MIN((rr*(0.25*ABS(uu)/rr)**2)/(ABS(tau)+1.e-9),1.)*tau !ust limited at maximum 0.25*u 

	IF (tau*uu<0.) THEN ! they have different signs --> tau enhancing the flow instead of slowing it down
		uu = uu - tau*scal*dt/dz 	!! only uu is adjusted
	ELSE 
		uu = uu / (1. + ABS(tau)*scal*dt/dz/MAX(ABS(uu),1.e-9)) 	!! only uu is adjusted ! implicit = more stable
	ENDIF 
	end
	
	subroutine wall_fun_rho_VD(uu,vv,rr,dz,scal,dt,kn,kappa,nu_mol,ust,tau)
		
	implicit none

c*************************************************************
c	Wall function based on Van Driest velocity profile  
c   Based on (Roulund et al. 2005) Numerical and experimental investigation of ﬂow and scour around a circular pile
c	Determine tau-wall and adapt uu	
c	Substract tau_wall*dt*dx*dy/(dx*dy*dz) from uu (kg/m3*m/s)
c*************************************************************
	real uu,vv,rr,dz,dt,ust,z0,kn,kappa,yplus,tau,nu_mol,absU
	real int1,zzz,AA,nu_tot,zplus,ddzz(1:50),zz(1:50),scal,fac,kplus,dzplus
	integer tel,tel0,nlayer
	
	nlayer = 50
	AA = 25.
	! Van Driest damping factor (Wang and Moin 2002) and VanBalen use value 19, others often use 26, Roulund et al. use 25. In my matlab tests 26 gives better velocity profiles, but 19 gives better (closer to log law) tau based on velocity at certain distance.	
	ddzz(1) = dz/100. !0.5*dz/50. -->input velocity uu,vv defined at 0.5*dz from wall
	do tel =2,nlayer ! apply 50 grid layers between wall z=0 and first velocity point at 0.5*dz with grow factor 
		ddzz(tel)=ddzz(tel-1)*1.1 
	enddo
	ddzz=0.5*dz/SUM(ddzz)*ddzz 
	zz(1) = 0.5*ddzz(1)
	do tel =2,nlayer 
		zz(tel)=zz(tel-1)+0.5*ddzz(tel-1)+0.5*ddzz(tel) 
	enddo	
	int1=0.
!	absU=sqrt((uu/rr)**2+(vv/rr)**2)
!	if (ust<0.0001*absU) then 		
!		ust=0.1*absU
!		if (kn>0.) then !walls with rougness (log-law used which can be hydr smooth or hydr rough):
!			do tel=1,10 ! 10 iter is more than enough
!				z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
!				ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
!			enddo
!		else
!			do tel=1,10 ! 10 iter is more than enough
!				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
!				ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
!			enddo
!			if (yplus<30.) then
!			  do tel=1,10 ! 10 iter is more than enough
!				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
!				ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
!			  enddo	
!			endif
!			if (yplus<5.) then !viscous sublayer uplus=yplus
!				ust=sqrt(absU*nu_mol/(0.5*dz))
!			endif
!		endif	  
!	endif 
	if (kn>0.) then
		z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
	else
		z0=0.
	endif
	kplus = MAX(5.,MIN(2000.,kn*ust/nu_mol)) 
	dzplus = 0.9*(sqrt(kplus)-kplus*exp(-kplus/6.))
	do tel=1,nlayer
		zplus=zz(tel)*ust/nu_mol+dzplus 
		int1=int1+(1./(1.+sqrt(1.+4.*kappa**2*zplus**2*(1.-exp(-zplus/AA))**2)))*ddzz(tel)*ust/nu_mol
	enddo
	tau = rr*((ABS(uu)/rr)/(2.*int1+1.e-9))**2 
	tau = MIN((rr*(0.25*ABS(uu)/rr)**2)/(ABS(tau)+1.e-9),1.)*tau !ust limited at maximum 0.25*u 
	uu = uu / (1. + tau*scal*dt/dz/MAX(ABS(uu),1.e-9)) 	!! only uu is adjusted ! implicit = more stable

	end
	
	

	subroutine wall_fun_TBL_Ploc(uu,vv,Fix,Fiy,rr,dz,kn,kappa,nu_mol,ust,ustnew)
		
	implicit none

c*************************************************************
c	Wall function based on simplified TBL with influence of dpdx 
c   Based on (Wang and Moin 2002) Dynamic wall modeling for large-eddy simulation of complex turbulent flows
c	Determine ust for Pressure location 
c*************************************************************
	real uu,vv,Fix,Fiy,rr,dz,kn,kappa,nu_mol,ust !in coming variable
	real ustnew !out going variable
	real int1,int2,zzz,AA,nu_tot,zplus,ddzz(1:50),zz(1:50),fac,kplus,z0,yplus,absU,taux,tauy !local variable
	integer tel,tel0,nlayer
	
	nlayer=50
	AA = 19. ! Van Driest damping factor (Wang and Moin 2002) and VanBalen use value 19, others often use 26. In my matlab tests 26 gives better velocity profiles, but 19 gives better (closer to log law) tau based on velocity at certain distance.
	ddzz(1) = dz/100. !0.5*dz/50. -->input velocity uu,vv defined at 0.5*dz from wall
	do tel =2,nlayer ! apply 50 grid layers between wall z=0 and first velocity point at 0.5*dz with grow factor 
		ddzz(tel)=ddzz(tel-1)*1.1 
	enddo
	ddzz=0.5*dz/SUM(ddzz)*ddzz 
	zz(1) = 0.5*ddzz(1)
	do tel =2,nlayer 
		zz(tel)=zz(tel-1)+0.5*ddzz(tel-1)+0.5*ddzz(tel) 
	enddo	
	
	int1=0.
	int2=0. 
!	absU=sqrt(uu**2+vv**2)
!	if (ust<0.0001*absU) then 		
!		ust=0.1*absU
!		! wall_fun_rho_TBL_Ploc is only used for sediment with kn>0
!!		if (kn>0.) then !walls with rougness (log-law used which can be hydr smooth or hydr rough):
!			do tel=1,10 ! 10 iter is more than enough
!				z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
!				ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
!			enddo
!!		else
!!			do tel=1,10 ! 10 iter is more than enough
!!				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
!!				ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
!!			enddo
!!			if (yplus<30.) then
!!			  do tel=1,10 ! 10 iter is more than enough
!!				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
!!				ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
!!			  enddo	
!!			endif
!!			if (yplus<5.) then !viscous sublayer uplus=yplus
!!				ust=sqrt(absU*nu_mol/(0.5*dz))
!!			endif
!!		endif	  
!	endif 
!	if (kn>0.) then
		z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
!	else
!		z0=0.
!	endif
	! I found a nice relation between 1/kappa*log(y/y0) for rough walls and VanDriest velocity profile with z0 added to z in zplus in combination with Crimaldi 2006 correction 
	kplus = kn*ust/nu_mol ! no limiting for negative kn needed 
	fac = MIN(MAX((60.-kplus)/55.,0.),1.) ! Crimaldi 2006 correction for rough surfaces
	do tel =1,nlayer 
		zplus = (zz(tel)+z0)*ust/nu_mol 		
		nu_tot = nu_mol+nu_mol*kappa*zplus*(1.-fac*exp(-zplus/AA))**2  
		int1 = int1+1./nu_tot*ddzz(tel)
		int2 = int2+zzz/nu_tot*ddzz(tel) 
	enddo
	taux = rr/int1*(uu-Fix*int2) 
!	if (taux*uu<0.) then ! they have different signs which would result in tau enhancing the flow instead of slowing it down
!		taux = ABS(rr/int1*(uu))
!	else 
!		taux = ABS(taux)
!	endif 
	tauy = rr/int1*(vv-Fiy*int2) 
!	if (tauy*vv<0.) then ! they have different signs which would result in tau enhancing the flow instead of slowing it down
!		tauy = ABS(rr/int1*(vv))
!	else 
!		tauy = ABS(tauy)
!	endif 
	ustnew=((taux/rr)**2+(tauy/rr)**2)**0.25
	ustnew = MIN(ustnew,0.25*sqrt(uu*uu+vv*vv)) !ust limited at maximum 0.25*u 
	
	end
	
	
		subroutine wall_fun_VD_Ploc(uu,vv,rr,dz,kn,kappa,nu_mol,ust,ustnew)
		
	implicit none

c*************************************************************
c	Wall function based on Van Driest velocity profile  
c   Based on (Roulund et al. 2005) Numerical and experimental investigation of ﬂow and scour around a circular pile
c	Determine ust for Pressure location 
c*************************************************************
	real uu,vv,rr,dz,kn,kappa,nu_mol,ust !in coming variable
	real ustnew !out going variable
	real int1,AA,nu_tot,zplus,ddzz(1:50),zz(1:50),kplus,z0,yplus,dzplus,absU !local variable
	integer tel,tel0,nlayer
	
	nlayer=50
	AA = 25.
	! Van Driest damping factor (Wang and Moin 2002) and VanBalen use value 19, others often use 26, Roulund et al. use 25. In my matlab tests 26 gives better velocity profiles, but 19 gives better (closer to log law) tau based on velocity at certain distance.	
	ddzz(1) = dz/100. !0.5*dz/50. -->input velocity uu,vv defined at 0.5*dz from wall
	do tel =2,nlayer ! apply 50 grid layers between wall z=0 and first velocity point at 0.5*dz with grow factor 
		ddzz(tel)=ddzz(tel-1)*1.1 
	enddo
	ddzz=0.5*dz/SUM(ddzz)*ddzz 
	zz(1) = 0.5*ddzz(1)
	do tel =2,nlayer 
		zz(tel)=zz(tel-1)+0.5*ddzz(tel-1)+0.5*ddzz(tel) 
	enddo	
	
	int1=0.
!	absU=sqrt(uu**2+vv**2)
!	if (ust<0.0001*absU) then 		
!		ust=0.1*absU
!		! wall_fun_rho_TBL_Ploc is only used for sediment with kn>0
!!		if (kn>0.) then !walls with rougness (log-law used which can be hydr smooth or hydr rough):
!			do tel=1,10 ! 10 iter is more than enough
!				z0=kn/30.+0.11*nu_mol/MAX(ust,1.e-9)
!				ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
!			enddo
!!		else
!!			do tel=1,10 ! 10 iter is more than enough
!!				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
!!				ust=absU/MAX((2.5*log(yplus)+5.5),2.) !ust maximal 0.5*absU			
!!			enddo
!!			if (yplus<30.) then
!!			  do tel=1,10 ! 10 iter is more than enough
!!				yplus=MAX(0.5*dz*ust/nu_mol,1e-12)
!!				ust=absU/MAX((5.*log(yplus)-3.05),2.) !ust maximal 0.5*absU			
!!			  enddo	
!!			endif
!!			if (yplus<5.) then !viscous sublayer uplus=yplus
!!				ust=sqrt(absU*nu_mol/(0.5*dz))
!!			endif
!!		endif	  
!	endif 
	kplus = MAX(5.,MIN(2000.,kn*ust/nu_mol)) 
	dzplus = 0.9*(sqrt(kplus)-kplus*exp(-kplus/6.))
	do tel=1,nlayer
		zplus=zz(tel)*ust/nu_mol+dzplus 
		int1=int1+(1./(1.+sqrt(1.+4.*kappa**2*zplus**2*(1.-exp(-zplus/AA))**2)))*ddzz(tel)*ust/nu_mol
	enddo
	ustnew=sqrt((uu/(2.*int1+1.e-9))**2+(vv/(2.*int1+1.e-9))**2)
	ustnew = MIN(ustnew,0.25*sqrt(uu*uu+vv*vv)) !ust limited at maximum 0.25*u
	
	end
	

	subroutine update_nvol_bedplume(tt)
	
	USE nlist
	
	implicit none
	
	real tt
	integer n2,inout,n,tel
	real xTSHD(4),yTSHD(4),phi,interpseries,xx,yy
	real fbx2,fbx,fby2,fby,fbz2
	real Propx_dummy(0:i1,0:px*jmax+1,1:kmax)
	real Propy_dummy(0:i1,0:px*jmax+1,1:kmax)
	real Propz_dummy(0:i1,0:px*jmax+1,1:kmax)

	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif	
!	DO n2=1,nbedplume !make all forces zero before new forces are applied for bedplume
!		IF (bp(n2)%u.ne.-99999.) THEN ! apply bedplume velocity boundary condition: !(bp(n2).h_tseriesfile.ne.''.or.bp(n2).zb_tseriesfile.ne.''.or.bp(n2).nmove>0).and.
!			IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end)) THEN
!				Propx_dummy=0.
!				Propy_dummy=0.
!				Propz_dummy=0.
!			ENDIF
!		ENDIF 
!	ENDDO	
!	DO n2=1,nbedplume !make all forces zero before new forces are applied for bedplume
!!		IF (bp(n2).h_tseriesfile.ne.''.or.bp(n2).zb_tseriesfile.ne.''.or.bp(n2).nmove>0) THEN ! apply bedplume velocity boundary condition:	
!		   IF (bp(n2)%u.ne.-99999.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end) THEN
!			do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz))! do k=k1,0,-1 !from top to bottom
!			  do i=0,i1  
!				do j=0,jmax*px+1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
!					xx=Rp(i)*cos_ut(j)-schuif_x
!					yy=Rp(i)*sin_ut(j)
!					xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
!					yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!					CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
!					if (inout.eq.1) then		   
!					  fbx2=ABS(bp(n2)%u)*bp(n2)%u/(bp(n2)%x(2)-bp(n2)%x(1))
!					  fby2=ABS(bp(n2)%v)*bp(n2)%v/(bp(n2)%y(3)-bp(n2)%y(2))
!					  fbz2=ABS(bp(n2)%w)*bp(n2)%w/(bp(n2)%height-bp(n2)%zbottom)
!					  fbx   =  fbx2 * cos(phi) - fby2 * sin(phi)
!					  fby   =  fbx2 * sin(phi) + fby2 * cos(phi)		
!					  Propx_dummy(i,j,k) = Propx_dummy(i,j,k)+0.5*(cos_ut(j)*fbx+sin_ut(j)*fby)
!					  Propx_dummy(MAX(i-1,0),j,k) = Propx_dummy(MAX(i-1,0),j,k) + 0.5*(cos_ut(j)*fbx+sin_ut(j)*fby)
!					  Propy_dummy(i,j,k) = Propy_dummy(i,j,k)+0.5*(-sin_vt(j)*fbx+cos_vt(j)*fby)
!					  Propy_dummy(i,MAX(j-1,0),k) = Propy_dummy(i,MAX(j-1,0),k) + 0.5*(-sin_vt(MAX(j-1,0))*fbx+cos_vt(MAX(j-1,0))*fby)
!					  Propz_dummy(i,j,k) = Propz_dummy(i,j,k)+0.5*fbz2
!					  Propz_dummy(i,j,MAX(k-1,0)) = Propz_dummy(i,j,MAX(k-1,0)) + 0.5*fbz2	
!					endif
!				enddo
!			  enddo
!			enddo
!		   ENDIF
!!		ENDIF 
!	ENDDO
!	DO n2=1,nbedplume !make all forces zero before new forces are applied for bedplume
!		IF (bp(n2)%u.ne.-99999.) THEN ! apply bedplume velocity boundary condition: !(bp(n2).h_tseriesfile.ne.''.or.bp(n2).zb_tseriesfile.ne.''.or.bp(n2).nmove>0).and.
!			IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end)) THEN	
!				Ppropx(0:i1,0:j1,1:kmax)=Propx_dummy(0:i1,rank*jmax:rank*jmax+j1,1:kmax)*rhU(0:i1,0:j1,1:kmax)
!				Ppropy(0:i1,0:j1,1:kmax)=Propy_dummy(0:i1,rank*jmax:rank*jmax+j1,1:kmax)*rhV(0:i1,0:j1,1:kmax)
!				Ppropz(0:i1,0:j1,1:kmax)=Propz_dummy(0:i1,rank*jmax:rank*jmax+j1,1:kmax)*rhW(0:i1,0:j1,1:kmax)
!		   ENDIF
!		ENDIF 
!	ENDDO
	
	DO n2=1,nbedplume !make all forces zero before new forces are applied for bedplume
		IF (bp(n2)%u.ne.-99999.) THEN ! apply bedplume velocity boundary condition: !(bp(n2).h_tseriesfile.ne.''.or.bp(n2).zb_tseriesfile.ne.''.or.bp(n2).nmove>0).and.
			IF ((bp(n2)%forever.eq.1.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end)) THEN
				Ppropx=0.
				Ppropy=0.
				Ppropz=0.
			ENDIF
		ENDIF 
	ENDDO	
	DO n2=1,nbedplume !make all forces zero before new forces are applied for bedplume
!		IF (bp(n2).h_tseriesfile.ne.''.or.bp(n2).zb_tseriesfile.ne.''.or.bp(n2).nmove>0) THEN ! apply bedplume velocity boundary condition:	
		   IF (bp(n2)%u.ne.-99999.and.time_np.gt.bp(n2)%t0.and.time_np.lt.bp(n2)%t_end) THEN
			do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz))! do k=k1,0,-1 !from top to bottom
			   do tel=1,bp(n2)%tmax 
				 i=bp(n2)%i(tel) 
				 j=bp(n2)%j(tel) 			
!!			  do i=0,i1  
!!				do j=0,j1 !jmax*px+1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
!!					xx=Rp(i)*cos_u(j)-schuif_x
!!					yy=Rp(i)*sin_u(j)
!!					xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
!!					yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
!!					CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
!!					if (inout.eq.1) then		   
					  fbx2=ABS(bp(n2)%u)*bp(n2)%u/(bp(n2)%x(2)-bp(n2)%x(1))
					  fby2=ABS(bp(n2)%v)*bp(n2)%v/(bp(n2)%y(3)-bp(n2)%y(2))
					  fbz2=ABS(bp(n2)%w)*bp(n2)%w/(bp(n2)%height-bp(n2)%zbottom)
					  fbx   =  fbx2 * cos(phi) - fby2 * sin(phi)
					  fby   =  fbx2 * sin(phi) + fby2 * cos(phi)		
					  Ppropx(i,j,k) = Ppropx(i,j,k)+0.5*(cos_u(j)*fbx+sin_u(j)*fby)*rhU(i,j,k)
					  Ppropx(MAX(i-1,0),j,k) = Ppropx(MAX(i-1,0),j,k) + 0.5*(cos_u(j)*fbx+sin_u(j)*fby)*rhU(MAX(i-1,0),j,k)
					  Ppropy(i,j,k) = Ppropy(i,j,k)+0.5*(-sin_v(j)*fbx+cos_v(j)*fby)*rhV(i,j,k)
					  Ppropy(i,MAX(j-1,0),k) = Ppropy(i,MAX(j-1,0),k) + 0.5*(-sin_v(MAX(j-1,0))*fbx+cos_v(MAX(j-1,0))*fby)*rhV(i,MAX(j-1,0),k)
					  Ppropz(i,j,k) = Ppropz(i,j,k)+0.5*fbz2*rhW(i,j,k)
					  Ppropz(i,j,MAX(k-1,0)) = Ppropz(i,j,MAX(k-1,0)) + 0.5*fbz2*rhW(i,j,MAX(k-1,0))
!!					endif
!!				enddo
			  enddo
			enddo
		   ENDIF
!		ENDIF 
	ENDDO
	
	DO n2=1,nbedplume
		IF (bp(n2).h_tseriesfile.ne.''.or.bp(n2).zb_tseriesfile.ne.''.or.bp(n2).nmove>0) THEN !necessary to update volncells!
		  IF (bp(n2).h_tseriesfile.ne.'') THEN
			bp(n2)%height=interpseries(bp(n2)%h_tseries,bp(n2)%h_series,bp(n2)%h_seriesloc,tt)
		  ENDIF
		  IF (bp(n2).zb_tseriesfile.ne.'') THEN
			bp(n2)%zbottom=interpseries(bp(n2)%zb_tseries,bp(n2)%zb_series,bp(n2)%zb_seriesloc,tt)
		  ENDIF	
		  bp(n2)%volncells=0.
		  do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
		   do i=1,imax 
			 do j=1,jmax*px        
		  xx=Rp(i)*cos_ut(j)-schuif_x !global xx over different processors
		  yy=Rp(i)*sin_ut(j)          !global yy over different processors
!		  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle: 
			xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
			yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
		    if (bp(n2)%radius.gt.0.) then
			  inout=0
		      IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			    inout=1
			  ENDIF
			else 
			  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
			endif 
!		  ELSE 
!			inout=0
!		  ENDIF
		  if (inout.eq.1) then
		   bp(n2)%volncells=bp(n2)%volncells+vol_V(i,j)*fc_global(i,j,k)
		   endif
		  enddo
		 enddo
		enddo	  
!		 if (bp(n2)%volncells.le.0.and.tt<bp(n2)%t_end) then
!		   write(*,*),'WARNING, no cells found for bedplume number: ',n2,bp(n2)%volncells,rank
!		   write(*,*),'In case a sedflux or Q has been defined the code will crash because of division by a zero volncells'
!		 endif
		ENDIF ! bedplume loop
	ENDDO	
	
	end
	
	subroutine update_his_bedplume(tt)
	
	USE nlist
	
	implicit none
	
	include 'mpif.h'
	real tt,globalsum
	integer n2,inout,n,ierr,tel
	real xTSHD(4),yTSHD(4),phi,xx,yy

	! rotation ship for ambient side current
	if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
	  phi=atan2(V_b,1.e-12)
	else
	  phi=atan2(V_b,(U_TSHD-U_b))
	endif	
	DO n2=1,nbedplume
		IF (bp(n2)%dt_history>0.AND.tt.ge.bp(n2)%t_bphis_output.AND.bp(n2)%t_bphis_output.le.te_output+1e-12
     &        .and.bp(n2)%istep_bphis_output.le.20000) THEN 
		bp(n2)%istep_bphis_output=bp(n2)%istep_bphis_output+1
		bp(n2)%t_bphis_output=bp(n2)%t_bphis_output+bp(n2)%dt_history	 
		  do k=MAX(1,CEILING(bp(n2)%zbottom/dz)),MIN(kmax,FLOOR(bp(n2)%height/dz)) ! 1,kmax
		   do i=1,imax 
			 do j=1,jmax        
		  xx=Rp(i)*cos_u(j)-schuif_x 
		  yy=Rp(i)*sin_u(j)          
!		  IF (k.le.FLOOR(bp(n2)%height/dz).and.k.ge.CEILING(bp(n2)%zbottom/dz)) THEN ! obstacle: 
			xTSHD(1:4)=bp(n2)%x*cos(phi)-bp(n2)%y*sin(phi)
			yTSHD(1:4)=bp(n2)%x*sin(phi)+bp(n2)%y*cos(phi)
		    if (bp(n2)%radius.gt.0.) then
			  inout=0
		      IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n2)%radius)**2) THEN
			    inout=1
			  ENDIF
			else 
			  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
			endif 
!		  ELSE 
!			inout=0
!		  ENDIF
		  if (inout.eq.1) then
		     DO n=1,nfrac
				Chisbp(n,n2,bp(n2)%istep_bphis_output)=Chisbp(n,n2,bp(n2)%istep_bphis_output)+
     &				vol_V(i,j+rank*jmax)*cnew(n,i,j,k)*fc_global(i,j+jmax*rank,k)		
		     ENDDO
				Uhisbp(n2,bp(n2)%istep_bphis_output)=Uhisbp(n2,bp(n2)%istep_bphis_output)+
     &				vol_V(i,j+rank*jmax)*Unew(i,j,k)*fc_global(i,j+jmax*rank,k)				 
				Vhisbp(n2,bp(n2)%istep_bphis_output)=Vhisbp(n2,bp(n2)%istep_bphis_output)+
     &				vol_V(i,j+rank*jmax)*Vnew(i,j,k)*fc_global(i,j+jmax*rank,k)				 
				Whisbp(n2,bp(n2)%istep_bphis_output)=Whisbp(n2,bp(n2)%istep_bphis_output)+
     &				vol_V(i,j+rank*jmax)*Wnew(i,j,k)*fc_global(i,j+jmax*rank,k)		
		   endif
		  enddo
		 enddo
		enddo
		DO n=1,nfrac
		  call mpi_allreduce(Chisbp(n,n2,bp(n2)%istep_bphis_output),globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		  Chisbp(n,n2,bp(n2)%istep_bphis_output)=globalsum/bp(n2)%volncells
		ENDDO 
		  call mpi_allreduce(Uhisbp(n2,bp(n2)%istep_bphis_output),globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		  Uhisbp(n2,bp(n2)%istep_bphis_output)=globalsum/bp(n2)%volncells		
		  call mpi_allreduce(Vhisbp(n2,bp(n2)%istep_bphis_output),globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		  Vhisbp(n2,bp(n2)%istep_bphis_output)=globalsum/bp(n2)%volncells		
		  call mpi_allreduce(Whisbp(n2,bp(n2)%istep_bphis_output),globalsum,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
		  Whisbp(n2,bp(n2)%istep_bphis_output)=globalsum/bp(n2)%volncells				  
		thisbp(n2,bp(n2)%istep_bphis_output)=tt 
		zhisbp(n2,bp(n2)%istep_bphis_output)=(MAX(1,CEILING(bp(n2)%zbottom/dz))+MIN(kmax,FLOOR(bp(n2)%height/dz)))/2.*dz+0.5*dz 			
		ENDIF ! bedplume loop

	ENDDO
	
	end
	
	
	subroutine update_QSc_bedplume(tt)
	
	USE nlist
	
	implicit none
	
	real tt
	integer n2,inout,n3
	real interpseries,Utot
	
	DO n2=1,nbedplume
		  IF (bp(n2).Q_tseriesfile.ne.'') THEN
			bp(n2)%Q=interpseries(bp(n2)%Q_tseries,bp(n2)%Q_series,bp(n2)%Q_seriesloc,tt)
		  ENDIF
		  IF (bp(n2).u_tseriesfile.ne.'') THEN
			bp(n2)%u=interpseries(bp(n2)%u_tseries,bp(n2)%u_series,bp(n2)%u_seriesloc,tt)
			bp(n2)%uinput=bp(n2)%u
		  ENDIF		  
		  IF (bp(n2).v_tseriesfile.ne.'') THEN
			bp(n2)%v=interpseries(bp(n2)%v_tseries,bp(n2)%v_series,bp(n2)%v_seriesloc,tt)
		  ENDIF			  
		  IF (bp(n2).w_tseriesfile.ne.'') THEN
			bp(n2)%w=interpseries(bp(n2)%w_tseries,bp(n2)%w_series,bp(n2)%w_seriesloc,tt)
		  ENDIF			  
		  IF (bp(n2).nmove.gt.0.and.bp(n2)%u.ne.-99999.) THEN 
			Utot = SQRT(bp(n2)%uinput**2)
			bp(n2)%u = Utot*bp(n2)%move_nx_series(MAX(bp(n2)%nmove_present,1))
			bp(n2)%v = Utot*bp(n2)%move_ny_series(MAX(bp(n2)%nmove_present,1))			
		  ENDIF
		  IF (bp(n2).S_tseriesfile.ne.'') THEN
			CALL interpseries3(bp(n2)%S_tseries,bp(n2)%S_series(1:nfrac,:),bp(n2)%S_seriesloc,tt,nfrac,bp(n2)%sedflux(1:nfrac))
		  ENDIF	
		  IF (bp(n2).c_tseriesfile.ne.'') THEN
			CALL interpseries3(bp(n2)%c_tseries,bp(n2)%c_series(1:nfrac,:),bp(n2)%c_seriesloc,tt,nfrac,bp(n2)%c(1:nfrac))
		  ENDIF			  
	ENDDO	
	
	end
	
	subroutine update_Qc_plume(tt)
	
	USE nlist
	
	implicit none
	
	real tt
	integer n2,inout,n3
	real interpseries,Qplume
	
		  IF (plumeQtseriesfile.ne.'') THEN
			Q_j=interpseries(plumeQtseries,plumeQseries,plumeQseriesloc,tt)
			W_j=-Q_j/Aplume
		  ENDIF
		  IF (plumectseriesfile.ne.'') THEN
			CALL interpseries3(plumectseries,plumecseries(1:nfrac,:),plumecseriesloc,tt,nfrac,frac(1:nfrac)%c)
		  ENDIF			  
	
	end	
	
	
	subroutine update_location_bedplume 
	
	USE nlist
	
	implicit none
	include 'mpif.h'
	
	integer inout,n,ierr,tel
	real xx,yy,phi,xTSHD(1:4),yTSHD(1:4),zbed_max,zbed_mean,Abed_bedplume,Abed_bedplume_tot,zbed_tot,zb,Utot
	
	DO n=1,nbedplume
	  IF ((bp(n)%forever.eq.1.and.time_np.gt.bp(n)%t0.and.time_np.lt.bp(n)%t_end).and.
     & bp(n)%move_zbed_criterium(MAX(bp(n)%nmove_present,1))<depth) THEN
		! rotation ship for ambient side current
		if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
		IF (bp(n)%move_zbed_type(MAX(bp(n)%nmove_present,1)).eq.2) THEN !check for avg bed level:
			zbed_mean=0.
			Abed_bedplume=0.
			do i=1,imax !0,i1  
			 do j=jmax,1,-1 !j1,0,-1       
				xx=Rp(i)*cos_u(j)-schuif_x
				yy=Rp(i)*sin_u(j)
				if (bp(n)%radius2.gt.0.) then !use bp()%x(1),y(1) and radius for zbed-check
				  xTSHD(1:4)=bp(n)%x2*cos(phi)-bp(n)%y2*sin(phi)
				  yTSHD(1:4)=bp(n)%x2*sin(phi)+bp(n)%y2*cos(phi)			
				  inout=0
				  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n)%radius2)**2) THEN
					inout=1
				  ENDIF
				else 
				  xTSHD(1:4)=bp(n)%x2*cos(phi)-bp(n)%y2*sin(phi)
				  yTSHD(1:4)=bp(n)%x2*sin(phi)+bp(n)%y2*cos(phi)			
				  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
				endif 			
				if (inout.eq.1) then
					zb=MAX(REAL(kbed(i,j))*bednotfixed(i,j,kbed(i,j)),0.)*dz !use zb = kbed*dz as this corresponds with bed cells that are blocked in flow, buffer layer is not real bed but bookkeeping
					! only count zb if top cell bed (kbed) is erodable				
					zbed_mean=zbed_mean+zb*vol_V(i,j+rank*jmax)*bednotfixed(i,j,kbed(i,j)) 
					Abed_bedplume=Abed_bedplume+vol_V(i,j+rank*jmax)*bednotfixed(i,j,kbed(i,j)) 
				endif
			 enddo
			enddo
			call mpi_allreduce(zbed_mean,zbed_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			call mpi_allreduce(Abed_bedplume,Abed_bedplume_tot,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
			zbed_tot=zbed_tot/Abed_bedplume_tot
		ELSE !check for maximum bed level:
			zbed_max=0.
			do i=1,imax !0,i1  
			 do j=jmax,1,-1 !j1,0,-1       
				xx=Rp(i)*cos_u(j)-schuif_x
				yy=Rp(i)*sin_u(j)
				if (bp(n)%radius2.gt.0.) then !use bp()%x(1),y(1) and radius for zbed-check
				  xTSHD(1:4)=bp(n)%x2*cos(phi)-bp(n)%y2*sin(phi)
				  yTSHD(1:4)=bp(n)%x2*sin(phi)+bp(n)%y2*cos(phi)			
				  inout=0
				  IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n)%radius2)**2) THEN
					inout=1
				  ENDIF
				else 
				  xTSHD(1:4)=bp(n)%x2*cos(phi)-bp(n)%y2*sin(phi)
				  yTSHD(1:4)=bp(n)%x2*sin(phi)+bp(n)%y2*cos(phi)			
				  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
				endif 			
				if (inout.eq.1) then
					!zb=REAL(MAX(kbed(i,j)-1,0))*dz+(SUM(dcdtbot(1:nfrac,i,j))+SUM(Clivebed(1:nfrac,i,j,kbed(i,j))))/cfixedbed*dz
					!zb=REAL(MAX(kbed(i,j),0))*dz !use zb = kbed*dz as this corresponds with bed cells that are blocked in flow, buffer layer is not real bed but bookkeeping
					zb=MAX(REAL(kbed(i,j))*bednotfixed(i,j,kbed(i,j)),0.)*dz !use zb = kbed*dz as this corresponds with bed cells that are blocked in flow, buffer layer is not real bed but bookkeeping
					! only count zb if top cell bed (kbed) is erodable				
					zbed_max=MAX(zbed_max,zb)
				endif
			 enddo
			enddo
			call mpi_allreduce(zbed_max,zbed_tot,1,mpi_real8,mpi_max,mpi_comm_world,ierr)
		ENDIF 
		IF(zbed_tot>bp(n)%move_zbed_criterium(MAX(bp(n)%nmove_present,1))) THEN
			bp(n)%nmove_present=bp(n)%nmove_present+1
			istep_output_bpmove = istep_output_bpmove +1
			IF ((bp(n)%nmove_present-bp(n)%nmove).eq.1) THEN
				IF (rank.eq.0) THEN
					write(*,'(a,i6,a)'),'* WARNING, bedplume : ',n,' should have been moved, because'
					write(*,*),' move_zbed_criterium has been met, but move_dx_series is not long enough'
					write(*,*),' simulation is continued with bedplume at last location, '
					write(*,*),' a mvbp output file is created at this moment now move_zbed_criterium has been reached'
				ENDIF
				call output_nc('mvbp3D_',istep_output_bpmove,time_np)
			ELSEIF (bp(n)%nmove_present>bp(n)%nmove) THEN
				! move_dx_series not long enough: do nothing and leave bp at present location
			ELSE
				bp(n)%x = bp(n)%x+bp(n)%move_dx_series(bp(n)%nmove_present)
				bp(n)%y = bp(n)%y+bp(n)%move_dy_series(bp(n)%nmove_present)
				bp(n)%height = bp(n)%height+bp(n)%move_dz_series(bp(n)%nmove_present)*bp(n)%move_dz_height_factor
				bp(n)%zbottom = bp(n)%zbottom+bp(n)%move_dz_series(bp(n)%nmove_present)*bp(n)%move_dz_zbottom_factor
				bp(n)%x2 = bp(n)%x2+bp(n)%move_dx2_series(bp(n)%nmove_present)
				bp(n)%y2 = bp(n)%y2+bp(n)%move_dy2_series(bp(n)%nmove_present)				
				IF (bp(n)%u.ne.-99999.) THEN
					Utot = SQRT(bp(n)%uinput**2)
					bp(n)%u = Utot*bp(n)%move_nx_series(bp(n)%nmove_present)
					bp(n)%v = Utot*bp(n)%move_ny_series(bp(n)%nmove_present)
				ENDIF
				
				IF (rank.eq.0) THEN
				  if (bp(n)%move_zbed_type(MAX(bp(n)%nmove_present-1,1)).eq.2) then 
				    write(*,'(a,i6,a,f8.4,f8.4,a,i4,a)'),'* Bedplume : ',n,' with avg bedlevel and criterium: ',zbed_tot,
     &			    bp(n)%move_zbed_criterium(MAX(bp(n)%nmove_present-1,1)),' has been moved for the ',bp(n)%nmove_present,' time.'
				  else 
					write(*,'(a,i6,a,f8.4,f8.4,a,i4,a)'),'* Bedplume : ',n,' with max bedlevel and criterium: ',zbed_tot,
     &			    bp(n)%move_zbed_criterium(MAX(bp(n)%nmove_present-1,1)),' has been moved for the ',bp(n)%nmove_present,' time.'				  
				  endif 
				  write(*,'(a,f11.2,f11.2,f11.2,f11.2)'),'Location after move      : x=',bp(n)%x
				  write(*,'(a,f11.2,f11.2,f11.2,f11.2)'),'Location after move      : y=',bp(n)%y
				  write(*,'(a,f11.2,f11.2,a,f8.4,a,i4)'),'Location after move      : z=',bp(n)%zbottom,bp(n)%height,
     &			  '; z-criterium after move=',bp(n)%move_zbed_criterium(MAX(bp(n)%nmove_present,1))
     & ,' and move_zbed_type after move=',bp(n)%move_zbed_type(MAX(bp(n)%nmove_present,1))
				  write(*,'(a,f11.2,f11.2,f11.2,f11.2)'),'Check bedlevel after move: x2=',bp(n)%x2
				  write(*,'(a,f11.2,f11.2,f11.2,f11.2)'),'Check bedlevel after move: y2=',bp(n)%y2
				ENDIF
				IF (bp(n)%move_outputfile_series(bp(n)%nmove_present).eq.1) THEN
					call output_nc('mvbp3D_',istep_output_bpmove,time_np)
				ENDIF
			ENDIF
		ENDIF
	  ENDIF
	ENDDO
	
	DO n=1,nbedplume
	  IF ((bp(n)%forever.eq.1.and.time_np.gt.bp(n)%t0.and.time_np.lt.bp(n)%t_end)) THEN
        bp(n)%x=bp(n)%x+dt*bp(n)%move_u 
		bp(n)%y=bp(n)%y+dt*bp(n)%move_v
        bp(n)%x2=bp(n)%x2+dt*bp(n)%move_u 
		bp(n)%y2=bp(n)%y2+dt*bp(n)%move_v		
		bp(n)%height=bp(n)%height+dt*bp(n)%move_w
		bp(n)%zbottom=bp(n)%zbottom+dt*bp(n)%move_w
	  ENDIF
	ENDDO
	
	DO n=1,nbedplume
	  IF ((bp(n)%forever.eq.1.and.time_np.gt.bp(n)%t0.and.time_np.lt.bp(n)%t_end).and.
     & (bp(n)%move_zbed_criterium(MAX(bp(n)%nmove_present,1))<depth.or.(ABS(bp(n)%move_u)+ABS(bp(n)%move_v).gt.0.))) THEN
	    bp(n)%tmax=0
		! rotation ship for ambient side current
		if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
		do i=0,i1  
		 do j=j1,0,-1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
			xx=Rp(i)*cos_u(j)-schuif_x
			yy=Rp(i)*sin_u(j)
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
			if (inout.eq.1) then
			  bp(n)%tmax = bp(n)%tmax+1
			  if (bp(n)%tmax.le.10000) then
				  bp(n)%i(bp(n)%tmax)=i 
				  bp(n)%j(bp(n)%tmax)=j
			  else 
				write(*,*),'Bedplume : ',n 
			    CALL writeerror(281)
			  endif 
			endif
		 enddo
		enddo
	  ENDIF
	ENDDO
	
	end	

	subroutine init_location_bedplume 
	
	USE nlist
	
	implicit none
	include 'mpif.h'
	
	integer inout,n,ierr
	real xx,yy,phi,xTSHD(1:4),yTSHD(1:4)
	
	DO n=1,nbedplume
	    bp(n)%tmax=0
		! rotation ship for ambient side current
		if ((U_TSHD-U_b).eq.0.or.LOA<0.) then
		  phi=atan2(V_b,1.e-12)
		else
		  phi=atan2(V_b,(U_TSHD-U_b))
		endif
		do i=0,i1  
		 do j=j1,0,-1       ! bedplume loop is only initial condition: do not bother to have U,V,W initial staggering perfect 
			xx=Rp(i)*cos_u(j)-schuif_x
			yy=Rp(i)*sin_u(j)
			xTSHD(1:4)=bp(n)%x*cos(phi)-bp(n)%y*sin(phi)
			yTSHD(1:4)=bp(n)%x*sin(phi)+bp(n)%y*cos(phi)			
		    if (bp(n)%radius.gt.0.) then
			  bp(n)%x(2)=bp(n)%x(1)+2.*bp(n)%radius !to get dpdx force correct for bp()%u
			  bp(n)%y(3)=bp(n)%y(2)+2.*bp(n)%radius !to get dpdy force correct for bp()%v			
			  inout=0
		      IF (((xx-xTSHD(1))**2+(yy-yTSHD(1))**2).lt.(bp(n)%radius)**2) THEN
			    inout=1
			  ENDIF
			else 
			  CALL PNPOLY (xx,yy, xTSHD(1:4), yTSHD(1:4), 4, inout ) 
			endif 
			if (inout.eq.1) then
			  bp(n)%tmax = bp(n)%tmax+1
			  if (bp(n)%tmax.le.10000) then
				  bp(n)%i(bp(n)%tmax)=i 
				  bp(n)%j(bp(n)%tmax)=j
			  else 
				write(*,*),'Bedplume : ',n 
			    CALL writeerror(281)
			  endif 
			endif
		 enddo
		enddo
	ENDDO
	
	end	
	
	subroutine update_fc_global 
	
	USE nlist
	
	implicit none
	include 'mpif.h'
	
	integer status4(MPI_STATUS_SIZE),t,r,j2,tel,ierr
	real*8 fc_local(1:imax,1:jmax,1:kmax),fc_local_vec(imax*jmax*kmax),fc_global_vec(imax*jmax*px*kmax)



	fc_local=1.
	  do t=1,tmax_inPpuntTSHD
 	    i=i_inPpuntTSHD(t)
 	    j=j_inPpuntTSHD(t)		
 	    k=k_inPpuntTSHD(t)		
		IF (i.ge.1.and.i.le.imax.and.j.ge.1.and.j.le.jmax.and.k.ge.1.and.k.le.kmax) THEN
			fc_local(i,j,k)=0.
		ENDIF
	  enddo
      
		do j=1,jmax
			do i=1,imax
				do k=1,kbed(i,j)
					fc_local(i,j,k)=0.
				enddo
			enddo
		enddo	

		do i=1,imax
			do j=1,jmax
				do k=1,kmax	
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
		fc_global(0:i1,0:jmax*px+1,0)=fc_global(0:i1,0:jmax*px+1,1)
		fc_global(0:i1,0:jmax*px+1,k1)=fc_global(0:i1,0:jmax*px+1,kmax)
	
	
	end	
	
	
	REAL function interpseries(tseries,series,sloc,tt)

	REAL tseries(1:10000),series(1:10000),tt
	INTEGER sloc

	!IF (tt>tseries(sloc+1)) THEN
	!	sloc=sloc+1
	!ENDIF	
	DO WHILE (tt>tseries(sloc+1)) 
		sloc=sloc+1
	ENDDO
	
!	write(*,*)'sloc,tseries(sloc),tseries(sloc+1)',sloc,tseries(sloc),tseries(sloc+1)
!	write(*,*)'series(sloc),series(sloc+1),interpval',series(sloc),series(sloc+1),
!     & (tt-tseries(sloc))/(tseries(sloc+1)-tseries(sloc))*(series(sloc+1)-series(sloc))+series(sloc)
	interpseries=(tt-tseries(sloc))/(tseries(sloc+1)-tseries(sloc))*(series(sloc+1)-series(sloc))+series(sloc)	

	end function interpseries


	subroutine interpseries3(tseries,series,sloc,tt,nfrac,output)
	implicit none

	REAL tseries(1:10000),series(1:nfrac,1:10000),tt,fac 
	REAL output(1:nfrac)
	INTEGER sloc,nfrac,n

	!IF (tt>tseries(sloc+1)) THEN
	!	sloc=sloc+1
	!ENDIF
	DO WHILE (tt>tseries(sloc+1)) 
		sloc=sloc+1
	ENDDO	
	fac=(tt-tseries(sloc))/(tseries(sloc+1)-tseries(sloc))
	DO n=1,nfrac
		output(n)=fac*(series(n,sloc+1)-series(n,sloc))+series(n,sloc)	
	ENDDO

	end	
	
	
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

	subroutine shiftb_k(UT,UP)

      USE nlist


      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf,ierr
!       parameter (ileng= (k1+1)*(i1+1) )
      include 'mpif.h'
      integer itag,status(MPI_STATUS_SIZE),l
      real*8 ut(0:j1,0:k1)
      real*8 up(0:k1),UTMP(0:K1)
	!real up(0:i1,0:k1),UTMP(0:I1,0:K1)
      
	 do k=0,k1
	  utmp(k) =UT(1,k)
          enddo
      
      itag = 10
      ileng = (k1+1)
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

	subroutine shiftf_k(UT,UP)
      USE nlist

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf
!       parameter (ileng= (k1+1)*(i1+1))
      include 'mpif.h'
      real*8 UT(0:j1,0:k1),UP(0:k1),UTMP(0:k1)
	!real UT(0:i1,0:j1,0:k1),UP(0:i1,0:k1),UTMP(0:i1,0:k1)
      integer  itag,status(MPI_STATUS_SIZE),l,ierr
      itag = 11
      ileng = (k1+1)
      do k=0,k1
	  UTMP(k) =UT(jmax,k)
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

	subroutine shiftb_intl(UT1,UP1)
      USE nlist


      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf,ierr
!       parameter (ileng= (k1+1)*(i1+1) )
      include 'mpif.h'
      integer itag1,status(MPI_STATUS_SIZE),l
      integer*8 ut1(0:i1,0:j1)
      integer*8 up1(0:i1),UTMP1(0:I1)
	!real up(0:i1,0:k1),UTMP(0:I1,0:K1)
      do i=0,i1
	  utmp1(i) =UT1(i,1)
      enddo
      itag1 = 121
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

      call mpi_sendrecv(utmp1 ,ileng,MPI_INTEGER8,rankb,itag1,
     $                  up1 ,ileng,MPI_INTEGER8,rankf,itag1, MPI_COMM_WORLD,status,ierr)

!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankb,itag,
!     $                  up ,ileng,MPI_REAL,rankf,itag, MPI_COMM_WORLD,status,ierr)
      end

	subroutine shiftf_intl(UT,UP)
      USE nlist

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf
!       parameter (ileng= (k1+1)*(i1+1))
      include 'mpif.h'
      integer*8 UT(0:i1,0:j1),UP(0:i1),UTMP(0:i1)
	!real UT(0:i1,0:j1,0:k1),UP(0:i1,0:k1),UTMP(0:i1,0:k1)
      integer  itag2,status(MPI_STATUS_SIZE),l,ierr
      itag2 = 131
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

      call mpi_sendrecv(utmp ,ileng,MPI_INTEGER8,rankf,itag2,
     $                  up   ,ileng,MPI_INTEGER8,rankb,itag2, MPI_COMM_WORLD,status,ierr)
!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankf,itag,
!     $                  up   ,ileng,MPI_REAL,rankb,itag, MPI_COMM_WORLD,status,ierr)

      end

	subroutine shiftb_lreverse(UT1,UP1)

      USE nlist

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf,ierr,n
!       parameter (ileng= (k1+1)*(i1+1) )
      include 'mpif.h'
      integer itag1,status(MPI_STATUS_SIZE),l
      real*8 ut1(1:nfrac,0:i1,0:j1)
      real*8 up1(1:nfrac,0:i1),UTMP1(1:nfrac,0:I1)
	!real up(0:i1,0:k1),UTMP(0:I1,0:K1)
	 do n=1,nfrac
      do i=0,i1
	  utmp1(n,i) =UT1(n,i,0) !UT1(i,1)
      enddo
	 enddo
      itag1 = 1112
      ileng = nfrac*(i1+1)
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
	 
      end

	subroutine shiftf_lreverse(UT,UP)
      USE nlist

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf,n
!       parameter (ileng= (k1+1)*(i1+1))
      include 'mpif.h'
      real*8 UT(1:nfrac,0:i1,0:j1),UP(1:nfrac,0:i1),UTMP(1:nfrac,0:i1)
	!real UT(0:i1,0:j1,0:k1),UP(0:i1,0:k1),UTMP(0:i1,0:k1)
      integer  itag2,status(MPI_STATUS_SIZE),l,ierr
      itag2 = 1113
      ileng = nfrac*(i1+1)
	 do n=1,nfrac	  
        do i=0,i1
	  UTMP(n,i) =UT(n,i,j1) !UT(i,jmax)
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

      call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankf,itag2,
     $                  up   ,ileng,MPI_REAL8,rankb,itag2, MPI_COMM_WORLD,status,ierr)
 
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

	subroutine shiftb_intl2(UT1,UP1)

      USE nlist


      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf,ierr
!       parameter (ileng= (k1+1)*(i1+1) )
      include 'mpif.h'
      integer itag1,status(MPI_STATUS_SIZE),l
      integer*8 ut1(0:i1,0:j1)
      integer*8 up1(0:i1),UTMP1(0:I1)
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

      call mpi_sendrecv(utmp1 ,ileng,MPI_INTEGER8,rankb,itag1,
     $                  up1 ,ileng,MPI_INTEGER8,rankf,itag1, MPI_COMM_WORLD,status,ierr)
!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankb,itag,
!     $                  up ,ileng,MPI_REAL,rankf,itag, MPI_COMM_WORLD,status,ierr)
      end

	subroutine shiftf_intl2(UT,UP)
      USE nlist

      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf
!       parameter (ileng= (k1+1)*(i1+1))
      include 'mpif.h'
      integer*8 UT(0:i1,0:j1),UP(0:i1),UTMP(0:i1)
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

      call mpi_sendrecv(utmp ,ileng,MPI_INTEGER8,rankf,itag2,
     $                  up   ,ileng,MPI_INTEGER8,rankb,itag2, MPI_COMM_WORLD,status,ierr)
!      call mpi_sendrecv(utmp ,ileng,MPI_REAL,rankf,itag,
!     $                  up   ,ileng,MPI_REAL,rankb,itag, MPI_COMM_WORLD,status,ierr)
	end 


	subroutine shiftb_T(UT,UP)

      USE nlist


      implicit none
!       include 'param.txt'
!       include 'common.txt'
      integer ileng,rankb,rankf,ierr
!       parameter (ileng= (k1+1)*(i1+1) )
      include 'mpif.h'
      integer itag,status(MPI_STATUS_SIZE),l
      real*8 ut(1:i1,0:jmax*px+1,1:kmax/px+1)
      real*8 up(1:i1,0:jmax*px+1),UTMP(1:I1,0:jmax*px+1)
	!real up(0:i1,0:k1),UTMP(0:I1,0:K1)
      do i=1,i1
	 do j=0,jmax*px+1
	  utmp(i,j) =UT(i,j,1)
          enddo
      enddo
      itag = 10
      ileng = (jmax*px+2)*(i1)
      rankf=rank+1
      rankb=rank-1
           if(rank.eq.px-1)rankf=MPI_PROC_NULL
           if(rank.eq.   0)rankb=MPI_PROC_NULL
      call mpi_sendrecv(utmp ,ileng,MPI_REAL8,rankb,itag,
     $                  up ,ileng,MPI_REAL8,rankf,itag, MPI_COMM_WORLD,status,ierr)

      end
	  