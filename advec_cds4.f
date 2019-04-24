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

      subroutine advecu_CDS4(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,wf,wd)
      implicit none
c
c********************************************************************
c
c     advecu calculates the advection of the u-velocity, which is
c     the velocity in the radial direction.
c
c     In formula:
c
c         1 d(ruu)     1 d(uv)     d(uw)     vv
c    - (  - ------  +  - -----  +  -----  -  --  )
c         r   dr       r  dphi      dz        r
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          Utmp              : contains velocity at oldest timestep
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1),numdif,dzi
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm

	real Ax,Bx,Cx,DDx
	real Ay,By,Cy,DDy
	real Az,Bz,Cz,DDz
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy,wd
	real*8 ubb(0:i1,0:k1),ubf(0:i1,0:k1),Uvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real wwdd(1:3),wf(0:i1,0:j1,0:k1),numdif2,u4mag
	real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1)
	real Vvel2(-2:i1+2,-2:j1+2,-2:k1+2),Wvel2(-2:i1+2,-2:j1+2,-2:k1+2)

	Uvel2(0:i1,0:j1,0:k1)=Uvel

!c get stuff from other CPU's
	  call shiftf2(Uvel,ubf)
	  call shiftb2(Uvel,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   Uvel2(i,-1,k) = Uvel(i,0,k)
		   Uvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   Uvel2(i,-1,k) = Ubf(i,k)
		   Uvel2(i,j1+1,k) =Uvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   Uvel2(i,-1,k) = Ubf(i,k)
		   Uvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   Uvel2(i,-1,k) = Ubf(i,k)
		   Uvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	  call shiftf3(Uvel,ubf)
	  call shiftb3(Uvel,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   Uvel2(i,-2,k) = Uvel(i,0,k)
		   Uvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   Uvel2(i,-2,k) = Ubf(i,k)
		   Uvel2(i,j1+2,k) =Uvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   Uvel2(i,-2,k) = Ubf(i,k)
		   Uvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   Uvel2(i,-2,k) = Ubf(i,k)
		   Uvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	if (periodicx.eq.0.or.periodicx.eq.2) then
		Uvel2(-2,0:j1,0:k1)=Uvel(0,0:j1,0:k1)
		Uvel2(-1,0:j1,0:k1)=Uvel(0,0:j1,0:k1)
		Uvel2(i1+1,0:j1,0:k1)=Uvel(ie,0:j1,0:k1)
		Uvel2(i1+2,0:j1,0:k1)=Uvel(ie,0:j1,0:k1)
	else 
		Uvel2(-2,0:j1,0:k1)=Uvel(ie-2,0:j1,0:k1)
		Uvel2(-1,0:j1,0:k1)=Uvel(ie-1,0:j1,0:k1)
		Uvel2(i1+1,0:j1,0:k1)=Uvel(2,0:j1,0:k1)
		Uvel2(i1+2,0:j1,0:k1)=Uvel(3,0:j1,0:k1)
	endif

	Uvel2(0:i1,0:j1,-2)=Uvel(0:i1,0:j1,0)
	Uvel2(0:i1,0:j1,-1)=Uvel(0:i1,0:j1,0)
	Uvel2(0:i1,0:j1,k1+1)=Uvel(0:i1,0:j1,ke)
	Uvel2(0:i1,0:j1,k1+2)=Uvel(0:i1,0:j1,ke)


	if (wd.eq.1) then
		Vvel2(0:i1,0:j1,0:k1)=Vvel
	!c get stuff from other CPU's
		  call shiftf2(Vvel,ubf)
		  call shiftb2(Vvel,ubb) 

		if (periodicy.eq.0.or.periodicy.eq.2) then
		  if (rank.eq.0) then
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Vvel2(i,j1+1,k) =Ubb(i,k)
			   Vvel2(i,-1,k) = Vvel(i,0,k)
			   enddo
			enddo
		  elseif (rank.eq.px-1) then
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Vvel2(i,j1,k) =Vvel(i,je,k)
			   Vvel2(i,j1+1,k) =Vvel(i,je,k)
			   Vvel2(i,-1,k) = Ubf(i,k)
			   enddo
			enddo
		  else 
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Vvel2(i,-1,k) = Ubf(i,k)
			   Vvel2(i,j1+1,k) =Ubb(i,k)
			   enddo
			enddo
		  endif
		else
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Vvel2(i,-1,k) = Ubf(i,k)
			   Vvel2(i,j1+1,k) =Ubb(i,k)
			   enddo
			enddo
		endif
		  call shiftf3(Vvel,ubf)
		  call shiftb3(Vvel,ubb) 

		if (periodicy.eq.0.or.periodicy.eq.2) then
		  if (rank.eq.0) then
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Vvel2(i,-2,k) = Vvel(i,0,k)
			   Vvel2(i,j1+2,k) =Ubb(i,k)
			   enddo
			enddo
		  elseif (rank.eq.px-1) then
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Vvel2(i,-2,k) = Ubf(i,k)
			   Vvel2(i,j1+2,k) =Vvel(i,je,k)
			   enddo
			enddo
		  else 
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Vvel2(i,-2,k) = Ubf(i,k)
			   Vvel2(i,j1+2,k) =Ubb(i,k)
			   enddo
			enddo
		  endif
		else 
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Vvel2(i,-2,k) = Ubf(i,k)
			   Vvel2(i,j1+2,k) =Ubb(i,k)
			   enddo
			enddo
		endif

		if (periodicx.eq.0.or.periodicx.eq.2) then
			Vvel2(-1,-2:j1+2,0:k1)=  Vvel2(0,-2:j1+2,0:k1)
			Vvel2(-2,-2:j1+2,0:k1)=  Vvel2(0,-2:j1+2,0:k1)
			Vvel2(i1+1,-2:j1+2,0:k1)=Vvel2(i1,-2:j1+2,0:k1)
			Vvel2(i1+2,-2:j1+2,0:k1)=Vvel2(i1,-2:j1+2,0:k1)
		else
			Vvel2(-1,-2:j1+2,0:k1)=  Vvel2(ie-1,-2:j1+2,0:k1)
			Vvel2(-2,-2:j1+2,0:k1)=  Vvel2(ie-2,-2:j1+2,0:k1)
			Vvel2(i1+1,-2:j1+2,0:k1)=Vvel2(2,-2:j1+2,0:k1)
			Vvel2(i1+2,-2:j1+2,0:k1)=Vvel2(3,-2:j1+2,0:k1)
		endif
	!	!22-12-2017: test bc bed Uvel=0 (if it works also at immersed boundary all is fine! and with U(i,j,-1)=U(i,j,0) wiggly time-avg results for periodic turb channel MKM near bed.
		Vvel2(-2:i1+2,-2:j1+2,-1)=  Vvel2(-2:i1+2,-2:j1+2,0)
		Vvel2(-2:i1+2,-2:j1+2,-2)=  Vvel2(-2:i1+2,-2:j1+2,0)
		Vvel2(-2:i1+2,-2:j1+2,k1+1)=Vvel2(-2:i1+2,-2:j1+2,k1)
		Vvel2(-2:i1+2,-2:j1+2,k1+2)=Vvel2(-2:i1+2,-2:j1+2,k1)

			Wvel2(0:i1,0:j1,0:k1)=Wvel

	!c get stuff from other CPU's
		  call shiftf2(Wvel,ubf)
		  call shiftb2(Wvel,ubb) 
		if (periodicy.eq.0.or.periodicy.eq.2) then
		  if (rank.eq.0) then
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Wvel2(i,j1+1,k) =Ubb(i,k)
			   Wvel2(i,-1,k) = Wvel(i,0,k)
			   enddo
			enddo
		  elseif (rank.eq.px-1) then
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Wvel2(i,j1+1,k) =Wvel(i,j1,k)
			   Wvel2(i,-1,k) = Ubf(i,k)
			   enddo
			enddo
		  else 
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Wvel2(i,-1,k) = Ubf(i,k)
			   Wvel2(i,j1+1,k) =Ubb(i,k)
			   enddo
			enddo
		  endif
		else 
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Wvel2(i,-1,k) = Ubf(i,k)
			   Wvel2(i,j1+1,k) =Ubb(i,k)
			   enddo
			enddo
		endif
		  call shiftf3(Wvel,ubf)
		  call shiftb3(Wvel,ubb) 

		if (periodicy.eq.0.or.periodicy.eq.2) then
		  if (rank.eq.0) then
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Wvel2(i,j1+2,k) =Ubb(i,k)
			   Wvel2(i,-2,k) = Wvel(i,0,k)
			   enddo
			enddo
		  elseif (rank.eq.px-1) then
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Wvel2(i,j1+2,k) =Wvel(i,j1,k)
			   Wvel2(i,-2,k) = Ubf(i,k)
			   enddo
			enddo
		  else 
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Wvel2(i,-2,k) = Ubf(i,k)
			   Wvel2(i,j1+2,k) =Ubb(i,k)
			   enddo
			enddo
		  endif
		else 
			do k=0,k1
			   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
			   Wvel2(i,-2,k) = Ubf(i,k)
			   Wvel2(i,j1+2,k) =Ubb(i,k)
			   enddo
			enddo
		endif

		if (periodicx.eq.0.or.periodicx.eq.2) then
			Wvel2(-1,-2:j1+2,0:k1)=  Wvel2(0,-2:j1+2,0:k1)
			Wvel2(-2,-2:j1+2,0:k1)=  Wvel2(0,-2:j1+2,0:k1)
			Wvel2(i1+1,-2:j1+2,0:k1)=Wvel2(i1,-2:j1+2,0:k1)
			Wvel2(i1+2,-2:j1+2,0:k1)=Wvel2(i1,-2:j1+2,0:k1)
		else
			Wvel2(-1,-2:j1+2,0:k1)=  Wvel2(ie-1,-2:j1+2,0:k1)
			Wvel2(-2,-2:j1+2,0:k1)=  Wvel2(ie-2,-2:j1+2,0:k1)
			Wvel2(i1+1,-2:j1+2,0:k1)=Wvel2(2,-2:j1+2,0:k1)
			Wvel2(i1+2,-2:j1+2,0:k1)=Wvel2(3,-2:j1+2,0:k1)
		endif
		Wvel2(-2:i1+2,-2:j1+2,-1)=  Wvel2(-2:i1+2,-2:j1+2,0)
		Wvel2(-2:i1+2,-2:j1+2,-2)=  Wvel2(-2:i1+2,-2:j1+2,0)
		Wvel2(-2:i1+2,-2:j1+2,k1+1)=Wvel2(-2:i1+2,-2:j1+2,k1)
		Wvel2(-2:i1+2,-2:j1+2,k1+2)=Wvel2(-2:i1+2,-2:j1+2,k1)

	
      do k=1,k1
        kp=k+1
        km=k-1
		kpp=k+2 !MIN(k+2,k1)
		kmm=k-2 !MAX(k-2,0)
        do j=1,j1
          jp=j+1
          jm=j-1
		  jpp=j+2
		  jmm=j-2
          do  i=1,i1
            ip=i+1
            im=i-1
		    ipp=i+2
		    imm=i-2
			!Wiggle detector Sagaut 2001/PhD Simon 2003:
			!wwdd only on own directions of U,V,W gave best results in combination with grow/decrease factor 0.001
			wwdd(1)=MIN((Uvel2(ip,j ,k )-Uvel2(i,j,k))*(Uvel2(i,j,k)-Uvel2(im,j ,k )),0.)
     & 				*MIN((Uvel2(i,j,k)-Uvel2(im,j ,k ))*(Uvel2(im,j ,k )-Uvel2(imm,j  ,k  )),0.)
			wwdd(2)=MIN((Vvel2(i ,jp,k )-Vvel2(i,j,k))*(Vvel2(i,j,k)-Vvel2(i ,jm,k )),0.)
     &				*MIN((Vvel2(i,j,k)-Vvel2(i ,jm,k ))*(Vvel2(i ,jm,k )-Vvel2(i  ,jmm,k  )),0.)
			wwdd(3)=MIN((Wvel2(i ,j ,kp)-Wvel2(i,j,k))*(Wvel2(i,j,k)-Wvel2(i ,j ,km)),0.)
     &				*MIN((Wvel2(i,j,k)-Wvel2(i ,j ,km))*(Wvel2(i ,j ,km)-Wvel2(i  ,j  ,kmm)),0.)	
	 
			u4mag=((0.5*(Uvel2(i,j,k)+Uvel2(im,j,k)))**2+(0.5*(Vvel2(i,j,k)+Vvel2(i,jm,k)))**2+
     &			(0.5*(Wvel2(i,j,k)+Wvel2(i,j,km)))**2)
			u4mag=u4mag+1e-12
			wwdd=wwdd/u4mag
			IF (MAXVAL(wwdd).gt.1.e-12) THEN
				wf(i,j,k)=MIN(1.,wf(i,j,k)+0.001) !wf(i,j,k)=1. 
			ELSE
				wf(i,j,k)=MAX(0.,wf(i,j,k)-0.001) !wf(i,j,k)=0.
			ENDIF
		  enddo
	    enddo
	  enddo
	endif ! wf=1. is already default when wiggle_detector==0

	
		  Ax=9.   !CDS4-Wicker: 7 !CDS4-Morinishi: 9
		  Bx=-1.  !CDS4-Wicker:-1 !CDS4-Morinishi:-1
		  Cx=0.
		  DDx=1./(Ax+Bx+Cx)
		  Ay=9.
		  By=-1.
		  Cy=0.
		  DDy=1./(Ay+By+Cy)
		  Az=9.
		  Bz=-1.
		  Cz=0.
		  DDz=1./(Az+Bz+Cz)
	dzi=1./dz 
      do k=kb,ke
      kp=k+1
      km=k-1
		  kpp=k+2
		  kppp=k+3
		  kmm=k-2
		  kmmm=k-3
        do j=jb,je
        jp=j+1
        jm=j-1
		  jpp=j+2
		  jppp=j+3
		  jmm=j-2
		  jmmm=j-3
          do  i=ib,ie
          ip=i+1
          im=i-1
		  ipp=i+2
		  ippp=i+3
		  imm=i-2
		  immm=i-3
          !!rhojp =0.25*(rho(i,j,k)+rho(i,jp,k)+rho(ip,j,k)+rho(ip,jp,k))
          !!rhojm =0.25*(rho(i,j,k)+rho(i,jm,k)+rho(ip,j,k)+rho(ip,jm,k))
          !!rhokp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
          !!rhokm =0.25*(rho(i,j,k)+rho(i,j,km)+rho(ip,j,k)+rho(ip,j,km))
      
      !!uuR=(Uvel(i,j,k)+Uvel(ip,j,k))*rho(ip,j,k)*Rp(ip)
      !!uuL=(Uvel(i,j,k)+Uvel(im,j,k))*rho(i,j,k)*Rp(i)
      !!vvR=(Vvel(i,j,k)+Vvel(ip,j,k))*rhojp
      !!vvL=(Vvel(i,jm,k)+Vvel(ip,jm,k))*rhojm
      !!wwR=(Wvel(i,j,k)+Wvel(ip,j,k))*rhokp
      !!wwL=(Wvel(i,j,km)+Wvel(ip,j,km))*rhokm
	  
            uuR=(Uvel(i,j,k)*rhU(i,j,k)+Uvel(ip,j,k)*rhU(ip,j,k))*Rp(ip)
            uuL=(Uvel(i,j,k)*rhU(i,j,k)+Uvel(im,j,k)*rhU(im,j,k))*Rp(i)
            vvR=(Vvel(i,j,k)*rhV(i,j,k)+Vvel(ip,j,k)*rhV(ip,j,k))
            vvL=(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(ip,jm,k)*rhV(ip,jm,k))
            wwR=(Wvel(i,j,k)*rhW(i,j,k)+Wvel(ip,j,k)*rhW(ip,j,k))
            wwL=(Wvel(i,j,km)*rhW(i,j,km)+Wvel(ip,j,km)*rhW(ip,j,km))	  
	  
		numdif2=numdif*(wf(i,j,k)+wf(i+1,j,k))*0.5
		
		putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (
     1 (uuR    * DDx*(Ax*(Uvel2(i,j,k)+Uvel2(ip,j,k))+Bx*(Uvel2(im,j,k)+Uvel2(ipp,j,k))+Cx*(Uvel2(imm,j,k)+Uvel2(ippp,j,k)))-
     1 numdif2*
     1 ABS(uuR)* Cx*(10.*(-Uvel2(i,j,k)+Uvel2(ip,j,k))-5.*(-Uvel2(im,j,k)+Uvel2(ipp,j,k))+(-Uvel2(imm,j,k)+Uvel2(ippp,j,k))) -
     1 (uuL    * DDx*(Ax*(Uvel2(im,j,k)+Uvel2(i,j,k))+Bx*(Uvel2(imm,j,k)+Uvel2(ip,j,k))+Cx*(Uvel2(immm,j,k)+Uvel2(ipp,j,k)))-
     1 numdif2*
     1 ABS(uuL)* Cx*(10.*(-Uvel2(im,j,k)+Uvel2(i,j,k))-5.*(-Uvel2(imm,j,k)+Uvel2(ip,j,k))+(-Uvel2(immm,j,k)+Uvel2(ipp,j,k)))) )
     1  / ( Ru(i) * ( Rp(ip)-Rp(i) ) )
     +                         + 
     2 (vvR    * DDy*(Ay*(Uvel2(i,j,k)+Uvel2(i,jp,k))+By*(Uvel2(i,jm,k)+Uvel2(i,jpp,k))+Cy*(Uvel2(i,jmm,k)+Uvel2(i,jppp,k)))-
     2 numdif2*
     2 ABS(vvR)* Cy*(10.*(-Uvel2(i,j,k)+Uvel2(i,jp,k))-5.*(-Uvel2(i,jm,k)+Uvel2(i,jpp,k))+(-Uvel2(i,jmm,k)+Uvel2(i,jppp,k))) -
     2 (vvL    * DDy*(Ay*(Uvel2(i,jm,k)+Uvel2(i,j,k))+By*(Uvel2(i,jmm,k)+Uvel2(i,jp,k))+Cy*(Uvel2(i,jmmm,k)+Uvel2(i,jpp,k)))-
     2 numdif2*
     2 ABS(vvL)* Cy*(10.*(-Uvel2(i,jm,k)+Uvel2(i,j,k))-5.*(-Uvel2(i,jmm,k)+Uvel2(i,jp,k))+(-Uvel2(i,jmmm,k)+Uvel2(i,jpp,k)))) )
     2  / ( Ru(i) * (phiv(j)-phiv(jm)) )
     +                         +
     3 (wwR    * DDz*(Az*(Uvel2(i,j,k)+Uvel2(i,j,kp))+Bz*(Uvel2(i,j,km)+Uvel2(i,j,kpp))+Cz*(Uvel2(i,j,kmm)+Uvel2(i,j,kppp)))-
     3 numdif2*
     3 ABS(wwR)* Cz*(10.*(-Uvel2(i,j,k)+Uvel2(i,j,kp))-5.*(-Uvel2(i,j,km)+Uvel2(i,j,kpp))+(-Uvel2(i,j,kmm)+Uvel2(i,j,kppp))) -
     3 (wwL    * DDz*(Az*(Uvel2(i,j,km)+Uvel2(i,j,k))+Bz*(Uvel2(i,j,kmm)+Uvel2(i,j,kp))+Cz*(Uvel2(i,j,kmmm)+Uvel2(i,j,kpp)))-
     3 numdif2*
     3 ABS(wwL)* Cz*(10.*(-Uvel2(i,j,km)+Uvel2(i,j,k))-5.*(-Uvel2(i,j,kmm)+Uvel2(i,j,kp))+(-Uvel2(i,j,kmmm)+Uvel2(i,j,kpp)))) )
     3  *dzi
     +                         -
     4 0.5*(rho(i,j,k)+rho(ip,j,k))*( Vvel(i,j,k) + Vvel(ip,j,k) + Vvel(i,jm,k) + Vvel(ip,jm,k) )**2
     4  / ( 4.0 * Ru(i) )
     +                         )
           enddo
        enddo
      enddo
      return
      end


      subroutine advecv_CDS4(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,wf,wd)
      implicit none
c
c********************************************************************
c
c     advecv calculates the advection of the v-velocity, which is
c     the velocity in the tangential direction.
c
c     In formula:
c
c         1 d(ruv)     1 d(vv)     d(vw)     uv
c    - (  - ------  +  - -----  +  -----  +  --  )
c         r   dr       r  dphi      dz        r
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          Vtmp              : contains velocity at oldest timestep
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1),dzi
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm
	real Ax,Bx,Cx,DDx
	real Ay,By,Cy,DDy,numdif
	real Az,Bz,Cz,DDz
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy,wd

	real*8 ubb(0:i1,0:k1),ubf(0:i1,0:k1),Vvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real numdif2,wf(0:i1,0:j1,0:k1)
	real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1)

	Vvel2(0:i1,0:j1,0:k1)=Vvel

!c get stuff from other CPU's
	  call shiftf2(Vvel,ubf)
	  call shiftb2(Vvel,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   Vvel2(i,j1+1,k) =Ubb(i,k)
		   Vvel2(i,-1,k) = Vvel(i,0,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   Vvel2(i,j1,k) =Vvel(i,je,k)
		   Vvel2(i,j1+1,k) =Vvel(i,je,k)
		   Vvel2(i,-1,k) = Ubf(i,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   Vvel2(i,-1,k) = Ubf(i,k)
		   Vvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else
		do k=1,ke
		   do i=1,ie
		   Vvel2(i,-1,k) = Ubf(i,k)
		   Vvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	  call shiftf3(Vvel,ubf)
	  call shiftb3(Vvel,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   Vvel2(i,-2,k) = Vvel(i,0,k)
		   Vvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   Vvel2(i,-2,k) = Ubf(i,k)
		   Vvel2(i,j1+2,k) =Vvel(i,je,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   Vvel2(i,-2,k) = Ubf(i,k)
		   Vvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   Vvel2(i,-2,k) = Ubf(i,k)
		   Vvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	endif

	if (periodicx.eq.0.or.periodicx.eq.2) then
		Vvel2(-2,0:j1,0:k1)=Vvel(0,0:j1,0:k1)
		Vvel2(-1,0:j1,0:k1)=Vvel(0,0:j1,0:k1)
		Vvel2(i1+1,0:j1,0:k1)=Vvel(ie,0:j1,0:k1)
		Vvel2(i1+2,0:j1,0:k1)=Vvel(ie,0:j1,0:k1)
	else 
		Vvel2(-2,0:j1,0:k1)=Vvel(ie-2,0:j1,0:k1)
		Vvel2(-1,0:j1,0:k1)=Vvel(ie-1,0:j1,0:k1)
		Vvel2(i1+1,0:j1,0:k1)=Vvel(2,0:j1,0:k1)
		Vvel2(i1+2,0:j1,0:k1)=Vvel(3,0:j1,0:k1)
	endif

	Vvel2(0:i1,0:j1,-2)=Vvel(0:i1,0:j1,0)
	Vvel2(0:i1,0:j1,-1)=Vvel(0:i1,0:j1,0)
	Vvel2(0:i1,0:j1,k1+1)=Vvel(0:i1,0:j1,ke)
	Vvel2(0:i1,0:j1,k1+2)=Vvel(0:i1,0:j1,ke)

		  Ax=9.   !CDS4-Wicker: 7 !CDS4-Morinishi: 9
		  Bx=-1.  !CDS4-Wicker:-1 !CDS4-Morinishi:-1
		  Cx=0.
		  DDx=1./(Ax+Bx+Cx)
		  Ay=9.
		  By=-1.
		  Cy=0.
		  DDy=1./(Ay+By+Cy)
		  Az=9.
		  Bz=-1.
		  Cz=0.
		  DDz=1./(Az+Bz+Cz)
	dzi=1./dz 
      do k=kb,ke
      kp=k+1
      km=k-1
		  kpp=k+2
		  kppp=k+3
		  kmm=k-2
		  kmmm=k-3
        do j=jb,je
        jp=j+1
        jm=j-1
		  jpp=j+2
		  jppp=j+3
		  jmm=j-2
		  jmmm=j-3
          do  i=ib,ie
          ip=i+1
          im=i-1
		  ipp=i+2
		  ippp=i+3
		  imm=i-2
		  immm=i-3
          !rhoip =0.25*(rho(i,j,k)+rho(ip,j,k)+rho(i,jp,k)+rho(ip,jp,k))
          !rhoim =0.25*(rho(i,j,k)+rho(im,j,k)+rho(i,jp,k)+rho(im,jp,k))
          !rhojp =rho(i,jp,k)
          !rhojm =rho(i,j,k )
          !rhokp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
          !rhokm =0.25*(rho(i,j,k)+rho(i,j,km)+rho(i,jp,k)+rho(i,jp,km))
      !uuR=(Uvel(i,j,k)+Uvel(i,jp,k))*rhoip*Ru(i) !*Ru(i)
      !uuL=(Uvel(im,j,k)+Uvel(im,jp,k))*rhoim*Ru(im) !*Ru(im)
      !vvR=(Vvel(i,j,k)+Vvel(i,jp,k))*rhojp
      !vvL=(Vvel(i,jm,k)+Vvel(i,j,k))*rhojm
      !wwR=(Wvel(i,j,k)+Wvel(i,jp,k))*rhokp
      !wwL=(Wvel(i,j,km)+Wvel(i,jp,km))*rhokm
           uuR=(Uvel(i,j,k)*rhU(i,j,k)+Uvel(i,jp,k)*rhU(i,jp,k))*Ru(i) !*Ru(i)
           uuL=(Uvel(im,j,k)*rhU(im,j,k)+Uvel(im,jp,k)*rhU(im,jp,k))*Ru(im) !*Ru(im)
     	    vvR=(Vvel(i,j,k)*rhV(i,j,k)+Vvel(i,jp,k)*rhV(i,jp,k))
           vvL=(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(i,j,k)*rhV(i,j,k))
           wwR=(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,jp,k)*rhW(i,jp,k))
           wwL=(Wvel(i,j,km)*rhW(i,j,km)+Wvel(i,jp,km)*rhW(i,jp,km))	  

		numdif2=numdif*(wf(i,j,k)+wf(i,j+1,k))*0.5
		
      putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (

     1 (uuR    * DDx*(Ax*(Vvel2(i,j,k)+Vvel2(ip,j,k))+Bx*(Vvel2(im,j,k)+Vvel2(ipp,j,k))+Cx*(Vvel2(imm,j,k)+Vvel2(ippp,j,k)))-
     1 numdif2*
     1 ABS(uuR)* Cx*(10.*(-Vvel2(i,j,k)+Vvel2(ip,j,k))-5.*(-Vvel2(im,j,k)+Vvel2(ipp,j,k))+(-Vvel2(imm,j,k)+Vvel2(ippp,j,k))) -
     1 (uuL    * DDx*(Ax*(Vvel2(im,j,k)+Vvel2(i,j,k))+Bx*(Vvel2(imm,j,k)+Vvel2(ip,j,k))+Cx*(Vvel2(immm,j,k)+Vvel2(ipp,j,k)))-
     1 numdif2*
     1 ABS(uuL)* Cx*(10.*(-Vvel2(im,j,k)+Vvel2(i,j,k))-5.*(-Vvel2(imm,j,k)+Vvel2(ip,j,k))+(-Vvel2(immm,j,k)+Vvel2(ipp,j,k)))) )
     1  / ( Rp(i)* dr(i) )   !/ ( Rp(i) * Rp(i)* dr(i) ) 
     +                         + 
     2 (vvR    * DDy*(Ay*(Vvel2(i,j,k)+Vvel2(i,jp,k))+By*(Vvel2(i,jm,k)+Vvel2(i,jpp,k))+Cy*(Vvel2(i,jmm,k)+Vvel2(i,jppp,k)))-
     2 numdif2*
     2 ABS(vvR)* Cy*(10.*(-Vvel2(i,j,k)+Vvel2(i,jp,k))-5.*(-Vvel2(i,jm,k)+Vvel2(i,jpp,k))+(-Vvel2(i,jmm,k)+Vvel2(i,jppp,k))) -
     2 (vvL    * DDy*(Ay*(Vvel2(i,jm,k)+Vvel2(i,j,k))+By*(Vvel2(i,jmm,k)+Vvel2(i,jp,k))+Cy*(Vvel2(i,jmmm,k)+Vvel2(i,jpp,k)))-
     2 numdif2*
     2 ABS(vvL)* Cy*(10.*(-Vvel2(i,jm,k)+Vvel2(i,j,k))-5.*(-Vvel2(i,jmm,k)+Vvel2(i,jp,k))+(-Vvel2(i,jmmm,k)+Vvel2(i,jpp,k)))) )
     2  / ( Rp(i) * (phiv(j)-phiv(jm)) )
     +                         +
     3 (wwR    * DDz*(Az*(Vvel2(i,j,k)+Vvel2(i,j,kp))+Bz*(Vvel2(i,j,km)+Vvel2(i,j,kpp))+Cz*(Vvel2(i,j,kmm)+Vvel2(i,j,kppp)))-
     3 numdif2*
     3 ABS(wwR)* Cz*(10.*(-Vvel2(i,j,k)+Vvel2(i,j,kp))-5.*(-Vvel2(i,j,km)+Vvel2(i,j,kpp))+(-Vvel2(i,j,kmm)+Vvel2(i,j,kppp))) -
     3 (wwL    * DDz*(Az*(Vvel2(i,j,km)+Vvel2(i,j,k))+Bz*(Vvel2(i,j,kmm)+Vvel2(i,j,kp))+Cz*(Vvel2(i,j,kmmm)+Vvel2(i,j,kpp)))-
     3 numdif2*
     3 ABS(wwL)* Cz*(10.*(-Vvel2(i,j,km)+Vvel2(i,j,k))-5.*(-Vvel2(i,j,kmm)+Vvel2(i,j,kp))+(-Vvel2(i,j,kmmm)+Vvel2(i,j,kpp)))) )
     3  *dzi

     +                         +
!     4   0.5*(rho(i,j,k)+rho(i,jp,k))*
!     4   ((Uvel(im,j,k)+Uvel(i,j,k))*(Vvel(i,jm,k)+Vvel(i,j,k))+(Uvel(im,jp,k)+Uvel(i,jp,k))*(Vvel(i,jp,k)+Vvel(i,j,k)))
!     4   / ( 2.0*Rp(i) )

     4   0.5*(rho(i,j,k)+rho(i,jp,k))*
     4   (Uvel(im,j,k)+Uvel(i,j,k)+Uvel(im,jp,k)+Uvel(i,jp,k))*Vvel(i,j,k)
     4   / ( Rp(i) )
     +                         )
           enddo
        enddo
      enddo
      return
      end

      subroutine advecw_CDS4(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,wf,wd)
      implicit none
c
c********************************************************************
c
c     advecw calculates the advection of the w-velocity, which is
c     the velocity in the axial direction.
c
c     In formula:
c
c         1 d(ruw)     1 d(wv)     d(ww)
c    - (  - ------  +  - -----  +  -----  )
c         r   dr       r  dphi      dz
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          Wtmp              : contains velocity at oldest timestep
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1),dzi
      real rhoip,rhoim,rhojp,rhojm
	real Ax,Bx,Cx,DDx
	real Ay,By,Cy,DDy
	real Az,Bz,Cz,DDz
	real uuR,uuL,vvR,vvL,wwR,wwL,numdif
	integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy,wd
	real*8 ubb(0:i1,0:k1),ubf(0:i1,0:k1),Wvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real numdif2,wf(0:i1,0:j1,0:k1)	
	real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1)

	Wvel2(0:i1,0:j1,0:k1)=Wvel

!c get stuff from other CPU's
	  call shiftf2(Wvel,ubf)
	  call shiftb2(Wvel,ubb) 
	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   Wvel2(i,j1+1,k) =Ubb(i,k)
		   Wvel2(i,-1,k) = Wvel(i,0,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   Wvel2(i,j1+1,k) =Wvel(i,j1,k)
		   Wvel2(i,-1,k) = Ubf(i,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   Wvel2(i,-1,k) = Ubf(i,k)
		   Wvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   Wvel2(i,-1,k) = Ubf(i,k)
		   Wvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	  call shiftf3(Wvel,ubf)
	  call shiftb3(Wvel,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   Wvel2(i,j1+2,k) =Ubb(i,k)
		   Wvel2(i,-2,k) = Wvel(i,0,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   Wvel2(i,j1+2,k) =Wvel(i,j1,k)
		   Wvel2(i,-2,k) = Ubf(i,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   Wvel2(i,-2,k) = Ubf(i,k)
		   Wvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   Wvel2(i,-2,k) = Ubf(i,k)
		   Wvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	endif

	if (periodicx.eq.0.or.periodicx.eq.2) then
		Wvel2(-2,0:j1,0:k1)=Wvel(0,0:j1,0:k1)
		Wvel2(-1,0:j1,0:k1)=Wvel(0,0:j1,0:k1)
		Wvel2(i1+1,0:j1,0:k1)=Wvel(ie,0:j1,0:k1)
		Wvel2(i1+2,0:j1,0:k1)=Wvel(ie,0:j1,0:k1)
	else 
		Wvel2(-2,0:j1,0:k1)=Wvel(ie-2,0:j1,0:k1)
		Wvel2(-1,0:j1,0:k1)=Wvel(ie-1,0:j1,0:k1)
		Wvel2(i1+1,0:j1,0:k1)=Wvel(2,0:j1,0:k1)
		Wvel2(i1+2,0:j1,0:k1)=Wvel(3,0:j1,0:k1)
	endif

	Wvel2(0:i1,0:j1,-2)=Wvel(0:i1,0:j1,0)
	Wvel2(0:i1,0:j1,-1)=Wvel(0:i1,0:j1,0)
	Wvel2(0:i1,0:j1,k1)=Wvel(0:i1,0:j1,ke)
	Wvel2(0:i1,0:j1,k1+1)=Wvel(0:i1,0:j1,ke)
	Wvel2(0:i1,0:j1,k1+2)=Wvel(0:i1,0:j1,ke)


	
	
	dzi=1./dz 


		  Ax=9.   !CDS4-Wicker: 7 !CDS4-Morinishi: 9
		  Bx=-1.  !CDS4-Wicker:-1 !CDS4-Morinishi:-1
		  Cx=0.
		  DDx=1./(Ax+Bx+Cx)
		  Ay=9.
		  By=-1.
		  Cy=0.
		  DDy=1./(Ay+By+Cy)
		  Az=9.
		  Bz=-1.
		  Cz=0.
		  DDz=1./(Az+Bz+Cz)
      do k=kb,ke
      kp=k+1
      km=k-1
		  kpp=k+2
		  kppp=k+3
		  kmm=k-2
		  kmmm=k-3
        do j=jb,je
        jp=j+1
        jm=j-1
		  jpp=j+2
		  jppp=j+3
		  jmm=j-2
		  jmmm=j-3
          do  i=ib,ie
          ip=i+1
          im=i-1
		  ipp=i+2
		  ippp=i+3
		  imm=i-2
		  immm=i-3
!!          rhoip =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
!!          rhoim =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(im,j,k)+rho(im,j,kp))
!!          rhojp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
!!          rhojm =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jm,k)+rho(i,jm,kp))
!!      uuR=(Uvel(i,j,k)+Uvel(i,j,kp))*rhoip*Ru(i)
!!      uuL=(Uvel(im,j,k)+Uvel(im,j,kp))*rhoim*Ru(im)
!!      vvR=(Vvel(i,j,k)+Vvel(i,j,kp))*rhojp
!!      vvL=(Vvel(i,jm,k)+Vvel(i,jm,kp))*rhojm
!!      wwR=(Wvel(i,j,k)+Wvel(i,j,kp))*rho(i,j,kp)
!!      wwL=(Wvel(i,j,k)+Wvel(i,j,km))*rho(i,j,k)

      	    uuR=(Uvel(i,j,k)*rhU(i,j,k)+Uvel(i,j,kp)*rhU(i,j,kp))*Ru(i)
      	    uuL=(Uvel(im,j,k)*rhU(im,j,k)+Uvel(im,j,kp)*rhU(im,j,kp))*Ru(im)
			vvR=(Vvel(i,j,k)*rhV(i,j,k)+Vvel(i,j,kp)*rhV(i,j,kp))
    	    vvL=(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(i,jm,kp)*rhV(i,jm,kp))
    	    wwR=(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,j,kp)*rhW(i,j,kp))
			wwL=(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,j,km)*rhW(i,j,km))  
			

	  numdif2=numdif*(wf(i,j,k)+wf(i,j,k+1))*0.5
	  
      putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (
     1 (uuR    * DDx*(Ax*(Wvel2(i,j,k)+Wvel2(ip,j,k))+Bx*(Wvel2(im,j,k)+Wvel2(ipp,j,k))+Cx*(Wvel2(imm,j,k)+Wvel2(ippp,j,k)))-
     1 numdif2*
     1 ABS(uuR)* Cx*(10.*(-Wvel2(i,j,k)+Wvel2(ip,j,k))-5.*(-Wvel2(im,j,k)+Wvel2(ipp,j,k))+(-Wvel2(imm,j,k)+Wvel2(ippp,j,k))) -
     1 (uuL    * DDx*(Ax*(Wvel2(im,j,k)+Wvel2(i,j,k))+Bx*(Wvel2(imm,j,k)+Wvel2(ip,j,k))+Cx*(Wvel2(immm,j,k)+Wvel2(ipp,j,k)))-
     1 numdif2*
     1 ABS(uuL)* Cx*(10.*(-Wvel2(im,j,k)+Wvel2(i,j,k))-5.*(-Wvel2(imm,j,k)+Wvel2(ip,j,k))+(-Wvel2(immm,j,k)+Wvel2(ipp,j,k)))) )
     1  / ( Rp(i)* dr(i) ) 
     +                         + 
     2 (vvR    * DDy*(Ay*(Wvel2(i,j,k)+Wvel2(i,jp,k))+By*(Wvel2(i,jm,k)+Wvel2(i,jpp,k))+Cy*(Wvel2(i,jmm,k)+Wvel2(i,jppp,k)))-
     2 numdif2*
     2 ABS(vvR)* Cy*(10.*(-Wvel2(i,j,k)+Wvel2(i,jp,k))-5.*(-Wvel2(i,jm,k)+Wvel2(i,jpp,k))+(-Wvel2(i,jmm,k)+Wvel2(i,jppp,k))) -
     2 (vvL    * DDy*(Ay*(Wvel2(i,jm,k)+Wvel2(i,j,k))+By*(Wvel2(i,jmm,k)+Wvel2(i,jp,k))+Cy*(Wvel2(i,jmmm,k)+Wvel2(i,jpp,k)))-
     2 numdif2*
     2 ABS(vvL)* Cy*(10.*(-Wvel2(i,jm,k)+Wvel2(i,j,k))-5.*(-Wvel2(i,jmm,k)+Wvel2(i,jp,k))+(-Wvel2(i,jmmm,k)+Wvel2(i,jpp,k)))) )
     2  / ( Rp(i) * (phiv(j)-phiv(jm)) )
     +                         +
     3 (wwR    * DDz*(Az*(Wvel2(i,j,k)+Wvel2(i,j,kp))+Bz*(Wvel2(i,j,km)+Wvel2(i,j,kpp))+Cz*(Wvel2(i,j,kmm)+Wvel2(i,j,kppp)))-
     3 numdif2*
     3 ABS(wwR)* Cz*(10.*(-Wvel2(i,j,k)+Wvel2(i,j,kp))-5.*(-Wvel2(i,j,km)+Wvel2(i,j,kpp))+(-Wvel2(i,j,kmm)+Wvel2(i,j,kppp))) -
     3 (wwL    * DDz*(Az*(Wvel2(i,j,km)+Wvel2(i,j,k))+Bz*(Wvel2(i,j,kmm)+Wvel2(i,j,kp))+Cz*(Wvel2(i,j,kmmm)+Wvel2(i,j,kpp)))-
     3 numdif2*
     3 ABS(wwL)* Cz*(10.*(-Wvel2(i,j,km)+Wvel2(i,j,k))-5.*(-Wvel2(i,j,kmm)+Wvel2(i,j,kp))+(-Wvel2(i,j,kmmm)+Wvel2(i,j,kpp)))) )
     3  *dzi
     +                         )
           enddo
         enddo
      enddo
      return
      end


