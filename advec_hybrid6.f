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


      subroutine advecu_HYB6(putout,Uvel,Vvel,Wvel,RHO,rhu,rhv,rhw,Ru,Rp,dr,phiv,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,kbed,wf,wd,nd2,fcg)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,kbed(0:i1,0:j1)
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1),numdif,nd2
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy,wd,jmax 
	real*8 ubb(0:i1,0:k1),ubf(0:i1,0:k1),Uvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1)
	real wwdd(1:3),wf(0:i1,0:j1,0:k1),numdif2,u4mag 
	real Vvel2(-2:i1+2,-2:j1+2,-2:k1+2),Wvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real fcg(0:i1,0:px*(j1-1)+1,0:k1),fcg_local
	
	Uvel2(0:i1,0:j1,0:k1)=Uvel
	jmax=j1-1 



!c get stuff from other CPU's
	  call shiftf2(Uvel,ubf)
	  call shiftb2(Uvel,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly 
		   Uvel2(i,-1,k) = Uvel(i,0,k)
		   Uvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   Uvel2(i,-1,k) = Ubf(i,k)
		   Uvel2(i,j1+1,k) =Uvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   Uvel2(i,-1,k) = Ubf(i,k)
		   Uvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   Uvel2(i,-1,k) = Ubf(i,k)
		   Uvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	  call shiftf3(Uvel,ubf)
	  call shiftb3(Uvel,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   Uvel2(i,-2,k) = Uvel(i,0,k)
		   Uvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   Uvel2(i,-2,k) = Ubf(i,k)
		   Uvel2(i,j1+2,k) =Uvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   Uvel2(i,-2,k) = Ubf(i,k)
		   Uvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   Uvel2(i,-2,k) = Ubf(i,k)
		   Uvel2(i,j1+2,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	if (periodicx.eq.0.or.periodicx.eq.2) then
		Uvel2(-1,-2:j1+2,0:k1)=  Uvel2(0,-2:j1+2,0:k1)
		Uvel2(-2,-2:j1+2,0:k1)=  Uvel2(0,-2:j1+2,0:k1)
		Uvel2(i1+1,-2:j1+2,0:k1)=Uvel2(i1,-2:j1+2,0:k1)
		Uvel2(i1+2,-2:j1+2,0:k1)=Uvel2(i1,-2:j1+2,0:k1)
	else
		Uvel2(-1,-2:j1+2,0:k1)=  Uvel2(ie-1,-2:j1+2,0:k1)
		Uvel2(-2,-2:j1+2,0:k1)=  Uvel2(ie-2,-2:j1+2,0:k1)
		Uvel2(i1+1,-2:j1+2,0:k1)=Uvel2(2,-2:j1+2,0:k1)
		Uvel2(i1+2,-2:j1+2,0:k1)=Uvel2(3,-2:j1+2,0:k1)
	endif
!	!22-12-2017: test bc bed Uvel=0 (if it works also at immersed boundary all is fine! and with U(i,j,-1)=U(i,j,0) wiggly time-avg results for periodic turb channel MKM near bed.
	Uvel2(-2:i1+2,-2:j1+2,-1)=  Uvel2(-2:i1+2,-2:j1+2,0)
	Uvel2(-2:i1+2,-2:j1+2,-2)=  Uvel2(-2:i1+2,-2:j1+2,0)
	Uvel2(-2:i1+2,-2:j1+2,k1+1)=Uvel2(-2:i1+2,-2:j1+2,k1)
	Uvel2(-2:i1+2,-2:j1+2,k1+2)=Uvel2(-2:i1+2,-2:j1+2,k1)


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
				wf(i,j,k)=MAX(nd2,wf(i,j,k)-0.001) !wf(i,j,k)=0.
			ENDIF
		  enddo
	    enddo
	  enddo
	elseif (wd.eq.2) then 
		do i=ib,ie 
		  do j=jb,je 
			do k=kb,ke 
			fcg_local=MIN(fcg(i,j+rank*jmax,k),fcg(i+1,j+rank*jmax,k),fcg(i-1,j+rank*jmax,k),
     &			fcg(i,j+rank*jmax+1,k),fcg(i,j+rank*jmax-1,k),fcg(i,j+rank*jmax,k+1),fcg(i,j+rank*jmax,k-1))
			  wf(i,j,k)=MAX((1.-fcg_local),nd2) !high diffusion cells direct neighbours of ibm and low diffusion otherwise
			enddo 
		  enddo 
		enddo
		do i=0,i1 
		  do j=0,j1 
			do k=0,k1 
			  wf(i,j,k)=MAX((1.-fcg(i,j+rank*jmax,k)),nd2) !high diffusion all cells inside ibm and low diffusion outside ibm 
			enddo 
		    do k=0,kbed(i,j)
			  wf(i,j,k)=1. !high diffusion inside bed and first cell above 
			enddo 			
		  enddo 
		enddo
	endif ! wf=1. is already default when wiggle_detector==0	
	
	  do  i=ib,ie
		ip=i+1
		im=i-1
		ipp=i+2
		ippp=i+3
		imm=i-2
		immm=i-3
        do j=jb,je
          jp=j+1
          jm=j-1
		  jpp=j+2
		  jppp=j+3
		  jmm=j-2
		  jmmm=j-3
		  do k=kb,ke 
			kp=k+1
			km=k-1
			kpp=k+2
			kppp=k+3
			kmm=k-2
			kmmm=k-3		

!            uuR=(Uvel(i,j,k)*rhU(i,j,k)+Uvel(ip,j,k)*rhU(ip,j,k))*Rp(ip)
!            uuL=(Uvel(i,j,k)*rhU(i,j,k)+Uvel(im,j,k)*rhU(im,j,k))*Rp(i)
!            vvR=(Vvel(i,j,k)*rhV(i,j,k)+Vvel(ip,j,k)*rhV(ip,j,k))
!            vvL=(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(ip,jm,k)*rhV(ip,jm,k))
!            wwR=(Wvel(i,j,k)*rhW(i,j,k)+Wvel(ip,j,k)*rhW(ip,j,k))
!            wwL=(Wvel(i,j,km)*rhW(i,j,km)+Wvel(ip,j,km)*rhW(ip,j,km))
			
            uuR=(Uvel(i,j,k)*rhU(i,j,k)+Uvel(ip,j,k)*rhU(ip,j,k))*Rp(ip)
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(ip,j+rank*jmax,k),fcg(MIN(i1,ipp),j+rank*jmax,k))
            uuL=(Uvel(i,j,k)*rhU(i,j,k)+Uvel(im,j,k)*rhU(im,j,k))*Rp(i)
     &		*MIN(fcg(im,j+rank*jmax,k),fcg(i,j+rank*jmax,k),fcg(ip,j+rank*jmax,k))			
            vvR=(Vvel(i,j,k)*rhV(i,j,k)+Vvel(ip,j,k)*rhV(ip,j,k))
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(i,jp+rank*jmax,k),fcg(ip,j+rank*jmax,k),fcg(ip,jp+rank*jmax,k))			
            vvL=(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(ip,jm,k)*rhV(ip,jm,k))
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(i,jm+rank*jmax,k),fcg(ip,j+rank*jmax,k),fcg(ip,jm+rank*jmax,k))			
            wwR=(Wvel(i,j,k)*rhW(i,j,k)+Wvel(ip,j,k)*rhW(ip,j,k))
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(i,j+rank*jmax,kp),fcg(ip,j+rank*jmax,k),fcg(ip,j+rank*jmax,kp))			
            wwL=(Wvel(i,j,km)*rhW(i,j,km)+Wvel(ip,j,km)*rhW(ip,j,km))			
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(i,j+rank*jmax,km),fcg(ip,j+rank*jmax,k),fcg(ip,j+rank*jmax,km))			
			
				!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))	

			
		numdif2=numdif*(wf(i,j,k)+wf(i+1,j,k))*0.5

      putout(i,j,k) = - 0.25 * (
     1 (uuR    * ((Uvel2(i,j,k)+Uvel2(ip,j,k)))-
     1 numdif2*
     1 ABS(uuR)* (10.*(-Uvel2(i,j,k)+Uvel2(ip,j,k))-5.*(-Uvel2(im,j,k)+Uvel2(ipp,j,k))+(-Uvel2(imm,j,k)+Uvel2(ippp,j,k))) -
     1 (uuL    * ((Uvel2(im,j,k)+Uvel2(i,j,k)))-
     1 numdif2*
     1 ABS(uuL)* (10.*(-Uvel2(im,j,k)+Uvel2(i,j,k))-5.*(-Uvel2(imm,j,k)+Uvel2(ip,j,k))+(-Uvel2(immm,j,k)+Uvel2(ipp,j,k)))) )
     1  / ( Ru(i) * ( Rp(ip)-Rp(i) ) )
     +                         + 
     2 (vvR    * ((Uvel2(i,j,k)+Uvel2(i,jp,k)))-
     2 numdif2*
     2 ABS(vvR)* (10.*(-Uvel2(i,j,k)+Uvel2(i,jp,k))-5.*(-Uvel2(i,jm,k)+Uvel2(i,jpp,k))+(-Uvel2(i,jmm,k)+Uvel2(i,jppp,k))) -
     2 (vvL    * ((Uvel2(i,jm,k)+Uvel2(i,j,k)))-
     2 numdif2*
     2 ABS(vvL)* (10.*(-Uvel2(i,jm,k)+Uvel2(i,j,k))-5.*(-Uvel2(i,jmm,k)+Uvel2(i,jp,k))+(-Uvel2(i,jmmm,k)+Uvel2(i,jpp,k)))) )
     2  / ( Ru(i) * (phiv(j)-phiv(jm)) )
     +                         +
     3 (wwR    * ((Uvel2(i,j,k)+Uvel2(i,j,kp)))-
     3 numdif2*
     3 ABS(wwR)* (10.*(-Uvel2(i,j,k)+Uvel2(i,j,kp))-5.*(-Uvel2(i,j,km)+Uvel2(i,j,kpp))+(-Uvel2(i,j,kmm)+Uvel2(i,j,kppp))) -
     3 (wwL    * ((Uvel2(i,j,km)+Uvel2(i,j,k)))-
     3 numdif2*
     3 ABS(wwL)* (10.*(-Uvel2(i,j,km)+Uvel2(i,j,k))-5.*(-Uvel2(i,j,kmm)+Uvel2(i,j,kp))+(-Uvel2(i,j,kmmm)+Uvel2(i,j,kpp)))) )
     3  / ( dz )
     +                         -
     4 0.5*(rho(i,j,k)+rho(ip,j,k))*( Vvel(i,j,k) + Vvel(ip,j,k) + Vvel(i,jm,k) + Vvel(ip,jm,k) )**2
     4  / ( 4.0 * Ru(i) )
     +                         )
           enddo
        enddo
      enddo
      return
      end


      subroutine advecv_HYB6(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,kbed,wf,wd,fcg)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,kbed(0:i1,0:j1)
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1)
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm
	real numdif
	real uuR,uuL,vvR,vvL,wwR,wwL,uuR2,uuL2
	integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy,wd,jmax

	real*8 ubb(0:i1,0:k1),ubf(0:i1,0:k1),Vvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1),numdif2,wf(0:i1,0:j1,0:k1)
	real fcg(0:i1,0:px*(j1-1)+1,0:k1)
	
	jmax=j1-1 
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

	  do  i=ib,ie
		ip=i+1
		im=i-1
		ipp=i+2
		ippp=i+3
		imm=i-2
		immm=i-3
		do j=jb,je
          jp=j+1
          jm=j-1
		  jpp=j+2
		  jppp=j+3
		  jmm=j-2
		  jmmm=j-3
		  do k=kb,ke 
			kp=k+1
			km=k-1
			kpp=k+2
			kppp=k+3
			kmm=k-2
			kmmm=k-3		
           uuR=(Uvel(i,j,k)*rhU(i,j,k)+Uvel(i,jp,k)*rhU(i,jp,k))*Ru(i) !*Ru(i)
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(ip,j+rank*jmax,k),fcg(i,jp+rank*jmax,k),fcg(ip,jp+rank*jmax,k))
           uuL=(Uvel(im,j,k)*rhU(im,j,k)+Uvel(im,jp,k)*rhU(im,jp,k))*Ru(im) !*Ru(im)
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(im,j+rank*jmax,k),fcg(i,jp+rank*jmax,k),fcg(im,jp+rank*jmax,k))		   
    	   vvR=(Vvel(i,j,k)*rhV(i,j,k)+Vvel(i,jp,k)*rhV(i,jp,k))
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(i,jp+rank*jmax,k),fcg(i,MIN(jpp+rank*jmax,px*jmax+1),k))		   		   
           vvL=(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(i,j,k)*rhV(i,j,k))
     &		*MIN(fcg(i,jm+rank*jmax,k),fcg(i,j+rank*jmax,k),fcg(i,jp+rank*jmax,k))		   
           wwR=(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,jp,k)*rhW(i,jp,k))
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(i,j+rank*jmax,kp),fcg(i,jp+rank*jmax,k),fcg(i,jp+rank*jmax,kp))		   
           wwL=(Wvel(i,j,km)*rhW(i,j,km)+Wvel(i,jp,km)*rhW(i,jp,km))		   
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(i,j+rank*jmax,km),fcg(i,jp+rank*jmax,k),fcg(i,jp+rank*jmax,km))		   		   
			   
		   
		numdif2=numdif*(wf(i,j,k)+wf(i,j+1,k))*0.5

      !putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (

     1 (uuR    * ((Vvel2(i,j,k)+Vvel2(ip,j,k)))-
     1 numdif2*
     1 ABS(uuR)* (10.*(-Vvel2(i,j,k)+Vvel2(ip,j,k))-5.*(-Vvel2(im,j,k)+Vvel2(ipp,j,k))+(-Vvel2(imm,j,k)+Vvel2(ippp,j,k))) -
     1 (uuL    * ((Vvel2(im,j,k)+Vvel2(i,j,k)))-
     1 numdif2*
     1 ABS(uuL)* (10.*(-Vvel2(im,j,k)+Vvel2(i,j,k))-5.*(-Vvel2(imm,j,k)+Vvel2(ip,j,k))+(-Vvel2(immm,j,k)+Vvel2(ipp,j,k)))) )
     1  / ( Rp(i)* dr(i) )   !/ ( Rp(i) * Rp(i)* dr(i) ) 
!     1  / ( Rp(i)*Rp(i)* dr(i) )   !/ ( Rp(i) * Rp(i)* dr(i) ) 

!     1 (uuR    * (Vvel2(i,j,k)+Vvel2(ip,j,k))-uuL    * (Vvel2(im,j,k)+Vvel2(i,j,k)))/ ( Rp(i)*Rp(i)* dr(i) )
!     1 - (
!     1 numdif*ABS(uuR2)* (10.*(-Vvel2(i,j,k)+Vvel2(ip,j,k))-5.*(-Vvel2(im,j,k)+Vvel2(ipp,j,k))+(-Vvel2(imm,j,k)+Vvel2(ippp,j,k))) -
!     1 numdif*ABS(uuL2)* (10.*(-Vvel2(im,j,k)+Vvel2(i,j,k))-5.*(-Vvel2(imm,j,k)+Vvel2(ip,j,k))+(-Vvel2(immm,j,k)+Vvel2(ipp,j,k))) )
!     1  / ( Rp(i)* dr(i) )   !/ ( Rp(i) * Rp(i)* dr(i) ) 
     +                         + 
     2 (vvR    * ((Vvel2(i,j,k)+Vvel2(i,jp,k)))-
     2 numdif2*
     2 ABS(vvR)* (10.*(-Vvel2(i,j,k)+Vvel2(i,jp,k))-5.*(-Vvel2(i,jm,k)+Vvel2(i,jpp,k))+(-Vvel2(i,jmm,k)+Vvel2(i,jppp,k))) -
     2 (vvL    * ((Vvel2(i,jm,k)+Vvel2(i,j,k)))-
     2 numdif2*
     2 ABS(vvL)* (10.*(-Vvel2(i,jm,k)+Vvel2(i,j,k))-5.*(-Vvel2(i,jmm,k)+Vvel2(i,jp,k))+(-Vvel2(i,jmmm,k)+Vvel2(i,jpp,k)))) )
     2  / ( Rp(i) * (phiv(j)-phiv(jm)) )
     +                         +
     3 (wwR    * ((Vvel2(i,j,k)+Vvel2(i,j,kp)))-
     3 numdif2*
     3 ABS(wwR)* (10.*(-Vvel2(i,j,k)+Vvel2(i,j,kp))-5.*(-Vvel2(i,j,km)+Vvel2(i,j,kpp))+(-Vvel2(i,j,kmm)+Vvel2(i,j,kppp))) -
     3 (wwL    * ((Vvel2(i,j,km)+Vvel2(i,j,k)))-
     3 numdif2*
     3 ABS(wwL)* (10.*(-Vvel2(i,j,km)+Vvel2(i,j,k))-5.*(-Vvel2(i,j,kmm)+Vvel2(i,j,kp))+(-Vvel2(i,j,kmmm)+Vvel2(i,j,kpp)))) )
     3  / ( dz )

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

      subroutine advecw_HYB6(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phiv,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,kbed,wf,wd,fcg)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,kbed(0:i1,0:j1)
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1)
      real rhoip,rhoim,rhojp,rhojm
	real uuR,uuL,vvR,vvL,wwR,wwL,numdif
	integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy,wd,jmax
	real*8 ubb(0:i1,0:k1),ubf(0:i1,0:k1),Wvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1),numdif2,wf(0:i1,0:j1,0:k1)
	real fcg(0:i1,0:px*(j1-1)+1,0:k1)
	
	jmax=j1-1 

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

	do  i=ib,ie
		ip=i+1
		im=i-1
		ipp=i+2
		ippp=i+3
		imm=i-2
		immm=i-3
        do j=jb,je
          jp=j+1
          jm=j-1
		  jpp=j+2
		  jppp=j+3
		  jmm=j-2
		  jmmm=j-3
		  do k=kb,ke 
			kp=k+1
			km=k-1
			kpp=k+2
			kppp=k+3
			kmm=k-2
			kmmm=k-3		  

      	    uuR=(Uvel(i,j,k)*rhU(i,j,k)+Uvel(i,j,kp)*rhU(i,j,kp))*Ru(i)
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(ip,j+rank*jmax,k),fcg(i,j+rank*jmax,kp),fcg(ip,j+rank*jmax,kp))			
      	    uuL=(Uvel(im,j,k)*rhU(im,j,k)+Uvel(im,j,kp)*rhU(im,j,kp))*Ru(im)
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(im,j+rank*jmax,k),fcg(i,j+rank*jmax,kp),fcg(im,j+rank*jmax,kp))			
			vvR=(Vvel(i,j,k)*rhV(i,j,k)+Vvel(i,j,kp)*rhV(i,j,kp))
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(i,jp+rank*jmax,k),fcg(i,j+rank*jmax,kp),fcg(i,jp+rank*jmax,kp))			
    	    vvL=(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(i,jm,kp)*rhV(i,jm,kp))
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(i,jm+rank*jmax,k),fcg(i,j+rank*jmax,kp),fcg(i,jm+rank*jmax,kp))						
    	    wwR=(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,j,kp)*rhW(i,j,kp))
     &		*MIN(fcg(i,j+rank*jmax,k),fcg(i,j+rank*jmax,kp),fcg(i,j+rank*jmax,MIN(k1,kpp)))						
			wwL=(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,j,km)*rhW(i,j,km))
     &		*MIN(fcg(i,j+rank*jmax,km),fcg(i,j+rank*jmax,k),fcg(i,j+rank*jmax,kp))						

			
			numdif2=numdif*(wf(i,j,k)+wf(i,j,k+1))*0.5

      !putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (
     1 (uuR    * ((Wvel2(i,j,k)+Wvel2(ip,j,k)))-
     1 numdif2*
     1 ABS(uuR)* (10.*(-Wvel2(i,j,k)+Wvel2(ip,j,k))-5.*(-Wvel2(im,j,k)+Wvel2(ipp,j,k))+(-Wvel2(imm,j,k)+Wvel2(ippp,j,k))) -
     1 (uuL    * ((Wvel2(im,j,k)+Wvel2(i,j,k)))-
     1 numdif2*
     1 ABS(uuL)* (10.*(-Wvel2(im,j,k)+Wvel2(i,j,k))-5.*(-Wvel2(imm,j,k)+Wvel2(ip,j,k))+(-Wvel2(immm,j,k)+Wvel2(ipp,j,k)))) )
     1  / ( Rp(i)* dr(i) ) 
     +                         + 
     2 (vvR    * ((Wvel2(i,j,k)+Wvel2(i,jp,k)))-
     2 numdif2*
     2 ABS(vvR)* (10.*(-Wvel2(i,j,k)+Wvel2(i,jp,k))-5.*(-Wvel2(i,jm,k)+Wvel2(i,jpp,k))+(-Wvel2(i,jmm,k)+Wvel2(i,jppp,k))) -
     2 (vvL    * ((Wvel2(i,jm,k)+Wvel2(i,j,k)))-
     2 numdif2*
     2 ABS(vvL)* (10.*(-Wvel2(i,jm,k)+Wvel2(i,j,k))-5.*(-Wvel2(i,jmm,k)+Wvel2(i,jp,k))+(-Wvel2(i,jmmm,k)+Wvel2(i,jpp,k)))) )
     2  / ( Rp(i) * (phiv(j)-phiv(jm)) )
     +                         +
     3 (wwR    * ((Wvel2(i,j,k)+Wvel2(i,j,kp)))-
     3 numdif2*
     3 ABS(wwR)* (10.*(-Wvel2(i,j,k)+Wvel2(i,j,kp))-5.*(-Wvel2(i,j,km)+Wvel2(i,j,kpp))+(-Wvel2(i,j,kmm)+Wvel2(i,j,kppp))) -
     3 (wwL    * ((Wvel2(i,j,km)+Wvel2(i,j,k)))-
     3 numdif2*
     3 ABS(wwL)* (10.*(-Wvel2(i,j,km)+Wvel2(i,j,k))-5.*(-Wvel2(i,j,kmm)+Wvel2(i,j,kp))+(-Wvel2(i,j,kmmm)+Wvel2(i,j,kpp)))) )
     3  / ( dz )

     +                         )
           enddo
         enddo
      enddo
      return
      end


      subroutine advecu_TVD(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt,periodicx,periodicy,itel)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,px,itel
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phip(0:j1),phipt(0:je*px+1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter 
      real rho(0:i1,0:j1,0:k1),cfll,limiter3,dt
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,Ru2(-1:i1+1)
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kmm,kpp,imm,ipp,jmm,jpp,rank,periodicx,periodicy
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-1:i1+1,-1:j1+1,0:k1)
      real 	dz_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)	

	Ru2(0:i1)=Ru
	Ru2(-1)=Ru(0)-(Ru(2)-Ru(1)) !needed for periodicx sim
 	Ru2(i1+1)=Ru(i1)+(Ru(i1)-Ru(i1-1)) !needed for periodicx sim	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim	

	putin2(0:i1,0:j1,0:k1)=Uvel
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=Uvel(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Uvel(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=Uvel(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Uvel(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
	  call shiftf2(Uvel,pbf)
	  call shiftb2(Uvel,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = Uvel(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =Uvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif
	
		dz_i=1./dz
      do k=kb,ke
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
        do j=jb,je
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
          do  i=ib,ie
			  ip=i+1
			  im=i-1
			  imm=im
			  ipp=ip
			  if (periodicx.eq.0.or.periodicx.eq.2) then
				if (i.eq.1) imm=i
				if (i.eq.ie) ipp=i
			  endif
			  varx_grid_Rpos=(Ru2(ip)-Ru2(i))/(Ru2(i)-Ru2(im))
			  varx_grid_Rneg=(Ru2(ipp)-Ru2(i))/(Ru2(ipp+1)-Ru2(ipp))
			  varx_grid_Lpos=(Ru2(i)-Ru2(imm))/(Ru2(imm)-Ru2(imm-1))
			  varx_grid_Lneg=(Ru2(i)-Ru2(im))/(Ru2(ip)-Ru2(i))

			if (mod(itel,2).eq.0) then
			! divergence version: d(ruu)dx
			rhoip = rho(ip,j,k)
			rhoim = rho(i ,j,k)			
            rhojp =0.25*(rho(i,j,k)+rho(i,jp,k)+rho(ip,j,k)+rho(ip,jp,k))
            rhojm =0.25*(rho(i,j,k)+rho(i,jm,k)+rho(ip,j,k)+rho(ip,jm,k))
            rhokp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
            rhokm =0.25*(rho(i,j,k)+rho(i,j,km)+rho(ip,j,k)+rho(ip,j,km))			
            uuR=(Uvel(i,j,k)+Uvel(ip,j,k))*0.5
            uuL=(Uvel(i,j,k)+Uvel(im,j,k))*0.5
            vvR=(Vvel(i,j,k)+Vvel(ip,j,k))*0.5
            vvL=(Vvel(i,jm,k)+Vvel(ip,jm,k))*0.5
            wwR=(Wvel(i,j,k)+Wvel(ip,j,k))*0.5
            wwL=(Wvel(i,j,km)+Wvel(ip,j,km))*0.5
			else

			! advective version: rud(u)dx
			rhoip = (rho(ip,j,k)+rho(i ,j,k))*0.5
			rhoim = (rho(ip,j,k)+rho(i ,j,k))*0.5		
            rhojp = (rho(ip,j,k)+rho(i ,j,k))*0.5
            rhojm = (rho(ip,j,k)+rho(i ,j,k))*0.5
            rhokp = (rho(ip,j,k)+rho(i ,j,k))*0.5
            rhokm = (rho(ip,j,k)+rho(i ,j,k))*0.5			
            uuR=(Uvel(ip,j,k)+Uvel(im,j,k)+2.*Uvel(i,j,k))*0.25 !Uvel(i,j,k)
            uuL=(Uvel(ip,j,k)+Uvel(im,j,k)+2.*Uvel(i,j,k))*0.25 !Uvel(i,j,k)
            vvR=(Vvel(i,j,k)+Vvel(ip,j,k)+Vvel(i,jm,k)+Vvel(ip,jm,k))*0.25
            vvL=(Vvel(i,j,k)+Vvel(ip,j,k)+Vvel(i,jm,k)+Vvel(ip,jm,k))*0.25
            wwR=(Wvel(i,j,k)+Wvel(ip,j,k)+Wvel(i,j,km)+Wvel(ip,j,km))*0.25
            wwL=(Wvel(i,j,k)+Wvel(ip,j,k)+Wvel(i,j,km)+Wvel(ip,j,km))*0.25	
			endif			
            putout(i,j,k) = 0.0
			IF (uuR.ge.0.) THEN
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuR)/(Ru(ip)-Ru(i ))
			  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter3(rRpos,cfll)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = - Rp(ip)*rhoip*uuR*cRpos/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ELSE
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuR)/(Ru(ip)-Ru(i ))
			  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
			  cRneg = putin2(ip,j,k) + 0.5*limiter3(rRneg,cfll)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfll)
			  putout(i,j,k) = - Rp(ip)*rhoip*uuR*cRneg/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ENDIF	
			IF (uuL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuL)/(Ru(i)-Ru(im))
			  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
			  cLpos = putin2(im,j,k) + 0.5*limiter3(rLpos,cfll)*(putin2( i,j,k) - putin2(im,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + Rp(i)*rhoim*uuL*cLpos/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ELSE
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuL)/(Ru(i)-Ru(im))
			  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
			  cLneg = putin2(i ,j,k) + 0.5*limiter3(rLneg,cfll)*(putin2(im,j,k) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + Rp(i)*rhoim*uuL*cLneg/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ENDIF	
			IF (vvR.ge.0.) THEN
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter3(rRpos,cfll)*(putin2(i,jp,k) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhojp*vvR*cRpos/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
			  cRneg = putin2(i,jp,k) + 0.5*limiter3(rRneg,cfll)*(putin2( i,j,k) - putin2(i,jp,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhojp*vvR*cRneg/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF		
			IF (vvL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
			  cLpos = putin2(i,jm,k) + 0.5*limiter3(rLpos,cfll)*(putin2( i,j,k) - putin2(i,jm,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhojm*vvL*cLpos/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
			  cLneg = putin2(i,j ,k) + 0.5*limiter3(rLneg,cfll)*(putin2(i,jm,k) - putin2(i,j ,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhojm*vvL*cLneg/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF	
			IF (wwR.ge.0.) THEN
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwR)*dz_i
			  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
			  cRpos = putin2(i ,j,k) + 0.5*limiter3(rRpos,cfll)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhokp*wwR*cRpos*dz_i 
			ELSE
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwR)*dz_i
			  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
			  cRneg = putin2(i,j,kp) + 0.5*limiter3(rRneg,cfll)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhokp*wwR*cRneg*dz_i 
			ENDIF
			IF (wwL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwL)*dz_i
			  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
			  cLpos = putin2(i,j,km) + 0.5*limiter3(rLpos,cfll)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhokm*wwL*cLpos*dz_i 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwL)*dz_i
			  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
			  cLneg = putin2(i,j,k ) + 0.5*limiter3(rLneg,cfll)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhokm*wwL*cLneg*dz_i 
			ENDIF
			putout(i,j,k) = putout(i,j,k) 
     & + (rho(i,j,k)+rho(ip,j,k))*(Vvel(i,j,k) + Vvel(ip,j,k) + Vvel(i,jm,k) + Vvel(ip,jm,k))**2/(32.0*Ru(i))
           enddo
        enddo
      enddo
      return
      end


      subroutine advecv_TVD(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt,periodicx,periodicy,itel)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,px,itel
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phip(0:j1),phipt(0:je*px+1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),	 
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter 
      real rho(0:i1,0:j1,0:k1),cfll,limiter3,dt
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,Rp2(-1:i1+1)
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kmm,kpp,imm,ipp,jmm,jpp,rank,periodicx,periodicy
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-1:i1+1,-1:j1+1,0:k1)
      real 	dz_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phivt2(-1:je*px+2)	

	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim	
	phivt2(0:je*px+1)=phivt
	phivt2(-1)=phivt(0)-(phivt(1)-phivt(0)) !needed for periodicy sim
 	phivt2(je*px+1+1)=phivt(je*px+1)+(phivt(je*px+1)-phivt(je*px+1-1)) !needed for periodicy sim	

	putin2(0:i1,0:j1,0:k1)=Vvel
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=Vvel(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Vvel(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=Vvel(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Vvel(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
	  call shiftf2(Vvel,pbf)
	  call shiftb2(Vvel,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = Vvel(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =Vvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif
	
		dz_i=1./dz
      do k=kb,ke
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
        do j=jb,je
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif
      vary_grid_Rpos=(phivt2(rank*je+jp) -phivt2(rank*je+j))  /(phivt2(rank*je+j)    -phivt2(rank*je+jm))
      vary_grid_Rneg=(phivt2(rank*je+jpp)-phivt2(rank*je+j))  /(phivt2(rank*je+jpp+1)-phivt2(rank*je+jpp))
      vary_grid_Lpos=(phivt2(rank*je+j)  -phivt2(rank*je+jmm))/(phivt2(rank*je+jmm)  -phivt2(rank*je+jmm-1))
      vary_grid_Lneg=(phivt2(rank*je+j)  -phivt2(rank*je+jm)) /(phivt2(rank*je+jp)   -phivt2(rank*je+j))	
	  
	  
          do  i=ib,ie
			  ip=i+1
			  im=i-1
			  imm=im
			  ipp=ip
			  if (periodicx.eq.0.or.periodicx.eq.2) then
				if (i.eq.1) imm=i
				if (i.eq.ie) ipp=i
			  endif
			  varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
			  varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
			  varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
			  varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))			  
			  
			if (mod(itel,2).eq.0) then
			! divergence version: d(ruu)dx  
			  rhoip =0.25*(rho(i,j,k)+rho(ip,j,k)+rho(i,jp,k)+rho(ip,jp,k))
			  rhoim =0.25*(rho(i,j,k)+rho(im,j,k)+rho(i,jp,k)+rho(im,jp,k))
			  rhojp =rho(i,jp,k)
			  rhojm =rho(i,j ,k)
			  rhokp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
			  rhokm =0.25*(rho(i,j,k)+rho(i,j,km)+rho(i,jp,k)+rho(i,jp,km))			
            uuR=(Uvel(i,j,k)+Uvel(i,jp,k))*0.5
            uuL=(Uvel(im,j,k)+Uvel(im,jp,k))*0.5
      	    vvR=(Vvel(i,j,k)+Vvel(i,jp,k))*0.5
            vvL=(Vvel(i,jm,k)+Vvel(i,j,k))*0.5
            wwR=(Wvel(i,j,k)+Wvel(i,jp,k))*0.5
            wwL=(Wvel(i,j,km)+Wvel(i,jp,km))*0.5
			else
			! advective version: rud(u)dx
			  rhoip =(rho(i,jp,k)+rho(i,j ,k))*0.5
			  rhoim =(rho(i,jp,k)+rho(i,j ,k))*0.5
			  rhojp =(rho(i,jp,k)+rho(i,j ,k))*0.5
			  rhojm =(rho(i,jp,k)+rho(i,j ,k))*0.5
			  rhokp =(rho(i,jp,k)+rho(i,j ,k))*0.5
			  rhokm =(rho(i,jp,k)+rho(i,j ,k))*0.5			
            uuR=(Uvel(i,j,k)+Uvel(i,jp,k)+Uvel(im,j,k)+Uvel(im,jp,k))*0.25
            uuL=(Uvel(i,j,k)+Uvel(i,jp,k)+Uvel(im,j,k)+Uvel(im,jp,k))*0.25
      	    vvR=(Vvel(i,jp,k)+Vvel(i,jm,k)+2.*Vvel(i,j,k))*0.25 !Vvel(i,j,k)
            vvL=(Vvel(i,jp,k)+Vvel(i,jm,k)+2.*Vvel(i,j,k))*0.25 !Vvel(i,j,k)
            wwR=(Wvel(i,j,k)+Wvel(i,jp,k)+Wvel(i,j,km)+Wvel(i,jp,km))*0.25
            wwL=(Wvel(i,j,k)+Wvel(i,jp,k)+Wvel(i,j,km)+Wvel(i,jp,km))*0.25		
			endif			
			putout(i,j,k) = 0.0
			IF (uuR.ge.0.) THEN
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter3(rRpos,cfll)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = - Ru(i)*rhoip*uuR*cRpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
			  cRneg = putin2(ip,j,k) + 0.5*limiter3(rRneg,cfll)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfll)
			  putout(i,j,k) = - Ru(i)*rhoip*uuR*cRneg/ ( Rp(i)* dr(i) )
			ENDIF	
			IF (uuL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
			  cLpos = putin2(im,j,k) + 0.5*limiter3(rLpos,cfll)*(putin2( i,j,k) - putin2(im,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + Ru(im)*rhoim*uuL*cLpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
			  cLneg = putin2(i ,j,k) + 0.5*limiter3(rLneg,cfll)*(putin2(im,j,k) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + Ru(im)*rhoim*uuL*cLneg/ ( Rp(i)* dr(i) )
			ENDIF	
			IF (vvR.ge.0.) THEN
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvR)/(Rp(i)*(phivt2(rank*je+jp)-phivt2(rank*je+j)))
			  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter3(rRpos,cfll)*(putin2(i,jp,k) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhojp*vvR*cRpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvR)/(Rp(i)*(phivt2(rank*je+jp)-phivt2(rank*je+j)))
			  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
			  cRneg = putin2(i,jp,k) + 0.5*limiter3(rRneg,cfll)*(putin2( i,j,k) - putin2(i,jp,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhojp*vvR*cRneg/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF		
			IF (vvL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvL)/(Rp(i)*(phivt2(rank*je+j)-phivt2(rank*je+jm)))
			  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
			  cLpos = putin2(i,jm,k) + 0.5*limiter3(rLpos,cfll)*(putin2( i,j,k) - putin2(i,jm,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhojm*vvL*cLpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvL)/(Rp(i)*(phivt2(rank*je+j)-phivt2(rank*je+jm)))
			  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
			  cLneg = putin2(i,j ,k) + 0.5*limiter3(rLneg,cfll)*(putin2(i,jm,k) - putin2(i,j ,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhojm*vvL*cLneg/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF	
			IF (wwR.ge.0.) THEN
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwR)*dz_i
			  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
			  cRpos = putin2(i ,j,k) + 0.5*limiter3(rRpos,cfll)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhokp*wwR*cRpos*dz_i 
			ELSE
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwR)*dz_i
			  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
			  cRneg = putin2(i,j,kp) + 0.5*limiter3(rRneg,cfll)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhokp*wwR*cRneg*dz_i 
			ENDIF
			IF (wwL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwL)*dz_i
			  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
			  cLpos = putin2(i,j,km) + 0.5*limiter3(rLpos,cfll)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhokm*wwL*cLpos*dz_i 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwL)*dz_i
			  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
			  cLneg = putin2(i,j,k ) + 0.5*limiter3(rLneg,cfll)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhokm*wwL*cLneg*dz_i 
			ENDIF
			putout(i,j,k) = putout(i,j,k) -
     4   (rho(i,j,k)+rho(i,jp,k))*(Uvel(im,j,k)+Uvel(i,j,k)+Uvel(im,jp,k)+Uvel(i,jp,k))*Vvel(i,j,k)
     4   / ( 8.*Rp(i) )

           enddo
        enddo
      enddo
      return
      end	  



      subroutine advecw_TVD(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt,periodicx,periodicy,itel)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,px,itel
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phip(0:j1),phipt(0:je*px+1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter 
      real rho(0:i1,0:j1,0:k1),cfll,limiter3,dt
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,Rp2(-1:i1+1)
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kmm,kpp,imm,ipp,jmm,jpp,rank,periodicx,periodicy
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-1:i1+1,-1:j1+1,0:k1)
      real 	dz_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)	

	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim	

	putin2(0:i1,0:j1,0:k1)=Wvel
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=Wvel(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Wvel(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=Wvel(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Wvel(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
	  call shiftf2(Wvel,pbf)
	  call shiftb2(Wvel,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = Wvel(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =Wvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif
	
		dz_i=1./dz
      do k=kb,ke
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
        do j=jb,je
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
          do  i=ib,ie
			  ip=i+1
			  im=i-1
			  imm=im
			  ipp=ip
			  if (periodicx.eq.0.or.periodicx.eq.2) then
				if (i.eq.1) imm=i
				if (i.eq.ie) ipp=i
			  endif
			  varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
			  varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
			  varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
			  varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))	
			if (mod(itel,2).eq.0) then
			! divergence version: d(ruu)dx   
			rhoip =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
            rhoim =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(im,j,k)+rho(im,j,kp))
            rhojp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
            rhojm =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jm,k)+rho(i,jm,kp))
			rhokp = rho(i,j,kp)
			rhokm = rho(i,j,k)
      	    uuR=(Uvel(i,j,k)+Uvel(i,j,kp))*0.5
      	    uuL=(Uvel(im,j,k)+Uvel(im,j,kp))*0.5
			vvR=(Vvel(i,j,k)+Vvel(i,j,kp))*0.5
    	    vvL=(Vvel(i,jm,k)+Vvel(i,jm,kp))*0.5
    	    wwR=(Wvel(i,j,k)+Wvel(i,j,kp))*0.5
			wwL=(Wvel(i,j,k)+Wvel(i,j,km))*0.5
			else
			! advective version: rud(u)dx   
			rhoip =(rho(i,j,k)+rho(i,j,kp))*0.5
            rhoim =(rho(i,j,k)+rho(i,j,kp))*0.5
            rhojp =(rho(i,j,k)+rho(i,j,kp))*0.5
            rhojm =(rho(i,j,k)+rho(i,j,kp))*0.5
			rhokp =(rho(i,j,k)+rho(i,j,kp))*0.5
			rhokm =(rho(i,j,k)+rho(i,j,kp))*0.5
      	    uuR=(Uvel(i,j,k)+Uvel(i,j,kp)+Uvel(im,j,k)+Uvel(im,j,kp))*0.25
      	    uuL=(Uvel(i,j,k)+Uvel(i,j,kp)+Uvel(im,j,k)+Uvel(im,j,kp))*0.25
			vvR=(Vvel(i,j,k)+Vvel(i,j,kp)+Vvel(i,jm,k)+Vvel(i,jm,kp))*0.25
    	    vvL=(Vvel(i,j,k)+Vvel(i,j,kp)+Vvel(i,jm,k)+Vvel(i,jm,kp))*0.25
    	    wwR=(Wvel(i,j,kp)+Wvel(i,j,km)+2.*Wvel(i,j,k))*0.25 !Wvel(i,j,k)
			wwL=(Wvel(i,j,kp)+Wvel(i,j,km)+2.*Wvel(i,j,k))*0.25 !Wvel(i,j,k)
			endif
            putout(i,j,k) = 0.0

			IF (uuR.ge.0.) THEN
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter3(rRpos,cfll)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = - Ru(i)*rhoip*uuR*cRpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
			  cRneg = putin2(ip,j,k) + 0.5*limiter3(rRneg,cfll)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfll)
			  putout(i,j,k) = - Ru(i)*rhoip*uuR*cRneg/( Rp(i)* dr(i) )
			ENDIF	
			IF (uuL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
			  cLpos = putin2(im,j,k) + 0.5*limiter3(rLpos,cfll)*(putin2( i,j,k) - putin2(im,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + Ru(im)*rhoim*uuL*cLpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
			  cLneg = putin2(i ,j,k) + 0.5*limiter3(rLneg,cfll)*(putin2(im,j,k) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + Ru(im)*rhoim*uuL*cLneg/ ( Rp(i)* dr(i) )
			ENDIF	
			IF (vvR.ge.0.) THEN
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter3(rRpos,cfll)*(putin2(i,jp,k) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhojp*vvR*cRpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
			  cRneg = putin2(i,jp,k) + 0.5*limiter3(rRneg,cfll)*(putin2( i,j,k) - putin2(i,jp,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhojp*vvR*cRneg/ ( Rp(i) * (phiv(j)-phiv(jm)) )
			ENDIF		
			IF (vvL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
			  cLpos = putin2(i,jm,k) + 0.5*limiter3(rLpos,cfll)*(putin2( i,j,k) - putin2(i,jm,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhojm*vvL*cLpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
			  cLneg = putin2(i,j ,k) + 0.5*limiter3(rLneg,cfll)*(putin2(i,jm,k) - putin2(i,j ,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhojm*vvL*cLneg/ ( Rp(i) * (phiv(j)-phiv(jm)) )
			ENDIF	
			IF (wwR.ge.0.) THEN
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwR)*dz_i
			  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
			  cRpos = putin2(i ,j,k) + 0.5*limiter3(rRpos,cfll)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhokp*wwR*cRpos*dz_i 
			ELSE
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwR)*dz_i
			  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
			  cRneg = putin2(i,j,kp) + 0.5*limiter3(rRneg,cfll)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - rhokp*wwR*cRneg*dz_i 
			ENDIF
			IF (wwL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwL)*dz_i
			  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
			  cLpos = putin2(i,j,km) + 0.5*limiter3(rLpos,cfll)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhokm*wwL*cLpos*dz_i 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  cfll=dt*ABS(wwL)*dz_i
			  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
			  cLneg = putin2(i,j,k ) + 0.5*limiter3(rLneg,cfll)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + rhokm*wwL*cLneg*dz_i 
			ENDIF
           enddo
        enddo
      enddo
      return
      end



      subroutine advecu_C2Blend(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt,periodicx,periodicy,numdif)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phip(0:j1),phipt(0:je*px+1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter 
      real rho(0:i1,0:j1,0:k1),cfll,limiter5,dt,numdif
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,Ru2(-1:i1+1)
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kmm,kpp,imm,ipp,jmm,jpp,rank,periodicx,periodicy
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-1:i1+1,-1:j1+1,0:k1)
      real 	dz_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)	
	  real Ax,Bx,DDx,blend,cLho,cRho
	  real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1)	  

	Ru2(0:i1)=Ru
	Ru2(-1)=Ru(0)-(Ru(2)-Ru(1)) !needed for periodicx sim
 	Ru2(i1+1)=Ru(i1)+(Ru(i1)-Ru(i1-1)) !needed for periodicx sim	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim	

	putin2(0:i1,0:j1,0:k1)=Uvel
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=Uvel(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Uvel(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=Uvel(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Uvel(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
	  call shiftf2(Uvel,pbf)
	  call shiftb2(Uvel,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = Uvel(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =Uvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif
	
		dz_i=1./dz
      do k=kb,ke
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
        do j=jb,je
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
          do  i=ib,ie
			  ip=i+1
			  im=i-1
			  imm=im
			  ipp=ip
			  if (periodicx.eq.0.or.periodicx.eq.2) then
				if (i.eq.1) imm=i
				if (i.eq.ie) ipp=i
			  endif
			  varx_grid_Rpos=(Ru2(ip)-Ru2(i))/(Ru2(i)-Ru2(im))
			  varx_grid_Rneg=(Ru2(ipp)-Ru2(i))/(Ru2(ipp+1)-Ru2(ipp))
			  varx_grid_Lpos=(Ru2(i)-Ru2(imm))/(Ru2(imm)-Ru2(imm-1))
			  varx_grid_Lneg=(Ru2(i)-Ru2(im))/(Ru2(ip)-Ru2(i))
			
            uuR=0.5*(Uvel(i,j,k)*rhU(i,j,k)+Uvel(ip,j,k)*rhU(ip,j,k))*Rp(ip)
            uuL=0.5*(Uvel(i,j,k)*rhU(i,j,k)+Uvel(im,j,k)*rhU(im,j,k))*Rp(i)
            vvR=0.5*(Vvel(i,j,k)*rhV(i,j,k)+Vvel(ip,j,k)*rhV(ip,j,k))
            vvL=0.5*(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(ip,jm,k)*rhV(ip,jm,k))
            wwR=0.5*(Wvel(i,j,k)*rhW(i,j,k)+Wvel(ip,j,k)*rhW(ip,j,k))
            wwL=0.5*(Wvel(i,j,km)*rhW(i,j,km)+Wvel(ip,j,km)*rhW(ip,j,km))			
			
            putout(i,j,k) = 0.0
			IF (uuR.ge.0.) THEN
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Ru(ip)-Ru(i ))
			  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter5(rRpos,numdif)*(putin2(ip,j,k) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = - uuR*cRpos/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ELSE
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Ru(ip)-Ru(i ))
			  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
			  cRneg = putin2(ip,j,k) + 0.5*limiter5(rRneg,numdif)*(putin2( i,j,k) - putin2(ip,j,k))!*(1.-cfll)
			  putout(i,j,k) = - uuR*cRneg/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ENDIF	
			IF (uuL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Ru(i)-Ru(im))
			  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
			  cLpos = putin2(im,j,k) + 0.5*limiter5(rLpos,numdif)*(putin2( i,j,k) - putin2(im,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + uuL*cLpos/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ELSE
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Ru(i)-Ru(im))
			  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
			  cLneg = putin2(i ,j,k) + 0.5*limiter5(rLneg,numdif)*(putin2(im,j,k) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + uuL*cLneg/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ENDIF	
			IF (vvR.ge.0.) THEN
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter5(rRpos,numdif)*(putin2(i,jp,k) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - vvR*cRpos/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
			  cRneg = putin2(i,jp,k) + 0.5*limiter5(rRneg,numdif)*(putin2( i,j,k) - putin2(i,jp,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - vvR*cRneg/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF		
			IF (vvL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
			  cLpos = putin2(i,jm,k) + 0.5*limiter5(rLpos,numdif)*(putin2( i,j,k) - putin2(i,jm,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + vvL*cLpos/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
			  cLneg = putin2(i,j ,k) + 0.5*limiter5(rLneg,numdif)*(putin2(i,jm,k) - putin2(i,j ,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + vvL*cLneg/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF	
			IF (wwR.ge.0.) THEN
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
			  cRpos = putin2(i ,j,k) + 0.5*limiter5(rRpos,numdif)*(putin2(i,j,kp) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - wwR*cRpos*dz_i 
			ELSE
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
			  cRneg = putin2(i,j,kp) + 0.5*limiter5(rRneg,numdif)*(putin2(i,j, k) - putin2(i,j,kp))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - wwR*cRneg*dz_i 
			ENDIF
			IF (wwL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
			  cLpos = putin2(i,j,km) + 0.5*limiter5(rLpos,numdif)*(putin2(i,j, k) - putin2(i,j,km))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + wwL*cLpos*dz_i 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
			  cLneg = putin2(i,j,k ) + 0.5*limiter5(rLneg,numdif)*(putin2(i,j,km) - putin2(i,j, k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + wwL*cLneg*dz_i 
			ENDIF
			putout(i,j,k) = putout(i,j,k) 
     & + (rho(i,j,k)+rho(ip,j,k))*(Vvel(i,j,k) + Vvel(ip,j,k) + Vvel(i,jm,k) + Vvel(ip,jm,k))**2/(32.0*Ru(i))
           enddo
        enddo
      enddo
      return
      end


      subroutine advecv_C2blend(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt,periodicx,periodicy,numdif)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phip(0:j1),phipt(0:je*px+1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),	 
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter 
      real rho(0:i1,0:j1,0:k1),cfll,limiter5,dt,numdif
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,Rp2(-1:i1+1)
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kmm,kpp,imm,ipp,jmm,jpp,rank,periodicx,periodicy
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-1:i1+1,-1:j1+1,0:k1)
      real 	dz_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phivt2(-1:je*px+2)	
	  real Ax,Bx,DDx,blend,cLho,cRho
	  real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1)	  

	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim	
	phivt2(0:je*px+1)=phivt
	phivt2(-1)=phivt(0)-(phivt(1)-phivt(0)) !needed for periodicy sim
 	phivt2(je*px+1+1)=phivt(je*px+1)+(phivt(je*px+1)-phivt(je*px+1-1)) !needed for periodicy sim	

	putin2(0:i1,0:j1,0:k1)=Vvel
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=Vvel(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Vvel(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=Vvel(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Vvel(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
	  call shiftf2(Vvel,pbf)
	  call shiftb2(Vvel,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = Vvel(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =Vvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif
	
		dz_i=1./dz
      do k=kb,ke
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
        do j=jb,je
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif
      vary_grid_Rpos=(phivt2(rank*je+jp) -phivt2(rank*je+j))  /(phivt2(rank*je+j)    -phivt2(rank*je+jm))
      vary_grid_Rneg=(phivt2(rank*je+jpp)-phivt2(rank*je+j))  /(phivt2(rank*je+jpp+1)-phivt2(rank*je+jpp))
      vary_grid_Lpos=(phivt2(rank*je+j)  -phivt2(rank*je+jmm))/(phivt2(rank*je+jmm)  -phivt2(rank*je+jmm-1))
      vary_grid_Lneg=(phivt2(rank*je+j)  -phivt2(rank*je+jm)) /(phivt2(rank*je+jp)   -phivt2(rank*je+j))	
	  
	  
          do  i=ib,ie
			  ip=i+1
			  im=i-1
			  imm=im
			  ipp=ip
			  if (periodicx.eq.0.or.periodicx.eq.2) then
				if (i.eq.1) imm=i
				if (i.eq.ie) ipp=i
			  endif
			  varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
			  varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
			  varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
			  varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))			  
			  
           uuR=0.5*(Uvel(i,j,k)*rhU(i,j,k)+Uvel(i,jp,k)*rhU(i,jp,k))*Ru(i) !*Ru(i)
           uuL=0.5*(Uvel(im,j,k)*rhU(im,j,k)+Uvel(im,jp,k)*rhU(im,jp,k))*Ru(im) !*Ru(im)
     	   vvR=0.5*(Vvel(i,j,k)*rhV(i,j,k)+Vvel(i,jp,k)*rhV(i,jp,k))
           vvL=0.5*(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(i,j,k)*rhV(i,j,k))
           wwR=0.5*(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,jp,k)*rhW(i,jp,k))
           wwL=0.5*(Wvel(i,j,km)*rhW(i,j,km)+Wvel(i,jp,km)*rhW(i,jp,km))			
			
			putout(i,j,k) = 0.0
			IF (uuR.ge.0.) THEN
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter5(rRpos,numdif)*(putin2(ip,j,k) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = - uuR*cRpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
			  cRneg = putin2(ip,j,k) + 0.5*limiter5(rRneg,numdif)*(putin2( i,j,k) - putin2(ip,j,k))!*(1.-cfll)
			  putout(i,j,k) = - uuR*cRneg/ ( Rp(i)* dr(i) )
			ENDIF	
			IF (uuL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
			  cLpos = putin2(im,j,k) + 0.5*limiter5(rLpos,numdif)*(putin2( i,j,k) - putin2(im,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + uuL*cLpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
			  cLneg = putin2(i ,j,k) + 0.5*limiter5(rLneg,numdif)*(putin2(im,j,k) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + uuL*cLneg/ ( Rp(i)* dr(i) )
			ENDIF	
			IF (vvR.ge.0.) THEN
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Rp(i)*(phivt2(rank*je+jp)-phivt2(rank*je+j)))
			  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter5(rRpos,numdif)*(putin2(i,jp,k) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - vvR*cRpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Rp(i)*(phivt2(rank*je+jp)-phivt2(rank*je+j)))
			  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
			  cRneg = putin2(i,jp,k) + 0.5*limiter5(rRneg,numdif)*(putin2( i,j,k) - putin2(i,jp,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - vvR*cRneg/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF		
			IF (vvL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Rp(i)*(phivt2(rank*je+j)-phivt2(rank*je+jm)))
			  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
			  cLpos = putin2(i,jm,k) + 0.5*limiter5(rLpos,numdif)*(putin2( i,j,k) - putin2(i,jm,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + vvL*cLpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Rp(i)*(phivt2(rank*je+j)-phivt2(rank*je+jm)))
			  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
			  cLneg = putin2(i,j ,k) + 0.5*limiter5(rLneg,numdif)*(putin2(i,jm,k) - putin2(i,j ,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + vvL*cLneg/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF	
			IF (wwR.ge.0.) THEN
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
			  cRpos = putin2(i ,j,k) + 0.5*limiter5(rRpos,numdif)*(putin2(i,j,kp) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - wwR*cRpos*dz_i 
			ELSE
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
			  cRneg = putin2(i,j,kp) + 0.5*limiter5(rRneg,numdif)*(putin2(i,j, k) - putin2(i,j,kp))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - wwR*cRneg*dz_i 
			ENDIF
			IF (wwL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
			  cLpos = putin2(i,j,km) + 0.5*limiter5(rLpos,numdif)*(putin2(i,j, k) - putin2(i,j,km))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + wwL*cLpos*dz_i 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
			  cLneg = putin2(i,j,k ) + 0.5*limiter5(rLneg,numdif)*(putin2(i,j,km) - putin2(i,j, k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + wwL*cLneg*dz_i 
			ENDIF
			putout(i,j,k) = putout(i,j,k) -
     4   (rho(i,j,k)+rho(i,jp,k))*(Uvel(im,j,k)+Uvel(i,j,k)+Uvel(im,jp,k)+Uvel(i,jp,k))*Vvel(i,j,k)
     4   / ( 8.*Rp(i) )

           enddo
        enddo
      enddo
      return
      end	  



      subroutine advecw_C2Blend(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt,periodicx,periodicy,numdif)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phip(0:j1),phipt(0:je*px+1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter 
      real rho(0:i1,0:j1,0:k1),cfll,limiter5,dt,numdif
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,Rp2(-1:i1+1)
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kmm,kpp,imm,ipp,jmm,jpp,rank,periodicx,periodicy
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-1:i1+1,-1:j1+1,0:k1)
      real 	dz_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	  	  real Ax,Bx,DDx,blend,cLho,cRho
	  real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1)

	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim	

	putin2(0:i1,0:j1,0:k1)=Wvel
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=Wvel(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Wvel(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=Wvel(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Wvel(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
	  call shiftf2(Wvel,pbf)
	  call shiftb2(Wvel,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = Wvel(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =Wvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif
	
		dz_i=1./dz
      do k=kb,ke
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
        do j=jb,je
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
          do  i=ib,ie
			  ip=i+1
			  im=i-1
			  imm=im
			  ipp=ip
			  if (periodicx.eq.0.or.periodicx.eq.2) then
				if (i.eq.1) imm=i
				if (i.eq.ie) ipp=i
			  endif
			  varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
			  varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
			  varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
			  varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))	
 
      	    uuR=0.5*(Uvel(i,j,k)*rhU(i,j,k)+Uvel(i,j,kp)*rhU(i,j,kp))*Ru(i)
      	    uuL=0.5*(Uvel(im,j,k)*rhU(im,j,k)+Uvel(im,j,kp)*rhU(im,j,kp))*Ru(im)
			vvR=0.5*(Vvel(i,j,k)*rhV(i,j,k)+Vvel(i,j,kp)*rhV(i,j,kp))
    	    vvL=0.5*(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(i,jm,kp)*rhV(i,jm,kp))
    	    wwR=0.5*(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,j,kp)*rhW(i,j,kp))
			wwL=0.5*(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,j,km)*rhW(i,j,km))			

            putout(i,j,k) = 0.0

			IF (uuR.ge.0.) THEN
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter5(rRpos,numdif)*(putin2(ip,j,k) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = - uuR*cRpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
			  cRneg = putin2(ip,j,k) + 0.5*limiter5(rRneg,numdif)*(putin2( i,j,k) - putin2(ip,j,k))!*(1.-cfll)
			  putout(i,j,k) = - uuR*cRneg/( Rp(i)* dr(i) )
			ENDIF	
			IF (uuL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
			  cLpos = putin2(im,j,k) + 0.5*limiter5(rLpos,numdif)*(putin2( i,j,k) - putin2(im,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + uuL*cLpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
			  cLneg = putin2(i ,j,k) + 0.5*limiter5(rLneg,numdif)*(putin2(im,j,k) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + uuL*cLneg/ ( Rp(i)* dr(i) )
			ENDIF	
			IF (vvR.ge.0.) THEN
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
			  cRpos = putin2(i ,j,k) + 0.5*limiter5(rRpos,numdif)*(putin2(i,jp,k) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - vvR*cRpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
			  cRneg = putin2(i,jp,k) + 0.5*limiter5(rRneg,numdif)*(putin2( i,j,k) - putin2(i,jp,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - vvR*cRneg/ ( Rp(i) * (phiv(j)-phiv(jm)) )
			ENDIF		
			IF (vvL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
			  cLpos = putin2(i,jm,k) + 0.5*limiter5(rLpos,numdif)*(putin2( i,j,k) - putin2(i,jm,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + vvL*cLpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
			  cLneg = putin2(i,j ,k) + 0.5*limiter5(rLneg,numdif)*(putin2(i,jm,k) - putin2(i,j ,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + vvL*cLneg/ ( Rp(i) * (phiv(j)-phiv(jm)) )
			ENDIF	
			IF (wwR.ge.0.) THEN
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
			  cRpos = putin2(i ,j,k) + 0.5*limiter5(rRpos,numdif)*(putin2(i,j,kp) - putin2(i ,j,k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - wwR*cRpos*dz_i 
			ELSE
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
			  cRneg = putin2(i,j,kp) + 0.5*limiter5(rRneg,numdif)*(putin2(i,j, k) - putin2(i,j,kp))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) - wwR*cRneg*dz_i 
			ENDIF
			IF (wwL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
			  cLpos = putin2(i,j,km) + 0.5*limiter5(rLpos,numdif)*(putin2(i,j, k) - putin2(i,j,km))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + wwL*cLpos*dz_i 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
			  cLneg = putin2(i,j,k ) + 0.5*limiter5(rLneg,numdif)*(putin2(i,j,km) - putin2(i,j, k))!*(1.-cfll)
			  putout(i,j,k) = putout(i,j,k) + wwL*cLneg*dz_i 
			ENDIF
           enddo
        enddo
      enddo
      return
      end




      subroutine advecu_C4Blend(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt,periodicx,periodicy,numdif)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phip(0:j1),phipt(0:je*px+1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter 
      real rho(0:i1,0:j1,0:k1),cfll,limiter3,dt,numdif
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,Ru2(-1:i1+1)
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kmm,kpp,imm,ipp,jmm,jpp,rank,periodicx,periodicy
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-1:i1+1,-1:j1+1,0:k1)
      real 	dz_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)	
	  real Ax,Bx,DDx,blend,cLho,cRho,limiter5
	  real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1)

	Ru2(0:i1)=Ru
	Ru2(-1)=Ru(0)-(Ru(2)-Ru(1)) !needed for periodicx sim
 	Ru2(i1+1)=Ru(i1)+(Ru(i1)-Ru(i1-1)) !needed for periodicx sim	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim	

	putin2(0:i1,0:j1,0:k1)=Uvel
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=Uvel(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Uvel(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=Uvel(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Uvel(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
	  call shiftf2(Uvel,pbf)
	  call shiftb2(Uvel,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = Uvel(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =Uvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif
	
		  Ax=7.   !CDS4-own minimal dispersion phase speed error:6  !CDS4-Wicker: 7 !CDS4-Morinishi: 9
		  Bx=-1.  !CDS4-own minimal dispersion phase speed error:-1 !CDS4-Wicker:-1 !CDS4-Morinishi:-1
		  DDx=0.5/(Ax+Bx)	
		dz_i=1./dz
      do k=kb,ke
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
        do j=jb,je
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
          do  i=ib,ie
			  ip=i+1
			  im=i-1
			  imm=im
			  ipp=ip
			  if (periodicx.eq.0.or.periodicx.eq.2) then
				if (i.eq.1) imm=i
				if (i.eq.ie) ipp=i
			  endif
			  varx_grid_Rpos=(Ru2(ip)-Ru2(i))/(Ru2(i)-Ru2(im))
			  varx_grid_Rneg=(Ru2(ipp)-Ru2(i))/(Ru2(ipp+1)-Ru2(ipp))
			  varx_grid_Lpos=(Ru2(i)-Ru2(imm))/(Ru2(imm)-Ru2(imm-1))
			  varx_grid_Lneg=(Ru2(i)-Ru2(im))/(Ru2(ip)-Ru2(i))
			
            uuR=0.5*(Uvel(i,j,k)*rhU(i,j,k)+Uvel(ip,j,k)*rhU(ip,j,k))*Rp(ip)
            uuL=0.5*(Uvel(i,j,k)*rhU(i,j,k)+Uvel(im,j,k)*rhU(im,j,k))*Rp(i)
            vvR=0.5*(Vvel(i,j,k)*rhV(i,j,k)+Vvel(ip,j,k)*rhV(ip,j,k))
            vvL=0.5*(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(ip,jm,k)*rhV(ip,jm,k))
            wwR=0.5*(Wvel(i,j,k)*rhW(i,j,k)+Wvel(ip,j,k)*rhW(ip,j,k))
            wwL=0.5*(Wvel(i,j,km)*rhW(i,j,km)+Wvel(ip,j,km)*rhW(ip,j,km))			
			
            putout(i,j,k) = 0.0
			cRho = DDx*(Ax*(putin2(i,j,k)+putin2(ip,j,k))+Bx*(putin2(im,j,k)+putin2(ip+1,j,k)))
			cLho = DDx*(Ax*(putin2(i,j,k)+putin2(im,j,k))+Bx*(putin2(im-1,j,k)+putin2(ip,j,k)))
			IF (uuR.ge.0.) THEN
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Ru(ip)-Ru(i ))
			  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
			  blend = limiter5(rRpos,numdif)
			  cRpos = (1.-blend)*putin2(i ,j,k) + blend * cRho 
			  putout(i,j,k) = - uuR*cRpos/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ELSE
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Ru(ip)-Ru(i ))
			  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
			  blend = limiter5(rRneg,numdif)
			  cRneg = (1.-blend)*putin2(ip,j,k) + blend * cRho
			  putout(i,j,k) = - uuR*cRneg/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ENDIF	
			IF (uuL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Ru(i)-Ru(im))
			  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
			  blend = limiter5(rLpos,numdif)
			  cLpos = (1.-blend)*putin2(im,j,k) + blend * cLho
			  putout(i,j,k) = putout(i,j,k) + uuL*cLpos/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ELSE
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Ru(i)-Ru(im))
			  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
			  blend = limiter5(rLneg,numdif)
			  cLneg = (1.-blend)*putin2(im,j,k) + blend * cLho
			  putout(i,j,k) = putout(i,j,k) + uuL*cLneg/ ( Ru(i) * ( Rp(ip)-Rp(i) ) )
			ENDIF	
			cRho = DDx*(Ax*(putin2(i,j,k)+putin2(i,jp,k))+Bx*(putin2(i,jm,k)+putin2(i,jp+1,k)))
			cLho = DDx*(Ax*(putin2(i,j,k)+putin2(i,jm,k))+Bx*(putin2(i,jm-1,k)+putin2(i,jp,k)))			
			IF (vvR.ge.0.) THEN
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
			  blend = limiter5(rRpos,numdif)
			  cRpos = (1.-blend)*putin2(i ,j,k) + blend * cRho
			  putout(i,j,k) = putout(i,j,k) - vvR*cRpos/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
			  blend = limiter5(rRneg,numdif)
			  cRneg = (1.-blend)*putin2(i,jp,k) + blend * cRho
			  putout(i,j,k) = putout(i,j,k) - vvR*cRneg/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF		
			IF (vvL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
			  blend = limiter5(rLpos,numdif)
			  cLpos = (1.-blend)*putin2(i,jm,k) + blend * cLho
			  putout(i,j,k) = putout(i,j,k) + vvL*cLpos/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
			  blend = limiter5(rLneg,numdif)
			  cLneg = (1.-blend)*putin2(i,j ,k) + blend * cLho
			  putout(i,j,k) = putout(i,j,k) + vvL*cLneg/ ( Ru(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF	
			cRho = DDx*(Ax*(putin2(i,j,k)+putin2(i,j,kp))+Bx*(putin2(i,j,km)+putin2(i,j,kp+1)))
			cLho = DDx*(Ax*(putin2(i,j,k)+putin2(i,j,km))+Bx*(putin2(i,j,km-1)+putin2(i,j,kp)))				
			IF (wwR.ge.0.) THEN
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
			  blend = limiter5(rRpos,numdif)
			  cRpos = (1.-blend)*putin2(i ,j,k) + blend * cRho
			  putout(i,j,k) = putout(i,j,k) - wwR*cRpos*dz_i 
			ELSE
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
			  blend = limiter5(rRneg,numdif)
			  cRneg = (1.-blend)*putin2(i,j,kp) + blend * cRho
			  putout(i,j,k) = putout(i,j,k) - wwR*cRneg*dz_i 
			ENDIF
			IF (wwL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
			  blend = limiter5(rLpos,numdif)
			  cLpos = (1.-blend)*putin2(i,j,km) + blend * cLho
			  putout(i,j,k) = putout(i,j,k) + wwL*cLpos*dz_i 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
			  blend = limiter5(rLneg,numdif)
			  cLneg = (1.-blend)*putin2(i,j,k ) + blend * cLho
			  putout(i,j,k) = putout(i,j,k) + wwL*cLneg*dz_i 
			ENDIF
			putout(i,j,k) = putout(i,j,k) 
     & + (rho(i,j,k)+rho(ip,j,k))*(Vvel(i,j,k) + Vvel(ip,j,k) + Vvel(i,jm,k) + Vvel(ip,jm,k))**2/(32.0*Ru(i))
           enddo
        enddo
      enddo
      return
      end


      subroutine advecv_C4blend(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt,periodicx,periodicy,numdif)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phip(0:j1),phipt(0:je*px+1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),	 
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter 
      real rho(0:i1,0:j1,0:k1),cfll,limiter3,dt,numdif
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,Rp2(-1:i1+1)
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kmm,kpp,imm,ipp,jmm,jpp,rank,periodicx,periodicy
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-1:i1+1,-1:j1+1,0:k1)
      real 	dz_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phivt2(-1:je*px+2)	
	  real Ax,Bx,DDx,blend,cLho,cRho,limiter5
		real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1)

	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim	
	phivt2(0:je*px+1)=phivt
	phivt2(-1)=phivt(0)-(phivt(1)-phivt(0)) !needed for periodicy sim
 	phivt2(je*px+1+1)=phivt(je*px+1)+(phivt(je*px+1)-phivt(je*px+1-1)) !needed for periodicy sim	

	putin2(0:i1,0:j1,0:k1)=Vvel
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=Vvel(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Vvel(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=Vvel(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Vvel(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
	  call shiftf2(Vvel,pbf)
	  call shiftb2(Vvel,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = Vvel(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =Vvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif
	
		  Ax=7.   !CDS4-own minimal dispersion phase speed error:6  !CDS4-Wicker: 7 !CDS4-Morinishi: 9
		  Bx=-1.  !CDS4-own minimal dispersion phase speed error:-1 !CDS4-Wicker:-1 !CDS4-Morinishi:-1
		  DDx=0.5/(Ax+Bx)	
		dz_i=1./dz
      do k=kb,ke
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
        do j=jb,je
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif
      vary_grid_Rpos=(phivt2(rank*je+jp) -phivt2(rank*je+j))  /(phivt2(rank*je+j)    -phivt2(rank*je+jm))
      vary_grid_Rneg=(phivt2(rank*je+jpp)-phivt2(rank*je+j))  /(phivt2(rank*je+jpp+1)-phivt2(rank*je+jpp))
      vary_grid_Lpos=(phivt2(rank*je+j)  -phivt2(rank*je+jmm))/(phivt2(rank*je+jmm)  -phivt2(rank*je+jmm-1))
      vary_grid_Lneg=(phivt2(rank*je+j)  -phivt2(rank*je+jm)) /(phivt2(rank*je+jp)   -phivt2(rank*je+j))	
	  
	  
          do  i=ib,ie
			  ip=i+1
			  im=i-1
			  imm=im
			  ipp=ip
			  if (periodicx.eq.0.or.periodicx.eq.2) then
				if (i.eq.1) imm=i
				if (i.eq.ie) ipp=i
			  endif
			  varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
			  varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
			  varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
			  varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))			  
			  
           uuR=0.5*(Uvel(i,j,k)*rhU(i,j,k)+Uvel(i,jp,k)*rhU(i,jp,k))*Ru(i) !*Ru(i)
           uuL=0.5*(Uvel(im,j,k)*rhU(im,j,k)+Uvel(im,jp,k)*rhU(im,jp,k))*Ru(im) !*Ru(im)
     	   vvR=0.5*(Vvel(i,j,k)*rhV(i,j,k)+Vvel(i,jp,k)*rhV(i,jp,k))
           vvL=0.5*(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(i,j,k)*rhV(i,j,k))
           wwR=0.5*(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,jp,k)*rhW(i,jp,k))
           wwL=0.5*(Wvel(i,j,km)*rhW(i,j,km)+Wvel(i,jp,km)*rhW(i,jp,km))			
			
			putout(i,j,k) = 0.0
			cRho = DDx*(Ax*(putin2(i,j,k)+putin2(ip,j,k))+Bx*(putin2(im,j,k)+putin2(ip+1,j,k)))
			cLho = DDx*(Ax*(putin2(i,j,k)+putin2(im,j,k))+Bx*(putin2(im-1,j,k)+putin2(ip,j,k)))			
			IF (uuR.ge.0.) THEN
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
			  blend = limiter5(rRpos,numdif)
			  cRpos = (1.-blend)*putin2(i ,j,k) + blend*cRho 
			  putout(i,j,k) = - uuR*cRpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
			  blend = limiter5(rRneg,numdif)
			  cRneg = (1.-blend)*putin2(ip,j,k) + blend*cRho 
			  putout(i,j,k) = - uuR*cRneg/ ( Rp(i)* dr(i) )
			ENDIF	
			IF (uuL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
			  blend = limiter5(rLpos,numdif)
			  cLpos = (1.-blend)*putin2(im,j,k) + blend*cLho 
			  putout(i,j,k) = putout(i,j,k) + uuL*cLpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
			  blend = limiter5(rLneg,numdif)
			  cLneg = (1.-blend)*putin2(i ,j,k) + blend*cLho 
			  putout(i,j,k) = putout(i,j,k) + uuL*cLneg/ ( Rp(i)* dr(i) )
			ENDIF	
			cRho = DDx*(Ax*(putin2(i,j,k)+putin2(i,jp,k))+Bx*(putin2(i,jm,k)+putin2(i,jp+1,k)))
			cLho = DDx*(Ax*(putin2(i,j,k)+putin2(i,jm,k))+Bx*(putin2(i,jm-1,k)+putin2(i,jp,k)))			
			IF (vvR.ge.0.) THEN
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Rp(i)*(phivt2(rank*je+jp)-phivt2(rank*je+j)))
			  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
			  blend = limiter5(rRpos,numdif)
			  cRpos = (1.-blend)*putin2(i ,j,k) + blend*cRho 
			  putout(i,j,k) = putout(i,j,k) - vvR*cRpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Rp(i)*(phivt2(rank*je+jp)-phivt2(rank*je+j)))
			  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
			  blend = limiter5(rRneg,numdif)
			  cRneg = (1.-blend)*putin2(i,jp,k) + blend*cRho 
			  putout(i,j,k) = putout(i,j,k) - vvR*cRneg/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF		
			IF (vvL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Rp(i)*(phivt2(rank*je+j)-phivt2(rank*je+jm)))
			  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
			  blend = limiter5(rLpos,numdif)
			  cLpos = (1.-blend)*putin2(i,jm,k) + blend*cLho 
			  putout(i,j,k) = putout(i,j,k) + vvL*cLpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Rp(i)*(phivt2(rank*je+j)-phivt2(rank*je+jm)))
			  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
			  blend = limiter5(rLneg,numdif)
			  cLneg = (1.-blend)*putin2(i,j ,k) + blend*cLho 
			  putout(i,j,k) = putout(i,j,k) + vvL*cLneg/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ENDIF
			cRho = DDx*(Ax*(putin2(i,j,k)+putin2(i,j,kp))+Bx*(putin2(i,j,km)+putin2(i,j,kp+1)))
			cLho = DDx*(Ax*(putin2(i,j,k)+putin2(i,j,km))+Bx*(putin2(i,j,km-1)+putin2(i,j,kp)))			
			IF (wwR.ge.0.) THEN
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
			  blend = limiter5(rRpos,numdif)
			  cRpos = (1.-blend)*putin2(i ,j,k) + blend*cRho 
			  putout(i,j,k) = putout(i,j,k) - wwR*cRpos*dz_i 
			ELSE
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
			  blend = limiter5(rRneg,numdif)
			  cRneg = (1.-blend)*putin2(i,j,kp) + blend*cRho 
			  putout(i,j,k) = putout(i,j,k) - wwR*cRneg*dz_i 
			ENDIF
			IF (wwL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
			  blend = limiter5(rLpos,numdif)
			  cLpos = (1.-blend)*putin2(i,j,km) + blend*cLho 
			  putout(i,j,k) = putout(i,j,k) + wwL*cLpos*dz_i 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
			  blend = limiter5(rLneg,numdif)
			  cLneg = (1.-blend)*putin2(i,j,k ) + blend*cLho 
			  putout(i,j,k) = putout(i,j,k) + wwL*cLneg*dz_i 
			ENDIF
			putout(i,j,k) = putout(i,j,k) -
     4   (rho(i,j,k)+rho(i,jp,k))*(Uvel(im,j,k)+Uvel(i,j,k)+Uvel(im,jp,k)+Uvel(i,jp,k))*Vvel(i,j,k)
     4   / ( 8.*Rp(i) )

           enddo
        enddo
      enddo
      return
      end	  



      subroutine advecw_C4Blend(putout,Uvel,Vvel,Wvel,RHO,rhU,rhV,rhW,Ru,Rp,dr,phip,phiv,phipt,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,dt,periodicx,periodicy,numdif)
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
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phip(0:j1),phipt(0:je*px+1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter 
      real rho(0:i1,0:j1,0:k1),cfll,limiter3,dt,numdif
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,Rp2(-1:i1+1)
	real uuR,uuL,vvR,vvL,wwR,wwL
	integer kmm,kpp,imm,ipp,jmm,jpp,rank,periodicx,periodicy
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-1:i1+1,-1:j1+1,0:k1)
      real 	dz_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)	
	  real Ax,Bx,DDx,blend,cLho,cRho,limiter5 
	  real rhU(0:i1,0:j1,0:k1),rhV(0:i1,0:j1,0:k1),rhW(0:i1,0:j1,0:k1)

	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim	

	putin2(0:i1,0:j1,0:k1)=Wvel
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=Wvel(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Wvel(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=Wvel(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=Wvel(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
	  call shiftf2(Wvel,pbf)
	  call shiftb2(Wvel,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = Wvel(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =Wvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif
		  Ax=7.   !CDS4-own minimal dispersion phase speed error:6  !CDS4-Wicker: 7 !CDS4-Morinishi: 9
		  Bx=-1.  !CDS4-own minimal dispersion phase speed error:-1 !CDS4-Wicker:-1 !CDS4-Morinishi:-1
		  DDx=0.5/(Ax+Bx)	
		dz_i=1./dz
      do k=kb,ke
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
        do j=jb,je
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
          do  i=ib,ie
			  ip=i+1
			  im=i-1
			  imm=im
			  ipp=ip
			  if (periodicx.eq.0.or.periodicx.eq.2) then
				if (i.eq.1) imm=i
				if (i.eq.ie) ipp=i
			  endif
			  varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
			  varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
			  varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
			  varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))	
 
      	    uuR=0.5*(Uvel(i,j,k)*rhU(i,j,k)+Uvel(i,j,kp)*rhU(i,j,kp))*Ru(i)
      	    uuL=0.5*(Uvel(im,j,k)*rhU(im,j,k)+Uvel(im,j,kp)*rhU(im,j,kp))*Ru(im)
			vvR=0.5*(Vvel(i,j,k)*rhV(i,j,k)+Vvel(i,j,kp)*rhV(i,j,kp))
    	    vvL=0.5*(Vvel(i,jm,k)*rhV(i,jm,k)+Vvel(i,jm,kp)*rhV(i,jm,kp))
    	    wwR=0.5*(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,j,kp)*rhW(i,j,kp))
			wwL=0.5*(Wvel(i,j,k)*rhW(i,j,k)+Wvel(i,j,km)*rhW(i,j,km))			

            putout(i,j,k) = 0.0
			cRho = DDx*(Ax*(putin2(i,j,k)+putin2(ip,j,k))+Bx*(putin2(im,j,k)+putin2(ip+1,j,k)))
			cLho = DDx*(Ax*(putin2(i,j,k)+putin2(im,j,k))+Bx*(putin2(im-1,j,k)+putin2(ip,j,k)))
			IF (uuR.ge.0.) THEN
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
			  blend = limiter5(rRpos,numdif)
			  cRpos = putin2(i ,j,k)*(1.-blend) + blend * cRho 
			  putout(i,j,k) = - uuR*cRpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(ip,j,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuR)/(Rp(ip)-Rp(i ))
			  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
			  blend = limiter5(rRneg,numdif)
			  cRneg = putin2(ip,j,k)*(1.-blend) + blend * cRho 
			  putout(i,j,k) = - uuR*cRneg/( Rp(i)* dr(i) )
			ENDIF	
			IF (uuL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
			  blend = limiter5(rLpos,numdif)
			  cLpos = putin2(im,j,k)*(1.-blend) + blend * cLho 
			  putout(i,j,k) = putout(i,j,k) + uuL*cLpos/ ( Rp(i)* dr(i) )
			ELSE
			  noemer = putin2(i,j,k)-putin2(im,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(uuL)/(Rp(i)-Rp(im))
			  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
			  blend = limiter5(rLneg,numdif)
			  cLneg = putin2(i ,j,k)*(1.-blend) + blend * cLho 
			  putout(i,j,k) = putout(i,j,k) + uuL*cLneg/ ( Rp(i)* dr(i) )
			ENDIF	
			cRho = DDx*(Ax*(putin2(i,j,k)+putin2(i,jp,k))+Bx*(putin2(i,jm,k)+putin2(i,jp+1,k)))
			cLho = DDx*(Ax*(putin2(i,j,k)+putin2(i,jm,k))+Bx*(putin2(i,jm-1,k)+putin2(i,jp,k)))			
			IF (vvR.ge.0.) THEN
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
			  blend = limiter5(rRpos,numdif)
			  cRpos = putin2(i ,j,k)*(1.-blend) + blend * cRho 
			  putout(i,j,k) = putout(i,j,k) - vvR*cRpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,jp,k)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvR)/(Ru(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
			  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
			  blend = limiter5(rRneg,numdif)
			  cRneg = putin2(i,jp,k)*(1.-blend) + blend * cRho 
			  putout(i,j,k) = putout(i,j,k) - vvR*cRneg/ ( Rp(i) * (phiv(j)-phiv(jm)) )
			ENDIF		
			IF (vvL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
			  blend = limiter5(rLpos,numdif)
			  cLpos = putin2(i,jm,k)*(1.-blend) + blend * cLho 
			  putout(i,j,k) = putout(i,j,k) + vvL*cLpos/ ( Rp(i) * (phiv(j)-phiv(jm)) ) 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,jm,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(vvL)/(Ru(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
			  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
			  blend = limiter5(rLneg,numdif)
			  cLneg = putin2(i,j ,k)*(1.-blend) + blend * cLho 
			  putout(i,j,k) = putout(i,j,k) + vvL*cLneg/ ( Rp(i) * (phiv(j)-phiv(jm)) )
			ENDIF	
			cRho = DDx*(Ax*(putin2(i,j,k)+putin2(i,j,kp))+Bx*(putin2(i,j,km)+putin2(i,j,kp+1)))
			cLho = DDx*(Ax*(putin2(i,j,k)+putin2(i,j,km))+Bx*(putin2(i,j,km-1)+putin2(i,j,kp)))			
			IF (wwR.ge.0.) THEN
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
			  blend = limiter5(rRpos,numdif)
			  cRpos = putin2(i ,j,k)*(1.-blend) + blend * cRho 
			  putout(i,j,k) = putout(i,j,k) - wwR*cRpos*dz_i 
			ELSE
			  noemer = putin2(i,j,kp)-putin2(i,j,k)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwR)*dz_i
			  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
			  blend = limiter5(rRneg,numdif)
			  cRneg = putin2(i,j,kp)*(1.-blend) + blend * cRho 
			  putout(i,j,k) = putout(i,j,k) - wwR*cRneg*dz_i 
			ENDIF
			IF (wwL.ge.0.) THEN
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
			  blend = limiter5(rLpos,numdif)
			  cLpos = putin2(i,j,km)*(1.-blend) + blend * cLho 
			  putout(i,j,k) = putout(i,j,k) + wwL*cLpos*dz_i 
			ELSE
			  noemer = putin2(i,j,k)-putin2(i,j,km)
			  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
			  !cfll=dt*ABS(wwL)*dz_i
			  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
			  blend = limiter5(rLneg,numdif)
			  cLneg = putin2(i,j,k )*(1.-blend) + blend * cLho 
			  putout(i,j,k) = putout(i,j,k) + wwL*cLneg*dz_i 
			ENDIF
           enddo
        enddo
      enddo
      return
      end


	  

      real function limiter3(r,cfl)
      real r,cfl,alpha,slp
! limiter = third order accurate (with cfl=0 equal to Koren limiter) CFL dependent limiter from Arora and Roe 1997
	!alpha=(1.+cfl)/3.
	!limiter2=MAX(0.,MIN(MIN(2./MAX(cfl,1e-12)*r,1.+alpha*(r-1.)),1.99/(1.-cfl))) !looser TVD contstraints gives in TUDflow3d e-7 undershoot, therefore not advised...
	!limiter2=MAX(0.,MIN(MIN(1.99/MAX(cfl,1e-12)*r,1.+alpha*(r-1.)),2.)) !looser TVD contstraints gives in TUDflow3d e-7 undershoot, therefore not advised...
!	limiter2=MAX(0.,MIN(MIN(2.*r,1.+alpha*(r-1.)),2.)) !advised for non-linear systems a more restrictive TVD constraints, gives no undershoot in matlab and TUDflow3d but more numerical diff... 18-11-17 SWITCHED THIS ON
	
  
	
	!limiter2=(r+ABS(r))/MAX(1.+r,1.) 	! Van Leer limiter
	!limiter2=MAX(0.,MIN(2.*r,1.),MIN(r,2.))	! Superbee limiter
	!limiter3=MAX(0.,MIN(2.*r/MAX(cfl,1.e-6),1.),MIN(r,2./(1.-MIN(cfl,0.999999))))	! Superbee limiter loosened CFL-dependent restrictions 
	limiter3=1. !LW
	!limiter3=MIN(MAX(1.,r),2./(1.-cfl)) !LW or anti diffusion
	!limiter3=1./(1.-MIN(cfl,0.999999)) ! CDS2
	!limiter3=MAX(0.,r/(ABS(r)+1.e-6)*1./(1.-MIN(cfl,0.999999))) ! CDS2 if r>0 else UPW1 --> much too dissipative!	
	!limiter3=MAX(r/(1.-MIN(cfl,0.999999)),MIN(2.*r/MAX(cfl,1.e-6),1.),MIN(r,2./(1.-MIN(cfl,0.999999))))	
	!limiter3=MAX(r/(1.-MIN(cfl,0.999999)),1.,MIN(r,2./(1.-MIN(cfl,0.999999))))	
!	limiter3=MAX(0.,MIN(2.*r/MAX(cfl,1.e-6),1.),MIN(r,2./(1.-MIN(cfl,0.999999))))	! Superbee limiter loosened CFL-dependent restrictions
!	limiter3=limiter3/(1.-MIN(cfl,0.999999))
	
	!limiter3=MAX(MAX(r/(1.-MIN(cfl,0.00),MIN(2.*r/MAX(cfl,1.e-6),1.)),MIN(r,2./(1.-MIN(cfl,0.999999))))	! Superbee limiter loosened CFL-dependent restrictions 
	!in combi with EE1: always apply LW or move towards CDS2 or towards Downwind1
!	IF (r<-1.) THEN
		!limiter3=1. !LW
!		limiter3=MAX(0.1*r+1./(1.-MIN(cfl,0.999999)),1.)
		!limiter3=MAX((r+1.)/(1.-MIN(cfl,0.999999)),1.)
		!limiter3=(r+3.+MAX(1.+r,0.))/(4.*(1.-MIN(cfl,0.999999))) ! QUICK geblend naar CDS2 bij r=0
!		limiter3=(r+3.+MAX(3.+r,0.))/(4.*(1.-MIN(cfl,0.999999))) ! QUICK geblend naar CDS2 bij r=-1
		
		
		!limiter3=MAX((r+3.+MAX(1.+r,0.))/(4.*(1.-MIN(cfl,0.999999))),1.) ! QUICK geblend naar CDS2 bij r=0 met limit op LW
		!alpha=1. !alpha=0 is limit op CDS2 en alpha=1 limit op LW
		!limiter3=MAX((r+3.+MAX(1.+r,0.))/(4.*(1.-MIN(cfl,0.999999))),(1.-alpha*cfl)/(1.-MIN(cfl,0.999999))) ! QUICK geblend naar CDS2 bij r=0 met limit tussen CDS2 en LW in 
		!limiter3=r/(1.-MIN(cfl,0.999999))   !UPW2
		!limiter3=(r+3.)/(4.*(1.-MIN(cfl,0.999999)))   !QUICK
!	ELSE
!		limiter3=1./(1.-MIN(cfl,0.999999)) ! CDS2
		
		
		!limiter3=MAX(0.,MIN(2.*r/MAX(cfl,1.e-6),1.),MIN(r,2./(1.-MIN(cfl,0.999999))))	! Superbee limiter loosened CFL-dependent restrictions 
		!limiter3=MAX(0.,MIN(2.*r/MAX(cfl,1.e-6),1./(1.-MIN(cfl,0.999999))),MIN(r,2./(1.-MIN(cfl,0.999999))))	! Superbee limiter loosened CFL-dependent restrictions 
		!limiter3=MAX(0.,MIN(2.*r/MAX(cfl,1.e-6),1.),MIN(r,2./(1.-MIN(cfl,0.999999))))	! Superbee limiter loosened CFL-dependent restrictions
		!limiter3=MAX(0.75,MIN(2.*r/MAX(cfl,1.e-6),1.),MIN(r,2./(1.-MIN(cfl,0.999999))))	! Superbee limiter loosened CFL-dependent restrictions
		!limiter3=limiter3/(1.-MIN(cfl,0.999999))		
!	ENDIF

! with ABv the following is stable without much wiggles for jicf 	
!	IF (r<0.) THEN
!		!limiter3=(r+3.+MAX(1.+r,0.))/(4.*(1.-MIN(cfl,0.999999))) ! QUICK geblend naar CDS2 bij r=0
!		limiter3=MAX((r+1.)/(1.-MIN(cfl,0.999999)),1.) !blend van CDS2 naar LW toe
!	ELSE
!		limiter3=MAX(1./(1.-MIN(cfl,0.999999)),MIN(r,2./(1.-MIN(cfl,0.999999)))) ! CDS2 and upper limit of SuperbeeCFL
!	ENDIF
	
!	limiter3=1./(1.-MIN(cfl,0.999999))*(1.-cfl)**0.1 ! CDS2 with minor LW influence
! with ABv the following is stable without much wiggles for jicf 	
!	IF (r<-1.) THEN
!		limiter3=(r+3.+MAX(3.+r,0.))/(4.*(1.-MIN(cfl,0.999999))) ! QUICK geblend naar CDS2 bij r=-1
!	ELSE
!		limiter3=1./(1.-MIN(cfl,0.999999)) ! CDS2
!	ENDIF
      return
      end
	  
	  
      real function limiter5(r,numdif)
      real r,numdif

	!limiter5= MIN(MAX(1.+numdif*r,0.),1.) !own blend which is CDS2 for r>0 and slowly goes towards upw1 for r<0 (slope phi-r line is numdif)
	limiter5= MIN(MAX(1.+numdif*r,1.-numdif),1.)

      return
      end