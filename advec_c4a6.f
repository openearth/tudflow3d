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

      subroutine advecu_C4A6(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,phiv)
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
	integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy	  
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1),numdif,dzi,phiv(0:j1)
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,rhoipp,rhoimm,rhojpp,rhojmm,rhokpp,rhokmm
	real Axa,Bxa,DDxa,Axb,Bxb,DDxb,Axa2,Bxa2,DDxa2,Axb2,Bxb2,DDxb2,facA,facB
	real uuRA,uuLA,vvRA,vvLA,wwRA,wwLA,uuRB,uuLB,vvRB,vvLB,wwRB,wwLB
	real uuXRA,uuXLA,uuYRA,uuYLA,uuZRA,uuZLA,uuXRB,uuXLB,uuYRB,uuYLB,uuZRB,uuZLB,order4yes

	real ubb(0:i1,0:k1),ubf(0:i1,0:k1),Uvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real Vvel2(-2:i1+2,-2:j1+2,-2:k1+2),Wvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real rho2(-1:i1+1,-1:j1+1,-1:k1+1),phivt2(-1:je*px+2),Rp2(-1:i1+1)

	
	
	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(1)-Rp(0)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim	
	phivt2(0:je*px+1)=phivt
	phivt2(-1)=phivt(0)-(phivt(1)-phivt(0)) !needed for periodicy sim
 	phivt2(je*px+1+1)=phivt(je*px+1)+(phivt(je*px+1)-phivt(je*px+1-1)) !needed for periodicy sim
	
	Uvel2(0:i1,0:j1,0:k1)=Uvel

c get stuff from other CPU's
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
	!Uvel2(-2:i1+2,-2:j1+2,0) =0.
	Uvel2(-2:i1+2,-2:j1+2,-1)=  Uvel2(-2:i1+2,-2:j1+2,0)
	Uvel2(-2:i1+2,-2:j1+2,-2)=  Uvel2(-2:i1+2,-2:j1+2,0)
	!Uvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Uvel2(-2:i1+2,-2:j1+2,1)+1./3.*Uvel2(-2:i1+2,-2:j1+2,2)
	!Uvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Uvel2(-2:i1+2,-2:j1+2,1)+2.*Uvel2(-2:i1+2,-2:j1+2,2)	
	Uvel2(-2:i1+2,-2:j1+2,k1+1)=Uvel2(-2:i1+2,-2:j1+2,k1)
	Uvel2(-2:i1+2,-2:j1+2,k1+2)=Uvel2(-2:i1+2,-2:j1+2,k1)

	
	Vvel2(0:i1,0:j1,0:k1)=Vvel

c get stuff from other CPU's
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
	!Vvel2(-2:i1+2,-2:j1+2,0) =0.
	Vvel2(-2:i1+2,-2:j1+2,-1)=  Vvel2(-2:i1+2,-2:j1+2,0)
	Vvel2(-2:i1+2,-2:j1+2,-2)=  Vvel2(-2:i1+2,-2:j1+2,0)
	!Vvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Vvel2(-2:i1+2,-2:j1+2,1)+1./3.*Vvel2(-2:i1+2,-2:j1+2,2)
	!Vvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Vvel2(-2:i1+2,-2:j1+2,1)+2.*Vvel2(-2:i1+2,-2:j1+2,2)
	Vvel2(-2:i1+2,-2:j1+2,k1+1)=Vvel2(-2:i1+2,-2:j1+2,k1)
	Vvel2(-2:i1+2,-2:j1+2,k1+2)=Vvel2(-2:i1+2,-2:j1+2,k1)

		Wvel2(0:i1,0:j1,0:k1)=Wvel

c get stuff from other CPU's
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
!	Wvel2(-2:i1+2,-2:j1+2,-1)=  -Wvel2(-2:i1+2,-2:j1+2,1)
!	Wvel2(-2:i1+2,-2:j1+2,-2)=  -Wvel2(-2:i1+2,-2:j1+2,2)	
	Wvel2(-2:i1+2,-2:j1+2,k1+1)=Wvel2(-2:i1+2,-2:j1+2,k1)
	Wvel2(-2:i1+2,-2:j1+2,k1+2)=Wvel2(-2:i1+2,-2:j1+2,k1)

	rho2(0:i1,0:j1,0:k1)=rho

c get stuff from other CPU's
	  call shiftf2(rho,ubf)
	  call shiftb2(rho,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = rho(i,0,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =rho(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	
	if (periodicx.eq.0.or.periodicx.eq.2) then
		rho2(-1,-1:j1+1,0:k1)=  rho2(1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(ie,-1:j1+1,0:k1)
	else
		rho2(-1,-1:j1+1,0:k1)=  rho2(ie-1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(2,-1:j1+1,0:k1)
	endif
	rho2(-1:i1+1,-1:j1+1,-1)=  rho2(-1:i1+1,-1:j1+1,0)
	rho2(-1:i1+1,-1:j1+1,k1+1)=rho2(-1:i1+1,-1:j1+1,k1)
	
	dzi=1./dz 
      do k=kb,ke
      kp=k+1
      km=k-1
		  kpp=k+2 !MIN(k+2,k1)
		  kppp=k+3 !MIN(k+3,k1)
		  kmm=k-2 !MAX(k-2,0)
		  kmmm=k-3 !MAX(k-3,0)
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
	  
	      rhoim = rho(i ,j,k)
		  rhoip = rho(ip,j,k) 
          rhojp =0.25*(rho(i,j,k)+rho(i,jp,k)+rho(ip,j,k)+rho(ip,jp,k))
          rhojm =0.25*(rho(i,j,k)+rho(i,jm,k)+rho(ip,j,k)+rho(ip,jm,k))
          rhokp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
          rhokm =0.25*(rho(i,j,k)+rho(i,j,km)+rho(ip,j,k)+rho(ip,j,km))
		  
	      rhoimm = rho2(im ,j,k)
		  rhoipp = rho2(ipp,j,k) 
          rhojpp =0.25*(rho2(i,jp,k )+rho2(i,jpp,k  )+rho2(ip,jp,k )+rho2(ip,jpp,k  ))
          rhojmm =0.25*(rho2(i,jm,k )+rho2(i,jmm,k  )+rho2(ip,jm,k )+rho2(ip,jmm,k  ))
          rhokpp =0.25*(rho2(i,j ,kp)+rho2(i,j  ,kpp)+rho2(ip,j ,kp)+rho2(ip,j  ,kpp))
          rhokmm =0.25*(rho2(i,j ,km)+rho2(i,j  ,kmm)+rho2(ip,j ,km)+rho2(ip,j  ,kmm))
		  
      putout(i,j,k) = 0.0
!      uuR=(Uvel(i,j,k)+Uvel(ip,j,k))*rho(ip,j,k)*Rp(ip)
!      uuL=(Uvel(i,j,k)+Uvel(im,j,k))*rho(i,j,k)*Rp(i)
!      vvR=(Vvel(i,j,k)+Vvel(ip,j,k))*rhojp
!      vvL=(Vvel(i,jm,k)+Vvel(ip,jm,k))*rhojm
!      wwR=(Wvel(i,j,k)+Wvel(ip,j,k))*rhokp
!      wwL=(Wvel(i,j,km)+Wvel(ip,j,km))*rhokm
		   if ((periodicx.eq.1.and.periodicy.eq.1.and.(k.le.2.or.k.ge.ke-1))
     &	   .or.((periodicx.eq.1.and.periodicy.ne.1).and.
     &           (k.le.2.or.k.ge.ke-1.or.(rank.eq.0.and.j.le.2).or.(rank.eq.px-1.and.j.ge.je-1))) 
     &     .or.((periodicy.eq.1.and.periodicx.ne.1).and.(i.le.2.or.i.ge.ie-1.or.k.le.2.or.k.ge.ke-1))
     &     .or.((periodicy.ne.1.and.periodicx.ne.1).and.
     &           (i.le.2.or.i.ge.ie-1.or.k.le.2.or.k.ge.ke-1.or.(rank.eq.0.and.j.le.2).or.(rank.eq.px-1.and.j.ge.je-1)))) then  !switch to cds2
!			if (.false.) then
			  Axa=1.
			  Bxa=0.
			  DDxa=1./(Axa+Bxa)
			  Axb=1.
			  Bxb=0.
			  DDxb=1./(Axb+Bxb)	
			  Axa2=1.
			  Bxa2=0.
			  DDxa2=1./(Axa2+Bxa2)
			  Axb2=1.
			  Bxb2=0.
			  DDxb2=1./(Axb2+Bxb2)	
			  facA=1.
			  facB=0.
		  else
			  Axa=9. !1. !9.
			  Bxa=-1. !0. !-1.
			  DDxa=1./(Axa+Bxa)
			  Axb=0. !9. !0. !9.
			  Bxb=1. !-1. !1. !-1.
			  DDxb=1./(Axb+Bxb)	
			  Axa2=9. !9. !1. 
			  Bxa2=-1. !-1. !0. 
			  DDxa2=1./(Axa2+Bxa2)
			  Axb2=0. 
			  Bxb2=1. 
			  DDxb2=1./(Axb2+Bxb2)				  
			  facA=1. !1.+MIN(1.,DBLE(k-1)/3.)*1./8. !9./8.
			  facB=0. !MIN(1.,DBLE(k-1)/3.)*-1./8.    !-1./8.					  
			  order4yes=0.
		  endif	

	  ! advection velocities:
      uuRA=DDxa*(Axa*(Uvel2(i ,j ,k )+Uvel2(ip,j ,k ))+Bxa*(Uvel2(im ,j ,k )+Uvel2(ipp,j ,k )))
      uuLA=DDxa*(Axa*(Uvel2(im,j ,k )+Uvel2(i ,j ,k ))+Bxa*(Uvel2(imm,j ,k )+Uvel2(ip ,j ,k )))
	  vvRA=DDxa*(Axa*(Vvel2(i ,j ,k )+Vvel2(ip,j ,k ))+Bxa*(Vvel2(im ,j ,k )+Vvel2(ipp,j ,k )))
      vvLA=DDxa*(Axa*(Vvel2(i ,jm,k )+Vvel2(ip,jm,k ))+Bxa*(Vvel2(im ,jm,k )+Vvel2(ipp,jm,k )))
      wwRA=DDxa*(Axa*(Wvel2(i ,j ,k )+Wvel2(ip,j ,k ))+Bxa*(Wvel2(im ,j ,k )+Wvel2(ipp,j ,k )))
      wwLA=DDxa*(Axa*(Wvel2(i ,j ,km)+Wvel2(ip,j ,km))+Bxa*(Wvel2(im ,j ,km)+Wvel2(ipp,j ,km)))
!      wwRA=1*(1*(Wvel2(i ,j ,k )+Wvel2(ip,j ,k )))
!      wwLA=1*(1*(Wvel2(i ,j ,km)+Wvel2(ip,j ,km)))	  
      uuRB=DDxb*(Axb*(Uvel2(ip ,j  ,k  )+Uvel2(ipp,j  ,k  ))+Bxb*(Uvel2(i   ,j  ,k  )+Uvel2(ippp,j  ,k  )))
      uuLB=DDxb*(Axb*(Uvel2(imm,j  ,k  )+Uvel2(im ,j  ,k  ))+Bxb*(Uvel2(immm,j  ,k  )+Uvel2(i   ,j  ,k  )))
	  vvRB=DDxb*(Axb*(Vvel2(i  ,jp ,k  )+Vvel2(ip ,jp ,k  ))+Bxb*(Vvel2(im  ,jp ,k  )+Vvel2(ipp ,jp ,k  )))
      vvLB=DDxb*(Axb*(Vvel2(i  ,jmm,k  )+Vvel2(ip ,jmm,k  ))+Bxb*(Vvel2(im  ,jmm,k  )+Vvel2(ipp ,jmm,k  )))
      wwRB=DDxb*(Axb*(Wvel2(i  ,j  ,kp )+Wvel2(ip ,j  ,kp ))+Bxb*(Wvel2(im  ,j  ,kp )+Wvel2(ipp ,j  ,kp )))
      wwLB=DDxb*(Axb*(Wvel2(i  ,j  ,kmm)+Wvel2(ip ,j  ,kmm))+Bxb*(Wvel2(im  ,j  ,kmm)+Wvel2(ipp ,j  ,kmm)))
	  ! advected velocities:
      uuXRA=DDxa2*(Axa2*(Uvel2(i ,j ,k )+Uvel2(ip,j ,k ))+Bxa2*(Uvel2(im ,j  ,k )+Uvel2(ipp,j  ,k  )))
      uuXLA=DDxa2*(Axa2*(Uvel2(im,j ,k )+Uvel2(i ,j ,k ))+Bxa2*(Uvel2(imm,j  ,k )+Uvel2(ip ,j  ,k  )))
	  uuYRA=DDxa2*(Axa2*(Uvel2(i ,j ,k )+Uvel2(i ,jp,k ))+Bxa2*(Uvel2(i  ,jm ,k )+Uvel2(i  ,jpp,k  )))
      uuYLA=DDxa2*(Axa2*(Uvel2(i ,jm,k )+Uvel2(i ,j ,k ))+Bxa2*(Uvel2(i  ,jmm,k )+Uvel2(i  ,jp ,k  )))
      uuZRA=DDxa2*(Axa2*(Uvel2(i ,j ,k )+Uvel2(i ,j ,kp))+Bxa2*(Uvel2(i  ,j  ,km )+Uvel2(i ,j  ,kpp)))
      uuZLA=DDxa2*(Axa2*(Uvel2(i ,j ,km)+Uvel2(i ,j ,k ))+Bxa2*(Uvel2(i  ,j  ,kmm)+Uvel2(i ,j  ,kp )))
      uuXRB=DDxb2*(Axb2*(Uvel2(ip  ,j   ,k   )+Uvel2(ipp ,j   ,k   ))+Bxb2*(Uvel2(i   ,j   ,k   )+Uvel2(ippp,j   ,k   )))
      uuXLB=DDxb2*(Axb2*(Uvel2(imm ,j   ,k   )+Uvel2(im  ,j   ,k   ))+Bxb2*(Uvel2(immm,j   ,k   )+Uvel2(i   ,j   ,k   )))
	  uuYRB=DDxb2*(Axb2*(Uvel2(i   ,jp  ,k   )+Uvel2(i   ,jpp ,k   ))+Bxb2*(Uvel2(i   ,j   ,k   )+Uvel2(i   ,jppp,k   )))
      uuYLB=DDxb2*(Axb2*(Uvel2(i   ,jmm ,k   )+Uvel2(i   ,jm  ,k   ))+Bxb2*(Uvel2(i   ,jmmm,k   )+Uvel2(i   ,j   ,k   )))
      uuZRB=DDxb2*(Axb2*(Uvel2(i   ,j   ,kp  )+Uvel2(i   ,j   ,kpp ))+Bxb2*(Uvel2(i   ,j   ,k   )+Uvel2(i   ,j   ,kppp)))
      uuZLB=DDxb2*(Axb2*(Uvel2(i   ,j   ,kmm )+Uvel2(i   ,j   ,km  ))+Bxb2*(Uvel2(i   ,j   ,kmmm)+Uvel2(i   ,j   ,k   )))
	  
!	  IF (ABS(Uvel2(i,j,k)-Uvel(i,j,k)).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,U,U2',rank,i,j,k,Uvel2(i,j,k),Uvel(i,j,k)
!	  ENDIF
!	  IF (ABS(Vvel2(i,j,k)-Vvel(i,j,k)).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,V,V2',rank,i,j,k,Vvel2(i,j,k),Vvel(i,j,k)
!	  ENDIF
!	  IF (ABS(Wvel2(i,j,k)-Wvel(i,j,k)).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,W,W2',rank,i,j,k,Wvel2(i,j,k),Wvel(i,j,k)
!	  ENDIF
!	  IF (ABS(Uvel2(im,jm,km)-Uvel(im,jm,km)).gt.1e-14) THEN
!		write(*,*),'rank,im,jm,km,U,U2',rank,im,jm,km,Uvel2(im,jm,km),Uvel(im,jm,km)
!	  ENDIF
!	  IF (ABS(Vvel2(im,jm,km)-Vvel(im,jm,km)).gt.1e-14) THEN
!		write(*,*),'rank,im,jm,km,V,V2',rank,im,jm,km,Vvel2(im,jm,km),Vvel(im,jm,km)
!	  ENDIF
!	  IF (ABS(Wvel2(im,jm,km)-Wvel(im,jm,km)).gt.1e-14) THEN
!		write(*,*),'rank,im,jm,km,W,W2',rank,im,jm,km,Wvel2(im,jm,km),Wvel(im,jm,km)
!	  ENDIF
!	  IF (ABS(Uvel2(ip,jp,kp)-Uvel(ip,jp,kp)).gt.1e-14) THEN
!		write(*,*),'rank,ip,jp,kp,U,U2',rank,ip,jp,kp,Uvel2(ip,jp,kp),Uvel(ip,jp,kp)
!	  ENDIF
!	  IF (ABS(Vvel2(ip,jp,kp)-Vvel(ip,jp,kp)).gt.1e-14) THEN
!		write(*,*),'rank,ip,jp,kp,V,V2',rank,ip,jp,kp,Vvel2(ip,jp,kp),Vvel(ip,jp,kp)
!	  ENDIF
!	  IF (ABS(Wvel2(ip,jp,kp)-Wvel(ip,jp,kp)).gt.1e-14) THEN
!		write(*,*),'rank,ip,jp,kp,W,W2',rank,ip,jp,kp,Wvel2(ip,jp,kp),Wvel(ip,jp,kp)
!	  ENDIF	  
!	  uuRB=DDxa*(Axa*(Uvel(i ,j ,k )+Uvel(ip,j ,k ))+Bxa*(Uvel2(im ,j ,k )+Uvel2(ipp,j ,k )))
!      uuLB=DDxa*(Axa*(Uvel(im,j ,k )+Uvel(i ,j ,k ))+Bxa*(Uvel2(imm,j ,k )+Uvel2(ip ,j ,k )))
!	  vvRB=DDxa*(Axa*(Vvel(i ,j ,k )+Vvel(ip,j ,k ))+Bxa*(Vvel2(im ,j ,k )+Vvel2(ipp,j ,k )))
!      vvLB=DDxa*(Axa*(Vvel(i ,jm,k )+Vvel(ip,jm,k ))+Bxa*(Vvel2(im ,jm,k )+Vvel2(ipp,jm,k )))
!      wwRB=DDxa*(Axa*(Wvel(i ,j ,k )+Wvel(ip,j ,k ))+Bxa*(Wvel2(im ,j ,k )+Wvel2(ipp,j ,k )))
!      wwLB=DDxa*(Axa*(Wvel(i ,j ,km)+Wvel(ip,j ,km))+Bxa*(Wvel2(im ,j ,km)+Wvel2(ipp,j ,km)))
!	  IF (ABS(uuRB-uuRA).gt.1e-14.or.ABS(uuLB-uuLA).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,uU',rank,i,j,k,uuRA,uuRB,uuLB,uuLA 
!	  ENDIF
!	  IF (ABS(vvRB-vvRA).gt.1e-14.or.ABS(vvLB-vvLA).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,uV',rank,i,j,k,vvRA,vvRB,vvLB,vvLA 
!	  ENDIF
!	  IF (ABS(wwRB-wwRA).gt.1e-14.or.ABS(wwLB-wwLA).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,uW',rank,i,j,k,wwRA,wwRB,wwLB,wwLA 
!	  ENDIF
	  
	  
      putout(i,j,k) = - 0.25 * ( 
     1 facA*(rhoip*uuRA*Rp(ip)*uuXRA - numdif/facA*ABS(rhoip*uuRA*Rp(ip)) * 
     1 (10.*(-Uvel2(i,j,k)+Uvel2(ip,j,k))-5.*(-Uvel2(im,j,k)+Uvel2(ipp,j,k))+(-Uvel2(imm,j,k)+Uvel2(ippp,j,k))) -
     1 (rhoim*uuLA*Rp(i)*uuXLA - numdif/facA*ABS(rhoim*uuLA*Rp(i)) * 
     1 (10.*(-Uvel2(im,j,k)+Uvel2(i,j,k))-5.*(-Uvel2(imm,j,k)+Uvel2(ip,j,k))+(-Uvel2(immm,j,k)+Uvel2(ipp,j,k)))) )
     1  / ( Ru(i) * ( Rp(ip)-Rp(i) ) )
     +                         + 
     1 facB*(rhoipp*uuRB*Rp2(ipp)*uuXRB - rhoimm*uuLB*Rp2(im)*uuXLB)
     1  / ( Ru(i) * ( Rp2(ipp)-Rp2(im) ) )
     +                         + 	 
     2 facA*(rhojp*vvRA*uuYRA - numdif/facA*
     2 ABS(rhojp*vvRA)* (10.*(-Uvel2(i,j,k)+Uvel2(i,jp,k))-5.*(-Uvel2(i,jm,k)+Uvel2(i,jpp,k))+(-Uvel2(i,jmm,k)+Uvel2(i,jppp,k))) -
     2 (rhojm*vvLA*uuYLA - numdif/facA*
     2 ABS(rhojm*vvLA)* (10.*(-Uvel2(i,jm,k)+Uvel2(i,j,k))-5.*(-Uvel2(i,jmm,k)+Uvel2(i,jp,k))+(-Uvel2(i,jmmm,k)+Uvel2(i,jpp,k)))) )
     2  / ( Ru(i) * (phiv(j)-phiv(jm)) ) !/ ( Ru(i) * (phivt2(rank*je+j)-phivt2(rank*je+jm)) )
     +                         + 	 
     2 facB*(rhojpp*vvRB*uuYRB - rhojmm*vvLB*uuYLB)   
     2  / ( Ru(i) * (phivt2(rank*je+jp)-phivt2(rank*je+jmm)) )	 
     +                         +
     3 facA*(rhokp*wwRA*uuZRA - numdif/facA*	 
     3 ABS(rhokp*wwRA)* (10.*(-Uvel2(i,j,k)+Uvel2(i,j,kp))-5.*(-Uvel2(i,j,km)+Uvel2(i,j,kpp))+(-Uvel2(i,j,kmm)+Uvel2(i,j,kppp))) -
     3 (rhokm*wwLA*uuZLA - numdif/facA*
     3 ABS(rhokm*wwLA)* (10.*(-Uvel2(i,j,km)+Uvel2(i,j,k))-5.*(-Uvel2(i,j,kmm)+Uvel2(i,j,kp))+(-Uvel2(i,j,kmmm)+Uvel2(i,j,kpp)))) )
     3  *dzi
     +                         +
     3 facB*(rhokpp*wwRB*uuZRB - rhokmm*wwLB*uuZLB)*dzi/3.	 
     +                         -	 
     4 order4yes*0.5*(rho(i,j,k)+rho(ip,j,k))*(-Vvel2(im,j ,k)+9.*Vvel2(i,j ,k)+9.*Vvel2(ip,j ,k)-Vvel2(ipp,j ,k)     
     4                                  -Vvel2(im,jm,k)+9.*Vvel2(i,jm,k)+9.*Vvel2(ip,jm,k)-Vvel2(ipp,jm,k))**2 !4th order in x-dir
     4  / ( 256.0 * Ru(i) )
     +                         -
     4 (1.-order4yes)*0.5*(rho(i,j,k)+rho(ip,j,k))*( Vvel(i,j,k) + Vvel(ip,j,k) + Vvel(i,jm,k) + Vvel(ip,jm,k) )**2   !2nd order
     4  / ( 4.0 * Ru(i) )
     +                         )	 
           enddo
        enddo
      enddo
      return
      end


      subroutine advecv_C4A6(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,phiv)
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
	  integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),Ru2(-1:i1+1)
      real rho(0:i1,0:j1,0:k1),dzi,phiv(0:j1)
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,rhoipp,rhoimm,rhojpp,rhojmm,rhokpp,rhokmm
	real Axa,Bxa,DDxa,Axb,Bxb,DDxb,Axa2,Bxa2,DDxa2,Axb2,Bxb2,DDxb2,numdif,facA,facB,order4yes
	

	real ubb(0:i1,0:k1),ubf(0:i1,0:k1),Vvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real Uvel2(-2:i1+2,-2:j1+2,-2:k1+2),Wvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real uuRA,uuLA,vvRA,vvLA,wwRA,wwLA,uuRB,uuLB,vvRB,vvLB,wwRB,wwLB
	real vvXRA,vvXLA,vvYRA,vvYLA,vvZRA,vvZLA,vvXRB,vvXLB,vvYRB,vvYLB,vvZRB,vvZLB
	real rho2(-1:i1+1,-1:j1+1,-1:k1+1),phivt2(-1:je*px+2)

	
	
	Ru2(0:i1)=Ru
	Ru2(-1)=Ru(0)-(Ru(1)-Ru(0)) !needed for periodicx sim
 	Ru2(i1+1)=Ru(i1)+(Ru(i1)-Ru(i1-1)) !needed for periodicx sim	
	phivt2(0:je*px+1)=phivt
	phivt2(-1)=phivt(0)-(phivt(1)-phivt(0)) !needed for periodicy sim
 	phivt2(je*px+1+1)=phivt(je*px+1)+(phivt(je*px+1)-phivt(je*px+1-1)) !needed for periodicy sim
	
	Uvel2(0:i1,0:j1,0:k1)=Uvel

c get stuff from other CPU's
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
!	Uvel2(-2:i1+2,-2:j1+2,0) =0.
	Uvel2(-2:i1+2,-2:j1+2,-1)=  Uvel2(-2:i1+2,-2:j1+2,0)
	Uvel2(-2:i1+2,-2:j1+2,-2)=  Uvel2(-2:i1+2,-2:j1+2,0)
	!Uvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Uvel2(-2:i1+2,-2:j1+2,1)+1./3.*Uvel2(-2:i1+2,-2:j1+2,2)
	!Uvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Uvel2(-2:i1+2,-2:j1+2,1)+2.*Uvel2(-2:i1+2,-2:j1+2,2)	
	Uvel2(-2:i1+2,-2:j1+2,k1+1)=Uvel2(-2:i1+2,-2:j1+2,k1)
	Uvel2(-2:i1+2,-2:j1+2,k1+2)=Uvel2(-2:i1+2,-2:j1+2,k1)

	
	Vvel2(0:i1,0:j1,0:k1)=Vvel

c get stuff from other CPU's
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
!	Vvel2(-2:i1+2,-2:j1+2,0) =0.
	Vvel2(-2:i1+2,-2:j1+2,-1)=  Vvel2(-2:i1+2,-2:j1+2,0)
	Vvel2(-2:i1+2,-2:j1+2,-2)=  Vvel2(-2:i1+2,-2:j1+2,0)
	!Vvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Vvel2(-2:i1+2,-2:j1+2,1)+1./3.*Vvel2(-2:i1+2,-2:j1+2,2)
	!Vvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Vvel2(-2:i1+2,-2:j1+2,1)+2.*Vvel2(-2:i1+2,-2:j1+2,2)	
	Vvel2(-2:i1+2,-2:j1+2,k1+1)=Vvel2(-2:i1+2,-2:j1+2,k1)
	Vvel2(-2:i1+2,-2:j1+2,k1+2)=Vvel2(-2:i1+2,-2:j1+2,k1)

		Wvel2(0:i1,0:j1,0:k1)=Wvel

c get stuff from other CPU's
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
!	Wvel2(-2:i1+2,-2:j1+2,-1)=  -Wvel2(-2:i1+2,-2:j1+2,1)
!	Wvel2(-2:i1+2,-2:j1+2,-2)=  -Wvel2(-2:i1+2,-2:j1+2,2)
	Wvel2(-2:i1+2,-2:j1+2,k1+1)=Wvel2(-2:i1+2,-2:j1+2,k1)
	Wvel2(-2:i1+2,-2:j1+2,k1+2)=Wvel2(-2:i1+2,-2:j1+2,k1)

	rho2(0:i1,0:j1,0:k1)=rho

c get stuff from other CPU's
	  call shiftf2(rho,ubf)
	  call shiftb2(rho,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = rho(i,0,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =rho(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	
	if (periodicx.eq.0.or.periodicx.eq.2) then
		rho2(-1,-1:j1+1,0:k1)=  rho2(1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(ie,-1:j1+1,0:k1)
	else
		rho2(-1,-1:j1+1,0:k1)=  rho2(ie-1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(2,-1:j1+1,0:k1)
	endif
	rho2(-1:i1+1,-1:j1+1,-1)=  rho2(-1:i1+1,-1:j1+1,0)
	rho2(-1:i1+1,-1:j1+1,k1+1)=rho2(-1:i1+1,-1:j1+1,k1)
	
	dzi=1./dz 
      do k=kb,ke
      kp=k+1
      km=k-1
		  kpp=k+2 !MIN(k+2,k1)
		  kppp=k+3 !MIN(k+3,k1)
		  kmm=k-2 !MAX(k-2,0)
		  kmmm=k-3 !MAX(k-3,0)
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
	  
          rhoip =0.25*(rho(i,j,k)+rho(ip,j,k)+rho(i,jp,k)+rho(ip,jp,k))
          rhoim =0.25*(rho(i,j,k)+rho(im,j,k)+rho(i,jp,k)+rho(im,jp,k))
          rhojp =rho(i,jp,k)
          rhojm =rho(i,j,k )
          rhokp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
          rhokm =0.25*(rho(i,j,k)+rho(i,j,km)+rho(i,jp,k)+rho(i,jp,km))
          rhoipp =0.25*(rho2(ip,j,k)+rho2(ipp,j,k)+rho2(ip,jp,k)+rho2(ipp,jp,k))
          rhoimm =0.25*(rho2(im,j,k)+rho2(imm,j,k)+rho2(im,jp,k)+rho2(imm,jp,k))
          rhojpp =rho2(i,jpp,k)
          rhojmm =rho2(i,jm ,k )
          rhokpp =0.25*(rho2(i,j,kp)+rho2(i,j,kpp)+rho2(i,jp,kp)+rho2(i,jp,kpp))
          rhokmm =0.25*(rho2(i,j,km)+rho2(i,j,kmm)+rho2(i,jp,km)+rho2(i,jp,kmm))
		  
!      uuR=(Uvel(i,j,k)+Uvel(i,jp,k))*rhoip*Ru(i) !*Ru(i)
!      uuL=(Uvel(im,j,k)+Uvel(im,jp,k))*rhoim*Ru(im) !*Ru(im)
!      vvR=(Vvel(i,j,k)+Vvel(i,jp,k))*rhojp
!      vvL=(Vvel(i,jm,k)+Vvel(i,j,k))*rhojm
!      wwR=(Wvel(i,j,k)+Wvel(i,jp,k))*rhokp
!      wwL=(Wvel(i,j,km)+Wvel(i,jp,km))*rhokm

		   if ((periodicx.eq.1.and.periodicy.eq.1.and.(k.le.2.or.k.ge.ke-1))
     &	   .or.((periodicx.eq.1.and.periodicy.ne.1).and.
     &           (k.le.2.or.k.ge.ke-1.or.(rank.eq.0.and.j.le.2).or.(rank.eq.px-1.and.j.ge.je-1))) 
     &     .or.((periodicy.eq.1.and.periodicx.ne.1).and.(i.le.2.or.i.ge.ie-1.or.k.le.2.or.k.ge.ke-1))
     &     .or.((periodicy.ne.1.and.periodicx.ne.1).and.
     &           (i.le.2.or.i.ge.ie-1.or.k.le.2.or.k.ge.ke-1.or.(rank.eq.0.and.j.le.2).or.(rank.eq.px-1.and.j.ge.je-1)))) then  !switch to cds2
!			if (.false.) then
			  Axa=1.
			  Bxa=0.
			  DDxa=1./(Axa+Bxa)
			  Axb=1.
			  Bxb=0.
			  DDxb=1./(Axb+Bxb)	
			  Axa2=1.
			  Bxa2=0.
			  DDxa2=1./(Axa2+Bxa2)
			  Axb2=1.
			  Bxb2=0.
			  DDxb2=1./(Axb2+Bxb2)	
			  facA=1.
			  facB=0.
		  else
			  Axa=9. !1. !9.
			  Bxa=-1. !0. !-1.
			  DDxa=1./(Axa+Bxa)
			  Axb=0. !9. !0. !9.
			  Bxb=1. !-1. !1. !-1.
			  DDxb=1./(Axb+Bxb)	
			  Axa2=9. !9. !1. 
			  Bxa2=-1. !-1. !0. 
			  DDxa2=1./(Axa2+Bxa2)
			  Axb2=0. 
			  Bxb2=1. 
			  DDxb2=1./(Axb2+Bxb2)				  
			  facA=1. !1.+MIN(1.,DBLE(k-1)/3.)*1./8. !9./8.
			  facB=0. !MIN(1.,DBLE(k-1)/3.)*-1./8.    !-1./8.							
			  order4yes=0.
		  endif	

	  ! advection velocities:
      uuRA=DDxa*(Axa*(Uvel2(i ,j ,k )+Uvel2(i ,jp,k ))+Bxa*(Uvel2(i ,jm ,k )+Uvel2(i ,jpp,k )))
      uuLA=DDxa*(Axa*(Uvel2(im,j ,k )+Uvel2(im,jp,k ))+Bxa*(Uvel2(im,jm ,k )+Uvel2(im,jpp,k )))
	  vvRA=DDxa*(Axa*(Vvel2(i ,j ,k )+Vvel2(i ,jp,k ))+Bxa*(Vvel2(i ,jm ,k )+Vvel2(i ,jpp,k )))
      vvLA=DDxa*(Axa*(Vvel2(i ,jm,k )+Vvel2(i ,j ,k ))+Bxa*(Vvel2(i ,jmm,k )+Vvel2(i ,jp ,k )))
      wwRA=DDxa*(Axa*(Wvel2(i ,j ,k )+Wvel2(i ,jp,k ))+Bxa*(Wvel2(i ,jm ,k )+Wvel2(i ,jpp,k )))
      wwLA=DDxa*(Axa*(Wvel2(i ,j ,km)+Wvel2(i ,jp,km))+Bxa*(Wvel2(i ,jm ,km)+Wvel2(i ,jpp,km)))
  
      uuRB=DDxb*(Axb*(Uvel2(ip ,j  ,k  )+Uvel2(ip ,jp ,k  ))+Bxb*(Uvel2(ip ,jm  ,k  )+Uvel2(ip ,jpp ,k  )))
      uuLB=DDxb*(Axb*(Uvel2(imm,j  ,k  )+Uvel2(imm,jp ,k  ))+Bxb*(Uvel2(imm,jm  ,k  )+Uvel2(imm,jpp ,k  )))
	  vvRB=DDxb*(Axb*(Vvel2(i  ,jp ,k  )+Vvel2(i  ,jpp,k  ))+Bxb*(Vvel2(i  ,j   ,k  )+Vvel2(i  ,jppp,k  )))
      vvLB=DDxb*(Axb*(Vvel2(i  ,jmm,k  )+Vvel2(i  ,jm ,k  ))+Bxb*(Vvel2(i  ,jmmm,k  )+Vvel2(i  ,j   ,k  )))
      wwRB=DDxb*(Axb*(Wvel2(i  ,j  ,kp )+Wvel2(i  ,jp ,kp ))+Bxb*(Wvel2(i  ,jm  ,kp )+Wvel2(i  ,jpp ,kp )))
      wwLB=DDxb*(Axb*(Wvel2(i  ,j  ,kmm)+Wvel2(i  ,jp ,kmm))+Bxb*(Wvel2(i  ,jm  ,kmm)+Wvel2(i  ,jpp ,kmm)))
	  ! advected velocities:
      vvXRA=DDxa2*(Axa2*(Vvel2(i ,j ,k )+Vvel2(ip,j ,k ))+Bxa2*(Vvel2(im ,j  ,k )+Vvel2(ipp,j  ,k  )))
      vvXLA=DDxa2*(Axa2*(Vvel2(im,j ,k )+Vvel2(i ,j ,k ))+Bxa2*(Vvel2(imm,j  ,k )+Vvel2(ip ,j  ,k  )))
	  vvYRA=DDxa2*(Axa2*(Vvel2(i ,j ,k )+Vvel2(i ,jp,k ))+Bxa2*(Vvel2(i  ,jm ,k )+Vvel2(i  ,jpp,k  )))
      vvYLA=DDxa2*(Axa2*(Vvel2(i ,jm,k )+Vvel2(i ,j ,k ))+Bxa2*(Vvel2(i  ,jmm,k )+Vvel2(i  ,jp ,k  )))
      vvZRA=DDxa2*(Axa2*(Vvel2(i ,j ,k )+Vvel2(i ,j ,kp))+Bxa2*(Vvel2(i  ,j  ,km )+Vvel2(i ,j  ,kpp)))
      vvZLA=DDxa2*(Axa2*(Vvel2(i ,j ,km)+Vvel2(i ,j ,k ))+Bxa2*(Vvel2(i  ,j  ,kmm)+Vvel2(i ,j  ,kp )))
      vvXRB=DDxb2*(Axb2*(Vvel2(ip  ,j   ,k   )+Vvel2(ipp ,j   ,k   ))+Bxb2*(Vvel2(i   ,j   ,k   )+Vvel2(ippp,j   ,k   )))
      vvXLB=DDxb2*(Axb2*(Vvel2(imm ,j   ,k   )+Vvel2(im  ,j   ,k   ))+Bxb2*(Vvel2(immm,j   ,k   )+Vvel2(i   ,j   ,k   )))
	  vvYRB=DDxb2*(Axb2*(Vvel2(i   ,jp  ,k   )+Vvel2(i   ,jpp ,k   ))+Bxb2*(Vvel2(i   ,j   ,k   )+Vvel2(i   ,jppp,k   )))
      vvYLB=DDxb2*(Axb2*(Vvel2(i   ,jmm ,k   )+Vvel2(i   ,jm  ,k   ))+Bxb2*(Vvel2(i   ,jmmm,k   )+Vvel2(i   ,j   ,k   )))
      vvZRB=DDxb2*(Axb2*(Vvel2(i   ,j   ,kp  )+Vvel2(i   ,j   ,kpp ))+Bxb2*(Vvel2(i   ,j   ,k   )+Vvel2(i   ,j   ,kppp)))
      vvZLB=DDxb2*(Axb2*(Vvel2(i   ,j   ,kmm )+Vvel2(i   ,j   ,km  ))+Bxb2*(Vvel2(i   ,j   ,kmmm)+Vvel2(i   ,j   ,k   )))

!	  IF (ABS(Uvel2(i,j,k)-Uvel(i,j,k)).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,U,U2',rank,i,j,k,Uvel2(i,j,k),Uvel(i,j,k)
!	  ENDIF
!	  IF (ABS(Vvel2(i,j,k)-Vvel(i,j,k)).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,V,V2',rank,i,j,k,Vvel2(i,j,k),Vvel(i,j,k)
!	  ENDIF
!	  IF (ABS(Wvel2(i,j,k)-Wvel(i,j,k)).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,W,W2',rank,i,j,k,Wvel2(i,j,k),Wvel(i,j,k)
!	  ENDIF
!	  IF (ABS(Uvel2(im,jm,km)-Uvel(im,jm,km)).gt.1e-14) THEN
!		write(*,*),'rank,im,jm,km,U,U2',rank,im,jm,km,Uvel2(im,jm,km),Uvel(im,jm,km)
!	  ENDIF
!	  IF (ABS(Vvel2(im,jm,km)-Vvel(im,jm,km)).gt.1e-14) THEN
!		write(*,*),'rank,im,jm,km,V,V2',rank,im,jm,km,Vvel2(im,jm,km),Vvel(im,jm,km)
!	  ENDIF
!	  IF (ABS(Wvel2(im,jm,km)-Wvel(im,jm,km)).gt.1e-14) THEN
!		write(*,*),'rank,im,jm,km,W,W2',rank,im,jm,km,Wvel2(im,jm,km),Wvel(im,jm,km)
!	  ENDIF
!	  IF (ABS(Uvel2(ip,jp,kp)-Uvel(ip,jp,kp)).gt.1e-14) THEN
!		write(*,*),'rank,ip,jp,kp,U,U2',rank,ip,jp,kp,Uvel2(ip,jp,kp),Uvel(ip,jp,kp)
!	  ENDIF
!	  IF (ABS(Vvel2(ip,jp,kp)-Vvel(ip,jp,kp)).gt.1e-14) THEN
!		write(*,*),'rank,ip,jp,kp,V,V2',rank,ip,jp,kp,Vvel2(ip,jp,kp),Vvel(ip,jp,kp)
!	  ENDIF
!	  IF (ABS(Wvel2(ip,jp,kp)-Wvel(ip,jp,kp)).gt.1e-14) THEN
!		write(*,*),'rank,ip,jp,kp,W,W2',rank,ip,jp,kp,Wvel2(ip,jp,kp),Wvel(ip,jp,kp)
!	  ENDIF	
!      uuRB=DDya*(Aya*(Uvel(i ,j ,k )+Uvel(i ,jp,k ))+Bya*(Uvel2(i ,jm ,k )+Uvel2(i ,jpp,k )))
!      uuLB=DDya*(Aya*(Uvel(im,j ,k )+Uvel(im,jp,k ))+Bya*(Uvel2(im,jm ,k )+Uvel2(im,jpp,k )))
!	  vvRB=DDya*(Aya*(Vvel(i ,j ,k )+Vvel(i ,jp,k ))+Bya*(Vvel2(i ,jm ,k )+Vvel2(i ,jpp,k )))
!      vvLB=DDya*(Aya*(Vvel(i ,jm,k )+Vvel(i ,j ,k ))+Bya*(Vvel2(i ,jmm,k )+Vvel2(i ,jp ,k )))
!      wwRB=DDya*(Aya*(Wvel(i ,j ,k )+Wvel(i ,jp,k ))+Bya*(Wvel2(i ,jm ,k )+Wvel2(i ,jpp,k )))
!      wwLB=DDya*(Aya*(Wvel(i ,j ,km)+Wvel(i ,jp,km))+Bya*(Wvel2(i ,jm ,km)+Wvel2(i ,jpp,km)))
!	  IF (ABS(uuRB-uuRA).gt.1e-14.or.ABS(uuLB-uuLA).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,uU',rank,i,j,k,uuRA,uuRB,uuLB,uuLA 
!	  ENDIF
!	  IF (ABS(vvRB-vvRA).gt.1e-14.or.ABS(vvLB-vvLA).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,uV',rank,i,j,k,vvRA,vvRB,vvLB,vvLA 
!	  ENDIF
!	  IF (ABS(wwRB-wwRA).gt.1e-14.or.ABS(wwLB-wwLA).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,uW',rank,i,j,k,wwRA,wwRB,wwLB,wwLA 
!	  ENDIF
	  
!      uuR=DDy*(Ay*(Uvel2(i,j,k)+Uvel2(i,jp,k))+By*(Uvel2(i,jm,k)+Uvel2(i,jpp,k))+Cy*(Uvel2(i,jmm,k)+Uvel2(i,jppp,k)))*rhoip*Ru(i) !*Ru(i)
!      uuL=DDy*(Ay*(Uvel2(im,j,k)+Uvel2(im,jp,k))+By*(Uvel2(im,jm,k)+Uvel2(im,jpp,k))+Cy*(Uvel2(im,jmm,k)+Uvel2(im,jppp,k)))
!     &	  *rhoim*Ru(im) !*Ru(im)
!      vvR=DDy*(Ay*(Vvel2(i,j,k)+Vvel2(i,jp,k))+By*(Vvel2(i,jm,k)+Vvel2(i,jpp,k))+Cy*(Vvel2(i,jmm,k)+Vvel2(i,jppp,k)))*rhojp
!      vvL=DDy*(Ay*(Vvel2(i,jm,k)+Vvel2(i,j,k))+By*(Vvel2(i,jmm,k)+Vvel2(i,jp,k))+Cy*(Vvel2(i,jmmm,k)+Vvel2(i,jpp,k)))*rhojm
!      wwR=DDy*(Ay*(Wvel2(i,j,k)+Wvel2(i,jp,k))+By*(Wvel2(i,jm,k)+Wvel2(i,jpp,k))+Cy*(Wvel2(i,jmm,k)+Wvel2(i,jppp,k)))*rhokp
!      wwL=DDy*(Ay*(Wvel2(i,j,km)+Wvel2(i,jp,km))+By*(Wvel2(i,jm,km)+Wvel2(i,jpp,km))+Cy*(Wvel2(i,jmm,km)+Wvel2(i,jppp,km)))*rhokm
  
      putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (

     1 facA*(uuRA*rhoip*Ru(i)*vvXRA-numdif/facA*ABS(uuRA*rhoip*Ru(i))*
     1 (10.*(-Vvel2(i,j,k)+Vvel2(ip,j,k))-5.*(-Vvel2(im,j,k)+Vvel2(ipp,j,k))+(-Vvel2(imm,j,k)+Vvel2(ippp,j,k))) -
     1 (uuLA*rhoim*Ru(im)*vvXLA-numdif/facA*ABS(uuLA*rhoim*Ru(im))*
     1 (10.*(-Vvel2(im,j,k)+Vvel2(i,j,k))-5.*(-Vvel2(imm,j,k)+Vvel2(ip,j,k))+(-Vvel2(immm,j,k)+Vvel2(ipp,j,k)))) )
     1  / ( Rp(i)* dr(i) )   !/ ( Rp(i) * Rp(i)* dr(i) ) 
     +                         +
     1 facB*(uuRB*rhoipp*Ru2(ip)*vvXRB-uuLB*rhoimm*Ru2(imm)*vvXLB)
     1  / ( Rp(i)* (Ru2(ip)-Ru2(imm)))   !/ ( Rp(i) * Rp(i)* dr(i) ) 
     +                         +
     2 facA*(vvRA*rhojp*vvYRA-numdif/facA*ABS(vvRA*rhojp)* 
     2 (10.*(-Vvel2(i,j,k)+Vvel2(i,jp,k))-5.*(-Vvel2(i,jm,k)+Vvel2(i,jpp,k))+(-Vvel2(i,jmm,k)+Vvel2(i,jppp,k))) -
     2 (vvLA*rhojm*vvYLA-numdif/facA*ABS(vvLA*rhojm)* 
     2 (10.*(-Vvel2(i,jm,k)+Vvel2(i,j,k))-5.*(-Vvel2(i,jmm,k)+Vvel2(i,jp,k))+(-Vvel2(i,jmmm,k)+Vvel2(i,jpp,k)))) )
     2  / ( Rp(i) * (phiv(j)-phiv(jm)) ) !/ ( Rp(i) * (phivt2(rank*je+j)-phivt2(rank*je+jm)) )
     +                         +
     2 facB*(vvRB*rhojpp*vvYRB-vvLB*rhojmm*vvYLB)
     2  / ( Rp(i) * (phivt2(rank*je+jp)-phivt2(rank*je+jmm)) )
     +                         +
     3 facA*(wwRA*rhokp*vvZRA-numdif/facA*ABS(wwRA*rhokp)*
     3 (10.*(-Vvel2(i,j,k)+Vvel2(i,j,kp))-5.*(-Vvel2(i,j,km)+Vvel2(i,j,kpp))+(-Vvel2(i,j,kmm)+Vvel2(i,j,kppp))) -
     3 (wwLA*rhokm*vvZLA-numdif/facA*ABS(wwLA*rhokm)*
     3 (10.*(-Vvel2(i,j,km)+Vvel2(i,j,k))-5.*(-Vvel2(i,j,kmm)+Vvel2(i,j,kp))+(-Vvel2(i,j,kmmm)+Vvel2(i,j,kpp)))) )
     3  *dzi
     +                         +
     3 facB*(wwRB*rhokpp*vvZRB-wwLB*rhokmm*vvZLB)*dzi/3.
     +                         +
     4   order4yes*0.5*(rho(i,j,k)+rho(i,jp,k))*										!4th order in y-dir
     4   (-Uvel2(im,jm,k)+9.*Uvel2(im,j,k)+9.*Uvel2(im,jp,k)-Uvel2(im,jpp,k)
     4	  -Uvel2(i ,jm,k)+9.*Uvel2(i ,j,k)+9.*Uvel2(i ,jp,k)-Uvel2(i ,jpp,k))*Vvel(i,j,k)
     4   / ( 8.0*Rp(i) )
     +                         +	 
     4   (1.-order4yes)*0.5*(rho(i,j,k)+rho(i,jp,k))*									!2nd order
     4   (Uvel(im,j,k)+Uvel(i,j,k)+Uvel(im,jp,k)+Uvel(i,jp,k))*Vvel(i,j,k)
     4   / ( Rp(i) )	 
     +                         )
           enddo
        enddo
      enddo
      return
      end

      subroutine advecw_C4A6(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,phiv)
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
	  integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),Ru2(-1:i1+1)
      real rho(0:i1,0:j1,0:k1),dzi,phiv(0:j1)
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,rhoipp,rhoimm,rhojpp,rhojmm,rhokpp,rhokmm
	real Axa,Bxa,DDxa,Axb,Bxb,DDxb,Axa2,Bxa2,DDxa2,Axb2,Bxb2,DDxb2,numdif,facA,facB
	
	real ubb(0:i1,0:k1),ubf(0:i1,0:k1),Wvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real Uvel2(-2:i1+2,-2:j1+2,-2:k1+2),Vvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real uuRA,uuLA,vvRA,vvLA,wwRA,wwLA,uuRB,uuLB,vvRB,vvLB,wwRB,wwLB
	real wwXRA,wwXLA,wwYRA,wwYLA,wwZRA,wwZLA,wwXRB,wwXLB,wwYRB,wwYLB,wwZRB,wwZLB
	real rho2(-1:i1+1,-1:j1+1,-1:k1+1),phivt2(-1:je*px+2)

	
	
	Ru2(0:i1)=Ru
	Ru2(-1)=Ru(0)-(Ru(1)-Ru(0)) !needed for periodicx sim
 	Ru2(i1+1)=Ru(i1)+(Ru(i1)-Ru(i1-1)) !needed for periodicx sim	
	phivt2(0:je*px+1)=phivt
	phivt2(-1)=phivt(0)-(phivt(1)-phivt(0)) !needed for periodicy sim
 	phivt2(je*px+1+1)=phivt(je*px+1)+(phivt(je*px+1)-phivt(je*px+1-1)) !needed for periodicy sim
	
	Uvel2(0:i1,0:j1,0:k1)=Uvel

c get stuff from other CPU's
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
!	Uvel2(-2:i1+2,-2:j1+2,0) =0.
	Uvel2(-2:i1+2,-2:j1+2,-1)=  Uvel2(-2:i1+2,-2:j1+2,0)
	Uvel2(-2:i1+2,-2:j1+2,-2)=  Uvel2(-2:i1+2,-2:j1+2,0)
	!Uvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Uvel2(-2:i1+2,-2:j1+2,1)+1./3.*Uvel2(-2:i1+2,-2:j1+2,2)
	!Uvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Uvel2(-2:i1+2,-2:j1+2,1)+2.*Uvel2(-2:i1+2,-2:j1+2,2)	
	Uvel2(-2:i1+2,-2:j1+2,k1+1)=Uvel2(-2:i1+2,-2:j1+2,k1)
	Uvel2(-2:i1+2,-2:j1+2,k1+2)=Uvel2(-2:i1+2,-2:j1+2,k1)

	
	Vvel2(0:i1,0:j1,0:k1)=Vvel

c get stuff from other CPU's
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
!	Vvel2(-2:i1+2,-2:j1+2,0) =0.
	Vvel2(-2:i1+2,-2:j1+2,-1)=  Vvel2(-2:i1+2,-2:j1+2,0)
	Vvel2(-2:i1+2,-2:j1+2,-2)=  Vvel2(-2:i1+2,-2:j1+2,0)
	!Vvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Vvel2(-2:i1+2,-2:j1+2,1)+1./3.*Vvel2(-2:i1+2,-2:j1+2,2)
	!Vvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Vvel2(-2:i1+2,-2:j1+2,1)+2.*Vvel2(-2:i1+2,-2:j1+2,2)	
	Vvel2(-2:i1+2,-2:j1+2,k1+1)=Vvel2(-2:i1+2,-2:j1+2,k1)
	Vvel2(-2:i1+2,-2:j1+2,k1+2)=Vvel2(-2:i1+2,-2:j1+2,k1)

		Wvel2(0:i1,0:j1,0:k1)=Wvel

c get stuff from other CPU's
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
!	Wvel2(-2:i1+2,-2:j1+2,-1)=  -Wvel2(-2:i1+2,-2:j1+2,1)
!	Wvel2(-2:i1+2,-2:j1+2,-2)=  -Wvel2(-2:i1+2,-2:j1+2,2)
	Wvel2(-2:i1+2,-2:j1+2,k1+1)=Wvel2(-2:i1+2,-2:j1+2,k1)
	Wvel2(-2:i1+2,-2:j1+2,k1+2)=Wvel2(-2:i1+2,-2:j1+2,k1)

	rho2(0:i1,0:j1,0:k1)=rho

c get stuff from other CPU's
	  call shiftf2(rho,ubf)
	  call shiftb2(rho,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = rho(i,0,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =rho(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	
	if (periodicx.eq.0.or.periodicx.eq.2) then
		rho2(-1,-1:j1+1,0:k1)=  rho2(1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(ie,-1:j1+1,0:k1)
	else
		rho2(-1,-1:j1+1,0:k1)=  rho2(ie-1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(2,-1:j1+1,0:k1)
	endif
	rho2(-1:i1+1,-1:j1+1,-1)=  rho2(-1:i1+1,-1:j1+1,0)
	rho2(-1:i1+1,-1:j1+1,k1+1)=rho2(-1:i1+1,-1:j1+1,k1)
	
	dzi=1./dz 
      do k=kb,ke
      kp=k+1
      km=k-1
		  kpp=k+2 !MIN(k+2,k1)
		  kppp=k+3 !MIN(k+3,k1)
		  kmm=k-2 !MAX(k-2,0)
		  kmmm=k-3 !MAX(k-3,0)
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
	  
          rhoip =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
          rhoim =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(im,j,k)+rho(im,j,kp))
          rhojp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
          rhojm =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jm,k)+rho(i,jm,kp))
		  rhokp =rho(i,j,kp)
		  rhokm =rho(i,j,k )
          rhoipp =0.25*(rho2(ip,j,k)+rho2(ip,j,kp)+rho2(ipp,j,k)+rho2(ipp,j,kp))
          rhoimm =0.25*(rho2(im,j,k)+rho2(im,j,kp)+rho2(imm,j,k)+rho2(imm,j,kp))
          rhojpp =0.25*(rho2(i,jp,k)+rho2(i,jp,kp)+rho2(i,jpp,k)+rho2(i,jpp,kp))
          rhojmm =0.25*(rho2(i,jm,k)+rho2(i,jm,kp)+rho2(i,jmm,k)+rho2(i,jmm,kp))
		  rhokpp =rho2(i,j,kpp)
		  rhokmm =rho2(i,j,km )


	  
!       uuR=DDz*(Az*(Uvel2(i,j,k)+Uvel2(i,j,kp))+Bz*(Uvel2(i,j,km)+Uvel2(i,j,kpp))+Cz*(Uvel2(i,j,kmm)+Uvel2(i,j,kppp)))*rhoip*Ru(i)
!       uuL=DDz*(Az*(Uvel2(im,j,k)+Uvel2(im,j,kp))+Bz*(Uvel2(im,j,km)+Uvel2(im,j,kpp))+Cz*(Uvel2(im,j,kmm)+Uvel2(im,j,kppp)))
!     &	  *rhoim*Ru(im)
!       vvR=DDz*(Az*(Vvel2(i,j,k)+Vvel2(i,j,kp))+Bz*(Vvel2(i,j,km)+Vvel2(i,j,kpp))+Cz*(Vvel2(i,j,kmm)+Vvel2(i,j,kppp)))*rhojp
!       vvL=DDz*(Az*(Vvel2(i,jm,k)+Vvel2(i,jm,kp))+Bz*(Vvel2(i,jm,km)+Vvel2(i,jm,kpp))+Cz*(Vvel2(i,jm,kmm)+Vvel2(i,jm,kppp)))*rhojm	   
!	     wwR=DDz*(Az*(Wvel2(i,j,k)+Wvel2(i,j,kp))+Bz*(Wvel2(i,j,km)+Wvel2(i,j,kpp))+Cz*(Wvel2(i,j,kmm)+Wvel2(i,j,kppp)))*rho(i,j,kp)
!		wwL=DDz*(Az*(Wvel2(i,j,km)+Wvel2(i,j,k))+Bz*(Wvel2(i,j,kmm)+Wvel2(i,j,kp))+Cz*(Wvel2(i,j,kmmm)+Wvel2(i,j,kpp)))*rho(i,j,k)
!      	    uuR=(Uvel(i,j,k)+Uvel(i,j,kp))*rhoip*Ru(i)
!      	    uuL=(Uvel(im,j,k)+Uvel(im,j,kp))*rhoim*Ru(im)
!   	    vvR=(Vvel(i,j,k)+Vvel(i,j,kp))*rhojp
!    	    vvL=(Vvel(i,jm,k)+Vvel(i,jm,kp))*rhojm
!    	    wwR=(Wvel(i,j,k)+Wvel(i,j,kp))*rho(i,j,kp)
!   	    wwL=(Wvel(i,j,k)+Wvel(i,j,km))*rho(i,j,k)	
		   if ((periodicx.eq.1.and.periodicy.eq.1.and.(k.le.2.or.k.ge.ke-1))
     &	   .or.((periodicx.eq.1.and.periodicy.ne.1).and.
     &           (k.le.2.or.k.ge.ke-1.or.(rank.eq.0.and.j.le.2).or.(rank.eq.px-1.and.j.ge.je-1))) 
     &     .or.((periodicy.eq.1.and.periodicx.ne.1).and.(i.le.2.or.i.ge.ie-1.or.k.le.2.or.k.ge.ke-1))
     &     .or.((periodicy.ne.1.and.periodicx.ne.1).and.
     &           (i.le.2.or.i.ge.ie-1.or.k.le.2.or.k.ge.ke-1.or.(rank.eq.0.and.j.le.2).or.(rank.eq.px-1.and.j.ge.je-1)))) then  !switch to cds2
!			if (.false.) then
			  Axa=1.
			  Bxa=0.
			  DDxa=1./(Axa+Bxa)
			  Axb=1.
			  Bxb=0.
			  DDxb=1./(Axb+Bxb)	
			  Axa2=1.
			  Bxa2=0.
			  DDxa2=1./(Axa2+Bxa2)
			  Axb2=1.
			  Bxb2=0.
			  DDxb2=1./(Axb2+Bxb2)	
			  facA=1.
			  facB=0.
		  else
			  Axa=9. !1. !9.
			  Bxa=-1. !0. !-1.
			  DDxa=1./(Axa+Bxa)
			  Axb=0. !9. !0. !9.
			  Bxb=1. !-1. !1. !-1.
			  DDxb=1./(Axb+Bxb)	
			  Axa2=9. !9. !1. 
			  Bxa2=-1. !-1. !0. 
			  DDxa2=1./(Axa2+Bxa2)
			  Axb2=0. 
			  Bxb2=1. 
			  DDxb2=1./(Axb2+Bxb2)				  
			  facA=1. !1.+MIN(1.,DBLE(k-1)/3.)*1./8. !9./8.
			  facB=0. !MIN(1.,DBLE(k-1)/3.)*-1./8.    !-1./8.				  
		  endif	

	  ! advection velocities:
      uuRA=DDxa*(Axa*(Uvel2(i ,j ,k )+Uvel2(i ,j ,kp))+Bxa*(Uvel2(i ,j ,km )+Uvel2(i ,j ,kpp)))
      uuLA=DDxa*(Axa*(Uvel2(im,j ,k )+Uvel2(im,j ,kp))+Bxa*(Uvel2(im,j ,km )+Uvel2(im,j ,kpp)))
	  vvRA=DDxa*(Axa*(Vvel2(i ,j ,k )+Vvel2(i ,j ,kp))+Bxa*(Vvel2(i ,j ,km )+Vvel2(i ,j ,kpp)))
      vvLA=DDxa*(Axa*(Vvel2(i ,jm,k )+Vvel2(i ,jm,kp))+Bxa*(Vvel2(i ,jm,km )+Vvel2(i ,jm,kpp)))
      wwRA=DDxa*(Axa*(Wvel2(i ,j ,k )+Wvel2(i ,j ,kp))+Bxa*(Wvel2(i ,j ,km )+Wvel2(i ,j ,kpp)))
      wwLA=DDxa*(Axa*(Wvel2(i ,j ,km)+Wvel2(i ,j ,k ))+Bxa*(Wvel2(i ,j ,kmm)+Wvel2(i ,j ,kp )))
	  
      uuRB=DDxb*(Axb*(Uvel2(ip ,j  ,k  )+Uvel2(ip ,j  ,kp ))+Bxb*(Uvel2(ip ,j  ,km  )+Uvel2(ip ,j  ,kpp )))
      uuLB=DDxb*(Axb*(Uvel2(imm,j  ,k  )+Uvel2(imm,j  ,kp ))+Bxb*(Uvel2(imm,j  ,km  )+Uvel2(imm,j  ,kpp )))
	  vvRB=DDxb*(Axb*(Vvel2(i  ,jp ,k  )+Vvel2(i  ,jp ,kp ))+Bxb*(Vvel2(i  ,jp ,km  )+Vvel2(i  ,jp ,kpp )))
      vvLB=DDxb*(Axb*(Vvel2(i  ,jmm,k  )+Vvel2(i  ,jmm,kp ))+Bxb*(Vvel2(i  ,jmm,km  )+Vvel2(i  ,jmm,kpp )))
      wwRB=DDxb*(Axb*(Wvel2(i  ,j  ,kp )+Wvel2(i  ,j  ,kpp))+Bxb*(Wvel2(i  ,j  ,k   )+Wvel2(i  ,j  ,kppp)))
      wwLB=DDxb*(Axb*(Wvel2(i  ,j  ,kmm)+Wvel2(i  ,j  ,km ))+Bxb*(Wvel2(i  ,j  ,kmmm)+Wvel2(i  ,j  ,k   )))
	  ! advected velocities:
      wwXRA=DDxa2*(Axa2*(Wvel2(i ,j ,k )+Wvel2(ip,j ,k ))+Bxa2*(Wvel2(im ,j  ,k  )+Wvel2(ipp,j  ,k  )))
      wwXLA=DDxa2*(Axa2*(Wvel2(im,j ,k )+Wvel2(i ,j ,k ))+Bxa2*(Wvel2(imm,j  ,k  )+Wvel2(ip ,j  ,k  )))
	  wwYRA=DDxa2*(Axa2*(Wvel2(i ,j ,k )+Wvel2(i ,jp,k ))+Bxa2*(Wvel2(i  ,jm ,k  )+Wvel2(i  ,jpp,k  )))
      wwYLA=DDxa2*(Axa2*(Wvel2(i ,jm,k )+Wvel2(i ,j ,k ))+Bxa2*(Wvel2(i  ,jmm,k  )+Wvel2(i  ,jp ,k  )))
      wwZRA=DDxa2*(Axa2*(Wvel2(i ,j ,k )+Wvel2(i ,j ,kp))+Bxa2*(Wvel2(i  ,j  ,km )+Wvel2(i  ,j  ,kpp)))
      wwZLA=DDxa2*(Axa2*(Wvel2(i ,j ,km)+Wvel2(i ,j ,k ))+Bxa2*(Wvel2(i  ,j  ,kmm)+Wvel2(i  ,j  ,kp )))
      wwXRB=DDxb2*(Axb2*(Wvel2(ip  ,j   ,k   )+Wvel2(ipp ,j   ,k   ))+Bxb2*(Wvel2(i   ,j   ,k   )+Wvel2(ippp,j   ,k   )))
      wwXLB=DDxb2*(Axb2*(Wvel2(imm ,j   ,k   )+Wvel2(im  ,j   ,k   ))+Bxb2*(Wvel2(immm,j   ,k   )+Wvel2(i   ,j   ,k   )))
	  wwYRB=DDxb2*(Axb2*(Wvel2(i   ,jp  ,k   )+Wvel2(i   ,jpp ,k   ))+Bxb2*(Wvel2(i   ,j   ,k   )+Wvel2(i   ,jppp,k   )))
      wwYLB=DDxb2*(Axb2*(Wvel2(i   ,jmm ,k   )+Wvel2(i   ,jm  ,k   ))+Bxb2*(Wvel2(i   ,jmmm,k   )+Wvel2(i   ,j   ,k   )))
      wwZRB=DDxb2*(Axb2*(Wvel2(i   ,j   ,kp  )+Wvel2(i   ,j   ,kpp ))+Bxb2*(Wvel2(i   ,j   ,k   )+Wvel2(i   ,j   ,kppp)))
      wwZLB=DDxb2*(Axb2*(Wvel2(i   ,j   ,kmm )+Wvel2(i   ,j   ,km  ))+Bxb2*(Wvel2(i   ,j   ,kmmm)+Wvel2(i   ,j   ,k   )))

!	  IF (ABS(Uvel2(i,j,k)-Uvel(i,j,k)).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,U,U2',rank,i,j,k,Uvel2(i,j,k),Uvel(i,j,k)
!	  ENDIF
!	  IF (ABS(Vvel2(i,j,k)-Vvel(i,j,k)).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,V,V2',rank,i,j,k,Vvel2(i,j,k),Vvel(i,j,k)
!	  ENDIF
!	  IF (ABS(Wvel2(i,j,k)-Wvel(i,j,k)).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,W,W2',rank,i,j,k,Wvel2(i,j,k),Wvel(i,j,k)
!	  ENDIF
!	  IF (ABS(Uvel2(im,jm,km)-Uvel(im,jm,km)).gt.1e-14) THEN
!		write(*,*),'rank,im,jm,km,U,U2',rank,im,jm,km,Uvel2(im,jm,km),Uvel(im,jm,km)
!	  ENDIF
!	  IF (ABS(Vvel2(im,jm,km)-Vvel(im,jm,km)).gt.1e-14) THEN
!		write(*,*),'rank,im,jm,km,V,V2',rank,im,jm,km,Vvel2(im,jm,km),Vvel(im,jm,km)
!	  ENDIF
!	  IF (ABS(Wvel2(im,jm,km)-Wvel(im,jm,km)).gt.1e-14) THEN
!		write(*,*),'rank,im,jm,km,W,W2',rank,im,jm,km,Wvel2(im,jm,km),Wvel(im,jm,km)
!	  ENDIF
!	  IF (ABS(Uvel2(ip,jp,kp)-Uvel(ip,jp,kp)).gt.1e-14) THEN
!		write(*,*),'rank,ip,jp,kp,U,U2',rank,ip,jp,kp,Uvel2(ip,jp,kp),Uvel(ip,jp,kp)
!	  ENDIF
!	  IF (ABS(Vvel2(ip,jp,kp)-Vvel(ip,jp,kp)).gt.1e-14) THEN
!		write(*,*),'rank,ip,jp,kp,V,V2',rank,ip,jp,kp,Vvel2(ip,jp,kp),Vvel(ip,jp,kp)
!	  ENDIF
!	  IF (ABS(Wvel2(ip,jp,kp)-Wvel(ip,jp,kp)).gt.1e-14) THEN
!		write(*,*),'rank,ip,jp,kp,W,W2',rank,ip,jp,kp,Wvel2(ip,jp,kp),Wvel(ip,jp,kp)
!	  ENDIF	
!      uuRB=DDza*(Aza*(Uvel(i ,j ,k )+Uvel(i ,j ,kp))+Bza*(Uvel2(i ,j ,km )+Uvel2(i ,j ,kpp)))
!      uuLB=DDza*(Aza*(Uvel(im,j ,k )+Uvel(im,j ,kp))+Bza*(Uvel2(im,j ,km )+Uvel2(im,j ,kpp)))
!	  vvRB=DDza*(Aza*(Vvel(i ,j ,k )+Vvel(i ,j ,kp))+Bza*(Vvel2(i ,j ,km )+Vvel2(i ,j ,kpp)))
!      vvLB=DDza*(Aza*(Vvel(i ,jm,k )+Vvel(i ,jm,kp))+Bza*(Vvel2(i ,jm,km )+Vvel2(i ,jm,kpp)))
!      wwRB=DDza*(Aza*(Wvel(i ,j ,k )+Wvel(i ,j ,kp))+Bza*(Wvel2(i ,j ,km )+Wvel2(i ,j ,kpp)))
!      wwLB=DDza*(Aza*(Wvel(i ,j ,km)+Wvel(i ,j ,k ))+Bza*(Wvel2(i ,j ,kmm)+Wvel2(i ,j ,kp )))
!	  IF (ABS(uuRB-uuRA).gt.1e-14.or.ABS(uuLB-uuLA).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,uU',rank,i,j,k,uuRA,uuRB,uuLB,uuLA 
!	  ENDIF
!	  IF (ABS(vvRB-vvRA).gt.1e-14.or.ABS(vvLB-vvLA).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,uV',rank,i,j,k,vvRA,vvRB,vvLB,vvLA 
!	  ENDIF
!	  IF (ABS(wwRB-wwRA).gt.1e-14.or.ABS(wwLB-wwLA).gt.1e-14) THEN
!		write(*,*),'rank,i,j,k,uW',rank,i,j,k,wwRA,wwRB,wwLB,wwLA 
!	  ENDIF
	  
	  
      putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (
     1 facA*(uuRA*rhoip*Ru(i)*wwXRA-numdif/facA*ABS(uuRA*rhoip*Ru(i))* 
     1 (10.*(-Wvel2(i,j,k)+Wvel2(ip,j,k))-5.*(-Wvel2(im,j,k)+Wvel2(ipp,j,k))+(-Wvel2(imm,j,k)+Wvel2(ippp,j,k))) -
     1 (uuLA*rhoim*Ru(im)*wwXLA-numdif/facA*ABS(uuLA*rhoim*Ru(im))*
     1  (10.*(-Wvel2(im,j,k)+Wvel2(i,j,k))-5.*(-Wvel2(imm,j,k)+Wvel2(ip,j,k))+(-Wvel2(immm,j,k)+Wvel2(ipp,j,k)))) )
     1  / ( Rp(i)* dr(i) )
     +                         + 	 
     1 facB*(uuRB*rhoipp*Ru2(ip)*wwXRB-uuLB*rhoimm*Ru2(imm)*wwXLB)
     1  / ( Rp(i)* (Ru2(ip)-Ru2(imm)))	 
     +                         + 
     2 facA*(vvRA*rhojp*wwYRA-numdif/facA*ABS(vvRA*rhojp)*
     2  (10.*(-Wvel2(i,j,k)+Wvel2(i,jp,k))-5.*(-Wvel2(i,jm,k)+Wvel2(i,jpp,k))+(-Wvel2(i,jmm,k)+Wvel2(i,jppp,k))) -
     2 (vvLA*rhojm*wwYLA-numdif/facA*ABS(vvLA*rhojm)*
     2  (10.*(-Wvel2(i,jm,k)+Wvel2(i,j,k))-5.*(-Wvel2(i,jmm,k)+Wvel2(i,jp,k))+(-Wvel2(i,jmmm,k)+Wvel2(i,jpp,k)))) )
     2  / ( Rp(i) * (phivt2(rank*je+j)-phivt2(rank*je+jm)) )
     +                         + 
     2 facB*(vvRB*rhojpp*wwYRB-vvLB*rhojmm*wwYLB)
     2  / ( Rp(i) * (phivt2(rank*je+jp)-phivt2(rank*je+jmm)) )	 
     +                         +
     3 facA*(wwRA*rhokp*wwZRA-numdif/facA*ABS(wwRA*rhokp)*
     3  (10.*(-Wvel2(i,j,k)+Wvel2(i,j,kp))-5.*(-Wvel2(i,j,km)+Wvel2(i,j,kpp))+(-Wvel2(i,j,kmm)+Wvel2(i,j,kppp))) -
     3 (wwLA*rhokm*wwZLA-numdif/facA*ABS(wwLA*rhokm)*
     3  (10.*(-Wvel2(i,j,km)+Wvel2(i,j,k))-5.*(-Wvel2(i,j,kmm)+Wvel2(i,j,kp))+(-Wvel2(i,j,kmmm)+Wvel2(i,j,kpp)))) )
     3  *dzi
     +                         +
     3 facB*(wwRB*rhokpp*wwZRB-wwLB*rhokmm*wwZLB)*dzi/3.	 
     +                         )
           enddo
         enddo
      enddo
      return
      end

	  
      subroutine advecu_C4A6old(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,phiv)
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
	integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy	  
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1),numdif,dzi,phiv(0:j1)
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,rhoipp,rhoimm,rhojpp,rhojmm,rhokpp,rhokmm
	real Axa,Bxa,DDxa,Axb,Bxb,DDxb,Axa2,Bxa2,DDxa2,Axb2,Bxb2,DDxb2,facA,facB
	real uuRA,uuLA,vvRA,vvLA,wwRA,wwLA,uuRB,uuLB,vvRB,vvLB,wwRB,wwLB
	real uuXRA,uuXLA,uuYRA,uuYLA,uuZRA,uuZLA,uuXRB,uuXLB,uuYRB,uuYLB,uuZRB,uuZLB,order4yes

	real ubb(0:i1,0:k1),ubf(0:i1,0:k1),Uvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real Vvel2(-2:i1+2,-2:j1+2,-2:k1+2),Wvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real rho2(-1:i1+1,-1:j1+1,-1:k1+1),phivt2(-1:je*px+2),Rp2(-1:i1+1)

	
	
	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(1)-Rp(0)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim	
	phivt2(0:je*px+1)=phivt
	phivt2(-1)=phivt(0)-(phivt(1)-phivt(0)) !needed for periodicy sim
 	phivt2(je*px+1+1)=phivt(je*px+1)+(phivt(je*px+1)-phivt(je*px+1-1)) !needed for periodicy sim
	
	Uvel2(0:i1,0:j1,0:k1)=Uvel

c get stuff from other CPU's
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
	!Uvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Uvel2(-2:i1+2,-2:j1+2,1)+1./3.*Uvel2(-2:i1+2,-2:j1+2,2)
	!Uvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Uvel2(-2:i1+2,-2:j1+2,1)+2.*Uvel2(-2:i1+2,-2:j1+2,2)	
	Uvel2(-2:i1+2,-2:j1+2,k1+1)=Uvel2(-2:i1+2,-2:j1+2,k1)
	Uvel2(-2:i1+2,-2:j1+2,k1+2)=Uvel2(-2:i1+2,-2:j1+2,k1)

	
	Vvel2(0:i1,0:j1,0:k1)=Vvel

c get stuff from other CPU's
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
	!Vvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Vvel2(-2:i1+2,-2:j1+2,1)+1./3.*Vvel2(-2:i1+2,-2:j1+2,2)
	!Vvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Vvel2(-2:i1+2,-2:j1+2,1)+2.*Vvel2(-2:i1+2,-2:j1+2,2)
	Vvel2(-2:i1+2,-2:j1+2,k1+1)=Vvel2(-2:i1+2,-2:j1+2,k1)
	Vvel2(-2:i1+2,-2:j1+2,k1+2)=Vvel2(-2:i1+2,-2:j1+2,k1)

		Wvel2(0:i1,0:j1,0:k1)=Wvel

c get stuff from other CPU's
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
!	Wvel2(-2:i1+2,-2:j1+2,-1)=  -Wvel2(-2:i1+2,-2:j1+2,1)
!	Wvel2(-2:i1+2,-2:j1+2,-2)=  -Wvel2(-2:i1+2,-2:j1+2,2)	
	Wvel2(-2:i1+2,-2:j1+2,k1+1)=Wvel2(-2:i1+2,-2:j1+2,k1)
	Wvel2(-2:i1+2,-2:j1+2,k1+2)=Wvel2(-2:i1+2,-2:j1+2,k1)

	rho2(0:i1,0:j1,0:k1)=rho

c get stuff from other CPU's
	  call shiftf2(rho,ubf)
	  call shiftb2(rho,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = rho(i,0,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =rho(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	
	if (periodicx.eq.0.or.periodicx.eq.2) then
		rho2(-1,-1:j1+1,0:k1)=  rho2(1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(ie,-1:j1+1,0:k1)
	else
		rho2(-1,-1:j1+1,0:k1)=  rho2(ie-1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(2,-1:j1+1,0:k1)
	endif
	rho2(-1:i1+1,-1:j1+1,-1)=  rho2(-1:i1+1,-1:j1+1,0)
	rho2(-1:i1+1,-1:j1+1,k1+1)=rho2(-1:i1+1,-1:j1+1,k1)
	
	facA=4./3.
	facB=-1./6.	
	dzi=1./dz 
      do k=kb,ke
        do j=jb,je
          do  i=ib,ie
		  ipp=i+2
          ip=i+1
          im=i-1
		  imm=i-2
		  jpp=j+2
          jp=j+1
          jm=j-1
		  jmm=j-2
		  kpp=k+2
          kp=k+1
          km=k-1
		  kmm=k-2
		  rhojpp =0.25*(rho2(i,jp,k)+rho2(i,jpp,k)+rho2(ip,jp,k)+rho2(ip,jpp,k))
          rhojp =0.25*(rho2(i,j,k)+rho2(i,jp,k)+rho2(ip,j,k)+rho2(ip,jp,k))
          rhojm =0.25*(rho2(i,j,k)+rho2(i,jm,k)+rho2(ip,j,k)+rho2(ip,jm,k))
		  rhojmm =0.25*(rho2(i,jm,k)+rho2(i,jmm,k)+rho2(ip,jm,k)+rho2(ip,jmm,k))
		  rhokpp =0.25*(rho2(i,j,kp)+rho2(i,j,kpp)+rho2(ip,j,kp)+rho2(ip,j,kpp))
          rhokp =0.25*(rho2(i,j,k)+rho2(i,j,kp)+rho2(ip,j,k)+rho2(ip,j,kp))
          rhokm =0.25*(rho2(i,j,k)+rho2(i,j,km)+rho2(ip,j,k)+rho2(ip,j,km))
		  rhokmm =0.25*(rho2(i,j,km)+rho2(i,j,kmm)+rho2(ip,j,km)+rho2(ip,j,kmm))
		  IF (MAX(rhojpp,rhojp,rhojm,rhojmm,rhokpp,rhokp,rhokm,rhokmm).gt.1002.5) THEN
			write (*,*),'error max rho',MAX(rhojpp,rhojp,rhojm,rhojmm,rhokpp,rhokp,rhokm,rhokmm)
		  ENDIF
		  IF (MIN(rhojpp,rhojp,rhojm,rhojmm,rhokpp,rhokp,rhokm,rhokmm).lt.999.9999999999) THEN
			write (*,*),'error min rho',MIN(rhojpp,rhojp,rhojm,rhojmm,rhokpp,rhokp,rhokm,rhokmm)
		  ENDIF			
      putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (
     1 facA*(Rp2(ip)*(Uvel2(i,j,k)+Uvel2(ip,j,k))*(Uvel2(i,j,k)+Uvel2(ip,j,k))*rho2(ip,j,k) -
     1  Rp2(i )*(Uvel2(i,j,k)+Uvel2(im,j,k))*(Uvel2(i,j,k)+Uvel2(im,j,k))*rho2(i ,j,k)  )
     1  / ( Ru(i) * ( Rp2(ip)-Rp2(i) ) )
     +                         + 
     1 facB*(Rp2(i)*(Uvel2(im,j,k)+Uvel2(i,j,k))*(Uvel2(im,j,k)+Uvel2(i,j,k))*rho2(i,j,k) -
     1  Rp2(im )*(Uvel2(im,j,k)+Uvel2(imm,j,k))*(Uvel2(im,j,k)+Uvel2(imm,j,k))*rho2(im ,j,k)  )
     1  / ( Ru(im) * ( Rp2(i)-Rp2(im) ) )
     +                         + 
     1 facB*(Rp2(ipp)*(Uvel2(ip,j,k)+Uvel2(ipp,j,k))*(Uvel2(ip,j,k)+Uvel2(ipp,j,k))*rho2(ipp,j,k) -
     1  Rp2(ip )*(Uvel2(ip,j,k)+Uvel2(i,j,k))*(Uvel2(ip,j,k)+Uvel2(i,j,k))*rho2(ip ,j,k)  )
     1  / ( Ru(ip) * ( Rp2(ipp)-Rp2(ip) ) )
     +                         + 	 
     2 facA*(      (Vvel2(i,j,k) +Vvel2(ip,j,k) )*(Uvel2(i,j,k)+Uvel2(i,jp,k))*rhojp -
     2        (Vvel2(i,jm,k)+Vvel2(ip,jm,k))*(Uvel2(i,j,k)+Uvel2(i,jm,k))*rhojm  )
     2  / ( Ru(i) * (phivt2(rank*je+j)-phivt2(rank*je+jm)) )
     +                         + 	 
     2 facB*(      (Vvel2(i,jm,k) +Vvel2(ip,jm,k) )*(Uvel2(i,jm,k)+Uvel2(i,j,k))*rhojm -
     2        (Vvel2(i,jmm,k)+Vvel2(ip,jmm,k))*(Uvel2(i,jm,k)+Uvel2(i,jmm,k))*rhojmm  )
     2  / ( Ru(i) * (phivt2(rank*je+jm)-phivt2(rank*je+jmm)) )
     +                         + 	 
     2 facB*(      (Vvel2(i,jp,k) +Vvel2(ip,jp,k) )*(Uvel2(i,jp,k)+Uvel2(i,jpp,k))*rhojpp -
     2        (Vvel2(i,j,k)+Vvel2(ip,j,k))*(Uvel2(i,jp,k)+Uvel2(i,j,k))*rhojp  )
     2  / ( Ru(i) * (phivt2(rank*je+j)-phivt2(rank*je+jm)) )	 
     +                         +
     3 facA*(      (Wvel2(i,j,k) +Wvel2(ip,j,k) )*(Uvel2(i,j,k)+Uvel2(i,j,kp)) *rhokp -
     3        (Wvel2(i,j,km)+Wvel2(ip,j,km))*(Uvel2(i,j,k)+Uvel2(i,j,km))*rhokm  )
     3  / ( dz )
     +                         +
     3 facB*(      (Wvel2(i,j,km) +Wvel2(ip,j,km) )*(Uvel2(i,j,km)+Uvel2(i,j,k))*rhokm -
     3        (Wvel2(i,j,kmm)+Wvel2(ip,j,kmm))*(Uvel2(i,j,km)+Uvel2(i,j,kmm))*rhokmm  )
     3  / ( dz )
     +                         +
     3 facB*(      (Wvel2(i,j,kp) +Wvel2(ip,j,kp) )*(Uvel2(i,j,kp)+Uvel2(i,j,kpp))*rhokpp -
     3        (Wvel2(i,j,k)+Wvel2(ip,j,k))*(Uvel2(i,j,kp)+Uvel2(i,j,k))*rhokp  )
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


      subroutine advecv_C4A6old(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,phiv)
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
	  integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),Ru2(-1:i1+1)
      real rho(0:i1,0:j1,0:k1),dzi,phiv(0:j1)
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,rhoipp,rhoimm,rhojpp,rhojmm,rhokpp,rhokmm
	real Axa,Bxa,DDxa,Axb,Bxb,DDxb,Axa2,Bxa2,DDxa2,Axb2,Bxb2,DDxb2,numdif,facA,facB,order4yes
	

	real ubb(0:i1,0:k1),ubf(0:i1,0:k1),Vvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real Uvel2(-2:i1+2,-2:j1+2,-2:k1+2),Wvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real uuRA,uuLA,vvRA,vvLA,wwRA,wwLA,uuRB,uuLB,vvRB,vvLB,wwRB,wwLB
	real vvXRA,vvXLA,vvYRA,vvYLA,vvZRA,vvZLA,vvXRB,vvXLB,vvYRB,vvYLB,vvZRB,vvZLB
	real rho2(-1:i1+1,-1:j1+1,-1:k1+1),phivt2(-1:je*px+2)

	
	
	Ru2(0:i1)=Ru
	Ru2(-1)=Ru(0)-(Ru(1)-Ru(0)) !needed for periodicx sim
 	Ru2(i1+1)=Ru(i1)+(Ru(i1)-Ru(i1-1)) !needed for periodicx sim	
	phivt2(0:je*px+1)=phivt
	phivt2(-1)=phivt(0)-(phivt(1)-phivt(0)) !needed for periodicy sim
 	phivt2(je*px+1+1)=phivt(je*px+1)+(phivt(je*px+1)-phivt(je*px+1-1)) !needed for periodicy sim
	
	Uvel2(0:i1,0:j1,0:k1)=Uvel

c get stuff from other CPU's
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
	!Uvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Uvel2(-2:i1+2,-2:j1+2,1)+1./3.*Uvel2(-2:i1+2,-2:j1+2,2)
	!Uvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Uvel2(-2:i1+2,-2:j1+2,1)+2.*Uvel2(-2:i1+2,-2:j1+2,2)	
	Uvel2(-2:i1+2,-2:j1+2,k1+1)=Uvel2(-2:i1+2,-2:j1+2,k1)
	Uvel2(-2:i1+2,-2:j1+2,k1+2)=Uvel2(-2:i1+2,-2:j1+2,k1)

	
	Vvel2(0:i1,0:j1,0:k1)=Vvel

c get stuff from other CPU's
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
	!Vvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Vvel2(-2:i1+2,-2:j1+2,1)+1./3.*Vvel2(-2:i1+2,-2:j1+2,2)
	!Vvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Vvel2(-2:i1+2,-2:j1+2,1)+2.*Vvel2(-2:i1+2,-2:j1+2,2)	
	Vvel2(-2:i1+2,-2:j1+2,k1+1)=Vvel2(-2:i1+2,-2:j1+2,k1)
	Vvel2(-2:i1+2,-2:j1+2,k1+2)=Vvel2(-2:i1+2,-2:j1+2,k1)

		Wvel2(0:i1,0:j1,0:k1)=Wvel

c get stuff from other CPU's
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
!	Wvel2(-2:i1+2,-2:j1+2,-1)=  -Wvel2(-2:i1+2,-2:j1+2,1)
!	Wvel2(-2:i1+2,-2:j1+2,-2)=  -Wvel2(-2:i1+2,-2:j1+2,2)
	Wvel2(-2:i1+2,-2:j1+2,k1+1)=Wvel2(-2:i1+2,-2:j1+2,k1)
	Wvel2(-2:i1+2,-2:j1+2,k1+2)=Wvel2(-2:i1+2,-2:j1+2,k1)

	rho2(0:i1,0:j1,0:k1)=rho

c get stuff from other CPU's
	  call shiftf2(rho,ubf)
	  call shiftb2(rho,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = rho(i,0,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =rho(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	
	if (periodicx.eq.0.or.periodicx.eq.2) then
		rho2(-1,-1:j1+1,0:k1)=  rho2(1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(ie,-1:j1+1,0:k1)
	else
		rho2(-1,-1:j1+1,0:k1)=  rho2(ie-1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(2,-1:j1+1,0:k1)
	endif
	rho2(-1:i1+1,-1:j1+1,-1)=  rho2(-1:i1+1,-1:j1+1,0)
	rho2(-1:i1+1,-1:j1+1,k1+1)=rho2(-1:i1+1,-1:j1+1,k1)
	
	facA=4./3.
	facB=-1./6.	
	dzi=1./dz 
      do k=kb,ke
        do j=jb,je
          do  i=ib,ie
		  ipp=i+2
          ip=i+1
          im=i-1
		  imm=i-2
		  jpp=j+2
          jp=j+1
          jm=j-1
		  jmm=j-2
		  kpp=k+2
          kp=k+1
          km=k-1
		  kmm=k-2
		  rhoipp =0.25*(rho2(ip,j,k)+rho2(ipp,j,k)+rho2(ip,jp,k)+rho2(ipp,jp,k))
          rhoip =0.25*(rho2(i,j,k)+rho2(ip,j,k)+rho2(i,jp,k)+rho2(ip,jp,k))
          rhoim =0.25*(rho2(i,j,k)+rho2(im,j,k)+rho2(i,jp,k)+rho2(im,jp,k))
		  rhoimm =0.25*(rho2(im,j,k)+rho2(imm,j,k)+rho2(im,jp,k)+rho2(imm,jp,k))
          rhokpp =0.25*(rho2(i,j,kp)+rho2(i,j,kpp)+rho2(i,jp,kp)+rho2(i,jp,kpp))		  
          rhokp =0.25*(rho2(i,j,k)+rho2(i,j,kp)+rho2(i,jp,k)+rho2(i,jp,kp))
          rhokm =0.25*(rho2(i,j,k)+rho2(i,j,km)+rho2(i,jp,k)+rho2(i,jp,km))
		  rhokmm =0.25*(rho2(i,j,km)+rho2(i,j,kmm)+rho2(i,jp,km)+rho2(i,jp,kmm))
      putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (
     1 facA*( Ru2(i ) * Ru2(i ) *rhoip*
     1     (Uvel2(i,j,k) +Uvel2(i,jp,k) )*(Vvel2(i,j,k)+Vvel2(ip,j,k)) -
     1   Ru2(im) * Ru2(im) *rhoim*
     1     (Uvel2(im,j,k)+Uvel2(im,jp,k))*(Vvel2(i,j,k)+Vvel2(im,j,k))  )
     1  / ( Rp(i) * Rp(i) * (Ru2(i)-Ru2(im)) )	 
     +                         + 
     1 facB*( Ru2(im ) * Ru2(im ) *rhoim*
     1     (Uvel2(im,j,k) +Uvel2(im,jp,k) )*(Vvel2(im,j,k)+Vvel2(i,j,k)) -
     1   Ru2(imm) * Ru2(imm) *rhoimm*
     1     (Uvel2(imm,j,k)+Uvel2(imm,jp,k))*(Vvel2(im,j,k)+Vvel2(imm,j,k))  )
     1  / ( Rp(im) * Rp(im) * (Ru2(im)-Ru2(imm)) )	 
     +                         + 
     1 facB*( Ru2(ip ) * Ru2(ip ) *rhoipp*
     1     (Uvel2(ip,j,k) +Uvel2(ip,jp,k) )*(Vvel2(ip,j,k)+Vvel2(ipp,j,k)) -
     1   Ru2(i) * Ru2(i) *rhoip*
     1     (Uvel2(i,j,k)+Uvel2(i,jp,k))*(Vvel2(ip,j,k)+Vvel2(i,j,k))  )
     1  / ( Rp(ip) * Rp(ip) * (Ru2(ip)-Ru2(i)) )	 
     +                         + 	 
     2 facA*(   (Vvel2(i,j,k) +Vvel2(i,jp,k) )*(Vvel2(i,j,k)+Vvel2(i,jp,k))*rho2(i,jp,k) -
     2     (Vvel2(i,jm,k)+Vvel2(i,j,k)  )*(Vvel2(i,j,k)+Vvel2(i,jm,k))*rho2(i,j,k )  )
     2  / ( Rp(i) * (phivt2(rank*je+j)-phivt2(rank*je+jm)) )	 
     +                         + 	 
     2 facB*(   (Vvel2(i,jm,k) +Vvel2(i,j,k) )*(Vvel2(i,jm,k)+Vvel2(i,j,k))*rho2(i,j,k) -
     2     (Vvel2(i,jmm,k)+Vvel2(i,jm,k)  )*(Vvel2(i,jm,k)+Vvel2(i,jmm,k))*rho2(i,jm,k )  )
     2  / ( Rp(i) * (phivt2(rank*je+jm)-phivt2(rank*je+jmm)) )	 
     +                         + 	 
     2 facB*(   (Vvel2(i,jp,k) +Vvel2(i,jpp,k) )*(Vvel2(i,jp,k)+Vvel2(i,jpp,k))*rho2(i,jpp,k) -
     2     (Vvel2(i,j,k)+Vvel2(i,jp,k)  )*(Vvel2(i,jp,k)+Vvel2(i,j,k))*rho2(i,jp,k )  )
     2  / ( Rp(i) * (phivt2(rank*je+jp)-phivt2(rank*je+j)) )	 	 
     +                         +
     3 facA*(   (Wvel2(i,j,k) +Wvel2(i,jp,k) )*(Vvel2(i,j,k)+Vvel2(i,j,kp))*rhokp -
     3     (Wvel2(i,j,km)+Wvel2(i,jp,km))*(Vvel2(i,j,k)+Vvel2(i,j,km))*rhokm  )
     3  / ( dz )
     +                         +
     3 facB*(   (Wvel2(i,j,km) +Wvel2(i,jp,km) )*(Vvel2(i,j,km)+Vvel2(i,j,k))*rhokm -
     3     (Wvel2(i,j,kmm)+Wvel2(i,jp,kmm))*(Vvel2(i,j,km)+Vvel2(i,j,kmm))*rhokmm  )
     3  / ( dz )
     +                         +
     3 facB*(   (Wvel2(i,j,kp) +Wvel2(i,jp,kp) )*(Vvel2(i,j,kp)+Vvel2(i,j,kpp))*rhokpp -
     3     (Wvel2(i,j,k)+Wvel2(i,jp,k))*(Vvel2(i,j,kp)+Vvel2(i,j,k))*rhokp  )
     3  / ( dz )	 
     +                         )
           enddo
        enddo
      enddo
      return
      end

      subroutine advecw_C4A6old(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phivt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,numdif,periodicx,periodicy,phiv)
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
	  integer kpp,kppp,kmm,kmmm,jpp,jppp,jmm,jmmm,ipp,ippp,imm,immm,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phivt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),Ru2(-1:i1+1)
      real rho(0:i1,0:j1,0:k1),dzi,phiv(0:j1)
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm,rhoipp,rhoimm,rhojpp,rhojmm,rhokpp,rhokmm
	real Axa,Bxa,DDxa,Axb,Bxb,DDxb,Axa2,Bxa2,DDxa2,Axb2,Bxb2,DDxb2,numdif,facA,facB
	
	real ubb(0:i1,0:k1),ubf(0:i1,0:k1),Wvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real Uvel2(-2:i1+2,-2:j1+2,-2:k1+2),Vvel2(-2:i1+2,-2:j1+2,-2:k1+2)
	real uuRA,uuLA,vvRA,vvLA,wwRA,wwLA,uuRB,uuLB,vvRB,vvLB,wwRB,wwLB
	real wwXRA,wwXLA,wwYRA,wwYLA,wwZRA,wwZLA,wwXRB,wwXLB,wwYRB,wwYLB,wwZRB,wwZLB
	real rho2(-1:i1+1,-1:j1+1,-1:k1+1),phivt2(-1:je*px+2)

	
	
	Ru2(0:i1)=Ru
	Ru2(-1)=Ru(0)-(Ru(1)-Ru(0)) !needed for periodicx sim
 	Ru2(i1+1)=Ru(i1)+(Ru(i1)-Ru(i1-1)) !needed for periodicx sim	
	phivt2(0:je*px+1)=phivt
	phivt2(-1)=phivt(0)-(phivt(1)-phivt(0)) !needed for periodicy sim
 	phivt2(je*px+1+1)=phivt(je*px+1)+(phivt(je*px+1)-phivt(je*px+1-1)) !needed for periodicy sim
	
	Uvel2(0:i1,0:j1,0:k1)=Uvel

c get stuff from other CPU's
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
	!Uvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Uvel2(-2:i1+2,-2:j1+2,1)+1./3.*Uvel2(-2:i1+2,-2:j1+2,2)
	!Uvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Uvel2(-2:i1+2,-2:j1+2,1)+2.*Uvel2(-2:i1+2,-2:j1+2,2)	
	Uvel2(-2:i1+2,-2:j1+2,k1+1)=Uvel2(-2:i1+2,-2:j1+2,k1)
	Uvel2(-2:i1+2,-2:j1+2,k1+2)=Uvel2(-2:i1+2,-2:j1+2,k1)

	
	Vvel2(0:i1,0:j1,0:k1)=Vvel

c get stuff from other CPU's
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
	!Vvel2(-2:i1+2,-2:j1+2,-1)=  2./3.*Vvel2(-2:i1+2,-2:j1+2,1)+1./3.*Vvel2(-2:i1+2,-2:j1+2,2)
	!Vvel2(-2:i1+2,-2:j1+2,-2)=  -1.*Vvel2(-2:i1+2,-2:j1+2,1)+2.*Vvel2(-2:i1+2,-2:j1+2,2)	
	Vvel2(-2:i1+2,-2:j1+2,k1+1)=Vvel2(-2:i1+2,-2:j1+2,k1)
	Vvel2(-2:i1+2,-2:j1+2,k1+2)=Vvel2(-2:i1+2,-2:j1+2,k1)

		Wvel2(0:i1,0:j1,0:k1)=Wvel

c get stuff from other CPU's
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
!	Wvel2(-2:i1+2,-2:j1+2,-1)=  -Wvel2(-2:i1+2,-2:j1+2,1)
!	Wvel2(-2:i1+2,-2:j1+2,-2)=  -Wvel2(-2:i1+2,-2:j1+2,2)
	Wvel2(-2:i1+2,-2:j1+2,k1+1)=Wvel2(-2:i1+2,-2:j1+2,k1)
	Wvel2(-2:i1+2,-2:j1+2,k1+2)=Wvel2(-2:i1+2,-2:j1+2,k1)

	rho2(0:i1,0:j1,0:k1)=rho

c get stuff from other CPU's
	  call shiftf2(rho,ubf)
	  call shiftb2(rho,ubb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = rho(i,0,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =rho(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=0,k1
		   do i=0,i1 !passing inflow bc at i=0 and i=i1 correctly
		   rho2(i,-1,k) = Ubf(i,k)
		   rho2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	endif
	
	if (periodicx.eq.0.or.periodicx.eq.2) then
		rho2(-1,-1:j1+1,0:k1)=  rho2(1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(ie,-1:j1+1,0:k1)
	else
		rho2(-1,-1:j1+1,0:k1)=  rho2(ie-1,-1:j1+1,0:k1)
		rho2(i1+1,-1:j1+1,0:k1)=rho2(2,-1:j1+1,0:k1)
	endif
	rho2(-1:i1+1,-1:j1+1,-1)=  rho2(-1:i1+1,-1:j1+1,0)
	rho2(-1:i1+1,-1:j1+1,k1+1)=rho2(-1:i1+1,-1:j1+1,k1)
	
	facA=4./3.
	facB=-1./6.	
	dzi=1./dz 
      do k=kb,ke
        do j=jb,je
          do i=ib,ie
		  ipp=i+2
          ip=i+1
          im=i-1
		  imm=i-2
		  jpp=j+2
          jp=j+1
          jm=j-1
		  jmm=j-2
		  kpp=k+2
          kp=k+1
          km=k-1
		  kmm=k-2
		  rhoipp =0.25*(rho2(ip,j,k)+rho2(ip,j,kp)+rho2(ipp,j,k)+rho2(ipp,j,kp))
          rhoip =0.25*(rho2(i,j,k)+rho2(i,j,kp)+rho2(ip,j,k)+rho2(ip,j,kp))
          rhoim =0.25*(rho2(i,j,k)+rho2(i,j,kp)+rho2(im,j,k)+rho2(im,j,kp))
		  rhoimm =0.25*(rho2(im,j,k)+rho2(im,j,kp)+rho2(imm,j,k)+rho2(imm,j,kp))
		  rhojpp =0.25*(rho2(i,jp,k)+rho2(i,jp,kp)+rho2(i,jpp,k)+rho2(i,jpp,kp))
          rhojp =0.25*(rho2(i,j,k)+rho2(i,j,kp)+rho2(i,jp,k)+rho2(i,jp,kp))
          rhojm =0.25*(rho2(i,j,k)+rho2(i,j,kp)+rho2(i,jm,k)+rho2(i,jm,kp))
		  rhojmm =0.25*(rho2(i,jm,k)+rho2(i,jm,kp)+rho2(i,jmm,k)+rho2(i,jmm,kp))
      putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (
     1 facA*(Ru2(i )*(Uvel2(i,j,k) +Uvel2(i,j,kp) )*(Wvel2(i,j,k)+Wvel2(ip,j,k))*rhoip -
     1  Ru2(im)*(Uvel2(im,j,k)+Uvel2(im,j,kp))*(Wvel2(i,j,k)+Wvel2(im,j,k))*rhoim )
     1  / ( Rp(i) * (Ru2(i)-Ru2(im)) )
     +                         + 	 
     1 facB*(Ru2(im )*(Uvel2(im,j,k) +Uvel2(im,j,kp) )*(Wvel2(im,j,k)+Wvel2(i,j,k))*rhoim -
     1  Ru2(imm)*(Uvel2(imm,j,k)+Uvel2(imm,j,kp))*(Wvel2(im,j,k)+Wvel2(imm,j,k))*rhoimm )
     1  / ( Rp(im) * (Ru2(im)-Ru2(imm)) )
     +                         + 	 
     1 facB*(Ru2(ip )*(Uvel2(ip,j,k) +Uvel2(ip,j,kp) )*(Wvel2(ip,j,k)+Wvel2(ipp,j,k))*rhoipp -
     1  Ru2(i )*(Uvel2(i,j,k)+Uvel2(i,j,kp))*(Wvel2(ip,j,k)+Wvel2(i,j,k))*rhoip )
     1  / ( Rp(ip) * (Ru2(ip)-Ru2(i)) )	 
     +                         + 
     2 facA*(   (Vvel2(i,j,k) +Vvel2(i,j,kp) )*(Wvel2(i,j,k)+Wvel2(i,jp,k))*rhojp -
     2     (Vvel2(i,jm,k)+Vvel2(i,jm,kp))*(Wvel2(i,j,k)+Wvel2(i,jm,k))*rhojm  )
     2  / ( Rp(i) * (phivt2(rank*je+j)-phivt2(rank*je+jm)) )
     +                         + 
     2 facB*(   (Vvel2(i,jm,k) +Vvel2(i,jm,kp) )*(Wvel2(i,jm,k)+Wvel2(i,j,k))*rhojm -
     2     (Vvel2(i,jmm,k)+Vvel2(i,jmm,kp))*(Wvel2(i,jm,k)+Wvel2(i,jmm,k))*rhojmm  )
     2  / ( Rp(i) * (phivt2(rank*je+jm)-phivt2(rank*je+jmm)) )
     +                         + 
     2 facB*(   (Vvel2(i,jp,k) +Vvel2(i,jp,kp) )*(Wvel2(i,jp,k)+Wvel2(i,jpp,k))*rhojpp -
     2     (Vvel2(i,j,k)+Vvel2(i,j,kp))*(Wvel2(i,jp,k)+Wvel2(i,j,k))*rhojp  )
     2  / ( Rp(i) * (phivt2(rank*je+jp)-phivt2(rank*je+j)) )
     +                         +
     3 facA*(   (Wvel2(i,j,k) +Wvel2(i,j,kp) )*(Wvel2(i,j,k)+Wvel2(i,j,kp))*rho2(i,j,kp) -
     3     (Wvel2(i,j,k) +Wvel2(i,j,km) )*(Wvel2(i,j,k)+Wvel2(i,j,km))*rho2(i,j,k )  )
     3  / ( dz )
     +                         +
     3 facB*(   (Wvel2(i,j,km) +Wvel2(i,j,k) )*(Wvel2(i,j,km)+Wvel2(i,j,k))*rho2(i,j,k) -
     3     (Wvel2(i,j,km) +Wvel2(i,j,kmm) )*(Wvel2(i,j,km)+Wvel2(i,j,kmm))*rho2(i,j,km )  )
     3  / ( dz )
     +                         +
     3 facB*(   (Wvel2(i,j,kp) +Wvel2(i,j,kpp) )*(Wvel2(i,j,kp)+Wvel2(i,j,kpp))*rho2(i,j,kpp) -
     3     (Wvel2(i,j,kp) +Wvel2(i,j,k) )*(Wvel2(i,j,kp)+Wvel2(i,j,k))*rho2(i,j,kp )  )
     3  / ( dz )	 
     +                         )
           enddo
         enddo
      enddo
      return
      end	  

