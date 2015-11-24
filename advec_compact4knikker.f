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



      subroutine advecu_COM4(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)

      implicit none

	include 'mpif.h'
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
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
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
      integer  rank,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1),uu(0:i1,0:j1,0:k1)

      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm
      real wall_value
!      real aax_0i1(0:i1),bbx_0i1(0:i1),ccx_0i1(0:i1),ddx_0i1(0:i1)
      real aax2_0ie(0:ie),bbx2_0ie(0:ie),ccx2_0ie(0:ie),ddx2_0ie(0:ie)
      real aax_0ie(0:ie),bbx_0ie(0:ie),ccx_0ie(0:ie),ddx_0ie(0:ie)
!      real aax_1i1(1:i1),bbx_1i1(1:i1),ccx_1i1(1:i1),ddx_1i1(1:i1)
      real aax_1ie(1:ie),bbx_1ie(1:ie),ccx_1ie(1:ie),ddx_1ie(1:ie)

!      real aax_0ie1(0:i1),bbx_0ie1(0:i1),ccx_0ie1(0:i1),ddx_0ie1(0:i1),dfx1(0:i1)

      real aaz_0ke(0:ke),bbz_0ke(0:ke),ccz_0ke(0:ke),ddz_0ke(0:ke)
      real aaz_1ke(1:ke),bbz_1ke(1:ke),ccz_1ke(1:ke),ddz_1ke(1:ke)

      real ruu(0:i1),u_r(0:i1),dfx(0:ie),ru_u(0:ie,0:j1,0:k1),w_II(0:ie,0:j1,0:k1)
	real ru_II(0:ke),ruw(0:ke),dfz(1:ke)
	real v_I(0:ie,0:j1,0:k1)
	real putout_T(1:ie,1:je*px,1:ke/px),v_I_T(1:ie,0:je*px+1,1:ke/px),ru_u_T(1:ie,0:je*px+1,1:ke/px)
	real putouty(1:ie,1:je,1:ke)

	real dfy_T(1:je*px),ru_I_T(0:je*px),ruv_I_T(0:je*px)
      real aay_0jepx(0:je*px),bby_0jepx(0:je*px),ccy_0jepx(0:je*px),ddy_0jepx(0:je*px)
      real aay_1jepx(1:je*px),bby_1jepx(1:je*px),ccy_1jepx(1:je*px),ddy_1jepx(1:je*px)
	
      integer ileng,ierr,itag,status(MPI_STATUS_SIZE)
	real const_1_6,const_1_22,const_12_11,const_2_3
	putout=0.

	const_1_6=1./6.
	const_1_22=1./22.
	const_12_11=12./11.
	const_2_3=2./3.

       	bbx_1ie=1.
	aax_1ie=const_1_6 !3./10.
      	aax_1ie(1)=0.
	aax_1ie(ie)=0. !0. !1. !0.
	ccx_1ie=const_1_6 !3./10.
	ccx_1ie(1)=0. !0. !0. !1. !0. !1.
	ccx_1ie(ie)=0.

      	bbx2_0ie=1.
	aax2_0ie=const_1_6 !3./10.
      	aax2_0ie(0)=0.
	aax2_0ie(ie)=0. !5. !0.
	ccx2_0ie=const_1_6 !3./10.
	ccx2_0ie(0)=0. !5. !0.
	ccx2_0ie(ie)=0.

	bbx_0ie=1.
	aax_0ie=const_1_22 !9./62.
	aax_0ie(0)=0.
	aax_0ie(ie)=3.!*3./8. !331./15. !23. !0. !const_1_22 !0. !1./23. !	1./23. !0. !2. !15.
	ccx_0ie=const_1_22 !9./62.
	ccx_0ie(0)=3.!*3./8. !331./15. !23. !const_1_22 !23. !0. !1./23. !1./23. !0. !2. !15.	
	ccx_0ie(ie)=0.	


      do k=0,k1 !kb,ke
        do j=0,j1
!		!! interpolate_x rho -> u-loc
		ddx2_0ie(0)=0.5*(rho(0,j,k)+rho(1,j,k)) 
		ddx2_0ie(ie)=0.5*(rho(ie,j,k)+rho(i1,j,k)) 
		do i=1,ie-1
			ddx2_0ie(i)=const_2_3*(rho(i+1,j,k)+rho(i,j,k))
		enddo
		CALL solve_tridiag(ru_u(0:ie,j,k),aax2_0ie,bbx2_0ie,ccx2_0ie,ddx2_0ie,i1)	! interpolate rho -> u-loc -> 0:i1
		ru_u(0:ie,j,k)=ru_u(0:ie,j,k)*Uvel(0:ie,j,k)  
	enddo
      enddo

	uu=Uvel
      do k=kb,ke
        do j=jb,je
!		!! interpolate_x u -> rho-loc
	       ddx_1ie(1)=0.5*(uu(0,j,k)+uu(1,j,k)) 
	       ddx_1ie(ie)=0.5*(uu(ie-1,j,k)+uu(ie,j,k))  
		do i=2,ie-1
			ddx_1ie(i)=const_2_3*(uu(i-1,j,k)+uu(i,j,k))
		enddo
		CALL solve_tridiag(u_r(1:ie),aax_1ie,bbx_1ie,ccx_1ie,ddx_1ie,ie)	
		! derivative d(Rrhouu)Rdr
		ruu(1:ie)=Rp(1:ie)*u_r(1:ie)*u_r(1:ie)*rho(1:ie,j,k) 

              ddx_0ie(0)=8./3.*(-Ru(0)*uu(0,j,k)*uu(0,j,k)*rho(0,j,k)+ruu(2))/(Ru(0)*(Rp(1)-Rp(0)))
              ddx_0ie(ie)=8./3.*(Ru(ie)*uu(ie,j,k)*uu(ie,j,k)*rho(i1,j,k)-ruu(ie-1))/(Ru(ie)*(Rp(i1)-Rp(ie)))
		do i=1,ie-1
			ddx_0ie(i)=const_12_11*(ruu(i+1)-ruu(i))/(Ru(i)*(Rp(i+1)-Rp(i)))
		enddo
		CALL solve_tridiag(dfx(0:ie),aax_0ie,bbx_0ie,ccx_0ie,ddx_0ie,ie+1)   ! df on Uvel(i,j,k) location	-> 1:ie 
		putout(1:ie,j,k)=putout(1:ie,j,k)-dfx(1:ie) ! derivative: d(ruu)/rdr
	enddo
      enddo

      do k=0,ke ! for correct bc on w_II also k=0 is included 
        do j=jb,je
		!! interpolate_x Wvel -> II-loc 
		ddx2_0ie(0)=0.5*(Wvel(0,j,k)+Wvel(1,j,k)) 
		ddx2_0ie(ie)=0.5*(Wvel(ie,j,k)+Wvel(i1,j,k)) 
		do i=1,ie-1
			ddx2_0ie(i)=const_2_3*(Wvel(i+1,j,k)+Wvel(i,j,k))
		enddo
		CALL solve_tridiag(w_II(0:ie,j,k),aax2_0ie,bbx2_0ie,ccx2_0ie,ddx2_0ie,i1)	! interpolate W -> II-loc -> 0:i1 (with i1 as not used value)
	enddo
      enddo


      	bbz_0ke=1.
	aaz_0ke=const_1_6 !3./10.
      	aaz_0ke(0)=0.
	aaz_0ke(ke)=0. !5. !0.
	ccz_0ke=const_1_6 !3./10.
	ccz_0ke(0)=0. !5. !0.
	ccz_0ke(ke)=0.
      	bbz_1ke=1.
	aaz_1ke=const_1_22 !9./62.
      	aaz_1ke(1)=0.
	aaz_1ke(ke)=-1. !-1.!*26./24. !0. !const_1_22 !1./11. 
	ccz_1ke=const_1_22 !9./62.
	ccz_1ke(1)=-1. !-1.!*26./24. !0. !const_1_22 !1./11.
	ccz_1ke(ke)=0.

      do i=ib,ie
        do j=jb,je
		! interpolate_z ru_u -> II-loc
		ddz_0ke(0)=0.5*(ru_u(i,j,0)+ru_u(i,j,1))
		ddz_0ke(ke)=0.5*(ru_u(i,j,ke)+ru_u(i,j,k1))
		do k=1,ke-1
			ddz_0ke(k)=const_2_3*(ru_u(i,j,k+1)+ru_u(i,j,k))
		enddo
		CALL solve_tridiag(ru_II(0:ke),aaz_0ke,bbz_0ke,ccz_0ke,ddz_0ke,ke+1)	! ru_II -> II-loc 0:ke
		ruw(0:ke)=ru_II(0:ke)*w_II(i,j,0:ke)

		ddz_1ke(1)=(-ruw(0)+2.*ruw(1)-ruw(2))/dz 
		ddz_1ke(ke)=(+ruw(ke)-2.*ruw(ke-1)+ruw(ke-2))/dz 
		do k=2,ke-1 
			ddz_1ke(k)=const_12_11*(ruw(k)-ruw(k-1))/dz
		enddo
		CALL solve_tridiag(dfz(1:ke),aaz_1ke,bbz_1ke,ccz_1ke,ddz_1ke,ke)   ! dfz on Uvel(i,j,k) location	-> 1:ke
		putout(i,j,1:ke)=putout(i,j,1:ke)-dfz(1:ke)
	enddo
      enddo

!!	interpolate_x Vvel -> I-loc
      do k=kb,ke ! for correct bc on v_I also j=0 is included
        do j=0,je
		!! interpolate_x Vvel -> I-loc		
		ddx2_0ie(0)=0.5*(Vvel(0,j,k)+Vvel(1,j,k)) 
		ddx2_0ie(ie)=0.5*(Vvel(ie,j,k)+Vvel(i1,j,k)) 
		do i=1,ie-1
			ddx2_0ie(i)=const_2_3*(Vvel(i+1,j,k)+Vvel(i,j,k))
		enddo
		CALL solve_tridiag(v_I(0:ie,j,k),aax2_0ie,bbx2_0ie,ccx2_0ie,ddx2_0ie,i1)	! interpolate_x V -> I-loc -> 0:i1 (with 0 and i1 as not used value)
	enddo
      enddo	

	call t2np(ru_u(1:ie,1:je,1:ke),ru_u_T(1:ie,1:je*px,1:ke/px))
	call t2np(v_I(1:ie,1:je,1:ke),v_I_T(1:ie,1:je*px,1:ke/px))

	!! also pass over boundaries at j=0 :
	IF (rank.eq.0) THEN
	  do i=1,px-1
      		call mpi_send(v_I(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i,MPI_COMM_WORLD,status,ierr)
      		call mpi_send(ru_u(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+100,MPI_COMM_WORLD,status,ierr)
	  enddo
                v_I_T(1:ie,0,1:ke/px)=v_I(1:ie,0,1:ke/px)
                ru_u_T(1:ie,0,1:ke/px)=ru_u(1:ie,0,1:ke/px)
	else
		call mpi_recv(v_I_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(ru_u_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+100,MPI_COMM_WORLD,status,ierr)
	endif
	!! also pass over boundaries at j=jmax+1 :
	IF (rank.eq.px-1) THEN
	  do i=0,px-2
      		call mpi_send(ru_u(1:ie,je+1,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+200,MPI_COMM_WORLD,status,ierr)
	  enddo
                ru_u_T(1:ie,je*px+1,1:ke/px)=ru_u(1:ie,je+1,rank*ke/px+1:(rank+1)*ke/px)
	else
		call mpi_recv(ru_u_T(1:ie,je*px+1,1:ke/px),ie*ke/px,MPI_REAL8,px-1,rank+200,MPI_COMM_WORLD,status,ierr)
	endif

      	bby_0jepx=1.
	aay_0jepx=const_1_6 !3./10.
      	aay_0jepx(0)=0.
	aay_0jepx(je*px)=0. !5. !0.
	ccy_0jepx=const_1_6 !3./10.
	ccy_0jepx(0)=0. !5. !0.
	ccy_0jepx(je*px)=0.
      	bby_1jepx=1.
	aay_1jepx=const_1_22 !9./62.
      	aay_1jepx(1)=0.
	aay_1jepx(je*px)=-1. !-1.!*26./24. !0. !const_1_22 !0. !const_1_22 !0.! -1.

	ccy_1jepx=const_1_22 !9./62.
	ccy_1jepx(1)=-1. !-1.!*26./24. !0. !const_1_22 !0. !const_1_22 !0. !-1.
	ccy_1jepx(je*px)=0.
	
	do i=ib,ie
	  do k=1,ke/px
		!! interpolate_y ru_u_T -> I-loc		
		ddy_0jepx(0)= 0.5*(ru_u_T(i,0,k)+ru_u_T(i,1,k)) 
		ddy_0jepx(je*px)=0.5*(ru_u_T(i,je*px,k)+ru_u_T(i,je*px+1,k)) 
		do j=1,je*px-1
			ddy_0jepx(j)=const_2_3*(ru_u_T(i,j+1,k)+ru_u_T(i,j,k))
		enddo
		CALL solve_tridiag(ru_I_T(0:je*px),aay_0jepx,bby_0jepx,ccy_0jepx,ddy_0jepx,je*px+1)	!interpolate ru_u_T -> I-loc -> 0:je*px

		ruv_I_T(0:je*px)=ru_I_T(0:je*px)*v_I_T(i,0:je*px,k)

		ddy_1jepx(1)=(-ruv_I_T(0)+2.*ruv_I_T(1)-ruv_I_T(2))/(Ru(i)*dphi) 
		ddy_1jepx(je*px)=(ruv_I_T(je*px)-2.*ruv_I_T(je*px-1)+ruv_I_T(je*px-2))/(Ru(i)*dphi)
		do j=2,je*px-1 
			ddy_1jepx(j)=const_12_11*(ruv_I_T(j)-ruv_I_T(j-1))/(Ru(i)*dphi)
		enddo
               CALL solve_tridiag(dfy_T(1:je*px),aay_1jepx,bby_1jepx,ccy_1jepx,ddy_1jepx,je*px)   !dfy on Uvel(i,j,k) location -> 1:je*px
               putout_T(i,1:je*px,k)=-dfy_T(1:je*px)

	  enddo
	enddo

	call t2fp(putout_T(ib:ie,1:je*px,kb:ke/px),putouty(ib:ie,jb:je,kb:ke))
	putout(ib:ie,jb:je,kb:ke)=putout(ib:ie,jb:je,kb:ke)+putouty(ib:ie,jb:je,kb:ke)

      do k=kb,ke
        do j=jb,je
        jp=j+1
        jm=j-1
          do  i=ib,ie
          ip=i+1
          im=i-1
      putout(i,j,k) = putout(i,j,k) + 0.25 * (
     4 0.5*(rho(i,j,k)+rho(ip,j,k))*( Vvel(i,j,k) + Vvel(ip,j,k) + Vvel(i,jm,k) + Vvel(ip,jm,k) )**2
     4  / ( 4.0 * Ru(i) )
     +                         )
           enddo
        enddo
      enddo

      return
      end

      subroutine advecv_COM4(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
      implicit none
      include 'mpif.h'
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
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
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
      integer  rank,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1)
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm

      real aax_0ie(0:ie),bbx_0ie(0:ie),ccx_0ie(0:ie),ddx_0ie(0:ie)
      real aax_1ie(1:ie),bbx_1ie(1:ie),ccx_1ie(1:ie),ddx_1ie(1:ie)

!      real ruu(0:i1),u_r(0:i1),ru_u(0:i1,0:j1,0:k1)
      real u_I(0:ie,1:je,1:ke)

      real v_I_T(0:ie,0:je*px+1,1:ke/px)
      real rvu(0:ie),dfx(1:ie),u_I_T(0:ie,0:je*px,1:ke)
      real rv_I(0:ie)
	
	real rho_T(0:ie,0:je*px+1,1:ke/px),rho_v_T(0:ie,0:je*px,1:ke/px),rv_v(0:ie,1:je,1:ke)
!	real rho_T(1:ie,0:je*px+1,0:ke/px+1),rho_v_T(1:ie,0:je*px,0:ke/px+1),rv_v(1:ie,1:je,0:ke+1)

	real u_T(0:ie,0:je*px+1,1:ke/px),v_T(1:ie,0:je*px,1:ke/px),w_T(1:ie,0:je*px+1,1:ke/px)
      real rvv_T(1:je*px),dfy_T(1:ie,0:je*px,1:ke/px),putouty(1:ie,1:je,1:ke)
      real aay_0jepx(0:je*px),bby_0jepx(0:je*px),ccy_0jepx(0:je*px),ddy_0jepx(0:je*px)
      real aay_1jepx(1:je*px),bby_1jepx(1:je*px),ccy_1jepx(1:je*px),ddy_1jepx(1:je*px)

      real w_II(1:ie,1:je,0:ke),dfz(1:ke),rvw(0:ke),rv_II(0:ke),w_II_T(1:ie,0:je*px,1:ke/px)

      real aaz_0ke(0:ke),bbz_0ke(0:ke),ccz_0ke(0:ke),ddz_0ke(0:ke)
      real aaz_1ke(1:ke),bbz_1ke(1:ke),ccz_1ke(1:ke),ddz_1ke(1:ke)

      integer ileng,ierr,itag,status(MPI_STATUS_SIZE)
	real const_1_6,const_1_22,const_12_11,const_2_3
	putout=0.

	const_1_6=1./6.
	const_1_22=1./22.
	const_12_11=12./11.
	const_2_3=2./3.

      	bbx_0ie=1.
	aax_0ie=const_1_6 !3./10.
      	aax_0ie(0)=0.
	aax_0ie(ie)=0.
	ccx_0ie=const_1_6 !3./10.
	ccx_0ie(0)=0.
	ccx_0ie(ie)=0.
	bbx_1ie=1.

	aax_1ie=const_1_22 !9./62.
	aax_1ie(1)=0.
	aax_1ie(ie)=-1. !-1.!*26./24. !-1. !const_1_22 !-1. !0. !const_1_22 !0. !	1./23. !0. !2. !15.
	ccx_1ie=const_1_22 !9./62.
	ccx_1ie(1)=-1. !-1.!*26./24. !-1. !const_1_22 !-1. !0. !const_1_22  !1./23. !0. !2. !15.	
	ccx_1ie(ie)=0.	

	!! swap rho -> rho_T & u -> u_T
	call t2np_0ie(rho(0:ie,1:je,1:ke),rho_T(0:ie,1:je*px,1:ke/px))
	call t2np_0ie(Uvel(0:ie,1:je,1:ke),u_T(0:ie,1:je*px,1:ke/px))
	!! also pass over boundaries at j=0 :
	IF (rank.eq.0) THEN
	  do i=1,px-1
      		call mpi_send( rho(0:ie,0,i*ke/px+1:(i+1)*ke/px),i1*ke/px,MPI_REAL8,i,i+100,MPI_COMM_WORLD,status,ierr)
      		call mpi_send(Uvel(0:ie,0,i*ke/px+1:(i+1)*ke/px),i1*ke/px,MPI_REAL8,i,i+200,MPI_COMM_WORLD,status,ierr)
	  enddo
                u_T  (0:ie,0,1:ke/px)=Uvel(0:ie,0,1:ke/px)
                rho_T(0:ie,0,1:ke/px)= rho(0:ie,0,1:ke/px)
	else
		call mpi_recv(rho_T(0:ie,0,1:ke/px),i1*ke/px,MPI_REAL8,0,rank+100,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(  u_T(0:ie,0,1:ke/px),i1*ke/px,MPI_REAL8,0,rank+200,MPI_COMM_WORLD,status,ierr)
	endif
	!! also pass over boundaries at j=jmax+1 :
	IF (rank.eq.px-1) THEN
	  do i=0,px-2
      		call mpi_send( rho(0:ie,je+1,i*ke/px+1:(i+1)*ke/px),i1*ke/px,MPI_REAL8,i,i+300,MPI_COMM_WORLD,status,ierr)
      		call mpi_send(Uvel(0:ie,je+1,i*ke/px+1:(i+1)*ke/px),i1*ke/px,MPI_REAL8,i,i+400,MPI_COMM_WORLD,status,ierr)
	  enddo
                rho_T(0:ie,je*px+1,1:ke/px)= rho(0:ie,je+1,rank*ke/px+1:(rank+1)*ke/px)
                  u_T(0:ie,je*px+1,1:ke/px)=Uvel(0:ie,je+1,rank*ke/px+1:(rank+1)*ke/px)

	else
		call mpi_recv(rho_T(0:ie,je*px+1,1:ke/px),i1*ke/px,MPI_REAL8,px-1,rank+300,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(  u_T(0:ie,je*px+1,1:ke/px),i1*ke/px,MPI_REAL8,px-1,rank+400,MPI_COMM_WORLD,status,ierr)
	endif

     	bby_0jepx=1.
	aay_0jepx=const_1_6 !3./10.
      	aay_0jepx(0)=0.
	aay_0jepx(je*px)=0.
	ccy_0jepx=const_1_6 !3./10.
	ccy_0jepx(0)=0.
	ccy_0jepx(je*px)=0.

	do i=0,ie
	  do k=1,ke/px !k=0,ke/px+1
		!! interpolate_y rho_T -> v-loc		
		ddy_0jepx(0)=0.5*(rho_T(i,0,k)+rho_T(i,1,k)) 
		ddy_0jepx(je*px)=0.5*(rho_T(i,je*px,k)+rho_T(i,je*px+1,k)) 
 		do j=1,je*px-1
			ddy_0jepx(j)=const_2_3*(rho_T(i,j+1,k)+rho_T(i,j,k))
		enddo
		CALL solve_tridiag(rho_v_T(i,0:je*px,k),aay_0jepx,bby_0jepx,ccy_0jepx,ddy_0jepx,je*px+1)	! interpolate rho_T -> v-loc -> 0:je*px
	  enddo
	enddo

	do i=0,ie
	  do k=1,ke/px
		!! interpolate_y u_T -> I-loc		
		ddy_0jepx(0)=0.5*(u_T(i,0,k)+u_T(i,1,k)) 
		ddy_0jepx(je*px)=0.5*(u_T(i,je*px,k)+u_T(i,je*px+1,k)) 
 		do j=1,je*px-1
			ddy_0jepx(j)=const_2_3*(u_T(i,j+1,k)+u_T(i,j,k))
		enddo
		CALL solve_tridiag(u_I_T(i,0:je*px,k),aay_0jepx,bby_0jepx,ccy_0jepx,ddy_0jepx,je*px+1)	! interpolate_y u_T -> I-loc -> 0:je*px
	  enddo
	enddo

	! swap back
	call t2fp_0ie(rho_v_T(0:ie,1:je*px,1:ke/px),rv_v(0:ie,1:je,1:ke))
	call t2fp_0ie(u_I_T(0:ie,1:je*px,kb:ke/px),u_I(0:ie,1:je,1:ke))

	rv_v(0:ie,1:je,1:ke)=rv_v(0:ie,1:je,1:ke)*Vvel(0:ie,1:je,1:ke)

      do k=kb,ke !kb,ke
        do j=jb,je
		!! interpolate_x rv_v -> I-loc
		ddx_0ie(0)=0.5*(rv_v(0,j,k)+rv_v(1,j,k)) 
		ddx_0ie(ie)=rv_v(ie,j,k) 
		do i=1,ie-1
			ddx_0ie(i)=const_2_3*(rv_v(i+1,j,k)+rv_v(i,j,k))
		enddo
		CALL solve_tridiag(rv_I(0:ie),aax_0ie,bbx_0ie,ccx_0ie,ddx_0ie,i1)	! interpolate_x rv -> I-loc -> 0:ie
		!! derivative d(Rrhovu)Rdr
		rvu(0:ie)=rv_I(0:ie)*u_I(0:ie,j,k)*Ru(0:ie) 
                ddx_1ie(1)=(2.*rvu(1)-rvu(0)-rvu(2))/ ( Rp(1) * dr(1) ) 
                ddx_1ie(ie)=(-2.*rvu(ie-1)+rvu(ie)+rvu(ie-2))/ ( Rp(ie) * dr(ie) ) 
		do i=2,ie-1 
			ddx_1ie(i)=const_12_11*(rvu(i)-rvu(i-1))/ ( Rp(i) * dr(i) )
		enddo
		CALL solve_tridiag(dfx(1:ie),aax_1ie,bbx_1ie,ccx_1ie,ddx_1ie,ie)   ! df on Uvel(i,j,k) location	-> 1:ie 
		putout(1:ie,j,k)=putout(1:ie,j,k)-dfx(1:ie) ! derivative: d(rvu)/rdr

	enddo
      enddo

	call t2np(Vvel(1:ie,1:je,1:ke),v_T(1:ie,1:je*px,1:ke/px))
	!! also pass over boundaries at j=0 :
	IF (rank.eq.0) THEN
	  do i=1,px-1
      		call mpi_send(Vvel(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+100,MPI_COMM_WORLD,status,ierr)
		v_T(1:ie,0,1:ke/px)=Vvel(1:ie,0,1:ke/px)
	  enddo
	else
		call mpi_recv(v_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+100,MPI_COMM_WORLD,status,ierr)
	endif

      	bby_0jepx=1.
	aay_0jepx=const_1_22 !9./62.
      	aay_0jepx(0)=0.
	aay_0jepx(je*px)=3.!*3./8. !0. !1./23. !2. !0. !2. !1./15. !1./23.
	ccy_0jepx=const_1_22 !9./62.
	ccy_0jepx(0)=3.!*3./8. !0. !1./23. !2. !0. !2. !1./15. !1./23.
	ccy_0jepx(je*px)=0.

      	bby_1jepx=1.
	aay_1jepx=const_1_6 !3./10.
      	aay_1jepx(1)=0.
	aay_1jepx(je*px)=0. !1. !0. !1.
	ccy_1jepx=const_1_6 !3./10.
	ccy_1jepx(je*px)=0.
	ccy_1jepx(1)=0. !1. !0. !1. 


	do i=ib,ie
	  do k=1,ke/px
		!! interpolate_y v_T -> rho-loc		
		ddy_1jepx(1)=0.5*(v_T(i,0,k)+v_T(i,1,k)) !
		ddy_1jepx(je*px)=0.5*(v_T(i,je*px,k)+v_T(i,je*px-1,k)) !
 		do j=2,je*px-1
			ddy_1jepx(j)=const_2_3*(v_T(i,j,k)+v_T(i,j-1,k))
		enddo
		CALL solve_tridiag(rvv_T(1:je*px),aay_1jepx,bby_1jepx,ccy_1jepx,ddy_1jepx,je*px)	! interpolate v_T -> rho-loc -> 1:je*px
		rvv_T(1:je*px)=rvv_T(1:je*px)*rvv_T(1:je*px)*rho_T(i,1:je*px,k)
!		ddy_0jepx(0)=0.
!		ddy_0jepx(je*px)=0.
!		ddy_0jepx(0)=(-25.*rvv_T(1)+26.*rvv_T(2)-rvv_T(3))/(Rp(i)*dphi)
!		ddy_0jepx(je*px)=(+25.*rvv_T(je*px)-26.*rvv_T(je*px-1)+rvv_T(je*px-2))/(Rp(i)*dphi)
		ddy_0jepx(0)=8./3.*(-v_T(i,0,k)**2*0.5*(rho_T(i,0,k)+rho_T(i,1,k))+rvv_T(2))/(Rp(i)*dphi)
		ddy_0jepx(je*px)=8./3.*(v_T(i,je*px,k)**2*0.5*(rho_T(i,je*px,k)+rho_T(i,je*px+1,k))-rvv_T(je*px-1))/(Rp(i)*dphi)
		do j=1,je*px-1
			ddy_0jepx(j)=const_12_11*(rvv_T(j+1)-rvv_T(j))/(Rp(i)*dphi)
		enddo
		CALL solve_tridiag(dfy_T(i,0:je*px,k),aay_0jepx,bby_0jepx,ccy_0jepx,ddy_0jepx,je*px+1)   ! dfy on Uvel(i,j,k) location	-> 1:je*px
	  enddo
	enddo	
	call t2fp(dfy_T(1:ie,1:je*px,1:ke/px),putouty(ib:ie,jb:je,kb:ke))

	putout(ib:ie,jb:je,kb:ke)=putout(ib:ie,jb:je,kb:ke)-putouty(ib:ie,jb:je,kb:ke)

!!!!!!!!
	call t2np(Wvel(1:ie,1:je,1:ke),w_T(1:ie,1:je*px,1:ke/px))
	!! also pass over boundaries at j=0 :
	IF (rank.eq.0) THEN
	  do i=1,px-1
      		call mpi_send(Wvel(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+200,MPI_COMM_WORLD,status,ierr)
	  enddo
                w_T(1:ie,0,1:ke/px)=Wvel(1:ie,0,1:ke/px)
	else
		call mpi_recv(w_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+200,MPI_COMM_WORLD,status,ierr)
	endif
	!! also pass over boundaries at j=jmax+1 :
	IF (rank.eq.px-1) THEN
	  do i=0,px-2
      		call mpi_send(Wvel(1:ie,je+1,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+400,MPI_COMM_WORLD,status,ierr)
	  enddo
                w_T(1:ie,je*px+1,1:ke/px)=Wvel(1:ie,je+1,rank*ke/px+1:(rank+1)*ke/px)
	else
		call mpi_recv(w_T(1:ie,je*px+1,1:ke/px),ie*ke/px,MPI_REAL8,px-1,rank+400,MPI_COMM_WORLD,status,ierr)
	endif

     	bby_0jepx=1.
	aay_0jepx=const_1_6 !3./10.
      	aay_0jepx(0)=0.
	aay_0jepx(je*px)=0.
	ccy_0jepx=const_1_6 !3./10.
	ccy_0jepx(0)=0.
	ccy_0jepx(je*px)=0.

	do i=ib,ie
	  do k=1,ke/px
		!! interpolate_y w_T -> II-loc		
		ddy_0jepx(0)=0.5*(w_T(i,0,k)+w_T(i,1,k)) 
		ddy_0jepx(je*px)=0.5*(w_T(i,je*px,k)+w_T(i,je*px+1,k)) 
 		do j=1,je*px-1
			ddy_0jepx(j)=const_2_3*(w_T(i,j+1,k)+w_T(i,j,k))
		enddo
		CALL solve_tridiag(w_II_T(i,0:je*px,k),aay_0jepx,bby_0jepx,ccy_0jepx,ddy_0jepx,je*px+1)	! interpolate w_T -> II-loc -> 0:je*px
	  enddo
	enddo
 
	call t2fp(w_II_T(1:ie,1:je*px,1:ke/px),w_II(1:ie,1:je,1:ke))
	w_II(ib:ie,jb:je,0)=0.5*(Wvel(ib:ie,jb:je,0)+Wvel(ib:ie,jb+1:je+1,0))  !non compact interpolation, no problem -> W(i,j,0)=0.

      	bbz_0ke=1.
	aaz_0ke=const_1_6 !3./10.
      	aaz_0ke(0)=0.
	aaz_0ke(ke)=0.
	ccz_0ke=const_1_6 !3./10.
	ccz_0ke(0)=0.
	ccz_0ke(ke)=0.

      	bbz_1ke=1.
	aaz_1ke=const_1_22 !9./62.
      	aaz_1ke(1)=0.
	aaz_1ke(ke)=-1. !-1.!*26./24. !0. !const_1_22 !1./11. !! met aanpassing voor k=0 en k=ke ontploft hij, zowel van Lele, als van Hokpunna 2010
	aaz_1ke(2)=const_1_22!*21./24.
	aaz_1ke(ke-1)=const_1_22!*21./24.

	ccz_1ke=const_1_22 !9./62.
	ccz_1ke(1)=-1. !-1.!*26./24. !0. !const_1_22 !1./11.
	ccz_1ke(ke)=0.
	ccz_1ke(2)=const_1_22!*21./24.
	ccz_1ke(ke-1)=const_1_22!*21./24.

      do i=ib,ie 
        do j=jb,je
		!! interpolate_z rv_v -> II-loc
!		ddz_0ke(0)=0.125*(Vvel(i,j,0)+Vvel(i,j,1))*(rho(i,j,0)+rho(i,j,1)+rho(i,j+1,0)+rho(i,j+1,1)) 
!		ddz_0ke(ke)=0.125*(Vvel(i,j,k1)+Vvel(i,j,ke))*(rho(i,j,k1)+rho(i,j,ke)+rho(i,j+1,k1)+rho(i,j+1,ke)) 
		ddz_0ke(0)=rv_v(i,j,1) !0.5*(rv_v(i,j,1)+rv_v(i,j,0))
		ddz_0ke(ke)=rv_v(i,j,ke) !0.5*(rv_v(i,j,ke+1)+rv_v(i,j,ke))
		do k=1,ke-1
			ddz_0ke(k)=const_2_3*(rv_v(i,j,k+1)+rv_v(i,j,k))
		enddo
		CALL solve_tridiag(rv_II(0:ke),aaz_0ke,bbz_0ke,ccz_0ke,ddz_0ke,k1)	! interpolate_z rv_v -> II-loc -> 0:ke

		!! derivative d(rhovw)dz
		rvw(0:ke)=rv_II(0:ke)*w_II(i,j,0:ke)

		ddz_1ke(1)=(-rvw(0)+2.*rvw(1)-rvw(2))/ dz 
		ddz_1ke(ke)=(+rvw(ke)-2.*rvw(ke-1)+rvw(ke-2))/ dz 
		do k=2,ke-1
			ddz_1ke(k)=const_12_11*(rvw(k)-rvw(k-1))/ dz
		enddo
		CALL solve_tridiag(dfz(1:ke),aaz_1ke,bbz_1ke,ccz_1ke,ddz_1ke,ke)   ! df on Uvel(i,j,k) location	-> 1:ke 
		putout(i,j,1:ke)=putout(i,j,1:ke)-dfz(1:ke) ! derivative: d(rvw)/dz
	enddo
      enddo



       do k=kb,ke
         do j=jb,je
         jp=j+1
         jm=j-1
           do  i=ib,ie
           ip=i+1
           im=i-1
       putout(i,j,k) = putout(i,j,k)- 0.25 * (
     4	0.5*(rho(i,j,k)+rho(i,jp,k))*((Uvel(i,j,k)+Uvel(i,jp,k)+Uvel(im,j,k)+Uvel(im,jp,k))*Vvel(i,j,k))/Rp(i)	 
     4  )
            enddo
         enddo
       enddo
       return
      end

      subroutine advecw_COM4(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
      implicit none

	include 'mpif.h'
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
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
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
      integer  rank,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1)
      real rhoip,rhoim,rhojp,rhojm

	real rw_w(0:ie,0:j1,0:ke),u_I(0:ie,jb:je,0:ke),rwu_I(0:ie),rw_I(0:ie)
	real dfx(1:ie),dfy_T(ib:ie,jb:je*px,kb:ke/px),putouty(ib:ie,jb:je,kb:ke),dfz(0:ke)
	real rw_w_T(ib:ie,0:je*px+1,1:ke/px),v_II(ib:ie,0:je,0:ke),v_II_T(ib:ie,0:je*px,1:ke/px)
	real rw_II_T(0:je*px),rwv_II_T(0:je*px),w_r(1:ke),rww_r(1:ke),rww_r2(0:k1)
	real ww(0:i1,0:j1,0:k1)

      real aax_0ie(0:ie),bbx_0ie(0:ie),ccx_0ie(0:ie),ddx_0ie(0:ie)
      real aax_1ie(1:ie),bbx_1ie(1:ie),ccx_1ie(1:ie),ddx_1ie(1:ie)

      real aay_0jepx(0:je*px),bby_0jepx(0:je*px),ccy_0jepx(0:je*px),ddy_0jepx(0:je*px)
      real aay_1jepx(1:je*px),bby_1jepx(1:je*px),ccy_1jepx(1:je*px),ddy_1jepx(1:je*px)

      real aaz_0ke(0:ke),bbz_0ke(0:ke),ccz_0ke(0:ke),ddz_0ke(0:ke)
      real aaz_1ke(1:ke),bbz_1ke(1:ke),ccz_1ke(1:ke),ddz_1ke(1:ke)

      integer ileng,ierr,itag,status(MPI_STATUS_SIZE)
	real const_1_6,const_1_22,const_12_11,const_2_3
	putout=0.

	const_1_6=1./6.
	const_1_22=1./22.
	const_12_11=12./11.
	const_2_3=2./3.


      	bbz_0ke=1.
	aaz_0ke=const_1_6 !3./10.
      	aaz_0ke(0)=0.
	aaz_0ke(ke)=0.
	ccz_0ke=const_1_6 !3./10.
	ccz_0ke(0)=0.
	ccz_0ke(ke)=0.

      	bbx_0ie=1.
	aax_0ie=const_1_6 !3./10.
      	aax_0ie(0)=0.
	aax_0ie(ie)=0.
	ccx_0ie=const_1_6 !3./10.
	ccx_0ie(0)=0.
	ccx_0ie(ie)=0.

	bbx_1ie=1.
	aax_1ie=const_1_22 !9./62.
	aax_1ie(1)=0.
	aax_1ie(ie)=-1. !-1.!*26./24. !-1. !const_1_22 !-1. !0. !const_1_22 !0. !	1./23. !0. !2. !15.
	ccx_1ie=const_1_22 !9./62.
	ccx_1ie(1)=-1. !0. !-1.!*26./24. !-1. !const_1_22 !-1. !0. !const_1_22  !1./23. !0. !2. !15.	
	ccx_1ie(ie)=0.	

      do i=0,ie
        do j=0,j1
		!! interpolate_z rho -> w-loc
		ddz_0ke(0)=0.5*(rho(i,j,0)+rho(i,j,1)) 
		ddz_0ke(ke)=0.5*(rho(i,j,ke)+rho(i,j,k1))
		do k=1,ke-1
			ddz_0ke(k)=const_2_3*(rho(i,j,k+1)+rho(i,j,k))
		enddo
		CALL solve_tridiag(rw_w(i,j,0:ke),aaz_0ke,bbz_0ke,ccz_0ke,ddz_0ke,k1)	! interpolate rho -> w-loc -> 0:ke
		rw_w(i,j,0:ke)=rw_w(i,j,0:ke)*Wvel(i,j,0:ke) 
	enddo
      enddo
      do i=0,ie
        do j=jb,je
		!! interpolate_z u -> I-loc
		ddz_0ke(0)=0.5*(Uvel(i,j,0)+Uvel(i,j,1)) 
		ddz_0ke(ke)=0.5*(Uvel(i,j,ke)+Uvel(i,j,k1))
		do k=1,ke-1
			ddz_0ke(k)=const_2_3*(Uvel(i,j,k+1)+Uvel(i,j,k))
		enddo
		CALL solve_tridiag(u_I(i,j,0:ke),aaz_0ke,bbz_0ke,ccz_0ke,ddz_0ke,k1)	! interpolate_z u -> I-loc -> 0:ke
	enddo
      enddo

      do k=kb,ke !kb,ke
        do j=jb,je
		!! interpolate_x rw_w -> I-loc
		ddx_0ie(0)=0.5*(rw_w(0,j,k)+rw_w(1,j,k))
		ddx_0ie(ie)=rw_w(ie,j,k) !uitstroom is upstream interpolatie 
		do i=1,ie-1
			ddx_0ie(i)=const_2_3*(rw_w(i+1,j,k)+rw_w(i,j,k))
		enddo
		CALL solve_tridiag(rw_I(0:ie),aax_0ie,bbx_0ie,ccx_0ie,ddx_0ie,i1)	! interpolate_x rw -> I-loc -> 0:ie

		!! derivative d(Rrhowu)Rdr
		rwu_I(0:ie)=rw_I(0:ie)*u_I(0:ie,j,k)*Ru(0:ie)
                ddx_1ie(1)=(-rwu_I(0)+2.*rwu_I(1)-rwu_I(2))/ ( Rp(1) * dr(1) )
                ddx_1ie(ie)=(rwu_I(ie)-2.*rwu_I(ie-1)+rwu_I(ie-2))/ ( Rp(ie) * dr(ie) )
		do i=2,ie-1 !
			ddx_1ie(i)=const_12_11*(rwu_I(i)-rwu_I(i-1))/ ( Rp(i) * dr(i) )
		enddo
		CALL solve_tridiag(dfx(1:ie),aax_1ie,bbx_1ie,ccx_1ie,ddx_1ie,ie)   ! dfx on Wvel(i,j,k) location	-> 1:ie 
		putout(1:ie,j,k)=putout(1:ie,j,k)-dfx(1:ie) ! derivative: d(rwu)/rdr

	enddo
      enddo

      do i=ib,ie
        do j=0,je
		!! interpolate_z v -> II-loc
		ddz_0ke(0)=0.5*(Vvel(i,j,0)+Vvel(i,j,1)) 
		ddz_0ke(ke)=0.5*(Vvel(i,j,ke)+Vvel(i,j,k1))
		do k=1,ke-1
			ddz_0ke(k)=const_2_3*(Vvel(i,j,k+1)+Vvel(i,j,k))
		enddo
		CALL solve_tridiag(v_II(i,j,0:ke),aaz_0ke,bbz_0ke,ccz_0ke,ddz_0ke,k1)	! interpolate_z v -> II-loc -> 0:ke
	enddo
      enddo

	call t2np(v_II(1:ie,1:je,1:ke),v_II_T(1:ie,1:je*px,1:ke/px))
	!! also pass over boundaries at j=0 :
	IF (rank.eq.0) THEN
	  do i=1,px-1
      		call mpi_send(v_II(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+201,MPI_COMM_WORLD,status,ierr)
	  enddo
                v_II_T(1:ie,0,1:ke/px)=v_II(1:ie,0,1:ke/px)
	else
		call mpi_recv(v_II_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+201,MPI_COMM_WORLD,status,ierr)
	endif
	call t2np(rw_w(1:ie,1:je,1:ke),rw_w_T(1:ie,1:je*px,1:ke/px))
	!! also pass over boundaries at j=0 :
	IF (rank.eq.0) THEN
	  do i=1,px-1
      		call mpi_send(rw_w(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+202,MPI_COMM_WORLD,status,ierr)
	  enddo
                rw_w_T(1:ie,0,1:ke/px)=rw_w(1:ie,0,1:ke/px)
	else
		call mpi_recv(rw_w_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+202,MPI_COMM_WORLD,status,ierr)
	endif

	!! also pass over boundaries at j=jmax+1 :
	IF (rank.eq.px-1) THEN
	  do i=0,px-2
      		call mpi_send(rw_w(1:ie,je+1,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+401,MPI_COMM_WORLD,status,ierr)
	  enddo
                rw_w_T(1:ie,je*px+1,1:ke/px)=rw_w(1:ie,je+1,rank*ke/px+1:(rank+1)*ke/px)
	else
		call mpi_recv(rw_w_T(1:ie,je*px+1,1:ke/px),ie*ke/px,MPI_REAL8,px-1,rank+401,MPI_COMM_WORLD,status,ierr)
	endif

      	bby_0jepx=1.
	aay_0jepx=const_1_6 !3./10.
      	aay_0jepx(0)=0.
	aay_0jepx(je*px)=0.
	ccy_0jepx=const_1_6 !3./10.
	ccy_0jepx(0)=0.
	ccy_0jepx(je*px)=0.
      	bby_1jepx=1.
	aay_1jepx=const_1_22 !9./62.
      	aay_1jepx(1)=0.
	aay_1jepx(je*px)=-1. !-1.!*26./24. !0. !const_1_22 !0. !const_1_22 !0.! -1.
	ccy_1jepx=const_1_22 !9./62.
	ccy_1jepx(1)=-1. !-1.!*26./24. !0. !const_1_22 !0. !const_1_22 !0. !-1.
	ccy_1jepx(je*px)=0.

	do i=ib,ie
	  do k=kb,ke/px
		!! interpolate_y rw_w_T -> II-loc		
		ddy_0jepx(0)=0.5*(rw_w_T(i,0,k)+rw_w_T(i,1,k)) 
		ddy_0jepx(je*px)=0.5*(rw_w_T(i,je*px,k)+rw_w_T(i,je*px+1,k)) 
 		do j=1,je*px-1
			ddy_0jepx(j)=const_2_3*(rw_w_T(i,j+1,k)+rw_w_T(i,j,k))
		enddo
		CALL solve_tridiag(rw_II_T(0:je*px),aay_0jepx,bby_0jepx,ccy_0jepx,ddy_0jepx,je*px+1)	! interpolate rw_w_T -> II-loc -> 0:je*px

		rwv_II_T(0:je*px)=rw_II_T(0:je*px)*v_II_T(i,0:je*px,k)
		!! diff_y drwv/(rdphi) 
		ddy_1jepx(1)=(-rwv_II_T(0)+2.*rwv_II_T(1)-rwv_II_T(2))/(Rp(i)*dphi)
		ddy_1jepx(je*px)=(rwv_II_T(je*px)-2.*rwv_II_T(je*px-1)+rwv_II_T(je*px-2))/(Rp(i)*dphi)
 		do j=2,je*px-1
			ddy_1jepx(j)=const_12_11*(rwv_II_T(j)-rwv_II_T(j-1))/(Rp(i)*dphi)
		enddo
		CALL solve_tridiag(dfy_T(i,1:je*px,k),aay_1jepx,bby_1jepx,ccy_1jepx,ddy_1jepx,je*px)	! diff_y  -> W-loc -> 1:je*px
	  enddo
	enddo
	call t2fp(dfy_T(1:ie,1:je*px,1:ke/px),putouty(ib:ie,jb:je,kb:ke))

	putout(ib:ie,jb:je,kb:ke)=putout(ib:ie,jb:je,kb:ke)-putouty(ib:ie,jb:je,kb:ke)

      	bbz_0ke=1.
	aaz_0ke=const_1_22 !9./62.
      	aaz_0ke(0)=0.
	aaz_0ke(ke)=3.!*3./8. !15. !23. !1./23. !1./15. !1./23.
	ccz_0ke=const_1_22 !9./62.
	ccz_0ke(0)=3.!*3./8. !15. !23. !0. !1./23.
	ccz_0ke(ke)=0.

      	bbz_1ke=1.
	aaz_1ke=const_1_6 !3./10.
      	aaz_1ke(1)=0.
	aaz_1ke(ke)=0. !1. !0. !1. !0. !1.
	ccz_1ke=const_1_6 !3./10.
	ccz_1ke(ke)=0.
	ccz_1ke(1)=0. !0. !1. !0. !1. !0. !1. 

	ww=Wvel
	do i=ib,ie
	  do j=jb,je
		!! interpolate_z Wvel -> rho-loc		
		ddz_1ke(1)=0.5*ww(i,j,0)+.5*ww(i,j,1)
		ddz_1ke(ke)=0.5*ww(i,j,ke)+.5*ww(i,j,ke-1)
 		do k=2,ke-1
			ddz_1ke(k)=const_2_3*(ww(i,j,k)+ww(i,j,k-1))
		enddo
		CALL solve_tridiag(w_r(1:ke),aaz_1ke,bbz_1ke,ccz_1ke,ddz_1ke,ke)	! interpolate_z Wvel -> rho-loc -> 1:ke
		rww_r(1:ke)=w_r(1:ke)*w_r(1:ke)*rho(i,j,1:ke)

		ddz_0ke(ke)=8./3.*(ww(i,j,ke)*ww(i,j,ke)*0.5*(rho(i,j,k1)+rho(i,j,ke))-rww_r(ke-1))/(dz)
		ddz_0ke(0)=8./3.*(-ww(i,j,0)*ww(i,j,0)*0.5*(rho(i,j,1)+rho(i,j,0))+rww_r(2))/(dz)
		do k=1,ke-1
			ddz_0ke(k)=const_12_11*(rww_r(k+1)-rww_r(k))/dz
		enddo
		CALL solve_tridiag(dfz(0:ke),aaz_0ke,bbz_0ke,ccz_0ke,ddz_0ke,ke+1)   ! dfz on Wvel(i,j,k) location	-> 0:ke
		putout(i,j,1:ke)=putout(i,j,1:ke)-dfz(1:ke)
	  enddo
	enddo	
      return
      end
  

      subroutine solve_tridiag(x,a,b,c,v,n)
      implicit none
!	 a - sub-diagonal (means it is the diagonal below the main diagonal) 
!	 [a(1) is on line 1 of the matrix system and does not exist [give dummy value], a(end) is on end line of matrix system]
!	 b - the main diagonal 
!	 c - sup-diagonal (means it is the diagonal above the main diagonal) 
!	[c(1) is on line 1 of the matrix system, c(end) is on end line of matrix system and does not exist [give dummy value]]
!	 v - right part
!	 x - the answer
!	 n - number of equations
 
        integer,intent(in) :: n
        real,dimension(n),intent(in) :: a,b,c,v
        real,dimension(n),intent(out) :: x
        real,dimension(n) :: bp,vp
        real :: m
        integer i
 
! Make copies of the b and v variables so that they are unaltered by this sub
        bp(1) = b(1)
        vp(1) = v(1)
 
        !The first pass (setting coefficients):
	do i = 2,n
         m = a(i)/bp(i-1)
         bp(i) = b(i) - m*c(i-1)
         vp(i) = v(i) - m*vp(i-1)
        end do 
 
         x(n) = vp(n)/bp(n)
        !The second pass (back-substition)
	do i = n-1, 1, -1
          x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
        end do 
 
      end subroutine solve_tridiag

   
      subroutine solve_tridiag_periodic(x,a,b,c,v,n)
      implicit none
!	 a - sub-diagonal (means it is the diagonal below the main diagonal)
!	 b - the main diagonal
!	 c - sup-diagonal (means it is the diagonal above the main diagonal)
!	 v - right part
!	 x - the answer
!	 n - number of equations
 
        integer,intent(in) :: n
        real,dimension(n),intent(in) :: a,b,c,v
        real,dimension(n),intent(out) :: x
        real,dimension(n) :: bp,vp,bper,u,y,q,up
        real :: m,corr
        integer i
 
! Make new slightly adjusted diagonal to solve periodic system:
	bper = b
	bper(1) = b(1) + c(1)
	bper(n) = b(n) + a(n)
! Make extra rhs:
	u = 0.
	u(1) = 1.
	u(n) = -1.
	

!! Solve first system: Aper*y=v
!! Make copies of the b and v variables so that they are unaltered by this sub
        bp(1) = bper(1)
        vp(1) = v(1)
 
        !The first pass (setting coefficients):
	do i = 2,n
         m = a(i)/bp(i-1)
         bp(i) = bper(i) - m*c(i-1)
         vp(i) = v(i) - m*vp(i-1)
        end do 
 
         y(n) = vp(n)/bp(n)
        !The second pass (back-substition)
	do i = n-1, 1, -1
          y(i) = (vp(i) - c(i)*y(i+1))/bp(i)
        end do 
!! Solve second system: Aper*q=u
!! Make copies of the b and v variables so that they are unaltered by this sub
        bp(1) = bper(1)
        up(1) = u(1)
 
        !The first pass (setting coefficients):
	do i = 2,n
         m = a(i)/bp(i-1)
         bp(i) = bper(i) - m*c(i-1)
         up(i) = u(i) - m*up(i-1)
        end do 
 
         q(n) = up(n)/bp(n)
        !The second pass (back-substition)
	do i = n-1, 1, -1
          q(i) = (up(i) - c(i)*q(i+1))/bp(i)
        end do  

!!	CALL solve_tridiag(y,a,bper,c,v,n)
!	CALL solve_tridiag(q,a,bper,c,u,n)

	!! correction: x=y + (vvT*y/(1-vvT*q))*q
	!! vvt=[c(n),0,0,...,0,-a(1)]
	corr=(c(n)*y(1)-a(1)*y(n))/(1.-c(n)*q(1)+a(1)*q(n))
	x=y+corr*q

      end subroutine solve_tridiag_periodic



