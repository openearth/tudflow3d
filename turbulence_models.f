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


	subroutine compact_filter(Uvel,Vvel,Wvel,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px,a,
     &   tmax_inPpunt,i_inPpunt,j_inPpunt,tmax_inUpunt,i_inUpunt,j_inUpunt,tmax_inVpunt,i_inVpunt,j_inVpunt,kjet)
	implicit none

	include 'mpif.h'

	real Uvel(0:i1,0:j1,0:k1),Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1)
      integer  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px
	real dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1),a
      	integer i,j,k
        integer ileng,ierr,itag,status(MPI_STATUS_SIZE)

	real aax_1ie(1:ie),bbx_1ie(1:ie),ccx_1ie(1:ie),ddx_1ie(1:ie)
	real aaz_1ke(1:ke),bbz_1ke(1:ke),ccz_1ke(1:ke),ddz_1ke(1:ke)
	real aay_1jepx(1:je*px),bby_1jepx(1:je*px),ccy_1jepx(1:je*px),ddy_1jepx(1:je*px)

	real var(0:i1,1:je,0:k1),filt_var(1:ie,1:je,1:ke)
	real Uvel_T(1:ie,0:je*px+1,1:ke/px),Vvel_T(1:ie,1:je*px,1:ke/px),Wvel_T(1:ie,0:je*px+1,1:ke/px)
	real var_T(1:ie,0:je*px+1,1:ke/px),filt_var_T(1:ie,1:je*px,1:ke/px)

	real Uvel2(0:i1,0:j1,0:k1),Vvel2(0:i1,0:j1,0:k1),Wvel2(0:i1,0:j1,0:k1)
	integer tmax_inPpunt,i_inPpunt(i1*j1),j_inPpunt(i1*j1)
	integer tmax_inVpunt,i_inVpunt(i1*j1),j_inVpunt(i1*j1)
	integer tmax_inUpunt,i_inUpunt(i1*j1),j_inUpunt(i1*j1)
	integer kjet,t

       	bbx_1ie=1.
	aax_1ie=a
      	aax_1ie(1)=0.
	aax_1ie(ie)=0.
	ccx_1ie=a
	ccx_1ie(1)=0.
	ccx_1ie(ie)=0.
       	bbz_1ke=1.
	aaz_1ke=a
      	aaz_1ke(1)=0.
      	aaz_1ke(ke)=0.
	ccz_1ke=a
	ccz_1ke(1)=0.
	ccz_1ke(ke)=0.

       	bby_1jepx=1.
	aay_1jepx=a
      	aay_1jepx(1)=0.
      	aay_1jepx(je*px)=0.
	ccy_1jepx=a
	ccy_1jepx(1)=0.
	ccy_1jepx(je*px)=0.

       !! Set boundary conditions jet in:
	Uvel2=Uvel
	Vvel2=Vvel
	Wvel2=Wvel

	var(0:i1,1:je,1:ke)=Uvel(0:i1,1:je,1:ke)
	do j=jb,je
	  do k=kb,ke
!	    ddx_1ie(1)=(0.0625*a+0.9375)*var(1,j,k)+(0.25+0.75*a)*var(2,j,k)
!     &                 -0.375*(1-a)*var(3,j,k)+0.25*(1-a)*var(4,j,k)-0.0625*(1-a)*var(5,j,k)
!	    ddx_1ie(ie)=(0.0625*a+0.9375)*var(ie,j,k)+(0.25+0.75*a)*var(ie-1,j,k)
!     &                 -0.375*(1-a)*var(ie-2,j,k)+0.25*(1-a)*var(ie-3,j,k)-0.0625*(1-a)*var(ie-4,j,k)
	    ddx_1ie(1)=var(1,j,k)
	    ddx_1ie(ie)=var(ie,j,k)
!	    ddx_1ie(2)=(0.0625+0.875*a)*var(1,j,k)+(0.75+0.5*a)*var(2,j,k)
!     &                 +(0.375+0.25*a)*var(3,j,k)-0.5*(0.5-a)*var(4,j,k)+0.125*(0.5-a)*var(5,j,k)
!	    ddx_1ie(ie-1)=(0.0625+0.875*a)*var(ie,j,k)+(0.75+0.5*a)*var(ie-1,j,k)
!     &                 +(0.375+0.25*a)*var(ie-2,j,k)-0.5*(0.5-a)*var(ie-3,j,k)+0.125*(0.5-a)*var(ie-4,j,k)
	    do i=2,ie-1
	      ddx_1ie(i)= ( 0.625+0.75*a )*var(i,j,k)
     &                   +( 0.25 +0.5 *a )*(var(i+1,j,k)+var(i-1,j,k))
     &                   +(-0.125+0.25*a )*(var(i+2,j,k)+var(i-2,j,k)) 
   	    enddo
	    CALL solve_tridiag(filt_var(1:ie,j,k),aax_1ie,bbx_1ie,ccx_1ie,ddx_1ie,ie)
	  enddo
	enddo
	Uvel(1:ie,1:je,1:ke)=filt_var(1:ie,1:je,1:ke)
	var(0:i1,1:je,1:ke)=Vvel(0:i1,1:je,1:ke)
	do j=jb,je
	  do k=kb,ke
!	    ddx_1ie(1)=(0.0625*a+0.9375)*var(1,j,k)+(0.25+0.75*a)*var(2,j,k)
!     &                 -0.375*(1-a)*var(3,j,k)+0.25*(1-a)*var(4,j,k)-0.0625*(1-a)*var(5,j,k)
!	    ddx_1ie(ie)=(0.0625*a+0.9375)*var(ie,j,k)+(0.25+0.75*a)*var(ie-1,j,k)
!     &                 -0.375*(1-a)*var(ie-2,j,k)+0.25*(1-a)*var(ie-3,j,k)-0.0625*(1-a)*var(ie-4,j,k)
	    ddx_1ie(1)=var(1,j,k)
	    ddx_1ie(ie)=var(ie,j,k)
!	    ddx_1ie(2)=(0.0625+0.875*a)*var(1,j,k)+(0.75+0.5*a)*var(2,j,k)
!     &                 +(0.375+0.25*a)*var(3,j,k)-0.5*(0.5-a)*var(4,j,k)+0.125*(0.5-a)*var(5,j,k)
!	    ddx_1ie(ie-1)=(0.0625+0.875*a)*var(ie,j,k)+(0.75+0.5*a)*var(ie-1,j,k)
!     &                 +(0.375+0.25*a)*var(ie-2,j,k)-0.5*(0.5-a)*var(ie-3,j,k)+0.125*(0.5-a)*var(ie-4,j,k)
	    do i=2,ie-1
	      ddx_1ie(i)= ( 0.625+0.75*a )*var(i,j,k)
     &                   +( 0.25 +0.5 *a )*(var(i+1,j,k)+var(i-1,j,k))
     &                   +(-0.125+0.25*a )*(var(i+2,j,k)+var(i-2,j,k)) 
   	    enddo
	    CALL solve_tridiag(filt_var(1:ie,j,k),aax_1ie,bbx_1ie,ccx_1ie,ddx_1ie,ie)
	  enddo
	enddo
	Vvel(1:ie,1:je,1:ke)=filt_var(1:ie,1:je,1:ke)          
	var(0:i1,1:je,1:ke)=Wvel(0:i1,1:je,1:ke)
	do j=jb,je
	  do k=kb,ke
!	    ddx_1ie(1)=(0.0625*a+0.9375)*var(1,j,k)+(0.25+0.75*a)*var(2,j,k)
!     &                 -0.375*(1-a)*var(3,j,k)+0.25*(1-a)*var(4,j,k)-0.0625*(1-a)*var(5,j,k)
!	    ddx_1ie(ie)=(0.0625*a+0.9375)*var(ie,j,k)+(0.25+0.75*a)*var(ie-1,j,k)
!     &                 -0.375*(1-a)*var(ie-2,j,k)+0.25*(1-a)*var(ie-3,j,k)-0.0625*(1-a)*var(ie-4,j,k)
	    ddx_1ie(1)=var(1,j,k)
	    ddx_1ie(ie)=var(ie,j,k)
!	    ddx_1ie(2)=(0.0625+0.875*a)*var(1,j,k)+(0.75+0.5*a)*var(2,j,k)
!     &                 +(0.375+0.25*a)*var(3,j,k)-0.5*(0.5-a)*var(4,j,k)+0.125*(0.5-a)*var(5,j,k)
!	    ddx_1ie(ie-1)=(0.0625+0.875*a)*var(ie,j,k)+(0.75+0.5*a)*var(ie-1,j,k)
!     &                 +(0.375+0.25*a)*var(ie-2,j,k)-0.5*(0.5-a)*var(ie-3,j,k)+0.125*(0.5-a)*var(ie-4,j,k)
	    do i=2,ie-1
	      ddx_1ie(i)= ( 0.625+0.75*a )*var(i,j,k)
     &                   +( 0.25 +0.5 *a )*(var(i+1,j,k)+var(i-1,j,k))
     &                   +(-0.125+0.25*a )*(var(i+2,j,k)+var(i-2,j,k)) 
   	    enddo
	    CALL solve_tridiag(filt_var(1:ie,j,k),aax_1ie,bbx_1ie,ccx_1ie,ddx_1ie,ie)
	  enddo
	enddo
	Wvel(1:ie,1:je,1:ke)=filt_var(1:ie,1:je,1:ke)

	var(1:ie,1:je,0:k1)=Uvel(1:ie,1:je,0:k1)
	do j=jb,je
	  do i=ib,ie
!	    ddz_1ke(1)=(0.0625*a+0.9375)*var(i,j,1)+(0.25+0.75*a)*var(i,j,2)
!     &                 -0.375*(1-a)*var(i,j,3)+0.25*(1-a)*var(i,j,4)-0.0625*(1-a)*var(i,j,5)
!	    ddz_1ke(ke)=(0.0625*a+0.9375)*var(i,j,ke)+(0.25+0.75*a)*var(i,j,ke-1)
!     &                 -0.375*(1-a)*var(i,j,ke-2)+0.25*(1-a)*var(i,j,ke-3)-0.0625*(1-a)*var(i,j,ke-4)
	    ddz_1ke(1)=var(i,j,1)
	    ddz_1ke(ke)=var(i,j,ke)
	    do k=2,ke-1
	      ddz_1ke(k)= ( 0.625+0.75*a )*var(i,j,k)
     &                   +( 0.25 +0.5 *a )*(var(i,j,k+1)+var(i,j,k-1))
     &                   +(-0.125+0.25*a )*(var(i,j,k+2)+var(i,j,k-2)) 
   	    enddo
	    CALL solve_tridiag(filt_var(i,j,1:ke),aaz_1ke,bbz_1ke,ccz_1ke,ddz_1ke,ke)
	  enddo
	enddo
	Uvel(1:ie,1:je,1:ke)=filt_var(1:ie,1:je,1:ke)
	var(1:ie,1:je,0:k1)=Vvel(1:ie,1:je,0:k1)
	do j=jb,je
	  do i=ib,ie
!	    ddz_1ke(1)=(0.0625*a+0.9375)*var(i,j,1)+(0.25+0.75*a)*var(i,j,2)
!     &                 -0.375*(1-a)*var(i,j,3)+0.25*(1-a)*var(i,j,4)-0.0625*(1-a)*var(i,j,5)
!	    ddz_1ke(ke)=(0.0625*a+0.9375)*var(i,j,ke)+(0.25+0.75*a)*var(i,j,ke-1)
!     &                 -0.375*(1-a)*var(i,j,ke-2)+0.25*(1-a)*var(i,j,ke-3)-0.0625*(1-a)*var(i,j,ke-4)
	    ddz_1ke(1)=var(i,j,1)
	    ddz_1ke(ke)=var(i,j,ke)
	    do k=2,ke-1
	      ddz_1ke(k)= ( 0.625+0.75*a )*var(i,j,k)
     &                   +( 0.25 +0.5 *a )*(var(i,j,k+1)+var(i,j,k-1))
     &                   +(-0.125+0.25*a )*(var(i,j,k+2)+var(i,j,k-2)) 
   	    enddo
	    CALL solve_tridiag(filt_var(i,j,1:ke),aaz_1ke,bbz_1ke,ccz_1ke,ddz_1ke,ke)
	  enddo
	enddo
	Vvel(1:ie,1:je,1:ke)=filt_var(1:ie,1:je,1:ke)
	var(1:ie,1:je,1:ke)=Wvel(1:ie,1:je,1:ke)
	do j=jb,je
	  do i=ib,ie
!	    ddz_1ke(1)=(0.0625*a+0.9375)*var(i,j,1)+(0.25+0.75*a)*var(i,j,2)
!     &                 -0.375*(1-a)*var(i,j,3)+0.25*(1-a)*var(i,j,4)-0.0625*(1-a)*var(i,j,5)
!	    ddz_1ke(ke)=(0.0625*a+0.9375)*var(i,j,ke)+(0.25+0.75*a)*var(i,j,ke-1)
!     &                 -0.375*(1-a)*var(i,j,ke-2)+0.25*(1-a)*var(i,j,ke-3)-0.0625*(1-a)*var(i,j,ke-4)
	    ddz_1ke(1)=var(i,j,1)
	    ddz_1ke(ke)=var(i,j,ke)
	    ddz_1ke(2)=(0.0625+0.875*a)*var(i,j,1)+(0.75+0.5*a)*var(i,j,2)
     &                 +(0.375+0.25*a)*var(i,j,3)-0.5*(0.5-a)*var(i,j,4)+0.125*(0.5-a)*var(i,j,5)
	    ddz_1ke(ke-1)=(0.0625+0.875*a)*var(i,j,ke)+(0.75+0.5*a)*var(i,j,ke-1)
     &                 +(0.375+0.25*a)*var(i,j,ke-2)-0.5*(0.5-a)*var(i,j,ke-3)+0.125*(0.5-a)*var(i,j,ke-4)
	    do k=3,ke-2
	      ddz_1ke(k)= ( 0.625+0.75*a )*var(i,j,k)
     &                   +( 0.25 +0.5 *a )*(var(i,j,k+1)+var(i,j,k-1))
     &                   +(-0.125+0.25*a )*(var(i,j,k+2)+var(i,j,k-2)) 
   	    enddo
	    CALL solve_tridiag(filt_var(i,j,1:ke),aaz_1ke,bbz_1ke,ccz_1ke,ddz_1ke,ke)
	  enddo
	enddo
	Wvel(1:ie,1:je,1:ke)=filt_var(1:ie,1:je,1:ke)

	!! first swap u,v,w -> u_T,v_T,w_T 
	call t2np(Uvel(1:ie,1:je,1:ke),Uvel_T(1:ie,1:je*px,1:ke/px))
	call t2np(Vvel(1:ie,1:je,1:ke),Vvel_T(1:ie,1:je*px,1:ke/px))
	call t2np(Wvel(1:ie,1:je,1:ke),Wvel_T(1:ie,1:je*px,1:ke/px))
	!! also pass over boundaries at j=0 :
	IF (rank.eq.0) THEN
	  do i=1,px-1
      		call mpi_send(Uvel(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+200,MPI_COMM_WORLD,status,ierr)
      		call mpi_send(Wvel(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+201,MPI_COMM_WORLD,status,ierr)
	  enddo
                Uvel_T(1:ie,0,1:ke/px)=Uvel(1:ie,0,1:ke/px)
                Wvel_T(1:ie,0,1:ke/px)=Wvel(1:ie,0,1:ke/px)
	else
		call mpi_recv(Uvel_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+200,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(Wvel_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+201,MPI_COMM_WORLD,status,ierr)
	endif
	!! also pass over boundaries at j=jmax+1 :
	IF (rank.eq.px-1) THEN
	  do i=0,px-2
      		call mpi_send(Uvel(1:ie,je+1,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+400,MPI_COMM_WORLD,status,ierr)
      		call mpi_send(Wvel(1:ie,je+1,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+401,MPI_COMM_WORLD,status,ierr)
	  enddo
                Uvel_T(1:ie,je*px+1,1:ke/px)=Uvel(1:ie,je+1,rank*ke/px+1:(rank+1)*ke/px)
                Wvel_T(1:ie,je*px+1,1:ke/px)=Wvel(1:ie,je+1,rank*ke/px+1:(rank+1)*ke/px)
	else
		call mpi_recv(Uvel_T(1:ie,je*px+1,1:ke/px),ie*ke/px,MPI_REAL8,px-1,rank+400,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(Wvel_T(1:ie,je*px+1,1:ke/px),ie*ke/px,MPI_REAL8,px-1,rank+401,MPI_COMM_WORLD,status,ierr)
	endif



	var_T(1:ie,0:je*px+1,1:ke/px)=Uvel_T(1:ie,0:je*px+1,1:ke/px)
	do i=ib,ie
	  do k=kb,ke/px
!	    ddy_1jepx(1)=(0.0625*a+0.9375)*var_T(i,1,k)+(0.25+0.75*a)*var_T(i,2,k)
!     &                 -0.375*(1-a)*var_T(i,3,k)+0.25*(1-a)*var_T(i,4,k)-0.0625*(1-a)*var_T(i,5,k)
!	    ddy_1jepx(je*px)=(0.0625*a+0.9375)*var_T(i,je*px,k)+(0.25+0.75*a)*var_T(i,je*px-1,k)
!     &                 -0.375*(1-a)*var_T(i,je*px-2,k)+0.25*(1-a)*var_T(i,je*px-3,k)-0.0625*(1-a)*var_T(i,je*px-4,k)
	    ddy_1jepx(1)=var_T(i,1,k)
	    ddy_1jepx(je*px)=var_T(i,je*px,k)
!	    ddy_1jepx(2)=(0.0625+0.875*a)*var_T(i,1,k)+(0.75+0.5*a)*var_T(i,2,k)
!     &                 +(0.375+0.25*a)*var_T(i,3,k)-0.5*(0.5-a)*var_T(i,4,k)+0.125*(0.5-a)*var_T(i,5,k)
!	    ddy_1jepx(je*px-1)=(0.0625+0.875*a)*var_T(i,je*px,k)+(0.75+0.5*a)*var_T(i,je*px-1,k)
!     &                 +(0.375+0.25*a)*var_T(i,je*px-2,k)-0.5*(0.5-a)*var_T(i,je*px-3,k)+0.125*(0.5-a)*var_T(i,je*px-4,k)
	    do j=2,je*px-1
	      ddy_1jepx(j)= ( 0.625+0.75*a )*var_T(i,j,k)
     &                   +( 0.25 +0.5 *a )*(var_T(i,j+1,k)+var_T(i,j-1,k))
     &                   +(-0.125+0.25*a )*(var_T(i,j+2,k)+var_T(i,j-2,k)) 
   	    enddo
            CALL solve_tridiag(filt_var_T(i,1:je*px,k),aay_1jepx,bby_1jepx,ccy_1jepx,ddy_1jepx,je*px)
	  enddo
	enddo
	Uvel_T(1:ie,1:je*px,1:ke/px)=filt_var_T(1:ie,1:je*px,1:ke/px)
	var_T(1:ie,1:je*px,1:ke/px)=Vvel_T(1:ie,1:je*px,1:ke/px)
	do i=ib,ie
	  do k=kb,ke/px
!	    ddy_1jepx(1)=(0.0625*a+0.9375)*var_T(i,1,k)+(0.25+0.75*a)*var_T(i,2,k)
!     &                 -0.375*(1-a)*var_T(i,3,k)+0.25*(1-a)*var_T(i,4,k)-0.0625*(1-a)*var_T(i,5,k)
	    ddy_1jepx(2)=(0.0625+0.875*a)*var_T(i,1,k)+(0.75+0.5*a)*var_T(i,2,k)
     &                 +(0.375+0.25*a)*var_T(i,3,k)-0.5*(0.5-a)*var_T(i,4,k)+0.125*(0.5-a)*var_T(i,5,k)
!	    ddy_1jepx(je*px)=(0.0625*a+0.9375)*var_T(i,je*px,k)+(0.25+0.75*a)*var_T(i,je*px-1,k)
!     &                 -0.375*(1-a)*var_T(i,je*px-2,k)+0.25*(1-a)*var_T(i,je*px-3,k)-0.0625*(1-a)*var_T(i,je*px-4,k)
	    ddy_1jepx(je*px-1)=(0.0625+0.875*a)*var_T(i,je*px,k)+(0.75+0.5*a)*var_T(i,je*px-1,k)
     &                 +(0.375+0.25*a)*var_T(i,je*px-2,k)-0.5*(0.5-a)*var_T(i,je*px-3,k)+0.125*(0.5-a)*var_T(i,je*px-4,k)
	    ddy_1jepx(1)=var_T(i,1,k)
	    ddy_1jepx(je*px)=var_T(i,je*px,k)

	    do j=3,je*px-2
	      ddy_1jepx(j)= ( 0.625+0.75*a )*var_T(i,j,k)
     &                   +( 0.25 +0.5 *a )*(var_T(i,j+1,k)+var_T(i,j-1,k))
     &                   +(-0.125+0.25*a )*(var_T(i,j+2,k)+var_T(i,j-2,k)) 
   	    enddo
            CALL solve_tridiag(filt_var_T(i,1:je*px,k),aay_1jepx,bby_1jepx,ccy_1jepx,ddy_1jepx,je*px)
	  enddo
	enddo
	Vvel_T(1:ie,1:je*px,1:ke/px)=filt_var_T(1:ie,1:je*px,1:ke/px)
	var_T(1:ie,0:je*px+1,1:ke/px)=Wvel_T(1:ie,0:je*px+1,1:ke/px)
	do i=ib,ie
	  do k=kb,ke/px
!	    ddy_1jepx(1)=(0.0625*a+0.9375)*var_T(i,1,k)+(0.25+0.75*a)*var_T(i,2,k)
!     &                 -0.375*(1-a)*var_T(i,3,k)+0.25*(1-a)*var_T(i,4,k)-0.0625*(1-a)*var_T(i,5,k)
!	    ddy_1jepx(je*px)=(0.0625*a+0.9375)*var_T(i,je*px,k)+(0.25+0.75*a)*var_T(i,je*px-1,k)
!     &                 -0.375*(1-a)*var_T(i,je*px-2,k)+0.25*(1-a)*var_T(i,je*px-3,k)-0.0625*(1-a)*var_T(i,je*px-4,k)
	    ddy_1jepx(1)=var_T(i,1,k)
	    ddy_1jepx(je*px)=var_T(i,je*px,k)

!	    ddy_1jepx(2)=(0.0625+0.875*a)*var_T(i,1,k)+(0.75+0.5*a)*var_T(i,2,k)
!     &                 +(0.375+0.25*a)*var_T(i,3,k)-0.5*(0.5-a)*var_T(i,4,k)+0.125*(0.5-a)*var_T(i,5,k)
!	    ddy_1jepx(je*px-1)=(0.0625+0.875*a)*var_T(i,je*px,k)+(0.75+0.5*a)*var_T(i,je*px-1,k)
!     &                 +(0.375+0.25*a)*var_T(i,je*px-2,k)-0.5*(0.5-a)*var_T(i,je*px-3,k)+0.125*(0.5-a)*var_T(i,je*px-4,k)
	    do j=2,je*px-1
	      ddy_1jepx(j)= ( 0.625+0.75*a )*var_T(i,j,k)
     &                   +( 0.25 +0.5 *a )*(var_T(i,j+1,k)+var_T(i,j-1,k))
     &                   +(-0.125+0.25*a )*(var_T(i,j+2,k)+var_T(i,j-2,k)) 
   	    enddo
            CALL solve_tridiag(filt_var_T(i,1:je*px,k),aay_1jepx,bby_1jepx,ccy_1jepx,ddy_1jepx,je*px)
	  enddo
	enddo
	Wvel_T(1:ie,1:je*px,1:ke/px)=filt_var_T(1:ie,1:je*px,1:ke/px)


	call t2fp(Uvel_T(1:ie,1:je*px,1:ke/px),Uvel(1:ie,1:je,1:ke))
	call t2fp(Vvel_T(1:ie,1:je*px,1:ke/px),Vvel(1:ie,1:je,1:ke))
	call t2fp(Wvel_T(1:ie,1:je*px,1:ke/px),Wvel(1:ie,1:je,1:ke))



	end


      subroutine LES_smagorinsky(Uvel,Vvel,Wvel,rr)
      USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      !include 'mpif.h'
	real shear
      real xx,yy,f,dzi,uu,vv,absU,ust,z0,yplus
      integer n,t,tel,kt
      real ebb(0:i1,0:k1)
      real ebf(0:i1,0:k1)
      real Lm2,Rudphi_i,Rumdphi_i,Rpdphi_i,dRpp_i,dRp_i,divergentie
      real     Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),rr(0:i1,0:j1,0:k1)

	real ekm2(0:i1,0:j1,0:k1)

!	real shr(1:imax,1:jmax,1:kmax),ekm1(1:imax,1:jmax,1:kmax),ekm2(1:imax,1:jmax,1:kmax)
!	real shr_T(1:imax,1:jmax*px,1:kmax/px),ekm2_T(1:imax,1:jmax*px,1:kmax/px),ekm2y(1:imax,1:jmax,1:kmax)

	shear = 0.

	dzi=1./dz
	do i=1,imax
	  Rudphi_i=1./(Ru(i)*dphi)
	  Rumdphi_i=1./(Ru(i-1)*dphi)
	  Rpdphi_i=1./(Rp(i)*dphi)
	  Lm2=Lmix(i)*Lmix(i)
	  dRpp_i=1./(Rp(i+1)-Rp(i))
	  dRp_i=1./(Rp(i)-Rp(i-1))
	  
	  
	  do j=1,jmax
	     do k=1,kmax

				shear = 2.0*(
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*dphi) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*dphi) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)) +	
     1			((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi) )
     
     
				shear = shear + 0.25*(
	1			((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))*Rudphi_i +
     2			( Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))*Rudphi_i +
     2			( Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))*Rumdphi_i +
     2			( Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2 +
	1			((Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))*Rumdphi_i +
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
     2			( Wvel(i  ,j+1,k  )-Wvel(i  ,j  ,k  ))*Rpdphi_i)**2 +
	1			((Vvel(i  ,j  ,k  )-Vvel(i  ,j  ,k-1))*dzi +
     2			( Wvel(i  ,j+1,k-1)-Wvel(i  ,j  ,k-1))*Rpdphi_i)**2 +
	1			((Vvel(i  ,j-1,k+1)-Vvel(i  ,j-1,k  ))*dzi +
     2			( Wvel(i  ,j  ,k  )-Wvel(i  ,j-1,k  ))*Rpdphi_i)**2 +
	1			((Vvel(i  ,j-1,k  )-Vvel(i  ,j-1,k-1))*dzi +
     2			( Wvel(i  ,j  ,k-1)-Wvel(i  ,j-1,k-1))*Rpdphi_i)**2 
     e				)

!! 
      divergentie= 2./3.*(( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) ) / ( dz ))   

		shear = shear 
     +                +1.5*divergentie*divergentie
     +                -divergentie*2.*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))
     +                -divergentie*2.*((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*dphi) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))
     +                -divergentie*2.*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)

				
				ekm(i,j,k) = rr(i,j,k) * Lm2 * sqrt(shear) ! Smagorinsky
     			enddo			
		enddo
	enddo
    

c*************************************************************
c	Van Driest damping of sgs viscosity
c	Determine ust, then yplus and adjust 
c	mu_t with mu_t = mu_t * (1 - exp(-y+/A))
c	A = 26 cf. Modeling and simulation of turbulent flows (Schiestel 2008) & Pope (2000) p.304
c*************************************************************

!       ust=0.1*U_b
!       do i=1,10
! 	z0=0.11*1e-6/ust+kn/30.
! 	ust=U_b*kappa/(log(((depth-bc_obst_h)-bc_obst_h)/z0)-1);
!       enddo   



 	if (slip_bot.eq.1) then
	ekm2=ekm
 	  do i=1,imax
 	    do j=1,jmax
  	      uu=0.5*(Uvel(i,j,1)+Uvel(i-1,j,1))
  	      vv=0.5*(Vvel(i,j,1)+Vvel(i,j-1,1))
  	      absU=sqrt(uu*uu+vv*vv)
  	      ust=0.1*absU
  	      do tel=1,10 ! 10 iter is more than enough
  		  z0=kn/30+0.11*nu_mol/MAX(ust,1.e-6)
  		  ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
  	      enddo
  	      yplus=0.5*dz*ust/nu_mol
  	      if (yplus<11.225) then	!viscous sublayer uplus=yplus
  		  ust=sqrt(absU*nu_mol/(0.5*dz))
  	      endif
  ! 	      write(*,*) 'ust van Driest damping',ust
 	      do k=1,10 ! first 10 grid layers from wall are adjusted
 		  yplus=(REAL(k)-0.5)*dz*ust/nu_mol
 		  ekm(i,j,k)=ekm(i,j,k) * (1. - exp(-yplus/26.))
 	      enddo
 	    enddo
 	  enddo
	       !! Set ekm jet back:
	      do t=1,tmax_inPpunt
		i=i_inPpunt(t)
		j=j_inPpunt(t)
		do k=k1-kjet+1,k1
		  ekm(i,j,k)=ekm2(i,j,k)
		enddo
	      enddo
 	endif

!!      Boundary conditions Neumann
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
	if (periodicx.eq.0) then
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
        ekm=ekm+ekm_mol
	ekm2=ekm
        IF (LOA>0.) THEN ! ship:
          do t=1,tmax_inPpuntTSHD
            i=i_inPpuntTSHD(t)
            j=j_inPpuntTSHD(t)
            k=k_inPpuntTSHD(t)
            ekm(i,j,k)=0.
          enddo
	ELSEIF(kjet>0) THEN
          ekm(0:i1,0:j1,kmax-kjet+1:k1)=0. !maak ekm in zone rondom buisje nul
	ENDIF

	if (LOA>0.and.kn_TSHD>0.) then
	  do t=1,tmax_inUpunt_tauTSHD
 	    i=i_inUpunt_tauTSHD(t)
 	    j=j_inUpunt_tauTSHD(t)		
 	    k=k_inUpunt_tauTSHD(t)		

		uu=0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))
  	      vv=0.5*(Vvel(i,j,k)+Vvel(i,j-1,k))
  	      absU=sqrt(uu*uu+vv*vv)
  	      ust=0.1*absU
  	      do tel=1,10 ! 10 iter is more than enough
  		  z0=kn/30+0.11*nu_mol/MAX(ust,1.e-6)
  		  ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
  	      enddo
  	      yplus=0.5*dz*ust/nu_mol
  	      if (yplus<11.225) then	!viscous sublayer uplus=yplus
  		  ust=sqrt(absU*nu_mol/(0.5*dz))
  	      endif
  ! 	      write(*,*) 'ust van Driest damping',ust
 	      do kt=k,k-10,-1 ! first 10 grid layers from wall are adjusted
 		  yplus=(REAL(kmax-kjet-kt)+0.5)*dz*ust/nu_mol
 		  ekm(i,j,kt)=ekm(i,j,kt) * (1. - exp(-yplus/26.))
 	      enddo


	  enddo
	endif



        do k=kmax-kjet-2,k1 ! laat ekm in buisje ongemoeid
            do t=1,tmax_inPpunt
              i=i_inPpunt(t)
              j=j_inPpunt(t)
              ekm(i,j,k)=ekm2(i,j,k)
            enddo
        enddo

        Diffcof=ekm/Sc/rr
        do k=kmax-kjet+1,k1 ! maak Diffcof rand buisje nul
          do t=1,tmax_inPpuntrand
            i=i_inPpuntrand(t)
            j=j_inPpuntrand(t)
            Diffcof(i,j,k)=0.
          enddo
        enddo
      end

      subroutine LES_mixinglengthdamped(Uvel,Vvel,Wvel,rr)
      USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      !include 'mpif.h'
	real shear
      real xx,yy,f,dzi,uu,vv,absU,ust,z0,yplus
      integer n,t,tel,kt
      real ebb(0:i1,0:k1)
      real ebf(0:i1,0:k1)
      real Lm2,Rudphi_i,Rumdphi_i,Rpdphi_i,dRpp_i,dRp_i,divergentie
      real     Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),rr(0:i1,0:j1,0:k1)

	real ekm2(0:i1,0:j1,0:k1)
	real z,Ri,dudz,dvdz,drdz,Diffcof2(0:i1,0:j1,0:k1),Lm2_Bak(1:kmax),MuAn_factor

!	real shr(1:imax,1:jmax,1:kmax),ekm1(1:imax,1:jmax,1:kmax),ekm2(1:imax,1:jmax,1:kmax)
!	real shr_T(1:imax,1:jmax*px,1:kmax/px),ekm2_T(1:imax,1:jmax*px,1:kmax/px),ekm2y(1:imax,1:jmax,1:kmax)

	shear = 0.
	MuAn_factor=0.
	if (damping_drho_dz.eq.'MuAn') then
	  MuAn_factor=1.
	endif

	do k=1,kmax
	  z=k*dz-0.5*dz 
	  Lm2_Bak(k)=(kappa*z)**2*(1.-z/depth)
	enddo 
		

	dzi=1./dz
	do i=1,imax
	  Rudphi_i=1./(Ru(i)*dphi)
	  Rumdphi_i=1./(Ru(i-1)*dphi)
	  Rpdphi_i=1./(Rp(i)*dphi)
	  Lm2=Lmix(i)*Lmix(i)
	  dRpp_i=1./(Rp(i+1)-Rp(i))
	  dRp_i=1./(Rp(i)-Rp(i-1))
	  
	  
	  do j=1,jmax
	     do k=1,kmax

				shear = 2.0*(
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*dphi) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*dphi) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)) +	
     1			((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi) )
     
     
				shear = shear + 0.25*(
	1			((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))*Rudphi_i +
     2			( Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))*Rudphi_i +
     2			( Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))*Rumdphi_i +
     2			( Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2 +
	1			((Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))*Rumdphi_i +
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
     2			( Wvel(i  ,j+1,k  )-Wvel(i  ,j  ,k  ))*Rpdphi_i)**2 +
	1			((Vvel(i  ,j  ,k  )-Vvel(i  ,j  ,k-1))*dzi +
     2			( Wvel(i  ,j+1,k-1)-Wvel(i  ,j  ,k-1))*Rpdphi_i)**2 +
	1			((Vvel(i  ,j-1,k+1)-Vvel(i  ,j-1,k  ))*dzi +
     2			( Wvel(i  ,j  ,k  )-Wvel(i  ,j-1,k  ))*Rpdphi_i)**2 +
	1			((Vvel(i  ,j-1,k  )-Vvel(i  ,j-1,k-1))*dzi +
     2			( Wvel(i  ,j  ,k-1)-Wvel(i  ,j-1,k-1))*Rpdphi_i)**2 
     e				)

!! 
      divergentie= 2./3.*(( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) ) / ( dz ))   

		shear = shear 
     +                +1.5*divergentie*divergentie
     +                -divergentie*2.*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))
     +                -divergentie*2.*((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*dphi) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))
     +                -divergentie*2.*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)

		dudz = ((Uvel(i  ,j  ,k+1)-Uvel(i  ,j  ,k  ))*dzi +
     &		        (Uvel(i  ,j  ,k  )-Uvel(i  ,j  ,k-1))*dzi +
     &		        (Uvel(i-1,j  ,k+1)-Uvel(i-1,j  ,k  ))*dzi +
     &		        (Uvel(i-1,j  ,k  )-Uvel(i-1,j  ,k-1))*dzi ) * 0.25
		dvdz = ((Vvel(i  ,j  ,k+1)-Vvel(i  ,j  ,k  ))*dzi +
     &		        (Vvel(i  ,j  ,k  )-Vvel(i  ,j  ,k-1))*dzi +
     &		        (Vvel(i  ,j-1,k+1)-Vvel(i  ,j-1,k  ))*dzi +
     &		        (Vvel(i  ,j-1,k  )-Vvel(i  ,j-1,k-1))*dzi ) * 0.25

		drdz=0.5*(rr(i,j,k+1  )-rr(i,j,k-1))*dzi
		Ri = -gz/rr(i,j,k)*drdz/(dudz**2+dvdz**2+1.e-18)
		Ri = MuAn_factor*MAX(0.,Ri)	
		ekm(i,j,k) = rr(i,j,k) * Lm2_Bak(k) * sqrt(shear) ! Mixing length model with Ri-damping

		ekm(i,j,k) = ekm(i,j,k) * (1.+damping_a1*Ri)**damping_b1
		Diffcof(i,j,k) = (ekm(i,j,k)+ekm_mol)/rr(i,j,k)/Sc * (1.+damping_a2*Ri)**damping_b2 ! damping diffusivity is by both damping functions
     			enddo			
		enddo
	enddo

!!      Boundary conditions Neumann
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
	if (periodicx.eq.0) then
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
        call shiftf(Diffcof,ebf) 
        call shiftb(Diffcof,ebb) 
	if (periodicy.eq.0.or.periodicy.eq.2) then
		if (rank.eq.0) then ! boundaries in j-direction
		        do k=1,kmax
		           do i=1,imax
		           Diffcof(i,0,k) = Diffcof(i,1,k) 
		           Diffcof(i,j1,k)= ebb(i,k) 
		           enddo
		        enddo
		elseif (rank.eq.px-1) then
		        do k=1,kmax
		           do i=1,imax
		           Diffcof(i,0,k) = ebf(i,k)
		           Diffcof(i,j1,k)= Diffcof(i,jmax,k) 
		           enddo
		        enddo   
		else
		        do k=1,kmax
		           do i=1,imax
		           Diffcof(i,0,k) = ebf(i,k)
		           Diffcof(i,j1,k) =ebb(i,k) 
		           enddo
		        enddo
		endif
	else
	        do k=1,kmax
	           do i=1,imax
	           Diffcof(i,0,k) = ebf(i,k)
	           Diffcof(i,j1,k) =ebb(i,k) 
	           enddo
	        enddo
	endif
	if (periodicx.eq.0) then
		do k=1,kmax ! boundaries in i-direction
		        do j=0,j1
		                Diffcof(0,j,k) = Diffcof(1,j,k)
		                Diffcof(i1,j,k) = Diffcof(imax,j,k)
		        enddo
		enddo
	else
		do k=1,kmax ! boundaries in i-direction
		        do j=0,j1
		                Diffcof(0,j,k) = Diffcof(imax,j,k)
		                Diffcof(i1,j,k) = Diffcof(1,j,k)
		        enddo
		enddo
	endif
        do i=0,i1 ! boundaries in k-direction
                do j=0,j1
                        Diffcof(i,j,0) = Diffcof(i,j,1)
                        Diffcof(i,j,k1) = Diffcof(i,j,kmax)
                enddo
        enddo
 
 
       ekm=ekm+ekm_mol
	ekm2=ekm
	Diffcof2=Diffcof
        IF (LOA>0.) THEN ! ship:
          do t=1,tmax_inPpuntTSHD
            i=i_inPpuntTSHD(t)
            j=j_inPpuntTSHD(t)
            k=k_inPpuntTSHD(t)
            ekm(i,j,k)=0.
	    Diffcof(i,j,k)=0.
          enddo
	ELSEIF(kjet>0) THEN
          ekm(0:i1,0:j1,kmax-kjet+1:k1)=0. !maak ekm in zone rondom buisje nul
          Diffcof(0:i1,0:j1,kmax-kjet+1:k1)=0. !maak ekm in zone rondom buisje nul
	ENDIF

	if (LOA>0.and.kn_TSHD>0.) then
	  do t=1,tmax_inUpunt_tauTSHD
 	    i=i_inUpunt_tauTSHD(t)
 	    j=j_inUpunt_tauTSHD(t)		
 	    k=k_inUpunt_tauTSHD(t)		

		uu=0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))
  	      vv=0.5*(Vvel(i,j,k)+Vvel(i,j-1,k))
  	      absU=sqrt(uu*uu+vv*vv)
  	      ust=0.1*absU
  	      do tel=1,10 ! 10 iter is more than enough
  		  z0=kn/30+0.11*nu_mol/MAX(ust,1.e-6)
  		  ust=absU/MAX(1./kappa*log(0.5*dz/z0),2.) !ust maximal 0.5*absU
  	      enddo
  	      yplus=0.5*dz*ust/nu_mol
  	      if (yplus<11.225) then	!viscous sublayer uplus=yplus
  		  ust=sqrt(absU*nu_mol/(0.5*dz))
  	      endif
  ! 	      write(*,*) 'ust van Driest damping',ust
 	      do kt=k,k-10,-1 ! first 10 grid layers from wall are adjusted
 		  yplus=(REAL(kmax-kjet-kt)+0.5)*dz*ust/nu_mol
 		  ekm(i,j,kt)=ekm(i,j,kt) * (1. - exp(-yplus/26.))
 		  Diffcof(i,j,kt)=Diffcof(i,j,kt) * (1. - exp(-yplus/26.))
 	      enddo
	  enddo
	endif
        do k=kmax-kjet-2,k1 ! laat ekm in buisje ongemoeid
            do t=1,tmax_inPpunt
              i=i_inPpunt(t)
              j=j_inPpunt(t)
              ekm(i,j,k)=ekm2(i,j,k)
	      Diffcof(i,j,k)=Diffcof2(i,j,k)
            enddo
        enddo
        do k=kmax-kjet+1,k1 ! maak Diffcof rand buisje nul
          do t=1,tmax_inPpuntrand
            i=i_inPpuntrand(t)
            j=j_inPpuntrand(t)
            Diffcof(i,j,k)=0.
          enddo
        enddo
      end


      subroutine LES_WALE(Uvel,Vvel,Wvel,rr)
      USE nlist
      implicit none
      !include 'mpif.h'
	real SijSij,SdijSdij
      real xx,yy,f,dzi,uu,vv,absU,ust,z0,yplus
      integer n,t,tel
      real ebb(0:i1,0:k1)
      real ebf(0:i1,0:k1)
      real Lm2,Rudphi_i,Rumdphi_i,Rpdphi_i,dRpp_i,dRp_i,divergentie
      real     Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),rr(0:i1,0:j1,0:k1)

	real ekm2(0:i1,0:j1,0:k1)
	real dudx,dudy,dudz
	real dvdx,dvdy,dvdz
	real dwdx,dwdy,dwdz
	real Sd11,Sd12,Sd13,Sd22,Sd23,Sd33
	real MuAn_factor,drdz,Ri 


	MuAn_factor=0.
	if (damping_drho_dz.eq.'MuAn') then
	  MuAn_factor=1.
	endif

	SijSij = 0.
	ekm=0.

	dzi=1./dz
	do i=1,imax
	  Rudphi_i=1./(Ru(i)*dphi)
	  Rumdphi_i=1./(Ru(i-1)*dphi)
	  Rpdphi_i=1./(Rp(i)*dphi)
	  Lm2=Lmix(i)*Lmix(i)
	  dRpp_i=1./(Rp(i+1)-Rp(i))
	  dRp_i=1./(Rp(i)-Rp(i-1))
	  
	  
	  do j=1,jmax
	     do k=1,kmax
		dudx = (Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i) 
		dvdy = (Vvel(i,j,k)-Vvel(i,j-1,k))*Rpdphi_i + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)
		dwdz = (Wvel(i,j,k)-Wvel(i,j,k-1))*dzi

		dudy = ((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))*Rudphi_i + 
     &		        (Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))*Rudphi_i +
     &		        (Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))*Rumdphi_i +
     &		        (Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))*Rumdphi_i ) * 0.25
		dudz = ((Uvel(i  ,j  ,k+1)-Uvel(i  ,j  ,k  ))*dzi +
     &		        (Uvel(i  ,j  ,k  )-Uvel(i  ,j  ,k-1))*dzi +
     &		        (Uvel(i-1,j  ,k+1)-Uvel(i-1,j  ,k  ))*dzi +
     &		        (Uvel(i-1,j  ,k  )-Uvel(i-1,j  ,k-1))*dzi ) * 0.25
		dvdx =(( Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) +
     &		       ( Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) + 
     &		       ( Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1)  +
     &		       ( Vvel(i  ,j-1,k  )/Rp(i)-Vvel(i-1,j-1,k  )/Rp(i-1))*dRp_i*Ru(i-1) ) * 0.25
		dvdz = ((Vvel(i  ,j  ,k+1)-Vvel(i  ,j  ,k  ))*dzi +
     &		        (Vvel(i  ,j  ,k  )-Vvel(i  ,j  ,k-1))*dzi +
     &		        (Vvel(i  ,j-1,k+1)-Vvel(i  ,j-1,k  ))*dzi +
     &		        (Vvel(i  ,j-1,k  )-Vvel(i  ,j-1,k-1))*dzi ) * 0.25
		dwdx =(( Wvel(i+1,j  ,k  )-Wvel(i  ,j  ,k  ))*dRpp_i + 
     &		       ( Wvel(i+1,j  ,k-1)-Wvel(i  ,j  ,k-1))*dRpp_i +
     &		       ( Wvel(i  ,j  ,k  )-Wvel(i-1,j  ,k  ))*dRp_i +
     &		       ( Wvel(i  ,j  ,k-1)-Wvel(i-1,j  ,k-1))*dRp_i ) * 0.25
		dwdy =(( Wvel(i  ,j+1,k  )-Wvel(i  ,j  ,k  ))*Rpdphi_i +
     &		       ( Wvel(i  ,j+1,k-1)-Wvel(i  ,j  ,k-1))*Rpdphi_i +
     &		       ( Wvel(i  ,j  ,k  )-Wvel(i  ,j-1,k  ))*Rpdphi_i +
     &		       ( Wvel(i  ,j  ,k-1)-Wvel(i  ,j-1,k-1))*Rpdphi_i ) * 0.25

		divergentie = 1./3.*(dudx**2+dvdy**2+dwdz**2)
		Sd11 =     (dudx*dudx + dudy*dvdx + dudz*dwdx) - divergentie 
                Sd22 =     (dvdx*dudy + dvdy*dvdy + dvdz*dwdy) - divergentie 
                Sd33 =     (dwdx*dudz + dwdy*dvdz + dwdz*dwdz) - divergentie
		Sd12 = 0.5*(dudx*dudy + dudy*dvdy + dudz*dwdy + dvdx*dudx + dvdy*dvdx + dvdz*dwdx)
                Sd13 = 0.5*(dudx*dudz + dudy*dvdz + dudz*dwdz + dwdx*dudx + dwdy*dvdx + dwdz*dwdx)
                Sd23 = 0.5*(dvdx*dudz + dvdy*dvdz + dvdz*dwdz + dwdx*dudy + dwdy*dvdy + dwdz*dwdy)

		SdijSdij = Sd11*Sd11 + Sd22*Sd22 + Sd33*Sd33 + 2.*Sd12*Sd12 + 2.*Sd13*Sd13 + 2.*Sd23*Sd23
 
!      divergentie= 1./3.*(( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
!     +              +
!     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) *Rpdphi_i 
!     +              +
!     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) ) *dzi)

		divergentie=0.

		Sd11 =     (dudx - divergentie)
		Sd22 =     (dvdy - divergentie)
		Sd33 =     (dwdz - divergentie)
		Sd12 = 0.5*(dvdx+dudy)
		Sd13 = 0.5*(dudz+dwdx)
		Sd23 = 0.5*(dvdz+dwdy)

		SijSij = Sd11*Sd11 + Sd22*Sd22 + Sd33*Sd33 + 2.*Sd12*Sd12 + 2.*Sd13*Sd13 + 2.*Sd23*Sd23

		ekm(i,j,k) = rr(i,j,k) * Lm2 * (SdijSdij**(1.5))/(1.e-12+SijSij**(2.5)+SdijSdij**(1.25))
		
!		drdz=0.5*(rr(i,j,k+1  )-rr(i,j,k-1))*dzi
!		Ri = -gz/rr(i,j,k)*drdz/(dudz**2+dvdz**2+1.e-18)
!		Ri = MuAn_factor*MAX(0.,Ri)	
!		ekm(i,j,k) = ekm(i,j,k) * (1.+damping_a1*Ri)**(-damping_b1)
     	    enddo			
	  enddo
	enddo

!!      Boundary conditions Neumann
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
	if (periodicx.eq.0) then
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


        ekm=ekm+ekm_mol
	ekm2=ekm
        IF (LOA>0.) THEN ! ship:
          do t=1,tmax_inPpuntTSHD
            i=i_inPpuntTSHD(t)
            j=j_inPpuntTSHD(t)
            k=k_inPpuntTSHD(t)
            ekm(i,j,k)=0.
          enddo
	ELSEIF (kjet>0) THEN
          ekm(0:i1,0:j1,kmax-kjet+1:k1)=0. !maak ekm in zone rondom buisje nul
	ENDIF

        do k=kmax-kjet+1,k1 ! laat ekm in buisje ongemoeid
            do t=1,tmax_inPpunt
              i=i_inPpunt(t)
              j=j_inPpunt(t)
              ekm(i,j,k)=ekm2(i,j,k)
            enddo
        enddo

        Diffcof=ekm/Sc/rr
        do k=kmax-kjet+1,k1 ! maak Diffcof rand buisje nul
          do t=1,tmax_inPpuntrand
            i=i_inPpuntrand(t)
            j=j_inPpuntrand(t)
            Diffcof(i,j,k)=0.
          enddo
        enddo

      end

      subroutine LES_Sigma(Uvel,Vvel,Wvel,rr)
      USE nlist
      implicit none
      !include 'mpif.h'
      real xx,yy,dzi
      integer n,t,tel
      real ebb(0:i1,0:k1)
      real ebf(0:i1,0:k1)
      real Lm2,Rudphi_i,Rumdphi_i,Rpdphi_i,dRpp_i,dRp_i,divergentie
      real     Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),rr(0:i1,0:j1,0:k1)

	real ekm2(0:i1,0:j1,0:k1)
	real dudx,dudy,dudz
	real dvdx,dvdy,dvdz
	real dwdx,dwdy,dwdz
	real G11,G22,G33,G12,G13,G23
	real s1,s2,s3,a1,a2,a3,II1,II2,II3
	real rrr,qqq,ppp,alpha,ss1,ss2,ss3,sm1,sm2


	dzi=1./dz
	do i=1,imax
	  Rudphi_i=1./(Ru(i)*dphi)
	  Rumdphi_i=1./(Ru(i-1)*dphi)
	  Rpdphi_i=1./(Rp(i)*dphi)
	  Lm2=Lmix(i)*Lmix(i)
	  dRpp_i=1./(Rp(i+1)-Rp(i))
	  dRp_i=1./(Rp(i)-Rp(i-1))
	  
	  
	  do j=1,jmax
	     do k=1,kmax
		dudx = (Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i) 
		dvdy = (Vvel(i,j,k)-Vvel(i,j-1,k))*Rpdphi_i + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)
		dwdz = (Wvel(i,j,k)-Wvel(i,j,k-1))*dzi

		dudy = ((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))*Rudphi_i + 
     &		        (Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))*Rudphi_i +
     &		        (Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))*Rumdphi_i +
     &		        (Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))*Rumdphi_i ) * 0.25
		dudz = ((Uvel(i  ,j  ,k+1)-Uvel(i  ,j  ,k  ))*dzi +
     &		        (Uvel(i  ,j  ,k  )-Uvel(i  ,j  ,k-1))*dzi +
     &		        (Uvel(i-1,j  ,k+1)-Uvel(i-1,j  ,k  ))*dzi +
     &		        (Uvel(i-1,j  ,k  )-Uvel(i-1,j  ,k-1))*dzi ) * 0.25
		dvdx =(( Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) +
     &		       ( Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) + 
     &		       ( Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1)  +
     &		       ( Vvel(i  ,j-1,k  )/Rp(i)-Vvel(i-1,j-1,k  )/Rp(i-1))*dRp_i*Ru(i-1) ) * 0.25
		dvdz = ((Vvel(i  ,j  ,k+1)-Vvel(i  ,j  ,k  ))*dzi +
     &		        (Vvel(i  ,j  ,k  )-Vvel(i  ,j  ,k-1))*dzi +
     &		        (Vvel(i  ,j-1,k+1)-Vvel(i  ,j-1,k  ))*dzi +
     &		        (Vvel(i  ,j-1,k  )-Vvel(i  ,j-1,k-1))*dzi ) * 0.25
		dwdx =(( Wvel(i+1,j  ,k  )-Wvel(i  ,j  ,k  ))*dRpp_i + 
     &		       ( Wvel(i+1,j  ,k-1)-Wvel(i  ,j  ,k-1))*dRpp_i +
     &		       ( Wvel(i  ,j  ,k  )-Wvel(i-1,j  ,k  ))*dRp_i +
     &		       ( Wvel(i  ,j  ,k-1)-Wvel(i-1,j  ,k-1))*dRp_i ) * 0.25
		dwdy =(( Wvel(i  ,j+1,k  )-Wvel(i  ,j  ,k  ))*Rpdphi_i +
     &		       ( Wvel(i  ,j+1,k-1)-Wvel(i  ,j  ,k-1))*Rpdphi_i +
     &		       ( Wvel(i  ,j  ,k  )-Wvel(i  ,j-1,k  ))*Rpdphi_i +
     &		       ( Wvel(i  ,j  ,k-1)-Wvel(i  ,j-1,k-1))*Rpdphi_i ) * 0.25


		! Sigma model from Nicoud et al. (2011) Using singular values to build a sgs model, Phys. Fluids 23
		! build G matrix Gij = gki*gkj
		! gki = du_k/dx_i

		G11 = dudx*dudx + dvdx*dvdx + dwdx*dwdx
		G22 = dudy*dudy + dvdy*dvdy + dwdy*dwdy
		G33 = dudz*dudz + dvdz*dvdz + dwdz*dwdz
		G12 = dudx*dudy + dvdx*dvdy + dwdx*dwdy
		G13 = dudx*dudz + dvdx*dvdz + dwdx*dwdz
		G23 = dudz*dudy + dvdz*dvdy + dwdz*dwdy

!!		METHOD FROM PAPER NICOUD ET AL (2011) TO GET EIGENVALUES IS NOT ROBUST ENOUGH (THEREFORE NOT USED):
!		ALGORITHM FROM WIKIPEDIA TO GET EIGENVALUES FROM 3X3 SYMMETRIC MATRIX IS ROBUST:
!		   p = A(1,2)^2 + A(1,3)^2 + A(2,3)^2
!		   q = trace(A)/3
!		   p = (A(1,1) - q)^2 + (A(2,2) - q)^2 + (A(3,3) - q)^2 + 2 * p
!		   p = sqrt(p / 6)
!		   B = (1 / p) * (A - q * I)       % I is the identity matrix
!		   r = det(B) / 2
!		 
!		   % In exact arithmetic for a symmetric matrix  -1 <= r <= 1
!		   % but computation error can leave it slightly outside this range.
!		   if (r <= -1) 
!		      phi = pi / 3
!		   elseif (r >= 1)
!		      phi = 0
!		   else
!		      phi = acos(r) / 3
!		   end
!		   % the eigenvalues satisfy eig3 <= eig2 <= eig1
!		   eig1 = q + 2 * p * cos(phi)
!		   eig3 = q + 2 * p * cos(phi + pi * (2/3))
!		   eig2 = 3 * q - eig1 - eig3     % since trace(A) = eig1 + eig2 + eig3
		qqq = (G11 + G22 + G33)/3.
		ppp = (G11-qqq)**2 + (G22-qqq)**2 + (G33-qqq)**2 + 2.*(G12**2+G13**2+G23**2)
		ppp = sqrt (ABS(ppp)/6.)
		ppp = MAX(ppp,1.e-25)

		G11 = G11-qqq
		G22 = G22-qqq
		G33 = G33-qqq
		G11 = G11/ppp
		G22 = G22/ppp
		G33 = G33/ppp
		G12 = G12/ppp
		G13 = G13/ppp
		G23 = G23/ppp
		rrr = 0.5*(G11*(G22*G33-G23*G23) - G12*(G12*G33-G23*G13) + G13*(G12*G23-G13*G22))

		! keep rrr between -1 and +1 to not get NaN from acos:
		rrr = MAX(rrr,-1.)
		rrr = MIN(rrr,1.)
		alpha = acos(rrr)/3.
		s1 = qqq + 2.*ppp*cos(alpha)
		s3 = qqq + 2.*ppp*cos(alpha+2.*pi/3)
		s2 = 3.*qqq - s1 - s3

		s1=sqrt(ABS(s1))
		s2=sqrt(ABS(s2))
		s3=sqrt(ABS(s3))
		! make sure that s1>s2>s3>0 otherwise negative ekm, therefore sort them,
		! Wiki algorithm fails in correct order when two eigenvalues are very small:
		sm1=MAX(s1,s2)
		sm2=MAX(s2,s3)
		ss1=MAX(sm1,sm2)
		ss2=MIN(sm1,sm2)
		ss3=MIN(s1,s2,s3) 

		ekm(i,j,k) = rr(i,j,k) * Lm2 * (ss3*(ss1-ss2)*(ss2-ss3))/(MAX(ss1*ss1,1.e-25))

    	    enddo			
	  enddo
	enddo
!!      Boundary conditions Neumann
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
	if (periodicx.eq.0) then
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

        ekm=ekm+ekm_mol
	ekm2=ekm
        IF (LOA>0.) THEN ! ship:
          do t=1,tmax_inPpuntTSHD
            i=i_inPpuntTSHD(t)
            j=j_inPpuntTSHD(t)
            k=k_inPpuntTSHD(t)
            ekm(i,j,k)=0.
          enddo
	ELSEIF (kjet>0) THEN
          ekm(0:i1,0:j1,kmax-kjet+1:k1)=0. !maak ekm in zone rondom buisje nul
	ENDIF

        do k=kmax-kjet+1,k1 ! laat ekm in buisje ongemoeid
            do t=1,tmax_inPpunt
              i=i_inPpunt(t)
              j=j_inPpunt(t)
              ekm(i,j,k)=ekm2(i,j,k)
            enddo
        enddo

        Diffcof=ekm/Sc/rr
        do k=kmax-kjet+1,k1 ! maak Diffcof rand buisje nul
          do t=1,tmax_inPpuntrand
            i=i_inPpuntrand(t)
            j=j_inPpuntrand(t)
            Diffcof(i,j,k)=0.
          enddo
        enddo

      end





      subroutine LES_filteredSmagorinsky(Uvel4,Vvel4,Wvel4,rr)
      USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      !include 'mpif.h'
	real shear
      real xx,yy,f,dzi,uu,vv,absU,ust,z0,yplus
      integer n,t,tel
      real ebb(0:i1,0:k1)
      real ebf(0:i1,0:k1)
      real Lm2,Rudphi_i,Rumdphi_i,Rpdphi_i,dRpp_i,dRp_i,divergentie
      real     Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),rr(0:i1,0:j1,0:k1)
	real ekm2(0:i1,0:j1,0:k1)

      real     Uvel4(0:i1,0:j1,0:k1),
     +         Vvel4(0:i1,0:j1,0:k1),Wvel4(0:i1,0:j1,0:k1)
	integer imm,jmm,kmm,ipp,jpp,kpp,ie,je,ke
	real*8 ubb(0:i1,0:k1),ubf(0:i1,0:k1)
	REAL Uvel2(-1:i1+1,-1:j1+1,-1:k1+1),Vvel2(-1:i1+1,-1:j1+1,-1:k1+1),Wvel2(-1:i1+1,-1:j1+1,-1:k1+1)

	Uvel=Uvel4
	Vvel=Vvel4
	Wvel=Wvel4

	!! HighPass filter according to 
	!! Ducros 1996 Large-eddy simulation of transition to turbulence in a boundary layer developing
	!! spatially over a flat plate
	!! Also used by Schluter 2000 in LES of JICF Re=82000 -> HP filtered Smag. better than standard Smag.
	!! HP filter is applied three times:

!	write(*,*),'Uvel,Vvel,Wvel(1,1,1):',Uvel(1,1,1),Vvel(1,1,1),Wvel(1,1,1)
!	write(*,*),'Uvel,Vvel,Wvel(0,0,0):',Uvel(0,0,0),Vvel(0,0,0),Wvel(0,0,0)
!	write(*,*),'Uvel,Vvel,Wvel(ie,je,ke):',Uvel(ie,je,ke),Vvel(ie,je,ke),Wvel(ie,je,ke)
	ie=imax
	je=jmax
	ke=kmax
	do t=1,nr_HPfilter
	Uvel2(0:ie,0:j1,0:k1)=Uvel(0:ie,0:j1,0:k1)
c get stuff from other CPU's
	  call shiftf2(Uvel,ubf)
	  call shiftb2(Uvel,ubb) 

	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1
		   Uvel2(i,-1,k) = Uvel(i,0,k)
		   Uvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1
		   Uvel2(i,-1,k) = Ubf(i,k)
		   Uvel2(i,j1+1,k) =Uvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1
		   Uvel2(i,-1,k) = Ubf(i,k)
		   Uvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	Uvel2(i1,0:j1,0:k1)=Uvel(ie,0:j1,0:k1)
	Uvel2(-1,0:j1,0:k1)=Uvel(0,0:j1,0:k1)
	Uvel2(i1+1,0:j1,0:k1)=Uvel(i1,0:j1,0:k1)
	Uvel2(0:i1,0:j1,-1)=Uvel(0:i1,0:j1,0)
	Uvel2(0:i1,0:j1,k1+1)=Uvel(0:i1,0:j1,k1)

	Vvel2(0:i1,0:j1,0:k1)=Vvel(0:i1,0:j1,0:k1)
c get stuff from other CPU's
	  call shiftf2(Vvel,ubf)
	  call shiftb2(Vvel,ubb) 

	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1
		   Vvel2(i,-1,k) = Vvel(i,0,k)
		   Vvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1
		   Vvel2(i,-1,k) = Ubf(i,k)
		   Vvel2(i,j1,k) =Vvel(i,je,k)
		   Vvel2(i,j1+1,k) =Vvel(i,je,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1
		   Vvel2(i,-1,k) = Ubf(i,k)
		   Vvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	Vvel2(-1,0:j1,0:k1)=Vvel(0,0:j1,0:k1)
	Vvel2(i1+1,0:j1,0:k1)=Vvel(i1,0:j1,0:k1)
	Vvel2(0:i1,0:j1,-1)=Vvel(0:i1,0:j1,0)
	Vvel2(0:i1,0:j1,k1+1)=Vvel(0:i1,0:j1,k1)

	Wvel2(0:i1,0:j1,0:k1)=Wvel(0:i1,0:j1,0:k1)
c get stuff from other CPU's
	  call shiftf2(Wvel,ubf)
	  call shiftb2(Wvel,ubb) 

	  if (rank.eq.0) then
		do k=0,k1
		   do i=0,i1
		   Wvel2(i,-1,k) = Wvel(i,0,k)
		   Wvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=0,k1
		   do i=0,i1
		   Wvel2(i,-1,k) = Ubf(i,k)
		   Wvel2(i,j1+1,k) =Wvel(i,j1,k)
		   enddo
		enddo
	  else 
		do k=0,k1
		   do i=0,i1
		   Wvel2(i,-1,k) = Ubf(i,k)
		   Wvel2(i,j1+1,k) =Ubb(i,k)
		   enddo
		enddo
	  endif
	Wvel2(-1,0:j1,0:k1)=Wvel(0,0:j1,0:k1)
	Wvel2(i1+1,0:j1,0:k1)=Wvel(i1,0:j1,0:k1)
	Wvel2(0:i1,0:j1,-1)=Wvel(0:i1,0:j1,0)
	Wvel2(0:i1,0:j1,k1)=Wvel(0:i1,0:j1,ke)
	Wvel2(0:i1,0:j1,k1+1)=Wvel(0:i1,0:j1,k1)
 	  do i=0,i1
!	    imm=MAX(i-2,-1)
!	    ipp=MIN(i+2,i1+1)
	    do j=0,j1
!	      jmm=MAX(j-2,-1)
!	      jpp=MIN(j+2,j1+1)
	      do k=0,k1
!	        kmm=MAX(k-2,-1)
!	        kpp=MIN(k+2,k1+1)
	!! 2nd diss: according to Sagaut, Comte, Ducros 2000 Filtered subgrid-scale models
		Uvel(i,j,k) = Uvel2(i-1,j,k)-2.*Uvel2(i,j,k)+Uvel2(i+1,j,k)+
     &		              Uvel2(i,j-1,k)-2.*Uvel2(i,j,k)+Uvel2(i,j+1,k)+ 
     &		              Uvel2(i,j,k-1)-2.*Uvel2(i,j,k)+Uvel2(i,j,k+1)
		Vvel(i,j,k) = Vvel2(i-1,j,k)-2.*Vvel2(i,j,k)+Vvel2(i+1,j,k)+
     &		              Vvel2(i,j-1,k)-2.*Vvel2(i,j,k)+Vvel2(i,j+1,k)+ 
     &		              Vvel2(i,j,k-1)-2.*Vvel2(i,j,k)+Vvel2(i,j,k+1)
		Wvel(i,j,k) = Wvel2(i-1,j,k)-2.*Wvel2(i,j,k)+Wvel2(i+1,j,k)+
     &		              Wvel2(i,j-1,k)-2.*Wvel2(i,j,k)+Wvel2(i,j+1,k)+ 
     &		              Wvel2(i,j,k-1)-2.*Wvel2(i,j,k)+Wvel2(i,j,k+1)

!		Uvel(i,j,k)=Uvel4(i,j,k)+0.25*Uvel(i,j,k)! low pass filtered 
!		Vvel(i,j,k)=Vvel4(i,j,k)+0.25*Vvel(i,j,k)! low pass filtered
!		Wvel(i,j,k)=Wvel4(i,j,k)+0.25*Wvel(i,j,k)! low pass filtered
		Uvel(i,j,k)=0.25*Uvel(i,j,k)! high pass filtered 
		Vvel(i,j,k)=0.25*Vvel(i,j,k)! high pass filtered
		Wvel(i,j,k)=0.25*Wvel(i,j,k)! high pass filtered

!		Uvel(i,j,k)=Uvel4(i,j,k)-Wvel(i,j,k)
!		Vvel(i,j,k)=Vvel4(i,j,k)-Wvel(i,j,k)
!		Wvel(i,j,k)=Wvel4(i,j,k)-Wvel(i,j,k)

!	!! 4th diss: according to Nicoud, Ducros 1999 sgs stress models based on the square of the velocity gradient tensor:
!		Uvel(i,j,k) = -Uvel2(imm,j,k)+4.*Uvel2(i-1,j,k)-6.*Uvel2(i,j,k)+4.*Uvel2(i+1,j,k)-Uvel2(ipp,j,k)
!     &		              -Uvel2(i,jmm,k)+4.*Uvel2(i,j-1,k)-6.*Uvel2(i,j,k)+4.*Uvel2(i,j+1,k)-Uvel2(i,jpp,k) 
!     &		              -Uvel2(i,j,kmm)+4.*Uvel2(i,j,k-1)-6.*Uvel2(i,j,k)+4.*Uvel2(i,j,k+1)-Uvel2(i,j,kpp)
!		Vvel(i,j,k) = -Vvel2(imm,j,k)+4.*Vvel2(i-1,j,k)-6.*Vvel2(i,j,k)+4.*Vvel2(i+1,j,k)-Vvel2(ipp,j,k)
!     &		              -Vvel2(i,jmm,k)+4.*Vvel2(i,j-1,k)-6.*Vvel2(i,j,k)+4.*Vvel2(i,j+1,k)-Vvel2(i,jpp,k) 
!     &		              -Vvel2(i,j,kmm)+4.*Vvel2(i,j,k-1)-6.*Vvel2(i,j,k)+4.*Vvel2(i,j,k+1)-Vvel2(i,j,kpp)
!		Wvel(i,j,k) = -Wvel2(imm,j,k)+4.*Wvel2(i-1,j,k)-6.*Wvel2(i,j,k)+4.*Wvel2(i+1,j,k)-Wvel2(ipp,j,k)
!     &		              -Wvel2(i,jmm,k)+4.*Wvel2(i,j-1,k)-6.*Wvel2(i,j,k)+4.*Wvel2(i,j+1,k)-Wvel2(i,jpp,k) 
!     &		              -Wvel2(i,j,kmm)+4.*Wvel2(i,j,k-1)-6.*Wvel2(i,j,k)+4.*Wvel2(i,j,k+1)-Wvel2(i,j,kpp)
	      enddo
	    enddo
	  enddo
	enddo

!	Uvel(:,:,1)=Uvel4(:,:,1)
!	Vvel(:,:,1)=Vvel4(:,:,1)
!	Wvel(:,:,1)=Wvel4(:,:,1)
!	Uvel(:,:,kmax)=Uvel4(:,:,kmax)
!	Vvel(:,:,kmax)=Vvel4(:,:,kmax)
!	Wvel(:,:,kmax)=Wvel4(:,:,kmax)
!	Uvel(1,:,:)=Uvel4(1,:,:)
!	Vvel(1,:,:)=Vvel4(1,:,:)
!	Wvel(1,:,:)=Wvel4(1,:,:)
!	Uvel(imax,:,:)=Uvel4(:,:,:)
!	Vvel(imax,:,:)=Vvel4(:,:,:)
!	Wvel(imax,:,:)=Wvel4(:,:,:)


	shear = 0.

	dzi=1./dz
	do i=1,imax
	  Rudphi_i=1./(Ru(i)*dphi)
	  Rumdphi_i=1./(Ru(i-1)*dphi)
	  Rpdphi_i=1./(Rp(i)*dphi)
	  Lm2=Lmix(i)*Lmix(i)
	  dRpp_i=1./(Rp(i+1)-Rp(i))
	  dRp_i=1./(Rp(i)-Rp(i-1))
	  
	  
	  do j=1,jmax
	     do k=1,kmax
				shear = 2.0*(
		!!	2S11S11+2S22S22+2S33S33	
     1			((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)) + 
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*dphi) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))*
     1			((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*dphi) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i)) +	
     1			((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi) )
     
     		!!	S12S12+S21S21=0.25*(4*S12S12)**2
				shear = shear + 0.25*(
	1			((Uvel(i  ,j+1,k  )-Uvel(i  ,j  ,k  ))*Rudphi_i +
     2			( Vvel(i+1,j  ,k  )/Rp(i+1)-Vvel(i  ,j  ,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i  ,j  ,k  )-Uvel(i  ,j-1,k  ))*Rudphi_i +
     2			( Vvel(i+1,j-1,k  )/Rp(i+1)-Vvel(i  ,j-1,k  )/Rp(i))*dRpp_i*Ru(i) )**2 +
	1			((Uvel(i-1,j+1,k  )-Uvel(i-1,j  ,k  ))*Rumdphi_i +
     2			( Vvel(i  ,j  ,k  )/Rp(i)-Vvel(i-1,j  ,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2 +
	1			((Uvel(i-1,j  ,k  )-Uvel(i-1,j-1,k  ))*Rumdphi_i +
     2			( Vvel(i  ,j-1,k  )/Rp(i)-Vvel(i-1,j-1,k  )/Rp(i-1))*dRp_i*Ru(i-1) )**2
     e				)
     		!!	S13S13+S31S31=0.25*(4*S31S31)**2
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
     		!!	S23S23+S32S32=0.25*(4*S23S23)**2
				shear = shear + 0.25*(
	1			((Vvel(i  ,j  ,k+1)-Vvel(i  ,j  ,k  ))*dzi +
     2			( Wvel(i  ,j+1,k  )-Wvel(i  ,j  ,k  ))*Rpdphi_i)**2 +
	1			((Vvel(i  ,j  ,k  )-Vvel(i  ,j  ,k-1))*dzi +
     2			( Wvel(i  ,j+1,k-1)-Wvel(i  ,j  ,k-1))*Rpdphi_i)**2 +
	1			((Vvel(i  ,j-1,k+1)-Vvel(i  ,j-1,k  ))*dzi +
     2			( Wvel(i  ,j  ,k  )-Wvel(i  ,j-1,k  ))*Rpdphi_i)**2 +
	1			((Vvel(i  ,j-1,k  )-Vvel(i  ,j-1,k-1))*dzi +
     2			( Wvel(i  ,j  ,k-1)-Wvel(i  ,j-1,k-1))*Rpdphi_i)**2 
     e				)

!! 
      divergentie= 2./3.*(( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) ) / ( dz ))   

		shear = shear 
     +                +1.5*divergentie*divergentie
     +                -divergentie*2.*((Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i))
     +                -divergentie*2.*((Vvel(i,j,k)-Vvel(i,j-1,k))/(Rp(i)*dphi) + 0.5*(Uvel(i,j,k)+Uvel(i-1,j,k))/Rp(i))
     +                -divergentie*2.*((Wvel(i,j,k)-Wvel(i,j,k-1))*dzi)

				
				ekm(i,j,k) = rr(i,j,k) * Lm2 * sqrt(shear) ! Smagorinsky
     			enddo			
		enddo
	enddo
    
!!	Boundary conditions Neumann
	call shiftf(ekm,ebf) 
	call shiftb(ekm,ebb) 
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
	do k=1,kmax ! boundaries in i-direction
		do j=0,j1
			ekm(0,j,k) = ekm(1,j,k)
			ekm(i1,j,k) = ekm(imax,j,k)
		enddo
	enddo
	do i=0,i1 ! boundaries in k-direction
		do j=0,j1
			ekm(i,j,0) = ekm(i,j,1)
			ekm(i,j,k1) = ekm(i,j,kmax)
		enddo
	enddo

        ekm=ekm+ekm_mol
	ekm2=ekm
        IF (LOA>0.) THEN ! ship:
          do t=1,tmax_inPpuntTSHD
            i=i_inPpuntTSHD(t)
            j=j_inPpuntTSHD(t)
            k=k_inPpuntTSHD(t)
            ekm(i,j,k)=0.
          enddo
	ELSEIF (kjet>0) THEN
          ekm(0:i1,0:j1,kmax-kjet+1:k1)=0. !maak ekm in zone rondom buisje nul
	ENDIF

        do k=kmax-kjet+1,k1 ! laat ekm in buisje ongemoeid
            do t=1,tmax_inPpunt
              i=i_inPpunt(t)
              j=j_inPpunt(t)
              ekm(i,j,k)=ekm2(i,j,k)
            enddo
        enddo

        Diffcof=ekm/Sc/rr
        do k=kmax-kjet+1,k1 ! maak Diffcof rand buisje nul
          do t=1,tmax_inPpuntrand
            i=i_inPpuntrand(t)
            j=j_inPpuntrand(t)
            Diffcof(i,j,k)=0.
          enddo
        enddo

      end


 
!       subroutine SEM_turb_bc(Ub1new,Vb1new,Wb1new,Ub2new,Vb2new,Wb2new,Ub1old,Vb1old,Wb1old,Ub2old,Vb2old,Wb2old)
      subroutine SEM_turb_bc
      USE nlist

      implicit none
!       include 'param.txt'
!       include 'common.txt'
!       include 'mpif.h'
!       integer modes
!       parameter (modes=150)
!       real Ub1new(0:i1,0:k1),Vb1new(0:i1,0:k1),Wb1new(0:i1,0:k1),Ub2new(0:j1,0:k1),Vb2new(0:j1+1,0:k1),Wb2new(0:j1,0:k1)
!       real Ub1old(0:i1,0:k1),Vb1old(0:i1,0:k1),Wb1old(0:i1,0:k1),Ub2old(0:j1,0:k1),Vb2old(0:j1+1,0:k1),Wb2old(0:j1,0:k1)
      integer ii,tt
      real fun,uu,vv,ww,fac1,fac2,x,y,z,yymin,yymax,Vbox1,Vbox2,boxside_x,boxside_y,phi
 	
      Ub1old=Ub1new
      Vb1old=Vb1new
      Wb1old=Wb1new      
      Ub2old=Ub2new
      Vb2old=Vb2new
      Wb2old=Wb2new 
!       fac=1./sqrt(REAL(nmax))*sqrt(1.5)*1.5
! 	  yymin=Rp(0)*sin((-0.5*jmax*px)*dphi)-boxsize
! 	  yymax=Rp(0)*sin((0.5*px*jmax)*dphi)+boxsize
!       Vbox=2.*boxsize*(yymax-yymin)*(depth-bc_obst_h)


      if (nmax2.gt.0) then
	  	yymin=Rp(0)*sin((-0.5*jmax*px)*dphi)-0.22*(depth-bc_obst_h) !boxsize
	  	yymax=Rp(0)*sin((0.5*px*jmax)*dphi)+0.22*(depth-bc_obst_h) !boxsize
      		Vbox2=1.2*(depth-bc_obst_h)*(yymax-yymin)*(depth-bc_obst_h)
      		fac2=1./sqrt(REAL(nmax2))*sqrt(1.5)*1.5

		phi=MAX(0.5*px*jmax*dphi,0.01/180.*pi) ! minimal 0.01 degrees
		boxside_x=0.6*(depth-bc_obst_h)+0.22*(depth-bc_obst_h)/tan(phi) !Lx_max=0.6 (depth-bc_obst_h), Ly_max=0.22 (depth-bc_obst_h)
		boxside_y=2*0.22*(depth-bc_obst_h)+(Rp(imax)-Rp(0))*sin(phi)
		Vbox1=2.*boxside_x*boxside_y*(depth-bc_obst_h)
		fac1=1./sqrt(REAL(nmax1))*sqrt(1.5)*1.5

      !! boundary at i=0
	      i=0
	      do k=kbed_bc+1,kmax
		z=(k-kbed_bc)*dz-0.5*dz
		do j=0,j1
		  x=Ru(0)*cos_u(j)-schuif_x
		  y=Ru(0)*sin_u(j)
		  uu=0.
		  vv=0.
		  ww=0.
		  do ii=1,llmax2(j,k)
		      tt=llist2(j,k,ii)
	! 	  do tt=1,nmax2
		      fun=sqrt(Vbox2/(lmxSEM2(tt)*lmySEM2(tt)*lmzSEM2(tt)))*(1.-MIN(ABS(x-xSEM2(tt))/lmxSEM2(tt),1.))
     & *(1.-MIN(ABS(y-ySEM2(tt))/lmySEM2(tt),1.))*(1.-MIN(ABS(z-zSEM2(tt))/lmzSEM2(tt),1.))
			uu=uu+epsSEM2(1,tt)*fun*fac2
			vv=vv+epsSEM2(2,tt)*fun*fac2
			ww=ww+epsSEM2(3,tt)*fun*fac2
		  enddo
		  do ii=1,llmax1(i,k)	!also interaction with SEM points at side
		      tt=llist1(i,k,ii)
		      fun=sqrt(Vbox1/(lmxSEM1(tt)*lmySEM1(tt)*lmzSEM1(tt)))*(1.-MIN(ABS(x-xSEM1(tt))/lmxSEM1(tt),1.))
     & *(1.-MIN(ABS(y-ySEM1(tt))/lmySEM1(tt),1.))*(1.-MIN(ABS(z-zSEM1(tt))/lmzSEM1(tt),1.))
			uu=uu+epsSEM1(1,tt)*fun*fac1
			vv=vv+epsSEM1(2,tt)*fun*fac1
			ww=ww+epsSEM1(3,tt)*fun*fac1
		  enddo

		    Ub2new(j,k)=uu*AA(1,1,k)+vv*AA(1,2,k)+ww*AA(1,3,k) 
		    Wb2new(j,k)=uu*AA(2,1,k)+vv*AA(2,2,k)+ww*AA(2,3,k)
		    Vb2new(j,k)=uu*AA(3,1,k)+vv*AA(3,2,k)+ww*AA(3,3,k)
		enddo
	      enddo
      endif


      if(nmax1.gt.0) then
	phi=MAX(0.5*px*jmax*dphi,0.01/180.*pi) ! minimal 0.01 degrees
	boxside_x=0.6*(depth-bc_obst_h)+0.22*(depth-bc_obst_h)/tan(phi) !Lx_max=0.6 (depth-bc_obst_h), Ly_max=0.22 (depth-bc_obst_h)
	boxside_y=2*0.22*(depth-bc_obst_h)+(Rp(imax)-Rp(0))*sin(phi)
        Vbox1=2.*boxside_x*boxside_y*(depth-bc_obst_h)
        fac1=1./sqrt(REAL(nmax1))*sqrt(1.5)*1.5

  	yymin=Rp(0)*sin((-0.5*jmax*px)*dphi)-0.22*(depth-bc_obst_h) !boxsize
  	yymax=Rp(0)*sin((0.5*px*jmax)*dphi)+0.22*(depth-bc_obst_h) !boxsize
	Vbox2=1.2*(depth-bc_obst_h)*(yymax-yymin)*(depth-bc_obst_h)
	fac2=1./sqrt(REAL(nmax2))*sqrt(1.5)*1.5
		

	       if (rank.eq.0) then	      !! boundary at j=0
	       j=0
	 	do k=kbed_bc+1,kmax
	 	  z=(k-kbed_bc)*dz-0.5*dz
	 	  do i=0,i1
	 	    x=Rp(i)*cos_v(0)-schuif_x
	 	    y=Rp(i)*sin_v(0)
	 	    uu=0.
	 	    vv=0.
	 	    ww=0.
	 	    do ii=1,llmax1(i,k)
	 		tt=llist1(i,k,ii)
		      fun=sqrt(Vbox1/(lmxSEM1(tt)*lmySEM1(tt)*lmzSEM1(tt)))*(1.-MIN(ABS(x-xSEM1(tt))/lmxSEM1(tt),1.))
     & *(1.-MIN(ABS(y-ySEM1(tt))/lmySEM1(tt),1.))*(1.-MIN(ABS(z-zSEM1(tt))/lmzSEM1(tt),1.))
	 		uu=uu+epsSEM1(1,tt)*fun*fac1
	 		vv=vv+epsSEM1(2,tt)*fun*fac1
	 		ww=ww+epsSEM1(3,tt)*fun*fac1
	 	    enddo
	 	    do ii=1,llmax2(j,k)  ! also interaction with SEM points at front inflow
	 		tt=llist2(j,k,ii)
		      fun=sqrt(Vbox2/(lmxSEM2(tt)*lmySEM2(tt)*lmzSEM2(tt)))*(1.-MIN(ABS(x-xSEM2(tt))/lmxSEM2(tt),1.))
     & *(1.-MIN(ABS(y-ySEM2(tt))/lmySEM2(tt),1.))*(1.-MIN(ABS(z-zSEM2(tt))/lmzSEM2(tt),1.))
	 		uu=uu+epsSEM2(1,tt)*fun*fac2
	 		vv=vv+epsSEM2(2,tt)*fun*fac2
	 		ww=ww+epsSEM2(3,tt)*fun*fac2
	 	    enddo
	 	    Ub1new(i,k)=uu*AA(1,1,k)+vv*AA(1,2,k)+ww*AA(1,3,k) 
	 	    Wb1new(i,k)=uu*AA(2,1,k)+vv*AA(2,2,k)+ww*AA(2,3,k)
	 	    Vb1new(i,k)=uu*AA(3,1,k)+vv*AA(3,2,k)+ww*AA(3,3,k)
	 	  enddo
	 	enddo   
	       ! boundary at j=px*jmax  
	       elseif (rank.eq.px-1) then
		j=jmax
	 	do k=kbed_bc+1,kmax
	 	  z=(k-kbed_bc)*dz-0.5*dz
	 	  do i=0,i1
	 	    x=Rp(i)*cos_v(jmax)-schuif_x
	 	    y=Rp(i)*sin_v(jmax)
	 	    uu=0.
	 	    vv=0.
	 	    ww=0.
	 	    do ii=1,llmax1(i,k)
	 		tt=llist1(i,k,ii)
		      fun=sqrt(Vbox1/(lmxSEM1(tt)*lmySEM1(tt)*lmzSEM1(tt)))*(1.-MIN(ABS(x-xSEM1(tt))/lmxSEM1(tt),1.))
     & *(1.-MIN(ABS(y-ySEM1(tt))/lmySEM1(tt),1.))*(1.-MIN(ABS(z-zSEM1(tt))/lmzSEM1(tt),1.))
	 		uu=uu+epsSEM1(1,tt)*fun*fac1
	 		vv=vv+epsSEM1(2,tt)*fun*fac1
	 		ww=ww+epsSEM1(3,tt)*fun*fac1
	 	    enddo
	 	    do ii=1,llmax2(j,k)  ! also interaction with SEM points at front inflow
	 		tt=llist2(j,k,ii)
		      fun=sqrt(Vbox2/(lmxSEM2(tt)*lmySEM2(tt)*lmzSEM2(tt)))*(1.-MIN(ABS(x-xSEM2(tt))/lmxSEM2(tt),1.))
     & *(1.-MIN(ABS(y-ySEM2(tt))/lmySEM2(tt),1.))*(1.-MIN(ABS(z-zSEM2(tt))/lmzSEM2(tt),1.))
	 		uu=uu+epsSEM2(1,tt)*fun*fac2
	 		vv=vv+epsSEM2(2,tt)*fun*fac2
	 		ww=ww+epsSEM2(3,tt)*fun*fac2
	 	    enddo
	 	    Ub1new(i,k)=uu*AA(1,1,k)+vv*AA(1,2,k)+ww*AA(1,3,k)
	 	    Wb1new(i,k)=uu*AA(2,1,k)+vv*AA(2,2,k)+ww*AA(2,3,k)
	 	    Vb1new(i,k)=uu*AA(3,1,k)+vv*AA(3,2,k)+ww*AA(3,3,k)
	 	  enddo
	 	enddo
	       endif
	endif
      end


      subroutine create_init_eddies_SEM
      USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      include 'mpif.h'

      integer m,n,ierr,clock,ii
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      real yy(nmax2),xx(nmax2),y,z
      real z0,xxmin,xxmax,yymin,yymax,zzmin,zzmax,x0,y0,phi,ust,fac,phi2
      character*60 fmatname
      real boxside_x

	CALL SYSTEM_CLOCK(COUNT=clock)
	CALL RANDOM_SEED(size = n)
	ALLOCATE(seed(n))
	CALL SYSTEM_CLOCK(COUNT=clock)
	seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	CALL RANDOM_SEED(PUT = seed)


	ust=sqrt(U_bSEM**2+V_bSEM**2)
	if (slip_bot.eq.1) then
	  do i=1,10
	    z0=0.11*nu_mol/ust+kn/30
	    ust=sqrt(U_bSEM**2+V_bSEM**2)*kappa/(log((depth-bc_obst_h)/z0)-1);
	  enddo
	endif
	! rotation ship for ambient side current
	if ((U_TSHD-U_bSEM).eq.0.or.LOA<0.) then
	  phi=atan2(V_bSEM,1.e-12)
	else
	  phi=atan2(V_bSEM,(U_TSHD-U_bSEM))
	endif


	if (rank.eq.1) then ! 2nd processor can generate SEM2
	  call random_number(zSEM2) ! uniform distribution 0,1
	  call random_number(ySEM2) ! uniform distribution 0,1
	  call random_number(xSEM2) ! uniform distribution 0,1
	  call random_number(epsSEM2) ! uniform distribution 0,1
	  epsSEM2=anint(epsSEM2)
	  epsSEM2=epsSEM2*2.-1. ! +1 or -1

	  zzmin=z0+1.e-6
	  zzmax=(depth-bc_obst_h)
	  zSEM2=(zzmax-zzmin)*zSEM2+zzmin
	  yymin=Rp(0)*sin((-0.5*jmax*px)*dphi)-0.22*(depth-bc_obst_h) !boxsize
	  yymax=Rp(0)*sin((0.5*px*jmax)*dphi)+0.22*(depth-bc_obst_h) !boxsize
	  ySEM2=(yymax-yymin)*ySEM2+yymin
	  xxmin=Rp(0)*cos(0.5*jmax*px*dphi)-schuif_x-0.6*(depth-bc_obst_h) !boxsize
	  xxmax=Rp(0)-schuif_x+0.6*(depth-bc_obst_h) !boxsize
	  xSEM2=(xxmax-xxmin)*xSEM2+xxmin
	  fac=kappa*(depth-bc_obst_h)
	  do i=1,nmax2
!	      uSEM2(i)=ust/kappa*log(zSEM2(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)  !correct sign with and without TSHD
	      uSEM2(i)=ust/kappa*log(zSEM2(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)+U_TSHD  !correct sign with and without TSHD
	      lmxSEM2(i)=0.5*(depth-bc_obst_h)/sqrt(U_bSEM**2+V_bSEM**2)*ust/kappa*log(zSEM2(i)/z0)
	      lmxSEM2(i)=max(lm_min,lmxSEM2(i))
	      lmySEM2(i)=0.22*(depth-bc_obst_h)
	      lmzSEM2(i)=fac*zSEM2(i)/(depth-bc_obst_h)*sqrt(1.-zSEM2(i)/(depth-bc_obst_h))
	      lmzSEM2(i)=max(lm_min,lmzSEM2(i))
	  enddo
	endif
	call mpi_bcast(xSEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(ySEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(zSEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(uSEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(epsSEM2,3*nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(lmxSEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(lmySEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(lmzSEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	
!!	make linked list on every node for SEM2:
	llmax2=0
 	do k=kbed_bc+1,kmax
 	  z=(k-kbed_bc)*dz-0.5*dz
	  do j=0,j1
	    y=Ru(0)*sin_u(j)
	    do ii=1,nmax2
	      if(ABS(y-ySEM2(ii))/lmySEM2(ii).lt.1..and.ABS(z-zSEM2(ii))/lmzSEM2(ii).lt.1.) then
		llmax2(j,k)=llmax2(j,k)+1
		llist2(j,k,llmax2(j,k))=ii
	      endif
	    enddo
	  enddo
	enddo

 	if (rank.eq.0) then       !! boundary at j=0
 	  call random_number(zSEM1) ! uniform distribution 0,1
 	  call random_number(xSEM1) ! uniform distribution 0,1
 	  call random_number(ySEM1) ! uniform distribution 0,1
 	  call random_number(epsSEM1) ! uniform distribution 0,1
 	  epsSEM1=anint(epsSEM1)
 	  epsSEM1=epsSEM1*2.-1. ! +1 or -1
 
 	  zzmin=z0+1.e-6
 	  zzmax=(depth-bc_obst_h)
 	  zSEM1=(zzmax-zzmin)*zSEM1+zzmin
	
	  phi2=MAX(0.5*px*jmax*dphi,0.01/180.*pi) ! minimal 0.01 degrees (phi2 is positive...)
	  yymin=Rp(0)*sin(phi2)-0.22*(depth-bc_obst_h)
	  yymax=Rp(imax)*sin(phi2)+0.22*(depth-bc_obst_h)
	  ySEM1=-(yymax-yymin)*ySEM1-yymin  !! boundary at j=0, so y is negative!!
	  boxside_x=0.6*(depth-bc_obst_h)+0.22*(depth-bc_obst_h)/tan(phi2) !Lx_max=0.6 (depth-bc_obst_h), Ly_max=0.22 (depth-bc_obst_h)
	  fac=kappa*(depth-bc_obst_h)
 	  do i=1,nmax1
	      xxmin=-ySEM1(i)/tan(phi2)-schuif_x-boxside_x ! xxmin is positive
	      xxmax=-ySEM1(i)/tan(phi2)-schuif_x+boxside_x ! xxmax is positive
	      xSEM1(i)=(xxmax-xxmin)*xSEM1(i)+xxmin  
!	      uSEM1(i)=ust/kappa*log(zSEM1(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)  !correct sign with and without TSHD
	      uSEM1(i)=ust/kappa*log(zSEM1(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)+U_TSHD  !correct sign with and without TSHD
	      lmxSEM1(i)=0.5*(depth-bc_obst_h)/sqrt(U_bSEM**2+V_bSEM**2)*ust/kappa*log(zSEM1(i)/z0)
	      lmxSEM1(i)=max(lm_min,lmxSEM1(i))
	      lmySEM1(i)=0.22*(depth-bc_obst_h)
	      lmzSEM1(i)=fac*zSEM1(i)/(depth-bc_obst_h)*sqrt(1.-zSEM1(i)/(depth-bc_obst_h))
	      lmzSEM1(i)=max(lm_min,lmzSEM1(i))
 	  enddo
 
   !!	make linked list for SEM1:
 	  llmax1=0
 	do k=kbed_bc+1,kmax
 	  z=(k-kbed_bc)*dz-0.5*dz
 	    do i=0,i1
 	      y=Rp(i)*sin_v(0)
 	      do ii=1,nmax1
 		if(ABS(y-ySEM1(ii))/lmySEM1(ii).lt.1..and.ABS(z-zSEM1(ii))/lmzSEM1(ii).lt.1.) then
 		  llmax1(i,k)=llmax1(i,k)+1
 		  llist1(i,k,llmax1(i,k))=ii
 		endif
 	      enddo
 	    enddo
 	  enddo
 	endif

 	if (rank.eq.px-1) then
 	  call random_number(zSEM1) ! uniform distribution 0,1
 	  call random_number(xSEM1) ! uniform distribution 0,1
 	  call random_number(ySEM1) ! uniform distribution 0,1
 	  call random_number(epsSEM1) ! uniform distribution 0,1
 	  epsSEM1=anint(epsSEM1)
 	  epsSEM1=epsSEM1*2.-1. ! +1 or -1
 
 	  zzmin=z0+1.e-6
 	  zzmax=(depth-bc_obst_h)
 	  zSEM1=(zzmax-zzmin)*zSEM1+zzmin
	
	  phi2=MAX(0.5*px*jmax*dphi,0.01/180.*pi) ! minimal 0.01 degrees (phi2 is positive...)
	  yymin=Rp(0)*sin(phi2)-0.22*(depth-bc_obst_h)
	  yymax=Rp(imax)*sin(phi2)+0.22*(depth-bc_obst_h)
	  ySEM1=(yymax-yymin)*ySEM1+yymin  !! boundary at j=jmax, so y is positive !!
	  boxside_x=0.6*(depth-bc_obst_h)+0.22*(depth-bc_obst_h)/tan(phi2) !Lx_max=0.6 (depth-bc_obst_h), Ly_max=0.22 (depth-bc_obst_h)
	  fac=kappa*(depth-bc_obst_h)
 	  do i=1,nmax1
	      xxmin=ySEM1(i)/tan(phi2)-schuif_x-boxside_x 
	      xxmax=ySEM1(i)/tan(phi2)-schuif_x+boxside_x
	      xSEM1(i)=(xxmax-xxmin)*xSEM1(i)+xxmin  
!	      uSEM1(i)=ust/kappa*log(zSEM1(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)  !correct sign with and without TSHD
	      uSEM1(i)=ust/kappa*log(zSEM1(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)+U_TSHD  !correct sign with and without TSHD
	      lmxSEM1(i)=0.5*(depth-bc_obst_h)/sqrt(U_bSEM**2+V_bSEM**2)*ust/kappa*log(zSEM1(i)/z0)
	      lmxSEM1(i)=max(lm_min,lmxSEM1(i))
	      lmySEM1(i)=0.22*(depth-bc_obst_h)
	      lmzSEM1(i)=fac*zSEM1(i)/(depth-bc_obst_h)*sqrt(1.-zSEM1(i)/(depth-bc_obst_h))
	      lmzSEM1(i)=max(lm_min,lmzSEM1(i))
 	  enddo
 
   !!	make linked list for SEM1:
 	  llmax1=0
 	do k=kbed_bc+1,kmax
 	  z=(k-kbed_bc)*dz-0.5*dz
 	    do i=0,i1
 	      y=Rp(i)*sin_v(jmax)
 	      do ii=1,nmax1
 		if(ABS(y-ySEM1(ii))/lmySEM1(ii).lt.1..and.ABS(z-zSEM1(ii))/lmzSEM1(ii).lt.1.) then
 		  llmax1(i,k)=llmax1(i,k)+1
 		  llist1(i,k,llmax1(i,k))=ii
 		endif
 	      enddo
 	    enddo
 	  enddo
 	endif


	!! create Cholesky decomposition of Re stress tensor
	if (LOA<0.and.kjet>0) then
        do k=1,kmax-kjet
	 if(k.gt.kbed_bc) then
	  IF (wallup.eq.1) THEN
	    z=depth-((k-kbed_bc)*dz-0.5*dz)
	  ELSE
	    z=(k-kbed_bc)*dz-0.5*dz
	  ENDIF
          call Chol_tensor_from_DNS_Ret395(z/(depth-kbed_bc*dz),AA(:,:,k))
          AA(:,:,k)=AA(:,:,k)*ust ! scale back with ust
	 else 
	  AA(:,:,k)=0.
	 endif
        enddo

	else

	do k=1,kmax
	 if(k.gt.kbed_bc) then
	  IF (wallup.eq.1) THEN
	    z=depth-((k-kbed_bc)*dz-0.5*dz)
	  ELSE
	    z=(k-kbed_bc)*dz-0.5*dz
	  ENDIF
	  call Chol_tensor_from_DNS_Ret395(z/(depth-kbed_bc*dz),AA(:,:,k))
	  AA(:,:,k)=AA(:,:,k)*ust ! scale back with ust
	 else 
	  AA(:,:,k)=0. 
	 endif
	enddo
	endif

	end


      subroutine move_eddies_SEM
      USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      include 'mpif.h'

      integer n,clock,ierr
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      real zz,eps(3),yy
      real z0,xxmin,xxmax,yymin,yymax,x0,y0,phi,ust,fac,phi2
      integer move
      character*60 fmatname
      integer jminSEM,jmaxSEM,kminSEM,kmaxSEM,iii,tel,ii
      real phiSEM,dphiSEM,y,z,zzmin,zzmax,lmSEM2old(1:nmax2),ySEM2old(1:nmax2),zSEM2old(1:nmax2),xSEM2old(1:nmax2)
      real lmxSEM2old(1:nmax2),lmySEM2old(1:nmax2),lmzSEM2old(1:nmax2)
      integer ind(nmax2),iimax
      real boxside_x
      

      CALL SYSTEM_CLOCK(COUNT=clock)
      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))
      CALL SYSTEM_CLOCK(COUNT=clock)
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)


      fac=kappa*(depth-bc_obst_h)

	ust=sqrt(U_bSEM**2+V_bSEM**2)
	if (slip_bot.eq.1) then
	  do i=1,10
	    z0=0.11*nu_mol/ust+kn/30
	    ust=sqrt(U_bSEM**2+V_bSEM**2)*kappa/(log((depth-bc_obst_h)/z0)-1);
	  enddo
	endif
	! rotation ship for ambient side current
	if ((U_TSHD-U_bSEM).eq.0.or.LOA<0.) then
	  phi=atan2(V_bSEM,1.e-12)
	else
	  phi=atan2(V_bSEM,(U_TSHD-U_bSEM))
	endif

      phi2=MAX(0.5*px*jmax*dphi,0.01/180.*pi) ! minimal 0.01 degrees (phi2 is positive...)
      boxside_x=0.6*(depth-bc_obst_h)+0.22*(depth-bc_obst_h)/tan(phi2) !Lx_max=0.6 (depth-bc_obst_h), Ly_max=0.22 (depth-bc_obst_h)	
      zzmin=z0+1.e-6
      zzmax=(depth-bc_obst_h)
      yymin=Rp(0)*sin(phi2)-0.22*(depth-bc_obst_h)
      yymax=Rp(imax)*sin(phi2)+0.22*(depth-bc_obst_h)
      if (rank.eq.0) then
   !! 	grens op j=0
 	do i=1,nmax1
 	  xSEM1(i)=xSEM1(i)+uSEM1(i)*dt
	  xxmax=-ySEM1(i)/tan(phi2)-schuif_x+boxside_x ! xxmax is positive
 	  if (xSEM1(i)>xxmax)  then  ! put on start line inflow Vbox:
		!! first remove from linked list:
		  kminSEM=INT(FLOOR((zSEM1(i)-lmzSEM1(i))/dz))+kbed_bc
		  kmaxSEM=INT(CEILING((zSEM1(i)+lmzSEM1(i))/dz))+kbed_bc
		  kminSEM=MAX(kminSEM,1)
		  kminSEM=MIN(kminSEM,kmax)
		  kmaxSEM=MAX(kmaxSEM,1)
		  kmaxSEM=MIN(kmaxSEM,kmax)
		  do k=kminSEM,kmaxSEM
	 	    z=(k-kbed_bc)*dz-0.5*dz
		    do ii=0,i1 
		      y=Rp(ii)*sin_v(0)	
		      if(ABS(y-ySEM1(i))/lmySEM1(i).lt.1..and.ABS(z-zSEM1(i))/lmzSEM1(i).lt.1.) then
			do tel=1,llmax1(ii,k)
			  if(llist1(ii,k,tel).eq.i) then
			    llist1(ii,k,tel)=llist1(ii,k,llmax1(ii,k))
			    llmax1(ii,k)=llmax1(ii,k)-1
			  endif
			enddo
		      endif
		    enddo
		  enddo
		! create new position/eps:
		  call random_number(zz) ! uniform distribution 0,1
	 	  call random_number(yy) ! uniform distribution 0,1
	  	  zSEM1(i)=(zzmax-zzmin)*zz+zzmin
		  ySEM1(i)=-(yymax-yymin)*yy-yymin  !! boundary at j=0, so y is negative!!
	 	  xSEM1(i)=-ySEM1(i)/tan(phi2)-schuif_x-boxside_x
	 	  call random_number(eps)
	 	  eps=anint(eps)
	 	  do j=1,3
	 	    epsSEM1(i,j)=eps(j)*2.-1. ! +1 or -1	
	 	  enddo
!	      	  uSEM1(i)=ust/kappa*log(zSEM1(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)  !correct sign with and without TSHD
	      	  uSEM1(i)=ust/kappa*log(zSEM1(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)+U_TSHD  !correct sign with and without TSHD
		  lmxSEM1(i)=0.5*(depth-bc_obst_h)/sqrt(U_bSEM**2+V_bSEM**2)*ust/kappa*log(zSEM1(i)/z0)
		  lmxSEM1(i)=max(lm_min,lmxSEM1(i))
		  lmySEM1(i)=0.22*(depth-bc_obst_h)
		  lmzSEM1(i)=fac*zSEM1(i)/(depth-bc_obst_h)*sqrt(1.-zSEM1(i)/(depth-bc_obst_h))
		  lmzSEM1(i)=max(lm_min,lmzSEM1(i))
		!! make linked list for new SEM point:
		  kminSEM=INT(FLOOR((zSEM1(i)-lmzSEM1(i))/dz))+kbed_bc
		  kmaxSEM=INT(CEILING((zSEM1(i)+lmzSEM1(i))/dz))+kbed_bc
		  kminSEM=MAX(kminSEM,1)
		  kminSEM=MIN(kminSEM,kmax)
		  kmaxSEM=MAX(kmaxSEM,1)
		  kmaxSEM=MIN(kmaxSEM,kmax)
	 	  do k=kminSEM,kmaxSEM
	 	    z=(k-kbed_bc)*dz-0.5*dz
	 	    do ii=0,i1
	 	      y=Rp(ii)*sin_v(0)
	 	      if(ABS(y-ySEM1(i))/lmySEM1(i).lt.1..and.ABS(z-zSEM1(i))/lmzSEM1(i).lt.1.) then
	 		  llmax1(ii,k)=llmax1(ii,k)+1
	 		  llist1(ii,k,llmax1(ii,k))=i
	 	      endif
	 	    enddo
	 	  enddo
 	  endif
 	enddo
       endif
       if (rank.eq.px-1) then
 !! 	grens op j=jmax
 	do i=1,nmax1
 	  xSEM1(i)=xSEM1(i)+uSEM1(i)*dt
	  xxmax=ySEM1(i)/tan(phi2)-schuif_x+boxside_x ! xxmax is positive
 	  if (xSEM1(i)>xxmax)  then  ! put on start line inflow Vbox:
		!! first remove from linked list:
		  kminSEM=INT(FLOOR((zSEM1(i)-lmzSEM1(i))/dz))+kbed_bc
		  kmaxSEM=INT(CEILING((zSEM1(i)+lmzSEM1(i))/dz))+kbed_bc
		  kminSEM=MAX(kminSEM,1)
		  kminSEM=MIN(kminSEM,kmax)
		  kmaxSEM=MAX(kmaxSEM,1)
		  kmaxSEM=MIN(kmaxSEM,kmax)
		  do k=kminSEM,kmaxSEM
 	 	    z=(k-kbed_bc)*dz-0.5*dz
		    do ii=0,i1 
		      y=Rp(ii)*sin_v(jmax)	
		      if(ABS(y-ySEM1(i))/lmySEM1(i).lt.1..and.ABS(z-zSEM1(i))/lmzSEM1(i).lt.1.) then
			do tel=1,llmax1(ii,k)
			  if(llist1(ii,k,tel).eq.i) then
			    llist1(ii,k,tel)=llist1(ii,k,llmax1(ii,k))
			    llmax1(ii,k)=llmax1(ii,k)-1
			  endif
			enddo
		      endif
		    enddo
		  enddo
		  call random_number(zz) ! uniform distribution 0,1
	 	  call random_number(yy) ! uniform distribution 0,1
	  	  zSEM1(i)=(zzmax-zzmin)*zz+zzmin
		  ySEM1(i)=(yymax-yymin)*yy+yymin  !! boundary at j=jmax, so y is positive!!
	 	  xSEM1(i)=ySEM1(i)/tan(phi2)-schuif_x-boxside_x
	 	  call random_number(eps)
	 	  eps=anint(eps)
	 	  do j=1,3
	 	    epsSEM1(i,j)=eps(j)*2.-1. ! +1 or -1	
	 	  enddo
!	      	  uSEM1(i)=ust/kappa*log(zSEM1(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)  !correct sign with and without TSHD
	      	  uSEM1(i)=ust/kappa*log(zSEM1(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)+U_TSHD  !correct sign with and without TSHD
		  lmxSEM1(i)=0.5*(depth-bc_obst_h)/sqrt(U_bSEM**2+V_bSEM**2)*ust/kappa*log(zSEM1(i)/z0)
		  lmxSEM1(i)=max(lm_min,lmxSEM1(i))
		  lmySEM1(i)=0.22*(depth-bc_obst_h)
		  lmzSEM1(i)=fac*zSEM1(i)/(depth-bc_obst_h)*sqrt(1.-zSEM1(i)/(depth-bc_obst_h))
		  lmzSEM1(i)=max(lm_min,lmzSEM1(i))
		!! make linked list for new SEM point:
		  kminSEM=INT(FLOOR((zSEM1(i)-lmzSEM1(i))/dz))+kbed_bc
		  kmaxSEM=INT(CEILING((zSEM1(i)+lmzSEM1(i))/dz))+kbed_bc
		  kminSEM=MAX(kminSEM,1)
		  kminSEM=MIN(kminSEM,kmax)
		  kmaxSEM=MAX(kmaxSEM,1)
		  kmaxSEM=MIN(kmaxSEM,kmax)
	 	  do k=kminSEM,kmaxSEM
	 	    z=(k-kbed_bc)*dz-0.5*dz
	 	    do ii=0,i1
	 	      y=Rp(ii)*sin_v(jmax)
	 	      if(ABS(y-ySEM1(i))/lmySEM1(i).lt.1..and.ABS(z-zSEM1(i))/lmzSEM1(i).lt.1.) then
	 		  llmax1(ii,k)=llmax1(ii,k)+1
	 		  llist1(ii,k,llmax1(ii,k))=i
	 	      endif
	 	    enddo
	 	  enddo
 	  endif
 	enddo
       endif


      move=0
      iimax=0
!! 	grens op i=0
      xxmin=Rp(0)*cos(0.5*jmax*px*dphi)-schuif_x-0.6*(depth-bc_obst_h) !boxsize
      xxmax=Rp(0)-schuif_x+ 0.6*(depth-bc_obst_h) !boxsize
      zzmin=z0+1.e-6
      zzmax=(depth-bc_obst_h)
      yymin=Rp(0)*sin((-0.5*jmax*px)*dphi)-0.22*(depth-bc_obst_h) !boxsize
      yymax=Rp(0)*sin((0.5*px*jmax)*dphi)+0.22*(depth-bc_obst_h) !boxsize	    


!       lmSEM2old=lmSEM2
      lmxSEM2old=lmxSEM2
      lmySEM2old=lmySEM2
      lmzSEM2old=lmzSEM2
      zSEM2old=zSEM2
      ySEM2old=ySEM2
      do i=1,nmax2
	xSEM2(i)=xSEM2(i)+uSEM2(i)*dt
	if (xSEM2(i)>xxmax)  then  ! just put 2*boxsizes* back in x-direction, keep y,z the same to maintain linked list, randomize eps
	  move=1
	  if (rank.eq.1) then
	    iimax=iimax+1
	    ind(iimax)=i
	    xSEM2(i)=xxmin
	    call random_number(yy)
	    ySEM2(i)=(yymax-yymin)*yy+yymin
	    call random_number(zz)
	    zSEM2(i)=(zzmax-zzmin)*zz+zzmin
	    call random_number(eps)
	    eps=anint(eps)
	    do j=1,3
	      epsSEM2(i,j)=eps(j)*2.-1. ! +1 or -1	
	    enddo	
!	      uSEM2(i)=ust/kappa*log(zSEM2(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)  !correct sign with and without TSHD)
	      uSEM2(i)=ust/kappa*log(zSEM2(i)/z0)*(-signU_bSEM*cos(phi)-signV_bSEM*sin(phi))*LOA/ABS(LOA)+U_TSHD !correct sign with and without TSHD)
	      lmxSEM2(i)=0.5*(depth-bc_obst_h)/sqrt(U_bSEM**2+V_bSEM**2)*ust/kappa*log(zSEM2(i)/z0)
	      lmxSEM2(i)=max(lm_min,lmxSEM2(i))
	      lmySEM2(i)=0.22*(depth-bc_obst_h)
	      lmzSEM2(i)=fac*zSEM2(i)/(depth-bc_obst_h)*sqrt(1.-zSEM2(i)/(depth-bc_obst_h))
	      lmzSEM2(i)=max(lm_min,lmzSEM2(i))
	  endif
	endif
      enddo	 
      if (move.eq.1) then
	call mpi_bcast(iimax,1,MPI_INTEGER,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(ind(1:iimax),iimax,MPI_INTEGER,1,MPI_COMM_WORLD,ierr)
	! remove old location moved items from linked list on every node for SEM2:
	do iii=1,iimax
	  ii=ind(iii)
	  kminSEM=INT(FLOOR((zSEM2old(ii)-lmzSEM2old(ii))/dz))+kbed_bc
	  kmaxSEM=INT(CEILING((zSEM2old(ii)+lmzSEM2old(ii))/dz))+kbed_bc
	  kminSEM=MAX(kminSEM,1)
	  kminSEM=MIN(kminSEM,kmax)
	  kmaxSEM=MAX(kmaxSEM,1)
	  kmaxSEM=MIN(kmaxSEM,kmax)

! 	  dphiSEM=1.1*lmSEM2old(ii)/Ru(0)
! 	  phiSEM=atan(ySEM2old(ii)/(xSEM2old(ii)+schuif_x))
! 	  jminSEM=INT(FLOOR((phiSEM-dphiSEM-(rank*jmax-0.5*jmax*px-0.5)*dphi)/dphi))
! 	  jminSEM=MAX(jminSEM,1)
! 	  jminSEM=MIN(jminSEM,j1)
! 	  jmaxSEM=INT(CEILING((phiSEM+dphiSEM-(rank*jmax-0.5*jmax*px-0.5)*dphi)/dphi))
! 	  jmaxSEM=MAX(jmaxSEM,1)
! 	  jmaxSEM=MIN(jmaxSEM,j1)

!           if (iii.eq.1) then
! 	    write(*,*)'****************************************'
! 	    write(*,*)' REMOVE rank,ii',rank,ii
! 	    write(*,*)'zSEM2(ii),xSEM2(ii),lmSEM2(ii)',zSEM2(ii),xSEM2(ii),lmSEM2(ii)
! 	    write(*,*)'phiSEM,dphiSEM',phiSEM,dphiSEM
! 	    write(*,*)'rank,jminSEM,jmaxSEM,kminSEM,kmaxSEM',rank,jminSEM,jmaxSEM,kminSEM,kmaxSEM
! 	    write(*,*)'****************************************'
!           endif

	  do k=kminSEM,kmaxSEM
	    z=(k-kbed_bc)*dz-0.5*dz
	    do j=0,j1 !j=jminSEM,jmaxSEM
	      y=Ru(0)*sin_u(j)	
	      if(ABS(y-ySEM2old(ii))/lmySEM2old(ii).lt.1..and.ABS(z-zSEM2old(ii))/lmzSEM2old(ii).lt.1.) then
		do tel=1,llmax2(j,k)
		  if(llist2(j,k,tel).eq.ii) then
		    llist2(j,k,tel)=llist2(j,k,llmax2(j,k))
		    llmax2(j,k)=llmax2(j,k)-1
		  endif
		enddo
	      endif
	    enddo
	  enddo
	enddo

	call mpi_bcast(xSEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(ySEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(zSEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(uSEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(epsSEM2,3*nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(lmxSEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(lmySEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)
	call mpi_bcast(lmzSEM2,nmax2,MPI_REAL8,1,MPI_COMM_WORLD,ierr)

	! add new location moved items in linked list on every node for SEM2:
	do iii=1,iimax
	  ii=ind(iii)
	  kminSEM=INT(FLOOR((zSEM2(ii)-lmzSEM2(ii))/dz))+kbed_bc
	  kmaxSEM=INT(CEILING((zSEM2(ii)+lmzSEM2(ii))/dz))+kbed_bc
	  kminSEM=MAX(kminSEM,1)
	  kminSEM=MIN(kminSEM,kmax)
	  kmaxSEM=MAX(kmaxSEM,1)
	  kmaxSEM=MIN(kmaxSEM,kmax)

! 	  dphiSEM=1.1*lmSEM2(ii)/Ru(0)
! 	  phiSEM=atan(ySEM2(ii)/(xSEM2(ii)+schuif_x))
! 	  jminSEM=INT(FLOOR((phiSEM-dphiSEM-(rank*jmax-0.5*jmax*px-0.5)*dphi)/dphi))
! 	  jminSEM=MAX(jminSEM,1)
! 	  jminSEM=MIN(jminSEM,j1)
! 	  jmaxSEM=INT(CEILING((phiSEM+dphiSEM-(rank*jmax-0.5*jmax*px-0.5)*dphi)/dphi))
! 	  jmaxSEM=MAX(jmaxSEM,1)
! 	  jmaxSEM=MIN(jmaxSEM,j1)
	  do k=kminSEM,kmaxSEM
 	    z=(k-kbed_bc)*dz-0.5*dz
	    do j=0,j1 !jminSEM,jmaxSEM
	      y=Ru(0)*sin_u(j)	
	      if(ABS(y-ySEM2(ii))/lmySEM2(ii).lt.1..and.ABS(z-zSEM2(ii))/lmzSEM2(ii).lt.1.) then
		llmax2(j,k)=llmax2(j,k)+1
		llist2(j,k,llmax2(j,k))=ii
	      endif
	    enddo
	  enddo
	enddo
      endif

!       if (rank.eq.1) then
! 	write(*,*)'rank,ind(3),ii',rank,ind(3),ii
!       endif
!       if (rank.eq.3) then
! 	write(*,*)'rank,ind(3),ii',rank,ind(3),ii
!       endif
! 

      end

      subroutine Chol_tensor_from_DNS_Ret395(y,A)
!       subroutine Chol_tensor_from_DNS_Ret395(y,A)
! 	note in this subroutine y,z is switched compared with Dflow3D, in this subroutine y=vertical (and R22=Ryy=vertical Re stress)
! 	input:
! 	y	vertical position divided by depth (-)
! 	output Cholesky decomposition tensor A(i,j) from target Re stress tensor from DNS with Ret=395:
! 	u_i =		sum_j Aij*ufluc_j
! 	Aij		3x3 tensor constructed from Reynolds stress tensor from DNS with Ret=395
! 	DNS Data from:
! 	# Authors: Moser, Kim & Mansour
! 	# Reference: DNS of Turbulent Channel Flow up to Re_tau=590, 1999,
! 	#            Physics of Fluids, vol 11, 943-945.
! 	# Numerical Method: Kim, Moin & Moser, 1987, J. Fluid Mech. vol 177, 133-166 
! 	# Re_tau = 392.24
! 	# Normalization: U_tau, h 
! 	# Description: Components of the Reynolds stress tensor 


      implicit none
	
      real ydns(130),R11dns(129),R22dns(129),R33dns(129),R21dns(129),R31dns(129),R32dns(129)
      real y,R11,R22,R33,R21,R31,R32,A(3,3)
      real factor,ydns1
      integer i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      ydns(:)=(/0.,7.5298e-05,0.00030118,0.00067762,0.0012045,0.0018819,0.0027095,0.0036874,0.0048153,
     & 0.006093,0.0075205,0.0090974,0.010823,0.012699,0.014722,0.016895,0.019215,0.021683,0.024298,0.02706,
     & 0.029969,0.033024,0.036224,0.039569,0.04306,0.046694,0.050472,0.054393,0.058456,0.062661,0.067007,
     & 0.071494,0.07612,0.080886,0.08579,0.090832,0.096011,0.10133,0.10678,0.11236,0.11808,0.12393,0.12991,
     & 0.13603,0.14227,0.14864,0.15515,0.16178,0.16853,0.17541,0.18242,0.18954,0.19679,0.20416,0.21165,0.21926,
     & 0.22699,0.23483,0.24279,0.25086,0.25905,0.26735,0.27575,0.28427,0.29289,0.30162,0.31046,0.3194,0.32844,
     & 0.33758,0.34683,0.35617,0.36561,0.37514,0.38477,0.39449,0.4043,0.4142,0.42419,0.43427,0.44443,0.45467,
     & 0.465,0.47541,0.4859,0.49646,0.5071,0.51782,0.5286,0.53946,0.55039,0.56138,0.57244,0.58357,0.59476,0.60601,
     & 0.61732,0.62868,0.6401,0.65158,0.66311,0.67469,0.68632,0.69799,0.70972,0.72148,0.73329,0.74513,0.75702,0.76894,
     & 0.7809,0.79289,0.80491,0.81696,0.82904,0.84114,0.85327,0.86542,0.87759,0.88978,0.90198,0.9142,0.92644,0.93868,
     & 0.95093,0.96319,0.97546,0.98773,1.,9999./) ! last value added to keep loop from running infinitely
      R11dns(:)=(/5.3376e-28,0.00013684,0.0021828,0.01102,0.034718,0.084454,0.17438,0.32123,0.54336,0.85851,1.2802,1.8132,2.449,
     & 3.1644,3.9229,4.6811,5.3953,6.0286,6.5558,6.9642,7.2524,7.4274,7.5016,7.4897,7.4069,7.2679,7.0865,6.875,6.6441,6.4025,6.1573
     & ,5.9141,5.6767,5.4477,5.2291,5.0223,4.8279,4.6462,4.4768,4.3194,4.173,4.037,3.9106,3.793,3.6834,3.5821,3.4896,3.4053,3.3278,
     & 3.2557,3.1886,3.1269,3.0702,3.0178,2.9678,2.9191,2.8712,2.8239,2.7776,2.7324,2.6883,2.6453,2.6031,2.5612,2.5196,2.4793,
     & 2.4411,2.4055,2.3714,2.338,2.3046,2.2708,2.2365,2.2023,2.1682,2.1335,2.0983,2.0631,2.0277,1.992,1.9568,1.9229,1.8896,1.8562,
     & 1.8232,1.7908,1.7591,1.7274,1.6948,1.661,1.6264,1.591,1.555,1.5184,1.4808,1.4425,1.4037,1.3647,1.326,1.2881,1.2513,1.2162,
     & 1.1827,1.15,1.1175,1.0843,1.0507,1.0173,0.98478,0.95377,0.92454,0.89706,0.87113,0.84645,0.82296,0.80062,0.77929,0.75889,
     & 0.73932,0.72051,0.70261,0.68576,0.67023,0.6562,0.64389,0.63375,0.62626,0.62174,0.62024/)
      R22dns(:)=(/4.9018e-29,9.4734e-11,2.3603e-08,5.7878e-07,5.44e-06,3.0014e-05,0.00011756,0.00036204,0.00093203,0.0020879,
     & 0.0041852,0.007661,0.013005,0.020724,0.0313,0.045157,0.062636,0.083971,0.10928,0.13854,0.17162,0.20824,0.24803,0.2905,
     & 0.3351,0.38121,0.42818,0.47535,0.52211,0.5679,0.61229,0.65493,0.69554,0.73388,0.76971,0.80285,0.83316,0.86061,0.8852,0.907,
     & 0.92607,0.94247,0.95626,0.96754,0.97651,0.98341,0.9885,0.99197,0.99396,0.99455,0.99379,0.99183,0.98887,0.98512,0.98075,
     & 0.97583,0.97044,0.96465,0.95846,0.95179,0.94459,0.93679,0.92845,0.91976,0.91087,0.90185,0.89268,0.88335,0.87385,0.86418,
     & 0.85436,0.84439,0.83421,0.82369,0.81285,0.80181,0.79074,0.77968,0.76856,0.75731,0.74588,0.73425,0.72249,0.71064,0.69862,
     & 0.68642,0.67417,0.66201,0.65007,0.63843,0.62712,0.61613,0.60545,0.59512,0.58523,0.57589,0.56715,0.55896,0.5512,0.54372,
     & 0.53642,0.52919,0.52197,0.51473,0.5075,0.50029,0.49312,0.48605,0.47922,0.47268,0.46644,0.46047,0.45478,0.44932,0.4441,
     & 0.43915,0.43458,0.43044,0.42671,0.42339,0.42046,0.41796,0.41591,0.41431,0.41313,0.41234,0.41188,0.41166,0.4116/)     
      R33dns(:)=(/5.9575e-28,5.2943e-05,0.00082715,0.0040392,0.012147,0.027852,0.053564,0.090947,0.14062,0.20209,0.27384,0.35369,
     & 0.43919,0.52804,0.61834,0.70871,0.79819,0.8861,0.97178,1.0545,1.1333,1.2075,1.2762,1.3387,1.3948,1.4444,1.4878,1.5251,1.5568
     & ,1.5833,1.6047,1.6216,1.6344,1.6439,1.6512,1.6567,1.6606,1.6623,1.6618,1.6594,1.6553,1.65,1.6439,1.6373,1.6299,1.6209,1.6105
     & ,1.5991,1.5872,1.5754,1.5637,1.551,1.5367,1.5212,1.5052,1.4889,1.4721,1.4551,1.4379,1.4206,1.4035,1.3867,1.3696,1.3515,
     & 1.3322,1.3118,1.291,1.27,1.2488,1.2276,1.2066,1.1855,1.1641,1.1423,1.1202,1.0983,1.0767,1.0551,1.0339,1.0131,0.99279,0.97254
     & ,0.95187,0.93057,0.90936,0.88926,0.87042,0.85243,0.83469,0.8169,0.79939,0.78241,0.76587,0.74931,0.73258,0.71572,0.69868,
     & 0.68159,0.66492,0.64916,0.63424,0.62001,0.6066,0.59377,0.58101,0.56815,0.55541,0.54304,0.53113,0.51971,0.50882,0.49853,
     & 0.48883,0.47982,0.47164,0.46392,0.45618,0.44842,0.44071,0.43345,0.42708,0.42169,0.41709,0.41319,0.40984,0.40689,0.40431,
     & 0.40247,0.40179/)
      R21dns(:)=(/1.1206e-30,-2.4414e-08,-1.5728e-06,-1.8113e-05,-0.00010324,-0.00040027,-0.0012143,-0.0031009,-0.0069486,-0.014015
     & ,-0.025858,-0.044135,-0.070276,-0.10513,-0.14869,-0.19999,-0.25724,-0.31808,-0.38001,-0.44073,-0.49837,-0.55162,-0.59969,
     & -0.64224,-0.67924,-0.7109,-0.73759,-0.7598,-0.7781,-0.79302,-0.80503,-0.81452,-0.82184,-0.82733,-0.83133,-0.83412,-0.8359,
     & -0.83683,-0.837,-0.83649,-0.83523,-0.83321,-0.83044,-0.82687,-0.82252,-0.81762,-0.81249,-0.80726,-0.8018,-0.79592,-0.78961,
     & -0.78296,-0.77613,-0.76929,-0.76243,-0.75552,-0.74851,-0.74146,-0.73445,-0.72737,-0.72012,-0.71271,-0.70521,-0.69761,
     & -0.68985,-0.68192,-0.67396,-0.66598,-0.65794,-0.64977,-0.64148,-0.63298,-0.62415,-0.615,-0.60551,-0.59561,-0.5853,-0.5746,
     & -0.5635,-0.55204,-0.5404,-0.52883,-0.51747,-0.50626,-0.49507,-0.48386,-0.47269,-0.46164,-0.45072,-0.43999,-0.42949,-0.4191,
     & -0.40876,-0.39853,-0.38843,-0.37835,-0.3681,-0.35761,-0.34695,-0.33621,-0.32544,-0.31468,-0.30393,-0.29315,-0.28231,-0.27133
     & ,-0.26019,-0.24893,-0.23753,-0.22605,-0.21462,-0.20326,-0.19188,-0.18043,-0.16898,-0.15752,-0.14601,-0.13442,-0.12272,
     & -0.11094,-0.09912,-0.087301,-0.075441,-0.063452,-0.051252,-0.038772,-0.026001,-0.013031,0./)
      R31dns(:)=(/-6.1282e-30,5.9464e-08,8.3228e-07,3.3616e-06,1.0732e-05,2.8297e-05,6.6865e-05,0.00014506,0.00029097,0.00054723,
     & 0.00097373,0.0016529,0.0026847,0.0041721,0.0061889,0.0087562,0.011834,0.015387,0.019449,0.024133,0.029501,0.035418,0.041511,
     & 0.04733,0.052512,0.056926,0.060659,0.063835,0.066418,0.068207,0.069001,0.068896,0.068285,0.067481,0.06645,0.064937,0.062875,
     & 0.06071,0.059182,0.058527,0.058437,0.058512,0.058376,0.058084,0.058125,0.058537,0.058809,0.058494,0.057876,0.057034,0.055746
     & ,0.054259,0.052957,0.051819,0.050872,0.050514,0.050564,0.050229,0.049171,0.047565,0.045716,0.043916,0.042134,0.040077,
     & 0.03752,0.034358,0.030858,0.027094,0.023371,0.0205,0.018795,0.017891,0.017572,0.017517,0.017095,0.015913,0.014306,0.012812,
     & 0.011907,0.01164,0.011503,0.011097,0.010218,0.0090004,0.0078753,0.0071617,0.0070121,0.0073802,0.0082299,0.0093954,0.010585,
     & 0.011834,0.013728,0.016563,0.019835,0.022555,0.023965,0.024099,0.023378,0.022236,0.020871,0.019224,0.017273,0.015116,
     & 0.013055,0.01136,0.0099787,0.0087677,0.0078868,0.0076156,0.0079424,0.0086379,0.009584,0.01054,0.011159,0.011209,0.010617,
     & 0.009526,0.0082998,0.0072287,0.0064502,0.0061696,0.0064277,0.0070238,0.0077818,0.0087489,0.0098627,0.010731,0.011032/)
      R32dns(:)=(/6.2793e-32,-2.8525e-11,-1.7141e-09,-1.9051e-08,-1.0411e-07,-3.8272e-07,-1.0969e-06,-2.6288e-06,-5.4596e-06,
     & -1.0091e-05,-1.7119e-05,-2.7764e-05,-4.4587e-05,-7.1044e-05,-0.00010893,-0.0001546,-0.00019552,-0.00021029,-0.00017262,
     & -6.049e-05,0.00013268,0.00039086,0.00067818,0.00095511,0.0011989,0.0014016,0.0015532,0.001636,0.0016465,0.0016026,0.0015073
     & ,0.0013268,0.0010303,0.00064119,0.00023642,-0.00013612,-0.00047569,-0.00078412,-0.0010539,-0.0012655,-0.0014126,-0.001493
     & ,-0.0015118,-0.0015196,-0.0016228,-0.0018295,-0.0020145,-0.0020826,-0.0021101,-0.0021628,-0.0022441,-0.002322,-0.0023519,
     & -0.0023018,-0.0020851,-0.0016209,-0.00085802,0.00011783,0.0010337,0.001644,0.0019851,0.0022291,0.0024143,0.0025128,
     & 0.0025006,0.0023726,0.0021952,0.0020517,0.0019999,0.00197,0.0018189,0.0014336,0.00075375,-7.929e-05,-0.0008372,-0.0013413
     & ,-0.0015771,-0.0015179,-0.0011762,-0.00070498,-0.00020689,0.00020529,0.00050971,0.00076211,0.00092292,0.00093017,0.00072814
     & ,0.00027352,-0.00029303,-0.00074565,-0.00095384,-0.00095389,-0.0010075,-0.001297,-0.001762,-0.0022305,-0.0025637,-0.0027183
     & ,-0.0027095,-0.0025722,-0.0024646,-0.0025472,-0.0027871,-0.0030094,-0.0031202,-0.0031772,-0.003241,-0.0033171,-0.0033738
     & ,-0.0034135,-0.0034897,-0.003672,-0.0039882,-0.0043896,-0.0047685,-0.0050417,-0.0051687,-0.005104,-0.0048287,-0.0044056,
     & -0.0038984,-0.0033462,-0.0027927,-0.0022376,-0.0016855,-0.0011858,-0.00076554,-0.00038835,0./)

      i=1
      do while (y.gt.ydns(i))
	i=i+1
      enddo
	if (i>129) then
	  write(*,*)'Error y>1 in Re_stress_tensor',i,ydns(i),y
	  stop
	endif

      if (i>1) then
	factor=(y-ydns(i-1))/(ydns(i)-ydns(i-1))
	R11=(1.-factor)*R11dns(i-1)+factor*R11dns(i)	
	R22=(1.-factor)*R22dns(i-1)+factor*R22dns(i)	
	R33=(1.-factor)*R33dns(i-1)+factor*R33dns(i)	
	R21=(1.-factor)*R21dns(i-1)+factor*R21dns(i)	
	R31=(1.-factor)*R31dns(i-1)+factor*R31dns(i)
	R32=(1.-factor)*R32dns(i-1)+factor*R32dns(i)
      else
	R11=R11dns(1)
    	R22=R22dns(1)
	R33=R33dns(1)
	R21=R21dns(1)
	R31=R31dns(1)
	R32=R32dns(1)
      endif
      
      A=0.
      A(1,1)=		sqrt(R11)
      A(2,1)=		R21/A(1,1)
      A(3,1)=		R31/A(1,1)
      A(2,2)=		sqrt(R22-A(2,1)*A(2,1))
      A(3,2)=		(R32-A(2,1)*A(3,1))/A(2,2)
      A(3,3)=		sqrt(R33-A(3,1)**2-A(3,2)**2)

      end
