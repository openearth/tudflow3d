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


      subroutine stress_terms(Uvel,Vvel,Wvel,ib,ie,jb,je,kb,ke)

      USE nlist

      implicit none

      include 'mpif.h' 

c
c*****************************************************************
!	determine stress terms:

!	Srr, Spr, Szr, Spp, Spz, Szz

!	with 

!	r   : 	r-direction
!	p   :   phi-direction
!	z   :   z-direction
c
c*****************************************************************

      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
      real     Uvel(0:i1,0:j1,0:k1),Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         eppo,epmo,epop,epom
	real divv

	real aax_1ie(1:ie),bbx_1ie(1:ie),ccx_1ie(1:ie),ddx_1ie(1:ie)
	real aaz_1ke(1:ke),bbz_1ke(1:ke),ccz_1ke(1:ke),ddz_1ke(1:ke)
	real aay_1jepx(1:je*px),bby_1jepx(1:je*px),ccy_1jepx(1:je*px),ddy_1jepx(1:je*px)
	real aax_0ie(0:ie),bbx_0ie(0:ie),ccx_0ie(0:ie),ddx_0ie(0:ie)
        real aay_0jepx(0:je*px),bby_0jepx(0:je*px),ccy_0jepx(0:je*px),ddy_0jepx(0:je*px)
	real aaz_0ke(0:ke),bbz_0ke(0:ke),ccz_0ke(0:ke),ddz_0ke(0:ke)

	real dvdr(0:ie),dwdr(0:ie),dudrphi_T(0:je*px),dwdrphi_T(0:je*px),dvdz(0:ke),dudz(0:ke)
	real Spp_T(1:ie,1:je*px,1:ke/px),Spr_T(0:ie,0:je*px,1:ke/px),Spz_T(1:ie,0:je*px,1:ke/px)  !Srr_T(1:ie,1:je*px,1:ke/px)
	real Uvel_T(0:ie,0:je*px+1,1:ke/px),Vvel_T(1:ie,0:je*px+1,1:ke/px),Wvel_T(1:ie,0:je*px+1,1:ke/px)

        integer ileng,ierr,itag,status(MPI_STATUS_SIZE)

	!! first swap u,v,w -> u_T,v_T,w_T 
	call t2np_0ie(Uvel(0:ie,1:je,1:ke),Uvel_T(0:ie,1:je*px,1:ke/px))
	call t2np(Vvel(1:ie,1:je,1:ke),Vvel_T(1:ie,1:je*px,1:ke/px))
	call t2np(Wvel(1:ie,1:je,1:ke),Wvel_T(1:ie,1:je*px,1:ke/px))
	!! also pass over boundaries at j=0 :
	IF (rank.eq.0) THEN
	  do i=1,px-1
      		call mpi_send(Uvel(0:ie,0,i*ke/px+1:(i+1)*ke/px),i1*ke/px,MPI_REAL8,i,i+100,MPI_COMM_WORLD,status,ierr)
      		call mpi_send(Vvel(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+101,MPI_COMM_WORLD,status,ierr)
      		call mpi_send(Wvel(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+102,MPI_COMM_WORLD,status,ierr)
	  enddo
                Uvel_T(0:ie,0,1:ke/px)=Uvel(0:ie,0,1:ke/px)
                Vvel_T(1:ie,0,1:ke/px)=Vvel(1:ie,0,1:ke/px)
                Wvel_T(1:ie,0,1:ke/px)=Wvel(1:ie,0,1:ke/px)
	else
		call mpi_recv(Uvel_T(0:ie,0,1:ke/px),i1*ke/px,MPI_REAL8,0,rank+100,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(Vvel_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+101,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(Wvel_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+102,MPI_COMM_WORLD,status,ierr)
	endif
	!! also pass over boundaries at j=jmax+1 :
	IF (rank.eq.px-1) THEN
	  do i=0,px-2
      		call mpi_send(Uvel(0:ie,je+1,i*ke/px+1:(i+1)*ke/px),i1*ke/px,MPI_REAL8,i,i+300,MPI_COMM_WORLD,status,ierr)
      		call mpi_send(Vvel(1:ie,je+1,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+301,MPI_COMM_WORLD,status,ierr)
      		call mpi_send(Wvel(1:ie,je+1,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+302,MPI_COMM_WORLD,status,ierr)
	  enddo
                Uvel_T(0:ie,je*px+1,1:ke/px)=Uvel(0:ie,je+1,rank*ke/px+1:(rank+1)*ke/px)
                Vvel_T(1:ie,je*px+1,1:ke/px)=Vvel(1:ie,je+1,rank*ke/px+1:(rank+1)*ke/px)
                Wvel_T(1:ie,je*px+1,1:ke/px)=Wvel(1:ie,je+1,rank*ke/px+1:(rank+1)*ke/px)
	else
		call mpi_recv(Uvel_T(0:ie,je*px+1,1:ke/px),i1*ke/px,MPI_REAL8,px-1,rank+300,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(Vvel_T(1:ie,je*px+1,1:ke/px),ie*ke/px,MPI_REAL8,px-1,rank+301,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(Wvel_T(1:ie,je*px+1,1:ke/px),ie*ke/px,MPI_REAL8,px-1,rank+302,MPI_COMM_WORLD,status,ierr)
	endif

       	bbx_1ie=1.
	aax_1ie=1./22. 
      	aax_1ie(1)=0.
	aax_1ie(ie)=0. 
	ccx_1ie=aax_1ie
       	bbz_1ke=1.
	aaz_1ke=1./22. 
      	aaz_1ke(1)=0.
	aaz_1ke(ke)=0. 
	ccz_1ke=aaz_1ke
       	bby_1jepx=1.
	aay_1jepx=1./22. 
      	aay_1jepx(1)=0.
	aay_1jepx(je*px)=0. 
	ccy_1jepx=aay_1jepx
       	bbx_0ie=1.
	aax_0ie=1./22. 
      	aax_0ie(0)=0.
	aax_0ie(ie)=0. 
	ccx_0ie=aax_0ie
       	bby_0jepx=1.
	aay_0jepx=1./22. 
      	aay_0jepx(0)=0.
	aay_0jepx(je*px)=0. 
	ccy_0jepx=aay_0jepx
       	bbz_0ke=1.
	aaz_0ke=1./22. 
      	aaz_0ke(0)=0.
	aaz_0ke(ke)=0. 
	ccz_0ke=aaz_0ke

	do j=jb,je
	  do k=kb,ke
	    ddx_1ie(1)=(Uvel(1,j,k)-Uvel(0,j,k))/dr(1)
	    ddx_1ie(ie)=(Uvel(ie,j,k)-Uvel(ie-1,j,k))/dr(ie)
	    do i=2,ie-1
	      ddx_1ie(i)=12./11.*(Uvel(i,j,k)-Uvel(i-1,j,k))/dr(i)
   	    enddo
            CALL solve_tridiag(Srr(1:ie,j,k),aax_1ie,bbx_1ie,ccx_1ie,ddx_1ie,ie)
	  enddo
	enddo
	Srr=Srr*ekm(1:ie,1:je,1:ke)*2.

	do i=ib,ie
	  do j=jb,je
	    ddz_1ke(1)=(Wvel(i,j,1)-Wvel(i,j,0))/dz
	    ddz_1ke(ke)=(Wvel(i,j,ke)-Wvel(i,j,ke-1))/dz
	    do k=2,ke-1
	      ddz_1ke(k)=12./11.*(Wvel(i,j,k)-Wvel(i,j,k-1))/dz
   	    enddo
            CALL solve_tridiag(Szz(i,j,1:ke),aaz_1ke,bbz_1ke,ccz_1ke,ddz_1ke,ke)
	  enddo
	enddo
	Szz=Szz*ekm(1:ie,1:je,1:ke)*2.

	do i=ib,ie
	  do k=kb,ke/px
	    ddy_1jepx(1)=(Vvel_T(i,1,k)-Vvel_T(i,0,k))/(Rp(i)*dphi)
	    ddy_1jepx(je*px)=(Vvel_T(i,je*px,k)-Vvel_T(i,je*px-1,k))/(Rp(i)*dphi)
	    do j=2,je*px-1
	      ddy_1jepx(j)=12./11.*(Vvel_T(i,j,k)-Vvel_T(i,j-1,k))/(Rp(i)*dphi)
   	    enddo
            CALL solve_tridiag(Spp_T(i,1:je*px,k),aay_1jepx,bby_1jepx,ccy_1jepx,ddy_1jepx,je*px)
	  enddo
	enddo
	call t2fp(Spp_T(1:ie,1:je*px,1:ke/px),Spp(1:ie,1:je,1:ke))
	do j=jb,je
	  do k=kb,ke
	    Spp(1:ie,j,k)=ekm(1:ie,j,k)*2.*(Spp(1:ie,j,k)+(Uvel(1:ie,j,k)+Uvel(0:ie-1,j,k))/(2.*Rp(1:ie)))
	  enddo	
	enddo

      do i=ib,ie
	do k=kb,ke
	  do j=jb,je
      divv= ( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) ) / ( dz )   
      divv=2./3.*ekm(i,j,k)*divv
      Srr(i,j,k) = Srr(i,j,k) - divv
      Spp(i,j,k) = Spp(i,j,k) - divv
      Szz(i,j,k) = Szz(i,j,k) - divv
            enddo
         enddo
      enddo

	do j=0,je !include j=0
	  do k=kb,ke
	    ddx_0ie(0)=(Vvel(1,j,k)-Vvel(0,j,k))/(Rp(1)-Rp(0))
	    ddx_0ie(ie)=(Vvel(i1,j,k)-Vvel(ie,j,k))/(Rp(i1)-Rp(ie))
	    do i=1,ie-1
	      ddx_0ie(i)=12./11.*(Vvel(i+1,j,k)-Vvel(i,j,k))/(Rp(i+1)-Rp(i))
   	    enddo
            CALL solve_tridiag(dvdr,aax_0ie,bbx_0ie,ccx_0ie,ddx_0ie,i1)
	    Spr(0:ie,j,k)=dvdr
	  enddo
	enddo
	do j=0,je !include j=0
	  do i=ib,ie 
	    ddz_0ke(0)=(Vvel(i,j,1)-Vvel(i,j,0))/dz
	    ddz_0ke(ke)=(Vvel(i,j,k1)-Vvel(i,j,ke))/dz
	    do k=1,ke-1
	      ddz_0ke(k)=12./11.*(Vvel(i,j,k+1)-Vvel(i,j,k))/dz
   	    enddo
            CALL solve_tridiag(dvdz,aaz_0ke,bbz_0ke,ccz_0ke,ddz_0ke,k1)
	    Spz(i,j,0:ke)=dvdz
	  enddo
	enddo

	call t2np_0ie(Spr(0:ie,1:je,1:ke),Spr_T(0:ie,1:je*px,1:ke/px))	
	call t2np(Spz(1:ie,1:je,1:ke),Spz_T(1:ie,1:je*px,1:ke/px))
	IF (rank.eq.0) THEN !! also pass over boundaries at j=0 :
	  do i=1,px-1
      		call mpi_send(Spz(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+100,MPI_COMM_WORLD,status,ierr)
      		call mpi_send(Spr(0:ie,0,i*ke/px+1:(i+1)*ke/px),i1*ke/px,MPI_REAL8,i,i+200,MPI_COMM_WORLD,status,ierr)
	  enddo
                Spz_T(1:ie,0,1:ke/px)=Spz(1:ie,0,1:ke/px)
                Spr_T(0:ie,0,1:ke/px)=Spr(0:ie,0,1:ke/px)
	else
		call mpi_recv(Spz_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+100,MPI_COMM_WORLD,status,ierr)
		call mpi_recv(Spr_T(0:ie,0,1:ke/px),i1*ke/px,MPI_REAL8,0,rank+200,MPI_COMM_WORLD,status,ierr)
	endif


	do i=0,ie !include i=0
	  do k=kb,ke/px
	    ddy_0jepx(0)=(Uvel_T(i,1,k)-Uvel_T(i,0,k))/(Ru(i)*dphi)
	    ddy_0jepx(je*px)=(Uvel_T(i,je*px+1,k)-Uvel_T(i,je*px,k))/(Ru(i)*dphi)
	    do j=1,je*px-1
	      ddy_0jepx(j)=12./11.*(Uvel_T(i,j+1,k)-Uvel_T(i,j,k))/(Ru(i)*dphi)
   	    enddo
            CALL solve_tridiag(dudrphi_T,aay_0jepx,bby_0jepx,ccy_0jepx,ddy_0jepx,je*px+1)
	    Spr_T(i,0:je*px,k)=Spr_T(i,0:je*px,k)+dudrphi_T
	  enddo
	enddo	
	do i=ib,ie
	  do k=kb,ke/px !should include k=0,but is zero because Wvel(i,j,0)=0, so left out
	    ddy_0jepx(0)=(Wvel_T(i,1,k)-Wvel_T(i,0,k))/(Rp(i)*dphi)
	    ddy_0jepx(je*px)=(Wvel_T(i,je*px+1,k)-Wvel_T(i,je*px,k))/(Rp(i)*dphi)
	    do j=1,je*px-1
	      ddy_0jepx(j)=12./11.*(Wvel_T(i,j+1,k)-Wvel_T(i,j,k))/(Rp(i)*dphi)
   	    enddo
            CALL solve_tridiag(dwdrphi_T,aay_0jepx,bby_0jepx,ccy_0jepx,ddy_0jepx,je*px+1)
	    Spz_T(i,0:je*px,k)=Spz_T(i,0:je*px,k)+dwdrphi_T
	  enddo
	enddo

	call t2fp_0ie(Spr_T(0:ie,1:je*px,1:ke/px),Spr(0:ie,1:je,1:ke))		
	call t2fp(Spz_T(1:ie,1:je*px,1:ke/px),Spz(1:ie,1:je,1:ke))	
	IF (rank.eq.0) THEN !! also pass over boundaries at j=0 :
 	  do i=1,px-1
      		call mpi_recv(Spz(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+100,MPI_COMM_WORLD,status,ierr)
      		call mpi_recv(Spr(0:ie,0,i*ke/px+1:(i+1)*ke/px),i1*ke/px,MPI_REAL8,i,i+200,MPI_COMM_WORLD,status,ierr)
	  enddo
                Spz(1:ie,0,1:ke/px)=Spz_T(1:ie,0,1:ke/px)
                Spr(0:ie,0,1:ke/px)=Spr_T(0:ie,0,1:ke/px)
	else
		call mpi_send(Spz_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+100,MPI_COMM_WORLD,status,ierr)
		call mpi_send(Spr_T(0:ie,0,1:ke/px),i1*ke/px,MPI_REAL8,0,rank+200,MPI_COMM_WORLD,status,ierr)
	endif
	do j=0,je
	  do k=kb,ke
	    Spr(0:ie,j,k)=(Spr(0:ie,j,k)-(Vvel(0:ie,j,k)+Vvel(1:i1,j,k))/(2.*Ru(0:ie)))*
     &      0.25*(ekm(0:ie,j,k)+ekm(1:i1,j,k)+ekm(0:ie,j+1,k)+ekm(1:i1,j+1,k))
	  enddo
	enddo
	Spz(1:ie,0:je,0:ke)=Spz(1:ie,0:je,0:ke)*
     &  0.25*(ekm(1:ie,0:je,0:ke)+ekm(1:ie,1:j1,0:ke)+ekm(1:ie,0:je,1:k1)+ekm(1:ie,1:j1,1:k1))

	do j=1,je 
	  do i=0,ie !include i=0 
	    ddz_0ke(0)=(Uvel(i,j,1)-Uvel(i,j,0))/dz
	    ddz_0ke(ke)=(Uvel(i,j,k1)-Uvel(i,j,ke))/dz
	    do k=1,ke-1
	      ddz_0ke(k)=12./11.*(Uvel(i,j,k+1)-Uvel(i,j,k))/dz
   	    enddo
            CALL solve_tridiag(dudz,aaz_0ke,bbz_0ke,ccz_0ke,ddz_0ke,k1)
	    Szr(i,j,0:ke)=dudz
	  enddo
	enddo
	do j=jb,je 
	  do k=0,ke !include k=0
	    ddx_0ie(0)=(Wvel(1,j,k)-Wvel(0,j,k))/(Rp(1)-Rp(0))
	    ddx_0ie(ie)=(Wvel(i1,j,k)-Wvel(ie,j,k))/(Rp(i1)-Rp(ie))
	    do i=1,ie-1
	      ddx_0ie(i)=12./11.*(Wvel(i+1,j,k)-Wvel(i,j,k))/(Rp(i+1)-Rp(i))
   	    enddo
            CALL solve_tridiag(dwdr,aax_0ie,bbx_0ie,ccx_0ie,ddx_0ie,i1)
	    Szr(0:ie,j,k)=Szr(0:ie,j,k)+dwdr
	  enddo
	enddo
	Szr=Szr*0.25*(ekm(0:ie,1:je,0:ke)+ekm(1:i1,1:je,0:ke)+ekm(0:ie,1:je,1:k1)+ekm(1:i1,1:je,1:k1))

	end

      subroutine diffu_com4(putout,Sigrr,Sigpr,Sigzr,Sigpp,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)

      implicit none
      include 'mpif.h'
c
c*****************************************************************
c
c      diffu calculates the diffusion of u-velocity, which is
c      the velocity in the radial direction.
c
c
c     In formula:  (4 terms)
c
c    1  d                  1  d                     d
c    - -- (r Sigma(r r)) + - ---- (Sigma(phi r)) + -- (Sigma(z r)) -
c    r dr                  r dphi                  dz
c
c
c    1
c    - Sigma(phi phi)
c    r
c
c        r   : direction  ---> explicit (subroutine diffu)
c        phi : direction  ---> implicit (subroutine predic)
c        z   : direction  ---> explicit (subroutine diffu)
c
c      on input :
c
c          putout            : advection part
c          Uvel,Vvel,Wvel    : contain velocities at n-1
c          ekm               : diffusion coefficients (for velocity) in
c                              center points of the grid cells
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              diffusion part has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection and diffusion part
c          other parameters  : all unchanged
c
c*****************************************************************
      integer  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px
      real     putout(0:i1,0:j1,0:k1),Sigrr(ib:ie,jb:je,kb:ke),
     +         Sigpr(0:ie,0:je,kb:ke),Sigzr(0:ie,1:je,0:ke),
     +         Sigpp(ib:ie,jb:je,kb:ke),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
      	integer i,j,k
        integer ileng,ierr,itag,status(MPI_STATUS_SIZE)

      	real dSdr(0:ie),dSdphi(1:je*px),dSdz(1:ke),putout_T(1:ie,1:je*px,1:ke/px)
      	real Sigpr_T(1:ie,0:je*px,1:ke/px)
        real aax_0ie(0:ie),bbx_0ie(0:ie),ccx_0ie(0:ie),ddx_0ie(0:ie)
	real aay_1jepx(1:je*px),bby_1jepx(1:je*px),ccy_1jepx(1:je*px),ddy_1jepx(1:je*px)
	real aaz_1ke(1:ke),bbz_1ke(1:ke),ccz_1ke(1:ke),ddz_1ke(1:ke)


       	bbx_0ie=1.
	aax_0ie=1./22. 
      	aax_0ie(0)=0.
	aax_0ie(ie)=23.
	ccx_0ie=1./22. 
      	ccx_0ie(0)=23.
	ccx_0ie(ie)=0.
       	bby_1jepx=1.
	aay_1jepx=1./22. 
      	aay_1jepx(1)=0.
	aay_1jepx(je*px)=0. !-1. 
	ccy_1jepx=1./22. 
      	ccy_1jepx(1)=0. !-1.
	ccy_1jepx(je*px)=0.
       	bbz_1ke=1.
	aaz_1ke=1./22. 
      	aaz_1ke(1)=0.
	aaz_1ke(ke)=0. !-1.
	ccz_1ke=1./22. 
      	ccz_1ke(1)=0. !-1.
	ccz_1ke(ke)=0.

	do j=jb,je
	  do k=kb,ke
	    Sigrr(1:ie,j,k)=Sigrr(1:ie,j,k)*Rp(1:ie)
	    ddx_0ie(0)=(-25.*Sigrr(1,j,k)+26.*Sigrr(2,j,k)-Sigrr(3,j,k))/(Ru(0)*(Rp(1)-Rp(0)))
	    ddx_0ie(ie)=(25.*Sigrr(ie,j,k)-26.*Sigrr(ie-1,j,k)+Sigrr(ie-2,j,k))/(Ru(ie)*(Rp(ie+1)-Rp(ie)))
	    do i=1,ie-1
	      ddx_0ie(i)=12./11.*(Sigrr(i+1,j,k)-Sigrr(i,j,k))/(Ru(i)*(Rp(i+1)-Rp(i)))
   	    enddo
            CALL solve_tridiag(dSdr(0:ie),aax_0ie,bbx_0ie,ccx_0ie,ddx_0ie,ie+1)
	    putout(1:ie,j,k)=putout(1:ie,j,k)+dSdr(1:ie)-Sigpp(1:ie,j,k)/Rp(1:ie)
	  enddo
	enddo    
	!! swap Spr -> Spr_T
	call t2np(Sigpr(1:ie,1:je,1:ke),Sigpr_T(1:ie,1:je*px,1:ke/px))
	call t2np(putout(1:ie,1:je,1:ke),putout_T(1:ie,1:je*px,1:ke/px))
	!! also pass over boundaries at j=0 :
	IF (rank.eq.0) THEN
	  do i=1,px-1
      		call mpi_send(Sigpr(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+100,MPI_COMM_WORLD,status,ierr)
	  enddo
                Sigpr_T(1:ie,0,1:ke/px)=Sigpr(1:ie,0,1:ke/px)
	else
		call mpi_recv(Sigpr_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+100,MPI_COMM_WORLD,status,ierr)
	endif

	do i=ib,ie
	  do k=kb,ke/px 
!	    ddy_1jepx(1)=(-Sigpr_T(i,0,k)+2.*Sigpr_T(i,1,k)-Sigpr_T(i,2,k))/(Ru(i)*dphi)
!	    ddy_1jepx(je*px)=(Sigpr_T(i,je*px,k)-2.*Sigpr_T(i,je*px-1,k)+Sigpr_T(i,je*px-2,k))/(Ru(i)*dphi)
	    ddy_1jepx(1)=(Sigpr_T(i,1,k)-Sigpr_T(i,0,k))/(Ru(i)*dphi)
	    ddy_1jepx(je*px)=(Sigpr_T(i,je*px,k)-Sigpr_T(i,je*px-1,k))/(Ru(i)*dphi)
	    do j=2,je*px-1
	      ddy_1jepx(j)=12./11.*(Sigpr_T(i,j,k)-Sigpr_T(i,j-1,k))/(Ru(i)*dphi)
   	    enddo
            CALL solve_tridiag(dSdphi(1:je*px),aay_1jepx,bby_1jepx,ccy_1jepx,ddy_1jepx,je*px)
	    putout_T(i,1:je*px,k)=putout_T(i,1:je*px,k)+dSdphi(1:je*px)
	  enddo
	enddo
	call t2fp(putout_T(1:ie,1:je*px,1:ke/px),putout(1:ie,1:je,1:ke))	

	do i=ib,ie
	  do j=jb,je
!	    ddz_1ke(1)=(-Sigzr(i,j,0)+2.*Sigzr(i,j,1)-Sigzr(i,j,2))/dz
!	    ddz_1ke(ke)=(Sigzr(i,j,ke)-2.*Sigzr(i,j,ke-1)+Sigzr(i,j,ke-2))/dz
	    ddz_1ke(1)=(Sigzr(i,j,1)-Sigzr(i,j,0))/dz
	    ddz_1ke(ke)=(Sigzr(i,j,ke)-Sigzr(i,j,ke-1))/dz
	    do k=2,ke-1
	      ddz_1ke(k)=12./11.*(Sigzr(i,j,k)-Sigzr(i,j,k-1))/dz
   	    enddo
            CALL solve_tridiag(dSdz(1:ke),aaz_1ke,bbz_1ke,ccz_1ke,ddz_1ke,ke)
	    putout(i,j,1:ke)=putout(i,j,1:ke)+dSdz(1:ke)
	  enddo
	enddo

      return
      end

      subroutine diffv_com4(putout,Sigpr,Sigpp,Sigpz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)

      implicit none

	include 'mpif.h'
c
c*****************************************************************
c
c      diffv calculates the diffusion of v-velocity, which is
c      the velocity in the tangential direction.
c
c
c     In formula:  (4 terms)
c
c  1  d                    1  d                       d
c  - -- (r Sigma(r phi)) + - ---- (Sigma(phi phi)) + -- (Sigma(z phi)) +
c  r dr                    r dphi                    dz
c
c
c    1
c    - Sigma(r phi)
c    r
c
c        r   : direction  ---> explicit (subroutine diffv)
c        phi : direction  ---> implicit (subroutine predic)
c        z   : direction  ---> explicit (subroutine diffv)
c
c      on input :
c
c          putout            : advection part
c          Uvel,Vvel,Wvel    : contain velocities at n-1
c          ekm               : diffusion coefficients (for velocity) in
c                              center points of the grid cells
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              diffusion part has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection and diffusion part
c          other parameters  : all unchanged
c
c*****************************************************************

      integer  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px
      real     putout(0:i1,0:j1,0:k1),Sigpp(ib:ie,jb:je,kb:ke),
     +         Sigpr(0:ie,0:je,kb:ke),Sigpz(1:ie,0:je,0:ke),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
      	integer i,j,k
        integer ileng,ierr,itag,status(MPI_STATUS_SIZE)

	real rSpr(0:ie),dSdr(1:ie),dSdphi(0:je*px),dSdz(1:ke)
      	real Sigpp_T(1:ie,1:je*px,1:ke/px),putout_T(1:ie,1:je*px,1:ke/px)
	real aax_1ie(1:ie),bbx_1ie(1:ie),ccx_1ie(1:ie),ddx_1ie(1:ie)
        real aay_0jepx(0:je*px),bby_0jepx(0:je*px),ccy_0jepx(0:je*px),ddy_0jepx(0:je*px)
	real aaz_1ke(1:ke),bbz_1ke(1:ke),ccz_1ke(1:ke),ddz_1ke(1:ke)

       	bbx_1ie=1.
	aax_1ie=1./22. 
      	aax_1ie(1)=0.
	aax_1ie(ie)=0. !-1.
	ccx_1ie=1./22. 
      	ccx_1ie(1)=0. !-1.
	ccx_1ie(ie)=0.
       	bby_0jepx=1.
	aay_0jepx=1./22. 
      	aay_0jepx(0)=0.
	aay_0jepx(je*px)=23. 
	ccy_0jepx=1./22. 
      	ccy_0jepx(0)=23.
	ccy_0jepx(je*px)=0.
       	bbz_1ke=1.
	aaz_1ke=1./22. 
      	aaz_1ke(1)=0.
	aaz_1ke(ke)=0. !-1.
	ccz_1ke=1./22. 
      	ccz_1ke(1)=0. !-1.
	ccz_1ke(ke)=0.

	do j=jb,je
	  do k=kb,ke
	    rSpr(0:ie)=Sigpr(0:ie,j,k)*Ru(0:ie)
!	    ddx_1ie(1)=(-rSpr(0)+2.*rSpr(1)-rSpr(2))/(Rp(1)*dr(1))
!	    ddx_1ie(ie)=(rSpr(ie)-2.*rSpr(ie-1)+rSpr(ie-2))/(Rp(ie)*dr(ie))
	    ddx_1ie(1)=(rSpr(1)-rSpr(0))/(Rp(1)*dr(1))
	    ddx_1ie(ie)=(rSpr(ie)-rSpr(ie-1))/(Rp(ie)*dr(ie))
	    do i=2,ie-1
	      ddx_1ie(i)=12./11.*(rSpr(i)-rSpr(i-1))/(Rp(i)*dr(i))
   	    enddo
            CALL solve_tridiag(dSdr(1:ie),aax_1ie,bbx_1ie,ccx_1ie,ddx_1ie,ie)
	    putout(1:ie,j,k)=putout(1:ie,j,k)+dSdr(1:ie)+(Sigpr(1:ie,j,k)+Sigpr(0:ie-1,j,k))/(2.*Rp(1:ie)) 
	  enddo
	enddo 
	!! swap Spp -> Spp_T
	call t2np(Sigpp(1:ie,1:je,1:ke),Sigpp_T(1:ie,1:je*px,1:ke/px))
	call t2np(putout(1:ie,1:je,1:ke),putout_T(1:ie,1:je*px,1:ke/px))

	do i=ib,ie
	  do k=kb,ke/px 
	    ddy_0jepx(0)=(-25.*Sigpp_T(i,1,k)+26.*Sigpp_T(i,2,k)-Sigpp_T(i,3,k))/(Rp(i)*dphi)
	    ddy_0jepx(je*px)=(25.*Sigpp_T(i,je*px,k)-26.*Sigpp_T(i,je*px-1,k)+Sigpp_T(i,je*px-2,k))/(Rp(i)*dphi)
	    do j=1,je*px-1
	      ddy_0jepx(j)=12./11.*(Sigpp_T(i,j+1,k)-Sigpp_T(i,j,k))/(Rp(i)*dphi)
   	    enddo
            CALL solve_tridiag(dSdphi(0:je*px),aay_0jepx,bby_0jepx,ccy_0jepx,ddy_0jepx,je*px+1)
	    putout_T(i,1:je*px,k)=putout_T(i,1:je*px,k)+dSdphi(1:je*px)
	  enddo
	enddo
	call t2fp(putout_T(1:ie,1:je*px,1:ke/px),putout(1:ie,1:je,1:ke))	
	do i=ib,ie
	  do j=jb,je
!	    ddz_1ke(1)=(-Sigpz(i,j,0)+2.*Sigpz(i,j,1)-Sigpz(i,j,2))/dz
!	    ddz_1ke(ke)=(Sigpz(i,j,ke)-2.*Sigpz(i,j,ke-1)+Sigpz(i,j,ke-2))/dz
	    ddz_1ke(1)=(Sigpz(i,j,1)-Sigpz(i,j,0))/dz
	    ddz_1ke(ke)=(Sigpz(i,j,ke)-Sigpz(i,j,ke-1))/dz
	    do k=2,ke-1
	      ddz_1ke(k)=12./11.*(Sigpz(i,j,k)-Sigpz(i,j,k-1))/dz
   	    enddo
            CALL solve_tridiag(dSdz(1:ke),aaz_1ke,bbz_1ke,ccz_1ke,ddz_1ke,ke)
	    putout(i,j,1:ke)=putout(i,j,1:ke)+dSdz(1:ke)
	  enddo
	enddo

      return
      end

      subroutine diffw_com4(putout,Sigzr,Sigpz,Sigzz,Ru,Rp,dr,dphi,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)

      implicit none

	include 'mpif.h'

      integer  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px
      real     putout(0:i1,0:j1,0:k1),Sigzz(ib:ie,jb:je,kb:ke),
     +         Sigzr(0:ie,1:je,0:ke),Sigpz(1:ie,0:je,0:ke),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
      	integer i,j,k
        integer ileng,ierr,itag,status(MPI_STATUS_SIZE)

	real dSdr(1:ie),dSdphi(1:je*px),dSdz(0:ke),putout_T(1:ie,1:je*px,1:ke/px)
	real Sigpz_T(1:ie,0:je*px,1:ke/px),rSzr(0:ie)
	real aax_1ie(1:ie),bbx_1ie(1:ie),ccx_1ie(1:ie),ddx_1ie(1:ie)
	real aay_1jepx(1:je*px),bby_1jepx(1:je*px),ccy_1jepx(1:je*px),ddy_1jepx(1:je*px)
        real aaz_0ke(0:ke),bbz_0ke(0:ke),ccz_0ke(0:ke),ddz_0ke(0:ke)

       	bbx_1ie=1.
	aax_1ie=1./22. 
      	aax_1ie(1)=0.
	aax_1ie(ie)=0. !-1.
	ccx_1ie=1./22. 
      	ccx_1ie(1)=0. !-1.
	ccx_1ie(ie)=0.
       	bby_1jepx=1.
	aay_1jepx=1./22. 
      	aay_1jepx(1)=0.
	aay_1jepx(je*px)=0. !-1. 
	ccy_1jepx=1./22. 
      	ccy_1jepx(1)=0. !-1.
	ccy_1jepx(je*px)=0.
       	bbz_0ke=1.
	aaz_0ke=1./22. 
      	aaz_0ke(0)=0.
	aaz_0ke(ke)=23.
	ccz_0ke=1./22. 
      	ccz_0ke(0)=23.
	ccz_0ke(ke)=0.

	do j=jb,je
	  do k=kb,ke
	    rSzr(0:ie)=Sigzr(0:ie,j,k)*Ru(0:ie)
!	    ddx_1ie(1)=(-rSzr(0)+2.*rSzr(1)-rSzr(2))/(Rp(1)*dr(1))
!	    ddx_1ie(ie)=(rSzr(ie)-2.*rSzr(ie-1)+rSzr(ie-2))/(Rp(ie)*dr(ie))
	    ddx_1ie(1)=(rSzr(1)-rSzr(0))/(Rp(1)*dr(1))
	    ddx_1ie(ie)=(rSzr(ie)-rSzr(ie-1))/(Rp(ie)*dr(ie))
	    do i=2,ie-1
	      ddx_1ie(i)=12./11.*(rSzr(i)-rSzr(i-1))/(Rp(i)*dr(i))
   	    enddo
            CALL solve_tridiag(dSdr(1:ie),aax_1ie,bbx_1ie,ccx_1ie,ddx_1ie,ie)
	    putout(1:ie,j,k)=putout(1:ie,j,k)+dSdr(1:ie)
	  enddo
	enddo 
	!! swap Spz -> Spz_T
	call t2np(Sigpz(1:ie,1:je,1:ke),Sigpz_T(1:ie,1:je*px,1:ke/px))
	call t2np(putout(1:ie,1:je,1:ke),putout_T(1:ie,1:je*px,1:ke/px))
	!! also pass over boundaries at j=0 :
	IF (rank.eq.0) THEN
	  do i=1,px-1
      		call mpi_send(Sigpz(1:ie,0,i*ke/px+1:(i+1)*ke/px),ie*ke/px,MPI_REAL8,i,i+100,MPI_COMM_WORLD,status,ierr)
	  enddo
                Sigpz_T(1:ie,0,1:ke/px)=Sigpz(1:ie,0,1:ke/px)
	else
		call mpi_recv(Sigpz_T(1:ie,0,1:ke/px),ie*ke/px,MPI_REAL8,0,rank+100,MPI_COMM_WORLD,status,ierr)
	endif
	do i=ib,ie
	  do k=kb,ke/px 
!	    ddy_1jepx(1)=(-Sigpz_T(i,0,k)+2.*Sigpz_T(i,1,k)-Sigpz_T(i,2,k))/(Rp(i)*dphi)
!	    ddy_1jepx(je*px)=(Sigpz_T(i,je*px,k)-2.*Sigpz_T(i,je*px-1,k)+Sigpz_T(i,je*px-2,k))/(Rp(i)*dphi)
	    ddy_1jepx(1)=(Sigpz_T(i,1,k)-Sigpz_T(i,0,k))/(Rp(i)*dphi)
	    ddy_1jepx(je*px)=(Sigpz_T(i,je*px,k)-Sigpz_T(i,je*px-1,k))/(Rp(i)*dphi)
	    do j=2,je*px-1
	      ddy_1jepx(j)=12./11.*(Sigpz_T(i,j,k)-Sigpz_T(i,j-1,k))/(Rp(i)*dphi)
   	    enddo
            CALL solve_tridiag(dSdphi(1:je*px),aay_1jepx,bby_1jepx,ccy_1jepx,ddy_1jepx,je*px)
	    putout_T(i,1:je*px,k)=putout_T(i,1:je*px,k)+dSdphi(1:je*px)
	  enddo
	enddo
	call t2fp(putout_T(1:ie,1:je*px,1:ke/px),putout(1:ie,1:je,1:ke))
	do i=ib,ie
	  do j=jb,je
	    ddz_0ke(0)=(-25.*Sigzz(i,j,1)+26.*Sigzz(i,j,2)-Sigzz(i,j,3))/dz
	    ddz_0ke(ke)=(25.*Sigzz(i,j,ke)-26.*Sigzz(i,j,ke-1)+Sigzz(i,j,ke-2))/dz
	    do k=1,ke-1
	      ddz_0ke(k)=12./11.*(Sigzz(i,j,k+1)-Sigzz(i,j,k))/dz
   	    enddo
            CALL solve_tridiag(dSdz(0:ke),aaz_0ke,bbz_0ke,ccz_0ke,ddz_0ke,ke+1)
	    putout(i,j,1:ke)=putout(i,j,1:ke)+dSdz(1:ke)
	  enddo
	enddo

      return
      end


      subroutine diffu_CDS2(putout,Uvel,Vvel,Wvel,ib,ie,jb,je,kb,ke)

      USE nlist

      implicit none
      
!       include 'param.txt'
!       include 'common.txt'


c
c*****************************************************************
c
c      diffu calculates the diffusion of u-velocity, which is
c      the velocity in the radial direction.
c
c
c     In formula:  (4 terms)
c
c    1  d                  1  d                     d
c    - -- (r Sigma(r r)) + - ---- (Sigma(phi r)) + -- (Sigma(z r)) -
c    r dr                  r dphi                  dz
c
c
c    1
c    - Sigma(phi phi)
c    r
c
c        r   : direction  ---> explicit (subroutine diffu)
c        phi : direction  ---> implicit (subroutine predic)
c        z   : direction  ---> explicit (subroutine diffu)
c
c      on input :
c
c          putout            : advection part
c          Uvel,Vvel,Wvel    : contain velocities at n-1
c          ekm               : diffusion coefficients (for velocity) in
c                              center points of the grid cells
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              diffusion part has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection and diffusion part
c          other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         eppo,epmo,epop,epom,drp,dzi,divergentie
      real theta_U,theta_V,xx,yy,f,fluc,Wjet
      real divvR,divvL,CNz
      integer n,t

	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSE
	  CNz=1.
	ENDIF
      dzi =1./dz
      do i=ib,ie
      ip=i+1
      im=i-1
      drp = Rp(ip)-Rp(i)
	do k=kb,ke
	kp=k+1
	km=k-1
	  do j=jb,je
	  jp=j+1
	  jm=j-1
c
      eppo = 0.25 * (
     +   ekm(i,j,k) + ekm(ip,j,k) + ekm(ip,jp,k) + ekm(i,jp,k)  )
      epmo = 0.25 * (
     +   ekm(i,j,k) + ekm(ip,j,k) + ekm(ip,jm,k) + ekm(i,jm,k)  )
      epop = 0.25 * (
     +   ekm(i,j,k) + ekm(ip,j,k) + ekm(ip,j,kp) + ekm(i,j,kp)  )
      epom = 0.25 * (
     +   ekm(i,j,k) + ekm(ip,j,k) + ekm(ip,j,km) + ekm(i,j,km)  )

      divvL= ( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) ) * dzi  
      divvR= ( Ru(ip)*Uvel(ip,j,k) - Ru(ip-1)*Uvel(ip-1,j,k) ) / ( Rp(ip)*dr(ip) )
     +              +
     2  (       Vvel(ip,j,k) -         Vvel(ip,j-1,k) ) / ( Rp(ip)*dphi )
     +              +
     3  (       Wvel(ip,j,k) -         Wvel(ip,j,k-1) ) * dzi   

      putout(i,j,k) = putout(i,j,k) +
     1 (  Rp(ip) * ekm(ip,j,k) *
     1            ((Uvel(ip,j,k) - Uvel(i,j,k) ) / ( dr(ip) ) - 1./3.*divvR) -
     1    Rp(i ) * ekm(i,j,k)  *
     1            ((Uvel(i,j,k)  - Uvel(im,j,k)) / ( dr(i)  ) - 1./3.*divvL) )  /
     1 ( 0.5 * Ru(i) * ( drp ) )
     +              +
     2 ( eppo * ( Ru(i) * ( Vvel(ip,j,k)/Rp(ip)  - Vvel(i,j,k)/Rp(i) ) /
     2                    ( Rp(ip) - Rp(i) )
     2            + (Uvel(i,jp,k)  - Uvel(i,j,k) ) / ( Ru(i) * dphi )
     2          )             -
     2   epmo * ( Ru(i) * ( Vvel(ip,jm,k)/Rp(ip) - Vvel(i,jm,k)/Rp(i))/
     2                    (drp )
     2            + (Uvel(i,j,k)   - Uvel(i,jm,k)) / ( Ru(i) * dphi )
     2          ) ) / ( Ru(i) * dphi )
     +              +
     3 ( epop * (   CNz*(Uvel(i,j,kp)  - Uvel(i,j,k) ) * dzi
     3            + (Wvel(ip,j,k)  - Wvel(i,j,k) ) / (Rp(ip) - Rp(i))
     3          )             -
     3   epom * (   CNz*(Uvel(i,j,k)   - Uvel(i,j,km)) * dzi 
     3            + (Wvel(ip,j,km) - Wvel(i,j,km)) / (Rp(ip) - Rp(i))
     3          ) ) * dzi
     +              -
     4   (ekm(i,j,k) + ekm(ip,j,k)) * ( Uvel(i,j,k) +
     4   (Vvel(ip,j,k) + Vvel(i,j,k) - Vvel(ip,jm,k) - Vvel(i,jm,k)) /
     4   (2.0 * dphi) )/ ( Ru(i) * Ru(i) )
            enddo
         enddo
      enddo
c
      return
      end

      subroutine diffv_CDS2(putout,Uvel,Vvel,Wvel,ib,ie,jb,je,kb,ke)

      USE nlist

      implicit none
  
    
!       include 'param.txt'
!       include 'common.txt'

c
c*****************************************************************
c
c      diffv calculates the diffusion of v-velocity, which is
c      the velocity in the tangential direction.
c
c
c     In formula:  (4 terms)
c
c  1  d                    1  d                       d
c  - -- (r Sigma(r phi)) + - ---- (Sigma(phi phi)) + -- (Sigma(z phi)) +
c  r dr                    r dphi                    dz
c
c
c    1
c    - Sigma(r phi)
c    r
c
c        r   : direction  ---> explicit (subroutine diffv)
c        phi : direction  ---> implicit (subroutine predic)
c        z   : direction  ---> explicit (subroutine diffv)
c
c      on input :
c
c          putout            : advection part
c          Uvel,Vvel,Wvel    : contain velocities at n-1
c          ekm               : diffusion coefficients (for velocity) in
c                              center points of the grid cells
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              diffusion part has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection and diffusion part
c          other parameters  : all unchanged
c
c      Note :
c
c          The 'source' term [ Sigma (r phi) ] / r has been
c          incorporated into the radial derivative to avoid
c          interpolation problems at the centerline.
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         eppo,empo,eopp,eopm,dzi
      real theta_U,theta_V,xx,yy,f,fluc,Wjet
      integer n,t
	real divvL,divvR,CNz

	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSE
	  CNz=1.
	ENDIF

      dzi=1./dz
	do i=ib,ie
	  ip=i+1
	  im=i-1
	      do k=kb,ke
	      kp=k+1
	      km=k-1
		do j=jb,je
		jp=j+1
		jm=j-1

      eppo = 0.25 * (
     +   ekm(i,j,k) + ekm(ip,j,k) + ekm(ip,jp,k) + ekm(i,jp,k)  )
      empo = 0.25 * (
     +   ekm(i,j,k) + ekm(im,j,k) + ekm(i,jp,k)  + ekm(im,jp,k) )
      eopp = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,kp) + ekm(i,jp,k)  + ekm(i,jp,kp) )
      eopm = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,km) + ekm(i,jp,k)  + ekm(i,jp,km) )

      divvL= ( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) )  * dzi   
      divvR= ( Ru(i)*Uvel(i,jp,k) - Ru(i-1)*Uvel(i-1,jp,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,jp,k) -         Vvel(i,jp-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       Wvel(i,jp,k) -         Wvel(i,jp,k-1) ) * dzi 

      putout(i,j,k) = putout(i,j,k) +
     1 ( eppo*( Ru(i )*Ru(i )*(Vvel(ip,j,k)/Rp(ip) - Vvel(i,j,k)/Rp(i))/
     1                      ( Rp(ip) - Rp(i) )
     1        + (Uvel(i,jp,k)  - Uvel(i,j,k) ) / (dphi)
     1        ) * Ru(i) -
     1   empo*( Ru(im)*Ru(im)*(Vvel(i,j,k)/Rp(i) - Vvel(im,j,k)/Rp(im))/
     1                      ( Rp(i) - Rp(im) )
     1        + (Uvel(im,jp,k) - Uvel(im,j,k)) / (dphi)
     1        ) * Ru(im) ) / ( Rp(i) * Rp(i) * dr(i) )
!     1 ( eppo*( Ru(i )*(Vvel(ip,j,k)/Rp(ip) - Vvel(i,j,k)/Rp(i))/
!     1                      ( Rp(ip) - Rp(i) )
!     1        + (Uvel(i,jp,k)  - Uvel(i,j,k) ) / (Ru(i)*dphi)
!     1        ) * Ru(i) -
!     1   empo*( Ru(im)*(Vvel(i,j,k)/Rp(i) - Vvel(im,j,k)/Rp(im))/
!     1                      ( Rp(i) - Rp(im) )
!     1        + (Uvel(im,jp,k) - Uvel(im,j,k)) / (Ru(im)*dphi)
!     1        ) * Ru(im) ) / ( Rp(i) * dr(i) )
     +              +
     2 ( ekm(i,jp,k) * (   (Uvel(i,jp,k) + Uvel(im,jp,k)) / 2.0
     2                   + (Vvel(i,jp,k) - Vvel(i,j,k)  ) / dphi
     2                   - 1./3.*divvR*Rp(i)
     2                 )             -
     2   ekm(i,j,k)  * (   (Uvel(i,j,k)  + Uvel(im,j,k) ) / 2.0
     2                   + (Vvel(i,j,k)  - Vvel(i,jm,k) ) / dphi
     2			 - 1./3.*divvL*Rp(i)
     2                 ) ) / ( 0.5 * Rp(i) * Rp(i) * dphi)
     +              +
     3 (   eopp * (  CNz*(Vvel(i,j,kp)  - Vvel(i,j,k) ) * dzi
     3              +(Wvel(i,jp,k)  - Wvel(i,j,k) ) / (Rp(i)*dphi)
     3            ) -
     3     eopm * (  CNz*(Vvel(i,j,k)   - Vvel(i,j,km)) * dzi
     3              +(Wvel(i,jp,km) - Wvel(i,j,km)) / (Rp(i)*dphi)
     3            )  ) * dzi

            enddo
         enddo
      enddo    
      return
      end

      subroutine diffw_CDS2(putout,Uvel,Vvel,Wvel,ib,ie,jb,je,kb,ke)
  
      USE nlist

      implicit none

!       include 'param.txt'
!       include 'common.txt'
c
c*****************************************************************
c
c      diffw calculates the diffusion of w-velocity, which is
c      the velocity in the axial direction.
c
c
c     In formula:  (3 terms)
c
c     1  d                  1  d                     d
c     - -- (r Sigma(r z)) + - ---- (Sigma(phi z)) + -- (Sigma(z z))
c     r dr                  r dphi                  dz
c
c        r   : direction  ---> explicit (subroutine diffw)
c        phi : direction  ---> implicit (subroutine predic)
c        z   : direction  ---> explicit (subroutine diffw)
c
c      on input :
c
c          putout            : advection part
c          Uvel,Vvel,Wvel    : contain velocities at n-1
c          ekm               : diffusion coefficients (for velocity) in
c                              center points of the grid cells
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              diffusion part has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection and diffusion part
c          other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         epop,emop,eopp,eomp
      real theta_U,theta_V,xx,yy,f,fluc,Wjet,dzi
      integer n,t
	real divvR,divvL,CNz

	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSE
	  CNz=1.
	ENDIF
	dzi=1./dz
        do i=ib,ie
        ip=i+1
        im=i-1
	      do k=kb,ke
	      kp=k+1
	      km=k-1
		do j=jb,je
		jp=j+1
		jm=j-1

      epop = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,kp) + ekm(ip,j,k) + ekm(ip,j,kp) )
      emop = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,kp) + ekm(im,j,k) + ekm(im,j,kp) )
      eopp = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,kp) + ekm(i,jp,k) + ekm(i,jp,kp) )
      eomp = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,kp) + ekm(i,jm,k) + ekm(i,jm,kp) )

      divvL= ( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) )   *dzi   
      divvR= ( Ru(i)*Uvel(i,j,kp) - Ru(i-1)*Uvel(i-1,j,kp) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,kp) -         Vvel(i,j-1,kp) ) / ( Rp(i)*dphi )
     +              +
     3  (       Wvel(i,j,kp) -         Wvel(i,j,kp-1) ) *dzi   

      putout(i,j,k) =  putout(i,j,k)+
     1 (Ru(i )*epop*( (Uvel(i,j,kp)  - Uvel(i,j,k)  ) *dzi
     1               +(Wvel(ip,j,k)  - Wvel(i,j,k)  ) / (Rp(ip)-Rp(i))
     1              ) -
     1  Ru(im)*emop*( (Uvel(im,j,kp) - Uvel(im,j,k) ) *dzi
     1               +(Wvel(i,j,k)   - Wvel(im,j,k) ) / (Rp(i)-Rp(im))
     1              ) ) / ( Rp(i) * dr(i) )
     +             +
     2 (  eopp * (  (Vvel(i,j,kp)  - Vvel(i,j,k)  ) *dzi
     2             +(Wvel(i,jp,k)  - Wvel(i,j,k)  ) /( Rp(i) * dphi )
     2           ) -
     2    eomp * (  (Vvel(i,jm,kp) - Vvel(i,jm,k) ) *dzi
     2             +(Wvel(i,j,k)   - Wvel(i,jm,k) )/( Rp(i) * dphi )
     2           ) ) / ( Rp(i) * dphi )
     +             +
     3 ( ekm(i,j,kp) * (CNz*(Wvel(i,j,kp) - Wvel(i,j,k ))*dzi - 1./3.*divvR ) -
     3   ekm(i,j,k ) * (CNz*(Wvel(i,j,k)  - Wvel(i,j,km))*dzi - 1./3.*divvL )
     3  	 ) *2.*dzi 
           enddo
        enddo
      enddo
      return
      end

      subroutine diffc_CDS2(putout,putin,ekh,
     +                 ib,ie,jb,je,kb,ke)

      USE nlist

      implicit none

!       include 'param.txt'
!       include 'common.txt'


c
c*****************************************************************
c
c      diffc calculates the diffusion of a scalar variable, which is
c      situated in the center point of the grid cell.
c
c
c     In formula:  (3 terms)
c
c        1  d      dC     1  d    K  dC       d    dC
c        - -- (r K -- ) + - ---- (- ---- ) + -- (K -- )
c        r dr      dr     r dphi  r dphi     dz    dz
c
c        r   : direction  ---> explicit (subroutine diffc)
c        phi : direction  ---> implicit (subroutine predic)
c        z   : direction  ---> explicit (subroutine diffc)
c
c      on input :
c
c          putout            : advection part
c          putin             : variable for which the diffusion has
c                              to be calculated at n-1
c          ekh               : diffusion coefficients (for scalar) in
c                              center points of the grid cells
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              diffusion part has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection and diffusion part
c          other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke,t,kplus
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         ekh(0:i1,0:j1,0:k1),putin2(0:i1,0:j1,0:k1)
      real	dzdz_i,Rpdr_i,Rpdphi2_i,CNz


	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSE
	  CNz=1.
	ENDIF
      putin2=putin
      do k=k1,k1 !-kjet,k1
	do t=1,tmax_inPpunt
	  i=i_inPpunt(t)
	  j=j_inPpunt(t)
!	  putin2(i,j,k) = 0. !! No diffusion in horizontal direction from jet to surroundings (now some material will diffuse into the jet, so it dissappears..)
	  putin2(i,j,k) = putin2(i,j,kmax) ! No diffusion over vertical inflow boundary, this line is needed for exact influx
	enddo
      enddo

	if (nobst>0) then
	 DO i=0,i1
	  DO j=0,j1
	    kplus = MIN(kbed(i,j)+1,k1)
	  	putin2(i,j,kbed(i,j))=putin2(i,j,kplus) ! apply neumann boundary over obstacles to get correct drift flux settling
	  ENDDO
	 ENDDO
	endif
	
	  
      dzdz_i=1./(dz*dz)
c
c     -------------------------------------------start i-loop
      do 100 i=ib,ie
c
      ip=i+1
      im=i-1
      Rpdr_i=1./(Rp(i)*dr(i))
      Rpdphi2_i=1./(Rp(i)*dphi)*1./(Rp(i)*dphi)
c
c     -------------------------------------------start j-loop
        do 200 j=jb,je
c
        jp=j+1
        jm=j-1

c
c     -------------------------------------------start k-loop
	do 300 k=kb,ke
c
	kp=k+1
	km=k-1

c
      putout(i,j,k) = putout(i,j,k) +  0.5 * (
     1 ( Ru(i)  * (ekh(i,j,k)+ekh(ip,j,k)) *
     1            (putin2(ip,j,k)-putin2(i,j,k) ) / (Rp(ip) - Rp(i))   -
     1   Ru(im) * (ekh(i,j,k)+ekh(im,j,k)) *
     1            (putin2(i,j,k) -putin2(im,j,k)) / (Rp(i) - Rp(im))    )
     1 * Rpdr_i !/ ( Rp(i) * dr(i) )
     +              +
     2 (    (ekh(i,j,k)+ekh(i,jp,k)) * (putin2(i,jp,k)-putin2(i,j,k) ) -
     2      (ekh(i,j,k)+ekh(i,jm,k)) * (putin2(i,j,k) -putin2(i,jm,k))  )
     2 * Rpdphi2_i !/ ( Rp(i) * dphi * Rp(i) * dphi )
     +              +
     3 (CNz*(ekh(i,j,k)+ekh(i,j,kp)) * (putin2(i,j,kp)-putin2(i,j,k) ) -
     3  CNz*(ekh(i,j,k)+ekh(i,j,km)) * (putin2(i,j,k) -putin2(i,j,km))  )
     3 *dzdz_i !/ ( dz * dz )
     +                                       )
c
300       continue
c     -------------------------------------------end i-loop
200     continue
c     -------------------------------------------end j-loop
100   continue
c     -------------------------------------------end k-loop
c
      return
      end
