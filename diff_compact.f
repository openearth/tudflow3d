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
	call t2np_i1(Uvel(0:ie,1:je,1:ke),Uvel_T(0:ie,1:je*px,1:ke/px))
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
	    ddy_1jepx(1)=(Vvel_T(i,1,k)-Vvel_T(i,0,k))/(Rp(i)*(phivt(1)-phivt(0)))
	    ddy_1jepx(je*px)=(Vvel_T(i,je*px,k)-Vvel_T(i,je*px-1,k))/(Rp(i)*(phivt(je*px)-phivt(je*px-1)))
	    do j=2,je*px-1
	      ddy_1jepx(j)=12./11.*(Vvel_T(i,j,k)-Vvel_T(i,j-1,k))/(Rp(i)*(phivt(j)-phivt(j-1)))
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
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
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

	call t2np_i1(Spr(0:ie,1:je,1:ke),Spr_T(0:ie,1:je*px,1:ke/px))	
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
	    ddy_0jepx(0)=(Uvel_T(i,1,k)-Uvel_T(i,0,k))/(Ru(i)*(phipt(1)-phipt(0)))
	    ddy_0jepx(je*px)=(Uvel_T(i,je*px+1,k)-Uvel_T(i,je*px,k))/(Ru(i)*(phipt(je*px+1)-phipt(je*px)))
	    do j=1,je*px-1
	      ddy_0jepx(j)=12./11.*(Uvel_T(i,j+1,k)-Uvel_T(i,j,k))/(Ru(i)*(phipt(j+1)-phipt(j)))
   	    enddo
            CALL solve_tridiag(dudrphi_T,aay_0jepx,bby_0jepx,ccy_0jepx,ddy_0jepx,je*px+1)
	    Spr_T(i,0:je*px,k)=Spr_T(i,0:je*px,k)+dudrphi_T
	  enddo
	enddo	
	do i=ib,ie
	  do k=kb,ke/px !should include k=0,but is zero because Wvel(i,j,0)=0, so left out
	    ddy_0jepx(0)=(Wvel_T(i,1,k)-Wvel_T(i,0,k))/(Rp(i)*(phipt(1)-phipt(0)))
	    ddy_0jepx(je*px)=(Wvel_T(i,je*px+1,k)-Wvel_T(i,je*px,k))/(Rp(i)*(phipt(je*px+1)-phipt(je*px)))
	    do j=1,je*px-1
	      ddy_0jepx(j)=12./11.*(Wvel_T(i,j+1,k)-Wvel_T(i,j,k))/(Rp(i)*(phipt(j+1)-phipt(j)))
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

      subroutine diffu_com4(putout,Sigrr,Sigpr,Sigzr,Sigpp,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)

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
     +         dr(0:i1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1)
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
	    ddy_1jepx(1)=(Sigpr_T(i,1,k)-Sigpr_T(i,0,k))/(Ru(i)*(phipt(1)-phipt(0)))
	    ddy_1jepx(je*px)=(Sigpr_T(i,je*px,k)-Sigpr_T(i,je*px-1,k))/(Ru(i)*(phipt(je*px+1)-phipt(je*px)))
	    do j=2,je*px-1
	      ddy_1jepx(j)=12./11.*(Sigpr_T(i,j,k)-Sigpr_T(i,j-1,k))/(Ru(i)*(phipt(j)-phipt(j-1)))
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

      subroutine diffv_com4(putout,Sigpr,Sigpp,Sigpz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)

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
     +         dr(0:i1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1)
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
	    ddy_0jepx(0)=(-25.*Sigpp_T(i,1,k)+26.*Sigpp_T(i,2,k)-Sigpp_T(i,3,k))/(Rp(i)*(phipt(1)-phipt(0)))
	    ddy_0jepx(je*px)=(25.*Sigpp_T(i,je*px,k)-26.*Sigpp_T(i,je*px-1,k)+Sigpp_T(i,je*px-2,k))/(Rp(i)*(phipt(je*px+1)-phipt(je*px)))
	    do j=1,je*px-1
	      ddy_0jepx(j)=12./11.*(Sigpp_T(i,j+1,k)-Sigpp_T(i,j,k))/(Rp(i)*(phipt(j+1)-phipt(j)))
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

      subroutine diffw_com4(putout,Sigzr,Sigpz,Sigzz,Ru,Rp,dr,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)

      implicit none

	include 'mpif.h'

      integer  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px
      real     putout(0:i1,0:j1,0:k1),Sigzz(ib:ie,jb:je,kb:ke),
     +         Sigzr(0:ie,1:je,0:ke),Sigpz(1:ie,0:je,0:ke),
     +         dr(0:i1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1)
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
	    ddy_1jepx(1)=(Sigpz_T(i,1,k)-Sigpz_T(i,0,k))/(Rp(i)*(phipt(1)-phipt(0)))
	    ddy_1jepx(je*px)=(Sigpz_T(i,je*px,k)-Sigpz_T(i,je*px-1,k))/(Rp(i)*(phipt(je*px)-phipt(je*px-1)))
	    do j=2,je*px-1
	      ddy_1jepx(j)=12./11.*(Sigpz_T(i,j,k)-Sigpz_T(i,j-1,k))/(Rp(i)*(phipt(j)-phipt(j-1)))
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
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke,km2
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         eppo,epmo,epop,epom,drp,dzi,divergentie
      real xx,yy,f,fluc,Wjet
      real divvR,divvL,CNz,CNx,CNy
      integer n,t
	  !real fcg2(0:i1,0:px*jmax+1,0:k1)
	  real fcg2(0:i1,0:j1,0:k1)
	  logical me2

	  if (momentum_exchange_obstacles.eq.100.or.momentum_exchange_obstacles.eq.110) then 
	    do j=0,j1 
		  fcg2(0:i1,j,0:k1)=fc_global(0:i1,j+rank*jmax,0:k1)
		enddo 
		!fcg2=fc_global
	  else 
		fcg2=1. !all momentum interactions are active
	  endif	  
	if (momentum_exchange_obstacles.eq.111) then 
		me2=.true. !	!111 means dUVdn = 0 for cells directly above bed
	else 
		me2=.false.
	endif 	  

	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSEIF (CNdiffz.eq.2.or.CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNz=0. 
	ELSE
	  CNz=1.
	ENDIF
	IF (CNdiffz.eq.11) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNx=0.5 !0.5*2. !0.45 !0.5
	  CNy=0.5 !0.5*2. !0.45 !0.5
	  CNz=0.5 !0.5*2.
	ELSEIF (CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNx=0. 
	  CNy=0. 
	ELSE
	  CNx=1.
	  CNy=1. 
	ENDIF	
      dzi =1./dz
      do i=ib,ie
      ip=i+1
      im=i-1
      drp = Rp(ip)-Rp(i)
	  do j=jb,je
	  jp=j+1
	  jm=j-1
		do k=kb,ke !k=MAX(kb,kbed(i,k)),ke !kb,ke
		kp=k+1
		km=k-1
  	    km2=km		  
		if (kbed(i,j).eq.km2.and.me2) km2=k !	--> dUVdn = 0 (like ordinary boundary)
		  
      eppo = 0.25 * (
     +   ekm(i,j,k) + ekm(ip,j,k) + ekm(ip,jp,k) + ekm(i,jp,k)  )
     + * MIN(fcg2(i,j,k),fcg2(ip,j,k),fcg2(ip,jp,k),fcg2(i,jp,k))	 
      epmo = 0.25 * (
     +   ekm(i,j,k) + ekm(ip,j,k) + ekm(ip,jm,k) + ekm(i,jm,k)  )
     + * MIN(fcg2(i,j,k),fcg2(ip,j,k),fcg2(ip,jm,k),fcg2(i,jm,k))	 
      epop = 0.25 * (
     +   ekm(i,j,k) + ekm(ip,j,k) + ekm(ip,j,kp) + ekm(i,j,kp)  )
     + * MIN(fcg2(i,j,k),fcg2(ip,j,k),fcg2(ip,j,kp),fcg2(i,j,kp))	 
      epom = 0.25 * (
     +   ekm(i,j,k) + ekm(ip,j,k) + ekm(ip,j,km) + ekm(i,j,km)  )
     + * MIN(fcg2(i,j,k),fcg2(ip,j,k),fcg2(ip,j,km),fcg2(i,j,km))	 

      divvL= ( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) ) * dzi  
      divvR= ( Ru(ip)*Uvel(ip,j,k) - Ru(ip-1)*Uvel(ip-1,j,k) ) / ( Rp(ip)*dr(ip) )
     +              +
     2  (       Vvel(ip,j,k) -         Vvel(ip,j-1,k) ) / ( Rp(ip)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       Wvel(ip,j,k) -         Wvel(ip,j,k-1) ) * dzi   

      putout(i,j,k) = putout(i,j,k) +
     1 (  Rp(ip) * ekm(ip,j,k) * fcg2(ip,j,k)*CNx *
!     1 (  Rp(ip) * ekm(ip,j,k)*CNx  *	 
     1            ((Uvel(ip,j,k) - Uvel(i,j,k) ) / ( dr(ip) ) - 1./3.*divvR) -
     1    Rp(i ) * ekm(i,j,k) * fcg2(i,j,k)*CNx *
!     1    Rp(i ) * ekm(i,j,k)*CNx * 	 
     1            ((Uvel(i,j,k)  - Uvel(im,j,k)) / ( dr(i)  ) - 1./3.*divvL) )  /
     1 ( 0.5 * Ru(i) * ( drp ) )
     +              +
     2 ( CNy*eppo * ( Ru(i) * ( Vvel(ip,j,k)/Rp(ip)  - Vvel(i,j,k)/Rp(i) ) /
     2                    ( Rp(ip) - Rp(i) )
     2            + (Uvel(i,jp,k)  - Uvel(i,j,k) ) / ( Ru(i) * (phip(jp)-phip(j)) )
     2          )             -
     2   CNy*epmo * ( Ru(i) * ( Vvel(ip,jm,k)/Rp(ip) - Vvel(i,jm,k)/Rp(i))/
     2                    (drp )
     2            + (Uvel(i,j,k)   - Uvel(i,jm,k)) / ( Ru(i) * (phip(j)-phip(jm)) )
     2          ) ) / ( Ru(i) * (phiv(j)-phiv(jm)) )
     +              +
     3 ( epop *CNz* ( (Uvel(i,j,kp)  - Uvel(i,j,k) ) * dzi
     3            + (Wvel(ip,j,k)  - Wvel(i,j,k) ) / (Rp(ip) - Rp(i))
     3          )             -
     3   epom *CNz* (   (Uvel(i,j,k)   - Uvel(i,j,km2)) * dzi 
     3            + (Wvel(ip,j,km) - Wvel(i,j,km)) / (Rp(ip) - Rp(i))
     3          ) ) * dzi
     +              -
     4   CNy*(ekm(i,j,k) + ekm(ip,j,k)) * MIN(fcg2(i,j,k),fcg2(ip,j,k)) * ( Uvel(i,j,k) +
!     4   CNy*(ekm(i,j,k) + ekm(ip,j,k)) * ( Uvel(i,j,k) +	 
     4   (Vvel(ip,j,k) + Vvel(i,j,k) - Vvel(ip,jm,k) - Vvel(i,jm,k)) /
     4   (2.0 * (phiv(j)-phiv(jm))) )/ ( Ru(i) * Ru(i) )
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
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke,km2
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         eppo,empo,eopp,eopm,dzi
      real xx,yy,f,fluc,Wjet
      integer n,t
	real divvL,divvR,CNz,CNx,CNy
	  !real fcg2(0:i1,0:px*jmax+1,0:k1)
	  real fcg2(0:i1,0:j1,0:k1)
	  logical me2

	  if (momentum_exchange_obstacles.eq.100.or.momentum_exchange_obstacles.eq.110) then 
	    do j=0,j1 
		  fcg2(0:i1,j,0:k1)=fc_global(0:i1,j+rank*jmax,0:k1)
		enddo 
		!fcg2=fc_global
	  else 
		fcg2=1. !all momentum interactions are active
	  endif	
	if (momentum_exchange_obstacles.eq.111) then 
		me2=.true. !	!111 means dUVdn = 0 for cells directly above bed
	else 
		me2=.false.
	endif 
	
	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSEIF (CNdiffz.eq.2.or.CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNz=0. 
	ELSE
	  CNz=1.
	ENDIF
	IF (CNdiffz.eq.11) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNx=0.5 !0.5*2. !0.45 !0.5
	  CNy=0.5 !0.5*2. !0.45 !0.5
	  CNz=0.5 !0.5*2.
	ELSEIF (CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNx=0. 
	  CNy=0. 
	ELSE
	  CNx=1.
	  CNy=1. 
	ENDIF	

      dzi=1./dz
	do i=ib,ie
	  ip=i+1
	  im=i-1
		do j=jb,je
		jp=j+1
		jm=j-1
	      do k=kb,ke !MAX(kb,kbed(i,k)),ke !do k=kb,ke
	      kp=k+1
	      km=k-1
		  km2=km		  
		  if (kbed(i,j).eq.km2.and.me2) km2=k !	--> dUVdn = 0 (like ordinary boundary)

      eppo = 0.25 * (
     +   ekm(i,j,k) + ekm(ip,j,k) + ekm(ip,jp,k) + ekm(i,jp,k)  )
     + * MIN(fcg2(i,j,k),fcg2(ip,j,k),fcg2(ip,jp,k),fcg2(i,jp,k))	 
      empo = 0.25 * (
     +   ekm(i,j,k) + ekm(im,j,k) + ekm(i,jp,k)  + ekm(im,jp,k) )
     + * MIN(fcg2(i,j,k),fcg2(im,j,k),fcg2(i,jp,k),fcg2(im,jp,k))	 
      eopp = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,kp) + ekm(i,jp,k)  + ekm(i,jp,kp) )
     + * MIN(fcg2(i,j,k),fcg2(i,j,kp),fcg2(i,jp,k),fcg2(i,jp,kp))	 
      eopm = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,km) + ekm(i,jp,k)  + ekm(i,jp,km) )
     + * MIN(fcg2(i,j,k),fcg2(i,j,km),fcg2(i,jp,k),fcg2(i,jp,km))	 

      divvL= ( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*(phiv(j)-phiv(j-1)) )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) )  * dzi   
      divvR= ( Ru(i)*Uvel(i,jp,k) - Ru(i-1)*Uvel(i-1,jp,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,jp,k) -         Vvel(i,jp-1,k) ) / ( Rp(i)*(phiv(jp)-phiv(jp-1)) )
     +              +
     3  (       Wvel(i,jp,k) -         Wvel(i,jp,k-1) ) * dzi 

      putout(i,j,k) = putout(i,j,k) +
     1 ( eppo*CNx*( Ru(i )*Ru(i )*(Vvel(ip,j,k)/Rp(ip) - Vvel(i,j,k)/Rp(i))/
     1                      ( Rp(ip) - Rp(i) )
     1        + (Uvel(i,jp,k)  - Uvel(i,j,k) ) / ((phip(jp)-phip(j)))
     1        ) * Ru(i) -
     1   empo*CNx*( Ru(im)*Ru(im)*(Vvel(i,j,k)/Rp(i) - Vvel(im,j,k)/Rp(im))/
     1                      ( Rp(i) - Rp(im) )
     1        + (Uvel(im,jp,k) - Uvel(im,j,k)) / ((phip(jp)-phip(j)))
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
     2 ( CNy*ekm(i,jp,k)*fcg2(i,jp,k) * (   (Uvel(i,jp,k) + Uvel(im,jp,k)) / 2.0
!     2 ( CNy*ekm(i,jp,k) * (   (Uvel(i,jp,k) + Uvel(im,jp,k)) / 2.0	 
     2                   + (Vvel(i,jp,k) - Vvel(i,j,k)  ) / (phiv(jp)-phiv(j))
     2                   - 1./3.*divvR*Rp(i)
     2                 )             -
     2   CNy*ekm(i,j,k)*fcg2(i,j,k)  * (   (Uvel(i,j,k)  + Uvel(im,j,k) ) / 2.0
!     2   CNy*ekm(i,j,k)  * (   (Uvel(i,j,k)  + Uvel(im,j,k) ) / 2.0	 
     2                   + (Vvel(i,j,k)  - Vvel(i,jm,k) ) / (phiv(j)-phiv(jm))
     2			 - 1./3.*divvL*Rp(i)
     2                 ) ) / ( 0.5 * Rp(i) * Rp(i) * (phip(jp)-phip(j)))
     +              +
     3 (   eopp * CNz*(  (Vvel(i,j,kp)  - Vvel(i,j,k) ) * dzi
     3              +(Wvel(i,jp,k)  - Wvel(i,j,k) ) / (Rp(i)*(phip(jp)-phip(j)))
     3            ) -
     3     eopm * CNz*(  (Vvel(i,j,k)   - Vvel(i,j,km2)) * dzi
     3              +(Wvel(i,jp,km) - Wvel(i,j,km)) / (Rp(i)*(phip(jp)-phip(j)))
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
      real xx,yy,f,fluc,Wjet,dzi
      integer n,t
	real divvR,divvL,CNz,CNx,CNy
	  !real fcg2(0:i1,0:px*jmax+1,0:k1)
	  real fcg2(0:i1,0:j1,0:k1)

	  if (momentum_exchange_obstacles.eq.100.or.momentum_exchange_obstacles.eq.110) then 
	    do j=0,j1 
		  fcg2(0:i1,j,0:k1)=fc_global(0:i1,j+rank*jmax,0:k1)
		enddo 
		!fcg2=fc_global
	  else 
		fcg2=1. !all momentum interactions are active
	  endif	

	  
	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSEIF (CNdiffz.eq.2.or.CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNz=0. 
	ELSE
	  CNz=1.
	ENDIF
	IF (CNdiffz.eq.11) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNx=0.5 !0.5*2. !0.45 !0.5
	  CNy=0.5 !0.5*2. !0.45 !0.5
	  CNz=0.5 !0.5*2.
	ELSEIF (CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNx=0. 
	  CNy=0. 
	ELSE
	  CNx=1.
	  CNy=1. 
	ENDIF	
	dzi=1./dz
        do i=ib,ie
        ip=i+1
        im=i-1
		do j=jb,je
		jp=j+1
		jm=j-1
	      do k=kb,ke !MAX(kb,kbed(i,k)),ke ! do k=kb,ke
	      kp=k+1
	      km=k-1
      epop = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,kp) + ekm(ip,j,k) + ekm(ip,j,kp) )
     + * MIN(fcg2(i,j,k),fcg2(i,j,kp),fcg2(ip,j,k),fcg2(ip,j,kp))
      emop = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,kp) + ekm(im,j,k) + ekm(im,j,kp) )
     + * MIN(fcg2(i,j,k),fcg2(i,j,kp),fcg2(im,j,k),fcg2(im,j,kp))	 
      eopp = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,kp) + ekm(i,jp,k) + ekm(i,jp,kp) )
     + * MIN(fcg2(i,j,k),fcg2(i,j,kp),fcg2(i,jp,k),fcg2(i,jp,kp))	 
      eomp = 0.25 * (
     +   ekm(i,j,k) + ekm(i,j,kp) + ekm(i,jm,k) + ekm(i,jm,kp) )
     + * MIN(fcg2(i,j,k),fcg2(i,j,kp),fcg2(i,jm,k),fcg2(i,jm,kp))	 

      divvL= ( Ru(i)*Uvel(i,j,k) - Ru(i-1)*Uvel(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,k) -         Vvel(i,j-1,k) ) / ( Rp(i)*(phiv(j+1)-phiv(j)) )
     +              +
     3  (       Wvel(i,j,k) -         Wvel(i,j,k-1) )   *dzi   
      divvR= ( Ru(i)*Uvel(i,j,kp) - Ru(i-1)*Uvel(i-1,j,kp) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vvel(i,j,kp) -         Vvel(i,j-1,kp) ) / ( Rp(i)*(phiv(j+1)-phiv(j)) )
     +              +
     3  (       Wvel(i,j,kp) -         Wvel(i,j,kp-1) ) *dzi   

      putout(i,j,k) =  putout(i,j,k)+
     1 (Ru(i )*epop*CNx*( (Uvel(i,j,kp)  - Uvel(i,j,k)  ) *dzi
     1               +(Wvel(ip,j,k)  - Wvel(i,j,k)  ) / (Rp(ip)-Rp(i))
     1              ) -
     1  Ru(im)*emop*CNx*( (Uvel(im,j,kp) - Uvel(im,j,k) ) *dzi
     1               +(Wvel(i,j,k)   - Wvel(im,j,k) ) / (Rp(i)-Rp(im))
     1              ) ) / ( Rp(i) * dr(i) )
     +             +
     2 (  eopp *CNy* (  (Vvel(i,j,kp)  - Vvel(i,j,k)  ) *dzi
     2             +(Wvel(i,jp,k)  - Wvel(i,j,k)  ) /( Rp(i) * (phip(jp)-phip(j)) )
     2           ) -
     2    eomp *CNy* (  (Vvel(i,jm,kp) - Vvel(i,jm,k) ) *dzi
     2             +(Wvel(i,j,k)   - Wvel(i,jm,k) )/( Rp(i) * (phip(j)-phip(jm)) )
     2           ) ) / ( Rp(i) * (phiv(j)-phiv(jm)) )
     +             +
     3 ( CNz*ekm(i,j,kp)*fcg2(i,j,kp) * (CNz*(Wvel(i,j,kp) - Wvel(i,j,k ))*dzi - 1./3.*divvR ) -
     3   CNz*ekm(i,j,k )*fcg2(i,j,k) * (CNz*(Wvel(i,j,k)  - Wvel(i,j,km))*dzi - 1./3.*divvL )	 
!     3 ( ekm(i,j,kp)* CNz*((Wvel(i,j,kp) - Wvel(i,j,k ))*dzi - 1./3.*divvR ) -
!     3   ekm(i,j,k )* CNz*((Wvel(i,j,k)  - Wvel(i,j,km))*dzi - 1./3.*divvL )
     3  	 ) *2.*dzi 
           enddo
        enddo
      enddo
      return
      end
	  
      subroutine diffu_CDS2_CNexpl(putout,Uvel,Vvel,Wvel,ib,ie,jb,je,kb,ke)

      USE nlist

      implicit none
      
!       include 'param.txt'
!       include 'common.txt'


c
c*****************************************************************
c
c      diffu calculates the diffusion of u-velocity, which is
c      the velocity in the radial direction.
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
c          putout            : diffusion part
c          other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         aaax,aaay,aaaz,bbbx,bbby,bbbz,cccx,cccy,cccz,ekm_min,ekm_plus
      real CNz,CNx,CNy,dzi 

	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSEIF (CNdiffz.eq.2.or.CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNz=0. 
	ELSE
	  CNz=1.
	ENDIF
	IF (CNdiffz.eq.11) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNx=0.5 !0.5*2. !0.45 !0.5
	  CNy=0.5 !0.5*2. !0.45 !0.5
	  CNz=0.5 !0.5*2.
	ELSEIF (CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNx=0. 
	  CNy=0. 
	ELSE
	  CNx=1.
	  CNy=1. 
	ENDIF	
      dzi =1./dz
      do i=ib,ie
      ip=i+1
      im=i-1
	  do j=jb,je
	  jp=j+1
	  jm=j-1
		do k=kb,ke !k=MAX(kb,kbed(i,k)),ke !kb,ke
			kp=k+1
			km=k-1
			! not *dt because this is done in solve.f 
			! d^2udx^2 
			ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
			ekm_plus=ekm(ip,j,k)!*fc_global(ip,j+rank*jmax,k)
			aaax=CNx*ekm_min*Rp(i)/(dr(i)*(Rp(ip)-Rp(i))*Ru(i))
			cccx=CNx*ekm_plus*Rp(ip)/(dr(ip)*(Rp(ip)-Rp(i))*Ru(i))
			bbbx=-aaax-cccx
			! d^2udy^2
			ekm_min= 0.25*(ekm(i,j,k)+ekm(i,jm,k)+ekm(i+1,j,k)+ekm(i+1,jm,k))
			ekm_plus=0.25*(ekm(i,j,k)+ekm(i,jp,k)+ekm(i+1,j,k)+ekm(i+1,jp,k))
			aaay=CNy*ekm_min/(Ru(i)*(phip(j)-phip(jm))*Ru(i)*dphi2(j))
			cccy=CNy*ekm_plus/(Ru(i)*(phip(jp)-phip(j))*Ru(i)*dphi2(j))
			bbby=-aaay-cccy		
			! d^2udz^2 
			ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))
			ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))
			aaaz=CNz*ekm_min/dz**2
			cccz=CNz*ekm_plus/dz**2
			bbbz=-aaaz-cccz
	 
			putout(i,j,k) =  !putout(i,j,k)+
     1      aaax*Uvel(im,j,k) + bbbx*Uvel(i,j,k) + cccx*Uvel(ip,j,k) +				
     2      aaay*Uvel(i,jm,k) + bbby*Uvel(i,j,k) + cccy*Uvel(i,jp,k) +
     3      aaaz*Uvel(i,j,km) + bbbz*Uvel(i,j,k) + cccz*Uvel(i,j,kp)		 
		  enddo
        enddo
      enddo
c
      return
      end	  
	  
      subroutine diffv_CDS2_CNexpl(putout,Uvel,Vvel,Wvel,ib,ie,jb,je,kb,ke)

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
c          putout            : diffusion part
c          other parameters  : all unchanged
c
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         aaax,aaay,aaaz,bbbx,bbby,bbbz,cccx,cccy,cccz,ekm_min,ekm_plus
	real CNz,CNx,CNy,dzi 

	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSEIF (CNdiffz.eq.2.or.CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNz=0. 
	ELSE
	  CNz=1.
	ENDIF
	IF (CNdiffz.eq.11) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNx=0.5 !0.5*2. !0.45 !0.5
	  CNy=0.5 !0.5*2. !0.45 !0.5
	  CNz=0.5 !0.5*2.
	ELSEIF (CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNx=0. 
	  CNy=0. 
	ELSE
	  CNx=1.
	  CNy=1. 
	ENDIF	

      dzi=1./dz
	do i=ib,ie
	  ip=i+1
	  im=i-1
		do j=jb,je
		jp=j+1
		jm=j-1
	      do k=kb,ke !MAX(kb,kbed(i,k)),ke !do k=kb,ke
			kp=k+1
			km=k-1		
			! not *dt because this is done in solve.f 
			! d^2vdx^2 
			ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j+1,k)+ekm(im,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(im,j+1+rank*jmax,k))		
			ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j+1,k)+ekm(ip,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(ip,j+1+rank*jmax,k))	
			aaax=CNx*ekm_min*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))
			cccx=CNx*ekm_plus*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))
			bbbx=-aaax-cccx
			! d^2vdy^2
			ekm_min= ekm(i,j,k)
			ekm_plus=ekm(i,jp,k)
			aaay=CNy*ekm_min/(Rp(i)*dphi2(j)*Rp(i)*(phip(jp)-phip(j)))
			cccy=CNy*ekm_plus/(Rp(i)*dphi2(jp)*Rp(i)*(phip(jp)-phip(j)))
			bbby=-aaay-cccy		
			! d^2vdz^2 
			ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))		
			ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))		
			aaaz=CNz*ekm_min/dz**2
			cccz=CNz*ekm_plus/dz**2
			bbbz=-aaaz-cccz
	 
			putout(i,j,k) = !putout(i,j,k)+	
     1      aaax*Vvel(im,j,k) + bbbx*Vvel(i,j,k) + cccx*Vvel(ip,j,k) +				
     2      aaay*Vvel(i,jm,k) + bbby*Vvel(i,j,k) + cccy*Vvel(i,jp,k) +
     3      aaaz*Vvel(i,j,km) + bbbz*Vvel(i,j,k) + cccz*Vvel(i,j,kp)	 

          enddo
        enddo
      enddo    
      return
      end	  
	  
      subroutine diffw_CDS2_CNexpl(putout,Uvel,Vvel,Wvel,ib,ie,jb,je,kb,ke)
  
      USE nlist

      implicit none

!       include 'param.txt'
!       include 'common.txt'
c
c*****************************************************************
c
c      diffw calculates the diffusion of w-velocity, which is
c      the velocity in the vertical direction.
c
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
c          putout            : diffusion part
c          other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         aaax,aaay,aaaz,bbbx,bbby,bbbz,cccx,cccy,cccz,ekm_min,ekm_plus
      real dzi
      integer n,t
	real CNz,CNx,CNy

	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSEIF (CNdiffz.eq.2.or.CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNz=0. 
	ELSE
	  CNz=1.
	ENDIF
	IF (CNdiffz.eq.11) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNx=0.5 !0.5*2. !0.45 !0.5
	  CNy=0.5 !0.5*2. !0.45 !0.5
	  CNz=0.5 !0.5*2.
	ELSEIF (CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNx=0. 
	  CNy=0. 
	ELSE
	  CNx=1.
	  CNy=1. 
	ENDIF	
	dzi=1./dz
        do i=ib,ie
        ip=i+1
        im=i-1
		do j=jb,je
		jp=j+1
		jm=j-1
	      do k=kb,ke !MAX(kb,kbed(i,k)),ke ! do k=kb,ke
			kp=k+1
			km=k-1
			! not *dt because this is done in solve.f 
			! d^2wdx^2 
			ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j,k+1)+ekm(im,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
	!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(im,j+rank*jmax,k+1))		
			ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j,k+1)+ekm(ip,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
	!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(ip,j+rank*jmax,k+1))	
			aaax=CNx*ekm_min*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))
			cccx=CNx*ekm_plus*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))
			bbbx=-aaax-cccx
			! d^2wdy^2
			ekm_min= 0.25*(ekm(i,j,k)+ekm(i,jm,k)+ekm(i,j,k+1)+ekm(i,jm,k+1))
			ekm_plus=0.25*(ekm(i,j,k)+ekm(i,jp,k)+ekm(i,j,k+1)+ekm(i,jp,k+1))
			aaay=CNy*ekm_min/(Rp(i)*(phip(j)-phip(jm))*Rp(i)*dphi2(j))
			cccy=CNy*ekm_plus/(Rp(i)*(phip(jp)-phip(j))*Rp(i)*dphi2(j))
			bbby=-aaay-cccy		
			! d^2wdz^2 
			ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
			ekm_plus=ekm(i,j,kp)!*fc_global(i,j+rank*jmax,kp)
			aaaz=CNz*ekm_min/dz**2
			cccz=CNz*ekm_plus/dz**2
			bbbz=-aaaz-cccz 				
			
			putout(i,j,k) = !putout(i,j,k)+ 	
     1      aaax*Wvel(im,j,k) + bbbx*Wvel(i,j,k) + cccx*Wvel(ip,j,k) +				
     2      aaay*Wvel(i,jm,k) + bbby*Wvel(i,j,k) + cccy*Wvel(i,jp,k) +
     3      aaaz*Wvel(i,j,km) + bbbz*Wvel(i,j,k) + cccz*Wvel(i,j,kp)	
           enddo
        enddo
      enddo
      return
      end	  

      subroutine diffuvw_xdir_CDS2_CNimpl

        USE nlist

        implicit none

        integer  im,ip,jm,jp,km,kp !,ib,ie,jb,je,kb,ke
        real CNz
		real ekm_min,ekm_plus,rho_min,rho_plus
		real aaax(0:i1),bbbx(0:i1),cccx(0:i1),rhssx(0:i1)
		integer perx,pery 


		IF (CNdiffz.eq.1) THEN 
			CNz=0.5 !0.55
		ELSEIF (CNdiffz.eq.11) THEN 
			CNz=CNdiff_factor !0.5 !0.55 
		ELSEIF (CNdiffz.eq.12) THEN 
			CNz=(3.-2.*CNdiff_factor)/3.
		ELSE 
			CNz=1.
		ENDIF 
		perx=0
		pery=0 		
		IF (periodicx.eq.1) perx=1
		IF (periodicy.eq.1) pery=1		
!		 if (.false.) then
!!		 if (CNdiffz.eq.11) then 
!			!Hirsch 2007 2nd edition Eq. 9.4.21 p.478 mentions that also explicit part of CN should build on for 3 subsequent directions
!			!here x-dir is done as 2nd step on velocity field with finished CN explicit and implicit actions in y-dir:
!			call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 	Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
!			call diffuvw_xdir_CDS2_CNexpl(dnew,dnew2,dold,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
!			dUdt = dUdt + dt*dnew/rhU  
!			dVdt = dVdt + dt*dnew2/rhV
!			dWdt = dWdt + dt*dold/rhW
!			call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 	Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)			
!		 endif  		
		 

		IF (periodicx.eq.0) THEN !inflow bc i=0 and Neumann outflow i=imax 
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax-1 
				!im=MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
				ekm_plus=ekm(ip,j,k)!*fc_global(ip,j+rank*jmax,k)
				aaax(i)=-CNz*ekm_min*Rp(i)*dt/(dr(i)*(Rp(ip)-Rp(i))*Ru(i))/rhU(i,j,k) !/drdt(i,j,k)
				cccx(i)=-CNz*ekm_plus*Rp(ip)*dt/(dr(ip)*(Rp(ip)-Rp(i))*Ru(i))/rhU(i,j,k) !/drdt(ip,j,k)
				bbbx(i)=1.-aaax(i)-cccx(i) 
            enddo
			aaax(0)=0.
			cccx(0)=0.
			bbbx(0)=1.
			bbbx(imax-1)=bbbx(imax-1)+cccx(imax-1) !d.dn=0; should be applied on imax-1 because staggered positino dUdt(imax) and also explicit bc is U(imax)=U(imax-1)
			cccx(imax-1)=0. 
			rhssx(0:imax-1)=dUdt(0:imax-1,j,k)
			CALL solve_tridiag_switchperiodic2(dUdt(0:imax-1,j,k),aaax(0:imax-1),bbbx(0:imax-1),
     &			cccx(0:imax-1),rhssx(0:imax-1),imax,perx) 
           enddo
		  enddo
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax !0,i1
				im=i-1 !MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j+1,k)+ekm(im,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(im,j+1+rank*jmax,k))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j+1,k)+ekm(ip,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(ip,j+1+rank*jmax,k))	
				rho_min=0.25*(drdt(i,j,k)+drdt(im,j,k)+drdt(i,j+1,k)+drdt(im,j+1,k))
				rho_plus=0.25*(drdt(i,j,k)+drdt(ip,j,k)+drdt(i,j+1,k)+drdt(ip,j+1,k))
				aaax(i)=-CNz*ekm_min*dt*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))/rhV(i,j,k) !/rho_min
				cccx(i)=-CNz*ekm_plus*dt*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))/rhV(i,j,k) !/rho_plus
				bbbx(i)=1.-aaax(i)-cccx(i)
			enddo
			aaax(0)=0.
			cccx(0)=0.
			bbbx(0)=1.
			bbbx(imax)=bbbx(imax)+cccx(imax)
			cccx(imax)=0.
!!!			rhssx(0:imax)=dVdt2(0:imax,j,k)
			rhssx(0:imax)=dVdt(0:imax,j,k)	
			CALL solve_tridiag_switchperiodic2(dVdt(0:imax,j,k),aaax(0:imax),bbbx(0:imax),
     &			cccx(0:imax),rhssx(0:imax),i1,perx)
           enddo
		  enddo
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax !0,i1
				im=i-1 !MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j,k+1)+ekm(im,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(im,j+rank*jmax,k+1))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j,k+1)+ekm(ip,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(ip,j+rank*jmax,k+1))	
				rho_min=0.25*(drdt(i,j,k)+drdt(im,j,k)+drdt(i,j,k+1)+drdt(im,j,k+1))
				rho_plus=0.25*(drdt(i,j,k)+drdt(ip,j,k)+drdt(i,j,k+1)+drdt(ip,j,k+1))
				aaax(i)=-CNz*ekm_min*dt*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))/rhW(i,j,k) !/rho_min
				cccx(i)=-CNz*ekm_plus*dt*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))/rhW(i,j,k) !/rho_plus
				bbbx(i)=1.-aaax(i)-cccx(i)				
            enddo
			aaax(0)=0.
			cccx(0)=0.
			bbbx(0)=1.
			bbbx(imax)=bbbx(imax)+cccx(imax)
			cccx(imax)=0.
!!!			rhssx(0:imax)=dWdt2(0:imax,j,k)
			rhssx(0:imax)=dWdt(0:imax,j,k)
			CALL solve_tridiag_switchperiodic2(dWdt(0:imax,j,k),aaax(0:imax),bbbx(0:imax),
     &			cccx(0:imax),rhssx(0:imax),i1,perx)
		   enddo
		  enddo	
		ELSEIF (periodicx.eq.1) THEN 
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax 
				!im=MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
				ekm_plus=ekm(ip,j,k)!*fc_global(ip,j+rank*jmax,k)
				aaax(i)=-CNz*ekm_min*Rp(i)*dt/(dr(i)*(Rp(ip)-Rp(i))*Ru(i))/rhU(i,j,k) !/drdt(i,j,k)
				cccx(i)=-CNz*ekm_plus*Rp(ip)*dt/(dr(ip)*(Rp(ip)-Rp(i))*Ru(i))/rhU(i,j,k) !/drdt(ip,j,k)
				bbbx(i)=1.-aaax(i)-cccx(i) 
            enddo
!!!			rhssx(1:imax)=dUdt2(1:imax,j,k)
			rhssx(1:imax)=dUdt(1:imax,j,k)
			CALL solve_tridiag_switchperiodic2(dUdt(1:imax,j,k),aaax(1:imax),bbbx(1:imax),
     &			cccx(1:imax),rhssx(1:imax),imax,perx) 
           enddo
		  enddo
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax !0,i1
				im=i-1 !MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j+1,k)+ekm(im,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(im,j+1+rank*jmax,k))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j+1,k)+ekm(ip,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(ip,j+1+rank*jmax,k))	
				rho_min=0.25*(drdt(i,j,k)+drdt(im,j,k)+drdt(i,j+1,k)+drdt(im,j+1,k))
				rho_plus=0.25*(drdt(i,j,k)+drdt(ip,j,k)+drdt(i,j+1,k)+drdt(ip,j+1,k))
				aaax(i)=-CNz*ekm_min*dt*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))/rhV(i,j,k) !/rho_min
				cccx(i)=-CNz*ekm_plus*dt*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))/rhV(i,j,k) !/rho_plus
				bbbx(i)=1.-aaax(i)-cccx(i)
			enddo
!!!			rhssx(1:imax)=dVdt2(1:imax,j,k)
			rhssx(1:imax)=dVdt(1:imax,j,k)
			CALL solve_tridiag_switchperiodic2(dVdt(1:imax,j,k),aaax(1:imax),bbbx(1:imax),
     &			cccx(1:imax),rhssx(1:imax),imax,perx)
           enddo
		  enddo
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax !0,i1
				im=i-1 !MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j,k+1)+ekm(im,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(im,j+rank*jmax,k+1))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j,k+1)+ekm(ip,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(ip,j+rank*jmax,k+1))	
				rho_min=0.25*(drdt(i,j,k)+drdt(im,j,k)+drdt(i,j,k+1)+drdt(im,j,k+1))
				rho_plus=0.25*(drdt(i,j,k)+drdt(ip,j,k)+drdt(i,j,k+1)+drdt(ip,j,k+1))
				aaax(i)=-CNz*ekm_min*dt*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))/rhW(i,j,k) !/rho_min
				cccx(i)=-CNz*ekm_plus*dt*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))/rhW(i,j,k) !/rho_plus
				bbbx(i)=1.-aaax(i)-cccx(i)				
            enddo
!!!			rhssx(1:imax)=dWdt2(1:imax,j,k)
			rhssx(1:imax)=dWdt(1:imax,j,k)
			CALL solve_tridiag_switchperiodic2(dWdt(1:imax,j,k),aaax(1:imax),bbbx(1:imax),
     &			cccx(1:imax),rhssx(1:imax),imax,perx)
		   enddo
		  enddo	
		ELSEIF(periodicx.eq.2) THEN !U=0 i=0 and i=imax + V=-V at i=0 and i=imax + W=-W at i=0 and W(i1)=2*Wbc-W(imax)
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax-1 
				!im=MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
				ekm_plus=ekm(ip,j,k)!*fc_global(ip,j+rank*jmax,k)
				aaax(i)=-CNz*ekm_min*Rp(i)*dt/(dr(i)*(Rp(ip)-Rp(i))*Ru(i))/rhU(i,j,k) !/drdt(i,j,k)
				cccx(i)=-CNz*ekm_plus*Rp(ip)*dt/(dr(ip)*(Rp(ip)-Rp(i))*Ru(i))/rhU(i,j,k) !/drdt(ip,j,k)
				bbbx(i)=1.-aaax(i)-cccx(i) 
            enddo
			aaax(0)=0.
			cccx(0)=0.
			bbbx(0)=1.
			bbbx(imax)=1. 
			aaax(imax)=0. 
			cccx(imax)=0. 
!!!			rhssx(0:imax)=dUdt2(0:imax,j,k)
			rhssx(0:imax)=dUdt(0:imax,j,k)
			rhssx(0)=0.
			rhssx(imax)=0.
			CALL solve_tridiag_switchperiodic2(dUdt(0:imax,j,k),aaax(0:imax),bbbx(0:imax),
     &			cccx(0:imax),rhssx(0:imax),i1,perx) 
           enddo
		  enddo
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax !0,i1
				im=i-1 !MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j+1,k)+ekm(im,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(im,j+1+rank*jmax,k))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j+1,k)+ekm(ip,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(ip,j+1+rank*jmax,k))	
				rho_min=0.25*(drdt(i,j,k)+drdt(im,j,k)+drdt(i,j+1,k)+drdt(im,j+1,k))
				rho_plus=0.25*(drdt(i,j,k)+drdt(ip,j,k)+drdt(i,j+1,k)+drdt(ip,j+1,k))
				aaax(i)=-CNz*ekm_min*dt*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))/rhV(i,j,k) !/rho_min
				cccx(i)=-CNz*ekm_plus*dt*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))/rhV(i,j,k) !/rho_plus
				bbbx(i)=1.-aaax(i)-cccx(i)
			enddo
			bbbx(1)=bbbx(1)-aaax(1)
			aaax(1)=0.
			bbbx(imax)=bbbx(imax)-cccx(imax)
			cccx(imax)=0.
!!!			rhssx(1:imax)=dVdt2(1:imax,j,k)
			rhssx(1:imax)=dVdt(1:imax,j,k)
			CALL solve_tridiag_switchperiodic2(dVdt(1:imax,j,k),aaax(1:imax),bbbx(1:imax),
     &			cccx(1:imax),rhssx(1:imax),imax,perx)
           enddo
		  enddo
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax !0,i1
				im=i-1 !MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j,k+1)+ekm(im,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(im,j+rank*jmax,k+1))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j,k+1)+ekm(ip,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(ip,j+rank*jmax,k+1))	
				rho_min=0.25*(drdt(i,j,k)+drdt(im,j,k)+drdt(i,j,k+1)+drdt(im,j,k+1))
				rho_plus=0.25*(drdt(i,j,k)+drdt(ip,j,k)+drdt(i,j,k+1)+drdt(ip,j,k+1))
				aaax(i)=-CNz*ekm_min*dt*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))/rhW(i,j,k) !/rho_min
				cccx(i)=-CNz*ekm_plus*dt*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))/rhW(i,j,k) !/rho_plus
				bbbx(i)=1.-aaax(i)-cccx(i)				
            enddo
			bbbx(1)=bbbx(1)-aaax(1)
			aaax(1)=0.
			bbbx(imax)=bbbx(imax)-cccx(imax)  			!force W(i1)=2*W_ox-W(imax) 
!!!			rhssx(1:imax)=dWdt2(1:imax,j,k)
			rhssx(1:imax)=dWdt(1:imax,j,k)
			rhssx(imax)=rhssx(imax)-2.*cccx(imax)*W_ox!*rho_b !force W(i1)=2*W_ox-W(imax) 
!			rhssx(imax)=rhssx(imax)-2.*cccx(imax)*W_ox 
!     & 			+cccx(imax)*(Wnew(imax,j,k)+Wnew(i1,j,k))*(12.-CNdiffz) !force W^n+1(i1)=2*W_ox-W^n+1(imax)
			!W^n+1(i1)=2*W_ox-W^n+1(imax) --> [W^n+1(i1)-W^n(i1)]=2*W_ox-[W^n+1(imax)-W^n(imax)]-W^n(imax)-W^n(i1) 
			cccx(imax)=0.			
			CALL solve_tridiag_switchperiodic2(dWdt(1:imax,j,k),aaax(1:imax),bbbx(1:imax),
     &			cccx(1:imax),rhssx(1:imax),imax,perx)
		   enddo
		  enddo	
		ENDIF 
		 
		return
		end	
	  
      subroutine diffuvw_ydir_CDS2_CNimpl

        USE nlist

        implicit none
		
		include 'mpif.h'

        integer  im,ip,jm,jp,km,kp !,ib,ie,jb,je,kb,ke
        real CNz,CNx,CNy,dzi 
		real ekm_min,ekm_plus,rho_min,rho_plus
		real aaay(0:jmax*px+1),bbby(0:jmax*px+1),cccy(0:jmax*px+1),rhssy(0:jmax*px+1)
		real ans_T(1:i1,0:jmax*px+1,1:kmax/px+1),ekm_T(1:i1,0:jmax*px+1,1:kmax/px+1)!,ekm2(0:i1,0:j1,0:k1) 
		real uu_T(1:imax,0:jmax*px+1,1:kmax/px),vv_T(1:imax,0:jmax*px+1,1:kmax/px),ww_T(1:imax,0:jmax*px+1,1:kmax/px)
		real rhU_T(1:imax,0:jmax*px+1,1:kmax/px),rhV_T(1:imax,0:jmax*px+1,1:kmax/px),rhW_T(1:imax,0:jmax*px+1,1:kmax/px)
		real interface_T(1:i1,0:jmax*px+1)
		integer ileng,ierr,itag,status(MPI_STATUS_SIZE),perx,pery 


		IF (CNdiffz.eq.1) THEN 
			CNz=0.5 !0.55
		ELSEIF (CNdiffz.eq.11) THEN 
			CNz=CNdiff_factor !0.5 !0.55 
		ELSEIF (CNdiffz.eq.12) THEN 
			CNz=(3.-2.*CNdiff_factor)/3.
		ELSE 
			CNz=1.
		ENDIF 
		perx=0
		pery=0 		
		IF (periodicx.eq.1) perx=1
		IF (periodicy.eq.1) pery=1		
		
!		if (CNdiffz.eq.11) then 
!			!Hirsch 2007 2nd edition Eq. 9.4.21 p.478 mentions that also explicit part of CN should build on for 3 subsequent directions
!			!here y-dir is done as first step on velocity field previous time step:
!			!call diffuvw_ydir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
!			call diffuvw_ydir_CDS2_CNexpl(dnew,dnew2,dold,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
!			dUdt = dUdt + dt*dnew  
!			dVdt = dVdt + dt*dnew2
!			dWdt = dWdt + dt*dold	
!		endif 		
		!ekm2=ekm/drdt   
		!call t2np_i1(ekm(1:i1,1:jmax,1:kmax),ekm_T(1:i1,1:jmax*px,1:kmax/px))
		!call t2np_i1(drdt(1:i1,1:jmax,1:kmax),ans_T(1:i1,1:jmax*px,1:kmax/px))
		call t2np_i1(ekm(1:i1,1:jmax,1:kmax),ans_T(1:i1,1:jmax*px,1:kmax/px))
		call t2np(dUdt(1:imax,1:jmax,1:kmax),uu_T(1:imax,1:jmax*px,1:kmax/px))
		call t2np(dVdt(1:imax,1:jmax,1:kmax),vv_T(1:imax,1:jmax*px,1:kmax/px))
		call t2np(dWdt(1:imax,1:jmax,1:kmax),ww_T(1:imax,1:jmax*px,1:kmax/px))
		call t2np(rhU(1:imax,1:jmax,1:kmax),rhU_T(1:imax,1:jmax*px,1:kmax/px))
		call t2np(rhV(1:imax,1:jmax,1:kmax),rhV_T(1:imax,1:jmax*px,1:kmax/px))
		call t2np(rhW(1:imax,1:jmax,1:kmax),rhW_T(1:imax,1:jmax*px,1:kmax/px))		  
		!! also pass over boundaries at j=0 :
		IF (rank.eq.0) THEN
			do i=1,px-1
			  !call mpi_send(ekm(1:i1,0,i*kmax/px+1:(i+1)*kmax/px),i1*kmax/px,MPI_REAL8,i,i+1000,MPI_COMM_WORLD,status,ierr)
			  call mpi_send(dUdt(1:imax,0,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+2000,MPI_COMM_WORLD,status,ierr)
			  call mpi_send(dVdt(1:imax,0,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+3000,MPI_COMM_WORLD,status,ierr)
			  call mpi_send(dWdt(1:imax,0,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+4000,MPI_COMM_WORLD,status,ierr)
			  call mpi_send(ekm(1:i1,0,i*kmax/px+1:(i+1)*kmax/px),i1*kmax/px,MPI_REAL8,i,i+5000,MPI_COMM_WORLD,status,ierr)
!			  call mpi_send(rhU(1:imax,0,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+6000,MPI_COMM_WORLD,status,ierr)
!			  call mpi_send(rhV(1:imax,0,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+7000,MPI_COMM_WORLD,status,ierr)
!			  call mpi_send(rhW(1:imax,0,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+8000,MPI_COMM_WORLD,status,ierr)			  
			enddo
			!ekm_T(1:i1,0,1:kmax/px)=ekm2(1:i1,0,1:kmax/px)
			uu_T(1:imax,0,1:kmax/px)=dUdt(1:imax,0,1:kmax/px)			
			vv_T(1:imax,0,1:kmax/px)=dVdt(1:imax,0,1:kmax/px)			
			ww_T(1:imax,0,1:kmax/px)=dWdt(1:imax,0,1:kmax/px)			
			ans_T(1:i1,0,1:kmax/px)=ekm(1:i1,0,1:kmax/px)
!			rhU_T(1:imax,0,1:kmax/px)=rhU(1:imax,0,1:kmax/px)			
!			rhV_T(1:imax,0,1:kmax/px)=rhV(1:imax,0,1:kmax/px)			
!			rhW_T(1:imax,0,1:kmax/px)=rhW(1:imax,0,1:kmax/px)				
		ELSE
			!call mpi_recv(ekm_T(1:i1,0,1:kmax/px),i1*kmax/px,MPI_REAL8,0,rank+1000,MPI_COMM_WORLD,status,ierr)
			call mpi_recv(uu_T(1:imax,0,1:kmax/px),imax*kmax/px,MPI_REAL8,0,rank+2000,MPI_COMM_WORLD,status,ierr)
			call mpi_recv(vv_T(1:imax,0,1:kmax/px),imax*kmax/px,MPI_REAL8,0,rank+3000,MPI_COMM_WORLD,status,ierr)
			call mpi_recv(ww_T(1:imax,0,1:kmax/px),imax*kmax/px,MPI_REAL8,0,rank+4000,MPI_COMM_WORLD,status,ierr)
			call mpi_recv(ans_T(1:i1,0,1:kmax/px),i1*kmax/px,MPI_REAL8,0,rank+5000,MPI_COMM_WORLD,status,ierr)
!			call mpi_recv(rhU_T(1:imax,0,1:kmax/px),imax*kmax/px,MPI_REAL8,0,rank+6000,MPI_COMM_WORLD,status,ierr)
!			call mpi_recv(rhV_T(1:imax,0,1:kmax/px),imax*kmax/px,MPI_REAL8,0,rank+7000,MPI_COMM_WORLD,status,ierr)
!			call mpi_recv(rhW_T(1:imax,0,1:kmax/px),imax*kmax/px,MPI_REAL8,0,rank+8000,MPI_COMM_WORLD,status,ierr)			
		ENDIF
		!! also pass over boundaries at j=jmax+1 :
		IF (rank.eq.px-1) THEN
		    do i=0,px-2
			 !call mpi_send(ekm (1:i1,jmax+1,i*kmax/px+1:(i+1)*kmax/px),i1*kmax/px,MPI_REAL8,i,i+1000,MPI_COMM_WORLD,status,ierr)
			 call mpi_send(dUdt(1:imax,jmax+1,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+2000,MPI_COMM_WORLD,status,ierr)	 
			 call mpi_send(dVdt(1:imax,jmax+1,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+3000,MPI_COMM_WORLD,status,ierr)
			 call mpi_send(dWdt(1:imax,jmax+1,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+4000,MPI_COMM_WORLD,status,ierr)
			 call mpi_send(ekm(1:i1,jmax+1,i*kmax/px+1:(i+1)*kmax/px),i1*kmax/px,MPI_REAL8,i,i+5000,MPI_COMM_WORLD,status,ierr)
	!		 call mpi_send(rhU(1:imax,jmax+1,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+6000,MPI_COMM_WORLD,status,ierr)	 
	!		 call mpi_send(rhV(1:imax,jmax+1,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+7000,MPI_COMM_WORLD,status,ierr)
	!		 call mpi_send(rhW(1:imax,jmax+1,i*kmax/px+1:(i+1)*kmax/px),imax*kmax/px,MPI_REAL8,i,i+8000,MPI_COMM_WORLD,status,ierr)		 
			enddo
			  !ekm_T(1:i1,jmax*px+1,1:kmax/px)=ekm(1:i1,jmax+1,rank*kmax/px+1:(rank+1)*kmax/px)
			  uu_T(1:imax,jmax*px+1,1:kmax/px)=dUdt(1:imax,jmax+1,rank*kmax/px+1:(rank+1)*kmax/px)
			  vv_T(1:imax,jmax*px+1,1:kmax/px)=dVdt(1:imax,jmax+1,rank*kmax/px+1:(rank+1)*kmax/px)
			  ww_T(1:imax,jmax*px+1,1:kmax/px)=dWdt(1:imax,jmax+1,rank*kmax/px+1:(rank+1)*kmax/px)
			  ans_T(1:i1,jmax*px+1,1:kmax/px)=ekm(1:i1,jmax+1,rank*kmax/px+1:(rank+1)*kmax/px)
!			  rhU_T(1:imax,jmax*px+1,1:kmax/px)=rhU(1:imax,jmax+1,rank*kmax/px+1:(rank+1)*kmax/px)
!			  rhV_T(1:imax,jmax*px+1,1:kmax/px)=rhV(1:imax,jmax+1,rank*kmax/px+1:(rank+1)*kmax/px)
!			  rhW_T(1:imax,jmax*px+1,1:kmax/px)=rhW(1:imax,jmax+1,rank*kmax/px+1:(rank+1)*kmax/px)			  
		ELSE
		    !call mpi_recv(ekm_T(1:i1,jmax*px+1,1:kmax/px),i1*kmax/px,MPI_REAL8,px-1,rank+1000,MPI_COMM_WORLD,status,ierr)
		    call mpi_recv(uu_T(1:imax,jmax*px+1,1:kmax/px),imax*kmax/px,MPI_REAL8,px-1,rank+2000,MPI_COMM_WORLD,status,ierr)
			call mpi_recv(vv_T(1:imax,jmax*px+1,1:kmax/px),imax*kmax/px,MPI_REAL8,px-1,rank+3000,MPI_COMM_WORLD,status,ierr)
			call mpi_recv(ww_T(1:imax,jmax*px+1,1:kmax/px),imax*kmax/px,MPI_REAL8,px-1,rank+4000,MPI_COMM_WORLD,status,ierr)
			call mpi_recv(ans_T(1:i1,jmax*px+1,1:kmax/px),i1*kmax/px,MPI_REAL8,px-1,rank+5000,MPI_COMM_WORLD,status,ierr)
!		    call mpi_recv(rhU_T(1:imax,jmax*px+1,1:kmax/px),imax*kmax/px,MPI_REAL8,px-1,rank+6000,MPI_COMM_WORLD,status,ierr)
!			call mpi_recv(rhV_T(1:imax,jmax*px+1,1:kmax/px),imax*kmax/px,MPI_REAL8,px-1,rank+7000,MPI_COMM_WORLD,status,ierr)
!			call mpi_recv(rhW_T(1:imax,jmax*px+1,1:kmax/px),imax*kmax/px,MPI_REAL8,px-1,rank+8000,MPI_COMM_WORLD,status,ierr)			
		ENDIF
		!! also pass over boundaries at k=kmax/px+1:
		!call shiftb_T(ekm_T,interface_T)
		!if (rank.eq.px-1) then 
		!ekm_T(1:i1,0:jmax*px+1,kmax/px+1)=ekm_T(1:i1,0:jmax*px+1,kmax/px)
		!else 
		!  ekm_T(1:i1,0:jmax*px+1,kmax/px+1)=interface_T(1:i1,0:jmax*px+1)
		!endif 
		call shiftb_T(ans_T,interface_T)
		if (rank.eq.px-1) then 
			ans_T(1:i1,0:jmax*px+1,kmax/px+1)=ans_T(1:i1,0:jmax*px+1,kmax/px)
		else 
			ans_T(1:i1,0:jmax*px+1,kmax/px+1)=interface_T(1:i1,0:jmax*px+1)
		endif 
		IF (periodicy.eq.0) THEN !bc defined at both lateral boundaries 
		  do k=1,kmax/px 
           do i=1,imax
            do j=1,px*jmax !0,px*jmax+1
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= 0.25*(ans_T(i,j,k)+ans_T(i,jm,k)+ans_T(i+1,j,k)+ans_T(i+1,jm,k))
				ekm_plus=0.25*(ans_T(i,j,k)+ans_T(i,jp,k)+ans_T(i+1,j,k)+ans_T(i+1,jp,k))
				aaay(j)=-CNz*ekm_min*dt/(Ru(i)*(phipt(j)-phipt(jm))*Ru(i)*dphi2t(j))/rhU_T(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Ru(i)*(phipt(jp)-phipt(j))*Ru(i)*dphi2t(j))/rhU_T(i,j,k)				
				bbby(j)=1.-aaay(j)-cccy(j) 
			enddo
			aaay(0)=0.
			aaay(px*jmax+1)=0.
			cccy(0)=0.
			cccy(px*jmax+1)=0.
			bbby(0)=1.
			bbby(px*jmax+1)=1.
			rhssy=uu_T(i,0:px*jmax+1,k)
			CALL solve_tridiag_switchperiodic2(uu_T(i,0:px*jmax+1,k),aaay,bbby,cccy,rhssy,px*jmax+2,pery)
           enddo
		  enddo
		  do k=1,kmax/px 
           do i=1,imax
            do j=1,px*jmax-1 !0,px*jmax+1
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= ans_T(i,j,k)
				ekm_plus=ans_T(i,jp,k)
				aaay(j)=-CNz*ekm_min*dt/(Rp(i)*dphi2t(j)*Rp(i)*(phipt(jp)-phipt(j)))/rhV_T(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Rp(i)*dphi2t(jp)*Rp(i)*(phipt(jp)-phipt(j)))/rhV_T(i,j,k)				
				bbby(j)=1.-aaay(j)-cccy(j) 
			enddo
			aaay(0)=0.
			cccy(0)=0.
			bbby(0)=1.
			aaay(px*jmax)=0. !bc dVdt in j-dir is on px*jmax
			bbby(px*jmax)=1. 
			cccy(px*jmax)=0. 
			rhssy(0:px*jmax)=vv_T(i,0:px*jmax,k)
			CALL solve_tridiag_switchperiodic2(vv_T(i,0:px*jmax,k),aaay(0:px*jmax),bbby(0:px*jmax),
     &			cccy(0:px*jmax),rhssy(0:px*jmax),px*jmax+1,pery)
           enddo
		  enddo		  
		  do k=1,kmax/px 
           do i=1,imax
            do j=1,px*jmax !0,px*jmax+1
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= 0.25*(ans_T(i,j,k)+ans_T(i,jm,k)+ans_T(i,j,k+1)+ans_T(i,jm,k+1))
				ekm_plus=0.25*(ans_T(i,j,k)+ans_T(i,jp,k)+ans_T(i,j,k+1)+ans_T(i,jp,k+1))				
				aaay(j)=-CNz*ekm_min*dt/(Rp(i)*(phipt(j)-phipt(jm))*Rp(i)*dphi2t(j))/rhW_T(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Rp(i)*(phipt(jp)-phipt(j))*Rp(i)*dphi2t(j))/rhW_T(i,j,k)
				bbby(j)=1.-aaay(j)-cccy(j) 				
			enddo
			aaay(0)=0.
			aaay(px*jmax+1)=0.
			cccy(0)=0.
			cccy(px*jmax+1)=0.
			bbby(0)=1.
			bbby(px*jmax+1)=1.
			rhssy=ww_T(i,0:px*jmax+1,k)
			CALL solve_tridiag_switchperiodic2(ww_T(i,0:px*jmax+1,k),aaay,bbby,cccy,rhssy,px*jmax+2,pery)
           enddo
		  enddo
		ELSEIF (periodicy.eq.1) THEN !periodic lateral boundaries 
		  do k=1,kmax/px 
           do i=1,imax
            do j=1,px*jmax !0,px*jmax+1
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= 0.25*(ans_T(i,j,k)+ans_T(i,jm,k)+ans_T(i+1,j,k)+ans_T(i+1,jm,k))
				ekm_plus=0.25*(ans_T(i,j,k)+ans_T(i,jp,k)+ans_T(i+1,j,k)+ans_T(i+1,jp,k))
				aaay(j)=-CNz*ekm_min*dt/(Ru(i)*(phipt(j)-phipt(jm))*Ru(i)*dphi2t(j))/rhU_T(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Ru(i)*(phipt(jp)-phipt(j))*Ru(i)*dphi2t(j))/rhU_T(i,j,k)		
				bbby(j)=1.-aaay(j)-cccy(j) 
			enddo
			rhssy(1:px*jmax)=uu_T(i,1:px*jmax,k)
			CALL solve_tridiag_switchperiodic2(uu_T(i,1:px*jmax,k),aaay(1:px*jmax),bbby(1:px*jmax),cccy(1:px*jmax),
     &			rhssy(1:px*jmax),px*jmax,pery)
           enddo
		  enddo
		  do k=1,kmax/px 
           do i=1,imax
            do j=1,px*jmax !0,px*jmax+1
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= ans_T(i,j,k)
				ekm_plus=ans_T(i,jp,k)
				aaay(j)=-CNz*ekm_min*dt/(Rp(i)*dphi2t(j)*Rp(i)*(phipt(jp)-phipt(j)))/rhV_T(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Rp(i)*dphi2t(jp)*Rp(i)*(phipt(jp)-phipt(j)))/rhV_T(i,j,k)
				bbby(j)=1.-aaay(j)-cccy(j) 
			enddo
			rhssy(1:px*jmax)=vv_T(i,1:px*jmax,k)
			CALL solve_tridiag_switchperiodic2(vv_T(i,1:px*jmax,k),aaay(1:px*jmax),bbby(1:px*jmax),
     &			cccy(1:px*jmax),rhssy(1:px*jmax),px*jmax,pery)
           enddo
		  enddo		  
		  do k=1,kmax/px 
           do i=1,imax
            do j=1,px*jmax !0,px*jmax+1
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= 0.25*(ans_T(i,j,k)+ans_T(i,jm,k)+ans_T(i,j,k+1)+ans_T(i,jm,k+1))
				ekm_plus=0.25*(ans_T(i,j,k)+ans_T(i,jp,k)+ans_T(i,j,k+1)+ans_T(i,jp,k+1))
				aaay(j)=-CNz*ekm_min*dt/(Rp(i)*(phipt(j)-phipt(jm))*Rp(i)*dphi2t(j))/rhW_T(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Rp(i)*(phipt(jp)-phipt(j))*Rp(i)*dphi2t(j))/rhW_T(i,j,k)				
				bbby(j)=1.-aaay(j)-cccy(j) 	
			enddo
			rhssy(1:px*jmax)=ww_T(i,1:px*jmax,k)
			CALL solve_tridiag_switchperiodic2(ww_T(i,1:px*jmax,k),aaay(1:px*jmax),bbby(1:px*jmax),
     &			cccy(1:px*jmax),rhssy(1:px*jmax),px*jmax,pery)
           enddo
		  enddo		
		ELSEIF (periodicy.eq.2) THEN !free slip lateral boundaries dUWdn=0 at both ends and V=0
		  do k=1,kmax/px 
           do i=1,imax
            do j=1,px*jmax !0,px*jmax+1
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= 0.25*(ans_T(i,j,k)+ans_T(i,jm,k)+ans_T(i+1,j,k)+ans_T(i+1,jm,k))
				ekm_plus=0.25*(ans_T(i,j,k)+ans_T(i,jp,k)+ans_T(i+1,j,k)+ans_T(i+1,jp,k))
				aaay(j)=-CNz*ekm_min*dt/(Ru(i)*(phipt(j)-phipt(jm))*Ru(i)*dphi2t(j))/rhU_T(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Ru(i)*(phipt(jp)-phipt(j))*Ru(i)*dphi2t(j))/rhU_T(i,j,k)	
				bbby(j)=1.-aaay(j)-cccy(j) 
			enddo
			bbby(px*jmax)=bbby(px*jmax)+cccy(px*jmax) !d.dn=0
			cccy(px*jmax)=0.
			bbby(1)=bbby(1)+aaay(1) !d.dn=0
			aaay(1)=0.
			rhssy(1:px*jmax)=uu_T(i,1:px*jmax,k)
			CALL solve_tridiag_switchperiodic2(uu_T(i,1:px*jmax,k),aaay(1:px*jmax),bbby(1:px*jmax),
     &			cccy(1:px*jmax),rhssy(1:px*jmax),px*jmax,pery)			
           enddo
		  enddo
		  do k=1,kmax/px 
           do i=1,imax
            do j=1,px*jmax-1 !0,px*jmax+1
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= ans_T(i,j,k)
				ekm_plus=ans_T(i,jp,k)
				aaay(j)=-CNz*ekm_min*dt/(Rp(i)*dphi2t(j)*Rp(i)*(phipt(jp)-phipt(j)))/rhV_T(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Rp(i)*dphi2t(jp)*Rp(i)*(phipt(jp)-phipt(j)))/rhV_T(i,j,k)	
				bbby(j)=1.-aaay(j)-cccy(j) 
			enddo
			aaay(0)=0.
			cccy(0)=0.
			bbby(0)=1.
			aaay(px*jmax)=0. !bc dVdt in j-dir is on px*jmax
			bbby(px*jmax)=1. 
			cccy(px*jmax)=0. 
			rhssy(0:px*jmax)=vv_T(i,0:px*jmax,k)
			rhssy(0)=0.
			rhssy(px*jmax)=0. 			
			CALL solve_tridiag_switchperiodic2(vv_T(i,0:px*jmax,k),aaay(0:px*jmax),bbby(0:px*jmax),
     &			cccy(0:px*jmax),rhssy(0:px*jmax),px*jmax+1,pery)
           enddo
		  enddo		  
		  do k=1,kmax/px 
           do i=1,imax
            do j=1,px*jmax !0,px*jmax+1
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= 0.25*(ans_T(i,j,k)+ans_T(i,jm,k)+ans_T(i,j,k+1)+ans_T(i,jm,k+1))
				ekm_plus=0.25*(ans_T(i,j,k)+ans_T(i,jp,k)+ans_T(i,j,k+1)+ans_T(i,jp,k+1))
				aaay(j)=-CNz*ekm_min*dt/(Rp(i)*(phipt(j)-phipt(jm))*Rp(i)*dphi2t(j))/rhW_T(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Rp(i)*(phipt(jp)-phipt(j))*Rp(i)*dphi2t(j))/rhW_T(i,j,k)
				bbby(j)=1.-aaay(j)-cccy(j) 				
			enddo
			bbby(px*jmax)=bbby(px*jmax)+cccy(px*jmax) !d.dn=0
			cccy(px*jmax)=0.
			bbby(1)=bbby(1)+aaay(1) !d.dn=0
			aaay(1)=0.
			rhssy(1:px*jmax)=ww_T(i,1:px*jmax,k)
			CALL solve_tridiag_switchperiodic2(ww_T(i,1:px*jmax,k),aaay(1:px*jmax),bbby(1:px*jmax),
     &			cccy(1:px*jmax),rhssy(1:px*jmax),px*jmax,pery)	
           enddo
		  enddo
		ENDIF 
		
		call t2fp(uu_T(1:imax,1:jmax*px,1:kmax/px),dUdt(1:imax,1:jmax,1:kmax))
		call t2fp(vv_T(1:imax,1:jmax*px,1:kmax/px),dVdt(1:imax,1:jmax,1:kmax))
		call t2fp(ww_T(1:imax,1:jmax*px,1:kmax/px),dWdt(1:imax,1:jmax,1:kmax))	
		 
		return
		end	

	  
      subroutine diffuvw_zdir_CDS2_CNimpl

        USE nlist

        implicit none

        integer  im,ip,jm,jp,km,kp !,ib,ie,jb,je,kb,ke
        real CNz,CNx,CNy,dzi 
		real aaa(0:k1),bbb(0:k1),ccc(0:k1),ekm_min,ekm_plus,rhss(0:k1),rho_min,rho_plus
		integer perx,pery 


		IF (CNdiffz.eq.1) THEN 
			CNz=0.5 !0.55
		ELSEIF (CNdiffz.eq.11) THEN 
			CNz=CNdiff_factor !0.5 !0.55 
		ELSEIF (CNdiffz.eq.12) THEN 
			CNz=(3.-2.*CNdiff_factor)/3.
		ELSE 
			CNz=1.
		ENDIF 

		
!		 if (.false.) then 
!		 !if (CNdiffz.eq.11) then 
!			!Hirsch 2007 2nd edition Eq. 9.4.21 p.478 mentions that also explicit part of CN should build on for 3 subsequent directions
!			!here z-dir is done as 3th step on velocity field with finished CN explicit and implicit actions in x and y-dir:
!			call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 	Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
!			call diffuvw_zdir_CDS2_CNexpl(dnew,dnew2,dold,dUdt,dVdt,dWdt,ib,ie,jb,je,kb,ke)
!			dUdt(1:imax,1:jmax,1:kmax) = dUdt(1:imax,1:jmax,1:kmax) + dt*dnew (1:imax,1:jmax,1:kmax)/rhU(1:imax,1:jmax,1:kmax)  
!			dVdt(1:imax,1:jmax,1:kmax) = dVdt(1:imax,1:jmax,1:kmax) + dt*dnew2(1:imax,1:jmax,1:kmax)/rhV(1:imax,1:jmax,1:kmax)
!			dWdt(1:imax,1:jmax,1:kmax) = dWdt(1:imax,1:jmax,1:kmax) + dt*dold (1:imax,1:jmax,1:kmax)/rhW(1:imax,1:jmax,1:kmax)
!			call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),0,time_np,Ub1new,Vb1new,
!     & 	Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)	
!		 endif 		
	   IF (slip_bot.eq.-1) THEN !no slip bottom, free slip Neumann UV top 
		do j=1,jmax
         do i=1,imax
            do k=1,kmax !0,k1
				km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))
				rho_min=0.25*(drdt(i,j,k)+drdt(i,j,km)+drdt(i+1,j,k)+drdt(i+1,j,km))
				rho_plus=0.25*(drdt(i,j,k)+drdt(i,j,kp)+drdt(i+1,j,k)+drdt(i+1,j,kp))
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhU(i,j,k) !/rho_min
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhU(i,j,k) !/rho_plus
				bbb(k)=1.-aaa(k)-ccc(k) 				
            enddo
			bbb(1)=bbb(1)-aaa(1) 
			aaa(1)=0.
			bbb(kmax)=bbb(kmax)+ccc(kmax)
			ccc(kmax)=0. 
!!!			rhss(1:kmax)=dUdt2(i,j,1:kmax)
			rhss(1:kmax)=dUdt(i,j,1:kmax)
			CALL solve_tridiag(dUdt(i,j,1:kmax),aaa(1:kmax),bbb(1:kmax),ccc(1:kmax),rhss(1:kmax),kmax) 
         enddo
		enddo
		do j=1,jmax
         do i=1,imax
            do k=1,kmax !0,k1
				km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))		
				rho_min=0.25*(drdt(i,j,k)+drdt(i,j,km)+drdt(i,j+1,k)+drdt(i,j+1,km))
				rho_plus=0.25*(drdt(i,j,k)+drdt(i,j,kp)+drdt(i,j+1,k)+drdt(i,j+1,kp))
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhV(i,j,k) !/rho_min
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhV(i,j,k) !/rho_plus
				bbb(k)=1.-aaa(k)-ccc(k) 				
			enddo
			bbb(1)=bbb(1)-aaa(1) 
			aaa(1)=0.
			bbb(kmax)=bbb(kmax)+ccc(kmax)
			ccc(kmax)=0. 
!!!			rhss(1:kmax)=dVdt2(i,j,1:kmax)
			rhss(1:kmax)=dVdt(i,j,1:kmax)
			CALL solve_tridiag(dVdt(i,j,1:kmax),aaa(1:kmax),bbb(1:kmax),ccc(1:kmax),rhss(1:kmax),kmax)
         enddo
		enddo
		do j=1,jmax
         do i=1,imax
            do k=1,kmax-1 !0,k1
				!km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
				ekm_plus=ekm(i,j,kp)!*fc_global(i,j+rank*jmax,kp)
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhW(i,j,k) !/drdt(i,j,k)
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhW(i,j,k) !/drdt(i,j,kp)
				bbb(k)=1.-aaa(k)-ccc(k) 
            enddo
			aaa(0)=0.
			ccc(0)=0.
			bbb(0)=1.
			aaa(kmax)=0. !bc dWdt in k-dir is on kmax 
			bbb(kmax)=1.
			ccc(kmax)=0.
!!!			rhss(0:kmax)=dWdt2(i,j,0:kmax)
			rhss(0:kmax)=dWdt(i,j,0:kmax)
			CALL solve_tridiag(dWdt(i,j,0:kmax),aaa(0:kmax),bbb(0:kmax),ccc(0:kmax),rhss(0:kmax),k1) 
         enddo
		enddo
	  ELSEIF (slip_bot.eq.-2) THEN !no slip top and bottom 
		do j=1,jmax
         do i=1,imax
            do k=1,kmax !0,k1
				km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))
				rho_min=0.25*(drdt(i,j,k)+drdt(i,j,km)+drdt(i+1,j,k)+drdt(i+1,j,km))
				rho_plus=0.25*(drdt(i,j,k)+drdt(i,j,kp)+drdt(i+1,j,k)+drdt(i+1,j,kp))
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhU(i,j,k) !/rho_min
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhU(i,j,k) !/rho_plus
				bbb(k)=1.-aaa(k)-ccc(k) 				
            enddo
			bbb(1)=bbb(1)-aaa(1) 
			aaa(1)=0.
			bbb(kmax)=bbb(kmax)-ccc(kmax)
			ccc(kmax)=0. 
!!!			rhss(1:kmax)=dUdt2(i,j,1:kmax)
			rhss(1:kmax)=dUdt(i,j,1:kmax)
			CALL solve_tridiag(dUdt(i,j,1:kmax),aaa(1:kmax),bbb(1:kmax),ccc(1:kmax),rhss(1:kmax),kmax) 
         enddo
		enddo
		do j=1,jmax
         do i=1,imax
            do k=1,kmax !0,k1
				km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))		
				rho_min=0.25*(drdt(i,j,k)+drdt(i,j,km)+drdt(i,j+1,k)+drdt(i,j+1,km))
				rho_plus=0.25*(drdt(i,j,k)+drdt(i,j,kp)+drdt(i,j+1,k)+drdt(i,j+1,kp))
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhV(i,j,k) !/rho_min
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhV(i,j,k) !/rho_plus
				bbb(k)=1.-aaa(k)-ccc(k) 				
			enddo
			bbb(1)=bbb(1)-aaa(1) 
			aaa(1)=0.
			bbb(kmax)=bbb(kmax)-ccc(kmax)
			ccc(kmax)=0. 
!!!			rhss(1:kmax)=dVdt2(i,j,1:kmax)
			rhss(1:kmax)=dVdt(i,j,1:kmax)
			CALL solve_tridiag(dVdt(i,j,1:kmax),aaa(1:kmax),bbb(1:kmax),ccc(1:kmax),rhss(1:kmax),kmax)
         enddo
		enddo
		do j=1,jmax
         do i=1,imax
            do k=1,kmax-1 !0,k1
				!km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
				ekm_plus=ekm(i,j,kp)!*fc_global(i,j+rank*jmax,kp)
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhW(i,j,k) !/drdt(i,j,k)
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhW(i,j,k) !/drdt(i,j,kp)
				bbb(k)=1.-aaa(k)-ccc(k) 
            enddo
			aaa(0)=0.
			ccc(0)=0.
			bbb(0)=1.
			aaa(kmax)=0. !bc dWdt in k-dir is on kmax 
			bbb(kmax)=1.
			ccc(kmax)=0.
!!!			rhss(0:kmax)=dWdt2(i,j,0:kmax)
			rhss(0:kmax)=dWdt(i,j,0:kmax)
			CALL solve_tridiag(dWdt(i,j,0:kmax),aaa(0:kmax),bbb(0:kmax),ccc(0:kmax),rhss(0:kmax),k1) 
         enddo
		enddo	  
	  ELSE !Neuman UV top and bottom
		do j=1,jmax
         do i=1,imax
            do k=1,kmax !0,k1
				km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))
				rho_min=0.25*(drdt(i,j,k)+drdt(i,j,km)+drdt(i+1,j,k)+drdt(i+1,j,km))
				rho_plus=0.25*(drdt(i,j,k)+drdt(i,j,kp)+drdt(i+1,j,k)+drdt(i+1,j,kp))
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhU(i,j,k) !/rho_min
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhU(i,j,k) !/rho_plus
				bbb(k)=1.-aaa(k)-ccc(k) 				
            enddo
			bbb(1)=bbb(1)+aaa(1) 
			aaa(1)=0.
			bbb(kmax)=bbb(kmax)+ccc(kmax)
			ccc(kmax)=0. 
!!!			rhss(1:kmax)=dUdt2(i,j,1:kmax)
			rhss(1:kmax)=dUdt(i,j,1:kmax)
			CALL solve_tridiag(dUdt(i,j,1:kmax),aaa(1:kmax),bbb(1:kmax),ccc(1:kmax),rhss(1:kmax),kmax) 
         enddo
		enddo
		do j=1,jmax
         do i=1,imax
            do k=1,kmax !0,k1
				km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))		
				rho_min=0.25*(drdt(i,j,k)+drdt(i,j,km)+drdt(i,j+1,k)+drdt(i,j+1,km))
				rho_plus=0.25*(drdt(i,j,k)+drdt(i,j,kp)+drdt(i,j+1,k)+drdt(i,j+1,kp))
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhV(i,j,k) !/rho_min
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhV(i,j,k) !/rho_plus
				bbb(k)=1.-aaa(k)-ccc(k) 				
			enddo
			bbb(1)=bbb(1)+aaa(1) 
			aaa(1)=0.
			bbb(kmax)=bbb(kmax)+ccc(kmax)
			ccc(kmax)=0. 
!!!			rhss(1:kmax)=dVdt2(i,j,1:kmax)
			rhss(1:kmax)=dVdt(i,j,1:kmax)
			CALL solve_tridiag(dVdt(i,j,1:kmax),aaa(1:kmax),bbb(1:kmax),ccc(1:kmax),rhss(1:kmax),kmax)
         enddo
		enddo
		do j=1,jmax
         do i=1,imax
            do k=1,kmax-1 !0,k1
				!km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
				ekm_plus=ekm(i,j,kp)!*fc_global(i,j+rank*jmax,kp)
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhW(i,j,k) !/drdt(i,j,k)
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhW(i,j,k) !/drdt(i,j,kp)
				bbb(k)=1.-aaa(k)-ccc(k) 
            enddo
			aaa(0)=0.
			ccc(0)=0.
			bbb(0)=1.
			aaa(kmax)=0. !bc dWdt in k-dir is on kmax 
			bbb(kmax)=1.
			ccc(kmax)=0.
!!!			rhss(0:kmax)=dWdt2(i,j,0:kmax)
			rhss(0:kmax)=dWdt(i,j,0:kmax)
			CALL solve_tridiag(dWdt(i,j,0:kmax),aaa(0:kmax),bbb(0:kmax),ccc(0:kmax),rhss(0:kmax),k1) 
         enddo
		enddo
	  ENDIF 
		 
		return
		end			


      subroutine diffuvw_CDS2_3DCNimpl

        USE nlist

        implicit none

        integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
        real CNz
		real ekm_min,ekm_plus,rho_min,rho_plus
		real aaax(0:i1),bbbx(0:i1),cccx(0:i1)
		real aaay(0:j1),bbby(0:j1),cccy(0:j1)
		real aaa(0:k1),bbb(0:k1),ccc(0:k1)
		real rhssx(0:i1),rhssy(0:jmax*px+1),rhss(0:k1)
		real Ax3D(0:i1,0:j1,0:k1),Ay3D(0:i1,0:j1,0:k1),Az3D(0:i1,0:j1,0:k1)
		real Cx3D(0:i1,0:j1,0:k1),Cy3D(0:i1,0:j1,0:k1),Cz3D(0:i1,0:j1,0:k1)
		real DD3D(0:i1,0:j1,0:k1),RHS3D(0:i1,0:j1,0:k1)
		real u0(0:i1,0:j1,0:k1),v0(0:i1,0:j1,0:k1),w0(0:i1,0:j1,0:k1)
		
		CNz = CNdiff_factor
		! determine start condition iteration:
		IF(CNdiff_ini.eq.1) THEN
			u0 = dUdt 
			v0 = dVdt 
			w0 = dWdt 
		ELSEIF (CNdiff_ini.eq.2) THEN  
			u0 = Unew 
			v0 = Vnew 
			w0 = Wnew 
		ELSEIF (CNdiff_ini.eq.3) THEN  !start iteration with velocity field determined with explicit diffusion terms:
			ib=1 
			ie=imax 
			jb=1 
			je=jmax 
			kb=1 
			ke=kmax		
			! to reduce memory usage re-use of existing variables: Cx3d==dnew; Cy3d==dnew2; Cz3d==dold 
			call diffuvw_ydir_CDS2_CNexpl(Cx3D,Cy3D,Cz3D,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			u0 = dUdt + CNdiff_factor*dt*Cx3D/rhU  
			v0 = dVdt + CNdiff_factor*dt*Cy3D/rhV
			w0 = dWdt + CNdiff_factor*dt*Cz3D/rhW				
			call diffuvw_zdir_CDS2_CNexpl(Cx3D,Cy3D,Cz3D,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			u0 = u0 + CNdiff_factor*dt*Cx3D/rhU  
			v0 = v0 + CNdiff_factor*dt*Cy3D/rhV
			w0 = w0 + CNdiff_factor*dt*Cz3D/rhW
			call diffuvw_xdir_CDS2_CNexpl(Cx3D,Cy3D,Cz3D,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			u0 = u0 + CNdiff_factor*dt*Cx3D/rhU  
			v0 = v0 + CNdiff_factor*dt*Cy3D/rhV
			w0 = w0 + CNdiff_factor*dt*Cz3D/rhW		
			
			call bound_incljet(u0,v0,w0,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 		Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
		ELSEIF (CNdiff_ini.eq.4) THEN  !start iteration with velocity field determined with implicit ADI diffusion terms:
			Ax3D = dUdt !backup of dUdt 
			Ay3D = dVdt !backup of dVdt 
			Az3D = dWdt !backup of dWdt		
			! to reduce memory usage re-use of existing variables: Cx3d==dnew; Cy3d==dnew2; Cz3d==dold  
			! CNdiff_factor acts as theta in CN --> 0.5 gives 2nd order dt Douglas-Gunn ADI and 1 gives 1st order Douglass-Rachford; factor 0.5-1 gives blend
			! Douglas-Gunn and Douglass-Rachford ADI first take 100% explicit diffusion of 2 dirs not implicit
			ib=1 
			ie=imax 
			jb=1 
			je=jmax 
			kb=1 
			ke=kmax
			call diffuvw_ydir_CDS2_CNexpl(Cx3D,Cy3D,Cz3D,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			dUdt = dUdt + (1.-CNdiff_factor)*dt*Cx3D/rhU  
			dVdt = dVdt + (1.-CNdiff_factor)*dt*Cy3D/rhV
			dWdt = dWdt + (1.-CNdiff_factor)*dt*Cz3D/rhW				
			call diffuvw_zdir_CDS2_CNexpl(Cx3D,Cy3D,Cz3D,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			dUdt = dUdt + dt*Cx3D/rhU  
			dVdt = dVdt + dt*Cy3D/rhV
			dWdt = dWdt + dt*Cz3D/rhW
			call diffuvw_xdir_CDS2_CNexpl(Cx3D,Cy3D,Cz3D,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			dUdt = dUdt + dt*Cx3D/rhU  
			dVdt = dVdt + dt*Cy3D/rhV
			dWdt = dWdt + dt*Cz3D/rhW				

			call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)	
			CNdiffz=11 !get right theta inside subroutine:
			CALL diffuvw_ydir_CDS2_CNimpl
			CNdiffz=31
			! Douglas-Gunn ADI after each direction 50% implicit is finished subtract 50% explicit 
			! diffuvw_xdir_CDS2_CNexpl was executed as last direction; therefore Cx3D,Cy3D,Cz3D are already correct for xdir_expl:
			!call diffuvw_xdir_CDS2_CNexpl(Cx3D,Cy3D,Cz3D,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke) 
			dUdt = dUdt - CNdiff_factor*dt*Cx3D/rhU  
			dVdt = dVdt - CNdiff_factor*dt*Cy3D/rhV
			dWdt = dWdt - CNdiff_factor*dt*Cz3D/rhW
			call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)
			CNdiffz=11 !get right theta inside subroutine:
			CALL diffuvw_xdir_CDS2_CNimpl 
			CNdiffz=31 
			! Douglas-Gunn ADI after each direction 50% implicit is finished subtract 50% explicit 
			call diffuvw_zdir_CDS2_CNexpl(Cx3D,Cy3D,Cz3D,Unew,Vnew,Wnew,ib,ie,jb,je,kb,ke)
			dUdt = dUdt - CNdiff_factor*dt*Cx3D/rhU  
			dVdt = dVdt - CNdiff_factor*dt*Cy3D/rhV
			dWdt = dWdt - CNdiff_factor*dt*Cz3D/rhW			
			call bound_incljet(dUdt,dVdt,dWdt,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 			Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)								
			CNdiffz=11 !get right theta inside subroutine:
			CALL diffuvw_zdir_CDS2_CNimpl		
			CNdiffz=31 
			u0 = dUdt 
			v0 = dVdt 
			w0 = dWdt	
			call bound_incljet(u0,v0,w0,dRdt,MIN(0,slip_bot),time_np,Ub1new,Vb1new,
     & 		Wb1new,Ub2new,Vb2new,Wb2new,Ub3new,Vb3new,Wb3new)			
			dUdt = Ax3D !put back original dUdt 
			dVdt = Ay3D 
			dWdt = Az3D	 
		ELSE
			u0 = 0.*dUdt 
			v0 = 0.*dVdt 
			w0 = 0.*dWdt 		
		ENDIF
		

		!! Prof. Vuik mentions in the CFD2 lecture notes that enforcing boundary condition in PCG can be done naturally via ghost-cells also for the pre-conditioner 
		! for Dirichlet boundaries this approach is used (otherwise symmetry would get destroyed)
		! for Neumann type boundares the diagonal and off-diag terms in the matrix A are corrected which leaves a symmetric matrix


		!! Fill off-diags and diagonal forming matrix A and RHS3D to solve implicit CN diffusion terms of the form A*velocity=RHS3D 
		!! U-velocity:
		Ax3D=0. !off-diagonal i-1 in x-dir 
		Cx3D=0. !off-diagonal i+1 in x-dir 
		Ay3D=0. !off-diagonal j-1 in y-dir 
		Cy3D=0. !off-diagonal j+1 in y-dir 
		Az3D=0. !off-diagonal k-1 in z-dir 
		Cz3D=0. !off-diagonal k+1 in z-dir 
		DD3D=1. !diagonal for xyz-dir combined
		RHS3D=dUdt 
		IF (Apvisc_interp.eq.1.or.Apvisc_interp.eq.3) THEN ! use interpolation for apparent viscosity and for moleculair and eddy viscosity:
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax 
				!im=MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
				ekm_plus=ekm(ip,j,k)!*fc_global(ip,j+rank*jmax,k)
				aaax(i)=-CNz*ekm_min*Rp(i)*dt/(dr(i)*(Rp(ip)-Rp(i))*Ru(i))/rhU(i,j,k) !/drdt(i,j,k)
				cccx(i)=-CNz*ekm_plus*Rp(ip)*dt/(dr(ip)*(Rp(ip)-Rp(i))*Ru(i))/rhU(i,j,k) !/drdt(ip,j,k)
				bbbx(i)=1.-aaax(i)-cccx(i) 
            enddo
			Ax3D(1:imax,j,k)=Ax3D(1:imax,j,k)+aaax(1:imax)
			Cx3D(1:imax,j,k)=Cx3D(1:imax,j,k)+cccx(1:imax)
			DD3D(1:imax,j,k)=DD3D(1:imax,j,k)+bbbx(1:imax)-1.
           enddo
		  enddo
		  do k=1,kmax 
           do i=1,imax 
            do j=1,jmax 
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= 0.25*(ekm(i,j,k)+ekm(i,jm,k)+ekm(i+1,j,k)+ekm(i+1,jm,k))
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,jp,k)+ekm(i+1,j,k)+ekm(i+1,jp,k))
				aaay(j)=-CNz*ekm_min*dt/(Ru(i)*(phip(j)-phip(jm))*Ru(i)*dphi2(j))/rhU(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Ru(i)*(phip(jp)-phip(j))*Ru(i)*dphi2(j))/rhU(i,j,k)				
				bbby(j)=1.-aaay(j)-cccy(j) 
			enddo
			Ay3D(i,1:jmax,k)=Ay3D(i,1:jmax,k)+aaay(1:jmax)
			Cy3D(i,1:jmax,k)=Cy3D(i,1:jmax,k)+cccy(1:jmax)
			DD3D(i,1:jmax,k)=DD3D(i,1:jmax,k)+bbby(1:jmax)-1.
           enddo
		  enddo
		  do j=1,jmax
           do i=1,imax
            do k=1,kmax !0,k1
				km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhU(i,j,k) !/rho_min
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhU(i,j,k) !/rho_plus
				bbb(k)=1.-aaa(k)-ccc(k) 				
            enddo
			Az3D(i,j,1:kmax)=Az3D(i,j,1:kmax)+aaa(1:kmax)
			Cz3D(i,j,1:kmax)=Cz3D(i,j,1:kmax)+ccc(1:kmax)
			DD3D(i,j,1:kmax)=DD3D(i,j,1:kmax)+bbb(1:kmax)-1.
           enddo
		  enddo
		ELSEIF (Apvisc_interp.eq.2.or.Apvisc_interp.eq.4) THEN !use maximum of neighbouring cells for apparent visocity and linear interpolation for mol and turb visc
		  ekm = ekm - muA 
		  !! Fill off-diags and diagonal forming matrix A and RHS3D to solve implicit CN diffusion terms of the form A*velocity=RHS3D 
		  !! U-velocity:
		  Ax3D=0. !off-diagonal i-1 in x-dir 
		  Cx3D=0. !off-diagonal i+1 in x-dir 
		  Ay3D=0. !off-diagonal j-1 in y-dir 
		  Cy3D=0. !off-diagonal j+1 in y-dir 
		  Az3D=0. !off-diagonal k-1 in z-dir 
		  Cz3D=0. !off-diagonal k+1 in z-dir 
		  DD3D=1. !diagonal for xyz-dir combined
		  RHS3D=dUdt 
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax 
				!im=MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=ekm(i,j,k) + muA(i,j,k) !*fc_global(i,j+rank*jmax,k)
				ekm_plus=ekm(ip,j,k) + muA(ip,j,k) !*fc_global(ip,j+rank*jmax,k)
				aaax(i)=-CNz*ekm_min*Rp(i)*dt/(dr(i)*(Rp(ip)-Rp(i))*Ru(i))/rhU(i,j,k) !/drdt(i,j,k)
				cccx(i)=-CNz*ekm_plus*Rp(ip)*dt/(dr(ip)*(Rp(ip)-Rp(i))*Ru(i))/rhU(i,j,k) !/drdt(ip,j,k)
				bbbx(i)=1.-aaax(i)-cccx(i) 
            enddo
			Ax3D(1:imax,j,k)=Ax3D(1:imax,j,k)+aaax(1:imax)
			Cx3D(1:imax,j,k)=Cx3D(1:imax,j,k)+cccx(1:imax)
			DD3D(1:imax,j,k)=DD3D(1:imax,j,k)+bbbx(1:imax)-1.
           enddo
		  enddo
		  do k=1,kmax 
           do i=1,imax 
            do j=1,jmax 
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= 0.25*(ekm(i,j,k)+ekm(i,jm,k)+ekm(i+1,j,k)+ekm(i+1,jm,k))
     & 			+ MAX(muA(i,j,k),muA(i,jm,k),muA(i+1,j,k),muA(i+1,jm,k))				
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,jp,k)+ekm(i+1,j,k)+ekm(i+1,jp,k))
     &			+ MAX(muA(i,j,k),muA(i,jp,k),muA(i+1,j,k),muA(i+1,jp,k))				
				aaay(j)=-CNz*ekm_min*dt/(Ru(i)*(phip(j)-phip(jm))*Ru(i)*dphi2(j))/rhU(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Ru(i)*(phip(jp)-phip(j))*Ru(i)*dphi2(j))/rhU(i,j,k)				
				bbby(j)=1.-aaay(j)-cccy(j) 
			enddo
			Ay3D(i,1:jmax,k)=Ay3D(i,1:jmax,k)+aaay(1:jmax)
			Cy3D(i,1:jmax,k)=Cy3D(i,1:jmax,k)+cccy(1:jmax)
			DD3D(i,1:jmax,k)=DD3D(i,1:jmax,k)+bbby(1:jmax)-1.
           enddo
		  enddo
		  do j=1,jmax
           do i=1,imax
            do k=1,kmax !0,k1
				km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))!*MIN(fc_global(i,j+rank*jmax,k),
     &			+ MAX(muA(i,j,k),muA(i,j,km),muA(i+1,j,k),muA(i+1,j,km))				
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))!*MIN(fc_global(i,j+rank*jmax,k),
     &			+ MAX(muA(i,j,k),muA(i,j,kp),muA(i+1,j,k),muA(i+1,j,kp))				
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhU(i,j,k) !/rho_min
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhU(i,j,k) !/rho_plus
				bbb(k)=1.-aaa(k)-ccc(k) 				
            enddo
			Az3D(i,j,1:kmax)=Az3D(i,j,1:kmax)+aaa(1:kmax)
			Cz3D(i,j,1:kmax)=Cz3D(i,j,1:kmax)+ccc(1:kmax)
			DD3D(i,j,1:kmax)=DD3D(i,j,1:kmax)+bbb(1:kmax)-1.
           enddo
		  enddo		
		ENDIF

		  IF (periodicy.eq.0) THEN !bc defined at both lateral boundaries 
			jb=1
			je=jmax			
			IF (rank.eq.0) THEN 
				Ax3D(0:i1,0,0:k1)=0.
				Cx3D(0:i1,0,0:k1)=0.
				Ay3D(0:i1,0,0:k1)=0.
				Cy3D(0:i1,0,0:k1)=0.
				Az3D(0:i1,0,0:k1)=0.
				Cz3D(0:i1,0,0:k1)=0.				
				DD3D(0:i1,0,0:k1)=1.
				jb=1 !leave symmetry and reduce amount of work !0 
				je=jmax
			ELSEIF (rank.eq.px-1) THEN 
				Ax3D(0:i1,j1,0:k1)=0.
				Cx3D(0:i1,j1,0:k1)=0.
				Ay3D(0:i1,j1,0:k1)=0.
				Cy3D(0:i1,j1,0:k1)=0.				
				Az3D(0:i1,j1,0:k1)=0.
				Cz3D(0:i1,j1,0:k1)=0.
				DD3D(0:i1,j1,0:k1)=1.
				jb=1
				je=jmax !leave symmetry and reduce amount of work !j1 
			ENDIF
		  ELSEIF (periodicy.eq.2) THEN !d.dn=0
			jb=1
			je=jmax	
			IF (rank.eq.0) THEN 
				DD3D(0:i1,1,0:k1)=DD3D(0:i1,1,0:k1)+Ay3D(0:i1,1,0:k1) !1.-Cy3D(0:i1,1,0:k1) !d.dn=0
				Ay3D(0:i1,1,0:k1)=0.
				!Ax3D(0:i1,1,0:k1)=0.
				!Cx3D(0:i1,1,0:k1)=0.
				!Az3D(0:i1,1,0:k1)=0.
				!Cz3D(0:i1,1,0:k1)=0.				
			ELSEIF (rank.eq.px-1) THEN 
				DD3D(0:i1,jmax,0:k1)=DD3D(0:i1,jmax,0:k1)+Cy3D(0:i1,jmax,0:k1) !1.-Ay3D(0:i1,jmax,0:k1) !d.dn=0
				Cy3D(0:i1,jmax,0:k1)=0.
				!Ax3D(0:i1,jmax,0:k1)=0.
				!Cx3D(0:i1,jmax,0:k1)=0.
				!Az3D(0:i1,jmax,0:k1)=0.
				!Cz3D(0:i1,jmax,0:k1)=0.				
			ENDIF
		  ELSEIF (periodicy.eq.1) THEN 
			jb=1
			je=jmax			
		  ENDIF 
		IF (slip_bot.eq.-1) THEN !no slip bottom, free slip Neumann UV top 
			DD3D(0:i1,0:j1,kmax)=DD3D(0:i1,0:j1,kmax)+Cz3D(0:i1,0:j1,kmax) !1.-Az3D(0:i1,0:j1,kmax)
			Cz3D(0:i1,0:j1,kmax)=0.	
			!Ax3D(0:i1,0:j1,kmax)=0.
			!Cx3D(0:i1,0:j1,kmax)=0.
			!Ay3D(0:i1,0:j1,kmax)=0.
			!Cy3D(0:i1,0:j1,kmax)=0.	
			DD3D(0:i1,0:j1,1)=DD3D(0:i1,0:j1,1)-Az3D(0:i1,0:j1,1) !1.-2.*Az3D(0:i1,0:j1,1)-Cz3D(0:i1,0:j1,1)
			Az3D(0:i1,0:j1,1)=0.		
			!Cz3D(0:i1,0:j1,1)=0.
			!!RHS3D(0:i1,0:j1,1)=0.
			!Ax3D(0:i1,0:j1,1)=0.
			!Cx3D(0:i1,0:j1,1)=0.
			!Ay3D(0:i1,0:j1,1)=0.
			!Cy3D(0:i1,0:j1,1)=0.
			kb=1
			ke=kmax
		ELSEIF (slip_bot.eq.-2) THEN !no slip top and bottom 
			DD3D(0:i1,0:j1,kmax)=DD3D(0:i1,0:j1,kmax)-Cz3D(0:i1,0:j1,kmax) !1.-Az3D(0:i1,0:j1,kmax)-2.*Cz3D(0:i1,0:j1,kmax)
			Cz3D(0:i1,0:j1,kmax)=0.		
			!Az3D(0:i1,0:j1,kmax)=0.
			!!RHS3D(0:i1,0:j1,kmax)=0.
!			Ax3D(0:i1,0:j1,kmax)=0.
!			Cx3D(0:i1,0:j1,kmax)=0.
!			Ay3D(0:i1,0:j1,kmax)=0.
!			Cy3D(0:i1,0:j1,kmax)=0.			
			DD3D(0:i1,0:j1,1)=DD3D(0:i1,0:j1,1)-Az3D(0:i1,0:j1,1) !1.-2.*Az3D(0:i1,0:j1,1)-Cz3D(0:i1,0:j1,1)
			Az3D(0:i1,0:j1,1)=0.		
			!Cz3D(0:i1,0:j1,1)=0.
			!!RHS3D(0:i1,0:j1,1)=0.
!			Ax3D(0:i1,0:j1,1)=0.
!			Cx3D(0:i1,0:j1,1)=0.
!			Ay3D(0:i1,0:j1,1)=0.
!			Cy3D(0:i1,0:j1,1)=0.						
			kb=1
			ke=kmax			
		ELSE !Neuman UV top and bottom
			DD3D(0:i1,0:j1,kmax)=DD3D(0:i1,0:j1,kmax)+Cz3D(0:i1,0:j1,kmax) !1.-Az3D(0:i1,0:j1,kmax)
			Cz3D(0:i1,0:j1,kmax)=0.	
			!Ax3D(0:i1,0:j1,kmax)=0.
			!Cx3D(0:i1,0:j1,kmax)=0.
			!Ay3D(0:i1,0:j1,kmax)=0.
			!Cy3D(0:i1,0:j1,kmax)=0.				
			DD3D(0:i1,0:j1,1)=DD3D(0:i1,0:j1,1)+Az3D(0:i1,0:j1,1) !1.-Cz3D(0:i1,0:j1,1)
			Az3D(0:i1,0:j1,1)=0.		
			!Ax3D(0:i1,0:j1,1)=0.
			!Cx3D(0:i1,0:j1,1)=0.
			!Ay3D(0:i1,0:j1,1)=0.
			!Cy3D(0:i1,0:j1,1)=0.			
			kb=1
			ke=kmax			
		ENDIF 
		IF (periodicx.eq.0) THEN
			!RHS3D(0,0:j1,0:k1)=0.
			DD3D(0,0:j1,0:k1)=1.
			Ax3D(0,0:j1,0:k1)=0.
			Cx3D(0,0:j1,0:k1)=0.
			Ay3D(0,0:j1,0:k1)=0.
			Cy3D(0,0:j1,0:k1)=0.
			Az3D(0,0:j1,0:k1)=0.
			Cz3D(0,0:j1,0:k1)=0.
			DD3D(imax-1,0:j1,0:k1)=DD3D(imax-1,0:j1,0:k1)+Cx3D(imax-1,0:j1,0:k1) !1.-Ax3D(imax-1,0:j1,0:k1)
			!Ax3D(imax-1,0:j1,0:k1)=0.
			Cx3D(imax-1,0:j1,0:k1)=0.
			!Ay3D(imax-1,0:j1,0:k1)=0.
			!Cy3D(imax-1,0:j1,0:k1)=0.
			!Az3D(imax-1,0:j1,0:k1)=0.
			!Cz3D(imax-1,0:j1,0:k1)=0.
!			DD3D(imax,0:j1,0:k1)=1.-Ax3D(imax,0:j1,0:k1)
!			!Ax3D(imax,0:j1,0:k1)=0.
!			Cx3D(imax,0:j1,0:k1)=0.
!			Ay3D(imax,0:j1,0:k1)=0.
!			Cy3D(imax,0:j1,0:k1)=0.
!			Az3D(imax,0:j1,0:k1)=0.
!			Cz3D(imax,0:j1,0:k1)=0.
			ib=1 !leave symmetry and reduce amount of work !0
			ie=imax-1			
		ELSEIF (periodicx.eq.2) THEN
			!RHS3D(0,0:j1,0:k1)=0.
			DD3D(0,0:j1,0:k1)=1.
			Ax3D(0,0:j1,0:k1)=0.
			Cx3D(0,0:j1,0:k1)=0.
			Ay3D(0,0:j1,0:k1)=0.
			Cy3D(0,0:j1,0:k1)=0.
			Az3D(0,0:j1,0:k1)=0.
			Cz3D(0,0:j1,0:k1)=0.
			!RHS3D(imax,0:j1,0:k1)=0.
			DD3D(imax,0:j1,0:k1)=1.
			Ax3D(imax,0:j1,0:k1)=0.
			Cx3D(imax,0:j1,0:k1)=0.
			Ay3D(imax,0:j1,0:k1)=0.
			Cy3D(imax,0:j1,0:k1)=0.
			Az3D(imax,0:j1,0:k1)=0.
			Cz3D(imax,0:j1,0:k1)=0.		
			ib=1 !leave symmetry and reduce amount of work !0
			ie=imax-1 !leave symmetry and reduce amount of work !imax
		ELSEIF (periodicx.eq.1) THEN 
			ib=1
			ie=imax		
		ENDIF	
		
		IF (CNdiff_pc.eq.0) THEN !no pre-conditioner
			call CN3Dcg(dUdt,u0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank) 
			!call CN3Dcg2(dUdt,u0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol)
		ELSEIF (CNdiff_pc.eq.1) THEN !diag pre-conditioner 
			call CN3Dpcg_d(dUdt,u0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)
		ELSEIF (CNdiff_pc.eq.2) THEN !incomplete Cholesky pre-conditioner 
			call CN3Dpcg_ic(dUdt,u0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.3) THEN !ILU(0) pre-conditioner 
			call CN3Dpcg_ilu(dUdt,u0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.4) THEN !pol pc 
			call CN3Dpcg_pol(dUdt,u0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.11) THEN !diag pre-conditioner directly imposed 
			call CN3Dpcg_d2(dUdt,u0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.41) THEN !Pol pc without diagional scaling 
			call CN3Dpcg_pol2(dUdt,u0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ENDIF 
	
		!! Fill off-diags and diagonal forming matrix A and RHS3D to solve implicit CN diffusion terms of the form A*velocity=RHS3D 
		!! V-velocity:
		Ax3D=0. !off-diagonal i-1 in x-dir 
		Cx3D=0. !off-diagonal i+1 in x-dir 
		Ay3D=0. !off-diagonal j-1 in y-dir 
		Cy3D=0. !off-diagonal j+1 in y-dir 
		Az3D=0. !off-diagonal k-1 in z-dir 
		Cz3D=0. !off-diagonal k+1 in z-dir 
		DD3D=1. !diagonal for xyz-dir combined
		RHS3D=dVdt 
		IF (Apvisc_interp.eq.1.or.Apvisc_interp.eq.3) THEN ! use interpolation for apparent viscosity and for moleculair and eddy viscosity:   
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax !0,i1
				im=i-1 !MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j+1,k)+ekm(im,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(im,j+1+rank*jmax,k))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j+1,k)+ekm(ip,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(ip,j+1+rank*jmax,k))	
				aaax(i)=-CNz*ekm_min*dt*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))/rhV(i,j,k) !/rho_min
				cccx(i)=-CNz*ekm_plus*dt*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))/rhV(i,j,k) !/rho_plus
				bbbx(i)=1.-aaax(i)-cccx(i)
			enddo
			Ax3D(1:imax,j,k)=Ax3D(1:imax,j,k)+aaax(1:imax)
			Cx3D(1:imax,j,k)=Cx3D(1:imax,j,k)+cccx(1:imax)
			DD3D(1:imax,j,k)=DD3D(1:imax,j,k)+bbbx(1:imax)-1.	
           enddo
		  enddo
		  do k=1,kmax 
           do i=1,imax
            do j=1,jmax
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= ekm(i,j,k)
				ekm_plus=ekm(i,jp,k)
				aaay(j)=-CNz*ekm_min*dt/(Rp(i)*dphi2(j)*Rp(i)*(phip(jp)-phip(j)))/rhV(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Rp(i)*dphi2(jp)*Rp(i)*(phip(jp)-phip(j)))/rhV(i,j,k)
				bbby(j)=1.-aaay(j)-cccy(j) 
			enddo
			Ay3D(i,1:jmax,k)=Ay3D(i,1:jmax,k)+aaay(1:jmax)
			Cy3D(i,1:jmax,k)=Cy3D(i,1:jmax,k)+cccy(1:jmax)
			DD3D(i,1:jmax,k)=DD3D(i,1:jmax,k)+bbby(1:jmax)-1.
           enddo
		  enddo		  
		  do j=1,jmax
           do i=1,imax
            do k=1,kmax !0,k1
				km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))		
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhV(i,j,k) !/rho_min
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhV(i,j,k) !/rho_plus
				bbb(k)=1.-aaa(k)-ccc(k) 				
			enddo
			Az3D(i,j,1:kmax)=Az3D(i,j,1:kmax)+aaa(1:kmax)
			Cz3D(i,j,1:kmax)=Cz3D(i,j,1:kmax)+ccc(1:kmax)
			DD3D(i,j,1:kmax)=DD3D(i,j,1:kmax)+bbb(1:kmax)-1.
           enddo
		  enddo
		ELSEIF (Apvisc_interp.eq.2.or.Apvisc_interp.eq.4) THEN !use maximum of neighbouring cells for apparent visocity and linear interpolation for mol and turb visc
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax !0,i1
				im=i-1 !MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j+1,k)+ekm(im,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
     & 			+ MAX(muA(i,j,k),muA(im,j,k),muA(i,j+1,k),muA(im,j+1,k))				
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(im,j+1+rank*jmax,k))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j+1,k)+ekm(ip,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
     &			+ MAX(muA(i,j,k),muA(ip,j,k),muA(i,j+1,k),muA(ip,j+1,k))				
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(ip,j+1+rank*jmax,k))	
				aaax(i)=-CNz*ekm_min*dt*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))/rhV(i,j,k) !/rho_min
				cccx(i)=-CNz*ekm_plus*dt*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))/rhV(i,j,k) !/rho_plus
				bbbx(i)=1.-aaax(i)-cccx(i)
			enddo
			Ax3D(1:imax,j,k)=Ax3D(1:imax,j,k)+aaax(1:imax)
			Cx3D(1:imax,j,k)=Cx3D(1:imax,j,k)+cccx(1:imax)
			DD3D(1:imax,j,k)=DD3D(1:imax,j,k)+bbbx(1:imax)-1.	
           enddo
		  enddo
		  do k=1,kmax 
           do i=1,imax
            do j=1,jmax
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= ekm(i,j,k)+muA(i,j,k)
				ekm_plus=ekm(i,jp,k)+muA(i,jp,k)
				aaay(j)=-CNz*ekm_min*dt/(Rp(i)*dphi2(j)*Rp(i)*(phip(jp)-phip(j)))/rhV(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Rp(i)*dphi2(jp)*Rp(i)*(phip(jp)-phip(j)))/rhV(i,j,k)
				bbby(j)=1.-aaay(j)-cccy(j) 
			enddo
			Ay3D(i,1:jmax,k)=Ay3D(i,1:jmax,k)+aaay(1:jmax)
			Cy3D(i,1:jmax,k)=Cy3D(i,1:jmax,k)+cccy(1:jmax)
			DD3D(i,1:jmax,k)=DD3D(i,1:jmax,k)+bbby(1:jmax)-1.
           enddo
		  enddo		  
		  do j=1,jmax
           do i=1,imax
            do k=1,kmax !0,k1
				km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))!*MIN(fc_global(i,j+rank*jmax,k),
     &			+ MAX(muA(i,j,k),muA(i,j,km),muA(i,j+1,k),muA(i,j+1,km))				
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))!*MIN(fc_global(i,j+rank*jmax,k),
     &			+ MAX(muA(i,j,k),muA(i,j,kp),muA(i,j+1,k),muA(i,j+1,kp))				
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))		
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhV(i,j,k) !/rho_min
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhV(i,j,k) !/rho_plus
				bbb(k)=1.-aaa(k)-ccc(k) 				
			enddo
			Az3D(i,j,1:kmax)=Az3D(i,j,1:kmax)+aaa(1:kmax)
			Cz3D(i,j,1:kmax)=Cz3D(i,j,1:kmax)+ccc(1:kmax)
			DD3D(i,j,1:kmax)=DD3D(i,j,1:kmax)+bbb(1:kmax)-1.
           enddo
		  enddo		
		ENDIF 


			kb=1 
			ke=kmax 
		IF (slip_bot.eq.-1) THEN !no slip bottom, free slip Neumann UV top 
			DD3D(0:i1,0:j1,kmax)=DD3D(0:i1,0:j1,kmax)+Cz3D(0:i1,0:j1,kmax) !1.-Az3D(0:i1,0:j1,kmax)
			Cz3D(0:i1,0:j1,kmax)=0.	
			!Ax3D(0:i1,0:j1,kmax)=0.
			!Cx3D(0:i1,0:j1,kmax)=0.
			!Ay3D(0:i1,0:j1,kmax)=0.
			!Cy3D(0:i1,0:j1,kmax)=0.	
			DD3D(0:i1,0:j1,1)=DD3D(0:i1,0:j1,1)-Az3D(0:i1,0:j1,1) !1.-2.*Az3D(0:i1,0:j1,1)-Cz3D(0:i1,0:j1,1)
			Az3D(0:i1,0:j1,1)=0.		
			!Cz3D(0:i1,0:j1,1)=0.
			!!RHS3D(0:i1,0:j1,1)=0.
			!Ax3D(0:i1,0:j1,1)=0.
			!Cx3D(0:i1,0:j1,1)=0.
			!Ay3D(0:i1,0:j1,1)=0.
			!Cy3D(0:i1,0:j1,1)=0.
			kb=1 
			ke=kmax 
		ELSEIF (slip_bot.eq.-2) THEN !no slip top and bottom 
			DD3D(0:i1,0:j1,kmax)=DD3D(0:i1,0:j1,kmax)-Cz3D(0:i1,0:j1,kmax) !1.-Az3D(0:i1,0:j1,kmax)-2.*Cz3D(0:i1,0:j1,kmax)
			Cz3D(0:i1,0:j1,kmax)=0.		
			!Az3D(0:i1,0:j1,kmax)=0.
			!!RHS3D(0:i1,0:j1,kmax)=0.
!			Ax3D(0:i1,0:j1,kmax)=0.
!			Cx3D(0:i1,0:j1,kmax)=0.
!			Ay3D(0:i1,0:j1,kmax)=0.
!			Cy3D(0:i1,0:j1,kmax)=0.			
			DD3D(0:i1,0:j1,1)=DD3D(0:i1,0:j1,1)-Az3D(0:i1,0:j1,1) !1.-2.*Az3D(0:i1,0:j1,1)-Cz3D(0:i1,0:j1,1)
			Az3D(0:i1,0:j1,1)=0.		
			!Cz3D(0:i1,0:j1,1)=0.
			!!RHS3D(0:i1,0:j1,1)=0.
!			Ax3D(0:i1,0:j1,1)=0.
!			Cx3D(0:i1,0:j1,1)=0.
!			Ay3D(0:i1,0:j1,1)=0.
!			Cy3D(0:i1,0:j1,1)=0.						
			kb=1 
			ke=kmax 			
		ELSE !Neuman UV top and bottom
			DD3D(0:i1,0:j1,kmax)=DD3D(0:i1,0:j1,kmax)+Cz3D(0:i1,0:j1,kmax) !1.-Az3D(0:i1,0:j1,kmax)
			Cz3D(0:i1,0:j1,kmax)=0.	
			!Ax3D(0:i1,0:j1,kmax)=0.
			!Cx3D(0:i1,0:j1,kmax)=0.
			!Ay3D(0:i1,0:j1,kmax)=0.
			!Cy3D(0:i1,0:j1,kmax)=0.				
			DD3D(0:i1,0:j1,1)=DD3D(0:i1,0:j1,1)+Az3D(0:i1,0:j1,1) !1.-Cz3D(0:i1,0:j1,1)
			Az3D(0:i1,0:j1,1)=0.		
			!Ax3D(0:i1,0:j1,1)=0.
			!Cx3D(0:i1,0:j1,1)=0.
			!Ay3D(0:i1,0:j1,1)=0.
			!Cy3D(0:i1,0:j1,1)=0.			
			kb=1 
			ke=kmax 			
		ENDIF 

	  IF (periodicy.eq.0.or.periodicy.eq.2) THEN !bc defined at both lateral boundaries, which is V=0 for periodicy.eq.2; and periodicy.eq.1 do nothing	
		jb=1 
		je=jmax		  
		IF (rank.eq.0) THEN
			DD3D(0:i1,0,0:k1)=1.
			Ax3D(0:i1,0,0:k1)=0.
			Cx3D(0:i1,0,0:k1)=0.
			Ay3D(0:i1,0,0:k1)=0.
			Cy3D(0:i1,0,0:k1)=0.			
			Az3D(0:i1,0,0:k1)=0.
			Cz3D(0:i1,0,0:k1)=0.
			jb=1 !leave symmetry and reduce amount of work !0 
			je=jmax 
		ELSEIF (rank.eq.px-1) THEN
			DD3D(0:i1,jmax,0:k1)=1.
			Ax3D(0:i1,jmax,0:k1)=0.
			Cx3D(0:i1,jmax,0:k1)=0.
			Ay3D(0:i1,jmax,0:k1)=0.
			Cy3D(0:i1,jmax,0:k1)=0.			
			Az3D(0:i1,jmax,0:k1)=0.
			Cz3D(0:i1,jmax,0:k1)=0.
			jb=1 
			je=jmax-1 !leave symmetry and reduce amount of work !jmax			
		ENDIF 
	  ELSEIF (periodicy.eq.1) THEN  
		jb=1 
		je=jmax		  
	  ENDIF	

		IF (periodicx.eq.0) THEN
			!RHS3D(0,0:j1,0:k1)=0.
			DD3D(0,0:j1,0:k1)=1.
			Ax3D(0,0:j1,0:k1)=0.
			Cx3D(0,0:j1,0:k1)=0.
			Ay3D(0,0:j1,0:k1)=0.
			Cy3D(0,0:j1,0:k1)=0.
			Az3D(0,0:j1,0:k1)=0.
			Cz3D(0,0:j1,0:k1)=0.
			DD3D(imax,0:j1,0:k1)=DD3D(imax,0:j1,0:k1)+Cx3D(imax,0:j1,0:k1) !1.-Ax3D(imax,0:j1,0:k1)
			!Ax3D(imax,0:j1,0:k1)=0.
			Cx3D(imax,0:j1,0:k1)=0.
			!Ay3D(imax,0:j1,0:k1)=0.
			!Cy3D(imax,0:j1,0:k1)=0.
			!Az3D(imax,0:j1,0:k1)=0.
			!Cz3D(imax,0:j1,0:k1)=0.	
			ib=1 !leave symmetry and reduce amount of work !0 
			ie=imax 
		ELSEIF (periodicx.eq.2) THEN 
			DD3D(1,0:j1,0:k1)=DD3D(1,0:j1,0:k1)-Ax3D(1,0:j1,0:k1) !1.-2.*Ax3D(1,0:j1,0:k1)-Cx3D(1,0:j1,0:k1)
			Ax3D(1,0:j1,0:k1)=0.
!			!Cx3D(1,0:j1,0:k1)=0.
!			Ay3D(1,0:j1,0:k1)=0.
!			Cy3D(1,0:j1,0:k1)=0.
!			Az3D(1,0:j1,0:k1)=0.
!			Cz3D(1,0:j1,0:k1)=0.	
			DD3D(imax,0:j1,0:k1)=DD3D(imax,0:j1,0:k1)-Cx3D(imax,0:j1,0:k1) !1.-Ax3D(imax,0:j1,0:k1)-2.*Cx3D(imax,0:j1,0:k1)
!			!Ax3D(imax,0:j1,0:k1)=0.
			Cx3D(imax,0:j1,0:k1)=0.
!			Ay3D(imax,0:j1,0:k1)=0.
!			Cy3D(imax,0:j1,0:k1)=0.
!			Az3D(imax,0:j1,0:k1)=0.
!			Cz3D(imax,0:j1,0:k1)=0.
			ib=1 
			ie=imax	
		ELSEIF (periodicx.eq.1) THEN 
			ib=1 
			ie=imax			
		ENDIF
		
		IF (CNdiff_pc.eq.0) THEN !no pre-conditioner
			call CN3Dcg(dVdt,v0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank) 
			!call CN3Dcg2(dVdt,v0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol)
		ELSEIF (CNdiff_pc.eq.1) THEN !diag pre-conditioner 
			call CN3Dpcg_d(dVdt,v0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)
		ELSEIF (CNdiff_pc.eq.2) THEN !incomplete Cholesky pre-conditioner 
			call CN3Dpcg_ic(dVdt,v0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.3) THEN !ILU(0) pre-conditioner 
			call CN3Dpcg_ilu(dVdt,v0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.4) THEN !pol pc 
			call CN3Dpcg_pol(dVdt,v0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.11) THEN !diag pre-conditioner directly imposed 
			call CN3Dpcg_d2(dVdt,v0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.41) THEN !Pol pc without diagional scaling 
			call CN3Dpcg_pol2(dVdt,v0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ENDIF 
	  
	  

		!! Fill off-diags and diagonal forming matrix A and RHS3D to solve implicit CN diffusion terms of the form A*velocity=RHS3D 
		!! W-velocity:
		Ax3D=0. !off-diagonal i-1 in x-dir 
		Cx3D=0. !off-diagonal i+1 in x-dir 
		Ay3D=0. !off-diagonal j-1 in y-dir 
		Cy3D=0. !off-diagonal j+1 in y-dir 
		Az3D=0. !off-diagonal k-1 in z-dir 
		Cz3D=0. !off-diagonal k+1 in z-dir 
		DD3D=1. !diagonal for xyz-dir combined
		RHS3D=dWdt 
		
		IF (Apvisc_interp.eq.1.or.Apvisc_interp.eq.3) THEN !use linear interpolation of neighbouring cells for apparent visocity and mol and turb visc
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax !0,i1
				im=i-1 !MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j,k+1)+ekm(im,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(im,j+rank*jmax,k+1))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j,k+1)+ekm(ip,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(ip,j+rank*jmax,k+1))	
				aaax(i)=-CNz*ekm_min*dt*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))/rhW(i,j,k) !/rho_min
				cccx(i)=-CNz*ekm_plus*dt*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))/rhW(i,j,k) !/rho_plus
				bbbx(i)=1.-aaax(i)-cccx(i)				
            enddo
			Ax3D(1:imax,j,k)=Ax3D(1:imax,j,k)+aaax(1:imax)
			Cx3D(1:imax,j,k)=Cx3D(1:imax,j,k)+cccx(1:imax)
			DD3D(1:imax,j,k)=DD3D(1:imax,j,k)+bbbx(1:imax)-1. 	 
		   enddo
		  enddo	
		  do k=1,kmax 
           do i=1,imax
            do j=1,jmax 
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= 0.25*(ekm(i,j,k)+ekm(i,jm,k)+ekm(i,j,k+1)+ekm(i,jm,k+1))
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,jp,k)+ekm(i,j,k+1)+ekm(i,jp,k+1))
				aaay(j)=-CNz*ekm_min*dt/(Rp(i)*(phip(j)-phip(jm))*Rp(i)*dphi2(j))/rhW(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Rp(i)*(phip(jp)-phip(j))*Rp(i)*dphi2(j))/rhW(i,j,k)				
				bbby(j)=1.-aaay(j)-cccy(j) 	
			enddo
			Ay3D(i,1:jmax,k)=Ay3D(i,1:jmax,k)+aaay(1:jmax)
			Cy3D(i,1:jmax,k)=Cy3D(i,1:jmax,k)+cccy(1:jmax)
			DD3D(i,1:jmax,k)=DD3D(i,1:jmax,k)+bbby(1:jmax)-1.
           enddo
		  enddo		
		  ! all slip_bot values give same boundary condition for W:
		  do j=1,jmax
           do i=1,imax
            do k=1,kmax!-1 !0,k1
				!km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
				ekm_plus=ekm(i,j,kp)!*fc_global(i,j+rank*jmax,kp)
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhW(i,j,k) !/drdt(i,j,k)
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhW(i,j,k) !/drdt(i,j,kp)
				bbb(k)=1.-aaa(k)-ccc(k) 
            enddo
			Az3D(i,j,1:kmax)=Az3D(i,j,1:kmax)+aaa(1:kmax)
			Cz3D(i,j,1:kmax)=Cz3D(i,j,1:kmax)+ccc(1:kmax)
			DD3D(i,j,1:kmax)=DD3D(i,j,1:kmax)+bbb(1:kmax)-1.
           enddo
		  enddo		
		ELSEIF (Apvisc_interp.eq.2.or.Apvisc_interp.eq.4) THEN !use maximum of neighbouring cells for apparent visocity and linear interpolation for mol and turb visc
		  do k=1,kmax 
           do j=1,jmax
            do i=1,imax !0,i1
				im=i-1 !MAX(0,i-1)
				ip=i+1 !MIN(i1,i+1)
				ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j,k+1)+ekm(im,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
     & 			+ MAX(muA(i,j,k),muA(im,j,k),muA(i,j,k+1),muA(im,j,k+1))				
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(im,j+rank*jmax,k+1))		
				ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j,k+1)+ekm(ip,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
     &			+ MAX(muA(i,j,k),muA(ip,j,k),muA(i,j,k+1),muA(ip,j,k+1))				
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(ip,j+rank*jmax,k+1))	
				aaax(i)=-CNz*ekm_min*dt*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))/rhW(i,j,k) !/rho_min
				cccx(i)=-CNz*ekm_plus*dt*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))/rhW(i,j,k) !/rho_plus
				bbbx(i)=1.-aaax(i)-cccx(i)				
            enddo
			Ax3D(1:imax,j,k)=Ax3D(1:imax,j,k)+aaax(1:imax)
			Cx3D(1:imax,j,k)=Cx3D(1:imax,j,k)+cccx(1:imax)
			DD3D(1:imax,j,k)=DD3D(1:imax,j,k)+bbbx(1:imax)-1. 	 
		   enddo
		  enddo	
		  do k=1,kmax 
           do i=1,imax
            do j=1,jmax 
				jm=j-1 !MAX(0,j-1)
				jp=j+1 !MIN(px*jmax+1,j+1)
				ekm_min= 0.25*(ekm(i,j,k)+ekm(i,jm,k)+ekm(i,j,k+1)+ekm(i,jm,k+1))
     &			+ MAX(muA(i,j,k),muA(i,jm,k),muA(i,j,k+1),muA(i,jm,k+1))				
				ekm_plus=0.25*(ekm(i,j,k)+ekm(i,jp,k)+ekm(i,j,k+1)+ekm(i,jp,k+1))
     &			+ MAX(muA(i,j,k),muA(i,jp,k),muA(i,j,k+1),muA(i,jp,k+1))				
				aaay(j)=-CNz*ekm_min*dt/(Rp(i)*(phip(j)-phip(jm))*Rp(i)*dphi2(j))/rhW(i,j,k)
				cccy(j)=-CNz*ekm_plus*dt/(Rp(i)*(phip(jp)-phip(j))*Rp(i)*dphi2(j))/rhW(i,j,k)				
				bbby(j)=1.-aaay(j)-cccy(j) 	
			enddo
			Ay3D(i,1:jmax,k)=Ay3D(i,1:jmax,k)+aaay(1:jmax)
			Cy3D(i,1:jmax,k)=Cy3D(i,1:jmax,k)+cccy(1:jmax)
			DD3D(i,1:jmax,k)=DD3D(i,1:jmax,k)+bbby(1:jmax)-1.
           enddo
		  enddo		
		  ! all slip_bot values give same boundary condition for W:
		  do j=1,jmax
           do i=1,imax
            do k=1,kmax!-1 !0,k1
				!km=k-1 !MAX(0,k-1)
				kp=k+1 !MIN(k1,k+1)
				ekm_min=ekm(i,j,k)+muA(i,j,k) !*fc_global(i,j+rank*jmax,k)
				ekm_plus=ekm(i,j,kp)+muA(i,j,kp) !*fc_global(i,j+rank*jmax,kp)
				aaa(k)=-CNz*ekm_min*dt/dz**2/rhW(i,j,k) !/drdt(i,j,k)
				ccc(k)=-CNz*ekm_plus*dt/dz**2/rhW(i,j,k) !/drdt(i,j,kp)
				bbb(k)=1.-aaa(k)-ccc(k) 
            enddo
			Az3D(i,j,1:kmax)=Az3D(i,j,1:kmax)+aaa(1:kmax)
			Cz3D(i,j,1:kmax)=Cz3D(i,j,1:kmax)+ccc(1:kmax)
			DD3D(i,j,1:kmax)=DD3D(i,j,1:kmax)+bbb(1:kmax)-1.
           enddo
		  enddo
		ekm = ekm + muA 
		ENDIF 
		IF (periodicy.eq.0) THEN !bc defined at both lateral boundaries 
			jb=1
			je=jmax			
			IF (rank.eq.0) THEN 
				Ax3D(0:i1,0,0:k1)=0.
				Cx3D(0:i1,0,0:k1)=0.
				Ay3D(0:i1,0,0:k1)=0.
				Cy3D(0:i1,0,0:k1)=0.
				Az3D(0:i1,0,0:k1)=0.
				Cz3D(0:i1,0,0:k1)=0.				
				DD3D(0:i1,0,0:k1)=1.
				jb=1 !leave symmetry and reduce amount of work !0 
				je=jmax
			ELSEIF (rank.eq.px-1) THEN 
				Ax3D(0:i1,j1,0:k1)=0.
				Cx3D(0:i1,j1,0:k1)=0.
				Ay3D(0:i1,j1,0:k1)=0.
				Cy3D(0:i1,j1,0:k1)=0.				
				Az3D(0:i1,j1,0:k1)=0.
				Cz3D(0:i1,j1,0:k1)=0.
				DD3D(0:i1,j1,0:k1)=1.
				jb=1
				je=jmax !leave symmetry and reduce amount of work !j1 
			ENDIF
		ELSEIF (periodicy.eq.2) THEN !d.dn=0
			jb=1
			je=jmax			
			IF (rank.eq.0) THEN 
				DD3D(0:i1,1,0:k1)=DD3D(0:i1,1,0:k1)+Ay3D(0:i1,1,0:k1) !1.-Cy3D(0:i1,1,0:k1) !d.dn=0
				Ay3D(0:i1,1,0:k1)=0.
				!Ax3D(0:i1,1,0:k1)=0.
				!Cx3D(0:i1,1,0:k1)=0.
				!Az3D(0:i1,1,0:k1)=0.
				!Cz3D(0:i1,1,0:k1)=0.								
				jb=1
				je=jmax					
			ELSEIF (rank.eq.px-1) THEN 
				DD3D(0:i1,jmax,0:k1)=DD3D(0:i1,jmax,0:k1)+Cy3D(0:i1,jmax,0:k1) !1.-Ay3D(0:i1,jmax,0:k1) !DD3D(0:i1,jmax,0:k1)+Cy3D(0:i1,jmax,0:k1) !d.dn=0
				Cy3D(0:i1,jmax,0:k1)=0.
				!Ax3D(0:i1,jmax,0:k1)=0.
				!Cx3D(0:i1,jmax,0:k1)=0.
				!Az3D(0:i1,jmax,0:k1)=0.
				!Cz3D(0:i1,jmax,0:k1)=0.				
				jb=1
				je=jmax					
			ENDIF
		ELSEIF (periodicy.eq.1) THEN 
			jb=1
			je=jmax			
		ENDIF !periodicy.eq.1 to nothing

		IF (periodicx.eq.0) THEN
			DD3D(0,0:j1,0:k1)=1.
			Ax3D(0,0:j1,0:k1)=0.
			Cx3D(0,0:j1,0:k1)=0.
			Ay3D(0,0:j1,0:k1)=0.
			Cy3D(0,0:j1,0:k1)=0.
			Az3D(0,0:j1,0:k1)=0.
			Cz3D(0,0:j1,0:k1)=0.
			DD3D(imax,0:j1,0:k1)=DD3D(imax,0:j1,0:k1)+Cx3D(imax,0:j1,0:k1) !1.-Ax3D(imax,0:j1,0:k1)
			!Ax3D(imax,0:j1,0:k1)=0.
			Cx3D(imax,0:j1,0:k1)=0.
!			Ay3D(imax,0:j1,0:k1)=0.
!			Cy3D(imax,0:j1,0:k1)=0.
!			Az3D(imax,0:j1,0:k1)=0.
!			Cz3D(imax,0:j1,0:k1)=0.	
			ib = 1 !leave symmetry and reduce amount of work !0 
			ie = imax 
		ELSEIF (periodicx.eq.2) THEN 
			DD3D(1,0:j1,0:k1)=DD3D(1,0:j1,0:k1)-Ax3D(1,0:j1,0:k1) !1.-2.*Ax3D(1,0:j1,0:k1)-Cx3D(1,0:j1,0:k1)
			Ax3D(1,0:j1,0:k1)=0.
!			!Cx3D(1,0:j1,0:k1)=0.
!			Ay3D(1,0:j1,0:k1)=0.
!			Cy3D(1,0:j1,0:k1)=0.
!			Az3D(1,0:j1,0:k1)=0.
!			Cz3D(1,0:j1,0:k1)=0.		
		!this is needed, when done inside loop DD3D still contains contributions from z-dir and y-dir and this gives wrong bc at imax
			RHS3D(imax,0:j1,0:k1)=dWdt(imax,0:j1,0:k1)-2.*Cx3D(imax,0:j1,0:k1)*W_ox 
!			!RHS3D(imax,0:j1,0:k1)=W_ox-2.*Cx3D(imax,0:j1,0:k1)*W_ox 
!			DD3D(imax,0:j1,0:k1)=1.-Ax3D(imax,0:j1,0:k1)-2.*Cx3D(imax,0:j1,0:k1)!-Az3D(imax,0:j1,0:k1)-Cz3D(imax,0:j1,0:k1)  
			DD3D(imax,0:j1,0:k1)=DD3D(imax,0:j1,0:k1)-Cx3D(imax,0:j1,0:k1) !+Ay3D(imax,0:j1,0:k1)+Cy3D(imax,0:j1,0:k1)
!			!Ax3D(imax,0:j1,0:k1)=0.
			Cx3D(imax,0:j1,0:k1)=0.			
!			Ay3D(imax,0:j1,0:k1)=0. !have to switch off diffusion in y-dir and z-dir for this boundary condition; otherwise bizar near wall W velocity sdc for high viscosity?
!			Cy3D(imax,0:j1,0:k1)=0.
!			Az3D(imax,0:j1,0:k1)=0.
!			Cz3D(imax,0:j1,0:k1)=0.
!			!!DD3D(imax,0:j1,1)=DD3D(imax,0:j1,1)+Az3D(imax,0:j1,1)
!			!!Az3D(imax,0:j1,1)=0.
!			!!DD3D(imax,0:j1,kmax-1)=DD3D(imax,0:j1,kmax-1)+Cz3D(imax,0:j1,kmax-1)
!			!!Cz3D(imax,0:j1,kmax-1)=0.

!			DD3D(imax,0:j1,0:k1)=DD3D(imax,0:j1,0:k1)-Cx3D(imax,0:j1,0:k1)  
!			RHS3D(imax,0:j1,0:k1)=dWdt(imax,0:j1,0:k1)-2.*Cx3D(imax,0:j1,0:k1)*W_ox
!			Cx3D(imax,0:j1,0:k1)=0.
			
!			RHS3D(imax,0:j1,0:k1)=W_ox-0.5*(Wnew(imax,0:j1,0:k1)-Wnew(imax-1,0:j1,0:k1)) !approximate bc which is stable for a not so important test case 
!			DD3D(imax,0:j1,0:k1)=1.
!			Ax3D(imax,0:j1,0:k1)=0.
!			Cx3D(imax,0:j1,0:k1)=0.			
!			Ay3D(imax,0:j1,0:k1)=0. !have to switch off diffusion in y-dir and z-dir for this boundary condition; otherwise bizar near wall W velocity sdc for high viscosity?
!			Cy3D(imax,0:j1,0:k1)=0.
!			Az3D(imax,0:j1,0:k1)=0.
!			Cz3D(imax,0:j1,0:k1)=0.
			
			
			ib = 1 
			ie = imax	
		ELSEIF (periodicx.eq.1) THEN 
			ib = 1 
			ie = imax		
		ENDIF	
		! for all slip_bot values W=0 top and bottom:
		DD3D(0:i1,0:j1,0)=1. 
		Ax3D(0:i1,0:j1,0)=0.
		Cx3D(0:i1,0:j1,0)=0.			
		Ay3D(0:i1,0:j1,0)=0.
		Cy3D(0:i1,0:j1,0)=0.
		Az3D(0:i1,0:j1,0)=0.
		Cz3D(0:i1,0:j1,0)=0.
		DD3D(0:i1,0:j1,kmax)=1. 
		Ax3D(0:i1,0:j1,kmax)=0.
		Cx3D(0:i1,0:j1,kmax)=0.			
		Ay3D(0:i1,0:j1,kmax)=0.
		Cy3D(0:i1,0:j1,kmax)=0.
		Az3D(0:i1,0:j1,kmax)=0.
		Cz3D(0:i1,0:j1,kmax)=0.
		kb = 1 !leave symmetry and reduce amount of work ! 0 
		ke = kmax-1 !leave symmetry and reduce amount of work !kmax 		
		
		IF (CNdiff_pc.eq.0) THEN !no pre-conditioner
			call CN3Dcg(dWdt,w0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank) 
			!call CN3Dcg2(dWdt,w0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol)
		ELSEIF (CNdiff_pc.eq.1) THEN !diag pre-conditioner 
			call CN3Dpcg_d(dWdt,w0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)
		ELSEIF (CNdiff_pc.eq.2) THEN !incomplete Cholesky pre-conditioner 
			call CN3Dpcg_ic(dWdt,w0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.3) THEN !ILU(0) pre-conditioner 
			call CN3Dpcg_ilu(dWdt,w0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.4) THEN !pol pc 
			call CN3Dpcg_pol(dWdt,w0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.11) THEN !diag pre-conditioner directly imposed 
			call CN3Dpcg_d2(dWdt,w0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ELSEIF (CNdiff_pc.eq.41) THEN !Pol pc without diagional scaling 
			call CN3Dpcg_pol2(dWdt,w0,DD3D,Ax3D,Cx3D,Ay3D,Cy3D,Az3D,Cz3D,RHS3D,i1,j1,k1,ib,ie,jb,je,kb,ke,CNdiff_maxi,CNdiff_tol,rank)			
		ENDIF  
		
		return
		end	
		



      subroutine diffuvw_xdir_CDS2_CNexpl(putoutu,putoutv,putoutw,Uvel,Vvel,Wvel,ib,ie,jb,je,kb,ke)

      USE nlist

      implicit none
      
!       include 'param.txt'
!       include 'common.txt'


c
c*****************************************************************
c
c      calculates the diffusion of uvw-velocities in r direction 
c      on input :
c
c          putoutu,v,w       : empty on input
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
c          putoutu,v,w       : diffusion part in r direction on u,v,w
c          other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
      real     putoutu(0:i1,0:j1,0:k1),putoutv(0:i1,0:j1,0:k1),putoutw(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         aaax,bbbx,cccx,ekm_min,ekm_plus
      real CNz,CNx,CNy,dzi 

	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSEIF (CNdiffz.eq.2.or.CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNz=0. 
	ELSE
	  CNz=1.
	ENDIF
	IF (CNdiffz.eq.11) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNx=1. !0.5 !0.5*2. !0.45 !0.5
	  CNy=1. !0.5 !0.5*2. !0.45 !0.5
	  CNz=1. !0.5 !0.5*2.
	ELSEIF (CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNx=0. 
	  CNy=0. 
	ELSE
	  CNx=1.
	  CNy=1. 
	ENDIF	
	IF (CNdiffz.eq.31) THEN 
	  CNx=1.
	  CNy=1. 
	  CNz=1. 
	ENDIF 	
      dzi =1./dz
      do i=ib,ie
      ip=i+1
      im=i-1
	  do j=jb,je
	  jp=j+1
	  jm=j-1
		do k=kb,ke !k=MAX(kb,kbed(i,k)),ke !kb,ke
			kp=k+1
			km=k-1
			! not *dt because this is done in solve.f 
			! d^2udx^2 
			ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
			ekm_plus=ekm(ip,j,k)!*fc_global(ip,j+rank*jmax,k)
			aaax=CNx*ekm_min*Rp(i)/(dr(i)*(Rp(ip)-Rp(i))*Ru(i))
			cccx=CNx*ekm_plus*Rp(ip)/(dr(ip)*(Rp(ip)-Rp(i))*Ru(i))
			bbbx=-aaax-cccx
			putoutu(i,j,k) = aaax*Uvel(im,j,k) + bbbx*Uvel(i,j,k) + cccx*Uvel(ip,j,k)
			! d^2vdx^2 
			ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j+1,k)+ekm(im,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(im,j+1+rank*jmax,k))		
			ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j+1,k)+ekm(ip,j+1,k))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+1+rank*jmax,k),fc_global(ip,j+1+rank*jmax,k))	
			aaax=CNx*ekm_min*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))
			cccx=CNx*ekm_plus*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))
			bbbx=-aaax-cccx
			putoutv(i,j,k) = aaax*Vvel(im,j,k) + bbbx*Vvel(i,j,k) + cccx*Vvel(ip,j,k)
			! d^2wdx^2 
			ekm_min=0.25*(ekm(i,j,k)+ekm(im,j,k)+ekm(i,j,k+1)+ekm(im,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
	!     & 			fc_global(im,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(im,j+rank*jmax,k+1))		
			ekm_plus=0.25*(ekm(i,j,k)+ekm(ip,j,k)+ekm(i,j,k+1)+ekm(ip,j,k+1))!*MIN(fc_global(i,j+rank*jmax,k),
	!     & 			fc_global(ip,j+rank*jmax,k),fc_global(i,j+rank*jmax,k+1),fc_global(ip,j+rank*jmax,k+1))	
			aaax=CNx*ekm_min*Ru(im)/((Rp(i)-Rp(im))*dr(i)*Rp(i))
			cccx=CNx*ekm_plus*Ru(i)/((Rp(ip)-Rp(i))*dr(i)*Rp(i))
			bbbx=-aaax-cccx
			putoutw(i,j,k) = aaax*Wvel(im,j,k) + bbbx*Wvel(i,j,k) + cccx*Wvel(ip,j,k)
		  enddo
        enddo
      enddo
c
      return
      end	
	  
      subroutine diffuvw_ydir_CDS2_CNexpl(putoutu,putoutv,putoutw,Uvel,Vvel,Wvel,ib,ie,jb,je,kb,ke)

      USE nlist

      implicit none
      
!       include 'param.txt'
!       include 'common.txt'


c
c*****************************************************************
c
c      calculates the diffusion of uvw-velocities in phi (or called y) direction 
c      on input :
c
c          putoutu,v,w       : empty on input
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
c          putoutu,v,w       : diffusion part in phi direction on u,v,w
c          other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
      real     putoutu(0:i1,0:j1,0:k1),putoutv(0:i1,0:j1,0:k1),putoutw(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         aaay,bbby,cccy,ekm_min,ekm_plus
      real CNz,CNx,CNy,dzi 

	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSEIF (CNdiffz.eq.2.or.CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNz=0. 
	ELSE
	  CNz=1.
	ENDIF
	IF (CNdiffz.eq.11) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNx=1. !0.5 !0.5*2. !0.45 !0.5
	  CNy=1. !0.5 !0.5*2. !0.45 !0.5
	  CNz=1. !0.5 !0.5*2.
	ELSEIF (CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNx=0. 
	  CNy=0. 
	ELSE
	  CNx=1.
	  CNy=1. 
	ENDIF	
	IF (CNdiffz.eq.31) THEN 
	  CNx=1.
	  CNy=1.
	  CNz=1.
	ENDIF 
      dzi =1./dz
      do i=ib,ie
      ip=i+1
      im=i-1
	  do j=jb,je
	  jp=j+1
	  jm=j-1
		do k=kb,ke !k=MAX(kb,kbed(i,k)),ke !kb,ke
			kp=k+1
			km=k-1
			! not *dt because this is done in solve.f 
			! d^2udy^2
			ekm_min= 0.25*(ekm(i,j,k)+ekm(i,jm,k)+ekm(i+1,j,k)+ekm(i+1,jm,k))
			ekm_plus=0.25*(ekm(i,j,k)+ekm(i,jp,k)+ekm(i+1,j,k)+ekm(i+1,jp,k))
			aaay=CNy*ekm_min/(Ru(i)*(phip(j)-phip(jm))*Ru(i)*dphi2(j))
			cccy=CNy*ekm_plus/(Ru(i)*(phip(jp)-phip(j))*Ru(i)*dphi2(j))
			bbby=-aaay-cccy	
			putoutu(i,j,k) = aaay*Uvel(i,jm,k) + bbby*Uvel(i,j,k) + cccy*Uvel(i,jp,k)
			! d^2vdy^2
			ekm_min= ekm(i,j,k)
			ekm_plus=ekm(i,jp,k)
			aaay=CNy*ekm_min/(Rp(i)*dphi2(j)*Rp(i)*(phip(jp)-phip(j)))
			cccy=CNy*ekm_plus/(Rp(i)*dphi2(jp)*Rp(i)*(phip(jp)-phip(j)))
			bbby=-aaay-cccy	
			putoutv(i,j,k) = aaay*Vvel(i,jm,k) + bbby*Vvel(i,j,k) + cccy*Vvel(i,jp,k)
			! d^2wdy^2
			ekm_min= 0.25*(ekm(i,j,k)+ekm(i,jm,k)+ekm(i,j,k+1)+ekm(i,jm,k+1))
			ekm_plus=0.25*(ekm(i,j,k)+ekm(i,jp,k)+ekm(i,j,k+1)+ekm(i,jp,k+1))
			aaay=CNy*ekm_min/(Rp(i)*(phip(j)-phip(jm))*Rp(i)*dphi2(j))
			cccy=CNy*ekm_plus/(Rp(i)*(phip(jp)-phip(j))*Rp(i)*dphi2(j))
			bbby=-aaay-cccy	
			putoutw(i,j,k) = aaay*Wvel(i,jm,k) + bbby*Wvel(i,j,k) + cccy*Wvel(i,jp,k)
		  enddo
        enddo
      enddo
c
      return
      end	
	  
	  
      subroutine diffuvw_zdir_CDS2_CNexpl(putoutu,putoutv,putoutw,Uvel,Vvel,Wvel,ib,ie,jb,je,kb,ke)

      USE nlist

      implicit none
      
!       include 'param.txt'
!       include 'common.txt'


c
c*****************************************************************
c
c      calculates the diffusion of uvw-velocities in z direction 
c      on input :
c
c          putoutu,v,w       : empty on input
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
c          putoutu,v,w       : diffusion part in z direction on u,v,w
c          other parameters  : all unchanged
c
c*****************************************************************
      integer  im,ip,jm,jp,km,kp,ib,ie,jb,je,kb,ke
      real     putoutu(0:i1,0:j1,0:k1),putoutv(0:i1,0:j1,0:k1),putoutw(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         aaaz,bbbz,cccz,ekm_min,ekm_plus
      real CNz,CNx,CNy,dzi 

	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSEIF (CNdiffz.eq.2.or.CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNz=0. 
	ELSE
	  CNz=1.
	ENDIF
	IF (CNdiffz.eq.11) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNx=1. !0.5 !0.5*2. !0.45 !0.5
	  CNy=1. !0.5 !0.5*2. !0.45 !0.5
	  CNz=1. !0.5 !0.5*2.
	ELSEIF (CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNx=0. 
	  CNy=0. 
	ELSE
	  CNx=1.
	  CNy=1. 
	ENDIF	
	IF (CNdiffz.eq.31) THEN 
	  CNx=1.
	  CNy=1. 
	  CNz=1. 
	ENDIF 	
      dzi =1./dz
      do i=ib,ie
      ip=i+1
      im=i-1
	  do j=jb,je
	  jp=j+1
	  jm=j-1
		do k=kb,ke !k=MAX(kb,kbed(i,k)),ke !kb,ke
			kp=k+1
			km=k-1
			! not *dt because this is done in solve.f 
			! d^2udz^2 
			ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i+1,j,k)+ekm(i+1,j,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,km))
			ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i+1,j,k)+ekm(i+1,j,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i+1,j+rank*jmax,k),fc_global(i+1,j+rank*jmax,kp))
			aaaz=CNz*ekm_min/dz**2
			cccz=CNz*ekm_plus/dz**2
			bbbz=-aaaz-cccz
			putoutu(i,j,k) = aaaz*Uvel(i,j,km) + bbbz*Uvel(i,j,k) + cccz*Uvel(i,j,kp)
			! d^2vdz^2 
			ekm_min=0.25*(ekm(i,j,k)+ekm(i,j,km)+ekm(i,j+1,k)+ekm(i,j+1,km))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,km),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,km))		
			ekm_plus=0.25*(ekm(i,j,k)+ekm(i,j,kp)+ekm(i,j+1,k)+ekm(i,j+1,kp))!*MIN(fc_global(i,j+rank*jmax,k),
!     & 			fc_global(i,j+rank*jmax,kp),fc_global(i,j+1+rank*jmax,k),fc_global(i,j+1+rank*jmax,kp))		
			aaaz=CNz*ekm_min/dz**2
			cccz=CNz*ekm_plus/dz**2
			bbbz=-aaaz-cccz
			putoutv(i,j,k) = aaaz*Vvel(i,j,km) + bbbz*Vvel(i,j,k) + cccz*Vvel(i,j,kp)
			! d^2wdz^2 
			ekm_min=ekm(i,j,k)!*fc_global(i,j+rank*jmax,k)
			ekm_plus=ekm(i,j,kp)!*fc_global(i,j+rank*jmax,kp)
			aaaz=CNz*ekm_min/dz**2
			cccz=CNz*ekm_plus/dz**2
			bbbz=-aaaz-cccz 
			putoutw(i,j,k) = aaaz*Wvel(i,j,km) + bbbz*Wvel(i,j,k) + cccz*Wvel(i,j,kp)
		  enddo
        enddo
      enddo
c
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
      real	dzdz_i,Rpdr_i,Rpdphi2_i,CNz,CNx,CNy


	IF (CNdiffz.eq.1) THEN !CN diff in z-dir is half old, half new timestep (this is half old timestep)
	  CNz=0.5
	ELSEIF (CNdiffz.eq.2.or.CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNz=0. 
	ELSE
	  CNz=1.
	ENDIF
	IF (CNdiffz.eq.11) THEN !CNdiffz 11 is 0.5*DIFF(C^n+1-C^n) + 1*DIFF(C^n) --> this is explicit last step
	  CNx=0.5 !0.5*2. !0.45 !0.5
	  CNy=0.5 !0.5*2. !0.45 !0.5
	  CNz=0.5 !0.5*2.	  
	ELSEIF (CNdiffz.eq.12) THEN !CN diff in z-dir is 100% new timestep (Euler backward)
	  CNx=0. 
	  CNy=0. 
	ELSE
	  CNx=1.
	  CNy=1. 
	ENDIF	
      putin2=putin
      do k=k1,k1 !-kjet,k1
	do t=1,tmax_inPpunt
	  i=i_inPpunt(t)
	  j=j_inPpunt(t)
	  putin2(i,j,k) = putin2(i,j,kmax) ! No diffusion over vertical inflow boundary, this line is needed for exact influx
	enddo
      enddo

	  
	 DO i=0,i1 ! not needed anymore because Diffcoff is made 0 around kbed
	  DO j=0,j1
	  	putin2(i,j,kbed(i,j))=putin2(i,j,kbed(i,j)+1) ! neumann boundary to prevent vertical diffusion into bed
	  ENDDO
	 ENDDO

	
	  
      dzdz_i=1./(dz*dz)
c
c     -------------------------------------------start i-loop
      do 100 i=ib,ie
c
      ip=i+1
      im=i-1
      Rpdr_i=1./(Rp(i)*dr(i))
      
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
     1 ( CNx*Ru(i)  * (ekh(i,j,k)+ekh(ip,j,k)) * MIN(fc_global(i,j+rank*jmax,k),fc_global(ip,j+rank*jmax,k)) *
     1            (putin2(ip,j,k)-putin2(i,j,k) ) / (Rp(ip) - Rp(i))   -
     1   CNx*Ru(im) * (ekh(i,j,k)+ekh(im,j,k)) * MIN(fc_global(i,j+rank*jmax,k),fc_global(im,j+rank*jmax,k)) *
     1            (putin2(i,j,k) -putin2(im,j,k)) / (Rp(i) - Rp(im))    )
     1 * Rpdr_i !/ ( Rp(i) * dr(i) )
     +              +
     2 (    CNy*(ekh(i,j,k)+ekh(i,jp,k)) * MIN(fc_global(i,j+rank*jmax,k),fc_global(i,jp+rank*jmax,k)) 
     2 * (putin2(i,jp,k)-putin2(i,j,k) )/(Rp(i)*(phip(jp)-phip(j ))) -
     2      CNy*(ekh(i,j,k)+ekh(i,jm,k)) * MIN(fc_global(i,j+rank*jmax,k),fc_global(i,jm+rank*jmax,k))
     2 * (putin2(i,j,k) -putin2(i,jm,k))/(Rp(i)*(phip(j )-phip(jm)))  )
     2 / ( Rp(i) *(phiv(j )-phiv(jm)) )
     +              +
     3 (CNz*(ekh(i,j,k)+ekh(i,j,kp)) * MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,kp)) 
     3 * (putin2(i,j,kp)-putin2(i,j,k) ) -
     3  CNz*(ekh(i,j,k)+ekh(i,j,km)) * MIN(fc_global(i,j+rank*jmax,k),fc_global(i,j+rank*jmax,km)) 
     3 * (putin2(i,j,k) -putin2(i,j,km))  )
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
