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

      subroutine chkdt

      USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      include 'mpif.h'
      integer ierr
      real  dr2,dz2,df,df2,kcoeff,tmp1,tmp2,tmp3,Courant,dtmp

	  
!		IF (nobst>0.and.bp(1)%forever.eq.0.and.time_np+dt.gt.bp(1)%t0.and.counter<10) THEN
!			counter=counter+1
!			Courant=0.1
!			IF (rank.eq.0) THEN
!				write(*,*),'Placement bedplume, CFL=0.1 for 10 time steps ',counter,dt
!			ENDIF
!		ELSE 
			Courant = CFL
!		ENDIF
		
      IF (istep.le.10) THEN
	dt = dt_ini
      ELSE
	dt = dt_max
      ENDIF
	  
      
      IF (CNdiffz.eq.1) THEN
	dz2 = 1.e9*dz*dz ! no dt restriction for vertical diff with CN implicit scheme
      ELSE
        dz2 = dz    * dz
      ENDIF
      do k=1,kmax
         do j=1,jmax
            do i=1,imax
            df = Rp(i)*dphi
      	    df2 = df * df
            dr2 = dr(i) * dr(i)
            kcoeff = ekm(i,j,k) /rnew(i,j,k)
            tmp1 = ( abs(Unew(i,j,k)) / ( Rp(i+1)-Rp(i) ) ) +
     &             ( abs(Vnew(i,j,k)) /         df        ) +
     &             ( abs(Wnew(i,j,k)) /         dz        )
            tmp2 = ( 1.0/dr2 + 1.0/df2 + 1.0/dz2 )
            tmp3 = 1.0 / ( 1.0 * tmp2 * kcoeff + tmp1 + 1.e-12 )
            tmp3 =Courant *tmp3 
            dt = min( dt , tmp3 )
            dtmp = dt
!		if (abs(dt)/abs(dt_max).le.0.7) then
!			write(*,*)'rank,i,j,k,dt=',rank,i,j,k,dt
!			write(*,*)'U,V,W:',Unew(i,j,k),Vnew(i,j,k),Wnew(i,j,k)
!			write(*,*)'*i-1: U,V,W:',Unew(i-1,j,k),Vnew(i-1,j,k),Wnew(i-1,j,k)
!			write(*,*)'*j-1:U,V,W:',Unew(i,j-1,k),Vnew(i,j-1,k),Wnew(i,j-1,k)
!			write(*,*)'*k-1:U,V,W:',Unew(i,j,k-1),Vnew(i,j,k-1),Wnew(i,j,k-1)
!			write(*,*)'*i+1: U,V,W:',Unew(i+1,j,k),Vnew(i+1,j,k),Wnew(i+1,j,k)
!			write(*,*)'*j+1:U,V,W:',Unew(i,j+1,k),Vnew(i,j+1,k),Wnew(i,j+1,k)
!			write(*,*)'*k+1:U,V,W:',Unew(i,j,k+1),Vnew(i,j,k+1),Wnew(i,j,k+1)
!			write(*,*)'Rp*dphi,Rp(i+1)-Rp(i),dz',df,Rp(i+1)-Rp(i),dz
!			write(*,*)'dyn viscosity nu,ekm,rho:',kcoeff,ekm(i,j,k),rnew(i,j,k)
!			write(*,*)'dt CFL,dt visc:',1./tmp1,1./(tmp2*kcoeff)
!			call output(99999,time_np)
!			CALL writeerror(10000)
!		endif
            enddo
         enddo
      enddo

	if (dt<1.e-12) then 
	   write(*,*),'rank,dt',rank,dt 
	endif
      
	call mpi_allreduce(dtmp,dt,1,mpi_real8,mpi_min,mpi_comm_world,ierr)
!	call mpi_allreduce(dtmp,dt,1,mpi_real,mpi_min,mpi_comm_world,ierr)


	if (isnan(dt)) stop 'ERROR, QUITING DFLOW3D: "dt" is a NaN'
	if (dt<1.e-12) stop 'ERROR, QUITING DFLOW3D: "dt" is smaller than 1e-12'

	if (time_int.eq.'AB2'.or.time_int.eq.'AB3') then 
	  if (dt<dt_max) then
		stop 'ERROR, QUITING DFLOW3D: "dt" is adjusted in AB2 or AB3, sign of instability try again with smaller dt'
	  endif
	endif
!! 	if AB2 or AB3 timestep must be fixed, if unfortunately initially a too large timestep has been chosen 
!! 	then a 25% reduction must lead to stable simulation from this point on (one change in timestep can be 
!!	stable, but subsequent changes in timestep are not stable for AB2 or AB3)
!	if (time_int.eq.'AB2'.or.time_int.eq.'AB3') then 
!	  if (dt<dt_max) then  
!		WRITE(*,'(a,a,a,f12.9,a,f12.9,a)') ' WARNING: Using ',time_int,' with dt too large, dt changed from : '
!     &           ,dt_max,' s, to : ',0.75*dt_max,' s'  
!	    dt_max=0.75*dt_max
!	    dt=dt_max
!	  endif
!	endif
      end

      subroutine chkdiv
      USE nlist
      implicit none
!       include 'param.txt'
!       include 'common.txt'
      include 'mpif.h'
      real   div1,divmax1,divbar1,divmax_tot1,divbar_tot1,rhoip,rhoim,rhojp,rhojm,rhokp,rhokm
      real   div2,divmax2,divbar2,divmax_tot2,divbar_tot2
      integer ierr
      divbar1 = 0.0
      divmax1 = 0.0
      divbar2 = 0.0
      divmax2 = 0.0
      do k=2,kmax-1
         do j=2,jmax-1
            do i=2,imax-1
!       rhoip =0.5*(drdt(i,j,k)+drdt(i+1,j,k))
!       rhoim =0.5*(drdt(i,j,k)+drdt(i-1,j,k))
!       rhojp =0.5*(drdt(i,j,k)+drdt(i,j+1,k))
!       rhojm =0.5*(drdt(i,j,k)+drdt(i,j-1,k))
!       rhokp =0.5*(drdt(i,j,k)+drdt(i,j,k+1))
!       rhokm =0.5*(drdt(i,j,k)+drdt(i,j,k-1))
!       div = (Ru(i)*Unew(i,j,k)*rhoip-Ru(i-1)*Unew(i-1,j,k)*rhoim)*dphi*dz   +
!      +      (      Vnew(i,j,k)*rhojp-        Vnew(i,j-1,k)*rhojm)*dr(i)*dz  +
!      +      (      Wnew(i,j,k)*rhokp-        Wnew(i,j,k-1)*rhokm)*Rp(i)*dphi*dr(i)+
!      +  ((3*drdt(i,j,k)-4*rnew(i,j,k)+rold(i,j,k))/(2*dt))*rp(i)*dr(i)*dphi*dz

!       rhoip =0.5*(drdt(i,j,k)+drdt(i+1,j,k))
!       rhoim =0.5*(drdt(i,j,k)+drdt(i-1,j,k))
!       rhojp =0.5*(drdt(i,j,k)+drdt(i,j+1,k))
!       rhojm =0.5*(drdt(i,j,k)+drdt(i,j-1,k))
!       rhokp =0.5*(drdt(i,j,k)+drdt(i,j,k+1))
!       rhokm =0.5*(drdt(i,j,k)+drdt(i,j,k-1))

!       div = (Ru(i)*Unew(i,j,k)*rhoip-Ru(i-1)*Unew(i-1,j,k)*rhoim)/(rp(i)*dr(i))   +
!      +      (      Vnew(i,j,k)*rhojp-        Vnew(i,j-1,k)*rhojm)/(rp(i)*dphi)  +
!      +      (      Wnew(i,j,k)*rhokp-        Wnew(i,j,k-1)*rhokm)/dz+
!      +  ((3*drdt(i,j,k)-4*rnew(i,j,k)+rold(i,j,k))/(2*dt))

	div1 =
     1  ( Ru(i)*dUdt(i,j,k) - Ru(i-1)*dUdt(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       dVdt(i,j,k) -         dVdt(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       dWdt(i,j,k) -         dWdt(i,j,k-1) ) / ( dz )
!     +              +  (3.*drdt(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/(2.*dt)
     +              +  (3.*drdt(i,j,k)-4.*rnew(i,j,k)+rold(i,j,k))/((3.*time_np-4.*time_n+time_nm))

	    div2= ( Ru(i)*Unew(i,j,k) - Ru(i-1)*Unew(i-1,j,k) ) / ( Rp(i)*dr(i) )
     +              +
     2  (       Vnew(i,j,k) -         Vnew(i,j-1,k) ) / ( Rp(i)*dphi )
     +              +
     3  (       Wnew(i,j,k) -         Wnew(i,j,k-1) ) / ( dz )  

      divbar1 = divbar1 + div1
      div1    = abs(div1)
      divmax1 = max( divmax1 , div1 )
      divbar2 = divbar2 + div2
      div2    = abs(div2)
      divmax2 = max( divmax2 , div2 )

           enddo
         enddo
      enddo
      call mpi_allreduce(divbar1,divbar_tot1,1,mpi_real8,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(divmax1,divmax_tot1,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(divbar2,divbar_tot2,1,mpi_real8,mpi_max,mpi_comm_world,ierr)
      call mpi_allreduce(divmax2,divmax_tot2,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
      !call mpi_allreduce(divbar,divbar_tot,1,mpi_real,mpi_max,mpi_comm_world,ierr)
      !call mpi_allreduce(divmax,divmax_tot,1,mpi_real,mpi_sum,mpi_comm_world,ierr)
      
	if (rank.eq.0) write(6,100)divbar_tot1,divmax_tot1   
100   format('Mass loss/gain : Tot = ',e13.6,
     +                      '  Max = ',e13.6)
	if (rank.eq.0) write(6,101)divbar_tot2,divmax_tot2  
101   format('Div(u) : Tot = ',e13.6,
     +                      '  Max = ',e13.6)
      end
      
      
!      subroutine calc_div
!      USE nlist
!      implicit none

!      real dz_i,Rpdr_i,Rpdphi_i

!      dz_i=1./dz
!      do i=1,imax
!	Rpdr_i=1./(Rp(i)*dr(i))
!	Rpdphi_i=1./(Rp(i)*dphi)
!	do j=1,jmax
!	  do k=1,kmax

!	    div(i,j,k)= ( Ru(i)*Unew(i,j,k) - Ru(i-1)*Unew(i-1,j,k) ) *Rpdr_i ! / ( Rp(i)*dr(i) )
!     +              +
!     2  (       Vnew(i,j,k) -         Vnew(i,j-1,k) ) * Rpdphi_i !/ ( Rp(i)*dphi )
!     +              +
!     3  (       Wnew(i,j,k) -         Wnew(i,j,k-1) ) * dz_i !/ ( dz )  ) 
!	  enddo
!	enddo
!      enddo
!      end

	MODULE work_array
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: Uavg,Vavg,Wavg,Ravg,Pavg,muavg
	  REAL, DIMENSION(:,:,:), ALLOCATABLE :: sigU2,sigV2,sigW2,sigR2,sigUV,sigUW,sigVW
	  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: sigC2,Cavg,sigUC,sigVC,sigWC,Cmax,Cmin
	  INTEGER stat_count
	REAL stat_time_count
!        REAL, DIMENSION(:,:,:), ALLOCATABLE :: fUavg,sigfU2
      END MODULE
      


	
	subroutine statistics(U,V,W,C,R)
	USE work_array
	USE nlist
	implicit none
! 	include 'param.txt'
! 	include 'common.txt'

	REAL U(0:i1,0:j1,0:k1),V(0:i1,0:j1,0:k1),W(0:i1,0:j1,0:k1)
	REAL R(0:i1,0:j1,0:k1),C(nfrac,0:i1,0:j1,0:k1)
	REAL uu,vv,ww
	INTEGER n

!	REAL, ALLOCATABLE :: Uavg(:,:,:),Vavg(:,:,:),Wavg(:,:,:),Cavg(:,:,:),Ravg(:,:,:)
!	REAL, ALLOCATABLE :: sigU2(:,:,:),sigV2(:,:,:),sigW2(:,:,:),sigC2(:,:,:),sigR2(:,:,:)  
!	REAL, ALLOCATABLE :: Urms(:,:,:),Vrms(:,:,:),Wrms(:,:,:),Crms(:,:,:),Rrms(:,:,:)
!	INTEGER :: stat_count

	
	IF (.not.ALLOCATED(Uavg)) then
		ALLOCATE (Uavg(1:imax,1:jmax,1:kmax),Vavg(1:imax,1:jmax,1:kmax),
     1                Wavg(1:imax,1:jmax,1:kmax),Cavg(nfrac,1:imax,1:jmax,1:kmax),
     1                Ravg(1:imax,1:jmax,1:kmax),
     1		      Pavg(1:imax,1:jmax,1:kmax),muavg(1:imax,1:jmax,1:kmax),						
     1		      sigU2(1:imax,1:jmax,1:kmax),sigV2(1:imax,1:jmax,1:kmax),
     1                sigW2(1:imax,1:jmax,1:kmax),sigC2(nfrac,1:imax,1:jmax,1:kmax),
     1                sigR2(1:imax,1:jmax,1:kmax),
     1		      sigUV(1:imax,1:jmax,1:kmax),sigVW(1:imax,1:jmax,1:kmax),
     1                sigUW(1:imax,1:jmax,1:kmax) ,
     1		      sigUC(nfrac,1:imax,1:jmax,1:kmax),sigVC(nfrac,1:imax,1:jmax,1:kmax),
     1                sigWC(nfrac,1:imax,1:jmax,1:kmax),Cmax(nfrac,1:imax,1:jmax,1:kmax),
     1                Cmin(nfrac,1:imax,1:jmax,1:kmax))
 !    1 			,fUavg(1:imax,1:jmax,1:kmax),sigfU2(1:imax,1:jmax,1:kmax))

		stat_count=0
		stat_time_count=0.
		sigU2=0.
		sigV2=0.
		sigW2=0.
		sigC2=0.
		sigR2=0.
		sigUV=0.
		sigVW=0.
		sigUW=0.
		sigUC=0.	
		sigVC=0.
		sigWC=0.

		Uavg=0.
		Vavg=0.
		Wavg=0.
		Cavg=0.
		Cmax=0.
		Cmin=1.e18
		Ravg=0.
		Pavg=0.
		muavg=0.
!		sigfU2=0.
!		fUavg=0.
	endif


	stat_count=stat_count+1
	stat_time_count=stat_time_count+dt
	
	do i=1,imax
	  do j=1,jmax
	    do k=1,kmax
		  uu=0.5*(U(i,j,k)+U(i-1,j,k))*cos_u(j)-0.5*(V(i,j,k)+V(i,j-1,k))*sin_u(j)
		  vv=0.5*(V(i,j,k)+V(i,j-1,k))*cos_u(j)+0.5*(U(i,j,k)+U(i-1,j,k))*sin_u(j)
		  ww=0.5*(W(i,j,k)+W(i,j,k-1))
                !  uu=(U(i,j,k))*cos_u(j)-(V(i,j,k))*sin_v(j)
                !  vv=(V(i,j,k))*cos_v(j)+(U(i,j,k))*sin_u(j)
                !  ww=W(i,j,k)

		  sigU2(i,j,k) = sigU2(i,j,k) + uu*uu*dt
		  sigV2(i,j,k) = sigV2(i,j,k) + vv*vv*dt
		  sigW2(i,j,k) = sigW2(i,j,k) + ww*ww*dt
		  sigR2(i,j,k) = sigR2(i,j,k) + R(i,j,k)*R(i,j,k)*dt
		  sigUV(i,j,k) = sigUV(i,j,k) + uu*vv*dt
		  sigUW(i,j,k) = sigUW(i,j,k) + uu*ww*dt
		  sigVW(i,j,k) = sigVW(i,j,k) + vv*ww*dt

		  Uavg(i,j,k)  = Uavg(i,j,k) + uu*dt
		  Vavg(i,j,k)  = Vavg(i,j,k) + vv*dt
		  Wavg(i,j,k)  = Wavg(i,j,k) + ww*dt
		  Ravg(i,j,k)  = Ravg(i,j,k) + R(i,j,k)*dt
		  Pavg(i,j,k)  = Pavg(i,j,k) + (p(i,j,k)+pold(i,j,k))*dt
		  muavg(i,j,k) = muavg(i,j,k) + ekm(i,j,k)*dt
		do n=1,nfrac
		  Cavg(n,i,j,k)  = Cavg(n,i,j,k) + C(n,i,j,k)*dt
		  sigC2(n,i,j,k) = sigC2(n,i,j,k) + C(n,i,j,k)*C(n,i,j,k)*dt
		  sigUC(n,i,j,k) = sigUC(n,i,j,k) + uu*C(n,i,j,k)*dt
		  sigVC(n,i,j,k) = sigVC(n,i,j,k) + vv*C(n,i,j,k)*dt
		  sigWC(n,i,j,k) = sigWC(n,i,j,k) + ww*C(n,i,j,k)*dt
		  Cmax(n,i,j,k) = MAX(Cmax(n,i,j,k),C(n,i,j,k))
		  Cmin(n,i,j,k) = MIN(Cmin(n,i,j,k),C(n,i,j,k))
		enddo
!!! Favre average results in differences of <1 promille for Umean and <1% for U' between Favre and Reynolds avg
!!! for CO2 plume simulation (Wang et al. 2008), this is not important for dredge plumes, thus only Reynolds avg is calculated from now on
!		  fUavg(i,j,k) = fUavg(i,j,k) + uu*R(i,j,k)
!		  sigfU2(i,j,k) = sigfU2(i,j,k) + uu*uu*R(i,j,k)
	    enddo
	  enddo
	enddo
!	write(*,*),'stat_count,sigU2(1,1,1)',stat_count,sigU2(1,1,1),Uavg(1,1,1)
	end

	MODULE history_array
	use nlist
	IMPLICIT NONE

	INTEGER,PARAMETER :: tdim=200000 ! length of history array in the time
	INTEGER :: nhispoint
	
	type history
	    integer ::	i,j,k
	    REAL ::	x,y,z
	    REAL*8, DIMENSION(tdim) :: U,V,W,P,RHO,t
	    REAL*8, DIMENSION(20,tdim) :: C ! 20 fractions in history, must be enough..
	end type history

	TYPE(history), DIMENSION(:), ALLOCATABLE :: his
      END MODULE	
	
	subroutine init_his
	USE history_array
	USE error_functions
	USE nlist
	implicit none
	
!       include 'param.txt'
!       include 'common.txt'

	integer ios,n

	type history_init
	  integer :: i,j,k
	end type history_init

	TYPE(history_init), DIMENSION(100) :: hist

	NAMELIST /histories/hist
	
	hist(:)%i=0
	OPEN(1,FILE=hisfile,IOSTAT=ios,ACTION='read')
	IF (ios/=0) CALL writeerror(100)

	READ (UNIT=1,NML=histories,IOSTAT=ios)

	nhispoint=0
	DO WHILE (hist(nhispoint+1)%i.NE.0)
	  nhispoint=nhispoint+1
	END DO

	!write(*,*)'hist1,2,3,4',hist(1)%i,hist(2)%i,hist(3)%i,hist(4)%i
	!write(*,*)'nhispoint',nhispoint
	CLOSE(1)

	ALLOCATE(his(nhispoint))
	!write(*,*)'his,2,3,4',his(1)%i,his(2)%i,his(3)%i,his(4)%i
	
	DO n=1,nhispoint
	  his(n)%i=hist(n)%i
	  his(n)%j=hist(n)%j
	  his(n)%k=hist(n)%k

	  IF (his(n)%i>imax) CALL writeerror(110)
	  IF (his(n)%i<0) CALL writeerror(111)
	  IF (his(n)%j>jmax*px) CALL writeerror(120)
	  IF (his(n)%j<0) CALL writeerror(121)
	  IF (his(n)%k>kmax) CALL writeerror(130)
	  IF (his(n)%k<0) CALL writeerror(131)

	  his(n)%x=Rp(his(n)%i)*cos_u(his(n)%j-rank*jmax)-schuif_x
	  his(n)%y=Rp(his(n)%i)*sin_u(his(n)%j-rank*jmax)
	  his(n)%z=his(n)%k*dz-0.5*dz
	ENDDO
	    

	END SUBROUTINE init_his

! 	subroutine appendhis(rank,istep,time,Unew,Vnew,Wnew,P,Cnew,Rnew)
	subroutine appendhis
	USE history_array
	USE nlist

	implicit none


!       include 'param.txt'
      !include 'common.txt'


! 	REAL Unew(0:i1,0:j1,0:k1),Vnew(0:i1,0:j1,0:k1),Wnew(0:i1,0:j1,0:k1)
! 	REAL Rnew(0:i1,0:j1,0:k1),Cnew(0:i1,0:j1,0:k1),P(0:i1,0:j1,0:k1)
! 	REAL time

	integer n,r,nf
	
	his(1)%t(istep)=time_np !Unew is already updated, therefore it is output at time_np
	DO n=1,nhispoint
	  r=INT((his(n)%j-1)/jmax)
	  IF (rank.eq.r) THEN
	    i=his(n)%i
	    j=his(n)%j-rank*jmax
	    k=his(n)%k
	    his(n)%U(istep)=Unew(i,j,k)*cos_u(j)-Vnew(i,j,k)*sin_v(j)
	    his(n)%V(istep)=Vnew(i,j,k)*cos_v(j)+Unew(i,j,k)*sin_u(j)
	    his(n)%W(istep)=Wnew(i,j,k)
	    his(n)%P(istep)=Pold(i,j,k)+P(i,j,k)
	    his(n)%RHO(istep)=Rnew(i,j,k)
	    DO nf=1,nfrac
	    	his(n)%C(nf,istep)=Cnew(nf,i,j,k)
	    ENDDO
	  ENDIF
	ENDDO
	end SUBROUTINE

! 	subroutine finalize_his(rank,istep)
	subroutine finalize_his
	USE history_array
	USE nlist
        USE netcdf
	implicit none

!       include 'param.txt'
      include 'mpif.h'
!	real UU(1:200000)
	integer ierr,n,r,tag,status(MPI_STATUS_SIZE),nf
	INTEGER (8) :: ps
	REAL :: Uhis(1:nhispoint,1:istep),Vhis(1:nhispoint,1:istep),Whis(1:nhispoint,1:istep)
	REAL :: Phis(1:nhispoint,1:istep),Chis(1:nfrac,1:nhispoint,1:istep),RHOhis(1:nhispoint,1:istep)
	REAL :: xhis(1:nhispoint),yhis(1:nhispoint),zhis(1:nhispoint)
	INTEGER :: ihis(1:nhispoint),jhis(1:nhispoint),khis(1:nhispoint)
	REAL :: this(1:istep)

       ! We are writing 3D, 2D and 1D data, a nhispoint x istep grid.
       integer, parameter :: NDIMS1 = 2
       integer, parameter :: NDIMS2 = 1
       integer, parameter :: NDIMS3 = 3
       integer, parameter :: NDIMS4 = 1
!       integer, parameter :: NX = imax, NY = jmax, NZ = kmax
     
       ! When we create netCDF files, variables and dimensions, we get back
       ! an ID for each one.
       integer :: ncid, varid1,varid2,varid3, varid4, varid5, varid6, varid7, varid8 
       integer :: varid9,varid10,varid11,varid12,varid13
       integer :: dimids1(NDIMS1), dimids2(NDIMS2),dimids3(NDIMS3),dimids4(NDIMS4)
       integer :: nhis_dimid,time_dimid,nfrac_dimid
	character(1024) :: svnversion
	character(1024) :: svnurl
      include 'version.inc'


	IF (rank>0) THEN
	  DO n=1,nhispoint
	    r=INT((his(n)%j-1)/jmax)
	    IF (rank.eq.r) THEN
	      tag=INT(r*n)
	      call mpi_sendrecv_replace(his(n)%U ,200000,MPI_REAL8,0,10,0,10, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%V ,200000,MPI_REAL8,0,11,0,11, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%W ,200000,MPI_REAL8,0,12,0,12, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%P ,200000,MPI_REAL8,0,13,0,13, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%C,20*200000,MPI_REAL8,0,14,0,14, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%RHO ,200000,MPI_REAL8,0,15,0,15, MPI_COMM_WORLD,status,ierr)
 	      call mpi_sendrecv_replace(his(n)%x ,1,MPI_REAL8,0,16,0,16, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%y ,1,MPI_REAL8,0,17,0,17, MPI_COMM_WORLD,status,ierr)	      
	      call mpi_sendrecv_replace(his(n)%z ,1,MPI_REAL8,0,18,0,18, MPI_COMM_WORLD,status,ierr)
	    ENDIF
	  ENDDO
	ELSE
	  DO n=1,nhispoint
	    r=INT((his(n)%j-1)/jmax)
	    IF (rank.ne.r) THEN
 	      call mpi_sendrecv_replace(his(n)%U ,200000,MPI_REAL8,r,10,r,10, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%V ,200000,MPI_REAL8,r,11,r,11, MPI_COMM_WORLD,status,ierr)	      
	      call mpi_sendrecv_replace(his(n)%W ,200000,MPI_REAL8,r,12,r,12, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%P ,200000,MPI_REAL8,r,13,r,13, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%C ,20*200000,MPI_REAL8,r,14,r,14, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%RHO ,200000,MPI_REAL8,r,15,r,15, MPI_COMM_WORLD,status,ierr)
 	      call mpi_sendrecv_replace(his(n)%x ,1,MPI_REAL8,r,16,r,16, MPI_COMM_WORLD,status,ierr)
	      call mpi_sendrecv_replace(his(n)%y ,1,MPI_REAL8,r,17,r,17, MPI_COMM_WORLD,status,ierr)	      
	      call mpi_sendrecv_replace(his(n)%z ,1,MPI_REAL8,r,18,r,18, MPI_COMM_WORLD,status,ierr)
	    ENDIF
	  ENDDO
	ENDIF

	IF (rank.eq.0) THEN
	  DO n=1,nhispoint
	    Uhis(n,1:istep)=his(n)%U(1:istep)
	    Vhis(n,1:istep)=his(n)%V(1:istep)
	    Whis(n,1:istep)=his(n)%W(1:istep)
	    DO k=1,nfrac
	      Chis(k,n,1:istep)=his(n)%C(k,1:istep)
	    ENDDO
	    Phis(n,1:istep)=his(n)%P(1:istep)
	    RHOhis(n,1:istep)=his(n)%RHO(1:istep)
		
	    xhis(n)=his(n)%x
	    yhis(n)=his(n)%y
	    zhis(n)=his(n)%z
	    ihis(n)=his(n)%i
	    jhis(n)=his(n)%j
	    khis(n)=his(n)%k

	  ENDDO
	this(1:istep)=his(1)%t(1:istep)

       ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
       ! overwrite this file, if it already exists.
       call check( nf90_create('history.nc', NF90_CLOBBER, ncid) )
     
       ! Define the dimensions. NetCDF will hand back an ID for each.
       call check( nf90_def_dim(ncid, "nhis_dim", nhispoint, nhis_dimid) )
       call check( nf90_def_dim(ncid, "time_dim", istep, time_dimid) )
	if (nfrac>0) then
       call check( nf90_def_dim(ncid, "nfrac_dim", nfrac, nfrac_dimid) )
	endif
     
       ! The dimids array is used to pass the IDs of the dimensions of
       ! the variables. Note that in fortran arrays are stored in
       ! column-major format.
       dimids1 =  (/ nhis_dimid, time_dimid/)
       dimids2 =  (/ nhis_dimid /)
	if (nfrac>0) then
       dimids3 =  (/ nfrac_dimid, nhis_dimid, time_dimid/)
	endif
       dimids4 =  (/ time_dimid/)
    
       ! Define the variable. The type of the variable in this case is
       ! NF90_DOUBLE (4-byte double).
       call check( nf90_def_var(ncid, "Uhis", NF90_DOUBLE, dimids1, varid1) )
       call check( nf90_put_att(ncid, varid1, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid1, 'long_name', 'Time history of U velocity at specific history location') )

       call check( nf90_def_var(ncid, "Vhis", NF90_DOUBLE, dimids1, varid2) )
       call check( nf90_put_att(ncid, varid2, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid2, 'long_name', 'Time history of V velocity at specific history location') )

       call check( nf90_def_var(ncid, "Whis", NF90_DOUBLE, dimids1, varid3) )
       call check( nf90_put_att(ncid, varid3, 'units', 'm/s') )
       call check( nf90_put_att(ncid, varid3, 'long_name', 'Time history of W velocity at specific history location') )

       call check( nf90_def_var(ncid, "RHOhis", NF90_DOUBLE, dimids1, varid4) )
       call check( nf90_put_att(ncid, varid4, 'units', 'kg/m3') )
       call check( nf90_put_att(ncid, varid4, 'long_name', 'Time history of density at specific history location') )

       call check( nf90_def_var(ncid, "Phis", NF90_DOUBLE, dimids1, varid5) )
       call check( nf90_put_att(ncid, varid5, 'units', 'Pa') )
       call check( nf90_put_att(ncid, varid5, 'long_name', 'Time history of pressure at specific history location') )

       if (nfrac>0) then
       call check( nf90_def_var(ncid, "Chis", NF90_DOUBLE, dimids3, varid6) )
       call check( nf90_put_att(ncid, varid6, 'units', '-') )
       call check( nf90_put_att(ncid, varid6, 'long_name', 'Time history of volume concentration 
     & of each fraction at specific history location') )
	endif

       call check( nf90_def_var(ncid, "xhis", NF90_DOUBLE, dimids2, varid7) )
       call check( nf90_put_att(ncid, varid7, 'units', 'm') )
       call check( nf90_put_att(ncid, varid7, 'long_name', 'X coordinate of specific history location') )

       call check( nf90_def_var(ncid, "yhis", NF90_DOUBLE, dimids2, varid8) )
       call check( nf90_put_att(ncid, varid8, 'units', 'm') )
       call check( nf90_put_att(ncid, varid8, 'long_name', 'Y coordinate of specific history location') )

       call check( nf90_def_var(ncid, "zhis", NF90_DOUBLE, dimids2, varid9) )
       call check( nf90_put_att(ncid, varid9, 'units', 'm') )
       call check( nf90_put_att(ncid, varid9, 'long_name', 'Z coordinate of specific history location') )

       call check( nf90_def_var(ncid, "ihis", NF90_INT, dimids2, varid10) )
       call check( nf90_put_att(ncid, varid10, 'units', '-') )
       call check( nf90_put_att(ncid, varid10, 'long_name', 'i coordinate of specific history location') )

       call check( nf90_def_var(ncid, "jhis", NF90_INT, dimids2, varid11) )
       call check( nf90_put_att(ncid, varid11, 'units', '-') )
       call check( nf90_put_att(ncid, varid11, 'long_name', 'j coordinate of specific history location') )

       call check( nf90_def_var(ncid, "khis", NF90_INT, dimids2, varid12) )
       call check( nf90_put_att(ncid, varid12, 'units', '-') )
       call check( nf90_put_att(ncid, varid12, 'long_name', 'k coordinate of specific history location') )

       call check( nf90_def_var(ncid, "time_his", NF90_DOUBLE, dimids4, varid13) )
       call check( nf90_put_att(ncid, varid13, 'units', 's') )
       call check( nf90_put_att(ncid, varid13, 'long_name', 'Time series of history') )

	! also add svn info in output files:
       CALL check( nf90_put_att(ncid,nf90_global, "svnversion", trim(svnversion)))
       CALL check( nf90_put_att(ncid,nf90_global, "svnurl", trim(svnurl)))

       ! End define mode. This tells netCDF we are done defining metadata.
       call check( nf90_enddef(ncid) )
     
       ! Write the pretend data to the file. Although netCDF supports
       ! reading and writing subsets of data, in this case we write all the
       ! data in one operation.
       call check( nf90_put_var(ncid, varid1, Uhis) )
       call check( nf90_put_var(ncid, varid2, Vhis) )
       call check( nf90_put_var(ncid, varid3, Whis) )
       call check( nf90_put_var(ncid, varid4, RHOhis) )
       call check( nf90_put_var(ncid, varid5, Phis) )
	if (nfrac>0) then
       call check( nf90_put_var(ncid, varid6, Chis) )
	endif
       call check( nf90_put_var(ncid, varid7, xhis) )
       call check( nf90_put_var(ncid, varid8, yhis) )
       call check( nf90_put_var(ncid, varid9, zhis) )
       call check( nf90_put_var(ncid, varid10, ihis) )
       call check( nf90_put_var(ncid, varid11, jhis) )
       call check( nf90_put_var(ncid, varid12, khis) )
       call check( nf90_put_var(ncid, varid13, this) )

       call check( nf90_close(ncid) )

	ENDIF

	END SUBROUTINE
    
