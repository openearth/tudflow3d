      MODULE sparse_matrix
!___________________________________________________________________________
!     Copyright 2001 Svasek Hydraulics BV, Rotterdam 
!     Copyright 2004 Technische Universiteit Delft
!         Robert Jan Labeur
!        
!         Environmental Hydraulics Section
!         Faculty of Civil Engineering and Geosciences
!         PO Box 5048
!         2600 GA Delft
!         The Netherlands
!        
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Lesser General Public
!     License as published by the Free Software Foundation; either
!     version 2.1 of the License, or (at your option) any later version.
!
!     This library is distributed in the hope that it will be useful,
!     but without any warranty; without even the implied warranty of
!     merchantability or fitness for a particular purpose.  See the GNU
!     Lesser General Public License for more details.
!
!     You should have received a copy of the GNU Lesser General Public
!     License along with this library; if not, please contact the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA
!___________________________________________________________________________

      ! List of modules used
      !USE finite_elements
      !USE matrix_functions
      USE nlist

      IMPLICIT NONE
	  
      INCLUDE 'dmumps_struc.h'	  	  

      SAVE

      ! Direct solver structure
      TYPE (DMUMPS_STRUC) id	  
      TYPE (DMUMPS_STRUC), ALLOCATABLE, DIMENSION(:) :: id2
	  

	  
      ! Local performance variables
      INTEGER          :: solver_flag,solver_niter
      DOUBLE PRECISION :: solver_runtime,solver_residual

      ! Interfaces
      INTERFACE sparse
          MODULE PROCEDURE real_sparse
      END INTERFACE sparse

      INTERFACE MATVEC
          MODULE PROCEDURE real_MATVEC
      END INTERFACE MATVEC

      INTERFACE MILU
          MODULE PROCEDURE real_MILU
      END INTERFACE MILU

      INTERFACE UPPINV
          MODULE PROCEDURE real_UPPINV
      END INTERFACE UPPINV

      INTERFACE LOWINV
          MODULE PROCEDURE real_LOWINV
      END INTERFACE LOWINV

      INTERFACE MATINV
          MODULE PROCEDURE real_MATINV
      END INTERFACE MATINV
  
      INTERFACE BiCGSTAB
          MODULE PROCEDURE REAL_BiCGSTAB
      END INTERFACE BiCGSTAB  

      ! Sparse matrix variables
!      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, TARGET :: B,LUB
!      INTEGER,          DIMENSION(:), ALLOCATABLE, TARGET :: beg,jco,di,iro

      ! Local module variables
!      INTEGER, PRIVATE :: i,ii1,i2,j,j0,jj1,j2,k,kk1,k2,n,m
      INTEGER, PRIVATE :: ii,ii1,i2,jj,j0,jj1,j2,kk,kk1,k2,n,m
      
	  
      CONTAINS

!      SUBROUTINE INIT_SPARSE
      ! Allocate sparse matrix system
!      ALLOCATE(B(nmat),LUB(0:nmat))
!      ALLOCATE(jco(nmat),beg(nnode+1),di(nnode),iro(nmat))
      ! Fill jco, beg and di
!      END SUBROUTINE INIT_SPARSE

      SUBROUTINE real_sparse(A,B,p)
      ! Puts matrix B into sparse matrix A with symmetric permutation vector p
      ! Dummy variables
      DOUBLE PRECISION,    DIMENSION(:),   INTENT(INOUT) :: A
      DOUBLE PRECISION,    DIMENSION(:,:), INTENT(IN)    :: B
      INTEGER, DIMENSION(:),   INTENT(IN)    :: p
      n=SIZE(p)
      DO ii=1,n
          jj=p(ii)
          j0=di(jj)
          jj1=beg(jj)
          j2=beg(jj+1) - 1
          DO kk=1,n
              SELECT CASE (p(kk)-jj)
              CASE(0)  
                  kk1=j0
                  k2=j0
              CASE(1:)                  
                  kk1=j0 + 1
                  k2=j2
              CASE(:-1) 
                  kk1=jj1
                  k2=j0 - 1 
              END SELECT
              WHERE (jco(kk1:k2)==p(kk)) A(kk1:k2)=A(kk1:k2) + B(ii,kk)
          END DO
      END DO
      END SUBROUTINE real_sparse

      FUNCTION real_MATVEC(A,x) RESULT (y)
      ! Matrix-vector multiply
      ! Dummy variables
      DOUBLE PRECISION, DIMENSION(:) :: A,x
      ! Result variable
      DOUBLE PRECISION, DIMENSION(SIZE(x)) :: y
      ! Local variables
      INTEGER, DIMENSION(:), POINTER :: p
      DO ii=1,SIZE(x)
          ii1=beg(ii)
          i2=beg(ii+1)-1
          p=>jco(ii1:i2)
          y(ii)=DOT_PRODUCT(x(p),A(ii1:i2))
      END DO
      END FUNCTION real_MATVEC

      FUNCTION real_MILU(A,modify) RESULT(LU) 
      ! ILU/MILU-decomposition of matrix A 
      ! Dummy variables
      DOUBLE PRECISION, DIMENSION(:) :: A
      LOGICAL, OPTIONAL :: modify
      ! Result variable
      DOUBLE PRECISION, DIMENSION(0:SIZE(A)) :: LU
      ! Local variables
      INTEGER, DIMENSION(:), ALLOCATABLE :: point
      INTEGER v,w,right,diag
      LOGICAL change_diag
      IF(PRESENT(modify)) THEN
          change_diag=modify
      ELSE
          change_diag=.FALSE.
      END IF
      ALLOCATE(point(SIZE(di)))
      LU(1:)=A
      point=0
      DO ii=2,SIZE(di)
          right=beg(ii+1) - 1
          diag=di(ii)
          LU(0)=0
          DO v=beg(ii)+1,right
              point(jco(v))=v
          END DO 
          DO v=beg(ii),diag-1
              jj=jco(v)
              LU(v)=LU(v)/LU(di(jj))
              DO w=di(jj)+1,beg(jj+1)-1
                  kk=point(jco(w))
                  LU(kk)=LU(kk) - LU(v)*LU(w)
              END DO
          END DO   
          IF (change_diag) LU(diag)=LU(diag) + LU(0)
          DO v=beg(ii)+1,right
              point(jco(v))=0
          END DO
      END DO 
      DEALLOCATE(point)
      END FUNCTION real_MILU

      FUNCTION real_LOWINV(LU,x) RESULT (y)
      ! Computes inverse of L
      ! Dummy variables
      DOUBLE PRECISION, DIMENSION(0:) :: LU
      DOUBLE PRECISION, DIMENSION(:)  :: x
      ! Result variable
      DOUBLE PRECISION, DIMENSION(SIZE(x)) :: y
      ! Local variables
      INTEGER, DIMENSION(:), POINTER :: p
      ! Forward elimination
      DO ii=1,SIZE(x)
          ii1=beg(ii)
          i2=di(ii)-1
          p=>jco(ii1:i2)
          y(ii)=x(ii) - SUM(y(p)*LU(ii1:i2))
      END DO
      END FUNCTION real_LOWINV

      FUNCTION real_UPPINV(LU,x) RESULT (y)
      ! Computes inverse of U
      ! Dummy variables
      DOUBLE PRECISION, DIMENSION(0:) :: LU
      DOUBLE PRECISION, DIMENSION(:)  :: x
      ! Result variable
      DOUBLE PRECISION, DIMENSION(SIZE(x)) :: y
      ! Back substitution
      DO ii=SIZE(x),1,-1
          ii1=di(ii) + 1
          i2=beg(ii+1)-1
          y(ii)=x(ii)
          DO jj=ii1,i2
              y(ii)=y(ii) - LU(jj)*y(jco(jj))
          END DO
          y(ii)=y(ii)/LU(di(ii))
      END DO
      END FUNCTION real_UPPINV

      FUNCTION real_MATINV(LU,x) RESULT (y)
      ! Computes inverse of LU 
      ! Dummy variables
      DOUBLE PRECISION, DIMENSION(0:) :: LU
      DOUBLE PRECISION, DIMENSION(:)  :: x
      ! Result variable
      DOUBLE PRECISION, DIMENSION(SIZE(x)) :: y
      ! Forward elimination
      y=LOWINV(LU,x)
      ! Back substitution
      y=UPPINV(LU,y)
      END FUNCTION real_MATINV

      FUNCTION L(x) RESULT(y)
      ! Function for BiCG iteration (pre-conditioning)
      ! Dummy variable
      DOUBLE PRECISION, DIMENSION(:) :: x
      ! Result variable
      DOUBLE PRECISION, DIMENSION(SIZE(x)) :: y
      y=MATVEC(LHS2,x)
      y=MATINV(LUB,y)
      END FUNCTION L

      FUNCTION Lnopre(x) RESULT(y)
      ! Function for BiCG iteration (pre-conditioning)
      ! Dummy variable
      DOUBLE PRECISION, DIMENSION(:) :: x
      ! Result variable
      DOUBLE PRECISION, DIMENSION(SIZE(x)) :: y
      y=MATVEC(LHS2,x)
!      y=MATINV(LUB,y)

      END FUNCTION Lnopre
	  
	  
      FUNCTION REAL_BiCGSTAB(A,b,convergence,max_iter,start) RESULT(x)
      ! Dummy variables
      INTERFACE
          DOUBLE PRECISION FUNCTION A(x)
              DOUBLE PRECISION, DIMENSION(:) :: x
              DIMENSION A(SIZE(x))
          END FUNCTION A
      END INTERFACE
      DOUBLE PRECISION, DIMENSION(:) :: b
      DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: start
      DOUBLE PRECISION, OPTIONAL :: convergence
      INTEGER, OPTIONAL :: max_iter
      ! Result variable
      DOUBLE PRECISION, DIMENSION(SIZE(b)) :: x
      ! Local variables
      DOUBLE PRECISION, DIMENSION(SIZE(b)) :: p,q,r,s,t,v
      DOUBLE PRECISION alpha,beta,rho,rho0,omega
      DOUBLE PRECISION residual0,accuracy,time0,time1
      INTEGER :: iter,niter
!      CALL CPU_TIME(time0)
      solver_flag     = 0
      solver_niter    = 0
      solver_residual = 0
      solver_runtime  = 0
      x=0
      accuracy=1D-6
      niter=SIZE(b)
      IF (PRESENT(start)) x=start 
      IF (PRESENT(convergence)) accuracy=convergence
      IF (PRESENT(max_iter)) niter=max_iter   
      ! Starting values
      r=b-A(x)
      p=r
      q=0
      v=0
      rho=1
      alpha=1
      omega=1
      residual0=SUM(r*r)
      IF (ABS(residual0)==0) THEN
          WRITE(*,*) 'BiCGSTAB not executed: initial residual is zero'
!          CALL CPU_TIME(time1)
          solver_runtime=time1-time0
          RETURN
      END IF
      ! Iteration
      DO iter=1,niter
          rho0=rho
          rho=SUM(p*r)
          beta=rho*alpha/(rho0*omega)
          q=r + beta*(q - omega*v)
          v=A(q)
          alpha=rho/SUM(p*v)
          x=x + alpha*q
          s=r - alpha*v
          rho0=SUM(s*s)
          IF (ABS(rho0/residual0)<accuracy*accuracy) EXIT
          t=A(s)
          omega=SUM(t*s)/SUM(t*t)
          x=x + omega*s
          r=s - omega*t
          rho0=SUM(r*r)
          IF (ABS(rho0/residual0)<accuracy*accuracy) EXIT
          !IF (MOD(iter,10)==0) THEN 
          !    CALL CPU_TIME(time1)
          !    WRITE(*,'(A11,A1,A2,F16.6,I6,3E16.6)') ' BiCGSTAB ',time1-time0,iter,ABS(residual0),ABS(rho0),ABS(rho0/residual0)
          !END IF
      END DO
!      CALL CPU_TIME(time1)
      IF (ABS(rho0/residual0)>=accuracy*accuracy) THEN
          WRITE(*,*) 'Warning: BiCGSTAB did not converge sufficiently'
          WRITE(*,'(A11,A1,A2,F16.6,I6,3E16.6)') ' BiCGSTAB ',time1-time0,iter,ABS(residual0),ABS(rho0),ABS(rho0/residual0)
          solver_flag = 1
      END IF
      solver_niter    = iter
      solver_residual = ABS(rho0/residual0)
      solver_runtime  = time1-time0
!      WRITE(*,'(A11,A1,A2,F16.6,I6,3E16.6)') ' BiCGSTAB ',time1-time0,iter,ABS(residual0),ABS(rho0),ABS(rho0/residual0)
      END FUNCTION REAL_BiCGSTAB   

      END MODULE sparse_matrix

