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

    MODULE Module_RwMatfiles

	! Construct filename
	!      mkmatname(5,12,FH,VEC,outfilename) => FH00012VEC.mat = outfilename
	! READ and WRITE matfiles (Matlab binary format) 
	! Write double, single, integer or strings 
	!      CALL wvar2matf(x,'xdouble.mat' ,'xdouble')
	!      CALL wvar2matf(x,'xsingle.mat' ,'xsingle')
	!      CALL wvar2matf('blabla','bla.mat' ,'string')
	! Read double, single, integer or strings 
	!      CALL rmatf2var(xr,'mvconst0000.mat' ,'xdouble')
	!      CALL rmatf2var(xr,'mvconst0000.mat' ,'xsingle')
	! Write Structure double, single, integer or strings 
	!      CALL Ini_Struc2matf('structure.mat','His','x',10)
	!      CALL wstruc2matf(x,'structure.mat','H','x',1)
	!
      ! interfaces

	! 
	!___________________________________________________________________________
	!     Copyright 2008 Svasek Hydraulics 
	! 	Permission to copy or distribute this software or documentation 
	!	in hard copy or soft copy granted only by written license 
	!	obtained from Svasek Hydraulics.
	!	All rights reserved. No part of this publication may be 
	!	reproduced, stored in a retrieval system (e.g., in memory, disk, 
	!	or core) or be transmitted by any means, electronic, mechanical, 
	!	photocopy, recording, or otherwise, without written permission 
	!	from the publisher.
	!___________________________________________________________________________
      ! B. Les 2008

     IMPLICIT NONE


     INTERFACE wvar2matf
       MODULE PROCEDURE doublemat3d_Wvar2Matf,singlemat3d_Wvar2Matf, &
					  doublemat_Wvar2Matf, singlemat_Wvar2Matf,&
					  doublevec_Wvar2Matf, singlevec_Wvar2Matf, &
						doublematint_Wvar2Matf,singlematint_Wvar2Matf, &
						doublereal_Wvar2Matf, &
	                    char_Wvar2Matf, singleint_wvar2matf, singlerealscal_wvar2matf
     END INTERFACE 
	INTERFACE wstruc2matf
	   MODULE PROCEDURE  doublemat_wstruc2matf, doublevec_wstruc2matf, &
                         singlemat_wstruc2matf, singlevec_wstruc2matf, &
                         integermat_wstruc2matf,integervec_wstruc2matf,integer_wstruc2matf, &
						 real_wstruc2matf,doublereal_wstruc2matf,char_wstruc2matf
     END INTERFACE 
	INTERFACE wstruc2ps
	   MODULE PROCEDURE  doublemat_wstruc2ps, doublevec_wstruc2ps, &
                         singlemat_wstruc2ps, singlevec_wstruc2ps, &
                         integermat_wstruc2ps,integervec_wstruc2ps,integer_wstruc2ps, &
						 real_wstruc2ps,doublereal_wstruc2ps,char_wstruc2ps
     END INTERFACE 
	INTERFACE rmatf2var
       MODULE PROCEDURE  doublemat3d_rmatf2var,doublemat_rmatf2var, doublevec_rmatf2var, doublerealscal_rmatf2var, &
	                     singlemat_rmatf2var, singlevec_rmatf2var, singlerealscal_rmatf2var, &
                         integermat_rmatf2var,integervec_rmatf2var 
     END INTERFACE 

     CONTAINS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE matgetdim(matname,keyw,dims)

    ! i/o variables
    INTEGER (4)                  ,   INTENT(OUT)    :: dims(2)
    CHARACTER(LEN=*)             ,   INTENT(IN)     :: matname
    CHARACTER(LEN=*)             ,   INTENT(IN)     :: keyw
    ! local variables         
    INTEGER (8) ::    mp, pa
    ! Matlab call variables
    INTEGER (4) ::    matGetVariable, mxIsNumeric,            &
                      mxGetNumberOfDimensions,                &
					  mxGetM,mxGetN,                          &
                      mxGetDimensions

    ! Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                               
    ! Create, write and free from memory  
    pa = matGetVariable( mp ,keyw)
    IF (mxIsNumeric(pa) .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable not numeric= ',keyw
       STOP
    ENDIF
    !Get dimensions 
    dims(1)=mxGetM(pa)
    dims(2)=mxGetN(pa)
	CALL mxDestroyArray(pa)
    !Close matfile
    CALL CloseMatFile(mp,matname)

    END SUBROUTINE matgetdim    


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    SUBROUTINE doublemat3d_Wvar2Matf(matrix,matname,keyw)
	! dummy variables
    REAL (8),      DIMENSION(:,:,:),   INTENT(IN)    :: matrix
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    ! local variables         
    INTEGER (4) ::       dims(3), mmax, nmax, status  
	INTEGER (8) :: pa, mp, mxCreateNumericArray, mxGetPr
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName, matPutVariable  

    !Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    !Create, write and free from memory  
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
    dims(3)=SIZE(matrix,DIM=3)
    pa=mxCreateNumericArray(3, dims,mxClassIDFromClassName('double'),0)
    CALL mxCopyReal8ToPtr(matrix, mxGetPr(pa), dims(1)*dims(2)*dims(3))
    status = matPutVariable(mp, keyw, pa)
    CALL mxDestroyArray(pa)
    !Close matfile
    CALL CloseMatFile(mp,matname)
    
    END SUBROUTINE doublemat3d_Wvar2Matf


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    SUBROUTINE singlemat3d_Wvar2Matf(matrix,matname,keyw)
	! dummy variables
    REAL (4),      DIMENSION(:,:,:),   INTENT(IN)    :: matrix
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    ! local variables         
    INTEGER (4) ::       dims(3), mmax, nmax, status  
	INTEGER (8) :: pa, mp, mxCreateNumericArray, mxGetPr
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName, matPutVariable  

    !Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    !Create, write and free from memory  
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
    dims(3)=SIZE(matrix,DIM=3)
    pa=mxCreateNumericArray(3, dims,mxClassIDFromClassName('single'),0)
    CALL mxCopyReal4ToPtr(matrix, mxGetPr(pa), dims(1)*dims(2)*dims(3))
    status = matPutVariable(mp, keyw, pa)
    CALL mxDestroyArray(pa)
   !Close matfile
    CALL CloseMatFile(mp,matname)
    
    END SUBROUTINE singlemat3d_Wvar2Matf


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    SUBROUTINE doublemat_Wvar2Matf(matrix,matname,keyw)
	! dummy variables
    REAL (8),      DIMENSION(:,:),   INTENT(IN)    :: matrix
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    ! local variables         
    INTEGER (4) ::       dims(2), mmax, nmax, status  
	INTEGER (8) :: pa, mp, mxCreateNumericArray, mxGetPr
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName, matPutVariable  

    !Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    !Create, write and free from memory  
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
    pa=mxCreateNumericArray(2, dims,mxClassIDFromClassName('double'),0)
	CALL mxCopyReal8ToPtr(matrix, mxGetPr(pa), dims(1)*dims(2))
    status = matPutVariable(mp, keyw, pa)
    CALL mxDestroyArray(pa)
    !Close matfile
    CALL CloseMatFile(mp,matname)
    
    END SUBROUTINE doublemat_Wvar2Matf


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    SUBROUTINE singlematint_Wvar2Matf(matrix,matname,keyw)
	! dummy variables
    INTEGER (4),      DIMENSION(:,:),   INTENT(IN)    :: matrix
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    ! local variables         
    INTEGER (4) ::       dims(2), mmax, nmax, status
	INTEGER (8) :: pa, mxGetPr,mxCreateNumericArray
    INTEGER (8) :: mp
	! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName, matPutVariable 

    !Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    !PRINT *,'single int matrix: ',keyw
    !Create, write and free from memory  
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
    !write(*,*) dims(1),dims(2),mxClassIDFromClassName('int32')
	pa=mxCreateNumericArray(2,dims,mxClassIDFromClassName('int32'),0)
    !write(*,*) 'createnummatrix',pa
	CALL mxCopyInteger4ToPtr(matrix, mxGetPr(pa), dims(1)*dims(2))
    !write(*,*) 'mxcopy'
    status = matPutVariable(mp, keyw, pa)
	!write(*,*) matput
    CALL mxDestroyArray(pa)
	!write(*,*) mxdestroy
    !Close matfile
    CALL CloseMatFile(mp,matname)
   
    END SUBROUTINE singlematint_Wvar2Matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    SUBROUTINE doublematint_Wvar2Matf(matrix,matname,keyw)
	! dummy variables
    INTEGER (8),      DIMENSION(:,:),   INTENT(IN)    :: matrix
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    ! local variables         
    INTEGER (4) ::       dims(2), mmax, nmax, status
	INTEGER (8) :: pa, mxGetPr,mxCreateNumericArray
    INTEGER (8) :: mp
	! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName, matPutVariable 

    !Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    !PRINT *,'single int matrix: ',keyw
    !Create, write and free from memory  
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
    !write(*,*) dims(1),dims(2),mxClassIDFromClassName('int32')
	pa=mxCreateNumericArray(2,dims,mxClassIDFromClassName('int64'),0)
    !write(*,*) 'createnummatrix',pa
	CALL mxCopyInteger4ToPtr(matrix, mxGetPr(pa), dims(1)*dims(2))
    !write(*,*) 'mxcopy'
    status = matPutVariable(mp, keyw, pa)
	!write(*,*) matput
    CALL mxDestroyArray(pa)
	!write(*,*) mxdestroy
    !Close matfile
    CALL CloseMatFile(mp,matname)
   
    END SUBROUTINE doublematint_Wvar2Matf


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    
    SUBROUTINE doublevec_Wvar2Matf(vector,matname,keyw)
	! dummy variables
    REAL (8),        DIMENSION(:),   INTENT(IN)    :: vector
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    ! local variables         
    INTEGER (4) ::       dims(2), mmax, status
	INTEGER (8) :: pa, mp, mxCreateNumericArray, mxGetPr
	! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName, matPutVariable 

	!Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    !Create, write and free from memory  
    dims(1)=SIZE(vector,DIM=1)
    dims(2)=1
    pa=mxCreateNumericArray(1, dims,mxClassIDFromClassName('double'),0)
	CALL mxCopyReal8ToPtr(vector, mxGetPr(pa), dims(1)*dims(2))
    status = matPutVariable(mp, keyw, pa)
    CALL mxDestroyArray(pa)
    !Close matfile
    CALL CloseMatFile(mp,matname)
  
    END SUBROUTINE doublevec_Wvar2Matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    SUBROUTINE doublereal_Wvar2Matf(vector,matname,keyw)
	! dummy variables
    REAL (8),   INTENT(IN)    :: vector
    CHARACTER(LEN=*)  ,   INTENT(IN)    :: matname,keyw
    ! local variables         
    INTEGER (4) ::       dims(2), mmax , status  
   	INTEGER (8) :: pa, mp, mxCreateNumericArray, mxGetPr
	! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName, matPutVariable

	!Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)
    !Create, write and free from memory  
    dims(1)=1
    dims(2)=1
    pa=mxCreateNumericArray(1, dims,mxClassIDFromClassName('double'),0)
	CALL mxCopyReal8ToPtr(vector, mxGetPr(pa), dims(1)*dims(2))
    status = matPutVariable(mp, keyw, pa)
    CALL mxDestroyArray(pa)
    !Close matfile
    CALL CloseMatFile(mp,matname)
    
    END SUBROUTINE doublereal_Wvar2Matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    SUBROUTINE singlemat_Wvar2Matf(matrix,matname,keyw)
    ! dummy variables
    REAL (4),      DIMENSION(:,:),   INTENT(IN)    :: matrix
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    ! local variables         
    INTEGER (4) ::       dims(2), status
    ! Matlab call variables
	INTEGER (8) :: pa, mp, mxGetPr,mxCreateNumericMatrix
	INTEGER (4) ::       mxClassIDFromClassName,   &
                         matPutVariable 
    
	!Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    !Create, write and free from memory  
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
	pa=mxCreateNumericMatrix(dims(1),dims(2),mxClassIDFromClassName('single'),0)
	CALL mxCopyReal4ToPtr(matrix, mxGetPr(pa), dims(1)*dims(2))
    status = matPutVariable(mp, keyw, pa)
    CALL mxDestroyArray(pa)

    !Close matfile
    CALL CloseMatFile(mp,matname)
    
    END SUBROUTINE singlemat_Wvar2Matf
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE singlevec_Wvar2Matf(vector,matname,keyw)
    ! dummy variables
    REAL (4),      DIMENSION(:)  ,   INTENT(IN)    :: vector
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    ! local variables         

	INTEGER (4) ::       dims(2), status
 	INTEGER (8) :: pa, mp, mxCreateNumericArray, mxGetPr
   ! Matlab call variables
    INTEGER (4) ::       mxCreateNumericMatrix,    &
	                     mxClassIDFromClassName,   &
                         matPutVariable

    !Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    !Create, write and free from memory  
    dims(1)=SIZE(vector,DIM=1)
    dims(2)=1
	pa=mxCreateNumericArray(2,dims,mxClassIDFromClassName('single'),0)
	CALL mxCopyReal4ToPtr(vector, mxGetPr(pa), dims(1)*dims(2))
    status = matPutVariable(mp, keyw, pa)
    CALL mxDestroyArray(pa)
    !Close matfile
    CALL CloseMatFile(mp,matname)
    
    END SUBROUTINE singlevec_Wvar2Matf
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE char_Wvar2Matf(string,matname,keyw)
    ! dummy variables
    CHARACTER(LEN=*)     ,INTENT(IN) :: string,matname,keyw
    ! local variables         
    INTEGER (4) ::       status
	INTEGER (8) :: pa, mp,mxCreateString
    ! Matlab call variables
    INTEGER (4) ::   	 matPutVariable 

	!Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    pa = mxCreateString(string)
    status = matPutVariable(mp, keyw, pa)
    if (status .ne. 0) then
         write(*,*) 'matPutVariable failed'
         stop
      end if
    CALL mxDestroyArray(pa)
    CALL CloseMatFile(mp,matname)
  
    END SUBROUTINE char_Wvar2Matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	SUBROUTINE singleint_wvar2matf(vector,matname,keyw)
    ! dummy variables
    INTEGER (4),  INTENT(IN)    :: vector
    CHARACTER(LEN=*), INTENT(IN)    :: matname,keyw
    ! local variables         
    INTEGER (4) ::       dims(2), status
   	INTEGER (8) :: pa, mp, mxCreateNumericArray, mxGetPr
	! Matlab call variables
    INTEGER (4) ::       mxCreateNumericMatrix,    &
	                     mxClassIDFromClassName,   &
                         matPutVariable 

	!Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    !Create, write and free from memory  
    dims(1)=1
    dims(2)=1
	pa=mxCreateNumericArray(2,dims,mxClassIDFromClassName('int32'),0)
	CALL mxCopyInteger4ToPtr(vector, mxGetPr(pa), dims(1)*dims(2))
    status = matPutVariable(mp, keyw, pa)
    CALL mxDestroyArray(pa)
    !Close matfile
    CALL CloseMatFile(mp,matname)

	END SUBROUTINE singleint_Wvar2Matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	SUBROUTINE singlerealscal_wvar2matf(scalar,matname,keyw)
    ! dummy variables
    REAL (4),  INTENT(IN)    :: scalar
    CHARACTER(LEN=*), INTENT(IN)    :: matname,keyw
    ! local variables         
    INTEGER (4) ::       dims(2), status
   	INTEGER (8) :: pa, mp, mxCreateNumericArray, mxGetPr
	! Matlab call variables
    INTEGER (4) ::       mxCreateNumericMatrix,    &
	                     mxClassIDFromClassName,   &
                         matPutVariable 
	
    !Check if matname exist and open new or existing file
    CALL openmatfile(mp,matname)                  
    !Create, write and free from memory  
    dims(1)=1
    dims(2)=1
	pa=mxCreateNumericArray(2,dims,mxClassIDFromClassName('single'),0)
	CALL mxCopyReal4ToPtr(scalar, mxGetPr(pa), dims(1)*dims(2))
    status = matPutVariable(mp, keyw, pa)
    CALL mxDestroyArray(pa)
    !Close matfile
    CALL CloseMatFile(mp,matname)

	END SUBROUTINE singlerealscal_Wvar2Matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE Ini_Struc2matf(matname,keyw,fieldname,indexfield)
    ! dummy variables
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw,fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield
    ! local variables         
    INTEGER(4) ::       dimfield(2) ,status 
    INTEGER(8) :: ps, mp, mxCreateStructArray
	! Matlab call variables
    INTEGER(4) ::       matGetVariable,   &
	                    matDeleteVariable
						 
    !Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)     
    ! check, initialize and make if ness. the structured array    
 	dimfield(1)=indexfield
	dimfield(2)=1
	ps=matGetVariable(mp, keyw) 
    IF (ps/=0) status=matDeleteVariable(ps,keyw)
	ps=mxCreateStructArray(2,dimfield,1,fieldname)
	CALL mxDestroyArray(ps)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE Ini_Struc2matf
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 SUBROUTINE char_wstruc2matf(string,matname,keyw,fieldname,indexfield,numberoffields)
    ! dummy variables
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: string,matname,keyw
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) :: dimfield(2),status,matPutVariable  
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString,matGetVariable
						 
    !Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)                  
    ! check and make if ness. the structured array    
    ps=matGetVariable(mp, keyw)
   	dimfield(1)=numberoffields
	dimfield(2)=1
    IF (ps==0) ps=mxCreateStructArray(2,dimfield,1,fieldname)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa = mxCreateString(string)
	CALL mxSetField(ps, indexfield, fieldname, pfa)
	status=matPutVariable(mp,keyw,ps)
    CALL mxDestroyArray(ps)
    !CALL mxDestroyArray(pf)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE char_wstruc2matf
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE doublemat_wstruc2matf(matrix,matname,keyw,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL (8),   DIMENSION(:,:)   ,   INTENT(IN)    :: matrix
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2), dimfield(2)  
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray, matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName,status,matPutVariable
						  !mxSetField
    
    !Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)                  
    ! get dimensions of matrix to write
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
    ! check and make if ness. the structured array    
    ps=matGetVariable(mp, keyw)     
   	dimfield(1)=numberoffields
	dimfield(2)=1
    IF (ps==0) ps=mxCreateStructArray(2,dimfield,1,fieldname)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(2,dims,mxClassIDFromClassName('double'),0)
    CALL mxCopyReal8ToPtr(matrix, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)
	status=matPutVariable(mp,keyw,ps)
    CALL mxDestroyArray(ps)
    !CALL mxDestroyArray(pf)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE doublemat_wstruc2matf
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE doublevec_wstruc2matf(vector,matname,keyw,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL (8),   DIMENSION(:)     ,   INTENT(IN)    :: vector
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2), dimfield(2)
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName,status,matPutVariable
						  !mxSetField

    !Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)                  
    ! get dimensions of matrix to write
    dims(1)=SIZE(vector,DIM=1)
    dims(2)=1
    ! check and make if ness. the structured array    
    ps=matGetVariable(mp, keyw)
   	dimfield(1)=numberoffields
	dimfield(2)=1
    IF (ps==0) ps=mxCreateStructArray(2,dimfield,1,fieldname)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('double'),0)
    CALL mxCopyReal8ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps,indexfield, fieldname, pfa)
	status=matPutVariable(mp,keyw,ps)
    CALL mxDestroyArray(ps)
    !CALL mxDestroyArray(pf)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE doublevec_wstruc2matf
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE singlemat_wstruc2matf(matrix,matname,keyw,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL (4),   DIMENSION(:,:)   ,   INTENT(IN)    :: matrix
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    
    ! local variables         
    INTEGER (4) ::       dims(2), dimfield(2)  
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName,status,matPutVariable
						  !mxSetField

	!Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)                  
    ! get dimensions of matrix to write
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
    ! check and make if ness. the structured array    
    ps=matGetVariable(mp, keyw)     
	dimfield(1)=numberoffields !indexfield GD veranderd 30-01-2007
	dimfield(2)=1
    IF (ps==0) ps=mxCreateStructArray(2,dimfield,1,fieldname)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(2,dims,mxClassIDFromClassName('single'),0)
    CALL mxCopyReal4ToPtr(matrix, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)
	status=matPutVariable(mp,keyw,ps)
    CALL mxDestroyArray(ps)
    !CALL mxDestroyArray(pf)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE singlemat_wstruc2matf
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE singlevec_wstruc2matf(vector,matname,keyw,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL (4),   DIMENSION(:)     ,   INTENT(IN)    :: vector
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    
    ! local variables         
    INTEGER (4) ::       dims(2), dimfield(2)  
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName,status,matPutVariable
						  !mxSetField

	!Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)                  
    ! get dimensions of matrix to write
    dims(1)=SIZE(vector,DIM=1)
    dims(2)=1
    ! check and make if ness. the structured array    
	ps=matGetVariable(mp, keyw)    
	dimfield(1)=numberoffields
	dimfield(2)=1
    IF (ps==0) ps=mxCreateStructArray(2,dimfield,1,fieldname)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('single'),0)
    CALL mxCopyReal4ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)
	status=matPutVariable(mp,keyw,ps)
    CALL mxDestroyArray(ps)
    !CALL mxDestroyArray(pf)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE singlevec_wstruc2matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE integermat_wstruc2matf(matrix,matname,keyw,fieldname,indexfield,numberoffields)
    ! dummy variables
    INTEGER (4), DIMENSION(:,:)  ,   INTENT(IN)    :: matrix
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2), dimfield(2)  
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName,status,matPutVariable
						  !mxSetField

	!Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)                  

    ! get dimensions of matrix to write
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
    ! check and make if ness. the structured array    
    ps=matGetVariable(mp, keyw)  
	dimfield(1)=numberoffields
	dimfield(2)=1
    IF (ps==0) ps=mxCreateStructArray(2,dimfield,1,fieldname)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('int32'),0)
    CALL mxCopyInteger4ToPtr(matrix, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)
	status=matPutVariable(mp,keyw,ps)
    CALL mxDestroyArray(ps)
	!CALL mxDestroyArray(pf)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE integermat_wstruc2matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE integervec_wstruc2matf(vector,matname,keyw,fieldname,indexfield,numberoffields)
    ! dummy variables
    INTEGER (2), DIMENSION(:)    ,   INTENT(IN)    :: vector
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2), dimfield(2)  
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName,status,matPutVariable
						  !mxSetField

	!Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)                  

    ! get dimensions of matrix to write
    dims(1)=SIZE(vector,DIM=1)
    dims(2)=1
    ! check and make if ness. the structured array  
	  
    ps=matGetVariable(mp, keyw)  
	dimfield(1)=numberoffields
	dimfield(2)=1
    IF (ps==0) ps=mxCreateStructArray(2,dimfield,1,fieldname)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('int16'),0)
    CALL mxCopyInteger2ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)
	status=matPutVariable(mp,keyw,ps)
    CALL mxDestroyArray(ps)
	!CALL mxDestroyArray(pf)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE integervec_wstruc2matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE integer_wstruc2matf(vector,matname,keyw,fieldname,indexfield,numberoffields)
    ! dummy variables
    INTEGER(4)                   ,   INTENT(IN)    :: vector
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2), dimfield(2)  
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName,status,matPutVariable
						  !mxSetField

	!Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)                  

    ! get dimensions of matrix to write
    dims(1)=1 !SIZE(vector,DIM=1)
    dims(2)=1
    ! check and make if ness. the structured array    
    ps=matGetVariable(mp, keyw)  
	dimfield(1)=numberoffields
	dimfield(2)=1
    IF (ps==0) ps=mxCreateStructArray(2,dimfield,1,fieldname)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('int32'),0)
    CALL mxCopyInteger4ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)
	status=matPutVariable(mp,keyw,ps)
    CALL mxDestroyArray(ps)
	!CALL mxDestroyArray(pf)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE integer_wstruc2matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE real_wstruc2matf(vector,matname,keyw,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL(4)                      ,   INTENT(IN)    :: vector
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2), dimfield(2)  
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName,status,matPutVariable
						  !mxSetField

	!Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)                  
    ! get dimensions of matrix to write
    dims(1)=1 !SIZE(vector,DIM=1)
    dims(2)=1
    ! check and make if ness. the structured array    
    ps=matGetVariable(mp, keyw)    
	dimfield(1)=numberoffields
	dimfield(2)=1
    IF (ps==0) ps=mxCreateStructArray(2,dimfield,1,fieldname)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('single'),0)
    CALL mxCopyReal4ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)
	status=matPutVariable(mp,keyw,ps)
    CALL mxDestroyArray(ps)
	!CALL mxDestroyArray(pf)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE real_wstruc2matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE doublereal_wstruc2matf(vector,matname,keyw,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL(8)                      ,   INTENT(IN)    :: vector
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2), dimfield(2)  
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName,status,matPutVariable
						  !mxSetField

	!Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)                  
    ! get dimensions of matrix to write
    dims(1)=1 !SIZE(vector,DIM=1)
    dims(2)=1
    ! check and make if ness. the structured array    
    ps=matGetVariable(mp, keyw)    
	dimfield(1)=numberoffields
	dimfield(2)=1
    IF (ps==0) ps=mxCreateStructArray(2,dimfield,1,fieldname)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('double'),0)
    CALL mxCopyReal8ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)
	status=matPutVariable(mp,keyw,ps)
    CALL mxDestroyArray(ps)
	!CALL mxDestroyArray(pf)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE doublereal_wstruc2matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE doublemat3d_rmatf2var(matrix,matname,keyw)
    IMPLICIT NONE
    ! dummy variables
    REAL (8),  DIMENSION(:,:,:), INTENT(OUT)    :: matrix
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: matname
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: keyw
    ! local variables         
    REAL (8)   , ALLOCATABLE    :: dummatrix(:,:,:)
    INTEGER (4)::      dims(3), mmax,    nmax, omax
    INTEGER (4)::      cmmax,   cnmax, comax   
    
	! local variables         	               
    ! Matlab call variables
	INTEGER (8)::      mp, pa, MxgetPr, matGetVariable,mxGetDimensions,mxGetNumberOfDimensions
    ! Matlab call variables
    INTEGER (4) ::     mxGetM, mxGetN, mxIsNumeric 

    ! Check if matname exist and open new or existing file
    CALL openmatfile_r(mp,matname)                                    
    ! Create, write and free from memory  
    pa = matGetVariable ( mp ,keyw)
    IF (pa .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable = ',keyw
       STOP
    ENDIF
    IF (mxIsNumeric(pa) .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable not numeric= ',keyw
       STOP
    ENDIF
    ! Get dimensions of mxArray, pm
    CALL mxCopyPtrToInteger4(mxGetDimensions(pa), dims,mxGetNumberOfDimensions(pa))
    CALL mxCopyPtrToReal8(MxgetPr(pa),matrix,dims(1)*dims(2)*dims(3))  
    
    !cmmax=mxGetM(pa)
    !cnmax=mxGetN(pa)
    !mmax=SIZE(matrix,DIM=1)
    !nmax=SIZE(matrix,DIM=2)
    !IF (mmax == cmmax .AND. nmax == cnmax) THEN
    !   CALL mxCopyPtrToReal8(MxgetPr(pa),matrix,mmax*nmax)  
    !ELSEIF (mmax == cnmax .AND. nmax == cmmax ) THEN
    !   ALLOCATE(dummatrix(nmax,mmax))
    !   CALL mxCopyPtrToReal8(MxgetPr(pa),dummatrix,mmax*nmax)  
    !   matrix=TRANSPOSE(dummatrix)
    !   DEALLOCATE(dummatrix)    
    !ELSE
    !   PRINT *,'ERROR READING Variable [',keyw,'], file = ', matname 
    !   PRINT *,'Dimension of m [',mmax,'] /= [',cmmax,'] in file'
    !   PRINT *,'Dimension of n [',nmax,'] /= [',cnmax,'] in file'
    !   STOP
    !ENDIF
    CALL mxDestroyArray(pa)
	!Close matfile    
	CALL CloseMatFile(mp,matname)
    RETURN

    END SUBROUTINE doublemat3d_rmatf2var
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE doublemat_rmatf2var(matrix,matname,keyw)
    IMPLICIT NONE
    ! dummy variables
    REAL (8),  DIMENSION(:,:),   INTENT(OUT)    :: matrix
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: matname
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: keyw
    ! local variables         
    REAL (8)   , ALLOCATABLE    :: dummatrix(:,:)
    INTEGER (4)::      mmax,    nmax
    INTEGER (4)::      cmmax,   cnmax   
    
	! local variables         	               
    ! Matlab call variables
	INTEGER (8)::      mp, pa, MxgetPr, matGetVariable
    ! Matlab call variables
    INTEGER (4) ::     mxGetM, mxGetN, mxIsNumeric 

    ! Check if matname exist and open new or existing file
    CALL openmatfile_r(mp,matname)                                    
    ! Create, write and free from memory  
    pa = matGetVariable ( mp ,keyw)
    IF (pa .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable = ',keyw
       STOP
    ENDIF
    IF (mxIsNumeric(pa) .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable not numeric= ',keyw
       STOP
    ENDIF
    cmmax=mxGetM(pa)
    cnmax=mxGetN(pa)
    mmax=SIZE(matrix,DIM=1)
    nmax=SIZE(matrix,DIM=2)
    IF (mmax == cmmax .AND. nmax == cnmax) THEN
       CALL mxCopyPtrToReal8(MxgetPr(pa),matrix,mmax*nmax)  
    ELSEIF (mmax == cnmax .AND. nmax == cmmax ) THEN
       ALLOCATE(dummatrix(nmax,mmax))
       CALL mxCopyPtrToReal8(MxgetPr(pa),dummatrix,mmax*nmax)  
       matrix=TRANSPOSE(dummatrix)
       DEALLOCATE(dummatrix)    
    ELSE
       PRINT *,'ERROR READING Variable [',keyw,'], file = ', matname 
       PRINT *,'Dimension of m [',mmax,'] /= [',cmmax,'] in file'
       PRINT *,'Dimension of n [',nmax,'] /= [',cnmax,'] in file'
       STOP
    ENDIF
    CALL mxDestroyArray(pa)
	!Close matfile    
	CALL CloseMatFile(mp,matname)
    RETURN

    END SUBROUTINE doublemat_rmatf2var
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE doublevec_rmatf2var(vector,matname,keyw)
    IMPLICIT NONE
    ! dummy variables
    REAL (8),  DIMENSION(:)  ,   INTENT(OUT)    :: vector
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: matname
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: keyw
    ! local variables         
    INTEGER (4)::      mmax,cmmax    
	INTEGER (8)::      mp, pa, MxgetPr, matGetVariable
    ! Matlab call variables
    INTEGER (4) ::     mxGetM, mxGetN, mxIsNumeric 
    
    ! Check if matname exist and open new or existing file
    CALL openmatfile_r(mp,matname)                                    
    ! Create, write and free from memory  
    pa = matGetVariable ( mp ,keyw)
    IF (mxIsNumeric(pa) .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable not numeric= ',keyw
       STOP
    ENDIF
    cmmax=mxGetM(pa)
    mmax=SIZE(vector)
    IF (mmax == cmmax ) THEN
       CALL mxCopyPtrToReal8(MxgetPr(pa),vector,mmax)  
    ELSE
       PRINT *,'ERROR READING Variable [',keyw,'], file = ', matname 
       PRINT *,'Dimension of m [',mmax,'] /= [',cmmax,'] in file'
       RETURN
    ENDIF
    CALL mxDestroyArray(pa)
    !Close matfile    
	CALL CloseMatFile(mp,matname)
    RETURN

    END SUBROUTINE doublevec_rmatf2var
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	SUBROUTINE doublerealscal_rmatf2var(scalar,matname,keyw)
    IMPLICIT NONE
    ! dummy variables
    REAL (8),     INTENT(OUT)    :: scalar
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: matname
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: keyw
	! local variables 
	REAL (8)   :: dum8scalar
	REAL (4)   :: dum4scalar 
	INTEGER (4)::      mmax,cmmax    
   	INTEGER (8)::      mp, pa, MxgetPr, matGetVariable
    ! Matlab call variables
    INTEGER (4) ::     mxGetM, mxGetN, mxIsNumeric 
	CHARACTER (LEN = 6)    ClassName,mxGetClassName

	! Check if matname exist and open new or existing file
    CALL openmatfile_r(mp,matname)        
	! Create, write and free from memory  
    pa = matGetVariable ( mp ,keyw)
	
	IF (mxIsNumeric(pa) .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable not numeric= ',keyw
       STOP
    ENDIF
	ClassName=mxGetClassName(pa)
    cmmax=mxGetM(pa)
    mmax=1 
    IF (mmax == cmmax ) THEN
	   IF (ClassName == 'single') THEN
           CALL mxCopyPtrToReal4(MxgetPr(pa),dum4scalar,mmax)  
		   scalar=DBLE(dum4scalar)
       ELSEIF (LLE(ClassName,'double')) THEN
	       CALL mxCopyPtrToReal8(MxgetPr(pa),scalar,mmax)  
	   ENDIF
       !CALL mxCopyPtrToReal4(MxgetPr(pa),scalar,mmax)  
    ELSE
       PRINT *,'ERROR READING Variable [',keyw,'], file = ', matname 
       PRINT *,'Dimension of scalar [',mmax,'] /= [',cmmax,'] in file'
       RETURN
    ENDIF
    CALL mxDestroyArray(pa)
    !Close matfile    
	CALL CloseMatFile(mp,matname)
    RETURN

    END SUBROUTINE doublerealscal_rmatf2var	    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE singlemat_rmatf2var(matrix,matname,keyw)
    IMPLICIT NONE
    ! dummy variables
    REAL (4),  DIMENSION(:,:),   INTENT(OUT)    :: matrix
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: matname
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: keyw
    ! local variables         
	REAL (8)   , ALLOCATABLE                    :: dum8matrix(:,:)
	REAL (4)   , ALLOCATABLE                    :: dum4matrix(:,:)
	INTEGER (4)::      mmax,    nmax
    INTEGER (4)::      cmmax,   cnmax   
   	INTEGER (8)::      mp, pa, MxgetPr, matGetVariable
	
	! Matlab call variables
    CHARACTER (LEN = 6)    ClassName,mxGetClassName
	INTEGER (4) ::     mxGetM, mxGetN, mxIsNumeric 
    ! Check if matname exist and open new or existing file
    CALL openmatfile_r(mp,matname)                                    
    ! Create, write and free from memory  
    pa = matGetVariable ( mp ,keyw)
    IF (pa .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable = ',keyw
       STOP
    ENDIF
    IF (mxIsNumeric(pa) .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable not numeric= ',keyw
       STOP
    ENDIF
    ClassName=mxGetClassName(pa)
	cmmax=mxGetM(pa)
    cnmax=mxGetN(pa)
    mmax=SIZE(matrix,DIM=1)
    nmax=SIZE(matrix,DIM=2)
    IF (mmax == cmmax .AND. nmax == cnmax) THEN
	   IF (ClassName == 'single') THEN
           CALL mxCopyPtrToReal4(MxgetPr(pa),matrix,mmax*nmax)  
       ELSEIF (LLE(ClassName,'double')) THEN
	       ALLOCATE(dum8matrix(mmax,nmax))
           CALL mxCopyPtrToReal8(MxgetPr(pa),dum8matrix,mmax*nmax)  
		   matrix=REAL(dum8matrix)
		   DEALLOCATE(dum8matrix)
	   ENDIF
    ELSEIF (mmax == cnmax .AND. nmax == cmmax ) THEN
	   IF (LLE(ClassName,'single')) THEN
	       ALLOCATE(dum4matrix(nmax,mmax))
           CALL mxCopyPtrToReal4(MxgetPr(pa),dum4matrix,mmax*nmax)  
		   matrix=TRANSPOSE(dum4matrix)
		   DEALLOCATE(dum4matrix)    
       ELSEIF (ClassName == 'double') THEN
	       ALLOCATE(dum8matrix(nmax,mmax))
           CALL mxCopyPtrToReal8(MxgetPr(pa),dum8matrix,mmax*nmax)  
	       matrix=REAL(TRANSPOSE(dum8matrix))
		   DEALLOCATE(dum8matrix)    
	   ENDIF
    ELSE
       PRINT *,'ERROR READING Variable [',keyw,'], file = ', matname 
       PRINT *,'Dimension of m [',mmax,'] /= [',cmmax,'] in file'
       PRINT *,'Dimension of n [',nmax,'] /= [',cnmax,'] in file'
       STOP
    ENDIF
    CALL mxDestroyArray(pa)
    !Close matfile    
	CALL CloseMatFile(mp,matname)
    RETURN

    END SUBROUTINE singlemat_rmatf2var
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE singlevec_rmatf2var(vector,matname,keyw)
    IMPLICIT NONE
    ! dummy variables
    REAL (4),  DIMENSION(:)  ,   INTENT(OUT)    :: vector
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: matname
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: keyw
    ! local variables 
	REAL (8)   , ALLOCATABLE                    :: dum8vector(:)
	REAL (4)   , ALLOCATABLE                    :: dum4vector(:)  
    INTEGER (4)::      mmax,cmmax    
   	INTEGER (8)::      mp, pa, MxgetPr, matGetVariable
	! Matlab call variables
    INTEGER (4) ::     mxGetM, mxGetN, mxIsNumeric 
	CHARACTER (LEN = 6)    ClassName,mxGetClassName
    
	! Check if matname exist and open new or existing file
    CALL openmatfile_r(mp,matname)     
	! Create, write and free from memory  
    pa = matGetVariable ( mp ,keyw)
	IF (mxIsNumeric(pa) .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable not numeric= ',keyw
       STOP
    ENDIF
	ClassName=mxGetClassName(pa)
	cmmax=mxGetM(pa)	
    mmax=SIZE(vector,DIM=1)
    IF (mmax == cmmax ) THEN	   
	   IF (ClassName == 'single') THEN
           CALL mxCopyPtrToReal4(MxgetPr(pa),vector,mmax)  
       ELSEIF (LLE(ClassName,'double')) THEN
	       ALLOCATE(dum8vector(mmax))
           CALL mxCopyPtrToReal8(MxgetPr(pa),dum8vector,mmax)  
		   vector=REAL(dum8vector)
		   DEALLOCATE(dum8vector)
	   ENDIF
       !CALL mxCopyPtrToReal4(MxgetPr(pa),vector,mmax)  
    ELSE
       PRINT *,'ERROR READING Single Vector Variable [',keyw,'], file = ', matname 
       PRINT *,'Dimension of m [',mmax,'] /= [',cmmax,'] in file'
       RETURN
    ENDIF
	CALL mxDestroyArray(pa)
	!Close matfile    
	CALL CloseMatFile(mp,matname)
	RETURN

    END SUBROUTINE singlevec_rmatf2var	
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE singlerealscal_rmatf2var(scalar,matname,keyw)
    IMPLICIT NONE
    ! dummy variables
    REAL (4),     INTENT(OUT)    :: scalar
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: matname
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: keyw
	! local variables 
	REAL (8)   :: dum8scalar
	REAL (4)   :: dum4scalar 
	INTEGER (4)::      mmax,cmmax    
   	INTEGER (8)::      mp, pa, MxgetPr, matGetVariable
    ! Matlab call variables
    INTEGER (4) ::     mxGetM, mxGetN, mxIsNumeric 
	CHARACTER (LEN = 6)    ClassName,mxGetClassName

	! Check if matname exist and open new or existing file
    CALL openmatfile_r(mp,matname)        
	! Create, write and free from memory  
    pa = matGetVariable ( mp ,keyw)
	
	IF (mxIsNumeric(pa) .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable not numeric= ',keyw
       STOP
    ENDIF
	ClassName=mxGetClassName(pa)
    cmmax=mxGetM(pa)
    mmax=1 
    IF (mmax == cmmax ) THEN
	   IF (ClassName == 'single') THEN
           CALL mxCopyPtrToReal4(MxgetPr(pa),scalar,mmax)  
       ELSEIF (LLE(ClassName,'double')) THEN
	       CALL mxCopyPtrToReal8(MxgetPr(pa),dum8scalar,mmax)  
		   scalar=REAL(dum8scalar)
	   ENDIF
       !CALL mxCopyPtrToReal4(MxgetPr(pa),scalar,mmax)  
    ELSE
       PRINT *,'ERROR READING Variable [',keyw,'], file = ', matname 
       PRINT *,'Dimension of scalar [',mmax,'] /= [',cmmax,'] in file'
       RETURN
    ENDIF
    CALL mxDestroyArray(pa)
    !Close matfile    
	CALL CloseMatFile(mp,matname)
    RETURN

    END SUBROUTINE singlerealscal_rmatf2var	    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE integermat_rmatf2var(matrix,matname,keyw)
    IMPLICIT NONE
    ! dummy variables
    INTEGER,   DIMENSION(:,:),  INTENT(OUT)     :: matrix
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: matname
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: keyw
    ! local variables         
    REAL (8)   , ALLOCATABLE    :: dummatrix(:,:)
    INTEGER (4) ::    mmax,    nmax
    INTEGER (4) ::    cmmax,   cnmax    
   	INTEGER (8)::      mp, pa, MxgetPr, matGetVariable
    ! Matlab call variables
    INTEGER (4) ::    mxGetM, mxGetN, mxIsNumeric 
    
    ! Check if matname exist and open new or existing file
    CALL openmatfile_r(mp,matname)                                    
    ! Create, write and free from memory  
    pa = matGetVariable ( mp ,keyw)
    IF (pa .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable = ',keyw
       STOP
    ENDIF
    IF (mxIsNumeric(pa) .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable not numeric= ',keyw
       STOP
    ENDIF
    cmmax=mxGetM(pa)
    cnmax=mxGetN(pa)
    mmax=SIZE(matrix,DIM=1)
    nmax=SIZE(matrix,DIM=2)
    IF (mmax == cmmax .AND. nmax == cnmax) THEN
       ALLOCATE(dummatrix(mmax,nmax))
       CALL mxCopyPtrToReal8(MxgetPr(pa),dummatrix,mmax*nmax)  
       matrix=IDINT(dummatrix)
       DEALLOCATE(dummatrix)
    ELSEIF (mmax == cnmax .AND. nmax == cmmax ) THEN
       ALLOCATE(dummatrix(nmax,mmax))
       CALL mxCopyPtrToReal8(MxgetPr(pa),dummatrix,mmax*nmax)  
       matrix=IDINT(TRANSPOSE(dummatrix))
       DEALLOCATE(dummatrix)    
    ELSE
       PRINT *,'ERROR READING Variable [',keyw,'], file = ', matname 
       PRINT *,'Dimension of m [',mmax,'] /= [',cmmax,'] in file'
       PRINT *,'Dimension of n [',nmax,'] /= [',cnmax,'] in file'
       STOP
    ENDIF
    CALL mxDestroyArray(pa)
    !Close matfile
    CALL CloseMatFile(mp,matname)
    RETURN

    END SUBROUTINE integermat_rmatf2var
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE integervec_rmatf2var(vector,matname,keyw)
    IMPLICIT NONE
    ! dummy variables
    INTEGER,   DIMENSION(:)  ,  INTENT(OUT)     :: vector
	
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: matname
    CHARACTER(LEN=*)         ,   INTENT(IN)     :: keyw
    ! local variables         
    REAL (8),    ALLOCATABLE                     :: dumvector(:)
	INTEGER (4) ::         mmax,cmmax    
   	INTEGER (8)::      mp, pa, MxgetPr, matGetVariable
    ! Matlab call variables
    INTEGER (4) ::   mxGetM, mxGetN, mxIsNumeric 
    
    ! Check if matname exist and open new or existing file
    CALL openmatfile_r(mp,matname)                  
    ! Create, write and free from memory  
    pa = matGetVariable ( mp ,keyw)
    IF (pa .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable = ',keyw
       STOP
    ENDIF
    IF (mxIsNumeric(pa) .eq. 0 ) THEN
       PRINT *,'Problem reading .mat file: ', matname 
       PRINT *,'Variable not numeric= ',keyw
       STOP
    ENDIF
    cmmax=mxGetM(pa)
    mmax=SIZE(vector)
    IF (mmax == cmmax ) THEN
       ALLOCATE(dumvector(mmax))
       CALL mxCopyPtrToReal8(MxgetPr(pa),dumvector,mmax)  
       vector=IDINT(dumvector)
       DEALLOCATE(dumvector)
    ELSE
       PRINT *,'ERROR READING Variable [',keyw,'], file = ', matname 
       PRINT *,'Dimension of m [',mmax,'] /= [',cmmax,'] in file'
       STOP
    ENDIF
    CALL mxDestroyArray(pa)
    !Close matfile
    CALL CloseMatFile(mp,matname)
    RETURN

    END SUBROUTINE integervec_rmatf2var

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE Ini_Struc2ps(fieldname,indexfield,ps)
    ! dummy variables
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield
    ! local variables         
    INTEGER(4) ::       dimfield(2) 
    INTEGER(8) :: ps, mp, mxCreateStructArray
	! Matlab call variables
    INTEGER(8) ::       matGetVariable,   &
	                    matDeleteVariable
						 
    !Check if matname exist and open new or existing file
	!CALL openmatfile(mp,matname)     
    ! check, initialize and make if ness. the structured array    
 	dimfield(1)=indexfield
	dimfield(2)=1
	!ps=matGetVariable(mp, keyw) 
    !IF (ps/=0) status=matDeleteVariable(ps,keyw)
	ps=mxCreateStructArray(2,dimfield,1,fieldname)
	!CALL mxDestroyArray(ps)
	!CALL CloseMatFile(mp,matname)

    END SUBROUTINE Ini_Struc2ps

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE wps2matf(matname,keyw,ps)
    ! dummy variables
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: matname,keyw
	INTEGER(8), INTENT(IN)						   :: ps
    ! local variables         
    INTEGER(8) :: mp
	! Matlab call variables
 !   INTEGER(4) ::       matGetVariable,   &
!	                    matDeleteVariable
      INTEGER(4) ::	status,matPutVariable
						 
    !Check if matname exist and open new or existing file
	CALL openmatfile(mp,matname)     
   
	status=matPutVariable(mp,keyw,ps)
	CALL mxDestroyArray(ps)
	CALL CloseMatFile(mp,matname)

    END SUBROUTINE wps2matf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 SUBROUTINE char_wstruc2ps(string,ps,fieldname,indexfield,numberoffields)
    ! dummy variables
	CHARACTER(LEN=*)             ,   INTENT(IN)    :: string
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
   	INTEGER (8) :: ps,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString,matGetVariable
						 
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa = mxCreateString(string)
	CALL mxSetField(ps, indexfield, fieldname, pfa)

    END SUBROUTINE char_wstruc2ps
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE doublemat_wstruc2ps(matrix,ps,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL (8),   DIMENSION(:,:)   ,   INTENT(IN)    :: matrix
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2)  
   	INTEGER (8) :: ps, pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray
    ! Matlab call variables
    INTEGER (4) ::       matGetVariable, mxClassIDFromClassName
						  !mxSetField

            
    ! get dimensions of matrix to write
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(2,dims,mxClassIDFromClassName('double'),0)
    CALL mxCopyReal8ToPtr(matrix, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)

    END SUBROUTINE doublemat_wstruc2ps
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE doublevec_wstruc2ps(vector,ps,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL (8),   DIMENSION(:)     ,   INTENT(IN)    :: vector
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2)
   	INTEGER (8) :: ps, pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName
						  !mxSetField
	    
    ! get dimensions of matrix to write
    dims(1)=SIZE(vector,DIM=1)
    dims(2)=1

    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('double'),0)
    CALL mxCopyReal8ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps,indexfield, fieldname, pfa)

    END SUBROUTINE doublevec_wstruc2ps
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE singlemat_wstruc2ps(matrix,ps,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL (4),   DIMENSION(:,:)   ,   INTENT(IN)    :: matrix
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    
    ! local variables         
    INTEGER (4) ::       dims(2) 
   	INTEGER (8) :: ps, pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName
						  !mxSetField
                
    ! get dimensions of matrix to write
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)
 
    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(2,dims,mxClassIDFromClassName('single'),0)
    CALL mxCopyReal4ToPtr(matrix, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)


    END SUBROUTINE singlemat_wstruc2ps
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE singlevec_wstruc2ps(vector,ps,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL (4),   DIMENSION(:)     ,   INTENT(IN)    :: vector
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    
    ! local variables         
    INTEGER (4) ::       dims(2) 
   	INTEGER (8) :: ps, pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName
						  !mxSetField

            
    ! get dimensions of matrix to write
    dims(1)=SIZE(vector,DIM=1)
    dims(2)=1

    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('single'),0)
    CALL mxCopyReal4ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)

    END SUBROUTINE singlevec_wstruc2ps

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE integermat_wstruc2ps(matrix,fieldname,indexfield,numberoffields)
    ! dummy variables
    INTEGER (4), DIMENSION(:,:)  ,   INTENT(IN)    :: matrix
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2) 
   	INTEGER (8) :: ps, pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName
						  !mxSetField

              

    ! get dimensions of matrix to write
    dims(1)=SIZE(matrix,DIM=1)
    dims(2)=SIZE(matrix,DIM=2)

    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('int32'),0)
    CALL mxCopyInteger4ToPtr(matrix, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)

    END SUBROUTINE integermat_wstruc2ps

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE integervec_wstruc2ps(vector,ps,fieldname,indexfield,numberoffields)
    ! dummy variables
    INTEGER (4), DIMENSION(:)    ,   INTENT(IN)    :: vector
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2)
   	INTEGER (8) :: ps, mp,pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName
						  !mxSetField

              

    ! get dimensions of matrix to write
    dims(1)=SIZE(vector,DIM=1)
    dims(2)=1

    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
    !pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('int16'),0)
    !CALL mxCopyInteger2ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('int32'),0)
    CALL mxCopyInteger4ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2)) ! changed GD, 13-03-2009
	CALL mxSetField(ps, indexfield, fieldname, pfa)

    END SUBROUTINE integervec_wstruc2ps

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE integer_wstruc2ps(vector,ps,fieldname,indexfield,numberoffields)
    ! dummy variables
    INTEGER(4)                   ,   INTENT(IN)    :: vector
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2)
   	INTEGER (8) :: ps, pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName
						  !mxSetField

             

    ! get dimensions of matrix to write
    dims(1)=1 !SIZE(vector,DIM=1)
    dims(2)=1

    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('int32'),0)
    CALL mxCopyInteger4ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)

    END SUBROUTINE integer_wstruc2ps

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE real_wstruc2ps(vector,ps,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL(4)                      ,   INTENT(IN)    :: vector
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2)
   	INTEGER (8) :: ps, pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName
						  !mxSetField

               
    ! get dimensions of matrix to write
    dims(1)=1 !SIZE(vector,DIM=1)
    dims(2)=1

    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('single'),0)
    CALL mxCopyReal4ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)

    END SUBROUTINE real_wstruc2ps

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE doublereal_wstruc2ps(vector,ps,fieldname,indexfield,numberoffields)
    ! dummy variables
    REAL(8)                      ,   INTENT(IN)    :: vector
    CHARACTER(LEN=*)             ,   INTENT(IN)    :: fieldname
    INTEGER(4)                   ,   INTENT(IN)    :: indexfield,numberoffields
    

    ! local variables         
    INTEGER (4) ::       dims(2) 
   	INTEGER (8) :: ps, pf,pfa, mxCreateStructArray, mxGetPr, &
	               mxGetField, mxAddField, mxCreateString, & 
				   mxCreateNumericArray,matGetVariable
    ! Matlab call variables
    INTEGER (4) ::       mxClassIDFromClassName
						  !mxSetField
             
    ! get dimensions of matrix to write
    dims(1)=1 !SIZE(vector,DIM=1)
    dims(2)=1

    pf=mxGetField(ps,indexfield, fieldname)
    IF (pf==0) pf=mxAddField(ps, fieldname)
	pfa=mxCreateNumericArray(1,dims,mxClassIDFromClassName('double'),0)
    CALL mxCopyReal8ToPtr(vector, mxGetPr(pfa), dims(1)*dims(2))
	CALL mxSetField(ps, indexfield, fieldname, pfa)

    END SUBROUTINE doublereal_wstruc2ps


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE openmatfile(mp,matname) 
	
	CHARACTER(LEN=*)      :: matname
    INTEGER (8)           :: mp
    INTEGER (8)           :: matOpen

	mp=matOpen(matname,'u')
    IF (mp .eq. 0) THEN
		mp=matOpen(matname,'w')
		IF (mp>0) THEN
! 			PRINT *, 'New File : ',matname
		ELSE
			PRINT *, 'File : ',matname
			PRINT *, 'Error opening MAT-file : ', matname
			STOP
		END IF
    ENDIF                  
    END SUBROUTINE openmatfile
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE openmatfile_r(mp,matname) 
	
	CHARACTER(LEN=*)      :: matname
    INTEGER (8)           :: mp
    INTEGER (8)           :: matOpen

	mp=matOpen(matname,'r')
    IF (mp .eq. 0) THEN
       PRINT *, 'Cannot read : ',matname
	   STOP
    ENDIF                  
    END SUBROUTINE openmatfile_r
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE CloseMatFile(mp,matname) 
	
	CHARACTER(LEN=*)      :: matname
	INTEGER (8)           :: mp
    INTEGER (4)           :: status
    INTEGER (4)	          :: matClose
	!Close matfile
    status = matClose(mp)
    IF (status .ne. 0) THEN
       PRINT *, 'Error closing MAT-file : ', matname
       STOP
    ENDIF
	!CALL mxdestroyarray(mp)
    END SUBROUTINE CloseMatFile 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SUBROUTINE mkmatname(nnum,i,charv,charn,outfilename)
!
!  Make a file name with filenumber
!
!	nnum		= number of places reserved for sequential numbers
!	i			= sequential number
!	charv		= char before number
!	charn		= char after number
!     outfilename = outputfile
!
!	mkmatname(4,45,'Hs','AVEC.mat',outfilename)
!	produces outfilename = Hs0045AVEC.mat
!	
	INTEGER            i,nnum
	CHARACTER (LEN=*)  charv,charn,outfilename
!
! process filename
!	
    SELECT CASE (nnum)	
    CASE (1)
	   !WRITE(outfilename,101) charv,i,charn
	   WRITE(outfilename,101) charv,i,charn
 	CASE (2)
	   WRITE(outfilename,102) charv,i,charn
	CASE (3)
	   WRITE(outfilename,103) charv,i,charn	   
    CASE (4)
	   WRITE(outfilename,104) charv,i,charn
	CASE (5)
	   WRITE(outfilename,105) charv,i,charn
    CASE (6)
	   WRITE(outfilename,106) charv,i,charn
    CASE DEFAULT 
	   WRITE(outfilename,107) charv,i,charn
    END SELECT
    101	FORMAT(a,i1.1,a) 
    102	FORMAT(a,i2.2,a) 
    103	FORMAT(a,i3.3,a)	
	104	FORMAT(a,i4.4,a)	
	105	FORMAT(a,i5.5,a)	
	106	FORMAT(a,i6.6,a)
	107	FORMAT(a,i7.7,a)	
	END SUBROUTINE mkmatname
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     FUNCTION DateNum(Y,M,D,hh,mm,ss) RESULT(MD)
! !
! ! CALL : Matlab_serial_time = datenum(YEAR,MONTH,DAY) or
! !        Matlab_serial_time = datenum(YEAR,MONTH,DAY,hour,min,sec) 
! !
! ! INPUT  : YEAR  INTEGER YEAR NUMBER
! ! INTEGER  MONTH INTEGER MONTH NUMBER ( 1...12 )
! !          DAY   INTEGER DAY NUMBER ( 1...31 )
! !          hh = hour
! !          mm = minutes
! !          ss = seconds
! !     output: MD: Matlab time (0 = 00-Jan-0000 )   DOUBLE PRECISION 
! !     VALID from 1900 to 2099
! !    
! 
!     INTEGER (4)  D,M,Y
!     INTEGER (4), OPTIONAL :: hh,mm,ss
!     INTEGER (4)  hhh,mmm,sss,JD
!     REAL    (8)  MD
! 
!     IF (NOT(PRESENT(hh))) THEN 
!       hhh=0 
!     ELSE 
!       hhh=hh 
!     ENDIF
!     IF (NOT(PRESENT(mm)))  THEN 
!       mmm=0 
!     ELSE 
!       mmm=mm
!     ENDIF
!     IF (NOT(PRESENT(ss)))  THEN 
!       sss=0 
!     ELSE 
!       sss=ss 
!     ENDIF
!    
!     JD = D-32075+1461*(Y+4800+(M-14)/12)/4+367*(M-2-(M-14)/12*12) &
!           /12-3*((Y+4900+(M-14)/12)/100)/4
! 
!     MD = DBLE(JD- 2415021 + 693962) +  &
!         DBLE(hhh)/24.0 + DBLE(mmm)/1440.0 +DBLE(sss)/86400.0
! 
!     END FUNCTION DateNum
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
END MODULE Module_RwMatfiles
