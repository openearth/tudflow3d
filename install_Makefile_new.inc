#
#  This file is part of MUMPS 4.10.0, built on Tue May 10 12:56:32 UTC 2011
#
#Begin orderings

# NOTE that PORD is distributed within MUMPS by default. If you would like to
# use other orderings, you need to obtain the corresponding package and modify
# the variables below accordingly.
#
# Metis/ParMetis and SCOTCH/PT-SCOTCH (ver 5.1 and later) orderings are now available for MUMPS.
#

# <<< Eddy: hier staan dus alleen de parallelle libs: pt-scotch
#===================================================================================
SCOTCHDIR  = /usr/local/scotch-5.1.12
ISCOTCH    = -I$(SCOTCHDIR)/include 
LSCOTCH    = -L$(SCOTCHDIR)/lib -lptesmumps -lptscotch -lptscotcherr
#===================================================================================


#===================================================================================
LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord
#===================================================================================

# <<< Eddy: hier staan dus beide libs: parmetis en metis
#===================================================================================
LMETISDIR = /usr/local/metis/parmetis-3.2.0
IMETIS    = -I$(LMETISDIR)/include
# <<< Eddy lib bijgevoegd
# <<< Eddy omdat blijkbaar ook libmetis nodig is, heb ik die uit 
# <<< /usr/local/metis/lib naar parmetis gecopieerd.
LMETIS    = -L$(LMETISDIR)/lib -lparmetis -lmetis
#===================================================================================

#===================================================================================
# The following variables will be used in the compilation process.
# Please note that -Dptscotch and -Dparmetis imply -Dscotch and -Dmetis respectively.
# <<< Eddy Klopt deze volgorde? Niet beter -Dptscotch -Dparmetis vooraan te zetten?
 
ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)
#===================================================================================



PLAT    =
LIBEXT  = .a
OUTC    = -o
OUTF    = -o
RM = /bin/rm -f
####CC = gcc
####FC = ifort
####FL = ifort
CC = mpicc
FC = mpif90
FL = mpif90
AR = ar vr 
#RANLIB = ranlib
RANLIB  = echo


SCALAP = /usr/local/intel/compiler11_official/11.1/072/mkl/lib/em64t/libmkl_blas95_lp64.a \
         /usr/local/intel/compiler11_official/11.1/072/mkl/lib/em64t/libmkl_lapack95_lp64.a \
         /usr/local/intel/compiler11_official/11.1/072/mkl/lib/em64t/libmkl_scalapack_lp64.a -Wl,--start-group \
         /usr/local/intel/compiler11_official/11.1/072/mkl/lib/em64t/libmkl_intel_lp64.a \
         /usr/local/intel/compiler11_official/11.1/072/mkl/lib/em64t/libmkl_core.a \
         /usr/local/intel/compiler11_official/11.1/072/mkl/lib/em64t/libmkl_sequential.a -Wl,--end-group \
         /usr/local/intel/compiler11_official/11.1/072/mkl/lib/em64t/libmkl_blacs_intelmpi_lp64.a -lpthread -lm 

MPI_ROOT = /usr/local/mpich2/1.4.1p1
LIBPAR = $(SCALAP) -L${MPI_ROOT}/lib -lmpich -lpthread -lrt -i-static
INCPAR = -I/usr/local/mpich2/1.4.1p1/include/
#MPI_ROOT ==/usr/local/openmpi-1.10.1
#LIBPAR = $(SCALAP) -L${MPI_ROOT}/lib -lpthread -lrt -i-static
#INCPAR = -I/usr/local/openmpi-1.10.1/include

INCSEQ = -I$(topdir)/libseq
LIBSEQ  =  -L$(topdir)/libseq -lmpiseq

LIBBLAS = /usr/local/intel/compiler11_official/11.1/072/mkl/lib/em64t/libmkl_blas95_lp64.a

# <<< Eddy Weet niet zeker of dit moet blijven staan !!!
CDEFS   = -DAdd_

#Begin Optimized options
OPTF    = -O3 -Dintel_ -DALLOW_NON_INIT -nofor_main -I/usr/local/intel/compiler11_official/11.1/072/mkl/include/em64t/lp64 -I/usr/local/intel/compiler11_official/11.1/072/mkl/include -132 -autodouble -traceback
OPTL    = -O3 -Dintel_ -nofor_main -I/usr/local/intel/compiler11_official/11.1/072/mkl/include/em64t/lp64 -I/usr/local/intel/compiler11_official/11.1/072/mkl/include -132 -autodouble -traceback
OPTC    = -O3 -132 -autodouble -traceback

INCS = $(INCPAR)
LIBS = $(LIBBLAS) $(LIBPAR)
LIBSEQNEEDED =
