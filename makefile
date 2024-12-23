F77 = mpif90 -autodouble -132 -traceback -O3 -heap-arrays 10 -qmkl -fstack-protector-all -march=x86-64-v3 # this line is 30% faster than line below with -Ofast on H7-16vcpu, on H7-4pcpu the performance is same for both lines, -march=x86-64-v4 does not work on 4pcpu older cores and is not tested 
#F77 = mpif90 -autodouble -132 -traceback -Ofast -heap-arrays 10 -qmkl -fstack-protector-all -march=x86-64-v3 # default option for Deltares runs 2024 #-axcore-avx2,avx ## dit was nodig, maar zonder malloc en heap-arrays opties 3% sneller
#F77 = mpif90 -autodouble -132 -traceback -O3 -heap-arrays 10 -qmkl -fstack-protector-all # default option for Deltares runs 2024 #-axcore-avx2,avx ## dit was nodig, maar zonder malloc en heap-arrays opties 3% sneller
#F77 = mpif90 -autodouble -132 -traceback -O0 -heap-arrays 10 -qmkl -fstack-protector-all # -axcore-avx2,avx ## dit was nodig, maar zonder malloc en heap-arrays opties 3% sneller
#F77 = mpiifort -autodouble -132 -traceback -O0 -heap-arrays 10 -qmkl # -axcore-avx2,avx ## dit was nodig, maar zonder malloc en heap-arrays opties 3% sneller
#F77 = mpif90 -autodouble -132 -traceback -O3 -qopt-malloc-options=4 -heap-arrays 10 # -axcore-avx2,avx ## dit was nodig, maar zonder malloc en heap-arrays opties 3% sneller
#F77 = mpif90 -autodouble -132 -traceback -O3 -axavx -qopt-malloc-options=4 -heap-arrays 10 # -axavx is even snel als axcor-avx2,avx
#F77 = mpif90 -autodouble -132 -traceback -O3 -qopt-malloc-options=4 -heap-arrays 10 -axcore-avx2,-axavx  #ax optie checkt bij runnen of betreffende processor optimization ondersteund--> 3% sneller op XEONv4 proc
## -qopt-malloc-options=4 -heap-arrays 10 are really needed otherwise spooky errors for some runs; --> ifort 2023 -qopt-malloc-options=4 not supported anymore
#F77 = mpif90 -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132 -fbacktrace -Ofast

FLAGS =  

## debug options two lines below:
#FLAGS = -diag-disable:7712 -diag-disable:5194 -g -fpe0 -ftrapuv -warn all,noexternal -debug extended -check all -check:noarg_temp_created
#FLAGS = -diag-disable:7712 -diag-disable:5194 -g -fpe0 -ftrapuv -warn all,noexternal -debug extended -check all -check:noarg_temp_created
#F77 = mpif90 -autodouble -132 -traceback -O0 -qopt-malloc-options=4 -heap-arrays 10 # -axcore-avx2,avx
#FLAGS = -diag-disable:7712 -diag-disable:5194 -g -fpe0 -ftrapuv -warn -debug extended -check all 
#FLAGS = -diag-disable:7712 -diag-disable:5194 -g -fpe0 -ftrapuv -check all 
#FLAGS =  -fpe3 -ftrapuv -debug extended -check bounds


#netcdfdir = /opt/netcdf/v4.6.2_v4.4.4_intel_18.0.3
#LIBS2 = -L/opt/netcdf/v4.6.2_v4.4.4_intel_18.0.3/lib -lnetcdff -lnetcdf  
#netcdfdir = /opt/apps/netcdf/v4.7.4_v4.5.3_intel18.0.3/
#LIBS2 = -L/opt/apps/netcdf/v4.7.4_v4.5.3_intel18.0.3/lib -lnetcdff -lnetcdf

netcdfdir = /opt/apps/netcdf/4.9.2_4.6.1_intel2023.1.0
LIBS2 = -L/opt/apps/netcdf/4.9.2_4.6.1_intel2023.1.0/lib -lnetcdff -lnetcdf

#netcdfdir= /opt/apps/netcdf-cxx/4.3.1_gcc12.2.0/
#LIBS2 = -L/opt/apps/netcdf-cxx/4.3.1_gcc12.2.0/lib -lnetcdff -lnetcdf

#netcdfdir = /opt/apps/netcdf-fortran/4.6.0_gcc12.2.0
#LIBS2 = -L/opt/apps/netcdf-fortran/4.6.0_gcc12.2.0/lib -lnetcdff -lnetcdf

#topdir = /usr/local/mumps/4.10.0
#relict from past when mumps was linked, now only needed for sparse_matrix.f

### below lines needed to use mumps solver:
###libdir = $(topdir)/lib
###include ./install_Makefile_new.inc
###LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)
###LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)
### end mumps part....

#######MKLPATH=$MKLROOT/lib/em64t
#####MKLROOT= /opt/intel/18.0.3/mkl
#####MKLPATH=$MKLROOT/lib/intel64_lin/
#####MKLINCLUDE=$MKLROOT/include
#MKLPATH=/opt/apps/intelmkl/2023.1.0/mkl/2023.1.0/lib 
#MKLINCLUDE=/opt/apps/intelmkl/2023.1.0/mkl/2023.1.0/include 
#LIBMKL = -L$MKLPATH -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread


FC = $(F77) $(FLAGS)
RM = rm -f

# commandos to obtain github info in fortran program and trace version of code for each simulation:
VR1 = $(shell echo "      gitversion = '`git rev-parse --short HEAD`'" > version.inc)
VR2 = $(shell echo "      url = 'https://github.com/openearth/tudflow3d'" >> version.inc )
VR3 = $(shell echo "      date_make = '`date`'" >> version.inc)

#VR1 = 
#VR2 = 
#VR3 = 


PROGRAM = TUDflow3d_vg_h7AL8.exe

SRCS    = error_functions.f nlist.f extra_functions.f sediment.f advec.f advec_compact4knikker.f advec_cds6.f advec_cds4.f advec_hybrid6.f advec_hybrid4.f advec_c4a6.f diff_compact.f bound.f init.f main.f solve.f turbulence_models.f vfft.f output.f rheology.f

OBJS    = error_functions.o nlist.o extra_functions.o sediment.o advec.o advec_compact4knikker.o advec_cds6.o advec_cds4.o advec_hybrid6.o advec_hybrid4.o advec_c4a6.o diff_compact.o bound.o init.o main.o solve.o turbulence_models.o vfft.o output.o rheology.o 

all: $(PROGRAM)

#$(PROGRAM): $(OBJS)
#	$(F77) $(FLAGS) -o $(PROGRAM) $(OBJS) $(LIBS2) $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBMETIS5) $(LIBOTHERS) $(LIBMKL)
#$(PROGRAM): $(OBJS)
#	$(F77) $(FLAGS) -o $(PROGRAM) $(OBJS) $(LIBS2) $(LIBMKL)	
$(PROGRAM): $(OBJS)
	$(F77) $(FLAGS) -o $(PROGRAM) $(OBJS) $(LIBS2)
	
###########################
main.o: main.f makefile
	$(F77) $(FLAGS) -c main.f

vfft.o: vfft.f makefile
	$(F77) $(FLAGS) -c vfft.f

sediment.o: sediment.f makefile
	$(F77) $(FLAGS) -c -I$(netcdfdir)/include sediment.f

error_functions.o: error_functions.f makefile
	$(F77) $(FLAGS) -c error_functions.f

turbulence_models.o: turbulence_models.f makefile
	$(F77) $(FLAGS) -c turbulence_models.f

nlist.o: nlist.f makefile	
	$(F77) $(FLAGS) -c nlist.f

advec.o: advec.f makefile	
	$(F77) $(FLAGS) -c advec.f

advec_compact4knikker.o: advec_compact4knikker.f makefile	
	$(F77) $(FLAGS) -c advec_compact4knikker.f

advec_cds6.o: advec_cds6.f makefile	
	$(F77) $(FLAGS) -c advec_cds6.f

advec_cds4.o: advec_cds4.f makefile	
	$(F77) $(FLAGS) -c advec_cds4.f

advec_hybrid6.o: advec_hybrid6.f makefile	
	$(F77) $(FLAGS) -c advec_hybrid6.f
	
advec_hybrid4.o: advec_hybrid4.f makefile	
	$(F77) $(FLAGS) -c advec_hybrid4.f	
	
advec_a4c6.o: advec_c4a6.f makefile	
	$(F77) $(FLAGS) -c advec_c4a6.f	

diff_compact.o: diff_compact.f makefile	
	$(F77) $(FLAGS) -c diff_compact.f

output.o: output.f makefile	
	$(F77) $(FLAGS) -c -I$(netcdfdir)/include output.f

extra_functions.o: extra_functions.f makefile
	$(F77) $(FLAGS) -c -I$(netcdfdir)/include extra_functions.f

bound.o: bound.f makefile	
	$(F77) $(FLAGS) -I$(netcdfdir)/include -c bound.f

init.o: init.f makefile	
	$(F77) $(FLAGS) -c -I$(netcdfdir)/include init.f

#sparse_matrix.o: sparse_matrix.f makefile	
#	$(F77) $(FLAGS) -c -I$(topdir)/include sparse_matrix.f
	
solve.o: solve.f makefile	
##	$(F77) $(FLAGS) -I$MKLINCLUDE -c solve.f
	$(F77) $(FLAGS) -c solve.f
	
rheology.o: rheology.f makefile	
	$(F77) $(FLAGS) -c rheology.f

version: $(VR1) $(VR2) $(VR3) 

clean:
	$(RM) a.out core *.mod $(PROGRAM) $(OBJS)

