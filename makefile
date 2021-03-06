F77 = mpif90 -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132 -fbacktrace -Ofast -ftree-loop-distribution -finit-local-zero 
#F77 = mpif90 -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-132 -fbacktrace -O0 

FLAGS =  
LIBS = -L/data/data_lwit/netcdf_libs/lib -lnetcdff 

FC = $(F77) $(FLAGS)
RM = rm -f

# commandos to obtain svn info in fortran program:
VR1 = $(shell echo "      svnversion = '`svnversion`'" > version.inc)
VR2 = $(shell echo "      svnurl = '`svn info | grep URL`'" >> version.inc )

PROGRAM = Dflow3d_vg.exe

SRCS    = error_functions.f nlist.f extra_functions.f sediment.f advec.f advec_compact4knikker.f advec_cds6.f advec_hybrid4.f advec_hybrid6.f advec_c4a6.f diff_compact.f bound.f init.f main.f solve.f turbulence_models.f vfft.f output.f

OBJS    = error_functions.o nlist.o extra_functions.o sediment.o advec.o advec_compact4knikker.o advec_cds6.o advec_hybrid4.o advec_hybrid6.o advec_c4a6.o diff_compact.o bound.o init.o main.o solve.o output.o turbulence_models.o vfft.o

all: $(VR1) $(VR2) $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(F77) $(FLAGS) -o $(PROGRAM) $(OBJS) $(LIBS)

###########################
main.o: main.f makefile
	$(F77) $(FLAGS) -c main.f

vfft.o: vfft.f makefile
	$(F77) $(FLAGS) -c vfft.f

sediment.o: sediment.f makefile
	$(F77) $(FLAGS) -c sediment.f

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

advec_hybrid4.o: advec_hybrid4.f makefile	
	$(F77) $(FLAGS) -c advec_hybrid4.f

advec_hybrid6.o: advec_hybrid6.f makefile	
	$(F77) $(FLAGS) -c advec_hybrid6.f

advec_a4c6.o: advec_c4a6.f makefile	
	$(F77) $(FLAGS) -c advec_c4a6.f

diff_compact.o: diff_compact.f makefile	
	$(F77) $(FLAGS) -c diff_compact.f

output.o: output.f makefile	
	$(F77) $(FLAGS) -c -I/data/data_lwit/netcdf_libs/include output.f

extra_functions.o: extra_functions.f makefile
	$(F77) $(FLAGS) -c -I/data/data_lwit/netcdf_libs/include extra_functions.f

bound.o: bound.f makefile	
	$(F77) $(FLAGS) -c bound.f

init.o: init.f makefile	
	$(F77) $(FLAGS) -c -I/data/data_lwit/netcdf_libs/include init.f

solve.o: solve.f makefile	
	$(F77) $(FLAGS) -c solve.f

clean:
	$(RM) *.mod *.o $(PROGRAM) $(OBJS)

version: $(VR1) $(VR2) 
