#FFLAGS= -g -traceback -fpp -mkl -openmp -O3 -mtune=native -march=native 
FFLAGS= -fpp -mkl -openmp -O3 -mtune=native -march=native 
#CCFLAGS = -g -fopenmp  -mcmodel=medium -Wall -pedantic -Ofast -funroll-loops -finline-functions -ftree-vectorize  
CCFLAGS = -fopenmp  -mcmodel=medium -Wall -pedantic -Ofast -funroll-loops -finline-functions -ftree-vectorize  

#LDFLAGS = -lstdc++ 
LDFLAGS = -lstdc++ -lgomp
LINKS=/Soft/intelf/composer_xe_2013.3.163/compiler/lib/intel64/libiomp5.a $(LDFLAGS)

FC=/Soft/intelf/composer_xe_2013.3.163/bin/intel64/ifort
CXX=/Soft/gcc/6.2.0/bin/g++

KIND_OBJS=  
F90_LEVELS_NEIG_OBJS= quicksort.o gauinv.o random.o glb.o  grid_info.o aniso_d.o graph_subs.o kdtree3.o cova_1D.o sortem.o ktsol.o exhaustive_srch.o sgs_lva_levels_neig.o Scaling-ISOMAP_landmark.o
FOR_OBJS=chknam.o  getindx.o
LAPACK_OBJS=LAPACK/eig.o LAPACK/solve.o

OBJS=   eig.o solve.o

CPP_OBJS= boostdijkstra.o boostdijkstra_capi.o Boostdijkstra.o

ALL_LEVELS_NEIG_OBJS= $(CPP_OBJS) $(FOR_OBJS) $(F90_LEVELS_NEIG_OBJS)
OUT_LEVELS_NEIG=SGS_LVA_levels_neig_opt05

#Scaling-ISOMAP_landmark.o : boostdijkstra_mod.o

default: $(KIND_OBJS) $(LAPACK_OBJS) $(ALL_LEVELS_NEIG_OBJS)
	$(FC) $(FFLAGS) $(LINKS) $(KIND_OBJS) $(OBJS) $(ALL_LEVELS_NEIG_OBJS) -o $(OUT_LEVELS_NEIG)
	cp $(OUT_LEVELS_NEIG) ../bin/

clean:
	rm *.mod *.o $(OUT_LEVELS_NEIG) ../bin/$(OUT_LEVELS_NEIG)

.SUFFIXES: .o .F90 .f90 .for .f .fpp .cpp

.F90.o : ; $(FC) -c $(FFLAGS) $*.F90
.f90.o : ; $(FC) -c $(FFLAGS) $*.f90
.for.o : ; $(FC) -c $(FFLAGS) $*.for
.f.o : ; $(FC) -c $(FFLAGS) $*.f
.fpp.o : ; $(FC) -c $(FFLAGS) $*.fpp
.cpp.o : ; ${CXX} -c ${CCFLAGS} $*.cpp
