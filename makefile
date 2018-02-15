# Directories
DIR=$(PWD)/
EXECDIR=$(DIR)
OBJDIR=$(DIR)
SRCDIR=$(DIR)

# GNU C compiler
MPICC=mpicxx 
CC=icc
CPP=icpc

#Intel MPI compiler for C++
MPICC=mpicxx

# Compiler flags: crashed when compiling with -O0
CFLAGS = -O2 -openmp  -I$(TACC_FFTW3_INC) -I$(TACC_MKL_INC)
FFTFLAGS = -L$(TACC_FFTW3_LIB) -lfftw3_threads -lfftw3 -lpthread -lm 
MKLFLAGS = -Wl,-rpath,$(TACC_MKL_LIB) -L$(TACC_MKL_LIB) -Wl,--start-group -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -Wl,--end-group -liomp5 -lpthread

# FEniCS Related Flags
BOOST_DIR = /opt/apps/intel17/boost/1.64/include/
EIGEN_DIR = $(WORK)/FEniCS/eigen-install/include/eigen3/
DOLFIN_INC_DIR=$(WORK)/FEniCS/dolfin-install/include/
MSHR_INC_DIR=$(WORK)/FEniCS/mshr-install/include/
UFC_DIR = $(WORK)/FEniCS/ffc-install/lib/python2.7/site-packages/ffc/backends/ufc/

FENICS_FLAGS=$(WORK)/FEniCS/dolfin-install/lib/libdolfin.so $(WORK)/FEniCS/mshr-install/lib/libmshr.so -I$(DOLFIN_INC_DIR) -I$(MSHR_INC_DIR) -I$(BOOST_DIR) -I$(EIGEN_DIR) -I$(UFC_DIR) -std=c++11

# Command definition
RM=rm -f

sources_FPL = FPL_main.cpp 
objects_FPL= $(sources_FPL:.c=.o)

sources_LP = LP_ompi_testscaling.cpp 
objects_LP= $(sources_LP:.c=.o)

sources_LP1 = LP_ompi_1.cpp 
objects_LP1= $(sources_LP1:.c=.o)

sources_LP2 = LP_ompi_2.cpp 
objects_LP2= $(sources_LP2:.c=.o)

sources_wt = WeightGenerator_mpi.cpp 
objects_wt= $(sources_wt:.c=.o)


FPL: $(objects_FPL)
	@echo "Building FPL solver"
	$(CPP) $(objects_FPL) $(CFLAGS)  -o fpl.out $(MKLFLAGS) $(FFTFLAGS)

LP: $(objects_LP)
	@echo "Building Landau-Poisson solver"
	@$(MPICC) $(CFLAGS) $(objects_LP) -o $(EXECDIR)lp_testscaling.out  $(FFTFLAGS) $(MKLFLAGS)

LP1: $(objects_LP1)
	@echo "Building Landau-Poisson solver1"
	@$(MPICC) $(CFLAGS) $(FENICS_FLAGS) $(objects_LP1) -o $(EXECDIR)LPsolver_nu0_2D_FEniCS_Test.out  $(FFTFLAGS) $(MKLFLAGS)

LP2: $(objects_LP2)
	@echo "Building Landau-Poisson solver1"
	@$(MPICC) $(CFLAGS) $(objects_LP2) -o $(EXECDIR)LPsolver_nu005_N8_Etest.out  $(FFTFLAGS) $(MKLFLAGS)

wts: $(objects_wt)
	@echo "Building mpi weights"
	@$(MPICC) -O2 -openmp $(objects_wt) -o $(EXECDIR)weight.out  

clean:
	$(RM) $(OBJDIR)*.o 
	$(RM) $(EXECDIR)*.out

# icpc -O2 -openmp LP_main.cpp -I$TACC_FFTW3_INC -I$TACC_MKL_INC -L$TACC_FFTW3_LIB -lfftw3_threads -lfftw3 -lpthread -lm -Wl,-rpath,$TACC_MKL_LIB -L$TACC_MKL_LIB -Wl,--start-group -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -Wl,--end-group -liomp5 -lpthread
