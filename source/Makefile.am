bin_PROGRAMS  = LPsolver_2D_NoField
AM_CPPFLAGS   = $(FFTW_CFLAGS) 
AM_CPPFLAGS  += -I$(OPENBLAS_INC)
LIBS          = $(BLAS_LIBS) $(MKL_LIBS)

h_sources   = EntropyCalculations.h LP_ompi.h NegativityChecks.h collisionRoutines_1.h \
	      MarginalCreation.h SetInit_1.h conservationRoutines.h \
	      MomentCalculations.h advection_1.h \
              PoissonFunctions.h PoissonPrint.h PoissonVariables.h

cpp_sources = EntropyCalculations.cpp LP_ompi.cpp NegativityChecks.cpp  collisionRoutines_1.cpp \
              MarginalCreation.cpp SetInit_1.cpp conservationRoutines.cpp \
	      MomentCalculations.cpp advection_1.cpp \
              PoissonFunctions.cpp PoissonPrint.cpp PoissonVariables.cpp

LPsolver_2D_NoField_SOURCES = $(cpp_sources) $(h_sources)
