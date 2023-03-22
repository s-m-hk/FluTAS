# Compiler definitions
#
# GNU (cpu only)
#
PREPROC   := -cpp 
CC        := mpicc
FC        := mpifort
ifeq ($(DO_DBG),1)
  DBG       := -O0 -g -fbacktrace -ffpe-trap=invalid -Wall -Wextra -pedantic -fcheck=all -finit-real=snan
else
  OPT       := -O3 -ffast-math -march=native
endif
PRE       := #-fdefault-real-8 -fdefault-double-8

# Take all the compiler flags together
FFLAGS  := $(OPT) $(DBG) $(PRE)
DFLAGS  := -D_TIMING -D_EPA2A -D_DECOMP_X #-D_TWOD
DFLAGS  += #-D_OVERWRITE -D_EVEN # FFLAGS_2DECOMP
LDFLAGS :=

# Architecture switches
USE_NVTX = 0

# Required for FFTW
LDFLAGS   += -L${FFTW_ROOT}/lib -lfftw3 -lfftw3_threads

# Required for INIT_MONTECARLO
#GSL_LIB   += -lgsl -lgslcblas -lm -lstdc++
#GSL_INC   += 

# Required for NVTX
# NVTX_LIB   +=
