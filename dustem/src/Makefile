
# Options g95
# FC = g95
#- Standard ----------------------------
# FFLAGS = -O2
#- develop -----------------------------
# FFLAGS = -g3 -fno-second-underscore -Wall -Wextra -fbounds-check -ftrace=full -Wno=136,161,102 -finteger=99999999 -flogical=none -freal=nan -fpointer=invalid -Wline-truncation -Wprecision-loss
# FFLAGS = -g3 -fno-second-underscore -Wall -Wextra -fbounds-check -ftrace=full -Wno=112,136,161,102,165 -finteger=99999999 -flogical=none -freal=nan -fpointer=invalid


# Options gfortran
  FC = gfortran
#- Standard ----------------------------
# FFLAGS = -O2 -fno-second-underscore
# --------------------------------------
# FFLAGS = -O2 -march=i686
#- develop -----------------------------
# FFLAGS = -g3 -O0 -fno-second-underscore -Wall -Wextra -fbounds-check -fimplicit-none -std=f2003 -pedantic -ffpe-trap=invalid,zero,overflow,denormal -gdwarf-2 #-fall-intrinsics -Wconversion -Wno-tabs -Waliasing
# FFLAGS = -g3 -O0 -fno-second-underscore -Wall -Wextra -fbounds-check -fimplicit-none -std=f2003 -pedantic -ffpe-trap=zero,overflow -gdwarf-2 -Wconversion -Waliasing
 FFLAGS = -O3 -fno-second-underscore -ftree-vectorize -falign-loops=16 -msse2 -ffast-math -fassociative-math	

# Options ifort
# FC = ifort
#- Standard ----------------------------
# FFLAGS = -O  -fp-model strict
# FFLAGS = -O2 -parallel
#- Intel Core2Duo ----------------------
# FFLAGS = -O2 -fp-model strict -axSSSE3 -m32 -mtune=core2 -parallel
# FFLAGS = -O2 -axSSSE3 -m32 -mtune=core2 -parallel
# LIBS = -lmkl_lapack -lmkl_intel -lmkl_core -lmkl_intel_thread -liomp5 -lpthread
#- Debug -------------------------------
# FFLAGS = -g3
# FFLAGS = -g3 -check all -warn all -traceback -debug all -diag-enable warn -fp-stack-check -mp1
# FFLAGS = -g3 -check all -warn all -traceback -debug all -ftrapuv -fpe0 -diag-enable warn -fp-stack-check -mp1 -stand f03

.SUFFIXES: .f90

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90

OBJ = \
      DM_constants.o \
      DM_utility.o \
      DM_dtls.o \
      DM_get_qext.o \
      DM_get_tdist.o \
      DM_zdist.o \
      DM_spin.o \
      DM_compute.o \
      DM_inout.o 

ALL: dustem
	echo "   That's it!"

dustem :: $(OBJ) DM_dustem.o
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) DM_dustem.o -o dustem

clean:
	/bin/rm -f *.o *.mod dustem *genmod* fort.*
