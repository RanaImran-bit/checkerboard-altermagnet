# ====== CPQMC Makefile (two-orbital layout) ======
# Compiler & flags
FC        = mpiifort
OPTIMIZE  = -O2 -g -traceback -qopenmp -check bounds -check all -fpe0 -traceback -init=snan,arrays -traceback -mcmodel=large
FFLAGS    = $(OPTIMIZE)
MKLROOT   = /opt/intel/oneapi/mkl/latest/
# Intel MKL (LP64) convenience: -qmkl is enough for mpiifort/ifort.
MKLFLAGS  = -mkl
INCLUDES  = -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LDFLAGS   = $(MKLFLAGS) $(INCLUDES) -fpe0 -ftrapuv
LIBS      =

# Sources (new structure)
SRCS = \
  cpCore.f90 \
  cpMain.f90 \
  cpOut.f90 \
  cpParallel.f90 \
  cpParameter.f90 \
  meas.f90 \
  t.f90 \
  Vee.f90 \
  Vph.f90 \
  z_energy.f90 \
  z_jiekou.f90 \
  z_libs.f90 \
  z_overlap.f90 \
  z_stabilize.f90 \
  z_walkers.f90

# A conservative build order (modules first)
OBJS = \
  cpParameter.o \
  z_libs.o \
  z_jiekou.o \
  t.o \
  Vee.o \
  Vph.o \
  z_overlap.o \
  z_stabilize.o \
  z_walkers.o \
  z_energy.o \
  meas.o \
  cpOut.o \
  cpParallel.o \
  cpCore.o \
  cpMain.o

# Default target
all: CPMC.exe

# Link
CPMC.exe: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
	@rm -f *.o *.mod

# Pattern rule to compile any .f90 into .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(MKLFLAGS)

# Phony helpers
.PHONY: clean veryclean print

print:
	@echo "FC       = $(FC)"
	@echo "FFLAGS   = $(FFLAGS)"
	@echo "LDFLAGS  = $(LDFLAGS)"
	@echo "LIBS     = $(LIBS)"
	@echo "SRCS     = $(SRCS)"
	@echo "OBJS     = $(OBJS)"

clean:
	@rm -f *.o *.mod

veryclean: clean
	@rm -f CPMC.exe *.exe out.dat time.dat
	@rm -rf dir-kVals
	@rm -f corr_00*
	@rm -rf dir-rVals
