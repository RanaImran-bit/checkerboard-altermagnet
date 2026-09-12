# ── GNU Fortran + OpenMPI + system LAPACK/BLAS ──────────────────────────────
# Use on clusters WITHOUT Intel OneAPI  (e.g. 113new)
# Requires:  module load openmpi/4.1.6   (done inside job_cluster.sh)
# Does NOT touch the original Makefile (Intel/MKL) used on 250-255.

FC       = /usr/local/openmpi-4.1.6/bin/mpif90
# -ffree-line-length-none: parameter.f90 has declaration lines past gfortran's
# 132-character free-form limit, which ifort on 250-255 accepts silently.
OPTIMIZE = -O2 -mcmodel=large -ffree-line-length-none
OBJECT1  = parameter.o jiekou.o mc2.o
OBJECT2  = libuph.o cpPara.o
LIBS     = -llapack -lblas -lpthread -lm

all: CPMC.exe

mc2.o: mc2duph.f90
	$(FC) $(OPTIMIZE) -g -fbacktrace -o mc2.o -c mc2duph.f90

cpOut.o: cpOut.f90
	$(FC) $(OPTIMIZE) -g -fbacktrace -o cpOut.o -c cpOut.f90

cp.o: cp.f90
	$(FC) $(OPTIMIZE) -g -fbacktrace -o cp.o -c cp.f90

libuph.o: libuph.f90
	$(FC) $(OPTIMIZE) -g -fbacktrace -o libuph.o -c libuph.f90

cpPara.o: cpPara.f90
	$(FC) $(OPTIMIZE) -g -fbacktrace -o cpPara.o -c cpPara.f90

jiekou.o: jiekou.f90
	$(FC) $(OPTIMIZE) -g -fbacktrace -o jiekou.o -c jiekou.f90

parameter.o: parameter.f90
	$(FC) $(OPTIMIZE) -g -fbacktrace -o parameter.o -c parameter.f90

CPMC.exe: $(OBJECT1) $(OBJECT2) cp.o cpOut.o
	$(FC) $(OPTIMIZE) -g -fbacktrace -o CPMC.exe \
	    cp.o cpOut.o $(OBJECT1) $(OBJECT2) $(LIBS)
	rm -f ./*.o ./*.mod

clear:
	rm -f ./*.o ./*.mod ./*.exe
