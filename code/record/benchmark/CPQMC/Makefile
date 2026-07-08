OPTIMIZE = -O2 -mcmodel=large
PARAMETER1=
OBJECT1 = parameter.o jiekou.o mc2.o 
OBJECT2 = libuph.o cpPara.o 
MKLROOT=/opt/intel/oneapi/mkl/2024.0/

all: CPMC.exe
mc2.o: mc2duph.f90
	mpiifort $(OPTIMIZE) -g -traceback -o mc2.o -c mc2duph.f90 -qmkl

cpOut.o: cpOut.f90
	mpiifort $(OPTIMIZE) -g -traceback -o cpOut.o -c cpOut.f90 -qmkl
	
cp.o: 	cp.f90
	mpiifort $(OPTIMIZE) -g -traceback -o cp.o -c cp.f90 -qmkl 

libuph.o: libuph.f90 
	mpiifort $(OPTIMIZE) -g -traceback -o libuph.o -c libuph.f90  -qmkl

cpPara.o: cpPara.f90
	mpiifort $(OPTIMIZE) -g -traceback -o cpPara.o -c cpPara.f90 -qmkl

jiekou.o: jiekou.f90
	mpiifort $(OPTIMIZE) -g -traceback -o jiekou.o -c jiekou.f90 -qmkl

parameter.o: parameter.f90
	mpiifort $(OPTIMIZE) -g -traceback -o parameter.o -c parameter.f90 -qmkl

CPMC.exe:  $(OBJECT1) $(OBJECT2) cp.o cpOut.o 
	mpiifort $(OPTIMIZE) -g -traceback -o CPMC.exe cp.o cpOut.o $(OBJECT1) $(OBJECT2) $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include 

	rm -f ./*.o ./*.mod
	
clear:
	rm -f ./*.o ./*.mod ./*.exe

