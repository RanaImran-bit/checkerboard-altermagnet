OPTIMIZE = -O2 
OBJECT1 = HCBSS.o HCBSSOutput.o AppHCBSS.o 
MKLROOT=/opt/intel/oneapi/mkl/2024.0
all: SOCh.exe
 
HCBSS.o: 	BSS.f90
	mpif90 $(OPTIMIZE) -o HCBSS.o -c BSS.f90

HCBSSOutput.o: 	BSSOutput.f90
	mpif90 $(OPTIMIZE) -o HCBSSOutput.o -c BSSOutput.f90	

AppHCBSS.o:	AppBSS.f90
	mpif90 $(OPTIMIZE) -o AppHCBSS.o -c AppBSS.f90

SOCh.exe:  $(OBJECT1)
	mpif90 $(OPTIMIZE) -o SOCh.exe $(OBJECT1) $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include 
	rm -f ./*.o ./*.mod
	
clear:
	rm -f ./*.o ./*.mod

