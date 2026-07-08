#!/bin/bash
# source /opt/intel/oneapi/setvars.sh
ulimit -s unlimited
make clear
make all
chmod +x SOCh.exe
mpirun -np 1 ./SOCh.exe > out.log