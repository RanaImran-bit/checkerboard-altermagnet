program AppBSS
use mbss
use link,only: DeallocatingVariables,ranGen
use link,only: cnfmake,cnfmeas
include 'mpif.h'
!implicit none
double precision::ranGenyu1,negs
integer::valu(8)
integer::i,j,nwarms,nmeas,nswps
character(12)::date1,clock1
real:: tStartSth,tEndSth
real:: tStartBin,tEndBin,binTime,runTime,runTimePerStep
!      fileList
!      100 real/momentum space definition
!      3,4 nonInteractive green function output for HFQMC
!      500 input
!      600 run overall details
!      700 <n>, other quantities
!      800 all calculated quantity
!      900 energy spectrum, initialization state
!       7  output for analysis

!initialize the time spent in each subroutine
!MPI feature
call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD,nsize,ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD,ocpu,ierror)
!ocpu=0
!nsize=1

DO iido=0,0
!finding random seeds from clock
call date_and_time(date1,clock1)
do i =1,ocpu+1
    call date_and_time(values=valu)
enddo
!   negative number to initialize generator
ISEED=(-valu(8)*(ocpu+1))
!	  iseed=123644332
!	  write(*,*) 'random No fixed'
do i =1,ocpu+1
    ranGenyu1=ranGen(ISEED)
enddo
!do iround=0,mpiround,40
!ncpu=ocpu+iround
ncpu=ocpu+1
write(900+ncpu+iido*1000,*) 'in mpi ncpu=',ncpu,'nsize=',nsize,'ierror',ierror
!!... prepare system for measurement
call sysdef (nwarms,nmeas,nswps)
write(900+ncpu+iido*1000,*) 'OK > sysdef:'
write(*,*) 'OK > sysdef:ncpu=',ncpu
!-------------------------------!
!testingHowLongTheRunWillLast   !
!-------------------------------!
!        write(900+ncpu,*)'Date', date1,'Time', clock1
!        call cpu_time(tStartSth)
!        do i = 1, 5
!            call cnfmake (accTmp,negTmp)
!        enddo
!        call cpu_time(tEndSth)
!        write(*,*) nwarms*(tEndSth-tStartSth)/5/60/60,'hrs needed for eqm'
!        write(*,*) (tEndSth-tStartSth)/5,' s per warmup step'
!
!        call cpu_time(tStartSth)
!        do i = 1, 5
!         call cnfmake (accTmp,negTmp)
!		 call cnfmake (accTmp,negTmp) !change line 696 about times of cnfmake flipping spin
!         call cnfmeas
!        enddo
!        call cpu_time(tEndSth)
!        write(*,*) nswps*(tEndSth-tStartSth)/5/60/60,'hrs needed for each Bin'
!        write(*,*) nswps*nmeas*(tEndSth-tStartSth)/5/60/60,'hrs needed for all Bin'
!        write(*,*) (tEndSth-tStartSth)/5,' s per measuring step'

!equilibrium runs
call cpu_time(tStartSth)
call sysequil (nwarms)
call cpu_time(tEndSth)
write(900+ncpu+iido*1000,*) 'OK > sysequil:',nwarms,'times'
write(*,*) 'OK > sysequil:ncpu=',ncpu
! Calculate the elapsed time in seconds: time for equilibrium
write(*,*)'EqmTime', (tEndSth-tStartSth)

!initialization
ram=0.0d0
rsign=0.0d0
ramc=cmplx(0.0d0,0.0d0)
! ... construct statistic averages ram(rd) for nmeas 'Monte Carlo runs'
! ... statistic average comes from repeatedly measure a physical quality
call cpu_time(tStartSth)
do i = 1, nmeas
    call cpu_time(tStartBin)
    call sysmeas (nswps,negs)
    call cpu_time(tEndBin)
    binTime = tEndBin-tStartBin
    call binOutput(ncpu+iido*1000,i)
!finding square and average value for different quanitites in each bin
    do j = 1, MEAS4
    ! mah is the result after nswps(sweeps) sampling in HS field in one Monte Carlo run
        ram(j,fAVG) = ram(j,fAVG) + mah(j)
        ram(j,fRMS) = ram(j,fRMS) + mah(j)*mah(j)

    end do

    rsign(1)=rsign(1)+mah(MEAS4+1)
    rsign(2)=rsign(2)+mah(MEAS4+1)*mah(MEAS4+1)
    do j = 1, MEAC3
        ramc(j,fAVG) = ramc(j,fAVG) + mahc(j)
        ramc(j,fRMS) = ramc(j,fRMS) + mahc(j)*mahc(j)
    end do
    write(*,*)'timing meas',i," ", binTime
end do
call cpu_time(tEndSth)
runTime = tEndSth - tStartSth
write(*,*)'MeasTime', runTime
runTimePerStep = runTime / nmeas / nswps
write(*,*)'timing meas', runTime

!!... analyze measurements
call measana(nmeas)

!!... output results
call system('mkdir -p PairKSpace')
call system('mkdir -p PairRSpace')

call analysisOutput
write(200+ncpu+iido*1000,*) 'Square'
write(200+ncpu+iido*1000,*) 'beta'
write(200+ncpu+iido*1000,*) beta
write(200+ncpu+iido*1000,*) 'NT'
write(200+ncpu+iido*1000,*) NT
write(200+ncpu+iido*1000,*) 'dt'
write(200+ncpu+iido*1000,*) dt
write(200+ncpu+iido*1000,*) 'NXNY'
write(200+ncpu+iido*1000,*) NX,NY
write(200+ncpu+iido*1000,*) 'mu'
write(200+ncpu+iido*1000,*) mu
write(200+ncpu+iido*1000,*) 'U'
write(200+ncpu+iido*1000,*) ue


call MPI_FINALIZE(ierror)

!release memory
call DeallocatingVariables
END DO
endprogram AppBSS
