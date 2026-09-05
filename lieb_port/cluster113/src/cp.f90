      program CPMCH
      use cpmc

!-------------------------HoKinPara-------------------------!
      include 'mpif.h'
      !use mpi
!-------------------------HoKinPara-------------------------!
      integer :: clock_start, clock_end, clock_rate
      real :: elapsed_time
!-------------------------HoKinPara-------------------------!


      integer::iblk,istp,iw,th,record
      real(sp)::ranGenyu1,rescale
      real(sp)::tt1,tt2
      integer::valu(8),icpu1,icpu2,icpu3 

!-------------------------HoKinPara-------------------------!
      integer:: clusterSize,cpuRank
      integer iTag,iErr,iProc,ISTATUS(MPI_STATUS_SIZE)
       
 
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,clusterSize,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,cpuRank,ierr) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      noOfProc=clusterSize
      myID=cpuRank
      write(*,*) 'in mpi ncpu=',cpuRank,'clusterSize=',clusterSize,'ierror',ierr
      write(*,*) 'myID=',myID
 
      call system_clock(clock_start, clock_rate) 
!-------------------------HoKinPara-------------------------!
  
      do i =1,cpuRank+1
        call date_and_time(values=valu)
      enddo
      ISEED=(-valu(8)*(cpuRank+1))

      do i =1,cpuRank+1
      ranGenyu1=ranGen(ISEED)
      enddo

      measl=0
      wbkwlkr=0.0
      nforward=1 
      open(5,file='time.dat',status='unknown')

     !!~~~~~INITALIZE~~~~~!!
         
      call SetUp
      call InitPop
  
!!~~~~~START QMC~~~~~!!
       !test delivery and reception
       call Scatter 
       call Gather  
      write(*,*) '53myID=',myID
      if (myID==0) then
       
      call InitEnergy
  
      
      call preprog         
      end if 

      call date_and_time(values=valu)
      write(5,*) 'Init relax:',valu(2),valu(3),valu(5),valu(6),valu(7)
      valu2=valu(2);valu3=valu(3);valu5=valu(5);valu6=valu(6);valu7=valu(7)  
        
!!~~~~~RELAXATION PHASE~~~~~!!
      
      do iblk=1,nblkeq
        do istp=1,nblkstps
        call Scatter
         
        call Step(istp)

        if(mod(istp,itvlorth) == 0) then
          call Stblz
        end if
       
       call Gather 

        if(mod(istp,itvlpceq) == 0) then   
        if (myID==0) then
          call Comb
          wgtwlkr=(real(NWLKRS)/sum(wgtwlkr))*wgtwlkr 
        end if 
        end if
  
        end do
 
        write(5,*) 'equilibrate:',iblk
        flush(5)
        if (myID==0) then
          write(*,*) 'equilibrate block',iblk,'of',nblkeq
          flush(6)
        end if

      end do 

      call date_and_time(values=valu)
      write(5,*) 'End relax and Init growth:',valu(3),valu(5),valu(6),valu(7)

      write(*,*) '100myID=',myID,etrial 

!!~~~ADJUST ETRIAL BY GROWTH ESTIMATE~~~!!

      
      if(nblkgr /= 0) then
         if (myID==0) then
         call EstEtrial
         
        end if 
      end if

      call MPI_Bcast(etrial,1,MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      write(*,*) '108myID=',myID,etrial 

      call date_and_time(values=valu)
      write(5,*) 'End growth and start Runmeas:',valu(3),valu(5),valu(6),valu(7)

!!~~~ MEASUREMENT PHASE ~~~!!

      !call Gather

      measl=1

    if(nblk /= 0) then

      if (myID==0) then  
      call InitRunMeas
      end if

      do iblk=1,nblk
        if (myID==0) then
        call InitBlkMeas
        end if

        do istp=1,nblkstps

          mstep=mod(istp,itvl_m)

          
          if( mstep == 1 ) then     ! save current wave functions
           !call Gather
           if (myID==0) then
           
           !! Adjust total weight to be real(NWLKS):
           rescale=real(NWLKRS)/sum(wgtwlkr)
           do th=1,NWLKRS
            wgtwlkr(th)=wgtwlkr(th)*rescale
           end do 

           call savephi 

           end if
           !call Scatter
          end if

         !call Gather 
!---------------------------------------------!  
  !      if (myID==0) then
   !      call Step1(istp)  
   !       if(mod(istp,itvlorth) == 0) then
   !         call Stblz1 
    !      end if 
    !    end if  

     !------------------para--------------------!  
        call Scatter
         
        call Step(istp)

        if(mod(istp,itvlorth) == 0) then
          call Stblz
        end if
       
       call Gather  
!-----------------------------------------------!        
        
          if(mod(istp,itvlmeas) == 0) then ! Mixed energy
          !call Gather
            if (myID==0) then
            call StpMeas 
            end if
          end if
 
           
          if( mstep == 0 ) then     ! Correlations 
          !call Gather
            if (myID==0) then
            call correl
            end if 
          end if

          if(mod(istp,itvlpc) == 0) then
           !call Gather
            if (myID==0) then
            call Comb
            end if
           !call Scatter
          end if

         end do                      ! Loop istp (within block)
         
         !call Gather
         if (myID==0) then
         call BlkMeas(iblk)
         write(5,*) 'back step:',iblk
         end if 
      end do                         ! Loop iblk
      
      if (myID==0) then
      call RunMeas
      end if

    end if                         ! Measurement Phase; nblk != 0

      call date_and_time(values=valu)
      write(5,*) 'End Runmeas:',valu(3),valu(5),valu(6),valu(7)
!-------------------------Time-------------------------! 
      call system_clock(clock_end)
      elapsed_time = real(clock_end - clock_start) / real(clock_rate)
      

      if (myid == 0) then
      write(*,*) 'Elapsed time for rank 0 (seconds):', elapsed_time,'second'
      write(5,*) 'Elapsed time for rank 0 (seconds):', elapsed_time,'second'
      valusec=(valu(2)-valu2)*30*86400+(valu(3)-valu3)*86400+(valu(5)-valu5)*3600+(valu(6)-valu6)*60+(valu(7)-valu7)
      write(*,*) 'time',valu(2)-valu2,'month',valu(3)-valu3,'day',valu(5)-valu5,'hour',valu(6)-valu6,'min',valu(7)-valu7,'second'
      write(5,*) 'time',valu(2)-valu2,'month',valu(3)-valu3,'day',valu(5)-valu5,'hour',valu(6)-valu6,'min',valu(7)-valu7,'second'
      write(*,*) 'time_sec', valusec,'second'
      write(5,*) 'time_sec', valusec,'second'
      end if

      close(5)
!-------------------------HoKinPara-------------------------!
      call MPI_FINALIZE(ierr)
!-------------------------HoKinPara-------------------------!  
      stop
      END program CPMCH



