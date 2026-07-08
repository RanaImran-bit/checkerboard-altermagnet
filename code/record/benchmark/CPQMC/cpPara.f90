!----------HoKin ParallelComputing--------------!
!For defining the index 
subroutine StartEnd(NPROC,IS1,IS2,ISTART,IEND,ICOUNT)
integer NPROC,IS1,IS2,ISTART(0:1000),IEND(0:1000),ICOUNT(0:1000),i,dataLength
dataLength=(IS2-IS1+1)/NPROC
lastDataLength=(IS2-IS1+1)-dataLength*NPROC
do i=0,NPROC-1
    if (i==NPROC-1) then
    ISTART(i)=1+i*dataLength
    IEND(i)=(i+1)*dataLength+lastDataLength
    else
    ISTART(i)=1+i*dataLength
    IEND(i)=(i+1)*dataLength
    endif
    ICOUNT(i)=IEND(i)-ISTART(i)+1
enddo
endsubroutine StartEnd



!--------------------------HoKin ParallelComputing--------------!
subroutine Scatter
use cpmc
use mpi 
integer iErr,ISTATUS(MPI_STATUS_SIZE) 
integer::kexpV_s1(NSTATES,M_MAX,NWLKRS),keph_s1(NSTATES,M_MAX,NWLKRS)
real(sp)::phi_up1(NSTATES,NUP,NWLKRS),phi_dn1(NSTATES,NDN,NWLKRS)
real(sp)::phi_cup1(NSTATES,NUP,NWLKRS),phi_cdn1(NSTATES,NDN,NWLKRS)
real(sp)::wgtwlkr1(NWLKRS),sgn1(NWLKRS),ovlpDET1(NWLKRS),wgt_c1(NWLKRS) 


      call MPI_BARRIER(MPI_COMM_WORLD,ierror) 

      call MPI_Scatter(kexpV_s(1,1,iwStart(myid)),iwCount(myid)*sizekexpV,MPI_INTEGER,kexpV_s1(1,1,iwStart(myid)),iwCount(myid)*sizekexpV,MPI_INTEGER,0,MPI_COMM_WORLD,iErr)

      call MPI_Scatter(keph_s(1,1,iwStart(myid)),iwCount(myid)*sizekexpV,MPI_INTEGER,keph_s1(1,1,iwStart(myid)),iwCount(myid)*sizekexpV,MPI_INTEGER,0,MPI_COMM_WORLD,iErr)  
       
      call MPI_Scatter(phi_up(1,1,iwStart(myid)),iwCount(myid)*sizePhiUp,MPI_REAL8,phi_up1(1,1,iwStart(myid)),iwCount(myid)*sizePhiUp,MPI_REAL8,0,MPI_COMM_WORLD,iErr)
      
      call MPI_Scatter(phi_dn(1,1,iwStart(myid)),iwCount(myid)*sizePhiDn,MPI_REAL8,phi_dn1(1,1,iwStart(myid)),iwCount(myid)*sizePhiDn,MPI_REAL8,0,MPI_COMM_WORLD,iErr)
       
      call MPI_Scatter(phi_cup(1,1,iwStart(myid)),iwCount(myid)*sizePhiUp,MPI_REAL8,phi_cup1(1,1,iwStart(myid)),iwCount(myid)*sizePhiUp,MPI_REAL8,0,MPI_COMM_WORLD,iErr)
      
      call MPI_Scatter(phi_cdn(1,1,iwStart(myid)),iwCount(myid)*sizePhiDn,MPI_REAL8,phi_cdn1(1,1,iwStart(myid)),iwCount(myid)*sizePhiDn,MPI_REAL8,0,MPI_COMM_WORLD,iErr)
       
      call MPI_Scatter(wgtwlkr(iwStart(myid)),iwCount(myid),MPI_REAL8,wgtwlkr1(iwStart(myid)),iwCount(myid),MPI_REAL8,0,MPI_COMM_WORLD,iErr)
       
      call MPI_Scatter(sgn(iwStart(myid)),iwCount(myid),MPI_REAL8,sgn1(iwStart(myid)),iwCount(myid),MPI_REAL8,0,MPI_COMM_WORLD,iErr)
      
      call MPI_Scatter(ovlpDET(iwStart(myid)),iwCount(myid),MPI_REAL8,ovlpDET1(iwStart(myid)),iwCount(myid),MPI_REAL8,0,MPI_COMM_WORLD,iErr)
      
      call MPI_Scatter(wgt_c(iwStart(myid)),iwCount(myid),MPI_REAL8,wgt_c1(iwStart(myid)),iwCount(myid),MPI_REAL8,0,MPI_COMM_WORLD,iErr)
       
      !call MPI_Scatter(etrial,1,MPI_REAL8,etrial1,1,MPI_REAL8,0,MPI_COMM_WORLD,iErr)
  
             

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       kexpV_s(:,:,:)=kexpV_s1(:,:,:)
       keph_s(:,:,:)=keph_s1(:,:,:)
       phi_up(:,:,:)=phi_up1(:,:,:)
       phi_dn(:,:,:)=phi_dn1(:,:,:)
       phi_cup(:,:,:)=phi_cup1(:,:,:)
       phi_cdn(:,:,:)=phi_cdn1(:,:,:)
       wgtwlkr(:)=wgtwlkr1(:)
       sgn(:)=sgn1(:)
       ovlpDET(:)=ovlpDET1(:)
       wgt_c(:)=wgt_c1(:)
       !etrial=etrial1

      
        
endsubroutine Scatter



subroutine Gather
use cpmc
use mpi 
integer iErr,ISTATUS(MPI_STATUS_SIZE) 
integer::kexpV_s1(NSTATES,M_MAX,NWLKRS),keph_s1(NSTATES,M_MAX,NWLKRS)
real(sp)::phi_up1(NSTATES,NUP,NWLKRS),phi_dn1(NSTATES,NDN,NWLKRS)
real(sp)::phi_cup1(NSTATES,NUP,NWLKRS),phi_cdn1(NSTATES,NDN,NWLKRS)
real(sp)::wgtwlkr1(NWLKRS),sgn1(NWLKRS),ovlpDET1(NWLKRS),wgt_c1(NWLKRS) 

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)


      call MPI_Gather(kexpV_s(1,1,iwStart(myid)),iwCount(myid)*sizekexpV,MPI_INTEGER,kexpV_s1(1,1,iwStart(myid)),iwCount(myid)*sizekexpV,MPI_INTEGER,0,MPI_COMM_WORLD,iErr) 

      call MPI_Gather(keph_s(1,1,iwStart(myid)),iwCount(myid)*sizekexpV,MPI_INTEGER,keph_s1(1,1,iwStart(myid)),iwCount(myid)*sizekexpV,MPI_INTEGER,0,MPI_COMM_WORLD,iErr)  
       
      call MPI_Gather(phi_up(1,1,iwStart(myid)),iwCount(myid)*sizePhiUp,MPI_REAL8,phi_up1(1,1,iwStart(myid)),iwCount(myid)*sizePhiUp,MPI_REAL8,0,MPI_COMM_WORLD,iErr)
      
      call MPI_Gather(phi_dn(1,1,iwStart(myid)),iwCount(myid)*sizePhiDn,MPI_REAL8,phi_dn1(1,1,iwStart(myid)),iwCount(myid)*sizePhiDn,MPI_REAL8,0,MPI_COMM_WORLD,iErr)
       
      call MPI_Gather(phi_cup(1,1,iwStart(myid)),iwCount(myid)*sizePhiUp,MPI_REAL8,phi_cup1(1,1,iwStart(myid)),iwCount(myid)*sizePhiUp,MPI_REAL8,0,MPI_COMM_WORLD,iErr)
      
      call MPI_Gather(phi_cdn(1,1,iwStart(myid)),iwCount(myid)*sizePhiDn,MPI_REAL8,phi_cdn1(1,1,iwStart(myid)),iwCount(myid)*sizePhiDn,MPI_REAL8,0,MPI_COMM_WORLD,iErr)
       
      call MPI_Gather(wgtwlkr(iwStart(myid)),iwCount(myid),MPI_REAL8,wgtwlkr1(iwStart(myid)),iwCount(myid),MPI_REAL8,0,MPI_COMM_WORLD,iErr)
       
      call MPI_Gather(sgn(iwStart(myid)),iwCount(myid),MPI_REAL8,sgn1(iwStart(myid)),iwCount(myid),MPI_REAL8,0,MPI_COMM_WORLD,iErr)
      
      call MPI_Gather(ovlpDET(iwStart(myid)),iwCount(myid),MPI_REAL8,ovlpDET1(iwStart(myid)),iwCount(myid),MPI_REAL8,0,MPI_COMM_WORLD,iErr)
      
      call MPI_Gather(wgt_c(iwStart(myid)),iwCount(myid),MPI_REAL8,wgt_c1(iwStart(myid)),iwCount(myid),MPI_REAL8,0,MPI_COMM_WORLD,iErr)
       
      !call MPI_Gather(etrial,1,MPI_REAL8,etrial1,1,MPI_REAL8,0,MPI_COMM_WORLD,iErr)
  
      

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       kexpV_s(:,:,:)=kexpV_s1(:,:,:)
       keph_s(:,:,:)=keph_s1(:,:,:)
       phi_up(:,:,:)=phi_up1(:,:,:)
       phi_dn(:,:,:)=phi_dn1(:,:,:)
       phi_cup(:,:,:)=phi_cup1(:,:,:)
       phi_cdn(:,:,:)=phi_cdn1(:,:,:)
       wgtwlkr(:)=wgtwlkr1(:)
       sgn(:)=sgn1(:)
       ovlpDET(:)=ovlpDET1(:)
       wgt_c(:)=wgt_c1(:)
       !etrial=etrial1

endsubroutine Gather 
!----------HoKin ParallelComputing--------------!