! Standalone replica of the tk construction in Checkerboard_Model/mc2duph.f90
! (GetMtrxs: coordinates, PBC, NN block with tam=0, checkerboard diagonal block).
! Compares nothing itself - just dumps tk(:,:,1) and tk(:,:,2) for Python to check.
program tk_check
  implicit none
  integer, parameter :: sp = kind(1.0d0)
  integer :: lx, ly, nsites_cu, nsites
  integer :: ixx, iyy, i, j, k, itmpx, itmpy
  integer, allocatable :: ixv(:), iyv(:)
  integer, allocatable :: iposit(:,:)
  real(sp), allocatable :: tk(:,:,:)
  real(sp) :: t0, t1, t2, tam, ttp, ttn, tp, tm
  character(len=64) :: arg

  call get_command_argument(1, arg); read(arg,*) lx
  ly = lx
  call get_command_argument(2, arg); read(arg,*) t0
  call get_command_argument(3, arg); read(arg,*) t1
  call get_command_argument(4, arg); read(arg,*) t2
  tam = 0.0_sp

  nsites_cu = lx*ly
  nsites = nsites_cu
  allocate(ixv(nsites), iyv(nsites))
  allocate(iposit(-lx:2*lx, -ly:2*ly))
  allocate(tk(nsites, nsites, 2))

  !! Set up coordinates (verbatim logic from mc2duph.f90)
  do ixx=1,lx
    do iyy=1,ly
      i=(ixx-1)*ly+iyy
      ixv(i)=ixx-1
      iyv(i)=iyy-1
      iposit( ixv(i),iyv(i) )=i
    end do
  end do

  !! Periodic boundary conditions (verbatim logic)
  do iyy=-ly,2*ly
    do ixx=-lx,2*lx
      itmpx=mod(ixx+lx,lx)
      itmpy=mod(iyy+ly,ly)
      iposit(ixx,iyy)=iposit(itmpx,itmpy)
    end do
  end do

  tk = 0.0_sp

  !! NN block (verbatim, with tam=0)
  ttp=t0+tam
  ttn=t0-tam
  do i=1,nsites_cu
    tk( i, iposit( ixv(i)+1,iyv(i) ),1 )=ttn
    tk( i, iposit( ixv(i)-1,iyv(i) ),1 )=ttn
    tk( i, iposit( ixv(i),iyv(i)+1 ),1 )=ttp
    tk( i, iposit( ixv(i),iyv(i)-1 ),1 )=ttp
    tk( i, iposit( ixv(i)+1,iyv(i) ),2 )=ttp
    tk( i, iposit( ixv(i)-1,iyv(i) ),2 )=ttp
    tk( i, iposit( ixv(i),iyv(i)+1 ),2 )=ttn
    tk( i, iposit( ixv(i),iyv(i)-1 ),2 )=ttn
  end do

  !! CHECKERBOARD diagonal block (verbatim from the patch; arXiv 2605.11669:
  !! every site keeps BOTH diagonals, strong orientation alternates A/B)
  tp = t1 + t2
  tm = t1 - t2
  do i=1,nsites_cu
    if ( mod( ixv(i)+iyv(i), 2 ) == 0 ) then
      tk( i, iposit( ixv(i)+1,iyv(i)+1 ),1 )=tp
      tk( i, iposit( ixv(i)-1,iyv(i)-1 ),1 )=tp
      tk( i, iposit( ixv(i)+1,iyv(i)-1 ),1 )=tm
      tk( i, iposit( ixv(i)-1,iyv(i)+1 ),1 )=tm
      tk( i, iposit( ixv(i)+1,iyv(i)+1 ),2 )=tp
      tk( i, iposit( ixv(i)-1,iyv(i)-1 ),2 )=tp
      tk( i, iposit( ixv(i)+1,iyv(i)-1 ),2 )=tm
      tk( i, iposit( ixv(i)-1,iyv(i)+1 ),2 )=tm
    else
      tk( i, iposit( ixv(i)+1,iyv(i)+1 ),1 )=tm
      tk( i, iposit( ixv(i)-1,iyv(i)-1 ),1 )=tm
      tk( i, iposit( ixv(i)+1,iyv(i)-1 ),1 )=tp
      tk( i, iposit( ixv(i)-1,iyv(i)+1 ),1 )=tp
      tk( i, iposit( ixv(i)+1,iyv(i)+1 ),2 )=tm
      tk( i, iposit( ixv(i)-1,iyv(i)-1 ),2 )=tm
      tk( i, iposit( ixv(i)+1,iyv(i)-1 ),2 )=tp
      tk( i, iposit( ixv(i)-1,iyv(i)+1 ),2 )=tp
    end if
  end do

  !! Dump: site coords then both spin matrices
  open(11, file='tk_dump.txt', status='replace')
  write(11,*) lx, ly
  do i=1,nsites
    write(11,*) ixv(i), iyv(i)
  end do
  do k=1,2
    do i=1,nsites
      do j=1,nsites
        if (abs(tk(i,j,k)) > 1.0e-12_sp) write(11,'(3i6,f20.12)') k, i, j, tk(i,j,k)
      end do
    end do
  end do
  close(11)
  write(*,*) 'dumped tk for lx=',lx,' t0=',t0,' t1=',t1,' t2=',t2
end program tk_check
