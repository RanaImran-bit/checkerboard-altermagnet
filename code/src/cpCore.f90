!====================FIRST PART==========================!

SUBROUTINE SetUp
use cpmc
use jiekou,only:InitPop,Initlatt,InitV,Initpair,Initk,MkExpT,MkExpV, initpara
integer :: ibas
external :: StartEnd
character*60 infile
!!!!!!READ INPUT PARAMETERS FROM in.dat.
open(8,file='in.dat',status='old')
read(8,*) nblkeq,nblkgr,nblk,nblkstps
read(8,*) itvlpceq,itvlpc,itvlorth,itvlmeas,itvl_m
read(8,*) deltau,etrial
read(8,*) infile
read(8,*) t1,t2,t3,t4,SOC
read(8,*) uxx,uxy,v
read(8,*) g_ph,w_ph
read(8,*) (cwfbas(ibas),ibas=1,NWFBAS)
read(8,*) (cwibas(ibas),ibas=1,NIW)
close(8)
hdeltau=half*deltau
call allocate_vars
call initout
!!!!!INPUT AND GENERATE VARIOUS MATRICIES
!     construct exp(-dt*T) and exp(-dt*V)
!      .. deltau/2 is used to symmetrize the kernel:
!!!!!exp(-0.5*deltau*T)*exp(-deltau*V)*exp(-0.5*deltau*T).
! call GetMtrxs
call Initlatt
call InitV
call Initpair
call Initk
call MkExpT(HALF*deltau)
! call MkInitOvlps
call MkExpV      ! Make e^V for each spin
call initpara
call InitPop

! return
END subroutine SetUp

SUBROUTINE Step(istp, mode)
  use cpmc
  use jiekou,only:HalfT
  integer::i,j,k,l,m,n,iselect,chan, mode
  integer::iw,istp,mx,iw1,iw2
  real(sp)::alpha_u
  ! external :: Fullph,VHubb,Vxy,Vznn
  iselect=mod(istp,itvlorth)
  if (mode==1) then
    iw1=iwStart(myID)
    iw2=iwEnd(myID)
  else
    iw1=1
    iw2=NWLKRS
  endif
!------------------HoKinPara----------------!
  do 100 iw=iw1,iw2
!------------------HoKinPara----------------!
    !! SKIP WALKERS WITH wgt=0
      if(wgtwlkr(iw) <= ZERO) go to 100
    !! APPLY growth control factor
      wgtwlkr(iw)=exp(deltau*etrial)*wgtwlkr(iw)

    !! APPLY EXP(-0.5*DelTau*T):
      if(iselect==1) then
      call HalfT(iw,1,istp)
      else
      call halfT(iw,2,istp)
      end if

      if(wgtwlkr(iw) <= ZERO) go to 100

    !! APPLY exp(-DelTau*V):
    if (abs(uxx)>0.01) then
      call VHubb(iw)
      if (wgtwlkr(iw) <= ZERO) go to 100
    endif
    if (abs(uxy)>0.01) then
      call Vxy(iw)
      if (wgtwlkr(iw) <= ZERO) go to 100
    endif
    if (abs(v)>0.01) then
      call Vznn(iw)
      if (wgtwlkr(iw) <= ZERO) go to 100
    endif

    !-------Zhongbing-phonon fields------!
      if(g_ph*w_ph>0.001) then
      call Fullph(iw)
      if(wgtwlkr(iw) <= ZERO) go to 100
      end if
    !------------------------------------!

      mx=mod(istp,itvl_m)
      if(mx==0) mx=itvl_m

      if(measl==1) then
        call rec_fields(iw,mx)
      end if
    !! APPLY exp(-0.5*DelTau*T):
      if(iselect==0) call HalfT(iw,1,istp)
100    continue
  return
END subroutine Step


SUBROUTINE StepMeas
  use cpmc
  use jiekou,only:ranGen
  implicit none
  integer::i,j,k,l,m,n
  integer::iw,ibas
  real(sp)::ekin,ecoul,ovlpNEW
  real(sp)::ewlkr(NWLKRS),ebas(NWFBAS),detwlkr(NWLKRS)
  real(sp)::detpbas(NWFBAS)
  external :: dgemm,dgedi,dgefa
  !! g = <a^{\dagger}a> = I-G = [R(LR)^{-1}L]/det(LR)
  estptop = ZERO
  estpbtm = ZERO
  ewlkr   = ZERO
  detwlkr = ZERO

  do iw = 1, NWLKRS
    if (wgtwlkr(iw) /= ZERO) then
      ! detpbas(:) = ZERO    
      ebas(:)    = ZERO
      do ibas = 1, NWFBAS
        call calgf(iw, ibas, 2, detpbas, NWFBAS)                
        call calenergy(ibas, ebas, NWFBAS)      
        ! ebas(ibas) = ebas(ibas) * detpbas(ibas)
      end do
      detwlkr(iw) = wgtwlkr(iw) * sum(detpbas) / ovlpDET(iw)
      ewlkr(iw)   = wgtwlkr(iw) * sum(ebas)    / ovlpDET(iw)
    end if
  end do

  estptop = estptop + sum(ewlkr)
  estpbtm = estpbtm + sum(detwlkr)
  nstp    = nstp + 1
END SUBROUTINE StepMeas

! subroutine rec_fields(iw,mx)
!   use cpmc
!   integer:: i, chan, iw, mx

!   ! Hubbard U
!   if (uxx>0) then
!     do i=1,nsites
!       if (ising_v(i,1) == 0) then
!           write(*,'(A)') '!!! ================================================='
!           write(*,'(A)') '!!! FATAL ERROR in rec_fields: ising_v is 0!'
!           write(*,'(A,4I8)') '!!! Walker(iw), Step(mx), Site(i), Chan:', iw, mx, i, 1
!           write(*,'(A)') '!!! ================================================='
!           stop "rec_fields: Uninitialized ising_v detected."
!       endif
!       kexpV_s(i,mx,1,iw)=ising_v(i,1)
!     enddo
!   endif

!   ! U'
!   if (abs(uxy)>0.01) then
!     do chan=2,5
!       do i=1,lxy
!         if (ising_v(i,chan) == 0) then
!             write(*,'(A)') '!!! ================================================='
!             write(*,'(A)') '!!! FATAL ERROR in rec_fields: ising_v is 0!'
!             write(*,'(A,4I8)') '!!! Walker(iw), Step(mx), Site(i), Chan:', iw, mx, i, chan
!             write(*,'(A)') '!!! ================================================='
!             stop "rec_fields: Uninitialized ising_v detected."
!         endif
!         kexpV_s(i,mx,chan,iw)=ising_v(i,chan)
!       enddo
!     enddo
!   endif

!   ! nn V
!   if (abs(uxy) > 0.01) then
!     do chan = 6, 37
!       if ( MOD((chan-6)/8, 2) == 0 ) then
!         do i = 1, lxy
!           if (ising_v(i,chan) == 0) then
!               write(*,'(A)') '!!! ================================================='
!               write(*,'(A)') '!!! FATAL ERROR in rec_fields: ising_v is 0!'
!               write(*,'(A,4I8)') '!!! Walker(iw), Step(mx), Site(i), Chan:', iw, mx, i, chan
!               write(*,'(A)') '!!! ================================================='
!               stop "rec_fields: Uninitialized ising_v detected."
!           endif
!           kexpV_s(i, mx, chan, iw) = ising_v(i, chan)
!         end do
!       else
!         do i = lxy+1, nsites
!           if (ising_v(i,chan) == 0) then
!               write(*,'(A)') '!!! ================================================='
!               write(*,'(A)') '!!! FATAL ERROR in rec_fields: ising_v is 0!'
!               write(*,'(A,4I8)') '!!! Walker(iw), Step(mx), Site(i), Chan:', iw, mx, i, chan
!               write(*,'(A)') '!!! ================================================='
!               stop "rec_fields: Uninitialized ising_v detected."
!           endif
!           kexpV_s(i, mx, chan, iw) = ising_v(i, chan)
!         end do
!       end if
!     end do
!   endif

!   if(g_ph*w_ph>0.001) then
!     do i=1,nsites          ! Saved phonon fields
!       keph_s(i,mx,iw)=ising_sum(i)
!     end do
!   endif
! endsubroutine rec_fields

subroutine rec_fields(iw,mx)
use cpmc
integer:: i, chan, iw, mx
  ! Hubbard U
  if (uxx>0) then
    do i=1,nsites
      kexpV_s(i,mx,1,iw)=ising_v(i,1)
    enddo
  endif
  ! U'
  if (abs(uxy)>0.01) then
    do chan=2,5
      do i=1,lxy
        kexpV_s(i,mx,chan,iw)=ising_v(i,chan)
      enddo
    enddo
  endif
  !nn V
  if (abs(v) > 0.01) then          ! was abs(uxy): this block records the v (neighbour) fields
    do i=1, lxy
      do j1=1,2
        do ixy=0,1 ! xy orbital of i site
          do jxy=0,1 ! xy orbital of j site
            
            isite = i + ixy*lxy 

            chan =  6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 0 )
            if (ising_v(isite, chan) /= 0) then
                kexpV_s(isite, mx, chan, iw) = ising_v(isite, chan)
            endif

            chan = 6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 1 )
            if (ising_v(isite, chan) /= 0) then
                kexpV_s(isite, mx, chan, iw) = ising_v(isite, chan)
            endif

            chan = 6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 2 )
            if (ising_v(isite, chan) /= 0) then
                kexpV_s(isite, mx, chan, iw) = ising_v(isite, chan)
            endif

            chan = 6 + ( (((j1-1)*2 + ixy)*2 + jxy) * 4 + 3 )
            if (ising_v(isite, chan) /= 0) then
                kexpV_s(isite, mx, chan, iw) = ising_v(isite, chan)
            endif

          end do
        end do
      enddo
    enddo
  endif

  ! if (abs(uxy) > 0.01) then
  !   do chan = 6, 37
  !     if ( MOD((chan-6)/8, 2) == 0 ) then
  !       do i = 1, lxy
  !         kexpV_s(i, mx, chan, iw) = ising_v(i, chan)
  !       end do
  !     else
  !       do i = lxy+1, nsites
  !         kexpV_s(i, mx, chan, iw) = ising_v(i, chan)
  !       end do
  !     end if
  !   end do
  ! end if

  if(g_ph*w_ph>0.001) then
    do i=1,nsites             ! Saved phonon fields
      keph_s(i,mx,iw)=ising_sum(i)
    end do
  end if
endsubroutine rec_fields

subroutine calgf(ibas, jbas, method, detpbas, length)
  use cpmc
  implicit none
  integer, intent(in) :: ibas, jbas, length
  real(sp) :: detpbas(length)

  integer, intent(in) :: method     ! 1: InitEnergy, 2: StpMeas
  integer :: i, info
  real(sp) :: det(2)

  ! ==== spinor（NSO×NE）====
  real(sp) :: pnew(NE,NE), pg(NE,NE)        ! pnew: inv(LR)，pg: inv或adj
  real(sp) :: work_ne(NE)
  integer  :: ipvt_ne(NE)

  ! ==== spinful（up: nsites×NUP；dn: nsites×NDN） ====
  real(sp) :: pnew_up(NUP,NUP),  pg_up(NUP,NUP)
  real(sp) :: pnew_dn(NDN,NDN),  pg_dn(NDN,NDN)
  real(sp) :: work_up(NUP),      work_dn(NDN)
  integer  :: ipvt_up(NUP),      ipvt_dn(NDN)

  external :: dgemm, dgefa, dgedi

  if (use_spinor) then
    ! =========================
    ! === Spinor===
    ! =========================
    detp(jbas) = 0.0_sp 

    if (method == 1) then
      associate(L => phiB(:,:,ibas), R => phiB(:,:,jbas))
        call dgemm('T','N', NE, NE, NSO, 1.0_8, L, NSO, R, NSO, 0.0_8, ovlpINV, NE)
      end associate
    else
      associate(L => phiZ(:,:,jbas), R => phi(:,:,ibas))
         call dgemm('N','N', NE, NE, NSO, 1.0_8, L, NSO, R, NSO, 0.0_8, ovlpINV, NE) ! Corrected 'T' to 'N' based on common usage for phiZ
      end associate
    end if

    pnew = ovlpINV
    call dgefa(pnew, NE, NE, ipvt_ne, info)
    if (info /= 0) then
      if (method == 1) then
          write(*,*) 'WARNING in calgf (spinor, method=1): dgefa failed with info =', info, ' for ibas=', ibas, ' jbas=', jbas
      else
          stop 'calgf(spinor): dgefa failed'
      endif

      pg = 0.0_sp
      gx = 0.0_sp
    else
      call dgedi(pnew, NE, NE, ipvt_ne, det, work_ne, 10_i4b)
      detp(jbas) = det(1) * TEN**det(2)

      if (method == 1) then
        pg = detp(jbas) * transpose(pnew)
        associate(L => phiB(:,:,ibas), R => phiB(:,:,jbas))
          call dgemm('N','N', NSO, NE, NE, 1.0_8, R, NSO, pg, NE, 0.0_8, tmp_ud, NSO)
          call dgemm('N','T', NSO, NSO, NE, -1.0_8, tmp_ud, NSO, L, NSO, 0.0_8, gx, NSO)
          do i = 1, NSO
             gx(i,i) = gx(i,i) + detp(jbas)
          end do
        end associate
      else ! method == 2
        pg = pnew
        associate(L => phiZ(:,:,jbas), R => phi(:,:,ibas))
           call dgemm('N','N', NSO, NE, NE, 1.0_8, R, NSO, pg, NE, 0.0_8, tmp_ud, NSO)
           ! Assuming phiZ (L) is NE x NSO, need 'T' and leading dim NSO
           call dgemm('N','T', NSO, NSO, NE, 1.0_8, tmp_ud, NSO, L, NSO, 0.0_8, gx, NSO)
        end associate
      endif
    endif

    detpbas(jbas) = cwibas(jbas) * detp(jbas)

  else
    ! =========================
    ! === Spinful (up/dn)  ===
    ! =========================
    detp_up(jbas) = 0.0_sp 
    detp_dn(jbas) = 0.0_sp 

    ! ---------- UP ----------
    if (method == 1) then
      associate(Lu => phiB_up(:,:,ibas), Ru => phiB_up(:,:,jbas))
        call dgemm('T','N', NUP, NUP, nsites, 1.0_8, Lu, nsites, Ru, nsites, 0.0_8, ovlpINV_up, NUP)
      end associate
    else
      associate(Lu => phiZ_up(:,:,jbas), Ru => phi_up(:,:,ibas))
         ! Assuming phiZ_up (Lu) is NUP x nsites, need 'N' and leading dim NUP
        call dgemm('N','N', NUP, NUP, nsites, 1.0_8, Lu, NUP, Ru, nsites, 0.0_8, ovlpINV_up, NUP)
      end associate
    end if

    pnew_up = ovlpINV_up
    call dgefa(pnew_up, NUP, NUP, ipvt_up, info)
    if (info /= 0) then
       if (method == 1) then
           write(*,*) 'WARNING in calgf (up, method=1): dgefa failed with info =', info, ' for ibas=', ibas, ' jbas=', jbas
      else
           stop 'calgf(up): dgefa failed'
      endif
      pg_up = 0.0_sp
      gx_up = 0.0_sp
    else
      call dgedi(pnew_up, NUP, NUP, ipvt_up, det, work_up, 10_i4b)
      detp_up(jbas) = det(1) * TEN**det(2)

      if (method == 1) then
        pg_up = detp_up(jbas) * transpose(pnew_up)
        associate(Lu => phiB_up(:,:,ibas), Ru => phiB_up(:,:,jbas))
          call dgemm('N','N', nsites, NUP, NUP, 1.0_8, Ru, nsites, pg_up, NUP, 0.0_8, tmp_up, nsites)
          call dgemm('N','T', nsites, nsites, NUP, -1.0_8, tmp_up, nsites, Lu, nsites, 0.0_8, gx_up, nsites)
          do i = 1, nsites
            gx_up(i,i) = gx_up(i,i) + detp_up(jbas)
          end do
        end associate
      else ! method == 2
        ! match method-1 convention: gx_up = detp_up * (I - R inv(ovlp) L)
        ! (calenergy expects the determinant-scaled I-G form; plain G here
        !  inverted the kinetic term and flipped the measured energy sign)
        pg_up = detp_up(jbas) * pnew_up
        associate(Lu => phiZ_up(:,:,jbas), Ru => phi_up(:,:,ibas))
          call dgemm('N','N', nsites, NUP, NUP, 1.0_8, Ru, nsites, pg_up, NUP, 0.0_8, tmp_up, nsites)
          call dgemm('N','N', nsites, nsites, NUP, -1.0_8, tmp_up, nsites, Lu, NUP, 0.0_8, gx_up, nsites)
          do i = 1, nsites
            gx_up(i,i) = gx_up(i,i) + detp_up(jbas)
          end do
        end associate
      endif
    endif

    ! ---------- DOWN ----------
    if (method == 1) then
      associate(Ld => phiB_dn(:,:,ibas), Rd => phiB_dn(:,:,jbas))
        call dgemm('T','N', NDN, NDN, nsites, 1.0_8, Ld, nsites, Rd, nsites, 0.0_8, ovlpINV_dn, NDN)
      end associate
    else
      associate(Ld => phiZ_dn(:,:,jbas), Rd => phi_dn(:,:,ibas))
        ! Assuming phiZ_dn (Ld) is NDN x nsites, need 'N' and leading dim NDN
         call dgemm('N','N', NDN, NDN, nsites, 1.0_8, Ld, NDN, Rd, nsites, 0.0_8, ovlpINV_dn, NDN)
      end associate
    end if

    pnew_dn = ovlpINV_dn
    call dgefa(pnew_dn, NDN, NDN, ipvt_dn, info)
    if (info /= 0) then
      if (method == 1) then
          write(*,*) 'WARNING in calgf (dn, method=1): dgefa failed with info =', info, ' for ibas=', ibas, ' jbas=', jbas
      else
          stop 'calgf(dn): dgefa failed'
      endif
      pg_dn = 0.0_sp
      gx_dn = 0.0_sp
    else
      call dgedi(pnew_dn, NDN, NDN, ipvt_dn, det, work_dn, 10_i4b)
      detp_dn(jbas) = det(1) * TEN**det(2)

      if (method == 1) then
        pg_dn = detp_dn(jbas) * transpose(pnew_dn)
        associate(Ld => phiB_dn(:,:,ibas), Rd => phiB_dn(:,:,jbas))
          call dgemm('N','N', nsites, NDN, NDN, 1.0_8, Rd, nsites, pg_dn, NDN, 0.0_8, tmp_dn, nsites)
          call dgemm('N','T', nsites, nsites, NDN, -1.0_8, tmp_dn, nsites, Ld, nsites, 0.0_8, gx_dn, nsites)
          do i = 1, nsites
            gx_dn(i,i) = gx_dn(i,i) + detp_dn(jbas)
          end do
        end associate
      else ! method == 2
        ! match method-1 convention: gx_dn = detp_dn * (I - R inv(ovlp) L)
        pg_dn = detp_dn(jbas) * pnew_dn
        associate(Ld => phiZ_dn(:,:,jbas), Rd => phi_dn(:,:,ibas))
          call dgemm('N','N', nsites, NDN, NDN, 1.0_8, Rd, nsites, pg_dn, NDN, 0.0_8, tmp_dn, nsites)
          call dgemm('N','N', nsites, nsites, NDN, -1.0_8, tmp_dn, nsites, Ld, NDN, 0.0_8, gx_dn, nsites)
          do i = 1, nsites
            gx_dn(i,i) = gx_dn(i,i) + detp_dn(jbas)
          end do
        end associate
      endif
    endif

    detpbas(jbas) = cwibas(jbas) * detp_up(jbas) * detp_dn(jbas)

  endif
end subroutine calgf
! subroutine calgf(ibas, jbas, method, detpbas)
!   use cpmc
!   implicit none
!   integer, intent(in) :: ibas, jbas
!   real(sp), intent(inout) :: detpbas(:)

!   integer, intent(in) :: method     ! 1: InitEnergy, 2: StpMeas
!   integer :: i, info
!   real(sp) :: det(2)

!   ! ==== spinor（NSO×NE）====
!   real(sp) :: pnew(NE,NE), pg(NE,NE)        ! pnew: inv(LR)，pg: inv或adj
!   real(sp) :: work_ne(NE)
!   integer  :: ipvt_ne(NE)

!   ! ==== spinful（up: nsites×NUP；dn: nsites×NDN） ====
!   real(sp) :: pnew_up(NUP,NUP),  pg_up(NUP,NUP)
!   real(sp) :: pnew_dn(NDN,NDN),  pg_dn(NDN,NDN)
!   real(sp) :: work_up(NUP),      work_dn(NDN)
!   integer  :: ipvt_up(NUP),      ipvt_dn(NDN)

!   external :: dgemm, dgefa, dgedi

!   if (use_spinor) then
!     ! =========================
!     ! === Spinor===
!     ! =========================
!     if (method == 1) then
!       ! L = phiB(:,:,ibas), R = phiB(:,:,jbas)
!       associate(L => phiB(:,:,ibas), R => phiB(:,:,jbas))
!         ! LR = L^T * R → 存入 ovlpINV
!         call dgemm('T','N', NE, NE, NSO, 1.0_8, L, NSO, R, NSO, 0.0_8, ovlpINV, NE)
!       end associate
!     else
!       ! L = phiZ(:,:,ibas), R = phi(:,:,jbas)
!       associate(L => phiZ(:,:,ibas), R => phi(:,:,jbas))
!         call dgemm('T','N', NE, NE, NSO, 1.0_8, L, NSO, R, NSO, 0.0_8, ovlpINV, NE)
!       end associate
!     end if

!     ! 求 inv(LR) 与 det(LR)
!     pnew = ovlpINV
!     call dgefa(pnew, NE, NE, ipvt_ne, info)
!     if (info /= 0) then
!       detp(jbas) = 0.0_sp
!       if (method == 2) stop 'calgf(spinor): dgefa failed'
!     else
!       call dgedi(pnew, NE, NE, ipvt_ne, det, work_ne, 10_i4b)
!       detp(jbas) = det(1) * TEN**det(2)
!     end if

!     if (method == 1) then
!       ! pg = adj(LR) = det(LR) * transpose(inv(LR))
!       pg = detp(jbas) * transpose(pnew)
!       ! gx = - R * pg * L^T + det(LR) * I
!       associate(L => merge(phiB(:,:,ibas), phiZ(:,:,ibas), .false.), &
!                R => merge(phiB(:,:,jbas), phi(:,:,jbas), .false.))
!         call dgemm('N','N', NSO, NE, NE, 1.0_8, R, NSO, pg, NE, 0.0_8, tmp_ud, NSO)
!         call dgemm('N','T', NSO, NSO, NE, -1.0_8, tmp_ud, NSO, L, NSO, 0.0_8, gx, NSO)
!         do i = 1, NSO
!           gx(i,i) = gx(i,i) + detp(jbas)
!         end do
!       end associate
!     else
!       ! pg = inv(LR)
!       pg = pnew
!       ! gx =  R * pg * L^T
!       associate(L => phiZ(:,:,ibas), R => phi(:,:,jbas))
!         call dgemm('N','N', NSO, NE, NE, 1.0_8, R, NSO, pg, NE, 0.0_8, tmp_ud, NSO)
!         call dgemm('N','T', NSO, NSO, NE, 1.0_8, tmp_ud, NSO, L, NSO, 0.0_8, gx, NSO)
!       end associate
!     end if
!     detpbas(jbas) = cwibas(jbas) * detp(jbas)

!   else
!     ! =========================
!     ! === Spinful (up/dn)  ===
!     ! =========================

!     ! ---------- UP ----------
!     if (method == 1) then
!       ! L_up = phiB_up(:,:,ibas), R_up = phiB_up(:,:,jbas)
!       associate(Lu => phiB_up(:,:,ibas), Ru => phiB_up(:,:,jbas))
!         call dgemm('T','N', NUP, NUP, nsites, 1.0_8, Lu, nsites, Ru, nsites, 0.0_8, ovlpINV_up, NUP)
!       end associate
!     else
!       ! L_up = phiZ_up(:,:,ibas), R_up = phi_up(:,:,jbas)
!       associate(Lu => phiZ_up(:,:,ibas), Ru => phi_up(:,:,jbas))
!         call dgemm('T','N', NUP, NUP, nsites, 1.0_8, Lu, nsites, Ru, nsites, 0.0_8, ovlpINV_up, NUP)
!       end associate
!     end if

!     pnew_up = ovlpINV_up
!     call dgefa(pnew_up, NUP, NUP, ipvt_up, info)
!     if (info /= 0) then
!       detp_up(jbas) = 0.0_sp
!       if (method == 2) stop 'calgf(up): dgefa failed'
!     else
!       call dgedi(pnew_up, NUP, NUP, ipvt_up, det, work_up, 10_i4b)
!       detp_up(jbas) = det(1) * TEN**det(2)
!     end if

!     if (method == 1) then
!       ! pg_up = adj(LR_up)
!       pg_up = detp_up(jbas) * transpose(pnew_up)
!       associate(Lu => phiB_up(:,:,ibas), Ru => phiB_up(:,:,jbas))
!         call dgemm('N','N', nsites, NUP, NUP, 1.0_8, Ru, nsites, pg_up, NUP, 0.0_8, tmp_up, nsites)
!         call dgemm('N','T', nsites, nsites, NUP, -1.0_8, tmp_up, nsites, Lu, nsites, 0.0_8, gx_up, nsites)
!         do i = 1, nsites
!           gx_up(i,i) = gx_up(i,i) + detp_up(jbas)
!         end do
!       end associate
!     else
!       ! pg_up = inv(LR_up)
!       pg_up = pnew_up
!       associate(Lu => phiZ_up(:,:,ibas), Ru => phi_up(:,:,jbas))
!         call dgemm('N','N', nsites, NUP, NUP, 1.0_8, Ru, nsites, pg_up, NUP, 0.0_8, tmp_up, nsites)
!         call dgemm('N','T', nsites, nsites, NUP, 1.0_8, tmp_up, nsites, Lu, nsites, 0.0_8, gx_up, nsites)
!       end associate
!     end if

!     ! ---------- DOWN ----------
!     if (method == 1) then
!       ! L_dn = phiB_dn(:,:,ibas), R_dn = phiB_dn(:,:,jbas)
!       associate(Ld => phiB_dn(:,:,ibas), Rd => phiB_dn(:,:,jbas))
!         call dgemm('T','N', NDN, NDN, nsites, 1.0_8, Ld, nsites, Rd, nsites, 0.0_8, ovlpINV_dn, NDN)
!       end associate
!     else
!       ! L_dn = phiZ_dn(:,:,ibas), R_dn = phi_dn(:,:,jbas)
!       associate(Ld => phiZ_dn(:,:,ibas), Rd => phi_dn(:,:,jbas))
!         call dgemm('T','N', NDN, NDN, nsites, 1.0_8, Ld, nsites, Rd, nsites, 0.0_8, ovlpINV_dn, NDN)
!       end associate
!     end if

!     pnew_dn = ovlpINV_dn
!     call dgefa(pnew_dn, NDN, NDN, ipvt_dn, info)
!     if (info /= 0) then
!       detp_dn(jbas) = 0.0_sp
!       if (method == 2) stop 'calgf(dn): dgefa failed'
!     else
!       call dgedi(pnew_dn, NDN, NDN, ipvt_dn, det, work_dn, 10_i4b)
!       detp_dn(jbas) = det(1) * TEN**det(2)
!     end if

!     if (method == 1) then
!       ! pg_dn = adj(LR_dn)
!       pg_dn = detp_dn(jbas) * transpose(pnew_dn)
!       associate(Ld => phiB_dn(:,:,ibas), Rd => phiB_dn(:,:,jbas))
!         call dgemm('N','N', nsites, NDN, NDN, 1.0_8, Rd, nsites, pg_dn, NDN, 0.0_8, tmp_dn, nsites)
!         call dgemm('N','T', nsites, nsites, NDN, -1.0_8, tmp_dn, nsites, Ld, nsites, 0.0_8, gx_dn, nsites)
!         do i = 1, nsites
!           gx_dn(i,i) = gx_dn(i,i) + detp_dn(jbas)
!         end do
!       end associate
!     else
!       ! pg_dn = inv(LR_dn)
!       pg_dn = pnew_dn
!       associate(Ld => phiZ_dn(:,:,ibas), Rd => phi_dn(:,:,jbas))
!         call dgemm('N','N', nsites, NDN, NDN, 1.0_8, Rd, nsites, pg_dn, NDN, 0.0_8, tmp_dn, nsites)
!         call dgemm('N','T', nsites, nsites, NDN, 1.0_8, tmp_dn, nsites, Ld, nsites, 0.0_8, gx_dn, nsites)
!       end associate
!     end if

!     detpbas(jbas) = cwibas(jbas) * detp_up(jbas) * detp_dn(jbas)

!   end if
! end subroutine calgf

subroutine backphi(iw, ibas, phip_up, phip_dn, phip)
  use cpmc
  use jiekou, only: stblzbk, stblzbk_spinor
  use, intrinsic :: ieee_arithmetic

  implicit none
  !-------------------- arguments --------------------!
  integer,  intent(in)  :: iw, ibas
  real(sp), intent(out) :: phip_up(NUP, nsites)
  real(sp), intent(out) :: phip_dn(NDN, nsites)
  real(sp), intent(out) :: phip(NE, NSO)

  !-------------------- locals -----------------------!
  integer  :: i, j, mx, ip, lsave, chan, count
  integer  :: n1, n2, n3, imax
  integer  :: j1, jnb, ixy, jxy, isite, jsite, ic
  ! spinor work (NE x NSO)
  real(sp) :: ud_tmp(NE, NSO, 2)
  ! spinful work
  real(sp) :: up_tmp(NUP, nsites, 2), dn_tmp(NDN, nsites, 2)

  ! external :: dgemm

  !===================== SPINOR =====================!
  if (use_spinor) then
    count = 0
    ud_tmp(:,:,1) = phiZ(:,:,ibas) 

    n1 = 1; n2 = 2
    do mx = itvl_m, 1, -1
      count = count + 1

      if (mod(count, itvlorth) == 1) then
        call dgemm('N','N', NE, NSO, NSO, 1.0_8, ud_tmp(:,:,n1), NE, expT(:,:),  NSO, 0.0_8, ud_tmp(:,:,n2), NE)
      else
        call dgemm('N','N', NE, NSO, NSO, 1.0_8, ud_tmp(:,:,n1), NE, exp2T(:,:), NSO, 0.0_8, ud_tmp(:,:,n2), NE)
      end if

      do chan = 1, channels
        imax = nsites
        if (chan >= 2 .and. chan <= 5) imax = lxy

        do ip = 1, nlsi(chan)
          do i = 1, imax
            lsave = kexpV_s(i, mx, chan, iw)
            if (lsave /= 0) then
              if (spinlsi(chan) == 1) then
                ud_tmp(ip, i, n2) = ud_tmp(ip, i, n2) * expV(i, lsave, 1, chan)
              else
                ud_tmp(ip, i, n2) = ud_tmp(ip, i, n2) * expV(i, lsave, 2, chan)
              end if
            end if
          end do
        end do

        do ip = 1, nlsj(chan)
          do j = 1, imax
            lsave = kexpV_s(j, mx, chan, iw)
            if (lsave /= 0) then
              if (spinlsj(chan) == 1) then
                ud_tmp(ip, j, n2) = ud_tmp(ip, j, n2) * expV(j, lsave, 1, chan)
              else
                ud_tmp(ip, j, n2) = ud_tmp(ip, j, n2) * expV(j, lsave, 2, chan)
              end if
            end if
          end do
        end do
      end do

      n3 = n2; n2 = n1; n1 = n3

      if (mod(count, itvlorth) == 0) then
        call dgemm('N','N', NE, NSO, NSO, 1.0_8, ud_tmp(:,:,n1), NE, expT(:,:), NSO, 0.0_8, ud_tmp(:,:,n2), NE) ! Exp(-K/2)
        n3 = n2; n2 = n1; n1 = n3
        call stblzbk_spinor(iw, ud_tmp(:,:,n1))
      end if
    end do

    phip = ud_tmp(:,:,n1)

  !===================== SPINFUL ====================!
  else
    count = 0
    up_tmp(:,:,1) =  phiZ_up(:,:,ibas) 
    dn_tmp(:,:,1) =  phiZ_dn(:,:,ibas) 

    n1 = 1; n2 = 2
    do mx = itvl_m, 1, -1
      count = count + 1

      if (mod(count, itvlorth) == 1) then
        call dgemm('N','N', NUP, nsites, nsites, 1.0_8, up_tmp(:,:,n1), NUP, expT_ud(:,:,1),  nsites, 0.0_8, up_tmp(:,:,n2), NUP)
        call dgemm('N','N', NDN, nsites, nsites, 1.0_8, dn_tmp(:,:,n1), NDN, expT_ud(:,:,2),  nsites, 0.0_8, dn_tmp(:,:,n2), NDN)
      else
        call dgemm('N','N', NUP, nsites, nsites, 1.0_8, up_tmp(:,:,n1), NUP, exp2T_ud(:,:,1), nsites, 0.0_8, up_tmp(:,:,n2), NUP)
        call dgemm('N','N', NDN, nsites, nsites, 1.0_8, dn_tmp(:,:,n1), NDN, exp2T_ud(:,:,2), nsites, 0.0_8, dn_tmp(:,:,n2), NDN)
      end if

      ! ---- replay recorded HS fields for slice mx (diagonal ops commute) ----
      ! Each interaction term V*n_{a,si}*n_{b,sj} scaled the FIRST site a by
      ! expV(a,field,1,chan) on spin si, and the PARTNER site b by
      ! expV(a,field,2,chan) on spin sj (field stored at a). backphi must mirror
      ! that bond structure -- the old loop applied the 2nd operator to site j
      ! (== a only for the same-site uxx channel), which is why uxx worked but
      ! uxy/v did not.
      ! -- uxx (chan 1): a = b = i --
      do i = 1, nsites
        lsave = kexpV_s(i, mx, 1, iw)
        if (lsave /= 0) then
          up_tmp(:, i, n2) = up_tmp(:, i, n2) * expV(i, lsave, 1, 1)   ! spinlsi(1)=up
          dn_tmp(:, i, n2) = dn_tmp(:, i, n2) * expV(i, lsave, 2, 1)   ! spinlsj(1)=dn
        end if
      end do
      ! -- uxy (chan 2-5): a = i (orbital 1), b = i+lxy (orbital 2) --
      do chan = 2, 5
        do i = 1, lxy
          lsave = kexpV_s(i, mx, chan, iw)
          if (lsave == 0) cycle
          if (spinlsi(chan) == 1) then
            up_tmp(:, i, n2) = up_tmp(:, i, n2) * expV(i, lsave, 1, chan)
          else
            dn_tmp(:, i, n2) = dn_tmp(:, i, n2) * expV(i, lsave, 1, chan)
          end if
          if (spinlsj(chan) == 1) then
            up_tmp(:, i+lxy, n2) = up_tmp(:, i+lxy, n2) * expV(i, lsave, 2, chan)
          else
            dn_tmp(:, i+lxy, n2) = dn_tmp(:, i+lxy, n2) * expV(i, lsave, 2, chan)
          end if
        end do
      end do
      ! -- v (chan 6-37): a = i+ixy*lxy, b = idis(i,j1)+jxy*lxy (mirror Vznn/rec_fields) --
      if (abs(v) > 0.01) then
        do i = 1, lxy
          do j1 = 1, 2
            jnb = idis(i, j1)
            do ixy = 0, 1
              do jxy = 0, 1
                isite = i   + ixy*lxy
                jsite = jnb + jxy*lxy
                chan  = 6 + ((((j1-1)*2 + ixy)*2 + jxy) * 4)
                do ic = 0, 3
                  lsave = kexpV_s(isite, mx, chan+ic, iw)
                  if (lsave == 0) cycle
                  if (spinlsi(chan+ic) == 1) then
                    up_tmp(:, isite, n2) = up_tmp(:, isite, n2) * expV(isite, lsave, 1, chan+ic)
                  else
                    dn_tmp(:, isite, n2) = dn_tmp(:, isite, n2) * expV(isite, lsave, 1, chan+ic)
                  end if
                  if (spinlsj(chan+ic) == 1) then
                    up_tmp(:, jsite, n2) = up_tmp(:, jsite, n2) * expV(isite, lsave, 2, chan+ic)
                  else
                    dn_tmp(:, jsite, n2) = dn_tmp(:, jsite, n2) * expV(isite, lsave, 2, chan+ic)
                  end if
                end do
              end do
            end do
          end do
        end do
      end if

      n3 = n2; n2 = n1; n1 = n3

      if (mod(count, itvlorth) == 0) then
        call dgemm('N','N', NUP, nsites, nsites, 1.0_8, up_tmp(:,:,n1), NUP, expT_ud(:,:,1), nsites, 0.0_8, up_tmp(:,:,n2), NUP)  ! Exp(-K/2)
        call dgemm('N','N', NDN, nsites, nsites, 1.0_8, dn_tmp(:,:,n1), NDN, expT_ud(:,:,2), nsites, 0.0_8, dn_tmp(:,:,n2), NDN)
        n3 = n2; n2 = n1; n1 = n3
        call stblzbk(iw, up_tmp(:,:,n1), dn_tmp(:,:,n1))
      end if
    end do

    phip_up = up_tmp(:,:,n1)
    phip_dn = dn_tmp(:,:,n1)
  end if
end subroutine backphi

! subroutine backphi(iw,ibas,phip_up,phip_dn,phip)
!   !iw is the walker indexuse cpmc
!   use cpmc
!   use jiekou,only:stblzbk
!   integer::i,j,k,ibas,m,count,chan
!   integer::n1,n2,n3
!   integer::iw,mx,ip,lsave
!   real(sp)::phip_up(NUP,nsites),phip_dn(NDN,nsites)
!   real(sp)::up_tmp(NUP,nsites,2),dn_tmp(NDN,nsites,2)
!   real(sp)::ud_tmp(NSO,NE,2)
!   external :: dgemm,dgedi,dgefa
!   integer :: i_lo,i_hi,j_lo,j_hi,t,ixy,jxy
  
!   !----------------------------------------------------------------!
!   ! < Phi(k',ibas) | = < Phi_T(ibas) | B(N+m,iw)B(N+m-1,iw) ... B(N+1,iw)!
!   !               = < Phi_T(ibas) | Exp(-K/2) V(N+m,iw)
!   !                              Exp(-K)   V(N+m-1,iw) ...
!   !                              Exp(-K)   V(N-1,iw)   Exp(-K/2)!
!   !-------------------------------------------------------------!
!   !-----------------------------------------------------------------!
!   ! For spin up
!   ! < Phi_T(ibas) | Exp(-K/2); (Np x nsites) * (nsites x nsites) ...!
!   !-----------------------------------------------------------------!
!   if (use_spinor) then    
!     count=0
!     ud_tmp(:,:,1)=phiZ(:,:,ibas)

!     n1=1
!     n2=2
!     do mx=itvl_m,1,-1
!       count=count+1
!       if(mod(count,itvlorth)==1) then
!         call dgemm('N', 'N', NE, NSO, NSO, 1.0_8, ud_tmp(:,:,n1), NE, expT(:,:,1), NSO, 0.0_8, ud_tmp(:,:,n2), NE) !cao
!       else
!         call dgemm('N', 'N', NE, NSO, NSO, 1.0_8, up_tmp(:,:,n1), NE, exp2T(:,:,1), NSO, 0.0_8, up_tmp(:,:,n2), NE) !cao
!       end if
!       !----------------------------------------------------------!
!       ! V(N+mx,iw); (Np x nsites) * (nsites x nsites) Diagonal!
!       !----------------------------------------------------------! 
!       do chan=1,channels
!         imax = nsites
!         if (chan >= 2 .and. chan <= 5) imax = lxy

!         ! every channcel different ni,nj, spini,spinj
!         do ip=1,nlsi(chan)
!           do i=1,imax
!             lsave=kexpV_s(i,mx,iw,chan)
!             if (spinlsi(chan)==1) then
!               ! write(*,*) 'lsave',lsave,i, chan
!               ud_tmp(ip,i,n2)=ud_tmp(ip,i,n2)*expV(i,lsave,1,chan)
!             else
!               ud_tmp(ip,i,n2)=ud_tmp(ip,i,n2)*expV(i,lsave,2,chan)
!             endif
!           end do
!         end do

!         do ip=1,nlsj(chan)
!           do j=1,imax
!           lsave=kexpV_s(j,mx,iw,chan)
!             if (spinlsj(chan)==1) then
!               ud_tmp(ip,j,n2)=ud_tmp(ip,j,n2)*expV(j,lsave,1,chan)
!             else
!               ud_tmp(ip,j,n2)=ud_tmp(ip,j,n2)*expV(j,lsave,2,chan)
!             endif
!           end do
!         enddo
!       enddo

!       !--------------Zhongbing-phonon fields-----------!
!       ! if(g_ph*w_ph>0.001) then
!       !   do 101 ip=1,NUP
!       !   do 101 j=1,nsites

!       !     lsave=keph_s(j,mx,iw)
!       !     up_tmp(ip,j,n2)=up_tmp(ip,j,n2)*eph(j,lsave,1)

!       !   101        continue

!       !         do 202 ip=1,NDN
!       !         do 202 j=1,nsites

!       !             lsave=keph_s(j,mx,iw)
!       !     dn_tmp(ip,j,n2)=dn_tmp(ip,j,n2)*eph(j,lsave,2)

!       !   202        continue
!       ! end if

!       n3=n2
!       n2=n1
!       n1=n3
!       if(mod(count,itvlorth)==0) then
!         call dgemm('N', 'N', NE, NSO, NSO, 1.0_8, ud_tmp(:,:,n1), NE, expT(:,:), NSO, 0.0_8, ud_tmp(:,:,n2), NE) !cao
!         !! Exp(-K/2)!!
!         n3=n2
!         n2=n1
!         n1=n3
!         call stblzbk_spinor(iw,ud_tmp(:,:,n1))
!       end if
!     end do
!     phip=ud_tmp(:,:,n1)

!   else 
!     count=0
!     up_tmp(:,:,1)=phiZ_up(:,:,ibas)
!     dn_tmp(:,:,1)=phiZ_dn(:,:,ibas)
!     n1=1
!     n2=2
!     do mx=itvl_m,1,-1
!       count=count+1
!       if(mod(count,itvlorth)==1) then
!         call dgemm('N', 'N', NUP, nsites, nsites, 1.0_8, up_tmp(:,:,n1), NUP, expT_ud(:,:,1), nsites, 0.0_8, up_tmp(:,:,n2), NUP) !cao
!         call dgemm('N', 'N', NDN, nsites, nsites, 1.0_8, dn_tmp(:,:,n1), NDN, expT_ud(:,:,2), nsites, 0.0_8, dn_tmp(:,:,n2), NDN) !cao
!         !up_tmp(:,:,n2)=matmul(up_tmp(:,:,n1),expT_ud(:,:,1))
!         !dn_tmp(:,:,n2)=matmul(dn_tmp(:,:,n1),expT_ud(:,:,2))
!       else
!         call dgemm('N', 'N', NUP, nsites, nsites, 1.0_8, up_tmp(:,:,n1), NUP, exp2T_ud(:,:,1), nsites, 0.0_8, up_tmp(:,:,n2), NUP) !cao
!         call dgemm('N', 'N', NDN, nsites, nsites, 1.0_8, dn_tmp(:,:,n1), NDN, exp2T_ud(:,:,2), nsites, 0.0_8, dn_tmp(:,:,n2), NDN) !cao
!         !up_tmp(:,:,n2)=matmul(up_tmp(:,:,n1),exp2T_ud(:,:,1))
!         !dn_tmp(:,:,n2)=matmul(dn_tmp(:,:,n1),exp2T_ud(:,:,2))
!       end if
!       !----------------------------------------------------------!
!       ! V(N+mx,iw); (Np x nsites) * (nsites x nsites) Diagonal!
!       !----------------------------------------------------------! 
!       do chan=1,channels
!         imax = nsites
!         if (chan >= 2 .and. chan <= 5) imax = lxy
!         ! every channcel different ni,nj, spini,spinj
!         do ip=1,nlsi(chan)
!           do i=1,imax
!             if (spinlsi(chan)==1) then
!               lsave=kexpV_s(i,mx,iw,chan)
!               ! write(*,*) 'lsave',lsave,i, chan
!               up_tmp(ip,i,n2)=up_tmp(ip,i,n2)*expV(i,lsave,1,chan)
!             else
!               lsave=kexpV_s(i,mx,iw,chan)
!               dn_tmp(ip,i,n2)=dn_tmp(ip,i,n2)*expV(i,lsave,2,chan)
!             endif
!           end do
!         end do

!         do ip=1,nlsj(chan)
!           do j=1,imax
!             if (spinlsj(chan)==1) then
!               lsave=kexpV_s(j,mx,iw,chan)
!               up_tmp(ip,j,n2)=up_tmp(ip,j,n2)*expV(j,lsave,1,chan)
!             else
!               lsave=kexpV_s(j,mx,iw,chan)
!               dn_tmp(ip,j,n2)=dn_tmp(ip,j,n2)*expV(j,lsave,2,chan)
!             endif
!           end do
!         enddo
!       enddo

!       !--------------Zhongbing-phonon fields-----------!
!       ! if(g_ph*w_ph>0.001) then
!       !   do 101 ip=1,NUP
!       !   do 101 j=1,nsites

!       !     lsave=keph_s(j,mx,iw)
!       !     up_tmp(ip,j,n2)=up_tmp(ip,j,n2)*eph(j,lsave,1)

!       !   101        continue

!       !         do 202 ip=1,NDN
!       !         do 202 j=1,nsites

!       !             lsave=keph_s(j,mx,iw)
!       !     dn_tmp(ip,j,n2)=dn_tmp(ip,j,n2)*eph(j,lsave,2)

!       !   202        continue
!       ! end if

!       n3=n2
!       n2=n1
!       n1=n3
!       if(mod(count,itvlorth)==0) then
!         call dgemm('N', 'N', NUP, nsites, nsites, 1.0_8, up_tmp(:,:,n1), NUP, expT_ud(:,:,1), nsites, 0.0_8, up_tmp(:,:,n2), NUP) !cao
!         call dgemm('N', 'N', NDN, nsites, nsites, 1.0_8, dn_tmp(:,:,n1), NDN, expT_ud(:,:,2), nsites, 0.0_8, dn_tmp(:,:,n2), NDN) !cao
!         !up_tmp(:,:,n2)=matmul(up_tmp(:,:,n1),expT_ud(:,:,1))
!         !dn_tmp(:,:,n2)=matmul(dn_tmp(:,:,n1),expT_ud(:,:,2))   !! Exp(-K/2)!!
!         n3=n2
!         n2=n1
!         n1=n3
!         call stblzbk(iw,up_tmp(:,:,n1),dn_tmp(:,:,n1))
!       end if
!     end do
!     phip_up=up_tmp(:,:,n1)
!     phip_dn=dn_tmp(:,:,n1)

!   end if
!   ! return
! end subroutine backphi