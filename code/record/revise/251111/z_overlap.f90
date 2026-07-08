! SUBROUTINE MkInitOvlps
!   use cpmc
!   use jiekou 
!   implicit none

!   integer :: i, j, k, l, m, n
!   integer :: ibas, jbas, info
!   real(sp) :: s(NWFBAS)
!   real(sp) :: det(2)

!   integer :: ipvt_e(NE),   ipvt_up(NUP),   ipvt_dn(NDN)
!   real(sp) :: work_e(NE),  work_up(NUP),   work_dn(NDN)

!   real(sp) :: temp(NE,NE)
!   real(sp) :: det_ud(NWFBAS)

!   real(sp) :: temp_up(NUP,NUP), temp_dn(NDN,NDN)
!   real(sp) :: det_up(NWFBAS), det_dn(NWFBAS)

!   external :: dgemm

!   do jbas = 1, NIW
!      s = 0.0

!      do ibas = 1, NWFBAS
!         if (use_spinor) then
!            ! temp = phiZ(:,:,ibas) * phiB(:,:,jbas)  ->  (NE x NSO) * (NSO x NE) = (NE x NE)
!            call dgemm('N','N', NE, NE, NSO, 1.0, phiZ(:,:,ibas), NE, &
!                       phiB(:,:,jbas), NSO, 0.0, temp, NE)

!            call dgefa(temp, NE, NE, ipvt_e, info)
!            if (info /= 0) stop 'Problem in dgefa (spinor)'

!            call sgedi(temp, NE, NE, ipvt_e, det, work_e, 10_i4b)
!            det_ud(ibas) = det(1) * TEN**det(2)

!            s(ibas) = cwfbas(ibas) * det_ud(ibas)

!         else
!            ! ----- up -----
!            ! temp_up = phiZ_up(:,:,ibas) * phiB_up(:,:,jbas)  -> (NUP x NSTATES) * (NSTATES x NUP) = (NUP x NUP)
!            call dgemm('N','N', NUP, NUP, NSTATES, 1.0, phiZ_up(:,:,ibas), NUP, &
!                       phiB_up(:,:,jbas), NSTATES, 0.0, temp_up, NUP)

!            call dgefa(temp_up, NUP, NUP, ipvt_up, info)
!            if (info /= 0) stop 'Problem in dgefa (up)'

!            call sgedi(temp_up, NUP, NUP, ipvt_up, det, work_up, 10_i4b)
!            det_up(ibas) = det(1) * TEN**det(2)

!            ! ----- down -----
!            ! temp_dn = phiZ_dn(:,:,ibas) * phiB_dn(:,:,jbas)  -> (NDN x NSTATES) * (NSTATES x NDN) = (NDN x NDN)
!            call dgemm('N','N', NDN, NDN, NSTATES, 1.0, phiZ_dn(:,:,ibas), NDN, &
!                       phiB_dn(:,:,jbas), NSTATES, 0.0, temp_dn, NDN)

!            call dgefa(temp_dn, NDN, NDN, ipvt_dn, info)
!            if (info /= 0) stop 'Problem in dgefa (dn)'

!            call sgedi(temp_dn, NDN, NDN, ipvt_dn, det, work_dn, 10_i4b)
!            det_dn(ibas) = det(1) * TEN**det(2)

!            s(ibas) = cwfbas(ibas) * det_up(ibas) * det_dn(ibas)
!         end if
!      end do

!      ovlpINIT(jbas) = sum(s)
!   end do

!   do ibas = 1, NIW
!      sgnINIT(ibas)  = sign(ONE, ovlpINIT(ibas))
!      ovlpINIT(ibas) = abs(ovlpINIT(ibas))
!   end do
! END SUBROUTINE MkInitOvlps

SUBROUTINE Ovlps(iw)
  use cpmc
!$ use omp_lib
  implicit none
  integer, intent(in) :: iw
  integer :: ibas, n, info, lwork
  integer, allocatable :: ipiv(:)
  real(sp), allocatable :: A(:,:), work(:), detlist(:)
  real(sp) :: ovlpNEW, detval
  real(sp), dimension(NWFBAS) :: det_up, det_dn

  ovlpNEW = ZERO

  if (use_spinor) then
    n = NE
    lwork = max(4*n, 64)
    allocate(ipiv(n), A(n,n), work(lwork), detlist(NWFBAS))

!$omp parallel do default(shared) private(ibas, A, ipiv, work, info, detval) reduction(+:ovlpNEW) if(n>=64)
    do ibas = 1, NWFBAS
      call dgemm('N','N', n, n, NSO, 1.0_8, phiZ(1,1,ibas), n, phi(1,1,iw), NSO, 0.0_8, A, n)
      call dgetrf(n, n, A, n, ipiv, info)
      if (info /= 0) stop 'Ovlps(spinor): dgetrf failed'
      detval = det_from_lu(A, ipiv, n)
      call dgetri(n, A, n, ipiv, work, lwork, info)
      if (info /= 0) stop 'Ovlps(spinor): dgetri failed'
!$omp critical
      g(:,:,ibas) = A
!$omp end critical
      ovlpNEW = ovlpNEW + cwfbas(ibas) * detval
    end do
!$omp end parallel do

    detbas(:) = cwfbas(:) * detlist(:)  ! optional if needed elsewhere

    deallocate(ipiv, A, work, detlist)

  else
    ! Legacy per-spin path
    ! up
    n = NUP
    lwork = max(4*n, 64)
    allocate(ipiv(n), A(n,n), work(lwork), detlist(NWFBAS))
!$omp parallel do default(shared) private(ibas, A, ipiv, work, info, detval) if(n>=64)
    do ibas = 1, NWFBAS
      call dgemm('N','N', n, n, nsites, 1.0_8, phiZ_up(1,1,ibas), n, phi_up(1,1,iw), nsites, 0.0_8, A, n)
      call dgetrf(n, n, A, n, ipiv, info)
      if (info /= 0) stop 'Ovlps(up): dgetrf failed'
      detval = det_from_lu(A, ipiv, n)
      call dgetri(n, A, n, ipiv, work, lwork, info)
      if (info /= 0) stop 'Ovlps(up): dgetri failed'
!$omp critical
      g_up(:,:,ibas) = A
      det_up(ibas) = detval
!$omp end critical
    end do
!$omp end parallel do
    deallocate(ipiv, A, work)

    ! down
    n = NDN
    lwork = max(4*n, 64)
    allocate(ipiv(n), A(n,n), work(lwork))
!$omp parallel do default(shared) private(ibas, A, ipiv, work, info, detval) if(n>=64)
    do ibas = 1, NWFBAS
      call dgemm('N','N', n, n, nsites, 1.0_8, phiZ_dn(1,1,ibas), n, phi_dn(1,1,iw), nsites, 0.0_8, A, n)
      call dgetrf(n, n, A, n, ipiv, info)
      if (info /= 0) stop 'Ovlps(dn): dgetrf failed'
      detval = det_from_lu(A, ipiv, n)
      call dgetri(n, A, n, ipiv, work, lwork, info)
      if (info /= 0) stop 'Ovlps(dn): dgetri failed'
!$omp critical
      g_dn(:,:,ibas) = A
      det_dn(ibas) = detval
      detbas(ibas)=cwfbas(ibas)*det_up(ibas)*det_dn(ibas)
!$omp end critical
    end do
!$omp end parallel do

    ovlpNEW = sum( cwfbas(:) * det_up(:) * det_dn(:) )

    deallocate(ipiv, A, work)
  end if

  if (ovlpNEW * sgn(iw) <= ZERO) then
    wgtwlkr(iw) = ZERO
  else
    wgtwlkr(iw) = (ovlpNEW/ovlpDET(iw)) * wgtwlkr(iw)
    ovlpDET(iw) = ovlpNEW
  end if

  return

contains
  pure real(sp) function det_from_lu(U, ipiv, n) result(detv)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: ipiv(n)
    real(sp), intent(in) :: U(n,n)
    integer :: i, sign_piv
    real(sp) :: prod
    prod = 1.0
    do i=1,n
      prod = prod * U(i,i)
    end do
    sign_piv = 1
    do i=1,n
      if (ipiv(i) /= i) sign_piv = -sign_piv
    end do
    detv = prod * real(sign_piv, kind=sp)
  end function det_from_lu
END SUBROUTINE Ovlps
