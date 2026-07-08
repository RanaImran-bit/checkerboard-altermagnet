!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

SUBROUTINE InitEnergy
use cpmc
use, intrinsic :: ieee_arithmetic
integer::i,j,k,l,m,n,ibas,jbas
! real(sp)::energy

real(sp)::ekin,ecoul,detu,detd
real(sp)::pnew_up(NUP,NUP),pnew_dn(NDN,NDN)
real(sp)::pg_up(NUP,NUP),pg_dn(NDN,NDN)
real(sp)::pg(NE,NE)
real(sp)::detpbas(NIW)
real(sp)::ewlkr(NIW),ebas(NIW),detwlkr(NIW)
real(sp)::det(2),work(nsites)
integer::info,ipvt(nsites)
external :: dgemm,dgefa,dgedi

!g = <a^{\dagger}a> = I-G = [R(LR)^{-1}L]/det(LR) /=0!
!g = computing adjoin matrix det(LR) =0 or /=0!

estptop=ZERO
estpbtm=ZERO
ewlkr=0.0
detwlkr=0.0
do 600 ibas=1,NIW

  do 700 jbas=1,NIW

   ebas(jbas)=ZERO
   if (use_spinor) then
     call dgemm('T', 'N', NE, NE, nsites, 1.0_8, phiB(:,:,ibas), NSO, phiB(:,:,jbas), NSO, 0.0_8, ovlpINV, NE) !li
     pnew_up=ovlpINV

     !compute determinant of overlap matrix LR!
     call dgefa(pnew_up,NUP,NUP,ipvt,info)
     if(info /= 0) then
     detp_up(jbas)=0.0
     else
     call dgedi(pnew_up,NUP,NUP,ipvt,det,work,10_i4b)
     detp_up(jbas)=det(1)*TEN**det(2)
     end if

     !compute g(R[(LR)^{-1}]L)!
     pg_up=0.0
     do i=1,NUP
       do j=1,NUP
       pnew_up=ovlpINV_up
       pnew_up(j,:)=0.0; pnew_up(:,i)=0.0; pnew_up(j,i)=1.0
       call dgefa(pnew_up,NUP,NUP,ipvt,info)
       if(info /= 0) then
        detu=0.0
       else
        call dgedi(pnew_up,NUP,NUP,ipvt,det,work,10_i4b)
        detu=det(1)*TEN**det(2)
       end if
       pg_up(i,j)=detu
       end do
     end do

     call dgemm('N', 'N', nsites, NUP, NUP, 1.0_8, phiB_up(:,:,jbas), nsites, pg_up, NUP, 0.0_8, tmp_up, nsites) !li
     call dgemm('N', 'T', nsites, nsites, NUP, -1.0_8, tmp_up, nsites, phiB_up(:,:,ibas), nsites, 0.0_8, gx_up(:,:), nsites)!li
     do i=1,nsites
     gx_up(i,i)=detp_up(jbas)+gx_up(i,i)
     end do

   else

   !----UP----compute overlap matrix LR, etc. for spin up!

   call dgemm('T', 'N', NUP, NUP, nsites, 1.0_8, phiB_up(:,:,ibas), nsites, phiB_up(:,:,jbas), nsites, 0.0_8, ovlpINV_up, NUP) !li

   pnew_up=ovlpINV_up

   !compute determinant of overlap matrix LR!
   call dgefa(pnew_up,NUP,NUP,ipvt,info)
   if(info /= 0) then
   detp_up(jbas)=0.0
   else
   call dgedi(pnew_up,NUP,NUP,ipvt,det,work,10_i4b)
   detp_up(jbas)=det(1)*TEN**det(2)
   end if

   !compute g(R[(LR)^{-1}]L)!
   pg_up=0.0
   do i=1,NUP
     do j=1,NUP
     pnew_up=ovlpINV_up
     pnew_up(j,:)=0.0; pnew_up(:,i)=0.0; pnew_up(j,i)=1.0
     call dgefa(pnew_up,NUP,NUP,ipvt,info)
     if(info /= 0) then
      detu=0.0
     else
      call dgedi(pnew_up,NUP,NUP,ipvt,det,work,10_i4b)
      detu=det(1)*TEN**det(2)
     end if
     pg_up(i,j)=detu
     end do
   end do

   call dgemm('N', 'N', nsites, NUP, NUP, 1.0_8, phiB_up(:,:,jbas), nsites, pg_up, NUP, 0.0_8, tmp_up, nsites) !li

   call dgemm('N', 'T', nsites, nsites, NUP, -1.0_8, tmp_up, nsites, phiB_up(:,:,ibas), nsites, 0.0_8, gx_up(:,:), nsites)!li

   do i=1,nsites
   gx_up(i,i)=detp_up(jbas)+gx_up(i,i)
   end do

   !---DOWN---compute overlap matrix LR, etc. for spin down!

   call dgemm('T', 'N', NDN, NDN, nsites, 1.0_8, phiB_dn(:,:,ibas), nsites, phiB_dn(:,:,jbas), nsites, 0.0_8, ovlpINV_dn, NDN)!li check

   !ovlpINV_dn=matmul(transpose(phiB_dn(:,:,ibas)),phiB_dn(:,:,jbas))
   pnew_dn=ovlpINV_dn

   !compute inverse and determinant of overlap matrix LR!
   call dgefa(pnew_dn,NDN,NDN,ipvt,info)
   if(info /= 0) then
    detp_dn(jbas)=0.0
   else
    call dgedi(pnew_dn,NDN,NDN,ipvt,det,work,10_i4b)
    detp_dn(jbas)=det(1)*TEN**det(2)
   end if

   !compute g (R[(LR)^1]L)!
   pg_dn=0.0
   do i=1,NDN
     do j=1,NDN
     pnew_dn=ovlpINV_dn
     pnew_dn(j,:)=0.0; pnew_dn(:,i)=0.0; pnew_dn(j,i)=1.0
     call dgefa(pnew_dn,NDN,NDN,ipvt,info)
     if(info /= 0) then
      detd=0.0
     else
      call dgedi(pnew_dn,NDN,NDN,ipvt,det,work,10_i4b)
      detd=det(1)*TEN**det(2)
     end if
     pg_dn(i,j)=detd
     end do
   end do

   call dgemm('N', 'N', nsites, NDN, NDN, 1.0_8, phiB_dn(:,:,jbas), nsites, pg_dn, NDN, 0.0_8, tmp_dn, nsites) !li

   !tmp_dn=matmul(phiB_dn(:,:,jbas),pg_dn)

   call dgemm('N', 'T', nsites, nsites, NDN, -1.0_8, tmp_dn, nsites, phiB_dn(:,:,ibas), nsites, 0.0_8, gx_dn(:,:), nsites) !li

   !gx_dn(:,:)=-matmul(tmp_dn,transpose(phiB_dn(:,:,ibas)))

   do i=1,nsites
    gx_dn(i,i)=detp_dn(jbas)+gx_dn(i,i)
   end do

   detpbas(jbas)=cwibas(jbas)*detp_up(jbas)*detp_dn(jbas)

 endif ! spinor

   call PairEnergy_Accumulate(ibas, jbas)

700      continue
   detwlkr(ibas)=cwibas(ibas)*sum(detpbas)
   ewlkr(ibas)=cwibas(ibas)*sum(ebas)
600   continue

estptop=sum(ewlkr)
estpbtm=sum(detwlkr)

etrial=estptop/(estpbtm)
e_var=etrial
write(*,*)estptop,estpbtm,etrial

write(*,*)'Initial Energy=',etrial
write(*,*)'Initial Energy Per N=',etrial/lxy
write(*,*)'###########################################'
write(*,*)

return
END subroutine InitEnergy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

SUBROUTINE EstEtrial
use cpmc
use jiekou,only:Step,Stblz,comb
! ,Step1,Stblz1
integer::i,j,k,l,m,n
integer::istp,iw
real(sp)::egrowth,w_up,w_down,wgtsum
external :: step1,stblz1
write(*,*) ' '
write(*,*) 'Growth Estimation:'
write(*,'((a),t11,(a),t28,(a),t38,(a))')&
      & ' Iblk','Egrowth', 'W_Up', 'W_Down'

egrowth=ZERO
do i=1,nblkgr*8
   wgtwlkr=(real(NWLKRS)/sum(wgtwlkr))*wgtwlkr
   w_up=sum(wgtwlkr)

   do istp=1,itvlorth
      call Step1(istp)
   end do

   w_down=sum(wgtwlkr)

   call Stblz1

   call comb

   egrowth=egrowth+w_up/w_down
   write(*,'(i5,3e12.4)') i,egrowth,w_up,w_down

end do

etrial=etrial+log(egrowth/real(nblkgr*8))/deltau/real(itvlorth)

write(*,*)
write(*,*)'Growth estimate of Etrial: ',etrial
write(*,*)'Growth estimate of Etrial Per N: ',etrial/lxy
write(*,*)

return
end subroutine EstEtrial


!===============================================================
! Compute and accumulate one- and two-body energy for (ibas,jbas)
! Adapts to spinor (use_spinor=.true.) and legacy up/down paths.
! This subroutine updates ebas(jbas) and detpbas(jbas) as needed.
!===============================================================
subroutine PairEnergy_Accumulate(ibas, jbas)
  use cpmc
  implicit none
  integer, intent(in) :: ibas, jbas
  integer :: i, j, info
  real(sp) :: ekin, ecoul

  if (use_spinor) then
     integer :: NSOloc, NEloc
     real(sp), allocatable :: A(:,:), invA(:,:), tmp(:,:), Gs(:,:)
     integer, allocatable :: ipiv(:)
     real(sp), allocatable :: work(:)
     real(sp) :: detval
     integer :: iu, ju, id, jd, lwork

     NSOloc = NSO
     NEloc  = NE

     allocate(A(NEloc,NEloc), invA(NEloc,NEloc), tmp(NSOloc,NEloc), Gs(NSOloc,NSOloc))
     allocate(ipiv(NEloc))
     lwork = max(4*NEloc, 64)
     allocate(work(lwork))

     ! A = L^T R = phiB(:,:,ibas)^T * phiB(:,:,jbas)
     call dgemm('T','N', NEloc, NEloc, NSOloc, 1.0_8, phiB(:,:,ibas), NSOloc, phiB(:,:,jbas), NSOloc, 0.0_8, A, NEloc)
     invA = A
     call dgetrf(NEloc, NEloc, invA, NEloc, ipiv, info)
     if (info /= 0) then
        detval = 0.0_sp
        Gs = 0.0_sp
     else
        ! determinant from LU (with pivot sign)
        detval = 1.0_sp
        do i=1, NEloc
           detval = detval * invA(i,i)
           if (ipiv(i) /= i) detval = -detval
        end do
        ! inverse
        call dgetri(NEloc, invA, NEloc, ipiv, work, lwork, info)
        ! Gs = R * (L^T R)^{-1} * L^T
        call dgemm('N','N', NSOloc, NEloc, NEloc, 1.0_8, phiB(:,:,jbas), NSOloc, invA, NEloc, 0.0_8, tmp, NSOloc)
        call dgemm('N','T', NSOloc, NSOloc, NEloc, 1.0_8, tmp, NSOloc, phiB(:,:,ibas), NSOloc, 0.0_8, Gs, NSOloc)
     end if

     detpbas(jbas) = cwibas(jbas) * detval


  end if

end subroutine PairEnergy_Accumulate
