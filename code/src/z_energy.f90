subroutine calenergy(jbas, ebas, length)
   use cpmc
   implicit none
   integer, intent(in) :: jbas, length
   integer::i,j
   real(sp)::ekin,ecoul,detu,detd
   real(sp):: ebas(length) 
if (use_spinor) then
   ! One-body energy: use block-diagonal tk; if spin-flip one-body exists, extend here.
   ekin = 0.0_sp
   do j=1, NSO
      do i=1, NSO
         ekin = ekin - tk(j,i)*gx(i,j)
      end do
   end do

   ! Two-body (Hubbard-U on-site) using spinor densities
   ecoul = 0.0_sp
   do i=1, nsites
      ecoul = ecoul + Vlist(i,1) * gx(i,i) * gx(i+nsites,i+nsites)
   end do

   ebas(jbas) = ebas(jbas) + ekin + ecoul
   ebas(jbas) = ebas(jbas) * cwibas(jbas)
   ! deallocate(A, invA, tmp, Gs, ipiv, work)

else
   ekin = 0.0_sp
   do j=1, nsites
      do i=1, nsites
         ekin = ekin - tk_ud(j,i,1)*gx_up(i,j)*detp_dn(jbas) - tk_ud(j,i,2)*gx_dn(i,j)*detp_up(jbas)
      end do
   end do
   ebas(jbas) = ebas(jbas) + ekin

   ecoul = 0.0_sp
   do i=1, nsites
      ecoul = ecoul + Vlist(i,1)*gx_up(i,i)*gx_dn(i,i) - &
                     Vlist(i,1)*detp_up(jbas)*gx_dn(i,i) - &
                     Vlist(i,1)*detp_dn(jbas)*gx_up(i,i) + &
                     Vlist(i,1)*detp_up(jbas)*detp_dn(jbas)
   end do
   ebas(jbas) = ebas(jbas) + ecoul
   ebas(jbas) = ebas(jbas) * cwibas(jbas)

endif
end subroutine calenergy

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

SUBROUTINE InitEnergy
use cpmc
use, intrinsic :: ieee_arithmetic
integer::i,j,k,l,m,n,ibas,jbas
! real(sp)::energy

real(sp)::ecoul,detu,detd

real(sp)::detpbas(NIW)
real(sp)::ewlkr(NIW),ebas(NIW),detwlkr(NIW)
real(sp)::det(2),work(nsites)
integer::info,ipvt(nsites)
external :: dgemm,dgefa,dgedi
logical :: has_nan_in_detpbas 

!g = <a^{\dagger}a> = I-G = [R(LR)^{-1}L]/det(LR) /=0!
!g = computing adjoin matrix det(LR) =0 or /=0!

estptop=ZERO
estpbtm=ZERO
ewlkr=0.0
detwlkr=0.0
do 600 ibas=1,NIW
   ! detpbas(:) = ZERO    
   ebas(:)    = ZERO
  do 700 jbas=1,NIW
   call calgf(ibas, jbas, 1, detpbas, NIW)
   call calenergy(jbas, ebas, NIW)
   ! ebas(jbas)=ebas(jbas)*cwibas(jbas)

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
integer::i,j,k,l,m,n
integer::istp,iw
real(sp)::egrowth,w_up,w_down,wgtsum
write(*,*) ' '
write(*,*) 'Growth Estimation:'
write(*,'((a),t11,(a),t28,(a),t38,(a))')&
      & ' Iblk','Egrowth', 'W_Up', 'W_Down'

egrowth=ZERO
do i=1,nblkgr*8
   wgtwlkr=(real(NWLKRS)/sum(wgtwlkr))*wgtwlkr
   w_up=sum(wgtwlkr)

   do istp=1,itvlorth
      call Step(istp,2)
   end do

   w_down=sum(wgtwlkr)

   call Stblz(2)

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
