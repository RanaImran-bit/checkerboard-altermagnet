
subroutine stblz(mode)
  use cpmc
  use jiekou,only:modgs
  implicit none

  integer::iw,iw1,iw2,mode
  real(sp)::rescale
  if (mode==1) then
    iw1=iwStart(myID)
    iw2=iwEnd(myID)
  else
    iw1=1
    iw2=NWLKRS
  endif
  !-----------------HoKinPara-------------------!
  do iw = iw1, iw2
  !-----------------HoKinPara-------------------!
    !!skip walkers with zero weight.
    if(wgtwlkr(iw)/=0.0) then
      if (use_spinor) then
        call modgs(phi(:,:,iw),NSO,NE,rescale)
        ovlpDET(iw)=rescale*ovlpDET(iw)
      else
        call modgs(phi_up(:,:,iw),nsites,NUP,rescale)
        ovlpDET(iw)=rescale*ovlpDET(iw)
        call modgs(phi_dn(:,:,iw),nsites,NDN,rescale)
        ovlpDET(iw)=rescale*ovlpDET(iw)
      endif
    end if
  end do
  return
end subroutine stblz


subroutine stblzbk(kp,phi_1,phi_2)
	use cpmc
	use jiekou,only:modgs
	implicit none
  integer::kp
  real(sp),dimension(NUP,nsites)::phi_1
  real(sp),dimension(NDN,nsites)::phi_2
  real(sp)::rescale,tmpup(nsites,NUP),tmpdn(nsites,NDN)

    tmpup=transpose(phi_1)
    tmpdn=transpose(phi_2)

    call modgs(tmpup,nsites,NUP,rescale)

    call modgs(tmpdn,nsites,NDN,rescale)

    phi_1=transpose(tmpup)
    phi_2=transpose(tmpdn)

	return
end subroutine stblzbk

subroutine stblzbk_spinor(kp,phi_1)
	use cpmc
	use jiekou,only:modgs
	implicit none
  integer::kp
  real(sp),dimension(NE,NSO)::phi_1
  real(sp)::rescale,tmpud(NSO,NE)
    tmpud=transpose(phi_1)
    call modgs(tmpud,NSO,NE,rescale)
    phi_1=transpose(tmpud)
	return
end subroutine stblzbk_spinor


subroutine modgs(phiin,nss,npp,rescale)
  use cpmc
  implicit none

  real(sp)::rescale,tmp,hld,sdot,dnrm2
  integer::nss,npp,ip,jp,n1

  real(sp),dimension(:,:)::phiin
  external :: dnrm2

  rescale=1.0_dp
  n1=1

  do ip=1,npp

     tmp=1.0/dnrm2(nss,phiin(:,ip),n1)
     rescale=rescale*tmp
     phiin(:,ip)=tmp*phiin(:,ip)

     do jp=ip+1,npp
     hld=dot_product(phiin(:,ip),phiin(:,jp))
     phiin(:,jp)=phiin(:,jp)-hld*phiin(:,ip)
     end do
  end do

  return
end subroutine modgs
