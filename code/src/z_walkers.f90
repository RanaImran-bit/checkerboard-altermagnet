!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

SUBROUTINE InitPop
use cpmc
integer::i,j,k,l,m,n
integer::ibas,istart,iend,istop,iw
real(sp)::cnorm,exch(nsites,1)

cnorm=0.0
do i=1,NIW
cnorm=cnorm+abs(cwibas(i))
end do

istart=1
do ibas=1,NIW
   istop=istart+nint(abs(cwibas(ibas)/cnorm)*NWLKRS)
   if(ibas == NIW) istop=NWLKRS
   iend=min(istop,NWLKRS)
   do iw=istart,iend

      if (use_spinor) then
         do l=1,NSO
            do k=1,NE
               phi(l,k,iw) = phiB(l,k,ibas)
            end do
         end do

         ! if (sgnINIT(ibas) < 0.0) then
         !    exch(:,1) = phi(:,1,iw)
         !    phi(:,1,iw) = phi(:,2,iw)
         !    phi(:,2,iw) = exch(:,1)
         ! end if

      else
        do l=1,nsites
           do k=1,NUP
              phi_up(l,k,iw)=phiB_up(l,k,ibas)
            end do
            do k=1,NDN
               phi_dn(l,k,iw)=phiB_dn(l,k,ibas)
            end do
      end do

      ! if(sgnINIT(ibas)<0.0) then
      ! exch(:,1)=phi_up(:,1,iw)
      ! phi_up(:,1,iw)=phi_up(:,2,iw)
      ! phi_up(:,2,iw)=exch(:,1)
      ! end if

      ovlpDET(iw)=ONE
      wgtwlkr(iw)=ONE
      sgn(iw)=ONE
    end if
   end do
   istart=iend+1
end do
return
END subroutine InitPop

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!


subroutine savephi
use cpmc
integer::iw,is,ip

  do 100 iw=1,NWLKRS
      if (use_spinor) then
         ! Save spinor walker matrices if cache exists
        phi_c(:,:,iw) = phi(:,:,iw)
      else
     phi_cup(:,:,iw) = phi_up(:,:,iw)
     phi_cdn(:,:,iw) = phi_dn(:,:,iw)
      end if
100     continue

return
end subroutine savephi

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

SUBROUTINE Comb
use cpmc
use jiekou,only:PopCopy,ranGen
integer::i,j,k,l,m,n
real(sp)::wgtTOT,wgtAVG,offset
real(sp)::psum(0:NWLKRS),wcomb(NWLKRS)
integer::iw
integer::imult(NWLKRS),iempty(NWLKRS)

!! form the partial sum over the weights
psum(0)=ZERO
psum(1)=wgtwlkr(1)
do iw=2,NWLKRS
   psum(iw)=psum(iw-1)+wgtwlkr(iw)
end do

wgtTOT=psum(NWLKRS)
wgtAVG=wgtTOT/real(NWLKRS)

!! make the comb and reset weights


do iw=1,NWLKRS
   offset=ranGen(ISEED)*wgtAVG
   wgtwlkr(iw)=wgtAVG
   wcomb(iw)=offset+real(iw-1)*wgtAVG
end do

!! make list marking walkers to be replicated or eliminated
k=1
do iw=1,NWLKRS
   imult(iw)=0
11      if(k <= NWLKRS) then
     if(psum(iw-1) < wcomb(k)) then
         if(wcomb(k) <= psum(iw)) then
            k=k+1
            imult(iw)=imult(iw)+1
            goto 11
         end if
      end if
   end if
end do

!! create a list marking memory gaps
k=0
do iw=1,NWLKRS
   if(imult(iw) == 0) then
      k=k+1
      iempty(k)=iw
   end if
end do

!! replicate walkers and fill in memory
k=0
do j=1,NWLKRS
   if(imult(j) > 1) then
      do l=2,imult(j)
         k=k+1
         call PopCopy(iempty(k),j)
      end do
   end if
end do

return
END subroutine Comb

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

SUBROUTINE PopCopy(iwempty,iw)
use cpmc
integer::i,j,k,l,m,n,chan
integer::iwempty,iw,my

  if (use_spinor) then
     phi(:,:,iwempty) = phi(:,:,iw)
  else
    phi_up(:,:,iwempty)=phi_up(:,:,iw)
    phi_dn(:,:,iwempty)=phi_dn(:,:,iw)
  end if

ovlpDET(iwempty)=ovlpDET(iw)
sgn(iwempty)=sgn(iw)

if(measl==1.and.mstep>0) then


   if (use_spinor) then
      phi_c(:,:,iwempty) = phi_c(:,:,iw)
   else
      phi_cup(:,:,iwempty)=phi_cup(:,:,iw)
       phi_cdn(:,:,iwempty)=phi_cdn(:,:,iw)
    endif

   do my=1,mstep       ! parent(mstep,iwempty)=parent(mstep,iw)
      do i=1,nsites
         do chan=1,channels
         kexpV_s(i,my,chan,iwempty)=kexpV_s(i,my,chan,iw)
         end do
      end do
   end do

   if(g_ph*w_ph>0.001) then
      do my=1,mstep       ! parent(mstep,iwempty)=parent(mstep,iw)
         do i=1,nsites
            keph_s(i,my,iwempty)=keph_s(i,my,iw)
         end do
      end do
   end if
end if

return
END subroutine PopCopy
