subroutine Initpair
  use cpmc
  !----------------------------------------------------------------!
  ! Set up pairing symmetry and relative position of nearest site:'!
  !----------------------------------------------------------------!

  sf(1)=1.0; sf(2)=1.0; sf(3)=1.0;  sf(4)=1.0
  df(1)=1.0; df(2)=1.0; df(3)=-1.0; df(4)=-1.0
  pf(1)=1.0; pf(2)=-1.0; pf(3)=0.0; pf(4)=0.0
  ppf(1)=0.0; ppf(2)=0.0; ppf(3)=1.0; ppf(4)=-1.
  ddf(1)=1.0;ddf(2)=1.0;ddf(3)=-1.0;ddf(4)=-1.0;

  ipx(1)=1.0; ipy(1)=0.0
  ipx(2)=-1.0;ipy(2)=0.0
  ipx(3)=0.0; ipy(3)=1.0
  ipx(4)=0.0; ipy(4)=-1.0

  ippx(1)=1.0; ippy(1)=1.0
  ippx(2)=-1.0;ippy(2)=-1.0
  ippx(3)=-1.0; ippy(3)=1.0
  ippx(4)=1.0; ippy(4)=-1.0

  idx(1)=2.0; idy(1)=0.0
  idx(2)=-2.0;idy(2)=0.0
  idx(3)=0.0; idy(3)=2.0
  idx(4)=0.0; idy(4)=-2.0


endsubroutine Initpair

subroutine Initk
  use cpmc
  integer::i,j,k

  do i=1,lx
  xk(i)=real(i-1)*2.0*PI/real(lx)
  end do

  do j=1,ly
  yk(j)=real(j-1)*2.0*PI/real(ly)
  end do

  k=0
  do i=1,lx
  do j=1,ly
  k=k+1
  kSet(k,1)=xk(i)
  kSet(k,2)=yk(j)
  end do
  end do
endsubroutine Initk


subroutine correl
use cpmc
use jiekou,only:pair,backphi
real(sp)::phip_up(NUP,nsites),phip_dn(NDN,nsites)
integer::i,j,k,l,m,n,iu,i1,i2,j1,dx,dy,iw,kp,i_a,i_b
integer::m1,m2,m3,m4,n1,n2,n3,n4,ibas,iave
integer::info,ipvt(nsites)
real(sp)::det(2),work(nsites)
external :: dgemm,dgefa,dgedi,FourierTransform

real(sp) :: oy_total = 0.0
integer :: oy_count = 0
!-------------------!
! Quantum Averages  !
!-------------------!

wgt_p=0.0
corwlkr=0.0
ave=0.0

do 1100 iw=1,NWLKRS
if(wgtwlkr(iw) /= zero ) then   ! sum over walkers.

wgt_c(iw)=wgtwlkr(iw)
kp=iw

do 1000 ibas=1,NWFBAS
!--------------------------------------------------------------------------!
! Get < Phi(k',l) |; up & down; phip_up( NUP x nsites ) Transposed already!
!--------------------------------------------------------------------------!
call backphi(kp,ibas,phip_up,phip_dn)

!----------------------------------------------!
! Get Green's Function and Overlap; up and down!
!----------------------------------------------!

!------------------------------------------------------------!
! For Spin Up                                                !
! L = < Phi(k',l) |; R = | Phi(k) > at step N(saved); k == iw!
! Compute overlap matrix LR, etc.                            !
!------------------------------------------------------------!

call dgemm('N', 'N', NUP, NUP, nsites, 1.0_8, phip_up, NUP, phi_cup(:,:,iw), nsites, 0.0_8, ovlpINV_up, NUP) !li

!ovlpINV_up=matmul(phip_up,phi_cup(:,:,iw))

!-----------------------------------------------------!
! compute inverse and determinant of overlap matrix LR!
!-----------------------------------------------------!

call dgefa(ovlpINV_up,NUP,NUP,ipvt,info)
if(info /= 0) stop 'Problem in dgefa routine'

call dgedi(ovlpINV_up,NUP,NUP,ipvt,det,work,11_i4b)
detp_up(ibas)=det(1)*TEN**det(2)

!--------------------------!
! compute g (R[(LR)^{-1}]L)!
!--------------------------!
call dgemm('N', 'N', nsites, NUP, NUP, 1.0_8, phi_cup(:,:,iw), nsites, ovlpINV_up, NUP, 0.0_8, tmp_up, nsites) !li
!tmp_up=matmul(phi_cup(:,:,iw),ovlpINV_up)

call dgemm('N', 'N', nsites, nsites, NUP, 1.0_8, tmp_up, nsites, phip_up, NUP, 0.0_8, gx_up, nsites)
!gx_up=matmul(tmp_up,phip_up)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! For Spin Down                                              !
! L = < Phi(k',l) |; R = | Phi(k) > at step N(saved); k == iw!
! Compute overlap matrix LR, etc.                            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

call dgemm('N', 'N', NDN, NDN, nsites, 1.0_8, phip_dn, NDN, phi_cdn(:,:,iw), nsites, 0.0_8, ovlpINV_dn, NDN) !li

!ovlpINV_dn=matmul(phip_dn,phi_cdn(:,:,iw))

!-----------------------------------------------------!
! compute inverse and determinant of overlap matrix LR!
!-----------------------------------------------------!

call dgefa(ovlpINV_dn,NDN,NDN,ipvt,info)
if(info /= 0) stop 'Problem in dgefa routine'
call dgedi(ovlpINV_dn,NDN,NDN,ipvt,det,work,11_i4b)
detp_dn(ibas)=det(1)*TEN**det(2)

!--------------------------!
! compute g (R[(LR)^{-1}]L)!
!--------------------------!

call dgemm('N', 'N', nsites, NDN, NDN, 1.0_8, phi_cdn(:,:,iw), nsites, ovlpINV_dn, NDN, 0.0_8, tmp_dn, nsites) !li

!tmp_dn=matmul(phi_cdn(:,:,iw),ovlpINV_dn)

call dgemm('N', 'N', nsites, nsites, NDN, 1.0_8, tmp_dn, nsites, phip_dn, NDN, 0.0_8, gx_dn, nsites)
!gx_dn=matmul(tmp_dn,phip_dn)

!------------------------------------------------!
! Here need change for more than one determinants!
!------------------------------------------------!

ave(ibas,0,:,:)=1.0

!--------------------------------------------------------!
! Averages, Wick's theorem, etc.  Take the following form!
!--------------------------------------------------------!
! < c_i^\dagger c_j > = g(j,i) * <Psi(l) | Phi> (overlap)!
!                    -> g(j,i)                           !
!--------------------------------------------------------!

do iave=1,NAVE
ave_cor(iave,:,:)=zero
end do
opn_up=0.0
opn_dn=0.0
cdw=0.0
pmdf=0.0
dsf=0.0
cdwAvg=0.0

have=0.0; hkin=0.0; hpot=0.0;
sdwz=0.0; sdwx=0.0; cdw_re=0.0
nk_up=0.0; nk_dn=0.0; nk_w=0.0
pmdf_re=0.0
dsf_re=0.0

!-----------------------------------------------------------------------!
! < n(i,s)n(j,s') >                                                     !
! < n(i,+)n(j,+) > = < c^\dagger(i,+) c(i,+) > < c^\dagger(j,+) c(j,+) >!
!                  + < c^\dagger(i,+) c(j,+) > < c(i,+) c^\dagger(j,+) >!
! < n(i,+)n(j,-) > = < c^\dagger(i,+) c(i,+) > < c^\dagger(j,-) c(j,-) >!
!                  + < c^\dagger(i,+) c(j,+) > < c(i,+) c^\dagger(j,+) >!
!-----------------------------------------------------------------------!

do i=1,nsites            ! diagonal
    opn_up(i)=opn_up(i)+gx_up(i,i)
    opn_dn(i)=opn_dn(i)+gx_dn(i,i)
end do

do i=1,nsites
    do j=1,nsites
        up_up=gx_up(i,i)*gx_up(j,j)-gx_up(j,i)*gx_up(i,j)
        dn_dn=gx_dn(i,i)*gx_dn(j,j)-gx_dn(j,i)*gx_dn(i,j)
        up_dn=gx_up(i,i)*gx_dn(j,j)+gx_up(j,j)*gx_dn(i,i) ! 2

        n_up(i,j)=gx_up(i,j)
        n_dn(i,j)=gx_dn(i,j)
        n_w(i,j)=gx_up(i,j)+gx_dn(i,j)

        if (i==j) then
            pmdf_re(i,j)=pmdf_re(i,j)+gx_up(j,i)*gx_dn(j,i)
            dsf_re(i,j)=dsf_re(i,j)+gx_up(j,i)*gx_dn(j,i)

            cdw_re(i,j)=opn_up(i)+opn_dn(i)+2.0*gx_up(i,i)*gx_dn(i,i)+&
            & dens*(dens-2.0*opn_up(i)-2.0*opn_dn(i))        ! CDW
            sdwz_re(i,j)=(opn_up(i)+opn_dn(i)-2.0*gx_up(i,i)*gx_dn(i,i))*0.25  ! SDW-z
            sdwx_re(i,j)=(opn_up(i)+opn_dn(i)-2.0*gx_up(i,i)*gx_dn(i,i))*0.25  ! SDW-x

        else
            pmdf_re(i,j)=pmdf_re(i,j)+gx_up(j,i)*gx_dn(j,i)
            dsf_re(i,j)=dsf_re(i,j)+(gx_up(i,i)*gx_up(j,j)+gx_up(j,i)*gx_up(i,j))*&
            &(gx_dn(i,i)*gx_dn(j,j)+gx_dn(j,i)*gx_dn(i,j))

            cdw_re(i,j)=up_up+dn_dn+up_dn+dens*(dens-opn_up(i)-opn_dn(i)-&
            & opn_up(j)-opn_dn(j))     ! CDW
            sdwz_re(i,j)=(up_up+dn_dn-up_dn)*0.25    ! SDW-z
            sdwx_re(i,j)=(-gx_up(j,i)*gx_dn(i,j)-&
            & gx_up(i,j)*gx_dn(j,i))*0.25 ! SDW-x

        endif
    end do
end do

cu_up=sum(opn_up(1:nsites))/real(nsites)
cu_dn=sum(opn_dn(1:nsites))/real(nsites)

z2_sum=0.0
do i=1,nsites
    z2_sum=z2_sum+opn_up(i)+opn_dn(i)-2.0*gx_up(i,i)*gx_dn(i,i)
end do
z2_sum=z2_sum/real(nsites)

do j=1,nsites
    do i=1,nsites
        if (use_spinor) then
            hkin=hkin+tk(j,i)*gx(i,j)
        else
            hkin=hkin+tk_ud(j,i,1)*gx_up(i,j)+tk_ud(j,i,2)*gx_dn(i,j)
        end if 
    end do
end do
do i=1,nsites
    hpot=hpot+Vlist(i,1)*gx_up(i,i)*gx_dn(i,i)
end do

if(abs(v)>0.001) then
    do i=1,nsites
        do j1=1,2
            j=idis(i,j1)
            hpot=hpot+v*(gx_up(i,i)*gx_up(j,j)-gx_up(i,j)*gx_up(j,i))&
            & +v*(gx_dn(i,i)*gx_dn(j,j)-gx_dn(i,j)*gx_dn(j,i))&
            & +v*(gx_up(i,i)*gx_dn(j,j)+gx_up(j,j)*gx_dn(i,i))
        end do
    end do
end if
have=hkin+hpot

cdwAvg=sum([(gx_up(i,i)*gx_dn(i,i),i=1,nsites)])/nsites
do m = 1, nsites
do n = 1, nsites
dsf_re(n,m)=dsf_re(n,m)-cdwAvg*(gx_up(n,n)*gx_dn(n,n)+gx_up(m,m)*gx_dn(m,m))+cdwAvg**2
end do
end do

!     dsf(k)=<(\rho_b(r)-\rho_avg)  (\rho(r')>-\rho_avg)>
!            =<\rho_b(r) \rho_b(r')> -\bar{\rho} (<\rho_b(r')> + <\rho_b(r)>) + \bar{\rho}^2
!       <\rho_b(r')>= double occupancy at r'= cdw(r')

!***************************************************!
!  compute pairing correlation functions!
!---------------------------------------------------!

call pair(gx_up,gx_dn)

call FourierTransform(gx_up,gxk_up)
call FourierTransform(gx_dn,gxk_dn)
call FourierTransform(n_up,nk_up)
call FourierTransform(n_dn,nk_dn)
call FourierTransform(n_w,nk_w)
call FourierTransform(cdw_re,cdw)
call FourierTransform(sdwz_re,sdwz)
call FourierTransform(sdwx_re,sdwx)
call FourierTransform(pmdf_re,pmdf)
call FourierTransform(dsf_re,dsf)

do i_a=1,NLA
    do i_b=1,NLA
        k=0
        do m=1,NSTATES
            k=k+1
            ave_cor(k+NSTATES*0,i_a,i_b)=nk_up(m,i_a,i_b)
            ave_cor(k+NSTATES*1,i_a,i_b)=nk_dn(m,i_a,i_b)
            ave_cor(k+NSTATES*2,i_a,i_b)=nk_w(m,i_a,i_b)
            ave_cor(k+NSTATES*3,i_a,i_b)=cdw(m,i_a,i_b)
            ave_cor(k+NSTATES*4,i_a,i_b)=sdwz(m,i_a,i_b)
            ave_cor(k+NSTATES*5,i_a,i_b)=sdwx(m,i_a,i_b)
            ave_cor(k+NSTATES*6,i_a,i_b)=pmdf(m,i_a,i_b)
            ave_cor(k+NSTATES*7,i_a,i_b)=dsf(m,i_a,i_b)
            ave_cor(k+NSTATES*8,i_a,i_b)=swave(m,i_a,i_b)
            ave_cor(k+NSTATES*9,i_a,i_b)=pwave(m,i_a,i_b)
            ave_cor(k+NSTATES*10,i_a,i_b)=dwave(m,i_a,i_b)
            ave_cor(k+NSTATES*11,i_a,i_b)=sowave(m,i_a,i_b)
            ave_cor(k+NSTATES*12,i_a,i_b)=sbwave(m,i_a,i_b)
            ave_cor(k+NSTATES*13,i_a,i_b)=pbwave(m,i_a,i_b)
            ave_cor(k+NSTATES*14,i_a,i_b)=dbwave(m,i_a,i_b)
            ave_cor(k+NSTATES*15,i_a,i_b)=pdsfbd1(m,i_a,i_b)
            ave_cor(k+NSTATES*16,i_a,i_b)=pdsfbd2(m,i_a,i_b)
            ave_cor(k+NSTATES*17,i_a,i_b)=pdsfbd12(m,i_a,i_b)
            ave_cor(k+NSTATES*18,i_a,i_b)=dd1wave(m,i_a,i_b)
            ave_cor(k+NSTATES*19,i_a,i_b)=dd2wave(m,i_a,i_b)
            ave_cor(k+NSTATES*20,i_a,i_b)=dd12wave(m,i_a,i_b)
            ave_cor(k+NSTATES*21,i_a,i_b)=puupxwave(m,i_a,i_b)
            ave_cor(k+NSTATES*22,i_a,i_b)=pddpxwave(m,i_a,i_b)
            ave_cor(k+NSTATES*23,i_a,i_b)=pudpxwave(m,i_a,i_b)
            ave_cor(k+NSTATES*24,i_a,i_b)=puupywave(m,i_a,i_b)
            ave_cor(k+NSTATES*25,i_a,i_b)=pddpywave(m,i_a,i_b)
            ave_cor(k+NSTATES*26,i_a,i_b)=pudpywave(m,i_a,i_b)
            ave_cor(k+NSTATES*27,i_a,i_b)=pbdx2y2d1(m,i_a,i_b)
            ave_cor(k+NSTATES*28,i_a,i_b)=pbdx2y2d2(m,i_a,i_b)
            ave_cor(k+NSTATES*29,i_a,i_b)=pbdx2y2d12(m,i_a,i_b)
        end do

        ave_cor(NAVE-5,i_a,i_b)=cu_up
        ave_cor(NAVE-4,i_a,i_b)=cu_dn
        ave_cor(NAVE-3,i_a,i_b)=z2_sum
        ave_cor(NAVE-2,i_a,i_b)=hkin
        ave_cor(NAVE-1,i_a,i_b)=hpot
        ave_cor(NAVE,i_a,i_b)=have

    end do
end do

do i_a=1,NLA
    do i_b=1,NLA
        k=NSTATES*ntypes
        do m=1,nsites
            do n=1,nsites
                if (sublatt(m)==i_a) then
                    if (sublatt(n)==i_b) then
                        k=k+1
                        ave_cor(k+NSTATES*NSTATES*0,sublatt(m),sublatt(n))=n_up(m,n)
                        ave_cor(k+NSTATES*NSTATES*1,sublatt(m),sublatt(n))=n_dn(m,n)
                        ave_cor(k+NSTATES*NSTATES*2,sublatt(m),sublatt(n))=n_w(m,n)
                        ave_cor(k+NSTATES*NSTATES*3,sublatt(m),sublatt(n))=cdw_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*4,sublatt(m),sublatt(n))=sdwz_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*5,sublatt(m),sublatt(n))=sdwx_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*6,sublatt(m),sublatt(n))=pmdf_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*7,sublatt(m),sublatt(n))=dsf_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*8,sublatt(m),sublatt(n))=swave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*9,sublatt(m),sublatt(n))=pwave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*10,sublatt(m),sublatt(n))=dwave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*11,sublatt(m),sublatt(n))=sowave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*12,sublatt(m),sublatt(n))=sbwave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*13,sublatt(m),sublatt(n))=pbwave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*14,sublatt(m),sublatt(n))=dbwave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*15,sublatt(m),sublatt(n))=pdsfbd1_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*16,sublatt(m),sublatt(n))=pdsfbd2_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*17,sublatt(m),sublatt(n))=pdsfbd12_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*18,sublatt(m),sublatt(n))=dd1wave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*19,sublatt(m),sublatt(n))=dd2wave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*20,sublatt(m),sublatt(n))=dd12wave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*21,sublatt(m),sublatt(n))=puupxwave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*22,sublatt(m),sublatt(n))=pddpxwave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*23,sublatt(m),sublatt(n))=pudpxwave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*24,sublatt(m),sublatt(n))=puupywave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*25,sublatt(m),sublatt(n))=pddpywave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*26,sublatt(m),sublatt(n))=pudpywave_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*27,sublatt(m),sublatt(n))=pbdx2y2d1_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*28,sublatt(m),sublatt(n))=pbdx2y2d2_re(m,n)
                        ave_cor(k+NSTATES*NSTATES*29,sublatt(m),sublatt(n))=pbdx2y2d12_re(m,n)
                    end if
                end if
            end do
        end do
    end do
end do

do i_a=1,NLA
    do i_b=1,NLA
        do iave=1,NAVE
            ave(ibas,iave,i_a,i_b)=ave_cor(iave,i_a,i_b)*ave(ibas,0,i_a,i_b)
        end do
    end do
end do

1000    continue                             ! loop over ibas

do i_a=1,NLA
    do i_b=1,NLA
        do iave=0,nave
            wgt_p(iave,i_a,i_b)=zero
            do ibas=1,NWFBAS
                wgt_p(iave,i_a,i_b)=wgt_p(iave,i_a,i_b)+ave(ibas,iave,i_a,i_b)
            end do                            ! loop over ibas & kp
        end do

        !*****************************************************************!
        !------------------------------------------------!
        ! Here need change for more than one determinants!
        !------------------------------------------------!

        do iave=0,nave
            corwlkr(iave,i_a,i_b)=wgt_c(iw)*wgt_p(iave,i_a,i_b)
        end do

        do iave=0,nave                         ! 0 for bottom
            corptop(iave,i_a,i_b)=corptop(iave,i_a,i_b)+corwlkr(iave,i_a,i_b)
        end do

    end do
end do
end if
1100    continue                                ! loop over k (==iw)

nstpcor=nstpcor+1

end subroutine correl


subroutine pair(gp_up, gp_dn)
use cpmc
implicit none
integer :: i1, i2, k, iu, m, n, m1, n1, mm, nn, ii, jj, i3, j3, m3, n3, i_a, i_b
real(sp), intent(in) :: gp_up(nsites,nsites), gp_dn(nsites,nsites)
external :: FourierTransform

do m=1,nsites
    do n=1,nsites

        swave_re(m,n)       = 0.0
        dwave_re(m,n)       = 0.0
        pwave_re(m,n)       = 0.0
        sowave_re(m,n)      = 0.0
        sbwave_re(m,n)      = 0.0
        dbwave_re(m,n)      = 0.0
        pbwave_re(m,n)      = 0.0
        puupxwave_re(m,n)   = 0.0
        pddpywave_re(m,n)   = 0.0
        pudpxwave_re(m,n)   = 0.0
        puupywave_re(m,n)   = 0.0
        pddpxwave_re(m,n)   = 0.0
        pudpywave_re(m,n)   = 0.0
        pdsfbd1_re(m,n)     = 0.0
        pdsfbd2_re(m,n)     = 0.0
        pdsfbd12_re(m,n)    = 0.0
        dd1wave_re(m,n)     = 0.0
        dd2wave_re(m,n)     = 0.0
        dd12wave_re(m,n)    = 0.0
        pbdx2y2d1_re(m,n)   = 0.0
        pbdx2y2d2_re(m,n)   = 0.0
        pbdx2y2d12_re(m,n)  = 0.0

        do i1=1,4
            m1=iposit(ixv(m)+ipx(i1),iyv(m)+ipy(i1),sublatt(m))

            do i2=1,4

                n1=iposit(ixv(n)+ipx(i2),iyv(n)+ipy(i2),sublatt(n))

                p_sdp=gp_up(m,n)*gp_dn(m1,n1)+gp_up(m1,n1)*gp_dn(m,n)
                swave_re(m,n)=swave_re(m,n)+sf(i1)*sf(i2)*p_sdp
                dwave_re(m,n)=dwave_re(m,n)+df(i1)*df(i2)*p_sdp
                pwave_re(m,n)=pwave_re(m,n)+pf(i1)*pf(i2)*p_sdp

                puu_sdp=2*(gp_up(m,n)*gp_up(m1,n1)-gp_up(m,n1)*gp_up(m1,n))
                pdd_sdp=2*(gp_dn(m,n)*gp_dn(m1,n1)-gp_dn(m,n1)*gp_dn(m1,n))
                puupxwave_re(m,n)=puupxwave_re(m,n)+pf(i1)*pf(i2)*puu_sdp
                pddpywave_re(m,n)=pddpywave_re(m,n)+ppf(i1)*ppf(i2)*pdd_sdp
                puupywave_re(m,n)=puupywave_re(m,n)+ppf(i1)*ppf(i2)*puu_sdp
                pddpxwave_re(m,n)=pddpxwave_re(m,n)+pf(i1)*pf(i2)*pdd_sdp
                pdsfb=2*(gp_up(n,m)*gp_up(n1,m1)-gp_up(n1,m)*gp_up(n,m1))*(gp_dn(n,m)*gp_dn(n1,m1)-gp_dn(n1,m)*gp_dn(n,m1))
                sbwave_re(m,n)=sbwave_re(m,n)+sf(i1)*sf(i2)*pdsfb
                dbwave_re(m,n)=dbwave_re(m,n)+df(i1)*df(i2)*pdsfb
                pbwave_re(m,n)=pbwave_re(m,n)+pf(i1)*pf(i2)*pdsfb
            end do
        end do
        sowave_re(m,n)=sowave_re(m,n)+gp_up(m,n)*gp_dn(m,n)

        !----------------------dsfb---------------!
        do ii=1,4
            do jj=1,4
                mm=iposit(ixv(m)+ippx(ii),iyv(m)+ippy(ii),sublatt(m))
                nn=iposit(ixv(n)+ippx(jj),iyv(n)+ippy(jj),sublatt(n))

                pdsfbd=2*(gp_up(n,m)*gp_up(nn,mm)-gp_up(nn,m)*gp_up(n,mm))*(gp_dn(n,m)*gp_dn(nn,mm)-gp_dn(nn,m)*gp_dn(n,mm))
                if(ddf(ii)*ddf(jj)==-1) then
                    pdsfbd1_re(m,n)=pdsfbd1_re(m,n)+ddf(ii)*ddf(jj)*pdsfbd
                else
                    pdsfbd2_re(m,n)=pdsfbd2_re(m,n)+ddf(ii)*ddf(jj)*pdsfbd
                end if

                p_ddp=gp_up(m,n)*gp_dn(mm,nn)+gp_up(mm,nn)*gp_dn(m,n)
                if(ddf(ii)*ddf(jj)==-1) then
                    dd1wave_re(m,n)=dd1wave_re(m,n)+ddf(ii)*ddf(jj)*p_ddp
                else
                    dd2wave_re(m,n)=dd2wave_re(m,n)+ddf(ii)*ddf(jj)*p_ddp
                end if
            end do
        end do

        !----------------------------dx2-y2-------!
        do i3=1,4
            do j3=1,4
                ! write(*,*) 'm,n,i3,j3',m,n,i3,j3
                ! write(*,*)'1,',ixv(m)+idx(i3),iyv(m)+idy(i3),sublatt(m)
                ! write(*,*)'2,',iposit(ixv(m)+idx(i3),iyv(m)+idy(i3),sublatt(m))

                m3=iposit(ixv(m)+idx(i3),iyv(m)+idy(i3),sublatt(m))
                n3=iposit(ixv(n)+idx(j3),iyv(n)+idy(j3),sublatt(n))
                pbdx2y2=2*(gp_up(n,m)*gp_up(n3,m3)-gp_up(n3,m)*gp_up(n,m3))*(gp_dn(n,m)*gp_dn(n3,m3)-gp_dn(n3,m)*gp_dn(n,m3))
                if(ddf(i3)*ddf(j3)==-1) then
                    pbdx2y2d1_re(m,n)=pbdx2y2d1_re(m,n)+ddf(i3)*ddf(j3)*pbdx2y2
                else
                    pbdx2y2d2_re(m,n)=pbdx2y2d2_re(m,n)+ddf(i3)*ddf(j3)*pbdx2y2
                end if
            end do
        end do

    end do
end do

! rVals->kVals
call FourierTransform(swave_re,swave)
call FourierTransform(dwave_re,dwave)
call FourierTransform(pwave_re,pwave)
call FourierTransform(sowave_re,sowave)
call FourierTransform(sbwave_re,sbwave)
call FourierTransform(dbwave_re,dbwave)
call FourierTransform(pbwave_re,pbwave)
call FourierTransform(puupxwave_re,puupxwave)
call FourierTransform(pddpywave_re,pddpywave)
call FourierTransform(puupywave_re,puupywave)
call FourierTransform(pddpxwave_re,pddpxwave)
call FourierTransform(pdsfbd1_re,pdsfbd1)
call FourierTransform(pdsfbd2_re,pdsfbd2)
call FourierTransform(dd1wave_re,dd1wave)
call FourierTransform(dd2wave_re,dd2wave)
call FourierTransform(pbdx2y2d1_re,pbdx2y2d1)
call FourierTransform(pbdx2y2d2_re,pbdx2y2d2)

do i_a=1,NLA
    do i_b=1,NLA
        do iu=1,NSTATES
            swave(iu,i_a,i_b)=swave(iu,i_a,i_b)/2
            dwave(iu,i_a,i_b)=dwave(iu,i_a,i_b)/2
            pwave(iu,i_a,i_b)=pwave(iu,i_a,i_b)/1.41421356
            sowave(iu,i_a,i_b)=sowave(iu,i_a,i_b)
            sbwave(iu,i_a,i_b)=sbwave(iu,i_a,i_b)/2
            dbwave(iu,i_a,i_b)=dbwave(iu,i_a,i_b)/2
            pbwave(iu,i_a,i_b)=pbwave(iu,i_a,i_b)/1.41421356
            puupxwave(iu,i_a,i_b)=puupxwave(iu,i_a,i_b)/1.41421356
            pddpywave(iu,i_a,i_b)=pddpywave(iu,i_a,i_b)/1.41421356
            puupywave(iu,i_a,i_b)=puupywave(iu,i_a,i_b)/1.41421356
            pddpxwave(iu,i_a,i_b)=pddpxwave(iu,i_a,i_b)/1.41421356
            pdsfbd1(iu,i_a,i_b)=pdsfbd1(iu,i_a,i_b)/2
            pdsfbd2(iu,i_a,i_b)=pdsfbd2(iu,i_a,i_b)/2
            dd1wave(iu,i_a,i_b)=dd1wave(iu,i_a,i_b)/2
            dd2wave(iu,i_a,i_b)=dd2wave(iu,i_a,i_b)/2
            pbdx2y2d1(iu,i_a,i_b)=pbdx2y2d1(iu,i_a,i_b)/2
            pbdx2y2d2(iu,i_a,i_b)=pbdx2y2d2(iu,i_a,i_b)/2
            pudpxwave(iu,i_a,i_b)=puupxwave(iu,i_a,i_b)+pddpxwave(iu,i_a,i_b)
            pudpywave(iu,i_a,i_b)=puupywave(iu,i_a,i_b)+pddpywave(iu,i_a,i_b)
            pdsfbd12(iu,i_a,i_b)=pdsfbd1(iu,i_a,i_b)+pdsfbd2(iu,i_a,i_b)
            dd12wave(iu,i_a,i_b)=dd1wave(iu,i_a,i_b)+dd2wave(iu,i_a,i_b)
            pbdx2y2d12(iu,i_a,i_b)=pbdx2y2d1(iu,i_a,i_b)+pbdx2y2d2(iu,i_a,i_b)
        end do
    end do
end do

return
end subroutine pair

subroutine FourierTransform(rval,kvalue)
use cpmc
real(sp), intent(in) :: rval(nsites,nsites)
real(sp), intent(inout) :: kvalue(NSTATES, NLA, NLA)
integer :: rxm, rym, rxn, ryn,iu,i,j,m,n,a,b
real(sp) :: phase,cc
integer :: dx(NLA), dy(NLA)
dx(1)= 0; dx(2)= 0.0
dy(1)= 0; dy(2)=0

iu = 0
kvalue = 0.0
do i = 1, lx
    do j = 1, ly
        iu = iu + 1
        do m = 1, nsites
            rxm = ixv(m) + dx(sublatt(m))
            rym = iyv(m) + dy(sublatt(m))
            b = sublatt(m)
            do n = 1, nsites
                rxn = ixv(n) + dx(sublatt(n))
                ryn = iyv(n) + dy(sublatt(n))
                a = sublatt(n)

                phase = xk(i) * (rxm - rxn) + yk(j) * (rym - ryn)
                cc = cos(phase)
                kvalue(iu, a, b) = kvalue(iu, a, b) + cc * rval(n, m)
            end do
        end do
    end do
end do

do a = 1, NLA
    do b = 1, NLA
        do iu = 1, NSTATES
            kvalue(iu, a, b) = kvalue(iu, a, b) / real(NSTATES, sp)
        end do
    end do
end do

end subroutine FourierTransform
