!====================FIRST PART==========================!

SUBROUTINE SetUp
use cpmc
use jiekou,only:GetMtrxs,MkExpT,MkExpV
character*60 infile

!!!!!!READ INPUT PARAMETERS FROM in.dat.

open(8,file='in.dat',status='old')

read(8,*) nblkeq,nblkgr,nblk,nblkstps
read(8,*) itvlpceq,itvlpc,itvlorth,itvlmeas,itvl_m
read(8,*) deltau,etrial
read(8,*) infile
read(8,*) t0,t1,t2,tam,alphat1
read(8,*) ud,vpd
read(8,*) g_ph,w_ph
read(8,*) (cwfbas(ibas),ibas=1,NWFBAS)
read(8,*) (cwibas(ibas),ibas=1,NIW)

close(8)

if(NSPIN == 2) then
NELECs(1)=NUP
NELECs(2)=NDN
end if
dens=real(NUP+NDN)/real(NSTATES)

!!!!!!WRITE CONTROL PARAMETERS TO STDOUT
write(*,*)
write(*,*)'=================================================='
write(*,*)' CPMC: one-band Hubbard model'
write(*,*)' singlet (closed shell) ground-state'
write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(*,*)
write(*,*)'Parameters: t0,t1,ud,vpd'
write(*,*) t0,t1,ud,vpd
write(*,*)
write(*,*)'# of states in SCF solution: NSTATES=',NSTATES
write(*,*)'# of up (or down) electrons: NELEC=',NELEC,NUP,NDN
write(*,*)'# of walkers (expected): NWLKRS=',NWLKRS
write(*,*)'# of blocks for relaxation: NblkEq=',nblkeq
write(*,*)'# of blocks for growth estimate: NblkGr=',nblkgr
write(*,*)'# of blocks after relaxation: Nblk=',nblk
write(*,*)'# of steps in each block: NblkStps=',nblkstps
write(*,*)
write(*,'((a),f7.4)')' Trotter step size: DeltaTau=',deltau
write(*,'((a),f16.8)')' Growth-control energy estimate=',etrial
write(*,*)
write(*,*)'Measurement interval: itvlmeas=',itvlmeas
write(*,*)'Orthonormalization interval: itvlorth=',itvlorth
write(*,*)'Pop control interval (relaxation phase):  ','itvlpceq=',itvlpceq
write(*,*)'Pop control interval (measurement phase): ','itvlpc=',itvlpc
write(*,*)
write(*,*)'Back Propagation steps: ','itvl_m=',itvl_m
write(*,*)
write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(*,*)

!!!!!INPUT AND GENERATE VARIOUS MATRICIES
!     construct exp(-dt*T) and exp(-dt*V)
!      .. deltau/2 is used to symmetrize the kernel:
!!!!!exp(-0.5*deltau*T)*exp(-deltau*V)*exp(-0.5*deltau*T).

call GetMtrxs

call MkExpT(HALF*deltau)

call MkExpV      ! Make e^V for each spin

!-------------------------HoKinPara-------------------------!
      allocate(iwStart(0:noOfProc-1),iwEnd(0:noOfProc-1),iwCount(0:noOfProc-1)) 
      call StartEnd(noOfProc,1,NWLKRS,iwStart,iwEnd,iwCount)
      sizekexpV=NSTATES*M_MAX
      sizePhiUp=NSTATES*NUP
      sizePhiDn=NSTATES*NDN
!-------------------------HoKinPara-------------------------!

return
END subroutine SetUp

!**************************************************!

SUBROUTINE GetMtrxs
use cpmc

!! Set up coordinates for Cu sites, lattice spacing=1
do 10 ixx=1,lx
do 10 iyy=1,ly
i=(ixx-1)*ly+iyy
ixv(i)=ixx-1
iyv(i)=iyy-1
iposit( ixv(i),iyv(i) )=i
10      continue
nsites_cu=lx*ly


nsites=nsites_cu
if(nsites /= nstates) stop 'nsites != nstates'

!! Making Periodic Boundary Conditions
do 30 iyy=-ly,2*ly
do 30 ixx=-lx,2*lx
itmpx=mod(ixx+lx,lx)
itmpy=mod(iyy+ly,ly)
iposit(ixx,iyy)=iposit(itmpx,itmpy)
30      continue

!! Hoppings( Assuming Spin-independent ).
do 110 i=1,nsites
do 110 j=1,nsites

do 110 k=1,NSPIN
tk(i,j,k)=zero
110     continue

!! Cu-Cu-t0
do 120 i=1,nsites_cu
!-----------------HKANHop------------------------!
ttp=t0+tam
ttn=t0-tam

! 1-up 2-down
tk( i, iposit( ixv(i)+1,iyv(i) ),1 )=ttn
tk( i, iposit( ixv(i)-1,iyv(i) ),1 )=ttn
tk( i, iposit( ixv(i),iyv(i)+1 ),1 )=ttp
tk( i, iposit( ixv(i),iyv(i)-1 ),1 )=ttp


tk( i, iposit( ixv(i)+1,iyv(i) ),2 )=ttp
tk( i, iposit( ixv(i)-1,iyv(i) ),2 )=ttp
tk( i, iposit( ixv(i),iyv(i)+1 ),2 )=ttn
tk( i, iposit( ixv(i),iyv(i)-1 ),2 )=ttn

!tk( i, iposit( ixv(i)+1,iyv(i) ),1 )=t0
!tk( i, iposit( ixv(i)-1,iyv(i) ),1 )=t0
!tk( i, iposit( ixv(i),iyv(i)+1 ),1 )=alpha*t0
!tk( i, iposit( ixv(i),iyv(i)-1 ),1 )=alpha*t0

!tk( i, iposit( ixv(i)+1,iyv(i) ),2 )=alpha*t0
!tk( i, iposit( ixv(i)-1,iyv(i) ),2 )=alpha*t0
!tk( i, iposit( ixv(i),iyv(i)+1 ),2 )=t0
!tk( i, iposit( ixv(i),iyv(i)-1 ),2 )=t0

!-----------------HKANHop------------------------!
120     continue
!-----------------HK output------------------!
do 121 i=1,nsites_cu
        tk( i, iposit( ixv(i)+1,iyv(i)+1 ),1 )=t1
        tk( i, iposit( ixv(i)-1,iyv(i)+1 ),1 )=alphat1*t1
        tk( i, iposit( ixv(i)-1,iyv(i)-1 ),1 )=t1
        tk( i, iposit( ixv(i)+1,iyv(i)-1 ),1 )=alphat1*t1

        tk( i, iposit( ixv(i)+1,iyv(i)+1 ),2 )=alphat1*t1
        tk( i, iposit( ixv(i)-1,iyv(i)+1 ),2 )=t1
        tk( i, iposit( ixv(i)-1,iyv(i)-1 ),2 )=alphat1*t1
        tk( i, iposit( ixv(i)+1,iyv(i)-1 ),2 )=t1
121     continue

! write(23,'(16 f 10.5)') tk

!---------------------------------------------!



!------------------HK disable for low demension case-----------!
! Cu_Cu-t1
!      do 130 i=1,nsites_cu
!         tk( i, iposit( ixv(i)+1,iyv(i)+1 ) )=t1
!         tk( i, iposit( ixv(i)-1,iyv(i)+1 ) )=t1
!         tk( i, iposit( ixv(i)-1,iyv(i)-1 ) )=t1
!         tk( i, iposit( ixv(i)+1,iyv(i)-1 ) )=t1
! 130     continue

!        !! Cu_Cu-t2
!     do 131 i=1,nsites_cu
!        tk( i, iposit( ixv(i)+2,iyv(i) ) )=t2
!        tk( i, iposit( ixv(i)-2,iyv(i) ) )=t2
!        tk( i, iposit( ixv(i),iyv(i)+2 ) )=t2
!        tk( i, iposit( ixv(i),iyv(i)-2 ) )=t2
!131     continue

!---------------------------------------------------------------!
!! U(i) and E(i)
do i=1,nsites_cu
hub_u(i)=ud
epsil(i)=zero
end do

k=0
do i=1,lx
do j=1,ly
k=k+1
icorx(k)=i-1
icory(k)=j-1
ncor(i,j) = k
end do
end do
if(k/=NSTATES) stop 'PAIR NUMBER ERROR!'

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

!----------------------------!
! set cu-cu relative position.!
!----------------------------!

do i=1,nsites
idis(i,1)=iposit(ixv(i)+1,iyv(i))
idis(i,2)=iposit(ixv(i),iyv(i)+1)
end do
!===================================================!

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

!! Delta Function
do i=1,nsites
do j=1,nsites
deltaf(i,j)=zero
end do
deltaf(i,i)=one
end do

return
END subroutine GetMtrxs

!**************************************************!

SUBROUTINE MkExpT(xdeltau)
use cpmc
use jiekou,only:tred2,tql2,MkInitOvlps
integer::ierr, iRead, ierr2
real(sp)::xdeltau,sum
real(sp)::d(NSTATES),e(NSTATES),z(NSTATES,NSTATES)

!----------------------HKAnHop------------------------!

real(sp)::d2(NSTATES),e2(NSTATES),z2(NSTATES,NSTATES)

call tred2(NSTATES,NSTATES,tk(:,:,1),d,e,z)
call tred2(NSTATES,NSTATES,tk(:,:,2),d2,e2,z2)
call tql2(NSTATES,NSTATES,d,e,z,ierr)
call tql2(NSTATES,NSTATES,d2,e2,z2,ierr2)
if(ierr /= 0) write(*,*) 'Ierr = ',ierr
if(ierr2 /= 0) write(*,*) 'Ierr = ',ierr2

!----------------------HKAnHop------------------------!
write(*,*) 'Eigenvalues of One-Body Matrix tk'

write(*,*) 'Up-spin'
write(*,'(4e14.6)') d

write(*,*) 'Dn-spin'

write(*,'(4e14.6)') d2
write(*,*)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
write(*,*) 'Make varitional Slater Determinent'
!---------HoKinReading----------!

iRead=0

if (iRead==1) then

open(22,file='wfup.txt',status='old')
read(22,*) phiT_up(:,:,1)
phiB_up=phiT_up
close(22)

open(22,file='wfdn.txt',status='old')
read(22,*) phiT_dn(:,:,1)
phiB_dn=phiT_dn
close(22)


!-------------------------------!

else
do ip=1,NUP
do is=1,NSTATES
phiT_up(is,ip,1)=z(is,ip)
end do
end do

do ip=1,NDN
do is=1,NSTATES
phiT_dn(is,ip,1)=z2(is,ip)
end do
end do

do ip=1,NUP
do is=1,NSTATES
phiB_up(is,ip,1)=z(is,ip)
end do
end do

do ip=1,NDN
do is=1,NSTATES
phiB_dn(is,ip,1)=z2(is,ip)
end do
end do
endif
 

do ibas=1,NWFBAS
phiZ_up(:,:,ibas)=transpose(phiT_up(:,:,ibas))
phiZ_dn(:,:,ibas)=transpose(phiT_dn(:,:,ibas))
end do

!-----------HK Output------------!
!            write(*,'(16 f 10.5)') phiT_up(:,:,1)
!            write(*,'(16 f 10.5)') phiT_dn(:,:,1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

do i=1,NSTATES
d(i)=exp(-xdeltau*d(i))
end do

do i=1,NSTATES
do j=i,NSTATES
sum=zero
do k=1,NSTATES
sum=sum+z(j,k)*d(k)*z(i,k)
end do
expT(j,i,1)=sum
expT(i,j,1)=sum
end do
end do 

call dgemm('N', 'N', NSTATES, NSTATES, NSTATES, 1.0_8, expT(:,:,1), NSTATES, expT(:,:,1), NSTATES, 0.0_8, exp2T(:,:,1), NSTATES) !li 
!write(*,*) exp2T(:,:,1)
!exp2T(:,:,1)=matmul(expT(:,:,1),expT(:,:,1))
!write(*,*) exp2T(:,:,1)
!write(*,*) '755'

do i=1,NSTATES
d2(i)=exp(-xdeltau*d2(i))
end do

do i=1,NSTATES
do j=i,NSTATES
sum=zero
do k=1,NSTATES
sum=sum+z2(j,k)*d2(k)*z2(i,k)
end do
expT(j,i,2)=sum
expT(i,j,2)=sum
end do
end do

call dgemm('N', 'N', NSTATES, NSTATES, NSTATES, 1.0_8, expT(:,:,2), NSTATES, expT(:,:,2), NSTATES, 0.0_8, exp2T(:,:,2), NSTATES) !li

!exp2T(:,:,2)=matmul(expT(:,:,2),expT(:,:,2))


call MkInitOvlps

return
ENDsubroutine MkExpT

!**************************************************!

!~~~~~Initialize Final Averages: cor_top/btm
!     1: averages; 2: averages square; 3: top-bottom cross term;
!~~~~~(**4: standard deviation; 5: covariance;**) 6: error.

SUBROUTINE InitRunMeas
use cpmc
! allocate(super_xo((ntypes-10) * NSTATES * NSTATES))
! allocate(superk_xo((ntypes-10)*NSTATES))
! allocate(super_xd((ntypes-10) * NSTATES * NSTATES))
! allocate(superk_xd((ntypes-10)*NSTATES))


allocate(rVals(ntypes))
allocate(kVals(ntypes))
do i = 1, ntypes
allocate(rVals(i)%arr(NSTATES,NSTATES))
allocate(kVals(i)%arr(NSTATES))
end do

allocate(super_xo((ntypes-10)*NSTATES*NSTATES,3))
allocate(super_xd((ntypes-10)*NSTATES*NSTATES,4))
allocate(superk_xo((ntypes-10)*NSTATES,3))
allocate(superk_xd((ntypes-10)*NSTATES,4))

! Correlations

do 10 j=1,6
do 10 iave=0,NAVE
cor_top(iave,j)=zero
10      continue
super_xo=0.0
super_xd=0.0
superk_xd=0.0
superk_xo=0.0

top=ZERO
hbin=ZERO    ! mean value of bin
btm=ZERO
hbinsq=ZERO  ! square of mean value
prd=ZERO

write(*,*)' IBLK      ENERGY       wgtTOT'

return
END subroutine InitRunMeas


!!!Initialize (Each) Block Averages: corptop/btm
SUBROUTINE InitBlkMeas
use cpmc

real(sp)::rescale


do iave=0,NAVE
corptop(iave)=zero
end do
nstpcor=0

estptop=ZERO
estpbtm=ZERO
nstp=0

return
END subroutine InitBlkMeas


!!!Measurements in Each Block
SUBROUTINE BlkMeas(iblk)
use cpmc
use jiekou,only:pair

integer::iblk
real(sp)::eblk,eblksq
real(sp)::gren_up(NSTATES,NSTATES),gren_dn(NSTATES,NSTATES),superk_tmp((ntypes-10)*NSTATES)
! real(sp)::swave(NSTATES),dwave(NSTATES),pwave(NSTATES),sowave(NSTATES)
! real(sp)::sbwave(NSTATES),dbwave(NSTATES),pbwave(NSTATES)
! real(sp)::puupxwave(NSTATES),pddpywave(NSTATES),pudpxwave(NSTATES)
! real(sp)::puupywave(NSTATES),pddpxwave(NSTATES),pudpywave(NSTATES)
! real(sp)::pdsfbd1(NSTATES),pdsfbd2(NSTATES),pdsfbd12(NSTATES)
! real(sp)::dd1wave(NSTATES),dd2wave(NSTATES),dd12wave(NSTATES)
! real(sp)::pbdx2y2d1(NSTATES),pbdx2y2d2(NSTATES),pbdx2y2d12(NSTATES)
save eblk,eblksq
! real(sp)::gx_up(NSTATES,NSTATES),gx_dn(NSTATES,NSTATES)
real(sp)::eavg

top   = top+(estptop)/real(nstp)
hbin = hbin+estptop/estpbtm
btm   = btm+(estpbtm)/real(nstp)
hbinsq = hbinsq+(estptop/estpbtm)**2
prd   = prd+(estptop)*(estpbtm)/real(nstp**2)

!-------------------------------------------!

do iave=0,NAVE
cor_top(iave,1) = cor_top(iave,1)+corptop(iave)/corptop(0)
cor_top(iave,2) = cor_top(iave,2)+(corptop(iave)/corptop(0))**2
end do

!~~~~~Covariance
do iave=1,NAVE
cor_top(iave,3) = cor_top(iave,3)&
& +( corptop(iave) )*( corptop(0) )/real(nstpcor**2)
end do

eavg=(estptop)/(estpbtm)
eblk=eblk+eavg
eblksq=eblksq+eavg*eavg

!*********************************************!
k=ntypes * NSTATES

do i=1,NSTATES
do j=1,NSTATES
k=k+1
gren_up(i,j)=corptop(k)/corptop(0)
end do
end do

do i=1,NSTATES
do j=1,NSTATES
k=k+1
gren_dn(i,j)=corptop(k)/corptop(0)
end do
end do

call pair(gx_up,gx_dn)

! real space
do i = 1, (ntypes-10)
super_xo((i-1)*NSTATES*NSTATES+1:i*NSTATES*NSTATES,1) = reshape(rVals(i+10)%arr, [NSTATES*NSTATES])
end do

m=ntypes*NSTATES+10*NSTATES*NSTATES
k=0
do i=1,(ntypes-10)*NSTATES*NSTATES
k=k+1; m=m+1
super_xd(k,1)=super_xd(k,1)+(corptop(m)/corptop(0)-super_xo(k,1))
super_xd(k,2)=super_xd(k,2)+(corptop(m)/corptop(0)-super_xo(k,1))**2
super_xd(k,4)=(corptop(m)/corptop(0)-super_xo(k,1))
end do

! moment space
do i = 1, (ntypes-10)
superk_xo((i-1)*NSTATES+1:i*NSTATES,1) = reshape(kVals(i+10)%arr, [NSTATES])
end do

m=10*NSTATES
k=0
do i=1,(ntypes-10)*NSTATES
k=k+1; m=m+1
superk_xd(k,1)=superk_xd(k,1)+(corptop(m)/corptop(0)-superk_xo(k,1))
superk_xd(k,2)=superk_xd(k,2)+(corptop(m)/corptop(0)-superk_xo(k,1))**2
superk_xd(k,4)=(corptop(m)/corptop(0)-superk_xo(k,1))
end do



!-------------------cao-r-k---------------------!
! superk_xd(k,2)=superk_xd(k,2)+cc*(super_xd(iu,4)**2)/real(nsites)                          !swave
! superk_xd(NSTATES+k,2)=superk_xd(NSTATES+k,2)+cc*(super_xd(NSTATES+iu,4)**2)/real(nsites)         !dwave
! superk_xd(2*NSTATES+k,2)=superk_xd(2*NSTATES+k,2)+cc*(super_xd(2*NSTATES+iu,4)**2)/real(nsites)   !pwave
! superk_xd(3*NSTATES+k,2)=superk_xd(3*NSTATES+k,2)+cc*(super_xd(3*NSTATES+iu,4)**2)/real(nsites)   !sowave

write(*,*) iblk,eavg,sum(wgtwlkr)

return
END subroutine BlkMeas


SUBROUTINE RunMeas
use cpmc

real(sp)::energy,stdtop,stdbtm,covar,arg,error



real(sp)::numberup,numberdn

top=top/real(nblk)
btm=btm/real(nblk)
prd=prd/real(nblk)
hbin=hbin/real(nblk)
hbinsq=hbinsq/real(nblk)

do j=1,3
do iave=0,NAVE
cor_top(iave,j) = cor_top(iave,j)/real(nblk)
end do
end do

energy=hbin

write(*,*) 'energy is the:',energy
write(*,*) 'energy Per N is the:',energy/lxy


do iave=1,NAVE
cor_final(iave) = cor_top(iave,1)
end do

do iave=1,NAVE         ! error
cor_top(iave,6)=sqrt(abs(cor_top(iave,2)-cor_top(iave,1)**2)&
& /real(nblk-1))
end do

do j=1,2
do i=1,(ntypes-10)*NSTATES*NSTATES
super_xo(i,j)=super_xo(i,j)/real(nblk)
super_xd(i,j)=super_xd(i,j)/real(nblk)
end do
end do

do j=1,2
do i=1,(ntypes-10)*NSTATES
superk_xo(i,j)=superk_xo(i,j)/real(nblk)
superk_xd(i,j)=superk_xd(i,j)/real(nblk)
end do
end do

do i=1,(ntypes-10)*NSTATES*NSTATES
super_xo(i,3)=sqrt(abs(super_xo(i,2)-super_xo(i,1)**2)/real(nblk-1))
super_xd(i,3)=sqrt(abs(super_xd(i,2)-super_xd(i,1)**2)/real(nblk-1))
end do

do i=1,(ntypes-10)*NSTATES
superk_xo(i,3)=sqrt(abs(superk_xo(i,2)-superk_xo(i,1)**2)/real(nblk-1))
superk_xd(i,3)=sqrt(abs(superk_xd(i,2)-superk_xd(i,1)**2)/real(nblk-1))
end do

numberup=0.0
numberdn=0.0
do j=1,NSITES
numberup=numberup+cor_final(j)
numberdn=numberdn+cor_final(nsites+j)
end do
write(*,*)
write(*,*)'Result:'
if(nblk <= 1)then
write(*,*) 'No Correlations'
write(*,'((a),f11.7)')'Energy=',energy
else 
call writeCorrDat
call writePairDat
call writeVertexDat
call writeUnpairDat
call cpOut
end if

RETURN
END subroutine RunMeas


subroutine writeCorrDat
use cpmc
open(unit=21,file='corr.dat',status='unknown')
write(21,*) '============= ONE BAND Hubbard Model========'
write(21,'((a20),i5,(a10),i5,(a10),i5,(a10),f5.3)') &
!&'lattice size=',nsites,'up=',NUP,'down=',NDN,'alpha=',alpha,'alphat1=',alphat1
&'lattice size=',nsites,'up=',NUP,'down=',NDN,'tam=',tam,'alphat1=',alphat1
write(21,'((a10),f10.5,(a10),f10.5,(a10),f10.5,(a10),f10.5)') &
&'Ud=',ud,'t=',t0,'t1=',t1,'vpd=',vpd
write(21,*) '--------------Parameters--------------------'
write(21,*) 'Time Slice is--------------------:',deltau
write(21,*) 'Variational energy---------------:',e_var
write(21,*) 'Block number of Equilibration----:',nblkeq
write(21,*) 'Block number of Growth Estimation:',nblkgr
write(21,*) 'Block number of Measurement------:',nblk
write(21,*) 'Steps of Each Block--------------:',nblkstps
write(21,*) 'Pop frequency of Equilibration---:',itvlpceq
write(21,*) 'Pop frequency of Measurement-----:',itvlpc
write(21,*) 'Orthogonalization frequency------:',itvlorth
write(21,*) 'STEPMEAS Energy steps------------:',itvlmeas
write(21,*) 'BACK PROPAGATION STEPS-----------:',itvl_m
write(21,*) '--------------------------------------------'

write(21,*) 'expection value of spin UP electrons=',numberup
write(21,*) 'expection value of spin DN electrons=',numberdn

error=sqrt(abs(hbinsq-hbin**2)/real(nblk-1))

write(21,'((a),f11.6,(a),f11.6)')&
& 'Mixed Energy=',energy,' with error',error

write(*,*) 'Correlation Functions in File corr.dat'

write(21,*) '********back propagation results:'
write(21,*) '                                 '
write(21,*) 'Charge__up:',&
&cor_final(NAVE-5),cor_top(NAVE-5,6)
write(21,*) 'Charge__dn:',&
&cor_final(NAVE-4),cor_top(NAVE-4,6)
write(21,*) 'magnetic moment:',&
&cor_final(NAVE-3),cor_top(NAVE-3,6)
write(21,*) 'kinetic energy:',&
&cor_final(NAVE-2),cor_top(NAVE-2,6)
write(21,*) 'potential energy:',&
&cor_final(NAVE-1),cor_top(NAVE-1,6)
write(21,*) 'total energy:',&
&cor_final(NAVE),cor_top(NAVE,6)

write(21,*) 'total energy per N:',&
&cor_final(NAVE)/lxy,cor_top(NAVE,6)/lxy

END subroutine writeCorrDat
!**************************************************!

subroutine writePairDat
use cpmc
integer :: i, j, k
real :: kx, ky

character(len=32), parameter :: kNames(32) = (/ &
    'gx_up','gx_dn','n_up','n_dn','n_w','cdw','sdwz','sdwx',&
    'pmdf','dsf','swave','pwave','dwave','sowave','sbwave','pbwave',&
    'dbwave','pdsfbd1','pdsfbd2','pdsfbd12','dd1wave','dd2wave','dd12wave',&
    'puupxwave','pddpxwave','pudpxwave','puupywave','pddpywave','pudpywave',&
    'pbdx2y2d1','pbdx2y2d2','pbdx2y2d12' /)

call system('mkdir -p ' // trim('dir-kVals/'))
do j = 1, ntypes
open(unit=100+j, file='dir-kVals/'//trim(kNames(j))//'.dat', status='unknown')
write(100+j,'(2A,2A)') 'k_x:','k_y:', trim(kNames(j))//':', trim(kNames(j))//'_error:'
do i = 1, NSTATES
k = (j-1)*NSTATES + i
ky =kSet(i,2)
kx =kSet(i,1)
if (kSet(i,1) > pi) then
    kx = kSet(i,1) - 2 * pi
end if
if (kSet(i,2) > pi) then
    ky = kSet(i,2) - 2 * pi
end if
write(100+j,'(2F15.9,2F15.9)') kx, ky, cor_final(k), cor_top(k,6)
end do
close(100+j)
end do

call system('mkdir -p ' // trim('dir-rVals/'))
do j = 1, ntypes
open(unit=200+j, file='dir-rVals/'//trim(kNames(j))//'.dat', status='unknown')
write(200+j,'(2A,2A)') 'r_x:','r_y:', trim(kNames(j))//':', trim(kNames(j))//'_error:'
do m = 1, NSTATES
do n = 1, NSTATES
k = (j-1)*NSTATES*NSTATES + (m-1)*NSTATES + n
write(200+j,'(2I8,1F20.10)') m, n, cor_final(k), cor_top(k,6)
end do
end do
close(200+j)
end do
END subroutine writePairDat

subroutine writeVertexDat
use cpmc
integer :: i, j, k
real :: kx, ky

character(len=32), parameter :: kNames(22) = (/ &
    'Vertex_swave','Vertex_pwave','Vertex_dwave','Vertex_sowave','Vertex_sbwave','Vertex_pbwave',&
    'Vertex_dbwave','Vertex_pdsfbd1','Vertex_pdsfbd2','Vertex_pdsfbd12','Vertex_dd1wave','Vertex_dd2wave','Vertex_dd12wave',&
    'Vertex_puupxwave','Vertex_pddpxwave','Vertex_pudpxwave','Vertex_puupywave','Vertex_pddpywave','Vertex_pudpywave',&
    'Vertex_pbdx2y2d1','Vertex_pbdx2y2d2','Vertex_pbdx2y2d12' /)
k=0
call system('mkdir -p ' // trim('dir-kVals/'))
do j = 1, (ntypes-10)
open(unit=300+j, file='dir-kVals/'//trim(kNames(j))//'.dat', status='unknown')
write(300+j,'(2A,2A)') 'k_x:','k_y:', trim(kNames(j))//':', trim(kNames(j))//'_error:'
do i = 1, NSTATES
k=k+1
ky =kSet(i,2)
kx =kSet(i,1)
if (kSet(i,1) > pi) then
    kx = kSet(i,1) - 2 * pi
end if
if (kSet(i,2) > pi) then
    ky = kSet(i,2) - 2 * pi
end if
write(300+j,'(2F15.9,2F15.9)') kx, ky, superk_xd(k,1),superk_xd(k,3)
end do
close(300+j)
end do

k=0
call system('mkdir -p ' // trim('dir-rVals/'))
do j = 1, (ntypes-10)
open(unit=400+j, file='dir-rVals/'//trim(kNames(j))//'.dat', status='unknown')
write(400+j,'(2A,2A)') 'r_x:','r_y:', trim(kNames(j))//':', trim(kNames(j))//'_error:'
do m = 1, NSTATES
do n = 1, NSTATES
k=k+1
write(400+j,'(2I8,1F20.10)') m, n, super_xd(k,1),super_xd(k,3)
end do
end do
close(400+j)
end do
END subroutine writeVertexDat

subroutine writeUnpairDat
use cpmc
integer :: i, j, k
real :: kx, ky 

character(len=32), parameter :: kNames(22) = (/ &
    'Unpair_swave','Unpair_pwave','Unpair_dwave','Unpair_sowave','Unpair_sbwave','Unpair_pbwave',&
    'Unpair_dbwave','Unpair_pdsfbd1','Unpair_pdsfbd2','Unpair_pdsfbd12','Unpair_dd1wave','Unpair_dd2wave','Unpair_dd12wave',&
    'Unpair_puupxwave','Unpair_pddpxwave','Unpair_pudpxwave','Unpair_puupywave','Unpair_pddpywave','Unpair_pudpywave',&
    'Unpair_pbdx2y2d1','Unpair_pbdx2y2d2','Unpair_pbdx2y2d12' /)
k=0
call system('mkdir -p ' // trim('dir-kVals/'))
do j = 1, (ntypes-10)
open(unit=300+j, file='dir-kVals/'//trim(kNames(j))//'.dat', status='unknown')
write(300+j,'(2A,2A)') 'k_x:','k_y:', trim(kNames(j))//':', trim(kNames(j))//'_error:'
do i = 1, NSTATES
k=k+1
ky =kSet(i,2)
kx =kSet(i,1)
if (kSet(i,1) > pi) then
    kx = kSet(i,1) - 2 * pi
end if
if (kSet(i,2) > pi) then
    ky = kSet(i,2) - 2 * pi
end if
write(300+j,'(2F15.9,2F15.9)') kx, ky, superk_xo(k,1),superk_xo(k,3)
end do
close(300+j)
end do

k=0
call system('mkdir -p ' // trim('dir-rVals/'))
do j = 1, (ntypes-10)
open(unit=400+j, file='dir-rVals/'//trim(kNames(j))//'.dat', status='unknown')
write(400+j,'(2A,2A)') 'r_x:','r_y:', trim(kNames(j))//':', trim(kNames(j))//'_error:'
do m = 1, NSTATES
do n = 1, NSTATES
k=k+1
write(400+i,'(2I8,1F20.10)') m, n, super_xo(k,1),super_xo(k,3)
end do
end do
close(400+j)
end do
END subroutine writeUnpairDat

subroutine FourierTransform(rval,kvalue)
use cpmc
real(sp), intent(in) :: rval(NSTATES,NSTATES)
real(sp), intent(inout) :: kvalue(NSTATES)

kvalue = 0.0
k = 0
do i = 1, lx
do j = 1, ly
k = k + 1
do m = 1, NSTATES
do n = 1, NSTATES
phase = xk(i)*(ixv(m) - ixv(n)) + yk(j)*(iyv(m) - iyv(n))
cc = cos(phase)
kvalue(k) = kvalue(k) + cc * rval(n,m)
end do
end do
kvalue(k) = kvalue(k)/real(NSTATES*NSTATES, sp)
end do
end do

end subroutine FourierTransform

subroutine correl
use cpmc
use jiekou,only:sgefa,sgedi,pair,backphi
real(sp)::phip_up(NUP,NSTATES),phip_dn(NDN,NSTATES)
integer::i,j,k,l,m,n,iu,i1,i2,dx,dy
integer::m1,m2,m3,m4,n1,n2,n3,n4
integer::info,ipvt(NSTATES)
real(sp)::det(2),work(NSTATES)

! type(wave_t), dimension(ntypes) :: kVals
kVals(1)%arr	=>	gxk_up
kVals(2)%arr	=>	gxk_dn
kVals(3)%arr	=>	nk_up
kVals(4)%arr	=>	nk_dn
kVals(5)%arr	=>	nk_w
kVals(6)%arr	=>	cdw
kVals(7)%arr	=>	sdwz
kVals(8)%arr	=>	sdwx
kVals(9)%arr	=>	pmdf
kVals(10)%arr	=>	dsf
kVals(11)%arr	=>	swave
kVals(12)%arr	=>	pwave
kVals(13)%arr	=>	dwave
kVals(14)%arr	=>	sowave
kVals(15)%arr	=>	sbwave
kVals(16)%arr	=>	pbwave
kVals(17)%arr	=>	dbwave
kVals(18)%arr	=>	pdsfbd1
kVals(19)%arr	=>	pdsfbd2
kVals(20)%arr	=>	pdsfbd12
kVals(21)%arr	=>	dd1wave
kVals(22)%arr	=>	dd2wave
kVals(23)%arr	=>	dd12wave
kVals(24)%arr	=>	puupxwave
kVals(25)%arr	=>	pddpxwave
kVals(26)%arr	=>	pudpxwave
kVals(27)%arr	=>	puupywave
kVals(28)%arr	=>	pddpywave
kVals(29)%arr	=>	pudpywave
kVals(30)%arr	=>	pbdx2y2d1
kVals(31)%arr	=>	pbdx2y2d2
kVals(32)%arr	=>	pbdx2y2d12


! type(wave_t), dimension(ntypes) :: rVals
rVals(1)%arr	=>	gx_up
rVals(2)%arr	=>	gx_dn
rVals(3)%arr	=>	n_up
rVals(4)%arr	=>	n_dn
rVals(5)%arr	=>	n_w
rVals(6)%arr	=>	cdw_re
rVals(7)%arr	=>	sdwz_re
rVals(8)%arr	=>	sdwx_re
rVals(9)%arr	=>	pmdf_re
rVals(10)%arr	=>	dsf_re
rVals(11)%arr	=>	swave_re
rVals(12)%arr	=>	pwave_re
rVals(13)%arr	=>	dwave_re
rVals(14)%arr	=>	sowave_re
rVals(15)%arr	=>	sbwave_re
rVals(16)%arr	=>	pbwave_re
rVals(17)%arr	=>	dbwave_re
rVals(18)%arr	=>	pdsfbd1_re
rVals(19)%arr	=>	pdsfbd2_re
rVals(20)%arr	=>	pdsfbd12_re
rVals(21)%arr	=>	dd1wave_re
rVals(22)%arr	=>	dd2wave_re
rVals(23)%arr	=>	dd12wave_re
rVals(24)%arr	=>	puupxwave_re
rVals(25)%arr	=>	pddpxwave_re
rVals(26)%arr	=>	pudpxwave_re
rVals(27)%arr	=>	puupywave_re
rVals(28)%arr	=>	pddpywave_re
rVals(29)%arr	=>	pudpywave_re
rVals(30)%arr	=>	pbdx2y2d1_re
rVals(31)%arr	=>	pbdx2y2d2_re
rVals(32)%arr	=>	pbdx2y2d12_re
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
! Get < Phi(k',l) |; up & down; phip_up( NUP x NSTATES ) Transposed already!
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

call dgemm('N', 'N', NUP, NUP, NSTATES, 1.0_8, phip_up, NUP, phi_cup(:,:,iw), NSTATES, 0.0_8, ovlpINV_up, NUP) !li
 
!ovlpINV_up=matmul(phip_up,phi_cup(:,:,iw)) 

!-----------------------------------------------------!
! compute inverse and determinant of overlap matrix LR!
!-----------------------------------------------------!

call sgefa(ovlpINV_up,NUP,NUP,ipvt,info)
if(info /= 0) stop 'Problem in sgefa routine'
call sgedi(ovlpINV_up,NUP,NUP,ipvt,det,work,11_i4b)
detp_up(ibas)=det(1)*TEN**det(2)

!--------------------------!
! compute g (R[(LR)^{-1}]L)!
!--------------------------!
call dgemm('N', 'N', NSTATES, NUP, NUP, 1.0_8, phi_cup(:,:,iw), NSTATES, ovlpINV_up, NUP, 0.0_8, tmp_up, NSTATES) !li

!tmp_up=matmul(phi_cup(:,:,iw),ovlpINV_up)

call dgemm('N', 'N', NSTATES, NSTATES, NUP, 1.0_8, tmp_up, NSTATES, phip_up, NUP, 0.0_8, gx_up, NSTATES) !li

!gx_up=matmul(tmp_up,phip_up)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! For Spin Down                                              !
! L = < Phi(k',l) |; R = | Phi(k) > at step N(saved); k == iw!
! Compute overlap matrix LR, etc.                            !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

call dgemm('N', 'N', NDN, NDN, NSTATES, 1.0_8, phip_dn, NDN, phi_cdn(:,:,iw), NSTATES, 0.0_8, ovlpINV_dn, NDN) !li

!ovlpINV_dn=matmul(phip_dn,phi_cdn(:,:,iw))

!-----------------------------------------------------!
! compute inverse and determinant of overlap matrix LR!
!-----------------------------------------------------!

call sgefa(ovlpINV_dn,NDN,NDN,ipvt,info)
if(info /= 0) stop 'Problem in sgefa routine'
call sgedi(ovlpINV_dn,NDN,NDN,ipvt,det,work,11_i4b)
detp_dn(ibas)=det(1)*TEN**det(2)

!--------------------------!
! compute g (R[(LR)^{-1}]L)!
!--------------------------!

call dgemm('N', 'N', NSTATES, NDN, NDN, 1.0_8, phi_cdn(:,:,iw), NSTATES, ovlpINV_dn, NDN, 0.0_8, tmp_dn, NSTATES) !li

!tmp_dn=matmul(phi_cdn(:,:,iw),ovlpINV_dn)

call dgemm('N', 'N', NSTATES, NSTATES, NDN, 1.0_8, tmp_dn, NSTATES, phip_dn, NDN, 0.0_8, gx_dn, NSTATES) !li

!gx_dn=matmul(tmp_dn,phip_dn)

!------------------------------------------------!
! Here need change for more than one determinants!
!------------------------------------------------!

ave(ibas,0)=1.0

!--------------------------------------------------------!
! Averages, Wick's theorem, etc.  Take the following form!
!--------------------------------------------------------!
! < c_i^\dagger c_j > = g(j,i) * <Psi(l) | Phi> (overlap)!
!                    -> g(j,i)                           !
!--------------------------------------------------------!
 
do iave=1,NAVE
ave_cor(iave)=zero
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

gxk_up=0.0
gxk_dn=0.0

call FourierTransform(gx_up,gxk_up)
call FourierTransform(gx_dn,gxk_dn)

do i=1,NSTATES            ! diagonal
opn_up(i)=gx_up(i,i)
opn_dn(i)=gx_dn(i,i)
enddo

do i=1,NSTATES
do j=1,NSTATES
up_up=gx_up(i,i)*gx_up(j,j)-gx_up(j,i)*gx_up(i,j)
dn_dn=gx_dn(i,i)*gx_dn(j,j)-gx_dn(j,i)*gx_dn(i,j)
up_dn=gx_up(i,i)*gx_dn(j,j)+gx_up(j,j)*gx_dn(i,i) ! 2

n_up(i,j)=gx_up(j,i)+gx_up(i,j)
n_dn(i,j)=gx_dn(j,i)+gx_dn(i,j)
n_w(i,j)=gx_up(j,i)+gx_dn(j,i)+gx_up(i,j)+gx_dn(i,j)

cdw_re(i,j)=up_up+dn_dn+up_dn+dens*(dens-opn_up(i)-opn_dn(i)-&
opn_up(j)-opn_dn(j))        ! CDW
sdwz_re(i,j)=up_up+dn_dn-up_dn       ! SDW-z
sdwx_re(i,j)=-gx_up(j,i)*gx_dn(i,j)-&
& gx_up(i,j)*gx_dn(j,i)   ! SDW-x

if (i==j) then
pmdf_re(i,j)=pmdf_re(i,j)+gx_up(j,i)*gx_dn(j,i)
dsf_re(i,j)=dsf_re(i,j)+gx_up(j,i)*gx_dn(j,i)
else
pmdf_re(i,j)=pmdf_re(i,j)+gx_up(j,i)*gx_dn(j,i)
dsf_re(i,j)=dsf_re(i,j)+(gx_up(i,i)*gx_up(j,j)+gx_up(j,i)*gx_up(i,j))*&
&(gx_dn(i,i)*gx_dn(j,j)+gx_dn(j,i)*gx_dn(i,j))
endif
enddo
enddo 

cu_up=sum(opn_up(1:nsites))/real(nsites)
cu_dn=sum(opn_dn(1:nsites))/real(nsites)

z2_sum=0.0
do i=1,nsites
z2_sum=z2_sum+opn_up(i)+opn_dn(i)-2.0*gx_up(i,i)*gx_dn(i,i)
end do
z2_sum=z2_sum/real(nsites)

do j=1,NSTATES
do i=1,NSTATES
hkin=hkin+tk(j,i,1)*gx_up(i,j)+tk(j,i,2)*gx_dn(i,j)
end do
end do
do i=1,NSTATES
hpot=hpot+hub_u(i)*gx_up(i,i)*gx_dn(i,i)
end do

if(abs(vpd)>0.001) then
do i=1,nsites
do j1=1,2
j=idis(i,j1)
hpot=hpot+vpd*(gx_up(i,i)*gx_up(j,j)-gx_up(i,j)*gx_up(j,i))&
& +vpd*(gx_dn(i,i)*gx_dn(j,j)-gx_dn(i,j)*gx_dn(j,i))&
& +vpd*(gx_up(i,i)*gx_dn(j,j)+gx_up(j,j)*gx_dn(i,i))
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

call FourierTransform(n_up,nk_up)
call FourierTransform(n_dn,nk_dn)
call FourierTransform(n_w,nk_w)
call FourierTransform(cdw_re*2.0,cdw)
call FourierTransform(sdwz_re*2.0,sdwz)
call FourierTransform(sdwx_re*2.0,sdwx)

!     dsf(k)=<(\rho_b(r)-\rho_avg)  (\rho(r')>-\rho_avg)>
!            =<\rho_b(r) \rho_b(r')> -\bar{\rho} (<\rho_b(r')> + <\rho_b(r)>) + \bar{\rho}^2
!       <\rho_b(r')>= double occupancy at r'= cdw(r')
call FourierTransform(pmdf_re,pmdf)
call FourierTransform(dsf_re,dsf)

sum_up = sum([(gx_up(i,i), i=1, nsites)])
sum_dn = sum([(gx_dn(i,i), i=1, nsites)])
sum_cdw = sum([(gx_up(i,i)*gx_dn(i,i), i=1, nsites)])
sum_op = sum( opn_up + opn_dn )              

term1=sum_op-2.0*sum_cdw
term2=(sum_op+2.0*sum_cdw)+dens*(real(nsites,sp)*dens-2.0*sum_op)

nk_up = nk_up + sum_up
nk_dn = nk_dn + sum_dn
nk_w  = nk_w  + sum_up + sum_dn
sdwz = sdwz + term1
sdwx = sdwx + term1
cdw = cdw + term2


!***************************************************!
!  compute pairing correlation functions!
!---------------------------------------------------!

call pair(gx_up,gx_dn)

do i = 1, ntypes
kValsArr((i-1)*NSTATES+1:i*NSTATES) = kVals(i)%arr(1:NSTATES)
end do
ave_cor(1:ntypes*NSTATES) = kValsArr(1:ntypes*NSTATES)

do i = 1, ntypes
rValsArr((i-1)*NSTATES*NSTATES+1:i*NSTATES*NSTATES) = reshape(rVals(i)%arr(1:NSTATES,1:NSTATES), [NSTATES*NSTATES])
end do
ave_cor(1:ntypes*NSTATES*NSTATES) = ave_cor(1:ntypes*NSTATES*NSTATES) + rValsArr(1:ntypes*NSTATES*NSTATES)

ave_cor(NAVE-5)=cu_up
ave_cor(NAVE-4)=cu_dn
ave_cor(NAVE-3)=z2_sum
ave_cor(NAVE-2)=hkin
ave_cor(NAVE-1)=hpot
ave_cor(NAVE)=have

do iave=1,NAVE
ave(ibas,iave)=ave_cor(iave)*ave(ibas,0)
end do

1000    continue                             ! loop over ibas


do iave=0,nave
wgt_p(iave)=zero
do ibas=1,NWFBAS
wgt_p(iave)=wgt_p(iave)+ave(ibas,iave)
end do                            ! loop over ibas & kp
end do

!*****************************************************************!
!------------------------------------------------!
! Here need change for more than one determinants!
!------------------------------------------------!

do iave=0,nave
corwlkr(iave)=wgt_c(iw)*wgt_p(iave)
end do

do iave=0,nave                         ! 0 for bottom
corptop(iave)=corptop(iave)+corwlkr(iave)
end do

end if

1100    continue                                ! loop over k (==iw)

nstpcor=nstpcor+1 

end subroutine correl


subroutine pair(gp_up, gp_dn)
use cpmc

integer :: i1, i2, k, iu, m, n, m1, n1, mm, nn, ii, jj, i3, j3, m3, n3
real(sp), intent(in) :: gp_up(:,:), gp_dn(:,:)
! real(sp), target :: p_sdp, puu_sdp, pdd_sdp, pdsfbd, p_ddp, pdsfb, pbdx2y2
! real(sp), target :: swave_re(NSTATES,NSTATES), dwave_re(NSTATES,NSTATES), pwave_re(NSTATES,NSTATES), sowave_re(NSTATES,NSTATES)
! real(sp), target :: sbwave_re(NSTATES,NSTATES), dbwave_re(NSTATES,NSTATES), pbwave_re(NSTATES,NSTATES)
! real(sp), target :: puupxwave_re(NSTATES,NSTATES), pddpywave_re(NSTATES,NSTATES), pudpxwave_re(NSTATES,NSTATES), puupywave_re(NSTATES,NSTATES)
! real(sp), target :: pddpxwave_re(NSTATES,NSTATES), pudpywave_re(NSTATES,NSTATES)
! real(sp), target :: pdsfbd1_re(NSTATES,NSTATES), pdsfbd2_re(NSTATES,NSTATES), pdsfbd12_re(NSTATES,NSTATES)
! real(sp), target :: dd1wave_re(NSTATES,NSTATES), dd2wave_re(NSTATES,NSTATES), dd12wave_re(NSTATES,NSTATES)
! real(sp), target :: pbdx2y2d1_re(NSTATES,NSTATES), pbdx2y2d2_re(NSTATES,NSTATES), pbdx2y2d12_re(NSTATES,NSTATES)

! real(sp), target :: swave(NSTATES), dwave(NSTATES), pwave(NSTATES), sowave(NSTATES)
! real(sp), target :: sbwave(NSTATES), dbwave(NSTATES), pbwave(NSTATES)
! real(sp), target :: puupxwave(NSTATES), pddpywave(NSTATES), puupywave(NSTATES), pddpxwave(NSTATES)
! real(sp), target :: pudpxwave(NSTATES), pudpywave(NSTATES)
! real(sp), target :: pdsfbd1(NSTATES), pdsfbd2(NSTATES), pdsfbd12(NSTATES)
! real(sp), target :: dd1wave(NSTATES), dd2wave(NSTATES), dd12wave(NSTATES)
! real(sp), target :: pbdx2y2d1(NSTATES), pbdx2y2d2(NSTATES), pbdx2y2d12(NSTATES)

swave_re       = 0.0
dwave_re       = 0.0
pwave_re       = 0.0
sowave_re      = 0.0
sbwave_re      = 0.0
dbwave_re      = 0.0
pbwave_re      = 0.0
puupxwave_re   = 0.0
pddpywave_re   = 0.0
pudpxwave_re   = 0.0
puupywave_re   = 0.0
pddpxwave_re   = 0.0
pudpywave_re   = 0.0
pdsfbd1_re     = 0.0
pdsfbd2_re     = 0.0
pdsfbd12_re    = 0.0
dd1wave_re     = 0.0
dd2wave_re     = 0.0
dd12wave_re    = 0.0
pbdx2y2d1_re   = 0.0
pbdx2y2d2_re   = 0.0
pbdx2y2d12_re  = 0.0

do m=1,nsites
do n=1,nsites

do i1=1,4
m1=iposit(ixv(m)+ipx(i1),iyv(m)+ipy(i1))
do i2=1,4
n1=iposit(ixv(n)+ipx(i2),iyv(n)+ipy(i2))

p_sdp=gp_up(m,n)*gp_dn(m1,n1)+gp_up(m1,n1)*gp_dn(m,n)
swave_re(m,n)=sf(i1)*sf(i2)*p_sdp
dwave_re(m,n)=df(i1)*df(i2)*p_sdp
pwave_re(m,n)=pf(i1)*pf(i2)*p_sdp 
puu_sdp=2*(gp_up(m,n)*gp_up(m1,n1)-gp_up(m,n1)*gp_up(m1,n))
pdd_sdp=2*(gp_dn(m,n)*gp_dn(m1,n1)-gp_dn(m,n1)*gp_dn(m1,n))
puupxwave_re(m,n)=pf(i1)*pf(i2)*puu_sdp
pddpywave_re(m,n)=ppf(i1)*ppf(i2)*pdd_sdp
puupywave_re(m,n)=ppf(i1)*ppf(i2)*puu_sdp
pddpxwave_re(m,n)=pf(i1)*pf(i2)*pdd_sdp  
pdsfb=2*(gp_up(n,m)*gp_up(n1,m1)-gp_up(n1,m)*gp_up(n,m1))*(gp_dn(n,m)*gp_dn(n1,m1)-gp_dn(n1,m)*gp_dn(n,m1))
sbwave_re(m,n)=sf(i1)*sf(i2)*pdsfb
dbwave_re(m,n)=df(i1)*df(i2)*pdsfb
pbwave_re(m,n)=pf(i1)*pf(i2)*pdsfb
end do
end do
sowave_re(m,n)=gp_up(m,n)*gp_dn(m,n)

!----------------------dsfb---------------!
do ii=1,4
do jj=1,4
mm=iposit(ixv(m)+ippx(ii),iyv(m)+ippy(ii)) 
nn=iposit(ixv(n)+ippx(jj),iyv(n)+ippy(jj))  
pdsfbd=2*(gp_up(n,m)*gp_up(nn,mm)-gp_up(nn,m)*gp_up(n,mm))*(gp_dn(n,m)*gp_dn(nn,mm)-gp_dn(nn,m)*gp_dn(n,mm)) 
if(ddf(ii)*ddf(jj)==-1) then 
pdsfbd1_re(m,n)=ddf(ii)*ddf(jj)*pdsfbd
else
pdsfbd2_re(m,n)=ddf(ii)*ddf(jj)*pdsfbd
end if

p_ddp=gp_up(m,n)*gp_dn(mm,nn)+gp_up(mm,nn)*gp_dn(m,n) 
if(ddf(ii)*ddf(jj)==-1) then 
dd1wave_re(m,n)=ddf(ii)*ddf(jj)*p_ddp
else
dd2wave_re(m,n)=ddf(ii)*ddf(jj)*p_ddp
end if
end do
end do 

!----------------------------dx2-y2-------!
do i3=1,4
do j3=1,4
m3=iposit(ixv(m)+idx(i3),iyv(m)+idy(i3)) 
n3=iposit(ixv(n)+idx(j3),iyv(n)+idy(j3))  
pbdx2y2=2*(gp_up(n,m)*gp_up(n3,m3)-gp_up(n3,m)*gp_up(n,m3))*(gp_dn(n,m)*gp_dn(n3,m3)-gp_dn(n3,m)*gp_dn(n,m3))
if(ddf(i3)*ddf(j3)==-1) then 
pbdx2y2d1_re(m,n)=+ddf(i3)*ddf(j3)*pbdx2y2
else
pbdx2y2d2_re(m,n)=ddf(i3)*ddf(j3)*pbdx2y2
end if
end do
end do 
end do
end do


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

do iu=1,NSTATES

swave(iu)=swave(iu)/2
dwave(iu)=dwave(iu)/2
pwave(iu)=pwave(iu)/1.41421356
sowave(iu)=sowave(iu)
sbwave(iu)=sbwave(iu)/2
dbwave(iu)=dbwave(iu)/2
pbwave(iu)=pbwave(iu)/1.41421356
puupxwave(iu)=puupxwave(iu)/1.41421356
pddpywave(iu)=pddpywave(iu)/1.41421356
puupywave(iu)=puupywave(iu)/1.41421356
pddpxwave(iu)=pddpxwave(iu)/1.41421356
pdsfbd1(iu)=pdsfbd1(iu)/2
pdsfbd2(iu)=pdsfbd2(iu)/2
dd1wave(iu)=dd1wave(iu)/2
dd2wave(iu)=dd2wave(iu)/2
pbdx2y2d1(iu)=pbdx2y2d1(iu)/2
pbdx2y2d2(iu)=pbdx2y2d2(iu)/2
pudpxwave(iu)=puupxwave(iu)+pddpxwave(iu) 
pudpywave(iu)=puupywave(iu)+pddpywave(iu) 
pdsfbd12(iu)=pdsfbd1(iu)+pdsfbd2(iu)
dd12wave(iu)=dd1wave(iu)+dd2wave(iu)
pbdx2y2d12(iu)=pbdx2y2d1(iu)+pbdx2y2d2(iu) 
end do

return
end subroutine pair


 
