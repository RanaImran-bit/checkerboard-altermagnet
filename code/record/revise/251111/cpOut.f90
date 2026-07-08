!**************************************************!

!~~~~~Initialize Final Averages: cor_top/btm
!     1: averages; 2: averages square; 3: top-bottom cross term;
!~~~~~(**4: standard deviation; 5: covariance;**) 6: error.

SUBROUTINE InitRunMeas
use cpmc
integer j,iave
! Correlations
do j=1,6
    do iave=0,NAVE
        cor_top(iave,j,:,:)=zero
    end do
end do

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
integer iave
real(sp)::rescale
do iave=0,NAVE
    corptop(iave,:,:)=zero
end do
nstpcor=0.0
estptop=ZERO
estpbtm=ZERO
eblk = 0.0
eblksq = 0.0
nstp=0.0
return
END subroutine InitBlkMeas


!!!Measurements in Each Block
SUBROUTINE BlkMeas(iblk)
use cpmc
use jiekou,only:pair
integer::iblk,k,i,j,m,n,iave,i_a,i_b
real(sp)::gren_up(nsites,nsites),gren_dn(nsites,nsites)
! real(sp)::superk_tmp((ntypes-10)*NSTATES)
! save eblk,eblksq
real(sp)::eavg

top   = top+(estptop)/real(nstp)
hbin = hbin+estptop/estpbtm
btm   = btm+(estpbtm)/real(nstp)
hbinsq = hbinsq+(estptop/estpbtm)**2
prd   = prd+(estptop)*(estpbtm)/real(nstp**2)

do i_a=1,NLA
    do i_b=1,NLA
        do iave=0,NAVE
            cor_top(iave,1,i_a,i_b) = cor_top(iave,1,i_a,i_b)&
            &+corptop(iave,i_a,i_b)/corptop(0,i_a,i_b)
            cor_top(iave,2,i_a,i_b) = cor_top(iave,2,i_a,i_b)&
            &+(corptop(iave,i_a,i_b)/corptop(0,i_a,i_b))**2
        end do

        !~~~~~Covariance
        do iave=1,NAVE
            cor_top(iave,3,i_a,i_b) = cor_top(iave,3,i_a,i_b)&
            & +( corptop(iave,i_a,i_b) )*( corptop(0,i_a,i_b) )/real(nstpcor**2)
        end do
    end do
end do

eavg=(estptop)/(estpbtm)
eblk=eblk+eavg
eblksq=eblksq+eavg*eavg

!*********************************************!

do i_a=1,NLA
    do i_b=1,NLA
        k=NSTATES*ntypes
        do i=1,NSTATES
            do j=1,NSTATES
                k=k+1
                gren_up(i+(i_a-1)*NSTATES,j+(i_b-1)*NSTATES)=corptop(k,i_a,i_b)/corptop(0,i_a,i_b)
            end do
        end do
    end do
end do

do i_a=1,NLA
    do i_b=1,NLA
        k=NSTATES*ntypes+NSTATES*NSTATES
        do i=1,NSTATES
            do j=1,NSTATES
                k=k+1
                gren_dn(i+(i_a-1)*NSTATES,j+(i_b-1)*NSTATES)=corptop(k,i_a,i_b)/corptop(0,i_a,i_b)
            end do
        end do
    end do
end do

call pair(gren_up,gren_dn)

do i_a=1,NLA
    do i_b=1,NLA
        ! real space
        k=0
        do m=1,nsites
            do n=1,nsites
                if (sublatt(m)==i_a) then
                    if (sublatt(n)==i_b) then
                        k=k+1
                        super_xo(k+NSTATES*NSTATES*0,1,sublatt(m),sublatt(n))=swave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*1,1,sublatt(m),sublatt(n))=pwave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*2,1,sublatt(m),sublatt(n))=dwave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*3,1,sublatt(m),sublatt(n))=sowave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*4,1,sublatt(m),sublatt(n))=sbwave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*5,1,sublatt(m),sublatt(n))=pbwave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*6,1,sublatt(m),sublatt(n))=dbwave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*7,1,sublatt(m),sublatt(n))=pdsfbd1_re(m,n)
                        super_xo(k+NSTATES*NSTATES*8,1,sublatt(m),sublatt(n))=pdsfbd2_re(m,n)
                        super_xo(k+NSTATES*NSTATES*9,1,sublatt(m),sublatt(n))=pdsfbd12_re(m,n)
                        super_xo(k+NSTATES*NSTATES*10,1,sublatt(m),sublatt(n))=dd1wave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*11,1,sublatt(m),sublatt(n))=dd2wave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*12,1,sublatt(m),sublatt(n))=dd12wave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*13,1,sublatt(m),sublatt(n))=puupxwave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*14,1,sublatt(m),sublatt(n))=pddpxwave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*15,1,sublatt(m),sublatt(n))=pudpxwave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*16,1,sublatt(m),sublatt(n))=puupywave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*17,1,sublatt(m),sublatt(n))=pddpywave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*18,1,sublatt(m),sublatt(n))=pudpywave_re(m,n)
                        super_xo(k+NSTATES*NSTATES*19,1,sublatt(m),sublatt(n))=pbdx2y2d1_re(m,n)
                        super_xo(k+NSTATES*NSTATES*20,1,sublatt(m),sublatt(n))=pbdx2y2d2_re(m,n)
                        super_xo(k+NSTATES*NSTATES*21,1,sublatt(m),sublatt(n))=pbdx2y2d12_re(m,n)
                    end if
                end if
            end do
        end do

        m=ntypes*NSTATES + 8*NSTATES*NSTATES

        k=0
        do i=1,(ntypes-8)*NSTATES*NSTATES
            k=k+1; m=m+1
            super_xd(k,1,i_a,i_b)=super_xd(k,1,i_a,i_b)+(corptop(m,i_a,i_b)/corptop(0,i_a,i_b)-super_xo(k,1,i_a,i_b))
            super_xd(k,2,i_a,i_b)=super_xd(k,2,i_a,i_b)+(corptop(m,i_a,i_b)/corptop(0,i_a,i_b)-super_xo(k,1,i_a,i_b))**2
            super_xd(k,4,i_a,i_b)=(corptop(m,i_a,i_b)/corptop(0,i_a,i_b)-super_xo(k,1,i_a,i_b))
        end do

        ! moment space
        k=0
        do m=1,NSTATES
            k=k+1
            superk_xo(k+NSTATES*0,1,i_a,i_b)=swave(m,i_a,i_b)
            superk_xo(k+NSTATES*1,1,i_a,i_b)=pwave(m,i_a,i_b)
            superk_xo(k+NSTATES*2,1,i_a,i_b)=dwave(m,i_a,i_b)
            superk_xo(k+NSTATES*3,1,i_a,i_b)=sowave(m,i_a,i_b)
            superk_xo(k+NSTATES*4,1,i_a,i_b)=sbwave(m,i_a,i_b)
            superk_xo(k+NSTATES*5,1,i_a,i_b)=pbwave(m,i_a,i_b)
            superk_xo(k+NSTATES*6,1,i_a,i_b)=dbwave(m,i_a,i_b)
            superk_xo(k+NSTATES*7,1,i_a,i_b)=pdsfbd1(m,i_a,i_b)
            superk_xo(k+NSTATES*8,1,i_a,i_b)=pdsfbd2(m,i_a,i_b)
            superk_xo(k+NSTATES*9,1,i_a,i_b)=pdsfbd12(m,i_a,i_b)
            superk_xo(k+NSTATES*10,1,i_a,i_b)=dd1wave(m,i_a,i_b)
            superk_xo(k+NSTATES*11,1,i_a,i_b)=dd2wave(m,i_a,i_b)
            superk_xo(k+NSTATES*12,1,i_a,i_b)=dd12wave(m,i_a,i_b)
            superk_xo(k+NSTATES*13,1,i_a,i_b)=puupxwave(m,i_a,i_b)
            superk_xo(k+NSTATES*14,1,i_a,i_b)=pddpxwave(m,i_a,i_b)
            superk_xo(k+NSTATES*15,1,i_a,i_b)=pudpxwave(m,i_a,i_b)
            superk_xo(k+NSTATES*16,1,i_a,i_b)=puupywave(m,i_a,i_b)
            superk_xo(k+NSTATES*17,1,i_a,i_b)=pddpywave(m,i_a,i_b)
            superk_xo(k+NSTATES*18,1,i_a,i_b)=pudpywave(m,i_a,i_b)
            superk_xo(k+NSTATES*19,1,i_a,i_b)=pbdx2y2d1(m,i_a,i_b)
            superk_xo(k+NSTATES*20,1,i_a,i_b)=pbdx2y2d2(m,i_a,i_b)
            superk_xo(k+NSTATES*21,1,i_a,i_b)=pbdx2y2d12(m,i_a,i_b)
        end do

        m=8*NSTATES
        k=0
        do i=1,(ntypes-8)*NSTATES
            k=k+1; m=m+1
            superk_xd(k,1,i_a,i_b)=superk_xd(k,1,i_a,i_b)&
            &+(corptop(m,i_a,i_b)/corptop(0,i_a,i_b)-superk_xo(k,1,i_a,i_b))

            superk_xd(k,2,i_a,i_b)=superk_xd(k,2,i_a,i_b)&
            &+(corptop(m,i_a,i_b)/corptop(0,i_a,i_b)-superk_xo(k,1,i_a,i_b))**2

            superk_xd(k,4,i_a,i_b)=(corptop(m,i_a,i_b)/corptop(0,i_a,i_b)&
            &-superk_xo(k,1,i_a,i_b))
        end do
    end do
end do

write(*,*) iblk,eavg,sum(wgtwlkr)

return
END subroutine BlkMeas


SUBROUTINE RunMeas
use cpmc
real(sp)::stdtop,stdbtm,covar,arg
integer::i,j,iave,i_a,i_b
external :: writeCorrDat,writePairDat,writeUnpairDat,writeVertexDat,cpOut
top=top/real(nblk)
btm=btm/real(nblk)
prd=prd/real(nblk)
hbin=hbin/real(nblk)
hbinsq=hbinsq/real(nblk)

do i_a=1,NLA
    do i_b=1,NLA
        do j=1,3
            do iave=0,NAVE
                cor_top(iave,j,i_a,i_b) = cor_top(iave,j,i_a,i_b)/real(nblk)
            end do
        end do
    end do
end do
energy=hbin

write(*,*) 'energy is the:',energy
write(*,*) 'energy Per N is the:',energy/lxy

do i_a=1,NLA
    do i_b=1,NLA

        do iave=1,NAVE
            cor_final(iave,i_a,i_b) = cor_top(iave,1,i_a,i_b)
        end do

        do iave=1,NAVE         ! error
            cor_top(iave,6,i_a,i_b)=sqrt(abs(cor_top(iave,2,i_a,i_b)&
            & -cor_top(iave,1,i_a,i_b)**2) /real(nblk-1))
        end do

        do j=1,2
            do i=1,(ntypes-8)*NSTATES*NSTATES
                super_xo(i,j,i_a,i_b)=super_xo(i,j,i_a,i_b)/real(nblk)
                super_xd(i,j,i_a,i_b)=super_xd(i,j,i_a,i_b)/real(nblk)
            end do
        end do

        do i=1,(ntypes-8)*NSTATES*NSTATES
            super_xo(i,3,i_a,i_b)=sqrt(abs(super_xo(i,2,i_a,i_b)-super_xo(i,1,i_a,i_b)**2)/real(nblk-1))
            super_xd(i,3,i_a,i_b)=sqrt(abs(super_xd(i,2,i_a,i_b)-super_xd(i,1,i_a,i_b)**2)/real(nblk-1))
        end do

        do j=1,2
            do i=1,(ntypes-8)*NSTATES
                superk_xo(i,j,i_a,i_b)=superk_xo(i,j,i_a,i_b)/real(nblk)
                superk_xd(i,j,i_a,i_b)=superk_xd(i,j,i_a,i_b)/real(nblk)
            end do
        end do

        do i=1,(ntypes-8)*NSTATES
            superk_xo(i,3,i_a,i_b)=sqrt(abs(superk_xo(i,2,i_a,i_b)-superk_xo(i,1,i_a,i_b)**2)/real(nblk-1))
            superk_xd(i,3,i_a,i_b)=sqrt(abs(superk_xd(i,2,i_a,i_b)-superk_xd(i,1,i_a,i_b)**2)/real(nblk-1))
        end do
    end do
end do

write(*,*)
write(*,*)'Result:'
if(nblk <= 1)then
    write(*,*) 'No Correlations'
    write(*,'((a),f11.7)')'Energy=',energy
else
    do i_a=1,NLA
        do i_b=1,NLA
            call writeCorrDat(i_a,i_b)
            call writePairDat(i_a,i_b)
            call writeVertexDat(i_a,i_b)
            call writeUnpairDat(i_a,i_b)
            call cpOut(i_a,i_b)
        end do
    end do
end if

RETURN
END subroutine RunMeas


subroutine writeCorrDat(idx_a, idx_b)
use cpmc
implicit none
integer, intent(in) :: idx_a, idx_b
integer :: j

numberup=0.0
numberdn=0.0
do j=1,nsites
    numberup=numberup+cor_final(j,idx_a, idx_b)
    numberdn=numberdn+cor_final(nsites+j,idx_a, idx_b)
end do
write(filename, '(A,I0.3,A,I0.3,A)') 'corr_', idx_a, '_', idx_b, '.dat'
open(unit=21, file=filename, status='unknown', action='write')
write(21,*) '============= ONE BAND Hubbard Model========'
write(21,'((a20),i5,(a10),i5,(a10),i5,(a10),f5.3)') &
!&'lattice size=',nsites,'up=',NUP,'down=',NDN,'alpha=',alpha,'alphat1=',alphat1
&'lattice size=',nsites,'up=',NUP,'down=',NDN
! write(21,'((a10),f10.5,(a10),f10.5,(a10),f10.5,(a10),f10.5,(a10),f10.5,(a10),f10.5)') &
! &'UdA=',udA,'UdB=',udB,'tpd=',tpd,'tpp=',tpp,'vpd=',vpd
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

! write(*,*) 'Correlation Functions in File corr.dat'

write(21,*) '********back propagation results:'
write(21,*) '                                 '
write(21,*) 'Charge__up:',&
&cor_final(NAVE-5,idx_a, idx_b),cor_top(NAVE-5,6,idx_a, idx_b)
write(21,*) 'Charge__dn:',&
&cor_final(NAVE-4,idx_a, idx_b),cor_top(NAVE-4,6,idx_a, idx_b)
write(21,*) 'magnetic moment:',&
&cor_final(NAVE-3,idx_a, idx_b),cor_top(NAVE-3,6,idx_a, idx_b)
write(21,*) 'kinetic energy:',&
&cor_final(NAVE-2,idx_a, idx_b),cor_top(NAVE-2,6,idx_a, idx_b)
write(21,*) 'potential energy:',&
&cor_final(NAVE-1,idx_a, idx_b),cor_top(NAVE-1,6,idx_a, idx_b)
write(21,*) 'total energy:',&
&cor_final(NAVE,idx_a, idx_b),cor_top(NAVE,6,idx_a, idx_b)

write(21,*) 'total energy per N:',&
&cor_final(NAVE,idx_a, idx_b)/lxy,cor_top(NAVE,6,idx_a, idx_b)/lxy

END subroutine writeCorrDat
!**************************************************!

subroutine writePairDat(idx_a, idx_b)
use cpmc
implicit none
integer, intent(in) :: idx_a, idx_b
integer :: i, j, k, m, n
real :: kx, ky
external :: system
character(len=30), parameter :: kNames(30) = (/ &
    'n_up','n_dn','n_w','cdw','sdwz','sdwx',&
    'pmdf','dsf','swave','pwave','dwave','sowave','sbwave','pbwave',&
    'dbwave','pdsfbd1','pdsfbd2','pdsfbd12','dd1wave','dd2wave','dd12wave',&
    'puupxwave','pddpxwave','pudpxwave','puupywave','pddpywave','pudpywave',&
    'pbdx2y2d1','pbdx2y2d2','pbdx2y2d12' /)

call system('mkdir -p ' // trim('dir-kVals/'))
do j = 1, ntypes
    write(filename, '(A,I0,A,I0,A,A)') 'dir-kVals/(', idx_a, ',', idx_b, ')'//trim(kNames(j))//'.dat'
    open(unit=100+j, file=trim(filename), status='unknown')
    ! open(unit=100+j, file='dir-kVals/'//trim((i_a,i_b)kNames(j))//'.dat', status='unknown')
    write(100+j,'(4a15)') 'k_x','k_y', 'value', 'error'
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
        write(100+j,'(2F15.9,2F15.9)') kx, ky, cor_final(k,idx_a, idx_b), cor_top(k,6,idx_a, idx_b)
        if (kx == pi) then
            write(100+j, '(2f15.9,2f15.9)') -kx, ky, cor_final(k,idx_a, idx_b), cor_top(k,6,idx_a, idx_b)
        end if
        if (ky == pi ) then
            write(100+j, '(2f15.9,2f15.9)') kx, -ky, cor_final(k,idx_a, idx_b), cor_top(k,6,idx_a, idx_b)
        end if
        if(kx == pi .and. ky == pi )  then
            write(100+j, '(2f15.9,2f15.9)') -kx, -ky, cor_final(k,idx_a, idx_b), cor_top(k,6,idx_a, idx_b)
        end if
    end do
    close(100+j)
end do

call system('mkdir -p ' // trim('dir-rVals/'))
do j = 1, ntypes
    write(filename, '(A,I0,A,I0,A,A)') 'dir-rVals/(', idx_a, ',', idx_b, ')'//trim(kNames(j))//'.dat'
    open(unit=200+j, file=trim(filename), status='unknown')
    ! open(unit=200+j, file='dir-rVals/'//trim(kNames(j))//'.dat', status='unknown')
    write(200+j,'(4a15)')  'r_x', 'r_y', 'value', 'error'
    do m = 1, NSTATES
        do n = 1, NSTATES
            k = (j-1)*NSTATES*NSTATES + (m-1)*NSTATES + n
            write(200+j,'(2I8,2F15.9)') m, n, cor_final(k,idx_a, idx_b), cor_top(k,6,idx_a, idx_b)
        end do
    end do
    close(200+j)
end do
END subroutine writePairDat

subroutine writeVertexDat(idx_a, idx_b)
use cpmc
implicit none
integer, intent(in) :: idx_a, idx_b
integer :: i, j, k, m, n
real :: kx, ky
external :: system

character(len=30), parameter :: kNames(22) = (/ &
    'Vertex_swave','Vertex_pwave','Vertex_dwave','Vertex_sowave','Vertex_sbwave','Vertex_pbwave',&
    'Vertex_dbwave','Vertex_pdsfbd1','Vertex_pdsfbd2','Vertex_pdsfbd12','Vertex_dd1wave','Vertex_dd2wave','Vertex_dd12wave',&
    'Vertex_puupxwave','Vertex_pddpxwave','Vertex_pudpxwave','Vertex_puupywave',&
    'Vertex_pddpywave','Vertex_pudpywave',&
    'Vertex_pbdx2y2d1','Vertex_pbdx2y2d2','Vertex_pbdx2y2d12' /)
k=0
call system('mkdir -p ' // trim('dir-kVals/'))
do j = 1, (ntypes-8)
    write(filename, '(A,I0,A,I0,A,A)') 'dir-kVals/(', idx_a, ',', idx_b, ')'//trim(kNames(j))//'.dat'
    open(unit=300 + j, file=trim(filename), status='unknown')
    ! open(unit=300+j, file='dir-kVals/'//trim(kNames(j))//'.dat', status='unknown')
    write(300+j,'(4a15)') 'k_x','k_y', 'value', 'error'
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
        write(300+j,'(2F15.9,2F15.9)') kx, ky, superk_xd(k,1,idx_a, idx_b),superk_xd(k,3,idx_a, idx_b)

        if (kx == pi) then
            write(300+j, '(2f15.9,2f15.9)') -kx, ky, superk_xd(k,1,idx_a, idx_b),superk_xd(k,3,idx_a, idx_b)
        end if
        if (ky == pi ) then
            write(300+j, '(2f15.9,2f15.9)') kx, -ky, superk_xd(k,1,idx_a, idx_b),superk_xd(k,3,idx_a, idx_b)
        end if
        if(kx == pi .and. ky == pi )  then
            write(300+j, '(2f15.9,2f15.9)') -kx, -ky, superk_xd(k,1,idx_a, idx_b),superk_xd(k,3,idx_a, idx_b)
        end if

    end do
    close(300+j)
end do

k=0
call system('mkdir -p ' // trim('dir-rVals/'))
do j = 1, (ntypes-8)
    write(filename, '(A,I0,A,I0,A,A)') 'dir-rVals/(', idx_a, ',', idx_b, ')'//trim(kNames(j))//'.dat'
    open(unit=400+j, file=trim(filename), status='unknown')
    ! open(unit=400+j, file='dir-rVals/'//trim(kNames(j))//'.dat', status='unknown')
    write(400+j,'(4a15)')  'r_x', 'r_y', 'value', 'error'
    do m = 1, NSTATES
        do n = 1, NSTATES
            k=k+1
            write(400+j,'(2I8,2F15.9)') m, n, super_xd(k,1,idx_a, idx_b),super_xd(k,3,idx_a, idx_b)
        end do
    end do
    close(400+j)
end do
END subroutine writeVertexDat

subroutine writeUnpairDat(idx_a, idx_b)
use cpmc
implicit none
integer, intent(in) :: idx_a, idx_b
integer :: i, j, k, m, n
real :: kx, ky
external :: system

character(len=30), parameter :: kNames(22) = (/ &
    'Unpair_swave','Unpair_pwave','Unpair_dwave','Unpair_sowave','Unpair_sbwave','Unpair_pbwave',&
    'Unpair_dbwave','Unpair_pdsfbd1','Unpair_pdsfbd2','Unpair_pdsfbd12','Unpair_dd1wave','Unpair_dd2wave','Unpair_dd12wave',&
    'Unpair_puupxwave','Unpair_pddpxwave','Unpair_pudpxwave','Unpair_puupywave','Unpair_pddpywave','Unpair_pudpywave',&
    'Unpair_pbdx2y2d1','Unpair_pbdx2y2d2','Unpair_pbdx2y2d12' /)
k=0
call system('mkdir -p ' // trim('dir-kVals/'))
do j = 1, (ntypes-8)
    write(filename, '(A,I0,A,I0,A,A)') 'dir-kVals/(', idx_a, ',', idx_b, ')'//trim(kNames(j))//'.dat'
    open(unit=500 + j, file=trim(filename), status='unknown')
    ! open(unit=500+j, file='dir-kVals/'//trim(kNames(j))//'.dat', status='unknown')
    write(500+j,'(4a15)') 'k_x','k_y', 'value', 'error'
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
        write(500+j,'(2F15.9,2F15.9)') kx, ky, superk_xo(k,1,idx_a, idx_b),superk_xo(k,3,idx_a, idx_b)

        if (kx == pi) then
            write(500+j, '(2f15.9,2f15.9)') -kx, ky, superk_xo(k,1,idx_a, idx_b),superk_xo(k,3,idx_a, idx_b)
        end if
        if (ky == pi ) then
            write(500+j, '(2f15.9,2f15.9)') kx, -ky, superk_xo(k,1,idx_a, idx_b),superk_xo(k,3,idx_a, idx_b)
        end if
        if(kx == pi .and. ky == pi )  then
            write(500+j, '(2f15.9,2f15.9)') -kx, -ky, superk_xo(k,1,idx_a, idx_b),superk_xo(k,3,idx_a, idx_b)
        end if

    end do
    close(500+j)
end do

k=0
call system('mkdir -p ' // trim('dir-rVals/'))
do j = 1, (ntypes-8)
    write(filename, '(A,I0,A,I0,A,A)') 'dir-rVals/(', idx_a, ',', idx_b, ')'//trim(kNames(j))//'.dat'
    open(unit=600+j, file=trim(filename), status='unknown')
    ! open(unit=600+j, file='dir-rVals/'//trim(kNames(j))//'.dat', status='unknown')
    write(600+j,'(4a15)')  'r_x', 'r_y', 'value', 'error'
        do m = 1, NSTATES
            do n = 1, NSTATES
                k=k+1
                write(600+j,'(2I8,2F15.9)') m, n, super_xo(k,1,idx_a, idx_b),super_xo(k,3,idx_a, idx_b)
            end do
        end do
    close(600+j)
end do
END subroutine writeUnpairDat

!-----------------HoKinAnHop------------------!
!Purpose: For convenient input in other plotting
      SUBROUTINE cpOut(idx_a, idx_b)
      use cpmc
      integer, intent(in) :: idx_a, idx_b
      integer:: i,j,k
        open(unit=50,file='out.dat',status='unknown')
        write(50,*) 'n_site'
        write(50,*) nsites
        write(50,*) 'n_e'
        write(50,*) ne
        write(50,*) 'filling'
        write(50,*) real(NE)/real(nsites)
        write(50,*) 'n_up'
        write(50,*) NUP
        write(50,*) 'n_dn'
        write(50,*) NDN
        ! write(50,*) 'Ud'
        ! write(50,*) ud
        ! write(50,*) 't0'
        ! write(50,*) t0
        ! write(50,*) 't1'
        ! write(50,*) t1
        ! write(50,*) 't2'
        ! write(50,*) t2
        ! write(50,*) 'Vpd'
        ! write(50,*) vpd
        !write(50,*) 'alpha'
        !write(50,*) alpha
        write(50,*) 'tam'
        write(50,*) tam
        write(50,*) 'alphat1'
        write(50,*) alphat1
        write(50,*) 'n_walker'
        write(50,*) NWLKRS
        write(50,*) 'dt'
        write(50,*) deltau
        write(50,*) 'bpStep'
        write(50,*) itvl_m
        write(50,*) 'ke'
        write(50,*) cor_top(NAVE-2,1,idx_a, idx_b),cor_top(NAVE-2,6,idx_a, idx_b)
        write(50,*) 'pe'
        write(50,*) cor_top(NAVE-1,1,idx_a, idx_b),cor_top(NAVE-1,6,idx_a, idx_b)
        write(50,*) 'totalEn'
        write(50,*) cor_top(NAVE,1,idx_a, idx_b),cor_top(NAVE,6,idx_a, idx_b)
        write(50,*) 'kSpace'
        do i=1,NSTATES
            write(50,'(2f15.9)') kSet(i,1),kSet(i,2)
        enddo
        ! write(50,*) 'Npair'
        ! k=(3*NSTATES+9)*NSTATES+5*NPAIR
        ! do i=1,NSTATES
        !     k=k+1
        !     write(50,'(2f15.9)')   cor_top(k,1,idx_a, idx_b),cor_top(k,6,idx_a, idx_b)
        ! enddo
        ! write(50,*) 'Dk'
        ! k=(3*NSTATES+9)*NSTATES+5*NPAIR
        ! do i=1,NSTATES
        !     k=k+1
        !     write(50,'(2f15.9)')   cor_top(k,1,idx_a, idx_b),cor_top(k,6,idx_a, idx_b)
        ! enddo
        close(50)
      endsubroutine
!-----------------HoKinAnHop------------------!

subroutine initout
  use cpmc
!!!!!!WRITE CONTROL PARAMETERS TO STDOUT
write(*,*)
write(*,*)'=================================================='
write(*,*)' CPMC: one-band Hubbard model'
write(*,*)' singlet (closed shell) ground-state'
write(*,*)'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write(*,*)
! write(*,*)'Parameters: tpd,tpp,udA,udb,vpd'
! write(*,*) tpd,tpp,udA,udB,vpd
write(*,*)
write(*,*)'# of states in SCF solution: NSTATES=',NSTATES
write(*,*)'# of up (or down) electrons: NE=',NE,NUP,NDN
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
endsubroutine
