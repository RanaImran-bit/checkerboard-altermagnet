module mbss
    !implicit none
    !use f90_kind
    ! ... for MPI environment
    !      include 'mpif.h'
    integer::nsize,ncpu,ocpu,ierror
    ! symbolic names for kind types of
    !  single- and double-precision reals:
    integer,parameter::k4b=selected_int_kind(9)
    integer,parameter::i4b=selected_int_kind(9)
    integer,parameter::i2b=selected_int_kind(2)
    integer,parameter::sp=kind(1.0d0)
    integer,parameter::dp=kind(1.0d0)
    integer,parameter::dpc=kind((1.0d0,1.0d0))
    !in visual fortranGen: k4b=4, i4b=4, i4b=1, sp=8, dpc=8
    
    !defining constant
    double precision::pi,zero,one,two,four,ten,half,fourth,sup,sdown
    double complex:: onec,zeroc
    parameter (PI=3.1415926535897932384626433)
    parameter (ZERO=0.0,FOURTH=0.25,HALF=0.5)
    parameter (ONE=1.0,TWO=2.0,FOUR=4.0,TEN=10.0)
    parameter (sup=1.0,sdown=-1.0)
    parameter (onec=dcmplx(1.0d0,0.0d0),zeroc=dcmplx(0.0d0,0.0d0))
    ! ... for ranGendom number generator
    integer::ISEED
    
    ! ... declarations
    integer::NX,NY,NCN,NNNC,NPAIR,NXY,NT,NSPIN,NXYS,lx,ly,NK,NLA
    integer::fAVG, fRMS, dAVG, dRMS, isudv,switchg
    integer::MEAS,MEAS1,MEAS2,MEAS3,MEAS4,MEAS5,nsign
    integer::MEAC1,MEAC2,MEAC3
    integer::MEAS11,MEAS12,MEAS13,MEAS14,MEAS21,MEAS22,NKK,NWM
    integer::MEAS_VTX_R, MEAS_VTX_K !vertex
    ! ... nsign govern the sign of a whole space-time lattice sweep, initialize in cnfinit, update in cnfmake.
    parameter (isudv=8)
    !       1-NX*NY   =   conduction band e-
    !       NX*NY+1   =   impurity
    parameter (NX=8, NY=8, NLA=1, NCN=4, NNNC=4, NK=NX*NY, NXY=NX*NY*NLA, NSPIN=2, NXYS=NXY*NSPIN,&
     NKK=NX*NY*NLA*NLA,NWM=5)
    !     1---NXY stand for f band, NXY+1---NXY*NB stand for d band
    !        NX -- number of lattice sites in x-direction
    !        NY -- number of lattice sites in y-direction
    !        NT -- maximum number of time steps
    !        NSPIN -- number of components of electron spin
    !        NCN -- coordinate number4 NCN for 2d, 1 NCN for two band, 6 for 1d triangle lattice, 3 for honeycomb
    ! ... all energies are in units of the half-band width
    
    ! ... parameters and derived constants
    double precision::dt
    double complex,dimension(:,:,:),allocatable::th2,th3
    double precision,dimension(:),allocatable::ue2,mu2
    ! independent site input from file
    !     th -- hopping parameter in x, y direction
    !     th1 -- hopping parameter in diaganol
    !     ue -- electrostatic energy
    !     mu(nspin) -- chemical potential
    !     dt -- spacing in the inverse temperature
    !     warms -- number of sweeps to achieve equilibrium!     runs -- number of Monte Carlo runs per job
    !     sweeps -- number of sweeps per Monte Carlo run
    !     sitefile -- character string, from x,y coordinates to site number
    !     pairfile -- from origin to all different points in the lattice considering symmetric
    !     kfile -- k vector of the k-space corresponding to the xspace
    double precision,dimension(:),allocatable::jhs
    double precision:: beta
    integer:: warms, sweeps, runs, ntimes
    
    !     1---NXY stand for f band, NXY+1---NXY*NB stand for d band
    !     emh -- exp(-dt*H), the H part of the B matrix,
    !            H has the td,tf,vh hopping terms and uf/2-muf, ud/2-mud terms,
    !            H does not have the HS term
    !     eph -- exp(+dt*H), the H part of the I/B matrix
    !     emv -- exp(-dt*V)
    !     epv -- exp(+dt*V)
    !     expdv -- the exponential of the change in V when
    !              a spin is flipped
    double complex,dimension(:,:),allocatable::emh,eph
    double complex,dimension(:,:),allocatable::emv, epv
    double complex,dimension(:,:),allocatable::expdv
    integer::nrolls, maxrolls
    ! ... the rolling at the end of each measurement equals to 8
    ! ... x-space connection
    integer,dimension(:,:,:),allocatable::iposit
    integer,dimension(:,:),allocatable::nextsite,nnnsite
    integer,dimension(:),allocatable::ixv,iyv
    integer,dimension(:),allocatable::ilav
    double precision,dimension(:),allocatable::xpos,ypos
    integer,dimension(:,:),allocatable::npx,npy
    ! ... k-space
    double precision,dimension(:),allocatable::kx,ky
    ! ... pair structure factor
    double precision,dimension(:),allocatable::sefac,dfac,pfac,pupfac
    ! ... Ising field at every time-space.
    integer,dimension(:,:),allocatable::spin
    !	  double precision,dimension(:,:),allocatable::spinaccu
    ! ... Green's function.
    !       g=<ci(tau)cj^+(0)>
    !       gt=<c c^+> g2t=<c^+ c>
    double complex,dimension(:,:),allocatable::g
    double complex,dimension(:,:,:),allocatable::gt
    double complex,dimension(:,:,:),allocatable::g2t
    !       gt1=<ci(tau)cj^+(0)>
    !       gt2=<ci^+(tau)cj(0)>
    double complex,dimension(:,:,:),allocatable::gt1,gt2
    ! ... UDV space
    double complex,dimension(:,:),allocatable::u, ui, v, vi
    double complex,dimension(:),allocatable::s
    double complex,dimension(:,:,:),allocatable::usar,vsar
    double complex,dimension(:,:),allocatable::ssar
    double complex,dimension(:,:,:),allocatable::usal,vsal
    double complex,dimension(:,:),allocatable::ssal
    
    !for [(u^-1_l)*(v^-1_r)+d_l*v_l*v_r*d_r]^-1
    double complex,dimension(:,:,:),allocatable::uprod,vprod
    double complex,dimension(:,:),allocatable::sprod
    ! ... configuration avgs (average over time(l) index)
    parameter (fAVG=1, fRMS=2, dAVG=3, dRMS=4)
    !Acc for accumulation during each run
    !role:  acc ==> mah
    !       kSpace ==> caverl
    !       final ==> ram
    double precision,dimension(:),allocatable::caverl
    ! ... results of average of nmeas  Monte Carlo runs.
    double precision,dimension(:,:),allocatable::ram
    double precision:: rsign(2)
    ! ... results of average over nswps  important sampling in HS field in one Monte Carlo run
    double precision,dimension(:),allocatable::mah
    ! complex no
    double complex,dimension(:),allocatable::caverlc
    double complex,dimension(:,:),allocatable::ramc
    double complex,dimension(:),allocatable::mahc
    
    ! ... negative U case
    logical:: negativeUSites(NXY)
    !1 - NXY : up
    !NXY+1 - 2*NXY : down
    double complex:: spinExtraRatio(-1:1,NXYS)
    !   input from fort.500 file
    double precision:: ue,mu
    double precision:: lamda,h,hx,tam,ttp,ttn
    !   the name for the variable
    character(10),dimension(:),allocatable:: name1,name2,name3
    
    !group of B
    double complex,dimension(:,:,:),allocatable:: Bg
    
    integer:: nb
    
    character(len=50):: NXYSstr
    character(len=50):: NXYSstrE
    end module mbss
    
    module link
    interface
    function ranGen(idum) !0
    use mbss
    double precision::ranGen
    integer,intent(inout)::idum
    end function ranGen
    
    subroutine AllocatingVariables!1.1.1
    use mbss
    endsubroutine AllocatingVariables
    
    subroutine AllocatingSaveMem!1.1.2
    use mbss
    endsubroutine AllocatingSaveMem
    
    subroutine DeallocatingVariables!1.2
    use mbss
    endsubroutine DeallocatingVariables
    
    subroutine sysdef (nwarms,nmeas,nswps) !1
    use mbss
    integer::nwarms,nmeas,nswps
    end subroutine sysdef
    
    subroutine sysinit !2
    use mbss
    integer:: j
    end subroutine sysinit
    
    subroutine cnfinit !3
    use mbss
    end subroutine cnfinit
    
    subroutine sysequil(nwarms) !4
    use mbss
    integer::nwarms
    end subroutine sysequil
    
    subroutine initMem(acceptc,negs)
    use mbss
    double precision::negs,acceptc
    endsubroutine initMem
    
    subroutine recMem
    use mbss
    endsubroutine recMem
    
    subroutine normMem(nswps)
    use mbss
    integer:: i,j,nswps
    endsubroutine normMem
    
    subroutine sysmeas(nswps,negs)!5
    use mbss
    integer::nswps
    double precision::negs
    end subroutine sysmeas
    
    subroutine FTSF(ns,cdwIn,sdwxyIn) !5.1.1
    use mbss
    double precision::ns(NXY),cdwIn(NXY,NXY),sdwxyIn(NXY,NXY)
    end subroutine FTSF
    
    subroutine FTSus(cdwtIn,sdwtIn,sdwxytIn) !5.1.2
    use mbss
    double precision::cdwtIn(NXY,NXY),sdwtIn(NXY,NXY),sdwxytIn(NXY,NXY)
    end subroutine FTSus
    
    subroutine FTSsupBar(gt1vt,gt2vt,g2In) !5.1.3
    use mbss
    double complex::gt2vt(NXYS,NXYS,0:NT), gt1vt(NXYS,NXYS,0:NT), g2In(NXYS,NXYS)
    end subroutine FTSsupBar
    
    subroutine sumSus(cdwtIn,sdwtIn,sdwxytIn) !5.1.4
    use mbss
    double precision::cdwtIn(NXY,NXY),sdwtIn(NXY,NXY),sdwxytIn(NXY,NXY)
    end subroutine sumSus
    
    subroutine FTDF(pmdf,dsf)
    use mbss
    double precision::pmdf(NXY,NXY),dsf(NXY,NXY)
    endsubroutine FTDF
    
    subroutine FTGt(gt2In,g1In,g2In)
    use mbss
    double complex:: gt2In(NXYS,NXYS,0:NT),g1In(NXYS,NXYS),g2In(NXYS,NXYS)
    endsubroutine
    
    subroutine FTPC
    endsubroutine FTPC
    
    subroutine FTCC
    endsubroutine FTCC
    
    subroutine CorSub
    endsubroutine CorSub
    
    subroutine measana(nmeas) !5.5
    use mbss
    integer::nmeas
    end subroutine measana
    
    subroutine cnfmake(acceptc,negs)!6
    use mbss
    double precision::acceptc,negs
    end subroutine cnfmake
    
    subroutine chkg(it)!7
    use mbss
    integer::it,ispin
    end subroutine chkg
    
    subroutine makeg(it)!8
    use mbss
    integer::it,ispin
    end subroutine makeg
    
    subroutine makeb !9
    use mbss
    end subroutine makeb
    
    subroutine udvb !9.05
    use mbss
    end subroutine udvb
    
    subroutine makebt !9.1
    use mbss
    end subroutine makebt
    
    subroutine udvbt
    use mbss
    end subroutine udvbt
    
    subroutine saveb8 !9.2
    use mbss
    end subroutine saveb8
    
    subroutine makeipb !10
    use mbss
    end subroutine makeipb
    
    subroutine matinv(a) !11
    use mbss
    double complex::a(NXYS,NXYS)
    end subroutine matinv
    
    subroutine cnfmeas !12
    use mbss
    end subroutine cnfmeas
    
    subroutine gtime(tgt) !13
    use mbss
    integer::tgt
    end subroutine gtime
    
    subroutine guet(im,clcp0,cplc0) !13.1
    use mbss
    integer::im
    double precision::clcp0(1:NXYS,1:NXYS),cplc0(1:NXYS,1:NXYS)
    end subroutine guet
    
    subroutine oneSiteMeas(nOut,mOut,n2Out,m2Out,upDnOut,vacOut,nsOut,msOut,upDnsOut, &!13.3
    upSiteOutput,dnsOut,n2sOut,m2sOut,enH,enMu,enU,enLan)
    use mbss
    !define local variable
    double precision::up,dn
    !Ahout energy
    double precision::enH,enMu,enU,enLan
    !about spin n/m
    double precision::nOut,mOut,n2Out,m2Out,upDnOut,vacOut
    double precision::nsOut(NXY),msOut(NXY),upDnsOut(NXY), &
    upsOut(NXY),dnsOut(NXY),n2sOut(NXY),m2sOut(NXY)
    end subroutine oneSiteMeas
    
    subroutine twoSiteMeas(cdwIn,sdwxyIn,pmdfIn,dsfIn) !13.4
    use mbss
    double precision::up1(0:NT),up2(0:NT)
    double precision::dn1(0:NT),dn2(0:NT)
    double precision::up_up(0:NT),dn_dn(0:NT),up_dn(0:NT)
    double precision::cdwIn(NXY,NXY),sdwxyIn(NXY,NXY)
    double precision::pmdfIn(NXY,NXY),dsfIn(NXY,NXY)
    endsubroutine twoSiteMeas
    
    subroutine twoSitetMeas(cdwtIn,sdwtIn,sdwxytIn) !13.5
    use mbss
    
    double precision::cdwtIn(NXY,NXY),sdwtIn(NXY,NXY),sdwxytIn(NXY,NXY)
    double precision::up_up1(0:NT),dn_dn1(0:NT),up_dn1(0:NT),dn_up1(0:NT)
    double precision::up1(0:NT),up2(0:NT)
    double precision::dn1(0:NT),dn2(0:NT)
    end subroutine twoSitetMeas
    
    subroutine udv(NA,uh,sh,vh)
    use mbss
    double complex::uh(NA,NA),sh(NA),vh(NA,NA)
    double complex::xn1(NA)
    double complex::temp(NA,NA)
    double precision::xnorm(NA),vhelp(NA),xmax
    double complex::v1(NA,NA)
    !double precision::test(NA,NA),test1(NA,NA)
    double complex::zeta(NA),work(4*NA)
    end subroutine udv
    
    subroutine bpmult(it,gh) !17
    use mbss
    integer::it
    double complex::gh(NXYS,NXYS)
    end subroutine bpmult
    
    subroutine bmmult(it,gh) !18
    use mbss
    integer::it
    double complex::gh(NXYS,NXYS)
    end subroutine bmmult
    
    subroutine multbm (it,gh) !19
    use mbss
    integer::it
    double complex::gh(NXYS,NXYS)
    end subroutine multbm
    
    subroutine multbp (it,gh) !20
    use mbss
    integer::it
    double complex::gh(NXYS,NXYS)
    end subroutine multbp
    
    subroutine analysisOutput
    use mbss
    endsubroutine analysisOutput
    
    function trans(NA,ma)
    integer::NA,i,j
    double complex::ma(NA,NA),trans(NA,NA)
    endfunction
    
    double complex function det(N,mat)
    integer:: N
    double complex,intent(inout),dimension(N,N):: mat
    end function det
    
    subroutine makeA (ib)
    integer:: ib
    endsubroutine makeA
    
    subroutine updateBg(ib)
    integer::ib
    endsubroutine updateBg
    
    subroutine makeBg
    endsubroutine makeBg
    
    subroutine makeANoCh (ib)
    integer:: ib,ispin
    endsubroutine makeANoCh
    
    subroutine udvBg
    endsubroutine udvBg
    end interface
    endmodule link
    
    !0~~~~~~~~~~~~random number generator from Numerical Recipe~~~~~~~~~~!
    function ranGen(idum)
    use mbss
    implicit none
    double precision::ranGen
    integer,intent(inout)::idum
    integer,parameter::ia=16807,im=2147483647,iq=127773,ir=2836
    double precision,save::am
    integer,save::ix=-1,iy=-1,kp
    if(idum<=0 .or. iy<0) then
        am=nearest(1.0,-1.0)/im
        iy=ior(ieor(888889999,abs(idum)),1)
        ix=ieor(777755555,abs(idum))
        idum=abs(idum)+1
    endif
    ix=ieor(ix,ishft(ix,13))
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    kp=iy/iq
    iy=ia*(iy-kp*iq)-ir*kp
    if(iy<0) iy=iy+im
    ranGen=am*ior(iand(im,ieor(ix,iy)),1)
    return
    end function ranGen
    
    !1.1.1 initialize allocatable variable
    subroutine AllocatingVariables
    use mbss
    implicit none
    !NXYS
    !1-NXY : up spin
    !NXY+1 - NXYS : dn spin
    allocate(th2(NXY,NSPIN,NCN))
    allocate(th3(NXY,NSPIN,NNNC))
    allocate(ue2(2*NXY),mu2(2*NXY))
    allocate(jhs(NXYS))
    allocate(emh(NXYS,NXYS))
    allocate(eph(NXYS,NXYS))
    allocate(emv(-1:1,NXYS),epv(-1:1,NXYS))
    allocate(expdv(-1:1,NXYS))
    allocate(iposit(-2*NX:2*NX,-2*NY:2*NY,NLA),nextsite(NCN,NXY),nnnsite(NNNC,NXY))
    allocate(ixv(NXY),iyv(NXY))
    allocate(xpos(NXY),ypos(NXY))
    allocate(ilav(NXY))
    allocate(sefac(NCN),dfac(NCN),pfac(NCN),pupfac(NCN))
    allocate(kx(NXY),ky(NXY))
    allocate(spin(1:NXY,1:NT))
    allocate(g(1:NXYS,1:NXYS),gt(1:NXYS,1:NXYS,0:NT),g2t(1:NXYS,1:NXYS,0:NT))
    allocate(gt1(NXYS,NXYS,0:NT),gt2(NXYS,NXYS,0:NT))
    allocate(u(1:NXYS,1:NXYS), ui(1:NXYS,1:NXYS), v(1:NXYS,1:NXYS), vi(1:NXYS,1:NXYS))
    allocate(s(1:NXYS))
    allocate(usar(1:NXYS,1:NXYS,0:NT/isudv),vsar(1:NXYS,1:NXYS,0:NT/isudv))
    allocate(ssar(1:NXYS,0:NT/isudv))
    allocate(usal(1:NXYS,1:NXYS,0:NT/isudv),vsal(1:NXYS,1:NXYS,0:NT/isudv))
    allocate(ssal(1:NXYS,0:NT/isudv))
    allocate(uprod(1:NXYS,1:NXYS,0:NT/isudv),vprod(1:NXYS,1:NXYS,0:NT/isudv))
    allocate(sprod(1:NXYS,0:NT/isudv))
    
    allocate(Bg(NXYS,NXYS,NT/isudv))
    endsubroutine AllocatingVariables
    
    !1.1.2 allocate memory for saving
    subroutine AllocatingSaveMem
    use mbss
    implicit none
    
    !---------------variable location-------------!
    !       overall property
    !        1	sn      <n>
    !        2	sm      <m>
    !        3	sn2     <n2>
    !        4	sm2     <m2>
    !        5	s2oc    <double occ>
    !        6	s0oc    <vancacy>
    !        7	eu      <PE>
    !        8	eh      <KE>
    !        9	emu     <CE>
    !        10	elan    <LE>
    !        11	totalE
    !        12	0.0d0;
    !        13 onpair  =>on-site pair corr
    !        14 pairlen => pair correlation length
    !        15 afmssf
    !        21 spsus
    !        22 sepsus
    !        23 dpsus
    !        24 ppsus
    !        25 puppsus
    !        26 spcor
    !        27 sepcor
    !        28 dpcor
    !        29 ppcor
    !        30 puppcor
    !        31 spsussp
    !        32 sepsussp
    !        33 dpsussp
    !        34 ppsussp
    !        35 puppsussp
    !---------------site property------------------!
    !        Meas+0*Nxy	sns
    !        Meas+1*Nxy	sms
    !        Meas+2*Nxy	s2ocs
    !        Meas+3*Nxy	sds
    !        Meas+4*Nxy	sn2s
    !        Meas+5*Nxy	sm2s
    !------------------k-space propery---------------------!
    !        MEAS13+0*NKK	    csf(q)
    !        MEAS13+1*NKK	    ssfz(q)
    !        MEAS13+2*NKK	    ssfxy(q)
    !        MEAS13+3*NKK
    !        MEAS13+4*NKK	    csup(q)
    !        MEAS13+5*NKK	    ssupz(q)
    !        MEAS13+6*NKK	    ssupxy(q)
    !        MEAS13+7*NKK
    !        MEAS13+8*NKK	    ssupz(q) => single particle contribution
    !        MEAS13+9*NKK	    ssupxy(q) => single particle contribution
    !        MEAS13+10*NKK
    !        MEAS13+11*NKK        fdf(q)
    !        MEAS13+12*NKK        fdfUp(q)
    !        MEAS13+13*NKK       fdfDn(q)
    !        MEAS13+14*NKK       pmdf(q)
    !        MEAS13+15*NKK       dsf(q)
    !---------------site-site correlation------------------!
    !        Meas1+0*Nxy*Nxy	csn(i,j), <ninj> for i>j, position q=(i-1)*NXY+j ,<mimj>z for j>i, position q=(j-1)*NXY+i
    !        Meas1+1*Nxy*Nxy	sdwxy(i,j) <mimj>xy for i>j, position q=(i-1)*NXY+j
    !        Meas1+2*Nxy*Nxy	0.d0
    !        Meas1+3*Nxy*Nxy	cdw1(i,j)
    !        Meas1+4*Nxy*Nxy	sdwz1(i,j)
    !        Meas1+5*Nxy*Nxy	sdwxy1(i,j)
    !        Meas1+6*Nxy*Nxy	pmdf(i,j)
    !        Meas1+7*Nxy*Nxy	dsf(i,j)
    !        Meas1+8*Nxy*Nxy    spl
    !        Meas1+9*Nxy*Nxy    sepl
    !        Meas1+10*Nxy*Nxy    dpl
    !        Meas1+11*Nxy*Nxy    ppl
    !        Meas1+12*Nxy*Nxy    puppl
    !        Meas1+13*Nxy*Nxy    sp
    !        Meas1+14*Nxy*Nxy    sep
    !        Meas1+15*Nxy*Nxy    dp
    !        Meas1+16*Nxy*Nxy    pp
    !        Meas1+17*Nxy*Nxy    pupp
    !---------------current-current correlation--------------!
    !        Meas21+0*NT*Nxy*Nxy	curcor(i,j,tau)
    !---------------current-current correlation--------------!
    !        Meas22+0*NWM*Nxy	    curqw(i,tau)
    !---------------pair correlation------------------------!
    
    !-------------------------------------------------------!
    !        Meas4+1 nsign
    !        Meas4+2 acc
    !        Meas4+3 negs
    
    !---------------complex no------------------------------!
    !        1+0*NXYS*NXYS      gt_vt(i,j), gvt(1)=sum(1,NT) g(t,1)
    !        1+1*NXYS*NXYS      gt2(i,j,0) <c^+ c>equal time
    !        1+2*NXYS*NXYS	    gt1(i,j,0) <c c^+>equal time
    !---------------k space G-----------------------------!
    !        MEAC1+0*NK         Gkupup <c^+ c>
    !        MEAC1+1*NK         Gkdndn <c^+ c>
    !        MEAC1+2*NK         Gkupdn <c^+dn cup>
    !        MEAC1+3*NK         Gkdnup <c^+up cdn>
    !        MEAC1+4*NK         Gkupup <c c^+>
    !        MEAC1+5*NK         Gkdndn <c c^+>
    !        MEAC1+6*NK         Gkupdn <c c^+>
    !        MEAC1+7*NK         Gkdnup <c c^+>
    !---------------------Green function------------------!
    !        MEAC2+0*(NT+1)*Nxys*Nxys	gt2(i,j,time)cplc0 <c^+ c>
    !        MEAC2+1*(NT+1)*Nxys*Nxys	gt1(i,j,time)clcp0 <c c^+>
    MEAS=37
    !some type of varibale in k-space only used NK space in NXY
    MEAS11=MEAS+8*NXY
    MEAS12=MEAS11
    MEAS13=MEAS12
    MEAS14 = MEAS13 + 16*NKK
    MEAS1=MEAS14+16*NKK
    MEAS21=MEAS1+18*NXY*NXY
    MEAS22=MEAS21+1*NXY*NXY*NT
    MEAS2=MEAS22+1*NXY*NWM
    MEAS3=MEAS2
    MEAS_VTX_R = MEAS3
    MEAS_VTX_K = MEAS_VTX_R + 5*NXY*NXY
    MEAS4 = MEAS_VTX_K + 5*NK
    ! MEAS4=MEAS3
    !MEAS4=MEAS3+2*(NT+1)*NXYS*NXYS
    MEAS5=MEAS4+3
    ! complex no
    MEAC1=3*NXYS*NXYS
    MEAC2=MEAC1+8*NK
    MEAC3=MEAC2+2*(NT+1)*NXYS*NXYS
    allocate(caverl(MEAS4))
    allocate(ram(MEAS4,2))
    allocate(mah(MEAS5))
    allocate(name1(MEAS5))
    allocate(caverlc(MEAC3))
    allocate(ramc(MEAC3,2))
    allocate(mahc(MEAC3))
    allocate(name2(MEAC3))
    allocate(name3(15))
    endsubroutine AllocatingSaveMem
    
    !1.2 deallocate variables
    subroutine DeallocatingVariables
    use mbss
    implicit none
    deallocate(th2)
    deallocate(ue2,mu2)
    deallocate(jhs)
    deallocate(emh, emv)
    deallocate(eph, epv)
    deallocate(expdv)
    deallocate(iposit,nextsite,nnnsite)
    deallocate(ixv,iyv)
    deallocate(xpos,ypos)
    deallocate(kx,ky)
    deallocate(spin)
    deallocate(g,gt)
    deallocate(gt1,gt2)
    deallocate(u, ui, v, vi)
    deallocate(s)
    deallocate(usar,vsar)
    deallocate(ssar)
    deallocate(usal,vsal)
    deallocate(ssal)
    deallocate(caverl)
    deallocate(ram)
    deallocate(mah)
    deallocate(caverlc)
    deallocate(ramc)
    deallocate(mahc)
    endsubroutine DeallocatingVariables
    
    !1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine sysdef (nwarms,nmeas,nswps)
    use mbss
    use link,only:sysinit,cnfinit,AllocatingVariables,AllocatingSaveMem
    implicit none
    integer::nwarms,nmeas,nswps
    character(len=64) :: input_file

    write(900+ncpu,*)'ncpu=',ncpu
    
    ! 显式打开输入文件，文件名格式为 input_XX.dat (例如 input_01.dat)
    write(input_file, '("fort.5", I2.2)') ncpu
    open(unit=500+ncpu, file=trim(input_file), status='old', action='read')

    read(500+ncpu,*) NT
    read(500+ncpu,*) dt
    read(500+ncpu,*) ue,mu
    read(500+ncpu,*) warms, runs, sweeps
    read(500+ncpu,*) lamda,tam
    read(500+ncpu,*) h,hx
    
    
    if (Mod(NT,8)/=0) then
        write(*,*) 'NT is not multiple of 8'
    endif
    !after reading==> allocating variable
    call AllocatingVariables
    call AllocatingSaveMem
    ntimes = NT
    nwarms = warms
    nmeas = runs
    nswps = sweeps
    ! ... check parameters:
    write(900+ncpu,*)'k4b',k4b,' i4b',i4b,' i2b',i2b,' sp',sp,' dpc',dpc,' dp',dp
    ! ... set misc physical parameters
    beta =dble(NT)*dt  !beta=16*0.0625=1
    
    ! ... for chkg: adaptive error processing
    nrolls=isudv   !isudv=8
    maxrolls=isudv
    
    ! ... construct hamiltonian and related matrices
    write(900+ncpu,*)'before sysinit'
    call sysinit
    
    ! ... initialize configuration dependent fields
    write(900+ncpu,*)'before cnfinit'
    call cnfinit
    write(900+ncpu,*)'after cnfinit,before output parameters'
    
    ! ... output system parameters
    write(600+ncpu,*) " "
    write(600+ncpu,*) 'Hubbard Model(one bands): TwoLegLadder, -ve U'
    write(600+ncpu,*) 'pivot udv-conditioned BSS method'
    write(600+ncpu,*) 'NX,NY,NB,NXY,NT-----system size            :',NX,NY,NB,NXY,NT
    write(600+ncpu,*) 'beta -- inverdipos(k,3)=se temperature                       :',beta
    write(600+ncpu,*) 'ue -- electrostatic energy                        :',ue
    write(600+ncpu,*) 'mu(ii)  -- chemical potential                     :',mu
    write(600+ncpu,*) 'dt -- spacing in the inverse temperature          :',dt
    write(600+ncpu,*) 'isudv -- matrix order between two stablizations   :',isudv
    write(600+ncpu,*) 'warms -- number of sweeps to equilibrate          :',nwarms
    write(600+ncpu,*) 'runs -- number of Monte Carlo runs per job        :',nmeas
    write(600+ncpu,*) 'sweeps -- number of sweeps per Monte Carlo run    :',nswps
    write(600+ncpu,*) 'iseed -- ranGendom number seed                    :',iseed
    return
    end subroutine sysdef
    
    !2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine sysinit
    use mbss
    use link,only: det
    implicit none
    integer::info,ipvt(NXYS)
    double complex:: work(NXYS)
    double complex::deter,sgn
    ! ... local varibles
    integer::i,j,k,nrt,ispin,ix,iy,ila
    double complex::htmp(NXYS,NXYS),tfm(NXYS,NXYS), tfminv(NXYS,NXYS)
    double complex::tmp(NXYS,NXYS)
    !for complex no calculation
    double complex::VL(NXYS,NXYS),VR(NXYS,NXYS),temp(NXYS,NXYS),temp2(NXYS,NXYS),temp3(NXYS,NXYS)
    double complex::eigI(NXYS)
    double complex:: workee(2*NXYS)
    double complex:: work2(4*NXYS)
    double precision:: rwork(2*NXYS)
    double precision:: calE(NK),Energy,DbEnergy,ZPartition
    write(NXYSstr,*) '(',NXYS,' f 20.15)'
    NXYSstr=trim(NXYSstr)
    nb=NT/isudv
    !ZGEEV( 'V', 'V', NXYS, htmp2, NXYS, eigval, VL, NXYS, VR, NXYS, work, 4*NXYS, rwork, INFO )
    !set x-space
    iposit=0
    k=0
    write(100+ncpu,*) 'rSpace'
    do ix=1,NX
        do iy=1,NY
            k=k+1
            iposit(ix,iy,1)=k
            ixv(k)=ix
            iyv(k)=iy
            ilav(k)=1
            xpos(k)=dble(ix)
            ypos(k)=dble(iy)
            write(100+ncpu,*) k,ix,iy,ilav(k),xpos(k),ypos(k)
        enddo
    enddo
    !for impurity located at (1,1), no iposit label
    !set k-space
    write(100+ncpu,*) 'kSpace'
    k=0
    do ix=0,NX-1
        do iy=0,NY-1
            k=k+1
            kx(k)=dble(ix-nx/2)*2.0d0*PI/dble(nx)
            ky(k)=dble(iy-ny/2)*2.0d0*PI/dble(ny)
            write(100+ncpu,*) k,kx(k),ky(k)
        enddo
    enddo
    !Energy dispersion
    do i=1,NK
        calE(i)=-2.0d0*(dcos(kx(i))+dcos(ky(i)))
    enddo
    write(160,*) 'expected E'
    write(160,'(10 f 10.5)') calE
    !Making Periodic Boundary Conditions
    do i=1,NXY
        iposit(ixv(i)+nx,iyv(i),ilav(i))=i     !!!!!----one-----!!!
        iposit(ixv(i)-nx,iyv(i),ilav(i))=i     !!!!!----two-----!!!
        iposit(ixv(i),iyv(i)+NY,ilav(i))=i
    
        iposit(ixv(i),iyv(i)-NY,ilav(i))=i
    
        iposit(ixv(i)+NX,iyv(i)+NY,ilav(i))=i
    
        iposit(ixv(i)+NX,iyv(i)-NY,ilav(i))=i
    
        iposit(ixv(i)-NX,iyv(i)+NY,ilav(i))=i
    
        iposit(ixv(i)-NX,iyv(i)-NY,ilav(i))=i
    end do
    
    ! nearest site and next-nearest site: nextsite(NCN,NXY),nnnsite(NNNC,NXY)
    !  4'   3    1'
    !  2 -- * -- 1
    !  2'   4    3'
    write(900+ncpu,*)'before set nextsite, nnnsite'
    write(100+ncpu,*) 'nextSite'
    do i=1,NXY
        nextsite(1,i)=iposit(ixv(i)+1,iyv(i),1)
    
        nextsite(2,i)=iposit(ixv(i)-1,iyv(i),1)
    
        nextsite(3,i)=iposit(ixv(i),iyv(i)+1,1)
    
        nextsite(4,i)=iposit(ixv(i),iyv(i)-1,1)
    
        nnnsite(1,i)=iposit(ixv(i)+1,iyv(i)+1,1)
        nnnsite(2,i)=iposit(ixv(i)-1,iyv(i)-1,1)
    
        nnnsite(3,i)=iposit(ixv(i)+1,iyv(i)-1,1)
        nnnsite(4,i)=iposit(ixv(i)-1,iyv(i)+1,1)
    
        write(100+ncpu,*) i
        write(100+ncpu,*) nextsite(:,i)
        write(100+ncpu,*) nnnsite(:,i)
    enddo
    !----------------------------------------------------------------!
    ! Set up pairing symmetry and relative position of nearest site:'!
    !----------------------------------------------------------------!
    !  4'   3    1'
    !  2 -- * -- 1
    !  2'   4    3'
    
    !       1
    !  1 -- * -- 1
    !       1
    sefac(1:4)=0.5d0
    !       -1
    !  1 -- * -- 1
    !       -1
    dfac(1)=0.5d0
    dfac(2)=0.5d0
    dfac(3)=-0.5d0
    dfac(4)=-0.5d0
    !        0
    !  -1 -- * -- 1
    !        0
    pfac(1)= 1.0d0 / dsqrt(2.0d0)  
    pfac(2)=-1.0d0 / dsqrt(2.0d0)  
    pfac(3)= 0.0d0 
    pfac(4)= 0.0d0 
    !        0
    !  0 -- * -- 1
    !        0
    pupfac(1)=1.0d0
    !set hopping term/ K-matrix
    !only impurity site have U
    ue2=ue
    mu2=mu
    ttp=1.0d0+tam
    ttn=1.0d0-tam
    do i=1,NXY
    !up spin case
        th2(i,1,1)=dcmplx(-ttn,0.0d0)
        th2(i,1,2)=dcmplx(-ttn,0.0d0)
        th2(i,1,3)=dcmplx(-ttp,0.0d0)
        th2(i,1,4)=dcmplx(-ttp,0.0d0)
    !down spin case
        th2(i,2,1)=dcmplx(-ttp,0.0d0)
        th2(i,2,2)=dcmplx(-ttp,0.0d0)
        th2(i,2,3)=dcmplx(-ttn,0.0d0)
        th2(i,2,4)=dcmplx(-ttn,0.0d0)
    !up-dn c_up^+ c_dn
        th3(i,1,1)=dcmplx(lamda,0.0d0)
        th3(i,1,2)=dcmplx(-lamda,0.0d0)
        th3(i,1,3)=dcmplx(0.0d0,lamda)
        th3(i,1,4)=dcmplx(0.0d0,-lamda)
    !dn-up c_dn^+ c_up
        th3(i,2,1)=dcmplx(-lamda,0.0d0)
        th3(i,2,2)=dcmplx(lamda,0.0d0)
        th3(i,2,3)=dcmplx(0.0d0,lamda)
        th3(i,2,4)=dcmplx(0.0d0,-lamda)
    enddo
    !Negative U sites picking up
    negativeUSites=.false.
    do i=1, NXY
        if (ue2(i)<0) then
            negativeUSites(i)=.true.
        endif
    enddo
    
    !set Hamiltionian H
    write(900+ncpu,*)'before set htmp' !nextsite contain the lattice information
    ! ... hamiltonian hopping terms in x-space: between nearest neighbour, on-site Coulomb interaction
    htmp=dcmplx(0.0d0,0.0d0)
    do i=1,NXY
    !up spin
        htmp(i,nextsite(1,i))=th2(i,1,1)
        htmp(i,nextsite(2,i))=th2(i,1,2)
        htmp(i,nextsite(3,i))=th2(i,1,3)
        htmp(i,nextsite(4,i))=th2(i,1,4)
        htmp(i,i)=dcmplx(ue2(i)/2.0d0-mu2(i),h)
    !dn spin
        htmp(NXY+i,NXY+nextsite(1,i))=th2(i,2,1)
        htmp(NXY+i,NXY+nextsite(2,i))=th2(i,2,2)
        htmp(NXY+i,NXY+nextsite(3,i))=th2(i,2,3)
        htmp(NXY+i,NXY+nextsite(4,i))=th2(i,2,4)
        htmp(NXY+i,NXY+i)=dcmplx(ue2(i)/2.0d0-mu2(i),-h)
    !dn-up flip spin c^+_dn c_up
        htmp(NXY+i,nextsite(1,i))=th3(i,2,1)
        htmp(NXY+i,nextsite(2,i))=th3(i,2,2)
        htmp(NXY+i,nextsite(3,i))=th3(i,2,3)
        htmp(NXY+i,nextsite(4,i))=th3(i,2,4)
        htmp(NXY+i,i)=dcmplx(hx,0.0d0)  !This is ihx*sigma_y
    
    !up-dn flip spin c^+_up c_dn
        htmp(i,NXY+nextsite(1,i))=th3(i,1,1)
        htmp(i,NXY+nextsite(2,i))=th3(i,1,2)
        htmp(i,NXY+nextsite(3,i))=th3(i,1,3)
        htmp(i,NXY+nextsite(4,i))=th3(i,1,4)
        htmp(i,NXY+i)=dcmplx(-hx,0.0d0)   !This is ihx*sigma_y
    
    enddo
    !write(180,*) 'real H'
    !write(180,'(16 f 7.3)') dreal(htmp)
    !write(180,*) 'img H'
    !write(180,'(16 f 7.3)') dimag(htmp)
    temp2=htmp
    call ZGEEV( 'V', 'V', NXYS, htmp, NXYS, eigI, VL, NXYS, VR, NXYS, workee, 2*NXYS,rwork, INFO )
    
    write(160,*) 'eigenvalue'
    write(160,"(10F10.5)") dreal(eigI)
    write(160,"(10F10.5)") dimag(eigI)
    
    
    !write(160,*) 'VL'
    !write(160,'(8 f 10.5)') real(VL)
    !write(160,*) 'VR'
    !write(160,'(8 f 10.5)') real(VR)
    !temp=dcmplx(0.0d0,0.0d0)
    !do i=1,NXYS
    !    temp(i,i)=eigval(i)
    !enddo
    !write(160,*) 'prod'
    !temp=matmul(VR,temp)
    !write(160,'(8 f 10.5)') real(temp)
    !write(160,*) 'prod2'
    !temp=dcmplx(0.0d0,0.0d0)
    !do i=1,NXYS
    !    temp(i,i)=eigval(i)
    !enddo
    !temp=matmul(temp2,VR)
    !write(160,'(8 f 10.5)') real(temp)
    ! find inverse of right vector
    temp=VR
    call ZGETRF(NXYS,NXYS,temp, NXYS, ipvt, info )
    call ZGETRI(NXYS,temp,NXYS,ipvt,work2,4*NXYS,info)
    temp2=dcmplx(0.0d0,0.0d0)
    do i=1,NXYS
        temp2(i,i)=zexp(-dt*eigI(i))
    enddo
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), VR, NXYS, temp2, NXYS, dcmplx(0.0d0,0.0d0), temp3, NXYS)
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), temp3, NXYS, temp, NXYS, dcmplx(0.0d0,0.0d0), temp2, NXYS)
    !temp2=matmul(VR,temp2)
    !temp2=matmul(temp2,temp)
    emh=temp2
    temp2=dcmplx(0.0d0,0.0d0)
    do i=1,NXYS
        temp2(i,i)=zexp(dt*eigI(i))
    enddo
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), VR, NXYS, temp2, NXYS, dcmplx(0.0d0,0.0d0), temp3, NXYS)
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), temp3, NXYS, temp, NXYS, dcmplx(0.0d0,0.0d0), temp2, NXYS)
    !temp2=matmul(VR,temp2)
    !temp2=matmul(temp2,temp)
    eph=temp2
    
    write(161,*) 'eph'
    write(161,'(32 f 20.15)') dreal(eph)
    write(161,'(32 f 20.15)') dimag(eph)
    write(161,*) 'emh'
    write(161,'(32 f 20.15)') dreal(eph)
    write(161,'(32 f 20.15)') dimag(eph)
    
    !Determinant calculation
    temp=emh(:,:)
    deter=det(NXYS,temp)
    write(900+ncpu,*)'in sysdef, emh deter(spin=1)=',deter
    
    temp2=emh(:,:)
    write(900+ncpu,*)'in sysdef, call tred2 tql2 to get eigenvalues of emh'
    
    call ZGEEV( 'V', 'V', NXYS, temp2, NXYS, eigI, VL, NXYS, VR, NXYS, workee, 2*NXYS,rwork, INFO )
    write(900+ncpu,*) 'eigenvalue'
    write(900+ncpu,'(10 f 10.5)') dreal(eigI)
    write(900+ncpu,'(10 f 10.5)') dimag(eigI)
    
    
    !      emh(ispin), eph(ispin) are used in [17-20 bpmult-multbp] to realize B_l^segma times g
    ! ... diagonal elements of exp(-dt*V) and their inverses
    ! ... a diagonal element equals exp(-sigma*J) for up electron spin (ispin=1)
    !     and exp(+sigma*J) for down electron spin (ispin=2)
    ! ... the change in exp of V due to spin flip. the HS spin argument
    !     is the sign of the current HS spin
    !      expdv is used in cnfmake to update g in itth time slice
    
    !   calculate lamda
    do j =1, NXY
        jhs(j) = dexp(+dt*dabs(ue2(j))/2.0d0)
        jhs(j) = dlog(jhs(j) + dsqrt(jhs(j)*jhs(j)-1.0d0))
        jhs(j+NXY) = dexp(+dt*dabs(ue2(j+NXY))/2.0d0)
        jhs(j+NXY) = dlog(jhs(j+NXY) + dsqrt(jhs(j+NXY)*jhs(j+NXY)-1.0d0))
    enddo
    expdv=dcmplx(0.0d0,0.0d0)
    emv=dcmplx(1.0d0,0.0d0)
    epv=dcmplx(1.0d0,0.0d0)
    spinExtraRatio=dcmplx(1.0d0,0.0d0)
    do i=1,NXY
        if (.NOT. negativeUSites(i)) then
            expdv(+1,i)=dcmplx(dexp(-jhs(i))*dexp(-jhs(i))-1.0d0,0.0d0)
            expdv(-1,i)=dcmplx(dexp(jhs(i))*dexp(jhs(i))-1.0d0,0.0d0)
            expdv(-1,i+NXY)=dcmplx(dexp(-jhs(i))*dexp(-jhs(i))-1.0d0,0.0d0)
            expdv(+1,i+NXY)=dcmplx(dexp(jhs(i))*dexp(jhs(i))-1.0d0,0.0d0)
            ! ... define map between electron and HS spins and the v matrix components
            !      emv,epv are used in [17-20 bpmult-multbp] to realize B_l^segma times g
            !     ... for B
            emv(+1,i)=dcmplx(dexp(jhs(i)),0.0d0)
            emv(-1,i)=dcmplx(dexp(-jhs(i)),0.0d0)
            emv(-1,i+NXY)=dcmplx(dexp(jhs(i)),0.0d0)
            emv(+1,i+NXY)=dcmplx(dexp(-jhs(i)),0.0d0)
            !     ... for B^-1
            epv(+1,i)=dcmplx(dexp(-jhs(i)),0.0d0)
            epv(-1,i)=dcmplx(dexp(jhs(i)),0.0d0)
            epv(-1,i+NXY)=dcmplx(dexp(-jhs(i)),0.0d0)
            epv(+1,i+NXY)=dcmplx(dexp(+jhs(i)),0.0d0)
        else
        ! for negative U, no sigma here
            expdv(+1,i)=dcmplx(dexp(-jhs(i))*dexp(-jhs(i))-1.0d0,0.0d0)
            expdv(-1,i)=dcmplx(dexp(jhs(i))*dexp(jhs(i))-1.0d0,0.0d0)
            expdv(+1,i+NXY)=dcmplx(dexp(-jhs(i))*dexp(-jhs(i))-1.0d0,0.0d0)
            expdv(-1,i+NXY)=dcmplx(dexp(jhs(i))*dexp(jhs(i))-1.0d0,0.0d0)
            ! ... define map between electron and HS spins and the v matrix components
            !      emv,epv are used in [17-20 bpmult-multbp] to realize B_l^segma times g
            !     ... for B
            emv(+1,i)=dcmplx(dexp(jhs(i)),0.0d0)
            emv(-1,i)=dcmplx(dexp(-jhs(i)),0.0d0)
            emv(+1,i+NXY)=dcmplx(dexp(jhs(i)),0.0d0)
            emv(-1,i+NXY)=dcmplx(dexp(-jhs(i)),0.0d0)
            !     ... for B^-1
            epv(+1,i)=dcmplx(dexp(-jhs(i)),0.0d0)
            epv(-1,i)=dcmplx(dexp(jhs(i)),0.0d0)
            epv(+1,i+NXY)=dcmplx(dexp(-jhs(i)),0.0d0)
            epv(-1,i+NXY)=dcmplx(dexp(jhs(i)),0.0d0)
            !    for ratio calculation -ve U
            spinExtraRatio(+1,i)=dcmplx(dexp(jhs(i)),0.0d0)
            spinExtraRatio(-1,i)=dcmplx(dexp(-jhs(i)),0.0d0)
            spinExtraRatio(+1,i+NXY)=dcmplx(dexp(jhs(i)),0.0d0)
            spinExtraRatio(-1,i+NXY)=dcmplx(dexp(-jhs(i)),0.0d0)
        endif
    enddo
    write(900+ncpu,*)'jhs(site 1)',jhs(1)
    write(900+ncpu,*)'expdv(site 1)',expdv(:,1)
    write(900+ncpu,*)'emv(site 1)',emv(:,1)
    write(900+ncpu,*)'epv(site 1)',epv(:,1)
    write(900+ncpu,*)'expdv(+1,site1,1) expdv(-1,site1,1)',expdv(+1,1),expdv(-1,1)
    write(900+ncpu,*)'epv(+1,site1,1) epv(+1,site1,2)',epv(+1,1),epv(+1,1+NXY)
    write(900+ncpu,*)'epv(-1,site1,1) epv(-1,site1,2)',epv(-1,1),epv(-1,1+NXY)
    
    return
    end
    
    !3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine cnfinit
    use mbss
    use link,only:makeg
    implicit none
    integer::info,ipvt(NXYS),ispin
    double complex::deterudv,sgn
    double precision::tmpsign
    double complex :: work(NXYS)
    double complex :: temp1(NXYS),temp2(NXYS,NXYS),temp(NXYS,NXYS)
    double precision :: ranGen
    integer::i,site, time
    nsign=-1
    !do while (nsign<0)
        ! ... initialize Hubbard-Stratonvich fields
        spin=0
        do site = 1, NXY
            do time = 1, NT
                if (ranGen(ISEED) > 0.5d0) then
                   spin(site,time) = 1
                else
                   spin(site,time) = -1
                endif
                write(130,*) ranGen(ISEED)
            end do
        end do
    
    !    spin=1
    !    spin(3,1)=-1
    !    write(31,'(16 I 5.2)') spin
        ! ... initialize the Green's function   (l=1,g(l,1:2))
        write(900+ncpu,*)'before makeg(1,1)'
        call makeg(1)
        write(180,*) 'real'
        write(180,NXYSstr) real(g(:,:))
        write(180,*) 'img'
        write(180,NXYSstr) aimag(g(:,:))
        ! ... get the sign of initial determinant.
        deterudv=dcmplx(1.0d0,0.0d0)
        temp=g(:,:)
        call zgetrf( NXYS, NXYS, temp, NXYS, ipvt, info )
        do i=1,NXYS
            deterudv=deterudv*temp(i,i)
        enddo
        sgn = dcmplx(1.0d0,0.0d0)
        do i = 1, NXYS
            if(ipvt(i) /= i) then
                sgn = -sgn
            endif
        enddo
        deterudv = sgn*deterudv
        if(dreal(deterudv) > 0.0d0) then
            nsign=1
        else
            nsign=-1
        endif
    !enddo
    write(900+ncpu,*)'initial success, nsign=',nsign
    write(900+ncpu,*) 'det',deterudv
    !call makeBg
    
    return
    end
    
    !3.1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine readSpin
    use mbss
    implicit none
    integer::NT3,NXY3,NX3,NY3,NLA3,warms3,i
    double precision:: dt3,beta3,ue3,mu3,h3,lamda3
    integer::info,ipvt(NXYS),ispin
    double complex::deterudv,sgn
    double precision::tmpsign
    double complex :: work(NXYS)
    double complex :: temp1(NXYS),temp2(NXYS,NXYS),temp(NXYS,NXYS)
    read(400+ncpu,*) NT3
    if (NT3/=NT) then
        write(*,*) '400: NT not match'
        stop 5555
    endif
    read(400+ncpu,*) NXY3
    if (NXY3/=NXY) then
        write(*,*) '400: NXY not match'
        stop 5555
    endif
    do i=1,NT
        read(400+ncpu,*) spin(1:NXY,i)
    enddo
    read(400+ncpu,*) dt3
    if (abs(dt3-dt)>0.0001d0) then
        write(*,*) '400: dt not match'
        stop 5555
    endif
    read(400+ncpu,*) beta3
    if (abs(beta3-beta)>0.0001d0) then
        write(*,*) '400: beta not match'
        stop 5555
    endif
    read(400+ncpu,*) ue3
    if (abs(ue3-ue)>0.0001d0) then
        write(*,*) '400: ue not match'
        stop 5555
    endif
    read(400+ncpu,*) mu3
    if (abs(mu3-mu)>0.0001d0) then
        write(*,*) '400: mu not match'
        stop 5555
    endif
    read(400+ncpu,*) warms3
    if (warms3/=warms) then
        write(*,*) '400: warms not match'
        stop 5555
    endif
    read(400+ncpu,*) lamda3
    if (abs(lamda3-lamda)>0.0001d0) then
        write(*,*) '400: lamda not match'
        stop 5555
    endif
    read(400+ncpu,*) NX3
    if (NX3/=NX) then
        write(*,*) '400: NX not match'
        stop 5555
    endif
    read(400+ncpu,*) NY3
    if (NY3/=NY) then
        write(*,*) '400: NY not match'
        stop 5555
    endif
    read(400+ncpu,*) NLA3
    if (NLA3/=NLA) then
        write(*,*) '400: NLA not match'
        stop 5555
    endif
    read(400+ncpu,*) h3
    if (abs(h3-h)>0.0001d0) then
        write(*,*) '400: h not match'
        stop 5555
    endif
    
    ! ... initialize the Green's function   (l=1,g(l,1:2))
    write(900+ncpu,*)'before makeg(1,1) in readspin'
        call makeg(1)
    !    write(180,*) 'real'
    !    write(180,'(48 f 10.5)') real(g(:,:))
    !    write(180,*) 'img'
    !    write(180,'(48 f 10.5)') aimag(g(:,:))
        ! ... get the sign of initial determinant.
        deterudv=dcmplx(1.0d0,0.0d0)
        temp=g(:,:)
        call zgetrf( NXYS, NXYS, temp, NXYS, ipvt, info )
        do i=1,NXYS
            deterudv=deterudv*temp(i,i)
        enddo
        sgn = dcmplx(1.0d0,0.0d0)
        do i = 1, NXYS
            if(ipvt(i) /= i) then
                sgn = -sgn
            endif
        enddo
        deterudv = sgn*deterudv
        if(dreal(deterudv) > 0.0d0) then
            nsign=1
        else
            nsign=-1
        endif
    !enddo
    write(900+ncpu,*)'initial success in readspin, nsign=',nsign
    endsubroutine readSpin
    
    
    !4~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine sysequil(nwarms)
    use mbss
    use link,only:cnfmake
    implicit none
    integer::i,j,nwarms
    double precision::acceptc,negs
    acceptc=0.0d0
    negs=0.0d0
    write(900+ncpu,*)'in sysequil,nwarms=',nwarms
    do i = 1, nwarms
        call cnfmake (acceptc,negs)
        if (Mod(i*10,nwarms)==0) then
            j=i*100/nwarms
            write(*,*) 'eqm',j,'% done'
    
        endif
    enddo
    return
    end
    
    !5.0.1
    subroutine initMem(acceptc,negs)
    use mbss
    implicit none
    double precision::negs,acceptc
    mah = 0.0d0
    caverl = 0.0d0
    !initialize accumulator
    negs = 0.0d0
    acceptc = 0.0d0
    mahc=dcmplx(0.0d0,0.0d0)
    caverlc=dcmplx(0.0d0,0.0d0)
    endsubroutine initMem
    
    !5.0.2
    subroutine recMem
    use mbss
    implicit none
    integer:: j
    !accumulate the value to accumulator
    do j = 1, MEAS4
        mah(j) = mah(j) + dble(nsign)*caverl(j)
    enddo
    !write(11,*) mah(1156)
    mah(MEAS4+1)=mah(MEAS4+1)+dble(nsign)
    do j = 1, MEAC3
        mahc(j) = mahc(j) + dcmplx(dble(nsign),0.0d0)*caverlc(j)
    enddo
    endsubroutine recMem
    
    !5.0.3
    subroutine normMem(nswps)
    use mbss
    implicit none
    integer:: i,j,nswps
    ! ... construct measurements from a bin including nswps samples
    do i = 1, MEAS5
        mah(i) = mah(i)/dble(nswps)
    enddo
    ! ... adjust the nswps sample bin measurement with a effective sample number: mah(MEAS4+1)
    do j = 1, MEAS4
    ! mah is the result after nswps(sweeps) sampling in HS field in one Monte Carlo run
    ! due to sign problem the effective sample time is: dble(nswps)*mah(MEAS4+1)
        mah(j) = mah(j)/mah(MEAS4+1)
    enddo
    do i = 1, MEAC3
        mahc(i) = mahc(i)/dble(nswps)
    enddo
    do j = 1, MEAC3
        mahc(j) = mahc(j)/dcmplx(mah(MEAS4+1),0.0d0)
    enddo
    endsubroutine normMem
    
    !5~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine sysmeas (nswps,negs)
    use mbss
    use link,only:cnfmake,cnfmeas,initMem
    use link,only:recMem,FTSF,FTSus,FTSsupBar,sumSus,FTDF,FTGt
    implicit none
    integer:: kpos
    double precision:: xdiff,ydiff,cckq
    integer::i,j,nswps
    double precision::negs,acceptc
    double precision::sns(NXY),csn(NXY,NXY),sdwxy(NXY,NXY),cdw1(NXY,NXY),sdwz1(NXY,NXY),sdwxy1(NXY,NXY)
    double precision::pmdf(NXY,NXY),dsf(NXY,NXY)
    integer::time,q,mpos
    !for vertex
    double complex::gt_vt(NXYS,NXYS)
    double complex::gtinv_vt(NXYS,NXYS)
    double complex::gt2vt(NXYS,NXYS,0:NT), gt1vt(NXYS,NXYS,0:NT)
    !for equal time g
    double complex:: g1(NXYS,NXYS),g2(NXYS,NXYS)

    double precision :: unpair_sw(NXY,NXY), full_sw(NXY,NXY), vertex_sw(NXY,NXY)
    double precision :: unpair_se(NXY,NXY), full_se(NXY,NXY), vertex_se(NXY,NXY)
    double precision :: unpair_d(NXY,NXY),  full_d(NXY,NXY),  vertex_d(NXY,NXY)
    double precision :: unpair_p(NXY,NXY),  full_p(NXY,NXY),  vertex_p(NXY,NXY)
    double precision :: unpair_pup(NXY,NXY),full_pup(NXY,NXY),vertex_pup(NXY,NXY)
    double precision :: p_sdp_bg
    integer :: nn1, nn2, inn, jnn

    call initMem(acceptc,negs)
    
    ! ... construct quantum averages over nswps 'space-time
    !     lattcie sweeps' in one Monte Carlo run.
    ! ... quantum averages comes from important sampling an auxiliary field.
    do i = 1,nswps
    !use cnfmake two times to lessen the correlation between data
        call cnfmake (acceptc,negs)
    !    call cnfmake (acceptc,negs) !change line 696 about times of cnfmake flipping spin
        call cnfmeas
        if (mod(i,100)==1) then
    !    write(30,*) caverl(11)
        endif
        call recMem
        !to know how many percentage done in src
        if (Mod(i*10,nswps)==0) then
            j=i*100/nswps
            write(*,*) j,'% done'
        endif
    end do
    !finish all sweep for measurement in one bind
    mah(MEAS4+2) = acceptc
    mah(MEAS4+3) = negs
    call normMem(nswps)
    ! ... acceptance ratio and negs
    !mah(MEAS4+2) = mah(MEAS4+2)/dble(NT*NXY*2)
    !mah(MEAS4+3) = mah(MEAS4+3)/dble(NT*NXY*2)
    mah(MEAS4+2) = mah(MEAS4+2)/dble(NT*NXY)
    mah(MEAS4+3) = mah(MEAS4+3)/dble(NT*NXY)
    !! Fourier tranform each bin
    ! preparing Fourier Transform
    !     use csn, sdwxy to calculate charge structure factor csf, spin structure factor csf, ssfz, ssfxy,
    !     use cdw1, sdwz1,sdwxy1 to calculate charge susceptbility, spin susceptbility z and xy
    !gt_vt=<c c+>
    !gtinv_vt=<c+ c>
    !<ci ci+>=1 - <ci+ ci>
    gt_vt=dcmplx(0.0d0,0.0d0)
    gtinv_vt=dcmplx(0.0d0,0.0d0)
    gt2vt=dcmplx(0.0d0,0.0d0)
    gt1vt=dcmplx(0.0d0,0.0d0)
    
    sns=mah(MEAS+1+0*NXY:MEAS+1*NXY)
    do i=1,NXY
        do j=1,NXY
            q=(i-1)*NXY+j
            csn(i,j)   =mah(MEAS1+q+0*NXY*NXY)
            sdwxy(i,j) =mah(MEAS1+q+1*NXY*NXY)
            cdw1(i,j)  =mah(MEAS1+q+3*NXY*NXY)
            sdwz1(i,j) =mah(MEAS1+q+4*NXY*NXY)
            sdwxy1(i,j)=mah(MEAS1+q+5*NXY*NXY)
            pmdf(i,j)=mah(MEAS1+q+6*NXY*NXY)
            dsf(i,j)=mah(MEAS1+q+7*NXY*NXY)
        enddo
    enddo
    do i=1,NXYS
        do j=1,NXYS
            q=(i-1)*NXYS+j
            gt_vt(i,j)=mahc(q)
        enddo
    enddo
    gtinv_vt=-gt_vt
    do i=1,NXYS
        gtinv_vt(i,i)=dcmplx(1.0d0,0.0d0)+gtinv_vt(i,i)
    enddo
    
    !unequal time green function for vertex
    do time=0,NT    ! <========== NT+1 time slice totally
        do i=1,NXYS
            do j=1,NXYS
                q=(i-1)*NXYS+j
                gt2vt(i,j,time)=mahc(MEAC2+time*NXYS*NXYS+q)
                gt1vt(i,j,time)=mahc(MEAC2+time*NXYS*NXYS+(NT+1)*NXYS*NXYS+q)
            enddo
        enddo
    enddo
    
    !equal time green function
    do i=1,NXYS
        do j=1,NXYS
            q=(i-1)*NXYS+j
            g2(i,j)=mahc(NXYS*NXYS+q)
            g1(i,j)=mahc(2*NXYS*NXYS+q)
        enddo
    enddo

    do i=1, NXY
        do j=1, NXY
            q = (i-1)*NXY + j
            
            ! 1. Full Correlation
            full_sw(i,j)  = mah(Meas1 + 13*NXY*NXY + q) ! sp    (local s-wave, equal time)
            full_se(i,j)  = mah(Meas1 + 14*NXY*NXY + q) ! sep   (extended s-wave, equal time)
            full_d(i,j)   = mah(Meas1 + 15*NXY*NXY + q) ! dp    (d-wave, equal time)
            full_p(i,j)   = mah(Meas1 + 16*NXY*NXY + q) ! pp    (p-wave, equal time)
            full_pup(i,j) = mah(Meas1 + 17*NXY*NXY + q) ! pupp  (p-up wave, equal time)

            ! Unpaired Background
            unpair_sw(i,j) = 2.0d0 * dreal(g2(j,i)*g2(j+NXY, i+NXY))
            unpair_se(i,j) = 0.0d0
            unpair_d(i,j)  = 0.0d0
            unpair_p(i,j)  = 0.0d0
            unpair_pup(i,j)= 0.0d0

            do nn1=1, NCN
                do nn2=1, NCN
                    inn = nextsite(nn1, i)
                    jnn = nextsite(nn2, j)
                    
                    p_sdp_bg = dreal(g2(j,i)*g2(jnn+NXY, inn+NXY) + g2(jnn, inn)*g2(j+NXY, i+NXY))
                    
                    unpair_se(i,j) = unpair_se(i,j) + sefac(nn1)*sefac(nn2) * p_sdp_bg
                    unpair_d(i,j)  = unpair_d(i,j)  + dfac(nn1)*dfac(nn2)  * p_sdp_bg
                    unpair_p(i,j)  = unpair_p(i,j)  + pfac(nn1)*pfac(nn2)  * p_sdp_bg
                    
                    unpair_pup(i,j) = unpair_pup(i,j) + pupfac(nn1)*pupfac(nn2) * &
                        dreal( g2(j,i)*g2(jnn, inn) - g2(jnn, i)*g2(j, inn) )
                enddo
            enddo

            vertex_sw(i,j)  = full_sw(i,j)  - unpair_sw(i,j)
            vertex_se(i,j)  = full_se(i,j)  - unpair_se(i,j)
            vertex_d(i,j)   = full_d(i,j)   - unpair_d(i,j)
            vertex_p(i,j)   = full_p(i,j)   - unpair_p(i,j)
            vertex_pup(i,j) = full_pup(i,j) - unpair_pup(i,j)

            mah(MEAS_VTX_R + 0*NXY*NXY + q) = vertex_sw(i,j)
            mah(MEAS_VTX_R + 1*NXY*NXY + q) = vertex_se(i,j)
            mah(MEAS_VTX_R + 2*NXY*NXY + q) = vertex_d(i,j)
            mah(MEAS_VTX_R + 3*NXY*NXY + q) = vertex_p(i,j)
            mah(MEAS_VTX_R + 4*NXY*NXY + q) = vertex_pup(i,j)

        enddo
    enddo

    do q = 1, NK
        do i = 1, NXY
            do j = 1, NXY
                kpos = (i-1)*NXY + j
                xdiff = dble(xpos(iposit(ixv(i),iyv(i),1)) - xpos(iposit(ixv(j),iyv(j),1)))
                ydiff = dble(ypos(iposit(ixv(i),iyv(i),1)) - ypos(iposit(ixv(j),iyv(j),1)))
                cckq = dcos(kx(q)*xdiff + ky(q)*ydiff)
                
                mah(MEAS_VTX_K + 0*NK + q) = mah(MEAS_VTX_K + 0*NK + q) + cckq * vertex_sw(i,j)
                mah(MEAS_VTX_K + 1*NK + q) = mah(MEAS_VTX_K + 1*NK + q) + cckq * vertex_se(i,j)
                mah(MEAS_VTX_K + 2*NK + q) = mah(MEAS_VTX_K + 2*NK + q) + cckq * vertex_d(i,j)
                mah(MEAS_VTX_K + 3*NK + q) = mah(MEAS_VTX_K + 3*NK + q) + cckq * vertex_p(i,j)
                mah(MEAS_VTX_K + 4*NK + q) = mah(MEAS_VTX_K + 4*NK + q) + cckq * vertex_pup(i,j)
            enddo
        enddo
    enddo
    mah(MEAS_VTX_K + 1 : MEAS_VTX_K + 5*NK) = mah(MEAS_VTX_K + 1 : MEAS_VTX_K + 5*NK) / dble(NK)

    !call CorSub
    
    !for FT of structure factor and ladder measurement
    call FTSF(sns,csn,sdwxy)
    
    !for FT of Susceptibility
    call FTSus(cdw1,sdwz1,sdwxy1)
    
    !for FT of Susceptibility bar==>for substracting single particle contribution
    call FTSsupBar(gt1vt,gt2vt,g2)
    
    !for FT of Real space Susceptibility
    !call sumSus(cdw1,sdwz1,sdwxy1)
    
    call FTDF(pmdf,dsf)
    
    call FTGt(gt2vt,g1,g2)
    
    call FTPC
    
    call FTCC
    
    !pair correlation
    !---------------pair correlation------------------------!
    !        Meas2+0*Nxy*Nxy*Nxy    <c1u^+ cjd^+ cnd cmu>
    !        Meas2+1*Nxy*Nxy*Nxy    <c1d^+ cju^+ cnu cmd>
    !        Meas2+2*Nxy*Nxy*Nxy    <c1u^+ cju^+ cnu cmu>
    !        Meas2+3*Nxy*Nxy*Nxy    <c1d^+ cjd^+ cnd cmd>
    !---------------pair property------------------!
    !        Meas11+0*Npair	    udPa
    !        Meas11+1*Npair	    uuPa
    !        Meas11+2*Npair	    ddPa
    mpos=NXY*NXY*NXY
    do i=1,Npair
        q=iposit(npx(i,2),npy(i,2),1)*NXY*NXY+iposit(npx(i,3),npy(i,3),1)*NXY+iposit(npx(i,4),npy(i,4),1)
        mah(Meas11+0*Npair+i)=mah(Meas2+q)
        mah(Meas11+1*Npair+i)=mah(Meas2+2*mpos+q)
        mah(Meas11+2*Npair+i)=mah(Meas2+3*mpos+q)
    enddo
    
    ! times beta for susceptibility like observable, for integration
    ! mah(MEAS13+4*NKK+1:MEAS13+7*NKK)=mah(MEAS13+4*NKK+1:MEAS13+7*NKK)*beta
    mah(Meas22+1:Meas22+NXY*NWM)=mah(Meas22+1:Meas22+NXY*NWM)*beta
    
    
    ! mah(21:25)=mah(21:25)*beta
    ! ! beta/NT comes in as d\tau means finish the intergrate for P^bar type variable
    ! mah(MEAS13+8*NKK+1:MEAS13+10*NKK)=mah(MEAS13+8*NKK+1:MEAS13+10*NKK)*beta/dble(NT)
    ! mah(31:35)=mah(31:35)*beta/dble(NT)
    
    mah(MEAS14+1 : MEAS14+5*NKK) = mah(MEAS14+1 : MEAS14+5*NKK) * beta
    ! beta/NT comes in as d\tau means finish the intergrate for P^bar type variable
    mah(MEAS13+8*NKK+1:MEAS13+10*NKK)=mah(MEAS13+8*NKK+1:MEAS13+10*NKK)*beta/dble(NT)
    mah(MEAS14+10*NKK+1 : MEAS14+15*NKK) = mah(MEAS14+10*NKK+1 : MEAS14+15*NKK) * beta/dble(NT)
    return
    endsubroutine sysmeas
    
    !5.1 Fourier TranGensform FT to k-space
    !5.1.1 csf/ssf
    !   FTSF   = fourier transform to structure factor
    !   ns     = n per site
    !   cdwIn  = charge density wave
    !   sdwxyIn= spin density wave XY
    subroutine FTSF(ns,cdwIn,sdwxyIn)
    use mbss
    implicit none
    !input: csn, sdwxy
    !output: csf, ssfz,ssfxy
    !csf(q)=(1/Ns)sum_ij(exp(i*q.(i-j))*cdw(i,j))
    !ssf(q)=(1/Ns)sum_ij(exp(i*q.(i-j))*sdw(i,j))
    double precision::ns(NXY),cdwIn(NXY,NXY),sdwxyIn(NXY,NXY)
    !real-space off-diagonal terms
    integer:: q,i,j,mpos,p1,p2
    double precision::cckq,xdiff,ydiff
    !AFM spin structure factor
    do i=2,NXY
        do j=1,i-1
            if (ilav(i)==ilav(j)) then
                cckq=2.0d0
            else
                cckq=-2.0d0
            endif
            mah(15)=mah(15)+cckq*cdwIn(j,i)
        enddo
    enddo
    do i=1,NXY
        mah(15)=mah(15)+ns(i)-cdwIn(i,i)
    enddo
    mah(15)=mah(15)/dble(NK)
    !charge/ spin structure factor
    do q=1,NK
        do i=1,NXY
            do j=1,NXY
                if (i/=j) then
                mpos=(ilav(i)-1)*NK+(ilav(j)-1)*NLA*NK+q
                xdiff=dble(xpos(iposit(ixv(i),iyv(i),1))-xpos(iposit(ixv(j),iyv(j),1)))
                ydiff=dble(ypos(iposit(ixv(i),iyv(i),1))-ypos(iposit(ixv(j),iyv(j),1)))
    !            cckq=2.0d0*dcos(kx(q)*(xpos(i)-xpos(j))+ky(q)*(ypos(i)-ypos(j)))
                cckq=dcos(kx(q)*xdiff+ky(q)*ydiff)
                if (i<j) then
                    p1=j
                    p2=i
                    else
                    p1=i
                    p2=j
                endif
                mah(MEAS13+mpos) = mah(MEAS13+mpos) + cckq * ( cdwIn(p1,p2) - ns(i)*ns(j) )  !csf(q)  use csn lower triangle
                mah(MEAS13+1*NKK+mpos)=mah(MEAS13+1*NKK+mpos)+cckq*cdwIn(p2,p1)  !ssfz(q)
                mah(MEAS13+2*NKK+mpos)=mah(MEAS13+2*NKK+mpos)+cckq*sdwxyIn(p1,p2) !ssfxy(q)
                endif
            enddo
        enddo
    enddo
    !real-space diagonal terms
    do i=1,NXY
    !lattice AA ==> 1:NK
    !lattice BB ==> NK*3+1:NK*4
        mpos=(ilav(i)-1)*(NKK-NK)
        mah(MEAS13+mpos+1:MEAS13+1*NK+mpos)=mah(MEAS13+mpos+1:MEAS13+1*NK+mpos)+( ns(i) + cdwIn(i,i) - ns(i)*ns(i) ) !csf(q)
        !ns(i) effect eliminated as opposed sign by cdwIn
        mah(MEAS13+mpos+1*NKK+1:MEAS13+mpos+1*NKK+NK)=mah(MEAS13+mpos+1*NKK+1:MEAS13+mpos+1*NKK+NK)+ns(i)-cdwIn(i,i) !ssfz(q)
        mah(MEAS13+mpos+2*NKK+1:MEAS13+mpos+2*NKK+NK)=mah(MEAS13+mpos+2*NKK+1:MEAS13+mpos+2*NKK+NK)+sdwxyIn(i,i) !ssfxy(q)
    end do
    
    mah(MEAS13+1:MEAS13+3*NKK)=mah(MEAS13+1:MEAS13+3*NKK)/dble(NK)
    !spin factor 1/2
    !mah(MEAS13+1*NKK+1:MEAS13+3*NKK)=mah(MEAS13+1*NKK+1:MEAS13+3*NKK)/4.0d0
    return
    end
    
    !5.1.2 csup/ssup
    !   FTSF   = fourier transform of susceptibility
    !   cdwtIn = charge density wave imaginary time average
    !   sdwtIn = spin density wave Z imaginary time average
    !   sdwxytIn=spin density wave XY imaginary time average
    subroutine FTSus(cdwtIn,sdwtIn,sdwxytIn)
    use mbss
    implicit none
    !input: cdw1,sdwz1,sdwxy1
    !output
    !csup(q)=(1/Ns)sum_ij(exp(i*q.(i-j))*integral_tau(<(n_iup+n_idn)_tau*(n_jup+n_jdn)_0>)
    !ssup(q)=(1/Ns)sum_ij(exp(i*q.(i-j))*integral_tau(<S_i(t)*S_j(0)>)
    double precision::cdwtIn(NXY,NXY),sdwtIn(NXY,NXY),sdwxytIn(NXY,NXY)
    integer:: q,i,j,mpos
    double precision::cckq,xdiff,ydiff
    do q=1,NK
        do i=1,NXY
            do j=1,NXY
                mpos=(ilav(i)-1)*NK+(ilav(j)-1)*NLA*NK+q
                xdiff=dble(xpos(iposit(ixv(i),iyv(i),1))-xpos(iposit(ixv(j),iyv(j),1)))
                ydiff=dble(ypos(iposit(ixv(i),iyv(i),1))-ypos(iposit(ixv(j),iyv(j),1)))
    !            cckq=dcos(kx(q)*(xpos(i)-xpos(j))+ky(q)*(ypos(i)-ypos(j)))
                cckq=dcos(kx(q)*xdiff+ky(q)*ydiff)
                mah(MEAS13+4*NKK+mpos)=mah(MEAS13+4*NKK+mpos)+cckq*cdwtIn(i,j) !csup(q)
                mah(MEAS13+5*NKK+mpos)=mah(MEAS13+5*NKK+mpos)+cckq*sdwtIn(i,j) !ssupz(q)
                mah(MEAS13+6*NKK+mpos)=mah(MEAS13+6*NKK+mpos)+cckq*sdwxytIn(i,j) !ssupxy(q)
            enddo
        enddo
    enddo
    mah(MEAS13+4*NKK+1:MEAS13+7*NKK)=mah(MEAS13+4*NKK+1:MEAS13+7*NKK)/dble(NK)
    return
    end
    
    !5.1.3
    !single particle contribution
    subroutine FTSsupBar(gt1vt,gt2vt,g2In)
    use mbss
    implicit none
    double complex::gt2vt(NXYS,NXYS,0:NT), gt1vt(NXYS,NXYS,0:NT), g2In(NXYS,NXYS)
    integer:: it,q,i,j,mpos,nn1,nn2,inn,jnn
    double precision::cckq,xdiff,ydiff,coef
    double precision::val_s,val_se,val_d,val_p,val_pup,p_sdp_bg, p_pup_bg
    
    !input: unequal it green function
    ! ssupz(q)^bar=(1/Ns)sum_ij(exp(i*q.(i-j))*integral_tau ( <G_up(\tau)>_ji*<G~_up(\tau)>_ij+<G_dn(\tau)>_ji*<G~_dn(\tau)>_ij )
    !ssupxy(q)^bar=(1/Ns)sum_ij(exp(i*q.(i-j))*integral_tau ( <G_up(\tau)>_ji*<G~_dn(\tau)>_ij+<G_dn(\tau)>_ji*<G~_up(\tau)>_ij )
    do it=0,NT-1
        do q=1,NK
            do i=1,NXY
                do j=1,NXY
                    mpos=(ilav(i)-1)*NK+(ilav(j)-1)*NLA*NK+q
                    xdiff=dble(xpos(iposit(ixv(i),iyv(i),1))-xpos(iposit(ixv(j),iyv(j),1)))
                    ydiff=dble(ypos(iposit(ixv(i),iyv(i),1))-ypos(iposit(ixv(j),iyv(j),1)))
    !            cckq=dcos(kx(q)*(xpos(i)-xpos(j))+ky(q)*(ypos(i)-ypos(j)))
                    cckq=dcos(kx(q)*xdiff+ky(q)*ydiff)
                    ! after substract gt and 1-gt term, cdw1 is the same as sdwz1
    !                mah(MEAS13+8*NKK+mpos)=mah(MEAS13+8*NKK+mpos)+cckq*dreal(gt2vt(j,i,it)*gt1vt(i,j,it)&
    !                + gt2vt(j+NXY,i+NXY,it)*gt1vt(i+NXY,j+NXY,it)) !ssupz(q)
    !                mah(MEAS13+9*NKK+mpos)=mah(MEAS13+9*NKK+mpos)+cckq*dreal(gt2vt(j,i,it)*gt1vt(i+NXY,j+NXY,it)&
    !                + gt2vt(j+NXY,i+NXY,it)*gt1vt(i,j,it)) !ssupxy(q)
                if (mod(it,4)==0) then
                    ! Spin susceptibility single-particle background
                    mah(MEAS13+8*NKK+mpos)=mah(MEAS13+8*NKK+mpos)+cckq*dreal(gt2vt(j,i,it)*gt1vt(i,j,it)&
                    + gt2vt(j+NXY,i+NXY,it)*gt1vt(i+NXY,j+NXY,it))*28.0d0/45.0d0 !ssupz(q)
                    mah(MEAS13+9*NKK+mpos)=mah(MEAS13+9*NKK+mpos)+cckq*dreal(gt2vt(j,i,it)*gt1vt(i+NXY,j+NXY,it)&
                    + gt2vt(j+NXY,i+NXY,it)*gt1vt(i,j,it))*28.0d0/45.0d0 !ssupxy(q)

                    ! Pairing susceptibility single-particle background (unequal-time)
                    ! Local s-wave
                    mah(MEAS14+10*NKK+mpos) = mah(MEAS14+10*NKK+mpos) + cckq * 2.0d0 * &
                        dreal(gt2vt(j,i,it)*gt2vt(j+NXY,i+NXY,it)) * 28.0d0/45.0d0

                    ! Extended pairing symmetries (unequal-time)
                    do nn1=1,NCN
                        do nn2=1,NCN
                            inn = nextsite(nn1, i)
                            jnn = nextsite(nn2, j)

                            ! Singlet background: <c_i+ c_j> * <c_i- c_j-> + <c_i- c_j-> * <c_i+ c_j+>
                            p_sdp_bg = dreal(gt2vt(j,i,it)*gt2vt(jnn+NXY,inn+NXY,it) + &
                                            gt2vt(jnn,inn,it)*gt2vt(j+NXY,i+NXY,it))

                            mah(MEAS14+11*NKK+mpos) = mah(MEAS14+11*NKK+mpos) + &
                                cckq * sefac(nn1)*sefac(nn2) * p_sdp_bg * 28.0d0/45.0d0
                            mah(MEAS14+12*NKK+mpos) = mah(MEAS14+12*NKK+mpos) + &
                                cckq * dfac(nn1)*dfac(nn2) * p_sdp_bg * 28.0d0/45.0d0
                            mah(MEAS14+13*NKK+mpos) = mah(MEAS14+13*NKK+mpos) + &
                                cckq * pfac(nn1)*pfac(nn2) * p_sdp_bg * 28.0d0/45.0d0

                            ! Triplet background: up-up
                            p_pup_bg = dreal(gt2vt(j,i,it)*gt2vt(jnn,inn,it) - &
                                            gt2vt(jnn,i,it)*gt2vt(j,inn,it))
                            mah(MEAS14+14*NKK+mpos) = mah(MEAS14+14*NKK+mpos) + &
                                cckq * pupfac(nn1)*pupfac(nn2) * p_pup_bg * 28.0d0/45.0d0
                        enddo
                    enddo
                elseif (mod(it,4)==1) then
                    ! Spin susceptibility single-particle background
                    mah(MEAS13+8*NKK+mpos)=mah(MEAS13+8*NKK+mpos)+cckq*dreal(gt2vt(j,i,it)*gt1vt(i,j,it)&
                    + gt2vt(j+NXY,i+NXY,it)*gt1vt(i+NXY,j+NXY,it))*64.0d0/45.0d0 !ssupz(q)
                    mah(MEAS13+9*NKK+mpos)=mah(MEAS13+9*NKK+mpos)+cckq*dreal(gt2vt(j,i,it)*gt1vt(i+NXY,j+NXY,it)&
                    + gt2vt(j+NXY,i+NXY,it)*gt1vt(i,j,it))*64.0d0/45.0d0 !ssupxy(q)

                    ! Pairing susceptibility single-particle background (unequal-time)
                    ! Local s-wave
                    mah(MEAS14+10*NKK+mpos) = mah(MEAS14+10*NKK+mpos) + cckq * 2.0d0 * &
                        dreal(gt2vt(j,i,it)*gt2vt(j+NXY,i+NXY,it)) * 64.0d0/45.0d0

                    ! Extended pairing symmetries (unequal-time)
                    do nn1=1,NCN
                        do nn2=1,NCN
                            inn = nextsite(nn1, i)
                            jnn = nextsite(nn2, j)

                            ! Singlet background: <c_i+ c_j> * <c_i- c_j-> + <c_i- c_j-> * <c_i+ c_j+>
                            p_sdp_bg = dreal(gt2vt(j,i,it)*gt2vt(jnn+NXY,inn+NXY,it) + &
                                            gt2vt(jnn,inn,it)*gt2vt(j+NXY,i+NXY,it))

                            mah(MEAS14+11*NKK+mpos) = mah(MEAS14+11*NKK+mpos) + &
                                cckq * sefac(nn1)*sefac(nn2) * p_sdp_bg * 64.0d0/45.0d0
                            mah(MEAS14+12*NKK+mpos) = mah(MEAS14+12*NKK+mpos) + &
                                cckq * dfac(nn1)*dfac(nn2) * p_sdp_bg * 64.0d0/45.0d0
                            mah(MEAS14+13*NKK+mpos) = mah(MEAS14+13*NKK+mpos) + &
                                cckq * pfac(nn1)*pfac(nn2) * p_sdp_bg * 64.0d0/45.0d0

                            ! Triplet background: up-up
                            p_pup_bg = dreal(gt2vt(j,i,it)*gt2vt(jnn,inn,it) - &
                                            gt2vt(jnn,i,it)*gt2vt(j,inn,it))
                            mah(MEAS14+14*NKK+mpos) = mah(MEAS14+14*NKK+mpos) + &
                                cckq * pupfac(nn1)*pupfac(nn2) * p_pup_bg * 64.0d0/45.0d0
                        enddo
                    enddo
                elseif (mod(it,4)==2) then
                    ! Spin susceptibility single-particle background
                    mah(MEAS13+8*NKK+mpos)=mah(MEAS13+8*NKK+mpos)+cckq*dreal(gt2vt(j,i,it)*gt1vt(i,j,it)&
                    + gt2vt(j+NXY,i+NXY,it)*gt1vt(i+NXY,j+NXY,it))*24.0d0/45.0d0 !ssupz(q)
                    mah(MEAS13+9*NKK+mpos)=mah(MEAS13+9*NKK+mpos)+cckq*dreal(gt2vt(j,i,it)*gt1vt(i+NXY,j+NXY,it)&
                    + gt2vt(j+NXY,i+NXY,it)*gt1vt(i,j,it))*24.0d0/45.0d0 !ssupxy(q)

                    ! Pairing susceptibility single-particle background (unequal-time)
                    ! Local s-wave
                    mah(MEAS14+10*NKK+mpos) = mah(MEAS14+10*NKK+mpos) + cckq * 2.0d0 * &
                        dreal(gt2vt(j,i,it)*gt2vt(j+NXY,i+NXY,it)) * 24.0d0/45.0d0

                    ! Extended pairing symmetries (unequal-time)
                    do nn1=1,NCN
                        do nn2=1,NCN
                            inn = nextsite(nn1, i)
                            jnn = nextsite(nn2, j)

                            ! Singlet background: <c_i+ c_j> * <c_i- c_j-> + <c_i- c_j-> * <c_i+ c_j+>
                            p_sdp_bg = dreal(gt2vt(j,i,it)*gt2vt(jnn+NXY,inn+NXY,it) + &
                                            gt2vt(jnn,inn,it)*gt2vt(j+NXY,i+NXY,it))

                            mah(MEAS14+11*NKK+mpos) = mah(MEAS14+11*NKK+mpos) + &
                                cckq * sefac(nn1)*sefac(nn2) * p_sdp_bg * 24.0d0/45.0d0
                            mah(MEAS14+12*NKK+mpos) = mah(MEAS14+12*NKK+mpos) + &
                                cckq * dfac(nn1)*dfac(nn2) * p_sdp_bg * 24.0d0/45.0d0
                            mah(MEAS14+13*NKK+mpos) = mah(MEAS14+13*NKK+mpos) + &
                                cckq * pfac(nn1)*pfac(nn2) * p_sdp_bg * 24.0d0/45.0d0

                            ! Triplet background: up-up
                            p_pup_bg = dreal(gt2vt(j,i,it)*gt2vt(jnn,inn,it) - &
                                            gt2vt(jnn,i,it)*gt2vt(j,inn,it))
                            mah(MEAS14+14*NKK+mpos) = mah(MEAS14+14*NKK+mpos) + &
                                cckq * pupfac(nn1)*pupfac(nn2) * p_pup_bg * 24.0d0/45.0d0
                        enddo
                    enddo
                elseif (mod(it,4)==3) then
                    ! Spin susceptibility single-particle background
                    mah(MEAS13+8*NKK+mpos)=mah(MEAS13+8*NKK+mpos)+cckq*dreal(gt2vt(j,i,it)*gt1vt(i,j,it)&
                    + gt2vt(j+NXY,i+NXY,it)*gt1vt(i+NXY,j+NXY,it))*64.0d0/45.0d0 !ssupz(q)
                    mah(MEAS13+9*NKK+mpos)=mah(MEAS13+9*NKK+mpos)+cckq*dreal(gt2vt(j,i,it)*gt1vt(i+NXY,j+NXY,it)&
                    + gt2vt(j+NXY,i+NXY,it)*gt1vt(i,j,it))*64.0d0/45.0d0 !ssupxy(q)

                    ! Pairing susceptibility single-particle background (unequal-time)
                    ! Local s-wave
                    mah(MEAS14+10*NKK+mpos) = mah(MEAS14+10*NKK+mpos) + cckq * 2.0d0 * &
                        dreal(gt2vt(j,i,it)*gt2vt(j+NXY,i+NXY,it)) * 64.0d0/45.0d0

                    ! Extended pairing symmetries (unequal-time)
                    do nn1=1,NCN
                        do nn2=1,NCN
                            inn = nextsite(nn1, i)
                            jnn = nextsite(nn2, j)

                            ! Singlet background: <c_i+ c_j> * <c_i- c_j-> + <c_i- c_j-> * <c_i+ c_j+>
                            p_sdp_bg = dreal(gt2vt(j,i,it)*gt2vt(jnn+NXY,inn+NXY,it) + &
                                            gt2vt(jnn,inn,it)*gt2vt(j+NXY,i+NXY,it))

                            mah(MEAS14+11*NKK+mpos) = mah(MEAS14+11*NKK+mpos) + &
                                cckq * sefac(nn1)*sefac(nn2) * p_sdp_bg * 64.0d0/45.0d0
                            mah(MEAS14+12*NKK+mpos) = mah(MEAS14+12*NKK+mpos) + &
                                cckq * dfac(nn1)*dfac(nn2) * p_sdp_bg * 64.0d0/45.0d0
                            mah(MEAS14+13*NKK+mpos) = mah(MEAS14+13*NKK+mpos) + &
                                cckq * pfac(nn1)*pfac(nn2) * p_sdp_bg * 64.0d0/45.0d0

                            ! Triplet background: up-up
                            p_pup_bg = dreal(gt2vt(j,i,it)*gt2vt(jnn,inn,it) - &
                                            gt2vt(jnn,i,it)*gt2vt(j,inn,it))
                            mah(MEAS14+14*NKK+mpos) = mah(MEAS14+14*NKK+mpos) + &
                                cckq * pupfac(nn1)*pupfac(nn2) * p_pup_bg * 64.0d0/45.0d0
                        enddo
                    enddo
                endif
                enddo
            enddo
        enddo
    enddo
    mah(MEAS13+8*NKK+1:MEAS13+10*NKK)=mah(MEAS13+8*NKK+1:MEAS13+10*NKK)/dble(NK)
    
    ! Pairing susceptibility single-particle background (unequal-time)
    ! Initialize arrays    
    mah(MEAS14+10*NKK+1 : MEAS14+15*NKK) = mah(MEAS14+10*NKK+1 : MEAS14+15*NKK)/dble(NK)
    endsubroutine 
    
    !5.1.4
    !   sumSus  = summation of susceptibility in real space
    subroutine sumSus(cdwtIn,sdwtIn,sdwxytIn)
    use mbss
    !input: cdw1,sdwz1,sdwxy1
    !output:
    !xcsup(i)=(1/Ns)sum_j*integral_tau(<(n_iup+n_idn)_tau*(n_jup+n_jdn)_0>)
    !xssup(i)=(1/Ns)sum_j*integral_tau(<S_i(t)*S_j(0)>)
    double precision::cdwtIn(NXY,NXY),sdwtIn(NXY,NXY),sdwxytIn(NXY,NXY)
    !do i=1,NPAIR
    !    do j=1,NXY
    !        n=iposit(ixv(j)+ipairx(i),iyv(j)+ipairy(i),1)
    !        mah(Meas11+21*Npair+i)=mah(Meas11+21*Npair+i)+cdwtIn(j,n) !xcsup(i)
    !        mah(Meas11+22*Npair+i)=mah(Meas11+22*Npair+i)+sdwtIn(j,n) !xssupz(i)
    !        mah(Meas11+23*Npair +i)=mah(Meas11+23*Npair +i)+sdwxytIn(j,n) !xssupxy(i)
    !    enddo
    !enddo
    !mah(Meas11+21*Npair+1:Meas11+24*Npair)=mah(Meas11+21*Npair+1:Meas11+24*Npair)/dble(NXY)
    return
    end
    
    !5.1.5
    subroutine FTDF(pmdf,dsf)
    use mbss
    implicit none
    double precision::pmdf(NXY,NXY),dsf(NXY,NXY)
    integer:: q,i,j,mpos
    double precision::cckq,xdiff,ydiff
    integer:: uppos,dnpos
    !        13 onpair  =>on-site pair corr
    !        14 pairlen => pair correlation length
    !        MEAS13+11*NKK        fdf(q)
    !        MEAS13+12*NKK        fdfUp(q)
    !        MEAS13+13*NKK       fdfDn(q)
    !        MEAS13+14*NKK       pmdf(q)
    !        MEAS13+15*NKK       dsf(q)
    !        Meas1+6*Nxy*Nxy	pmdf(i,j)
    !        Meas1+7*Nxy*Nxy	dsf(i,j)
    do i=1,NXY
        do j=1,NXY
            mah(13)=mah(13)+pmdf(i,j)
            !cckq here= distance betw
            xdiff=abs(xpos(i)-xpos(j))
            ydiff=abs(ypos(i)-ypos(j))
            if (xdiff>dble(NX)/2.0d0) then
                xdiff=dble(NX)-xdiff
            endif
            if (ydiff>dble(NY)/2.0d0) then
                ydiff=dble(NY)-ydiff
            endif
            mah(14)=mah(14)+pmdf(i,j)*(xdiff*xdiff+ydiff*ydiff)
        enddo
    enddo
    mah(13)=mah(13)/dble(NXY)
    mah(14)=mah(14)/dble(NXY)
    
    !        MEAS13+14*NKK       pmdf(q)
    !        MEAS13+15*NKK       dsf(q)
    !        Meas1+6*Nxy*Nxy	pmdf(i,j)
    !        Meas1+7*Nxy*Nxy	dsf(i,j)
    do q=1,NK
        do i=1,NXY
            do j=1,NXY
                mpos=(ilav(i)-1)*NK+(ilav(j)-1)*NLA*NK+q
                xdiff=dble(xpos(iposit(ixv(i),iyv(i),1))-xpos(iposit(ixv(j),iyv(j),1)))
                ydiff=dble(ypos(iposit(ixv(i),iyv(i),1))-ypos(iposit(ixv(j),iyv(j),1)))
        !            cckq=dcos(kx(q)*(xpos(i)-xpos(j))+ky(q)*(ypos(i)-ypos(j)))
                cckq=dcos(kx(q)*xdiff+ky(q)*ydiff)
            !x2 becos it count both (i,j) term and (j,i) term
            if (i==j) then
                mah(MEAS13+14*NKK+mpos)=mah(MEAS13+14*NKK+mpos)+pmdf(i,i)
                mah(MEAS13+15*NKK+mpos)=mah(MEAS13+15*NKK+mpos)+(dsf(i,i)-mah(5)*(mah(MEAS+i+2*NXY)&
                +mah(MEAS+i+2*NXY))+mah(5)*mah(5))
            !            dsf_kSpace(:)=dsf_kSpace(:)+DSF(i,i)
                else
                mah(MEAS13+14*NKK+mpos)=mah(MEAS13+14*NKK+mpos)+cckq*pmdf(i,j)
                mah(MEAS13+15*NKK+mpos)=mah(MEAS13+15*NKK+mpos)+cckq*(dsf(i,j)- mah(5)*(mah(MEAS+i+2*NXY)+mah(MEAS+j+2*NXY))&
                +mah(5)*mah(5))
            !	        dsf_kSpace(q)=dsf_kSpace(q)+cckq*DSF(i,j)
            endif
            enddo
        enddo
    enddo
    mah(MEAS13+14*NKK+1:MEAS13+16*NKK)=mah(MEAS13+14*NKK+1:MEAS13+16*NKK)/dble(NK)
    
    !        MEAS13+11*NKK        fdf(q)
    !        MEAS13+12*NKK        fdfUp(q)
    !        MEAS13+13*NKK       fdfDn(q)
    !        1+1*NXYS*NXYS      gt2(i,j,0) <c^+ c>equal time
    !        1+2*NXYS*NXYS	    gt1(i,j,0) <c c^+>equal time
    
    do q=1,NK
        do i=1,NXY
            do j=1,NXY
                uppos=1*NXYS*NXYS+(i-1)*NXYS+j
                dnpos=1*NXYS*NXYS +NXY*NXYS+(i-1)*NXYS+NXY+j
                mpos=(ilav(i)-1)*NK+(ilav(j)-1)*NLA*NK+q
                xdiff=dble(xpos(iposit(ixv(i),iyv(i),1))-xpos(iposit(ixv(j),iyv(j),1)))
                ydiff=dble(ypos(iposit(ixv(i),iyv(i),1))-ypos(iposit(ixv(j),iyv(j),1)))
        !            cckq=dcos(kx(q)*(xpos(i)-xpos(j))+ky(q)*(ypos(i)-ypos(j)))
                cckq=dcos(kx(q)*xdiff+ky(q)*ydiff)
    
                mah(MEAS13+11*NKK+mpos)=mah(MEAS13+11*NKK+mpos)+dreal(cckq*(mahc(uppos)+mahc(dnpos)))
                mah(MEAS13+12*NKK+mpos)=mah(MEAS13+12*NKK+mpos)+dreal(cckq*(mahc(uppos)))
                mah(MEAS13+13*NKK+mpos)=mah(MEAS13+13*NKK+mpos)+dreal(cckq*(mahc(dnpos)))
            enddo
        enddo
    enddo
    mah(MEAS13+11*NKK+1:MEAS13+14*NKK)=mah(MEAS13+11*NKK+1:MEAS13+14*NKK)/dble(NK)
    endsubroutine FTDF
    
    subroutine FTGt(gt2In,g1In,g2In)
    use mbss
    implicit none
    double complex:: gt2In(NXYS,NXYS,0:NT),g1In(NXYS,NXYS),g2In(NXYS,NXYS)
    integer:: i,j,q,it,mpos
    double complex:: cckq
    
    mahc(MEAC1+1:MEAC1+8*NK)=dcmplx(0.0d0,0.0d0)
    do q=1,NK
        do i=1,NXY
            do j=1,NXY
            cckq=zexp(dcmplx(0.0d0,kx(q)*(xpos(i)-xpos(j))+ky(q)*(ypos(i)-ypos(j))))
    !       <c^+up cup>
            mahc(MEAC1+0*NK+q)=mahc(MEAC1+0*NK+q)+cckq*g2In(i,j)
    !       <c^+dn cdn>
            mahc(MEAC1+1*NK+q)=mahc(MEAC1+1*NK+q)+cckq*g2In(i+NXY,j+NXY)
    !       <c^+dn cup>
            mahc(MEAC1+2*NK+q)=mahc(MEAC1+2*NK+q)+cckq*g2In(i,j+NXY)
    !       <c^+up cdn>
            mahc(MEAC1+3*NK+q)=mahc(MEAC1+3*NK+q)+cckq*g2In(i+NXY,j)
    !       <cup c^+up>
            mahc(MEAC1+4*NK+q)=mahc(MEAC1+4*NK+q)+cckq*g1In(i,j)
    !       <cdn c^+dn>
            mahc(MEAC1+5*NK+q)=mahc(MEAC1+5*NK+q)+cckq*g1In(i+NXY,j+NXY)
    !       <cup c^+dn>
            mahc(MEAC1+6*NK+q)=mahc(MEAC1+6*NK+q)+cckq*g1In(i,j+NXY)
    !       <cdn c^+up>
            mahc(MEAC1+7*NK+q)=mahc(MEAC1+7*NK+q)+cckq*g1In(i+NXY,j)
            enddo
        enddo
    enddo
    mahc(MEAC1+1:MEAC1+8*NK)=mahc(MEAC1+1:MEAC1+8*NK)/dcmplx(dble(NXY),0.0d0)
    
    endsubroutine
    
    subroutine FTPC
    use mbss
    implicit none
    integer:: i,j,q,kpos
    double precision:: xdiff,ydiff,cckq
    
    do q=1,NK
        do i=1,NXY
            do j=1,NXY
                kpos=(i-1)*NXY+j
                xdiff=xpos(i)-xpos(j)
                ydiff=ypos(i)-ypos(j)
                cckq=dcos(kx(q)*xdiff+ky(q)*ydiff)
                
                ! Susceptibility
                mah(MEAS14+0*NKK+q) = mah(MEAS14+0*NKK+q) + cckq*mah(Meas1+8*NXY*NXY+kpos)
                mah(MEAS14+1*NKK+q) = mah(MEAS14+1*NKK+q) + cckq*mah(Meas1+9*NXY*NXY+kpos)
                mah(MEAS14+2*NKK+q) = mah(MEAS14+2*NKK+q) + cckq*mah(Meas1+10*NXY*NXY+kpos)
                mah(MEAS14+3*NKK+q) = mah(MEAS14+3*NKK+q) + cckq*mah(Meas1+11*NXY*NXY+kpos)
                mah(MEAS14+4*NKK+q) = mah(MEAS14+4*NKK+q) + cckq*mah(Meas1+12*NXY*NXY+kpos)
                
                ! Correlation
                mah(MEAS14+5*NKK+q) = mah(MEAS14+5*NKK+q) + cckq*mah(Meas1+13*NXY*NXY+kpos)
                mah(MEAS14+6*NKK+q) = mah(MEAS14+6*NKK+q) + cckq*mah(Meas1+14*NXY*NXY+kpos)
                mah(MEAS14+7*NKK+q) = mah(MEAS14+7*NKK+q) + cckq*mah(Meas1+15*NXY*NXY+kpos)
                mah(MEAS14+8*NKK+q) = mah(MEAS14+8*NKK+q) + cckq*mah(Meas1+16*NXY*NXY+kpos)
                mah(MEAS14+9*NKK+q) = mah(MEAS14+9*NKK+q) + cckq*mah(Meas1+17*NXY*NXY+kpos)
            enddo
        enddo
    enddo
    
    mah(MEAS14+1 : MEAS14+10*NKK) = mah(MEAS14+1 : MEAS14+10*NKK)/dble(NK)
    endsubroutine
    
    ! !q=0 pairing susceptibility (just sum up)
    ! !        21 spsus
    ! !        22 sepsus
    ! !        23 dpsus
    ! !        24 ppsus
    ! !        25 puppsus
    ! !        26 spcor
    ! !        27 sepcor
    ! !        28 dpcor
    ! !        29 ppcor
    ! !        30 puppcor
    ! !        31 spsussp
    ! !        32 sepsussp
    ! !        33 dpsussp
    ! !        34 ppsussp
    ! !        35 puppsussp
    ! !From
    ! !        Meas1+8*Nxy*Nxy    spl
    ! !        Meas1+9*Nxy*Nxy    sepl
    ! !        Meas1+10*Nxy*Nxy    dpl
    ! !        Meas1+11*Nxy*Nxy    ppl
    ! !        Meas1+12*Nxy*Nxy    puppl
    ! !        Meas1+13*Nxy*Nxy    sp
    ! !        Meas1+14*Nxy*Nxy    sep
    ! !        Meas1+15*Nxy*Nxy    dp
    ! !        Meas1+16*Nxy*Nxy    pp
    ! !        Meas1+17*Nxy*Nxy    pupp
    ! do i=1,NXY
    !     do j=1,NXY
    !         q=(i-1)*NXY+j
    !         mah(21)=mah(21)+mah(Meas1+8*Nxy*Nxy+q)
    !         mah(22)=mah(22)+mah(Meas1+9*Nxy*Nxy+q)
    !         mah(23)=mah(23)+mah(Meas1+10*Nxy*Nxy+q)
    !         mah(24)=mah(24)+mah(Meas1+11*Nxy*Nxy+q)
    !         mah(25)=mah(25)+mah(Meas1+12*Nxy*Nxy+q)
    !         mah(26)=mah(26)+mah(Meas1+13*Nxy*Nxy+q)
    !         mah(27)=mah(27)+mah(Meas1+14*Nxy*Nxy+q)
    !         mah(28)=mah(28)+mah(Meas1+15*Nxy*Nxy+q)
    !         mah(29)=mah(29)+mah(Meas1+16*Nxy*Nxy+q)
    !         mah(30)=mah(30)+mah(Meas1+17*Nxy*Nxy+q)
    !     enddo
    ! enddo
    ! mah(21:30)=mah(21:30)/dble(NK)
    ! endsubroutine
    
    subroutine FTCC
    use mbss
    implicit none
    integer::i,j,k,q,w,it,q2
    double precision::coe,coe2,xdiff,ydiff
    !---------------current-current correlation--------------!
    !        Meas21+0*NT*Nxy*Nxy	curcor(i,j,tau)
    !---------------current-current correlation--------------!
    !        Meas22+0*NT*Nxy	    curqw(i,tau)
    mah(Meas22+1:Meas22+NXY*NWM)=0.0d0
    do k=1,NK
        do w=1,NWM
            q=(k-1)*NWM+w
            do i=1,NXY
                do j=1,NXY
                    do it=0,NT-1
    !                    if ((it>0) .or. (i/=j)) then
                    !for omega-t part
                    !wm=2pi*m/beta
                    !tau=it*dt
                    coe=dcos(2*pi*dble(w-1)/beta*dble(it)*dt)
                    !for k-real part
                    xdiff=xpos(i)-xpos(j)
                    ydiff=ypos(i)-ypos(j)
                    coe2=dcos(kx(k)*xdiff+ky(k)*ydiff)/dble(NK)
                    q2=it*NXY*NXY+(i-1)*NXY+j
                    !this is trapezadol rule
    !                mah(Meas22+q)=mah(Meas22+q)+coe*coe2*mah(Meas21+q2)
                    !try simpson rule
    !                if (it==0) then
    !                    mah(Meas22+q)=mah(Meas22+q)+coe*coe2*mah(Meas21+q2)*2.0d0/3.0d0
    !                elseif (mod(it,2)==0) then
    !                    mah(Meas22+q)=mah(Meas22+q)+coe*coe2*mah(Meas21+q2)*2.0d0/3.0d0
    !                elseif (mod(it,2)==1) then
    !                    mah(Meas22+q)=mah(Meas22+q)+coe*coe2*mah(Meas21+q2)*4.0d0/3.0d0
    !                endif
                    !try boole's rule
                    if (mod(it,4)==0) then
                        mah(Meas22+q)=mah(Meas22+q)+coe*coe2*mah(Meas21+q2)*28.0d0/45.0d0
                    elseif (mod(it,4)==1) then
                        mah(Meas22+q)=mah(Meas22+q)+coe*coe2*mah(Meas21+q2)*64.0d0/45.0d0
                    elseif (mod(it,4)==2) then
                        mah(Meas22+q)=mah(Meas22+q)+coe*coe2*mah(Meas21+q2)*24.0d0/45.0d0
                    elseif (mod(it,4)==3) then
                        mah(Meas22+q)=mah(Meas22+q)+coe*coe2*mah(Meas21+q2)*64.0d0/45.0d0
                    endif
    !                endif
                    enddo
                enddo
            enddo
        enddo
    enddo
    
    endsubroutine FTCC
    
    ! 4-operators correlation subtraction of background
    subroutine CorSub
    use mbss
    implicit none
    integer:: j,m,n,k,q,q2
    integer:: mp1,mp2,mp3,mp4
    !condensation fraction
    !        MEAS13+16*NKK       cfsin(q)
    !        MEAS13+17*NKK       cftri(q)
    !        MEAS13+18*NKK       cfup(q)
    !        MEAS13+19*NKK      cfdn(q)
    !write(41,*) 'before'
    !write(41,'(1 f10.5)') mah(Meas2+0*NXY*NXY*NXY+1:Meas2+4*NXY*NXY*NXY)
    !q=NXY*NXY*NXY
    !mpos=(j-1)*NXY*NXY+(m-1)*NXY+n
    ! !        <c1u^+ cju^+ cmu cnu> = <c1u^+ cnu> <cju^+ cnm> - <c1u^+ cmu> <cju^+ cnu>
    !                    upup=dreal(g2t(n,i,time)*g2t(m,j,time)-g2t(m,i,time)*g2t(n,j,time))
    !!        <c1d^+ cjd^+ cnd cmd> = <c1d^+ cmd> <cjd^+ cnd> - <c1d^+ cnd> <cjd^+ cmd>
    !                    dndn=dreal(g2t(n+NXY,i+NXY,time)*g2t(m+NXY,j+NXY,time)-g2t(m+NXY,i+NXY,time)*g2t(n+NXY,j+NXY,time))
    !                    tmp5=dreal(g2t(n,i,time)*g2t(m+NXY,j+NXY,time)-g2t(m+NXY,i,time)*g2t(n,j+NXY,time))
    !                    tmp6=dreal(g2t(n+NXY,i+NXY,time)*g2t(m,j,time)-g2t(m,i+NXY,time)*g2t(n+NXY,j,time))
    !!                   tmp7=dreal(g2t(m,i+NXY,time)*g2t(n+NXY,j,time)-g2t(n+NXY,i+NXY,time)*g2t(m,j,time))
    !!                   tmp8=dreal(g2t(m+NXY,i,time)*g2t(n,j+NXY,time)-g2t(n,i,time)*g2t(m+NXY,j+NXY,time))
    !                    caverl(Meas2+0*q+mpos)=caverl(Meas2+0*q+mpos)+tmp5
    !                    caverl(Meas2+1*q+mpos)=caverl(Meas2+1*q+mpos)+tmp6
    !                    caverl(Meas2+2*q+mpos)=caverl(Meas2+2*q+mpos)+upup
    !                    caverl(Meas2+3*q+mpos)=caverl(Meas2+3*q+mpos)+dndn
    q=NXY*NXY*NXY
    q2=NXYS*NXYS
    !mah(Meas2+0*q+1:Meas2+4*q)=0.0d0
    do j=1,NXY
        do m=1,NXY
            do n=1,NXY
                k=(j-1)*NXY*NXY+(m-1)*NXY+n
    !---------------pair correlation------------------------!
    !        Meas2+0*Nxy*Nxy*Nxy    <c1u^+ cjd^+ cnd cmu>
    !        Meas2+1*Nxy*Nxy*Nxy    <c1d^+ cju^+ cnu cmd>
    !        Meas2+2*Nxy*Nxy*Nxy    <c1u^+ cju^+ cnu cmu>
    !        Meas2+3*Nxy*Nxy*Nxy    <c1d^+ cjd^+ cnd cmd>
    !        <c1u^+ cjd^+ cmd cnu> = <c1u^+ cnu> <cjd^+ cmd> - <c1u^+ cmd> <cjd^+ cnu>
                mp1=(n-1)*NXYS+1                    !<c1u^+ cnu>
                mp2=(m-1+NXY)*NXYS+j+NXY          !<cjd^+ cmd>
                mp3=(m-1+NXY)*NXYS+1              !<c1u^+ cmd>
                mp4=(n-1)*NXYS+j+NXY                !<cjd^+ cnu>
                mah(Meas2+0*q+k)=mah(Meas2+0*q+k)-dreal(mahc(mp1+q2)*mahc(mp2+q2)-mahc(mp3+q2)*mahc(mp4+q2))
    !        <c1d^+ cju^+ cmu cnd> = <c1d^+ cnd> <cju^+ cmu> - <c1d^+ cmu> <cju^+ cnd>
                mp1=(n-1+NXY)*NXYS+1+NXY          !<c1d^+ cnd>
                mp2=(m-1)*NXYS+j                  !<cju^+ cmu>
                mp3=(m-1)*NXYS+1+NXY                !<c1d^+ cmu>
                mp4=(n-1+NXY)*NXYS+j              !<cju^+ cnd>
                mah(Meas2+1*q+k)=mah(Meas2+1*q+k)-dreal(mahc(mp1+q2)*mahc(mp2+q2)-mahc(mp3+q2)*mahc(mp4+q2))
    !        <c1u^+ cju^+ cmu cnu> = <c1u^+ cnu> <cju^+ cmu> - <c1u^+ cmu> <cju^+ cnu>
                mp1=(n-1)*NXYS+1                    !<c1u^+ cnu>
                mp2=(m-1)*NXYS+j                    !<cju^+ cmu>
                mp3=(m-1)*NXYS+1                    !<c1u^+ cmu>
                mp4=(n-1)*NXYS+j                    !<cju^+ cnu>
                mah(Meas2+2*q+k)=mah(Meas2+2*q+k)-dreal(mahc(mp1+q2)*mahc(mp2+q2)-mahc(mp3+q2)*mahc(mp4+q2))
    !        <c1d^+ cjd^+ cmd cnd> = <c1d^+ cnd> <cjd^+ cmd> - <c1d^+ cmd> <cjd^+ cnd>
                mp1=(n-1+NXY)*NXYS+1+NXY                     !<c1d^+ cnd>
                mp2=(m-1+NXY)*NXYS+j+NXY                     !<cjd^+ cmd>
                mp3=(m-1+NXY)*NXYS+1+NXY                     !<c1d^+ cmd>
                mp4=(n-1+NXY)*NXYS+j+NXY                     !<cjd^+ cnd>
                mah(Meas2+3*q+k)=mah(Meas2+3*q+k)-dreal(mahc(mp1+q2)*mahc(mp2+q2)-mahc(mp3+q2)*mahc(mp4+q2))
            enddo
        enddo
    enddo
    !write(41,*) 'after'
    !write(41,'(1 f10.5)') mah(Meas2+0*NXY*NXY*NXY+1:Meas2+4*NXY*NXY*NXY)
    endsubroutine CorSub
    
    !5.5~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !   nmt =   no of measurement
    subroutine measana(nmt)
    use mbss
    implicit none
    double precision::w,ww
    integer::i,nmt
    !   w=1/N
    w = 1.0d0/dble(nmt)
    !   w=1/(N-1)
    ww = 1.0d0/dble(nmt-1)
    !find average and error bar
    do i = 1, MEAS4
        ram(i,fAVG) = w*ram(i,fAVG)
        ram(i,fRMS) = dsqrt(max(0.0d0,ww*(w*ram(i,fRMS)-ram(i,fAVG)*ram(i,fAVG))))
    end do
    do i = 1, MEAC3
        ramc(i,fAVG) = w*ramc(i,fAVG)
        ramc(i,fRMS) = dsqrt(max(0.0d0,ww*(w*dreal(ramc(i,fRMS))-dreal(ramc(i,fAVG)*dconjg(ramc(i,fAVG))))))
    end do
    rsign(1)= w*rsign(1)
    rsign(2)= dsqrt(max(0.0d0,ww*(w*rsign(2)-rsign(1)*rsign(1))))
    return
    end
    
    !6~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine cnfmake(acceptc,negs)
    use mbss
    use link,only:chkg,bpmult,multbm
    implicit none
    integer::i,j,k,l,it,ix,ispin,north
    double precision::acceptc,negs,ranGen
    double complex::part(NSPIN),hold(NSPIN),ratio
    double complex::t(NSPIN),temp1(NXYS,NSPIN),temp2(NXYS,NSPIN)
    double complex:: gtmp
    double precision:: numErr
    double complex:: gErr(NXYS,NXYS),gDiff(NXYS,NXYS)
    north = min(NT,maxrolls)
    ! ... sweep the space-time lattice
    !     ... for each time
    do it=1,NT
        if (mod(it,north) == 1) then
            gErr(:,:)=g(:,:)
            if (it>1) then
    !        call updateBg(it/north)
            endif
    
    
            call makeA(it/north+1)
            call makeipb_inv(it/north+1,g(:,:))
            gDiff(:,:)=g(:,:)-gErr(:,:)
            numErr=0.0d0
            do i=1,NXYS
            do j=1,NXYS
                numErr=numErr+cdabs(gDiff(i,j))
            enddo
            enddo
    !            write(*,*) 'nE',it,numErr/dble(NXYS)/dble(NXYS)
    !        if (nb-(it/north+1)<5) then
    !        call makeg(it)
    !        gDiff(:,:)=g(:,:)-gErr(:,:)
    !        numErr=0.0d0
    !        do i=1,NXYS
    !        do j=1,NXYS
    !            numErr=numErr+cdabs(gDiff(i,j))
    !        enddo
    !        enddo
    !        write(*,*) 'makeg-nE',it,numErr/dble(NXYS)/dble(NXYS)
    !        endif
    !        call makeipb
    !        call matinv (g(:,:))
        endif
    !    write(51,*) it
    !    write(51,NXYSstr) dreal(g(:,:))
    !        ... for each HS site
        do ix=1,NXY
            hold(1)=expdv(spin(ix,it),ix)
            part(1)=dcmplx(1.0d0,0.0d0)+(dcmplx(1.0d0,0.0d0)-g(ix,ix))*hold(1)
    
                temp1(ix+NXY,1)=dcmplx(0.0d0,0.0d0)
                temp2(ix+NXY,1)=dcmplx(0.0d0,0.0d0)
                t(1)=hold(1)/part(1)
                temp1(ix+NXY,1)=g(ix+NXY,ix)
                temp2(ix+NXY,1)=t(1)*g(ix,ix+NXY)
    !           gtmp is g(ix+NXY,ix+NXY)
                gtmp=g(ix+NXY,ix+NXY)+temp1(ix+NXY,1)*temp2(ix+NXY,1)
    
            hold(2)=expdv(spin(ix,it),ix+NXY)
            part(2)=dcmplx(1.0d0,0.0d0)+(dcmplx(1.0d0,0.0d0)-gtmp)*hold(2)
            ratio=part(1)*part(2)*spinExtraRatio(spin(ix,it),ix)*spinExtraRatio(spin(ix,it),ix+NXY)
    !        ratio=zeroc
    !                write(20,*) part(1),'*',part(2),'=r',part(1)*part(2)
            if (ranGen(ISEED) < dabs(dreal(ratio))/(1.0d0+dabs(dreal(ratio)))) then
                if (dreal(ratio) < 0.0d0) nsign = nsign*(-1)
    !                nsign=1
                acceptc = acceptc + 1.0d0
                spin(ix,it)=-spin(ix,it)
                temp1=dcmplx(0.0d0,0.0d0)
                temp2=dcmplx(0.0d0,0.0d0)
                t(1)=hold(1)/part(1)
                do l=1,NXYS
                    temp1(l,1)=g(l,ix)
                    temp2(l,1)=t(1)*g(ix,l)
                enddo
                temp1(ix,1)=temp1(ix,1)-dcmplx(1.0d0,0.0d0)
                do k=1,NXYS
                    do j=1,NXYS
                    g(j,k)=g(j,k)+temp1(j,1)*temp2(k,1)
                    enddo
                enddo
    
                t(2)=hold(2)/part(2)
                do l=1,NXYS
                    temp1(l,2)=g(l,ix+NXY)
                    temp2(l,2)=t(2)*g(ix+NXY,l)
                enddo
                temp1(ix+NXY,2)=temp1(ix+NXY,2)-dcmplx(1.0d0,0.0d0)
                do k=1,NXYS
                    do j=1,NXYS
                    g(j,k)=g(j,k)+temp1(j,2)*temp2(k,2)
                    enddo
                enddo
            endif
    !        write(181,*) nsign
            if(nsign==(-1)) negs=negs+1.0d0
        enddo
    !        ... roll ahead to next step and check g
    !        ... if it=NT, then roll back to g(1).
    
    !        call chkg(it)
    ! just wait makeA in next time slice
    !    if (mod(it,north)>0) then
            call bpmult (it,g(:,:))
            call multbm (it,g(:,:))
    !    endif
    enddo
    !call updateBg(nb)
    return
    end
    
    !7~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine chkg (it)
    use mbss
    implicit none
    ! ... passed and local variables
    integer::it
    !      integer::i,j,k,loc
    !      double precision::gg (NXY,NXY), temp(NXY,NXY), ttemp(NXY,NXY)
    !      double precision::enorm,gnorm,relerr
    if(mod(it,nrolls) == 0) then
    ! ... if nrolls.gt.maxrolls, then compute the new g and measure the
    !        difference between it and the old one
    !      loc=NXY*NXY*(ispin-1)
    !      gg=g(:,:,ispin)
    !      ttemp(:,:)=gg
    !write(21,*) 'renew g'
    call makeg (it)
    !      temp(:,:)=g(:,:,ispin)-ttemp(:,:)
    !      enorm=0.0
    !      gnorm=0.0
    !      do i=1,NXY
    !         enorm=max(enorm,sum(temp(:,i)))
    !         gnorm=max(gnorm,sum(g(:,i,ispin)))
    !      end do
    ! ... if the error is small then wait a little longer before
    !        checking again
    !      relerr = 1.e-6
    !      if (enorm <= relerr*gnorm) then
    !         maxrolls(ispin)=maxrolls(ispin)+1
    !      else
    !     ... check again a little sooner
    !         maxrolls(ispin)=max(maxrolls(ispin)-2,10)
    !      endif
    endif
    return
    end
    
    !8~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine makeg (it)     !the number "l" of B_l multiply before B_M is it
    use mbss
    use link,only:bpmult,makeb,makeipb,matinv
    implicit none
    integer::i,j,it,north,kount
    double complex:: temp(NXYS,NXYS)
    ! ... initialize the product of the B's by defining the first B
    !     to be the identity factorized as U*W*V (u*s*v)
    do i = 1, NXYS
        do j = 1, NXYS
            v(i,j) = 0.0d0
            u(i,j) = 0.0d0
        enddo
        s(i) = 1.0d0
        v(i,i) = 1.0d0
        u(i,i) = 1.0d0
    enddo
    ! ... form the product of a B times the MGS for the previous products
    kount = 0
    !        north =number of times the B's are multiplied before
    !        the product is orthonormalized
    north = min(NT,maxrolls)
    !	  write(*,*)'before loop bpmult, north=',north,'it=',it,'dp1=',dp1
    do i=it,NT
    !    write(180,*) 'real u'
    !write(180,'(8 f 10.5)') real(u)
    !        ... multply B times U and store the result in U
        call bpmult (i,u)
        kount = kount + 1
    !        ... accumulate only maxrolls products at a time before MGSing
        if (mod(kount,north) == 0) then
        !        ... finish the product and leave the result in MGS form
            call makeb
        endif
    
    end do
    !do i=1,NXYS
    !    do j=1,NXYS
    !        temp(i,j)=u(i,j)*s(j)
    !    enddo
    !enddo
    !temp=matmul(temp,v)
    !write(180,*) 'prod udv 2'
    !write(180,'(8 f 10.5)') real(temp)
    
    do i=1,it-1
    call bpmult (i,u)
    kount = kount + 1
        if (mod(kount,north) == 0) then
            call makeb
        endif
    end do
    !do i=1,NXYS
    !    do j=1,NXYS
    !        temp(i,j)=u(i,j)*s(j)
    !    enddo
    !enddo
    !temp=matmul(temp,v)
    !write(180,*) 'prod udv 1'
    !write(180,'(8 f 10.5)') real(temp)
    !write(180,*) 'prod udv 1 img'
    !write(180,'(8 f 10.5)') aimag(temp)
    !write(180,*) 'u aft'
    !write(180,'(16 f 10.5)') real(u)
    !write(180,*) 'u aft'
    !write(180,'(16 f 10.5)') aimag(u)
    
    ! ... leave the result in  form
    if (mod(kount,north) /= 0) call makeb
    ! ... form I+product of B's
    call makeipb
    !(U*Ui)*Di*(Vi*V)
    !temp=matmul(u,ui)
    !do i=1,NXYS
    !    do j=1,NXYS
    !        temp(i,j)=temp(i,j)*s(j)
    !    enddo
    !enddo
    !temp=matmul(temp,vi)
    !temp=matmul(temp,v)
    !write(180,*) 'prod udv 3'
    !write(180,'(8 f 10.5)') real(temp)
    !write(180,*) 'prod udv 3 img'
    !write(180,'(8 f 10.5)') aimag(temp)
    
    ! ... invert the matrix I+B*B*B...  (matrix is stored in the form of u*s*v)
    call matinv (g(:,:))
    !write(180,*) 'g'
    !write(180,'(8 f 10.5)') real(g(:,:))
    !write(180,'(8 f 10.5)') aimag(g(:,:))
    !do i=1,NXYS
    !    do j=1,NXYS
    !        temp(i,j)=u(i,j)*s(j)
    !    enddo
    !enddo
    !temp=matmul(temp,v)
    !write(180,*) 'prod udv 4'
    !write(180,'(8 f 10.5)') real(temp)
    return
    end
    
    !9~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine makeb
    use mbss
    use link,only:udv
    implicit none
    integer::i,j
    double complex::temp(1:NXYS,1:NXYS)
    ! ... U comes in as B*U          u=(b^north)*u
    ! ... form U*W and save the current V in TEMP
    do j = 1, NXYS
        do i = 1, NXYS
            u(i,j) = u(i,j)*s(j)          !u=(b^north)*u*s temp=v
            temp(i,j) = v(i,j)
        end do
    enddo
    !write(180,*) 'real u1'
    !write(180,'(16 f 10.5)') real(u)
    !write(180,*) 'img u'
    !write(180,'(16 f 10.5)') aimag(u)
    call udv(NXYS,u,s,vi)
    ! ... form Vi*TEMP, TEMP contains the old V
    !     ... V and Vi are upper triangular, unit diagonal (mgs case)
    !     ... V and Vi are well condition (udv case)
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), vi, NXYS, temp, NXYS, dcmplx(0.0d0,0.0d0), v, NXYS)
    !v=matmul(vi,temp)
    !    write(180,*) 'real u2'
    !write(180,'(16 f 10.5)') real(u)
    !write(180,*) 'img u'
    !write(180,'(16 f 10.5)') aimag(u)
    return
    end
    
    !9.05~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine udvb
    use mbss
    use link,only:udv
    implicit none
    integer::i,j
    double complex::temp(NXYS,NXYS)
    ! ... U comes in as B*U          u=(b^north)*u
    ! ... form U*W and save the current V in TEMP
    do j = 1, NXYS
        do i = 1, NXYS
            u(i,j) = u(i,j)*s(j)          !u=(b^north)*u*s temp=v
            temp(i,j) = v(i,j)
        end do
    end do
    call udv(NXYS,u,s,vi)
    ! ... form Vi*TEMP, TEMP contains the old V
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), vi, NXYS, temp, NXYS, dcmplx(0.0d0,0.0d0), v, NXYS)
    !v=matmul(vi,temp)
    return
    end
    
    !9.1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine makebt
    use mbss
    use link,only:udv,trans
    implicit none
    integer::i,j
    double complex::temp(NXYS,NXYS)
    ! ... V*[W*U*(B^north)]=temp*[U'']=temp*[(U''^T)^T]=temp*[(U*W*Vi)^T]=[temp*Vi^T]*W*U^T=V_new*W_new*U_new
    ! ... U comes in as U*B          u=u*(b^north)
    ! ... form W*U and save the current V (V now is lower triangle matrix) in TEMP
    do i = 1, NXYS
        do j = 1, NXYS
            u(i,j) = s(i)*u(i,j)          !   u''=W*[u*(b^north)]
            temp(i,j) = v(i,j)            !   temp=v
        end do
    end do
    u=dconjg(trans(NXYS,u))
    !call mgs(NXY,u,s,vi)
    call udv(NXYS,u,s,vi)
    ! ... form TEMP*Vi^T, TEMP contains the old V (V is lower triangle matrix)
    vi=dconjg(trans(NXYS,vi))
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), temp, NXYS, vi, NXYS, dcmplx(0.0d0,0.0d0), v, NXYS)
    !v=matmul(temp,vi)
    !write(*,*)'v',v
    ! ...  U_new=U^T,  W_new=s
    u=dconjg(trans(NXYS,u))
    return
    end
    
    !9.15~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine udvbt
    use mbss
    use link,only:udv,trans
    implicit none
    integer::i,j
    double complex::temp(NXYS,NXYS)
    ! ... V*[W*U*(B^north)]=temp*[U'']=temp*[(U''^T)^T]=temp*[(U*W*Vi)^T]=[temp*Vi^T]*W*U^T=V_new*W_new*U_new
    ! ... U comes in as U*B          u=u*(b^north)
    ! ... form W*U and save the current V (V now is lower triangle matrix) in TEMP
    do i = 1,NXYS
        do j = 1,NXYS
            u(i,j) = s(i)*u(i,j)          !   u''=W*[u*(b^north)]
            temp(i,j) = v(i,j)            !   temp=v
        enddo
    enddo
    u=dconjg(trans(NXYS,u))
    call udv(NXYS,u,s,vi)
    ! ... form TEMP*Vi^T, TEMP contains the old V
    vi=dconjg(trans(NXYS,vi))
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), temp, NXYS, vi, NXYS, dcmplx(0.0d0,0.0d0), v, NXYS)
    !v=matmul(temp,vi)
    ! ...  U_new=U^T,  W_new=s
    u=dconjg(trans(NXYS,u))
    return
    end
    
    !9.2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine saveb8
    use mbss
    use link,only:bpmult,multbp,makeb,makebt
    implicit none
    integer::i,j,it,ispin,north
    ! ... initialize
    north=isudv
    usar=0.0d0;vsar=0.0d0
    usal=0.0d0;vsal=0.0d0
    do i=1,NXYS
        usar(i,i,:)=1.0d0;ssar(i,:)=1.0d0;vsar(i,i,:)=1.0d0
        usal(i,i,:)=1.0d0;ssal(i,:)=1.0d0;vsal(i,i,:)=1.0d0
    enddo
    !.... usar*ssar*vsar=b_l*....*b_1
        do i = 1, NXYS
            do j = 1, NXYS
            v(i,j) = 0.0d0;u(i,j) = 0.0d0
            end do
            s(i) = 1.0d0;v(i,i) = 1.0d0;u(i,i) = 1.0d0
        end do
        do it=1,NT
            call bpmult(it,u)
            if(mod(it,north)==0)then
                call udvb
                usar(:,:,it/north)=u(:,:)
                ssar(:,it/north)=s(:)
                vsar(:,:,it/north)=v(:,:)
            endif
        enddo
        !.... vsal*ssal*usal=b_L*...*b_l
        do i = 1, NXYS
            do j = 1, NXYS
            v(i,j) = 0.0d0;u(i,j) = 0.0d0
            end do
            s(i) = 1.0d0;v(i,i) = 1.0d0;u(i,i) = 1.0d0
        end do
        do it=NT,1,-1
            call multbp(it,u)
            if(mod(it-1,north)==0)then
                call udvbt
                usal(:,:,(it-1)/north)=u(:,:)
                ssal(:,(it-1)/north)=s(:)
                vsal(:,:,(it-1)/north)=v(:,:)
            endif
        enddo
    return
    end
    
    subroutine makeipb_inv(ib,a)
    use mbss
    use link,only:udv,trans
    integer::ib
    integer::i,j,info,ipiv(1:NXYS)
    double complex::temp(NXYS,NXYS),temp2(NXYS,NXYS),tempv(NXYS,NXYS),work(1:4*NXYS),a(NXYS,NXYS)
    if (ib<nb/2) then
    ! I+ ul dl vl vr dr ur
    !=ul(ul^-1 ur^-1 + dl vl vr dr) ur
    !   UL*DL*VL * VR*DR*UR
    !   temp=(DL*VL * VR*DR)
    !VL * VR
    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, vsal(:,:,ib-1), NXYS, vsar(:,:,nb-ib+1), NXYS, zeroc, u(:,:), NXYS)
    !DL*
    do i=1,NXYS
        do j=1,NXYS
            u(i,j)=u(i,j)*ssal(i,ib-1)
        enddo
    enddo
    
    !*DR
    do i=1,NXYS
        do j=1,NXYS
            u(i,j)=u(i,j)*ssar(j,nb-ib+1)
        enddo
    enddo
    
    !ul^-1 ur^-1
    temp=usar(:,:,nb-ib+1)
    !call zgetrf(NXYS,NXYS,temp,NXYS,ipiv,info)
    !call zgetri(NXYS,temp,NXYS,ipiv,work,NXYS,info)
    !call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS,temp(:,:) , NXYS, zeroc, temp2(:,:), NXYS)
    call zgemm('C', 'C', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS,temp(:,:) , NXYS, zeroc, temp2(:,:), NXYS)
    !(ul^-1 ur^-1 + dl vl vr dr)
    u(:,:)=u(:,:)+temp2(:,:)
    !(ul^-1 ur^-1 + dl vl vr dr)=u s v
    call udv(NXYS,u,s,v)
    
    !v'^-1
    temp2=v(:,:)
    call zgetrf(NXYS,NXYS,temp2,NXYS,ipiv,info)
    call zgetri(NXYS,temp2,NXYS,ipiv,work,NXYS,info)
    
    !!vr^-1 v'^-1
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, temp2(:,:) , NXYS, zeroc, tempv(:,:), NXYS)
    !ur^-1 v'^-1
    call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, temp2(:,:) , NXYS, zeroc, tempv(:,:), NXYS)
    
    
    !ul u'
    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS, u(:,:) , NXYS, zeroc, temp(:,:), NXYS)
    
    !vr^-1 v'^-1 d'^-1
    !ur^-1 v'^-1 d'^-1
    do i=1,NXYS
        do j=1,NXYS
            tempv(i,j)=tempv(i,j)/s(j)
        enddo
    enddo
    
    !ur^-1 v'^-1 d'^-1 u'^-1 ul^-1
    call zgemm('N', 'C', NXYS, NXYS, NXYS, onec, tempv(:,:), NXYS, temp(:,:) , NXYS, zeroc, a(:,:), NXYS)
    
    else
    ! I+ ul dl vl vr dr ur
    !=ul(ul^-1 ur^-1 + dl vl vr dr) ur
    !   UL*DL*VL * VR*DR*UR
    !   temp=(DL*VL * VR*DR)
    !VL * VR
    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, vsal(:,:,ib-1), NXYS, vsar(:,:,nb-ib+1), NXYS, zeroc, u(:,:), NXYS)
    !DL*
    do i=1,NXYS
        do j=1,NXYS
            u(i,j)=u(i,j)*ssal(i,ib-1)
        enddo
    enddo
    
    !*DR
    do i=1,NXYS
        do j=1,NXYS
            u(i,j)=u(i,j)*ssar(j,nb-ib+1)
        enddo
    enddo
    
    !ul^-1 ur^-1
    temp=usar(:,:,nb-ib+1)
    !call zgetrf(NXYS,NXYS,temp,NXYS,ipiv,info)
    !call zgetri(NXYS,temp,NXYS,ipiv,work,NXYS,info)
    !call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS,temp(:,:) , NXYS, zeroc, temp2(:,:), NXYS)
    call zgemm('C', 'C', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS,temp(:,:) , NXYS, zeroc, temp2(:,:), NXYS)
    !(ul^-1 ur^-1 + dl vl vr dr)
    u(:,:)=u(:,:)+temp2(:,:)
    u=dconjg(trans(NXYS,u))
    !(ul^-1 ur^-1 + dl vl vr dr)=v^T s^T u^T
    call udv(NXYS,u,s,v)
    
    !v'^-1
    !temp2=v(:,:)
    !v'^T^-1
    temp2=dconjg(trans(NXYS,v(:,:)))
    call zgetrf(NXYS,NXYS,temp2,NXYS,ipiv,info)
    call zgetri(NXYS,temp2,NXYS,ipiv,work,NXYS,info)
    
    !!vr^-1 v'^-1
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, temp2(:,:) , NXYS, zeroc, tempv(:,:), NXYS)
    !ur^-1 v'^T^-1
    !call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, temp2(:,:) , NXYS, zeroc, tempv(:,:), NXYS)
    !v^T^-1 u_l^-1
    call zgemm('N', 'C', NXYS, NXYS, NXYS, onec, temp2(:,:) , NXYS,usal(:,:,ib-1) , NXYS, zeroc, tempv(:,:), NXYS)
    
    !u_r^-1 u'^T^-1 = u_r^-1 u'
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS, u(:,:) , NXYS, zeroc, temp(:,:), NXYS)
    call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, usar(:,:,nb-ib+1), NXYS, u(:,:) , NXYS, zeroc, temp(:,:), NXYS)
    
    !vr^-1 v'^-1 d'^-1
    !ur^-1 u' d'^-1
    do i=1,NXYS
        do j=1,NXYS
    !        tempv(i,j)=tempv(i,j)/s(j)
            temp(i,j)=temp(i,j)/s(j)
        enddo
    enddo
    
    !ur^-1 v'^-1 d'^-1 u'^-1 ul^-1
    !call zgemm('N', 'C', NXYS, NXYS, NXYS, onec, tempv(:,:), NXYS, temp(:,:) , NXYS, zeroc, a(:,:), NXYS)
    !ur^-1 u'^T^-1 d'^-1 v'^T^-1 ul^-1 = ur^-1 u' d'^-1 v'^T^-1 ul^-1
    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, tempv(:,:) , NXYS, zeroc, a(:,:), NXYS)
    
    
    !! (I+ ur^T dr^T vr^T vl^T dl^T ul^T)^T
    !!=(ur^T (ur^T^-1 ul^T^-1 + dr vr^T vl^T dl) ul^T)^T
    !!   UL*DL*VL * VR*DR*UR
    !!   temp=(DL*VL * VR*DR)
    !!VL * VR
    !call zgemm('C', 'C', NXYS, NXYS, NXYS, onec,vsar(:,:,nb-ib+1) , NXYS, vsal(:,:,ib-1), NXYS, zeroc, u(:,:), NXYS)
    !!DL*
    !do i=1,NXYS
    !    do j=1,NXYS
    !        u(i,j)=u(i,j)*ssar(j,nb-ib+1)
    !    enddo
    !enddo
    !
    !!*DR
    !do i=1,NXYS
    !    do j=1,NXYS
    !        u(i,j)=u(i,j)*ssal(i,ib-1)
    !    enddo
    !enddo
    !
    !!ul^-1 ur^-1
    !!temp=usar(:,:,nb-ib+1)
    !!call zgetrf(NXYS,NXYS,temp,NXYS,ipiv,info)
    !!call zgetri(NXYS,temp,NXYS,ipiv,work,NXYS,info)
    !!call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS,temp(:,:) , NXYS, zeroc, temp2(:,:), NXYS)
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, usar(:,:,nb-ib+1), NXYS,usal(:,:,ib-1) , NXYS, zeroc, temp2(:,:), NXYS)
    !!(ur^T^-1 ul^T^-1 + dr vr^T vl^T dl)
    !u(:,:)=u(:,:)+temp2(:,:)
    !!(ur^T^-1 ul^T^-1 + dr vr^T vl^T dl)=u s v
    !call udv(NXYS,u,s,v)
    !
    !!v'^-1
    !temp2=v(:,:)
    !call zgetrf(NXYS,NXYS,temp2,NXYS,ipiv,info)
    !call zgetri(NXYS,temp2,NXYS,ipiv,work,NXYS,info)
    !
    !!!vr^-1 v'^-1
    !!call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, temp2(:,:) , NXYS, zeroc, tempv(:,:), NXYS)
    !!ul v'^-1
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS, temp2(:,:) , NXYS, zeroc, tempv(:,:), NXYS)
    !
    !
    !!ur u'
    !call zgemm('C', 'N', NXYS, NXYS, NXYS, onec,usar(:,:,nb-ib+1) , NXYS, u(:,:) , NXYS, zeroc, temp(:,:), NXYS)
    !
    !!vr^-1 v'^-1 d'^-1
    !!ur^-1 v'^-1 d'^-1
    !do i=1,NXYS
    !    do j=1,NXYS
    !        tempv(i,j)=tempv(i,j)/s(j)
    !    enddo
    !enddo
    !
    !!ur^-1 v'^-1 d'^-1 u'^-1 ul^-1
    !call zgemm('N', 'C', NXYS, NXYS, NXYS, onec, tempv(:,:), NXYS, temp(:,:) , NXYS, zeroc, a(:,:), NXYS)
    !
    !a=dconjg(trans(NXYS,a(:,:)))
    endif
    endsubroutine makeipb_inv
    
    !for unequal time green function record
    
    ! uprod=u^-1 sprod=s^-1 vprod=v^-1
    subroutine makeipb_invRec(ib,a)
    use mbss
    use link,only:udv,trans
    integer::ib
    integer::i,j,info,ipiv(1:NXYS)
    double complex::temp(NXYS,NXYS),temp2(NXYS,NXYS),tempv(NXYS,NXYS),work(1:4*NXYS),a(NXYS,NXYS)
    if (ib<nb/2) then
    ! I+ ul dl vl vr dr ur
    !=ul(ul^-1 ur^-1 + dl vl vr dr) ur
    !   UL*DL*VL * VR*DR*UR
    !   temp=(DL*VL * VR*DR)
    !VL * VR
    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, vsal(:,:,ib-1), NXYS, vsar(:,:,nb-ib+1), NXYS, zeroc, u(:,:), NXYS)
    !DL*
    do i=1,NXYS
        do j=1,NXYS
            u(i,j)=u(i,j)*ssal(i,ib-1)
        enddo
    enddo
    
    !*DR
    do i=1,NXYS
        do j=1,NXYS
            u(i,j)=u(i,j)*ssar(j,nb-ib+1)
        enddo
    enddo
    
    !ul^-1 ur^-1
    temp=usar(:,:,nb-ib+1)
    !call zgetrf(NXYS,NXYS,temp,NXYS,ipiv,info)
    !call zgetri(NXYS,temp,NXYS,ipiv,work,NXYS,info)
    !call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS,temp(:,:) , NXYS, zeroc, temp2(:,:), NXYS)
    call zgemm('C', 'C', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS,temp(:,:) , NXYS, zeroc, temp2(:,:), NXYS)
    !(ul^-1 ur^-1 + dl vl vr dr)
    u(:,:)=u(:,:)+temp2(:,:)
    !(ul^-1 ur^-1 + dl vl vr dr)=u s v
    call udv(NXYS,u,s,v)
    
    !uprod=u^-1 sprod=s^-1 vprod=v^-1
    uprod(:,:,ib-1)=dconjg(trans(NXYS,u))
    
    temp2=v(:,:)
    call zgetrf(NXYS,NXYS,temp2,NXYS,ipiv,info)
    call zgetri(NXYS,temp2,NXYS,ipiv,work,NXYS,info)
    vprod(:,:,ib-1)=temp2
    do i=1,NXYS
        sprod(i,ib-1)=onec/s(i)
    enddo
    
    !!vr^-1 v'^-1
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, temp2(:,:) , NXYS, zeroc, tempv(:,:), NXYS)
    !ur^-1 v'^-1
    call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, temp2(:,:) , NXYS, zeroc, tempv(:,:), NXYS)
    
    
    !ul u'
    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS, u(:,:) , NXYS, zeroc, temp(:,:), NXYS)
    
    !vr^-1 v'^-1 d'^-1
    !ur^-1 v'^-1 d'^-1
    do i=1,NXYS
        do j=1,NXYS
            tempv(i,j)=tempv(i,j)/s(j)
        enddo
    enddo
    
    !ur^-1 v'^-1 d'^-1 u'^-1 ul^-1
    call zgemm('N', 'C', NXYS, NXYS, NXYS, onec, tempv(:,:), NXYS, temp(:,:) , NXYS, zeroc, a(:,:), NXYS)
    
    else
    
    ! I+ ul dl vl vr dr ur
    !=ul(ul^-1 ur^-1 + dl vl vr dr) ur
    !   UL*DL*VL * VR*DR*UR
    !   temp=(DL*VL * VR*DR)
    !VL * VR
    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, vsal(:,:,ib-1), NXYS, vsar(:,:,nb-ib+1), NXYS, zeroc, u(:,:), NXYS)
    !DL*
    do i=1,NXYS
        do j=1,NXYS
            u(i,j)=u(i,j)*ssal(i,ib-1)
        enddo
    enddo
    
    !*DR
    do i=1,NXYS
        do j=1,NXYS
            u(i,j)=u(i,j)*ssar(j,nb-ib+1)
        enddo
    enddo
    
    !ul^-1 ur^-1
    temp=usar(:,:,nb-ib+1)
    !call zgetrf(NXYS,NXYS,temp,NXYS,ipiv,info)
    !call zgetri(NXYS,temp,NXYS,ipiv,work,NXYS,info)
    !call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS,temp(:,:) , NXYS, zeroc, temp2(:,:), NXYS)
    call zgemm('C', 'C', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS,temp(:,:) , NXYS, zeroc, temp2(:,:), NXYS)
    !(ul^-1 ur^-1 + dl vl vr dr)
    u(:,:)=u(:,:)+temp2(:,:)
    u=dconjg(trans(NXYS,u))
    !(ul^-1 ur^-1 + dl vl vr dr)=v^T s^T u^T
    call udv(NXYS,u,s,v)
    
    !uprod=v^T^-1 sprod=s^-1 vprod=u^T^-1=u
    vprod(:,:,ib-1)=u
    
    temp2=dconjg(trans(NXYS,v(:,:)))
    call zgetrf(NXYS,NXYS,temp2,NXYS,ipiv,info)
    call zgetri(NXYS,temp2,NXYS,ipiv,work,NXYS,info)
    uprod(:,:,ib-1)=temp2
    do i=1,NXYS
        sprod(i,ib-1)=onec/s(i)
    enddo
    !!vr^-1 v'^-1
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, temp2(:,:) , NXYS, zeroc, tempv(:,:), NXYS)
    !ur^-1 v'^T^-1
    !call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, temp2(:,:) , NXYS, zeroc, tempv(:,:), NXYS)
    !v^T^-1 u_l^-1
    call zgemm('N', 'C', NXYS, NXYS, NXYS, onec, temp2(:,:) , NXYS,usal(:,:,ib-1) , NXYS, zeroc, tempv(:,:), NXYS)
    
    !u_r^-1 u'^T^-1 = u_r^-1 u'
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS, u(:,:) , NXYS, zeroc, temp(:,:), NXYS)
    call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, usar(:,:,nb-ib+1), NXYS, u(:,:) , NXYS, zeroc, temp(:,:), NXYS)
    
    !vr^-1 v'^-1 d'^-1
    !ur^-1 u' d'^-1
    do i=1,NXYS
        do j=1,NXYS
    !        tempv(i,j)=tempv(i,j)/s(j)
            temp(i,j)=temp(i,j)/s(j)
        enddo
    enddo
    
    !ur^-1 v'^-1 d'^-1 u'^-1 ul^-1
    !call zgemm('N', 'C', NXYS, NXYS, NXYS, onec, tempv(:,:), NXYS, temp(:,:) , NXYS, zeroc, a(:,:), NXYS)
    !ur^-1 u'^T^-1 d'^-1 v'^T^-1 ul^-1 = ur^-1 u' d'^-1 v'^T^-1 ul^-1
    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, temp(:,:), NXYS, tempv(:,:) , NXYS, zeroc, a(:,:), NXYS)
    
    endif
    endsubroutine makeipb_invRec
    
    !10~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine makeipb
    use mbss
    use link,only:udv,trans
    implicit none
    integer::i,info,ipiv(1:NXYS)
    double complex::temp(NXYS,NXYS),work(1:4*NXYS),tempU(NXYS,NXYS)
    ! ... makes I + BB...B
    !        on entry BB...B is UWV.  So I+BB...B=I+UDV
    !        which equals U*(U^T*(I/V)+D)*V
    ! ... The sum is returned in MGS form as (U*Ui)*Di*(Vi*V)
    ! ... form U^T*(I/V)
    temp(:,:)=v(:,:)
    tempU(:,:)=u(:,:)
    !inverse the triangular matrix
    call zgetrf(NXYS,NXYS,temp,NXYS,ipiv,info)
    call zgetri(NXYS,temp,NXYS,ipiv,work,4*NXYS,info)
    !call zgemm('C', 'N', NXYS, NXYS, NXYS, onec, u, NXYS, temp, NXYS, zeroc, ui, NXYS)
    !tempU(:,:)=dconjg(trans(NXYS,tempU(:,:)))
    call zgetrf(NXYS,NXYS,tempU,NXYS,ipiv,info)
    call zgetri(NXYS,tempU,NXYS,ipiv,work,4*NXYS,info)
    !do i=1,NXYS
    !    do j=1,NXYS
    !        if (cdabs(tempU(i,j)-tempU2(i,j))>0.001d0) then
    !        write(*,*) 'not equal'
    !        endif
    !    enddo
    !enddo
    !write(60,*) 'tempv-re'
    !write(60,'(18 f 10.5)') dreal(temp)
    !write(60,*) 'tempv-img'
    !write(60,'(18 f 10.5)') dimag(temp)
    !write(60,*) 'v-re'
    !write(60,'(18 f 10.5)') dreal(v)
    !write(60,*) 'v-img'
    !write(60,'(18 f 10.5)') dimag(v)
    
    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, tempU, NXYS, temp, NXYS, zeroc, ui, NXYS)
    !ui=matmul(transpose(u),temp)        !ui=u^T*(I/V)
    ! ... form U^T*(I/V)+D
    do i = 1, NXYS
        ui(i,i) = ui(i,i) + s(i)
    end do
    ! ... find MGS of this matrix
    !recUDVRat=.true.
    call udv(NXYS,ui,s,vi)
    !recUDVRat=.false.
    return
    end
    
    !11~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine matinv (a)
    use mbss
    use link,only: trans
    implicit none
    integer::info,ipiv(NXYS)
    double complex::work(NXYS)
    ! ... passed and local variables
    integer::i,j
    double complex::a(NXYS,NXYS), tempv(NXYS,NXYS), tempu(NXYS,NXYS)
    !     use after I+BBBB, Vi,Ui obtained
    ! ... inverts (U*Ui)*Di*(Vi*V) ---input from u,ui,s,v,vi
    !     the result is [I/(Vi*V)]*(I/Di)*(U*Ui)^T
    !     ... form U*Ui
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), u, NXYS, ui, NXYS, dcmplx(0.0d0,0.0d0), tempu, NXYS)
    !tempu=matmul(u,ui)
    !     ... form Vi*V
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), vi, NXYS, v, NXYS, dcmplx(0.0d0,0.0d0), tempv, NXYS)
    !tempv=matmul(vi,v)
    !     ... invert the product
    !call triinv (NXY,tempv)
    call zgetrf(NXYS,NXYS,tempv,NXYS,ipiv,info)
    call zgetri(NXYS,tempv,NXYS,ipiv,work,NXYS,info)
    !     ... form (I/HOLD)/D    [I/(Vi*V)]/D
    do j = 1, NXYS
        do i = 1, NXYS
        tempv(i,j) = tempv(i,j)/s(j)
        end do
    end do
    !     ... form TEMPV*TEMPU[T]  [I/(Vi*V)]/D*([u*ui]^T)
    !a=matmul(tempv,trans(tempu))
    ![I/(Vi*V)]/D*(1/[u*ui])
    !call zgetrf(NXYS,NXYS,tempu,NXYS,ipiv,info)
    !call zgetri(NXYS,tempu,NXYS,ipiv,work,NXYS,info)
    call zgemm('N', 'C', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), tempv, NXYS, tempu, NXYS, dcmplx(0.0d0,0.0d0), a, NXYS)
    !a=matmul(tempv,dconjg(trans(NXYS,tempu)))
    return
    end
    
    
    subroutine cnfmeas
    use mbss
    use link,only:gtime
    implicit none
    ! ... local array declarations,
    double precision::eh,emu,eu,elan
    double precision::ehs(NXY),emus(NXY),eus(NXY)
    double precision::sn,sm,sn2,sm2,s2oc,s0oc
    double precision::sns(NXY), sms(NXY), sn2s(NXY), sm2s(NXY),s2ocs(NXY), sus(NXY), sds(NXY)
    double precision::csf(NXY), ssfz(NXY), ssfxy(NXY)
    double precision::csup(NXY),ssupz(NXY),ssupxy(NXY)
    double precision::csn(NXY,NXY),sdwxy(NXY,NXY),opn(NXY)
    double precision::cdw1(NXY,NXY),sdwz1(NXY,NXY),sdwxy1(NXY,NXY)
    double precision::pmdf(NXY,NXY),dsf(NXY,NXY)
    ! ... local varibles
    integer::i,j,ispin,time,q
    double precision::swv(NPAIR),sewv(NPAIR),pxwv(NPAIR),pywv(NPAIR),dxywv(NPAIR),dx2wv(NPAIR),fwv(NPAIR)
    ! one body green function expectation for vertex calculation
    double complex::gt_vt(NXYS,NXYS)
    double precision::Energy,spinx,spiny,spinz,DbEnergyC,DbspinxC,DbspinyC,DbspinzC
    character(len=50):: NTstr
    write(NTstr,*) '(',NT,' i 5.1)'
    NTstr=trim(NTstr)
    !  MEAS=2,MEAS1=2+6*NXY
    ! ... local designation of measured quanties (f band)
    !        sn  -- average occupancy
    !        sm  -- average moment
    !        sn2 -- average squared-occupancy
    !        sm2 -- average squared-moment
    !        s2oc -- average double occupancies
    !        sns -- site occupancies
    !        sms -- site moments
    !        s2ocs -- site double occupancies
    !        sus -- site up number
    !        sds -- site down number
    !        sn2s -- site occupancies squared
    !        sm2s -- site moments squared
    !        csf -- charge structure factor
    !        ssf -- spin structure factor
    !        csup -- charge susceptibility
    !        ssupz -- spin susceptibility z component
    !        ssupxy-- spin susceptibility xy component
    !        csn -- charge density wave, spin density wave, particle number product, at equal time
    !        sdwxy -- spin density wave xy component at equal time
    !        sdwz -- spin density wave z component at equal time
    !        cdw1 -- charge density wave at unequal time
    !        sdwz1 -- spin density wave at unequal time z component
    !        sdwxy1 -- spin density wave at unequal time xy component
    ! ... compute the equal-time Green's function for all time steps
    call gtime(0)
    !write(50,*) 'spin'
    !write(50,NTstr) spin
    !write(50,*) 'gt1(t=0)'
    !write(50,NXYSstr) dreal(gt1(:,:,0))
    !write(50,*) 'gt1(t=1)'
    !write(50,NXYSstr) dreal(gt1(:,:,1))
    !write(50,*) 'gt1(t=2)'
    !write(50,NXYSstr) dreal(gt1(:,:,2))
    !write(50,*) 'gt1(t=4)'
    !write(50,NXYSstr) dreal(gt1(:,:,4))
    !write(50,*) 'gt1(t=5)'
    !write(50,NXYSstr) dreal(gt1(:,:,5))
    !write(50,*) 'gt1(t=20)'
    !write(50,NXYSstr) dreal(gt1(:,:,20))
    !write(50,*) 'gt1(t=30)'
    !write(50,NXYSstr) dreal(gt1(:,:,30))
    !write(50,*) 'gt1(t=40)'
    !write(50,NXYSstr) dreal(gt1(:,:,40))
    !write(50,*) 'gt1(t=50)'
    !write(50,NXYSstr) dreal(gt1(:,:,50))
    !write(50,*) 'gt1(t=55)'
    !write(50,NXYSstr) dreal(gt1(:,:,55))
    !write(50,*) 'gt2(t=0)'
    !write(50,NXYSstr) dreal(gt2(:,:,0))
    !write(50,*) 'gt2(t=1)'
    !write(50,NXYSstr) dreal(gt2(:,:,1))
    !write(50,*) 'gt2(t=4)'
    !write(50,NXYSstr) dreal(gt2(:,:,4))
    !write(50,*) 'gt2(t=5)'
    !write(50,NXYSstr) dreal(gt2(:,:,5))
    !!write(20+ncpu,'(16 f 10.5)') real(gt(:,:,1))
    !write(50,*) 'eqgt(t=2)'
    !write(50,NXYSstr) dreal(gt(:,:,2))
    !write(50,*) 'eqgt(t=4)'
    !write(50,NXYSstr) dreal(gt(:,:,4))
    !write(50,*) 'eqgt(t=5)'
    !write(50,NXYSstr) dreal(gt(:,:,5))
    !write(50,*) 'gt1(t=0)'
    !write(50,'(32 f 10.5)') gt1(:,:,0)
    !write(50,*) 'gt1(t=1)'
    !write(50,'(32 f 10.5)') gt1(:,:,1)
    !write(50,*) 'gt1(t=2)'
    !write(50,'(32 f 10.5)') gt1(:,:,2)
    !write(20+ncpu,'(16 f 10.5)') real(gt(:,:,1))
    eh=0.0d0
    emu=0.0d0
    eu=0.0d0
    ehs=0.0d0
    emus=0.0d0
    eus=0.0d0;
    ! ... initialize sums
    sn=0.0d0
    sm=0.0d0
    sn2=0.0d0
    sm2=0.0d0
    s2oc=0.0d0
    sns=0.0d0
    sms=0.0d0
    sn2s=0.0d0
    sm2s=0.0d0
    s2ocs=0.0d0
    sus=0.0d0
    sds=0.0d0
    s0oc=0.0d0
    csf=0.0d0
    ssfxy=0.0d0
    ssfz=0.0d0
    csup=0.0d0
    ssupz=0.0d0
    ssupxy=0.0d0
    opn=0.0d0
    csn=0.0d0
    sdwxy=0.0d0 !sdwz=0.0
    cdw1=0.0d0
    sdwz1=0.0d0
    sdwxy1=0.0d0
    
    pmdf=0.0d0
    dsf=0.0d0
    
    
    Energy=0.0d0
    spinx=0.0d0
    spiny=0.0d0
    spinz=0.0d0
    DbEnergyC=0.0d0
    DbspinxC=0.0d0
    DbspinyC=0.0d0
    DbspinzC=0.0d0
    
    !  about N : sn, sn2, sns(p), sn2s(p), csf(q), csup(q)
    !        M   sm, sm2, sms(p), sm2s(p), ssf(q), ssup(q)
    !sn=ave(sum[1-g^sgm_ii(l)])
    !equal time g2t <c^+ c>
    g2t=-gt
    
    do time=0,NT
        do i=1,NXYS
            g2t(i,i,time)=dcmplx(1.0d0,0.0d0)+g2t(i,i,time)
        enddo
    enddo
    
    
    !calculate one-body equal green function expectation for P^bar type calculation for recording in caverl
    gt_vt=dcmplx(0.0d0,0.0d0)
    do i=1,NXYS
        do j=1,NXYS
            do time=0,NT-1
                gt_vt(i,j)=gt_vt(i,j)+gt(i,j,time)
            enddo
        enddo
    enddo
    call oneSiteMeas(sn,sm,sn2,sm2,s2oc,s0oc,sns,sms,s2ocs,sus,sds,sn2s,sm2s,eh,emu,eu,elan)
    
    call twoSiteMeas(csn,sdwxy,pmdf,dsf)
    
    call twoSitetMeas(cdw1,sdwz1,sdwxy1)
    
    
    
    !-----------------------------------------!
    !    store in caverl array
    caverl(1)=sn
    caverl(2)=sm
    caverl(3)=sn2
    caverl(4)=sm2
    caverl(5)=s2oc
    caverl(6)=s0oc
    caverl(7)=eu
    caverl(8)=eh
    caverl(9)=emu
    caverl(10)=elan
    caverl(11)=eu+eh+emu+elan
    
    
    
    
    
    caverl(MEAS+1+0*NXY:MEAS+1*NXY)=sns
    caverl(MEAS+1+1*NXY:MEAS+2*NXY)=sms
    caverl(MEAS+1+2*NXY:MEAS+3*NXY)=s2ocs
    caverl(MEAS+1+3*NXY:MEAS+4*NXY)=sds
    caverl(MEAS+1+4*NXY:MEAS+5*NXY)=sn2s
    caverl(MEAS+1+5*NXY:MEAS+6*NXY)=sm2s
    
    
    do i=1,NXY
        do j=1,NXY
            q=(i-1)*NXY+j
            caverl(MEAS1+q+0*NXY*NXY)=csn(i,j)
            caverl(MEAS1+q+1*NXY*NXY)=sdwxy(i,j)
            caverl(MEAS1+q+2*NXY*NXY)=0.0d0
            caverl(MEAS1+q+3*NXY*NXY)=cdw1(i,j)
            caverl(MEAS1+q+4*NXY*NXY)=sdwz1(i,j)
            caverl(MEAS1+q+5*NXY*NXY)=sdwxy1(i,j)
    
            caverl(MEAS1+q+6*NXY*NXY)=pmdf(i,j)
            caverl(MEAS1+q+7*NXY*NXY)=dsf(i,j)
        enddo
    enddo
    ! green function for vertex
    ! store G(\tau) for noninteracting correlation
    ! equal time  ( WITHOUT time index )
    do i=1,NXYS
        do j=1,NXYS
            q=(i-1)*NXYS+j
            caverlc(q)=gt_vt(i,j)
        enddo
    enddo
    !Equal time green function
    caverlc(1*NXYS*NXYS+1:3*NXYS*NXYS)=dcmplx(0.0d0,0.0d0)
    do i=1,NXYS
        do j=1,NXYS
            do time=0,NT-1
                q=(i-1)*NXYS+j
                caverlc(q+1*NXYS*NXYS)=caverlc(q+1*NXYS*NXYS)+g2t(i,j,time)
                caverlc(q+2*NXYS*NXYS)=caverlc(q+2*NXYS*NXYS)+gt(i,j,time)
            enddo
        enddo
    enddo
    
    !Unequal time green function
    do time=0,NT    ! <========== NT+1 time slice totally
        do i=1,NXYS
            do j=1,NXYS
                q=(i-1)*NXYS+j
                caverlc(MEAC2+q+time*NXYS*NXYS)=gt2(i,j,time)
                caverlc(MEAC2+q+time*NXYS*NXYS+(NT+1)*NXYS*NXYS)=gt1(i,j,time)
            enddo
        enddo
    enddo
    
    ! ... time-slice averaged    ! one body unequal time green functions remain time index
    caverl(1:MEAS3)=caverl(1:MEAS3)/dble(NT)
    ! ... site averaged
    caverl(1:6)=caverl(1:6)/dble(NXY)
    caverlc(1:3*NXYS*NXYS)=caverlc(1:3*NXYS*NXYS)/dcmplx(dble(NT),0.0d0)
    
    return
    end
    
    !13---------------------------------------------------------------
    subroutine gtime(tgt)
    use mbss
    implicit none
    ! ... local declaration
    integer::i,j
    integer::tgt,it,ispin,istep,itpi,north
    double complex::h1(NXYS,NXYS), h2(NXYS,NXYS)
    double complex:: temp1(NXYS,NXYS), temp2(NXYS,NXYS)
    ! ... for each spin, compute all equal-time components of the
    !        Green's function
    write(NXYSstrE,*) '(',NXYS,' e 25.15)'
    NXYSstrE=trim(NXYSstrE)
    north = min(NT,maxrolls)
    istep = maxrolls
        do it=1,NT
        if (mod(it,north) == 1) then
    !        ... finish the product and leave the result in MGS form
            call makeANoCh(it/north+1)
    !! ... form I+product of B's
    !        call makeipb
    !
    !! ... invert the matrix I+B*B*B...  (matrix is stored in the form of u*s*v)
    !        call matinv (g(:,:))
            call makeipb_invRec(it/north+1,g(:,:))
        endif
    !    write(702,*) it,'gup-real'
    !    write(702,strfor) dreal(g(:,:,1))
    !    write(702,*) it,'gup-img'
    !    write(702,strfor) dimag(g(:,:,1))
        gt(:,:,it-1) = g(:,:)
        if (mod(it,north)>0) then
         call bpmult (it,g(:,:))
         call multbm (it,g(:,:))
         endif
        enddo
    gt(:,:,NT)=gt(:,:,0)
    g(:,:)=gt(:,:,tgt)
    
    !gt(:,:,:,2)=gt(:,:,:,1)
    !do ispin = 1, NSPIN
    !    north = min(NT,maxrolls)
    !    !        ... determine the numer of times G needs to be made, and
    !    !        ... the number of backward and forward rolls
    !    istep = maxrolls
    !    do it = 1, NT, istep
    !        call makeg(it,ispin)
    !        gt(:,:,it-1,ispin) = g(:,:,ispin)
    !        do i = 1, istep-1
    !            !     ... roll forward
    !            !     in ths case G[m+1]=B[m]*G[m]*(I/B[m])
    !            itpi = it + i
    !            call bpmult (itpi-1,ispin,g(:,:,ispin))
    !            call multbm (itpi-1,ispin,g(:,:,ispin))
    !            gt(:,:,itpi-1,ispin) = g(:,:,ispin)
    !        end do
    !    end do
    !end do
    !g(:,:,1) = gt(:,:,tgt,1)
    !g(:,:,2) = gt(:,:,tgt,2)
    !gt(:,:,NT,1)=gt(:,:,0,1)
    !gt(:,:,NT,2)=gt(:,:,0,2)
    
        gt1(:,:,0)=gt(:,:,0)
        gt2(:,:,0)=-gt(:,:,0)
        !<cj^+(t)ci(t)>
        g2t(:,:,:)=-gt(:,:,:)
        do i =1,NXYS
            gt2(i,i,0)=onec-gt(i,i,0)
            do it=0,NT
            g2t(i,i,it)=onec-gt(i,i,it)
            enddo
        enddo
    
    
    !gt1=Gij=<ci(l1)cj^+(l2)> l1>l2
    !gt2=Gtij=<ci^+(l1)cj(l2)> l1>l2
    
        h1(:,:)=gt1(:,:,0)
        h2(:,:)=gt2(:,:,0)
        !it=0/it=NT is equal time estimated
        do it=1,NT-1
        ! take out it=0,8,16,...
            if (mod(it,north) == 0) then
                call guetRec(it/north+1,gt1(:,:,it),gt2(:,:,it))
                h1(:,:)=gt1(:,:,it)
                h2(:,:)=gt2(:,:,it)
            else
                call bpmult(it,h1(:,:))
                call multbm(it,h2(:,:))
                gt1(:,:,it) = h1(:,:)
                gt2(:,:,it) = h2(:,:)
            endif
        enddo
    
    gt1(:,:,NT) = gt1(:,:,0)
    gt2(:,:,NT) = gt2(:,:,0)
    
    return
    end
    
    ! ur dr vr recorded in nb-ib+1, BN BN-1 ...Bib
    ! ul dl vl recorded in ib-1, Bib-1 Bib-2 ... B1
    !gt1=<c(t)c^+(0)>=((B_l1...B1)^-1+BN BN-1 ... B_{l1+1})^-1)=((vr dr ur) + (ul dl vl)^-1)^-1
    !gt2=<c+(t)c(0)>=((B_l1...B1)+(BN BN-1 ... B_{l1+1})^-1)^-1)=((vr dr ur)^-1 + ul dl vl)^-1
    subroutine guetRec(ib,clcp0,cplc0)
    use mbss
    use link,only:trans
    implicit none
    double complex::clcp0(1:NXYS,1:NXYS),cplc0(1:NXYS,1:NXYS)
    integer:: ib,i,j
    double complex::temp1(1:NXYS,1:NXYS),temp2(1:NXYS,1:NXYS),temp3(1:NXYS,1:NXYS)
    ! gt1=<c(t)c^+(0)>
    !uprod*dl=um
    do i=1,NXYS
        temp1(:,i)=uprod(:,i,ib-1)*ssal(i,ib-1)
    enddo
    !(d'^-1)*(u'^-1*dr)=temp1
    do i=1,NXYS
        temp1(i,:)=sprod(i,ib-1)*temp1(i,:)
    enddo
    
    !temp2=ur^-1 vprod
    call zgemm('C', 'N', NXYS,NXYS,NXYS,onec,usar(:,:,nb-ib+1),NXYS, vprod(:,:,ib-1),NXYS, zeroc,temp2,NXYS)
    
    ! gt1=<c(t)c^+(0)>=temp2*temp1*vl
    call zgemm('N', 'N', NXYS,NXYS,NXYS,onec,temp2,NXYS, temp1,NXYS, zeroc,temp3,NXYS)
    call zgemm('N', 'N', NXYS,NXYS,NXYS,onec,temp3,NXYS, vsal(:,:,ib-1),NXYS, zeroc,clcp0(:,:),NXYS)
    
    
    ! gt2=<c^+(t)c(0)>
    !uprod=um
    !temp1=uprod*(ul^-1)
    call zgemm('N', 'C', NXYS, NXYS, NXYS, dcmplx(1.0,0.0), uprod(:,:,ib-1), &
           NXYS, usal(:,:,ib-1), NXYS, dcmplx(0.0,0.0), temp1, NXYS)
    ! call zgemm('N', 'C', NXYS, NXYS, NXYS, dcmplx(1.0,0.0), uprod(:,:,ib-1), NXYS, usal(:,:,ib-1), NXYS, dcmplx(0.0,0.0), temp1, NXYS)
    
    !temp2=dr*vprod*dprod
    do i=1,NXYS
        temp2(:,i)=vprod(:,i,ib-1)*sprod(i,ib-1)
    enddo
    do i=1,NXYS
        temp2(i,:)=ssar(i,nb-ib+1)*temp2(i,:)
    enddo
    
    ! gt2=<c^+(t)c(0)>=vr*temp2*temp1
    call zgemm('N', 'N', NXYS,NXYS,NXYS,onec,vsar(:,:,nb-ib+1),NXYS, temp2,NXYS, zeroc,temp3,NXYS)
    call zgemm('N', 'N', NXYS,NXYS,NXYS,onec,temp3,NXYS,temp1,NXYS, zeroc,cplc0(:,:),NXYS)
    
    
    !
    !!temp1=vr^-1 v'^-1
    !call zgemm('N', 'N', NXYS,NXYS,NXYS,onec,vrinv(:,:,ib-1),NXYS, vprod(:,:,ib-1),NXYS, zeroc,temp1,NXYS)
    !
    !!temp2=d'^-1 u'^-1 dl
    !do i=1,NXYS
    !    do j=1,NXYS
    !        temp2(i,j)=uprod(i,j,ib-1)*ssal(j,ib-1)
    !    enddo
    !enddo
    !
    !do i=1,NXYS
    !    do j=1,NXYS
    !        temp2(i,j)=temp2(i,j)*sprod(i,ib-1)
    !    enddo
    !enddo
    !
    !!temp3=vr^-1 v'^-1 d'^-1 u'^-1 dl
    !call zgemm('N', 'N', NXYS,NXYS,NXYS,onec,temp1,NXYS, temp2,NXYS, zeroc,temp3,NXYS)
    !
    !! vr^-1 v'^-1 d'^-1 u'^-1 dl vl
    !call zgemm('N', 'N', NXYS,NXYS,NXYS,onec,temp3,NXYS, vsal(:,:,ib-1),NXYS, zeroc,clcp0,NXYS)
    !
    !! gt2=<c^+(t)c(0)>
    !!temp1=u'^-1 ul^-1
    !call zgemm('N', 'C', NXYS,NXYS,NXYS,onec,uprod(:,:,ib-1),NXYS, usal(:,:,ib-1),NXYS, zeroc,temp1,NXYS)
    !
    !!temp2=dr v'^-1 d'^-1
    !do i=1,NXYS
    !    do j=1,NXYS
    !        temp2(i,j)=vprod(i,j,ib-1)*sprod(j,ib-1)
    !    enddo
    !enddo
    !
    !do i=1,NXYS
    !    do j=1,NXYS
    !        temp2(i,j)=temp2(i,j)*ssar(i,nb-ib+1)
    !    enddo
    !enddo
    !
    !! ur dr v'^-1 d'^-1
    !call zgemm('N', 'N', NXYS,NXYS,NXYS,onec,usar(:,:,nb-ib+1),NXYS, temp2,NXYS, zeroc,temp3,NXYS)
    !! ur dr v'^-1 d'^-1 u'^-1 ul^-1
    !call zgemm('N', 'N', NXYS,NXYS,NXYS,onec,temp3,NXYS, temp1,NXYS, zeroc,cplc0,NXYS)
    endsubroutine guetRec
    
    
    !13.1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine guet(im,clcp0,cplc0) !clcp0=GM1=gt1=<c(t)c^+(0)>, cplc0=GM2=gt2=<c^+(t)c(0)>
    use mbss
    use link,only: trans
    implicit none
    double complex::clcp0(1:NXYS,1:NXYS),cplc0(1:NXYS,1:NXYS)
    double complex::ur(1:NXYS,1:NXYS),dr(1:NXYS),vr(1:NXYS,1:NXYS)
    double complex::ul(1:NXYS,1:NXYS),dl(1:NXYS),vl(1:NXYS,1:NXYS)
    double complex::prod2(1:NXYS,1:NXYS),prod3(1:NXYS,1:NXYS)
    double complex::vt(NXYS,NXYS),v1(1:NXYS,1:NXYS),sd(1:NXYS)
    double complex::um(1:NXYS,1:NXYS),v2(1:NXYS,1:NXYS)
    double complex::work(1:NXYS)
    double complex::temp(1:NXYS,1:NXYS)
    integer::i,im,ir,il,ispin,iblock,info,ipiv(1:NXYS)
    ! ... initialize local matrix
    prod2=dcmplx(0.0d0,0.0d0)
    prod3=dcmplx(0.0d0,0.0d0)
    vt=dcmplx(0.0d0,0.0d0)
    v1=dcmplx(0.0d0,0.0d0)
    sd=dcmplx(0.0d0,0.0d0)
    um=dcmplx(0.0d0,0.0d0)
    v2=dcmplx(0.0d0,0.0d0)
    
    !r/l for right/left
    iblock=im/isudv;ir=iblock;il=iblock
        ur=usar(:,:,ir)
        dr=ssar(:,ir)
        vr=vsar(:,:,ir)
        vl=vsal(:,:,il)
        dl=ssal(:,il)
        ul=usal(:,:,il)
        !  clcp0=[(v^-1_r)*(d^-1_r)*(u^-1_r)+v_l*d_l*u_l]^-1
        !       =(u^-1_l)*[(u^-1_r)*(u^-1_l)+d_r*v_r*v_l*d_l]^-1*d_r*v_r
        !       =(u^-1_l)*(case)^-1*d_r*v_r
        !
        !  cplc0=[(u^-1_l)*(d^-1_l)*(v^-1_l)+u_r*d_r*v_r]^-1
        !       =v_l*d_l*[(u^-1_r)*(u^-1_l)+d_r*v_r*v_l*d_l]*(u^-1_r)
        !       =v_l*d_l*(case)^-1**(u^-1_r)
        !
        !  im<beta/2  case=u'd'v'   im>beta/2  case=v'd'u'
        !
        !... calculate value of case
        ! ul*ur=um
        call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), ul, NXYS, ur, NXYS, dcmplx(0.0d0,0.0d0), um, NXYS)
    !    um=matmul(ul,ur)
        ! vr*vl=vt
        call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), vr, NXYS, vl, NXYS, dcmplx(0.0d0,0.0d0), vt, NXYS)
    !    vt=matmul(vr,vl)
        ! dr*(vr*vl)=dr*vt=prod2
        do i=1,NXYS
            prod2(i,:)=dr(i)*vt(i,:)
        enddo
        !(ul*ur)^-1+(dr*vr*vl)*dl=(ut)^t+prod2*dl=prod2
        do i=1,NXYS
            prod2(:,i)=prod2(:,i)*dl(i)
        enddo
    !    temp=um
    !call zgetrf(NXYS,NXYS,temp,NXYS,ipiv,info)
    !call zgetri(NXYS,temp,NXYS,ipiv,work,NXYS,info)
    !    prod2=prod2+temp
        prod2=prod2+dconjg(trans(NXYS,um))
        ! ... case=prod2
        if(im<=NT/2)then
            ! tau<beta/2  udv case
            call udv(NXYS,prod2,sd,v1)
            !clcp0=[(u^-1_l)*(v'^-1)]*[(d'^-1)*(u'^-1)*(d_r)]*v_r
            !u'^-1=um
            um=dconjg(trans(NXYS,prod2))
    !        temp=prod2
    !        call zgetrf(NXYS,NXYS,temp,NXYS,ipiv,info)
    !        call zgetri(NXYS,temp,NXYS,ipiv,work,NXYS,info)
    !        um=temp
            !(u'^-1)*dr=um
            do i=1,NXYS
                um(:,i)=um(:,i)*dr(i)
            enddo
            !(d'^-1)*(u'^-1*dr)=prod3
            do i=1,NXYS
                prod3(i,:)=(1.0d0/sd(i))*um(i,:)
            enddo
            !v'^-1=vt
            vt=v1
            call zgetrf(NXYS,NXYS,vt,NXYS,ipiv,info)
            call zgetri(NXYS,vt,NXYS,ipiv,work,NXYS,info)
            !(ul^-1)*(v'^-1)=ul^T*vt=vt
        call zgemm('C', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), ul, NXYS, vt, NXYS, dcmplx(0.0d0,0.0d0), temp, NXYS)
            vt=temp
    !        vt=matmul(dconjg(trans(NXYS,ul)),vt)
    !        temp=ul
    !        call zgetrf(NXYS,NXYS,temp,NXYS,ipiv,info)
    !        call zgetri(NXYS,temp,NXYS,ipiv,work,NXYS,info)
    !        vt=matmul(temp,vt)
            !cplc0=v_l*[d_l*v'^-1*d'^-1]*[u'^-1*(u^-1_r)]
            !v'^-1=um
            um=v1
            !tri matrix inverse
            call zgetrf(NXYS,NXYS,um,NXYS,ipiv,info)
            call zgetri(NXYS,um,NXYS,ipiv,work,NXYS,info)
            !v'^-1*d'^-1=um
            do i=1,NXYS
                um(:,i)=um(:,i)/sd(i)
            enddo
            !u'^-1*(ur^-1)=(ur*u')^-1=(ur*u')^T=v1
            call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), ur, NXYS, prod2, NXYS, dcmplx(0.0d0,0.0d0), v1, NXYS)
    !        v1=matmul(ur,prod2)
    !        call zgetrf(NXYS,NXYS,v1,NXYS,ipiv,info)
    !        call zgetri(NXYS,v1,NXYS,ipiv,work,NXYS,info)
            v1=dconjg(trans(NXYS,v1))
            !d_l*[(v'^-1)*(d'^-1)]=dl*um=prod2
            do i=1,NXYS
                prod2(i,:)=dl(i)*um(i,:)
            enddo
            else
            ! tau>beta/2  vdu case  initially case is in prod2
            ! trans prod2 to get vdu suited for left side vl,dl,ul
            prod3=dconjg(trans(NXYS,prod2))
            call udv(NXYS,prod3,sd,v2)
            prod2=prod3
            !clcp0=[(u^-1_l)*u']*[(d'^-1)*(v'^T^-1)*(d_r)]*v_r        [v'^T*d'*u'^T]^-1=u'*d'^-1*(v'^T)^-1
            !ul^-1=ul
            ul=dconjg(trans(NXYS,ul))
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), ul, NXYS, prod2, NXYS, dcmplx(0.0d0,0.0d0), vt, NXYS)
            !ul^-1*u'=vt
    !        vt=matmul(ul,prod2)
            !v'^T^-1=um^-1=((um^T)^-1)^T=um
            um=v2
            !tri matrix inverse
            call zgetrf(NXYS,NXYS,um,NXYS,ipiv,info)
            call zgetri(NXYS,um,NXYS,ipiv,work,NXYS,info)
            um=dconjg(trans(NXYS,um))
            !v'^T^-1*ur^-1=v'^-1^T*ur^-1=um*ur^-1=v1
    call zgemm('N', 'C', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), um, NXYS, ur, NXYS, dcmplx(0.0d0,0.0d0), v1, NXYS)
    !        v1=matmul(um,dconjg(trans(NXYS,ur)))
            !v'^T^-1*dr=um*dr=prod3
            do i=1,NXYS
                prod3(:,i)=um(:,i)*dr(i)
            enddo
            !d'^-1*(v'^T^-1*dr)=d'^-1*prod3=prod3
            do i=1,NXYS
                prod3(i,:)=(1.0d0/sd(i))*prod3(i,:)
            enddo
            !cplc0=v_l*(d_l*u'*d'^-1)*(v'^T^-1*u^-1_r)      [v'^T*d'*u'^T]^-1=u'*d'^-1*(v'^T)^-1
            !v'^T^-1*ur^-1=v1
            !u'*d'^-1=prod2*(1/s)=prod2
            do i=1,NXYS
                prod2(:,i)=prod2(:,i)/sd(i)
            enddo
            !dl*(u'*d'^-1)=dl*prod2=prod2
            do i=1,NXYS
                prod2(i,:)=dl(i)*prod2(i,:)
            enddo
        endif
        ! gt1=clcp0= uv dud v=vt*prod3*vr
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), vt, NXYS, prod3, NXYS, dcmplx(0.0d0,0.0d0), um, NXYS)
    !    um=matmul(vt,prod3)
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), um, NXYS, vr, NXYS, dcmplx(0.0d0,0.0d0), prod3, NXYS)
    !    prod3=matmul(um,vr)
        ! gt2=cplc0= v dud vu=vl*prod2*v1
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), vl, NXYS, prod2, NXYS, dcmplx(0.0d0,0.0d0), um, NXYS)
    !    um=matmul(vl,prod2)
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), um, NXYS, v1, NXYS, dcmplx(0.0d0,0.0d0), prod2, NXYS)
    !    prod2=matmul(um,v1)
        clcp0(:,:)=prod3
        cplc0(:,:)=prod2
    return
    end
    
    
    !13.3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !   oneSiteMeas     = two operators terms
    !   with s          = per site variable
    !   en              = energy
    subroutine oneSiteMeas(nOut,mOut,n2Out,m2Out,upDnOut,vacOut,nsOut,msOut,upDnsOut, &
    upsOut,dnsOut,n2sOut,m2sOut,enH,enMu,enU,enLan)
    use mbss
    implicit none
    !Input: Unequal time green function, directly use global gt
    !        ETGreenFunction: equal time green function(if can be used)
    !Output:
    !        sn  -- average occupancy
    !        sm  -- average moment
    !        sn2 -- average squared-occupancy
    !        sm2 -- average squared-moment
    !        s2oc -- average double occupancies
    !       s0oc-average vacancy
    !        sns -- site occupancies
    !        sms -- site moments
    !        s2ocs -- site double occupancies
    !        sus -- site up number
    !        sds -- site down number
    !        sn2s -- site occupancies squared
    !        sm2s -- site moments squared
    !define local variable
    double precision::up,dn,updnupdn
    !Ahout energy
    double precision::enH,enMu,enU,enLan
    double precision::enHSites(NXY),enMuSites(NXY),enUSites(NXY),enLanSites(NXY)
    !about spin n/m
    double precision::nOut,mOut,n2Out,m2Out,upDnOut,vacOut
    double precision::nsOut(NXY),msOut(NXY),upDnsOut(NXY),upsOut(NXY),dnsOut(NXY),n2sOut(NXY),m2sOut(NXY)
    integer:: site,time,j,ispin,spos
    nOut=0.0d0
    mOut=0.0d0
    n2Out=0.0d0
    m2Out=0.0d0
    upDnOut=0.0d0
    vacOut=0.0d0
    nsOut=0.0d0
    msOut=0.0d0
    upDnsOut=0.0d0
    upsOut=0.0d0
    dnsOut=0.0d0
    n2sOut=0.0d0
    m2sOut=0.0d0
    enH=0.0d0
    enMu=0.0d0
    enU=0.0d0
    enHSites=0.0d0
    enMuSites=0.0d0
    enUSites=0.0d0
    enLan=0.0d0
    enLanSites=0.0d0
    ! ... equal-spin correlations sums,equal-time averages,equal-site sum,
    do site = 1, NXY
        do time=0,NT-1
    !    write(180,*) 'real',time
    !    write(180,'(8 f 10.5)') real(gt(:,:,time))
    !    write(180,*) 'img',time
    !        write(180,'(8 f 10.5)') aimag(gt(:,:,time))
            up = 1.0d0 - dreal(gt(site,site,time))
            dn = 1.0d0 - dreal(gt(site+NXY,site+NXY,time))
            updnupdn=dreal((1.0d0 - gt(site,site,time))*(1.0d0 - gt(site+NXY,site+NXY,time)))&
            -dreal(gt(site,site+NXY,time)*gt(site+NXY,site,time))
            nsOut(site) = nsOut(site) + (up + dn) !site occupancies, averages over time
            msOut(site) = msOut(site) + (up - dn)
            upsOut(site) =upsOut(site) + (up)
            dnsOut(site) = dnsOut(site) + (dn)
            upDnsOut(site)=upDnsOut(site) + updnupdn
            n2sOut(site)=n2sOut(site) + (up + dn + 2.0*updnupdn)
            m2sOut(site)=m2sOut(site) + (up + dn - 2.0*updnupdn)
            !calculate energy (two-operator terms)
            do j=1,NCN
            !KE= sum t*<ci_up^+ cj_up> + t*<ci_dn^+ cj_dn>
                enHSites(site)=enHSites(site)-dreal(gt(nextsite(j,site),site,time)*th2(site,1,j))
                enHSites(site)=enHSites(site)-dreal(gt(nextsite(j,site)+NXY,site+NXY,time)*th2(site,2,j))
            enddo
    
    !!dn-up flip spin c^+_dn c_up
    !    htmp(NXY+i,nextsite(1,i))=th3(i,2,1)
    !    htmp(NXY+i,nextsite(2,i))=th3(i,2,2)
    !    htmp(NXY+i,nextsite(3,i))=th3(i,2,3)
    !    htmp(NXY+i,nextsite(4,i))=th3(i,2,4)
    !!up-dn flip spin c^+_up c_dn
    !    htmp(i,NXY+nextsite(1,i))=th3(i,1,1)
    !    htmp(i,NXY+nextsite(2,i))=th3(i,1,2)
    !    htmp(i,NXY+nextsite(3,i))=th3(i,1,3)
    !    htmp(i,NXY+nextsite(4,i))=th3(i,1,4)
            do j=1,NCN
                !<c^+_dn c_up>=-<c_up c^+_dn>
                enLanSites(site)=enLanSites(site)-dreal(gt(nextsite(j,site),NXY+site,time)*th3(site,2,j))
                !<c^+_up c_dn>=-<c_dn c^+_up>
                enLanSites(site)=enLanSites(site)-dreal(gt(NXY+nextsite(j,site),site,time)*th3(site,1,j))
            enddo
    !        enMuSites(site)=enMuSites(site)+ (up+dn)*(-mu)+(up-dn)*h   !This is for h*sigma_z
             enMuSites(site)=enMuSites(site)+ (up+dn)*(-mu)-h*dreal(gt(NXY+site,site,time)+gt(site,site+NXY,time))  !This is for h*sigma_x
            enUSites(site) =enUSites(site) + updnupdn
        enddo
    enddo
    do site=1,NXY
        nOut=nOut+nsOut(site)
        mOut=mOut+msOut(site)
        n2Out=n2Out+n2sOut(site)
        m2Out=m2Out+m2sOut(site)
        upDnOut=upDnOut+upDnsOut(site)
        vacOut=vacOut+upDnsOut(site)-nsOut(site)+dble(NT)
        enH=enH+enHSites(site)
        enMu=enMu+enMuSites(site)
        enU=enU+enUSites(site)
        enLan=enLan+enLanSites(site)
    enddo
    enH=enH
    enMu=enMu
    enU=enU*ue
    return
    end
    
    !13.4~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine twoSiteMeas(cdwIn,sdwxyIn,pmdfIn,dsfIn)
    use mbss
    implicit none
    !Input: Unequal time green function, directly use global gt
    !        ETGreenFunction: equal time green function(if can be used)
    !Output:
    !real space,only two obtained
    !        csn -- charge density wave, spin density wave, particle number product, at equal time
    !        sdwxy -- spin density wave xy component at equal time
    !(not obtained here)
    !k space, after collecting all data, doing fourier tranGensform
    !        sdwz -- spin density wave z component at equal time (from csn)
    !        csf -- charge structure factor (from csn)
    !        ssf -- spin structure factor (xy part from sdwxy)
    !Details Method
    !cdw(i,j)=<(n_iup+n_idn)(n_jup+n_jdn)>
    ! i!= j
    !        = (1-gup_ii)(1-gup_jj)-gup_ji*gup_ij+(1-gup_ii)(1-gdn_jj)
    !         +(1-gdn_ii)(1-gup_jj)+(1-gdn_ii)(1-gdn_jj)-gdn_ji*gdn_ij
    ! i = j
    !        = 2(1-gup_ii)(1-gdn_ii)+(1-gup_ii)+(1-gdn_ii)
    !
    !sdwz(i,j)=<m^z_i*m^z_j>=(nup_i-ndn_i)(nup_j-ndn_j)
    ! i!= j
    !         = (1-gup_ii)(1-gup_jj)-gup_ji*gup_ij-(1-gup_ii)(1-gdn_jj)
    !          -(1-gdn_ii)(1-gup_jj)+(1-gdn_ii)(1-gdn_jj)-gdn_ji*gdn_ij
    ! i = j
    !         =-2(1-gup_ii)(1-gdn_ii)+(1-gup_ii)+(1-gdn_ii)
    !sdwxy(i,j)=<m^x_i*m^x_j>=<(C^+_iup*C_idn+C^+_idn*C_iup)(C^+_jup*C_jdn+C^+_jdn*C_jup)>
    ! i!= j
    !          =-gup_ji*gdn_ij-gdn_ji*gup_ij
    ! i = j
    !          = (1-gup_ii)*gdn_ii+(1-gdn_ii)*gup_ii
    
    double precision::up1,dn1,up2,dn2
    double precision::cdwIn(NXY,NXY),sdwxyIn(NXY,NXY),pmdfIn(NXY,NXY),dsfIn(NXY,NXY)
    double precision:: upup,updn,dnup,dndn
    integer:: it,i,j,is1,is2
    double precision::tmp5,tmp6
    double complex:: g2(NXYS,NXYS)
    double precision:: cckq
    integer:: m,n,q,mpos,mp1,mp2,mp3,mp4,nn1,nn2,inn,jnn
    double precision::xdiff,ydiff
    double precision :: p_sdp, p_pup
    do i=2,NXY
        do j=1,i-1
            do it=0,NT-1
            upup=1.0d0+dreal(-gt(i,i,it)-gt(j,j,it)+gt(i,i,it)*gt(j,j,it)-gt(j,i,it)*gt(i,j,it))
            dndn=1.0d0+dreal(-gt(i+NXY,i+NXY,it)-gt(j+NXY,j+NXY,it)+gt(i+NXY,i+NXY,it)*gt(j+NXY,j+NXY,it)&
            -gt(j+NXY,i+NXY,it)*gt(i+NXY,j+NXY,it))
            updn=1.0d0+dreal(-gt(i,i,it)-gt(j+NXY,j+NXY,it)+gt(i,i,it)*gt(j+NXY,j+NXY,it)&
            -gt(j+NXY,i,it)*gt(i,j+NXY,it))
            dnup=1.0d0+dreal(-gt(i+NXY,i+NXY,it)-gt(j,j,it)+gt(i+NXY,i+NXY,it)*gt(j,j,it)&
            -gt(j,i+NXY,it)*gt(i+NXY,j,it))
                cdwIn(i,j)=cdwIn(i,j)+upup+dndn+updn+dnup
            !sdw-z
                cdwIn(j,i)=cdwIn(j,i)+upup+dndn-updn-dnup
    !       sdw-y dir
                sdwxyIn(i,j)=sdwxyIn(i,j)+dreal(-gt(i+NXY,i,it)*gt(j+NXY,j,it)+gt(j+NXY,i,it)*gt(i+NXY,j,it)&
                -gt(i,i+NXY,it)*gt(j,j+NXY,it)+gt(j,i+NXY,it)*gt(i,j+NXY,it)&
                +gt(i+NXY,i,it)*gt(j,j+NXY,it)-gt(j,i,it)*gt(i+NXY,j+NXY,it)&
                +gt(i,i+NXY,it)*gt(j+NXY,j,it)-gt(j+NXY,i+NXY,it)*gt(i,j,it))
            enddo
        enddo
    enddo
    
    do i=1,NXY            ! calculate the diagonal terms
        do it=0,NT-1
            up1 = 1.0d0 - dreal(gt(i,i,it))
            dn1 = 1.0d0 - dreal(gt(i+NXY,i+NXY,it))
            !cdwIn lack <c^+ c> + <c^+ c>, would be added in Fourier transformation
            cdwIn(i,i)=cdwIn(i,i)+2.0d0*up1*dn1-2.0d0*dreal(gt(i,i+NXY,it)*(gt(i+NXY,i,it)))
            sdwxyIn(i,i)=sdwxyIn(i,i)+2.0d0*dreal(gt(i,i+NXY,it)*gt(i+NXY,i,it))&
                    +up1*dreal(gt(i+NXY,i+NXY,it))+dn1*dreal(gt(i,i,it))
        enddo
    enddo
    
    do i=1,NXY
        do j=1,NXY
            do it=0,NT-1
            up1 = 1.0d0 - dreal(gt(i,i,it))
            dn1 = 1.0d0 - dreal(gt(i+NXY,i+NXY,it))
            up2 = 1.0d0 - dreal(gt(j,j,it))
            dn2 = 1.0d0 - dreal(gt(j+NXY,j+NXY,it))
    
            upup=up1*up2-gt(j,i,it)*gt(i,j,it)
            dndn=dn1*dn2-gt(j+NXY,i+NXY,it)*gt(i+NXY,j+NXY,it)
            updn=up1*dn2+up2*dn1
            if (i==j) then
                pmdfIn(i,i)=pmdfIn(i,i) + 1.0d0+dreal(-gt(i,i,it)-gt(i+NXY,i+NXY,it)&
                +gt(i,i,it)*gt(i+NXY,i+NXY,it)-gt(i+NXY,i,it)*gt(i,i+NXY,it))
                dsfIn(i,i)= dsfIn(i,i)+ up1*dn1
            else
                pmdfIn(i,j)= pmdfIn(i,j)+ dreal(gt(j,i,it)*gt(j+NXY,i+NXY,it)&
                -gt(j+NXY,i,it)*gt(j,i+NXY,it))
                dsfIn(i,j)= dsfIn(i,j)+upup*dndn
            endif
            enddo
        enddo
    enddo
    !write(120,'(48 f 10.5)') pmdfIn
    
    !        Meas1+13*Nxy*Nxy    sp
    !        Meas1+14*Nxy*Nxy    sep
    !        Meas1+15*Nxy*Nxy    dp
    !        Meas1+16*Nxy*Nxy    pp
    !        Meas1+17*Nxy*Nxy    pupp
    caverl(Meas1+13*Nxy*Nxy+1:Meas1+18*Nxy*Nxy)=0.0d0
    mpos=Meas1+13*Nxy*Nxy
    do i=1,NXY
        do j=1,NXY
            q=(i-1)*NXY+j
            do it=0,NT-1
                caverl(mpos+q)=caverl(mpos+q)+2.0d0*dreal(g2t(j,i,it)*g2t(j+NXY,i+NXY,it))
            enddo
        enddo
    enddo
    
    mp1=Meas1+14*Nxy*Nxy
    mp2=Meas1+15*Nxy*Nxy
    mp3=Meas1+16*Nxy*Nxy
    mp4=Meas1+17*Nxy*Nxy
    do i=1,NXY
        do j=1,NXY
            q=(i-1)*NXY+j
            do it=0,NT-1
                do nn1=1,NCN
                    do nn2=1,NCN
                    inn=nextsite(nn1,i)
                    jnn=nextsite(nn2,j)
                    p_sdp = dreal(g2t(j,i,it)*g2t(jnn+NXY,inn+NXY,it) + g2t(jnn,inn,it)*g2t(j+NXY,i+NXY,it))
                    p_pup = dreal(g2t(j,i,it)*g2t(jnn,inn,it) - g2t(jnn,i,it)*g2t(j,inn,it))                    
                    caverl(mp1+q) = caverl(mp1+q) + sefac(nn1)*sefac(nn2)*p_sdp
                    caverl(mp2+q) = caverl(mp2+q) + dfac(nn1)*dfac(nn2)*p_sdp
                    caverl(mp3+q) = caverl(mp3+q) + pfac(nn1)*pfac(nn2)*p_sdp
                    caverl(mp4+q) = caverl(mp4+q) + pupfac(nn1)*pupfac(nn2)*p_pup
                    enddo
                enddo
            enddo
        enddo
    enddo
    
    return
    end
    
    !13.5~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine twoSitetMeas(cdwtIn,sdwtIn,sdwxytIn)
    use mbss
    implicit none
    !Input: Unequal time green function, directly use global gt
    
    !Output:
    !real space,only two obtained
    !        cdw1 -- charge density wave at unequal time
    !        sdwz1 -- spin density wave at unequal time z component
    !        sdwxy1 -- spin density wave at unequal time xy component
    !(not obtained here)
    !k space, after collecting all data, doing fourier tranGensform
    !       cdw1  -->„»csup(q),xcsup(i)
    !       sdwz1 --?„»ssupz(q),xssup(i)
    !       sdwxy1-->„»ssupxy(q),xssupxy(i)
    
    !Details Method
    !G~^rho_ji(l;0)=gt2(j,i,l,rho)=<C^+_irho(l)C_jrho(0)>
    ! G^rho_ij(l;0)=gt1(i,j,l,rho)=<C_irho(l)C^+_jrho(0)>
    !cdwt(i,j)=<(n_iup(l)+n_idn(l))(n_jup(0)+n_jdn(0))>
    !         = (1-gup_ii(l))(1-gup_jj(0))+G~^up_ji(l;0)*G^up_ij(l;0)+(1-gup_ii(l))(1-gdn_jj(0))
    !          +(1-gdn_ii(l))(1-gup_jj(0))+(1-gdn_ii(l))(1-gdn_jj(0))+G~^dn_ji(l;0)*G^dn_ij(l;0)
    !
    !sdwzt(i,j)=<m^z_i(l)*m^z_j(0)>=(nup_i(l)-ndn_i(l))(nup_j(0)-ndn_j(0))
    !          = (1-gup_ii(l))(1-gup_jj(0))+G~^up_ji(l;0)*G^up_ij(l;0)-(1-gup_ii(l))(1-gdn_jj(0))
    !          -(1-gdn_ii(l))(1-gup_jj(0))+(1-gdn_ii(l))(1-gdn_jj(0))+G~^dn_ji(l;0)*G^dn_ij(l;0)
    !
    !sdwxyt(i,j)=<m^x_i(l)*m^x_j(0)>=<(C^+_iup(l)*C_idn(l)+C^+_idn(l)*C_iup(l))(C^+_jup(0)*C_jdn(0)+C^+_jdn(0)*C_jup(0))>
    !           = G~^up_ji(l;0)*G^up_ij(l;0)+G~^dn_ji(l;0)*G^dn_ij(l;0)
    
    double precision::cdwtIn(NXY,NXY),sdwtIn(NXY,NXY),sdwxytIn(NXY,NXY)
    double precision::upup,dndn,updn,dnup
    double precision::up1,up2,dn1,dn2

    double precision :: p_sdp, p_pup
    double precision::tmp5,tmp6,tmp7,tmp8
    double precision::xdiff,ydiff
    integer:: i,j,it,m,n,q,mpos
    double precision:: cckq
    integer:: mp1,mp2,mp3,mp4,mp5,mp6,mp7,mp8
    integer:: s1,s2,s3,s4
    integer:: fac(0:1)
    integer:: inn,jnn,nn1,nn2
    !write(11,*) cdwtIn(1,1)
    do i=1,NXY
        do j=1,NXY
            do it=0,NT-1
                up1 = 1.0d0 - dreal(gt(i,i,it))
                up2 = 1.0d0 - dreal(gt(j,j,0))
                dn1 = 1.0d0 - dreal(gt(i+NXY,i+NXY,it))
                dn2 = 1.0d0 - dreal(gt(j+NXY,j+NXY,0))
    
                !up_up1(time)=gt2(j,i,time,1)*gt1(i,j,time,1)
                !dn_dn1(time)=gt2(j,i,time,2)*gt1(i,j,time,2)
                upup=up1*up2+dreal(gt2(j,i,it)*gt1(i,j,it))
                dndn=dn1*dn2+dreal(gt2(j+NXY,i+NXY,it)*gt1(i+NXY,j+NXY,it))
                updn=up1*dn2+dreal(gt2(j+NXY,i,it)*gt1(i,j+NXY,it))
                dnup=dn1*up2+dreal(gt2(j,i+NXY,it)*gt1(i+NXY,j,it))
                tmp5=dreal(gt(i+NXY,i,it)*gt(j+NXY,j,0)+gt2(j+NXY,i,it)*gt1(i+NXY,j,it))
                tmp6=dreal(gt(i,i+NXY,it)*gt(j,j+NXY,0)+gt2(j,i+NXY,it)*gt1(i,j+NXY,it))
                tmp7=dreal(gt(i+NXY,i,it)*gt(j,j+NXY,0)+gt2(j,i,it)*gt1(i+NXY,j+NXY,it))
                tmp8=dreal(gt(i,i+NXY,it)*gt(j+NXY,j,0)+gt2(j+NXY,i+NXY,it)*gt1(i,j,it))
    !            write(11,*) upup
    !            cdwtIn(i,j)=  cdwtIn(i,j)+upup+dndn+updn+dnup
    !            !cdw1, averages over time
    !            sdwtIn(i,j)= sdwtIn(i,j)+upup+dndn-updn-dnup
    !            !for spin orbit coupling -- <my_i(tau)*my_j(0)>
    !            sdwxytIn(i,j)=sdwxytIn(i,j)-tmp5-tmp6+tmp7+tmp8
                    !try boole's rule
                if (mod(it,4)==0) then
                    cdwtIn(i,j)=  cdwtIn(i,j)+(upup+dndn+updn+dnup)*28.0d0/45.0d0
                    sdwtIn(i,j)= sdwtIn(i,j)+(upup+dndn-updn-dnup)*28.0d0/45.0d0
                    sdwxytIn(i,j)=sdwxytIn(i,j)+(-tmp5-tmp6+tmp7+tmp8)*28.0d0/45.0d0
                elseif (mod(it,4)==1) then
                    cdwtIn(i,j)=  cdwtIn(i,j)+(upup+dndn+updn+dnup)*64.0d0/45.0d0
                    sdwtIn(i,j)= sdwtIn(i,j)+(upup+dndn-updn-dnup)*64.0d0/45.0d0
                    sdwxytIn(i,j)=sdwxytIn(i,j)+(-tmp5-tmp6+tmp7+tmp8)*64.0d0/45.0d0
                elseif (mod(it,4)==2) then
                    cdwtIn(i,j)=  cdwtIn(i,j)+(upup+dndn+updn+dnup)*24.0d0/45.0d0
                    sdwtIn(i,j)= sdwtIn(i,j)+(upup+dndn-updn-dnup)*24.0d0/45.0d0
                    sdwxytIn(i,j)=sdwxytIn(i,j)+(-tmp5-tmp6+tmp7+tmp8)*24.0d0/45.0d0
                elseif (mod(it,4)==3) then
                    cdwtIn(i,j)=  cdwtIn(i,j)+(upup+dndn+updn+dnup)*64.0d0/45.0d0
                    sdwtIn(i,j)= sdwtIn(i,j)+(upup+dndn-updn-dnup)*64.0d0/45.0d0
                    sdwxytIn(i,j)=sdwxytIn(i,j)+(-tmp5-tmp6+tmp7+tmp8)*64.0d0/45.0d0
                endif
            enddo
        enddo
    enddo
    ! normalize time-sum -> approximate integral: multiply by delta_tau = beta / NT
    do i=1,NXY
      do j=1,NXY
        cdwtIn(i,j) = cdwtIn(i,j) * beta 
        sdwtIn(i,j) = sdwtIn(i,j) * beta
        sdwxytIn(i,j) = sdwxytIn(i,j) * beta
      end do
    end do
    
    !write(11,*) cdwtIn(1,1)
    ! current current correlation
    caverl(MEAS21+1:MEAS21+NXY*NXY*NT)=0.0d0
    !<jxp(i,tau)*jxp(j,0)>
    do i=1,NXY
        do j=1,NXY
            do s1=0,1
                do s2=0,1
                    mp1=i+s1*NXY
                    mp2=nextsite(1,i)+s1*NXY
                    mp3=j+s2*NXY
                    mp4=nextsite(1,j)+s2*NXY
            do it=0,NT-1
                q=it*NXY*NXY+(i-1)*NXY+j
    !<ct i+x s1^+ ct i s1 c0 j+x s2^+ c0 j s2>
    caverl(MEAS21+q)=caverl(MEAS21+q)-dreal(g2t(mp1,mp2,it)*g2t(mp3,mp4,0)+gt2(mp3,mp2,it)*gt1(mp1,mp4,it))
    !<ct i s1^+ ct i+x s1 c0 j+x s2^+ c0 j s2>
    caverl(MEAS21+q)=caverl(MEAS21+q)+dreal(g2t(mp2,mp1,it)*g2t(mp3,mp4,0)+gt2(mp3,mp1,it)*gt1(mp2,mp4,it))
    !<ct i+x s1^+ ct i s1 c0 j s2^+ c0 j+x s2>
    caverl(MEAS21+q)=caverl(MEAS21+q)+dreal(g2t(mp1,mp2,it)*g2t(mp4,mp3,0)+gt2(mp4,mp2,it)*gt1(mp1,mp3,it))
    !<ct i s1^+ ct i+x s1 c0 j s2^+ c0 j+x s2>
    caverl(MEAS21+q)=caverl(MEAS21+q)-dreal(g2t(mp2,mp1,it)*g2t(mp4,mp3,0)+gt2(mp4,mp1,it)*gt1(mp2,mp3,it))
                    enddo
                enddo
            enddo
        enddo
    enddo
    
    fac(0)=1
    fac(1)=-1
    !<jxp(i,tau)*jxl(j,0)>
    do i=1,NXY
        do j=1,NXY
            do s1=0,1
                do s2=0,1
                    s3=-s2+1
                    mp1=i+s1*NXY
                    mp2=nextsite(1,i)+s1*NXY
                    mp3=j+s2*NXY
                    mp4=j+s3*NXY
                    mp5=nextsite(1,j)+s2*NXY
                    mp6=nextsite(1,j)+s3*NXY
            do it=0,NT-1
                q=it*NXY*NXY+(i-1)*NXY+j
    !<ct i+x s1^+ ct i s1 c0 j+x s2^+ c0 j s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)-lamda*fac(s2)*dreal(g2t(mp1,mp2,it)*g2t(mp4,mp5,0)&
                                    +gt2(mp4,mp2,it)*gt1(mp1,mp5,it))
    !<ct i+x s1^+ ct i s1 c0 j s2^+ c0 j+x s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)-lamda*fac(s2)*dreal(g2t(mp1,mp2,it)*g2t(mp6,mp3,0)&
                                    +gt2(mp6,mp2,it)*gt1(mp1,mp3,it))
    !<ct i s1^+ ct i+x s1 c0 j+x s2^+ c0 j s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)+lamda*fac(s2)*dreal(g2t(mp2,mp1,it)*g2t(mp4,mp5,0)&
                                    +gt2(mp4,mp1,it)*gt1(mp2,mp5,it))
    !<ct i s1^+ ct i+x s1 c0 j s2^+ c0 j+x s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)+lamda*fac(s2)*dreal(g2t(mp2,mp1,it)*g2t(mp6,mp3,0)&
                                    +gt2(mp6,mp1,it)*gt1(mp2,mp3,it))
                    enddo
                enddo
            enddo
        enddo
    enddo
    !<jxl(i,tau)*jxp(j,0)>
    do i=1,NXY
        do j=1,NXY
            do s1=0,1
                do s2=0,1
                    s3=-s2+1
                    mp1=j+s1*NXY
                    mp2=nextsite(1,j)+s1*NXY
                    mp3=i+s2*NXY
                    mp4=i+s3*NXY
                    mp5=nextsite(1,i)+s2*NXY
                    mp6=nextsite(1,i)+s3*NXY
            do it=0,NT-1
                q=it*NXY*NXY+(i-1)*NXY+j
    !<ct i+x s1^+ ct i s1 c0 j+x s2^+ c0 j s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)-lamda*fac(s2)*dreal(g2t(mp4,mp5,it)*g2t(mp1,mp2,0)&
                                    +gt2(mp1,mp5,it)*gt1(mp4,mp2,it))
    !<ct i+x s1^+ ct i s1 c0 j s2^+ c0 j+x s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)-lamda*fac(s2)*dreal(g2t(mp6,mp3,it)*g2t(mp1,mp2,0)&
                                    +gt2(mp1,mp3,it)*gt1(mp6,mp2,it))
    !<ct i s1^+ ct i+x s1 c0 j+x s2^+ c0 j s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)+lamda*fac(s2)*dreal(g2t(mp4,mp5,it)*g2t(mp2,mp1,0)&
                                    +gt2(mp2,mp5,it)*gt1(mp4,mp1,it))
    !<ct i s1^+ ct i+x s1 c0 j s2^+ c0 j+x s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)+lamda*fac(s2)*dreal(g2t(mp6,mp3,it)*g2t(mp2,mp1,0)&
                                    +gt2(mp2,mp3,it)*gt1(mp6,mp1,it))
                    enddo
                enddo
            enddo
        enddo
    enddo
    !<jxl(i,tau)*jxp(j,0)>
    do i=1,NXY
        do j=1,NXY
            do s1=0,1
                do s3=0,1
                    s2=-s1+1
                    s4=-s3+1
                    mp1=i+s1*NXY
                    mp2=i+s2*NXY
                    mp3=nextsite(1,i)+s1*NXY
                    mp4=nextsite(1,i)+s2*NXY
                    mp5=j+s3*NXY
                    mp6=j+s4*NXY
                    mp7=nextsite(1,j)+s3*NXY
                    mp8=nextsite(1,j)+s4*NXY
            do it=0,NT-1
                q=it*NXY*NXY+(i-1)*NXY+j
    !<ct i+x s1^+ ct i s1 c0 j+x s2^+ c0 j s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)-lamda*lamda*fac(s1)*fac(s3)*dreal(g2t(mp2,mp3,it)*g2t(mp6,mp7,0)&
                                    +gt2(mp6,mp3,it)*gt1(mp2,mp7,it))
    !<ct i+x s1^+ ct i s1 c0 j s2^+ c0 j+x s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)-lamda*lamda*fac(s1)*fac(s3)*dreal(g2t(mp2,mp3,it)*g2t(mp8,mp5,0)&
                                    +gt2(mp8,mp3,it)*gt1(mp2,mp5,it))
    !<ct i s1^+ ct i+x s1 c0 j+x s2^+ c0 j s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)-lamda*lamda*fac(s1)*fac(s3)*dreal(g2t(mp4,mp1,it)*g2t(mp6,mp7,0)&
                                    +gt2(mp6,mp1,it)*gt1(mp4,mp7,it))
    !<ct i s1^+ ct i+x s1 c0 j s2^+ c0 j+x s3>
    caverl(MEAS21+q)=caverl(MEAS21+q)-lamda*lamda*fac(s1)*fac(s3)*dreal(g2t(mp4,mp1,it)*g2t(mp8,mp5,0)&
                                    +gt2(mp8,mp1,it)*gt1(mp4,mp5,it))
                    enddo
                enddo
            enddo
        enddo
    enddo
    
    
    !        Meas1+8*Nxy*Nxy    spl
    !        Meas1+9*Nxy*Nxy    sepl
    !        Meas1+10*Nxy*Nxy    dpl
    !        Meas1+11*Nxy*Nxy    ppl
    !        Meas1+12*Nxy*Nxy    puppl
    caverl(Meas1+8*Nxy*Nxy+1:Meas1+13*Nxy*Nxy)=0.0d0
    mpos=Meas1+8*Nxy*Nxy
    do i=1,NXY
        do j=1,NXY
            q=(i-1)*NXY+j
            do it=0,NT-1
                ! up-dn + dn-up
                if (mod(it,4)==0) then
                    caverl(mpos+q) = caverl(mpos+q) + 28.0d0/45.0d0 * &
                        dreal(gt2(j,i,it)*gt2(j+NXY,i+NXY,it) + gt2(j+NXY,i+NXY,it)*gt2(j,i,it))
                elseif (mod(it,4)==1) then
                    caverl(mpos+q) = caverl(mpos+q) + 64.0d0/45.0d0 * &
                        dreal(gt2(j,i,it)*gt2(j+NXY,i+NXY,it) + gt2(j+NXY,i+NXY,it)*gt2(j,i,it))
                elseif (mod(it,4)==2) then
                    caverl(mpos+q) = caverl(mpos+q) + 24.0d0/45.0d0 * &
                        dreal(gt2(j,i,it)*gt2(j+NXY,i+NXY,it) + gt2(j+NXY,i+NXY,it)*gt2(j,i,it))
                elseif (mod(it,4)==3) then
                    caverl(mpos+q) = caverl(mpos+q) + 64.0d0/45.0d0 * &
                        dreal(gt2(j,i,it)*gt2(j+NXY,i+NXY,it) + gt2(j+NXY,i+NXY,it)*gt2(j,i,it))
                endif
            enddo
        enddo
    enddo
    
    mp1=Meas1+9*Nxy*Nxy
    mp2=Meas1+10*Nxy*Nxy
    mp3=Meas1+11*Nxy*Nxy
    mp4=Meas1+12*Nxy*Nxy
    do i=1,NXY
        do j=1,NXY
            q=(i-1)*NXY+j
            do it=0,NT-1
                do nn1=1,NCN
                    do nn2=1,NCN
                    inn=nextsite(nn1,i)
                    jnn=nextsite(nn2,j)
                    ! Singlet (up-dn + dn-up)
                    p_sdp = dreal(gt2(j,i,it)*gt2(jnn+NXY,inn+NXY,it) + gt2(jnn,inn,it)*gt2(j+NXY,i+NXY,it))
                        
                    ! Triplet (up-up)
                    p_pup = dreal(gt2(j,i,it)*gt2(jnn,inn,it) - gt2(jnn,i,it)*gt2(j,inn,it))
                    
                    if (mod(it,4)==0) then
                        caverl(mp1+q)=caverl(mp1+q)+28.0d0/45.0d0*sefac(nn1)*sefac(nn2)*p_sdp
                        caverl(mp2+q)=caverl(mp2+q)+28.0d0/45.0d0*dfac(nn1)*dfac(nn2)*p_sdp
                        caverl(mp3+q)=caverl(mp3+q)+28.0d0/45.0d0*pfac(nn1)*pfac(nn2)*p_sdp
                        caverl(mp4+q)=caverl(mp4+q)+28.0d0/45.0d0*pupfac(nn1)*pupfac(nn2)*p_pup
                    elseif (mod(it,4)==1) then
                        caverl(mp1+q)=caverl(mp1+q)+64.0d0/45.0d0*sefac(nn1)*sefac(nn2)*p_sdp
                        caverl(mp2+q)=caverl(mp2+q)+64.0d0/45.0d0*dfac(nn1)*dfac(nn2)*p_sdp
                        caverl(mp3+q)=caverl(mp3+q)+64.0d0/45.0d0*pfac(nn1)*pfac(nn2)*p_sdp
                        caverl(mp4+q)=caverl(mp4+q)+64.0d0/45.0d0*pupfac(nn1)*pupfac(nn2)*p_pup
                    elseif (mod(it,4)==2) then
                        caverl(mp1+q)=caverl(mp1+q)+24.0d0/45.0d0*sefac(nn1)*sefac(nn2)*p_sdp
                        caverl(mp2+q)=caverl(mp2+q)+24.0d0/45.0d0*dfac(nn1)*dfac(nn2)*p_sdp
                        caverl(mp3+q)=caverl(mp3+q)+24.0d0/45.0d0*pfac(nn1)*pfac(nn2)*p_sdp
                        caverl(mp4+q)=caverl(mp4+q)+24.0d0/45.0d0*pupfac(nn1)*pupfac(nn2)*p_pup
                    elseif (mod(it,4)==3) then
                        caverl(mp1+q)=caverl(mp1+q)+64.0d0/45.0d0*sefac(nn1)*sefac(nn2)*p_sdp
                        caverl(mp2+q)=caverl(mp2+q)+64.0d0/45.0d0*dfac(nn1)*dfac(nn2)*p_sdp
                        caverl(mp3+q)=caverl(mp3+q)+64.0d0/45.0d0*pfac(nn1)*pfac(nn2)*p_sdp
                        caverl(mp4+q)=caverl(mp4+q)+64.0d0/45.0d0*pupfac(nn1)*pupfac(nn2)*p_pup
                    endif
                    enddo
                enddo
            enddo
        enddo
    enddo
    
    return
    end
    
    !17~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine bpmult (it,gh)
    use mbss
    implicit none
    ! ...     exp(-dt*H)*exp(-dt*V)*G
    integer::i,it,j
    double complex::gh(NXYS,NXYS), temp(NXYS,NXYS), scal
    do i=1,NXY
    !	   scal=emv(spin(i,it),i,ispin)
        scal=emv(spin(i,it),i)
        temp(i,:)=scal*gh(i,:)
        j=i+NXY
        scal=emv(spin(i,it),j)
        temp(j,:)=gh(j,:)*scal
    enddo
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), emh, NXYS, temp, NXYS, dcmplx(0.0d0,0.0d0), gh, NXYS)
    !gh=matmul(emh(:,:),temp)
    return
    end
    
    !18~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine bmmult (it,gh)
    use mbss
    implicit none
    ! ...    exp(+dt*V)*exp(+dt*H)*G
    integer::i,it,j
    double complex::gh(NXYS,NXYS), temp(NXYS,NXYS), scal
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), eph, NXYS, gh, NXYS, dcmplx(0.0d0,0.0d0), temp, NXYS)
    !temp=matmul(eph(:,:),gh)
    do i=1,NXY
    !        scal=epv(spin(i,it),ispin)
        scal=epv(spin(i,it),i)
        gh(i,:)=scal*temp(i,:)
        j=i+NXY
        scal=epv(spin(i,it),j)
        gh(j,:)=scal*temp(j,:)
    enddo
    return
    end
    
    !19~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine multbm (it,gh)
    use mbss
    implicit none
    ! ...    G*exp(+dt*V)*exp(+dt*H)
    integer::i,it,spos,j
    double complex::gh(NXYS,NXYS), temp(NXYS,NXYS), scal
    do i=1,NXY
    !         scal=epv(spin(i,it),ispin)
        scal=epv(spin(i,it),i)
        temp(:,i)=gh(:,i)*scal
        j=i+NXY
        scal=epv(spin(i,it),j)
        temp(:,j)=gh(:,j)*scal
    enddo
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), temp, NXYS, eph, NXYS, dcmplx(0.0d0,0.0d0), gh, NXYS)
    !gh=matmul(temp,eph(:,:))
    return
    end
    
    !20~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    subroutine multbp (it,gh)
    use mbss
    implicit none
    ! ...    G*exp(-dt*H)*exp(-dt*V)
    integer::i,it,j
    double complex::gh(NXYS,NXYS), temp(NXYS,NXYS), scal
    call zgemm('N', 'N', NXYS, NXYS, NXYS, dcmplx(1.0d0,0.0d0), gh, NXYS, emh, NXYS, dcmplx(0.0d0,0.0d0), temp, NXYS)
    !temp=matmul(gh,emh(:,:))
    do i=1,NXY
    !         scal=emv(spin(i,it),ispin)
        scal=emv(spin(i,it),i)
        gh(:,i)=temp(:,i)*scal
        j=i+NXY
        scal=emv(spin(i,it),j)
        gh(:,j)=temp(:,j)*scal
    enddo
    return
    end
    
    subroutine udv(NA,uh,sh,vh)
    use mbss
    implicit none
    integer::NA
    double complex::uh(NA,NA),sh(NA),vh(NA,NA)
    double complex::xn1(NA)
    double complex::temp(NA,NA)
    double precision::xnorm(NA),vhelp(NA),xmax
    double complex::v1(NA,NA)
    !double precision::test(NA,NA),test1(NA,NA)
    double complex::zeta(NA),work(4*NA)
    double complex::x1,shi
    integer::i,nr,j,imax,iv,iinfo,j1
    integer::ivpt(NA),ivptm1(NA)
    !.. Normalize vectors:
    !write(180,*) 'ini udv'
    !write(180,'(8 f 10.5)') dreal(uh)
    do i=1,NA
        xnorm(i)=0.0d0
        do nr=1,NA
            xnorm(i)=xnorm(i)+dreal(uh(nr,i))*dreal(uh(nr,i))+dimag(uh(nr,i))*dimag(uh(nr,i))
        enddo
        xnorm(i)=dsqrt(xnorm(i))
        vhelp(i)=xnorm(i)
    enddo
    !write(180,*) vhelp
    do i=1,NA
        xmax=0.0d0
        do j=1,NA
            if(vhelp(j)>xmax)then
                imax=j
                xmax=vhelp(j)
            endif
        enddo
        vhelp(imax)=-1.0d0
        ivpt(i)=imax
        ivptm1(imax)=i
    enddo
    !write(180,*) ivpt
    !do i=NA/2+1,NA
    !    xnorm(i)=0.0d0
    !    do nr=1,NA
    !        xnorm(i)=xnorm(i)+real(uh(nr,i))*real(uh(nr,i))+aimag(uh(nr,i))*aimag(uh(nr,i))
    !    enddo
    !    xnorm(i)=dsqrt(xnorm(i))
    !    vhelp(i)=xnorm(i)
    !enddo
    !write(180,*) vhelp
    !do i=NA/2+1,NA
    !    xmax=0.0d0
    !    do j=NA/2+1,NA
    !        if(vhelp(j)>xmax)then
    !            imax=j
    !            xmax=vhelp(j)
    !        endif
    !    enddo
    !    vhelp(imax)=-1.0d0
    !    ivpt(i)=imax
    !    ivptm1(imax)=i
    !enddo
    !write(180,*) ivpt
    
    do i=1,NA
        iv=ivpt(i)
        x1=dcmplx(xnorm(iv),0.0d0)
        do nr=1,NA
            temp(nr,i)=uh(nr,iv)/x1
            v1(i,nr)=0.0d0
        enddo
    enddo
    !    write(180,*) 'real temp b'
    !write(180,'(8 f 10.5)') real(temp)
    !write(180,*) 'img  temp b'
    !write(180,'(8 f 10.5)') aimag(temp)
    !
    !  TEMP.Nag routines
    call ZGEQRF(NA,NA,temp,NA,zeta,work,4*NA,iinfo)
    !write(180,*) 'ZGEQRF ierror',iinfo
    !    write(180,*) 'real temp'
    !write(180,'(8 f 10.5)') real(temp)
    !write(180,*) 'img temp'
    !write(180,'(8 f 10.5)') aimag(temp)
    !  Scale v1 to a unit triangluar matrix
    do i=1,NA
        xn1(i)=dcmplx(xnorm(ivpt(i)),0.0d0)
        sh(i)=temp(i,i)
        shi=sh(i)
        sh(i)=xn1(i)*shi
        do j=i,NA
            v1(i,j)=temp(i,j)/shi
        enddo
    enddo
    !    write(180,*) 'real v1'
    !write(180,'(8 f 10.5)') real(v1)
    !write(180,*) 'img uhv1'
    !write(180,'(8 f 10.5)') aimag(v1)
    !  finish the pivoting.
    do i=1,NA-1
        x1=xn1(i)
        do j=i+1,NA
            v1(i,j)=v1(i,j)*xn1(j)/x1
        enddo
    enddo
    !    write(180,*) 'real v1-rearrange'
    !write(180,'(8 f 10.5)') real(v1)
    !write(180,*) 'img v1-rearrange'
    !write(180,'(8 f 10.5)') aimag(v1)
    do j=1,NA
        j1=ivptm1(j)
        vh(:,j)=v1(:,j1)
    enddo
    !    write(180,*) 'real vh-rearrange'
    !write(180,'(8 f 10.5)') real(vh)
    !write(180,*) 'img vh-rearrange'
    !write(180,'(8 f 10.5)') aimag(vh)
    !  Compute U -> TEMP
    call ZUNGQR(NA,NA,NA,temp,NA,zeta,work,4*NA,iinfo)
    !write(180,*) 'ZUNGQR ierror',iinfo
    do j=1,NA
        do i=1,NA
            uh(i,j)=temp(i,j)
        enddo
    enddo
    !    write(180,*) 'uh'
    !write(180,'(8 f 10.5)') real(uh)
    !write(180,*) 'uh'
    !write(180,'(8 f 10.5)') aimag(uh)
    !        do i=1,NA
    !            do j=1,NA
    !                temp(i,j)=uh(i,j)*sh(j)
    !            enddo
    !        enddo
    !temp=matmul(temp,vh)
    !write(180,*) 'prod udv'
    !write(180,'(8 f 10.5)') real(temp)
    return
    endsubroutine
    
    function trans(NA,ma)
    implicit none
    integer::NA,i,j
    double complex::ma(NA,NA),trans(NA,NA)
    do i=1,NA
        do j=1,NA
            trans(i,j)=ma(j,i)
        enddo
    enddo
    endfunction
    
    double complex function det(N, mat)
    implicit none
    integer:: N
    double complex, intent(inout), dimension(N,N) :: mat
    integer(kind=8) :: i, info
    integer, allocatable :: ipiv(:)
    double precision :: sgn
    allocate(ipiv(N))
    ipiv = 0
    call zgetrf(N, N, mat, N, ipiv, info)
    det = 1.0d0
    do i = 1, N
     det = det*mat(i, i)
    enddo
    sgn = 1.0d0
    do i = 1, N
    if(ipiv(i) /= i) then
        sgn = -sgn
      end if
    end do
    det = sgn*det
    end function det
    
    ! produce u d v of A=BBBBB..
    ! usar,ssar,vsar: RHS. 1. B6 2. B6*B5 3. B6*B5*B4 ...
    ! usal,ssal,vsal: LHS  1. B1 2 B2*B1  3. B3*B2*B1 ...
    subroutine makeA (ib)     !the number "l" of B_l multiply before B_M is it
    use mbss
    use link,only:bpmult,makeb,makeipb,matinv
    implicit none
    integer::i,j,ib,k,it
    double complex:: temp(NXYS,NXYS)
    ! ... initialize the product of the B's by defining the first B
    !     to be the identity factorized as U*W*V (u*s*v)
    if (ib==1) then
        usar(:,:,:)=dcmplx(0.0d0,0.0d0)
        ssar(:,:)=dcmplx(0.0d0,0.0d0)
        vsar(:,:,:)=dcmplx(0.0d0,0.0d0)
        usal(:,:,:)=dcmplx(0.0d0,0.0d0)
        ssal(:,:)=dcmplx(0.0d0,0.0d0)
        vsal(:,:,:)=dcmplx(0.0d0,0.0d0)
        !initialize LHS first matrix
        do i = 1, NXYS
        ssal(i,0) = onec
        vsal(i,i,0) = onec
        usal(i,i,0) = onec
        enddo
    
        do i = 1, NXYS
        do j = 1, NXYS
            temp(i,j)=zeroc
            v(i,j) = zeroc
            u(i,j) = zeroc
        enddo
        temp(i,i)=onec
        s(i) = onec
        v(i,i) = onec
        u(i,i) = onec
        enddo
        do i=1,nb
    !        ... multply Bg
    !        call zgemm('N', 'N', NXYS, NXYS, NXYS, onec,temp(:,:) , NXYS, Bg(:,:,nb-i+1), NXYS, zeroc, v(:,:), NXYS)
    !        call udvBg
    !                call zgemm('N', 'N', NXYS, NXYS, NXYS, onec,temp(:,:) , NXYS,Bg(:,:,nb-i+1) , NXYS, zeroc,u(:,:) , NXYS)
            do j=1,isudv
                it=NT-(i-1)*isudv-j+1
                call multbp(it,u)
            enddo
            call makebt
                do j=1,NXYS
                    do k=1,NXYS
                        usar(j,k,i)=u(j,k)
                        vsar(j,k,i)=v(j,k)
                    enddo
                    ssar(j,i)=s(j)
                enddo
            temp(:,:)=u(:,:)
        enddo
    !for first green function calculation g=(udv)^-1
        u(:,:)=v(:,:)
        v(:,:)=temp(:,:)
    
    !ib>1
        else
    !    do i = 1, NXY
    !    do j = 1, NXY
    !        v(i,j) = zeroc
    !        u(i,j) = zeroc
    !    enddo
    !    s(i) = onec
    !    v(i,i) = onec
    !    u(i,i) = onec
    !    enddo
        do j=1,NXYS
            do k=1,NXYS
                v(j,k)=vsal(j,k,ib-2)
                u(j,k)=usal(j,k,ib-2)
            enddo
            s(j)=ssal(j,ib-2)
        enddo
        do j=1,isudv
            it=j+(ib-2)*isudv
            call bpmult(it,u)
        enddo
    !    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, Bg(:,:,ib-1), NXYS, usal(:,:,ib-2), NXYS, zeroc, u(:,:), NXYS)
    
        call makeb
        do j=1,NXYS
            do k=1,NXYS
                usal(j,k,ib-1)=u(j,k)
                vsal(j,k,ib-1)=v(j,k)
            enddo
            ssal(j,ib-1)=s(j)
        enddo
    
    !!   UL*DL*VL * UR*DR*VR
    !!   temp=(DL*VL*UR*DR)
    !!VL*UR
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, vsal(:,:,ib-1), NXYS, usar(:,:,nb-ib+1), NXYS, zeroc, u(:,:), NXYS)
    !!DL*
    !do i=1,NXYS
    !    do j=1,NXYS
    !        u(i,j)=u(i,j)*ssal(i,ib-1)
    !    enddo
    !enddo
    !
    !!*DR
    !do i=1,NXYS
    !    do j=1,NXYS
    !        u(i,j)=u(i,j)*ssar(j,nb-ib+1)
    !    enddo
    !enddo
    !
    !call udv(NXYS,u,s,v)
    !
    !! u=uL*u v=v*vR
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS, u(:,:), NXYS, zeroc, temp(:,:), NXYS)
    !do i=1,NXYS
    !    do j=1,NXYS
    !        u(i,j)=temp(i,j)
    !    enddo
    !enddo
    !
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, v(:,:), NXYS, vsar(:,:,nb-ib+1), NXYS, zeroc, temp(:,:), NXYS)
    !do i=1,NXYS
    !    do j=1,NXYS
    !        v(i,j)=temp(i,j)
    !    enddo
    !enddo
    
    endif
    !write(60,*) 'B prod',ib
    !    do j=1,NXY
    !        do k=1,NXY
    !            temp(j,k)=u(j,k)*s(k)
    !        enddo
    !    enddo
    !        temp(:,:)=matmul(temp(:,:),v(:,:))
    !    write(60,*) 'udv'
    !    write(60,strfor) dreal(temp(:,:))
    !    write(60,strfor) dimag(temp(:,:))
    !
    !write(60,*) 'makeA',ib
    !do i=1,nb
    !    write(60,*) i
    !    write(60,*) 'ul'
    !    write(60,strfor) dreal(usal(:,:,i,ispin))
    !    write(60,strfor) dimag(usal(:,:,i,ispin))
    !        write(60,*) 'dl'
    !    write(60,strfor) dreal(ssal(:,i,ispin))
    !    write(60,strfor) dimag(ssal(:,i,ispin))
    !            write(60,*) 'vl'
    !    write(60,strfor) dreal(vsal(:,:,i,ispin))
    !    write(60,strfor) dimag(vsal(:,:,i,ispin))
    !    do j=1,NXY
    !        do k=1,NXY
    !            temp(j,k)=usal(j,k,i,ispin)*ssal(k,i,ispin)
    !        enddo
    !    enddo
    !        temp(:,:)=matmul(temp(:,:),vsal(:,:,i,ispin))
    !    write(60,*) 'udv'
    !    write(60,strfor) dreal(temp(:,:))
    !    write(60,strfor) dimag(temp(:,:))
    !
    !enddo
    !do i=1,nb
    !        write(60,*) 'ur'
    !    write(60,strfor) dreal(usar(:,:,i,ispin))
    !    write(60,strfor) dimag(usar(:,:,i,ispin))
    !        write(60,*) 'dr'
    !    write(60,strfor) dreal(ssar(:,i,ispin))
    !    write(60,strfor) dimag(ssar(:,i,ispin))
    !            write(60,*) 'vr'
    !    write(60,strfor) dreal(vsar(:,:,i,ispin))
    !    write(60,strfor) dimag(vsar(:,:,i,ispin))
    !
    !    do j=1,NXY
    !        do k=1,NXY
    !            temp(j,k)=usar(j,k,i,ispin)*ssar(k,i,ispin)
    !        enddo
    !    enddo
    !    temp(:,:)=matmul(temp(:,:),vsar(:,:,i,ispin))
    !    write(60,*) 'udv'
    !    write(60,strfor) dreal(temp(:,:))
    !    write(60,strfor) dimag(temp(:,:))
    !    enddo
    return
    endsubroutine makeA
    
    subroutine makeBg
    use mbss
    implicit none
    integer:: i,j,k,ispin
    do i=1,NXYS
        do j=1,NXYS
            do k=1,nb
            Bg(i,j,k)=dcmplx(0.0d0,0.0d0)
            enddo
        enddo
    enddo
    
    do i=1,NXYS
        do k=1,nb
        Bg(i,i,k)=dcmplx(1.0d0,0.0d0)
        enddo
    enddo
    
    do k=1,nb
    !        ... multply B times U and store the result in U
    !    do i=isudv,1,-1
    do i=1,isudv
            j=(k-1)*isudv+i
            call bpmult (j,Bg(:,:,k))
        enddo
    enddo
    endsubroutine makeBg
    
    subroutine updateBg(ib)
    use mbss
    implicit none
    integer::ib,i,j
    do i=1,NXYS
        do j=1,NXYS
            Bg(i,j,ib)=dcmplx(0.0d0,0.0d0)
        enddo
    enddo
    
    do i=1,NXYS
        Bg(i,i,ib)=dcmplx(1.0d0,0.0d0)
    enddo
    do i=1,isudv
    !do i=isudv,1,-1
        j=(ib-1)*isudv+i
        call bpmult (j,Bg(:,:,ib))
    enddo
    endsubroutine updateBg
    
    !for LHS already calculated, after the cnfmake
    subroutine makeANoCh (ib)     !the number "l" of B_l multiply before B_M is it
    use mbss
    use link,only:bpmult,makeb,makeipb,matinv,trans
    implicit none
    integer::i,j,ib,k,it
    double complex:: temp(NXYS,NXYS)
    ! ... initialize the product of the B's by defining the first B
    !     to be the identity factorized as U*W*V (u*s*v)
    if (ib==1) then
        usar(:,:,:)=dcmplx(0.0d0,0.0d0)
        ssar(:,:)=dcmplx(0.0d0,0.0d0)
        vsar(:,:,:)=dcmplx(0.0d0,0.0d0)
    !    usal(:,:,:,:)=dcmplx(0.0d0,0.0d0)
    !    ssal(:,:,:)=dcmplx(0.0d0,0.0d0)
    !    vsal(:,:,:,:)=dcmplx(0.0d0,0.0d0)
        !initialize LHS first matrix
    !    do i = 1, NXY
    !    ssal(i,1,ispin) = onec
    !    vsal(i,i,1,ispin) = onec
    !    usal(i,i,1,ispin) = onec
    !    enddo
    
        do i = 1, NXYS
        do j = 1, NXYS
            temp(i,j)=zeroc
            v(i,j) = zeroc
            u(i,j) = zeroc
        enddo
        temp(i,i)=onec
        s(i) = onec
        v(i,i) = onec
        u(i,i) = onec
        enddo
        do i=1,nb
    !        ... multply Bg
    !        call zgemm('N', 'N', NXYS, NXYS, NXYS, onec,temp(:,:) , NXYS, Bg(:,:,nb-i+1), NXYS, zeroc, v, NXYS)
    !        call udvBg
    !        call zgemm('N', 'N', NXYS, NXYS, NXYS, onec,temp(:,:) , NXYS,Bg(:,:,nb-i+1) , NXYS, zeroc,u(:,:) , NXYS)
            do j=1,isudv
                it=NT-(i-1)*isudv-j+1
                call multbp(it,u)
            enddo
            call makebt
                do j=1,NXYS
                    do k=1,NXYS
                        usar(j,k,i)=u(j,k)
                        vsar(j,k,i)=v(j,k)
                    enddo
                    ssar(j,i)=s(j)
                enddo
            temp(:,:)=u(:,:)
        enddo
    !for first green function calculation g=(udv)^-1
        u(:,:)=v(:,:)
        v(:,:)=temp(:,:)
    !ib>1
        else
    !!VL*UR
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, vsal(:,:,ib-1), NXYS, usar(:,:,nb-ib+1), NXYS, zeroc, u(:,:), NXYS)
    !!DL*
    !do i=1,NXYS
    !    do j=1,NXYS
    !        u(i,j)=u(i,j)*ssal(i,ib-1)
    !    enddo
    !enddo
    !
    !!*DR
    !do i=1,NXYS
    !    do j=1,NXYS
    !        u(i,j)=u(i,j)*ssar(j,nb-ib+1)
    !    enddo
    !enddo
    !
    !call udv(NXYS,u,s,v)
    !
    !! u=uL*u v=v*vR
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, usal(:,:,ib-1), NXYS, u(:,:), NXYS, zeroc, temp(:,:), NXYS)
    !do i=1,NXYS
    !    do j=1,NXYS
    !        u(i,j)=temp(i,j)
    !    enddo
    !enddo
    !
    !call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, v(:,:), NXYS, vsar(:,:,nb-ib+1), NXYS, zeroc, temp(:,:), NXYS)
    !do i=1,NXYS
    !    do j=1,NXYS
    !        v(i,j)=temp(i,j)
    !    enddo
    !enddo
    endif
    
    
    return
    endsubroutine makeANoCh
    
    !9~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !similar to makeb but it multiplies v
    subroutine udvBg
    use mbss
    use link,only:udv
    implicit none
    integer::i,j
    double complex::temp(NXYS,NXYS)
    ! ... V comes in as V*B          v=v*(b^north)
    ! ... save the current U in TEMP
    do j = 1, NXYS
        do i = 1, NXYS
            v(i,j) = v(i,j)*s(i)          !u=(b^north)*u*s temp=v
            temp(i,j) = u(i,j)
    !        vi(i,j)=zeroc
        end do
    end do
    call udv(NXYS,v,s,vi)
    ! ... form Vi*TEMP, TEMP contains the old V
    !     ... V and Vi are upper triangular, unit diagonal (mgs case)
    !     ... V and Vi are well condition (udv case)
    call zgemm('N', 'N', NXYS, NXYS, NXYS, onec, temp, NXYS, v, NXYS, zeroc, u, NXYS)
    v(:,:)=vi(:,:)
    !v=matmul(vi,temp)
    return
    endsubroutine    
    !-----------------Chiral spin density wave------------