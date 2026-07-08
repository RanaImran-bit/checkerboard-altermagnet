module cpmc
implicit none

integer,parameter::k4b=selected_int_kind(9)
integer,parameter::i4b=selected_int_kind(9)
integer,parameter::i2b=selected_int_kind(2)

! symbolic names for kind types of
!  single- and double-precision reals:

integer,parameter::sp=kind(1.0d0)
integer,parameter::dp=kind(1.0d0)

real(sp), parameter :: PI=3.14159265359
real(sp), parameter :: ZERO=0.0, FOURTH=0.25, HALF=0.5
real(sp), parameter :: ONE=1.0, TWO=2.0, FOUR=4.0, TEN=10.0

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
integer, parameter :: ntypes = 32

integer::lx,ly,lxy
parameter (lx=4,ly=4,lxy=lx*ly) !gedian

integer::NSTATES,NHS
parameter (NSTATES=lxy,NHS=NSTATES*NSTATES)

!        NSTATES  : number of states in SCF solution
!        NHS      : number of HS fields



integer::NUP,NDN,NELEC,NSPIN,NWFBAS,NIW,NWLKRS,M_MAX,NAVE
parameter (NUP=1,NDN=1,NELEC=NUP+NDN,NSPIN=2)
parameter (NWFBAS=1,NIW=1,NWLKRS=1000,M_MAX=125)
parameter (NAVE=(ntypes*NSTATES+ntypes)*NSTATES+6)

!        NUP      : number of up   electrons; Assumed >= NDN
!        NDN      : number of down electrons
!        NELEC    : number of up + down electrons
!        NSPIN    : number of possible spins
!        NWFBAS   : number of determinants in initial wavefunction
!        NWLKRS   : number of walkers
!        M_MAX    : maximum number of back-propagation
!        NAVE     : number of correlations
integer::NELECs(NSPIN)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

integer::ispin,ising

!SetUp

real(sp)::deltau,etrial,t0,t1,t2,ud,vpd,g_ph,w_ph

!        deltau   : delta tau, the Trotter step size
!        etrial   : growth control factor
!        t0       : t_0, hopping integral between NN sites
!        t1       : t_1, hopping integral between NNN sites
!        t2       : t_2, hopping integral between NNNN sites
!        ud       : U_d, Hubbard U on each site
!        vpd      : V_pd, Coulomb interaction between NN sites
!        g_ph     : g_ph, electron-phonon coupling strength
!        w_ph     : w_ph, phonon tunelling frequency

integer::itvlpceq,itvlpcgr,itvlpc,itvlmeas,itvlorth,itvl_m

!        itvlpceq : population control frequency in relaxation phase
!        itvlpcgr : population control frequency in growth estimator
!        itvlpc   : population control frequency in measurement phase
!        itvlmeas : measurement frequency
!        itvlorth : orthogonalization frequency
!        itvl_m   : back-propagator

integer::nblkeq,nblkgr,nblk,nblkstps

!        nblkeq   : number of equilibrium blocks
!        nblkgr   : number of growth blocks
!        nblk     : number of measurement blocks
!        nblkstps : number of steps in a block

integer(k4b)::nseed(4),ISEED

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!GetMtrxs(indexs)!##############################!

integer::nsites,nsites_cu,nsites_o

integer::ixx,iyy
integer::itmp,itmpx,itmpy
integer::itmpv(NSTATES)
real(sp)::tmpv(NSTATES)
real(sp)::tmpx,tmpy

integer::ixv(NSTATES),iyv(NSTATES)
integer::iposit(-lx:2*lx,-ly:2*ly),idis(lxy,2)
real(sp)::hub_u(NSTATES),epsil(NSTATES)

real(sp)::tk(NSTATES,NSTATES,NSPIN)

real(sp)::deltaf(NSTATES,NSTATES)

integer::icorx(NSTATES),icory(NSTATES),ipx(4),ipy(4),ippx(4),ippy(4),idx(4),idy(4),ncor(lx,ly) 

!        tk       : single-electron matrix (electronc kinetic)
!        hub_u    : Hubbard U at each site
!        epsil    : Chemical potential at each site
!        ixv      : x coordinates of each site
!        iyv      : y coordinates of each site
!        iposit   : (x,y) coordinates ==> site label

!GetPhit!########################################!

real(sp)::phiT_up(NSTATES,NUP,NWFBAS),phiT_dn(NSTATES,NDN,NWFBAS)
real(sp)::phiZ_up(NUP,NSTATES,NWFBAS),phiZ_dn(NDN,NSTATES,NWFBAS)
real(sp)::phiB_up(NSTATES,NUP,NIW),phiB_dn(NSTATES,NDN,NIW)
real(sp)::cwfbas(NWFBAS),cwibas(NIW),ranGenNYU
real(sp)::coeff(NSTATES,-1:1),coeffv(NSTATES,-1:1)
!-----------Zhongbing-phonon fields---------------!
real(sp)::coeffph(NSTATES,-2:2)
!-------------------------------------------------!

!        phiT_up  : trial wave function of a max NWFBAS Slater determinants
!        phiT_dn  : trial wave function of a max NWFBAS Slater determinants
!        cwfbas   : coefficients of phiT

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!Step!*******************************************!

real(sp)::phi_up(NSTATES,NUP,NWLKRS),phi_dn(NSTATES,NDN,NWLKRS)
real(sp)::wgtwlkr(NWLKRS),sgn(NWLKRS),ovlpDET(NWLKRS)
real(sp)::g_up(NUP,NUP,NWFBAS),g_dn(NDN,NDN,NWFBAS)
real(sp)::detbas(NWFBAS)
real(sp)::phi_cup(NSTATES,NUP,NWLKRS)
real(sp)::phi_cdn(NSTATES,NDN,NWLKRS)
real(sp)::wgt_c(NWLKRS)

!        phi_up   : rectangular matrix representing each walker with spin
!        wgtwlkr  : weight (multiplicity) of walker
!        ovlpDET  : overlap determinant of new to current phi
!        sgn      : sign of the overlap of phiT with initial phi's
!        g_up(dn) : Green's function for each walker with spin
!        phi_cup  : saved phi_up; for correlations

!---------------Zhongbing-phonon fields--------------------!
integer::ising_u(NSTATES),ising_v(lxy,2,4)
integer::ising_ph(NSTATES,NWLKRS),ising_sum(NSTATES)
integer::kexpV_s(NSTATES,M_MAX,NWLKRS)
integer::keVpd_s(lxy,M_MAX,NWLKRS,2,4)
integer::keph_s(NSTATES,M_MAX,NWLKRS)

real(sp)::expV(NSTATES,-1:1,NSPIN),DeltaV(NSTATES,-1:1,NSPIN)
real(sp)::eVpd(lxy,-1:1,NSPIN),DeVpd(lxy,-1:1,NSPIN)
real(sp)::eph(NSTATES,-2:2,NSPIN),Deltaph(NSTATES,-2:2,NSPIN)
!-----------------------------------------------------------!
!---------------------HoKinAnHop----------------------------!
real(sp)::expT(NSTATES,NSTATES,NSPIN),exp2T(NSTATES,NSTATES,NSPIN)
!---------------------HoKinAnHop----------------------------!
integer::mstep

integer::nstp,nstpcor,measl,mgreen
integer::nforward,nbakward
real(sp)::ovlpINIT(NIW),sgnINIT(NIW)
real(sp)::top,hbin,btm,hbinsq,prd,dens
real(sp)::cor_top(0:NAVE,6)
real(sp)::estptop,estpbtm

! real(sp)::super_xo((ntypes-10)*NSTATES*NSTATES,3),superk_xo((ntypes-10)*NSTATES,3)
! real(sp)::super_xd((ntypes-10)*NSTATES*NSTATES,4),superk_xd((ntypes-10)*NSTATES,4)



real(sp)::corptop(0:NAVE),xk(lx),yk(ly),e_var
real(sp)::sf(4),df(4),pf(4),ppf(4),ddf(4)


!------------HoKinAnHop-------------!

real(sp)::alpha,alphat1,ttp,ttn,tam

!for specify the k-point 1-kx 2-ky

real(sp)::kSet(NSTATES,2)

!------------HoKinAnHop-------------!

!------------HoKinPara-------------!
integer::noOfProc
integer::myID
!For each processor =>idStart= walker to start, idEnd=walker to end, iwCount=noOf walkers carried
!divide noOfWalkers
integer,allocatable::iwStart(:),iwEnd(:),iwCount(:)

!size of the matrix(For transferring process)
integer::sizekexpV !NSTATES*M_MAX
integer::sizePhiUp   !NSTATES*NUP
integer::sizePhiDn   !NSTATES*NDN
!------------HoKinPara-------------!
real(sp)::cor_final(NAVE)

type :: wave_r
real(sp), pointer :: arr(:,:)
end type wave_r

type :: wave_k
real(sp), pointer :: arr(:)
end type wave_k

type(wave_r), allocatable :: rVals(:)
real(sp), allocatable :: super_xo(:,:), super_xd(:,:)

type(wave_k), allocatable :: kVals(:)
real(sp), allocatable :: superk_xo(:,:), superk_xd(:,:)


real(sp)::ovlpINV_up(NUP,NUP),tmp_up(NSTATES,NUP),detp_up(NWFBAS)
real(sp)::ovlpINV_dn(NDN,NDN),tmp_dn(NSTATES,NDN),detp_dn(NWFBAS)
real(sp)::opn_up(NSTATES),opn_dn(NSTATES)
real(sp)::up_up,dn_dn,up_dn,hkin,hpot,have

real(sp), target::nk_up(NSTATES),nk_dn(NSTATES),nk_w(NSTATES)
real(sp), target::cdw(NSTATES),sdwz(NSTATES),sdwx(NSTATES)
real(sp), target::pmdf(NSTATES),dsf(NSTATES)

real(sp), target::gx_up(NSTATES,NSTATES),gx_dn(NSTATES,NSTATES)
real(sp), target::n_up(NSTATES,NSTATES),n_dn(NSTATES,NSTATES),n_w(NSTATES,NSTATES)
real(sp), target::cdw_re(NSTATES,NSTATES),sdwz_re(NSTATES,NSTATES),sdwx_re(NSTATES,NSTATES)
real(sp), target::pmdf_re(NSTATES,NSTATES),dsf_re(NSTATES,NSTATES)
real(sp), target::gxk_up(NSTATES),gxk_dn(NSTATES)

real(sp), target :: p_sdp, puu_sdp, pdd_sdp, pdsfbd, p_ddp, pdsfb, pbdx2y2
real(sp), target :: swave_re(NSTATES,NSTATES), dwave_re(NSTATES,NSTATES), pwave_re(NSTATES,NSTATES), sowave_re(NSTATES,NSTATES)
real(sp), target :: sbwave_re(NSTATES,NSTATES), dbwave_re(NSTATES,NSTATES), pbwave_re(NSTATES,NSTATES)
real(sp), target :: puupxwave_re(NSTATES,NSTATES), pddpywave_re(NSTATES,NSTATES), pudpxwave_re(NSTATES,NSTATES), puupywave_re(NSTATES,NSTATES)
real(sp), target :: pddpxwave_re(NSTATES,NSTATES), pudpywave_re(NSTATES,NSTATES)
real(sp), target :: pdsfbd1_re(NSTATES,NSTATES), pdsfbd2_re(NSTATES,NSTATES), pdsfbd12_re(NSTATES,NSTATES)
real(sp), target :: dd1wave_re(NSTATES,NSTATES), dd2wave_re(NSTATES,NSTATES), dd12wave_re(NSTATES,NSTATES)
real(sp), target :: pbdx2y2d1_re(NSTATES,NSTATES), pbdx2y2d2_re(NSTATES,NSTATES), pbdx2y2d12_re(NSTATES,NSTATES)

real(sp), target :: swave(NSTATES), dwave(NSTATES), pwave(NSTATES), sowave(NSTATES)
real(sp), target :: sbwave(NSTATES), dbwave(NSTATES), pbwave(NSTATES)
real(sp), target :: puupxwave(NSTATES), pddpywave(NSTATES), puupywave(NSTATES), pddpxwave(NSTATES)
real(sp), target :: pudpxwave(NSTATES), pudpywave(NSTATES)
real(sp), target :: pdsfbd1(NSTATES), pdsfbd2(NSTATES), pdsfbd12(NSTATES)
real(sp), target :: dd1wave(NSTATES), dd2wave(NSTATES), dd12wave(NSTATES)
real(sp), target :: pbdx2y2d1(NSTATES), pbdx2y2d2(NSTATES), pbdx2y2d12(NSTATES)
real(sp):: kValsArr(ntypes * NSTATES)
real(sp):: rValsArr(ntypes * NSTATES* NSTATES)

real(sp)::z2_sum,cu_up,cu_dn

real(sp)::wgt_p(0:NAVE), ave(NWFBAS,0:NAVE)
real(sp)::corwlkr(0:NAVE)
real(sp)::ave_cor(NAVE)
real(sp)::cdwAvg
real(sp)::sum_up, sum_dn
end module cpmc
    