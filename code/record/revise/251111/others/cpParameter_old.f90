module cpmc
implicit none

integer,parameter::k4b=selected_int_kind(9)
integer,parameter::i4b=selected_int_kind(9)
integer,parameter::i2b=selected_int_kind(2)


! parameter (NAVE=(ntypes*nsites+ntypes)*nsites+6)
!        NUP      : number of up   electrons; Assumed >= NDN
!        NDN      : number of down electrons
!        NE    : number of up + down electrons
!        NSPIN    : number of possible spins
!        NWFBAS   : number of determinants in initial wavefunction
!        NWLKRS   : number of walkers
!        M_MAX    : maximum number of back-propagation
!        NAVE     : number of correlations

! symbolic names for kind types of
!  single- and double-precision reals:

integer,parameter::sp=kind(1.0d0)
integer,parameter::dp=kind(1.0d0)

real(sp), parameter :: PI = 3.14159265359
real(sp), parameter :: ZERO=0.0, FOURTH=0.25, HALF=0.5
real(sp), parameter :: ONE=1.0, TWO=2.0, FOUR=4.0, TEN=10.0

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
integer, parameter :: lx=4, ly=4, NLA=2, lxy=lx*ly
integer, parameter :: nsites = lxy * NLA
integer, parameter :: NSTATES = lxy
integer, parameter :: NHS = NSTATES*NSTATES

!        NSTATES  : number of states in SCF solution
!        NHS      : number of HS fields

integer, parameter :: NUP=24, NDN=24, NSPIN=2

! ===== Spin-mixed (spinor) representation for Fock/SOC =====
integer, parameter :: NSO = 2*nsites          ! spin-orbital dimension
integer, parameter :: NE  = NUP + NDN         ! total electrons
logical :: use_spinor = .false.               ! set .true. to enable spin-mixed trial/walkers
real(sp), parameter :: dens = real(NUP+NDN)/real(nsites)
! Combined spinor Slater structures (allocated at runtime when use_spinor=.true.)
! yes spinor
real(sp), allocatable :: phi(:,:,:), phiT(:,:,:), phiZ(:,:,:), phiB(:,:,:)
real(sp), allocatable :: phi_c(:,:,:), g(:,:,:), gx(:,:)
real(sp), allocatable :: ovlpINV(:,:), tmp_ud(:,:), detp(:)

! no spinor
real(sp), allocatable :: phi_up(:,:,:), phi_dn(:,:,:)
real(sp), allocatable :: phiT_up(:,:,:), phiT_dn(:,:,:)
real(sp), allocatable :: phiZ_up(:,:,:), phiZ_dn(:,:,:)
real(sp), allocatable :: phiB_up(:,:,:), phiB_dn(:,:,:)
real(sp), allocatable :: phi_cup(:,:,:), phi_cdn(:,:,:)
real(sp), allocatable :: g_up(:,:,:), g_dn(:,:,:)
real(sp), allocatable :: gx_up(:,:), gx_dn(:,:)

real(sp), allocatable :: ovlpINV_up(:,:), tmp_up(:,:), detp_up(:)
real(sp), allocatable :: ovlpINV_dn(:,:), tmp_dn(:,:), detp_dn(:)

! different dimensions depends on spinor
real(sp), allocatable :: tk(:,:,:)
real(sp), allocatable :: expT(:,:,:), exp2T(:,:,:)

!same dimensions
real(sp)::expV(nsites,-1:1,NSPIN),DeltaV(nsites,-1:1,NSPIN)
real(sp)::eVpd(nsites,-1:1,NSPIN),DeVpd(nsites,-1:1,NSPIN)
real(sp)::eph(nsites,-2:2,NSPIN),Deltaph(nsites,-2:2,NSPIN)
real(sp):: gxk_up(NSTATES,NLA,NLA), gxk_dn(NSTATES,NLA,NLA)


integer, parameter :: NWFBAS=1, NIW=1, NWLKRS=1000, M_MAX=125
integer, parameter :: ntypes = 30
integer, parameter :: NAVE = ntypes*NSTATES + ntypes*NSTATES*NSTATES + 6



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

integer::ispin,ising

!SetUp

real(sp)::deltau,etrial,tpd,tpp,t1,t2,udA,udB,vpd,g_ph,w_ph

!        deltau   : delta tau, the Trotter step size
!        etrial   : growth control factor
!        tpd       : t_0A, hopping integral between NN sites for lattice A
!        tpp       : t_0B, hopping integral between NN sites for lattice B
!        t1       : t_1, hopping integral between NNN sites
!        t2       : t_2, hopping integral between NNNN sites
!        udA       : U_d, Hubbard U on each site for O
!        udB       : U_d, Hubbard U on each site for Cu
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


integer::ixx,iyy,i_a,i_b
integer::itmp,itmpx,itmpy
integer::itmpv(nsites)
real(sp)::tmpv(nsites)
real(sp)::tmpx,tmpy

integer::ixv(nsites),iyv(nsites),sublatt(nsites)
integer::iposit(-lx:2*lx,-ly:2*ly,NLA),idis(nsites,2)
real(sp)::hub_u(nsites),epsil(nsites)



real(sp)::deltaf(nsites,nsites)

integer::icorx(nsites),icory(nsites),ipx(4),ipy(4),ippx(4),ippy(4),idx(4),idy(4),ncor(lx,ly)

!        tk       : single-electron matrix (electronc kinetic)
!        hub_u    : Hubbard U at each site
!        epsil    : Chemical potential at each site
!        ixv      : x coordinates of each site
!        iyv      : y coordinates of each site
!        iposit   : (x,y) coordinates ==> site label

!GetPhit!########################################!


real(sp)::cwfbas(NWFBAS),cwibas(NIW),ranGenNYU
real(sp)::coeff(nsites,-1:1),coeffv(nsites,-1:1)
!-----------Zhongbing-phonon fields---------------!
real(sp)::coeffph(nsites,-2:2)
!-------------------------------------------------!

!        phiT_up  : trial wave function of a max NWFBAS Slater determinants
!        phiT_dn  : trial wave function of a max NWFBAS Slater determinants
!        cwfbas   : coefficients of phiT

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!Step!*******************************************!


real(sp)::wgtwlkr(NWLKRS),sgn(NWLKRS),ovlpDET(NWLKRS)

real(sp)::detbas(NWFBAS)

real(sp)::wgt_c(NWLKRS)

!        phi_up   : rectangular matrix representing each walker with spin
!        wgtwlkr  : weight (multiplicity) of walker
!        ovlpDET  : overlap determinant of new to current phi
!        sgn      : sign of the overlap of phiT with initial phi's
!        g_up(dn) : Green's function for each walker with spin
!        phi_cup  : saved phi_up; for correlations

!---------------Zhongbing-phonon fields--------------------!
integer::ising_u(nsites),ising_v(nsites,2,4)
integer::ising_ph(nsites,NWLKRS),ising_sum(nsites)
integer::kexpV_s(nsites,M_MAX,NWLKRS)
integer::keVpd_s(nsites,M_MAX,NWLKRS,2,4)
integer::keph_s(nsites,M_MAX,NWLKRS)


integer::mstep

integer::nstp,nstpcor,measl,mgreen
integer::nforward,nbakward
real(sp)::ovlpINIT(NIW),sgnINIT(NIW)
real(sp)::top,hbin,btm,hbinsq,prd
real(sp)::cor_top(0:NAVE,6,NLA,NLA)
real(sp)::estptop,estpbtm

! real(sp)::super_xo((ntypes-10)*nsites*nsites,3),superk_xo((ntypes-10)*nsites,3)
! real(sp)::super_xd((ntypes-10)*nsites*nsites,4),superk_xd((ntypes-10)*nsites,4)

real(sp)::corptop(0:NAVE,NLA,NLA),xk(lx),yk(ly),e_var
real(sp)::sf(4),df(4),pf(4),ppf(4),ddf(4)


real(sp)::alpha,alphat1,ttp,ttn,tam

!for specify the k-point 1-kx 2-ky

real(sp)::kSet(NSTATES,2)


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
integer:: sizePhi !Nsites*NE
!------------HoKinPara-------------!
real(sp)::cor_final(NAVE,NLA,NLA)

! type(wave_r), allocatable :: rVals(:)
real(sp):: super_xo((ntypes-8)*NSTATES*NSTATES,3,NLA,NLA)
real(sp):: super_xd((ntypes-8)*NSTATES*NSTATES,4,NLA,NLA)

! type(wave_k), allocatable :: kVals(:)
real(sp):: superk_xo((ntypes-8)*NSTATES,3,NLA,NLA)
real(sp):: superk_xd((ntypes-8)*NSTATES,4,NLA,NLA)


character(len=100) :: filename

real(sp)::opn_up(nsites),opn_dn(nsites)
real(sp)::up_up,dn_dn,up_dn,hkin,hpot,have,error,energy

real(sp)::nk_up(NSTATES,NLA,NLA),nk_dn(NSTATES,NLA,NLA),nk_w(NSTATES,NLA,NLA)
real(sp)::cdw(NSTATES,NLA,NLA),sdwz(NSTATES,NLA,NLA),sdwx(NSTATES,NLA,NLA)
real(sp)::pmdf(NSTATES,NLA,NLA),dsf(NSTATES,NLA,NLA)


real(sp)::n_up(nsites,nsites),n_dn(nsites,nsites),n_w(nsites,nsites)
real(sp)::cdw_re(nsites,nsites),sdwz_re(nsites,nsites),sdwx_re(nsites,nsites)
real(sp)::pmdf_re(nsites,nsites),dsf_re(nsites,nsites)


real(sp):: p_sdp, puu_sdp, pdd_sdp, pdsfbd, p_ddp, pdsfb, pbdx2y2
real(sp):: swave_re(nsites,nsites), dwave_re(nsites,nsites), pwave_re(nsites,nsites), sowave_re(nsites,nsites)
real(sp):: sbwave_re(nsites,nsites), dbwave_re(nsites,nsites), pbwave_re(nsites,nsites)
real(sp):: puupxwave_re(nsites,nsites), pddpywave_re(nsites,nsites), pudpxwave_re(nsites,nsites), puupywave_re(nsites,nsites)
real(sp):: pddpxwave_re(nsites,nsites), pudpywave_re(nsites,nsites)
real(sp):: pdsfbd1_re(nsites,nsites), pdsfbd2_re(nsites,nsites), pdsfbd12_re(nsites,nsites)
real(sp):: dd1wave_re(nsites,nsites), dd2wave_re(nsites,nsites), dd12wave_re(nsites,nsites)
real(sp):: pbdx2y2d1_re(nsites,nsites), pbdx2y2d2_re(nsites,nsites), pbdx2y2d12_re(nsites,nsites)

real(sp):: swave(NSTATES,NLA,NLA), dwave(NSTATES,NLA,NLA), pwave(NSTATES,NLA,NLA), sowave(NSTATES,NLA,NLA)
real(sp):: sbwave(NSTATES,NLA,NLA), dbwave(NSTATES,NLA,NLA), pbwave(NSTATES,NLA,NLA)
real(sp):: puupxwave(NSTATES,NLA,NLA), pddpywave(NSTATES,NLA,NLA), puupywave(NSTATES,NLA,NLA), pddpxwave(NSTATES,NLA,NLA)
real(sp):: pudpxwave(NSTATES,NLA,NLA), pudpywave(NSTATES,NLA,NLA)
real(sp):: pdsfbd1(NSTATES,NLA,NLA), pdsfbd2(NSTATES,NLA,NLA), pdsfbd12(NSTATES,NLA,NLA)
real(sp):: dd1wave(NSTATES,NLA,NLA), dd2wave(NSTATES,NLA,NLA), dd12wave(NSTATES,NLA,NLA)
real(sp):: pbdx2y2d1(NSTATES,NLA,NLA), pbdx2y2d2(NSTATES,NLA,NLA), pbdx2y2d12(NSTATES,NLA,NLA)
real(sp)::numberup,numberdn
real(sp)::z2_sum,cu_up,cu_dn

real(sp)::wgt_p(0:NAVE,NLA,NLA), ave(NWFBAS,0:NAVE,NLA,NLA)
real(sp)::corwlkr(0:NAVE,NLA,NLA)
real(sp)::ave_cor(NAVE,NLA,NLA)
real(sp)::cdwAvg
real(sp)::sum_up, sum_dn
end module cpmc

subroutine spinor_alloc
  use cpmc
  implicit none
  if (use_spinor) then
     ! Allocate only spinor (mixed-spin) structures
     allocate(phi(NSO,NE,NWLKRS))
     allocate(phiT(NSO,NE,NWFBAS))
     allocate(phiZ(NE,NSO,NWFBAS))
     allocate(phiB(NSO,NE,NIW))
     allocate(phi_cup(NSO,NWLKRS))
     allocate(g(NE,NE,NWFBAS))
     allocate(gx(NSO,NSO))

     allocate(tk(NSO,NSO,1))
     allocate(expT(NSO,NSO,1))
     allocate(exp2T(NSO,NSO,1))

     allocate(ovlpINV(NE,NE)),
     allocate(tmp_ud(NSO,NE))
     allocate(detp(NWFBAS))
  else
     ! Allocate only legacy per-spin structures
     allocate(phi_up(nsites,NUP,NWLKRS))
     allocate(phi_dn(nsites,NDN,NWLKRS))
     allocate(phiT_up(nsites,NUP,NWFBAS))
     allocate(phiT_dn(nsites,NDN,NWFBAS))
     allocate(phiZ_up(NUP,nsites,NWFBAS))
     allocate(phiZ_dn(NDN,nsites,NWFBAS))
     allocate(phiB_up(nsites,NUP,NIW))
     allocate(phiB_dn(nsites,NDN,NIW))
     allocate(phi_cup(nsites,NUP,NWLKRS))
     allocate(phi_cdn(nsites,NDN,NWLKRS))

     allocate(g_up(NUP,NUP,NWFBAS))
     allocate(g_dn(NDN,NDN,NWFBAS))

     allocate(tk(nsites,nsites,NSPIN))
     allocate(expT(nsites,nsites,NSPIN))
     allocate(exp2T(nsites,nsites,NSPIN))

     allocate(gx_up(nsites,nsites))
     allocate(gx_dn(nsites,nsites))

     allocate(ovlpINV_up(NUP,NUP))
     allocate(ovlpINV_dn(NDN,NDN))
     allocate(tmp_up(nsites,NUP))
     allocate(tmp_dn(nsites,NDN))
     allocate(detp_up(NWFBAS))
     allocate(detp_dn(NWFBAS))

  end if
endsubroutine spinor_alloc
