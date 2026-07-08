module cpmc
implicit none

integer,parameter::k4b=selected_int_kind(9)
integer,parameter::i4b=selected_int_kind(9)
integer,parameter::i2b=selected_int_kind(2)
integer,parameter::sp=kind(1.0d0)
integer,parameter::dp=kind(1.0d0)
real(sp), parameter :: PI = 3.14159265359
real(sp), parameter :: ZERO=0.0d0, FOURTH=0.25, HALF=0.5d0
real(sp), parameter :: ONE=1.0, TWO=2.0, FOUR=4.0, TEN=10.0
integer, parameter :: iRead=0

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! parameter (NAVE=(ntypes*nsites+ntypes)*nsites+6)
!        NSTATES  : number of states in SCF solution
!        NHS      : number of HS fields
!        NUP      : number of up   electrons; Assumed >= NDN
!        NDN      : number of down electrons
!        NE    : number of up + down electrons
!        NSPIN    : number of possible spins
!        NWFBAS   : number of determinants in initial wavefunction
!        NWLKRS   : number of walkers
!        M_MAX    : maximum number of back-propagation
!        NAVE     : number of correlations
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
!        itvlpceq : population control frequency in relaxation phase
!        itvlpcgr : population control frequency in growth estimator
!        itvlpc   : population control frequency in measurement phase
!        itvlmeas : measurement frequency
!        itvlorth : orthogonalization frequency
!        itvl_m   : back-propagator
!        nblkeq   : number of equilibrium blocks
!        nblkgr   : number of growth blocks
!        nblk     : number of measurement blocks
!        nblkstps : number of steps in a block
!        tk       : single-electron matrix (electronc kinetic+chemical potential)
!        ixv      : x coordinates of each site
!        iyv      : y coordinates of each site
!        iposit   : (x,y) coordinates ==> site label
!        idis     : nearest neighbor site index
!        phiT_up  : trial wave function of a max NWFBAS Slater determinants
!        phiT_dn  : trial wave function of a max NWFBAS Slater determinants
!        cwfbas   : coefficients of phiT
!        phi_up   : rectangular matrix representing each walker with spin
!        wgtwlkr  : weight (multiplicity) of walker
!        ovlpDET  : overlap determinant of new to current phi
!        sgn      : sign of the overlap of phiT with initial phi's
!        g_up(dn) : Green's function for each walker with spin
!        phi_cup  : saved phi_up; for correlations
logical :: use_spinor = .false.               ! set .true. to enable spin-mixed trial/walkers
integer, parameter :: lx=4, ly=4, NLA=2, lxy=lx*ly
integer, parameter :: nsites = lxy * NLA
integer, parameter :: NSTATES = lxy
integer, parameter :: NHS = NSTATES*NSTATES
integer, parameter :: NUP=16, NDN=16, NELEC=NUP+NDN, NSPIN=2
integer, parameter :: NWFBAS=1, NIW=1, NWLKRS=1000, M_MAX=125
integer, parameter :: ntypes = 30
integer, parameter :: NAVE = ntypes*NSTATES + ntypes*NSTATES*NSTATES + 6
integer, parameter :: channels = 37
integer, parameter :: NE  = NUP + NDN         ! total electrons
real(sp), parameter :: dens = real(NUP+NDN)/real(nsites)
real(sp)::deltau,hdeltau,etrial,t1,t2,t3,t4,uxx,uxy,v,g_ph,w_ph
integer::itvlpceq,itvlpcgr,itvlpc,itvlmeas,itvlorth,itvl_m
integer::nblkeq,nblkgr,nblk,nblkstps
integer(k4b)::nseed(4),ISEED
! ===== Spin-mixed (spinor) representation for Fock/SOC =====
integer, parameter :: NSO = 2*nsites          ! spin-orbital dimension
integer::ixv(nsites),iyv(nsites),sublatt(nsites)
integer::iposit(-2*lx:2*lx,-2*ly:2*ly,NLA),idis(nsites,2)
integer::icorx(nsites),icory(nsites)
integer::ipx(4),ipy(4),ippx(4),ippy(4),idx(4),idy(4)
real(sp)::cwfbas(NWFBAS),cwibas(NIW),ranGenNYU
real(sp)::wgtwlkr(NWLKRS),sgn(NWLKRS),ovlpDET(NWLKRS)
real(sp)::detbas(NWFBAS)
real(sp)::wgt_c(NWLKRS)
real(sp)::ovlpINIT(NIW),sgnINIT(NIW)
! Combined spinor Slater structures (allocated at runtime when use_spinor=.true.)
! spinor
real(sp), allocatable :: phi(:,:,:), phiT(:,:,:), phiZ(:,:,:), phiB(:,:,:)
real(sp), allocatable :: phi_c(:,:,:), g(:,:,:), gx(:,:)
real(sp), allocatable :: ovlpINV(:,:), tmp_ud(:,:), detp(:)
real(sp), allocatable :: tk(:,:)
real(sp), allocatable :: expT(:,:), exp2T(:,:)
real(sp), allocatable :: d(:), e(:), z(:,:), T(:,:)

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
real(sp), allocatable :: tk_ud(:,:,:)
real(sp), allocatable :: expT_ud(:,:,:), exp2T_ud(:,:,:)
real(sp), allocatable :: d_ud(:,:), e_ud(:,:), z_ud(:,:,:)

!------ ee interaction term ------!
real(sp)::coeffv(nsites,-1:1,channels) ! coefficents associated with +ve or -ve U
integer::ising_v(nsites,channels)
real(sp)::expV(nsites,-1:1,NSPIN,channels),DeltaV(nsites,-1:1,NSPIN,channels)
integer::kexpV_s(nsites,M_MAX,channels,NWLKRS)
! define the n and spin associated with each channel
integer::nlsi(channels),nlsj(channels), spinlsi(channels),spinlsj(channels)
integer:: orblsi(channels), orblsj(channels)
!--------phonon fields--------------------!
integer::ising_ph(nsites,NWLKRS),ising_sum(nsites)
real(sp)::eph(nsites,-2:2,NSPIN),Deltaph(nsites,-2:2,NSPIN)
real(sp)::coeffph(nsites,-2:2)
integer::keph_s(nsites,M_MAX,NWLKRS)

!-------global parameter---------!
integer::mstep ! parameter
integer::nstp,nstpcor,measl
real(sp)::top,hbin,btm,hbinsq,prd
real(sp)::cor_top(0:NAVE,6,NLA,NLA)
real(sp)::estptop,estpbtm
real(sp)::eblk,eblksq

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
integer, parameter ::sizekexpV = nsites * M_MAX * channels !NSTATES*M_MAX
integer, parameter ::sizePhiUp = nsites * NUP  !NSTATES*NUP
integer, parameter ::sizePhiDn = nsites * NDN  !NSTATES*NDN
integer, parameter ::sizeKeph  = nsites * M_MAX
integer, parameter ::sizePhi = NSO * NE
integer, allocatable :: counts(:), displs(:)

!------------HoKinPara-------------!
real(sp)::cor_final(NAVE,NLA,NLA)

! type(wave_r), allocatable :: rVals(:)
real(sp):: super_xo((ntypes-8)*NSTATES*NSTATES,3,NLA,NLA)
real(sp):: super_xd((ntypes-8)*NSTATES*NSTATES,4,NLA,NLA)

! type(wave_k), allocatable :: kVals(:)
real(sp):: superk_xo((ntypes-8)*NSTATES,3,NLA,NLA)
real(sp):: superk_xd((ntypes-8)*NSTATES,4,NLA,NLA)


character(len=100) :: filename

! real(sp)::ovlpINV_up(NUP,NUP),tmp_up(nsites,NUP),detp_up(NWFBAS)
! real(sp)::ovlpINV_dn(NDN,NDN),tmp_dn(nsites,NDN),detp_dn(NWFBAS)
real(sp)::opn_up(nsites),opn_dn(nsites)
real(sp)::up_up,dn_dn,up_dn,hkin,hpot,have,error,energy

real(sp)::nk_up(NSTATES,NLA,NLA),nk_dn(NSTATES,NLA,NLA),nk_w(NSTATES,NLA,NLA)
real(sp)::cdw(NSTATES,NLA,NLA),sdwz(NSTATES,NLA,NLA),sdwx(NSTATES,NLA,NLA)
real(sp)::pmdf(NSTATES,NLA,NLA),dsf(NSTATES,NLA,NLA)

! real(sp)::gx_up(nsites,nsites),gx_dn(nsites,nsites)
real(sp)::n_up(nsites,nsites),n_dn(nsites,nsites),n_w(nsites,nsites)
real(sp)::cdw_re(nsites,nsites),sdwz_re(nsites,nsites),sdwx_re(nsites,nsites)
real(sp)::pmdf_re(nsites,nsites),dsf_re(nsites,nsites)
real(sp)::gxk_up(NSTATES,NLA,NLA),gxk_dn(NSTATES,NLA,NLA)

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

real(sp)::Vlist(nsites,37)

integer :: kexpV_s1(nsites, M_MAX, channels, NWLKRS)
integer :: keph_s1(nsites, M_MAX, NWLKRS)
real(sp) :: phi_up1(nsites, NUP, NWLKRS), phi_dn1(nsites, NDN, NWLKRS)
real(sp) :: phi_cup1(nsites, NUP, NWLKRS), phi_cdn1(nsites, NDN, NWLKRS)
real(sp) :: phi1(NSO, NE, NWLKRS), phi_c1(NSO, NE, NWLKRS)
real(sp) :: wgtwlkr1(NWLKRS), sgn1(NWLKRS), ovlpDET1(NWLKRS), wgt_c1(NWLKRS), etrial1

end module cpmc

subroutine allocate_vars
  use cpmc
  allocate( tk_ud(nsites,nsites,NSPIN))

  if (use_spinor) then
      ! spinor
      allocate( phi(NSO,NE,NWFBAS), phiT(NSO,NE,NWFBAS), phiZ(NSO,NE,NWFBAS), phiB(NSO,NE,NWFBAS))
      allocate( phi_c(NSO,NE,NWFBAS), g(NE,NE,NWFBAS), gx(NSO,NSO))
      allocate( ovlpINV(NE,NE), tmp_ud(NSO,NE), detp(NWFBAS))
      allocate( tk(NSO,NSO))
      allocate( expT(NSO,NSO), exp2T(NSO,NSO))
      allocate( d(NSO), e(NSO), z(NSO,NSO), T(NSO,NSO))
  else
      ! no spinor
      allocate( phiT_up(nsites,NUP,NWFBAS), phiT_dn(nsites,NDN,NWFBAS))
      allocate( phiZ_up(NUP,nsites,NWFBAS),phiZ_dn(NDN,nsites,NWFBAS))
      allocate( phiB_up(nsites,NUP,NIW),phiB_dn(nsites,NDN,NIW))
      allocate( phi_up(nsites,NUP,NWLKRS),phi_dn(nsites,NDN,NWLKRS))
      allocate( phi_cup(nsites,NUP,NWLKRS),phi_cdn(nsites,NDN,NWLKRS))
      allocate( g_up(NUP,NUP,NWFBAS),g_dn(NDN,NDN,NWFBAS))
      allocate( gx_up(nsites,nsites),gx_dn(nsites,nsites))
      allocate( ovlpINV_up(NUP,NUP),tmp_up(nsites,NUP),detp_up(NWFBAS))
      allocate( ovlpINV_dn(NDN,NDN),tmp_dn(nsites,NDN),detp_dn(NWFBAS))
      allocate( expT_ud(nsites,nsites,NSPIN),exp2T_ud(nsites,nsites,NSPIN))
      allocate( d_ud(nsites,2), e_ud(nsites,2), z_ud(nsites,nsites,2))
  end if
end subroutine
