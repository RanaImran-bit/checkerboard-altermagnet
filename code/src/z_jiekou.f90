module jiekou
interface
function ranGen(idum)
use cpmc
implicit none
real(sp)::ranGen
integer(k4b),intent(inout)::idum
end function ranGen

subroutine stblz(mode)
use cpmc
implicit none
integer::mode
end subroutine stblz

subroutine modgs(phiin,nss,npp,rescale)
use cpmc
implicit none
real(sp)::rescale,tmp,hld,sdot,dnrm2
integer::nss,npp,ip,jp,n1
real(sp),dimension(:,:)::phiin
end subroutine modgs


subroutine tred2(nm,n,a,d,e,z)            ! tred2
use cpmc, only: sp, use_spinor, nsites, NSO
implicit none
integer::n,nm
real(sp),dimension(:,:)::a,z
real(sp),dimension(:)::d,e
end subroutine tred2

SUBROUTINE TQL2(NM,N,D,E,Z,IERR)          ! tql2
use cpmc, only: sp, use_spinor, nsites, NSO
implicit none
integer::N,NM,IERR
REAL(sp),dimension(:)::D,E
real(sp),dimension(:,:)::Z
end subroutine tql2

function pythag(a,b)                       ! pythag
use cpmc
implicit none
real(sp)::pythag
real(sp)::a,b
end function pythag

SUBROUTINE SetUp
use cpmc
end subroutine SetUp

SUBROUTINE Initlatt
use cpmc
end subroutine Initlatt

SUBROUTINE InitV
use cpmc
end subroutine InitV

SUBROUTINE Initpair
use cpmc
end subroutine Initpair

SUBROUTINE Initk
use cpmc
end subroutine Initk

subroutine initpara
  use cpmc
endsubroutine initpara

subroutine initout
  use cpmc
endsubroutine initout

SUBROUTINE MkExpT(xdeltau)
use cpmc
real(sp)::xdeltau
end subroutine MkExpT

SUBROUTINE InitPop
use cpmc
end subroutine InitPop

SUBROUTINE Comb
use cpmc
end subroutine Comb

SUBROUTINE PopCopy(iwempty,iw)
use cpmc
integer::iwempty,iw
end subroutine PopCopy

SUBROUTINE Step(istp, mode)
use cpmc
integer::istp,mode
end subroutine Step

SUBROUTINE HalfT(iw,icase,istp)
use cpmc
integer::iw,icase,istp
end subroutine HalfT

SUBROUTINE Ovlps(iw)
use cpmc
integer, intent(in) :: iw
end subroutine Ovlps

SUBROUTINE Fullph(iw)
use cpmc
integer::iw
end subroutine Fullph

! SUBROUTINE Vee(iw, isite, jsite, ispin, jspin, chan)
! use cpmc
! integer::iw, isite, jsite, ispin, jspin, chan
! end subroutine Vee

SUBROUTINE Vee(iw, isite, jsite,  chan)
  use cpmc
  integer::iw, isite, jsite, chan
end subroutine Vee

SUBROUTINE VHubb(iw)
use cpmc
integer::iw
end subroutine VHubb

SUBROUTINE Vznn(iw)
use cpmc
integer::iw
end subroutine Vznn

SUBROUTINE Vxy(iw)
use cpmc
integer::iw
end subroutine Vxy

SUBROUTINE Initph
use cpmc
end SUBROUTINE Initph

SUBROUTINE MkExpV
use cpmc
end subroutine MkExpV

SUBROUTINE invcosh(x,a)
use cpmc
implicit none
real(sp)::x,a
end subroutine invcosh

subroutine savephi
use cpmc
end subroutine savephi

SUBROUTINE EstEtrial
use cpmc
end subroutine EstEtrial

SUBROUTINE InitEnergy
use cpmc
end subroutine InitEnergy

SUBROUTINE StepMeas
use cpmc
end subroutine StepMeas

SUBROUTINE InitRunMeas
use cpmc
end subroutine InitRunMeas

SUBROUTINE InitBlkMeas
use cpmc
end subroutine InitBlkMeas

SUBROUTINE RunMeas
use cpmc
end subroutine RunMeas

subroutine correl
use cpmc
end subroutine correl

subroutine pair(gp_up, gp_dn)
use cpmc
real(sp), intent(in)::gp_up(nsites,nsites), gp_dn(nsites,nsites)
end subroutine pair

subroutine FourierTransform(rval,kvalue)
use cpmc
real(sp), intent(in) :: rval(nsites,nsites)
real(sp), intent(inout) :: kvalue(NSTATES, NLA, NLA)
integer :: rxm, rym, rxn, ryn
integer :: dx(NLA), dy(NLA)
end subroutine FourierTransform

subroutine backphi(kp,L,phip_up,phip_dn)
use cpmc
integer::kp,L
real(sp)::phip_up(NUP,NSTATES),phip_dn(NDN,NSTATES)
end subroutine backphi

subroutine Scatter
  use cpmc
  use mpi
  implicit none  
end subroutine Scatter

subroutine Gather
  use cpmc
  use mpi
  implicit none  
end subroutine Gather

subroutine stblzbk(kp,phi_1,phi_2)
use cpmc
implicit none
integer::kp
real(sp),dimension(NUP,NSTATES)::phi_1
real(sp),dimension(NDN,NSTATES)::phi_2
end subroutine stblzbk

subroutine stblzbk_spinor(kp,phi_1)
	use cpmc
	implicit none
  integer::kp
  real(sp),dimension(NE,NSO)::phi_1
end subroutine stblzbk_spinor

subroutine dgefa(a,lda,n,ipvt,info)
    integer lda,n,ipvt(1),info
    double precision a(lda,1)
    double precision t
    integer idamax,j,k,kp1,l,nm1
end

subroutine dgedi(a,lda,n,ipvt,det,work,job)
    integer lda,n,ipvt(1),job
    double precision a(lda,1),det(2),work(1)
    double precision t
    double precision ten
    integer i,j,k,kb,kp1,l,nm1
end

end interface

end module jiekou
