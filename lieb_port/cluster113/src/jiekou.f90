module jiekou

interface
function ranGen(idum)
use cpmc
implicit none
real(sp)::ranGen
integer(k4b),intent(inout)::idum
end function ranGen

subroutine stblz
use cpmc
implicit none
end subroutine stblz

subroutine stblz1
use cpmc
implicit none
end subroutine stblz1

subroutine modgs(phi,nss,npp,rescale)
use cpmc
implicit none
real(sp)::rescale,phi(:,:)
integer::nss,npp
end subroutine modgs

SUBROUTINE SGEFA(A,N,N1,IPVT,INFO)
use cpmc
implicit none
integer::N,N1,INFO
integer,dimension(:)::IPVT
REAL(sp),dimension(:,:)::A
end subroutine sgefa

SUBROUTINE SGEDI(A,N,N1,IPVT,DET,WORK,JOB)
use cpmc
implicit none
integer::N,N1,JOB
integer,dimension(:)::IPVT
REAL(sp),dimension(:)::DET,WORK
real(sp),dimension(:,:)::A
end subroutine sgedi

FUNCTION SNRM2(N,SX,INCX)
use cpmc
implicit none
real(sp)::SNRM2
real(sp),dimension(:)::SX
integer::N,INCX
END FUNCTION SNRM2

subroutine tred2(nm,n,a,d,e,z)            ! tred2
use cpmc
implicit none
integer::n,nm
real(sp),dimension(:,:)::a,z
real(sp),dimension(:)::d,e
end subroutine tred2

SUBROUTINE TQL2(NM,N,D,E,Z,IERR)          ! tql2
use cpmc
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

SUBROUTINE GetMtrxs
use cpmc
end subroutine GetMtrxs

SUBROUTINE MkExpT(xdeltau)
use cpmc
real(sp)::xdeltau
end subroutine MkExpT

SUBROUTINE MkInitOvlps
use cpmc
end subroutine MkInitOvlps

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

SUBROUTINE Step(istp)
use cpmc
integer::istp
end subroutine Step

SUBROUTINE Step1(istp)
use cpmc
integer::istp
end subroutine Step1

SUBROUTINE HalfT(iw,icase)
use cpmc
integer::iw,icase
end subroutine HalfT

SUBROUTINE Ovlps(iw)
use cpmc
integer::iw
end subroutine Ovlps

SUBROUTINE Fullph(iw)
use cpmc
integer::iw
end subroutine Fullph

SUBROUTINE FullV(iw)
use cpmc
integer::iw
end subroutine FullV

SUBROUTINE preprog
use cpmc
end SUBROUTINE preprog

SUBROUTINE MkExpV
use cpmc
end subroutine MkExpV

SUBROUTINE Coupling(x,a)
use cpmc
implicit none
real(sp)::x,a
end subroutine Coupling

subroutine savephi
use cpmc
end subroutine savephi

SUBROUTINE EstEtrial
use cpmc
end subroutine EstEtrial

SUBROUTINE InitEnergy
use cpmc
end subroutine InitEnergy

SUBROUTINE StpMeas
use cpmc
end subroutine StpMeas

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
real(sp), intent(in)::gp_up(:,:),gp_dn(:,:)
end subroutine pair

subroutine FourierTransform(rval,kvalue)
use cpmc
real(sp), intent(in) :: rval(:,:)
real(sp), allocatable :: kvalue(:)
end subroutine FourierTransform

subroutine backphi(kp,L,phip_up,phip_dn)
use cpmc
integer::kp,L
real(sp)::phip_up(NUP,NSTATES),phip_dn(NDN,NSTATES)
end subroutine backphi

subroutine stblzbk(kp,phi1,phi2)
use cpmc
implicit none
integer::kp
real(sp),dimension(NUP,NSTATES)::phi1
real(sp),dimension(NDN,NSTATES)::phi2
end subroutine stblzbk


end interface

end module jiekou
    