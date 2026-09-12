# Decouple the k-point count from the site count in the Lieb port.
#
# The group's code assumes #k-points == #sites == NSTATES, which holds on any
# lattice with one site per cell. The Lieb lattice has three sites per cell and
# we build it on a fine grid of spacing 1/2, so the site count is 3*(lx/2)*(ly/2)
# while the natural k mesh is still the fine grid's lx*ly points. Every array and
# stride that belongs to the k side therefore becomes NKPTS = lx*ly, and only the
# site side keeps NSTATES.
#
# This also repairs a genuine indexing bug carried by the published source: the
# dir-rVals writer omits the ntypes*NKPTS offset that separates the k block from
# the r block inside cor_final, so its "real space" files actually begin with the
# k-space values. Verified on a finished checkerboard run: the first rows of
# dir-rVals/gx_up.dat reproduce dir-kVals/gx_up.dat exactly. None of our own
# analysis reads dir-rVals, so no published checkerboard number is affected.
#
# Run from the directory holding mc2duph.f90 and parameter.f90.
import re, sys

mc = open("mc2duph.f90").read()
pa = open("parameter.f90").read()
lb = open("libuph.f90").read()
orig_mc, orig_pa = mc, pa

def sub1(s, old, new, tag):
    n = s.count(old)
    if n != 1:
        sys.exit(f"[{tag}] expected 1 occurrence, found {n}")
    return s.replace(old, new)

# ---------------- parameter.f90 ----------------
# NKPTS must be declared after lx and ly, not before them
pa = sub1(pa, "parameter (lx=4,ly=4,lxy=lx*ly) !gedian",
          "parameter (lx=4,ly=4,lxy=lx*ly) !gedian\n"
          "! LIEB: k points live on the fine grid, sites do not. Keep them apart.\n"
          "integer, parameter :: NKPTS = lx*ly", "NKPTS")

pa = sub1(pa, "parameter (NAVE=(ntypes*NSTATES+ntypes)*NSTATES+6)",
          "parameter (NAVE=ntypes*NSTATES*NSTATES+ntypes*NKPTS+6)", "NAVE")

pa = sub1(pa, "real(sp)::kSet(NSTATES,2)", "real(sp)::kSet(NKPTS,2)", "kSet")

for old, tag in [
    ("real(sp), target::nk_up(NSTATES),nk_dn(NSTATES),nk_w(NSTATES)", "nk"),
    ("real(sp), target::cdw(NSTATES),sdwz(NSTATES),sdwx(NSTATES)", "cdw"),
    ("real(sp), target::pmdf(NSTATES),dsf(NSTATES)", "pmdf"),
    ("real(sp), target::gxk_up(NSTATES),gxk_dn(NSTATES)", "gxk"),
    ("real(sp), target :: swave(NSTATES), dwave(NSTATES), pwave(NSTATES), sowave(NSTATES)", "swave"),
    ("real(sp), target :: sbwave(NSTATES), dbwave(NSTATES), pbwave(NSTATES)", "sbwave"),
    ("real(sp), target :: puupxwave(NSTATES), pddpywave(NSTATES), puupywave(NSTATES), pddpxwave(NSTATES)", "puupx"),
    ("real(sp), target :: pudpxwave(NSTATES), pudpywave(NSTATES)", "pudpx"),
    ("real(sp), target :: pdsfbd1(NSTATES), pdsfbd2(NSTATES), pdsfbd12(NSTATES)", "pdsfbd"),
    ("real(sp), target :: dd1wave(NSTATES), dd2wave(NSTATES), dd12wave(NSTATES)", "ddwave"),
    ("real(sp), target :: pbdx2y2d1(NSTATES), pbdx2y2d2(NSTATES), pbdx2y2d12(NSTATES)", "pbdx2y2"),
    ("real(sp):: kValsArr(ntypes * NSTATES)", "kValsArr"),
]:
    pa = sub1(pa, old, old.replace("NSTATES", "NKPTS"), tag)

# ---------------- mc2duph.f90 ----------------
# offset past the k block, in BlkMeas and in the measurement packer
mc = sub1(mc, "k=ntypes * NSTATES", "k=ntypes * NKPTS", "blk offset")
mc = sub1(mc, "k=NSTATES*ntypes", "k=NKPTS*ntypes", "meas offset")
mc = sub1(mc, "m=ntypes*NSTATES+10*NSTATES*NSTATES",
              "m=ntypes*NKPTS+10*NSTATES*NSTATES", "super_xd offset")
mc = sub1(mc, "m=10*NSTATES", "m=10*NKPTS", "superk_xd offset")

# the two k-space accumulator blocks: loop bound and stride
mc = sub1(mc, "k=0\ndo m=1,NSTATES\n    k=k+1\n    superk_xo(",
              "k=0\ndo m=1,NKPTS\n    k=k+1\n    superk_xo(", "superk loop")
mc = sub1(mc, "k=0\ndo m=1,NSTATES\n    k=k+1\n    ave_cor(",
              "k=0\ndo m=1,NKPTS\n    k=k+1\n    ave_cor(", "ave_cor loop")
mc, n = re.subn(r"superk_xo\(k\+NSTATES\*", "superk_xo(k+NKPTS*", mc)
assert n == 22, f"superk_xo strides: {n}"
mc, n = re.subn(r"ave_cor\(k\+NSTATES\*(\d+)\)", r"ave_cor(k+NKPTS*\1)", mc)
assert n == 32, f"ave_cor k strides: {n}"

# superk_xo / superk_xd length, in the allocate, the block loops and the scratch
mc, n = re.subn(r"\(ntypes-10\)\*NSTATES(?!\*NSTATES)", "(ntypes-10)*NKPTS", mc)
assert n == 8, f"superk lengths: {n}"

# electron counts read the two gxk blocks, not the sites
mc = sub1(mc, "do j=1,NSTATES\nnumberup=numberup+cor_final(j)\n"
              "numberdn=numberdn+cor_final(nsites+j)",
              "do j=1,NKPTS\nnumberup=numberup+cor_final(j)\n"
              "numberdn=numberdn+cor_final(NKPTS+j)", "numberup")

# writePairDat: k side
mc = sub1(mc, "do i = 1, NSTATES\nk = (j-1)*NSTATES + i",
              "do i = 1, NKPTS\nk = (j-1)*NKPTS + i", "writePairDat k")
# writePairDat: r side, with the missing block offset restored
mc = sub1(mc, "k = (j-1)*NSTATES*NSTATES + (m-1)*NSTATES + n",
              "k = ntypes*NKPTS + (j-1)*NSTATES*NSTATES + (m-1)*NSTATES + n",
              "writePairDat r")

# FourierTransform: lx*ly outputs, still normalised by the number of sites
mc = sub1(mc, "real(sp), intent(inout) :: kvalue(NSTATES)",
              "real(sp), intent(inout) :: kvalue(NKPTS)", "kvalue")
mc = sub1(mc, "do iu = 1, NSTATES\n    kvalue(iu) = kvalue(iu) / real(NSTATES, sp)",
              "do iu = 1, NKPTS\n    kvalue(iu) = kvalue(iu) / real(NSTATES, sp)", "FT norm")

# icorx/icory/ncor enumerate the fine grid, not the sites. They are written and
# never read, but the consistency check at the end of the loop aborts the run.
pa = sub1(pa, "integer::icorx(NSTATES),icory(NSTATES),",
              "integer::icorx(NKPTS),icory(NKPTS),", "icorx")
mc = sub1(mc, "if(k/=NSTATES) stop 'PAIR NUMBER ERROR!'",
              "if(k/=NKPTS) stop 'PAIR NUMBER ERROR!'", "pair number")

# "per N" must divide by the number of sites, which is no longer the grid size.
mc, n = re.subn(r"/lxy", "/real(nsites)", mc)
assert n == 3, f"per-N divisors: {n}"

open("mc2duph.f90", "w").write(mc)
open("parameter.f90", "w").write(pa)
# the startup banners divide by the grid size too
lb, n = re.subn(r"etrial/lxy", "etrial/real(nsites)", lb)
assert n == 2, f"banner divisors: {n}"
open("libuph.f90", "w").write(lb)
print(f"patched: mc2duph.f90 {len(orig_mc)} -> {len(mc)} bytes, "
      f"parameter.f90 {len(orig_pa)} -> {len(pa)} bytes")
