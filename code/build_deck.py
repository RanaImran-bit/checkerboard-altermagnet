"""Code-change tracking deck: Combined_Model (spin-dependent NNN) -> Checkerboard_Model.

Layout on every code slide: original on the left, modified on the right, reason
underneath. Code is set as real text in Courier New rather than as screenshots, so
it stays selectable, scales cleanly on a projector, and can be corrected in place.
"""
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

NAVY  = RGBColor(0x1E, 0x27, 0x61)
OLD   = RGBColor(0xD5, 0x5E, 0x00)     # vermillion  = original code
NEW   = RGBColor(0x00, 0x72, 0xB2)     # blue        = modified code
INK   = RGBColor(0x21, 0x21, 0x21)
MUTE  = RGBColor(0x6B, 0x6B, 0x6B)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
OLDBG = RGBColor(0xFD, 0xF1, 0xE8)
NEWBG = RGBColor(0xEA, 0xF3, 0xF9)
PANEL = RGBColor(0xF4, 0xF5, 0xF7)

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BLANK = prs.slide_layouts[6]

def tb(slide, x, y, w, h, text, size=14, bold=False, color=INK, font="Calibri",
       align=PP_ALIGN.LEFT, fill=None, line=None, anchor=MSO_ANCHOR.TOP, space=0):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = box.text_frame; tf.word_wrap = True; tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Inches(0.10)
    tf.margin_top = tf.margin_bottom = Inches(0.06)
    for i, ln in enumerate(str(text).split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(space)
        r = p.add_run(); r.text = ln
        r.font.size = Pt(size); r.font.bold = bold
        r.font.color.rgb = color; r.font.name = font
    if fill is not None:
        box.fill.solid(); box.fill.fore_color.rgb = fill
    else:
        box.fill.background()
    if line is not None:
        box.line.color.rgb = line; box.line.width = Pt(1.0)
    else:
        box.line.fill.background()
    return box

def title(slide, t, sub=None):
    tb(slide, 0.45, 0.28, 12.4, 0.62, t, size=30, bold=True, color=NAVY)
    if sub:
        tb(slide, 0.45, 0.92, 12.4, 0.34, sub, size=12, color=MUTE, font="Courier New")

def codeslide(t, sub, lcap, lcode, rcap, rcode, why, csize=8.5):
    s = prs.slides.add_slide(BLANK)
    title(s, t, sub)
    tb(s, 0.45, 1.28, 6.10, 0.34, lcap, size=12, bold=True, color=OLD)
    tb(s, 6.78, 1.28, 6.10, 0.34, rcap, size=12, bold=True, color=NEW)
    tb(s, 0.45, 1.62, 6.10, 3.28, lcode, size=csize, font="Courier New",
       color=INK, fill=OLDBG, line=OLD, space=0)
    tb(s, 6.78, 1.62, 6.10, 3.28, rcode, size=csize, font="Courier New",
       color=INK, fill=NEWBG, line=NEW, space=0)
    tb(s, 0.45, 5.06, 12.43, 2.04, why, size=13, color=INK, fill=PANEL, space=5)
    return s

# ---------------------------------------------------------------- 1. title
s = prs.slides.add_slide(BLANK)
bg = s.shapes.add_shape(1, 0, 0, prs.slide_width, prs.slide_height)
bg.fill.solid(); bg.fill.fore_color.rgb = NAVY; bg.line.fill.background()
tb(s, 0.9, 2.15, 11.5, 1.60, "From the spin-dependent NNN model\nto the checkerboard model",
   size=38, bold=True, color=WHITE, space=4)
tb(s, 0.9, 3.85, 11.5, 0.9,
   "Record of every code change, with the reason for each",
   size=19, color=RGBColor(0xCA, 0xDC, 0xFC))
tb(s, 0.9, 5.35, 11.5, 1.0,
   "Original   /home/phd25imran/CPQMC/Imran/Combined_Model\n"
   "Modified   /home/phd25imran/Checkerboard_Model",
   size=12, color=RGBColor(0x9F, 0xB0, 0xD8), font="Courier New", space=3)

# ---------------------------------------------------------------- 2. summary
s = prs.slides.add_slide(BLANK)
title(s, "Three files changed out of ten")
tb(s, 0.45, 1.15, 3.9, 1.55, "3", size=72, bold=True, color=NEW,
   fill=NEWBG, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
tb(s, 0.45, 2.78, 3.9, 0.42, "files modified", size=13, color=MUTE, align=PP_ALIGN.CENTER)
tb(s, 4.72, 1.15, 3.9, 1.55, "7", size=72, bold=True, color=MUTE,
   fill=PANEL, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
tb(s, 4.72, 2.78, 3.9, 0.42, "files byte-identical", size=13, color=MUTE, align=PP_ALIGN.CENTER)
tb(s, 8.99, 1.15, 3.89, 1.55, "0", size=72, bold=True, color=OLD,
   fill=OLDBG, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
tb(s, 8.99, 2.78, 3.89, 0.42, "changes to the QMC engine", size=13, color=MUTE, align=PP_ALIGN.CENTER)
tb(s, 0.45, 3.35, 6.1, 2.05,
   "MODIFIED\n"
   "  mc2duph.f90      +40 / -10   kinetic matrix tk\n"
   "  parameter.f90     +1         declares tp, tm\n"
   "  wf_unified.py   +103 / -30   trial wavefunction",
   size=12, font="Courier New", color=INK, fill=NEWBG, line=NEW, space=3)
tb(s, 6.78, 3.35, 6.10, 2.05,
   "UNCHANGED\n"
   "  libuph.f90    2699 lines     QMC engine\n"
   "  cp.f90  cpOut.f90  cpPara.f90  jiekou.f90\n"
   "  in.dat        Makefile",
   size=12, font="Courier New", color=INK, fill=PANEL, space=3)
tb(s, 0.45, 5.55, 12.43, 1.5,
   "The Monte Carlo engine was never touched. libuph.f90 is 2699 lines and is "
   "byte-identical between the two versions, as are the driver, the I/O, the "
   "parallel layer and the build files.\n"
   "Only the definition of the model changed. Any difference in the results "
   "therefore comes from the physics of the Hamiltonian and not from the sampler.",
   size=13, color=INK, fill=PANEL, space=6)

# ---------------------------------------------------------------- 3. the idea
s = prs.slides.add_slide(BLANK)
title(s, "What the change does, in one line")
tb(s, 0.45, 1.28, 6.10, 0.34, "BEFORE", size=12, bold=True, color=OLD)
tb(s, 6.78, 1.28, 6.10, 0.34, "AFTER", size=12, bold=True, color=NEW)
tb(s, 0.45, 1.62, 6.10, 1.70,
   "The two spin species hop differently.\n\n"
   "Spin up sees  +t1  on one diagonal,\n"
   "spin down sees  -t1  on the same bond.\n"
   "Identical at every site.",
   size=15, color=INK, fill=OLDBG, line=OLD, space=4)
tb(s, 6.78, 1.62, 6.10, 1.70,
   "The two spin species hop identically.\n\n"
   "Both spins see the same amplitude.\n"
   "What alternates is the GEOMETRY:\n"
   "A and B sublattices are 90 deg rotated.",
   size=15, color=INK, fill=NEWBG, line=NEW, space=4)
tb(s, 0.45, 3.50, 6.10, 1.35,
   "Spin splitting is BUILT IN.\nIt exists already at U = 0.",
   size=16, bold=True, color=OLD, fill=OLDBG, align=PP_ALIGN.CENTER,
   anchor=MSO_ANCHOR.MIDDLE, space=4)
tb(s, 6.78, 3.50, 6.10, 1.35,
   "Spin splitting must EMERGE.\nIt is exactly zero at U = 0.",
   size=16, bold=True, color=NEW, fill=NEWBG, align=PP_ALIGN.CENTER,
   anchor=MSO_ANCHOR.MIDDLE, space=4)
tb(s, 0.45, 5.06, 12.43, 1.9,
   "This is the whole point of the new model. In the old code the altermagnetic "
   "splitting is put in by hand through spin-dependent hopping, so the model can "
   "describe the consequences of altermagnetism but not its origin.\n"
   "In the new code nothing in the Hamiltonian distinguishes up from down. Any "
   "splitting that appears has to be generated by the on-site repulsion U, which "
   "is what makes the result an emergent altermagnet.",
   size=13, color=INK, fill=PANEL, space=6)

# ---------------------------------------------------------------- 4. mc2duph
codeslide(
  "mc2duph.f90  —  the kinetic matrix",
  "subroutine that fills tk(:,:,spin)     +40 / -10 lines     THE physics change",
  "BEFORE   Combined_Model",
"""do 121 i=1,nsites_cu
   tk(i, iposit(ixv(i)+1,iyv(i)+1), 1)= t1
   tk(i, iposit(ixv(i)-1,iyv(i)+1), 1)=-t1
   tk(i, iposit(ixv(i)-1,iyv(i)-1), 1)= t1
   tk(i, iposit(ixv(i)+1,iyv(i)-1), 1)=-t1

   tk(i, iposit(ixv(i)+1,iyv(i)+1), 2)=-t1
   tk(i, iposit(ixv(i)-1,iyv(i)+1), 2)= t1
   tk(i, iposit(ixv(i)-1,iyv(i)-1), 2)=-t1
   tk(i, iposit(ixv(i)+1,iyv(i)-1), 2)= t1
121 continue

   spin 1 and spin 2 get OPPOSITE signs
   no dependence on the site index""",
  "AFTER   Checkerboard_Model",
"""tp = t1 + t2
tm = t1 - t2
do 121 i=1,nsites_cu
 if ( mod( ixv(i)+iyv(i), 2 ) == 0 ) then
   ! A site:  '/' = tp ,  '\\' = tm
   tk(i, iposit(ixv(i)+1,iyv(i)+1), 1)=tp
   tk(i, iposit(ixv(i)-1,iyv(i)-1), 1)=tp
   tk(i, iposit(ixv(i)+1,iyv(i)-1), 1)=tm
   tk(i, iposit(ixv(i)-1,iyv(i)+1), 1)=tm
   tk(i, iposit(ixv(i)+1,iyv(i)+1), 2)=tp
   tk(i, iposit(ixv(i)-1,iyv(i)-1), 2)=tp
   tk(i, iposit(ixv(i)+1,iyv(i)-1), 2)=tm
   tk(i, iposit(ixv(i)-1,iyv(i)+1), 2)=tm
 else
   ! B site: same 8 lines, tp and tm swapped
 end if
121 continue""",
  "The spin index is the whole story. On the left, tk(...,1) and tk(...,2) carry "
  "opposite signs, so the two spin species see different Hamiltonians before the "
  "interaction is switched on.\n"
  "On the right they carry the same value on every bond. What alternates instead "
  "is the site parity mod(ixv+iyv,2), which makes the A and B sublattices differ "
  "by a 90 degree rotation. That is the checkerboard geometry, and it is the "
  "condition an altermagnet needs. The strengths are tp = t1+t2 and tm = t1-t2, "
  "so the anisotropy is set by t2 = -delta.")

# ---------------------------------------------------------------- 5. verification
s = prs.slides.add_slide(BLANK)
title(s, "Verification: the splitting really is gone at U = 0",
      "both kinetic matrices rebuilt exactly as mc2duph.f90 fills them, then diagonalised")
rows = [("", "max |K_up - K_dn|", "splitting at fixed k", ""),
        ("BEFORE", "0.6000", "2.4000", "spin split"),
        ("AFTER",  "0.0000", "0.0000", "degenerate")]
xs, ws = [0.45, 3.05, 6.55, 10.35], [2.5, 3.4, 3.7, 2.5]
for r, row in enumerate(rows):
    for c, val in enumerate(row):
        head = (r == 0)
        col = MUTE if head else (OLD if r == 1 else NEW)
        f = None if head else (OLDBG if r == 1 else NEWBG)
        tb(s, xs[c], 1.35 + r*0.72, ws[c], 0.62, val,
           size=11 if head else 15, bold=not head, color=col, fill=f,
           font="Calibri" if (head or c == 0 or c == 3) else "Courier New",
           align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
tb(s, 0.45, 3.72, 12.43, 1.25,
   "In the modified code K_up and K_dn are bit-for-bit identical, so no spin "
   "splitting is possible at any k before U is switched on.\n"
   "In the original code they differ by 0.6 and the splitting reaches 8|t1| = 2.4, "
   "with the d_xy form  E_up(k) - E_dn(k) = -8 t1 sin(kx) sin(ky), reproduced to 9e-16.",
   size=13, color=INK, fill=PANEL, space=6)
tb(s, 0.45, 5.15, 12.43, 1.85,
   "CAREFUL: comparing the sorted eigenvalue lists gives ZERO for both models.\n"
   "In the original model the two spin Hamiltonians are related by a 90 degree "
   "rotation, so their eigenvalue sets are identical even though the bands are "
   "split at each k. That degeneracy of the sets IS the altermagnetic "
   "compensation: equal density of states for both spins and zero net moment.\n"
   "This is why the paper measures the k-resolved polarization and not the total "
   "magnetization, which vanishes in both models and cannot tell them apart.",
   size=12.5, color=INK, fill=OLDBG, line=OLD, space=5)

# ---------------------------------------------------------------- 6. parameter.f90
codeslide(
  "parameter.f90  —  one declaration",
  "module variables     +1 line",
  "BEFORE   Combined_Model",
"""!------------HoKinAnHop-------------!

real(sp)::alpha,alphat1,ttp,ttn,tam




!for specify the k-point 1-kx 2-ky""",
  "AFTER   Checkerboard_Model",
"""!------------HoKinAnHop-------------!

real(sp)::alpha,alphat1,ttp,ttn,tam
real(sp)::tp,tm   ! checkerboard
                  ! anisotropic diagonals:
                  ! tp = t1+t2, tm = t1-t2

!for specify the k-point 1-kx 2-ky""",
  "The smallest change in the project, and it exists only because Fortran requires "
  "every variable to be declared before use.\n"
  "tp and tm hold the two diagonal strengths used by mc2duph.f90. Nothing else in "
  "the module was touched, so no existing variable changed meaning and no routine "
  "that used the old declarations needed recompiling differently.", csize=10)

# ---------------------------------------------------------------- 7. wf GetK
codeslide(
  "wf_unified.py  —  trial wavefunction, kinetic matrix",
  "def GetK(...)     the Hartree-Fock trial state handed to the QMC",
  "BEFORE   Combined_Model",
"""def GetK(Lx, Ly, t, tA, tt):

    # NN with spin-dependent anisotropy
    K_mat[i][down]  -= t - tA
    K_mat[i][right] -= t + tA

    # NNN, uniform over the lattice
    K_mat[i][up_right]   += tt
    K_mat[i][down_right] -= tt

    return K_mat""",
  "AFTER   Checkerboard_Model",
"""def GetK(Lx, Ly, t, tprime, delta):

    tp = tprime + delta
    tm = tprime - delta

    # NN uniform, spin independent
    K_mat[i][down]  -= t
    K_mat[i][right] -= t

    # diagonals alternate by sublattice
    if (ix + iy) % 2 == 0:
        ep, em = tp, tm     # A site
    else:
        ep, em = tm, tp     # B site

    return K_mat""",
  "The trial wavefunction has to describe the same model as the Fortran, otherwise "
  "the constrained path would be defined by a Hamiltonian the simulation is not "
  "solving.\n"
  "The signature changed from (t, tA, tt) to (t, tprime, delta), the spin-dependent "
  "nearest-neighbour term t +/- tA became a uniform t, and the uniform diagonals "
  "became the same parity branch used in mc2duph.f90. The parity rule was matched "
  "element by element against the Fortran tk matrix.")

# ---------------------------------------------------------------- 8. neel_seed
codeslide(
  "wf_unified.py  —  new starting point for the iteration",
  "def neel_seed(...)     new function, plus a seed argument on Iteration()",
  "BEFORE   Combined_Model",
"""def Iteration(Lx, Ly, U, tA, tt,
              NUP=None, NDN=None):

    # start from a random density
    # the spin-dependent hopping already
    # breaks the symmetry, so the
    # iteration has a direction to fall into""",
  "AFTER   Checkerboard_Model",
"""def neel_seed(Lx, Ly, NUP, NDN):
    for i in range(N):
        ix, iy = cart_coord(i + 1, Lx)
        if (ix + iy) % 2 == 0:
            nup[i], ndn[i] = 0.9, 0.1
        else:
            nup[i], ndn[i] = 0.1, 0.9
    return nup, ndn

def Iteration(Lx, Ly, U, tA, tt,
              NUP=None, NDN=None,
              seed='random'):""",
  "This addition follows directly from removing the spin dependence. In the old "
  "model the hopping itself broke the spin symmetry, so a random start would fall "
  "into the polarised solution on its own.\n"
  "In the new model the Hamiltonian is symmetric in spin, and a symmetric starting "
  "density stays symmetric under Hartree-Fock iteration. neel_seed puts 0.9 up on "
  "the A sublattice and 0.9 down on B, giving the iteration a broken configuration "
  "to converge from. The parity convention is the same one used in GetK.")

# ---------------------------------------------------------------- 9. closing
s = prs.slides.add_slide(BLANK)
title(s, "Summary of the change")
items = [
 ("mc2duph.f90", "Diagonal hopping made spin independent and alternated by "
                 "sublattice parity. This is the physics change.", NEW, NEWBG),
 ("parameter.f90", "Declares tp and tm. Bookkeeping only.", NEW, NEWBG),
 ("wf_unified.py", "Trial wavefunction rebuilt for the same model, and a "
                   "staggered seed added so the Hartree-Fock iteration can break "
                   "the spin symmetry.", NEW, NEWBG),
 ("everything else", "Untouched, including the 2699-line QMC engine libuph.f90.",
                     MUTE, PANEL),
]
y = 1.25
for name, txt, col, f in items:
    tb(s, 0.45, y, 2.75, 0.94, name, size=13, bold=True, color=col, fill=f,
       font="Courier New", anchor=MSO_ANCHOR.MIDDLE)
    tb(s, 3.40, y, 9.48, 0.94, txt, size=13, color=INK, anchor=MSO_ANCHOR.MIDDLE)
    y += 1.06
tb(s, 0.45, 5.60, 12.43, 1.45,
   "The result is that the old model imposes the altermagnetic splitting and the "
   "new one lets it emerge from U alone. Everything downstream of the Hamiltonian, "
   "including the sampling, the measurement and the analysis, is unchanged, so the "
   "two models can be compared directly.",
   size=13.5, color=INK, fill=PANEL, space=6)


# ================= PART 2: the Python susceptibility analysis =================

# ------------------------------------------------- 10. two kinds of change
s = prs.slides.add_slide(BLANK)
title(s, "The Python side changed in the opposite way")
tb(s, 0.45, 1.22, 6.10, 0.32, "FORTRAN   modified in place", size=12, bold=True, color=OLD)
tb(s, 6.78, 1.22, 6.10, 0.32, "PYTHON   extended, never modified", size=12, bold=True, color=NEW)
tb(s, 0.45, 1.58, 6.10, 1.74,
   "3 files modified\n0 files added\n+144 / -40 lines\n\n"
   "mc2duph.f90 builds tk internally,\nso the lattice is hard-coded and the\n"
   "source had to be edited.",
   size=13, color=INK, fill=OLDBG, line=OLD, space=3)
tb(s, 6.78, 1.58, 6.10, 1.74,
   "0 files modified\n13 files added\n+939 / -0 lines\n\n"
   "CPMC(..., K=..., K_dn=None) takes the\nhopping as an argument, so a new\n"
   "lattice is only a new matrix.",
   size=13, color=INK, fill=NEWBG, line=NEW, space=3)
tb(s, 0.45, 3.50, 12.43, 1.52,
   "3304b77   checkerboard model + d_xy pairing channel        569 lines, 5 files\n"
   "f7656fc   chi_zz(q) driver and three plot scripts          146 lines, 4 files\n"
   "57c7c8b   U-scan grid, timing, two more plot scripts       224 lines, 4 files",
   size=12.5, font="Courier New", color=INK, fill=PANEL, space=5)
tb(s, 0.45, 5.16, 12.43, 1.92,
   "The susceptibility engine cpqmc.py was never touched by any of this work. Its "
   "last change is commit f37d541, which predates the checkerboard entirely.\n"
   "So the two codebases give the same guarantee by different routes. In Fortran "
   "the sampler was left alone and only the model edited. In Python nothing "
   "existing was edited at all, and the checkerboard was added alongside as new "
   "files. Either way, results from the two lattices are produced by identical "
   "machinery.",
   size=13, color=INK, fill=PANEL, space=6)

# ------------------------------------------------- 11. checkerboard.py hopping
codeslide(
  "checkerboard.py  —  the same model, in Python",
  "pyqmc/checkerboard.py     new file, 166 lines     must agree with the Fortran exactly",
  "FORTRAN   mc2duph.f90 lines 155-178",
"""tp = t1 + t2
tm = t1 - t2
do 121 i=1,nsites_cu
 if ( mod( ixv(i)+iyv(i), 2 ) == 0 ) then
   ! A site
   tk(...,1)=tp     ! '/'
   tk(...,1)=tm     ! '\\'
 else
   ! B site: tp and tm swapped
 end if
121 continue""",
  "PYTHON   checkerboard_hopping()",
"""def checkerboard_hopping(lx, ly, t0=-1.0,
                         t1=0.3, t2=0.0):
    tp = t1 + t2
    tm = t1 - t2
    for x in range(lx):
      for y in range(ly):
        i = idx(x, y)
        K[i, idx(x+1, y)] = t0     # NN
        ...
        a, b = (tp, tm) if (x+y) % 2 == 0 \\
               else (tm, tp)
        K[i, idx(x+1, y+1)] = a    # '/'
        K[i, idx(x-1, y+1)] = b    # '\\'
    return 0.5 * (K + K.T)""",
  "The Python model has to be the same Hamiltonian as the Fortran, otherwise the "
  "two sets of results could not be compared and the trial wavefunction would not "
  "match the simulation.\n"
  "The parity test is identical: mod(ixv+iyv,2) in Fortran, (x+y) % 2 in Python, "
  "with the same tp and tm. This is not left to inspection. Gate 1 of the "
  "validation checks the two matrices element by element and requires max|diff| = 0.")

# ------------------------------------------------- 12. the dxy channel
codeslide(
  "diag_bond_factors()  —  the new pairing channel",
  "pyqmc/checkerboard.py     the d_xy form factor, which did not exist before",
  "BEFORE   square-lattice code",
"""# available pairing channels:
#
#   on-site s      (i == j)
#   extended s     NN bonds, all +1
#   d_x2-y2        NN bonds, +1 on x
#                            -1 on y
#
# all three live on NEAREST
# NEIGHBOUR bonds only.
#
# there is no form factor on the
# diagonal bonds, so d_xy cannot
# be measured at all.""",
  "AFTER   checkerboard code",
"""def diag_bond_factors(lx, ly):
    \"\"\"dxy form factor on the DIAGONAL
    bonds: +1 on '/', -1 on '\\'
    -> f(k) ~ sin kx sin ky\"\"\"
    for x in range(lx):
      for y in range(ly):
        m = idx(x, y)
        for (j, f) in [
          (idx(x+1, y+1), +1.0),
          (idx(x-1, y-1), +1.0),
          (idx(x-1, y+1), -1.0),
          (idx(x+1, y-1), -1.0)]:
            Fdxy[m, j] += f
    return Fdxy""",
  "This is new physics rather than a port. The square-lattice code carried s, "
  "extended s and d_x2-y2, all defined on nearest-neighbour bonds. d_xy has the "
  "form sin(kx)sin(ky), which is supported on the DIAGONAL bonds, so it had no "
  "form factor in the old code and simply could not be measured.\n"
  "Since d_xy pairing is the main result of the paper, the new lattice alone would "
  "not have been enough. The observable had to be added as well.")

# ------------------------------------------------- 13. validation
s = prs.slides.add_slide(BLANK)
title(s, "The new code was gated before it was used",
      "pyqmc/validate_checkerboard*.py     342 lines across three scripts")
gates = [
 ("GATE 1", "hopping bit-for-bit against the Fortran tk and against the "
            "independent ED builder", "max|diff| = 0 required"),
 ("GATE 2", "U = 0 CPMC energy against the free-fermion result",
            "checks the engine, not the model"),
 ("GATE 3", "U = 0 connected vertex is zero in every channel: s, d_x2-y2, d_xy",
            "the machinery gate for the new channel"),
 ("GATE 4", "interacting energy against QuSpin exact diagonalization on a small "
            "cluster", "optional, run with --ed"),
]
y = 1.42
for g, what, note in gates:
    tb(s, 0.45, y, 1.55, 0.92, g, size=13, bold=True, color=NEW, fill=NEWBG,
       font="Courier New", align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
    tb(s, 2.20, y, 7.05, 0.92, what, size=13, color=INK, anchor=MSO_ANCHOR.MIDDLE)
    tb(s, 9.40, y, 3.48, 0.92, note, size=11.5, color=MUTE, anchor=MSO_ANCHOR.MIDDLE)
    y += 1.02
tb(s, 0.45, 5.62, 12.43, 1.42,
   "Gate 3 is the one that matters for the new channel. A connected vertex must "
   "vanish at U = 0 by construction, so a nonzero value there would mean the d_xy "
   "form factor or the contraction was wrong rather than that pairing had been "
   "found. Passing it means the channel measures what it claims to measure.\n"
   "A second validator repeats the gates without QuSpin, so the checks can be run "
   "on machines where that library is not installed.",
   size=13, color=INK, fill=PANEL, space=6)


# ------------------------------------------------- 14. summary + next week
s = prs.slides.add_slide(BLANK)
title(s, "Summary")
bullets = [
 ("Emergent altermagnetism in a spin-independent model.",
  "The hopping carries no spin index, so the splitting is exactly zero at U = 0 "
  "and everything above it comes from the interaction."),
 ("The splitting is d_xy at every parameter set examined.",
  "126 of 126 interacting cells across L = 8, 10 and 12, classified by symmetry "
  "alone rather than by projection onto one harmonic."),
 ("Anisotropy ENHANCES d_xy pairing and suppresses d_x2-y2.",
  "The d_xy vertex leads for delta >= 0.3 at every U and grows by 2.3 to 3.2 times "
  "between L = 8 and L = 12. The earlier opposite conclusion came from L = 6."),
 ("Magnetism and pairing select the SAME representation.",
  "Both land in d_xy, and the channel is fixed by the lattice geometry rather than "
  "by the interaction."),
 ("Long-range order is not established at these sizes.",
  "S(pi,pi)/N falls as 1/N at every filling. The h -> 0 intercept is finite at each "
  "L but decreases with L, and L = 16 is running to settle the extrapolation."),
]
y = 1.12
for head, sub in bullets:
    tb(s, 0.45, y, 12.43, 0.36, head, size=14, bold=True, color=NAVY)
    tb(s, 0.75, y + 0.37, 12.13, 0.56, sub, size=12, color=MUTE, space=0)
    y += 0.96
tb(s, 0.45, 5.95, 6.10, 1.22,
   "NEXT WEEK\n\nFinish writing the manuscript.",
   size=15, bold=True, color=NEW, fill=NEWBG, line=NEW, space=6)
tb(s, 6.78, 5.95, 6.10, 1.22,
   "Title and abstract are settled and the introduction is drafted.\n"
   "Remaining: results and discussion sections, and Fig. 6 once L = 16 lands.",
   size=12.5, color=INK, fill=PANEL, space=5)

out = "code_changes_checkerboard.pptx"
prs.save(out)
print(f"wrote {out}  ({len(prs.slides.__iter__.__self__._sldIdLst)} slides)")
