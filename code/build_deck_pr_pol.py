#!/usr/bin/env python3
"""Deck: real-space P(R) across all fillings, plus the delta-polarisation plots."""
import os
from PIL import Image
from pptx import Presentation
from pptx.util import Inches as I, Pt
from pptx.dml.color import RGBColor as C
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

HERE = os.path.dirname(os.path.abspath(__file__)); IMG = os.path.join(HERE, "img")
NAVY = C(0x1E, 0x27, 0x61); ICE = C(0xCA, 0xDC, 0xFC); WHITE = C(0xFF, 0xFF, 0xFF)
INK = C(0x13, 0x18, 0x2F); MUTE = C(0x5A, 0x64, 0x84); RED = C(0xB2, 0x18, 0x2B)
BLUE = C(0x21, 0x66, 0xAC); GREEN = C(0x1B, 0x78, 0x37)
CARD = C(0xF4, 0xF6, 0xFB); CARDL = C(0xE2, 0xE7, 0xF2)
WARM = C(0xFD, 0xF3, 0xF2); WARML = C(0xF0, 0xD5, 0xD2)
COOL = C(0xF1, 0xF5, 0xFA); COOLL = C(0xD8, 0xE2, 0xEE)
HF, BF = "Cambria", "Calibri"
W, H, M = 13.333, 7.5, 0.55

prs = Presentation(); prs.slide_width, prs.slide_height = I(W), I(H)
BLANK = prs.slide_layouts[6]


def slide(bg=WHITE):
    s = prs.slides.add_slide(BLANK)
    f = s.background.fill; f.solid(); f.fore_color.rgb = bg
    return s


def tb(s, x, y, w, h, text, size=13, font=BF, color=INK, bold=False,
       italic=False, align=PP_ALIGN.LEFT):
    box = s.shapes.add_textbox(I(x), I(y), I(w), I(h)); tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        r = p.add_run(); r.text = line
        fo = r.font; fo.name = font; fo.size = Pt(size)
        fo.bold = bold; fo.italic = italic; fo.color.rgb = color
    return box


def bullets(s, x, y, w, h, items, size=12.5, color=INK, space=7):
    box = s.shapes.add_textbox(I(x), I(y), I(w), I(h)); tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, it in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(space)
        r = p.add_run(); r.text = "•  " + it
        fo = r.font; fo.name = BF; fo.size = Pt(size); fo.color.rgb = color


def card(s, x, y, w, h, fill=CARD, line=CARDL):
    sh = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, I(x), I(y), I(w), I(h))
    sh.fill.solid(); sh.fill.fore_color.rgb = fill
    sh.line.color.rgb = line; sh.line.width = Pt(1); sh.shadow.inherit = False
    sh.text_frame.text = ""
    return sh


def head(s, title, kicker=None):
    tb(s, M, 0.28, W - 2 * M, 0.7, title, size=28, font=HF, bold=True)
    if kicker:
        tb(s, M, 0.98, W - 2 * M, 0.34, kicker, size=12.5, color=MUTE, italic=True)


def pic_fit(s, name, y, maxw=W - 2 * M, maxh=4.3):
    """Size by BOTH dimensions. add_picture with width only auto-scales height,
    which previously pushed two-row figures through the text below them."""
    p = os.path.join(IMG, name)
    iw, ih = Image.open(p).size; ar = iw / ih
    w = min(maxw, maxh * ar); h = w / ar
    return s.shapes.add_picture(p, I((W - w) / 2), I(y), width=I(w), height=I(h))


def table(s, x, y, w, rows, colw, rowh=0.34, fs=11):
    shp = s.shapes.add_table(len(rows), len(rows[0]), I(x), I(y), I(w), I(rowh * len(rows)))
    t = shp.table
    for j, cw in enumerate(colw): t.columns[j].width = I(cw)
    for i in range(len(rows)):
        t.rows[i].height = I(rowh)
        for j in range(len(rows[0])):
            cell = t.cell(i, j); v = rows[i][j]
            txt, col, bold = (v if isinstance(v, tuple) else (v, INK, False))
            cell.text = ""; cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            cell.margin_left = cell.margin_right = I(0.05)
            p = cell.text_frame.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
            r = p.add_run(); r.text = str(txt)
            f = r.font; f.name = BF; f.size = Pt(fs)
            if i == 0:
                f.bold = True; f.color.rgb = WHITE
                cell.fill.solid(); cell.fill.fore_color.rgb = NAVY
            else:
                f.bold = bold; f.color.rgb = col
                cell.fill.solid()
                cell.fill.fore_color.rgb = WHITE if i % 2 else C(0xF7, 0xF9, 0xFD)


def notes(s, t): s.notes_slide.notes_text_frame.text = t


# ---------------- 1 title ----------------
s = slide(NAVY)
tb(s, M, 2.15, W - 2 * M, 0.9, "Real-space pairing and the delta polarisation",
   size=38, font=HF, bold=True, color=WHITE)
tb(s, M, 3.15, W - 2 * M, 0.6, "P(R) decomposition at every filling, and pairing "
   "against the magnetic polarisation", size=19, font=HF, color=ICE)
tb(s, M, 4.15, W - 2 * M, 0.4, "Checkerboard Hubbard  ·  CPQMC at T = 0  ·  "
   "L = 12, six fillings, five U, eight anisotropies", size=14, color=ICE)
tb(s, M, 4.62, W - 2 * M, 0.35, "14 August 2026", size=13, color=C(0x8F, 0xA3, 0xD9))
tb(s, M, 5.6, 11.6, 0.6,
   "Headline: the d_xy channel wins at long range only at half filling. "
   "Doping kills it, and it turns repulsive.", size=15, color=WHITE, italic=True)
notes(s, "Full P(R) cube is now complete: 40 chunks, 155520 rows, U in 0,2,4,6,8 "
         "times eight anisotropies times six fillings times six seeds.")

# ---------------- 2 method ----------------
s = slide(); head(s, "What P(R) measures, and why it was needed",
                  "the q = 0 susceptibility adds local and long-range pairing together, and they disagree")
bullets(s, M, 1.5, W - 2 * M, 2.5, [
    "The pairing susceptibility sums over ALL site pairs (i, j). Keeping the "
    "separation R instead of summing gives P(R), and the vertex v(R) = P(R) - bubble(R).",
    "Huang, Lin and Gubernatis (PRB 64, 205101) showed that extended s wins LOCALLY "
    "while d wins at LONG RANGE. A q = 0 sum mixes the two, which is why our earlier "
    "numbers put extended s ahead of d and looked like they contradicted White.",
    "The driver stores a SUM over the pairs at each separation, and the multiplicity "
    "runs from 144 to 1728, so everything here is divided by the pair count. "
    "Without that division the distant shells look artificially large.",
    "Long range means R > 2, following Huang."], size=12.5)
card(s, M, 4.5, 6.05, 1.5, COOL, COOLL)
tb(s, M + 0.24, 4.66, 5.6, 0.3, "Validation: the sum rule", size=13, bold=True, color=BLUE)
tb(s, M + 0.24, 5.0, 5.6, 0.9,
   "Summing v(R) over all 27 separations must return the original vertex from the "
   "q = 0 driver. It does, to 0.00% on all four channels.", size=12)
card(s, M + 6.3, 4.5, 6.05, 1.5, WARM, WARML)
tb(s, M + 6.54, 4.66, 5.6, 0.3, "Scope", size=13, bold=True, color=RED)
tb(s, M + 6.54, 5.0, 5.6, 0.9,
   "L = 12, six fillings, U = 0, 2, 4, 6, 8, eight anisotropies, six seeds. "
   "40 of 40 chunks complete, 155,520 rows.", size=12)
notes(s, "The sum rule is the thing to quote if anyone asks whether the R-resolved "
         "code is doing the same calculation as the original driver. It is exact.")

# ---------------- 3 profiles ----------------
s = slide(); head(s, "The vertex resolved by separation",
                  "U = 4, half filling: R = 1 is the bond, R = sqrt(2) the plaquette diagonal that the anisotropy acts on")
pic_fit(s, "fig_pr_profiles.png", 1.42, maxh=3.5)
card(s, M, 5.15, 12.23, 1.75, CARD, CARDL)
tb(s, M + 0.25, 5.32, 11.8, 0.3, "Reading the panels", size=13, bold=True, color=NAVY)
tb(s, M + 0.25, 5.68, 11.8, 1.1,
   "At delta = 0 extended s owns the bond and d_x2-y2 carries the tail, exactly the "
   "Huang result. Turning on delta lights up a large positive d_xy spike at R = sqrt(2), "
   "the diagonal bond the anisotropy strengthens, and d_x2-y2 goes negative there. "
   "By delta = 0.4 d_xy leads across the whole shaded long-range region. At delta = 0.7 "
   "it is weakening again.", size=12)
notes(s, "R = 0 is left off this figure because it is roughly twenty times larger than "
         "everything else; it has its own panel in the next slide.")

# ---------------- 4 local vs long ----------------
s = slide(); head(s, "Local and long-range pairing separate cleanly",
                  "U = 4, half filling: the two panels answer different questions")
pic_fit(s, "fig_pr_local_long.png", 1.4, maxh=3.7)
card(s, M, 5.3, 6.05, 1.6, COOL, COOLL)
tb(s, M + 0.24, 5.46, 5.6, 0.3, "Local  v(R = 0)", size=13, bold=True, color=BLUE)
tb(s, M + 0.24, 5.8, 5.6, 1.0,
   "Extended s leads at every delta. d_xy starts repulsive, crosses zero near "
   "delta = 0.15 and overtakes d_x2-y2 by delta = 0.3. On-site s is strongly "
   "repulsive throughout, as it must be.", size=12)
card(s, M + 6.3, 5.3, 6.05, 1.6, WARM, WARML)
tb(s, M + 6.54, 5.46, 5.6, 0.3, "Long range  R > 2", size=13, bold=True, color=RED)
tb(s, M + 6.54, 5.8, 5.6, 1.0,
   "d_x2-y2 leads at delta = 0 and 0.1, then collapses. d_xy takes over from "
   "delta = 0.2, peaks near delta = 0.3 and decays. Extended s is never the "
   "long-range winner.", size=12)
notes(s, "This is the slide that answers the supervisor's question about where each "
         "channel actually dominates: extended s locally, d or d_xy at long range.")

# ---------------- 5 all U ----------------
s = slide(); head(s, "The long-range takeover happens at every interaction strength",
                  "half filling; on-site s omitted because it is strongly repulsive and not a competitor")
pic_fit(s, "fig_pr_longrange_U.png", 1.42, maxh=3.5)
t = [["U", "long-range winner at delta = 0", "winner by delta = 0.4", "d_xy crossing"],
     ["2", "d_x2-y2", ("d_xy", RED, True), "0.3"],
     ["4", "d_x2-y2", ("d_xy", RED, True), "0.2"],
     ["6", "extended s", ("d_xy", RED, True), "0.2"],
     ["8", "extended s", ("d_xy", RED, True), "0.2"]]
table(s, M, 5.15, 7.1, t, [0.9, 2.6, 2.1, 1.5])
card(s, 8.1, 5.15, 4.68, 1.75, CARD, CARDL)
tb(s, 8.34, 5.32, 4.3, 0.3, "What stays the same", size=13, bold=True, color=NAVY)
tb(s, 8.34, 5.68, 4.3, 1.1,
   "d_xy ends up the long-range leader at every U. The signal shrinks steadily "
   "with U, so by U = 8 all three channels are pressed close to zero and the "
   "margin is small.", size=12)
notes(s, "Note the vertical scale is shared. The U = 8 panel is not empty, it is "
         "genuinely small.")

# ---------------- 6 fillings, the key new result ----------------
s = slide(); head(s, "Away from half filling the d_xy channel dies",
                  "long-range vertex against anisotropy at all six fillings, U = 4")
pic_fit(s, "fig_pr_fillings.png", 1.35, maxh=3.9)
card(s, M, 5.45, 12.23, 1.5, WARM, WARML)
tb(s, M + 0.25, 5.6, 11.8, 0.3, "This is the most important slide in the deck",
   size=13, bold=True, color=RED)
tb(s, M + 0.25, 5.96, 11.8, 1.0,
   "Only at n = 1.000 does d_xy rise above zero and take the long-range lead. At every "
   "doped filling it is negative and falls further as delta grows, reaching -0.9 at "
   "n = 0.500. d_x2-y2 is the long-range winner everywhere except half filling. "
   "The channel swap is a half-filling phenomenon, not a generic one.", size=12)
notes(s, "If asked why: the d_xy form factor lives on the plaquette diagonals, and the "
         "diagonal correlations that support it are tied to the commensurate structure "
         "at half filling. Doping destroys that and the channel goes repulsive.")

# ---------------- 7 vs filling ----------------
s = slide(); head(s, "The same result read along the filling axis",
                  "delta = 0.4, U = 4: local on the left, long range on the right")
pic_fit(s, "fig_pr_vs_filling.png", 1.4, maxh=3.7)
card(s, M, 5.3, 12.23, 1.6, CARD, CARDL)
tb(s, M + 0.25, 5.46, 11.8, 0.3, "Reading it", size=13, bold=True, color=NAVY)
tb(s, M + 0.25, 5.82, 11.8, 1.0,
   "Locally the picture is smooth: extended s and d_xy both strengthen towards half "
   "filling. At long range the three channels are all negative or near zero across the "
   "doped region and only separate in the last step to n = 1. The honest statement is "
   "that our long-range pairing signal exists near half filling and nowhere else.", size=12)

# ---------------- 8 polarisation definition ----------------
s = slide(); head(s, "The delta polarisation, and what plotting against it tests",
                  "replacing the m_AM = delta m times delta construction")
bullets(s, M, 1.5, W - 2 * M, 2.4, [
    "S_zz(pi, pi) is the spin structure factor at the Neel wavevector. "
    "m = sqrt( S_zz(pi,pi) / N ) is the staggered polarisation; the 1/N inside the "
    "root is what makes m comparable across system sizes.",
    "Delta m = m - m(U = 0) removes the free-fermion background. That background is "
    "not small: at L = 12 it is 60 to 93 per cent of m, so delta m is a difference "
    "between two comparable numbers and should be quoted with that caveat.",
    "Putting delta m on the x-axis collapses the two knobs U and delta into one "
    "question: is the pairing controlled by the magnetic polarisation? If every "
    "(U, delta) point fell on one curve the answer would be yes."], size=12.5)
card(s, M, 4.35, 6.05, 1.5, COOL, COOLL)
tb(s, M + 0.24, 4.5, 5.6, 0.3, "Why drop the times delta factor", size=13, bold=True, color=BLUE)
tb(s, M + 0.24, 4.85, 5.6, 1.0,
   "delta already enters through the suppression of m, so multiplying by it again "
   "double counts. arXiv 2607.06106 defines the order parameter on this lattice as "
   "the plain staggered magnetisation.", size=12)
card(s, M + 6.3, 4.35, 6.05, 1.5, WARM, WARML)
tb(s, M + 6.54, 4.5, 5.6, 0.3, "The answer, in one line", size=13, bold=True, color=RED)
tb(s, M + 6.54, 4.85, 5.6, 1.0,
   "Within each U the pairing tracks delta m tightly. Across U the curves are "
   "offset and do not collapse, so delta m alone does not control the pairing.", size=12)
notes(s, "m(U=0) is finite because free fermions already have spin correlations from "
         "Pauli statistics. It is a background, not order.")

# ---------------- 9 polarisation all fillings ----------------
s = slide(); head(s, "Pairing against the delta polarisation, every filling",
                  "each point is one (U, delta); lines join constant U; r is the within-U correlation")
pic_fit(s, "fig_pol_fillings.png", 1.35, maxh=3.9)
t2 = [["n", "extended s", "d_x2-y2", "d_xy"],
      ["1.000", ("+0.98", GREEN, True), ("+0.93", BLUE, True), ("-0.79", RED, True)],
      ["0.889", ("+0.92", GREEN, True), ("+0.93", BLUE, True), ("-0.83", RED, True)],
      ["0.778", "+0.47", ("+0.83", BLUE, True), "-0.25"],
      ["0.667", "+0.34", "+0.77", "+0.75"],
      ["0.556", "-0.17", "+0.06", "+0.07"],
      ["0.500", "-0.53", "-0.62", "-0.46"]]
table(s, M, 5.42, 5.4, t2, [1.0, 1.6, 1.4, 1.4], rowh=0.21, fs=10)
card(s, 6.35, 5.42, 6.43, 1.5, WARM, WARML)
tb(s, 6.59, 5.56, 6.0, 0.3, "The anti-correlation is a near-half-filling effect",
   size=12.5, bold=True, color=RED)
tb(s, 6.59, 5.9, 6.0, 1.0,
   "At n = 1.000 and 0.889 extended s and d_x2-y2 track the polarisation at r above "
   "+0.9 while d_xy runs against it at about -0.8. By n = 0.667 d_xy has flipped to "
   "+0.75, and below that the correlations are weak and unstable. Same boundary as "
   "the P(R) result.", size=11.5)

# ---------------- 10 pairing phase diagram ----------------
s = slide(); head(s, "Single-panel pairing phase diagram in the (n, delta) plane",
                  "built to mirror Fig. 4 of PRB 113, 134443, with pairing channels in place of ordering vectors")
pic_fit(s, "fig_pairing_phase.png", 1.3, maxh=4.05)
card(s, M, 5.5, 12.23, 1.45, WARM, WARML)
tb(s, M + 0.25, 5.64, 11.8, 0.3, "How to read it", size=13, bold=True, color=RED)
tb(s, M + 0.25, 5.98, 11.8, 1.0,
   "Colour is the difference between the two long-range vertices, so blue means "
   "d_x2-y2 leads and red means d_xy leads. The scale is symmetric log because the "
   "repulsive corner at low filling is ten times larger than the d_xy signal. Hatching "
   "marks where the two channels sit within one combined error bar. Extended s never "
   "leads anywhere on this grid, which is why the panel needs only two colours.", size=12)
notes(s, "The analogue of their red line is the d_x2-y2 to d_xy boundary. Their Q "
         "coordinate said which magnetic ordering vector wins; this difference says "
         "which pairing channel wins, and it is a single continuous number for the "
         "same reason: only two channels are ever in contention.")

# ---------------- 11 where extended s lives ----------------
s = slide(); head(s, "Where extended s dominates, and why the previous panel had only two colours",
                  "same plane, dominant channel by category, both pieces of the vertex")
pic_fit(s, "fig_pairing_phase_2panel.png", 1.3, maxh=3.9)
card(s, M, 5.35, 6.05, 1.55, COOL, COOLL)
tb(s, M + 0.24, 5.5, 5.6, 0.3, "Panel (a), local", size=13, bold=True, color=BLUE)
tb(s, M + 0.24, 5.86, 5.6, 1.0,
   "Extended s wins 17 of 48 cells, in a solid block at n = 0.889 and 1.000 for "
   "delta up to 0.5. d_xy takes 22 cells in the middle. On-site s wins nowhere, "
   "as expected for a repulsive channel.", size=12)
card(s, M + 6.3, 5.35, 6.05, 1.55, WARM, WARML)
tb(s, M + 6.54, 5.5, 5.6, 0.3, "Panel (b), long range", size=13, bold=True, color=RED)
tb(s, M + 6.54, 5.86, 5.6, 1.0,
   "Extended s wins 0 of 48. d_x2-y2 takes 37 cells and d_xy 11, all of them at "
   "high filling. That zero is why the continuous panel needed only two colours.", size=12)
notes(s, "This is the Huang result stated as a map rather than a pair of curves: "
         "extended s is a local channel and d is a long-range one. Hatched cells are "
         "where the top two are within one combined standard error.")

# ---------------- 12 summary ----------------
s = slide(); head(s, "Where this leaves us", "two independent measurements, one consistent boundary")
bullets(s, M, 1.45, W - 2 * M, 3.5, [
    "P(R) reproduces Huang at delta = 0: extended s wins locally, d_x2-y2 wins at "
    "long range. Our earlier q = 0 result was not in conflict with White, it was "
    "summing two competing pieces.",
    "Anisotropy hands long-range leadership to d_xy at every U, entering through a "
    "large positive spike on the plaquette diagonal at R = sqrt(2).",
    "That takeover exists only near half filling. At n <= 0.778 the long-range d_xy "
    "vertex is negative and falls further with delta.",
    "The delta-polarisation plots draw the same boundary independently: d_xy is "
    "anti-correlated with the polarisation at n = 1.000 and 0.889, and stops being "
    "so below that.",
    "Two different observables agreeing on the same filling boundary is the "
    "strongest evidence in the deck.",
    "The phase diagram puts the d_xy region at n above roughly 0.88 and delta above "
    "roughly 0.2, with a competition pocket at low delta and intermediate filling.",
    "Extended s dominates the LOCAL vertex over 17 of 48 cells but never leads at "
    "long range, which is the Huang local-versus-long-range split stated as a map.",
    "Still open: pinning the delta at which d_xy changes sign. Twist-averaged runs "
    "on a finer grid at U = 2 and 4 are in flight."], size=12.5)
card(s, M, 5.5, 12.23, 1.4, COOL, COOLL)
tb(s, M + 0.25, 5.66, 11.8, 0.3, "What we cannot claim", size=13, bold=True, color=BLUE)
tb(s, M + 0.25, 6.0, 11.8, 0.9,
   "v(R) decays into noise past R = 3.5 and chi/N is flat in system size, so there is "
   "no long-range order here. The defensible claim is a change in WHICH channel is "
   "attractive, not superconductivity.", size=12)

out = os.path.join(HERE, "PR_and_polarisation_2026-08-14.pptx")
prs.save(out); print("saved", out, len(prs.slides.__iter__.__self__._sldIdLst), "slides")
