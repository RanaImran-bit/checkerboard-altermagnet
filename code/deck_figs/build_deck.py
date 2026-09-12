#!/usr/bin/env python3
"""Weekly progress deck, built with python-pptx (no node in this environment)."""
import os
from pptx import Presentation
from pptx.util import Inches as I, Pt
from pptx.dml.color import RGBColor as C
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, "img")

NAVY = C(0x1E, 0x27, 0x61); ICE = C(0xCA, 0xDC, 0xFC); WHITE = C(0xFF, 0xFF, 0xFF)
INK = C(0x13, 0x18, 0x2F); MUTE = C(0x5A, 0x64, 0x84); RED = C(0xC0, 0x39, 0x2B)
TEAL = C(0x1C, 0x72, 0x93); GOLD = C(0x8A, 0x6D, 0x1F); GREEN = C(0x1E, 0x7A, 0x46)
CARD = C(0xF4, 0xF6, 0xFB); CARDL = C(0xE2, 0xE7, 0xF2)
GRN2 = C(0x2C, 0x7A, 0x2E); BLU2 = C(0x1F, 0x5C, 0x93)
WARM = C(0xFD, 0xF3, 0xF2); WARML = C(0xF0, 0xD5, 0xD2)
COOL = C(0xF1, 0xF5, 0xFA); COOLL = C(0xD8, 0xE2, 0xEE)
AMB = C(0xFF, 0xF8, 0xE6); AMBL = C(0xEF, 0xE0, 0xBC)
HF, BF = "Cambria", "Calibri"
W, H, M = 13.333, 7.5, 0.55

prs = Presentation()
prs.slide_width, prs.slide_height = I(W), I(H)
BLANK = prs.slide_layouts[6]


def slide(bg=WHITE):
    s = prs.slides.add_slide(BLANK)
    r = s.background.fill; r.solid(); r.fore_color.rgb = bg
    return s


def tb(s, x, y, w, h, text, size=13, font=BF, color=INK, bold=False,
       italic=False, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, space=0):
    box = s.shapes.add_textbox(I(x), I(y), I(w), I(h))
    tf = box.text_frame; tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = anchor
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(space)
        run = p.add_run(); run.text = line
        f = run.font
        f.name = font; f.size = Pt(size); f.bold = bold; f.italic = italic
        f.color.rgb = color
    return box


def bullets(s, x, y, w, h, items, size=12.5, color=INK, space=5):
    box = s.shapes.add_textbox(I(x), I(y), I(w), I(h))
    tf = box.text_frame; tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, it in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(space)
        run = p.add_run(); run.text = "•  " + it
        f = run.font; f.name = BF; f.size = Pt(size); f.color.rgb = color
    return box


def card(s, x, y, w, h, fill=CARD, line=CARDL):
    sh = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, I(x), I(y), I(w), I(h))
    sh.fill.solid(); sh.fill.fore_color.rgb = fill
    sh.line.color.rgb = line; sh.line.width = Pt(1)
    sh.shadow.inherit = False
    sh.text_frame.text = ""
    return sh


def badge(s, n, x, y, fill=NAVY, d=0.42):
    sh = s.shapes.add_shape(MSO_SHAPE.OVAL, I(x), I(y), I(d), I(d))
    sh.fill.solid(); sh.fill.fore_color.rgb = fill
    sh.line.fill.background(); sh.shadow.inherit = False
    tf = sh.text_frame; tf.word_wrap = False
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    run = p.add_run(); run.text = str(n)
    f = run.font; f.name = BF; f.size = Pt(13); f.bold = True; f.color.rgb = WHITE
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE


def head(s, title, kicker=None):
    tb(s, M, 0.30, W - 2 * M, 0.74, title, size=31, font=HF, bold=True)
    if kicker:
        tb(s, M, 1.03, W - 2 * M, 0.34, kicker, size=13, color=MUTE, italic=True)


def table(s, x, y, w, rows, colw, rowh=0.40, fs=11.5, hdr_fill=NAVY):
    nr, nc = len(rows), len(rows[0])
    shp = s.shapes.add_table(nr, nc, I(x), I(y), I(w), I(rowh * nr))
    t = shp.table
    for j, cw in enumerate(colw):
        t.columns[j].width = I(cw)
    for i in range(nr):
        t.rows[i].height = I(rowh)
        for j in range(nc):
            cell = t.cell(i, j)
            val = rows[i][j]
            txt, col, bold = (val if isinstance(val, tuple) else (val, INK, False))
            cell.text = ""
            cell.margin_left = cell.margin_right = I(0.06)
            cell.margin_top = cell.margin_bottom = I(0.02)
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            p = cell.text_frame.paragraphs[0]
            p.alignment = PP_ALIGN.CENTER
            r = p.add_run(); r.text = str(txt)
            f = r.font; f.name = BF; f.size = Pt(fs)
            if i == 0:
                f.bold = True; f.color.rgb = WHITE
                cell.fill.solid(); cell.fill.fore_color.rgb = hdr_fill
            else:
                f.bold = bold; f.color.rgb = col
                cell.fill.solid()
                cell.fill.fore_color.rgb = WHITE if i % 2 else C(0xF7, 0xF9, 0xFD)
    return shp


def notes(s, txt):
    s.notes_slide.notes_text_frame.text = txt


def pic_h(s, name, x, y, h):
    """Place by HEIGHT. add_picture with width only auto-scales the height, which
    silently overran the text on the two-row figures (aspect 2.24 -> 5.45in tall)."""
    return s.shapes.add_picture(os.path.join(IMG, name), I(x), I(y), height=I(h))


# ============ 1. TITLE ============
s = slide(NAVY)
tb(s, M, 2.05, W - 2 * M, 0.9, "Emergent altermagnetism and pairing",
   size=40, font=HF, bold=True, color=WHITE)
tb(s, M, 2.92, W - 2 * M, 0.75, "in the checkerboard Hubbard model",
   size=32, font=HF, color=ICE)
tb(s, M, 3.86, W - 2 * M, 0.4,
   "Constrained-path quantum Monte Carlo   ·   weekly progress", size=15, color=ICE)
tb(s, M, 4.32, W - 2 * M, 0.35, "12 August 2026", size=13, color=C(0x8F, 0xA3, 0xD9))
tb(s, M, 5.6, 10.5, 0.5,
   "Headline: anisotropy selects the d_xy pairing channel at half filling",
   size=15, color=WHITE, italic=True)
notes(s, "Two cubes finished this week. Lead with the pairing channel swap: size independent, "
         "holds at every U, 27 sigma. Magnetism is the supporting explanation, not the lead.")

# ============ 2. WEEK AT A GLANCE ============
s = slide(); head(s, "This week", "what ran, what it answered, and what is still in flight")
stats = [("9,840", "measurements across\nthe two complete cubes", NAVY),
         ("~300", "node-hours of CPQMC,\nsix nodes in parallel", RED),
         ("1e-13", "U = 0 null-test residual\nin all four channels", TEAL),
         ("0 / 48", "L = 12 fillings that\nare closed shell", GOLD)]
for i, (big, lab, col) in enumerate(stats):
    x = M + i * 3.1
    card(s, x, 1.5, 2.95, 1.5)
    tb(s, x + 0.12, 1.62, 2.71, 0.72, big, size=30, font=HF, bold=True,
       color=col, align=PP_ALIGN.CENTER)
    tb(s, x + 0.12, 2.36, 2.71, 0.58, lab, size=11, color=MUTE, align=PP_ALIGN.CENTER)

paras = [
    ("This week we ran roughly 300 node-hours of CPQMC across six cluster nodes in parallel.", INK, False),
    ("Both cubes are now closed. L = 8, 10, 12 at U = 0, 2, 4, 6, 8, fillings from n = 0.5 to 1, "
     "eight anisotropies δ = 0 to 0.7, six independent seeds per point. That is 9,840 measurements "
     "with no gaps and no NaN, and the U = 0 rows pass the null test to 1e-13. We also measured the "
     "full pair-field susceptibility at δ = 0, so the vertex can now be compared against the "
     "uncorrelated susceptibility the way White et al. do.", INK, False),
    ("Last week we derived that interactions drive Néel order near half filling, and that with "
     "δ ≠ 0 that order becomes an altermagnet of d_xy symmetry. We also found the equal-time vertex "
     "is CP-biased, so we moved to the integrated vertex. This week we tried to answer: does the "
     "pairing follow the magnetic order, or just the anisotropy?", INK, False),
    ("The answer is that it follows both, and the two are not separable. The d_xy channel turns "
     "attractive only where the anisotropy is large AND the filling is above n ≈ 0.85 — the same "
     "region where the moment is strong. Away from half filling the same anisotropy drives d_xy "
     "further into repulsion.", RED, True),
    ("Currently twist-averaged boundary conditions are running at L = 12, to test whether the "
     "open-shell trial state is distorting the result. Those land tomorrow.", MUTE, False),
]
y = 3.15
for txt, col, bold in paras:
    lines = max(1, -(-len(txt) // 130))          # ceil, ~130 chars per line at 12.5pt
    h = 0.24 * lines + 0.06
    tb(s, M, y, W - 2 * M, h, txt, size=12.5, color=col, bold=bold)
    y += h + 0.12
notes(s, "Narrative version. The key sentence is the fourth paragraph - pairing follows both the "
         "magnetic order and the anisotropy, and they cannot be separated. Everything after this "
         "slide is evidence for that claim.")

# ============ 3. METHOD ============
s = slide(); head(s, "Method and observables",
                  "CPQMC, back-propagated, six independent seeds per point")
card(s, M, 1.55, 5.6, 0.72, COOL, COOLL)
tb(s, M + 0.26, 1.72, 5.1, 0.42, "χ  =  P_full  −  P_bub",
   size=19, font=HF, bold=True, color=TEAL)
tb(s, 6.4, 1.66, 6.4, 0.55,
   "P_bub is the uncorrelated susceptibility — the same diagram with the interaction vertex "
   "removed. Written P̄ in the figures.", size=11.5, color=MUTE)

bl = [
    "Constrained-path QMC at T = 0, back-propagated, six independent seeds per point.",
    "The vertex is full minus bubble, integrated over imaginary time. Positive means the "
    "interaction is attractive in that channel; negative means it opposes pairing.",
    "The equal-time vertex is CP-biased against exact diagonalisation, so we report the "
    "integrated form throughout.",
    "Four channels measured in a single back-propagation pass at q = 0 — on-site s, extended s, "
    "d_x²−y² and d_xy — differing only in the bond pattern the form factor weights.",
    "χ is an unnormalised double site sum and therefore extensive, so every comparison across "
    "lattice sizes uses χ/N.",
    "Magnetic moment m = √( S(q*) / N ), with Δm = m − m(U=0) isolating the interaction part, "
    "since S(q) is finite even for free fermions.",
    "q* is the measured peak of S(q), not assumed to be (π,π) — see the correction later.",
]
bullets(s, M, 2.58, W - 2 * M, 2.8, bl, size=13, space=7)

card(s, M, 5.55, W - 2 * M, 1.35, WARM, WARML)
tb(s, M + 0.28, 5.68, 7.0, 0.3, "Validation that the machinery is sound",
   size=14, bold=True, color=RED)
tb(s, M + 0.28, 6.0, W - 2 * M - 0.56, 0.8,
   "At U = 0 the vertex must vanish identically and the full susceptibility must equal the bubble. "
   "Measured residual: 1.7e-13 to 7.6e-13 across all four channels and all three lattices. "
   "On-site s is repulsive in all 576 measured points, as a positive U requires.", size=12.5)
notes(s, "The U=0 null test is the strongest single argument that the full-minus-bubble "
         "decomposition is implemented correctly.")

# ============ 5. DATASET ============
s = slide(); head(s, "Both cubes are now complete", "no gaps, no NaN, six seeds everywhere")
t1 = [["", "U = 0", "U = 2", "U = 4", "U = 6", "U = 8"],
      [("L = 8", INK, True), "✓", "✓", "✓", "✓", "✓"],
      [("L = 10", INK, True), "✓", "✓", "✓", "✓", "✓"],
      [("L = 12", INK, True), "✓", "✓", "✓", "✓", "✓"]]
table(s, M, 1.72, 6.0, t1, [1.3, 0.94, 0.94, 0.94, 0.94, 0.94], rowh=0.42)
tb(s, M, 3.5, 6.0, 0.3, "✓ = all eight anisotropies δ = 0 … 0.7, six fillings, six seeds",
   size=11, color=MUTE, italic=True)
t2 = [["dataset", "rows", "status"],
      ["magnetic_master", "4,320", ("complete", GREEN, True)],
      ["pairing_master", "5,520", ("complete", GREEN, True)],
      ["pairfull (δ = 0)", "540", ("complete", GREEN, True)],
      ["S(q) grids", "15 files", ("collected", GREEN, True)],
      ["twist average", "in flight", ("tomorrow", GOLD, True)]]
table(s, 7.05, 1.72, 5.73, t2, [2.5, 1.5, 1.73], rowh=0.40)
card(s, M, 4.35, W - 2 * M, 1.6)
tb(s, M + 0.28, 4.5, 7.0, 0.3, "Two data problems found and fixed", size=14, bold=True, color=NAVY)
bullets(s, M + 0.28, 4.85, W - 2 * M - 0.56, 1.0, [
    "χ is an unnormalised double site sum, so it is extensive — every cross-size comparison must "
    "use χ/N. Raw ratios matched N ratios to under 1%.",
    "The L = 12 S(q) grids had never been collected; 33 files were scattered over six nodes in four "
    "directories, now merged into one file per (L, U)."], size=12)
notes(s, "9840 QMC rows total, every point six seeds. Roughly 300 node-hours this week across five nodes.")

# ============ 6. HEADLINE ============
s = slide(); head(s, "Headline: anisotropy selects the d_xy channel",
                  "per-site vertex, half filling, U = 4 — identical at all three lattice sizes")
s.shapes.add_picture(os.path.join(IMG, "collapse4.png"), I(M), I(1.5), width=I(12.2))
pts = [("on-site s", "flat", "−0.031 → −0.031", MUTE),
       ("extended s", "decays", "+0.036 → +0.020", MUTE),
       ("d_x²−y²", "decays", "+0.028 → +0.007", MUTE),
       ("d_xy", "CHANGES SIGN", "−0.010 → +0.026", RED)]
for i, (a, b, c, col) in enumerate(pts):
    x = M + i * 3.07
    card(s, x, 4.72, 2.98, 1.3, WARM if col == RED else CARD, WARML if col == RED else CARDL)
    tb(s, x + 0.16, 4.84, 2.66, 0.3, a, size=13, bold=True)
    tb(s, x + 0.16, 5.15, 2.66, 0.32, b, size=13.5, bold=True, color=col)
    tb(s, x + 0.16, 5.52, 2.66, 0.3, c, size=11, color=MUTE)
tb(s, M, 6.18, 12.2, 0.34,
   "δ = 0 → 0.7 at U = 4.  Three channels are dull; only d_xy crosses zero, at δ ≈ 0.15, peaking near δ = 0.4.",
   size=12.5, italic=True)
notes(s, "The paper's lead result. Three lattices lie on top of each other once chi is divided by N "
         "- 1% spread at U=4. dxy is the channel the checkerboard symmetry predicts.")

# ============ 7. REORDERING ============
s = slide(); head(s, "The reordering happens between δ = 0.1 and 0.4",
                  "pairing vertex vs U, one panel per anisotropy, half filling")
pic_h(s, "G1.png", 0.5, 1.5, 4.05)          # aspect 2.238 -> 9.07 wide
card(s, 9.75, 1.5, 3.05, 4.05)
tb(s, 9.97, 1.66, 2.65, 0.32, "Reading the panels", size=13, bold=True, color=NAVY)
bullets(s, 9.97, 2.05, 2.65, 3.35, [
    "δ = 0, 0.1 — d_xy (red) sits below zero: the worst of the four",
    "δ = 0.2 — it crosses zero",
    "δ = 0.3 — it passes d_x²−y²",
    "δ ≥ 0.4 — d_xy is the top curve, the strongest attractive channel",
    "every curve starts at exactly zero at U = 0",
    "all peak near U = 4–6, not at the strongest coupling"], size=11, space=7)
tb(s, 0.5, 5.75, 9.07, 0.4,
   "The whole reordering happens over a narrow window in anisotropy.", size=12.5, italic=True)
notes(s, "Point at the red curve rising through the panels. The peak at U=4-6 rather than U=8 is "
         "worth a sentence: attraction is strongest at intermediate coupling.")

# ============ 8. REVERSAL ============
s = slide(); head(s, "Anisotropy does opposite things at the two fillings",
                  "pairing vertex vs filling, one panel per anisotropy, U = 4")
pic_h(s, "G2.png", 0.5, 1.5, 4.05)          # aspect 2.238 -> 9.07 wide
card(s, 9.75, 1.5, 3.05, 1.92, WARM, WARML)
tb(s, 9.97, 1.64, 2.65, 0.3, "At half filling — δ helps", size=12.5, bold=True, color=RED)
tb(s, 9.97, 2.0, 2.65, 1.35,
   "d_xy vertex/N goes −0.010 → +0.039 as δ goes 0 → 0.4. Repulsive becomes strongly attractive.",
   size=11)
card(s, 9.75, 3.63, 3.05, 1.92, COOL, COOLL)
tb(s, 9.97, 3.77, 2.65, 0.3, "When doped — δ hurts", size=12.5, bold=True, color=TEAL)
tb(s, 9.97, 4.13, 2.65, 1.35,
   "At n = 0.778 the same anisotropy drives d_xy from −0.012 to −0.042, and the trough reaches "
   "−0.23 at n = 0.5, δ = 0.7.", size=11)
tb(s, 0.5, 5.75, 9.07, 0.4,
   "Same knob, opposite sign at the two ends of the filling axis.", size=12.5, italic=True)
notes(s, "Direct answer to whether dxy is suppressed away from half filling. Yes - and it gets "
         "worse with delta, not better.")

# ============ 8b. EXTENDED S AND d_x2-y2 ============
s = slide(); head(s, "The other two channels: extended s and d_x²−y²",
                  "both attractive near half filling, both weakened by anisotropy — but not equally")
pic_h(s, "sd.png", 0.48, 1.52, 2.80)          # aspect 4.42 -> 12.37 wide, fits
cc = [("Both decay with δ", "extended s falls ×1.8 from δ = 0 to 0.7 (0.0364 → 0.0201). "
                            "d_x²−y² falls ×3.8 (0.0277 → 0.0072) — twice as fast.", GRN2),
      ("Both need half filling", "extended s turns attractive above n = 0.733, d_x²−y² above "
                                 "n = 0.530. Anisotropy pushes both windows toward n = 1.", BLU2),
      ("Extended s leads — but only at n = 1", "the gap χ_s-ext − χ_d is +0.009 at δ = 0 and grows "
                                               "to +0.015. Away from half filling d_x²−y² leads.", NAVY)]
x = M
for a, b, col in cc:
    card(s, x, 4.72, 4.05, 1.5)
    tb(s, x + 0.24, 4.86, 3.6, 0.34, a, size=12.5, bold=True, color=col)
    tb(s, x + 0.24, 5.24, 3.6, 0.92, b, size=11, color=INK)
    x += 4.19
tb(s, M, 6.42, 12.2, 0.4,
   "By P/P_bub at δ = 0 the two cross near U = 6: extended s leads at weak coupling, d_x²−y² at strong "
   "(1.21 vs 1.26 at U = 8).", size=12, italic=True, color=MUTE)
notes(s, "These are the two conventional channels. Both behave the way White's square-lattice work "
         "predicts. The point of the slide is that anisotropy hurts dx2-y2 twice as fast as extended s, "
         "which is what clears the way for dxy to take over.")

# ============ 8c. WHY EXTENDED S LEADS ============
s = slide(); head(s, "Why extended s leads, and where to be careful",
                  "the one place our δ = 0 result departs from White — and it has an explanation")
card(s, M, 1.6, 6.05, 2.35, WARM, WARML)
tb(s, M + 0.24, 1.74, 5.6, 0.32, "The apparent disagreement", size=13.5, bold=True, color=RED)
tb(s, M + 0.24, 2.12, 5.6, 1.7,
   "White finds d_x²−y² the most attractive channel on the square lattice. We reproduce all four of "
   "his channel SIGNS, but our ordering differs: extended s is larger (0.0364 against 0.0277 at "
   "U = 4, δ = 0).", size=12)
card(s, M + 6.3, 1.6, 6.05, 2.35, COOL, COOLL)
tb(s, M + 6.54, 1.74, 5.6, 0.32, "The explanation, from Huang et al.", size=13.5, bold=True, color=TEAL)
tb(s, M + 6.54, 2.12, 5.6, 1.7,
   "They show extended s dominates the LOCAL pairing while d_x²−y² dominates the LONG-RANGE part. "
   "Our χ integrates over all separations, so the local piece wins on magnitude. Separating them "
   "needs the real-space P(R) decomposition we have not built yet.", size=12)
tb(s, M, 4.15, 12.2, 0.34, "Reliability of these two channels — both need care",
   size=14, bold=True, color=NAVY)
t = [["channel", "seed noise", "size collapse", "shell sensitivity", "verdict"],
     ["on-site s", "0.3%", "3%", "0.2× (none)", ("solid", GREEN, False)],
     ["extended s", "1.3%", "3%", ("16× at boundary", RED, True), ("check twist", GOLD, False)],
     [("d_x²−y²", INK, True), ("20.3%", RED, True), ("13–58%", RED, True), ("3.5×", GOLD, False),
      ("needs seeds", RED, True)],
     ["d_xy", "5.2%", "22%", "1.0× (none)", ("solid", GREEN, False)]]
table(s, M, 4.55, 12.2, t, [2.3, 2.3, 2.5, 2.9, 2.2], rowh=0.40, fs=11.5)
tb(s, M, 6.55, 12.2, 0.4,
   "d_x²−y² is the weakest channel in the dataset on every measure. Its sign and its decay with δ "
   "are safe; precise values are not.", size=12, italic=True, color=MUTE)
notes(s, "Be upfront that dx2-y2 is our noisiest channel. The headline dxy result is the cleanest "
         "one, which is fortunate. If asked what would fix it: more seeds, and the real-space "
         "decomposition to separate local from long-range.")

# ============ 9. FILLING BOUNDARY ============
s = slide(); head(s, "The effect switches off below n ≈ 0.85",
                  "d_xy vertex at L = 12, U = 4, as δ goes 0 → 0.7")
t1 = [["filling n", "d_xy with δ", "verdict"],
      ["0.500", "+0.024 → −0.228", ("reversed", TEAL, False)],
      ["0.556", "+0.018 → −0.180", ("reversed", TEAL, False)],
      ["0.667", "−0.010 → −0.082", ("no swap", MUTE, False)],
      ["0.778", "−0.012 → −0.013", ("no swap", MUTE, False)],
      [("0.889", INK, True), ("−0.013 → +0.021", INK, True), ("SWAP", RED, True)],
      [("1.000", INK, True), ("−0.010 → +0.026", INK, True), ("SWAP", RED, True)]]
table(s, M, 1.7, 6.0, t1, [1.5, 2.7, 1.8], rowh=0.40)
tb(s, 7.1, 1.7, 5.7, 0.3, "Filling at which each channel turns attractive",
   size=13.5, bold=True, color=NAVY)
t2 = [["δ", "extended s", "d_x²−y²", "d_xy"],
      ["0.0", "0.733", "0.530", "0.628"],
      ["0.2", "0.712", "0.570", "0.942"],
      ["0.4", "0.716", "0.729", "0.852"],
      ["0.7", "0.790", "0.752", "0.820"]]
table(s, 7.1, 2.1, 5.7, t2, [1.05, 1.65, 1.5, 1.5], rowh=0.40)
tb(s, 7.1, 4.05, 5.7, 0.85,
   "Anisotropy pushes every channel's attractive window toward half filling. d_x²−y² turns "
   "attractive from n = 0.53 at δ = 0, but only from n = 0.75 at δ = 0.7.", size=12, color=MUTE)
card(s, M, 4.9, W - 2 * M, 1.15)
tb(s, M + 0.28, 5.06, W - 2 * M - 0.56, 0.85,
   "The pairing boundary sits between n = 0.778 and n = 0.889 — the same place the magnetic moment "
   "becomes strong. The two hang together: the altermagnetic correlations that drive the d_xy "
   "attraction only exist near half filling. The claim should be scoped to n ≥ 0.85 explicitly, "
   "with the doped panels shown as the control.", size=12.5)
notes(s, "Scoping to near half filling is a strength - a falsifiable prediction with a control "
         "region where the effect vanishes.")

# ============ 9b. PHASE DIAGRAM ============
s = slide(); head(s, "(n, δ) phase diagram — every channel, every interaction",
                  "L = 12, one shared colour scale across all twenty panels")
pic_h(s, "pd.png", 0.35, 1.42, 5.85)         # aspect 0.956 -> 5.59 wide
cc = [("Read the columns", "Each column is one channel. Red is attractive, blue repulsive, "
                           "white is no effect. The whole U = 0 row is white — the null test, "
                           "seen directly.", NAVY),
      ("Read the rows", "The attractive region grows and pushes leftward as U increases, then "
                        "recedes slightly at U = 8. All four channels are repulsive on the "
                        "dilute side.", TEAL),
      ("The d_xy corner", "Top-right of the fourth column: attractive only where δ is large AND "
                          "n ≥ 0.89. That corner is the paper's result, and nothing else in the "
                          "figure looks like it.", RED)]
y = 1.5
for a, b, col in cc:
    card(s, 6.35, y, 6.45, 1.62)
    tb(s, 6.6, y + 0.14, 6.0, 0.32, a, size=13.5, bold=True, color=col)
    tb(s, 6.6, y + 0.5, 6.0, 1.0, b, size=12)
    y += 1.74
card(s, 6.35, 6.72, 6.45, 0.62, AMB, AMBL)
tb(s, 6.6, 6.83, 6.0, 0.44,
   "Colour scale is symmetric-log: on-site s reaches −0.59 while the d channels live at ±0.04, "
   "so a linear shared scale would flatten them to white.", size=10.5, color=GOLD)
notes(s, "This replaces the old scatter version and answers meeting point 5. One shared scale was "
         "requested; it needs to be symlog because the channels differ by more than a decade in "
         "magnitude. If asked, the linear version is available but the d channels vanish in it.")

# ============ 10. SUPERVISOR QUESTION ============
s = slide(); head(s, "Your question: does d_x²−y² dominate over extended s?",
                  "full P against the uncorrelated P_bub — the White et al. comparison, at δ = 0")
pic_h(s, "fullsus.png", 0.5, 1.48, 3.35)    # aspect 2.143 -> 7.18 wide
card(s, 8.0, 1.48, 4.8, 3.35, WARM, WARML)
tb(s, 8.24, 1.62, 4.35, 0.3, "The answer has a twist", size=13.5, bold=True, color=RED)
tb(s, 8.24, 1.98, 4.35, 2.7,
   "d_xy has the LARGEST full susceptibility of the four — yet its vertex is negative.\n\n"
   "Its bubble alone is 47.67, so it looks dominant for purely kinematic reasons while the "
   "interaction pushes against it.\n\n"
   "Ranking channels by P alone inverts the physical answer. That is exactly why White "
   "introduced the uncorrelated susceptibility.", size=11.5)
t = [["channel", "full P", "bubble P_bub", "vertex", "interaction"],
     ["on-site s", "3.85", "8.34", "−4.49", ("repulsive", TEAL, False)],
     ["extended s", "30.83", "25.58", "+5.25", ("ATTRACTIVE", RED, True)],
     ["d_x²−y²", "35.00", "31.02", "+3.98", ("ATTRACTIVE", RED, True)],
     [("d_xy", INK, True), ("46.25", INK, True), "47.67", "−1.42", ("repulsive", TEAL, False)]]
table(s, M, 5.05, 6.4, t, [1.5, 1.25, 1.25, 1.2, 1.2], rowh=0.36, fs=11)
tb(s, 7.35, 5.12, 5.45, 1.5,
   "By vertex the ordering is stable at every U: extended s and d_x²−y² attractive, on-site s and "
   "d_xy repulsive — reproducing White's channel signs on a lattice he never studied.",
   size=12, color=MUTE)
notes(s, "By vertex the ordering is stable at every U: extended s and dx2-y2 attractive, on-site s "
         "and dxy repulsive - reproducing White's signs on a lattice he never studied.")

# ============ 11. RATIO ============
s = slide(); head(s, "All four channels on one scale",
                  "P / P_bub divides out the kinematic bubble; above 1 the interaction is attractive")
s.shapes.add_picture(os.path.join(IMG, "ratio.png"), I(0.9), I(1.6), width=I(11.5))
tb(s, M, 6.35, 12.2, 0.62,
   "All four curves start at exactly 1.000 at U = 0 — the null test, made visible.  At U = 8 the "
   "interaction enhances d_x²−y² by 26% and extended s by 21%, while suppressing d_xy by 8% and "
   "on-site s by 46%.  Extended s crosses the line at n = 0.733.", size=12.5)
notes(s, "The U=0 flat line at exactly 1 is worth pausing on - the null test as a visual anchor.")

# ============ 12. MAGNETISM ============
s = slide(); head(s, "Meeting point 7: m_AM against U and δ",
                  "the two contributions separated — they pull in opposite directions")
s.shapes.add_picture(os.path.join(IMG, "mam.png"), I(0.95), I(1.5), width=I(7.5))
cc = [("U builds the moment", "Δm rises roughly linearly with U at every anisotropy: 0.0055 at "
                              "U = 2 to 0.0257 at U = 8.", RED),
      ("δ destroys it", "At U = 8 the same Δm falls by a factor of three, 0.0257 → 0.0086, as δ "
                        "goes 0 → 0.7.", TEAL),
      ("The product plateaus", "m_AM = Δm·δ climbs to δ ≈ 0.3 then flattens — two effects "
                               "cancelling, not saturation.", NAVY)]
y = 1.6
for a, b, col in cc:
    card(s, 8.7, y, 4.08, 1.5)
    tb(s, 8.92, y + 0.13, 3.7, 0.3, a, size=13, bold=True, color=col)
    tb(s, 8.92, y + 0.46, 3.7, 0.95, b, size=11.5)
    y += 1.62
card(s, 8.7, 6.5, 4.08, 0.62, AMB, AMBL)
tb(s, 8.92, 6.62, 3.7, 0.42,
   "Caveat: U = 6 and U = 8 are within 1σ of each other at δ ≥ 0.6.", size=11, color=GOLD)
notes(s, "Be honest: the U dependence is established up to delta about 0.4 and unresolved beyond, "
         "where the top two curves touch. Delta m is also still shrinking with lattice size.")

# ============ 13. Q* ============
s = slide(); head(s, "A correction we found ourselves",
                  "the magnetic peak leaves (π,π) — the moment was being read at the wrong wavevector")
s.shapes.add_picture(os.path.join(IMG, "sqmap3.png"), I(M), I(1.5), width=I(12.2))
s.shapes.add_picture(os.path.join(IMG, "qstar.png"), I(M), I(3.9), width=I(8.3))
t = [["δ", "peak q*", "S(π,π)/S(q*)"],
     ["0.0", "(π, 0.83π)", "0.927"],
     ["0.3", "(π, π)", "0.981"],
     ["0.5", "(π, π)", "0.944"],
     ["0.7", "(0, π)", ("0.804", RED, True)]]
table(s, 9.1, 3.95, 3.68, t, [0.88, 1.5, 1.3], rowh=0.36, fs=11)
tb(s, 9.1, 5.75, 3.68, 1.3,
   "Incommensurate at small δ, commensurate Néel near δ ≈ 0.35, stripe-like at large δ. "
   "Recomputing Δm at the true peak leaves the suppression intact — factor 2.3 instead of 3.0 — "
   "so the headline survives. Both limits match published work.", size=10.5)
notes(s, "Stress that we found this by checking our own order parameter, not because a referee "
         "asked - and the conclusion survived.")

# ============ 14. OPEN SHELL ============
s = slide(); head(s, "Meeting point 11: closed-shell filling",
                  "answered — and it turned out to matter")
t = [["lattice", "fillings closed-shell at every δ"],
     ["L = 8", "only n = 1.72, 1.84, 1.97"],
     [("L = 10", INK, True), ("n = 0.500, 0.740, 1.020", INK, True)],
     ["L = 12", ("none below n = 1.65", RED, True)]]
table(s, M, 1.7, 6.0, t, [1.7, 4.3], rowh=0.42)
tb(s, M, 3.5, 6.0, 0.6,
   "At L = 12 the open-shell problem cannot be fixed by choosing a different filling — zero of its "
   "48 (filling, δ) combinations is closed shell.", size=12, color=MUTE)
card(s, 7.05, 1.7, 5.73, 2.45, WARM, WARML)
tb(s, 7.3, 1.82, 5.25, 0.34, "Does it actually matter?  Yes, in two channels",
   size=13, bold=True, color=RED)
tb(s, 7.3, 2.22, 5.25, 1.8,
   "Testing at fixed L and fixed n, with the shell character flipping as δ varies:\n\n"
   "•  extended s jumps 16× its typical step at the boundary\n"
   "•  d_x²−y² jumps 3.5×\n"
   "•  on-site s (0.2×) and d_xy (1.0×) show nothing", size=12)
card(s, M, 4.35, W - 2 * M, 2.1)
tb(s, M + 0.28, 4.5, 9.0, 0.32,
   "The fix: twist-averaged boundary conditions — running now, results tomorrow",
   size=14, bold=True, color=NAVY)
bullets(s, M + 0.28, 4.88, W - 2 * M - 0.56, 1.5, [
    "A twist is a boundary phase ψ(r + L x̂) = e^{iθx} ψ(r), shifting the k grid off the awkward points.",
    "Only θ = 0 and π keep the hopping real, which the engine requires — so no rewrite is needed.",
    "The two mixed twists are degenerate by the lattice's x ↔ y symmetry (eigenvalues agree to "
    "1e-14), so three runs with weights 1 : 2 : 1 give the four-point average.",
    "Our headline channel d_xy shows no shell sensitivity, so the main claim is already protected."],
    size=11.5, space=3)
notes(s, "The empirical test matters: we did not assume open shells were harmless, we measured a "
         "case where the shell character flips at fixed L and n, and found real jumps in two of four channels.")

# ============ 15. LITERATURE ============
s = slide(); head(s, "Literature check: the field has moved",
                  "two 2025–26 papers already claim altermagnetism on this lattice family")
t = [["paper", "model / method", "overlap with us"],
     ["arXiv 2607.06106  (Jul 2026)", "checkerboard Hubbard, mean field + RPA",
      ("same lattice, weak coupling", RED, False)],
     ["arXiv 2503.08362", "t–t′–δ Hubbard, mean field + DMRG",
      ("same anisotropy; no pairing", GOLD, False)],
     ["arXiv 2505.12342  (our group)", "spin-anisotropic Hubbard, CPQMC", "anisotropy suppresses AFM"],
     ["arXiv 2603.29377  (our group)", "anisotropic t–J, CPQMC", "d + p mixed pairing"],
     ["White et al., PRB 39 839", "square Hubbard, DQMC", "our δ = 0 benchmark"],
     ["Huang, Lin, Gubernatis, PRB 64 205101", "t–t′–U, CPQMC at T = 0", "local vs long-range vertex"]]
table(s, M, 1.7, 12.2, t, [4.4, 4.2, 3.6], rowh=0.42, fs=11)
card(s, M, 4.8, 6.05, 1.9, COOL, COOLL)
tb(s, M + 0.24, 4.93, 5.6, 0.3, "What this costs us", size=13.5, bold=True, color=TEAL)
tb(s, M + 0.24, 5.26, 5.6, 1.3,
   "“Checkerboard hosts altermagnetism” is no longer available as a headline — 2607.06106 is five "
   "weeks old and says exactly that. Interaction-driven altermagnetism is an established category now.",
   size=12)
card(s, M + 6.3, 4.8, 6.05, 1.9, WARM, WARML)
tb(s, M + 6.54, 4.93, 5.6, 0.3, "What it gives us", size=13.5, bold=True, color=RED)
tb(s, M + 6.54, 5.26, 5.6, 1.3,
   "The t–t′–δ paper explicitly declines to compute pairing. Our d_xy channel swap is unclaimed, "
   "it is CPQMC at intermediate coupling where the competitors are mean field, and we have the "
   "doping axis nobody else has.", size=12)
notes(s, "Recommend reframing: lead with pairing, use magnetism as the explanation. Stronger "
         "position than the original plan, but it means rewriting the introduction.")

# ============ 16. MEETING POINTS ============
s = slide(); head(s, "Status against the eleven meeting points",
                  "four done, one newly answered, three substantially done, three outstanding")
left = [("1", "δ to 0.7", "DONE", GREEN),
        ("4", "s-channels: trivial or not", "DONE", GREEN),
        ("6", "spectral function in CPQMC", "DONE", GREEN),
        ("7", "m_AM vs U and δ", "DONE", GREEN),
        ("11", "closed-shell filling", "ANSWERED", GREEN),
        ("5", "heatmaps not scatter", "mostly", GOLD)]
right = [("9", "equal-time vertex", "resolved, unwritten", GOLD),
         ("10", "strengthen the novelty", "answer changed", GOLD),
         ("2", "twist-averaged BC", "running now", GOLD),
         ("3", "periodic BC, TRIM", "partial", GOLD),
         ("8", "write the manuscript", "fragments only", RED),
         ("—", "real-space P(R)", "not started", RED)]
for rows, x0 in ((left, M), (right, M + 6.3)):
    y = 1.68
    for num, lab, st, col in rows:
        card(s, x0, y, 6.05, 0.64, C(0xF7, 0xF9, 0xFD), CARDL)
        badge(s, num, x0 + 0.14, y + 0.14, col, d=0.36)
        tb(s, x0 + 0.62, y + 0.17, 3.28, 0.32, lab, size=12.5)
        tb(s, x0 + 4.02, y + 0.17, 1.88, 0.32, st, size=11.5, bold=True, color=col,
           align=PP_ALIGN.RIGHT)
        y += 0.74
tb(s, M, 6.35, 12.2, 0.6,
   "Point 6 closed with a decision rather than a measurement: G(k,τ) is measurable, but the usable "
   "window is τ ≈ 0.4–0.8, too short for analytic continuation. DQMC would cost 1–2 months and "
   "inherits the sign problem CPQMC was chosen to avoid.", size=11.5, color=MUTE, italic=True)
notes(s, "Point 6: propose Hartree-Fock band structure as the cheap way to display the altermagnetic "
         "spin splitting, since neither CPQMC nor DQMC shows it in a spin-symmetric simulation.")

# ============ 17. NEXT ============
s = slide(NAVY)
tb(s, M, 0.62, 6.0, 0.8, "Next", size=36, font=HF, bold=True, color=WHITE)
items = [("Tomorrow", "Twist average completes — confirm the d_xy swap is twist-robust", ICE),
         ("This week", "Real-space pair correlations P(R): the code does not exist, and without it "
                       "we cannot speak to long-range order at all", ICE),
         ("This week", "More seeds for d_x²−y² — 20% seed noise against 0.3–5% elsewhere", ICE),
         ("Decision", "Reframe the paper around pairing rather than magnetism?", C(0xFF, 0xD9, 0x8A))]
y = 1.85
for i, (a, b, col) in enumerate(items):
    badge(s, i + 1, M, y, GOLD if i == 3 else C(0x3A, 0x4C, 0x8F))
    tb(s, M + 0.62, y + 0.02, 2.5, 0.32, a, size=13, bold=True, color=col)
    tb(s, M + 3.2, y - 0.02, 9.0, 0.66, b, size=13.5, color=WHITE)
    y += 0.95
card(s, M, 5.75, W - 2 * M, 1.2, C(0x16, 0x20, 0x4A), C(0x3A, 0x4C, 0x8F))
tb(s, M + 0.28, 5.88, 8.0, 0.3, "Honest limitation to carry into the manuscript",
   size=13, bold=True, color=C(0xFF, 0xD9, 0x8A))
tb(s, M + 0.28, 6.2, W - 2 * M - 0.56, 0.68,
   "χ/N is flat in system size, so there is no tendency toward long-range pairing order. What we "
   "can defend is an interaction-driven change in which channel is attractive — not superconductivity. "
   "Huang, Lin and Gubernatis reach the same conclusion for t–t′–U.", size=12.5, color=ICE)
notes(s, "End on the limitation deliberately. Saying plainly what the data cannot support makes "
         "the positive claims more credible.")

out = os.path.join(HERE, "weekly_progress.pptx")
prs.save(out)
print("written:", out, os.path.getsize(out) // 1024, "KB,", len(prs.slides.__iter__.__self__._sldIdLst), "slides")
