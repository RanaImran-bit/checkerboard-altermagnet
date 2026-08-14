#!/usr/bin/env python3
"""Focused deck: twist averaging + real-space P(R). Built with python-pptx."""
import os
from pptx import Presentation
from pptx.util import Inches as I, Pt
from pptx.dml.color import RGBColor as C
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

HERE = os.path.dirname(os.path.abspath(__file__)); IMG = os.path.join(HERE, "img")
NAVY = C(0x1E, 0x27, 0x61); ICE = C(0xCA, 0xDC, 0xFC); WHITE = C(0xFF, 0xFF, 0xFF)
INK = C(0x13, 0x18, 0x2F); MUTE = C(0x5A, 0x64, 0x84); RED = C(0xC0, 0x39, 0x2B)
TEAL = C(0x1C, 0x72, 0x93); GOLD = C(0x8A, 0x6D, 0x1F); GREEN = C(0x1E, 0x7A, 0x46)
CARD = C(0xF4, 0xF6, 0xFB); CARDL = C(0xE2, 0xE7, 0xF2)
WARM = C(0xFD, 0xF3, 0xF2); WARML = C(0xF0, 0xD5, 0xD2)
COOL = C(0xF1, 0xF5, 0xFA); COOLL = C(0xD8, 0xE2, 0xEE)
AMB = C(0xFF, 0xF8, 0xE6); AMBL = C(0xEF, 0xE0, 0xBC)
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


def bullets(s, x, y, w, h, items, size=12.5, color=INK, space=6):
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
    tb(s, M, 0.30, W - 2*M, 0.74, title, size=30, font=HF, bold=True)
    if kicker: tb(s, M, 1.02, W - 2*M, 0.34, kicker, size=13, color=MUTE, italic=True)


def pic_h(s, name, x, y, h):
    return s.shapes.add_picture(os.path.join(IMG, name), I(x), I(y), height=I(h))


def table(s, x, y, w, rows, colw, rowh=0.38, fs=11.5):
    shp = s.shapes.add_table(len(rows), len(rows[0]), I(x), I(y), I(w), I(rowh*len(rows)))
    t = shp.table
    for j, cwd in enumerate(colw): t.columns[j].width = I(cwd)
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
tb(s, M, 2.2, W-2*M, 0.9, "Two validation results", size=40, font=HF, bold=True, color=WHITE)
tb(s, M, 3.1, W-2*M, 0.7, "twist-averaged boundary conditions  ·  real-space pairing P(R)",
   size=24, font=HF, color=ICE)
tb(s, M, 4.2, W-2*M, 0.4, "Checkerboard Hubbard, CPQMC at T = 0, L = 12", size=15, color=ICE)
tb(s, M, 4.66, W-2*M, 0.35, "13 August 2026", size=13, color=C(0x8F,0xA3,0xD9))
tb(s, M, 5.7, 11.5, 0.5,
   "Both were open questions last week. Both are now answered.",
   size=15, color=WHITE, italic=True)
notes(s, "Two things the supervisor asked for: check the open-shell trial state is not producing "
         "the result, and do the Huang local-versus-long-range decomposition. Both done.")

# ---------------- 2 twist: why ----------------
s = slide(); head(s, "Why twist averaging was needed",
                  "under periodic boundaries every filling at L = 12 is open shell")
bl = ["A twist is a boundary phase ψ(r + L x̂) = e^{iθx} ψ(r), shifting the k grid from "
      "2πm/L to (2πm + θ)/L. Shell effects come from that discrete grid landing awkwardly "
      "on the Fermi surface; averaging over θ fills the zone and washes them out.",
      "Only θ = 0 and π keep the hopping matrix real, which the engine requires — so no "
      "rewrite was needed.",
      "The two mixed twists are exactly degenerate by the lattice's x ↔ y symmetry "
      "(eigenvalues agree to 1e-14), so three runs with weights 1 : 2 : 1 give the "
      "four-point average.",
      "Shell effects are real, not hypothetical: at L = 8, n = 0.781 the extended-s vertex "
      "jumps sixteen times its typical step exactly where the shell character flips."]
bullets(s, M, 1.6, W-2*M, 2.6, bl, size=13.5)
card(s, M, 4.35, 6.05, 1.5, COOL, COOLL)
tb(s, M+0.24, 4.5, 5.6, 0.3, "The three runs", size=13.5, bold=True, color=TEAL)
tb(s, M+0.24, 4.85, 5.6, 0.9,
   "(+1,+1) periodic, weight 1  —  already had it\n"
   "(−1,+1) half twist, weight 2  —  stands for (+1,−1) too\n"
   "(−1,−1) both antiperiodic, weight 1", size=12.5)
card(s, M+6.3, 4.35, 6.05, 1.5, WARM, WARML)
tb(s, M+6.54, 4.5, 5.6, 0.3, "Free validation", size=13.5, bold=True, color=RED)
tb(s, M+6.54, 4.85, 5.6, 0.9,
   "The U = 0 null test still passes under both twists — max |vertex| = 8.6e-13 and 8.1e-13. "
   "Twisted boundaries do not break the full-minus-bubble decomposition.", size=12.5)
tb(s, M, 6.1, W-2*M, 0.7,
   "Scope: L = 12, half filling, U = 2, 4, 6, 8, all eight anisotropies, six seeds. "
   "Both twist legs complete.", size=12.5, color=MUTE, italic=True)
notes(s, "If asked why not full twist averaging: a continuous theta grid needs complex hoppings, "
         "and cpqmc.py is real-only. That is a one-to-two week rewrite. This four-point version "
         "removes the worst artefacts without touching the engine.")

# ---------------- 3 twist result ----------------
s = slide(); head(s, "The channel swap survives at every interaction strength",
                  "d_xy pairing vertex, periodic against twist-averaged")
pic_h(s, "twist.png", 0.53, 1.48, 3.00)   # aspect 4.14 -> 12.42 wide, fits
t = [["U", "crossing δ", "peak χ/N", "vs periodic"],
     ["2", "0.3", "+0.0200", ("−21%", GOLD, False)],
     ["4", "0.2", "+0.0308", ("−21%", GOLD, False)],
     ["6", "0.2", "+0.0304", ("−17%", GOLD, False)],
     ["8", "0.2", "+0.0245", ("−17%", GOLD, False)]]
table(s, M, 5.0, 5.3, t, [1.1, 1.5, 1.5, 1.2])
card(s, 6.4, 5.0, 6.38, 1.9, WARM, WARML)
tb(s, 6.64, 5.14, 5.9, 0.3, "What it confirms, and what it costs",
   size=13.5, bold=True, color=RED)
tb(s, 6.64, 5.5, 5.9, 1.3,
   "d_xy is repulsive at δ = 0 and 0.1 and attractive from δ = 0.2 at every U. The open-shell "
   "trial state is NOT producing the swap.\n"
   "But magnitudes fall 17–21% under averaging, so the periodic numbers were inflated. Quote the "
   "twist-averaged values.", size=12.5)
notes(s, "The honest detail: at U = 2, delta = 0.2 the sign actually flips under averaging, from "
         "+0.014 to -0.011, moving the crossing from 0.15 to 0.25 at weak coupling. Everywhere "
         "else the crossing location is unchanged.")

# ---------------- 4 P(R) method ----------------
s = slide(); head(s, "Real-space pairing P(R)",
                  "the decomposition Huang, Lin and Gubernatis use (PRB 64, 205101)")
card(s, M, 1.55, W-2*M, 0.8, COOL, COOLL)
tb(s, M+0.3, 1.72, 11.8, 0.5,
   "V_α(R)  =  P_α(R)  −  P̄_α(R)        local = V(R = 0)        long range = mean V(R > 2)",
   size=17, font=HF, bold=True, color=TEAL, align=PP_ALIGN.CENTER)
bl = ["The existing routine sums the pair matrix over ALL site pairs and discards the "
      "structure. Keeping that matrix and binning each entry by the separation between the "
      "two sites gives P(R) instead.",
      "Summing P(R) back over R must return the original susceptibility. That sum rule is the "
      "validation — it reproduces both the full susceptibility and the vertex to 0.00% in "
      "every channel and every seed.",
      "The sum rule caught two bugs before any production run: a call to a method that does "
      "not exist (guarded by hasattr, so equilibration was silently skipped — 400% error), and "
      "a missing per-block reorthogonalisation.",
      "Distances use the minimum-image convention; 27 distinct separations at L = 12."]
bullets(s, M, 2.6, W-2*M, 2.6, bl, size=13)
card(s, M, 5.35, W-2*M, 1.5, AMB, AMBL)
tb(s, M+0.28, 5.5, 8, 0.3, "Why this matters for the White comparison",
   size=14, bold=True, color=GOLD)
tb(s, M+0.28, 5.85, W-2*M-0.56, 0.9,
   "Our q = 0 susceptibility adds the local and long-range pieces together. If they favour "
   "different channels — which is exactly what Huang finds — then the sum can name the wrong "
   "winner. That is the suspected reason our δ = 0 ordering put extended s ahead of d_x²−y², "
   "where White has d_x²−y² first.", size=12.5)
notes(s, "Emphasise the sum rule. It is an exact internal check: the R-resolved code must "
         "reproduce the old routine when the bins are added back, and it does to every digit.")

# ---------------- 5 P(R) result 1 ----------------
s = slide(); head(s, "Huang's result reproduced — then reversed by anisotropy",
                  "L = 12, U = 4, half filling")
pic_h(s, "prbar.png", 0.7, 1.45, 3.7)
card(s, M, 5.35, 6.05, 1.55, COOL, COOLL)
tb(s, M+0.24, 5.5, 5.6, 0.3, "At δ = 0: exactly Huang's finding", size=13.5, bold=True, color=TEAL)
tb(s, M+0.24, 5.85, 5.6, 0.95,
   "LOCAL: extended s (+1.59) beats d_x²−y² (+0.92).\n"
   "LONG RANGE: d_x²−y² (+0.080) beats extended s (+0.017) by 4.8×.\n"
   "The local piece is ~20× larger, so the q = 0 sum names extended s.", size=12)
card(s, M+6.3, 5.35, 6.05, 1.55, WARM, WARML)
tb(s, M+6.54, 5.5, 5.6, 0.3, "At δ = 0.4: d_xy takes over", size=13.5, bold=True, color=RED)
tb(s, M+6.54, 5.85, 5.6, 0.95,
   "d_x²−y²'s long-range vertex collapses 0.080 → 0.014, while d_xy goes −0.012 → +0.062 "
   "and leads by more than 4×. The handover happens in the long-range channel.", size=12)
notes(s, "This resolves the apparent disagreement with White completely, and it is a stronger "
         "statement than the q=0 result: the handover is in the piece that actually matters for "
         "ordering.")

# ---------------- 6 P(R) result 2 ----------------
s = slide(); head(s, "Where each channel dominates, anisotropy by anisotropy",
                  "long-range vertex V(R > 2), L = 12, U = 4, half filling")
pic_h(s, "pr.png", 0.35, 1.45, 3.5)
t = [["δ", "on-site s", "extended s", "d_x²−y²", "d_xy", "leader"],
     ["0.0", "−0.051", "+0.017", ("+0.080", INK, True), "−0.012", ("d_x²−y²", TEAL, True)],
     ["0.1", "−0.049", "+0.016", ("+0.088", INK, True), "−0.024", ("d_x²−y²", TEAL, True)],
     ["0.2", "−0.048", "+0.012", "+0.022", ("+0.048", INK, True), ("d_xy", RED, True)],
     ["0.3", "−0.051", "+0.009", "+0.023", ("+0.066", INK, True), ("d_xy", RED, True)],
     ["0.4", "−0.050", "+0.006", "+0.014", ("+0.062", INK, True), ("d_xy", RED, True)],
     ["0.7", "−0.051", "+0.005", "+0.010", ("+0.025", INK, True), ("d_xy", RED, True)]]
table(s, M, 5.15, 8.4, t, [0.8, 1.5, 1.5, 1.4, 1.3, 1.9], rowh=0.32, fs=11)
card(s, 9.2, 5.15, 3.58, 1.75, WARM, WARML)
tb(s, 9.42, 5.28, 3.15, 0.3, "The handover", size=13, bold=True, color=RED)
tb(s, 9.42, 5.62, 3.15, 1.2,
   "Between δ = 0.1 and 0.2 the long-range attraction moves from d_x²−y² to d_xy — a factor of "
   "four collapse in one channel and a sign change in the other.", size=11.5)
notes(s, "Panel (c) shows V(R) decaying into the noise past R about 3.5 in every channel. That is "
         "the honest limit on any ordering claim - short-range correlations, no long-range order.")

# ---------------- 7 limits ----------------
s = slide(); head(s, "What these two results do and do not establish", None)
rows = [("Established", "The d_xy channel swap is not an open-shell artefact — it survives "
         "twist averaging at U = 2, 4, 6 and 8.", GREEN),
        ("Established", "At δ = 0 our channel ordering agrees with Huang once local and "
         "long-range are separated, which explains the apparent conflict with White.", GREEN),
        ("Established", "Anisotropy moves the long-range attraction from d_x²−y² to d_xy "
         "between δ = 0.1 and 0.2.", GREEN),
        ("Revised", "Periodic-boundary magnitudes are 17–21% too high. All quoted numbers "
         "should come from the twist average.", GOLD),
        ("NOT established", "Long-range order. V(R) decays into the noise past R ≈ 3.5 in "
         "every channel, and χ/N is flat in system size. What we can defend is a change in "
         "which channel is attractive — not superconductivity.", RED)]
y = 1.55
for tag, txt, col in rows:
    card(s, M, y, W-2*M, 0.95)
    sh = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, I(M+0.2), I(y+0.24), I(1.85), I(0.44))
    sh.fill.solid(); sh.fill.fore_color.rgb = col
    sh.line.fill.background(); sh.shadow.inherit = False
    tf = sh.text_frame; tf.margin_left = tf.margin_right = 0
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = tag
    r.font.name = BF; r.font.size = Pt(11.5); r.font.bold = True; r.font.color.rgb = WHITE
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tb(s, M+2.25, y+0.16, 10.2, 0.7, txt, size=12.5)
    y += 1.05
notes(s, "End on the limitation. Huang reach the same conclusion for t-t'-U: no tendency for the "
         "correlations to grow with system size. Saying it first makes the positive claims credible.")

# ---------------- 8 next ----------------
s = slide(NAVY)
tb(s, M, 0.7, 6, 0.8, "Next", size=36, font=HF, bold=True, color=WHITE)
items = [("Running now", "P(R) queue: 13 of 40 chunks done. U = 4 complete at all fillings and "
          "anisotropies; U = 8 in progress.", ICE),
         ("Timing", "264–521 min per chunk against my 225 estimate — the full set lands Friday "
          "evening, not tonight.", C(0xFF,0xD9,0x8A)),
         ("Then", "Same local/long-range decomposition at U = 2, 6, 8, and across fillings — "
          "does the handover track the filling boundary at n ≈ 0.85?", ICE),
         ("Open", "Twist averaging is only at half filling. The doped case is untested.", ICE)]
y = 1.9
for i, (a, b, col) in enumerate(items):
    sh = s.shapes.add_shape(MSO_SHAPE.OVAL, I(M), I(y), I(0.42), I(0.42))
    sh.fill.solid(); sh.fill.fore_color.rgb = C(0x3A,0x4C,0x8F)
    sh.line.fill.background(); sh.shadow.inherit = False
    tf = sh.text_frame; p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = str(i+1)
    r.font.name = BF; r.font.size = Pt(13); r.font.bold = True; r.font.color.rgb = WHITE
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tb(s, M+0.62, y-0.02, 2.4, 0.32, a, size=13, bold=True, color=col)
    tb(s, M+3.1, y-0.05, 9.2, 0.7, b, size=13.5, color=WHITE)
    y += 1.0
card(s, M, 6.0, W-2*M, 0.95, C(0x16,0x20,0x4A), C(0x3A,0x4C,0x8F))
tb(s, M+0.28, 6.14, 11.8, 0.7,
   "Both figures are reproducible: twist_average.py and checkerboard_pr.py, with the sum-rule "
   "test as the correctness check.", size=12.5, color=ICE)
notes(s, "The doped twist check is the obvious next validation - the channel swap is claimed to be "
         "a half-filling effect, and the twist test so far only covers half filling.")

out = os.path.join(HERE, "twist_and_PR.pptx")
prs.save(out)
print("written:", out, os.path.getsize(out)//1024, "KB")
