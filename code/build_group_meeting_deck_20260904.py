"""Group-meeting deck: three altermagnet projects (NNN paper, combined model, checkerboard).
python-pptx port of build.js (no Node on this machine)."""
import os, sys
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from lxml import etree
from PIL import Image

S = os.path.dirname(os.path.abspath(__file__))
P1 = lambda f: os.path.join(S, "p1", f)
P2 = lambda f: os.path.join(S, "p2", "figure", f)
P3 = lambda f: os.path.join(S, "p3", f)
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(S, "Imran_Group_Meeting_04_09.pptx")

C = dict(
    ink="0F2B2E", teal="1F5F5B", tealSoft="E6F0EF", orange="E07A2F", orangeSoft="FBEADF",
    text="1A1A1A", muted="5C6B6A", white="FFFFFF", p1="6C8E9B", p2="1F5F5B", p3="E07A2F",
    red="B03A2E", green="2E7D5B", pink="F6E9E7", grey="EEF2F3", darkcard="1A3D40",
    lightteal="9FC2BE", paleteal="D5E3E1", palegreen="E3EEEC",
)
def rgb(h): return RGBColor.from_string(h)
HFONT, BFONT = "Cambria", "Calibri"
W, H = 13.333, 7.5

prs = Presentation()
prs.slide_width, prs.slide_height = Inches(W), Inches(H)
BLANK = prs.slide_layouts[6]

# ---------------- helpers ----------------
def new_slide(bg=None):
    s = prs.slides.add_slide(BLANK)
    if bg:
        f = s.background.fill; f.solid(); f.fore_color.rgb = rgb(bg)
    return s

def shape(s, kind, x, y, w, h, fill, radius=None, line=None):
    sh = s.shapes.add_shape(kind, Inches(x), Inches(y), Inches(w), Inches(h))
    sh.fill.solid(); sh.fill.fore_color.rgb = rgb(fill)
    if line: sh.line.color.rgb = rgb(line); sh.line.width = Pt(0.75)
    else: sh.line.fill.background()
    sh.shadow.inherit = False
    if radius is not None and kind == MSO_SHAPE.ROUNDED_RECTANGLE:
        sh.adjustments[0] = radius
    sh.text_frame.text = ""
    return sh

def card(s, x, y, w, h, fill=None):
    return shape(s, MSO_SHAPE.ROUNDED_RECTANGLE, x, y, w, h, fill or C["tealSoft"], radius=0.06)

def _set_bullet(p, level_indent=0.22):
    pPr = p._p.get_or_add_pPr()
    pPr.set("marL", str(int(Inches(level_indent))))
    pPr.set("indent", str(-int(Inches(level_indent))))
    for tag in ("a:buNone", "a:buChar", "a:buAutoNum"):
        for e in pPr.findall(qn(tag)): pPr.remove(e)
    bu = etree.SubElement(pPr, qn("a:buChar")); bu.set("char", "•")

def text(s, runs, x, y, w, h, size=14, font=BFONT, color=None, bold=False, italic=False,
         align="left", valign="top", bullets=False, para_after=6, line_spacing=None):
    """runs: str | list of paragraphs; each paragraph is str or list of (text, opts) tuples."""
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame; tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = {"top": MSO_ANCHOR.TOP, "middle": MSO_ANCHOR.MIDDLE, "bottom": MSO_ANCHOR.BOTTOM}[valign]
    paras = runs if isinstance(runs, list) else [runs]
    for i, para in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = {"left": PP_ALIGN.LEFT, "center": PP_ALIGN.CENTER, "right": PP_ALIGN.RIGHT}[align]
        if bullets: _set_bullet(p)
        if para_after and len(paras) > 1: p.space_after = Pt(para_after)
        if line_spacing: p.line_spacing = line_spacing
        segs = para if isinstance(para, list) else [(para, {})]
        for t, o in segs:
            r = p.add_run(); r.text = t
            f = r.font; f.name = o.get("font", font); f.size = Pt(o.get("size", size))
            f.bold = o.get("bold", bold); f.italic = o.get("italic", italic)
            f.color.rgb = rgb(o.get("color", color or C["text"]))
    return tb

def bullets(s, items, x, y, w, h, size=14, color=None, head_color=None):
    paras = []
    for it in items:
        if isinstance(it, str): paras.append([(it, {})])
        else:
            head, body = it
            paras.append([(head + " ", {"bold": True, "color": head_color or C["teal"]}), (body, {})])
    return text(s, paras, x, y, w, h, size=size, color=color, bullets=True, para_after=6)

def img_fit(s, path, x, y, w, h):
    iw, ih = Image.open(path).size
    ar = iw / ih
    fw, fh = w, w / ar
    if fh > h: fh, fw = h, h * ar
    fx, fy = x + (w - fw) / 2, y + (h - fh) / 2
    s.shapes.add_picture(path, Inches(fx), Inches(fy), Inches(fw), Inches(fh))
    return fx, fy, fw, fh

def badge(s, n, color, label=None):
    shape(s, MSO_SHAPE.OVAL, W - 0.95, 0.32, 0.5, 0.5, color)
    text(s, str(n), W - 0.95, 0.32, 0.5, 0.5, size=18, font=HFONT, bold=True, color=C["white"], align="center", valign="middle")
    if label:
        text(s, label, W - 3.0, 0.32, 1.95, 0.5, size=11, color=C["muted"], align="right", valign="middle")

def title(s, t, dark=False):
    text(s, t, 0.6, 0.25, 9.6, 0.95, size=26, font=HFONT, bold=True, color=C["white"] if dark else C["teal"], valign="middle")

def caption(s, t, x, y, w, h=0.35):
    text(s, t, x, y, w, h, size=10.5, italic=True, color=C["muted"])

def stat(s, big, small, x, y, w, color=None, big_size=40):
    text(s, big, x, y, w, 0.85, size=big_size, font=HFONT, bold=True, color=color or C["teal"], align="center", valign="bottom")
    text(s, small, x, y + 0.88, w, 0.6, size=12, color=C["muted"], align="center")

def notes(s, t):
    s.notes_slide.notes_text_frame.text = t

def section(num, color, kicker, heading, sub, note):
    s = new_slide(C["ink"])
    shape(s, MSO_SHAPE.OVAL, 0.9, 2.3, 1.5, 1.5, color)
    text(s, str(num), 0.9, 2.3, 1.5, 1.5, size=48, font=HFONT, bold=True, color=C["white"], align="center", valign="middle")
    text(s, kicker, 2.9, 2.1, 9.5, 0.5, size=16, color=C["lightteal"])
    text(s, heading, 2.9, 2.6, 9.8, 1.9, size=30, font=HFONT, bold=True, color=C["white"])
    text(s, sub, 2.9, 4.6, 9.5, 1.4, size=15, color=C["paleteal"])
    notes(s, note)
    return s

def num_circle(s, n, col, x, y, d=0.55, size=18):
    shape(s, MSO_SHAPE.OVAL, x, y, d, d, col)
    text(s, n, x, y, d, d, size=size, font=HFONT, bold=True, color=C["white"], align="center", valign="middle")

def lead(head, body, hc):
    return [(head + "  ", {"bold": True, "color": hc}), (body, {})]

# ================= 1. Title =================
s = new_slide(C["ink"])
text(s, "From imposed to emergent altermagnetism", 0.8, 1.5, 11.5, 1.5, size=42, font=HFONT, bold=True, color=C["white"], valign="bottom")
text(s, "Three quantum Monte Carlo studies of spin-split magnetism in Hubbard models", 0.8, 3.1, 11.5, 0.6, size=20, color=C["paleteal"])
for i, (n, col, a, b) in enumerate([("1", C["p1"], "Spin-dependent NNN hopping", "published, PRB 113, 134443"),
                                     ("2", C["p2"], "Combined NN + NNN model", "major revision, lead author"),
                                     ("3", C["p3"], "Checkerboard lattice", "manuscript in preparation, lead author")]):
    x = 0.8 + i * 4.0
    num_circle(s, n, col, x, 4.3, 0.6, 20)
    text(s, a, x + 0.75, 4.25, 3.2, 0.4, size=14, bold=True, color=C["white"])
    text(s, b, x + 0.75, 4.62, 3.2, 0.4, size=12, color=C["lightteal"])
text(s, "Rana Imran Mushtaq   |   Group meeting, 4 September 2026", 0.8, 6.3, 11.5, 0.5, size=14, color=C["lightteal"])
notes(s, "Opening: three projects, one thread. Paper 1 is where I started when I joined the group. Papers 2 and 3 are the ones I lead, so most of the talk is on those.")

# ================= 2. Roadmap =================
s = new_slide()
title(s, "One thread through three projects")
text(s, "How does spin-dependent band structure reshape correlated magnetism, and can the splitting emerge from the interaction alone?", 0.6, 1.1, 12, 0.5, size=15, color=C["muted"])
cols = [
    ("1", C["p1"], "Spin-dependent NNN hopping t'", "Joined the group, late 2025", "Published 28 Apr 2026", C["green"], C["grey"],
     ["Square lattice, altermagnetic splitting imposed at U = 0", "CPQMC: Neel to spiral to stripe as t' and doping grow", "Joint first author with Y. Li"]),
    ("2", C["p2"], "Combined t_A + t' model", "Submitted 29 Jun 2026", "Major revision, 27 Aug 2026", C["orange"], C["tealSoft"],
     ["NN anisotropy t_A and NNN t' act together", "Four correlation types; Lifshitz condition t_A^2 + 4t'^2 = 1", "Lead author; response in progress"]),
    ("3", C["p3"], "Checkerboard Hubbard model", "Started Jul 2026", "Manuscript for PRB, writing", C["teal"], C["orangeSoft"],
     ["Spin-independent hopping: splitting must come from U", "d_xy altermagnet and d_xy pairing, three lattice sizes", "Lead author; intro and abstract written"]),
]
for i, (n, col, head, when, st, stc, fill, body) in enumerate(cols):
    x, y, w, h = 0.6 + i * 4.15, 1.8, 3.9, 4.9
    card(s, x, y, w, h, fill)
    num_circle(s, n, col, x + 0.25, y + 0.25)
    text(s, head, x + 0.95, y + 0.2, w - 1.15, 0.7, size=16, font=HFONT, bold=True, valign="middle")
    text(s, when, x + 0.25, y + 1.05, w - 0.5, 0.35, size=12, color=C["muted"])
    text(s, st, x + 0.25, y + 1.4, w - 0.5, 0.4, size=14, bold=True, color=stc)
    bullets(s, body, x + 0.25, y + 2.0, w - 0.5, h - 2.2, size=13)
notes(s, "Timeline in one slide. The common question: what does a spin-split band do to the magnetism of a Hubbard model? Papers 1 and 2 impose the splitting through spin-dependent hopping. Paper 3 asks whether the interaction can generate it on its own.")

# ================= 3. Section 1 =================
section(1, C["p1"], "Where I started", "Magnetic correlations with spin-dependent next-nearest-neighbor hopping",
        "Y. Li, R. I. Mushtaq, J. Liu, W. C. Yu, X. Yang, C.-T. Yip, H.-K. Tang, Phys. Rev. B 113, 134443 (2026)",
        "Paper 1. This is the project I joined when I entered the group. I will keep it short because it is published; the point is what it set up for the next two.")

# ================= 4. Paper 1 model =================
s = new_slide(); title(s, "Paper 1: the spin-dependent t' Hubbard model"); badge(s, 1, C["p1"], "NNN paper, published")
img_fit(s, P1("img-001.png"), 0.6, 1.3, 6.6, 4.6)
caption(s, "Fig. 1 of the paper. (a) Sign-alternating NNN hopping, opposite for the two spins. (b) t' = 0: spin degenerate. (c) t' = 0.5: the term 4t' sin kx sin ky splits the bands.", 0.6, 6.0, 6.6, 0.7)
bullets(s, [("Model.", "Square-lattice Hubbard, U = 4t, plus a spin-dependent, diagonally anisotropic NNN hopping. Zero net moment, momentum-dependent splitting."),
            ("Splitting at U = 0.", "The altermagnetic form factor sin kx sin ky is put in by hand. It is d_xy in symmetry."),
            ("Method.", "Constrained-path QMC; the momentum-resolved spin structure factor S^z(k) picks the dominant ordering vector Q."),
            ("Scan.", "Filling n from 0.5 to 1.0 and t' from 0 to 1, lattices L = 8 to 18."),
            ("My part.", "Joined here: learned the CPQMC pipeline, ran the scans and built the momentum-space analysis.")],
        7.6, 1.35, 5.2, 5.3, size=13.5)
notes(s, "Model of paper 1. Emphasise that the splitting exists already at U = 0: the hopping is spin dependent. This is the design choice paper 3 removes. Adjust the 'my part' line to what you want to say about your role.")

# ================= 5. Paper 1 results =================
s = new_slide(); title(s, "Paper 1: Neel, spiral, then stripe as t' and doping grow"); badge(s, 1, C["p1"], "NNN paper, published")
img_fit(s, P1("img-004.png"), 0.6, 1.3, 6.3, 3.7)
caption(s, "Fig. 4: peak coordinate Q of S^z(k) in the (n, t') plane, L = 16. Q = 1 is Neel (pi, pi), Q = 0 is stripe (pi, 0).", 0.6, 5.05, 6.3, 0.5)
img_fit(s, P1("img-003.png"), 7.2, 1.3, 5.6, 3.4)
caption(s, "Fig. 3: density of states. The single van Hove peak at small t' splits as the saddle points move into the zone interior.", 7.2, 4.75, 5.6, 0.5)
bullets(s, [("Half filling, small t':", "(pi, pi) Neel correlations."),
            ("Increase t' or dope:", "the peak moves to a spiral (pi, q), because t' degrades the (pi, pi) nesting."),
            ("Large t', low n:", "short-range regime where stripe (pi, 0) and spiral coexist.")], 0.6, 5.65, 12.2, 1.4, size=13.5)
notes(s, "Main result of paper 1: a continuous shift of the ordering vector, driven by nesting. The DOS panel shows the van Hove mechanism behind it.")

# ================= 6. Paper 1 real space =================
s = new_slide(); title(s, "Paper 1: real-space check and what it taught me"); badge(s, 1, C["p1"], "NNN paper, published")
img_fit(s, P1("img-009.png"), 0.6, 1.3, 5.6, 4.9)
caption(s, "Fig. 9: real-space spin correlation on 16 x 16. (a) Neel checkerboard at n = 0.945, t' = 0.1; (b)-(d) modulated textures on doping or at larger t'.", 0.6, 6.25, 5.6, 0.6)
bullets(s, [("Finite size mattered.", "The phase diagram was checked on L = 8, 12, 14 and 18 before the boundaries were trusted."),
            ("Long range vs short range.", "A peak in S^z(k) is not order. The paper separates the two with a sharpness index and real-space decay."),
            ("Referee round.", "One major revision. The reply added the size series and the U = 0 versus U = 4 comparison."),
            ("Toolkit carried forward.", "S^z(k) classification, the sharpness index R_p, and the Fermi-surface / van Hove analysis are reused in papers 2 and 3.")],
        6.6, 1.35, 6.2, 5.2, size=13.5)
notes(s, "Two lessons from paper 1 that shape everything after: check finite size early, and never call a structure-factor peak 'order' without a scaling argument. Both came back as referee points on paper 2.")

# ================= 7. Section 2 =================
section(2, C["p2"], "Lead author, under major revision", "Competing magnetic correlations with spin-dependent NN and NNN hopping",
        "R. I. Mushtaq, J. Liu, Y. Li, W. C. Yu, X. Yang, H.-K. Tang. Submitted to PRB 29 Jun 2026; major revision requested 27 Aug 2026 (BT15718)",
        "Paper 2, the combined model. This is the first project I lead. Submitted end of June, three referee reports back last week.")

# ================= 8. Paper 2 model =================
s = new_slide(); title(s, "Paper 2: two altermagnetic hoppings in one model"); badge(s, 2, C["p2"], "combined model")
img_fit(s, P2("1_Model.png"), 0.6, 1.25, 12.1, 3.4)
caption(s, "Fig. 1: (a) NN anisotropy t_A (thick vs thin bonds, opposite for the two spins) and NNN hopping t'. (b)-(d) Spin-resolved Fermi surfaces at half filling for t_A only, t' only, and both.", 0.6, 4.65, 12.1, 0.5)
y = 5.3
card(s, 0.6, y, 3.9, 1.6); text(s, [lead("t_A", "spin-dependent NN anisotropy. d_x2-y2 splitting, van Hove points stay at (pi, 0) and (0, pi), (pi, pi) nesting kept.", C["teal"])], 0.8, y + 0.1, 3.5, 1.4, size=12.5)
card(s, 4.7, y, 3.9, 1.6); text(s, [lead("t'", "spin-dependent NNN hopping. d_xy splitting, breaks (pi, pi) nesting, moves the van Hove points (paper 1).", C["teal"])], 4.9, y + 0.1, 3.5, 1.4, size=12.5)
card(s, 8.8, y, 3.9, 1.6, C["orangeSoft"]); text(s, [lead("Together", "U = 4t, L = 8, 10, 14, CPQMC. Band extremum at M turns into a saddle on the locus t_A^2 + 4t'^2 = 1, which couples the two.", C["orange"])], 9.0, y + 0.1, 3.5, 1.4, size=12.5)
notes(s, "The model combines the NN-only limit and the NNN-only limit of paper 1. The quadratic Lifshitz condition is the algebraic reason the two parameters cannot be treated separately; it reduces to t' = 1/2 at t_A = 0 and is never reached for t' = 0.")

# ================= 9. Paper 2 bands =================
s = new_slide(); title(s, "Paper 2: the band structure reorganizes across (t_A, t')"); badge(s, 2, C["p2"], "combined model")
img_fit(s, P2("2_VHS.png"), 0.6, 1.25, 5.6, 5.5)
caption(s, "Fig. 2: spin-up (red) and spin-down (blue) Fermi surfaces at n = 0.5 and 0.9; stars are van Hove points. Rows t_A = 0.1, 0.5, 0.9; columns t' = 0.1, 0.5, 0.9.", 0.6, 6.8, 5.6, 0.5)
img_fit(s, P2("3_DOS.png"), 6.5, 1.25, 6.3, 2.6)
caption(s, "Fig. 3: density of states. t' splits the single van Hove peak; t_A pushes weight away from E = 0 (spin-split van Hove energies at +/- 4 t_A).", 6.5, 3.85, 6.3, 0.5)
bullets(s, [("Small t':", "one logarithmic van Hove peak at the band centre, Fermi surfaces nest at (pi, pi)."),
            ("t' near 0.5:", "Lifshitz-like transition, saddle points migrate into the zone interior."),
            ("Large t_A:", "the two spin Fermi surfaces become quasi one-dimensional, one along x, one along y."),
            ("Why it matters:", "each regime opens a different particle-hole scattering channel, which is what the interacting S^z(k) then selects.")],
        6.5, 4.5, 6.3, 2.5, size=13)
notes(s, "Single-particle picture before the QMC. The referee asked us to draw the Lifshitz line on the phase diagrams; this is the figure it comes from.")

# ================= 10. Paper 2 main =================
s = new_slide(); title(s, "Paper 2: four correlations, selected jointly by n, t_A and t'"); badge(s, 2, C["p2"], "combined model")
img_fit(s, P2("5_Combined_Phase-Diagram.png"), 0.6, 1.2, 8.3, 4.9)
caption(s, "Fig. 5: dominant S^z(k) peak in the (n, t_A) plane at U = 4t, L = 14, for t' = 0 to 0.9. Circles Neel (pi, pi), stars diagonal (q, q), squares stripe (pi, 0), triangles spiral (pi, q).", 0.6, 6.15, 8.3, 0.6)
card(s, 9.2, 1.2, 3.6, 5.5)
bullets(s, [("t_A drives diagonal (q, q).", "Present only at finite t_A; absent at t_A = 0 in every panel."),
            ("t' selects stripe (pi, 0).", "Appears once t' passes about 0.5, through directional nesting."),
            ("n tunes the balance.", "Neel survives near half filling for all t_A; doping hands over to spiral or diagonal."),
            ("Joint effect.", "Diagonal is suppressed by t', stripe by t_A. The boundaries follow the single quadratic locus, not either parameter alone.")],
        9.45, 1.4, 3.15, 5.2, size=12.5)
notes(s, "The headline figure of paper 2. Six panels in t'. Note the top row of each panel: diagonal stars at large t_A vanish once t' is large, and stripes at large t' vanish once t_A is large. That mutual suppression is the non-additive part.")

# ================= 11. Paper 2 sharpness =================
s = new_slide(); title(s, "Paper 2: correlations versus order, and finite size"); badge(s, 2, C["p2"], "combined model")
img_fit(s, P2("6_Sharpness_Factor.png"), 0.6, 1.25, 6.0, 3.6)
caption(s, "Fig. 6: sharpness index R_p of the S^z(k) peak. Red is sharp (long-range tendency), blue is broad. Sharp only near half filling at small t'.", 0.6, 4.85, 6.0, 0.5)
img_fit(s, P2("10_Phase_Diagram_L8.png"), 6.9, 1.25, 5.9, 3.6)
caption(s, "Fig. 10: the same classification on L = 10 (top) and L = 8 (bottom). The four regimes persist; boundaries shift by a grid step.", 6.9, 4.85, 5.9, 0.5)
bullets(s, [("Honest wording.", "Away from half filling and at large t' the peaks are broad: these are correlations, not established order."),
            ("Size check.", "L = 8, 10, 14 agree on which regime dominates where; L = 8 alone has 1171 clean runs."),
            ("Referee 3 asked for more:", "S(Q)/N scaling and the bare susceptibility chi_0(q) next to S^z(k). Both are being added.")],
        0.6, 5.5, 12.2, 1.5, size=13)
notes(s, "Two supporting figures. The sharpness index is the tool from paper 1. The size comparison is what the referees want extended into a proper S(Q)/N scaling.")

# ================= 12. Referees =================
s = new_slide(); title(s, "Paper 2: what the three referees asked for"); badge(s, 2, C["p2"], "major revision, 27 Aug")
refs = [("Referee 1", "constructive", C["green"], C["tealSoft"],
         ["Name 2-3 joint effects that cannot be read off the single-parameter limits", "Discuss CPQMC systematic error; the sawtooth boundaries in Figs. 4 and 5", "U dependence: only U = 4t so far", "Material parameters unknown: weaken the experimental link or estimate them"]),
        ("Referee 2", "negative", C["red"], C["pink"],
         ["Sees the work as paper 1 plus t_A, same method, and does not recommend publication", "Microscopic origin of spin-dependent hopping?", "Trial-wavefunction details and sensitivity", "Rewrite the abstract: motivation and flow"]),
        ("Referee 3", "positive", C["green"], C["tealSoft"],
         ["Compare S^z(q) with the noninteracting chi_0(q)", "Separate correlations from long-range order: S(Q)/N scaling", "Draw the Lifshitz line on the phase diagrams", "Six references on altermagnetism, topology and transport; comment on large t_A, t'"])]
for i, (hd, tone, tc, fill, items) in enumerate(refs):
    x, y, w, h = 0.6 + i * 4.15, 1.25, 3.9, 5.2
    card(s, x, y, w, h, fill)
    text(s, hd, x + 0.25, y + 0.15, 2.2, 0.45, size=18, font=HFONT, bold=True)
    text(s, tone, x + w - 1.6, y + 0.2, 1.35, 0.35, size=12, bold=True, color=tc, align="right")
    bullets(s, items, x + 0.25, y + 0.8, w - 0.5, h - 0.9, size=14)
notes(s, "The reports. Referee 2 is the one that decides the outcome: the objection is novelty relative to paper 1. Everything in the revision is aimed at making the joint, non-additive effects explicit, with numbers, at two lattice sizes.")

# ================= 13. Revision status =================
s = new_slide(); title(s, "Paper 2: revision so far"); badge(s, 2, C["p2"], "major revision, 27 Aug")
card(s, 0.6, 1.25, 6.0, 2.55, C["pink"])
text(s, "Data audit found a real problem", 0.85, 1.35, 5.5, 0.4, size=16, font=HFONT, bold=True, color=C["red"])
bullets(s, ["65 of 800 L = 14 folders were generated with the hopping sign inverted (positive energies). They feed Figs. 5 and 10 and are being rerun.",
            "The 'Neel synergy' table drafted for the reply (up to 89x enhancement) came from those runs. Withdrawn: on the full L = 8 grid the maximum ratio is 1.2x.",
            "L = 8 (1171 runs) and L = 10 (416 runs) are clean."], 0.85, 1.8, 5.5, 1.95, size=12)
card(s, 0.6, 3.95, 6.0, 2.85)
text(s, "Joint effects that hold at both sizes", 0.85, 4.05, 5.5, 0.4, size=16, font=HFONT, bold=True, color=C["teal"])
bullets(s, [("Lifshitz locus separates regimes.", "Stripe fraction 0.6% below vs 11.9% above the line at L = 14 (2.1% vs 13.3% at L = 8); diagonal the reverse."),
            ("Diagonal needs t_A.", "Diagonal fraction rises monotonically with t_A at both sizes, 1% at t_A = 0 to 19% at 0.9 (L = 14)."),
            ("Algebraic non-additivity.", "det H at M = -4(t_A^2 + 4t'^2 - 1), verified symbolically.")], 0.85, 4.5, 5.5, 2.25, size=12)
stat(s, "824", "new CPQMC runs completed, 0 failures", 6.9, 1.35, 2.9)
stat(s, "U = 2, 8", "full grids added to answer the U-dependence point", 9.9, 1.35, 2.9)
stat(s, "48", "trial-wavefunction variants for Referee 2", 6.9, 3.3, 2.9, C["orange"])
stat(s, "7 / 7", "published points reproduced by the new pipeline, peaks within 0.2 to 1.3%", 9.9, 3.3, 2.9, C["orange"])
card(s, 6.9, 5.25, 5.9, 1.55, C["orangeSoft"])
text(s, [lead("Still to do:", "finish the 65 L = 14 reruns, S(Q)/N scaling and chi_0(q) figures, Lifshitz line on the phase diagrams, rewritten abstract, trial-sensitivity appendix, then the point-by-point reply.", C["orange"])], 7.1, 5.35, 5.5, 1.4, size=12.5)
notes(s, "Be upfront about the audit: it was found by us, before resubmission, and the affected cells are being regenerated with the correct convention. The withdrawn synergy claim is replaced by three effects that replicate at both L = 8 and L = 14. The right column is the compute already done for the reply.")

# ================= 14. Section 3 =================
section(3, C["p3"], "Lead author, manuscript in preparation for PRB", "Emergent altermagnetic correlations and d_xy pairing in the checkerboard Hubbard model",
        "Spin-independent hopping. Any spin splitting has to be generated by U. Constrained-path QMC on L = 8, 10, 12 (plus L = 14 to 18 for scaling), benchmarked against exact diagonalization.",
        "Paper 3, the checkerboard. This is the one with the new physics: the splitting is not put in, it emerges. Title is locked; abstract, introduction and section skeleton are written.")

# ================= 15. Paper 3 model =================
s = new_slide(); title(s, "Paper 3: take the spin dependence out of the hopping"); badge(s, 3, C["p3"], "checkerboard")
img_fit(s, P3("1_dnk.png"), 0.6, 1.25, 12.1, 3.3)
caption(s, "Fig. 1: (a) checkerboard lattice, NN hopping t and diagonal hopping t' +/- delta on alternating plaquettes; both spins see the same hopping. (b)-(d) Mean-field spin polarization Delta n(k): zero at delta = 0, a d_xy pattern that grows with delta once a moment M forms.", 0.6, 4.55, 12.1, 0.6)
y = 5.3
card(s, 0.6, y, 3.9, 1.6, C["orangeSoft"]); text(s, [lead("The idea.", "The two magnetic sublattices are related by a 90-degree rotation, not a translation. So a Neel moment M plus the anisotropy delta gives a compensated splitting, m(k) ~ M delta sin kx sin ky.", C["orange"])], 0.8, y + 0.1, 3.5, 1.45, size=12)
card(s, 4.7, y, 3.9, 1.6); text(s, [lead("Versus papers 1 and 2.", "There the splitting exists at U = 0. Here it vanishes if either U or delta is zero. Symmetry fixes the form factor to d_xy, orthogonal to the group's square-lattice d_x2-y2 pairing.", C["teal"])], 4.9, y + 0.1, 3.5, 1.45, size=12)
card(s, 8.8, y, 3.9, 1.6); text(s, [lead("Parameters.", "t = 1, t' = -0.3, delta from 0.1 to 0.7, U from 0 to 5, half filling for the main grid. Observables: Delta_tot = sum_k |n_up(k) - n_dn(k)|, S^z(q), and the four-channel pairing vertex.", C["teal"])], 9.0, y + 0.1, 3.5, 1.45, size=12)
notes(s, "The model in one line: the Hamiltonian is spin independent, so the constrained-path calculation cannot inherit a splitting from the band. Whatever polarisation we measure is interaction driven. Delta_tot is the same order parameter the group's published checkerboard PRB used.")

# ================= 16. Paper 3 context =================
s = new_slide(); title(s, "Paper 3: single-particle context and the doped phase diagram"); badge(s, 3, C["p3"], "checkerboard")
img_fit(s, P3("2_vhs_dos.png"), 0.6, 1.25, 6.4, 2.3)
caption(s, "Fig. 2: Fermi surface at n = 0.8 for delta = 0 and 0.4, and the DOS. delta moves the van Hove points to incommensurate momenta near the Fermi level.", 0.6, 3.6, 6.4, 0.5)
img_fit(s, P3("3_phase_diagram.png"), 7.2, 1.25, 5.6, 3.6)
caption(s, "Fig. 3: dominant S^z(k) peak in the (n, delta) plane, L = 14, U = 4. Neel near half filling; spiral and diagonal on doping.", 7.2, 4.9, 5.6, 0.5)
bullets(s, [("Same toolkit as papers 1 and 2:", "S^z(k) classification and the sharpness index, now on the checkerboard."),
            ("Doped side is short ranged.", "Only the half-filling Neel region carries sharp peaks, so the emergent-altermagnet study is done at n = 1."),
            ("Role in the paper:", "supporting material. The van Hove picture is band context, not the driver of the altermagnetism.")],
        0.6, 4.3, 6.4, 2.5, size=13)
notes(s, "This slide connects to the earlier two papers: the same magnetic classification, now on the checkerboard. The main narrative moves to half filling because that is where the moment lives.")

# ================= 17. Delta_tot =================
s = new_slide(); title(s, "Paper 3: the spin splitting switches on with U"); badge(s, 3, C["p3"], "checkerboard")
img_fit(s, P3("fig03_dtot_emergence.png"), 0.6, 1.25, 7.9, 4.2)
caption(s, "Fig. 3 of the draft: Delta_tot / N over the (U, delta) grid at half filling, L = 12. (a) map; the U = 0 column is exactly zero. (b) cuts at fixed U.", 0.6, 5.5, 7.9, 0.5)
stat(s, "0", "Delta_tot at U = 0, exactly, at every delta (49 cells)", 8.8, 1.3, 4.0)
stat(s, "delta = 0.4", "anisotropy that maximizes the splitting at U = 4.5 to 5", 8.8, 3.0, 4.0, C["orange"], big_size=32)
bullets(s, ["Splitting needs both the interaction and the anisotropy; it is zero if either is switched off.",
            "The optimum runs diagonally: delta as large as the Neel order can tolerate, moving outward as U rises (ridge tracks the order boundary, r = +0.98)."],
        8.8, 4.85, 4.0, 2.0, size=12.5)
notes(s, "First result. The U = 0 column is an exact zero, not a small number: with spin-independent hopping there is nothing to split. The ridge is the interplay: anisotropy is needed for the form factor, but too much anisotropy frustrates the order that carries it.")

# ================= 18. Symmetry + null =================
s = new_slide(); title(s, "Paper 3: the splitting is d_xy, and the delta = 0 control sits at the noise floor"); badge(s, 3, C["p3"], "checkerboard")
img_fit(s, P3("fig04_symmetry.png"), 0.6, 1.25, 7.6, 2.9)
caption(s, "Fig. 4: (a) every interacting cell at L = 8, 10, 12 sits at the d_xy corner (odd under C4, even under the diagonal mirror). (b)-(d) Delta n(k), its C4 sum (zero) and mirror sum (nonzero) at L = 12, U = 4, delta = 0.4.", 0.6, 4.15, 7.6, 0.6)
img_fit(s, P3("fig05_null_control.png"), 8.5, 1.25, 4.3, 3.4)
caption(s, "Fig. 5: delta = 0 (square lattice, splitting forbidden) sets the floor.", 8.5, 4.65, 4.3, 0.4)
stat(s, "126 / 126", "cells classify as d_xy across three sizes", 0.6, 5.0, 3.9)
stat(s, "6.3e-4", "per-site noise floor from the delta = 0 control", 4.7, 5.0, 3.9, C["orange"])
stat(s, "23 to 114x", "signal above the floor in the interacting cells", 8.8, 5.0, 3.9)
notes(s, "Two rigor points the paper leans on. Symmetry: the polarisation changes sign under a 90-degree rotation and is even under the diagonal mirror, which is d_xy and nothing else; this is checked cell by cell. Control: at delta = 0 the model is the square lattice and the splitting is forbidden, so that run measures our statistical floor.")

# ================= 19. Pairing =================
s = new_slide(); title(s, "Paper 3: d_xy pairing leads, and it grows with system size"); badge(s, 3, C["p3"], "checkerboard")
img_fit(s, P3("fig09_pairing_fss.png"), 0.6, 1.25, 8.2, 3.6)
caption(s, "Fig. 9: time-integrated pairing vertex at half filling, U = 4. (a) four channels at L = 12; (b) d_xy versus delta for L = 8, 10, 12; (c) d_xy against 1/N.", 0.6, 4.9, 8.2, 0.5)
bullets(s, [("Zero at U = 0", "in all four channels; the vertex is interaction generated."),
            ("d_xy leads above delta_c = 0.19,", "flat in U (spread 0.03 at L = 12), peak at delta = 0.4."),
            ("Grows with size:", "2.2 to 2.4x from L = 8 to L = 12 at every U. d_x2-y2 is never repulsive."),
            ("Tracks the magnetism:", "with delta and U removed, d_xy correlates with Delta_tot at +0.57 and d_x2-y2 at -0.51.")],
        9.1, 1.3, 3.7, 4.2, size=12.5)
card(s, 0.6, 5.55, 12.2, 1.3, C["orangeSoft"])
text(s, [lead("The claim:", "magnetism and pairing select the same d_xy representation, in a model where neither was put in by hand. On-site s is repulsive (U) and extended s is the large conventional channel, so the statement is about the unconventional, sign-changing channels.", C["orange"])], 0.85, 5.65, 11.7, 1.1, size=13)
notes(s, "The pairing result. Time-integrated vertex, which is the standard measure of a pairing instability and is free of the equal-time discrepancy discussed later. The size growth is the part that makes it more than a finite-cluster curiosity.")

# ================= 20. Crossover maps =================
s = new_slide(); title(s, "Paper 3: the d_xy / d_x2-y2 crossover converges with size"); badge(s, 3, C["p3"], "checkerboard")
img_fit(s, P3("fig_ddiff_vertex_3L.png"), 0.6, 1.25, 12.1, 4.0)
caption(s, "Time-integrated vertex difference d_xy minus d_x2-y2 over (U, delta) at half filling; black line is where the two channels are equal.", 0.6, 5.3, 12.1, 0.4)
stat(s, "0.235 / 0.197 / 0.192", "mean crossover delta_c at L = 8, 10, 12", 0.6, 5.65, 6.0, big_size=32)
stat(s, "0.148 / 0.049 / 0.028", "spread of delta_c across U at L = 8, 10, 12", 6.8, 5.65, 6.0, C["orange"], big_size=32)
notes(s, "Where the leading channel switches. As L grows the line straightens and its U dependence collapses, so the crossover is set by geometry, delta near 0.19, with the interaction driving the amplitude.")

# ================= 21. Finite size + open =================
s = new_slide(); title(s, "Paper 3: finite-size picture of the magnetism, and what is still open"); badge(s, 3, C["p3"], "checkerboard")
img_fit(s, P3("fig_dtot_6delta.png"), 0.6, 1.25, 6.8, 4.9)
caption(s, "Delta_tot / N against U, one curve per size, at six values of delta. Sizes collapse for delta >= 0.4; every size disagreement sits at delta <= 0.3, inside the ordered wedge.", 0.6, 6.2, 6.8, 0.6)
card(s, 7.7, 1.25, 5.1, 2.55)
text(s, "Neel order at half filling", 7.95, 1.35, 4.7, 0.4, size=15, font=HFONT, bold=True, color=C["teal"])
bullets(s, ["S(pi,pi)/N is flat across L = 12, 14, 16, 18 at delta = 0.2, U = 4 (0.074, 0.074, 0.072, 0.071), giving m = 0.27.",
            "Order occupies a wedge that widens with U: none at U <= 2, delta up to 0.1 at U = 3, 0.2 at U = 4, 0.4 at U = 5.",
            "U = 0 control decays as 1/N at every delta, as it must."], 7.95, 1.8, 4.7, 1.95, size=11.5)
card(s, 7.7, 3.95, 5.1, 2.85, C["pink"])
text(s, "Still open before Section III", 7.95, 4.05, 4.7, 0.4, size=15, font=HFONT, bold=True, color=C["red"])
bullets(s, [("Trial dependence.", "An audit of 304 half-filling runs found 28 where the trial used had a higher energy and a different order than the best seed. Those cells (a diagonal in the U-delta plane) are being rerun."),
            ("Equal-time vs time-integrated.", "Fortran and Python equal-time vertices disagree at U > 2; suspect the Trotter step (0.05 vs 0.01). One Delta-tau = 0.01 column at U = 5 will settle it."),
            ("Delta_tot at small U", "tracks the trial moment, not U. Quote only the U >= 3.5 plateau for now.")], 7.95, 4.5, 4.7, 2.3, size=11, head_color=C["red"])
notes(s, "The honest slide. Neel order at half filling is established by the structure-factor scaling to L = 18. The polarisation Delta_tot is a different quantity and is trial-sensitive at small U, which the audit has now traced to a wrong trial selection in 28 runs. Nothing in the headline claims uses those cells, but they are being rerun before the paper goes out.")

# ================= 22. Manuscript status =================
s = new_slide(); title(s, "Paper 3: where the manuscript stands"); badge(s, 3, C["p3"], "checkerboard")
rows = [("Title, abstract (no numbers, PRB style)", "done"), ("Introduction, five paragraphs plus roadmap", "done"),
        ("Section skeleton, 21 of ~70 references collected", "done"), ("Figs. 3, 4, 5, 9 and the two-row equal-time / vertex figure", "built"),
        ("Fig. 6, h -> 0 and finite size (L = 16 data landed)", "to build"), ("Section II, model and method", "next"),
        ("Sections III to VII", "not written"), ("Delta-tau test and the 28 trial reruns", "queued")]
colr = {"done": C["green"], "built": C["green"], "to build": C["orange"], "next": C["orange"], "not written": C["muted"], "queued": C["orange"]}
for i, (a, b) in enumerate(rows):
    y = 1.35 + i * 0.6
    shape(s, MSO_SHAPE.ROUNDED_RECTANGLE, 0.6, y, 7.4, 0.5, "F4F7F7" if i % 2 else C["white"], radius=0.1, line="E1E8E7")
    text(s, a, 0.8, y, 5.3, 0.5, size=13, valign="middle")
    text(s, b, 6.1, y, 1.8, 0.5, size=12.5, bold=True, color=colr[b], align="right", valign="middle")
card(s, 8.4, 1.35, 4.4, 4.75, C["orangeSoft"])
text(s, "Rules the draft follows", 8.65, 1.45, 4.0, 0.4, size=15, font=HFONT, bold=True, color=C["orange"])
bullets(s, ["Only L = 8, 10, 12 (and larger) are quoted. L = 6 once gave a 5-sigma wrong sign and is never used.",
            "No slope quoted at delta = 0.3, the crossover where periodic and antiperiodic boundaries disagree.",
            "Every figure script is in the repo and rerun from the CSVs; parameter maps use bilinear interpolation on the measured grid, never Delaunay.",
            "Symmetry ratios and the delta = 0 null control appear in every claim about the splitting."], 8.65, 1.95, 4.0, 4.1, size=11.5)
notes(s, "Status of the writing. Introduction and abstract are settled; Section II is the next writing task. The queued cluster jobs are small: one column of the vertex grid at the finer Trotter step, and the reruns from the trial audit.")

# ================= 23. Summary =================
s = new_slide(C["ink"]); title(s, "Summary and the next four weeks", dark=True)
cols = [("1", C["p1"], "NNN paper", ["Published, PRB 113, 134443.", "Its toolkit (S^z(k) classification, sharpness index, van Hove analysis) carries the other two."]),
        ("2", C["p2"], "Combined model", ["Finish the 65 L = 14 reruns.", "S(Q)/N scaling, chi_0(q), Lifshitz line on the diagrams.", "Trial-sensitivity appendix, rewritten abstract.", "Point-by-point reply, resubmit."]),
        ("3", C["p3"], "Checkerboard", ["Delta-tau test and the 28 trial reruns.", "Build Fig. 6 (finite size).", "Write Section II, then III.", "Target: full draft to the group for comments."])]
for i, (n, col, hd, lines) in enumerate(cols):
    x, y, w, h = 0.6 + i * 4.15, 1.4, 3.9, 5.1
    card(s, x, y, w, h, C["darkcard"])
    num_circle(s, n, col, x + 0.25, y + 0.25)
    text(s, hd, x + 0.95, y + 0.25, w - 1.15, 0.55, size=18, font=HFONT, bold=True, color=C["white"], valign="middle")
    bullets(s, lines, x + 0.25, y + 1.05, w - 0.5, h - 1.25, size=13, color=C["palegreen"])
text(s, "One thread: papers 1 and 2 impose the spin splitting through the hopping; paper 3 shows it can emerge from U and lattice geometry alone, with d_xy symmetry shared by magnetism and pairing.", 0.6, 6.65, 12.1, 0.6, size=13, italic=True, color=C["lightteal"])
notes(s, "Close on the thread and the concrete next steps. Ask for input on two things: the strongest way to answer Referee 2's novelty objection, and whether the equal-time panel should stay in paper 3.")

prs.save(OUT)
print("wrote", OUT, "slides:", len(prs.slides))
