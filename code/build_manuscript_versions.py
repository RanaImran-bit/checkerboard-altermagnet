#!/usr/bin/env python -u
"""Build the three manuscript versions from one shared body of text.

  v1  pairing AND magnetism   the combined story, the fullest paper
  v2  pairing only            the pairing vertex as the result
  v3  magnetism only          the emergent altermagnetic correlations

Written as a builder rather than three hand-maintained files so that a number or
a sentence corrected once propagates to every version that uses it. The section
text comes from manuscript/checkerboard_PRB.tex, which is split on its own
section headers; only the abstract, the closing paragraph of the introduction,
the figure set and the conclusion differ between versions.

House style enforced here: no em-dashes, no semicolons in prose, d_xy and
d_x2-y2 rather than the group-theory labels, constrained-path results never
called exact.

Figures are referenced WITHOUT a file extension so that graphicx takes the PDF
where one exists and the PNG otherwise. Eight of the figures are currently PDF
and the rest are PNG.
"""
import os, re, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(REPO, "manuscript", "checkerboard_PRB.tex")
OUT  = os.path.join(REPO, "manuscript")

# --------------------------------------------------------------------------
# figure captions. Each names what is plotted, the lattice size, and the point
# the figure is in the paper to make.
# --------------------------------------------------------------------------
FIG = {
"fig03_dtot_emergence": (
 r"Momentum-space spin polarization $\Delta_{\mathrm{tot}}/N$ over the "
 r"$(U,\delta)$ plane at half filling, $L=12$. The $U=0$ column is identically "
 r"zero at every anisotropy, which is required because the hopping in "
 r"Eq.~(\ref{eq:H}) is spin independent. The signal appears only once the "
 r"interaction is turned on, which is the emergent route rather than the "
 r"built-in one."),
"fig04_symmetry": (
 r"Symmetry classification of $\Delta n(\bm{k})$ for all 126 interacting "
 r"parameter sets at $L=8$, $10$ and $12$. The ratios $A_{\mathrm{odd}}$ of "
 r"Eq.~(\ref{eq:aodd}) and $M_{\mathrm{odd}}$ of Eq.~(\ref{eq:modd}) place every "
 r"cell at $A_{\mathrm{odd}}\simeq0$ and $M_{\mathrm{odd}}\simeq2$, which is odd "
 r"under $C_4$ and even under the diagonal mirror and identifies the $d_{xy}$ "
 r"channel. No cell is classified as $d_{x^2-y^2}$."),
"fig05_null_control": (
 r"Null control at vanishing anisotropy, logarithmic axis. At $\delta=0$ the two "
 r"sublattices are related by a translation and $\Delta_{\mathrm{tot}}$ must "
 r"vanish by symmetry at any $U$, so what is measured there is the statistical "
 r"floor of the estimator at production settings. The floor is "
 r"$6.3(2)\times10^{-4}$ per site and the measured signal at $\delta=0.3$ sits "
 r"between 23 and 114 times above it."),
"fig_dtot_6delta": (
 r"$\Delta_{\mathrm{tot}}/N$ against $U$ at six anisotropies from $\delta=0.1$ "
 r"to $0.6$, for $L=8$, $10$, $12$ and $14$. The background shades the value of "
 r"$\Delta_{\mathrm{tot}}/N$ on the same scale as the heat maps so that the two "
 r"presentations can be read against each other. $L=14$ forms a curve only for "
 r"$\delta\leq0.3$."),
"fig_dtot_maps_with_cuts": (
 r"$\Delta_{\mathrm{tot}}/N$ over the $(U,\delta)$ plane with the cuts used in "
 r"the preceding figure drawn on the map, $L=12$."),
"fig_smallU_L8to14_bg": (
 r"Behavior at small $U$ for $\delta=0.1$, $0.2$ and $0.3$ across $L=8$ to "
 r"$L=14$. The feature near $U=2$ does not survive to $L=14$ at $\delta=0.1$, "
 r"where $\Delta_{\mathrm{tot}}/N$ falls to about twice the $\delta=0$ floor and "
 r"the correlation with $\sin k_x \sin k_y$ is lost. At $\delta=0.3$ the decay "
 r"with size is steady. The $\delta=0.2$ column persists at every size."),
"fig06_fss_hscan": (
 r"(a) $\Delta_{\mathrm{tot}}/N$ against the symmetry-breaking field $h$ for "
 r"$L=8$ to $16$ at $U=4$ and $\delta=0.3$, with the straight-line fits over the "
 r"three smallest fields that define the $h\to0$ intercept. (b) those intercepts "
 r"against $1/L$. The two curves are fits to the $L\leq14$ points alone, linear "
 r"in $1/L$ and in $1/L^2$, extended to $L=16$. They predict $0.01325$ and "
 r"$0.01457$ there, and the measured value is $0.00733 \pm 0.00106$. The drop "
 r"from $L=14$ is a factor of $0.43$ where a $1/L$ law gives $0.88$, so the "
 r"nonzero intercept seen at each individual size is a finite-size effect."),
"fig_overlay_4L": (
 r"The four lattice sizes overlaid, so that the size dependence of "
 r"$\Delta_{\mathrm{tot}}/N$ can be read at fixed $U$ and $\delta$."),
"fig09_pairing_fss": (
 r"Time-integrated pairing vertex at half filling. (a) the four channels at "
 r"$L=12$ and $U=4$, (b) the $d_{xy}$ vertex against $\delta$ for $L=8$, $10$ "
 r"and $12$, and (c) the same quantity against $1/N$. The vertex is the "
 r"interaction-induced part of the pair-field correlator and vanishes "
 r"identically at $U=0$ in all four channels. $d_{xy}$ grows with system size at "
 r"every interaction strength studied."),
"fig_2row_eqtime_vertex_3L": (
 r"Equal-time pair correlator (top row) and time-integrated vertex (bottom row) "
 r"for $L=8$, $10$ and $12$. The channel crossover straightens and its spread "
 r"across $U$ narrows along the vertex row, which is the behavior expected of a "
 r"quantity fixed by the lattice geometry rather than by the coupling."),
"fig_ddiff_vertex_3L": (
 r"Difference between the $d_{xy}$ and $d_{x^2-y^2}$ vertices over the "
 r"$(U,\delta)$ plane at $L=8$, $10$ and $12$. The zero contour is the crossover "
 r"anisotropy $\delta_c$, which is nearly flat in $U$ and whose spread across "
 r"the interaction strengths studied falls from $0.148$ to $0.049$ to $0.028$ as "
 r"the lattice grows."),
"fig_dsum_L12": (
 r"Sum of the $d_{xy}$ and $d_{x^2-y^2}$ vertices at $L=12$, showing the total "
 r"weight carried by the two sign-changing channels."),
"fig_ddiff_L12": (
 r"Difference of the $d_{xy}$ and $d_{x^2-y^2}$ vertices at $L=12$, on the same "
 r"colour scale as the preceding figure."),
"fig_dchannels_eqtime_vs_vertex": (
 r"The two $d$ channels measured equal time (left) and time integrated (right) "
 r"at $L=12$. The equal-time correlator carries a large uncorrelated part, so "
 r"the vertex is the measure used for the claims in the text."),
"fig_d_channels_together": (
 r"Both $d$ channels in the $(U,\delta)$ plane at $L=12$. $d_{x^2-y^2}$ is "
 r"positive everywhere, so it is never repulsive even where it is subdominant, "
 r"and $d_{xy}$ is the larger of the two above $\delta\simeq0.19$."),
"fig_dchannels_over_dtot": (
 r"Contours of the two $d$ vertices drawn over the measured "
 r"$\Delta_{\mathrm{tot}}/N$. The pairing channel that leads changes across the "
 r"same region of the plane in which the magnetic signal changes character, so "
 r"two independently measured quantities turn over together."),
"fig_pairing_on_dtot_3L": (
 r"Pairing contours on the $\Delta_{\mathrm{tot}}/N$ background for $L=8$, $10$ "
 r"and $12$, showing that the relation between the two survives the change of "
 r"lattice size."),
}

def figure(key, width=r"\columnwidth", star=False, pos=None):
    # figure* is restricted to page tops in two-column revtex, so give the wide
    # ones a float page as well or they queue up and LaTeX reports them stuck
    if pos is None:
        pos = "tp" if star else "htbp"
    env = "figure*" if star else "figure"
    return ("\n\\begin{%s}[%s]\n\\includegraphics[width=%s]{figures/%s}\n"
            "\\caption{%s}\n\\label{fig:%s}\n\\end{%s}\n"
            % (env, pos, width, key, FIG[key], key, env))

# --------------------------------------------------------------------------
def load_blocks():
    s = open(SRC).read()
    pre, rest = s.split("% ---------------------------------------------------------------------------\n\\section{Introduction}", 1)
    body = "\\section{Introduction}" + rest
    body = body.split("\\begin{acknowledgments}")[0]
    tail = s[s.index("\\begin{acknowledgments}"):]
    parts = re.split(r"(?m)^% -{20,}\n\\section\{([^}]*)\}\n", "% " + "-"*75 + "\n" + body)
    blocks, i = {}, 1
    while i < len(parts):
        name = re.sub(r"\s+", " ", parts[i]).strip()
        blocks[name] = parts[i+1]
        i += 2
    return pre, blocks, tail

# --------------------------------------------------------------------------
# version-specific front matter and closing text
# --------------------------------------------------------------------------
ABSTRACT = {}
ABSTRACT["v1"] = r"""
We study the half-filled Hubbard model on the checkerboard lattice by
constrained-path quantum Monte Carlo. The Hamiltonian is spin independent, so
the two magnetic sublattices are related by a fourfold rotation rather than by a
translation, and any spin splitting of the momentum distribution must be
generated by the interaction. We find compensated spin-split correlations above
a finite interaction strength whose momentum-space polarization is odd under
$C_4$ and even under the diagonal mirror at all 126 interacting parameter sets
studied, which identifies the $d_{xy}$ channel. A control calculation at
vanishing lattice anisotropy, where the model reduces to the square lattice and
the splitting is forbidden by symmetry, fixes the statistical floor of the order
parameter at $6.3(2)\times10^{-4}$ per site, two orders of magnitude below the
measured signal. The interaction-induced pairing vertex vanishes identically in
the noninteracting limit and is largest in the $d_{xy}$ channel above an
anisotropy $\delta_c \simeq 0.19$ whose value converges with system size.
Magnetism and pairing therefore select the same representation, which
distinguishes this state from the $d_{x^2-y^2}$ pairing reported for the
checkerboard lattice in the regime where the two models coincide. Finite-size
analysis up to $16\times16$ sites places the extrapolated staggered moment below
the resolution of the calculation, so the correlations remain short ranged on
the lattices studied.
"""
ABSTRACT["v2"] = r"""
We compute the pair-field vertex of the half-filled Hubbard model on the
checkerboard lattice by constrained-path quantum Monte Carlo on lattices of
$8\times8$, $10\times10$ and $12\times12$ sites. Separating the
interaction-induced vertex from the uncorrelated part of the correlator, we find
that the vertex vanishes identically in all four channels at $U=0$, so every
feature reported here is generated by the interaction. The on-site $s$ vertex is
negative throughout, as a repulsive interaction requires. The leading
sign-changing channel changes with the lattice anisotropy $\delta$, from
extended $s$ at small $\delta$ to $d_{xy}$ above a crossover $\delta_c$. That
crossover is nearly independent of the interaction strength and its spread
across $U$ narrows from $0.148$ to $0.028$ between $L=8$ and $L=12$, converging
on $\delta_c \simeq 0.19$, which identifies it as a property of the lattice
geometry. The $d_{xy}$ vertex grows with system size by a factor between $2.18$
and $2.38$ over the same range, and $d_{x^2-y^2}$ is positive in all 126 cells
examined. The dominant channel lies outside the pairing basis used in earlier
work on this lattice.
"""
ABSTRACT["v3"] = r"""
We study magnetic correlations in the half-filled Hubbard model on the
checkerboard lattice by constrained-path quantum Monte Carlo. The hopping is
spin independent, so the momentum distributions of the two spin species coincide
at $U=0$ and any splitting must be generated by the interaction. The two
magnetic sublattices are exchanged by a fourfold rotation and by no translation,
which is the symmetry condition for a compensated spin splitting. We find such a
splitting above a finite interaction strength, with a momentum-space
polarization that is odd under $C_4$ and even under the diagonal mirror at all
126 interacting parameter sets studied, identifying the $d_{xy}$ channel. A
control at vanishing anisotropy, where the model reduces to the square lattice
and the splitting is forbidden, fixes the statistical floor of the estimator at
$6.3(2)\times10^{-4}$ per site, between 23 and 114 times below the measured
signal. Finite-size analysis on lattices from $8\times8$ to $16\times16$ sites
shows the equal-time structure factor falling as $1/N$ and the extrapolated
staggered moment dropping below the resolution of the calculation, so the
correlations remain short ranged on the lattices reached here.
"""

TITLE = {
 "v1": r"Emergent altermagnetic correlations and $d_{xy}$ pairing in the" "\n" r"checkerboard Hubbard model",
 "v2": r"Interaction-induced $d_{xy}$ pairing in the checkerboard Hubbard model",
 "v3": r"Emergent altermagnetic correlations in the checkerboard Hubbard model",
}

ROADMAP = {}
ROADMAP["v1"] = r"""
We report constrained-path quantum Monte Carlo (CPQMC) calculations for this
model at half filling on lattices of $8\times8$ up to $16\times16$ sites. Three
results follow. The momentum-resolved spin polarization is nonzero above a
finite $U$ and has $d_{xy}$ symmetry at every interacting parameter set we
examined. A null calculation at vanishing anisotropy, where symmetry forbids the
signal, places the statistical floor two orders of magnitude below the measured
values. The interaction-induced pairing vertex is largest in the $d_{xy}$
channel over the same anisotropy range, so the magnetic and pairing
instabilities lie in the same representation.
"""
ROADMAP["v2"] = r"""
We report constrained-path quantum Monte Carlo (CPQMC) calculations of the
pair-field vertex for this model at half filling on lattices of $8\times8$,
$10\times10$ and $12\times12$ sites. The vertex is the part of the pair
correlator that the interaction builds, and it is the quantity that distinguishes
a pairing tendency from the large positive correlator that free electrons already
carry. Three results follow. The vertex is zero in every channel at $U=0$. The
leading sign-changing channel is $d_{xy}$ above an anisotropy $\delta_c$ whose
value converges with system size. The $d_{xy}$ vertex grows with system size at
every interaction strength studied, while $d_{x^2-y^2}$ remains positive
throughout and is never repulsive.
"""
ROADMAP["v3"] = r"""
We report constrained-path quantum Monte Carlo (CPQMC) calculations for this
model at half filling on lattices of $8\times8$ up to $16\times16$ sites. Three
results follow. The momentum-resolved spin polarization is nonzero above a
finite $U$ and has $d_{xy}$ symmetry at every interacting parameter set we
examined. A null calculation at vanishing anisotropy, where symmetry forbids the
signal, places the statistical floor two orders of magnitude below the measured
values. Finite-size analysis to $16\times16$ sites shows that the correlations
stay short ranged on the lattices reached, so we describe them as correlations
rather than as order.
"""

CONCLUSION = {}
CONCLUSION["v1"] = r"""
The half-filled checkerboard Hubbard model supports compensated spin-split
correlations that are absent at $U=0$ and are $d_{xy}$ in symmetry at every
parameter set examined. A control at vanishing anisotropy, where the splitting is
forbidden, places the statistical floor two orders of magnitude below the signal.
The interaction-induced pairing vertex is largest in the same $d_{xy}$ channel,
and the anisotropy at which the pairing channel changes coincides with the
anisotropy at which the magnetic response changes sign.

The mechanism uses no spin-dependent hopping and no spin-orbit coupling. It
requires only a lattice whose magnetic sublattices are related by a rotation and
a uniform Hubbard repulsion, which suggests that the same construction should
apply to other lattices with the same sublattice relation. The correlations
remain short ranged on the lattices reached here, and the constrained-path bias
remains.
"""
CONCLUSION["v2"] = r"""
The interaction-induced pairing vertex of the half-filled checkerboard Hubbard
model is zero in every channel at $U=0$ and is largest in the $d_{xy}$ channel
above an anisotropy $\delta_c \simeq 0.19$. The crossover is nearly independent
of the interaction strength and converges with system size, which identifies it
with the lattice geometry rather than with the coupling. The $d_{xy}$ vertex
grows with system size at every interaction strength studied, and $d_{x^2-y^2}$
is positive in all 126 cells examined and therefore never repulsive.

The dominant channel lies outside the basis used in the closest published
calculation on this lattice, so the comparison is one of basis rather than of
disagreement. Whether the $d_{xy}$ tendency survives to the thermodynamic limit
is not settled by the sizes reached here, and the constrained-path bias remains.
"""
CONCLUSION["v3"] = r"""
The half-filled checkerboard Hubbard model supports compensated spin-split
correlations that are absent at $U=0$ and are $d_{xy}$ in symmetry at every
parameter set examined. A control at vanishing anisotropy, where the splitting is
forbidden by symmetry, places the statistical floor two orders of magnitude below
the signal, which establishes that the measured polarization is not an artifact
of the estimator or of the symmetry-broken trial state.

The mechanism uses no spin-dependent hopping and no spin-orbit coupling. It
requires only a lattice whose magnetic sublattices are related by a rotation and
a uniform Hubbard repulsion. The correlations remain short ranged on the
lattices reached here, so we do not claim long-range order, and the
constrained-path bias remains.
"""

# --------------------------------------------------------------------------
# which sections and which figures each version carries
# --------------------------------------------------------------------------
MAG_SECS  = ["Emergent altermagnetic correlations",
             "Null control at vanishing anisotropy",
             "Finite-size behavior of the magnetic correlations"]
PAIR_SECS = ["Pairing in the same representation",
             "Size dependence of the pairing vertex"]

PLAN = {
 "v1": dict(secs = MAG_SECS + PAIR_SECS + ["Robustness and comparison"],
            figs = {"Emergent altermagnetic correlations":
                        [("fig03_dtot_emergence", False), ("fig04_symmetry", True),
                         ("fig_dtot_6delta", True), ("fig_dtot_maps_with_cuts", False)],
                    "Null control at vanishing anisotropy": [("fig05_null_control", False)],
                    "Finite-size behavior of the magnetic correlations":
                        [("fig06_fss_hscan", True), ("fig_smallU_L8to14_bg", True),
                         ("fig_overlay_4L", False)],
                    "Pairing in the same representation":
                        [("fig09_pairing_fss", True), ("fig_dchannels_eqtime_vs_vertex", True),
                         ("fig_d_channels_together", False)],
                    "Size dependence of the pairing vertex":
                        [("fig_2row_eqtime_vertex_3L", True), ("fig_ddiff_vertex_3L", True),
                         ("fig_dchannels_over_dtot", True), ("fig_pairing_on_dtot_3L", True)]},
            apps = ["Trial wave function and the selection of the constrained path",
                    "Convergence in projection length", "Boundary conditions"]),
 "v2": dict(secs = PAIR_SECS + ["Robustness and comparison"],
            figs = {"Pairing in the same representation":
                        [("fig09_pairing_fss", True), ("fig_dchannels_eqtime_vs_vertex", True),
                         ("fig_d_channels_together", False), ("fig_dsum_L12", False),
                         ("fig_ddiff_L12", False)],
                    "Size dependence of the pairing vertex":
                        [("fig_2row_eqtime_vertex_3L", True), ("fig_ddiff_vertex_3L", True)]},
            apps = ["Trial wave function and the selection of the constrained path",
                    "Convergence in projection length", "Boundary conditions"]),
 "v3": dict(secs = MAG_SECS + ["Robustness and comparison"],
            figs = {"Emergent altermagnetic correlations":
                        [("fig03_dtot_emergence", False), ("fig04_symmetry", True),
                         ("fig_dtot_6delta", True), ("fig_dtot_maps_with_cuts", False),
                         ("fig_overlay_4L", False)],
                    "Null control at vanishing anisotropy": [("fig05_null_control", False)],
                    "Finite-size behavior of the magnetic correlations":
                        [("fig06_fss_hscan", True), ("fig_smallU_L8to14_bg", True)]},
            apps = ["Trial wave function and the selection of the constrained path",
                    "Convergence in projection length", "Boundary conditions"]),
}

# sentences that only make sense when the other half of the story is present
CROSSLINK = [
 (r"""The magnetic and pairing instabilities therefore select the same
representation. The polarization is $d_{xy}$ at all 126 parameter sets, and the
pairing vertex is largest in $d_{xy}$ over the anisotropy range where that
classification holds most strongly.""",
  r"""The $d_{xy}$ channel therefore leads over the whole range of anisotropy in
which the vertex is appreciable, and it does so at every interaction strength
studied."""),
 (r"""Magnetism and pairing move together rather than merely coexisting. With $\delta$
and $U$ both partialled out, the correlation between the pairing vertex and
$\Delta_{\mathrm{tot}}$ is $+0.57$ in the $d_{xy}$ channel and $-0.51$ in the
$d_{x^2-y^2}$ channel. The two channels respond to the magnetic signal with
opposite sign at fixed anisotropy and fixed interaction.""",
  r"""The two $d$ channels respond with opposite sign to the same change in
anisotropy at fixed interaction, so they are not simply rising and falling
together with the overall pairing scale."""),
]

# --------------------------------------------------------------------------
# Version-specific trims. A pairing-only paper that defines an order parameter
# it never measures invites the obvious referee question, and the same holds for
# a magnetism-only paper that carries a pairing comparison.
# --------------------------------------------------------------------------
PAIR_METHOD = r"""
The pair-field correlator is measured in four channels, on-site $s$, extended
$s$, $d_{x^2-y^2}$ and $d_{xy}$, each defined by a form factor
$f_{\bm{\delta}}$ on the bonds it connects. The extended $s$ and $d_{x^2-y^2}$
channels use the four nearest-neighbor bonds with $f$ uniform and with
$f = \pm 1$ on the $x$ and $y$ bonds respectively, and the $d_{xy}$ channel uses
the four diagonal bonds with $f = \pm 1$ on the two diagonal orientations. We
report the interaction-induced vertex, which is the correlator with its
uncorrelated Wick contribution subtracted. The subtraction matters because the
full correlator is large and positive already for free electrons, so only the
vertex distinguishes a pairing tendency built by the interaction from the
kinematic background.
"""

V2_BC = r"""
The sensitivity of the calculation to the momentum mesh was characterized on the
magnetic response rather than on the vertex. Its slope in $U$ agrees between the
two boundary conditions at $\delta = 0.1$ and at $\delta = 0.5$ and is ambiguous
at $\delta = 0.3$, where that slope passes through zero and its sign is most
sensitive to the mesh. We have not repeated the vertex calculation with
antiperiodic boundaries, so the mesh dependence of the pairing results is not
quantified here.
"""

TRIM = {
 "v2": {"Model and method": [
     # the order parameter, the field that induces it, and the symmetry ratios
     # all belong to the magnetic half of the story
     ("The order parameter is the momentum-space spin polarization", "unlike a projection onto a single lattice harmonic they are not sensitive to\nhow the weight is distributed among higher harmonics of the same symmetry.\n", PAIR_METHOD)],
    "Robustness and comparison": [
     ("The sign of $\\partial \\Delta_{\\mathrm{tot}} / \\partial U$ agrees", "% data/bc_test, analyze_bc_test.py\n", "")],
    "Boundary conditions": [
     ("The sign of $\\partial \\Delta_{\\mathrm{tot}} / \\partial U$ agrees", "a slope there.\n", V2_BC)]},
 "v3": {"Robustness and comparison": [
     ("The closest published work is a determinant quantum Monte Carlo study", "not that the two calculations disagree.\n", "")]},
}

def apply_trim(tag, name, body):
    for start, end, repl in TRIM.get(tag, {}).get(name, []):
        i = body.find(start)
        if i < 0:
            print("  WARN trim start not found in %s/%s" % (tag, name)); continue
        j = body.find(end, i)
        if j < 0:
            print("  WARN trim end not found in %s/%s" % (tag, name)); continue
        body = body[:i] + repl.strip() + ("\n" if repl.strip() else "") + body[j+len(end):]
    return body

def build(tag):
    pre, blocks, tail = load_blocks()
    plan = PLAN[tag]
    # front matter
    # lambda replacements: the LaTeX bodies contain backslash sequences that
    # re.sub would otherwise try to interpret as group references
    # brace-matched title replacement. A non-greedy .*?\} stops at the first
    # closing brace, which in this title is the one inside $d_{xy}$, so it
    # truncates the command and leaves the tail behind as stray text.
    i = pre.index("\\title{")
    j, depth = i + len("\\title{"), 1
    while depth:
        if pre[j] == "{": depth += 1
        elif pre[j] == "}": depth -= 1
        j += 1
    pre = pre[:i] + "\\title{%s}" % TITLE[tag] + pre[j:]
    pre = re.sub(r"\\begin\{abstract\}.*?\\end\{abstract\}",
                 lambda m: "\\begin{abstract}%s\\end{abstract}" % ABSTRACT[tag],
                 pre, flags=re.S)
    # loosen the float parameters: the defaults refuse a page that is mostly
    # figures, which is exactly what a 15-figure article needs
    pre = pre.replace("\\begin{document}",
        "\\renewcommand{\\topfraction}{0.9}\n"
        "\\renewcommand{\\bottomfraction}{0.8}\n"
        "\\renewcommand{\\textfraction}{0.07}\n"
        "\\renewcommand{\\floatpagefraction}{0.7}\n"
        "\\setcounter{topnumber}{3}\n"
        "\\setcounter{totalnumber}{4}\n\n\\begin{document}")
    pre = pre.replace("%  DRAFT v1 -- Physical Review B",
                      "%%  Physical Review B -- version %s" % tag)

    out = [pre]
    intro = blocks["Introduction"]
    old_road = re.search(r"We report constrained-path quantum Monte Carlo.*?representation\.\n",
                         intro, re.S)
    if old_road:
        intro = intro.replace(old_road.group(0), ROADMAP[tag].strip() + "\n")
    out.append(sec("Introduction", intro))
    out.append(sec("Model and method", apply_trim(tag, "Model and method",
                                                  blocks["Model and method"])))

    for name in plan["secs"]:
        body = apply_trim(tag, name, blocks[name])
        if tag != "v1":
            for combined, alone in CROSSLINK:
                body = body.replace(combined, alone)
        for key, star in plan["figs"].get(name, []):
            body += figure(key, r"\textwidth" if star else r"\columnwidth", star)
        out.append(sec(name, body))

    out.append(sec("Conclusion", CONCLUSION[tag].strip() + "\n"))
    out.append("\n\\appendix\n")
    for name in plan["apps"]:
        out.append(sec(name, apply_trim(tag, name, blocks[name])))
    txt = "".join(out) + tail if False else "".join(out)
    t = tail
    if tag == "v3":
        i = t.find("\\bibitem{Pan2024}")
        if i >= 0:
            j = t.find("\\bibitem", i + 10)
            t = t[:i] + (t[j:] if j > 0 else t[t.find("\\end{thebibliography}", i):])
    out.append(t)
    return "".join(out)

def sec(name, body):
    return ("\n%% %s\n\\section{%s}\n%s" % ("-"*73, name, body))

if __name__ == "__main__":
    for tag in ("v1", "v2", "v3"):
        txt = build(tag)
        path = os.path.join(OUT, "%s_%s.tex" % (tag, {"v1":"pairing_and_magnetism",
                                                      "v2":"pairing", "v3":"magnetism"}[tag]))
        open(path, "w").write(txt)
        nfig = txt.count("\\includegraphics")
        nsec = len(re.findall(r"(?m)^\\section\{", txt))
        print("%s  %-32s %5d lines  %2d sections  %2d figures"
              % (tag, os.path.basename(path), len(txt.splitlines()), nsec, nfig))
