# Title, abstract, and figure captions (Overleaf-ready)

## Title
Interaction-driven altermagnetism in the checkerboard Hubbard model

## Abstract
```latex
\begin{abstract}
Altermagnets combine antiferromagnetic compensation with a momentum-dependent, spin-split band
structure of definite $d$-wave symmetry, yet most microscopic studies impose this splitting through
spin-dependent hopping or treat it within mean-field theory. We study altermagnetism that emerges from
electronic correlations in a spin-independent single-band Hubbard model on the anisotropic checkerboard
(planar-pyrochlore) lattice, using constrained-path and determinant quantum Monte Carlo benchmarked
against exact diagonalization. The geometric diagonal anisotropy $\delta$ turns the interaction-driven
N\'eel order into an altermagnet whose spin splitting has $d_{xy}$ symmetry and scales as the product
$M\delta$ of the ordered moment $M$ and the anisotropy, vanishing when either is absent. In the
filling--anisotropy plane, long-range N\'eel order is confined to the vicinity of half-filling, while
doping yields short-range incommensurate correlations. Because the altermagnetic form factor is fixed by
the lattice geometry, it is $d_{xy}$ here rather than the $d_{x^2-y^2}$ of square-lattice models, which
implies a distinct symmetry channel for interaction-driven pairing. These results identify the
checkerboard Hubbard model as a minimal, unbiased setting for interaction-driven altermagnetism.
\end{abstract}
```

## Fig. 1 — model + Δn(k)  (fig1_model_dnk.png)
```latex
\caption{(a) Schematic diagram of the checkerboard Hubbard model with isotropic nearest-neighbor hopping
$t$ and anisotropic diagonal hopping $t_\pm = t' \pm \delta$ on alternating plaquettes. (b)--(d) Mean-field
momentum-resolved spin polarization $\Delta n(\mathbf{k})$ for $\delta = 0,\,0.2,\,0.4$. The splitting
vanishes at $\delta = 0$ and develops a $d_{xy}$ pattern that grows with $\delta$; the state is an
altermagnet with zero net moment. Parameters $t = 1$, $t' = -0.3$; the ordered moment $M$ and chemical
potential are representative mean-field values.}
```

## Fig. 2 — Fermi surface + DOS  (fig2_fs_dos.png)
```latex
\caption{Non-interacting single-particle structure of the checkerboard model ($t=1$, $t'=-0.3$, $n=0.8$).
(a),(b) Fermi surface for $\delta=0$ and $\delta=0.4$. White is the occupied Fermi sea, light blue the
empty states, the red curve is the Fermi surface, and gold stars mark the van Hove saddle points.
Increasing $\delta$ drives a Lifshitz reconstruction and moves the van Hove saddle points to incommensurate
momenta, with energy close to the Fermi level. (c) Density of states for the same two values of $\delta$,
with dashed lines at the Fermi level. Increasing $\delta$ shifts the lower-band peak from $E=4t'$ to
$E=-4\delta$ and raises the density of states near the Fermi level.}
```

## Fig. 3 — magnetic phase diagram  (3_phase_diagram.png)
Dominant S^z(k) peak in the (filling n, anisotropy δ) plane, L=14, U=4. Néel (π,π) near half-filling at all δ;
doping drives spiral (π,q) and diagonal (q,q) correlations. Spiral and diagonal share the q/π colorbar
(triangle vs star); no stripe point appears. Written in the group's PRB voice.
```latex
\begin{figure}[t]
  \centering
  \includegraphics[width=\columnwidth]{figures/3_phase_diagram.png}
  \caption{Distribution of the dominant $S^{z}(\mathbf{k})$ peak location in the $(n,\delta)$ plane for an
  $L=14$ lattice at $U=4$ and $t'=-0.3$. Each point is labeled by the ordering vector $\mathbf{Q}$. Three
  correlations are identified as N\'eel order at $(\pi,\pi)$ (red circles), spiral correlation at
  $(\pi,q)/(q,\pi)$ (triangles), and diagonal correlation at $(q,q)$ (stars). The color encodes the
  wavevector $q/\pi$. N\'eel correlations dominate near half filling for all $\delta$, while doping drives
  the system toward spiral and diagonal correlations.}
  \label{fig:phase_diagram}
\end{figure}
```

## Fig. 4 — R_p sharpness  (4_rp.png)
Sharpness index R_p of the dominant S^z(k) peak in the (n, δ) plane, L=14, U=4. Large R_p = sharp peak =
long-range order; small R_p = broad peak = short-range. Sharpest near half-filling (Néel), falls off with doping.
```latex
\begin{figure}[t]
  \centering
  \includegraphics[width=\columnwidth]{figures/4_rp.png}
  \caption{Sharpness index $R_{p}$ of the dominant $S^{z}(\mathbf{k})$ peak in the $(n,\delta)$ plane for an
  $L=14$ lattice at $U=4$ and $t'=-0.3$. A large $R_{p}$ marks a sharp peak and long-range magnetic order,
  while a small $R_{p}$ marks a broad peak and short-range correlation. The sharpest peaks appear near half
  filling, consistent with the N\'eel order, and $R_{p}$ decreases upon doping.}
  \label{fig:rp}
\end{figure}
```

## Fig. 5 — SDW momentum + real space  (5_sdw.png)
Wide 4-panel (figure*, 0.8\textwidth). (a)-(c) S^z(k) maps at fixed δ=0.20 for n=0.929 (Néel), 0.847
(spiral), 0.745 (diagonal); dominant peak walks (π,π) → (π,q)/(q,π) → (q,q) as n drops. (d) C_SDW(i) for
the same three n; Néel staggered/long-range, doped points decay fast.
INTERNAL NOTE: panel (d) is the x-direction cut only (C_SDW(i,0)); define C_SDW(i) as the x-correlation in
text, and do NOT describe (d) as "the decay of the diagonal order" (diagonal's true range is along x+y).
```latex
\begin{figure*}[t]
  \centering
  \includegraphics[width=0.8\textwidth]{figures/5_sdw.png}
  \caption{Momentum-space and real-space spin correlation at fixed anisotropy $\delta=0.20$ for three fillings. (a)--(c) Spin structure factor $S^{z}(\mathbf{k})$ for $n=0.929$, $0.847$, and $0.745$. As the filling decreases, the dominant peak moves from $(\pi,\pi)$ to $(\pi,q)/(q,\pi)$ and then to $(q,q)$. The peak intensity is largest at $n=0.929$, where $S^{z}(\pi,\pi)$ is more than an order of magnitude stronger than at the doped points. (d) Real-space correlation $C_{\mathrm{SDW}}(i)$ for the same three fillings. At $n=0.929$ the correlation shows a robust staggered structure that persists over the accessible distances, while at lower filling it decays rapidly.}
  \label{fig:sdw}
\end{figure*}
```

## Fig. 6 — pairing vertex crossover  (6_pairing.png)
Two-panel (single column). (a) all pairing channels N_ζ^Vertex vs δ at n=0.622 (on-site s repulsive/negative
from U, ext-s largest, d-waves and p-waves). (b) zoom on d_{x2-y2} vs d_{xy} for n=0.622 (solid) and n=0.582
(dashed); d_{x2-y2} falls, d_{xy} rises, cross near δ≈0.1–0.2. HONEST CLAIM = d-wave crossover among the
UNCONVENTIONAL channels, NOT global dxy dominance (ext-s larger = conventional/AFM-tied). EQUAL-TIME vertex
(not susceptibility) — CP-bias caveat pending the ED gate + 6×6/8×8 susceptibility.
```latex
\begin{figure}[t]
  \centering
  \includegraphics[width=\columnwidth]{figures/6_pairing.png}
  \caption{Equal-time connected pairing vertex $N^{\mathrm{Vertex}}_{\zeta}$ as a function of the anisotropy $\delta$ at $U=4$, $L=14$. (a) All pairing channels at fixed filling $n=0.622$. The on-site $s$ channel is strongly negative, reflecting the on-site Coulomb repulsion, and the extended $s$ channel is the largest. (b) The two $d$-wave channels for $n=0.622$ (solid) and $n=0.582$ (dashed). As $\delta$ increases the $d_{x^2-y^2}$ vertex decreases and the $d_{xy}$ vertex increases, and the two cross near $\delta\approx0.1$ to $0.2$. Among the unconventional pairing channels the checkerboard anisotropy selects $d_{xy}$.}
  \label{fig:pairing_crossover}
\end{figure}
```

## Fig. 7 — pairing phase diagram over (n, δ)  (7_pairing_pd.png)
3-panel scatter over (n,δ): (a) N^Vertex_dx2-y2 (bright small δ), (b) N^Vertex_dxy (bright large δ),
(c) difference (red=dxy dominant, blue=dx2-y2). Whole-plane version of the Fig 6 crossover; analog of paper-2
Fig 2(b,c). Frame: in the spin-INDEPENDENT checkerboard the crossover is geometry-driven, vs paper-2's imposed
spin-dependent t'. [cite paper 2]. figure*, 0.8\textwidth.
```latex
\begin{figure*}[t]
  \centering
  \includegraphics[width=0.8\textwidth]{figures/7_pairing_pd.png}
  \caption{Equal-time connected pairing vertex in the filling--anisotropy $(n,\delta)$ plane at $U=4$, $L=14$. (a) The $d_{x^2-y^2}$ vertex is largest at small $\delta$. (b) The $d_{xy}$ vertex is largest at large $\delta$. (c) The difference $N_{d_{xy}}-N_{d_{x^2-y^2}}$, positive (red) where $d_{xy}$ dominates and negative (blue) where $d_{x^2-y^2}$ dominates. The dominant $d$-wave pairing symmetry changes from $d_{x^2-y^2}$ to $d_{xy}$ across the plane.}
  \label{fig:pairing_pd}
\end{figure*}
```

## Fig. 8 — bare bubble (Unpair) vs connected vertex  (8_d_crossover.png)
Kills the kinematic objection. (a) uncorrelated/bubble (Unpair): dx2-y2 & dxy nearly flat, cross δ≈0.32.
(b) connected vertex: dx2-y2 collapses, dxy rises, cross δ≈0.22. => geometry sets the DIRECTION, correlations
DRIVE & SHARPEN. Honest: NOT purely correlation (bare drifts too) but interaction-dominated. This is the old
"Fig A / bare-vs-vertex". The k-space vertex maps figure (old Fig 8) is CUT from this set.
```latex
\begin{figure}[t]
  \centering
  \includegraphics[width=\columnwidth]{figures/8_d_crossover.png}
  \caption{Decomposition of the $d$-wave pairing correlations at $n=0.622$ into (a) Unpair and (b) Vertex. The unpair varies weakly with $\delta$ and the two $d$-wave channels cross near $\delta\approx0.32$, while the connected vertex shows a sharper crossover near $\delta\approx0.22$ with the $d_{x^2-y^2}$ channel strongly suppressed. The lattice geometry sets the direction of the crossover and the on-site repulsion drives and sharpens it.}
  \label{fig:bare_vs_vertex}
\end{figure}
```

## PAIRING SUSCEPTIBILITY (Figs 9-11) — the CORRECT observable, replaces Figs 6-8
Figs 6-8 use the EQUAL-TIME vertex, which ED proved is CP-bias sign-flipped for d-wave. Figs 9-11 use the
DYNAMIC (τ-integrated) susceptibility, ED-validated sign-faithful. HEADLINE RESULT (channel swap): the
anisotropy SUPPRESSES d_xy pairing and ENHANCES the conventional d_{x2-y2} channel (at n≈0.78, χ_d grows ~2×
from δ=0→0.4). Emergent AM decouples pairing symmetry from magnetic symmetry (opposite to imposed-AM papers).
All from the qmc-platform CP-AFQMC engine (code mirrored in cb_model/code; runs on 251/qmc48).

### Fig. 9 — pairing χ vs δ  (9_chi_delta.png)  [n=0.778, 6x6]
χ_d (dx2-y2, blue) attractive and GROWS with δ (+0.46→+0.93); χ_dxy (green) suppressed, driven negative.
```latex
\begin{figure}[t]
  \centering
  \includegraphics[width=\columnwidth]{figures/9_chi_delta.png}
  \caption{Connected-vertex pairing susceptibility $\chi^{\mathrm{vertex}}_{\zeta}$ versus anisotropy $\delta$
  at $n=0.78$ ($6\times6$, $U=4$). The $d_{x^2-y^2}$ vertex is attractive and grows with $\delta$, while the
  $d_{xy}$ vertex is repulsive and is driven increasingly negative. The geometric anisotropy enhances
  $d_{x^2-y^2}$ pairing and suppresses $d_{xy}$ pairing.}
  \label{fig:chi_delta}
\end{figure}
```

### Fig. 10 — pairing χ vs n  (10_chi_n.png)  [δ=0.4, 6x6]
χ_dxy repulsive across the doped range, positive only at half filling (magnetic region); χ_d attractive in
the doped window n≈0.67–0.89.
```latex
\begin{figure}[t]
  \centering
  \includegraphics[width=\columnwidth]{figures/10_chi_n.png}
  \caption{Connected-vertex pairing susceptibility versus filling $n$ at fixed anisotropy $\delta=0.4$
  ($6\times6$, $U=4$). The $d_{xy}$ vertex is repulsive throughout the doped range and turns weakly positive
  only near half filling, where the system is magnetically ordered. The $d_{x^2-y^2}$ vertex is weakly
  attractive at intermediate doping.}
  \label{fig:chi_n}
\end{figure}
```

### Fig. 11 — magnetic χ_zz(q)  (11_chizz_q.png)  [6x6]
Momentum-resolved magnetic spin susceptibility. Half filling peaks at (π,π) (Néel/AM); doped is weaker and
shifts off (π,π). The momentum-space magnetic susceptibility.
```latex
\begin{figure*}[t]
  \centering
  \includegraphics[width=0.85\textwidth]{figures/11_chizz_q.png}
  \caption{Momentum-resolved magnetic susceptibility $\chi_{zz}(\mathbf{q})$ over the Brillouin zone
  ($6\times6$, $U=4$). Near half filling the susceptibility peaks at $(\pi,\pi)$, signalling the N\'eel
  altermagnetic order. Upon doping the peak weakens and moves away from $(\pi,\pi)$.}
  \label{fig:chizz}
\end{figure*}
```

### Fig. 12 — staggered magnetic susceptibility vs n  (12_chizz_vs_n.png)
χ_zz(π,π) vs filling at δ=0.4; rises toward half filling (0.019→0.117, peak ~n=0.89). The magnetic-side
scalar susceptibility.
```latex
\begin{figure}[t]
  \centering
  \includegraphics[width=\columnwidth]{figures/12_chizz_vs_n.png}
  \caption{Staggered magnetic susceptibility $\chi_{zz}(\pi,\pi)$ versus filling $n$ at $\delta=0.4$
  ($6\times6$, $U=4$). It increases toward half filling, tracking the growth of the $(\pi,\pi)$
  antiferromagnetic (altermagnetic) correlations.}
  \label{fig:chizz_n}
\end{figure}
```

### Fig. 13 — magnetism vs pairing overlay  (13_mag_vs_pair.png)
Twin-axis: χ_zz(π,π) (magnetic) and χ_{dx2-y2} (pairing) vs n at δ=0.4. dx2-y2 turns attractive in the same
doped window where χ_zz(π,π) rises => dx2-y2 pairing TRACKS the (π,π) AFM susceptibility (AFM-mediated d-wave,
cuprate-like), while dxy is suppressed. The mechanism figure.
```latex
\begin{figure}[t]
  \centering
  \includegraphics[width=\columnwidth]{figures/13_mag_vs_pair.png}
  \caption{Magnetic and pairing susceptibilities versus filling at $\delta=0.4$ ($6\times6$, $U=4$). The
  $d_{x^2-y^2}$ pairing vertex (right axis) becomes attractive in the same doped window where the staggered
  magnetic susceptibility $\chi_{zz}(\pi,\pi)$ (left axis) grows, indicating that the $d_{x^2-y^2}$ pairing is
  driven by the $(\pi,\pi)$ antiferromagnetic correlations. The $d_{xy}$ channel, which does not couple to
  these correlations, is suppressed.}
  \label{fig:mag_vs_pair}
\end{figure}
```

### Reframed pairing paragraph (manuscript-ready)
> To determine whether the emergent altermagnetism drives an unconventional pairing instability, we compute
> the singlet pair-field susceptibility $\chi_\alpha = \int_0^\beta \langle \Delta_\alpha(\tau)\Delta_\alpha^\dagger(0)\rangle\, d\tau$
> and its connected (vertex) part for the extended-$s$, $d_{x^2-y^2}$, and $d_{xy}$ channels, benchmarked
> against exact diagonalization. This dynamic, imaginary-time-integrated susceptibility is the appropriate
> diagnostic. The equal-time pairing vertex is instead subject to a constrained-path sign bias in this regime:
> on the clusters where exact diagonalization is available, the connected $d$-wave equal-time vertex has the
> opposite sign to the exact result, whereas the sign of the $\tau$-integrated susceptibility is reproduced
> faithfully. We find that the connected $d_{xy}$ vertex is repulsive across the doped range and is driven
> increasingly negative by the anisotropy, so the geometry suppresses $d_{xy}$ pairing. The conventional
> $d_{x^2-y^2}$ channel behaves oppositely: at intermediate doping ($n\approx0.78$) its connected vertex is
> attractive and is enhanced by the anisotropy, growing by roughly a factor of two as $\delta$ increases from
> $0$ to $0.4$. The anisotropy therefore enhances $d$-wave pairing correlations, but in the $d_{x^2-y^2}$
> channel rather than the $d_{xy}$ channel. This distinguishes the emergent altermagnet from models in which
> altermagnetism is imposed through spin-dependent hopping, where the enhanced pairing shares the symmetry of
> the imposed order. In the present spin-independent model the pairing enhancement occurs in the channel
> orthogonal to the $d_{xy}$ symmetry of the altermagnet, so the magnetic and pairing symmetries decouple.

## Relation-to-prior-work paragraph
See the Research Plan docx (Stage 3) and the drafted paragraph positioning against Giuli et al.
PRB 111 L020401 (2025), arXiv:2507.00837, and the group's imposed-AM papers.
