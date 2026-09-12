# Reading list (benchmark literature)

## A. Closest — interaction-driven altermagnetism (position against these)
- Giuli, Mejuto-Zaera, Capone, "Altermagnetism from interaction-driven itinerant magnetism,"
  PRB 111, L020401 (2025). arXiv:2410.00909. **The key competitor.** Two-orbital, van-Hove-driven,
  gRISB (variational). We differ: single-orbital, geometry-driven, Néel near half-filling, unbiased QMC.
- "Spontaneous emergence of altermagnetism in the single-orbital extended Hubbard model," arXiv:2507.00837.
  Single-orbital square lattice + non-local V, mean-field, doped. We differ: checkerboard, on-site U only,
  geometry-driven, QMC.
- "Stripe-Ordered Altermagnetism from Correlation-Driven SDW Instability," arXiv:2607.11532.
- V. Leeb, A. Mook, L. Šmejkal, J. Knolle, PRL 132, 236701 (2024) — altermagnetism from orbital ordering.
- "Dirac points and topological phases in correlated altermagnets," PRR (doi:10.1103/7nvm-s225).

## B. van Hove / t-t' Hubbard benchmark (for the FS/VHS/DOS figures)
- H. Q. Lin, J. E. Hirsch, "2D Hubbard model with NN and NNN hopping," PRB 35, 3359 (1987).
  Establishes the −4t' cos kx cos ky term and the VHS at (π,0)/(0,π) at energy 4t'. Direct sanity check.
- "Weak ferromagnetism and instabilities of the 2D t-t' Hubbard model at Van Hove fillings," cond-mat/0306296.
- "Magnetic phase diagram of the Hubbard model with NNN hopping," arXiv:0904.4755.
- "Effects of Van Hove singularities on magnetism and superconductivity in the t-t' Hubbard model (parquet)."
- "Fermi condensation near van Hove singularities (triangular lattice)," PRL 112, 070403 (2014); arXiv:1308.0812.
  Good statement of the VHS <-> Lifshitz correspondence.

## C. Altermagnet Fermi surface / d-wave spin splitting
- Šmejkal, Sinova, Jungwirth, PRX 12, 040501 (2022) — foundational, d-wave spin-split Fermi surface.
- "Unconventional p-wave and finite-momentum SC induced by altermagnetism (Bogoliubov FS)," PRB 111, 054501 (2025);
  arXiv:2407.02059.

## D. Checkerboard / planar-pyrochlore electronic structure
- "Hidden mechanism for realizing flat bands: embedding Lieb, kagome, checkerboard" (checkerboard tight-binding).
- "The planar pyrochlore: a Valence Bond Crystal" (lattice / frustration background).

## Group's own work (cite as the imposed-AM lineage)
- J. Li, J. Liu, X. Yang, H.-K. Tang, "Enhancement of d-wave pairing in strongly correlated altermagnet" (2026).
- "Competing magnetic correlations in the square-lattice Hubbard model with spin-dependent NN and NNN hopping" (Paper 1).
- "Evolution of magnetic correlation in the doped Hubbard model with altermagnetic spin splitting,"
  PRB 113, 134443 (2026) (Paper 2).

## Methods
- Scalapino, White, Zhang, PRB 47, 7995 (1993) — pairing / pair-field susceptibility.
- White et al., PRB 40, 506 (1989) — determinant QMC.

## Quick benchmark you can quote
At δ=0 the model reduces to the square lattice with NNN hopping, whose VHS sits at (π,0)/(0,π) at energy
4t' = −1.2 (with t'=−0.3). Our band code reproduces this exactly, and δ>0 pushes the VHS to (±π/2,±π/2)
at energy −4δ (dxy). This matches the t-t' literature (group B).
