# Upsilon invariant-mass fit model

> **Validation status:** the archived fit numbers and uncertainties are **unvalidated**. The saved minimizer output contains a call-limit termination and non-positive-definite covariance matrices. The notebook establishes code and figure provenance, not an accepted numerical result.

## Source of record

The retained implementation is `notebooks/01_mass_fits_inclusive_and_rapidity.ipynb`, an exact source copy of the original `MassFit.ipynb`.

| Notebook cell | Analysis | Input | Saved output |
| --- | --- | --- | --- |
| 0 | Inclusive fit to an existing ROOT histogram | `~/merged_elecs_test.root`, object `h_invMass_Reco_All` | `Fit_Inclusive_Upsilon_DSCB_Scaled.pdf` |
| 1 | Inclusive, weighted tuple fit | `~/merged_elecs_test_10t100_10x100.root:UpsilonTree` | `Fit_Inclusive_Upsilon_DSCB_Nice.png` |
| 2 | Two rapidity-bin, weighted tuple fits | `~/merged_elecs_test_tuples_fix2.root:UpsilonTree` | `Fit_Inclusive_Upsilon_RapidityBins_{STATE_NORM_MODE}_{BR_WEIGHT_MODE}.png` |

Cells are numbered from zero, as in the notebook JSON. Cell 1 selects reconstructed mass from 8 to 12 GeV and `reco_W >= 0`; cell 2 additionally requires `abs(sigma_epz - 36) <= 3.6` and fits `-1 <= reco_y < 0` and `0 <= reco_y < 1`. Both tuple implementations use 200 bins over 8--12 GeV (20 MeV per bin).

## Mathematical definition

The fitted function has three Double-Sided Crystal Ball (DSCB) components plus an exponential background:

\[
F(m;\boldsymbol\theta) =
D_1(m) + D_2(m) + D_3(m) + B(m).
\]

For a state with parameters \((N,\mu,\sigma)\) and the shared tail parameters \((\alpha_L,n_L,\alpha_H,n_H)\), define

\[
t = \frac{m-\mu}{\sigma},\qquad
A_i = \left(\frac{n_i}{\alpha_i}\right)^{n_i}e^{-\alpha_i^2/2},\qquad
C_i = \frac{n_i}{\alpha_i}-\alpha_i.
\]

The notebook implements

\[
D(m) = N\begin{cases}
A_L(C_L-t)^{-n_L}, & t < -\alpha_L,\\
e^{-t^2/2}, & -\alpha_L \leq t \leq \alpha_H,\\
A_H(C_H+t)^{-n_H}, & t > \alpha_H.
\end{cases}
\]

The background is

\[
B(m) = N_{\mathrm{bkg}}e^{k_{\mathrm{bkg}}m}.
\]

The DSCB `Norm` parameters are peak amplitudes, not normalized event yields. The notebook obtains a component yield by numerically integrating its fitted DSCB over 8--12 GeV and dividing by the histogram bin width.

## Exact 15-parameter map

There is no `FixParameter` call in any fit cell. **All 15 parameters are free.** Limits constrain some parameters but do not fix them.

| Index | ROOT name | Role | Shared? | Limits in cells 0/1 | Limits in cell 2 |
| ---: | --- | --- | --- | --- | --- |
| 0 | `Norm 1S` | 1S amplitude | No | None | None |
| 1 | `Mean 1S` | 1S mean, GeV | No | 9.3--9.6 | 9.3--9.6 |
| 2 | `Sigma 1S` | 1S width, GeV | No | 0.01--0.15 | 0.01--0.15 |
| 3 | `Alpha L` | Low-side transition | Yes, all states | 0.1--5 | 0.1--5 |
| 4 | `n L` | Low-side power | Yes, all states | 1--20 | 1--20 |
| 5 | `Alpha H` | High-side transition | Yes, all states | 0.1--5 | 0.1--5 |
| 6 | `n H` | High-side power | Yes, all states | 1--20 | 1--20 |
| 7 | `Norm 2S` | 2S amplitude | No | None | None |
| 8 | `Mean 2S` | 2S mean, GeV | No | 9.9--10.15 | 9.9--10.15 |
| 9 | `Sigma 2S` | 2S width, GeV | No | 0.01--0.15 | 0.01--0.15 |
| 10 | `Norm 3S` | 3S amplitude | No | None | None |
| 11 | `Mean 3S` | 3S mean, GeV | No | 10.2--10.5 | 10.2--10.5 |
| 12 | `Sigma 3S` | 3S width, GeV | No | 0.01--0.15 | 0.01--0.15 |
| 13 | `Bkg Norm` | Exponential amplitude | No | None | None |
| 14 | `Bkg Slope` | Exponential slope | No | None | -10--1 |

The three states have independent amplitudes, means, and widths. Parameters 3--6 are reused for the 2S and 3S components, so the tail shape is common to all states. Signal amplitudes and background amplitude are not constrained to be non-negative. Cell 2 adds numerical guards (`abs`/minimum values for width and transition parameters, power-domain checks, and overflow handling); those guards do not make any parameter fixed.

## Weight and normalization modes

- Cell 1 currently sets `BR_WEIGHT_MODE = "visible_decay"`; the histogram weight is therefore the tuple `weight`. The alternative `inclusive_upsilon` mode multiplies each truth state by the inverse dielectron branching ratio.
- Cell 2 currently sets `STATE_NORM_MODE = "raw"` and `BR_WEIGHT_MODE = "visible_decay"`; its `fit_weight` is the raw tuple `weight`. The alternative state mode rescales truth-state weight sums to configured visible cross sections before fitting, and the alternative branching-ratio mode unfolds the dielectron branching ratio.
- These modes change the fitted histogram normalization. They must be recorded in the output filename and run metadata; results from different modes are not directly interchangeable.

## Post-fit STARlight scaling is not a fit constraint

Cell 0 declares `R2S_STARLIGHT = 0.445` and `R3S_STARLIGHT = 0.315`. It first fits independent 1S, 2S, and 3S amplitudes, then computes

\[
Y_i = \frac{1}{\Delta m}\int_{8\,\mathrm{GeV}}^{12\,\mathrm{GeV}}D_i(m)\,dm,
\]

and only afterward prints \(0.445Y_{2S}\) and \(0.315Y_{3S}\). The ratios do **not** tie `Norm 2S` or `Norm 3S` to `Norm 1S`, do not enter the objective function, and do not alter the drawn fitted components. These values must be described as post-fit scaled yields, never as constrained-fit yields.

## Cached minimizer evidence

The outputs saved inside the notebook record the following failures:

| Cell | Cached ROOT status | Additional evidence |
| ---: | --- | --- |
| 0 | `FCN=150.802 FROM HESSE STATUS=NOT POSDEF` | `EDM=8.97974e-07`; `ERR MATRIX NOT POS-DEF` |
| 1 | `FCN=109.858 FROM MIGRAD STATUS=CALL LIMIT 3651 CALLS 3652 TOTAL` | `EDM=0.128021`; `ERROR MATRIX UNCERTAINTY 15.0 per cent`; `Warning in <Fit>: Abnormal termination of minimization.` |
| 2, rapidity -1 to 0 | `FCN=114.943 FROM HESSE STATUS=NOT POSDEF` | `EDM=2.99153`; `ERR MATRIX NOT POS-DEF` |
| 2, rapidity 0 to 1 | `FCN=137.769 FROM HESSE STATUS=NOT POSDEF` | `EDM=6.26884e-08`; `ERR MATRIX NOT POS-DEF` |

Cell 1 also cached `Norm 2S = -2.52946e-02`, `Bkg Norm = -7.57927`, and integrated yields `348.75`, `-0.04`, and `0.94` for 1S, 2S, and 3S. The first rapidity fit cached a negative background amplitude. These are direct consequences of an unconstrained, non-converged fit and are not physical measurements.

## Acceptance requirements for automated fits

An automated workflow must treat a plot as failed, and must not publish its fitted numbers, unless all of the following are recorded and pass:

1. The minimizer terminates normally and does not report a call limit.
2. The covariance matrix is positive definite and has an acceptable ROOT covariance status.
3. EDM and parameter-at-bound checks meet explicit thresholds.
4. Signal and background components satisfy the chosen physical-positivity policy.
5. Fit quality, residuals or pulls, and stability under initialization changes are checked.
6. Dataset identity, selection, weight mode, normalization mode, code revision, and output path are captured in run metadata.

Until the model is rerun under those gates, **all masses, widths, yields, parameter errors, and covariance-derived quantities cached in the notebook remain unvalidated**.
