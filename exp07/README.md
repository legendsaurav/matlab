# EE208 Experiment 7 — MATLAB package

Discrete root-locus design of a proportional controller (5 % overshoot, 2 % settling time ≤ 9 s)
for the lab-sheet plant at T = 1 s

    G(z) = 0.32 (z + 0.62) / [ (z - 1)(z - 0.15) ]

and for the sampling-time sweep T = 0.1 … 1.0 s.

## Plant model used for the T sweep

The sweep needs G(z) at other sampling times, so the script recovers the underlying
continuous plant G(s) = K / [s(s+a)] from the T = 1 s sheet:

* a = -ln(0.15)/1 = 1.8971
* K = 0.32·a² / (a - 1 + 0.15) = 1.0999

**Important:** the sheet's G(z) is *not* an exact zero-order-hold equivalent of this G(s).
The exact ZOH model has its zero at -0.540, the sheet says -0.62 (14.8 % difference), so the
model gives Kp,max = 0.8524 at T = 1 s versus 0.7845 for the sheet coefficients.
The script prints this mismatch every run. The switch

    USE_GIVEN_AT_T1 = true;   % section 0, main_experiment7.m

makes every T = 1 s section use the exact sheet G(z) (default); set it to `false` to use the
recovered model everywhere (then the T = 1 s row is consistent with the rest of the sweep).

## Files

| File | Purpose |
|---|---|
| `main_experiment7.m` | **Main script** – runs the whole experiment, prints all tables, draws Figures 1–16 |
| `verify_report_numbers.m` | Pass/fail check of the function library (run first) |
| `zoh_coeffs.m` | Closed-form ZOH coefficients b1, b0, p for any T |
| `effective_zeta.m` | T-independent damping ratio of a z-plane pole |
| `breakaway_gain.m` | Exact ζ = 1 gain (discriminant of the characteristic equation) |
| `zeta_error.m` / `zeta_crossing_gain.m` | Gain where the locus crosses a ζ contour (fzero) |
| `evaluate_design.m` | Simulated %OS and 2 % settling time |
| `marginal_stability_gain.m` | Exact stability limit Kp = (1-p)/b0 |
| `fit_inverse_law.m` | Fit Kp(T) ≈ Kp0 / (1 + cT) |
| `pole_sensitivity.m` | dz/dKp by implicit differentiation |

## Run

1. Put all `.m` files in one folder and `cd` into it.
2. `verify_report_numbers` → should end with `ALL CHECKS PASSED`.
3. `main_experiment7` (or run it section by section with Ctrl+Enter).

Each run creates `assets/` and `artifacts/` beside the scripts. Generated
figures are saved as numbered PNG files in `assets/`; tables, sweep data,
heatmap matrices, and robustness/lead-search results are saved as CSV files
in `artifacts/`.

## Expected key numbers

| Quantity | Sheet G(z), T = 1 s | Recovered model, T = 1 s |
|---|---|---|
| Kp,min (ζ = 1, breakaway) | 0.4883 | 0.5261 |
| Kp,max (ζ = 0.6901, 5 % OS) | 0.7845 | 0.8524 |
| Kp,stab (marginal stability) | 4.284 | 4.919 |

Sweep (model): Kp,max falls from 1.5613 (T = 0.1 s) to 0.8524 (T = 1 s);
Kp,min from 0.7803 to 0.5261. A design at Kp,min (zeta = 1) has ts = 9.0 s at T = 1 s (zero margin).
Sec. 11 first scans the proportional window: ts = 5 s is already reachable for Kp of about
0.635-0.678 (OS about 1.5 %), so the lead compensator (zero cancels the plant pole z = 0.15)
is only accepted if it beats that: the script picks C(z) = 1.5922 (z - 0.15)/(z + 0.20),
OS = 4.91 %, ts = 4 s on the sheet plant (`ts_lead = 5`, `lead_margin = 1`).

## Requirements

MATLAB R2018b or newer (`yline`/`xline`; local functions in scripts need R2016b).
Control System Toolbox for `tf, c2d, ss, feedback, stepinfo, margin, bandwidth,
nyquistplot, rlocus, zgrid`. `fzero` and `fminsearch` are base MATLAB.

## Notes

* The code was written and cross-checked against an independent Python implementation,
  but it has not been run in MATLAB itself. Run `verify_report_numbers` first.
* Section 9.2 (Control System Designer) is a GUI step; the call is commented out.
* The heatmaps and the nomograph simulate a few thousand closed loops (1–3 min).
