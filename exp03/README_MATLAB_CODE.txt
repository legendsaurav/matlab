MATLAB CODE -- Controller Design on MATLAB Platform Using Analog Root Locus
============================================================================
Group 3C, IIT Ropar EE. Companion code to the Extended Analysis Report.

WHAT THIS IS
------------
Complete, runnable MATLAB scripts that regenerate every figure and table
in the extended report (Figures 1-22, Tables 1-12). Each of the five
"Part" scripts corresponds to one section of the report's analysis.
Every script was actually executed (via GNU Octave 8.4 + the Octave
control package, which mirrors MATLAB's Control System Toolbox syntax)
before delivery -- not just written and assumed correct.

FILES
-----
Run these five, in any order (each is self-contained):
  matlab_partA_rootlocus.m                    Figures 2-7,  Tables 1-2
  matlab_partB_damping_stability_sensitivity.m Figures 8,9,21, Tables 3-5,10
  matlab_partC_overshoot_settling_step.m       Figures 10-13,  Tables 6-7
  matlab_partD_extended_analysis.m             Figures 14-15,  Table 8
  matlab_partE_multi_a_and_heatmaps.m          Figures 16-20,  Table 12  [NEW]

Shared helper functions (keep these on the MATLAB path -- same folder is
simplest):
  charEqCoeffs.m          s^3+a*s^2+K*s+K coefficients
  dominantComplexPole.m   robust selector for the oscillatory pole
  realClosedLoopPole.m    selector for the real (non-oscillatory) pole
  robustStepResponse.m    true step response with guaranteed peak resolution
  breakawayPoints.m       breakaway/break-in solver for a given a
  viridisLike.m / plasmaLike.m / magmaLike.m   dependency-free colormaps

Figure 1 (block diagram) and Figure 22 (op-amp circuit) are reproduced
images from the original submission, not code-generated -- there is no
script for those two.

REQUIREMENTS
------------
MATLAB with the Control System Toolbox (for tf, rlocus, step). All five
scripts use only that one toolbox; no other add-ons are needed.

Each script creates its own output folder (figures_partA/ ... _partE/)
next to itself and writes PNGs at 150 dpi via print(...,'-dpng','-r150').
Figures are created with 'visible','off' so the scripts run unattended;
delete that name-value pair on any figure(...) call if you'd rather watch
them appear.

A NOTE ON ACCURACY -- PLEASE READ
----------------------------------
While testing this code I found that the "true simulated" overshoot
numbers in Table 6 of the Word report are slightly UNDER-estimated for
the lightly-damped, small-a cases (worst case: a=2, K=10000, where the
report says 97.70% but the correct value is 98.46% -- about 0.76 points
low). The cause: the original Python analysis behind the report gave
scipy's step-response solver only a moderately fine, uniformly-spaced
time grid, which under-resolves a sharp early peak when damping is very
light (zeta ~ 0.002-0.005). I verified this against a 3-million-point
reference grid.

This MATLAB code does NOT have that problem -- robustStepResponse.m
explicitly guarantees at least 200 samples before the analytically
predicted first-peak time, and its output matches the high-precision
reference to within 0.01 points across every case I checked. So the
numbers this code produces are the corrected ones, matching the current
(also corrected) Word report exactly.

SECOND FIX -- robustStepResponse.m (updated since first delivery)
-------------------------------------------------------------------
A second, unrelated bug was found and fixed after Parts A-D were first
delivered: the original robustStepResponse.m estimated its simulation
time window from Tpeak = pi/(wn*sqrt(1-zeta^2)), which explodes as
zeta -> 1. That happens for a > 9 at gains INSIDE the breakaway/break-in
loop of Sec. 5.3 (all three closed-loop poles real, so zeta from the
"dominant complex pole" selector is ~1, not a meaningful damping ratio).
The result: a multi-thousand-second time window spread over only 4000
points, so coarse that the true transient -- including the peak -- was
missed entirely, silently returning ~0% overshoot instead of the correct
value (example: a=20, K=90 returned 0.00% instead of 15.23%).

The fixed version (included in this package) instead sizes the window
from the raw pole magnitudes (slowest pole -> settling time, fastest
pole -> resolution needed), which has no equivalent failure mode for
real poles. Verified against fine-grid reference simulations across
lightly-damped, heavily-damped, and all-real-pole regimes.

Practical impact on what you already have: NONE of the figures in the
originally-delivered Parts A-D are affected -- I checked every call site.
Part C's grid uses K >= 10,000, always well outside any a's loop window;
Part D's ringing demo passes an explicit TendOverride, which happens to
keep even the old buggy code correct for that specific case. The bug
only mattered for Part E's new continuous K-sweep (Figure 17), which is
built with the fixed function and has been verified correct.

PART E -- NEW SINCE THE ORIGINAL DELIVERY
--------------------------------------------
matlab_partE_multi_a_and_heatmaps.m reproduces the report's newer
Sections 10.3-10.4:
  - Figure 16: 2x2 multi-a step-response waterfall (a = 2, 5, 9, 20)
  - Figure 17 + Table 12: true simulated PO(K) for the same four a,
    on one continuous K sweep
  - Figures 18-20: three consistent heatmaps on the identical (K,a)
    grid -- damping ratio zeta, peak overshoot PO, and 5% settling
    time Ts -- all derived from a single pole computation per grid
    point
  - A numerical check of the closed-form asymptote Ts -> 6/(a-1) as
    K -> infinity (derived from the report's own Sec. 6.2 high-gain
    pole asymptotics; confirmed to 4 significant figures)

It only depends on functions already in this package (charEqCoeffs,
dominantComplexPole, robustStepResponse, viridisLike) -- no new helpers
needed. Uses MATLAB's built-in 'parula' and 'hot' colormaps (both
standard since early MATLAB releases).
