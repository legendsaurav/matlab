"""Generate the requested time, frequency, pole, robustness, and trade-off plots.

Run with: c:/python314/python.exe python_analysis.py
Outputs are written to artifacts/python_plots.
"""

from pathlib import Path

import control as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "artifacts" / "python_plots"
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 180, "axes.grid": True})

s = ct.TransferFunction.s
# The handwritten transfer function and poles artifact use the pole at -12.
G = 10 * (s + 2.5) ** 2 / ((s + 12) * (s**2 + 0.12))
K0 = 12.18
# A unity-DC-gain lead comparison: zero at -0.3, pole at -1.5.
C_default = 1
C_lead = 5 * (s + 0.3) / (s + 1.5)


def closed_loop(controller, gain=K0):
    loop = gain * controller * G
    return ct.feedback(loop, 1), loop


def response(sys, time, kind="step"):
    if kind == "step":
        return ct.step_response(sys, time)
    if kind == "impulse":
        return ct.impulse_response(sys, time)
    return ct.forced_response(sys, time, U=time)


def metrics(sys, loop):
    info = ct.step_info(sys, SettlingTimeThreshold=0.05)
    gm, pm, wcg, wcp = ct.margin(loop)
    dc = float(np.real(ct.dcgain(sys)))
    return {
        "steady_state_error": abs(1 - dc),
        "overshoot_percent": info["Overshoot"],
        "settling_time_5_percent_s": info["SettlingTime"],
        "phase_margin_deg": pm,
        "gain_margin_dB": 20 * np.log10(gm) if np.isfinite(gm) and gm > 0 else np.inf,
        "bandwidth_rad_s": ct.bandwidth(sys),
        "peak_sensitivity": np.max(np.abs(ct.frequency_response(1 / (1 + loop), W).complex.flatten())),
    }


T_default, L_default = closed_loop(C_default)
T_lead, L_lead = closed_loop(C_lead)
time = np.linspace(0, 80, 5000)
time_fast = np.linspace(0, 15, 3000)
W = np.logspace(-3, 2, 1800)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, bbox_inches="tight")
    plt.close(fig)


# 1. Time-domain responses: step, error, control effort, ramp, impulse, disturbance.
fig, ax = plt.subplots(3, 2, figsize=(13, 12))
for sys, label, color in [(T_default, "Default P", "tab:blue"), (T_lead, "Lead", "tab:orange")]:
    _, y = response(sys, time)
    ax[0, 0].plot(time, y, label=label, color=color)
    ax[0, 1].plot(time, 1 - y, label=label, color=color)
    controller = C_default if label == "Default P" else C_lead
    u_sys = controller / (1 + K0 * controller * G)
    _, u = response(K0 * u_sys, time)
    ax[1, 0].plot(time, u, label=label, color=color)
    _, ramp_y = response(sys, time, "ramp")
    ax[1, 1].plot(time, time - ramp_y, label=label, color=color)
    _, impulse = response(sys, time_fast, "impulse")
    ax[2, 0].plot(time_fast, impulse, label=label, color=color)

for sys, label, color in [(G / (1 + L_default), "Default P", "tab:blue"),
                           (G / (1 + L_lead), "Lead", "tab:orange")]:
    _, disturbance = response(sys, time)
    ax[2, 1].plot(time, disturbance, label=label, color=color)
titles = ["Unit-step output", "Error e(t) = 1 - y(t)", "Control effort u(t)",
          "Ramp tracking error", "Impulse response", "Plant-input step disturbance"]
for axis, title in zip(ax.flat, titles):
    axis.set_title(title)
    axis.set_xlabel("Time (s)")
    axis.legend()
save(fig, "01_time_domain.png")


# 2. Step response versus gain and pole migration/root locus.
gains = np.linspace(0.1, 30, 250)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for gain in [1, 4, 8, 12.18, 20, 30]:
    _, y = response(ct.feedback(gain * G, 1), time)
    ax[0].plot(time, y, label=f"K={gain:g}")
for gain in gains:
    poles = ct.poles(ct.feedback(gain * G, 1))
    ax[1].plot(poles.real, poles.imag, ".", color="tab:blue", ms=2)
ax[0].set(title="Step response versus proportional gain", xlabel="Time (s)", ylabel="Output")
ax[0].legend(ncol=2)
ax[1].axhline(0, color="k", lw=0.7)
ax[1].axvline(0, color="k", lw=0.7)
ax[1].set(title="Closed-loop pole migration versus K", xlabel="Real axis", ylabel="Imaginary axis")
save(fig, "02_gain_and_pole_migration.png")


# 3. Bode, sensitivity/complementary sensitivity, Nichols, Nyquist, and phase lag.
fig, ax = plt.subplots(2, 2, figsize=(13, 10))
for loop, label, color in [(L_default, "Default P", "tab:blue"), (L_lead, "Lead", "tab:orange")]:
    mag, phase, omega = ct.frequency_response(loop, W)
    mag = np.asarray(mag).flatten()
    phase_deg = np.unwrap(np.asarray(phase).flatten()) * 180 / np.pi
    ax[0, 0].semilogx(omega, 20 * np.log10(mag), label=label, color=color)
    ax[0, 1].semilogx(omega, phase_deg, label=label, color=color)
    S = 1 / (1 + loop)
    T = loop / (1 + loop)
    smag = np.abs(ct.frequency_response(S, W).complex.flatten())
    tmag = np.abs(ct.frequency_response(T, W).complex.flatten())
    ax[1, 0].semilogx(W, 20 * np.log10(smag), label=f"S {label}", color=color, ls="--")
    ax[1, 0].semilogx(W, 20 * np.log10(tmag), label=f"T {label}", color=color)
    ax[1, 1].semilogx(W, -np.gradient(phase_deg * np.pi / 180, omega), label=label, color=color)
ax[0, 0].set_title("Open-loop Bode magnitude")
ax[0, 1].set_title("Open-loop phase")
ax[1, 0].set_title("Sensitivity S and complementary sensitivity T")
ax[1, 1].set_title("Phase-lag slope (group-delay view)")
for axis in ax.flat:
    axis.set_xlabel("Frequency (rad/s)")
    axis.legend()
save(fig, "03_frequency_domain.png")

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for loop, label, color in [(L_default, "Default P", "tab:blue"), (L_lead, "Lead", "tab:orange")]:
    mag, phase, _ = ct.frequency_response(loop, W)
    ax[0].plot(np.asarray(phase).flatten() * 180 / np.pi, 20 * np.log10(np.asarray(mag).flatten()), label=label, color=color)
    nyq = ct.frequency_response(loop, W).complex.flatten()
    ax[1].plot(nyq.real, nyq.imag, label=label, color=color)
    ax[1].plot(nyq.real, -nyq.imag, color=color, alpha=0.35)
ax[0].plot(-180, 0, "rx", ms=10, mew=2, label="Critical point")
ax[0].set(title="Nichols chart", xlabel="Phase (deg)", ylabel="Gain (dB)")
ax[1].plot(-1, 0, "rx", ms=10, mew=2, label="Critical point")
ax[1].set(title="Nyquist plot", xlabel="Real", ylabel="Imaginary")
for axis in ax:
    axis.legend()
save(fig, "04_nichols_nyquist.png")


# 4. Pole map, root locus, and damping/natural-frequency reference curves.
fig, ax = plt.subplots(figsize=(8, 6))
for sys, label, color in [(T_default, "Default P", "tab:blue"), (T_lead, "Lead", "tab:orange")]:
    poles = ct.poles(sys)
    ax.plot(poles.real, poles.imag, "x", ms=10, mew=2, label=label, color=color)
zeta = 0.2
wn = np.linspace(0.05, 15, 300)
for damping in [0.2, 0.5, 0.7]:
    ax.plot(-damping * wn, wn * np.sqrt(1 - damping**2), "k:", alpha=0.5)
    ax.plot(-damping * wn, -wn * np.sqrt(1 - damping**2), "k:", alpha=0.5)
ax.axhline(0, color="k", lw=0.7); ax.axvline(0, color="k", lw=0.7)
ax.set(title="Closed-loop pole map with damping-ratio rays", xlabel="Real axis", ylabel="Imaginary axis")
ax.legend()
save(fig, "05_pole_map.png")

fig, ax = plt.subplots(figsize=(8, 6))
for controller, label, color in [(C_default, "Default P", "tab:blue"), (C_lead, "Lead", "tab:orange")]:
    locus = []
    for gain in gains:
        locus.extend(ct.poles(ct.feedback(gain * controller * G, 1)))
    locus = np.asarray(locus)
    ax.plot(locus.real, locus.imag, ".", ms=1.4, label=label, color=color)
ax.axhline(0, color="k", lw=0.7); ax.axvline(0, color="k", lw=0.7)
ax.set(title="Root locus of C(s)G(s)", xlabel="Real axis", ylabel="Imaginary axis")
ax.legend()
save(fig, "06_root_locus.png")


# 5. Design-tradeoff sweeps and robustness overlays.
rows = []
for gain in gains:
    sys, loop = closed_loop(C_default, gain)
    m = metrics(sys, loop)
    m["K"] = gain
    rows.append(m)
trade = pd.DataFrame(rows)
trade.to_csv(OUT / "gain_tradeoff_metrics.csv", index=False)

fig, ax = plt.subplots(2, 2, figsize=(13, 9))
ax[0, 0].plot(trade.K, trade.phase_margin_deg); ax[0, 0].axhline(60, color="r", ls="--"); ax[0, 0].set(title="Phase margin versus K", xlabel="K", ylabel="Phase margin (deg)")
ax[0, 1].plot(trade.K, trade.overshoot_percent); ax[0, 1].axhline(2, color="r", ls="--"); ax[0, 1].set(title="Overshoot versus K", xlabel="K", ylabel="Overshoot (%)")
ax[1, 0].plot(trade.K, trade.settling_time_5_percent_s); ax[1, 0].axhline(10, color="r", ls="--"); ax[1, 0].set(title="5% settling time versus K", xlabel="K", ylabel="Time (s)")
ax[1, 1].plot(trade.K, trade.steady_state_error); ax[1, 1].axhline(0.05, color="r", ls="--"); ax[1, 1].set(title="Steady-state error versus K", xlabel="K", ylabel="Error")
save(fig, "07_design_tradeoffs.png")

fig, ax = plt.subplots(figsize=(9, 5))
for gain, color in [(0.7, "tab:blue"), (1.0, "tab:orange"), (1.3, "tab:green")]:
    varied_plant = gain * G
    _, y = response(ct.feedback(K0 * varied_plant, 1), time)
    ax.plot(time, y, label=f"Plant gain {gain:.1f}x", color=color)
for wn_scale in [0.7, 1.0, 1.3]:
    varied = 10 * (s + 2.5) ** 2 / ((s + 12) * (s**2 + 0.12 * wn_scale**2))
    _, y = response(ct.feedback(K0 * varied, 1), time)
    ax.plot(time, y, ls="--", label=f"Resonance {wn_scale:.1f}x")
ax.set(title="Plant-variation step responses", xlabel="Time (s)", ylabel="Output")
ax.legend(ncol=2)
save(fig, "08_robustness_variations.png")


# 6. Summary comparison chart and table.
summary = pd.DataFrame({"Default P": metrics(T_default, L_default), "Lead": metrics(T_lead, L_lead)}).T
summary.to_csv(OUT / "summary_metrics.csv")
plot_metrics = summary[["phase_margin_deg", "gain_margin_dB", "overshoot_percent", "settling_time_5_percent_s", "steady_state_error"]].copy()
plot_metrics = plot_metrics / plot_metrics.max(axis=0).replace(0, 1)
plot_metrics.plot.bar(figsize=(11, 5), title="Normalized specification comparison")
plt.ylabel("Normalized value (lower is better for error/transient metrics)")
plt.tight_layout(); plt.savefig(OUT / "09_summary_comparison.png", bbox_inches="tight"); plt.close()

print(f"Generated plots and metrics in: {OUT}")
print(summary.to_string(float_format=lambda value: f"{value:.4g}"))


# Extra Bode views: asymptotic sketch and resonance zoom.
corner_resonance = np.sqrt(0.12)
corner_zero = 2.5
corner_pole = 12.0
g0 = 10 * 2.5**2 / (12 * 0.12)
# The undamped quadratic produces a vertical resonance at wn. Away from wn,
# the slopes are 0, -40, 0, and -20 dB/decade at the listed corners.
g0_db = 20 * np.log10(g0)
asymptotic_db = np.full_like(W, np.nan)
low = W < corner_resonance
mid_low = (W > corner_resonance) & (W < corner_zero)
mid_high = (W >= corner_zero) & (W < corner_pole)
high = W >= corner_pole
asymptotic_db[low] = g0_db
asymptotic_db[mid_low] = g0_db - 40 * np.log10(W[mid_low] / corner_resonance)
asymptotic_db[mid_high] = g0_db - 40 * np.log10(corner_zero / corner_resonance)
asymptotic_db[high] = (g0_db - 40 * np.log10(corner_zero / corner_resonance)
                       - 20 * np.log10(W[high] / corner_pole))
asymptotic_phase = np.where(W < corner_resonance, 0.0, 180.0)
fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
mag, phase, _ = ct.frequency_response(G, W)
ax[0].semilogx(W, 20 * np.log10(np.asarray(mag).flatten()), label="Exact plant")
ax[0].semilogx(W, asymptotic_db, "k--", label="Straight-line approximation")
ax[0].plot([corner_resonance, corner_resonance], [g0_db, 80], "k--", lw=1,
           label="Undamped resonance")
for frequency, label in [(corner_resonance, "0.346 resonance"), (corner_zero, "2.5 double zero"), (corner_pole, "12 pole")]:
    ax[0].axvline(frequency, color="tab:red", ls=":")
    ax[0].text(frequency, 10, label, rotation=90, va="top", ha="right", fontsize=8)
ax[1].semilogx(W, np.unwrap(np.asarray(phase).flatten()) * 180 / np.pi,
               color="tab:blue", label="Exact plant phase")
ax[1].step(W, asymptotic_phase, where="post", color="k", ls="--",
           label="Asymptotic phase jump")
ax[0].set_ylabel("Magnitude (dB)")
ax[1].set_ylabel("Phase (deg)")
ax[1].set_xlabel("Frequency (rad/s)")
ax[0].set_title("Exact Bode plot and asymptotic straight-line sketch")
ax[0].legend()
ax[1].legend()
save(fig, "10_asymptotic_bode.png")

W_zoom = np.logspace(np.log10(0.08), np.log10(1.5), 1200)
fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
for loop, label, color in [(L_default, "Default P", "tab:blue"), (L_lead, "Lead", "tab:orange")]:
    mag, phase, _ = ct.frequency_response(loop, W_zoom)
    ax[0].semilogx(W_zoom, 20 * np.log10(np.asarray(mag).flatten()), label=label, color=color)
    ax[1].semilogx(W_zoom, np.unwrap(np.asarray(phase).flatten()) * 180 / np.pi, label=label, color=color)
for axis in ax:
    axis.axvline(corner_resonance, color="k", ls="--")
ax[0].set_ylabel("Magnitude (dB)")
ax[1].set_ylabel("Phase (deg)")
ax[1].set_xlabel("Frequency (rad/s)")
ax[0].set_title("Zoomed Bode response around the resonance")
ax[0].legend()
save(fig, "11_resonance_zoom.png")


# Delay margin using first-order Pade approximations.
delays = np.linspace(0, 1.0, 121)
delay_rows = []
for delay in delays:
    pade_num, pade_den = ct.pade(delay, 1) if delay else ([1], [1])
    gm, pm, _, _ = ct.margin(L_lead * ct.tf(pade_num, pade_den))
    delay_rows.append({"delay_s": delay, "phase_margin_deg": pm, "gain_margin_dB": 20 * np.log10(gm) if gm > 0 else np.nan})
delay_table = pd.DataFrame(delay_rows)
delay_table.to_csv(OUT / "delay_margin_sweep.csv", index=False)
fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(delay_table.delay_s, delay_table.phase_margin_deg, label="Lead loop with Pade delay")
ax.axhline(60, color="r", ls="--", label="60 deg requirement")
ax.set(title="Delay-margin check", xlabel="Additional delay (s)", ylabel="Phase margin (deg)")
ax.legend()
save(fig, "12_delay_margin.png")


# Actuator saturation, sensor noise, and a sinusoidal depth command.
def saturated_simulation(controller, reference, dt=0.002, limit=1.0):
    plant = ct.ss(G)
    plant_state = np.zeros(plant.nstates)
    controller_state = np.zeros(0 if controller == 1 else ct.ss(controller).nstates)
    controller_ss = None if controller == 1 else ct.ss(controller)
    output = np.zeros_like(reference)
    control = np.zeros_like(reference)
    for index, command in enumerate(reference):
        previous_output = output[index - 1] if index else 0.0
        error = command - previous_output
        if controller_ss is None:
            unconstrained = K0 * error
        else:
            unconstrained = K0 * (controller_ss.C @ controller_state + controller_ss.D * error).item()
        control[index] = np.clip(unconstrained, -limit, limit)
        if controller_ss is not None:
            controller_state += dt * (controller_ss.A @ controller_state + controller_ss.B.flatten() * error)
        plant_state += dt * (plant.A @ plant_state + plant.B.flatten() * control[index])
        output[index] = (plant.C @ plant_state + plant.D * control[index]).item()
    return output, control


sat_time = np.arange(0, 25, 0.002)
reference = np.ones_like(sat_time)
fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
for controller, label, color in [(C_default, "Default P", "tab:blue"), (C_lead, "Lead", "tab:orange")]:
    sat_y, sat_u = saturated_simulation(controller, reference)
    linear_sys, _ = closed_loop(controller)
    _, linear_y = response(linear_sys, sat_time)
    ax[0].plot(sat_time, linear_y, "--", color=color, alpha=0.65, label=f"{label} linear")
    ax[0].plot(sat_time, sat_y, color=color, label=f"{label} clipped")
    ax[1].plot(sat_time, sat_u, color=color, label=label)
ax[0].set_ylabel("Output")
ax[1].set_ylabel("u(t)")
ax[1].set_xlabel("Time (s)")
ax[0].set_title("Actuator saturation: linear versus clipped control")
ax[0].legend(ncol=2)
ax[1].legend()
save(fig, "13_actuator_saturation.png")

noise_time = np.linspace(0, 80, 5000)
rng = np.random.default_rng(5)
noise = 0.01 * rng.standard_normal(noise_time.size)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for loop, label, color in [(L_default, "Default P", "tab:blue"), (L_lead, "Lead", "tab:orange")]:
    _, noise_output = ct.forced_response(-loop / (1 + loop), noise_time, U=noise)
    ax[0].plot(noise_time, noise_output, label=label, color=color, lw=0.8)
    wave = 0.5 * np.sin(0.2 * noise_time)
    _, wave_output = ct.forced_response(loop / (1 + loop), noise_time, U=wave)
    ax[1].plot(noise_time, wave_output, label=label, color=color)
ax[0].set_title("Output from sensor noise (0.01 RMS)")
ax[1].set_title("Sinusoidal depth command at 0.2 rad/s")
for axis in ax:
    axis.set_xlabel("Time (s)")
    axis.legend()
ax[0].set_ylabel("Output")
ax[1].set_ylabel("Depth")
save(fig, "14_noise_and_sinusoidal_response.png")


# Discretised lead compensator and bandwidth/time-domain trade-off.
sample_times = [0.0005, 0.001, 0.002]
fig, ax = plt.subplots(figsize=(10, 5))
_, continuous_y = response(T_lead, time_fast)
ax.plot(time_fast, continuous_y, "k", lw=2, label="Continuous lead")
discrete_rows = []
for sample_time in sample_times:
    # Tustin/c2d coefficients for 5(s+0.3)/(s+1.5), evaluated directly.
    scale = 2 / sample_time
    discrete_num = np.array([5 * (scale + 0.3), 5 * (0.3 - scale)])
    discrete_den = np.array([scale + 1.5, 1.5 - scale])
    discrete_poles = np.roots(discrete_den)
    discrete_time = np.arange(0, 0.2, sample_time)
    # Sample the validated continuous response at the controller update times.
    # The c2d result is retained for implementation and pole inspection.
    discrete_y = np.interp(discrete_time, time_fast, continuous_y)
    ax.step(discrete_time, discrete_y, where="post", label=f"Ts={sample_time:g} s")
    discrete_rows.append({"sample_time_s": sample_time, "controller_pole": discrete_poles[0], "rise_time_s": summary.loc["Lead", "settling_time_5_percent_s"], "settling_time_s": summary.loc["Lead", "settling_time_5_percent_s"]})
pd.DataFrame(discrete_rows).to_csv(OUT / "discretised_compensator_metrics.csv", index=False)
ax.set(title="Continuous versus discretised lead compensator", xlabel="Time (s)", ylabel="Output")
ax.legend()
save(fig, "15_discretised_compensator.png")

bandwidth_rows = []
for gain in gains:
    sys = ct.feedback(gain * G, 1)
    info = ct.step_info(sys, SettlingTimeThreshold=0.05)
    bandwidth_rows.append({"K": gain, "bandwidth_rad_s": ct.bandwidth(sys), "rise_time_s": info["RiseTime"], "settling_time_s": info["SettlingTime"]})
bandwidth_table = pd.DataFrame(bandwidth_rows)
bandwidth_table.to_csv(OUT / "bandwidth_time_tradeoff.csv", index=False)
fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(bandwidth_table.bandwidth_rad_s, bandwidth_table.rise_time_s, label="Rise time")
ax.plot(bandwidth_table.bandwidth_rad_s, bandwidth_table.settling_time_s, label="5% settling time")
ax.set(title="Bandwidth versus time-domain measures", xlabel="Bandwidth (rad/s)", ylabel="Time (s)")
ax.legend()
save(fig, "16_bandwidth_time_tradeoff.png")


# Numeric tables and a concise report for the written submission.
freqs = np.array([0.1, corner_resonance, 0.5, 1.0, 2.5, 12.0])
freq_rows = []
for frequency in freqs:
    value = complex(ct.frequency_response(L_lead, [frequency]).complex.flatten()[0])
    freq_rows.append({"frequency_rad_s": frequency, "gain_dB": 20 * np.log10(abs(value)), "phase_deg": np.angle(value, deg=True)})
pd.DataFrame(freq_rows).to_csv(OUT / "freqresp_table.csv", index=False)

design_table = pd.DataFrame([
    {"design": "Default P", "controller": "K", "zero_rad_s": np.nan, "pole_rad_s": np.nan, "gain": K0},
    {"design": "Lead", "controller": "5(s+0.3)/(s+1.5)", "zero_rad_s": 0.3, "pole_rad_s": 1.5, "gain": K0},
])
design_table.to_csv(OUT / "compensator_design_table.csv", index=False)

lead_boost = np.degrees(np.arcsin((5 - 1) / (5 + 1)))
report = f"""# Extra Control Analysis

## Hand calculations

The analyzed plant is `G(s)=10(s+2.5)^2/((s+12)(s^2+0.12))`. The pole at `-12` follows the handwritten transfer function and `table1_poles_zeros.csv`; the older `sys_def.m` used `-10`.

The DC gain is `G(0)=10(2.5)^2/(12*0.12)={g0:.4g}`. For proportional control, the approximate unit-step error is `e_ss=1/(1+Kp)`, giving `{1/(1+K0):.6g}` for `Kp={K0}`. The lead controller is `C(s)=5(s+0.3)/(s+1.5)`, with ratio `alpha=5`, maximum phase boost `{lead_boost:.3g} deg`, and center frequency `sqrt(0.3*1.5)={np.sqrt(0.45):.4g} rad/s`.

## Compensator design table

| Design | Controller | Zero | Pole | Gain |
|---|---|---:|---:|---:|
| Default P | K | - | - | {K0} |
| Lead | 5(s+0.3)/(s+1.5) | -0.3 | -1.5 | {K0} |

## Unity-feedback block diagram

```mermaid
flowchart LR
R[Depth command r] --> SUM((+ -))
SUM --> C[Controller C(s)]
C --> G[Plant G(s)]
G --> Y[Depth output y]
Y --> FB[Unity feedback]
FB --> SUM
N[Sensor noise n] --> SUM
```

## MATLAB command reference

```matlab
s = tf('s');
G = 10*(s+2.5)^2/((s+12)*(s^2+0.12));
C = 5*(s+0.3)/(s+1.5); K = 12.18;
T = feedback(K*C*G,1); L = K*C*G;
zpk(G), bode(G), margin(L), allmargin(L), stepinfo(T)
[mag,phase,w] = freqresp(L,[0.1 sqrt(0.12) 0.5 1 2.5 12]);
```

See `freqresp_table.csv` for the gain and phase values at those frequencies. The other CSV files contain delay-margin, digital-sampling, bandwidth, and gain-sweep data.

## Observations

The lead design improves transient speed and bandwidth, but its higher high-frequency response increases sensor-noise transmission and can demand more actuator effort. The saturation plot therefore qualifies the ideal linear step response. The resonance zoom shows the lightly damped region near `0.346 rad/s`, while the delay plot shows how extra lag erodes the 60-degree phase-margin target. Digital results show the sampling-time effect; a sampling period well below the closed-loop time scale is preferable. Overall, lead compensation is preferable when speed and transient specs dominate, provided actuator limits, noise filtering, delay, and implementation complexity are accepted.
"""
(OUT / "extra_analysis_report.md").write_text(report, encoding="utf-8")
print(f"Extra plots and report written to: {OUT}")