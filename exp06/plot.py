"""Plot the measured temperature readings from the lab notebook."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OPEN_TIME = [0, 0.5, 10, 15, 20, 25, 30, 35, 40]
OPEN = [27.9, 28.1, 28.4, 28.9, 29.6, 30.1, 31.4, 32.5, 33.7]
CLOSED_LOOP_TIME = [0, 0.5, 10, 15, 20, 25, 30]
CLOSED_LOOP = [27.5, 28.1, 29.2, 29.8, 30.6, 32.5, 34.3]
PARTIALLY_OPEN_TIME = [0, 0.5, 10, 15, 20, 25, 30, 35]
PARTIALLY_OPEN = [27.9, 28.0, 28.5, 28.9, 30.2, 31.2, 32.0, 33.5]


def make_plot(time, temperature, title, color, output_path):
	"""Create and save one temperature-vs-time graph."""
	fig, ax = plt.subplots(figsize=(8, 5))
	ax.plot(time, temperature, marker="o", linewidth=2, color=color)
	ax.set_title(title)
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Temperature (deg C)")
	ax.set_xticks(time)
	ax.grid(True, alpha=0.3)
	fig.tight_layout()
	fig.savefig(output_path, dpi=200)
	return fig


def make_tangent_plot(output_path):
	"""Plot a fixed-slope line shifted until the upper and lower areas match."""
	cutoff_temperature = 33.0
	time = np.asarray(PARTIALLY_OPEN_TIME, dtype=float)
	temperature = np.asarray(PARTIALLY_OPEN, dtype=float)
	coefficients = np.polyfit(time, temperature, 3)
	inflection_candidates = np.roots(np.polyder(coefficients, 2))
	inflection_time = next(
		float(root.real)
		for root in inflection_candidates
		if abs(root.imag) < 1e-8 and time.min() <= root.real <= time.max()
	)
	inflection_temperature = float(np.polyval(coefficients, inflection_time))
	tangent_slope = float(np.polyval(np.polyder(coefficients), inflection_time))

	# Shift the line vertically so the signed area between it and the curve is zero.
	curve_integral = np.polyval(np.polyint(coefficients), time.max()) - np.polyval(
		np.polyint(coefficients), time.min()
	)
	time_integral = (time.max() ** 2 - time.min() ** 2) / 2
	tangent_intercept = (curve_integral - tangent_slope * time_integral) / (
		time.max() - time.min()
	)

	curve_time = np.linspace(time.min(), time.max(), 300)
	curve_temperature = np.polyval(coefficients, curve_time)
	tangent_temperature = tangent_slope * curve_time + tangent_intercept
	difference = curve_temperature - tangent_temperature
	below_cutoff = (curve_temperature <= cutoff_temperature) & (
		tangent_temperature <= cutoff_temperature
	)
	upper_difference = np.where(below_cutoff, np.maximum(difference, 0), 0)
	lower_difference = np.where(below_cutoff, np.maximum(-difference, 0), 0)
	upper_area = np.trapz(upper_difference, curve_time)
	lower_area = np.trapz(lower_difference, curve_time)
	cutoff_roots = np.roots(coefficients - np.array([0, 0, 0, cutoff_temperature]))
	cutoff_time = next(
		float(root.real)
		for root in cutoff_roots
		if abs(root.imag) < 1e-8 and time.min() <= root.real <= time.max()
	)

	fig, ax = plt.subplots(figsize=(8, 5))
	ax.plot(time, temperature, "o", color="tab:orange", label="Measured data")
	ax.plot(curve_time, curve_temperature, color="tab:orange", alpha=0.7, label="Cubic fit")
	ax.plot(
		curve_time,
		tangent_temperature,
		"--",
		color="tab:blue",
		label=f"Shifted tangent, slope = {tangent_slope:.4f} deg C/s",
	)
	ax.fill_between(curve_time, curve_temperature, tangent_temperature,
		where=below_cutoff & (difference >= 0),
		color="tab:green", alpha=0.25, label=f"Upper area = {upper_area:.2f}")
	ax.fill_between(curve_time, curve_temperature, tangent_temperature,
		where=below_cutoff & (difference < 0),
		color="tab:red", alpha=0.2, label=f"Lower area = {lower_area:.2f}")
	ax.axhline(cutoff_temperature, color="magenta", linestyle="-.",
		label=f"Top cutoff = {cutoff_temperature:.0f} deg C")
	ax.scatter([cutoff_time], [cutoff_temperature], color="magenta", zorder=4,
		label=f"Cutoff reached at {cutoff_time:.2f} s")
	ax.scatter([inflection_time], [inflection_temperature], color="black", zorder=3,
		label=f"Slope reference ({inflection_time:.2f} s, {inflection_temperature:.2f} deg C)")
	ax.set_title("Partially Open Curve: Equal-Area Region Below 33 deg C")
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Temperature (deg C)")
	ax.set_xticks(PARTIALLY_OPEN_TIME)
	ax.grid(True, alpha=0.3)
	ax.legend()
	fig.tight_layout()
	fig.savefig(output_path, dpi=200)

	print(f"Slope reference: t = {inflection_time:.2f} s, T = {inflection_temperature:.2f} deg C")
	print(f"Shifted tangent equation: T = {tangent_slope:.4f} t + {tangent_intercept:.4f}")
	print(f"Upper area = {upper_area:.4f}, lower area = {lower_area:.4f}")
	print(f"33 deg C cutoff reached at t = {cutoff_time:.2f} s")


def main():
	output_dir = Path(__file__).parent / "plots"
	output_dir.mkdir(exist_ok=True)

	make_plot(
		OPEN_TIME,
		OPEN,
		"Temperature Response - Open",
		"tab:green",
		output_dir / "open_temperature.png",
	)
	make_plot(
		CLOSED_LOOP_TIME,
		CLOSED_LOOP,
		"Temperature Response - Closed Loop",
		"tab:blue",
		output_dir / "closed_loop_temperature.png",
	)
	make_plot(
		PARTIALLY_OPEN_TIME,
		PARTIALLY_OPEN,
		"Temperature Response - Partially Open",
		"tab:orange",
		output_dir / "partially_open_temperature.png",
	)

	# Show both readings together for direct comparison when run interactively.
	fig, ax = plt.subplots(figsize=(8, 5))
	ax.plot(OPEN_TIME, OPEN, "o-", label="Open", color="tab:green")
	ax.plot(CLOSED_LOOP_TIME, CLOSED_LOOP, "o-", label="Closed loop", color="tab:blue")
	ax.plot(PARTIALLY_OPEN_TIME, PARTIALLY_OPEN, "o-", label="Partially open", color="tab:orange")
	ax.set_title("Temperature Comparison (SP = 30 deg C)")
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Temperature (deg C)")
	ax.grid(True, alpha=0.3)
	ax.legend()
	fig.tight_layout()
	fig.savefig(output_dir / "temperature_comparison.png", dpi=200)
	make_tangent_plot(output_dir / "partially_open_tangent.png")


if __name__ == "__main__":
	main()
