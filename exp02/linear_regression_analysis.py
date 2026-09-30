from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

base_dir = Path(__file__).resolve().parent
artifacts_dir = base_dir / "artifacts"
assets_dir = base_dir / "assets"
assets_dir.mkdir(exist_ok=True)


def load_experiment_data():
    mode_files = {
        "Manual CW": [
            artifacts_dir / "01_manual_mode_CW.csv",
        ],
        "Auto CW": [
            artifacts_dir / "02_auto_mode_CW.csv",
        ],
        "Auto ACW": [
            artifacts_dir / "03_auto_mode_ACW.csv",
        ],
    }

    mode_frames = {}
    for mode_name, file_list in mode_files.items():
        frames = []
        for file_path in file_list:
            df = pd.read_csv(file_path)
            df["Source"] = file_path.stem
            frames.append(df)

        combined = pd.concat(frames, ignore_index=True)
        combined["Frequency_Hz"] = pd.to_numeric(combined["Frequency_Hz"], errors="coerce")
        combined["Inverter_Voltage_V"] = pd.to_numeric(combined["Inverter_Voltage_V"], errors="coerce")
        combined["Motor_Speed_RPM"] = pd.to_numeric(combined["Motor_Speed_RPM"], errors="coerce")

        combined = combined.dropna(
            subset=["Frequency_Hz", "Inverter_Voltage_V", "Motor_Speed_RPM"]
        ).reset_index(drop=True)

        combined["Mode"] = mode_name
        mode_frames[mode_name] = combined

    return mode_frames


def linear_regression_summary(x, y, x_name, y_name, y_units):
    slope, intercept = np.polyfit(x, y, 1)
    predicted = intercept + slope * x
    residuals = y - predicted

    mae = float(np.mean(np.abs(residuals)))
    rmse = float(np.sqrt(np.mean(residuals ** 2)))
    r2 = float(1 - np.sum(residuals ** 2) / np.sum((y - np.mean(y)) ** 2))
    max_abs_dev = float(np.max(np.abs(residuals)))
    std_dev = float(np.std(residuals, ddof=1)) if len(residuals) > 1 else 0.0
    mean_bias = float(np.mean(residuals))

    stats = {
        "x_variable": x_name,
        "y_variable": y_name,
        "units": y_units,
        "slope": slope,
        "intercept": intercept,
        "r_squared": r2,
        "mae": mae,
        "rmse": rmse,
        "std_residual": std_dev,
        "mean_bias": mean_bias,
        "max_abs_deviation": max_abs_dev,
        "n_points": len(x),
    }

    summary_df = pd.DataFrame(
        [
            {
                "x_variable": x_name,
                "y_variable": y_name,
                "slope": slope,
                "intercept": intercept,
                "r_squared": r2,
                "mae": mae,
                "rmse": rmse,
                "std_residual": std_dev,
                "mean_bias": mean_bias,
                "max_abs_deviation": max_abs_dev,
                "n_points": len(x),
            }
        ]
    )
    return stats, summary_df, predicted, residuals


def make_figure(x, y, x_label, y_label, fit_label, y_units, title, output_name, plot_color="tab:blue"):
    slope, intercept = np.polyfit(x, y, 1)
    predicted = intercept + slope * x
    residuals = y - predicted
    abs_dev = np.abs(residuals)
    mae = float(np.mean(abs_dev))
    rmse = float(np.sqrt(np.mean(residuals ** 2)))
    r2 = float(1 - np.sum(residuals ** 2) / np.sum((y - np.mean(y)) ** 2))
    max_abs_dev = float(np.max(abs_dev))
    std_dev = float(np.std(residuals, ddof=1)) if len(residuals) > 1 else 0.0
    mean_bias = float(np.mean(residuals))

    sorted_idx = np.argsort(x)
    x_sorted = x[sorted_idx]
    y_sorted = predicted[sorted_idx]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)

    ax1.plot(x_sorted, y_sorted, color="tab:red", linewidth=2.5, linestyle="--", label=fit_label)
    ax1.scatter(x, y, s=35, color=plot_color, alpha=0.8, label="Measured data")
    ax1.set_title(title, fontsize=18, fontweight="bold")
    ax1.set_xlabel(x_label, fontsize=14)
    ax1.set_ylabel(y_label, fontsize=14)
    ax1.grid(True, linestyle="--", linewidth=0.8, alpha=0.5)
    ax1.tick_params(labelsize=12)
    ax1.legend(loc="best", frameon=True, fontsize=10)
    ax1.text(
        0.02,
        0.96,
        f"Fit: y = {slope:.3f}x + {intercept:.3f}\nR² = {r2:.4f}\nMAE = {mae:.4f} {y_units}\nRMSE = {rmse:.4f} {y_units}",
        transform=ax1.transAxes,
        va="top",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.7", alpha=0.9),
        fontsize=10,
    )

    dev_max = max(1.1 * max_abs_dev, 1e-6)
    ax2.axhline(0, color="black", linewidth=1.0, linestyle="--", alpha=0.7)
    ax2.plot(x, residuals, color=plot_color, linewidth=2.0, marker="o", markersize=6, label="Deviation")
    ax2.fill_between(x, residuals, 0, where=(residuals >= 0), color=plot_color, alpha=0.15)
    ax2.fill_between(x, residuals, 0, where=(residuals < 0), color="tab:red", alpha=0.12)
    ax2.set_title(f"Deviation from Linear Approximation", fontsize=18, fontweight="bold")
    ax2.set_xlabel(x_label, fontsize=14)
    ax2.set_ylabel(f"Deviation ({y_units})", fontsize=14)
    ax2.grid(True, linestyle="--", linewidth=0.8, alpha=0.5)
    ax2.tick_params(labelsize=12)
    ax2.set_ylim(-dev_max, dev_max)
    ax2.legend(loc="best", frameon=True, fontsize=10)
    ax2.text(
        0.02,
        0.96,
        f"Std residual = {std_dev:.4f} {y_units}\nMean bias = {mean_bias:.4f} {y_units}\nMax abs deviation = {max_abs_dev:.4f} {y_units}",
        transform=ax2.transAxes,
        va="top",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.7", alpha=0.9),
        fontsize=10,
    )

    fig.savefig(assets_dir / output_name, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r2,
        "mae": mae,
        "rmse": rmse,
        "std_residual": std_dev,
        "mean_bias": mean_bias,
        "max_abs_deviation": max_abs_dev,
        "n_points": len(x),
    }


def main():
    mode_data = load_experiment_data()
    summary_rows = []
    mode_plot_map = {
        "Manual CW": ("12_frequency_voltage_regression.png", "13_frequency_speed_regression.png"),
        "Auto CW": ("14_auto_cw_frequency_voltage_regression.png", "15_auto_cw_frequency_speed_regression.png"),
        "Auto ACW": ("16_auto_acw_frequency_voltage_regression.png", "17_auto_acw_frequency_speed_regression.png"),
    }

    for mode_name, data in mode_data.items():
        voltage_stats, voltage_summary, _, _ = linear_regression_summary(
            data["Frequency_Hz"].to_numpy(),
            data["Inverter_Voltage_V"].to_numpy(),
            "Frequency",
            "Inverter Voltage",
            "V",
        )

        speed_stats, speed_summary, _, _ = linear_regression_summary(
            data["Frequency_Hz"].to_numpy(),
            data["Motor_Speed_RPM"].to_numpy(),
            "Frequency",
            "Motor Speed",
            "RPM",
        )

        voltage_summary["Mode"] = mode_name
        speed_summary["Mode"] = mode_name
        summary_rows.extend([voltage_summary.iloc[0].to_dict(), speed_summary.iloc[0].to_dict()])

        voltage_fit = make_figure(
            data["Frequency_Hz"].to_numpy(),
            data["Inverter_Voltage_V"].to_numpy(),
            "Frequency (Hz)",
            "Inverter Voltage (V)",
            f"Voltage fit: V = {voltage_stats['slope']:.3f}*f + {voltage_stats['intercept']:.3f}",
            "V",
            f"{mode_name} - Frequency vs Voltage Regression",
            mode_plot_map[mode_name][0],
            plot_color="tab:blue",
        )

        speed_fit = make_figure(
            data["Frequency_Hz"].to_numpy(),
            data["Motor_Speed_RPM"].to_numpy(),
            "Frequency (Hz)",
            "Motor Speed (RPM)",
            f"Speed fit: RPM = {speed_stats['slope']:.3f}*f + {speed_stats['intercept']:.3f}",
            "RPM",
            f"{mode_name} - Frequency vs Speed Regression",
            mode_plot_map[mode_name][1],
            plot_color="tab:green",
        )

        print(f"\n{mode_name}")
        print("-" * 70)
        print("Voltage regression parameters:")
        print(voltage_fit)
        print("Speed regression parameters:")
        print(speed_fit)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(artifacts_dir / "linear_regression_summary.csv", index=False)

    print("\nLinear Regression Analysis for Experiment 02")
    print("-" * 70)
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
