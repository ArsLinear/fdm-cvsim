import csv
from pathlib import Path

import matplotlib.pyplot as plt

import ftacvsim


def main():
    # ftacvsim exposes the time array returned by the simulator.
    time = ftacvsim.time_array
    current_density = ftacvsim.current

    output_path = Path(__file__).with_name("ftacvsim.csv")
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow([
            "Time (s)",
            "Potential (V vs Ag/AgCl)",
            "Current density (A/cm^2)",
        ])
        writer.writerows(
            (float(t), float(potential), float(current))
            for t, potential, current in zip(
                time, ftacvsim.waveform_array, current_density, strict=True
            )
        )
    print(f"Saved {len(current_density)} data points to {output_path}")

    fig, ax = plt.subplots()
    ax.plot(time, current_density * 1e3, linewidth=0.8)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Current density (mA/cm²)")
    ax.set_title("Simulated FTACV current-time response")
    ax.axhline(0, color="black", linewidth=0.7)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
