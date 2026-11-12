import csv
from pathlib import Path
import matplotlib.pyplot as plt
import ftacvsim


def main():
    time = ftacvsim.time_array
    current_density = ftacvsim.current

    output_path = Path("./result").with_name("ftacvsim.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow([
            "Time (s)",
            "Potential (V)",
            "Current density (A/cm^2)",
        ])
        writer.writerows(
            (float(t), float(potential), float(current))
            for t, potential, current in zip(
                time, ftacvsim.waveform_array, current_density, strict=True
            )
        )
    print(f"Saved {len(current_density)} data points to {output_path}")


if __name__ == "__main__":
    main()
