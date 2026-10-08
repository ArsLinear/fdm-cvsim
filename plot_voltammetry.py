import csv
from pathlib import Path

import ftacvsim


def main():
    """Export the simulated CV potential and current density to a CSV file."""
    output_path = Path(__file__).with_name("voltammogram.csv")

    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Potential (V vs Ag/AgCl)", "Current density (A/cm^2)"])
        writer.writerows(
            (float(potential), float(current))
            for potential, current in zip(
                ftacvsim.waveform_array,
                ftacvsim.current,
                strict=True,
            )
        )

    print(f"Saved {len(ftacvsim.current)} data points to {output_path}")


if __name__ == "__main__":
    main()
