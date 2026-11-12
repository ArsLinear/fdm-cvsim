import csv
from pathlib import Path

import ftacvsim


def main():
    output_path = Path("./result").with_name("voltammogram.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["Potential (V)", "Current density (A/cm^2)"])
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
