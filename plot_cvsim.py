"""Run cvsim.py's example simulation and plot the cyclic voltammogram."""

import matplotlib.pyplot as plt

import cvsim


def main():
    # Importing cvsim runs its example simulation and exposes these arrays.
    potential = cvsim.waveform
    current_density = cvsim.current

    fig, ax = plt.subplots()
    ax.plot(potential, current_density * 1e3)
    ax.set_xlabel("Potential (V vs Ag/AgCl)")
    ax.set_ylabel("Current density (mA/cm²)")
    ax.set_title("Simulated cyclic voltammogram")
    ax.axhline(0, color="black", linewidth=0.7)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
