import matplotlib.pyplot as plt
import numpy as np

import ftacvsim


def main():
    # ftacvsim runs its example simulation on import and exposes the total
    # current density as `current`.
    time = np.arange(len(ftacvsim.current)) * ftacvsim.SAMPLE_TIME
    current_density = ftacvsim.current

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
