import matplotlib.pyplot as plt

import ftacvsim


def main():
    # ftacvsim exposes the time array returned by the simulator.
    time = ftacvsim.time_array
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
