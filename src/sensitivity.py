import matplotlib.pyplot as plt
import numpy as np


def plot_sensitivity(path_to_output):
    battery, transmission = np.meshgrid(np.linspace(0, 500, 100), np.linspace(0, 0.250, 100))
    z = 0.01 + (battery * 0.0008 + transmission * 0.4 + battery * transmission * 0.0006) * 0.25
    z = z[:-1, :-1]
    z_min, z_max = 0.01, np.abs(z).max()

    fig, ax = plt.subplots()
    c = ax.pcolormesh(battery, transmission, z, cmap='viridis', vmin=z_min, vmax=z_max)
    ax.set_xlabel('battery costs [€/kW]')
    ax.set_ylabel('transmission grid costs [€/kW/km]')
    # set the limits of the plot to the limits of the data
    ax.axis([battery.min(), battery.max(), transmission.min(), transmission.max()])
    fig.colorbar(c, ax=ax)

    fig.savefig(path_to_output, dpi=300, transparent=True)


if __name__ == "__main__":
    plot_sensitivity(path_to_output=snakemake.output[0])
