from pathlib import Path

import numpy as np
import geopandas as gpd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

RED = "#A01914"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

PATH_TO_FONT_AWESOME = Path(__file__).parent / 'fonts' / 'fa-solid-900.ttf'
LAYER_UNICODE = "\uf5fd"
MONEY_UNICODE = "\uf3d1"


def plot_map(path_to_continental_shape, path_to_national_shapes, path_to_regional_shapes, path_to_plot):
    """Plot maps of results."""
    np.random.seed(123456789)
    fig = plt.figure(figsize=(8, 8), constrained_layout=True)
    axes = fig.subplots(2, 2).flatten()
    norm = matplotlib.colors.Normalize(vmin=0, vmax=0.25)
    cmap = sns.light_palette(sns.desaturate(RED, 0.85), reverse=False, as_cmap=True)
    continental = gpd.read_file(path_to_continental_shape).to_crs(EPSG_3035_PROJ4)
    national = gpd.read_file(path_to_national_shapes).to_crs(EPSG_3035_PROJ4)
    regional = gpd.read_file(path_to_regional_shapes).to_crs(EPSG_3035_PROJ4)
    regional_relaxed = regional.copy()

    continental["cost"] = 0.1
    national["cost"] = np.random.normal(loc=0.14, scale=0.04, size=len(national))
    regional["cost"] = np.random.normal(loc=0.18, scale=0.04, size=len(regional))
    regional_relaxed["cost"] = np.random.normal(loc=0.12, scale=0.04, size=len(regional_relaxed))

    _plot_layer(continental, "continental", norm, cmap, axes[0])
    _plot_layer(national, "national", norm, cmap, axes[1])
    _plot_layer(regional, "regional", norm, cmap, axes[2])
    _plot_layer(regional_relaxed, "regional with trade", norm, cmap, axes[3])

    _plot_colorbar(fig, axes, norm, cmap)
    fig.savefig(path_to_plot, dpi=300, transparent=True)


def _plot_layer(units, layer_name, norm, cmap, ax):
    ax.set_aspect('equal')
    units.plot(
        linewidth=0.1,
        column="cost",
        vmin=norm.vmin,
        vmax=norm.vmax,
        cmap=cmap,
        ax=ax
    )
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)

    ax.annotate(
        f"{LAYER_UNICODE} ",
        xy=[0.10, 0.90],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color="black"
    )
    ax.annotate(
        f"{MONEY_UNICODE} ",
        xy=[0.10, 0.85],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color=sns.desaturate(RED, 0.85)
    )
    ax.annotate(layer_name, xy=[0.17, 0.90], xycoords='axes fraction')
    ax.annotate(f"{units.cost.mean():.2f}€/kWh", xy=[0.17, 0.85], xycoords='axes fraction')


def _plot_colorbar(fig, axes, norm, cmap):
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=axes, fraction=1, aspect=35, shrink=0.65)
    cbar.set_ticks(cbar.get_ticks())
    cbar.set_ticklabels(["{:.2f}€/kWh".format(tick)
                         for tick in cbar.get_ticks()])
    cbar.outline.set_linewidth(0)


if __name__ == "__main__":
    plot_map(
        path_to_continental_shape=snakemake.input.continental_shape,
        path_to_national_shapes=snakemake.input.national_shapes,
        path_to_regional_shapes=snakemake.input.regional_shapes,
        path_to_plot=snakemake.output[0]
    )
