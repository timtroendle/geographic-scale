from itertools import chain

import numpy as np
import xarray as xr
import geopandas as gpd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

RED = "#A01914"
BLUE = "#4F6DB8"
RED_TO_BLUE = [ # from https://gka.github.io using lightness correction
    '#002d6e', '#375aa2', '#6f8ad1', '#a7bffa',
    '#f5f5f5', '#fdad97', '#e36b55', '#b23125', '#720000'
]
CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("signature-BlRd", RED_TO_BLUE)
NORM = matplotlib.colors.Normalize(vmin=-1, vmax=3)
LAYER_FONT_WEIGHT = "medium"
EDGE_WIDTH = 0.06
EDGE_COLOR = "white"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

CONTINENTAL_SCENARIO = "continental-autarky-100-continental-grid"
NATIONAL_SCENARIO = "national-autarky-100-national-grid"
REGIONAL_SCENARIO = "regional-autarky-100-regional-grid"
DEMAND_TECH = "demand_elec"


def plot_map(path_to_national_shapes, path_to_regional_shapes, connected_regions,
             path_to_aggregated_results, path_to_shapes, path_to_plot):
    """Plot maps of results."""
    fig = plt.figure(figsize=(4.41, 2.2), constrained_layout=True)
    axes = fig.subplots(1, 3, gridspec_kw={'width_ratios': [30, 30, 1]}).flatten()
    shapes = (gpd.read_file(path_to_shapes)
                 .to_crs(EPSG_3035_PROJ4)
                 .rename(columns={"id": "locs"})
                 .set_index("locs")
                 .rename(index=lambda idx: idx.replace(".", "-")))

    national = gpd.read_file(path_to_national_shapes).to_crs(EPSG_3035_PROJ4).set_index("country_code")
    regional = (gpd.read_file(path_to_regional_shapes)
                   .to_crs(EPSG_3035_PROJ4)
                   .set_index("id")
                   .rename(index=lambda idx: idx.replace(".", "-")))

    cost_reader = CostReader(
        path_to_aggregated_results=path_to_aggregated_results,
        connected_regions=connected_regions,
        country_codes=shapes.country_code.to_xarray()
    )
    national["cost"] = cost_reader.read(scenario=NATIONAL_SCENARIO, level="national")
    regional["cost"] = cost_reader.read(scenario=REGIONAL_SCENARIO, level="regional")

    base_costs = cost_reader.read(scenario=CONTINENTAL_SCENARIO, level="continental")
    national["cost"] = national["cost"] / base_costs
    regional["cost"] = regional["cost"] / base_costs

    total_costs_national = cost_reader.read(scenario=NATIONAL_SCENARIO, level="continental") / base_costs
    total_costs_regional = cost_reader.read(scenario=REGIONAL_SCENARIO, level="continental") / base_costs

    _plot_layer(national, total_costs_national, "National supply and balancing", "A", axes[0])
    _plot_layer(regional, total_costs_regional, "Regional supply and balancing", "B", axes[1])

    _plot_colorbar(fig, axes[2].inset_axes([0, 0.175, 1, 0.65]))
    axes[2].axis("off")

    fig.savefig(path_to_plot, pil_kwargs={"compression": "tiff_lzw"})


def _plot_layer(units, total_cost, layer_name, panel_id, ax):
    ax.set_aspect('equal')
    units.plot(
        linewidth=EDGE_WIDTH,
        edgecolor=EDGE_COLOR,
        column="cost",
        vmin=NORM.vmin,
        vmax=NORM.vmax,
        cmap=CMAP,
        ax=ax
    )
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, top=True, bottom=True, left=True, right=True)

    ax.annotate(
        layer_name,
        xy=[0.5, -0.05],
        xycoords='axes fraction',
        ha='center',
        va='center',
        weight=LAYER_FONT_WEIGHT
    )
    ax.set_title(panel_id, loc="left")


def _plot_colorbar(fig, axis):
    s_m = matplotlib.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    cmap = s_m.get_cmap()
    colors = cmap(np.linspace(1 / 4, 1, cmap.N))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cut_jet', colors)
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(vmin=0, vmax=3))
    s_m.set_array([])
    cbar = fig.colorbar(s_m, cax=axis)
    cbar.set_ticks([0, 1, 2, 3])
    cbar.set_ticklabels(["0", "1", "2", "≥3"])
    cbar.outline.set_linewidth(0)
    cbar.ax.set_ylabel('Relative system cost', rotation=90)


class CostReader:

    def __init__(self, path_to_aggregated_results, country_codes, connected_regions):
        self.__results = xr.open_dataset(path_to_aggregated_results)
        self.__country_codes = country_codes
        self.__connected_regions = connected_regions

    def read(self, scenario, level):
        assert level in ["continental", "national", "regional"]
        cost = self.__results.cost.sel(scenario=scenario).sum("techs")
        demand = -self.__results.energy_cap.sel(techs=DEMAND_TECH, scenario=scenario)
        if level == "continental":
            cost = cost.sum(dim="locs").item()
            demand = demand.sum(dim="locs").item()
        elif level == "national":
            cost = (
                cost
                .groupby(self.__country_codes.sortby(cost.locs))
                .sum("locs")
                .to_series()
            )
            demand = (
                demand
                .groupby(self.__country_codes.sortby(demand.locs))
                .sum("locs")
                .to_series()
            )
        elif level == "regional":
            cost = self._postprocess_combined_regions(cost).to_series()
            demand = self._postprocess_combined_regions(demand).to_series()
        return (cost / demand)

    def _postprocess_combined_regions(self, da):
        region_pairs = list(self.__connected_regions.items())
        regions = list(chain(*region_pairs))
        if all([region in da.locs for region in regions]):
            # config is valid, and resolution is regional. apply config
            for region1, region2 in region_pairs:
                sum_data = da.sel(locs=region1) + da.sel(locs=region2)
                da.loc[dict(locs=region1)] = sum_data
                da.loc[dict(locs=region2)] = sum_data
            return da
        elif all([region not in da.locs for region in regions]):
            # config is valid, but resolution is not regional. do nothing
            return da
        else:
            raise ValueError("Config of connected regions is invalid.")


if __name__ == "__main__":
    plot_map(
        path_to_shapes=snakemake.input.shapes,
        path_to_national_shapes=snakemake.input.national_shapes,
        path_to_regional_shapes=snakemake.input.regional_shapes,
        path_to_aggregated_results=snakemake.input.results,
        connected_regions=snakemake.params.connected_regions,
        path_to_plot=snakemake.output[0]
    )
