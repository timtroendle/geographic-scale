from pathlib import Path

import calliope
import numpy as np
import geopandas as gpd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

RED = "#A01914"
BLUE = "#4F6DB8"
CMAP = 'RdBu_r'
NORM = matplotlib.colors.Normalize(vmin=-1, vmax=3)
PANEL_FONT_SIZE = 10
PANEL_FONT_WEIGHT = "bold"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

PATH_TO_FONT_AWESOME = Path(__file__).parent.parent / 'fonts' / 'fa-solid-900.ttf'
LAYER_UNICODE = "\uf5fd"
MONEY_UNICODE = "\uf3d1"

ALL_TECHS = ['wind_onshore_monopoly', 'demand_elec', 'biofuel',
             'roof_mounted_pv', 'hydro_reservoir', 'wind_offshore',
             'wind_onshore_competing', 'battery', 'hydrogen', 'free_transmission',
             'ac_transmission', 'open_field_pv', 'hydro_run_of_river',
             'pumped_hydro']
SUPPLY_TECHS = ['wind_onshore_monopoly', 'biofuel',
                'roof_mounted_pv', 'hydro_reservoir', 'wind_offshore',
                'wind_onshore_competing', 'open_field_pv', 'hydro_run_of_river']


def plot_map(path_to_continental_shape, path_to_national_shapes, path_to_regional_shapes,
             path_to_continental_result, path_to_national_result, path_to_regional_result,
             path_to_shapes, path_to_plot):
    """Plot maps of results."""
    fig = plt.figure(figsize=(8, 4), constrained_layout=True)
    axes = fig.subplots(1, 2).flatten()
    shapes = (gpd.read_file(path_to_shapes)
                 .to_crs(EPSG_3035_PROJ4)
                 .rename(columns={"id": "locs"})
                 .set_index("locs")
                 .rename(index=lambda idx: idx.replace(".", "-")))

    continental_model = calliope.read_netcdf(path_to_continental_result)
    national_model = calliope.read_netcdf(path_to_national_result)
    regional_model = calliope.read_netcdf(path_to_regional_result)

    continental = gpd.read_file(path_to_continental_shape).to_crs(EPSG_3035_PROJ4).set_index("id")
    national = gpd.read_file(path_to_national_shapes).to_crs(EPSG_3035_PROJ4).set_index("country_code")
    regional = (gpd.read_file(path_to_regional_shapes)
                   .to_crs(EPSG_3035_PROJ4)
                   .set_index("id")
                   .rename(index=lambda idx: idx.replace(".", "-")))

    continental["cost"] = _read_cost(shapes, continental_model, "continental")
    national["cost"] = _read_cost(shapes, national_model, "national")
    regional["cost"] = _read_cost(shapes, regional_model, "regional")

    base_costs = continental["cost"].iloc[0]
    national["cost"] = national["cost"] / base_costs
    regional["cost"] = regional["cost"] / base_costs

    total_costs_national = _read_total_system_costs(national_model) / base_costs
    total_costs_regional = _read_total_system_costs(regional_model) / base_costs

    _plot_layer(national, total_costs_national, "National", "a", axes[0])
    _plot_layer(regional, total_costs_regional, "Regional", "b", axes[1])

    _plot_colorbar(fig, axes)
    fig.savefig(path_to_plot, dpi=600, transparent=False)


def _plot_layer(units, total_cost, layer_name, panel_id, ax):
    ax.set_aspect('equal')
    units.plot(
        linewidth=0.1,
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
        f"{LAYER_UNICODE} ",
        xy=[0.19, 0.95],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color="black"
    )
    ax.annotate(
        f"{MONEY_UNICODE} ",
        xy=[0.19, 0.90],
        xycoords='axes fraction',
        fontproperties=matplotlib.font_manager.FontProperties(fname=PATH_TO_FONT_AWESOME.as_posix()),
        color=sns.desaturate(RED, 0.85)
    )
    ax.annotate(layer_name, xy=[0.26, 0.95], xycoords='axes fraction')
    ax.annotate(f"{total_cost:.1f}", xy=[0.26, 0.90], xycoords='axes fraction')
    ax.annotate(panel_id, xy=[0.10, 0.95], xycoords='axes fraction',
                fontsize=PANEL_FONT_SIZE, weight=PANEL_FONT_WEIGHT)


def _plot_colorbar(fig, axes):
    s_m = matplotlib.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    cmap = s_m.get_cmap()
    colors = cmap(np.linspace(1 / 4, 1, cmap.N))
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('cut_jet', colors)
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(vmin=0, vmax=3))
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=axes, fraction=1, aspect=35, shrink=0.65)
    cbar.set_ticks([0, 1, 2, 3])
    cbar.set_ticklabels(["0", "1", "2", "≥3"])
    cbar.outline.set_linewidth(0)
    cbar.ax.set_ylabel('Relative system cost', rotation=90)


def _read_cost(shapes, model, level):
    assert level in ["continental", "national", "regional"]
    assert set(model.inputs.techs.values) == set(ALL_TECHS)
    cost = model.get_formatted_array('cost').squeeze('costs').sum("techs")
    carrier_prod = (model
                    .get_formatted_array('carrier_prod')
                    .squeeze(['carriers']))
    carrier_prod = (carrier_prod
                    .sel(techs=SUPPLY_TECHS)
                    .sum(dim=['techs', 'timesteps']))
    if level == "continental":
        cost = cost.sum(dim="locs").item()
        carrier_prod = carrier_prod.sum(dim="locs").item()
    elif level == "national":
        cost = (cost.groupby(shapes.country_code.to_xarray().sortby(cost.locs))
                    .sum("locs")
                    .to_series())
        carrier_prod = (carrier_prod.groupby(shapes.country_code.to_xarray().sortby(carrier_prod.locs))
                                    .sum("locs")
                                    .to_series())
    elif level == "regional":
        cost = cost.to_series()
        carrier_prod.to_series()
    return (cost / carrier_prod)


def _read_total_system_costs(model):
    return model.get_formatted_array("total_levelised_cost").item()


if __name__ == "__main__":
    plot_map(
        path_to_shapes=snakemake.input.shapes,
        path_to_continental_shape=snakemake.input.continental_shape,
        path_to_national_shapes=snakemake.input.national_shapes,
        path_to_regional_shapes=snakemake.input.regional_shapes,
        path_to_continental_result=snakemake.input.continental_result,
        path_to_national_result=snakemake.input.national_result,
        path_to_regional_result=snakemake.input.regional_result,
        path_to_plot=snakemake.output[0]
    )
