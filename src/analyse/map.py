from pathlib import Path

import pandas as pd
import geopandas as gpd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import networkx as nx
import shapely

RED = "#A01914"
BLUE = "#4F6DB8"

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "

MAP_MIN_X = 2500000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6600000
MAP_MAX_Y = 5500000

PATH_TO_FONT_AWESOME = Path(__file__).parent.parent / 'fonts' / 'fa-solid-900.ttf'
LAYER_UNICODE = "\uf5fd"
MONEY_UNICODE = "\uf3d1"


def plot_map(path_to_continental_shape, path_to_national_shapes, path_to_regional_shapes,
             path_to_continental_result, path_to_national_result, path_to_regional_result,
             path_to_shapes, path_to_plot, scaling_factor_cost):
    """Plot maps of results."""
    fig = plt.figure(figsize=(8, 8), constrained_layout=True)
    axes = fig.subplots(2, 2).flatten()
    norm = matplotlib.colors.Normalize(vmin=0, vmax=6)
    cmap = sns.light_palette(sns.desaturate(RED, 0.85), reverse=False, as_cmap=True)
    shapes = (gpd.read_file(path_to_shapes)
                 .to_crs(EPSG_3035_PROJ4)
                 .rename(columns={"id": "locs"})
                 .set_index("locs")
                 .rename(index=lambda idx: idx.replace(".", "-")))
    continental = gpd.read_file(path_to_continental_shape).to_crs(EPSG_3035_PROJ4).set_index("id")
    national = gpd.read_file(path_to_national_shapes).to_crs(EPSG_3035_PROJ4).set_index("country_code")
    regional = (gpd.read_file(path_to_regional_shapes)
                   .to_crs(EPSG_3035_PROJ4)
                   .set_index("id")
                   .rename(index=lambda idx: idx.replace(".", "-")))

    continental["cost"] = _read_cost(shapes, path_to_continental_result, scaling_factor_cost, "continental")
    national["cost"] = _read_cost(shapes, path_to_national_result, scaling_factor_cost, "national")
    regional["cost"] = _read_cost(shapes, path_to_regional_result, scaling_factor_cost, "regional")

    base_costs = continental["cost"].iloc[0]
    continental["cost"] = continental["cost"] / base_costs
    national["cost"] = national["cost"] / base_costs
    regional["cost"] = regional["cost"] / base_costs

    total_costs_continental = _read_total_system_costs(path_to_continental_result, scaling_factor_cost) / base_costs
    total_costs_national = _read_total_system_costs(path_to_national_result, scaling_factor_cost) / base_costs
    total_costs_regional = _read_total_system_costs(path_to_regional_result, scaling_factor_cost) / base_costs

    network_continental = _read_network_graph(shapes, path_to_continental_result)
    network_national = _read_network_graph(shapes, path_to_national_result)
    network_regional = _read_network_graph(shapes, path_to_regional_result)

    _plot_layer(continental, network_continental, total_costs_continental, "continental", norm, cmap, axes[0])
    _plot_layer(national, network_national, total_costs_national, "national", norm, cmap, axes[2])
    _plot_layer(regional, network_regional, total_costs_regional, "regional", norm, cmap, axes[3])
    sns.despine(ax=axes[1], top=True, bottom=True, left=True, right=True)
    axes[1].set_xticks([])
    axes[1].set_yticks([])

    _plot_colorbar(fig, axes, norm, cmap)
    fig.savefig(path_to_plot, dpi=300, transparent=True)


def _plot_layer(units, network_graph, total_cost, layer_name, norm, cmap, ax):
    ax.set_aspect('equal')
    units.plot(
        linewidth=0.1,
        column="cost",
        vmin=norm.vmin,
        vmax=norm.vmax,
        cmap=cmap,
        ax=ax
    )
    gpd.GeoSeries([
        network_graph.edges[region, neighbour]["line"]
        for region, neighbour in network_graph.edges
    ], crs=units.crs).plot(ax=ax, linewidth=0.25, color=BLUE)
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
    ax.annotate(f"{total_cost:.1f}", xy=[0.17, 0.85], xycoords='axes fraction')


def _plot_colorbar(fig, axes, norm, cmap):
    s_m = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    cbar = fig.colorbar(s_m, ax=axes, fraction=1, aspect=35, shrink=0.65)
    cbar.set_ticks(cbar.get_ticks())
    cbar.set_ticklabels(["{:.0f}".format(tick)
                         for tick in cbar.get_ticks()])
    cbar.outline.set_linewidth(0)


def _read_cost(shapes, path_to_result, scaling_factor_cost, level):
    assert level in ["continental", "national", "regional"]
    results = xr.open_dataset(path_to_result)
    cost = results['cost']
    cost = split_loc_techs(cost).sum(dim=['techs']).squeeze('costs')
    carrier_prod = results['carrier_prod']
    supply_only_carrier_prod = carrier_prod.sel(
        loc_tech_carriers_prod=list(results.loc_tech_carriers_supply_all.values)
    )
    carrier_prod = (
        split_loc_techs(supply_only_carrier_prod)
        .sum(dim=['timesteps', 'techs'])
        .squeeze(['carriers'])
    )
    if level == "continental":
        cost = cost.sum(dim="locs").item()
        carrier_prod = carrier_prod.sum(dim="locs").item()
    elif level == "national":
        cost = (cost.groupby(shapes.country_code.to_xarray().sortby(cost.locs))
                    .sum(xr.ALL_DIMS)
                    .to_series())
        carrier_prod = (carrier_prod.groupby(shapes.country_code.to_xarray().sortby(carrier_prod.locs))
                                    .sum(xr.ALL_DIMS)
                                    .to_series())
    elif level == "regional":
        cost = cost.to_series()
        carrier_prod.to_series()
    return (cost / carrier_prod) * scaling_factor_cost / 1e1 # scale from €/MWh to €ct/kWh


def _read_total_system_costs(path_to_results, scaling_factor_cost):
    return (xr.open_dataset(path_to_results)["total_levelised_cost"]
              .squeeze(["costs", "carriers"])
              .item()) * scaling_factor_cost / 1e1 # scale from €/MWh to €ct/kWh


def _read_network_graph(shapes, path_to_results):
    g = nx.Graph()
    g.add_nodes_from((k, {"centroid": v}) for (k, v) in shapes.centroid.to_dict().items())

    model = xr.open_dataset(path_to_results)
    try:
        lines = [(line.split(":")[0], line.split(":")[-1]) for line in model.loc_techs_transmission.values]
        lines = [(a, b, {"line": shapely.geometry.LineString([g.nodes[a]["centroid"], g.nodes[b]["centroid"]])})
                 for a, b in lines]
        g.add_edges_from(lines)
    except AttributeError:
        pass # there are no lines
    return g


def reorganise_xarray_dimensions(data):
    """
    Reorganise Dataset or DataArray dimensions to be alphabetical *except*
    `timesteps`, which must always come last in any DataArray's dimensions
    """

    if not (isinstance(data, xr.Dataset) or isinstance(data, xr.DataArray)):
        raise TypeError('Must provide either xarray Dataset or DataArray to be reorganised')

    steps = [i for i in ['datesteps', 'timesteps'] if i in data.dims]

    if isinstance(data, xr.Dataset):
        new_dims = (
            sorted(list(set(data.dims.keys()) - set(steps)))
        ) + steps
    elif isinstance(data, xr.DataArray):
        new_dims = (
            sorted(list(set(data.dims) - set(steps)))
        ) + steps

    updated_data = data.transpose(*new_dims).reindex({k: data[k] for k in new_dims})

    return updated_data


def split_loc_techs(data_var, as_='DataArray'):
    """
    Get a DataArray with locations technologies, and possibly carriers
    split into separate coordinates.
    Parameters
    ----------
    data_var : xarray DataArray
        Variable from Calliope model_data, to split loc_techs dimension
    as_ : string
        'DataArray' to return xarray DataArray or 'Series' to return pandas
        Series with dimensions as a MultiIndex
    Returns
    -------
    updated_data_var : xarray DataArray of pandas Series
    """

    # Separately find the loc_techs(_carriers) dimension and all other dimensions
    loc_tech_dim = [i for i in data_var.dims if 'loc_tech' in i]
    if not loc_tech_dim:
        loc_tech_dim = [i for i in data_var.dims if 'loc_carrier' in i]

    if not loc_tech_dim:
        if as_ == 'Series':
            return data_var.to_series()
        elif as_ == 'DataArray':
            return data_var
        else:
            raise ValueError('`as_` must be `DataArray` or `Series`, '
                             'but `{}` given'.format(as_))

    loc_tech_dim = loc_tech_dim[0]
    # xr.Datarray -> pd.Series allows for string operations
    data_var_idx = data_var[loc_tech_dim].to_index()
    index_list = data_var_idx.str.split('::').tolist()

    # carrier_prod, carrier_con, and carrier_export will return an index_list
    # of size 3, all others will be an index list of size 2
    possible_names = ['loc', 'tech', 'carrier']
    names = [i + 's' for i in possible_names if i in loc_tech_dim]

    data_var_midx = pd.MultiIndex.from_tuples(index_list, names=names)

    # Replace the Datarray loc_tech_dim with this new MultiIndex
    updated_data_var = data_var.copy()
    updated_data_var.coords[loc_tech_dim] = data_var_midx
    updated_data_var = reorganise_xarray_dimensions(updated_data_var.unstack())

    if as_ == "Series":
        return updated_data_var.to_series()

    elif as_ == "DataArray":
        return updated_data_var

    else:
        raise ValueError('`as_` must be `DataArray` or `Series`, '
                         'but `{}` given'.format(as_))


if __name__ == "__main__":
    plot_map(
        path_to_shapes=snakemake.input.shapes,
        path_to_continental_shape=snakemake.input.continental_shape,
        path_to_national_shapes=snakemake.input.national_shapes,
        path_to_regional_shapes=snakemake.input.regional_shapes,
        path_to_continental_result=snakemake.input.continental_result,
        path_to_national_result=snakemake.input.national_result,
        path_to_regional_result=snakemake.input.regional_result,
        path_to_plot=snakemake.output[0],
        scaling_factor_cost=snakemake.params.scaling_factor_cost
    )
