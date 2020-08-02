import numpy as np
import shapely
import geopandas as gpd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns


EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "
GREY = "#C0C0C0"
BLUE = "#4F6DB8"
YELLOW = "#FABC3C"
SUPPLY_TECHS = [
    "hydro_reservoir", "hydro_run_of_river", "open_field_pv",
    "roof_mounted_pv", "wind_offshore", "wind_onshore_competing",
    "wind_onshore_monopoly"
]
DEMAND_TECH = "demand_elec"

MAP_MIN_X = 2200000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6300000
MAP_MAX_Y = 5500000


def bubble_map(path_to_shapes, path_to_continent_shape, scenario, resolution_km, colour, markersize,
               path_to_results, path_to_output):
    colour = {"yellow": YELLOW, "blue": BLUE}[colour]
    continent = (
        gpd
        .read_file(path_to_continent_shape)
        .to_crs(EPSG_3035_PROJ4)
        .rename(columns={"id": "locs"})
        .set_index("locs")
        .rename(index=lambda idx: idx.replace(".", "-"))
    )
    shapes = read_shapes(path_to_shapes, path_to_results, scenario)
    points = points_on_shape(continent.geometry.iloc[0], resolution_km2=resolution_km)
    points = generation_per_point(points, shapes)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.subplots(1, 1)
    continent.plot(ax=ax, color=GREY, alpha=0.2)
    points.plot(ax=ax, color=colour, markersize=points["generation"] if markersize == "gen" else int(markersize))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)
    sns.despine(fig=fig, top=True, bottom=True, left=True, right=True)

    fig.savefig(path_to_output)


def read_shapes(path_to_shapes, path_to_results, scenario):
    shapes = (
        gpd
        .read_file(path_to_shapes)
        .to_crs(EPSG_3035_PROJ4)
        .rename(columns={"id": "locs"})
        .set_index("locs")
        .rename(index=lambda idx: idx.replace(".", "-"))
    )
    ds = xr.open_dataset(path_to_results)
    demand_twh = (
        ds
        .carrier_con
        .sel(techs=DEMAND_TECH, scenario=scenario)
        .to_series()
        .reindex(shapes.index)
        .div(1e6)
        .mul(-1)
    )
    generation_twh = (
        ds
        .carrier_prod
        .sel(techs=SUPPLY_TECHS, scenario=scenario)
        .sum("techs")
        .to_series()
        .reindex(shapes.index)
        .div(1e6)
    )
    shapes["generation"] = generation_twh / demand_twh
    return shapes


def generation_per_point(points, shapes):
    points = gpd.sjoin(
        gpd.GeoDataFrame(geometry=points),
        shapes,
        how="left",
        op="within"
    )
    points.generation.fillna(value=0, inplace=True)
    points.index_right.fillna(value=0, inplace=True)
    points["generation"] = points.groupby("index_right").generation.transform(lambda x: x / x.count())

    max_value = 100
    points["generation"] = points["generation"] * 10
    points["generation"].where(points["generation"] < max_value, max_value, inplace=True)
    return points


def points_on_shape(shape_3035, resolution_km2):
    x_min, y_min, x_max, y_max = shape_3035.bounds
    all_points = [
        shapely.geometry.Point(x, y)
        for x in np.arange(start=x_min, stop=x_max, step=resolution_km2 * 1000)
        for y in np.arange(start=y_min, stop=y_max, step=resolution_km2 * 1000)
    ]
    simplification_strength = resolution_km2 * 1000 / 20
    surface_area = (
        shape_3035
        .simplify(simplification_strength)
    )
    prepared_shape = shapely.prepared.prep(surface_area)
    return gpd.GeoSeries(
        list(filter(
            lambda point: prepared_shape.intersects(point),
            all_points
        )),
        crs=EPSG_3035_PROJ4
    )


if __name__ == "__main__":
    bubble_map(
        path_to_shapes=snakemake.input.shapes,
        path_to_continent_shape=snakemake.input.continent_shape,
        scenario=snakemake.wildcards.scenario,
        colour=snakemake.wildcards.colour,
        markersize=snakemake.wildcards.markersize,
        resolution_km=snakemake.params.resolution_km,
        path_to_results=snakemake.input.results,
        path_to_output=snakemake.output[0]
    )
