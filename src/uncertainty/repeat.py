import pandas as pd
import calendar


def repeat_timeseries(path_to_timeseries, start_year, path_to_output):
    """Take a timeseries of one year and repeat it over more years."""
    ts = pd.read_csv(path_to_timeseries, index_col=0, parse_dates=True)
    assert len(ts.index) in [8760, 8784] # assume one year in hours
    final_year = ts.index[0].year
    if len(ts.index) == 8760:
        regular_year = ts
        leap_year = leap_from_regular(ts)
    else:
        leap_year = ts
        regular_year = regular_from_leap(ts)
    frankenstein = pd.concat(
        [leap_year if calendar.isleap(year) else regular_year
         for year in range(start_year, final_year + 1)]
    )
    frankenstein.index = pd.date_range(
        start=f"{start_year}-01-01 00:00",
        end=f"{final_year + 1}-01-01 00:00",
        closed="left",
        freq="H"
    )
    frankenstein.to_csv(path_to_output, header=True, index=True)


def leap_from_regular(ts):
    year = ts.index[0].year
    february = f"{year}-02"
    february28 = f"{year}-02-28"
    march = f"{year}-03"
    return pd.concat([ts[:february], ts[february28], ts[march:]])


def regular_from_leap(ts):
    year = ts.index[0].year
    february28 = f"{year}-02-28"
    march = f"{year}-03"
    return pd.concat([ts[:february28], ts[march:]])


if __name__ == "__main__":
    repeat_timeseries(
        path_to_timeseries=snakemake.input.timeseries,
        start_year=snakemake.params.start,
        path_to_output=snakemake.output[0]
    )
