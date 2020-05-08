import numpy as np
import pandas as pd


def merge(path_to_mf_lf_1, path_to_mf_lf_2, path_to_mf_hf_1, path_to_mf_hf_2, path_to_sf,
          path_to_mf_lf, path_to_mf_hf, path_to_sf_copy):
    """Merge multi fidelity results to be able to assess differences between scenarios."""
    mf_lf = merge_one_case(pd.read_csv(path_to_mf_lf_1, index_col=0), pd.read_csv(path_to_mf_lf_2, index_col=0))
    mf_hf = merge_one_case(pd.read_csv(path_to_mf_hf_1, index_col=0), pd.read_csv(path_to_mf_hf_2, index_col=0))
    sf = pd.read_csv(path_to_sf, index_col=0) # only bypassed
    mf_lf.to_csv(path_to_mf_lf, index=True, header=True)
    mf_hf.to_csv(path_to_mf_hf, index=True, header=True)
    sf.to_csv(path_to_sf_copy, index=True, header=True)


def merge_one_case(mf1, mf2):
    for df in [mf1, mf2]:
        df["y-supply-gw"] = (
            df["y-wind-gw"]
            + df["y-pv-gw"]
            + df["y-hydro-gw"]
        )
        df["y-balancing-gw"] = (
            df["y-biofuel-gw"]
            + df["y-storage-gw"]
        )
    return pd.DataFrame(
        index=mf1.index,
        data={
            'd_r': mf1['d_r'],
            'c_util': mf1['c_util'],
            'c_roof': mf1['c_roof'],
            'c_wind': mf1['c_wind'],
            'c_offshore': mf1['c_offshore'],
            'c_sts_p': mf1['c_sts_p'],
            'c_sts_e': mf1['c_sts_e'],
            'c_lts_p': mf1['c_lts_p'],
            'c_lts_e': mf1['c_lts_e'],
            'c_ntc': mf1['c_ntc'],
            'c_bio': mf1['c_bio'],
            'ac_bio': mf1['ac_bio'],
            "y-large-scale-cost-eur": mf1["y-cost-eur"],
            "y-small-scale-cost-eur": mf2["y-cost-eur"],
            "y-cost-diff-eur": mf2["y-cost-eur"] - mf1["y-cost-eur"],
            "y-cost-diff-relative": relative_measure(mf1, mf2, "y-cost-eur"),
            "y-supply-diff-relative": relative_measure(mf1, mf2, "y-supply-gw"),
            "y-wind-diff-relative": relative_measure(mf1, mf2, "y-wind-gw"),
            "y-balancing-diff-relative": relative_measure(mf1, mf2, "y-balancing-gw"),
            "y-large-scale-pv-gw": mf1["y-pv-gw"],
            "y-small-scale-pv-gw": mf2["y-pv-gw"],
            "y-large-scale-wind-gw": mf1["y-wind-gw"],
            "y-small-scale-wind-gw": mf2["y-wind-gw"],
            "y-large-scale-hydro-gw": mf1["y-hydro-gw"],
            "y-small-scale-hydro-gw": mf2["y-hydro-gw"],
            "y-large-scale-biofuel-gw": mf1["y-biofuel-gw"],
            "y-small-scale-biofuel-gw": mf2["y-biofuel-gw"],
            "y-large-scale-storage-gw": mf1["y-storage-gw"],
            "y-small-scale-storage-gw": mf2["y-storage-gw"],
            "y-large-scale-storage-gwh": mf1["y-storage-gwh"],
            "y-small-scale-storage-gwh": mf2["y-storage-gwh"],
            "y-large-scale-transmission-gwkm": mf1["y-transmission-gwkm"],
            "y-small-scale-transmission-gwkm": mf2["y-transmission-gwkm"],
        }
    )


def relative_measure(mf1, mf2, measure):
    diff = (mf2[measure] - mf1[measure]) / mf1[measure]
    return diff.fillna(0).replace(np.inf, 0) # special case treatment not necessary for final runs


if __name__ == "__main__":
    merge(
        path_to_mf_lf_1=snakemake.input.mf_lf_1,
        path_to_mf_lf_2=snakemake.input.mf_lf_2,
        path_to_mf_hf_1=snakemake.input.mf_hf_1,
        path_to_mf_hf_2=snakemake.input.mf_hf_2,
        path_to_sf=snakemake.input.sf,
        path_to_mf_lf=snakemake.output.mf_lf,
        path_to_mf_hf=snakemake.output.mf_hf,
        path_to_sf_copy=snakemake.output.sf
    )
