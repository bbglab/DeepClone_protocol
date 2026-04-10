# --- Utils to generate plots for extended figure 2 ---

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def add_variant_id(df: pd.DataFrame) -> pd.DataFrame:
    """Create a unique variant identifier as CHROM:POS_REF>ALT for merging across callers."""
    df["VARIANT_ID"] = (
        df["CHROM"].astype(str) + ":"
        + df["POS"].astype(str) + "_"
        + df["REF"].astype(str) + ">"
        + df["ALT"].astype(str)
    )
    return df

def build_merged_variants(
    sample_deepumi_df: pd.DataFrame,
    sample_dupcaller_df: pd.DataFrame,
    deepumi_pass_only: bool = False,
    dupcaller_pass_only: bool = False,
) -> pd.DataFrame:
    """Outer-merge deepUMI and DupCaller variant tables and label each row.

    Parameters
    ----------
    sample_deepumi_df : DataFrame
        deepUMIcaller variants (must contain VARIANT_ID, SAMPLE_ID, FILTER, VAF,
        TYPE, DEPTH, ALT_DEPTH, and the boolean ``is_pass`` column).
    sample_dupcaller_df : DataFrame
        DupCaller variants (must contain VARIANT_ID, SAMPLE_ID, FILTER, LR,
        TYPE, DEPTH, ALT_DEPTH, INFO).
    deepumi_pass_only : bool
        If True, keep only deepUMI variants where ``is_pass == True``.
    dupcaller_pass_only : bool
        If True, keep only DupCaller variants where ``FILTER == "PASS"``.

    Returns
    -------
    DataFrame
        Merged table with a ``status`` column: "Common", "Only deepUMI", or
        "Only DupCaller".
    """
    # Apply PASS filters when requested
    if deepumi_pass_only:
        deepumi_subset = sample_deepumi_df[sample_deepumi_df["is_pass"] == True].copy()
    else:
        deepumi_subset = sample_deepumi_df.copy()

    if dupcaller_pass_only:
        dupcaller_subset = sample_dupcaller_df[sample_dupcaller_df["FILTER"] == "PASS"].copy()
    else:
        dupcaller_subset = sample_dupcaller_df.copy()

    # Select and rename columns so they don't collide after the merge
    deepumi_subset = deepumi_subset[
        ["VARIANT_ID", "SAMPLE_ID", "FILTER", "VAF", "TYPE", "DEPTH", "ALT_DEPTH"]
    ].rename(columns={
        "FILTER": "FILTER_deepumi",
        "VAF": "VAF_deepumi",
        "TYPE": "TYPE_deepumi",
        "DEPTH": "DEPTH_deepumi",
        "ALT_DEPTH": "ALT_DEPTH_deepumi",
    })

    dupcaller_subset = dupcaller_subset[
        ["SAMPLE_ID", "VARIANT_ID", "FILTER", "LR", "TYPE", "DEPTH", "ALT_DEPTH", "INFO"]
    ].rename(columns={
        "FILTER": "FILTER_dupcaller",
        "LR": "LR_dupcaller",
        "TYPE": "TYPE_dupcaller",
        "DEPTH": "DEPTH_dupcaller",
        "ALT_DEPTH": "ALT_DEPTH_dupcaller",
        "INFO": "INFO_dupcaller",
    })

    # Remove the "_dupcaller" suffix so SAMPLE_IDs match
    dupcaller_subset["SAMPLE_ID"] = dupcaller_subset["SAMPLE_ID"].str.replace(
        "_dupcaller", "", regex=False
    )

    # Outer merge with indicator
    merged_df = pd.merge(
        deepumi_subset, dupcaller_subset,
        on=["VARIANT_ID", "SAMPLE_ID"],
        how="outer",
        indicator=True,
    )

    status_map = {
        "both": "Common",
        "left_only": "Only deepUMI",
        "right_only": "Only DupCaller",
    }
    merged_df["status"] = merged_df["_merge"].map(status_map)
    return merged_df.drop(columns=["_merge"])


def status_counts(merged_df: pd.DataFrame) -> tuple[int, int, int]:
    """Return (only_deepumi, only_dupcaller, common) counts from a merged table."""
    only_deepumi = int((merged_df["status"] == "Only deepUMI").sum())
    only_dupcaller = int((merged_df["status"] == "Only DupCaller").sum())
    common = int((merged_df["status"] == "Common").sum())
    return only_deepumi, only_dupcaller, common

def save_figure(fig: plt.Figure, *path_parts: str, output_base: Path) -> None:
    """Save a figure under OUTPUT_BASE, creating parent directories as needed."""
    save_path = output_base / Path(*path_parts)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path)
    plt.close(fig)