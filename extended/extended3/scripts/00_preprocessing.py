#!/usr/bin/env python
# This script performs the preprocessing of the data for the extended analysis 2.
# It includes loading the data, filtering, and preparing it for the comparison of mutations between deepUMI and DupCaller.
# The input data includes the variant calls from both tools (via deepCSA in the deepUMIcaller case),
# and the output is a merged dataset that can be used for further analysis and plotting.
#
# Usage:
#   python 00_preprocessing.py \
#       --cohort bladder \
#       --deepcsa-path /path/to/deepcsa/mutations \
#       --dupcaller-path /path/to/all_samples_info.maf
#
# Arguments:
#   --cohort         Cohort name. Supported values: "bladder".
#   --deepcsa-path   Path to the deepCSA output directory (mutations folder).
#                    Expected to contain:
#                      germline_somatic/all_samples.filtered.tsv.gz
#                      germline_somatic/CallerDeepumicaller.filtered.tsv.gz
#                      clean_somatic/CallerDeepumicaller.somatic.mutations.tsv
#   --dupcaller-path Path to the DupCaller MAF file (all_samples_info.maf).
#
# Output:
#   ../data/{cohort}/comparison_deepUMI_dupCaller.mutations.tsv
#   A merged table of all variants (all mutation types) called by deepUMI and/or DupCaller,
#   with a "status" column indicating: "Common", "Only deepUMI", or "Only DupCaller".

# --- 1. Setup: imports -------------------------------------------
import click
import pandas as pd
import numpy as np
from utils import add_variant_id, build_merged_variants


@click.command()
@click.option("--cohort", required=True, type=click.Choice(["bladder"]), help="Cohort name.")
@click.option(
    "--deepcsa-path",
    required=True,
    type=click.Path(exists=True),
    help="Path to the deepCSA output directory (mutations folder).",
)
@click.option(
    "--dupcaller-path",
    required=True,
    type=click.Path(exists=True),
    help="Path to the DupCaller MAF file.",
)
def main(cohort: str, deepcsa_path: str, dupcaller_path: str) -> None:
    """Merge deepUMI and DupCaller variant calls and save the result as a TSV."""

    # --- Cohort-specific sample lists -------------------------------------------
    if cohort == "bladder":
        samples = [
            "P19_0050_BDO_01", "P19_0050_BTR_01",
            "P19_0051_BDO_01", "P19_0051_BTR_01",
            "P19_0052_BDO_01", "P19_0052_BTR_01",
            "P19_0053_BDO_01", "P19_0053_BTR_01",
        ]

    dupcaller_samples = [f"{sample}_dupcaller" for sample in samples]

    # --- Derived paths (shared structure for both cohorts) ----------------------
    all_both = f"{deepcsa_path}/germline_somatic/all_samples.filtered.tsv.gz"
    all_deepumi = f"{deepcsa_path}/germline_somatic/CallerDeepumicaller.filtered.tsv.gz"
    deepumi_clean = f"{deepcsa_path}/clean_somatic/CallerDeepumicaller.somatic.mutations.tsv"

    print(f"Cohort: {cohort}  |  Samples: {len(samples)}")

    # --- Read raw tables --------------------------------------------------------
    dupcaller_df = pd.read_csv(dupcaller_path, sep="\t", low_memory=False)
    deepumi_df = pd.read_csv(all_deepumi, sep="\t", low_memory=False)
    deepumi_clean_df = pd.read_csv(deepumi_clean, sep="\t", low_memory=False)
    mutations_table = pd.read_csv(all_both, sep="\t", low_memory=False, compression="gzip")

    # --- Restrict to specified samples ------------------------------------------
    dupcaller_df = dupcaller_df[dupcaller_df["SAMPLE_ID"].isin(dupcaller_samples)]
    deepumi_df = deepumi_df[deepumi_df["SAMPLE_ID"].isin(samples)]
    deepumi_clean_df = deepumi_clean_df[deepumi_clean_df["SAMPLE_ID"].isin(samples)]

    # --- Add VARIANT_ID to both caller tables -----------------------------------
    dupcaller_df = add_variant_id(dupcaller_df)
    deepumi_df = add_variant_id(deepumi_df)

    # --- Mark PASS status for deepUMI using the clean somatic set ---------------
    deepumi_df["is_pass"] = deepumi_df["VARIANT_ID"].isin(deepumi_clean_df["MUT_ID"])

    # --- Harmonise mutation types (deepUMI reports INSERTION/DELETION separately)
    deepumi_df["TYPE"] = deepumi_df["TYPE"].replace({"INSERTION": "INDEL", "DELETION": "INDEL"})

    print(f"deepUMI variants: {len(deepumi_df):,}  |  DupCaller variants: {len(dupcaller_df):,}")
    print(f"Mutations table rows: {len(mutations_table):,}")

    # --- Build merged table of variants passing both callers --------------------
    merged_pass_both_df = build_merged_variants(
        deepumi_df, dupcaller_df,
        deepumi_pass_only=True, dupcaller_pass_only=True,
    )

    # --- Add TYPE column: prefer deepUMI type for Common/Only deepUMI rows ------
    merged_pass_both_df["TYPE"] = np.where(
        (merged_pass_both_df["status"] == "Only deepUMI") | (merged_pass_both_df["status"] == "Common"),
        merged_pass_both_df["TYPE_deepumi"],
        np.where(
            merged_pass_both_df["status"] == "Only DupCaller",
            merged_pass_both_df["TYPE_dupcaller"],
            "Unknown",
        ),
    )
    
    # --- Generate TSV output for plotting ---------------------------------------
    output_path = f"../data/{cohort}/comparison_deepUMI_dupCaller.mutations.tsv"
    merged_pass_both_df.to_csv(output_path, sep="\t", index=False)
    print(f"Output saved to: {output_path}")


if __name__ == "__main__":
    main()