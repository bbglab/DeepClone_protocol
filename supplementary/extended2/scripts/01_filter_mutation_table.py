#!/usr/bin/env python
# This script performs the filtering of the germline/somatic mutation table for the extended analysis 2.
# The table is obtained from the deepCSA output runned using deepUMIcaller and DupCaller variant calls for 8 bladder samples.
# The path inside deepCSA output is: mutations/germline_somatic/all_samples.filtered.tsv.gz
# The filtering is done to keep only somatic variants taking into account the definition of each caller. For that,
# this script will use the comparison_deepUMI_dupCaller.mutations.tsv obtained in the preprocessing step (00_preprocessing.py).
# We do this filtering to avoid publishing germline variants.
# Usage:
#   python 01_filter_mutation_table.py \
#       --cohort bladder \
#       --comparison-file ../data/bladder/comparison_deepUMI_dupCaller.mutations.tsv \
#       --mutations-table /path/to/deepcsa/mutations/germline_somatic/all_samples.filtered.tsv.gz
#
# Arguments:
#   --cohort:          Cohort name (used to name the output file).
#   --comparison-file: Path to the comparison file between deepUMI and DupCaller variant calls,
#                      obtained in the preprocessing step (00_preprocessing.py).
#   --mutations-table: Path to the deepCSA all_samples.filtered.tsv.gz file
#                      (mutations/germline_somatic/all_samples.filtered.tsv.gz).
#
# Output:
#   ../data/{cohort}/all_samples.filtered.somatic.tsv.gz:
#   A filtered table of somatic variants called by deepUMI and/or DupCaller,
#   with the three nucleotide context to plot the mutational profile.

# --- 1. Setup: imports -------------------------------------------
import click
import pandas as pd


@click.command()
@click.option("--cohort", required=True, type=click.Choice(["bladder"]), help="Cohort name.")
@click.option(
    "--comparison-file",
    required=True,
    type=click.Path(exists=True),
    help="Path to the comparison TSV file (comparison_deepUMI_dupCaller.mutations.tsv).",
)
@click.option(
    "--mutations-table",
    required=True,
    type=click.Path(exists=True),
    help="Path to the deepCSA all_samples.filtered.tsv.gz file.",
)
def main(cohort: str, comparison_file: str, mutations_table: str) -> None:
    """Filter the germline/somatic mutation table to keep only somatic variants."""

    comparison_df = pd.read_csv(comparison_file, sep="\t", low_memory=False)
    mutations_df = pd.read_csv(mutations_table, sep="\t", low_memory=False)

    mutations_table_snvs = mutations_df[mutations_df["TYPE"] == "SNV"].copy()
    mutations_table_snvs["SAMPLE_BASE"] = mutations_table_snvs["SAMPLE_ID"].str.replace(
        "_dupcaller", "", regex=False
    )

    # Filter the mutations table to keep only somatic variants according to the comparison file
    somatic_variants = mutations_table_snvs.merge(
        comparison_df[["VARIANT_ID", "SAMPLE_ID", "status"]],
        left_on=["MUT_ID", "SAMPLE_BASE"],
        right_on=["VARIANT_ID", "SAMPLE_ID"],
        how="inner",
    )
    # Deduplicate Common variants (appear once per caller in the mutations table)
    somatic_variants = somatic_variants.drop_duplicates(subset=["MUT_ID", "SAMPLE_BASE"])
    # Remove duplicated SAMPLE_ID column from the merge
    somatic_variants = somatic_variants.drop(columns=["SAMPLE_ID_y"]).rename(columns={"SAMPLE_ID_x": "SAMPLE_ID"})

    # Save the filtered table of somatic variants
    output_path = f"../data/{cohort}/all_samples.filtered.somatic.tsv.gz"
    somatic_variants.to_csv(output_path, sep="\t", index=False, compression="gzip")
    print(f"Output saved to: {output_path}")


if __name__ == "__main__":
    main()