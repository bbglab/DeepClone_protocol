
import pandas as pd


all_counts_df = pd.DataFrame()
# for file in os.listdir(f"{run}/processing_files/sortbamamfiltered/"):
run = "/data/bbg/nobackup2/prominent/duplex_seq_tests/error_rate/cord_blood/bbg/2026-03-17_tws_idt"

for sample in ["SC001_B1_1_H_2", "SC002_B1_1_H_2", "SC003_B1_1_H_2"]:
# SC001_B1_1_H_2,SC002_B1_1_H_2
    sample_counts = pd.read_csv(f"{run}/processing_files/sortbamamfiltered/{sample}.sorted.read_coords.grouped_cuts_freq.tsv.gz", sep="\t")
    sample_counts["sample"] = sample
    # sample_counts["expected_poisson_cuts"] = (sample_counts["relative_frequency"] * sample_counts["repeats"]).sum()
    all_counts_df = pd.concat([all_counts_df, sample_counts], ignore_index=True)

all_counts_df_combined = all_counts_df.groupby(by = ['CHROM', 'START_first', 'END_first', 'START_last', 'END_last'])["count"].sum().reset_index()
# all_counts_df_combined["count"]
all_counts_df_combined["count"].value_counts().to_frame(
    name = 'frequency').reset_index().rename(
        columns={"count": "repeats"}).to_csv(
            f"{run}/processing_files/sortbamamfiltered/IDT_combined.sorted.cuts_freq.tsv.gz",
            sep="\t", index=False, compression="gzip")

# .to_csv(f"{run}/processing_files/sortbamamfiltered/combined_cuts_freq.tsv", sep="\t", index=False)