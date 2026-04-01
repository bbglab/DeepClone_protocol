
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


cord_blood_samples_list = [# "SC003_B1_1_H_1", "SC002_B1_1_H_1", "SC001_B1_1_H_1",
                           "SC003_B1_1_H_2", "SC002_B1_1_H_2", "SC001_B1_1_H_2"]

sample = sys.argv[1]
run = sys.argv[2]


data = pd.read_csv(f"{run}/processing_files/sortbamamfiltered/{sample}.sorted.read_coords.tsv.gz",
                   sep="\t",
                   header=None)

data.columns = ["READ_NAME", "CHROM", "START", "END"]

merged_cuts = (
    data
    .groupby(["READ_NAME", "CHROM"])
    .agg(
        START_first=("START", "min"),
        END_first=("END", "min"),
        START_last=("START", "max"),
        END_last=("END", "max"),
    )
    .reset_index()
)

# # %%
# first_in_pair = data.sort_values(by = ["READ_NAME", "CHROM", "START"]).drop_duplicates(subset=["READ_NAME"], keep='first')
# last_in_pair = data.sort_values(by = ["READ_NAME", "CHROM", "START"]).drop_duplicates(subset=["READ_NAME"], keep='last')

# # %%
# merged_cuts = first_in_pair.merge(last_in_pair, on=["READ_NAME", 'CHROM'], suffixes=("_first", "_last"))
# merged_cuts

identical_cuts_size = merged_cuts.groupby(["CHROM", "START_first", "END_first", "START_last", "END_last"]).size()
frequency_identical_cuts_size = identical_cuts_size.to_frame().reset_index().rename(columns={0: "count"})
frequency_identical_cuts_size.sort_values(by = 'count', ascending=False).to_csv(
    f"{run}/processing_files/sortbamamfiltered/{sample}.sorted.read_coords.grouped_cuts_freq.tsv.gz",
    sep="\t", index=False, compression="gzip")

frequency_identical_cuts_size["count"].value_counts().to_frame(
    name = 'frequency').reset_index().rename(
        columns={"count": "repeats"}).to_csv(
            f"{run}/processing_files/sortbamamfiltered/{sample}.sorted.read_coords.cuts_freq.tsv.gz",
            sep="\t", index=False, compression="gzip")