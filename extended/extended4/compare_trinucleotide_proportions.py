#!/usr/bin/env python


import sys
import click
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sys.path.append('../../') 
from consensus_variables import *


def plot_trinucleotide_proportions(wgs_counts_file):
    """  
    This function generates two scatter plots showing how trinucleotide contexts are distributed across  
    WGS versus four consensus panels: all, non_protein_affecting, introns_intergenic, and exons_splice_sites.  

    Args:  
        wgs_counts_file (str): Path to tab-separated file containing WGS trinucleotide context counts.
    """  
    wgs_counts = pd.read_table(wgs_counts_file)
    wgs_counts.columns = ['CONTEXT', 'COUNT_WGS']


    counts_all = wgs_counts.copy()
    for cnsq in ["all", "non_protein_affecting", "introns_intergenic", "exons_splice_sites"]:
        try:
            wgs_counts_cnsq = pd.read_table(f"data/consensus.{cnsq}.tsv",
                                            header = 0,
                                            sep = '\t',
                                            usecols = ["CHROM", "POS", "CONTEXT_MUT", "CONTEXT"]
                                        )
            wgs_counts_cnsq = wgs_counts_cnsq.drop_duplicates()
            counts_panel_cnsq = wgs_counts_cnsq["CONTEXT"].value_counts().to_frame(name = f'COUNT_{cnsq}').reset_index()
            counts_all = counts_all.merge(counts_panel_cnsq, on = 'CONTEXT')

        except FileNotFoundError:
            print(f"File not found for consensus panel: {cnsq}")

    counts_all = counts_all.set_index("CONTEXT")
    proportions_all = counts_all / counts_all.sum()



    proportions_all_plot = proportions_all.copy()

    fig, axs = plt.subplots(1, 2, figsize=(4.2, 2))
    axs = axs.flatten()  # Flatten for easy iteration

    for i, cnsq in enumerate(["all", "non_protein_affecting"]):
        ax = axs[i]
        try :
            rmse = np.sqrt(((proportions_all_plot["COUNT_WGS"] - proportions_all_plot[f"COUNT_{cnsq}"])**2).mean())

        except KeyError:
            print(f"KeyError: 'COUNT_{cnsq}' not found in proportions_all_plot. Skipping RMSE calculation for this panel.")
            continue
        
        # Scatter plot
        sns.scatterplot(data=proportions_all_plot,
                        x="COUNT_WGS",
                        y=f"COUNT_{cnsq}",
                        hue="CONTEXT",
                        legend=False,
                        s = 10,
                        ax=ax)
        
        # Identity line (x = y)
        lims = [
            min(proportions_all_plot["COUNT_WGS"].min(), proportions_all_plot[f"COUNT_{cnsq}"].min()),
            max(proportions_all_plot["COUNT_WGS"].max(), proportions_all_plot[f"COUNT_{cnsq}"].max()),
        ]
        ax.plot(lims, lims, '--', color='gray')
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        
        # Annotate points with CONTEXT
        for namee, row in proportions_all_plot.iterrows():
            ax.text(row["COUNT_WGS"], row[f"COUNT_{cnsq}"], namee,
                    fontsize=3.5, alpha=0.7)
        
        ax.set_xlabel("Proportion WGS")
        ax.set_ylabel(f"Proportion {cnsq}")
        ax.set_title(f"{cnsq} (RMSE: {rmse:.4f})")

    fig.suptitle("Proportions of trinucleotide contexts in WGS vs. consensus panels", fontsize=7)
    plt.tight_layout()
    plt.savefig("Trinucleotide_content_comparisons.proportions.pdf", bbox_inches='tight', dpi=300)
    plt.show()



    # ## Count based
    # counts_all_plot = counts_all.copy()
    # fig, axs = plt.subplots(2, 2, figsize=(10, 9))
    # axs = axs.flatten()  # Flatten to make iteration easier

    # for i, cnsq in enumerate(["all", "non_protein_affecting", "introns_intergenic", "exons_splice_sites"]):
    #     ax = axs[i]

    #     # rmse = np.sqrt(((counts_all_plot["COUNT_WGS"] - counts_all_plot[f"COUNT_{cnsq}"])**2).mean())
    #     try:
    #         sns.scatterplot(data=counts_all_plot,
    #                         x="COUNT_WGS",
    #                         y=f"COUNT_{cnsq}",
    #                         hue="CONTEXT",
    #                         legend=False,
    #                         ax=ax)
    #     except KeyError:
    #         print(f"KeyError: 'COUNT_{cnsq}' not found in counts_all_plot. Skipping scatter plot for this panel.")
    #         continue
    
    #     # Annotate points with CONTEXT
    #     for namee, row in counts_all_plot.iterrows():
    #         ax.text(row["COUNT_WGS"], row[f"COUNT_{cnsq}"], namee,
    #                 fontsize=6, alpha=0.7)

    #     ax.set_xlabel("Counts WGS")
    #     ax.set_ylabel(f"Counts {cnsq}")
    #     ax.set_title(f"{cnsq}")
    # fig.suptitle("Counts of trinucleotide contexts in WGS vs. consensus panels", fontsize=16)
    # plt.tight_layout()
    # plt.savefig("Trinucleotide_content_comparisons.counts.pdf", dpi=300)
    # plt.show()


def main():
    click.echo("Comparing the trinucleotide proportions...")
    plot_trinucleotide_proportions("data/trinuc_counts.homo_sapiens.tsv")

if __name__ == '__main__':
    main()

