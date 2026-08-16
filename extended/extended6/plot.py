#!/usr/bin/env python


import sys
import click
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('../../') 
from consensus_variables import *

def plot_similarity_heatmaps(cosine_similarity_df, mode, selected_cols, label):

    # Plot clustermap
    clustergrid = sns.clustermap(cosine_similarity_df,
                    annot=False,
                    # fmt=".2f",
                    cmap='coolwarm',
                    cbar_kws={'label': 'Cosine Similarity'},
                    figsize=(max(4, 0.11*len(selected_cols)), max(3, 0.08*len(selected_cols))),
                    )
    plt.suptitle(f'Cosine Similarity Clustermap ({label})')
    clustergrid.savefig(f"{mode}.cosine_similarity_clustermap_{label}.pdf", dpi=300, bbox_inches='tight')
    plt.close()




def main():
    click.echo("Combining signature probabilities of all samples and groups...")
    cosine_similarity_df = pd.read_table("all.cosine_similarity_all.tsv", header=0, index_col=0, sep="\t")
    all_groups = cosine_similarity_df.index.tolist()

    samples = [g for g in all_groups if g.startswith("P19")]
    groups = [g for g in all_groups if not g.startswith("P19")]

    plot_similarity_heatmaps(cosine_similarity_df.loc[samples, samples], "all", samples, "samples")
    plot_similarity_heatmaps(cosine_similarity_df.loc[groups, groups], "all", groups, "groups")

if __name__ == '__main__':
    main()
