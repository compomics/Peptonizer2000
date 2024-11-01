import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from .taxon_manager import TaxonManager


def plot_peptonizer_results(input_file: str, output_file: str, number_of_taxa: int = 25):
    """
    Read the results of a Peptonizer run from a CSV-file (denoted by the input_file argument) and write bar charts
    representing these results to a PNG-file.
    """
    assert input_file.lower().endswith(".csv"), "Input file should be a CSV."
    assert output_file.lower().endswith(".png"), "Output file should be a PNG."


    # read csv using pandas
    ids = pd.read_csv(input_file, names=["ID", "score", "type"])
    tax_ids = ids.loc[ids["type"] == "taxon"]
    tax_ids = tax_ids.dropna()

    tax_ids.loc[:, "score"] = pd.to_numeric(tax_ids["score"], downcast="float")
    tax_ids = tax_ids.sort_values("score")
    taxa_check = tax_ids.ID.tolist()

    # translate taxids to scientific names
    taxa_name_dict = TaxonManager.get_names_for_taxa([int(x) for x in tax_ids["ID"]])
    taxa_names = [taxa_name_dict[int(tax)] for tax in taxa_check]
    scores = tax_ids["score"]

    # make the barplot
    fig, ax = plt.subplots()
    fig.set_size_inches(30, 15)
    bars = ax.barh(
        range(len(taxa_names[-number_of_taxa:])),
        scores[-number_of_taxa:],
        color="royalblue",
    )

    ax.set_yticks(range(len(taxa_names[-number_of_taxa:])))
    ax.set_yticklabels(taxa_names[-number_of_taxa:], fontsize=25)
    plt.xlim((0, 1))
    plt.xlabel("Posterior probability", fontsize=35)
    ax.xaxis.set_ticks(np.arange(0, 1.2, 0.2))
    ax.xaxis.set_ticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    bar_color = bars[0].get_facecolor()

    for bar in bars:
        ax.text(
            bar.get_width() + 0.05,
            bar.get_y() + bar.get_height() / 5,
            round(bar.get_width(), 3),
            fontsize=35,
            horizontalalignment="center",
            color=bar_color,
            weight="bold",
        )

    fig.tight_layout()

    plt.savefig(output_file)
    plt.close()
