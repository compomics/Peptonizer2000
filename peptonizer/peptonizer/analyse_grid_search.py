import pandas as pd
import rbo

from scipy.stats import entropy
from typing import List, NamedTuple, Tuple


class ParameterSet(NamedTuple):
    """
    Represents a set of parameters that have been used for a Peptonizer analysis to tweak the behaviour of the belief
    propagation algorithm.
    """
    alpha: float
    beta: float
    prior: float


def compute_goodness(df: pd.DataFrame, taxid_weights: pd.DataFrame):
    tax_ids = df.loc[df['type'] == 'taxon'].copy()
    tax_ids.loc[:, 'score'] = pd.to_numeric(tax_ids['score'], downcast = 'float')
    tax_ids.sort_values('score', ascending = False, inplace = True)
    tax_ids.dropna(inplace = True)

    # Compute entropy of the posterior probability distribution.
    computed_entropy = entropy(tax_ids['score'].tolist())

    # Compute the rank based similarity between the weight sorted taxa and the score sorted ID results.

    return rbo.RankingSimilarity(
        taxid_weights['HigherTaxa'].values,
        [int(i) for i in tax_ids['ID'].values]
    ).rbo() * (1 / computed_entropy ** 2)


def find_best_parameters(results: List[Tuple[pd.DataFrame, ParameterSet]], taxid_weights: pd.DataFrame):
    """
    Given the dataframes that have been run through the Belief Propagation Algorithm before and the matching parameter
    sets, compute a goodness metric for each of these dataframes and returns the ParameterSet that resulted in the
    highest goodness value.

    :param results: A list of tuples each holding two things:
        1. A dataframe containing taxa and their associated scores after running the belief propagation algorithm
        2. The parameter values that where used during the belief propagation algorithm for this set of taxa
    :param taxid_weights: A dataframe containing taxa and their corresponding 'scaled weights', as computed by the]
    weight_taxa step.
    """
    params = []
    goodness_list = []

    for df, param_set in results:
        goodness_list.append(compute_goodness(df, taxid_weights))
        params.append(param_set)

    metrics_params = zip(goodness_list, params)
    sorted_metric_param_pairs = sorted(metrics_params, reverse = True)

    # Return the ParameterSet that's associated with the highest computed goodness metric.
    return sorted_metric_param_pairs[0][1]
