import os
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, gaussian_kde, kendalltau
import csv
from collections import defaultdict, Counter
import seaborn as sns
from itertools import combinations
import sys
from sklearn.metrics import normalized_mutual_info_score
import math, statsmodels
from scipy.spatial.distance import pdist, squareform
import itertools
from statsmodels.stats.multitest import multipletests
import argparse

np.seterr(all='ignore')

def parse_args():
    parser = argparse.ArgumentParser(
        description="GRN analysis with TF-DNA, PPI priors, and RNA-seq data"
    )

    parser.add_argument(
        "-p", "--p_value_cutoff",
        type=float,
        default=0.01,
        help="P-value cutoff (e.g., 0.01)"
    )

    parser.add_argument(
        "-g", "--gene_of_interest",
        type=str,
        required=True,
        help="Path to gene ID file (e.g., ./data/All_gene_ID.txt)"
    )

    parser.add_argument(
        "-r", "--RNAseq_data",
        type=str,
        required=True,
        help="Path to RNA-seq TMM-normalized CPM file"
    )

    parser.add_argument(
        "-t", "--TF_DNA_prior",
        type=str,
        required=True,
        help="Path to TF–DNA prior JSON file"
    )

    parser.add_argument(
        "-c", "--protein_protein_colocalization_prior",
        type=str,
        required=True,
        help="Path to protein–protein colocalization prior JSON file"
    )

    parser.add_argument(
        "-a", "--annotation",
        type=str,
        required=True,
        help="Path to genome annotation file"
    )

    parser.add_argument(
        "-b", "--YEP_replicate_ID",
        type=str,
        required=True,
        help="Path to YEP replicate ID file"
    )

    return parser.parse_args()

args = parse_args()

def kde_likelihood_empirical_p(A, B, bw_method=None, n_permutations=2000,
                               alternative='greater', eps=1e-300,
                               random_seed=None, return_details=False):
    """
    Compute KDE-based geometric-mean score for B under KDE(A) and convert it to an
    empirical p-value via permutation.

    Parameters
    ----------
    A, B : array-like
        Samples. A is used to build the KDE, B are the test points.
    bw_method : str or scalar or callable, optional
        Passed to scipy.stats.gaussian_kde.
    n_permutations : int
        Number of permutations for the empirical null (must be >= 1).
    alternative : {'greater','less','two-sided'}
        How to treat extremeness:
          - 'greater'  : larger log-gm (i.e. larger gm) is considered more extreme
          - 'less'     : smaller log-gm is considered more extreme
          - 'two-sided': two-sided test (double smaller tail)
    eps : float
        Small value added inside log to avoid log(0).
    random_seed : int or None
        Seed for reproducibility.
    return_details : bool
        If True, return a dict with extra information:
          {'p_emp', 'obs_loggm', 'obs_gm', 'perm_loggms', 'perm_gms'}.

    Returns
    -------
    p_emp : float
        Empirical p-value (in (0,1]).
    or
    details : dict
        If return_details True, see keys above.
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1")

    nA = len(A)
    nB = len(B)
    if nA < 2:
        raise ValueError("A must contain at least 2 points for KDE")

    # helper: compute log geometric mean of KDE(Atrain) evaluated at Btest
    def log_geo_mean_for_split(Atrain, Btest):
        kde = gaussian_kde(Atrain, bw_method=bw_method)
        dens = kde(Btest)
        logdens = np.log(dens + eps)
        return float(np.mean(logdens)), np.exp(np.mean(logdens))  # (log-gm, gm)

    # observed score (A as training, B as test)
    obs_loggm, obs_gm = log_geo_mean_for_split(A, B)

    # permutation test: shuffle combined and recompute log-gm for each split
    combined = np.concatenate([A, B])
    perm_loggms = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        perm = np.random.permutation(combined)
        Aperm = perm[:nA]
        Bperm = perm[nA:]
        lgm, _ = log_geo_mean_for_split(Aperm, Bperm)
        perm_loggms[i] = lgm

    # compute empirical p-value according to alternative
    if alternative == 'greater':
        # larger log-gm (less negative) than null is evidence
        p_emp = (np.sum(perm_loggms >= obs_loggm) + 1) / (n_permutations + 1)
    elif alternative == 'less':
        p_emp = (np.sum(perm_loggms <= obs_loggm) + 1) / (n_permutations + 1)
    elif alternative == 'two-sided':
        greater = (np.sum(perm_loggms >= obs_loggm) + 1) / (n_permutations + 1)
        less    = (np.sum(perm_loggms <= obs_loggm) + 1) / (n_permutations + 1)
        p_emp = 2.0 * min(greater, less)
        p_emp = min(p_emp, 1.0)
    else:
        raise ValueError("alternative must be 'greater', 'less', or 'two-sided'")

    if return_details:
        return {
            'p_emp': float(p_emp),
            'obs_loggm': float(obs_loggm),
            'obs_gm': float(obs_gm),
            'perm_loggms': perm_loggms,
            'perm_gms': np.exp(perm_loggms)
        }
    return float(1-p_emp)

def teset_and_merge_welch(optimal_clusters, optimal_mapping, replicate_variability=None, alpha=args.p_value_cutoff):
    """
    Merge clusters greedily using Welch's t-test:
      - compute pairwise Welch t-test p-values,
      - find the pair with the largest p-value,
      - if largest p-value >= alpha, merge that pair (j into i) and repeat.

    Args:
        optimal_clusters: list of lists (clusters), each cluster contains numeric values
        optimal_mapping : dict-like mapping old cluster idx -> list of original indices
                          (it is only used to seed the mapping; keys should correspond to
                           positions in optimal_clusters). Fallback behavior preserved.
        replicate_variability: kept in signature for compatibility but not used in Welch merging.
        alpha: significance threshold for merging (default 0.01). If max pairwise p >= alpha -> merge.

    Returns:
        merged_clusters : list of lists (clusters after merging)
        merged_mapping  : dict with consecutive integer keys 0..(m-1) mapping to lists
                          of original indices that went into each merged cluster
    """
    # helper: average pairwise Euclidean distance inside a sublist (kept for compatibility)
    def avg_pairwise_dist(sub):
        arr = np.asarray(sub, dtype=float)
        arr = arr[~np.isnan(arr)]
        if arr.size < 2:
            return 0.0
        diffs = np.abs(arr[:, None] - arr)
        triu = diffs[np.triu_indices(arr.size, k=1)]
        return float(np.mean(triu))

    # Step 0: defensive conversions
    merged = [list(sub) for sub in optimal_clusters]

    # Build an aligned list of mapping entries. Prefer sequential indices if present.
    try:
        merged_map = [list(optimal_mapping[i]) for i in range(len(optimal_clusters))]
    except Exception:
        keys_sorted = sorted(optimal_mapping.keys())
        merged_map = [list(optimal_mapping[k]) for k in keys_sorted]
        if len(merged_map) < len(merged):
            start = max(keys_sorted) + 1 if keys_sorted else 0
            for idx in range(len(merged_map), len(merged)):
                merged_map.append([start + (idx - len(merged_map))])

    # Step 1: compute mean of within-cluster average pairwise distances (kept for completeness)
    within_dists = [avg_pairwise_dist(sub) for sub in merged]
    mean_within = float(np.mean(within_dists)) if within_dists else 0.0

    # Step 2: iterative greedy merging using Welch's t-test
    while True:
        n = len(merged)
        if n <= 1:
            break

        # build pairwise p-value matrix (initialize to -inf so invalid pairs are ignored)
        p_mat = np.full((n, n), -np.inf, dtype=float)

        # compute p-values for i < j
        for i in range(n):
            for j in range(i + 1, n):
                a = np.asarray(merged[i], dtype=float)
                b = np.asarray(merged[j], dtype=float)
                # drop NaNs
                a = a[~np.isnan(a)]
                b = b[~np.isnan(b)]

                # require at least 2 observations in each group for Welch's t-test
                if a.size < 2 or b.size < 2:
                    # cannot compute a reliable Welch t-test; leave as -inf (not mergeable)
                    continue

                # Welch's t-test
                try:
                    t_stat, p_val = ttest_ind(a, b, equal_var=False)
                    if len(a) < len(b):
                        gm = kde_likelihood_empirical_p(b, a, bw_method='scott', n_permutations=2000, alternative='greater', random_seed=0)
                    else:
                        gm = kde_likelihood_empirical_p(a, b, bw_method='scott', n_permutations=2000, alternative='greater', random_seed=0)
                    p_val = max(p_val, gm)
                except Exception:
                    # unexpected numeric error: mark as not mergeable
                    continue

                # treat NaN p-values as not mergeable
                if np.isnan(p_val):
                    continue

                p_mat[i, j] = p_val
                p_mat[j, i] = p_val

        # find the pair with the largest p-value
        max_p = np.max(p_mat)
        if not np.isfinite(max_p) or max_p < alpha:
            # no pair non-significant at alpha -> stop
            break

        # get indices of the maximum p-value (first occurrence)
        flat_idx = np.argmax(p_mat)
        i, j = divmod(flat_idx, n)
        # ensure i < j (we filled only upper triangle, but argmax might pick either)
        if i == j:
            # shouldn't happen because diagonal is -inf, but guard anyway
            break
        # If i > j swap to keep merge consistent (merge j into i as original)
        if i > j:
            i, j = j, i

        # Merge cluster j into i (keep order i then j)
        merged[i].extend(merged[j])
        merged_map[i].extend(merged_map[j])
        # Remove the j-th entries (list deletion keeps alignment)
        del merged[j]
        del merged_map[j]
        # continue loop to recompute pairwise p-values on new merged lists

    # build final mapping dict with consecutive keys 0..k-1
    merged_mapping = {k: merged_map[k] for k in range(len(merged_map))}
    return merged, merged_mapping

def find_elbow_idx_by_cutoff(seq, cutoff):
    """
    seq : sequence (list/tuple/1D-array) of numbers, length >= 1
    cutoff : numeric cutoff

    Returns:
      index (int) of the 'after' value where the selected drop occurs (0-based),
      or None if no such drop is found.

    Behavior:
      - Compute drops: drop[i] = seq[i] - seq[i+1]
      - Consider only drops > 0 (actual decreases)
      - Order drop sizes descending; for equal sizes scan indices left->right
      - For each candidate drop index i, check after_value = seq[i+1]
        if after_value <= cutoff: return i+1
    """
    if len(seq) < 2:
        return None

    freq = Counter(seq)
    most_val = max(freq.items(), key=lambda kv: (kv[1], kv[0]))[0]
    
    # compute drops
    drops = [seq[i] - seq[i + 1] for i in range(len(seq) - 1)]

    # collect unique positive drop sizes, sorted descending
    pos_sizes = sorted({d for d in drops if d > 0}, reverse=True)
    if not pos_sizes:
        return most_val, None

    # iterate drop sizes from largest to smallest
    for size in pos_sizes:
        # left-to-right indices where drop equals this size
        for i, d in enumerate(drops):
            if d == size:
                after_val = seq[i + 1]
                if after_val <= cutoff:
                    return most_val, i + 1

def composite_score(x, y):
    """
    Composite = |δ| * (|HL| / hl_scale), clipped to [0,1].
    hl_scale: a characteristic shift you consider “large” (e.g. max observed HL or a domain threshold).
    """
    if len(x) == 0 or len(y) == 0:
        return math.nan
    δ = abs(cliff_delta(x, y))         # in [0,1]
    hl = abs(hodges_lehmann(x, y))     # in data units
    return δ * hl

def hodges_lehmann(x, y):
    """Return the median of all x_i - y_j differences."""
    diffs = [xi - yj for xi in x for yj in y]
    return abs(np.median(diffs))

def cliff_delta(g1, g2):
    nx, ny = len(g1), len(g2)
    greater = sum(x > y for x in g1 for y in g2)
    less    = sum(x < y for x in g1 for y in g2)
    delta = (greater - less) / (nx*ny)
    return abs(delta)

def safe_kde(x, **kwargs):
    x = np.asarray(x)
    # need at least two unique points
    if x.size < 2 or np.all(x == x.flat[0]):
        return None
    return gaussian_kde(x, **kwargs)

def compute_all_ad_pairs(groups):
    """
    Compute t-test p-values for every unique pair.
    Returns a DataFrame with columns ['group1','group2','p_value'].
    """
    records = []
    for (i, g1), (j, g2) in combinations(enumerate(groups), 2):
        p = composite_score(g1, g2)
        records.append({'group1': i, 'group2': j, 'p_value': p})
    return pd.DataFrame(records)

def merge_with_threshold(groups, df_pairs, alpha):
    n = len(groups)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    # Merge groups based on threshold
    for _, row in df_pairs.iterrows():
        if row.p_value <= alpha:
            union(int(row.group1), int(row.group2))

    # Build connected components
    comps = defaultdict(list)
    for i in range(n):
        comps[find(i)].append(i)

    # Merge actual lists + track mapping
    merged = []
    merged_map = {}  # key: merged group index, value: original group indices

    for merged_idx, idxs in enumerate(comps.values()):
        buf = []
        for idx in idxs:
            buf.extend(groups[idx])
        merged.append(buf)
        merged_map[merged_idx] = idxs  # store original indices

    return merged, merged_map

def get_edges_from_json(json_file):
    # Dictionary to store edges
    edge_dict = {}
    # Load the JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)
    # Assuming edges are stored under the 'edges' key
    edges = data.get('edges', [])
    # Iterate over the edges and add them to the dictionary
    for edge in edges:
        source = edge['source']
        target = edge['target']
        if source in edge_dict:
            edge_dict[source].append(target)
        else:
            edge_dict[source] = [target]
    return edge_dict

def cluster_numbers(arr, threshold):
    # Sort the numbers
    arr.sort()
    # Initialize the first group and the list of clusters
    clusters = []
    current_cluster = [arr[0]]
    # Iterate through the sorted array and form clusters
    for i in range(1, len(arr)):
        # If the difference is within the threshold, add to the current cluster
        if arr[i] - arr[i - 1] <= threshold:
            current_cluster.append(arr[i])
        else:
            # Otherwise, finish the current cluster and start a new one
            clusters.append(current_cluster)
            current_cluster = [arr[i]]
    # Append the last cluster
    clusters.append(current_cluster)
    return clusters

def find_corresponding_letters(numbers, A, B):
    result = []
    for number in numbers:
        if number in B:
            index = B.index(number)
            result.append(A[index])
        else:
            result.append(None)
    return result

def get_edges_from_json(json_file):
    # Dictionary to store edges
    edge_dict = {}

    # Load the JSON data
    with open(json_file, 'r') as file:
        data = json.load(file)

    # Assuming edges are stored under the 'edges' key
    edges = data.get('edges', [])

    # Iterate over the edges and add them to the dictionary
    for edge in edges:
        source = edge['source']
        target = edge['target']
        if source in edge_dict:
            edge_dict[source].append(target)
        else:
            edge_dict[source] = [target]  # Store the edge as a tuple key

    return edge_dict

def assign_component(x, groups, mapping, original_comp):
    """
    Assign a batch of sample values x to the best reference group in `groups`.
    """
    if any(set(sublist) == set(x) for sublist in original_comp) and len(x) > 1:
        for keys in mapping:
            if [i for i, sublist in enumerate(original_comp) if Counter(sublist) == Counter(x)][0] in mapping[keys]:
                return keys
            else:
                pass
    else:
        distance_x_avg = []
        for each_subgroup in groups:
            distance_x_avg.append(abs(np.mean(each_subgroup) - np.mean(x)))
        idx_best = distance_x_avg.index(min(distance_x_avg))
    return idx_best

def select_and_convert_gmm_aicc(Samples_Dic, TPM_values_for_gene, Batch_values_for_gene, Factor_name_, max_components=8):
    values = [x for sublist in TPM_values_for_gene for x in sublist]
    batch_labels = [x for sublist in Batch_values_for_gene for x in sublist]

    ######################################## Calculate the variability between replicates ########################################
    replicate_variability = []
    for each_replicate_values in TPM_values_for_gene:
        replicate_variability.append(0 if len(each_replicate_values) < 2 else sum(abs(x-y) for i,x in enumerate(each_replicate_values) for y in each_replicate_values[i+1:]) / (len(each_replicate_values)*(len(each_replicate_values)-1)/2))
    ######################################## Calculate the variability between replicates ########################################

    ########################################## Use elbow method on the Cliff delta ###############################################
    clustering_sensitivity = 100 # the stepsize by which the group will be dissected 
    TPM_values_for_gene_cleaned = [vector for vector in TPM_values_for_gene if len(vector) > 1]
    df_pairs = compute_all_ad_pairs(TPM_values_for_gene_cleaned)
    
    #print(Factor_name_, flush=True)
    #print(df_pairs)
    number_of_clusters = [len(merge_with_threshold(TPM_values_for_gene_cleaned, df_pairs, test_alpha)[0]) for test_alpha in np.linspace(df_pairs['p_value'].min(), df_pairs['p_value'].max(), clustering_sensitivity, endpoint=False)]
    #print('num of clusters: ', number_of_clusters, flush=True)
    most_val, idx_of_elbow = find_elbow_idx_by_cutoff(number_of_clusters, 8)
    print('idx_of_elbow: ', idx_of_elbow, 'most_val: ', most_val, 'linspace: ', np.linspace(df_pairs['p_value'].min(), df_pairs['p_value'].max(), clustering_sensitivity, endpoint=False), flush=True)
    if idx_of_elbow == None:
        optimal_alpha = df_pairs['p_value'].min() # no elbow point, separate everything.
    elif most_val == 1:
        if number_of_clusters.count(1) >= len(number_of_clusters)-1:
            optimal_alpha = np.linspace(df_pairs['p_value'].min(), df_pairs['p_value'].max(), clustering_sensitivity, endpoint=False)[number_of_clusters.index(1)]
        else:
            #optimal_alpha = np.linspace(df_pairs['p_value'].min(), df_pairs['p_value'].max(), 100, endpoint=False)[idx_of_first_8]
            optimal_alpha = np.linspace(df_pairs['p_value'].min(), df_pairs['p_value'].max(), clustering_sensitivity, endpoint=False)[idx_of_elbow]
    else:
        #optimal_alpha = np.linspace(df_pairs['p_value'].min(), df_pairs['p_value'].max(), 100, endpoint=False)[idx_of_first_8]
        optimal_alpha = np.linspace(df_pairs['p_value'].min(), df_pairs['p_value'].max(), clustering_sensitivity, endpoint=False)[idx_of_elbow]
    optimal_clusters, optimal_mapping = merge_with_threshold(TPM_values_for_gene_cleaned, df_pairs, optimal_alpha)
    #print('after hodges_lehmann len(optimal_clusters): ', len(optimal_clusters), '\n', flush=True)
    filtered_clusters = []
    filtered_mapping = {}
    for new_idx, old_idx in enumerate(range(len(optimal_clusters))):
        cluster = optimal_clusters[old_idx]
        if len(cluster) > 1:
            filtered_clusters.append(cluster)
            filtered_mapping[len(filtered_clusters)-1] = optimal_mapping[old_idx]
    optimal_clusters = filtered_clusters
    optimal_mapping = {new_idx: filtered_mapping[old_key] for new_idx, old_key in enumerate(sorted(filtered_mapping.keys()))}
    optimal_clusters, optimal_mapping = teset_and_merge_welch(optimal_clusters, optimal_mapping, sum(sorted(replicate_variability)[-3:]) / 3)
    optimal_mapping = {new_idx: optimal_mapping[old_key] for new_idx, old_key in enumerate(sorted(optimal_mapping.keys()))}
    ########################################## Use elbow method on the Cliff delta ###############################################
    
    best_n = len(optimal_clusters)
    medians = [np.median(c) for c in optimal_clusters]
    means = [np.mean(c) for c in optimal_clusters]
    covs     = [np.std(c)    for c in optimal_clusters]
    weights  = [len(c)/len(values) for c in optimal_clusters]
    order = np.argsort(means)
    optimal_clusters = [optimal_clusters[i] for i in order]
    optimal_mapping = {new_idx: optimal_mapping[old_idx] for new_idx, old_idx in enumerate(order)}
    medians = [medians[i]          for i in order]
    covs = [covs[i]             for i in order]
    weights = [weights[i]          for i in order]
    means = [means[i]          for i in order]
    # map each value to its component’s mean
    outfile = open('./result/Steady_state_count.txt', 'a')
    outfile.write(Factor_name_+'\t'+str(best_n)+'\t')
    converted_GMM_TPM = []
    converted_GMM_std = []
    for idx, each_TPMs in enumerate(TPM_values_for_gene):
        if each_TPMs == []:
            converted_GMM_TPM.append([0])
            converted_GMM_std.append([0])
        else:
            component_label = assign_component(each_TPMs, optimal_clusters, optimal_mapping, TPM_values_for_gene_cleaned)
            converted_GMM_TPM.append(medians[component_label])
            converted_GMM_std.append(covs[component_label])
            outfile.write(sorted(Samples_Dic.keys())[idx] + ':' + str(component_label) + '\t')
    outfile.write('\n')
    # 5) Make plots with a biological replicates shown
    replicate_index = {}
    replicates_to_plot = ['SFL1', 'MET32', 'UME6', 'INO4', 'PHD1', 'PUT3', 'AFT2', 'SKN7', 'WT']
    replicates_to_plot = list(Samples_Dic.keys())
    for each in replicates_to_plot:
        replicate_index[each] = get_sample_index_in_TPM_list(Samples_Dic, TPM_values_for_gene, each)
    #print('result: ', best_n, '\n')

    pairs = list(itertools.combinations(range(len(optimal_clusters)), 2))
    pvals = []
    results = []

    replicate_test_p = []
    if len(optimal_clusters) > 1:
        for i, j in pairs:
            t_stat, p_val = ttest_ind(optimal_clusters[i], optimal_clusters[j], equal_var=False)  # Welch's t-test
            pvals.append(p_val)
            results.append((str(i), str(j), t_stat, p_val))

        # Correct for multiple comparisons (Holm recommended)
        reject, pvals_corr, _, _ = multipletests(pvals, method="holm")

        # Print results
        for (g1, g2, t, p), p_corr, r in zip(results, pvals_corr, reject):
            replicate_test_p.append(p_corr)
    else:
        pass
    outfile.close()
    
    plot_GMM_distribution(np.array(values), batch_labels, optimal_clusters, optimal_mapping, TPM_values_for_gene_cleaned, best_n, replicate_index, float(max(replicate_test_p)) if replicate_test_p else math.nan, Factor_name=Factor_name_)

    for each_i in range(0, len(converted_GMM_TPM)):
        if isinstance(converted_GMM_TPM[each_i], list):
            converted_GMM_TPM[each_i] = converted_GMM_TPM[each_i][0]
            converted_GMM_std[each_i] = converted_GMM_std[each_i][0]
        elif isinstance(converted_GMM_TPM[each_i], np.ndarray):
            converted_GMM_TPM[each_i] = converted_GMM_TPM[each_i].item()
            converted_GMM_std[each_i] = converted_GMM_std[each_i].item()
        else:
            pass

    return best_n, np.array(converted_GMM_TPM), np.array(converted_GMM_std)

def plot_GMM_distribution(values_, batch_labels, merged_clusters_, mapping_, original_comp_, best_n_, replicate_index, max_p_value, Factor_name='', outname='GMM'):
    fig, axes = plt.subplots(2, 1, figsize=(9, 10), sharex=True)
    Hexcode_colors = ['#FF0000', '#FF7F00', '#FFFF00', '#00FF00', '#0000FF', '#4B0082', '#8B00FF', '#FF00FF']
    axes[0].hist(values_, bins=100, density=True, alpha=0.3, color='gray')
    # Overlay each component PDF
    for cluster_i, cluster_data in enumerate(merged_clusters_):
        cluster_data = np.array(cluster_data)
        kde = safe_kde(cluster_data)
        if kde is None:
            axes[0].axhline(cluster_data.flat[0], color=Hexcode_colors[cluster_i])
        else:
            xs = np.linspace(cluster_data.min(), cluster_data.max(), 2*len(cluster_data))
            weight = len(cluster_data)/len(values_)
            axes[0].plot(xs, weight*kde(xs), color=Hexcode_colors[cluster_i])
    axes[0].text(0.95, 0.95, "Histogram of {} mRNA-seq samples\nMax ttest pvalue: {}".format(len(values_), format(float(max_p_value), ".2e")), horizontalalignment="right", verticalalignment="top", transform=axes[0].transAxes, fontsize=10, color="black")
    axes[0].tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    axes[0].set_ylabel("Normalized sample freq")
    axes[0].set_ylim()
    axes[0].set_title(fr"Bottom Up Model Fit for $\it{{{Factor_name}}}$ (n_component={best_n_})")
    df_combined = []
    df_all = []
    df_mean = []
    if len(replicate_index) > 0:
        for key in replicate_index:
            safe_key = key.replace('_', r'\_')
            x_scatter = []
            batch_scatter = []
            for j in range(replicate_index[key][0], replicate_index[key][1]):
                x_scatter.append(values_[j])
                batch_scatter.append(batch_labels[j])
            if len(x_scatter) != 0:
                df_temp = pd.DataFrame({
                    'value': x_scatter,
                    'batch' : batch_scatter,
                    'group': [fr"$\it{{{safe_key}}}$"] * len(x_scatter)
                })
                df_all.append(df_temp)
                df_mean.append(np.mean(x_scatter))
            else:
                pass
        df_all = [b for _, b in sorted(zip(df_mean, df_all), key=lambda x: x[0], reverse=False)]
    else:
        pass
    df_combined = pd.concat(df_all)
    unique_groups = df_combined['group'].unique()

    # Violin plot on the second subplot (with x ticks)
    for each_df_i in range(0, len(df_all)):
        component_label_ = assign_component(df_all[each_df_i]['value'].tolist(), merged_clusters_, mapping_, original_comp_)
        #sns.violinplot(x='value', y='group', data=df_all[each_df_i], split=True, ax=axes[1], orient='h', inner=None, linewidth=0, color=Hexcode_colors[component_label_], alpha=0.5)
    # Loop through groups and plot dots just below
    for i, group in enumerate(unique_groups):
        subset = df_combined[df_combined['group'] == group]
        y_pos = np.full(len(subset), i + 1)  # offset to bottom
        component_label_ = assign_component(subset['value'].tolist(), merged_clusters_, mapping_, original_comp_)
        axes[1].scatter(subset['value'], y_pos, color=Hexcode_colors[component_label_], s=45, alpha=1)
        axes[1].hlines(y=y_pos, xmin=min(subset['value']), xmax=max(subset['value']), colors=Hexcode_colors[component_label_], linestyles='-', linewidth=0.5)

    axes[1].tick_params(axis='y', which='both', left=False, labelleft=False)
    axes[1].set_ylabel("Normalized sample freq")
    axes[1].set_xlabel("TMM normalized read counts")
    axes[1].set_ylabel("Condition (replicates shown)")
    plt.tight_layout()
    os.makedirs('./result/GMM_figures', exist_ok=True)
    plt.savefig('./result/GMM_figures/AIC/{}_{}.jpg'.format(Factor_name, outname), dpi=300)
    plt.close()
    return

def get_sample_index_in_TPM_list(Samples_Dic, TPM_values_for_gene, Factor):
    sample_index = 0
    for i in range(0, sorted(Samples_Dic.keys()).index(Factor)):
        sample_index = sample_index + len(TPM_values_for_gene[i])
    sample_span = [sample_index, sample_index+len(TPM_values_for_gene[sorted(Samples_Dic.keys()).index(Factor)])]
    return sample_span

'''Step 1: obtain mRNA steady states'''
########################################################################################################################################################################
infile = open(args.gene_of_interest, 'r')

GRN_genes = []
ssTFs_len = {}
ssTFs_SGDID = {}
TMM_reads_df = pd.read_csv(args.RNAseq_data, sep='\t', index_col=0)
for line in infile:
    if line.split()[1] in TMM_reads_df.columns:
        ssTFs_len[line.split()[0]] = line.split()[2]
        ssTFs_SGDID[line.split()[1]] = line.split()[0]
        GRN_genes.append(line.split()[0])
    else:
        continue
infile.close()

mRNA_steady_states = {}
TMM_reads_selected = TMM_reads_df[[k for k in ssTFs_SGDID.keys() if k in TMM_reads_df.columns]]
for row_name, row_ in TMM_reads_selected.iterrows():
    dic_for_the_row = row_.to_dict()
    temp_converted_dic = {}
    for keys_ in dic_for_the_row:
        temp_converted_dic[ssTFs_SGDID[keys_]] = dic_for_the_row[keys_]
    mRNA_steady_states[row_name] = temp_converted_dic

Samples_Dic = {}
for each in list({item.split('_')[0] for item in list(mRNA_steady_states.keys())}):
    Samples_Dic[each] = []
    for each_sample in list(mRNA_steady_states.keys()):
        if each_sample.split('_')[0] == each:
            Samples_Dic[each].append(each_sample)
        else:
            pass
print('Sample_Dic: ', Samples_Dic, flush=True)
Column_order = GRN_genes

# output the Salmon TPMs average for each genotype.
outfile = open('./data/GRN_ssTFs_Salmon_SteadyStates_2025.txt', 'a')
for each_sample in sorted(Samples_Dic.keys()):
    line_to_write = ''
    if each_sample in Column_order:
        line_to_write = line_to_write + '{}\t'.format(Column_order.index(each_sample))
    else:
        line_to_write = line_to_write + '-1\t'
    for each_column in Column_order:
        temp_TPM_list = []
        for each_replicate in Samples_Dic[each_sample]:
            if each_column not in mRNA_steady_states[each_replicate]:
                temp_TPM_list.append(0)
            else:
                temp_TPM_list.append(float(mRNA_steady_states[each_replicate][each_column]))
        line_to_write = line_to_write + '{}\t'.format('_'.join(map(str,temp_TPM_list)))
    line_to_write = line_to_write[:-1] + '\n'
    outfile.write(line_to_write)
outfile.close()

TPM_matrix_T = []
std_matrix_T = []

row_header = []
for each_sample in sorted(Samples_Dic.keys()):
    if each_sample in Column_order:
        row_header.append(Column_order.index(each_sample))
    else:
        row_header.append(-1)
TPM_matrix_T.append(row_header)
std_matrix_T.append(row_header)

for each_column in Column_order:
    TPM_values_for_gene = []
    Batch_values_for_gene = []
    #print(each_column)
    for each_sample in sorted(Samples_Dic.keys()): # this is the row order
        TPM_values_for_replicates = []
        Batch_values_for_replicates = []
        if each_column == each_sample: # The gene has been knock-outed
            pass
        else:
            for each_replicate in Samples_Dic[each_sample]:
                TPM_values_for_replicates.append(float(mRNA_steady_states[each_replicate][each_column]))
                Batch_values_for_replicates.append(each_replicate.split('_')[2])
        TPM_values_for_gene.append(TPM_values_for_replicates)
        Batch_values_for_gene.append(Batch_values_for_replicates)
    #print('===>', each_column)
    best_n, convertion_Dic, std_Dic = select_and_convert_gmm_aicc(Samples_Dic, TPM_values_for_gene, Batch_values_for_gene, each_column)
    TPM_matrix_T.append(list(convertion_Dic))
    std_matrix_T.append(list(std_Dic))
    #print(len(convertion_Dic), 'convertion_Dic: ', convertion_Dic)
    #print('\n')

TPM_matrix = [ list(row) for row in zip(*TPM_matrix_T) ]
std_matrix = [ list(row) for row in zip(*std_matrix_T) ]
#print(TPM_matrix)
with open("./data/GRN_ssTFs_Salmon_SteadyStates_2025_discrete.txt", "a", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(TPM_matrix)
with open("./data/GRN_ssTFs_Salmon_SteadyStates_2025_std.txt", "a", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(std_matrix)
########################################################################################################################################################################


'''Step 2: obtain gene length'''
########################################################################################################################################################################
outfile = open('./data/GRN_ssTFs_Sc_gene_length.txt', 'a')
for each in Column_order:
    if (1+Column_order.index(each)) == len(Column_order):
        outfile.write(ssTFs_len[each])
    else:
        outfile.write(ssTFs_len[each]+'\t')
outfile.close()
########################################################################################################################################################################


'''Step 3: obtain the row and column headers'''
########################################################################################################################################################################
outfile = open('./data/GRN_ssTFs_row_names.txt', 'a')
for each in sorted(Samples_Dic.keys()):
    outfile.write(each+'\t')
outfile.close()

outfile = open('./data/GRN_ssTFs_column_names.txt', 'a')
for each in Column_order:
    outfile.write(each+'\t')
outfile.close()
########################################################################################################################################################################


'''Step 4: obtain TF-DNA'''
########################################################################################################################################################################
PPI_dic = get_edges_from_json(args.TF_DNA_prior)
PPI_matrix_TF_DNA = np.zeros((len(Column_order), len(Column_order)))
for i in range(0, len(Column_order)):
    for j in range(0, len(Column_order)):
        if Column_order[i] in PPI_dic and Column_order[j] in PPI_dic[Column_order[i]]:
            PPI_matrix_TF_DNA[i,j] = 1
        else:
            continue

np.savetxt('./data/GRN_ssTFs_Sc_TF_DNA.txt', PPI_matrix_TF_DNA, delimiter='\t', fmt='%d')
########################################################################################################################################################################


'''Step 5: obtain the union of Motif-based PPI and all-binding sites PPI'''
########################################################################################################################################################################
TF_DNA_all = get_edges_from_json('./data/Rossi_Ruihao_TF_DNA_union_motif_based.json')
TF_DNA = {}
for each_source in TF_DNA_all:
    for each_target in TF_DNA_all[each_source]:
        if each_source in ssTFs_len and each_target in ssTFs_len:
            if each_source not in TF_DNA:
                TF_DNA[each_source] = [each_target]
            else:
                TF_DNA[each_source].append(each_target)
        else:
            continue

PPI_dic = get_edges_from_json(args.protein_protein_colocalization_prior)
Motif_PPI_dic = get_edges_from_json(args.protein_protein_colocalization_prior)
#print(PPI_dic)
#print(Motif_PPI_dic, '\n\n')
for keys in Motif_PPI_dic:
    if keys not in PPI_dic:
        PPI_dic[keys] = Motif_PPI_dic[keys]
    else:
        for factors in Motif_PPI_dic[keys]:
            if factors not in PPI_dic[keys]:
                PPI_dic[keys].append(factors)
            else:
                pass
#print(PPI_dic)

PPI_matrix = np.zeros((len(Column_order), len(Column_order)))
for i in range(0, len(Column_order)):
    for j in range(0, len(Column_order)):
        if Column_order[i] in PPI_dic and Column_order[j] in PPI_dic[Column_order[i]]:
            PPI_matrix[i,j] = 1
        else:
            continue

# obtain genome annotation
inputfile = open(args.annotation, 'r')
Genome_anno = {}
NDR_or_NFR = {}
for line in inputfile:
    if line.split()[1] == '+':
        Genome_anno[line.split()[4]] = [line.split()[0], int(line.split()[2]), int(line.split()[3])]
        if line.split()[4] in ssTFs_SGDID:
            NDR_or_NFR[ssTFs_SGDID[line.split()[4]]] = [int(line.split()[-2]), int(line.split()[-1])]
    elif line.split()[1] == '-':
        Genome_anno[line.split()[4]] = [line.split()[0], int(line.split()[3]), int(line.split()[2])]
        if line.split()[4] in ssTFs_SGDID:
            NDR_or_NFR[ssTFs_SGDID[line.split()[4]]] = [int(line.split()[-2]), int(line.split()[-1])]
    if line.split()[4] in ssTFs_SGDID and ssTFs_SGDID[line.split()[4]] in NDR_or_NFR:
        NDR_or_NFR[ssTFs_SGDID[line.split()[4]]].append(line.split()[0])
inputfile.close()

infile = open(args.YEP_replicate_ID, 'r')
YEP_best_rep = {}
for line in infile:
    YEP_best_rep[line.split()[0].upper()] = line.split()[1]
    if 'RSC1' == line.split()[0].upper():
        YEP_best_rep['AFT1'] = line.split()[1]
    else:
        pass
infile.close()

Promoter_binding_summary = []
for each_source in TF_DNA:
    for each_target in TF_DNA[each_source]:
        if each_target == each_source:
            continue
        else:
            infile = open('./data/YEP_bed/{}_chexmix_filtered_peaks.bed'.format(YEP_best_rep[each_source]), 'r')
            dist_list = {}
            for line in infile:
                if line.split()[0] == NDR_or_NFR[each_target][2]:
                    dist_list[abs(int(line.split()[1])-0.5*abs(NDR_or_NFR[each_target][0] + NDR_or_NFR[each_target][1]))] = int(line.split()[1])
                else:
                    pass
            infile.close()
            if dist_list == {}:
                pass
            else:
                Promoter_binding_summary.append([each_target, each_source, dist_list[min(list(dist_list.keys()))]])
#print(Promoter_binding_summary)

Complexes = {}
for each in Promoter_binding_summary:
    if each[0] not in Complexes:
        Complexes[each[0]] = [[each[1]], [each[2]]]
    else:
        Complexes[each[0]][0].append(each[1])
        Complexes[each[0]][1].append(each[2])

Factor_binding_distance = 50

Complexes_names = {}
for each in Complexes:
    Complexes_names[each] = []
    for i in range(0, len(cluster_numbers(Complexes[each][1], Factor_binding_distance))):
        Complexes_names[each].append(find_corresponding_letters(cluster_numbers(Complexes[each][1], Factor_binding_distance)[i], Complexes[each][0], Complexes[each][1]))
#print(Complexes_names)

LG_String = ''
for each in Column_order:
    LG_String_temp = ''
    if each not in Complexes_names:
        LG_String_temp = ','.join(map(str, range(len(Column_order))))
    else:
        index_converter = {}
        for each_complex in Complexes_names[each]:
            if len(each_complex) == 1:
                pass
            else:
                indexes_for_the_complex = []
                for each_component in each_complex:
                    if each_component in Column_order:
                        indexes_for_the_complex.append(Column_order.index(each_component))
                    else:
                        pass
                for each_index in indexes_for_the_complex:
                    index_converter[each_index] = min(indexes_for_the_complex)
        for i in range(0, len(Column_order)):
            if i not in index_converter:
                LG_String_temp = LG_String_temp + str(i) + ','
            else:
                LG_String_temp = LG_String_temp + str(index_converter[i]) + ','
        LG_String_temp = LG_String_temp[:-1]
    LG_String = LG_String + ',' + LG_String_temp

#print('\n', LG_String[1:], '\n')
outfile = open('./data/GRN_ssTFs_Sc_LG.txt', 'a')
outfile.write(LG_String[1:])
outfile.close()
########################################################################################################################################################################


'''Step 6: obtain promoter strength'''
########################################################################################################################################################################
df_YEP = pd.read_excel('./data/41586_2021_3314_MOESM3_ESM.xlsx', sheet_name='Supplementary Data 1')

selected_columns = df_YEP[['Systematic ID', 'Chrom', 'Experiment_Left', 'Experiment_Right']]
ssTFs_NFRs_df = selected_columns[df_YEP['Systematic ID'].isin(ssTFs_SGDID.keys())]

ssTFs_NFRs_dic = {}
chr_convert = {'chr1': 'chrI', 'chr2': 'chrII', 'chr3': 'chrIII', 'chr4': 'chrIV', 'chr9': 'chrIX', 'chr5': 'chrV', 'chr6': 'chrVI',
              'chr7': 'chrVII', 'chr8': 'chrVIII', 'chr10': 'chrX', 'chr11': 'chrXI', 'chr12': 'chrXII', 'chr13': 'chrXIII', 'chr14': 'chrXIV',
              'chr15': 'chrXV', 'chr16': 'chrXVI'}
for index, row in ssTFs_NFRs_df.iterrows():
    ssTFs_NFRs_dic[row.values.tolist()[0]] = [chr_convert[row.values.tolist()[1]], int(row.values.tolist()[2]), int(row.values.tolist()[3])]

df_CAGE = pd.read_excel('./data/Supplemental_Data_S5_S8.xlsx', header=12, sheet_name='Data S5')
selected_columns = df_CAGE[['gene/transcript', 'chr', 'start', 'end', 'YPD.tpm']]

ssTFs_CAGE_df = selected_columns[df_CAGE['gene/transcript'].isin(ssTFs_SGDID.keys())].dropna(subset=['YPD.tpm'])

ssTFs_CAGE_dic = {}
for index, row in ssTFs_CAGE_df.iterrows():
    if row.values.tolist()[0] not in ssTFs_CAGE_dic:
        ssTFs_CAGE_dic[row.values.tolist()[0]] = [row.values.tolist()[1], int(row.values.tolist()[2]), int(row.values.tolist()[3]), float(row.values.tolist()[4])]
    else:
        if ssTFs_CAGE_dic[row.values.tolist()[0]][3] > float(row.values.tolist()[4]):
            pass
        else:
            ssTFs_CAGE_dic[row.values.tolist()[0]] = [row.values.tolist()[1], int(row.values.tolist()[2]), int(row.values.tolist()[3]), float(row.values.tolist()[4])]

outfile = open('./data/GRN_ssTFs_Sc_promoter_strength.txt', 'a')
for each_sample in sorted(Samples_Dic.keys()):
    for each in Column_order:
        if (1+Column_order.index(each)) == len(Column_order):
            #print(each_sample, each, ssTFs_CAGE_dic[[key for key, value in ssTFs_SGDID.items() if value == each][0]][-1])
            outfile.write(str(ssTFs_CAGE_dic[[key for key, value in ssTFs_SGDID.items() if value == each][0]][-1]))
        else:
            #print(each_sample, each, ssTFs_CAGE_dic[[key for key, value in ssTFs_SGDID.items() if value == each][0]][-1])
            outfile.write(str(ssTFs_CAGE_dic[[key for key, value in ssTFs_SGDID.items() if value == each][0]][-1])+'\t')
    outfile.write('\n')
outfile.close()
########################################################################################################################################################################


'''Step 7: obtain initial GRN from TF-DNA and TPM'''
########################################################################################################################################################################
GMM_TPM_levels = np.array(TPM_matrix)
KO_index_GMM = GMM_TPM_levels[:, 0]

TF_DNA_all = get_edges_from_json(args.TF_DNA_prior)
PPI_dic = {}
for each_source in TF_DNA_all:
    for each_target in TF_DNA_all[each_source]:
        if each_source in ssTFs_len and each_target in ssTFs_len:
            if each_source not in PPI_dic:
                PPI_dic[each_source] = [each_target]
            else:
                PPI_dic[each_source].append(each_target)
        else:
            continue
        
PPI_matrix_TF_DNA = np.zeros((len(Column_order), len(Column_order)))
for i in range(0, len(Column_order)):
    for j in range(0, len(Column_order)):
        if Column_order[i] in PPI_dic and Column_order[j] in PPI_dic[Column_order[i]]:
            PPI_matrix_TF_DNA[i,j] = 1
        else:
            continue

expression_data = GMM_TPM_levels[:, 1:].T
num_arrays = expression_data.shape[0]
nmi_matrix = np.zeros((num_arrays, num_arrays))
TPM_driven_matrix = np.zeros((num_arrays, num_arrays))
union_matrix = np.zeros((num_arrays, num_arrays))
nmi_cutoff = 0.5

for i in range(num_arrays):
    for j in range(num_arrays):
        print('expression_data[i]: ', expression_data[i])
        nmi_matrix[i, j] = normalized_mutual_info_score(list(map(str, expression_data[i])), list(map(str, expression_data[j])), average_method='arithmetic')
        if nmi_matrix[i, j] < 0.5 and PPI_matrix_TF_DNA[i, j] == 0:
            TPM_driven_matrix[i, j] = 0
        else:
            tau, p_ = kendalltau(expression_data[i], expression_data[j])
            if tau >= 0.25:
                TPM_driven_matrix[i, j] = 1
            elif tau <= -0.25:
                TPM_driven_matrix[i, j] = 2
            else:
                TPM_driven_matrix[i, j] = 0

for i in range(num_arrays):
    for j in range(num_arrays):
        if int(PPI_matrix_TF_DNA[i, j]) == 0 and int(TPM_driven_matrix[i, j]) == 0:
            pass
        else:
            union_matrix[i, j] = 1
np.savetxt('./data/GRN_ssTFs_Sc_TF_DNA_TPM_union.txt', union_matrix, delimiter='\t', fmt='%d')
            
TPM_driven_matrix = TPM_driven_matrix.astype(int)
String_TF_DNA = ''
for i in range(0, len(TPM_driven_matrix)):
    for j in range(0, len(TPM_driven_matrix[i])):
        String_TF_DNA = String_TF_DNA + str(TPM_driven_matrix[i][j])
outfile = open('./data/GRN_ssTFs_Sc_initial_condition.txt', 'a')
outfile.write(String_TF_DNA)
outfile.close()
########################################################################################################################################################################
