from utility_functions import *
from distance_functions import *
import numpy as np
import random
import itertools
import time

def find_duplicate_indices(lst):
    # Create a dictionary to store the indices for each number.
    indices = {}
    for idx, num in enumerate(lst):
        if num in indices:
            indices[num].append(idx)
        else:
            indices[num] = [idx]

    # Filter out numbers that appear only once.
    duplicates = {num: idxs for num, idxs in indices.items() if len(idxs) > 1}
    return duplicates

def biased_random(min_val, max_val, lambda_value=0.5):
    """
    Returns a random number between min_val and max_val (inclusive),
    with min_val having the highest probability and max_val the lowest.
    Probability decreases exponentially.
    """
    # Create range of numbers
    nums = np.arange(min_val, max_val + 1)

    # Compute exponential weights
    weights = np.exp(-lambda_value * (nums - min_val))  # Decay from min_val

    # Normalize to sum to 1 (convert weights to probabilities)
    probabilities = weights / weights.sum()

    # Sample from the distribution
    return np.random.choice(nums, p=probabilities)

def Mutation_LG_Simple(_HammingDistance, LG_Matrix, index_of_diff_gene):
    """ Randomly break a large complex into two smaller complexes """
    num_of_complexes = 0
    complexes_list = find_duplicate_indices(LG_Matrix)
    for each in complexes_list:
        num_of_complexes = num_of_complexes + 1
    _HammingDistance = biased_random(min(_HammingDistance, num_of_complexes), num_of_complexes)
    while _HammingDistance > 0:
        complexes_on_gene = find_duplicate_indices(LG_Matrix)
        if complexes_on_gene == {}:
            pass
        else:
            complex_to_mutate = np.random.choice(list(complexes_on_gene.keys()))
            subunit_to_mutate = np.random.choice(complexes_on_gene[complex_to_mutate])
            LG_Matrix[subunit_to_mutate] = 1 + max(LG_Matrix)
            _HammingDistance = _HammingDistance - 1
    return Sort_LG([LG_Matrix])[0]

def HammingMutation_1(n, strings, TF_DNA, diff_genes, max_attempts=1000):
    """
    Find a random valid string that obeys TF_DNA and is not in the tested strings.
    Instead of generating all possibilities, it tries up to max_attempts times,
    and uses frequency information from the tested strings to bias candidate generation.
    """
    # Directly returns False if all possibile candidates already exist in strings
    if 2**sum(TF_DNA) < len(strings):
        return False
    else:
        pass

    # For each position i, if TF_DNA[i] is 1, allowed digits are '1' or '2'; if 0, only '0'.
    # allowed_choices = [("12" if tf == 1 else "0") for tf in TF_DNA]
    allowed_choices = [
        ("12" if tf == 1 else "0") if i in diff_genes else "0"
        for i, tf in enumerate(TF_DNA)
    ]
    
    # Compute frequency counts for each allowed digit in each position.
    # Initialize a list of dictionaries, one for each position.
    freq = [{} for _ in range(n)]
    for s in strings:
        for i, ch in enumerate(s):
            if ch in allowed_choices[i]:
                freq[i][ch] = freq[i].get(ch, 0) + 1
    
    # For each position, create a biased probability distribution.
    # We want to favor digits that appear less frequently.
    biased_choices = []
    for i in range(n):
        allowed = list(allowed_choices[i])
        # Count frequency for each allowed digit (defaulting to 0 if not seen)
        counts = [freq[i].get(d, 0) for d in allowed]
        # Inverse frequency (adding 1 to avoid division by zero)
        inv_freq = [1.0 / (count + 1) for count in counts]
        total = sum(inv_freq)
        probabilities = [w / total for w in inv_freq]
        biased_choices.append((allowed, probabilities))
    
    for _ in range(max_attempts):
        # Generate a candidate string using the biased choices
        candidate_chars = [
            random.choices(allowed, weights=probs)[0] 
            for allowed, probs in biased_choices
        ]
        candidate = ''.join(candidate_chars)
        if candidate not in strings:
            return candidate
    return False

def HammingMutation_2(n, strings, TF_DNA, diff_genes, max_attempts=1000):
    '''
    Find a random valid string that partially obeys TF_DNA and is not in the tested strings.
    Instead of generating all possibilities, it tries up to max_attempts times,
    and uses frequency information from the tested strings to bias candidate generation.
    '''
    # Directly returns False if all possibile candidates already exist in strings
    if 3**sum(TF_DNA) < len(strings):
        return False
    else:
        pass
    
    # For each position i, if TF_DNA[i] is 1, allowed digits are '1' or '2'; if 0, only '0'.
    #allowed_choices = [("120" if tf == 1 else "0") for tf in TF_DNA]
    allowed_choices = [
        ("120" if tf == 1 else "0") if i in diff_genes else "0"
        for i, tf in enumerate(TF_DNA)
    ]
    
    # Compute frequency counts for each allowed digit in each position.
    # Initialize a list of dictionaries, one for each position.
    freq = [{} for _ in range(n)]
    for s in strings:
        for i, ch in enumerate(s):
            if ch in allowed_choices[i]:
                freq[i][ch] = freq[i].get(ch, 0) + 1
    
    # For each position, create a biased probability distribution.
    # We want to favor digits that appear less frequently.
    biased_choices = []
    for i in range(n):
        allowed = list(allowed_choices[i])
        # Count frequency for each allowed digit (defaulting to 0 if not seen)
        counts = [freq[i].get(d, 0) for d in allowed]
        # Inverse frequency (adding 1 to avoid division by zero)
        inv_freq = [1.0 / (count + 1) for count in counts]
        total = sum(inv_freq)
        probabilities = [w / total for w in inv_freq]
        biased_choices.append((allowed, probabilities))
    
    for _ in range(max_attempts):
        # Generate a candidate string using the biased choices
        candidate_chars = [
            random.choices(allowed, weights=probs)[0] 
            for allowed, probs in biased_choices
        ]
        candidate = ''.join(candidate_chars)
        if candidate not in strings:
            return candidate
    return False

def HammingMutation_3(n, strings, diff_genes, max_attempts=1000, network_sparsity=0.2):
    '''
    Find a random valid string that is not in the tested strings.
    Instead of generating all possibilities, it tries up to max_attempts times,
    and uses frequency information from the tested strings to bias candidate generation.
    '''
    allowed_choices = [
        ("12"+"0"*int(2/network_sparsity-2)) if i in diff_genes else "0"
        for i in range(0, n)
    ]

    for _ in range(max_attempts):
        # Generate a candidate string using the biased choices
        candidate_chars = [
            random.choice(allowed_strings)
            for allowed_strings in allowed_choices
        ]
        candidate = ''.join(candidate_chars)
        if candidate not in strings:
            return candidate
    return False
