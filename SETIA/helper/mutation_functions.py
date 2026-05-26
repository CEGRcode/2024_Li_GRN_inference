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

def HammingMutation_3(n, strings, TF_DNA, diff_genes, max_attempts=1000, network_sparsity=0.2):
    '''
    Find a random valid string that is not in the tested strings.
    Instead of generating all possibilities, it tries up to max_attempts times,
    and uses frequency information from the tested strings to bias candidate generation.
    '''

    if TF_DNA == []:
        TF_DNA = [0] * n

    # Precompute allowed choices per position.
    # For each position i: if TF_DNA[i] is 1, allowed digits are "1" or "2" or "0"; if 0, biased toward "0".
    allowed_choices = [
        ("120" if tf == 1 else "0") if i in diff_genes else "0"
        for i, tf in enumerate(TF_DNA)
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


def HammingMutation_local(current_config, n, strings, TF_DNA, diff_genes,
                           max_attempts=1000, hamming_dist=2, only_TF_DNA=True):
    """
    在当前config的局部邻域搜索未测试过的新config。
    
    与HammingMutation_1/2/3不同：不从头随机生成，而是在current_config
    基础上做hamming距离为1~hamming_dist的局部扰动。
    
    每次mutation：随机选1~hamming_dist个可变位置，改变它们的值。
    ChIP=1的位置可以是0/1/2（允许删边）；ChIP=0的位置保持0。
    
    返回：新config字符串，或False（邻域已穷尽）
    """
    if len(current_config) != n:
        return False

    # 可以被改变的位置：only_TF_DNA=True时只允许ChIP=1，False时允许所有diff_genes
    mutable = [i for i in range(n)
               if i in diff_genes and (not only_TF_DNA or (i < len(TF_DNA) and TF_DNA[i] == 1) or TF_DNA == [])]

    if not mutable:
        return False

    # 对每个可变位置，定义它能变成的值（ChIP=1 → 0/1/2）
    choices = {}
    for i in mutable:
        current_val = current_config[i]
        # 可以变成0（删边）、1（激活）、2（抑制），排除当前值
        choices[i] = [v for v in ['0', '1', '2'] if v != current_val]

    for _ in range(max_attempts):
        # 随机选1~hamming_dist个位置改变
        k = random.randint(1, min(hamming_dist, len(mutable)))
        positions = random.sample(mutable, k)

        candidate = list(current_config)
        for pos in positions:
            candidate[pos] = random.choice(choices[pos])
        candidate_str = ''.join(candidate)

        if candidate_str not in strings:
            return candidate_str

    return False  # 局部邻域已被穷尽，fallback到全局搜索

def HammingMutation_guided(current_config, n, strings, TF_DNA, diff_genes,
                            residuals, expression_data, t1_cache,
                            max_attempts=1000, hamming_dist=2, k=30, only_TF_DNA=True):
    """
    有方向性的局部mutation。
    
    residuals: list of (ss - target) for each condition
               正值 → ss > target → 这个条件需要更多抑制
               负值 → ss < target → 这个条件需要更多激活
    expression_data: 2D array [n_conds x n_genes]，原始表达量
    t1_cache: dict {gene_index: threshold t1}
    
    逻辑：
    1. 找误差大的条件，判断需要抑制还是激活
    2. 找在"需要抑制"条件里Hill高的TF作为repressor候选
       找在"需要激活"条件里Hill高的TF作为activator候选
    3. 从候选里随机选一个，在current_config里加/删/改这条边
    """
    import math

    def hill(X, gi_reg):
        t1 = t1_cache.get(gi_reg, 0)
        if t1 <= 0 or X <= 0: return 0.0
        try:
            r = X / t1
            return r**k / (r**k + 1)
        except:
            return 0.0

    n_conds = len(residuals)
    eps = 0.05  # 误差阈值，小于这个认为该条件已经拟合好

    # 按residual大小找问题条件
    need_rep = [c for c in range(n_conds) if residuals[c] >  eps]  # ss>target，需要抑制
    need_act = [c for c in range(n_conds) if residuals[c] < -eps]  # ss<target，需要激活

    # 可操作的位置：only_TF_DNA=True时只允许ChIP=1，False时允许所有diff_genes
    mutable = [i for i in range(n)
               if i in diff_genes and (not only_TF_DNA or (i < len(TF_DNA) and TF_DNA[i] == 1))]

    if not mutable:
        return False

    # 对每个候选TF，计算它在问题条件里的Hill值
    rep_candidates = []
    act_candidates = []

    for gi_reg in mutable:
        if need_rep:
            h_rep = sum(hill(expression_data[c][gi_reg], gi_reg) for c in need_rep) / len(need_rep)
        else:
            h_rep = 0.0
        if need_act:
            h_act = sum(hill(expression_data[c][gi_reg], gi_reg) for c in need_act) / len(need_act)
        else:
            h_act = 0.0

        # rep候选：在need_rep条件里Hill高，且比在need_act条件里高
        if h_rep > 0.2 and h_rep > h_act:
            rep_candidates.append((gi_reg, h_rep))
        # act候选：在need_act条件里Hill高，且比在need_rep条件里高
        if h_act > 0.2 and h_act > h_rep:
            act_candidates.append((gi_reg, h_act))

    # 按Hill值降序，给高质量候选更高的采样权重
    rep_candidates.sort(key=lambda x: -x[1])
    act_candidates.sort(key=lambda x: -x[1])

    # 如果没有明确方向的候选，退化为random local mutation
    if not rep_candidates and not act_candidates:
        return HammingMutation_local(current_config, n, strings, TF_DNA, diff_genes,
                                     max_attempts, hamming_dist)

    for _ in range(max_attempts):
        candidate = list(current_config)

        # 随机决定这次mutation的操作
        # 如果同时有两个方向的需求，各50%概率选
        # 按residual严重度加权选择action：
        # need_rep严重度高 → 大概率选add_rep；need_act严重度高 → 大概率选add_act
        rep_severity = sum(r for r in residuals if r > eps)
        act_severity = sum(-r for r in residuals if r < -eps)
        total_severity = rep_severity + act_severity + 1e-9

        action_pool = []
        action_weights = []
        if rep_candidates:
            action_pool.append('add_rep')
            action_weights.append(rep_severity / total_severity)
        if act_candidates:
            action_pool.append('add_act')
            action_weights.append(act_severity / total_severity)
        # 也允许删除当前的边（防止过度堆积），固定权重0.1
        current_edges = [(i, current_config[i]) for i in mutable if current_config[i] != '0']
        if current_edges:
            action_pool.append('remove_edge')
            action_weights.append(0.1)

        if not action_pool:
            continue

        # 归一化权重
        total_w = sum(action_weights)
        action_weights = [w / total_w for w in action_weights]
        action = random.choices(action_pool, weights=action_weights)[0]

        if action == 'add_rep' and rep_candidates:
            # 按Hill值加权采样，优先选Hill值高的repressor
            weights = [x[1] for x in rep_candidates]
            gi_reg = random.choices([x[0] for x in rep_candidates], weights=weights)[0]
            candidate[gi_reg] = '2'  # 设为repressor

        elif action == 'add_act' and act_candidates:
            weights = [x[1] for x in act_candidates]
            gi_reg = random.choices([x[0] for x in act_candidates], weights=weights)[0]
            candidate[gi_reg] = '1'  # 设为activator

        elif action == 'remove_edge' and current_edges:
            # 随机删一条现有的边
            gi_reg, _ = random.choice(current_edges)
            candidate[gi_reg] = '0'

        candidate_str = ''.join(candidate)
        if candidate_str not in strings and candidate_str != current_config:
            return candidate_str

    return False
