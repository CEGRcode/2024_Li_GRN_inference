from utility_functions import *

def GetAttractorDistance(mRNA1, WT_mRNA, transcriptionalprofilemax, transcriptionalprofilemin):
    '''Compute Hamming distances for matrices with flattening (for logic gate matrices)'''
    outdistance = 0.0
    for y in range(0, len(WT_mRNA)):
        outdistance = outdistance + abs((mRNA1[y]-WT_mRNA[y])/(transcriptionalprofilemax[y]-transcriptionalprofilemin[y]))
    return outdistance

def GetHammingDistance(_MatrixA, _MatrixB):
    '''Compute Hamming distances for matrices without flattening (for adjacency matrices)'''
    HammingDistance = 0
    if _MatrixA.shape == _MatrixB.shape:
        for i in range(0, _MatrixA.shape[0]):
            for j in range(0, _MatrixA.shape[1]):
                for z in range(0, _MatrixA.shape[2]):
                    if _MatrixA[i, j, z] != _MatrixB[i, j, z]:
                        HammingDistance = HammingDistance + 1
                    else:
                        pass
    else:
        raise Exception('Shapes don\'t match!')
    return HammingDistance


def GetHammingDistance_LG(_MA, _MB):
    '''Compute Hamming distances for matrices with flattening (for logic gate matrices)'''
    HammingDis = 0
    for i in range(0, len(_MA)):
        for j in range(0, len(_MA[i])):
            if _MA[i][j] == _MB[i][j]:
                pass
            else:
                HammingDis = HammingDis + 1
    return HammingDis

def WTTP_weights(WTTP, round_decimals=None, inplace=True, as_list=False):
    """
    Convert WTTP profiles to normalized weights highlighting non-background states.

    WTTP: dict like {'0': [[2], array([...])], '1': [[3], array([...])], ...}
    round_decimals: None or int. If int, round values before counting frequencies.
    inplace: if True, overwrite WTTP and return it; if False, return a new dict copy.
    as_list: if True, store weights as Python lists; otherwise as numpy arrays.

    Behavior:
      - Determine background value(s) per column as the most frequent value(s) across profiles.
      - For each profile: build a vector where non-background -> length_of_vector, background -> 1.
      - Normalize the vector by dividing by its sum (row-wise).
      - Put the normalized vector back into WTTP[key][1], keeping WTTP[key][0] as-is.
    """
    keys = list(WTTP.keys())
    if not keys:
        return {} if not inplace else WTTP

    # Extract arrays and ensure consistent shape
    arrays = [np.ravel(np.asarray(WTTP[k][1])) for k in keys]
    lengths = [a.size for a in arrays]
    if len(set(lengths)) != 1:
        raise ValueError(f"All profile arrays must have same length. Found lengths: {set(lengths)}")
    n_pos = lengths[0]
    n_profiles = len(arrays)

    mat = np.vstack(arrays)    # shape (n_profiles, n_pos)

    # Optionally round to group near-equal floats
    if round_decimals is not None:
        mat_count = np.round(mat, round_decimals)
    else:
        mat_count = mat

    # Find background value(s) per column (most frequent; ties all treated as background)
    background_sets = []
    for j in range(n_pos):
        col = mat_count[:, j]
        vals, counts = np.unique(col, return_counts=True)
        maxc = counts.max()
        bg_vals = set(vals[counts == maxc].tolist())
        background_sets.append(bg_vals)

    # Build binary matrix: 0 if value in background_sets[j], else 1
    bin_mat = np.zeros_like(mat_count, dtype=int)
    for j, bg in enumerate(background_sets):
        col = mat_count[:, j]
        # membership check per element
        bin_col = np.array([0 if v in bg else 1 for v in col], dtype=int)
        bin_mat[:, j] = bin_col

    # Map 1 -> n_pos, 0 -> 1, then normalize each row
    transformed = np.where(bin_mat == 1, n_pos, 1).astype(float)  # shape (n_profiles, n_pos)
    row_sums = transformed.sum(axis=1, keepdims=True)             # shape (n_profiles, 1)
    # Normalize (safe because row_sums > 0 since entries are at least 1)
    weights = transformed / row_sums

    # Prepare output dict
    out = WTTP if inplace else {k: [list(v[0]), np.ravel(np.asarray(v[1])).copy()] for k, v in WTTP.items()}

    # Put weights back, preserving the first element
    for i, k in enumerate(keys):
        w = weights[i]
        out[k][1] = w.tolist() if as_list else w

    return out
