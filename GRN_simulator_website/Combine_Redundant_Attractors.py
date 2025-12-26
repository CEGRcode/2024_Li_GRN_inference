import numpy as np
from distance_functions import *

def Combine_Redundant_Attractors(arrays, cutoff=0.05):
    """
    Combine arrays that are close to each other.

    Parameters:
      arrays: a 2D NumPy array of shape (n, m), where each row is an array.
      cutoff: a float; two arrays are "close" if, for every column i,
              abs(a[i] - b[i]) / max_i < cutoff,
              where max_i is the maximum value among all arrays in column i.

    Returns:
      combined: a 2D NumPy array where each row is the average of a group of close arrays.
      groups: a list of groups (each group is a list of indices from the input arrays).
    """
    # Compute the maximum per column across all arrays

    print('arrays', arrays)

    max_values = np.max(arrays, axis=0)

    n_arrays = len(arrays)
    used = np.full(n_arrays, False)
    groups = []
    print('Combining redundant attractor profiles...')
    # Loop over each array to group similar ones together.
    for i in range(n_arrays):
        if used[i]:
            continue  # Skip arrays already grouped
        # Start a new group with the i-th array as reference.
        group = [i]
        used[i] = True

        for j in range(i + 1, n_arrays):
            if used[j]:
                continue
            # Calculate the normalized absolute difference element-wise
            diff = np.abs(arrays[i] - arrays[j]) / max_values
            # Check if all differences are below the cutoff
            if np.all(diff < cutoff):
                group.append(j)
                used[j] = True

        groups.append(group)

    # Combine the arrays in each group by taking the mean
    combined = np.array([np.mean(arrays[group], axis=0) for group in groups])

    return combined

