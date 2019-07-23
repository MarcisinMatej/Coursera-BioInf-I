def find_clump(sequence, k, window_size, treshold):
    """

    :param sequence:
    :param window_size:
    :param k: k-mer length
    :param treshold: min number of matches of k-mer in window
    :return:
    """
    clump_candidates = set()
    # sliding of the window
    window = sequence[0:window_size]
    k_mers = get_k_mers_for_window(window, k)
    for i in range(len(sequence) + 1 - window_size):
        window = sequence[i:i + window_size]
        # decrease first k-mer count
        key = window[:k]
        k_mers[key] = k_mers[key] - 1
        # add new k-mer
        key = window[-k:]
        if key in k_mers.keys():
            k_mers[key] = k_mers[key] + 1
        else:
            k_mers[key] = 1
        clump_candidates.update(get_k_mers_over_threshold(k_mers, treshold))
    return clump_candidates


def get_k_mers_for_window(window, k):
    res = {}
    for i in range(len(window) - k + 1):
        key = window[i:i+k]
        if key in res.keys():
            res[key] = res[key] + 1
        else:
            res[key] = 1
    return res


def get_k_mers_over_threshold(k_mers, l):
    candidates = []
    for key in k_mers.keys():
        if k_mers[key] >= l:
            candidates.append(key)
    return candidates

