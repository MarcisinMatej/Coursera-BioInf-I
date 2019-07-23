from week_2_find_origin.hamming_tools import hamming_dist


def approximate_pattern_match(gene, pattern, max_d):
    """

    :param gene:
    :param pattern: pattern to be searched for
    :param max_d: max hamming distance acceptable
    :return: number of approximate matches
    """
    cnt = 0
    l = len(pattern)
    for i in range(len(gene)-l + 1):
        if hamming_dist(gene[i:i + l], pattern) <= max_d:
            cnt +=1
    return cnt


def approximate_pattern_match_indexes(gene, pattern, max_d):
    """
    :param gene:
    :param pattern: pattern to be searched for
    :param max_d: max hamming distance acceptable
    :return: number of approximate matches
    """
    cnt = 0
    l = len(pattern)
    for i in range(len(gene) - l + 1):
        if hamming_dist(gene[i:i + l], pattern) <= max_d:
            print(i, end=" ")
    return cnt
