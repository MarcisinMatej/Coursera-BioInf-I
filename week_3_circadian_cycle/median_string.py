import math

from support import read_input_median_file, read_input_median_file_distance
from week_2_find_origin.hamming_tools import hamming_dist, hamming_ball


def sum_distance_from_dna_list(pattern, dna_list):
    """
    sum of min_distances
    :param pattern:
    :param dna_list:
    :return:
    """
    res = 0
    for dna_sample in dna_list:
        res+=min_distance(pattern, dna_sample)
    return res

def min_distance(pattern, dna_string):
    """
     minimum Hamming distance between Pattern and any k-mer in Text
    :param pattern:
    :param dna_string:
    :return:
    """
    min_dist = math.inf
    k = len(pattern)
    for i in range(len(dna_string) - k + 1):
        sub_pattern = dna_string[i:i+k]
        if hamming_dist(sub_pattern,pattern) < min_dist:
            min_dist = hamming_dist(sub_pattern,pattern)
    return min_dist


def all_patterns(k):
    """
    creates all k-mer
    :param k:
    :return:
    """
    return hamming_ball('A'*k, k)


def MedianString(k, dna_list):
    """
    Input: An integer k, followed by a collection of strings Dna.
     Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers. (If there are
    multiple such strings Pattern, then you may return any one.)
    :param k: length of motif
    :param dna_list: list of dna strings
    :return: median string
    """

    min_dist = math.inf
    median_str = ""
    for pattern in all_patterns(k):
        tmp_dist = sum_distance_from_dna_list(pattern, dna_list)
        if tmp_dist <= min_dist:
            print(min_dist, pattern)
            min_dist = tmp_dist
            median_str = pattern
    return median_str


def run1():
    path = "/Users/marcisin/Downloads/dataset_158_9.txt"
    res = MedianString(*read_input_median_file(path))
    print(res)

def run2():
    path = "/Users/marcisin/Downloads/dataset_5164_1 (2).txt"
    res = sum_distance_from_dna_list(*read_input_median_file_distance(path))
    print(res)

def run3():
    li = ["CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC", "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC","GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"]
    res = MedianString(7,li)
    print(res)

if __name__ == '__main__':
    # run1()
    # run2()
    run3()