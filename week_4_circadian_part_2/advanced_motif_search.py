import math
import random

from support import read_input_motif_file
from week_3_circadian_cycle.profile_matrix import *


def random_k_mer(sequence, k):
    start = random.randint(0, len(sequence) - k)
    return sequence[start:start+k]


def find_best_random_motifs(Dna, k):
    best_motifs = [random_k_mer(sequence, k) for sequence in Dna]
    best_score = scoreMotifMatrix(best_motifs)
    tmp_motifs = best_motifs
    while True:
        # form profile matrix
        profile = formProfileFromMotifsPseudo(tmp_motifs)
        # find new motifs as most probable k-mer from each sequence
        tmp_motifs = [most_probable_kmer_from_profile(k, dna_seq, profile) for dna_seq in Dna]

        if scoreMotifMatrix(tmp_motifs) < best_score:
            best_motifs = tmp_motifs
            best_score = scoreMotifMatrix(best_motifs)
        else:
            break
    return best_motifs, best_score



def RandomizedMotifSearch(Dna, k, t=0, cyc_cnt = 1):
    the_best_motifs = []
    the_best_score = math.inf
    for i in range(cyc_cnt):
        motifs, score = find_best_random_motifs(Dna,k)
        if the_best_score > score:
            the_best_score = score
            the_best_motifs = motifs
            print("-------------", str(the_best_score))
    return the_best_motifs


def run1():
    k = 8
    t = 5
    dna_list = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
                "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
                "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
                "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
                "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
    print("\n".join(RandomizedMotifSearch(dna_list,k,t,1000)))


def run_tst():
    path = "/Users/marcisin/Downloads/dataset_161_5 (2).txt"

    res = RandomizedMotifSearch(*read_input_motif_file(path),1000)
    print("\n".join(res))


def run_prob_comp():
    """
    Compute the probability that ten randomly selected 15-mers from the ten 600-nucleotide long strings in
    the Subtle Motif Problem capture at least one implanted 15-mer. (Allowable error: 0.000001)
    prob(at least one) = 1 - prob(none)
    """
    p1 = 1 - (1 / (600-14))#(600-15)/(600-15+1)
    r = 1
    for i in range(10):
        r*=p1
    print(1-r)



if __name__ == '__main__':
    # run1()
    # run_tst()
    run_prob_comp()