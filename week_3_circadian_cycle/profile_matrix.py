from support import read_input_prob_k_mer, read_input_greed_motif_search
import numpy as np

ALPHABET = ['A','C','G','T']
letter_to_index = {'A':0, 'C':1, 'G':2, 'T':3}



def get_kmer_prob(k_mer, profile):
    res = 1
    for char, probs in zip(k_mer, profile):
        res *= probs[letter_to_index[char]]
    return res


def most_probable_kmer_from_profile(k, dna_string, profile):
    """
    Input: A string Text, an integer k, and a 4 × k matrix Profile.
    :param k: lenght of kmer
    :param dna_string:
    :param profile: transpose profile matrix, so i-th row is probabilities for i-th char in k-mer
    :return: A Profile-most probable k-mer in Text.
    """
    highest_prob = -1
    most_prob_kmer = ""

    for i in range(len(dna_string) - k + 1):
        k_mer = dna_string[i:i+k]
        prob = get_kmer_prob(k_mer, profile)
        if prob > highest_prob:
            highest_prob = prob
            most_prob_kmer = k_mer

    return most_prob_kmer


def scoreMotifMatrix(motif_list):
    """
    Scores list of motifs by creating profile matrix from candidate motifs.
    Score is sum of hamming distances from motifs to most probable string.
    :param motif_list:
    :return:
    """
    score = 0
    profile = formProfileFromMotifsPseudo(motif_list)
    for mot in motif_list:
        for i, char in enumerate(mot):
            if letter_to_index[char] != np.argmax(profile[i]):
                score += 1
    return score

def formProfileFromMotifsDummy(motifs):
    """
    creates profile matrix from a list of motifs in dummy way. We start with profile matrix full of 0.
    :param motifs:
    :return:
    """
    profile = np.ones((len(motifs[0]), 4))
    for motif in motifs:
        for i, char in enumerate(motif):
            profile[i, letter_to_index[char]] += 1
    return np.true_divide(profile, len(motifs)*2)


def formProfileFromMotifsPseudo(motifs):
    """
    Uses pseudocount, so there are no characters with 0 probability
    creates profile matrix from a list of motifs
    :param motifs:
    :return:
    """
    profile = np.ones((len(motifs[0]), 4))
    for motif in motifs:
        for i, char in enumerate(motif):
            profile[i, letter_to_index[char]] += 1
    return np.true_divide(profile, len(motifs)*2)


def greedyMotifSearch(dna_list, k, t=0):
    """
    Input: Integers k and t, followed by a collection of strings Dna.
    Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t).
    :param dna_list: list of dna sequences
    :param k: length of k-mer
    :param t: length of single dna-sequence
    :return:
    """
    # initial motifs as first k-mers from each
    best_motifs = [line[:k] for line in dna_list]
    for i in range(len(dna_list[0]) - k + 1):
        k_mer = dna_list[0][i:i+k]
        motifs_candidate = [k_mer]
        for dna_string in dna_list[1:]:
            profile = formProfileFromMotifsPseudo(motifs_candidate)
            motifs_candidate.append(most_probable_kmer_from_profile(k, dna_string, profile))
        if scoreMotifMatrix(best_motifs) > scoreMotifMatrix(motifs_candidate):
            best_motifs = motifs_candidate
    return best_motifs


def run1():
    path = "/Users/marcisin/Downloads/dataset_159_3 (1).txt"

    res = most_probable_kmer_from_profile(*read_input_prob_k_mer(path))
    print("".join(res))


def run2():
    k = 3
    t = 5
    li = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
    res = greedyMotifSearch(li, k,t)
    path = "/Users/marcisin/Downloads/dataset_160_9.txt"
    res = greedyMotifSearch(*read_input_greed_motif_search(path))
    print("\n".join(res))


def run3():
    prof = ""

def run_ttt():
    dna = ["TGACGTTC", "TAAGAGTT","GGACGAAA","CTGTTCGC"]
    m1 = ["TGA","GTT","GAA","TGT"]


    # Motifs(Profile(Motifs), Dna)
    p1 = formProfileFromMotifsPseudo(m1)
    m2 = []
    for dna_string in dna:
        profile = formProfileFromMotifsDummy(m1)
        m2.append(most_probable_kmer_from_profile(3, dna_string, profile))
    print(*m2,sep=" ")

if __name__ == '__main__':
    # run1()
    # run2()
    run_ttt()