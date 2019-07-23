from support import read_input_gibs_samp_file
from week_4_circadian_part_2.advanced_motif_search import *
from random import choices


def prob_of_k_mer(kmer, profile):
    res = 1
    for char, prob in zip(kmer, profile):
        res *= prob[letter_to_index[char]]
    # print(res)
    return res


def profile_randomly_generated_kmer(profile, sequence, k):
    """
    Firstly we compute probability of each k-mer in the sequence based on the profile matrix
    Than based of we "randomly" choose k mer with probability computed in previous step. So most likely
    we are to choose the most likely kmer based on the profile matrix. Downside of this approach is that
    we can get stuck in local minimum based on random assignemnt of original motifs which we selected
    and consequently compute profile matrix from them and later to pick k-mer which are close to original selection
    thus we are most likely to choose k-mer which is close to original motifs -> not exploring enough.
    Very sensitive to initial random choice.
    :param profile:
    :param sequence:
    :param k:
    :return:
    """
    population = []
    weights = []
    for i in range(len(sequence) - k + 1):
        population.append(i)
        weights.append(prob_of_k_mer(sequence[i:i+k], profile))
    r = choices(population, weights)[0]
    return sequence[r:r+k]


def Gibs_sampler_runner(dna,k,t,n,runs):
    """
    Repeated runs of gibs-sampling search. We return the best motifs from all the runs
    :param dna: list of genome sequences
    :param k: k-mer length
    :param t: length of single sequence
    :param n: number of sampling cycles for gibs sampler
    :param runs: number of gibs samplers runs
    :return:
    """
    best_motifs, best_score = GibbsSampler(dna,k,t,n)
    print(best_score)
    for i in range(runs - 1):
        run_motifs, run_score = GibbsSampler(dna,k,t,n)
        if run_score < best_score:
            best_score = run_score
            best_motifs = run_motifs
            print(best_score)
    return best_motifs


def GibbsSampler(Dna, k, t, N):
    """
    randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
        BestMotifs ← Motifs
        for j ← 1 to N
            i ← Random(t)
            Profile ← profile matrix constructed from all strings in Motifs except for Motifi
            Motifi ← Profile-randomly generated k-mer in the i-th sequence
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs
    :param Dna: list of genome sequences
    :param k: length of implanted(searched) k-mer
    :param t: length of single sequence
    :param N: number of repeated sampling cycles
    :return: motifs and score with smallest hamming distance from the sequences
    """
    best_motifs = [random_k_mer(sequence, k) for sequence in Dna]
    best_score = scoreMotifMatrix(best_motifs)
    tmp_motifs = best_motifs
    for cycle_index in range(N):
        i = random.randint(0, t-1)
        profile = formProfileFromMotifsPseudo([mot for j, mot in enumerate(tmp_motifs) if j != i])
        tmp_motifs[i] = profile_randomly_generated_kmer(profile,Dna[i],k)
        tmp_score = scoreMotifMatrix(tmp_motifs)
        if tmp_score < best_score:
            best_motifs = tmp_motifs
            best_score = tmp_score
    return best_motifs, best_score


def run_1():
    k = 8
    t = 5
    N =100
    dna_list = ["CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA",
                "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
                "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
                "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
                "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
    print("\n".join(Gibs_sampler_runner(dna_list,k,t,N,20)))


def run_eval():
    path = "/Users/marcisin/Downloads/dataset_163_4 (2).txt"

    res = Gibs_sampler_runner(*read_input_gibs_samp_file(path) ,20)
    print("\n".join(res))

if __name__ == '__main__':
    # run_1()
    run_eval()