from support import read_input_motif_file
from week_2_find_origin.approximate_matching import approximate_pattern_match
from week_2_find_origin.hamming_tools import hamming_ball


def MotifEnumeration(Dna, k, d):
    """
    Input: Integers k and d, followed by a collection of strings Dna.
     Output: All (k, d)-motifs in Dna.
    :param Dna: list of strings
    :param k: length of motif
    :param d: max distance
    :return:
    """
    motifs = {}
    for dna_part in Dna:
        for i in range(len(dna_part) -k + 1):
            pattern = dna_part[i:i+k]
            neighbourhood = hamming_ball(pattern,d)
            for neighbour in neighbourhood:
                add_flg = True
                for gene in Dna:
                    if approximate_pattern_match(gene,neighbour,d) == 0:
                        add_flg = False
                        break
                if add_flg:
                    if neighbour in motifs.keys():
                        motifs[neighbour] = motifs[neighbour] + 1
                    else:
                        motifs[neighbour] = 1
    return motifs.keys()


def run_motif1():
    tst_in = ["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"]
    res = MotifEnumeration(tst_in, 3, 1)
    print(" ".join(res))
    path = "/Users/marcisin/Downloads/dataset_156_8.txt"

    res = MotifEnumeration(*read_input_motif_file(path))
    print(" ".join(res))


if __name__ == '__main__':
    run_motif1()
