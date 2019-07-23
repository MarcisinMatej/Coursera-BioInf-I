

def read_input_file(path):
    with open(path) as f:
        return [line.rstrip('\n') for line in f]

def read_input_file_kokotno(path):
    with open(path) as f:
        return [line for line in f]


def read_input_motif_file(path):
    lines = read_input_file(path)
    k,d = map(int, lines[0].split(' '))
    return lines[1:] ,k ,d


def read_input_median_file(path):
    lines = read_input_file(path)
    k = int(lines[0])
    return k, lines[1:]


def read_input_median_file_distance(path):
    """

    :param path:
    :return: pattern, dna_list
    """
    lines = read_input_file(path)
    return lines[0], lines[1].split()

def read_matrix(sub_lines):
    tmp = []
    for l in sub_lines:
        tmp.append(list(map(float, l.split(' '))))
    return [[tmp[j][i] for j in range(len(tmp))] for i in range(len(tmp[0]))]


def read_input_prob_k_mer(path):
    """
    file in following form string Text, an integer k, and a 4 × k matrix Profile.
    ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
    5
    0.2 0.2 0.3 0.2 0.3
    0.4 0.3 0.1 0.5 0.1
    0.3 0.3 0.5 0.2 0.4
    0.1 0.2 0.1 0.1 0.2
    """
    lines = read_input_file(path)

    matrix = read_matrix(lines[2:])
    return int(lines[1]), lines[0], matrix


def read_input_greed_motif_search(path):
    lines = read_input_file(path)
    k, t = map(int, lines[0].split(' '))
    return lines[1:], k, t


def read_input_gibs_samp_file(path):
    lines = read_input_file(path)
    k, t, n = map(int, lines[0].split(' '))
    return lines[1:], k, t, n