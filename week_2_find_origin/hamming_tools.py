import itertools

def hamming_circle(s, n, alphabet):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.

    >>> sorted(hamming_circle('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_circle('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_circle('aaa', 2, 'ab'))
    ['abb', 'bab', 'bba']

    """
    for positions in itertools.combinations(range(len(s)), n):
        for replacements in itertools.product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin)


def hamming_ball(s, n, alphabet=['A', 'C', 'G', 'T']):
    """Generate strings over alphabet whose Hamming distance from s is
    less than or equal to n.

    >>> sorted(hamming_ball('abc', 0, 'abc'))
    ['abc']
    >>> sorted(hamming_ball('abc', 1, 'abc'))
    ['aac', 'aba', 'abb', 'abc', 'acc', 'bbc', 'cbc']
    >>> sorted(hamming_ball('aaa', 2, 'ab'))
    ['aaa', 'aab', 'aba', 'abb', 'baa', 'bab', 'bba']

    """
    return itertools.chain.from_iterable(hamming_circle(s, i, alphabet)
                                         for i in range(n + 1))


def hamming_dist(gene_1, gene_2):
    """
    Computes hamming distance between 2 strings
    :param gene_1:
    :param gene_2:
    :return:
    """
    ham_dist = 0
    for c1, c2 in zip(gene_1, gene_2):
        if c1 != c2:
            ham_dist += 1
    return ham_dist
