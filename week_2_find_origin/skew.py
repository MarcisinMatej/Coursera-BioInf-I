
def skew_sequence(genome):
    skew = 0
    tmp = [0]
    for nucleotide in genome.upper():
        if nucleotide is 'G':
            skew += 1
        if nucleotide is 'C':
            skew -= 1
        tmp.append(skew)
    print(tmp)
    return tmp


def MinimumSkew(Genome):
    """
    returns list of indeces where we have minimal skew
    :param Genome:
    :return:
    """
    skew = 0
    tmp = [0]
    for nucleotide in Genome.upper():
        if nucleotide is 'G':
            skew += 1
        if nucleotide is 'C':
            skew -= 1
        tmp.append(skew)
    indeces = []
    m_v = min(tmp)
    for i,val in enumerate(tmp):
        if val == m_v:
            indeces.append(i)
    return indeces



def Skew(gene, k):
    """
    Computes difference between G concrentratin and C concentration for k first nucleaotides
    :param gene: sequence of ATGC nucleotides
    :param k:
    :return: difference G count - C count
    """
    substr = gene[:k].upper()
    return substr.count('G') - substr.count('C')


def find_k_for_min_skew(gene):
    """
    the skew should achieve a minimum at the position where the reverse half-strand ends and the forward half-strand begins
    :param gene:
    :return:
    """
    # min_k = -1
    # min_val = 9999999
    # for i in range(len(gene)):
    #     if Skew(gene,i) < min_val:
    #         min_val = Skew(gene,i)
    #         min_k = i
    skew_list = skew_sequence(gene)
    min_val = min(skew_list)
    min_k = [i for i,val in enumerate(skew_list) if val==min_val]
    print("Min value: ", str(min_val), " Min k:", str(min_k))


def find_k_for_max_skew(gene):
    max_k = -1
    max_val = -1
    for k in range(len(gene)):
        if Skew(gene,k) > max_val:
            max_val = Skew(gene,k)
            max_k = k
    print("Maxn value: ", str(max_val), " Min k:", str(max_k))
