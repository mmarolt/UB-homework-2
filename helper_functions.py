def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("the brown cat", "these brownies", lambda x, y: [-1, 1][x == y])
    ('----the brown cat', 'thes--e brownies-', 5.0)
    Other alignments are also possible.
    """
    
    
    M = [[0] * (len(seq2)+1) for i in range(len(seq1)+1)]
    
    for i in range(len(seq1)+1):
        for j in range(len(seq2)+1):
            if i == j == 0:
                M[i][j] = 0.0
            elif i == 0:
                M[i][j] = M[i][j-1] + scoring_function('-', seq2[j-1])
            elif j == 0:
                M[i][j] = M[i-1][j] + scoring_function(seq1[i-1], '-')
            else:
                Mi_1j = M[i-1][j] + scoring_function(seq1[i-1], '-')
                Mij_1 = M[i][j-1] + scoring_function('-', seq2[j-1])
                Mi_1j_1 = M[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
                M[i][j] = max(Mi_1j, Mij_1, Mi_1j_1)
                
    
    i, j = len(seq1), len(seq2)
    aligned1, aligned2 = "", ""
    while i != 0 or j != 0:
        if i == 0:
            aligned1 = ('-' * j) + aligned1
            aligned2 = seq2[:j] + aligned2
            j = 0
        elif j == 0:
            aligned1 = seq1[:i] + aligned1
            aligned2 = ('-' * i) + aligned2
            i = 0
        else:
            Mi_1j = M[i-1][j] + scoring_function(seq1[i-1], '-')
            Mij_1 = M[i][j-1] + scoring_function('-', seq2[j-1])
            Mi_1j_1 = M[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            if M[i][j] == Mi_1j:
                aligned1 = seq1[i-1] + aligned1
                aligned2 = '-' + aligned2
                i -= 1
            elif M[i][j] == Mij_1:
                aligned1 = '-' + aligned1
                aligned2 = seq2[j-1] + aligned2
                j -= 1
            else:
                aligned1 = seq1[i-1] + aligned1
                aligned2 = seq2[j-1] + aligned2
                i -= 1
                j -= 1
            
    return aligned1, aligned2, M[len(seq1)][len(seq2)]


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("the brown cat", "these brownies", lambda x, y: [-1, 1][x == y])
    ('the-- brown', 'these brown', 7.0)

    Other alignments are also possible.

    """
    M = [[0] * (len(seq2)+1) for i in range(len(seq1)+1)]
    
    max_score = 0, 0, 0
    for i in range(len(seq1)+1):
        for j in range(len(seq2)+1):
            if i == j == 0:
                M[i][j] = 0.0
            elif i == 0:
                M[i][j] = max(0, M[i][j-1] + scoring_function('-', seq2[j-1]))
            elif j == 0:
                M[i][j] = max(0, M[i-1][j] + scoring_function(seq1[i-1], '-'))
            else:
                Mi_1j = M[i-1][j] + scoring_function(seq1[i-1], '-')
                Mij_1 = M[i][j-1] + scoring_function('-', seq2[j-1])
                Mi_1j_1 = M[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
                M[i][j] = max(0, Mi_1j, Mij_1, Mi_1j_1)
            if M[i][j] > max_score[0]:
                max_score = M[i][j], i, j
                
    
    i, j = len(seq1), len(seq2)
    aligned1, aligned2 = "", ""
    score, i, j = max_score
    while M[i][j] > 0:
        if i == 0:
            aligned1 = '-' + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1
        elif j == 0:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = '-' + aligned2
            i -= 1
        else:
            Mi_1j = M[i-1][j] + scoring_function(seq1[i-1], '-')
            Mij_1 = M[i][j-1] + scoring_function('-', seq2[j-1])
            Mi_1j_1 = M[i-1][j-1] + scoring_function(seq1[i-1], seq2[j-1])
            if M[i][j] == Mi_1j:
                aligned1 = seq1[i-1] + aligned1
                aligned2 = '-' + aligned2
                i -= 1
            elif M[i][j] == Mij_1:
                aligned1 = '-' + aligned1
                aligned2 = seq2[j-1] + aligned2
                j -= 1
            else:
                aligned1 = seq1[i-1] + aligned1
                aligned2 = seq2[j-1] + aligned2
                i -= 1
                j -= 1
            
    return aligned1, aligned2, score


def translate_to_protein(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }
    protein = ""
    for i in range(0, len(seq) - len(seq)%3, 3):
        protein += table[seq[i:i+3]]
    return protein