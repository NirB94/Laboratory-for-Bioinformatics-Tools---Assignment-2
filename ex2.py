# Bioinfotools exercise 2
# Put your code instead of the 'pass' statements
import random
import re


# Auxiliary function to read a FASTA file. Returns a string
def read_fasta(filename):
    f = open(filename)
    header = f.readline()
    x = f.read()
    x = x.replace("\n", "")
    return x


# Part A
# Scoring regime
def score(a, b):
    if (a == b):
        return 3
    elif ((a == "-" and b != "-") or (a != "-" and b == "-")):
        return -2
    else:
        return -3


def local_pwalignment(S, T, score=score):
    n = len(S)
    m = len(T)
    align_matrix = matrix_constructor(n, m)
    max_cell = (0, 0)
    max_val = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            align_matrix[i][j] = check_max(S, T, align_matrix, i, j)
            if align_matrix[i][j][0] > max_val:
                max_val = align_matrix[i][j][0]
                max_cell = (i, j)
    seq1 = ''
    seq2 = ''
    align_score = align_matrix[max_cell[0]][max_cell[1]][0]
    curr_cell = max_cell
    while curr_cell is not None:
        i2 = curr_cell[0]
        j2 = curr_cell[1]
        cell = align_matrix[i2][j2]
        if cell[1] is not None:
            father_i = cell[1][0]
            father_j = cell[1][1]
            if (i2 > father_i) and (j2 > father_j):  # go diagonally
                seq1 = S[i2 - 1] + seq1
                seq2 = T[j2 - 1] + seq2
            elif (i2 > father_i) and (j2 == father_j):  # go up
                seq1 = S[i2 - 1] + seq1
                seq2 = '-' + seq2
            elif (i2 == father_i) and (j2 > father_j):  # go left
                seq1 = '-' + seq1
                seq2 = T[j2 - 1] + seq2
        curr_cell = cell[1]
    return align_score, seq1, seq2


def matrix_constructor(n, m):
    return [[(0, None) for j in range(m + 1)] for i in range(n + 1)]


def check_max(seq1, seq2, mat, i, j, score=score):
    prev_diag = mat[i - 1][j - 1][0] + score(seq1[i - 1], seq2[j - 1])
    prev_row = mat[i - 1][j][0] + score(seq1[i - 1], "-")
    prev_col = mat[i][j - 1][0] + score("-", seq2[j - 1])
    maximum = max(prev_col, prev_diag, prev_row, 0)
    if maximum == prev_diag:
        return maximum, (i - 1, j - 1)
    if maximum == prev_row:
        return maximum, (i - 1, j)
    if maximum == prev_col:
        return maximum, (i, j - 1)
    if maximum == 0:
        return 0, None


# Part B
def find_strs(S, s, r):
    return len(re.findall('('+s+'){'+str(r)+',}', S))


def find_strs3(S, r):
    all_strs = generate_strs()
    total_str_num = 0
    for strs in all_strs:
        total_str_num += find_strs(S, strs, r)
    return total_str_num


def generate_strs():
    bases = ['G', 'C', 'A', 'T']
    str_list = []
    for base1st in bases:
        for base2nd in bases:
            for base3rd in bases:
                if base3rd == base2nd == base1st:
                    continue
                else:
                    str_list.append(base1st + base2nd + base3rd)
    return str_list


def permutation_test(S, r):
    perm_score_list = []
    original_score = find_strs3(S, r)
    for i in range(100):
        perm = generate_permutation(S)
        perm_score_list.append(find_strs3(perm, r))
    sum_true = sum([(original_score < perm_score) for perm_score in perm_score_list])
    return (sum_true / 100) < 0.05


def generate_permutation(S):
    return ''.join(random.sample(S, len(S)))


if __name__ == '__main__':
    ## Part a
    foxp1_human_path = 'Foxp1_Homo_sapiens.fasta'
    foxp1_mus_path = 'Foxp1_Mus_musculus.fasta'
    foxp1_bos_path = 'Foxp1_Bos_taurus.fasta'
    foxp1_human = read_fasta(foxp1_human_path)
    foxp1_mus = read_fasta(foxp1_mus_path)
    foxp1_bos = read_fasta(foxp1_bos_path)
    print(local_pwalignment("GGTTGACTA", "TGTTACGG", score=score)[0] == 13)
    print(local_pwalignment(foxp1_human, foxp1_bos))
    print(local_pwalignment(foxp1_bos, foxp1_human))


    ## Part b

    print(find_strs("AAGAGAGTTAGAGTCAGC", "AG", 2) == 2)
    print(find_strs("AAGAGAGTTAGAGTCAGC", "AG", 3) == 1)
    print(find_strs3("AAAGGAGGTGTTCGGTCGTCGTC", 2) == 4)
    print(find_strs3("AAAGGAGGTGTTCGGTCGTCGTC", 3) == 1)
    genome_path = "genome.fasta"
    pass