import numpy as np
from Bio import SeqIO

def Needleman_Wunsh_matrix(SEQ_list,GP=-2,DIFF=-1,SAME=2,MAX_SEQ_LENGTH=20):
    MatrixOfDistance = [[0 for x in range(len(SEQ_list))] for y in range(len(SEQ_list))]
    for i in range(len(SEQ_list)-1):
        for j in range(i+1,len(SEQ_list)):
            MatrixOfDistance[i][j]=Needleman_Wunsh(SEQ_list[i],SEQ_list[j],GP,DIFF,SAME,MAX_SEQ_LENGTH)
            MatrixOfDistance[j][i]=MatrixOfDistance[i][j]
    return MatrixOfDistance


def Needleman_Wunsh(SEQ1="", SEQ2="", GP=-2, DIFF=-5, SAME=5, MAX_SEQ_LENGTH=-1):
    def calculate_values(MatrixOfValues, SEQ1, SEQ2, GP, DIFF, SAME):
        for x in range(len(SEQ2)):
            for y in range(len(SEQ1)):
                left = MatrixOfValues[x + 1][y] + GP
                up = MatrixOfValues[x][y + 1] + GP
                if (SEQ2[x] == SEQ1[y]):
                    slant = MatrixOfValues[x][y] + SAME
                else:
                    slant = MatrixOfValues[x][y] + DIFF

                newValue = max(left, up, slant)
                MatrixOfValues[x + 1][y + 1] = newValue

    if (MAX_SEQ_LENGTH < 1):
        print(MAX_SEQ_LENGTH)
        print("Set MAX_SEQ_LENGTH to positive integer")
        return
    elif (MAX_SEQ_LENGTH < max(len(SEQ1), len(SEQ2))):
        print("MAX_SEQ_LENGTH reached")
        return

    # setting up needed matrices
    MatrixOfValues = [[0 for x in range(len(SEQ1) + 1)] for y in range(len(SEQ2) + 1)]
    for i in range(len(SEQ1) + 1):
        MatrixOfValues[0][i] = i * GP

    for i in range(len(SEQ2) + 1):
        MatrixOfValues[i][0] = i * GP

    # calculating values and updating directions
    calculate_values(MatrixOfValues, SEQ1, SEQ2, GP, DIFF, SAME)
    # creating solution based on MatrixOfDirections
    return MatrixOfValues[len(SEQ2)][len(SEQ1)]

def find_min_score_cell(table):
    search_matrix = np.diag(np.repeat(2*np.max(table), table.shape[0]))+table
    x, y = np.unravel_index(search_matrix.argmin(), search_matrix.shape)
    return x, y


def join_and_weigh_labels(weighted_labels, a, b):
    if b < a:
        a, b = b, a

    weighted_labels[a] = (weighted_labels[a][0], weighted_labels[b][0]), weighted_labels[a][1]+weighted_labels[b][1]

    weighted_labels.pop(b)

    return weighted_labels

def join_and_rescore_table(table, weighted_labels, a, b):
    if b < a:
        a, b = b, a

    new_table = table.copy()

    for i in range(0, a):
        new_table[a][i] = (table[a][i]*weighted_labels[a][1] + table[b][i]*weighted_labels[b][1]) / (weighted_labels[a][1]+weighted_labels[b][1])
        new_table[i][a] = (table[a][i]*weighted_labels[a][1] + table[b][i]*weighted_labels[b][1]) / (weighted_labels[a][1]+weighted_labels[b][1])

    for i in range(a + 1, b):
        new_table[a][i] = (table[a][i] * weighted_labels[a][1] + table[b][i] * weighted_labels[b][1]) / (weighted_labels[a][1] + weighted_labels[b][1])
        new_table[i][a] = (table[a][i] * weighted_labels[a][1] + table[b][i] * weighted_labels[b][1]) / (weighted_labels[a][1] + weighted_labels[b][1])


    for i in range(b + 1, len(weighted_labels)):
        new_table[a][i] = (table[a][i] * weighted_labels[a][1] + table[b][i] * weighted_labels[b][1]) / (weighted_labels[a][1] + weighted_labels[b][1])
        new_table[i][a] = (table[a][i] * weighted_labels[a][1] + table[b][i] * weighted_labels[b][1]) / (weighted_labels[a][1] + weighted_labels[b][1])

    new_table = np.c_[new_table[:, :b], new_table[:, (b+1):]]
    new_table = np.r_[new_table[:b, :], new_table[(b+1):, :]]
    return new_table

def UPGMA(matrix, labels):

    weighted_labels = list(zip(labels, [1]*len(labels)))
    current_table = matrix.copy()
    while current_table.shape[0] > 1:
        x, y = find_min_score_cell(current_table)

        current_table = join_and_rescore_table(current_table, weighted_labels, x, y)

        weighted_labels = join_and_weigh_labels(weighted_labels, x, y)

    return weighted_labels[0][0]


def fitch_margoliash(dist):
    cluster = []
    tree = {}

    n = dist.shape[0]
    n_ = n
    names = [i for i in range(n)]
    while len(cluster) < n - 1:
        i = dist.reshape(-1, ).argmin()

        row = i // n_
        col = i % n_

        row_name = names[row]
        col_name = names[col]

        if isinstance(row_name, int):
            cluster.append(row_name)

        if isinstance(col_name, int):
            cluster.append(col_name)

        ind = np.delete(np.arange(n_, dtype=int), [row, col])
        dist_row = dist[row, ind].mean()
        dist_col = dist[ind, col].mean()

        new_name = (row_name, col_name)

        tree[(row_name, new_name)] = (dist_row + dist[row, col] - dist_col) / 2
        tree[(col_name, new_name)] = (dist_col + dist[row, col] - dist_row) / 2

        if len(cluster) == n - 1:
            tree[(names[0], new_name)] = (dist_col - dist[row, col] + dist_row) / 2
            break

        names.remove(row_name)
        names.remove(col_name)
        names.append(new_name)

        dist = np.hstack((dist, np.zeros((n_, 1))))
        dist = np.vstack((dist, np.zeros((1, n_ + 1))))
        dist[-1, -1] = np.inf

        for i in range(n_ + 1):
            dist[i, -1] = (dist[i, row] + dist[i, col]) / 2
            dist[-1, i] = dist[i, -1]

        dist = np.delete(dist, [row, col], axis=1)
        dist = np.delete(dist, [row, col], axis=0)

        n_ -= 1

    return tree


class SmithWaterman:
    def __init__(self, sub=None, w=1):
        if sub is not None:
            self.sub = sub
        else:
            self.sub = lambda a, b: 1 if a == b else -1
        self.w = w

    def align(self, a, b):
        n = len(a)
        m = len(b)

        H = np.zeros((n + 1, m + 1))

        for i, aa in enumerate(a):
            for j, bb in enumerate(b):
                v1 = H[i, j] + self.sub(aa, bb)
                v2 = (H[:i + 1, j + 1] - self.w * np.arange(i + 1, 0, -1)).max()
                v3 = (H[i + 1, :j + 1] - self.w * np.arange(j + 1, 0, -1)).max()
                scores = [v1, v2, v3, 0]
                m_ = np.argmax(scores)
                H[i + 1, j + 1] = scores[m_]

        i = np.argmax(H.reshape(-1, ))
        j = i % (m + 1)
        i //= m + 1
        return H[i, j]

def SmithWaterman_matrix(SEQ_list):
    MatrixOfDistance = [[0 for x in range(len(SEQ_list))] for y in range(len(SEQ_list))]
    sw = SmithWaterman()
    for i in range(len(SEQ_list)-1):
        for j in range(i+1,len(SEQ_list)):
            MatrixOfDistance[i][j]=sw.align(SEQ_list[i], SEQ_list[j])
            MatrixOfDistance[j][i]=MatrixOfDistance[i][j]
    return MatrixOfDistance

human = SeqIO.read("P04439.fasta", "fasta").seq
zebrafish=SeqIO.read("Q2WGN8.fasta", "fasta").seq
mouse=SeqIO.read("Q7ZZ93.fasta", "fasta").seq
rat=SeqIO.read("Q9ESL4.fasta", "fasta").seq
fly=SeqIO.read("Q00900.fasta", "fasta").seq

all_seq=[human,zebrafish,mouse,rat,fly]
dist=SmithWaterman_matrix(all_seq)
dist=np.array(dist)
dist2=dist
np.fill_diagonal(dist, np.inf)
result=fitch_margoliash(dist)
print(result)
result2=UPGMA(dist2, ["H","Z","M","R","F"])
print(result2)