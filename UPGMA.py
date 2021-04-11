import numpy as np


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


matrix = np.array(
    [[0.0, 17.0, 21.0, 31.0, 23.0],
     [17.0, 0.0, 30.0, 34.0, 21.0],
     [21.0, 30.0, 0.0, 28.0, 39.0],
     [31.0, 34.0, 28.0, 0.0, 43.0],
     [23.0, 21.0, 39.0, 43.0, 0.0]])

print(UPGMA(matrix, ["A","B","C","D","E"]))