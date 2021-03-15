import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from Bio import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC


def gc_content(sequence):
    return sum([c in "GC" for c in sequence])


def gc_content2(sequence):
    return GC(sequence)


my_seq = Seq.Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
print(gc_content(my_seq))
print(gc_content2(my_seq))


def sequencer(record, window_size, step, cur=0):
    while cur + window_size < len(record.seq):
        y = 100 * gc_content(record.seq[cur:cur + window_size]) / window_size
        x = cur
        cur += step
        yield x, y


record = SeqIO.read("fasta_test.fa", "fasta")
window_size = 4000
step = 1000

xs = []
ys = []
print(f"Expected number of iterations ~{len(record.seq) // step}")
for x, y in tqdm(sequencer(record, window_size, step)):
    xs.append(x)
    ys.append(y)

plt.plot(xs, ys)
plt.xlim(0, xs[-1])
plt.ylim(0, 100)
plt.show()
