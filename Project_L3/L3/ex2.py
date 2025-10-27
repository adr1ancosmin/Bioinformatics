import math
import matplotlib.pyplot as plt

s = ""
with open("sequence.fasta") as f:
    for line in f:
        if not line.startswith(">"):
            s += line.strip()

window = 9
na = 0.05

def tm1(seq):
    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")
    return 4*(g+c) + 2*(a+t)

def tm2(seq):
    gc = seq.count("G") + seq.count("C")
    gc_percent = 100 * gc / len(seq)
    return 81.5 + 16.6*math.log10(na) + 0.41*gc_percent - 600/len(seq)

x = []
y1 = []
y2 = []

for i in range(len(s)-window+1):
    w = s[i:i+window]
    x.append(i+1+(window//2))
    y1.append(tm1(w))
    y2.append(tm2(w))

plt.plot(x, y1, label="tm1")
plt.plot(x, y2, label="tm2")
plt.xlabel("Position")
plt.ylabel("Tm (°C)")
plt.title("Sliding Window Tm (window=9)")
plt.legend()
plt.show()