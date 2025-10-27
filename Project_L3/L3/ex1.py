import math

s = "ATGCGACTGCGT"

a = s.count("A")
t = s.count("T")
g = s.count("G")
c = s.count("C")

tm = 4*(g+c) + 2*(a+t)

na = 0.05
gc_percent = 100*(g+c)/len(s)

tm2 = 81.5 + 16.6*math.log10(na) + 0.41*gc_percent - 600/len(s)

print("Tm1 =", tm)
print("Tm2 =", tm2)
