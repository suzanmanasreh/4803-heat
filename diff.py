import numpy as np
from numpy.linalg import norm

seq_file = open("2d/ex_0100.vtk", "r")
seq_nums = []
for line in seq_file.readlines():
    if line[0].isnumeric():
        seq_nums.append(float(line.strip()))

seq = np.array(seq_nums)

par_file = open("2dp/ex_0100.vtk", "r")
par_nums = []
for line in par_file.readlines():
    if line[0].isnumeric():
        par_nums.append(float(line.strip()))

par = np.array(par_nums)

diff = par - seq


count = 0
for num in diff:
    if num != 0:
        count += 1

# print(diff)

# print(norm(diff))

l2_norm = round(norm(diff), 5)

# i am not sure what a suitable upper bound on the L2 norm is

total = np.size(diff)
percent_diff = round((count / total) * 100, 2)

avg_diff = round(l2_norm / total, 7)

print(f"{percent_diff}% of values are different with a total difference of {l2_norm} and an average difference of {avg_diff}")
