import matplotlib.pyplot as plt

# Number of MPI Procs
x = [1, 4, 8, 16, 32]

# OMP_NUM_THREADS = 8
y8 = [0.116129, 0.135512, 0.200354, 0.181900, 0.264603]

y82 = []
for elem in y8:
    y82.append(y8[0] / elem)

print(y82)

# OMP_NUM_THREADS = 4
y4 = [0.215618, 0.234197, 0.258074, 0.284630, 0.370073]

y42 = []
for elem in y4:
    y42.append(y4[0] / elem)

print(y42)

# OMP_NUM_THREADS = 2
y2 = [0.424908, 0.444028, 0.460460, 0.496379, 0.587605]

y22 = []
for elem in y2:
    y22.append(y2[0] / elem)

print(y22)

# OMP_NUM_THREADS = 1
y1 = [0.843153, 0.863554, 0.881471, 0.920642, 1.018942]

y12 = []
for elem in y1:
    y12.append(y1[0] / elem)

print(y12)


# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(x, y8, marker='o', label='OMP_NUM_THREADS = 8')
plt.plot(x, y4, marker='o', label='OMP_NUM_THREADS = 4')
plt.plot(x, y2, marker='o', label='OMP_NUM_THREADS = 2')
plt.plot(x, y1, marker='o', label='OMP_NUM_THREADS = 1')

# plt.yscale('log')

# Add labels and title
plt.xlabel('Number of MPI Processes')
plt.ylabel('Execution Time (MPI Time)')
plt.title('Execution Time vs Number of MPI Processes for Different OMP_NUM_THREADS')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Show the plot
# plt.show()
plt.savefig('weak_scaling.png')
