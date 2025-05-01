import matplotlib.pyplot as plt

# Number of MPI Procs
x = [1, 4, 8, 16, 32]

# OMP_NUM_THREADS = 8
y8 = [10.894041, 3.083596, 1.751325, 1.077019, 0.745260]

y82 = []
for elem in y8:
    y82.append(y8[0] / elem)

print(y82)

# OMP_NUM_THREADS = 4
y4 = [21.306747, 5.669962, 3.058492, 1.736548, 1.094117]

y42 = []
for elem in y4:
    y42.append(y4[0] / elem)

print(y42)

# OMP_NUM_THREADS = 2
y2 = [42.142973, 10.920386, 5.675510, 3.053143, 1.895220]

y22 = []
for elem in y2:
    y22.append(y2[0] / elem)

print(y22)

# OMP_NUM_THREADS = 1
y1 = [83.815933, 21.398151, 10.942583, 5.706526, 3.091201]

y12 = []
for elem in y1:
    y12.append(y1[0] / elem)

print(y12)


# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(x, y82, marker='o', label='OMP_NUM_THREADS = 8')
plt.plot(x, y42, marker='o', label='OMP_NUM_THREADS = 4')
plt.plot(x, y22, marker='o', label='OMP_NUM_THREADS = 2')
plt.plot(x, y12, marker='o', label='OMP_NUM_THREADS = 1')

# plt.yscale('log')

# Add labels and title
plt.xlabel('Number of MPI Processes')
plt.ylabel('Speedup')
plt.title('Speedup vs Number of MPI Processes for Different OMP_NUM_THREADS')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Show the plot
# plt.show()
plt.savefig('scaling.png')
