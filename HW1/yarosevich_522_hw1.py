import numpy as np

n_0 = np.array([[100, 100, 100, 100]])

fecun = np.array([0, 1, 5, .5])
p_i = np.diag([.5, .9, .95, 0])

A = np.row_stack([fecun,p_i])
A = np.delete(A, 4, 0)

t_max = 50

n_t = np.zeros((4, t_max))
n_t[:, 0] = n_0

for t in range(1, t_max - 1):
    n_t[:, t] = np.dot(A, n_t[:, t-1])

