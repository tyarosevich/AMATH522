import numpy as np
import matplotlib.pyplot as plt
import random as rd

numsteps = 2000
dwell_list = []

A = np.array([[.98, .1, 0], [.02, .7, .05], [0, .2, .95]])
c1_sort = np.argsort(A[:,0])
c2_sort = np.argsort(A[:,1])
c3_sort = np.argsort(A[:,2])
sort_list = [c1_sort, c2_sort, c3_sort]

states_start = np.zeros(numsteps, dtype = int)
states = states_start

states[0] = 0

for i in range(0,numsteps):
    count = 0
    for k in np.arange(0, 500):
        count = count + 1
        rand = rd.random()
        if rand < A[sort_list[states[k]][0], states[k]]:
            states[k+1] = sort_list[states[k]][0]
        elif rand < A[sort_list[states[k]][1], states[k]]:
            states[k+1] = sort_list[states[k]][1]
        else:
            states[k+1] = sort_list[states[k]][2]
        if states[k+1] != states[k]:
            break
    dwell_list.append(count)

reduced_states = states
for i in range(0, len(reduced_states)):
    if reduced_states[i] == 1:
        reduced_states[i] = 0
for i in range(0, len(reduced_states)):
    if reduced_states[i] == 2:
        reduced_states[i] = 1


# n_bins = 20
# n, bins, patches = plt.hist(dwell_list, n_bins, facecolor='blue', alpha=0.5)
# plt.show()