import numpy as np
import matplotlib as plt
import random as rd

numsteps = 2000

A = np.array([[.98, .1, 0], [.02, .7, .05], [0, .2, .95]])
c1_sort = np.argsort(A[:,0])
c2_sort = np.argsort(A[:,1])
c3_sort = np.argsort(A[:,2])
sort_list = [c1_sort, c2_sort, c3_sort]

states = np.zeros(numsteps, dtype = int)

states[0] = 0

for k in np.arange(0, numsteps-1):
    rand = rd.random()
    if rand < A[sort_list[states[k]][0], states[k]]:
        states[k+1] = sort_list[states[k]][0]
    elif rand < A[sort_list[states[k]][1], states[k]]:
        states[k+1] = sort_list[states[k]][1]
    else:
        states[k+1] = sort_list[states[k]][2]

if 1 in states:
    print("There is a S2")
