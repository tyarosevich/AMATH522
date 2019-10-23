#%%
import numpy as np
import matplotlib.pyplot as plt
import random as rd

numsteps = 100000
dwell_list = []

A = np.array([[.98, .1, 0], [.02, .7, .05], [0, .2, .95]])
c1_sort = np.argsort(A[:,0])
c2_sort = np.argsort(A[:,1])
c3_sort = np.argsort(A[:,2])
sort_list = [c1_sort, c2_sort, c3_sort]


states = np.zeros(numsteps, dtype = 'int')
states[0] = 0


#for k in np.arange(0, numsteps - 1):
#    j = states[k]
#    rand = rd.random()
#    if j == 0:
#        if rand < .02:
#            states[k+1] = 1
#        else:
#            states[k+1] = 0
#    elif j == 1:
#        if rand < .1:
#            states[k+1] = 0
#        elif rand < .2:
#            states[k+1] = 2
#        else:
#            states[k+1] = 1
#    elif j == 2:
#        if rand < .05:
#            states[k+1] = 1
#        else:
#            states[k+1] = 2
for k in np.arange(0, numsteps-1):
    rand = rd.random()
    if rand < A[sort_list[states[k]][0], states[k]]:
        states[k+1] = sort_list[states[k]][0]
    elif rand < A[sort_list[states[k]][1], states[k]]:
        states[k+1] = sort_list[states[k]][1]
    else:
        states[k+1] = sort_list[states[k]][2]


#reduced_states = states[:]
#for i in range(0, len(reduced_states)):
#    if reduced_states[i] == 1:
#        reduced_states[i] = 0
#for i in range(0, len(reduced_states)):
#    if reduced_states[i] == 2:
#        reduced_states[i] = 1
#
#last_dwell = 0
#count = 1
#for i in range(1, len(reduced_states)):
#    if reduced_states[i] != reduced_states[i-1]:
#        dwell_list.append(count)
#        count = 1
#    count = count +1
print("Sim done")

#%%
unique, counts = np.unique(states, return_counts=True)
dict(zip(unique, counts))2


#%%

n_bins = 30
n, bins, patches = plt.hist(dwell_list, n_bins, facecolor='blue', alpha=0.5)
plt.show()
#%%
log_states = np.log(dwell_list)


n_bins = 30
n, bins, patches = plt.hist(dwell_list, n_bins, facecolor='blue', alpha=0.5)
plt.yscale('log', nonposy = 'clip')
plt.show()








