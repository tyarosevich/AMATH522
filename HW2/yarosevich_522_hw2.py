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

# This seems to be correct. There was something wrong with the last one.
for k in np.arange(0, numsteps - 1):
    if states[k] == 0:
        states[k+1] = np.random.choice(3, 1, p = [.98, .02, 0])
    if states[k] == 1:
        states[k+1] = np.random.choice(3, 1, p = [.1, .7, .2])
    if states[k] == 2:
        states[k+1] = np.random.choice(3, 1, p = [0, .05, .95])

#for k in np.arange(0, numsteps-1):
#    rand = rd.random()
#    if rand < A[sort_list[states[k]][0], states[k]]:
#        states[k+1] = sort_list[states[k]][0]
#    elif rand < A[sort_list[states[k]][1], states[k]]:
#        states[k+1] = sort_list[states[k]][1]
#    else:
#        states[k+1] = sort_list[states[k]][2]


reduced_states = states[:]
for i in range(0, len(reduced_states)):
    if reduced_states[i] == 1:
        reduced_states[i] = 0
for i in range(0, len(reduced_states)):
    if reduced_states[i] == 2:
        reduced_states[i] = 1

last_dwell = 0
count = 1
for i in range(1, len(reduced_states)):
    if reduced_states[i] != reduced_states[i-1]:
        dwell_list.append(count)
        count = 1
    count = count +1
print("Sim done")

#%%
#unique, counts = np.unique(states, return_counts=True)
#dict(zip(unique, counts))


#%%

n_bins = 30
n, bins, patches = plt.hist(dwell_list, n_bins, facecolor='blue', alpha=0.5)
plt.show()
#%%
log_states = np.log(dwell_list)


n_bins = 30
n, bins, patches = plt.hist(dwell_list, n_bins, facecolor='blue', alpha=0.5)
plt.yscale('log')
plt.xlabel('Dwell Bins')
plt.ylabel('Dwells per Bin')
plt.show()

#%%
# Part 4.)
import numpy as np
import matplotlib.pyplot as plt
import random as rd


N_in = 1000
N_out = 500
p_in = .4
p_out = .6
numsteps = 1000
T = np.arange(300)
p_T = np.zeros(len(T))
current_list = np.zeros(numsteps)
mean = N_in*p_in - N_out*p_out

# Simulate realizations compared against each threshhold T
for k in np.arange(0, len(T)):
    spike = 0.0
    
    # Simulate lots of realizations of N
    for i in np.arange(0, numsteps):
        # Generates a coin toss vector for the inward and outward channels
        # in which 'open' corresponds to 1. Sums the vector to ascertain
        # the number of open channels in a realization.
        N_in_open = np.sum(np.random.choice([1, 0], size = N_in, p = [p_in, 1 - p_in]))
        N_out_open = np.sum(np.random.choice([1,0], size = N_out, p = [p_out, 1 - p_out]))
        # The net current, open inward channels minus open outward.
        net_current = N_in_open - N_out_open
        # IF the net current is more than T, a spike occurs
        if net_current > k:
            spike = spike + 1.0
        if k == 1:
            current_list[i] = net_current
            
    # Once numsteps realizations have been simulated, divides the number
    # of spikes by the total number of realizations to get the probability
    # that a spike will occur with a given threshhold T
    p_T[k] = spike/numsteps
print("Sim done")
    
#%%

plt.figure(1)
plt.plot(T, p_T)
plt.xlabel('T*E(I)')
plt.ylabel('Probability of Spike')
plt.axis([0, 5*mean, 0, 1])

plt.show()