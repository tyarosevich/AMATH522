

#%% Part II
import numpy as np
import matplotlib.pyplot as plt

# Declared col. vector of initial populations
n_0 = np.array([[100, 100, 100, 100]])

# Vector of fecundity values, intermediate array
# for survival probabilities
fecun = np.array([0, 1, 5, .5])
p_i = np.diag([.5, .9, .95, 0])

# Adding fecundity onto the survival array and
#removing excess row.
A = np.row_stack([fecun,p_i])
A = np.delete(A, 4, 0)

t_max = 50
t_mesh = np.arange(1, t_max +1)

# Declaring population array and inputting initial condition.
n_t = np.zeros((4, t_max))
n_t[:, 0] = n_0

# Time-steps through the t-mesh starting at t=1
for t in range(1, t_max):
    n_t[:, t] = np.dot(A, n_t[:, t-1])

# Calculates the total population size, i.e.
# the sum of each age-population. Then declares
# an array representing the fraction of each age
# population of the total population at time t.
N_t = np.sum(n_t, 0)
w_t = n_t / N_t


plt.figure(1)
plt.plot(t_mesh, np.log(N_t), 'r--')
plt.xlabel('t')
plt.ylabel('log N(t)')

plt.figure(2)
plt.plot(t_mesh, w_t[0,:], 'r--', label = 'a = 0')
plt.plot(t_mesh, w_t[1,:], 'g--', label = 'a = 1')
plt.plot(t_mesh, w_t[2,:], 'b--', label = 'a = 2')
plt.plot(t_mesh, w_t[3,:], 'y--', label = 'a = 3')
plt.xlabel('t')
plt.ylabel('n_a(t) / N(t)')
plt.legend()

plt.show()

# Polyfits a 1st degree polynomial to the data and
# declares the estimate of lambda.
poly_coeff = np.polyfit(t_mesh, np.log(N_t),1)
lambda_est = np.exp(poly_coeff[0])

# Declares a variable holding the maximum mathematically
# calculated eigenvalue of A.
lambda_dom = np.amax(np.linalg.eigvals(A))
# #plt.legend()

#%% Part III
import numpy as np
import matplotlib.pyplot as plt

# b.) and c.)
#Constructing A matrix and inputing values
fecun_vec = .24 * np.ones(50)
fecun_vec[0:3] = 0
surv_mat = np.diag(.952 * np.ones(50))
surv_mat[0, 0] = 1
surv_mat[1,1,] = 1
surv_mat[2,2] = .0722

A = np.row_stack([fecun_vec,surv_mat])
A = np.delete(A, 50, 0)

# Setting up discrete time values and associated
# population value array with the n_0 values.
t_max = 100
t_mesh = np.arange(1, t_max + 1)
n_0 = 10 * np.ones(50)
n_t = np.zeros((50, t_max))
n_t[:,0] = n_0

# Iterating through the time steps.
for t in range(1, t_max):
    n_t[:, t] = np.dot(A, n_t[:, t-1])

# Total popluation for each time step.
N_t = np.sum(n_t, 0)
# plt.figure(3)
# plt.plot(t_mesh, np.log(N_t), 'r--')
# plt.xlabel('t')
# plt.ylabel('log N(t)')

# Estimating dominant eigenvalue and calculating directly.
poly_coeff = np.polyfit(t_mesh, np.log(N_t),1)
lambda_est = np.exp(poly_coeff[0])
lambda_dom = np.amax(np.linalg.eigvals(A))

#part d.)

# Storing the left and right eigenvectors, indexing the max
# eigenvalue. Taking the real values of the eigenvectors since the matrix
# A is power-positive.
w, vr = np.linalg.eig(A)
w2, vl = np.linalg.eig(A.T)
max_eig_index = np.argmax(w)
vr = vr[:,max_eig_index].real
vl = vl[:,max_eig_index].real

# Using an outer product of l/r eigenvectors and dividing
# by their dot product (scalar) to calculate sensitivity matrix.
# Weighting each entry to calculate the elasticity.
S_ij = (np.outer(vl, vr) / np.dot(vl,vr)).real
e_ij =  np.multiply( A / w[max_eig_index].real, S_ij)


#%%
# Part IV For Juvinilles
import numpy as np
import matplotlib.pyplot as plt

# Declaring the Leslie Matrix
fecun = np.array([0, .0043, .1132, 0])
prob = np.array([[.9775, .9111, 0, 0], [0, .0736, .9534, 0], [0, 0, .0452, .9804]])
A = np.row_stack([fecun,prob])

# Declaring eigenvalue and eigenvectors in orer to find
# the dominant eigenvectors in order use the ratios therein
# to find a stable population distribution with a total of 250
# individuals.
w, v = np.linalg.eig(A)
max_eig_index = w.argmax()
eigvec_stable = v[:,max_eig_index]
eig_vec_sum = np.sum(eigvec_stable)
n_0 = (eigvec_stable / eig_vec_sum) * 250

t_max = 500
t_mesh = np.arange(1, t_max +1)

# Declaring population array and inputting initial condition. Putting
# a placeholder value in the total population array in order to start
# the while loop. Holder variable for the sustainable h value, as well
# as the initial population.
n_t = np.zeros((4, t_max))
n_t[:, 0] = n_0
n_t_start = n_t
N_t = np.sum(n_t, 0)
N_t[-1] = 2
h_sustainable = 0
h_list = np.zeros(4)

# Note the while loop condition is that there be greater than or
# equal to 1 individual, since presumable 99% of an Orca is not an Orca.
while N_t[-1] >= 1:
    h_sustainable = h_list[1]
    h_list[1] = h_list[1] + 1
    n_t = n_t_start
    for t in range(1, t_max):
        n_t[:, t] = np.dot(A, n_t[:, t - 1]) - h_list
        if n_t[1, t] < 0:
            n_t[1, t] = 0
    N_t = np.sum(n_t, 0)

#%%
# Part IV For reproductive adults
import numpy as np
import matplotlib.pyplot as plt

fecun = np.array([0, .0043, .1132, 0])
prob = np.array([[.9775, .9111, 0, 0], [0, .0736, .9534, 0], [0, 0, .0452, .9804]])
A = np.row_stack([fecun,prob])

w, v = np.linalg.eig(A)
max_eig_index = w.argmax()
eigvec_stable = v[:,max_eig_index]
eig_vec_sum = np.sum(eigvec_stable)
n_0 = (eigvec_stable / eig_vec_sum) * 250

t_max = 500
t_mesh = np.arange(1, t_max +1)

# Declaring population array and inputting initial condition.
n_t = np.zeros((4, t_max))
n_t[:, 0] = n_0
n_t_start = n_t
N_t = np.sum(n_t, 0)
N_t[-1] = 2
h_sustainable = 0
h_list = np.zeros(4)

while N_t[-1] >= 1:
    h_sustainable = h_list[2]
    h_list[2] = h_list[2] + 1
    n_t = n_t_start
    for t in range(1, t_max):
        n_t[:, t] = np.dot(A, n_t[:, t - 1]) - h_list
        if n_t[2, t] < 0:
            n_t[2, t] = 0
    N_t = np.sum(n_t, 0)
