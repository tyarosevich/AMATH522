

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

#plt.show()

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

fecun_vec = .24 * np.ones(50)
fecun_vec[0:3] = 0
surv_mat = np.diag(.952 * np.ones(50))
surv_mat[0, 0] = 1
surv_mat[1,1,] = 1
surv_mat[2,2] = .0722

A = np.row_stack([fecun_vec,surv_mat])
A = np.delete(A, 50, 0)
t_max = 100
t_mesh = np.arange(1, t_max + 1)
n_0 = 10 * np.ones(50)
n_t = np.zeros((50, t_max))
n_t[:,0] = n_0

for t in range(1, t_max):
    n_t[:, t] = np.dot(A, n_t[:, t-1])

N_t = np.sum(n_t, 0)
plt.figure(3)
plt.plot(t_mesh, np.log(N_t), 'r--')
plt.xlabel('t')
plt.ylabel('log N(t)')


poly_coeff = np.polyfit(t_mesh, np.log(N_t),1)
lambda_est = np.exp(poly_coeff[0])

lambda_dom = np.amax(np.linalg.eigvals(A))

#part d.)

w, vr = np.linalg.eig(A)
vr = vr[:,1]

w, vl = np.linalg.eig(A.T)
vl = vl[:,1]

S_ij = np.outer(vl, vr) / np.dot(vl,vr)