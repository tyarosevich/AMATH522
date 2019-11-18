#%% Homework 3 Part 1
import numpy as np
import matplotlib.pyplot as plt
import random as rd
from scipy.integrate import solve_ivp


def repress(t, y, p):
   dy = np.zeros(6)
   dy[0] = -y[0] + p[0]/[1.+y[5]**p[3]]+ p[1]
   dy[1] = -y[1] + p[0]/[1.+y[3]**p[3]]+ p[1]
   dy[2] = -y[2] + p[0]/[1.+y[4]**p[3]]+ p[1]
   dy[3] = -p[2]*[y[3]-y[0]]
   dy[4] = -p[2]*[y[4]-y[1]]
   dy[5] = -p[2]*[y[5]-y[2]]
   return dy

alpha = 0
alpha0 = .5
beta = 1
n=2

p = np.array([alpha, alpha0, beta, n])
Tmax = 100
y0 = 30*np.random.uniform(0,1,6)

#%%
sol = solve_ivp(fun = lambda t,y: repress(t,y,p), [0, Tmax], y0, method = 'RK45')


# %%
