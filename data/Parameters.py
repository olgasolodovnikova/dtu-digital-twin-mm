## Parameters for simulation

#  --- parameters ---

# Model parameters
A = 300         # m^2
a = 3.0         # m^2
g = 0.981       # kg/s^2
rho = 998       # kg/m^3
powCoeff = 10.0 # W/(m^3/s)]

sigma = 1    # [m^3/s -> kg]
nW    = 1    # Number of noise sources --> only 1
seed  = 12   # Seed for random number generator

# Initial conditions
h0    = 10           # Water level in tank [m]
x0    = [h0*(rho*A)] # Water mass in tank [kg]
d0    = [20]         # Inlet flow [m^3/s]
u0    = [1]          # Valve position [-]

# Parameter vector
p = {'A': A, 'a': a, 'g': g, 'rho': rho, 'powCoeff': powCoeff, 'sigma': sigma}

# Simulation parameters
#Nsim = 200           # Number of steps for solvers

Nsim = 10       # Number of steps for solvers
opts = {'N': Nsim}