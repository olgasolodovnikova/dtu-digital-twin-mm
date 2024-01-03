######################################################################
# Stochastic simulations of water dam 2 power
#
# Author: Anders Hilmar Damm Andersen (ahda@dtu.dk)
# 28-11-2023
################################
import numpy as np
import matplotlib.pyplot as plt

from models.model import *
from solvers.solvers import *
from data.Parameters import *
from data.Deterministic_dist_profile import *

if __name__ == '__main__':
    
    # Convert Ws to MWh
    conv_unit_energy = 2.7778e-10
    # Convert seconds to hours
    s2h              = 1/(60*60)
    
    # Inputs
    u = np.array(u0)

    # Disturbance
    d = np.array(d0) #cm^3/s
    
    # Initial conditions

    x0 = np.array(x0)
    # Simulate 
    t0 = 0.0 #seconds
    tf = 60*60*10 #seconds

    ##### Deterministic simualtions #################################

    # Simulate system from t0 to tf with u and d as inputs
    (T, X) = explicitEulerFixedStepSize(driftModel, t0, tf, x0, u, d, p, opts)
    
    # Generate output sequence
    Z = np.zeros((1, Nsim+1))
    for i in range(0, Nsim+1):
        Z[:, i] = output(T[i], X[:,i], u, d, p)

    plt.figure(1)


    ax1 = plt.subplot(2, 1, 1)
    plt.plot(T*s2h, X[0, :]/(rho*A))
    plt.ylabel('Height in tank [m]')


    # equivalent but more general
    ax1 = plt.subplot(2, 1, 2)

    plt.plot(T*s2h, Z[0, :])
    plt.xlabel('Time [h]')
    plt.ylabel('Power [W]')


    ##### Stochastic simualtions #################################

    # Create Wiener process
    nW = 1
    Ns = 1
    seed = 123
    (W, Tw, dW) = std_wiener_process((tf-t0), Nsim+1, nW, seed)


    # Simulate system from t0 to tf with u and d as inputs
    (T, X) = explicitExplicit(driftModel, diffModel, t0, tf, x0, u, d, p, dW, opts)
    
    # Generate output sequence
    Z = np.zeros((1, Nsim+1))
    for i in range(0, Nsim+1):
        Z[:, i] = output(T[i], X[:,i], u, d, p)

    plt.figure(2)

    ax1 = plt.subplot(2, 1, 1)
    plt.plot(T*s2h, X[0, :]/(rho*A))
    plt.ylabel('Height in tank [m]')

    # equivalent but more general
    ax1 = plt.subplot(2, 1, 2)

    plt.plot(T*s2h, Z[0, :])
    plt.xlabel('Time [h]')
    plt.ylabel('Power [W]')
    plt.tight_layout()
    plt.show()


