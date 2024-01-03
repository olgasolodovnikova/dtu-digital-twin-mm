######################################################################
# Task b:
# Manually set valve configuration profile --> How can you reach the high demand?
#
# 10-11-2023
################################
import numpy as np
import matplotlib.pyplot as plt

from models.model import *
from solvers.solvers import *
from data.Parameters import *
from data.Deterministic_dist_profile import *
from data.Power_goal import *

 # Convert Ws to MWh
conv_unit_energy = 2.7778e-10

# Convert seconds to hours
s2h              = 1/(60*60)


if __name__ == '__main__':

    # ----------- Task: choose sequence of valve positions (between 0 and 1) -----------
    # Insert you sequence of value configurations (You need to specify Nd numbers of values)
    u0 = [0.3,0.2,0.2,0.25,0.35,0.4,0.4,0.41,0.48,0.5] # <---- insert here
    u0 = [1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ]
    u = np.array(u0)

    # Initial conditions
    x0 = np.array(x0)

    # Numpy matrices for storing data
    T = np.zeros(Nd*(Nsim))
    X = np.zeros(Nd*(Nsim))
    Z = np.zeros(Nd*(Nsim))
    Y = np.zeros(Nd*(Nsim))
    D = np.zeros(Nd*(Nsim))

    # *** Simulation ***
    # Generate stochastics data solving the system of stochastic differential equations (SDEs).
    # We split the total simulation into Nd sub-simulations. The inlet flow is constant in each sub-simulation.

    # Generate stochastic data from a standard Wiener process
    (W, Tw, dW) = std_wiener_process((tf-t0), (Nsim+1)*Nd, nW, seed)

    # Measurement noise
    Rk = 100000 # Maybe we need to change this
    vk = np.sqrt(Rk) * np.random.randn(1, Nd*(Nsim))
    # Set initial state for each sub-simulation
    xk = x0

    for i in range(0, Nd):

        # Compute the start and end indices for the current block of simulation data
        start_idx = i * (Nsim)
        end_idx = (i + 1) * (Nsim)

        # Solve the SDE with inflet flow D_det[i] (constant) and with initial condition xk.
        (t, x) = explicitExplicit(driftModel, diffModel, T_det[i], T_det[i+1], xk, u[i], D_det[i], p, dW[:,i * (Nsim+1):(i + 1) * (Nsim+1)], opts)
        # Update initial conditions for next simulation
        xk = x[:,-1]
        
         # Insert simulation data t and x into the corresponding block in T and X
        T[start_idx:end_idx] = t[:-1]
        X[start_idx:end_idx] = x[:, :-1]
        D[start_idx:end_idx] = np.tile(D_det[i], (Nsim)) + sigma*dW[:,i * (Nsim):(i + 1) * (Nsim)][0]
        
        # Generate output sequence (electrical power [W]).
        z = np.zeros((1, Nsim))
        for j in range(0, Nsim):
            z[:,j] = output(t[j], x[:, j], u[i], D_det[i], p)
        Z[start_idx:end_idx] = z
        Y[start_idx:end_idx] = z
    
    # Add measurement noise to Y
    Y =  Y + vk[0]

    # *** Plot data ***
    
    # Subplot for Tank water levels
    ax1 = plt.subplot(2, 2, 1)
    plt.plot(T*s2h, X/(rho*A))
    plt.ylabel('Height in tank [m]')
    plt.xlabel('Time [h]')
    plt.title('Dam water levels')

    # Subplot for disturbances
    ax2 = plt.subplot(2, 2, 2)
    plt.plot(T*s2h, D)  # Cumulative energy using np.cumsum
    plt.xlabel('Time [h]')
    plt.ylabel('Flow [m^3]/s')
    plt.title('Inflow')

    # Subplot for produced power and minimum required power.
    ax3 = plt.subplot(2, 2, 3)
    plt.plot(T*s2h, Z, label='Produced')
    plt.step(T_det*s2h, np.append(power_goal, power_goal[-1]), where='post', label='Min.')
    plt.xlabel('Time [h]')
    plt.ylabel('Power [W]')
    plt.title('Power from flow')
    ax3.legend()

    # Subplot for valve configurations
    ax4 = plt.subplot(2, 2, 4)
    plt.step(T_det*s2h, np.append(u, u[-1]), where='post')  # Cumulative energy using np.cumsum
    plt.xlabel('Time [h]')
    plt.ylabel('Valve configurations')
    plt.ylim([-0.9, 1.1])
    plt.title('Input sequence')


    # Adjust layout for better spacing
    plt.tight_layout()


    plt.show()
