######################################################################
# Task 1: Electricity production 
# What the the maximum amount of power generated for a given inflow profile?
#
# 08-11-2023
################################
import numpy as np
import matplotlib.pyplot as plt

from models.model import *
from solvers.solvers import *
from data.Parameters import *
from data.Deterministic_dist_profile import *

# *** Constants ***

# Convert Ws to MWh
conv_unit_energy = 2.7778e-10

# Convert seconds to hours
s2h              = 1/(60*60)


if __name__ == '__main__':

    # *** Valve settings ***
    # Insert you value for u0
    u0 = 0.5 # <---- substitude with your choice
    print(f'Constant valve setting u = {u0}')
    u = np.array(u0)

    # *** Simulation data ***
    # Initial conditions
    x0 = np.array(x0)

    # Numpy matrices for storing data
    T = np.zeros(Nd*(Nsim)) # time
    X = np.zeros(Nd*(Nsim)) # water mass in tank
    Z = np.zeros(Nd*(Nsim)) # generated power
    Y = np.zeros(Nd*(Nsim)) # measured power affected by measurement noise in discrete time - why not D affected by noise? 
    D = np.zeros(Nd*(Nsim)) # disturbance (inlet flow)

    D_det_plot = np.zeros(Nd*(Nsim)) #for plotting expected inflow in fig1

    # *** Simulation ***
    # Generate stochastics data solving the system of stochastic differential equations (SDEs).
    # We split the total simulation into Nd sub-simulations. The inlet flow is constant in each sub-simulation.

    # Generate stochastic data from a standard Wiener process
    (W, Tw, dW) = std_wiener_process((tf-t0), (Nsim+1)*Nd, nW, seed) 

    # Measurement noise
    Rk = 100000 # OBS: Maybe we need to change this
    
    vk = np.sqrt(Rk) * np.random.randn(1, Nd*(Nsim))
    # Set initial state for each sub-simulation
    xk = x0
    for i in range(0, Nd):

        # Compute the start and end indices for the current block of simulation data
        start_idx = i * (Nsim)
        end_idx = (i + 1) * (Nsim)

        D[start_idx:end_idx] = np.tile(D_det[i], (Nsim)) + sigma * dW[:,i * (Nsim):(i + 1) * (Nsim)][0]

        # Solve the SDE with inflet flow D_det[i] (constant) and with initial condition xk.
        (t, x) = explicitExplicit(driftModel, diffModel, T_det[i], T_det[i+1], xk, u, D_det[i], p, dW[:,i * (Nsim+1):(i + 1) * (Nsim+1)], opts)
        

        #Solve the SDE with inflet flow D[i] (stochastic) and with initial condition xk.
        #(t, x) = explicitExplicit(driftModel, diffModel, Tw[i], Tw[i+1], xk, u, D[i], p, dW[:,i * (Nsim+1):(i + 1) * (Nsim+1)], opts)
        
        # Update initial conditions for next sub-simulation
        xk = x[:,-1]
        
        # Insert simulation data t and x into the corresponding block in T and X
        T[start_idx:end_idx] = t[:-1]
        X[start_idx:end_idx] = x[:, :-1]
        #D[start_idx:end_idx] = np.tile(D_det[i], (Nsim)) + dW[:,i * (Nsim):(i + 1) * (Nsim)][0]

        D_det_plot[start_idx:end_idx] = np.tile(D_det[i], (Nsim))
        
        # Generate output sequence (electrical power [W]).
        z = np.zeros((1, Nsim))
        for j in range(0, Nsim):
            z[:,j] = output(t[j], x[:, j], u, D_det[i], p) 
            #z[:,j] = output(t[j], x[:, j], u, D[i], p) #Use the stochastic inflow
        Z[start_idx:end_idx] = z
        Y[start_idx:end_idx] = z

    # Add measurement noise to Y

    Y =  Y + vk[0] #OBS: why is it added to Y, not to D? 
    
    # *** Compute the total amount of energy produced ***
    # We use the total simulation data of Z to compute the accumulated power.
    # To convert this power into total energy by multiplying with the time.
    # The time-intervals are of equal length. We compute one time interval as:
    dt = T[1] - T[0]

    # Sum all simulation data in Z and multiply by the time-interval. We convert the result from Ws to MWh
    energy = np.sum(Z) * dt *conv_unit_energy, #OBS: why not use Y to calculate stochastic power?
    energy_Y = np.sum(Y) * dt *conv_unit_energy
    # Print total electrical energy produced
    print("Total electrical energy produced after {:.1f} hours is {:.3f} MWh".format((tf)*s2h, energy[0]))

    print(f'Deterministic energy: {energy}, stochastic energy: {energy_Y}')
    

    # *** Plot data ***
    plt.figure(1)

    # Subplot for disturbances
    ax2 = plt.subplot(2, 2, 1)
    plt.plot(T*s2h, D) 
    plt.plot(T*s2h, D_det_plot, '--', label ='Predicted inflow')
    plt.xlabel('Time [h]')
    plt.ylabel('Flow [m^3/s]')
    plt.title('Inflow')
    plt.legend()

    # Subplot for Tank water levels
    ax1 = plt.subplot(2, 2, 2)
    plt.plot(T*s2h, X/(rho*A))
    plt.ylabel('Height in tank [m]')
    plt.xlabel('Time [h]')
    plt.title('Dam water levels')

    

    # Subplot for electrical power
    ax3 = plt.subplot(2, 2, 3)
    plt.plot(T*s2h, Y/1e6) #scale to MW
    plt.xlabel('Time [h]')
    plt.ylabel('Power [MW]')
    plt.title(r'Generated power $z(t_k)$')

    # Subplot for accumulated energy
    ax4 = plt.subplot(2, 2, 4)
    plt.plot(T*s2h, np.cumsum(Z) * dt *conv_unit_energy,label='Sum')  # Cumulative energy using np.cumsum
    plt.axhline(y=energy, color='r', linestyle='--', label='Total produced energy')
    plt.xlabel('Time [h]')
    plt.ylabel('Energy [MWh]')
    ax4.legend(loc='lower left')
    plt.title('Energy production')

    # Adjust layout for better spacing
    plt.tight_layout()

    # Show plot
    plt.show()
