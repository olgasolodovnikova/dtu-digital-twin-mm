######################################################################
# Task b:
# 10-11-2023
################################
import numpy as np
import matplotlib.pyplot as plt

from models.model import *
from solvers.solvers import *
from data.Parameters import *
from data.Deterministic_dist_profile import *
from data.Power_goal import *

if __name__ == '__main__':
    # *** Valve settings ***
    # Insert you sequence of value configurations (You need to specify Nd = 10 values)
    
    u0 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] # <---- insert here
    u0  = [0.28, 0.44, 0.67, 0.67 ,0.6,  0.76, 0.93, 0.94 ,0.91 ,0.93]
    # *** Constants ***
     # Convert Ws to MWh
    conv_unit_energy = 2.7778e-10

    # Convert seconds to hours
    s2h              = 1/(60*60)

    # *** Simulation ***
    # Generate stochastics data solving the system of stochastic differential equations (SDEs).
    # We split the total simulation into Nd sub-simulations. The inlet flow is constant in each sub-simulation.

    u = np.array(u0)

    # Initial conditions
    x0 = np.array(x0)

    # Numpy matrices for storing data
    T = np.zeros(Nd*(Nsim)) # time
    X = np.zeros(Nd*(Nsim)) # water mass in tank
    Z = np.zeros(Nd*(Nsim)) # generated power
    D = np.zeros(Nd*(Nsim)) # disturbance (inlet flow)
    D_det_plot = np.zeros(Nd*(Nsim)) #for plotting expected inflow in fig1

    # Generate stochastic data from a standard Wiener process
    (W, Tw, dW) = std_wiener_process((tf-t0), (Nsim+1)*Nd, nW, seed)

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
        
        # Generate output sequence (electrical power [W]).
        z = np.zeros((1, Nsim))
        for j in range(0, Nsim):
            z[:,j] = output(t[j], x[:, j], u[i], D_det[i], p)
        Z[start_idx:end_idx] = z

        # Save deterministic (predected) stochastic (actual) inlet flow 
        D[start_idx:end_idx] = np.tile(D_det[i], (Nsim)) + sigma*dW[:,i * (Nsim):(i + 1) * (Nsim)][0]
        D_det_plot[start_idx:end_idx] = np.tile(D_det[i], (Nsim))
    
    # Print output
    print("\nSequence of valve configurations: ", u)
    print("Sequence of water levels:         ", np.round(X[0:10:]/(rho*A),2))
    print("Power requirement satisfied:      ",all(np.greater(np.max(np.reshape(Z,(Nd,Nsim)),axis=1),power_goal)))

    # *** Plot data ***
    plt.figure(1)

    # Subplot for disturbances
    ax2 = plt.subplot(2, 2, 1)
    plt.plot(T*s2h, D, 'k',label = 'Measured') 
    plt.plot(T*s2h, D_det_plot, 'r--', label = 'Predicted')
    plt.xlabel('Time [h]')
    plt.ylabel('Flow [m$^3$]/s')
    plt.title('Inflow')
    plt.legend()

    # Subplot for Tank water levels
    ax1 = plt.subplot(2, 2, 2)
    plt.plot(T*s2h, X/(rho*A),'k')
    plt.ylabel('Height [m]')
    plt.xlabel('Time [h]')
    plt.title('Dam water level $h(t_k)$')

    # Subplot for electrical power
    ax3 = plt.subplot(2, 2, 3)
    plt.plot(T*s2h, Z/1e3, 'k',label='Measured')
    plt.step(T_det*s2h, np.append(power_goal, power_goal[-1])/1e3,'-.b', where='post', label='Required')
    plt.xlabel('Time [h]')
    plt.ylabel('Power [kW]')
    plt.title('Generated power $z(t_k)$')
    ax3.legend()

    # Subplot for valve configurations
    ax4 = plt.subplot(2, 2, 4)
    plt.step(T_det*s2h, np.append(u, u[-1]),'k', where='post')  
    plt.xlabel('Time [h]')
    plt.ylabel('Valve configuration [-]')
    plt.ylim([-0.1, 1.1])
    plt.title('Input sequence $u(t_k)$')
 
    # Adjust layout for better spacing
    plt.tight_layout()

    # Show plot
    plt.show()
