######################################################################
# Task c:
# 18-12-2023
################################
import numpy as np
import matplotlib.pyplot as plt
from casadi import* 

from models.model import *
from solvers.solvers import *
from data.Parameters import *
from data.Deterministic_dist_profile import *
from data.Power_goal import *

if __name__ == '__main__':
     # Convert Ws to MWh
    conv_unit_energy = 2.7778e-10
    # Convert seconds to hours
    s2h              = 1/(60*60)

    # Parameters
    A   = p['A']
    a   = p['a']
    g   = p['g']
    rho = p['rho']

    Thor = tf # Time horizon
    N = Nd # number of control intervals

    # Declare model variables
    x = MX.sym('x')
    u = MX.sym('u')
    d = MX.sym('d')
    
    # Inlet flow is the disturbance
    qin = d
    # Compute water levels
    h = x/(rho*A)

    # Potential Outlet flow (if valve is completely open)
    potentialFlow = a*np.sqrt(2*g*h)
        
    # Outlet flow of dam
    qout = potentialFlow*u # The valve between [0, 1] 

    # ODE
    xdot = -rho*qout + rho*qin

    # Weight matrix
    Qu = 1  # Penalization of U
    Qs = 100
    Qz = 10
    # Fixed step Runge-Kutta 4 integrator
    M = 20 # RK4 steps per interval
    DT = Thor/N/M
    f  = Function('f', [x, u, d], [xdot])
    X0 = MX.sym('X0', 1)
    U = MX.sym('U')
    Pslack = MX.sym('Pslack')
    D = MX.sym('D')
    X = X0
    Q = 0
    Uold = 0
    for j in range(M):
        k1 = f(X, U, D)
        k2 = f(X + DT/2 * k1, U, D)
        k3 = f(X + DT/2 * k2, U, D)
        k4 = f(X + DT * k3, U, D)
        X=X+DT/6*(k1 +2*k2 +2*k3 +k4)
    F = Function('F', [X0, U, D], [X],['x0','u', 'd'],['xf'])

    # Evaluate at a test point
    Fk = F(x0=[x0[0]], u = u0[0], d = D_det[0])
    
    # Start with an empty NLP
    w   = []
    w0  = []
    lbw = []
    ubw = []
    phi = 0
    g_MS   = []
    lbg = []
    ubg = []

    # "Lift" initial conditions
    Xk = MX.sym('X0', 1)
    w += [Xk]
    lbw += [x0[0]]
    ubw += [x0[0]]
    w0 += [x0[0]]
    Uold = 0
    # Formulate the NLP
    for k in range(N):
        # New NLP variable for the control
        Uk = MX.sym('U_' + str(k))
        w   += [Uk]
        lbw += [0]
        ubw += [1]
        w0  += [u0[0]]

        # New NLP variable for slack variable
        Sk = MX.sym('S_' + str(k))
        w   += [Sk]
        lbw += [0]  # Slack variable is non-negative
        ubw += [inf]
        w0  += [0]  # Initial guess for the slack variable
        # Get disturbance
        Dk = D_det[k]

        # Integrate till the end of the interval
        Fk = F(x0=Xk, u=Uk, d = Dk)
        Xk_end = Fk['xf']

        # Compute water levels
        h = Xk/(rho*A)

        #   Potential Outlet flow (if valve is completely open)
        potentialFlow = a*np.sqrt(2*g*h)
            
        # Outlet flow of dam
        qout = rho*potentialFlow*Uk  # The valve between [0, 1] 

        # Conversion of outlet flow to power
        Zk =  powCoeff*qout
        Zk_bar = power_goal[k]
        phi = phi + Qu*Uk + Qs * Sk + Qz*Zk
        Uold = Uk

        # New NLP variable for state at end of interval
        Xk = MX.sym('X_' + str(k+1), 1)
        w   += [Xk]
        lbw += [-inf]
        ubw += [inf]
        w0  += [x0[0]]

        # Add equality constraint
        g_MS   += [Xk_end-Xk]
        lbg += [0]
        ubg += [0]

        buffer = 0 # <-- Update this to increase power (input in W)

        # Add inequality constraint: Zk - Zk_bar >= -Sk
        tolerance = inf # Upper tolerance
        g_MS += [Zk - (Zk_bar+buffer)]
        lbg   += [0]
        ubg   += [tolerance]  # Zk - Zk_bar >= -Sk

    # Create an NLP solver
    prob = {'f': phi, 'x': vertcat(*w), 'g': vertcat(*g_MS)}
    solver = nlpsol('solver', 'ipopt', prob, {'ipopt.print_level':0})

    # Solve the NLP
    sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    w_opt = sol['x'].full().flatten()

    # Collect output parameters
    x_opt = np.array(w_opt[0::3])
    u_opt = np.array(w_opt[1::3])
    u_opt = np.ceil(u_opt*100)/100 # Round up to 2 decimals

    Z = np.zeros((1, N))
    for i in range(0, N):
        Z[:, i] = output(0, x_opt[i], u_opt[i], D_det[i], p)

    # Print output
    print("\nSequence of valve configurations: ", u_opt)
    print("Sequence of water levels:         ", np.round(x_opt/(rho*A),2))
    print("Power requirement satisfied:      ",all(np.greater(Z[0, :],power_goal)))

    # *** Plot data ***
    plt.figure(1)
    ax1 = plt.subplot(2, 2, 2)
    plt.plot(T_det*s2h, x_opt/(rho*A),'r--')
    plt.xlabel('Time [h]')
    plt.ylabel('Height [m]')
    plt.title('Predicted water level')

    ax2 = plt.subplot(2, 2, 1)
    plt.step(T_det*s2h, np.append(D_det, D_det[-1]),'r--', where='post')  # Cumulative energy using np.cumsum
    # plt.plot(T_det*s2h, D)  # Cumulative energy using np.cumsum
    plt.xlabel('Time [h]')
    plt.ylabel('Flow [m$^3$/s]')
    plt.title('Predicted inflow')

    ax3 = plt.subplot(2, 2, 3)
    
    #plt.step(T_det*s2h, np.append(power_goal, power_goal[-1])/1e3, 'b--.',where='post', label='Required')
    plt.plot(T_det*s2h, np.append(power_goal, power_goal[-1])/1e3, 'b.',linestyle='dotted', label='Required')
    plt.plot(T_det*s2h, np.append(Z[0, :], Z[0, -1])/1e3,'r--',label='Predicted')
    plt.xlabel('Time [h]')
    plt.ylabel('Power [kW]')
    # plt.xlim([9.01,10])
    plt.title('Generated power')
    ax3.legend()

    ax4 = plt.subplot(2, 2, 4)
    plt.step(T_det*s2h, np.append(u_opt, u_opt[-1]),'k', where='post')
    plt.xlabel('Time [h]')
    plt.ylabel('Valve configuration [-]')
    plt.ylim([-0.1, 1.1])
    plt.title('Input sequence $u(t_k)$')

    # Adjust layout for better spacing
    plt.tight_layout()
    plt.show()