################################
# ODE solvers
#
# Author: Anders Hilmar Damm Andersen (ahda@dtu.dk)
# 28-11-2023
################################

import numpy as np

def explicitEulerFixedStepSize(fun, t0, tf, x0, u, d, p, opts):

    # Extract solver options
    N = opts['N']

    # Stepsize h
    h = (tf-t0)/N

    # Assign initial conditions
    tk = t0
    xk = x0

    # Dimensions
    nx = len(x0)

    # Output data
    T = np.zeros(N+1)       # Initialize vector of zeros for time
    X = np.zeros((nx, N+1))  # Initialize matrix of zeros for states
    
    # Update data
    T[0]    = tk  # Store the current state in the matrix
    X[:, 0] = xk  # Store the current state in the matrix
    
    for k in range(1, N+1):
        
        xk = xk + h*fun(tk, xk, u, d, p)
        tk = tk + h
        T[k]    = tk  
        X[:, k] = xk  

    return T, X


def std_wiener_process(T, N, nW, seed=None):
    if seed is not None:
        np.random.seed(seed)

    dt = T / N #kind of a big number for our D_det values
    dW = np.sqrt(dt) * np.random.randn(nW, N) # took the sqrt of sqrt(dt)
    
    Tw = np.arange(0, T + dt, dt)
    W = np.concatenate((np.zeros((nW, 1)), np.cumsum(dW, axis=1)), axis=1)

    return W, Tw, dW

def explicitExplicit(fun, sigmafun, t0, tf, x0, u, d, p, dW, opts):
    ### SDE solver for single realizations ###

    # Extract solver options
    N = opts['N']

    # Stepsize h
    h = (tf-t0)/N

    # Assign initial conditions
    tk = t0
    xk = x0

    # Dimensions
    nx = len(x0)

    # Output data
    T = np.zeros(N+1)       # Initialize vector of zeros for time
    X = np.zeros((nx, N+1))  # Initialize matrix of zeros for states
    
    # Update data
    T[0]    = tk  # Store the current state in the matrix
    X[:, 0] = xk  # Store the current state in the matrix
    
    for k in range(1, N+1):
        
        f = fun(tk, xk, u, d, p)
        g = sigmafun(tk, xk, u, d, p)
        phi = xk + g*dW[:, k]
        xk = phi + f*h
        tk = tk + h
        T[k]    = tk  
        X[:, k] = xk  

    return T, X