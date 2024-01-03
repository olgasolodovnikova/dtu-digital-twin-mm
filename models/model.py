################################
# Model of water damm 2 power
#
# Author: Anders Hilmar Damm Andersen (ahda@dtu.dk)
# 28-11-2023
################################

import numpy as np

def driftModel(t, x, u, d, p):

    # Parameters
    A   = p['A']
    a   = p['a']
    g   = p['g']
    rho = p['rho']

    # Inlet flow is the disturbance
    qin = d

    # Compute water levels
    h = x/(rho*A)

    # Potential Outlet flow (if valve is completely open)
    potentialFlow = a*np.sqrt(2*g*h)

    # Constrain valve
    if(u > 1.0):
        u = 1.0
    if(u < 0.0):
        u = 0.0
        
    # Outlet flow of dam
    qout = potentialFlow*u # The valve between [0, 1] 

    # ODE
    dxdt = -rho*qout + rho*qin
    return dxdt

def diffModel(t, x, u, d, p):
    # Parameters
    rho    =  p['rho']
    sigma  =  p['sigma']

    return rho*sigma


def measurement(t, x, u, d, p):

    # Parameters
    A          = p['A']
    a          = p['a']
    g          = p['g']
    rho        = p['rho']
    powCoeff   = p['powCoeff']

    # Compute water levels
    h = x/(rho*A)

    #   Potential Outlet flow (if valve is completely open)
    potentialFlow = a*np.sqrt(2*g*h)

    # Constrain valve
    if(u > 1.0):
        u = 1.0
    if(u < 0.0):
        u = 0.0
        
    # Outlet flow of dam
    qout = potentialFlow*u # The valve between [0, 1] 

    # Conversion of outlet flow to power
    return powCoeff*rho*qout

def output(t, x, u, d, p):

    # Parameters
    A          = p['A']
    a          = p['a']
    g          = p['g']
    rho        = p['rho']
    powCoeff   = p['powCoeff']

    # Compute water levels
    h = x/(rho*A)

    #   Potential Outlet flow (if valve is completely open)
    potentialFlow = a*np.sqrt(2*g*h)

    # Constrain valve
    if(u > 1.0):
        u = 1.0
    if(u < 0.0):
        u = 0.0
        
    # Outlet flow of dam
    qout = potentialFlow*u # The valve between [0, 1] 

    # Conversion of outlet flow to power
    return powCoeff*rho*qout
    