######################################################################
# Deterministic disturbance profile
# 08-11-2023
################################
import numpy as np

# Time data
t0 = 0.0
tf = 60*60*10

# Sequence of inlet flows [m^3/s]
D_det = np.array([17, 17, 18, 20, 25, 22, 21, 22, 18, 17]) # m^3/s

# Length of sequence
Nd    = len(D_det)

# Time intervals for inlet flow sequence
T_det = np.linspace(t0, tf, num=Nd+1)
