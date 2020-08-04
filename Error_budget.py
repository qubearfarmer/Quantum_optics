import numpy as np

#decoherence error
t_1 = 9.7e-6
t_2 = 7.1e-6
t_gate = 40e-9
t_phi = (t_2**-1 - (2*t_1)**-1)**-1
t_error = (t_1**-1 + t_phi**-1)**-1*3
error = t_gate/t_error

Qin = 2000
Qout = 6000
Q = (Qin**-1 + Qout**-1)**-1
kappa = 7.5e3/Q
print (kappa)