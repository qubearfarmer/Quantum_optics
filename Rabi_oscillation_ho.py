import matplotlib.pyplot as plt
import numpy as np
from qutip import *
#'''
def ho_Rabi(N, wa, wr, psi0, g, gamma1, tlist, solver):
    a = destroy(N)
    H = (wa-wr)*a.dag()*a + g*(a.dag() + a)
    c_ops = []
    if gamma1 > 0.0:
        c_ops.append(np.sqrt(gamma1)*sigmam())
    q_ops = [a.dag()*a]
    if solver == "me":
        output = mesolve(H, psi0, tlist, c_ops, q_ops)
    elif solver == "es":
        output = essolve(H, psi0, tlist, c_ops, q_ops)
    elif solver == "mc":
        ntraj = 250
        output = mcsolve(H, psi0, tlist, ntraj, c_ops, q_ops)
    else:
        raise ValueError("unknown solver")
    return output. expect[0]

N=2
psi0 = basis(N,0)
time = np.linspace(0,10,100)
gamma1 = 0
solver = "me"

wa = 10
wr = 10
g = 1
excite = ho_Rabi(N, wa, wr, psi0, g, gamma1, time, solver)
plt.plot(time,excite, linewidth = '2')
# wa = 10
# wr = 9
# g = 1
# excite = ho_Rabi(N, wa, wr, psi0, g, gamma1, time, solver)
# plt.plot(time,excite, linewidth = '2')
# wa = 10
# wr = 8
# g = 1
# excite = ho_Rabi(N, wa, wr, psi0, g, gamma1, time, solver)
# plt.plot(time,excite, linewidth = '2')

plt.xlabel('Time', size = 18)
plt.ylabel('Probability of having n excitation', size = 18)
plt.tick_params(labelsize=18)
plt.tick_params(labelsize=18)
plt.show()