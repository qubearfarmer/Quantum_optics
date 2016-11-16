import matplotlib.pyplot as plt
import numpy as np
from qutip import *

def thls_Rabi(w02, w01, w1, w2, psi0, g1, g2, gamma1, tlist, solver):
    state0 = basis(3,0)
    state1 = basis(3,1)
    state2 = basis(3,2)
    delta1 = w02 - w1
    w12 = w02-w01
    delta2 = w12 - w2
    delta = delta1 - delta2
    H = delta1*state2*state2.dag() + delta*state1*state1.dag() + 0.5*g1*(state0*state2.dag() + state2*state0.dag()) + 0.5*g2*(state1*state2.dag() + state2*state1.dag())
    c_ops = []
    if gamma1 > 0.0:
        c_ops.append(np.sqrt(gamma1)*(state0*state1.dag()))
    q_ops = [state0*state0.dag(),state1*state1.dag(), state2*state2.dag()]
    if solver == "me":
        output = mesolve(H, psi0, tlist, c_ops, q_ops)
    elif solver == "es":
        output = essolve(H, psi0, tlist, c_ops, q_ops)
    elif solver == "mc":
        ntraj = 250
        output = mcsolve(H, psi0, tlist, ntraj, c_ops, q_ops)
    else:
        raise ValueError("unknown solver")
    return output. expect[0], output.expect[1], output.expect[2]

w02 = 10
w01 = 2
#w12 = 8
w1 = 5
w2 = 3
g1 = 0.5
g2 = 0.5
gamma1 = 0.0
psi0 = basis(3,0)
time = np.linspace(0,500,1000)

pop_0, pop_1, pop_2 = thls_Rabi(w02, w01, w1, w2, psi0, g1, g2, gamma1, time, 'me')
plt.plot(time, pop_0, linewidth = '2', linestyle = '-', color = 'blue')
plt.plot(time, pop_1, linewidth = '2', linestyle = '-', color = 'green')
plt.plot(time, pop_2, linewidth = '2', linestyle = '-', color = 'red')

delta1 = w02-w1
omega = g1*g2/(4.0*delta1)
c_anal = (np.cos(omega*time))**2
plt.plot(time, c_anal, linewidth = '2', linestyle = '--', color = 'black')
plt.ylabel('Probability', size = '18')
plt.xlabel('Time', size = '18')
plt.tick_params(labelsize=18)
plt.show()


