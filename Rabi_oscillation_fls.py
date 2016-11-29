import matplotlib.pyplot as plt
import numpy as np
from qutip import *

def fls_Rabi(w03, w02, w01, w1, w2, w3, psi0, g1, g2, g3, gamma1, tlist, solver):
    state0 = basis(4,0)
    state1 = basis(4,1)
    state2 = basis(4,2)
    state3 = basis(4,3)
    w12 = w02-w01
    w23 = w03-w02

    delta1 = w03 - w1
    delta2 = w12 - w2
    delta3 = w23 - w3
    delta = delta1 - delta2
    H = delta1*state3*state3.dag() + (delta1-delta3)*state2*state2.dag() + delta*state1*state1.dag() \
        + 0.5*g1*(state0*state3.dag() + state3*state0.dag()) \
        + 0.5*g2*(state1*state2.dag() + state2*state1.dag()) \
        + 0.5*g3*(state3 * state2.dag() + state2 * state3.dag())
    c_ops = []
    if gamma1 > 0.0:
        c_ops.append(np.sqrt(gamma1)*(state0*state1.dag()))
    q_ops = [state0*state0.dag(),state1*state1.dag(), state2*state2.dag(), state3*state3.dag()]
    if solver == "me":
        output = mesolve(H, psi0, tlist, c_ops, q_ops)
    elif solver == "es":
        output = essolve(H, psi0, tlist, c_ops, q_ops)
    elif solver == "mc":
        ntraj = 250
        output = mcsolve(H, psi0, tlist, ntraj, c_ops, q_ops)
    else:
        raise ValueError("unknown solver")
    return output. expect[0], output.expect[1], output.expect[2], output.expect[3]

w03 = 16
w02 = 14
w01 = 2
#so w23 = 2, w12 = 12
w1 = 5
w2 = 1
w3 = 2
g1 = 1
g2 = 1
g3 = 1
gamma1 = 0.0
psi0 = basis(4,0)
time = np.linspace(0,5000,10000)
solver = 'me'

pop_0, pop_1, pop_2, pop_3 = fls_Rabi(w03, w02, w01, w1, w2, w3, psi0, g1, g2, g3, gamma1, time, solver)
plt.plot(time, pop_0, linewidth = '1', linestyle = '-', color = 'blue')
plt.plot(time, pop_1, linewidth = '1', linestyle = '-', color = 'green')
plt.plot(time, pop_2, linewidth = '1', linestyle = '-', color = 'red')
plt.plot(time, pop_3, linewidth = '1', linestyle = '-', color = 'yellow')

delta1 = w03-w1
omega = 0.5*g1*g2*g3/(4.0*delta1**2 - g3**2)
c_anal = (np.cos(omega*time))**2
plt.plot(time, c_anal, linewidth = '2', linestyle = '--', color = 'black')
plt.ylabel('Probability', size = '18')
plt.xlabel('Time', size = '18')
plt.tick_params(labelsize=18)
plt.show()
plt.show()