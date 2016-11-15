import matplotlib.pyplot as plt
import numpy as np
from qutip import *
#'''
def qubit_Rabi(wa, wr, psi0, g, gamma1, tlist, solver):
    H = -(wa-wr)/2.0*sigmaz() + g/2.0*(sigmay())
    c_ops = []
    if gamma1 > 0.0:
        c_ops.append(np.sqrt(gamma1)*sigmam())
    q_ops = [sigmax(),sigmay(),sigmaz()]
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

def tls_Rabi(wa, wr, psi0, g, gamma1, tlist, solver):
    state0 = basis(2,0)
    state1 = basis(2,1)
    H = (wa-wr)*state1*state1.dag() + g*(state0*state1.dag() + state1*state0.dag())
    c_ops = []
    if gamma1 > 0.0:
        c_ops.append(np.sqrt(gamma1)*(state0*state1.dag()))
    q_ops = [state0*state0.dag(),state1*state1.dag()]
    if solver == "me":
        output = mesolve(H, psi0, tlist, c_ops, q_ops)
    elif solver == "es":
        output = essolve(H, psi0, tlist, c_ops, q_ops)
    elif solver == "mc":
        ntraj = 250
        output = mcsolve(H, psi0, tlist, ntraj, c_ops, q_ops)
    else:
        raise ValueError("unknown solver")
    return output. expect[0], output.expect[1]

wa = 10
wr = 10
g = 1
gamma1 = 0.0
psi0 = basis(2,0)
time = np.linspace(0,10,1000)
sx, sy, sz = qubit_Rabi(wr, wa, psi0, g, gamma1, time, "me")
# bSphere = Bloch()
# theta = np.pi/4.0
# phi= np.pi/2.0
# idx = int(len(time)-1)
# arVec = [sx[idx],sy[idx],sz[idx]]
# bSphere.add_vectors(arVec)
# bSphere.add_points([sx,sy,sz])
# bSphere.make_sphere()



pop1, pop2 = tls_Rabi(wr, wa, psi0, g, gamma1, time, "me")
plt.plot(time, pop1, linewidth = '2')
# plt.plot(time, pop2, linewidth = '2')
plt.plot(time, sz, linewidth = '2')
plt.show()

# print sigmay()