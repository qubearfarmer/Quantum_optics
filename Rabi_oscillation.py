import matplotlib.pyplot as plt
import numpy as np
from qutip import *

def qubit_Rabi(wa, wr, psi0, g, gamma1, tlist, solver):
    H = (wa-wr)*sigmaz() - g*sigmax()
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

wa = 10
wr = 10
g = 2
gamma1 = 2
psi0 = basis(2,0)
time = np.linspace(0,20,1000)
sx, sy, sz = qubit_Rabi(wr, wa, psi0, g, gamma1, time, "me")
bSphere = Bloch()
theta = np.pi/4.0
phi= np.pi/2.0
idx = int(len(time)-1)
arVec = [sx[idx],sy[idx],sz[idx]]
bSphere.add_vectors(arVec)
bSphere.add_points([sx,sy,sz])
bSphere.make_sphere()
plt.show()
# plt.plot(time, sz)
# plt.show()