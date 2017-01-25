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
    H = (wa-wr)*state1*state1.dag() + 0.5*g*(state0*state1.dag() + state1*state0.dag())
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

def tls_analytical(wa, wr, g, c_o_o, c_1_o, t):
    delta = wa-wr
    omega = np.sqrt(delta**2.0 + g**2.0)
    c_o = np.exp(-1j*delta*t)*(c_o_o*np.cos(0.5*omega*t)+1j/omega*(delta*c_o_o-g*c_1_o)*np.sin(0.5*omega*t))
    c_1 = np.exp(-1j * delta * t) * (
    c_1_o * np.cos(0.5 * omega * t) - 1j / omega * (delta * c_1_o + g * c_o_o) * np.sin(0.5 * omega * t))
    return abs(c_o)**2, abs(c_1)**2

def dress_state(wa, wr, g):
    delta = wa - wr
    omega = np.sqrt(g**2 + delta**2)
    E0 = 0.5*(delta - omega)
    E1 = 0.5 * (delta + omega)
    return E0, E1
wa = 10
wr = 8
g = 5
gamma1 = 0.5
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


pop_a_0, pop_a_1 = tls_analytical(wr,wa,g,1,0,time)
pop_0, pop_1 = tls_Rabi(wr, wa, psi0, g, gamma1, time, "me")
plt.plot(time, pop_0, linewidth = '2', linestyle = '-', color = 'blue')
# plt.plot(time, pop_a_1, linewidth = '2', linestyle = '--', color = 'red')
plt.ylabel('Probability')
plt.xlabel('Time')
plt.tick_params(labelsize=18)


# wr_array = np.linspace(1,20,100)
# E0,E1 = dress_state(wa, wr_array, 2)
# plt.plot(wa - wr_array,E0,linewidth = '3', linestyle = '-',color = 'blue')
# plt.plot(wa - wr_array,E1,linewidth = '3', linestyle = '-',color = 'red')
# E0,E1 = dress_state(wa, wr_array, 0)
# plt.plot(wa - wr_array,E0,linewidth = '1', linestyle = '--',color = 'black')
# plt.plot(wa - wr_array,E1,linewidth = '1', linestyle = '--',color = 'black')
# plt.tick_params(
#     axis='y',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom='off',      # ticks along the bottom edge are off
#     top='off',         # ticks along the top edge are off
#     labelbottom='off') # labels along the bottom edge are off
# plt.axis('off')
# plt.tick_params(
#     axis='y',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom='off',      # ticks along the bottom edge are off
#     top='off',         # ticks along the top edge are off
#     labelbottom='off') # labels along the bottom edge are off

plt.show()

# print sigmay()