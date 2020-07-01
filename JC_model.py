import matplotlib.pyplot as plt
import numpy as np
from qutip import *

def JC_hamiltonian(N, wr, wa,g):

    a = tensor(destroy(N), qeye(2))
    smZ = tensor(qeye(N), sigmaz())
    sm = tensor(qeye(N), sigmam())
    nr = a.dag() * a
    Hr = wr*(a.dag()*a + 1.0/2)
    Ha = wa*smZ/2.0
    H = Hr + Ha -0.5 * g * (a.dag() + a) * (sm.dag() + sm)
    return H

def JC_hamiltonian_Hc(N, g):
    a = tensor(destroy(N), qeye(2))  # resonator operator
    sm = tensor(qeye(N), sigmam())  # atom operator
    H = -0.5 * g * (a.dag() + a) * (sm.dag() + sm)
    return H


#No coupling, resonance
# N = 10
# wr = 10
# wa = 10
# g = 1
# level_num = 7
# plt.figure(1)
# Ho = JC_hamiltonian_Ho(N, wr, wa)
# Hc = JC_hamiltonian_Hc(N, g)
# plot_energy_levels([Ho,Hc], 7)

#No coupling, resonance
# N = 10
# wr = 10
# wa = 11
# g = 0
# level_num = 7
# plt.figure(2)
# Ho = JC_hamiltonian_Ho(N, wr, wa)
# Hc = JC_hamiltonian_Hc(N, g)
# plot_energy_levels([Ho,Hc], 7)

#
# #Coupled, off resonance
N = 10
wr = 10
g = 1
level_num = 3
fig4 = plt.figure(4)
dtune = np.linspace(-3,3,101) #delta / g
dummy_y = np.zeros(len(dtune))
for idx in range(level_num):
    for idy, dw  in enumerate(dtune):
        wa = wr - dw*g
        dummy_y[idy] = JC_hamiltonian(N, wr, wa, g).eigenenergies()[idx]
    plt.plot(dtune/g, (dummy_y-wr)/wr, linewidth = '3')
    #No coupling case
    for idy, dw  in enumerate(dtune):
        wa = wr - dw*g
        dummy_y[idy] = JC_hamiltonian(N, wr, wa, 0).eigenenergies()[idx]
    plt.plot(dtune/g, (dummy_y-wr)/wr, linestyle = '--')
# plt.xlabel('Coupled, detuned')
plt.xlim([-3.5,3.5])
plt.ylim([-0.2,0.2])
plt.yticks([-0.2,-0.1,0,0.1,0.2])
plt.tick_params(labelsize = 15.0)
path = 'C:\\Users\\nguyen89\Google Drive\Research\Illustration\Thesis\Chapter 2\JC_resonant.pdf'
plt.savefig(path, dpi=300)
# plt.yticks([])


plt.show()

