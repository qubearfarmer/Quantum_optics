import matplotlib.pyplot as plt
import numpy as np
from qutip import *

def JC_hamiltonian_Ho(N, wr, wa):

    a = tensor(destroy(N), qeye(2))
    smZ = tensor(qeye(N), sigmaz())
    nr = a.dag() * a
    Hr = wr*(a.dag()*a + 1.0/2)
    Ha = wa*smZ/2.0
    H = Hr + Ha
    return H

def JC_hamiltonian_Hc(N, g):
    a = tensor(destroy(N), qeye(2))  # resonator operator
    sm = tensor(qeye(N), sigmam())  # atom operator
    H = -0.5 * g * (a.dag() + a) * (sm.dag() + sm)
    return H

dummy_x = np.linspace(3,8,6)
dummy_y = np.zeros(len(dummy_x))

#No coupling, resonance
N = 10
wr = 10
wa = 10
g = 1
level_num = 7
plt.figure(1)
Ho = JC_hamiltonian_Ho(N, wr, wa)
Hc = JC_hamiltonian_Hc(N, g)
plot_energy_levels([Ho,Hc], 7)

#No coupling, resonance
N = 10
wr = 10
wa = 11
g = 0
level_num = 7
plt.figure(2)
Ho = JC_hamiltonian_Ho(N, wr, wa)
Hc = JC_hamiltonian_Hc(N, g)
plot_energy_levels([Ho,Hc], 7)

#
# #Coupled, off resonance
# N = 10
# wr = 10
# g = 1
# level_num = 3
# fig4 = plt.figure(4)
# dtune = np.linspace(-3,3,101) #delta / g
# dummy_y = np.zeros(len(dtune))
# for idx in range(level_num):
#     for idy, dw  in enumerate(dtune):
#         wa = wr - dw*g
#         dummy_y[idy] = JC_hamiltonian(N, wr, wa, g).eigenenergies()[idx]
#     plt.plot(dtune, dummy_y, linewidth = '2', color = 'blue')
#     #No coupling case
#     for idy, dw  in enumerate(dtune):
#         wa = wr - dw*g
#         dummy_y[idy] = JC_hamiltonian(N, wr, wa, 0).eigenenergies()[idx]
#     plt.plot(dtune, dummy_y, linestyle = '--', color = 'green')
# plt.xlabel('Coupled, detuned')
# plt.xlim([-3.5,3.5])
# plt.ylim([8,12])


plt.show()

