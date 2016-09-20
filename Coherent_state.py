import matplotlib.pyplot as plt
import numpy as np
from qutip import *
N = 15

#coherent field
rho_coherent = coherent(N, 1.5)
plot_wigner_fock_distribution(rho_coherent)

#thermal state
rho_thermal = thermal_dm(N, 2)
plot_wigner_fock_distribution(rho_thermal)
plt.show()