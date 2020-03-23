import matplotlib.pyplot as plt
import numpy as np
from qutip import *
N = 100
xvec = np.linspace(-5,5,200)
#coherent field
rho_coherent = coherent(10, 1)
# plot_fock_distribution(rho_coherent)

# plt.xlim([-1,8])
# plt.ylim([0,0.5])

# plot_fock_distribution(rho_coherent)
# plt.tick_params(labelsize = 14.0)
# plt.ylim([0,0.2])
# plt.xlim([0,21])

plt.figure(figsize = [6,6])
W = wigner(rho_coherent, xvec, xvec)
cmap = wigner_cmap(W)
X, Y = np.meshgrid(xvec, xvec)
plt.contourf(X, Y, W, 50, cmap=cmap)
plt.tick_params(labelsize = 15.0)
path = 'C:\\Users\\nguyen89\Google Drive\Research\Illustration\Thesis\Chapter 2\coherent3.pdf'
plt.savefig(path, dpi=300)


# Q = qfunc(rho_coherent,xvec,xvec)
# qplt = plt.contourf(xvec, xvec, np.real(Q), 100, cmap = 'GnBu')
# plt.tick_params(labelsize = 15.0)
# path = 'C:\\Users\\nguyen89\Google Drive\Research\Illustration\Thesis\Chapter 2\coherent3.pdf'
# plt.savefig(path, dpi=300)

#thermal state
# rho_thermal = thermal_dm(N, 2)
# plot_wigner_fock_distribution(rho_thermal)
plt.show()