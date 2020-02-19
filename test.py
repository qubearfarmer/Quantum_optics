import matplotlib.pyplot as plt
import numpy as np
from qutip import*
IX2 = tensor(qeye(2),rx(phi=np.pi/2))
X2I = tensor(rx(phi=np.pi/2),qeye(2))
IY2 = tensor(qeye(2),ry(phi=np.pi/2))
Y2I = tensor(ry(phi=np.pi/2),qeye(2))

IX = tensor(qeye(2),rx(phi=np.pi))
XI = tensor(rx(phi=np.pi),qeye(2))
IY = tensor(qeye(2),ry(phi=np.pi))
IZ = tensor(qeye(2),sigmaz())
ZI = tensor(sigmaz(),qeye(2))
gnd = tensor(basis(2,0), basis(2,0))

XX = tensor(sigmax(),sigmax())
XZ = tensor(sigmax(),sigmaz())
ZX = tensor(sigmaz(),sigmax())
ZZ = tensor(sigmaz(),sigmaz())
XY = tensor(sigmax(),sigmay())
YX = tensor(sigmay(),sigmax())
YY = tensor(sigmay(),sigmay())

wA = 2*np.pi*72.4e6
wB = 2*np.pi*136.3e6
t=2400e-9
initial_state = tensor(basis(2,1), basis(2,1))
zA = tensor(rz(phi=wA*t),qeye(2))
zB = tensor(qeye(2), rz(phi=(wB-wA)*t))

t2 = (2*np.pi - (wB-wA)*t)/(wB-wA)
zB2 = tensor(qeye(2), rz(phi=(wB-wA)*t))
zA2 = tensor(rz(phi=(wB-wA)*t2), qeye(2))
rho = ket2dm(zA2*X2I*csign()*IX2*X2I*initial_state)
rho = ket2dm(IX2*csign()*IX2*X2I*initial_state)
matrix_histogram_complex(rho)
# plt.title(t)