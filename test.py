import matplotlib.pyplot as plt
import numpy as np
from qutip import*
IX2 = tensor(qeye(2),rx(phi=np.pi/2))
X2I = tensor(rx(phi=np.pi/2),qeye(2))
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

# rho = ket2dm(tensor(rx(phi=np.pi/2)*basis(2,0), rx(phi=np.pi/2)*basis(2,0)))
# matrix_histogram_complex(rho)
# plt.title('x/2, x/2')
# rho = ket2dm(cphase(theta=np.pi*0.7)*tensor(rx(phi=np.pi/2)*basis(2,0), rx(phi=np.pi/2)*basis(2,0)))
# matrix_histogram_complex(rho)
# plt.title('after CZ')
rho = ket2dm(IX2*cphase(theta=np.pi)*tensor(rx(phi=np.pi/2)*basis(2,0), rx(phi=np.pi/2)*basis(2,0)))
# matrix_histogram_complex(rho)
# plt.title('Bell')
# rho = ket2dm(IX2*IX*cphase(theta=np.pi)*tensor(rx(phi=np.pi/2)*basis(2,0), rx(phi=np.pi/2)*basis(2,0)))
# matrix_histogram_complex(rho)
# plt.title('Bell')
# print (csign()*tensor(qeye(2),rx(phi=np.pi*2)))
# plt.show()


print (expect(XX,bell_state('10')) - expect(XY,bell_state('10')) + expect(YX,bell_state('10')) + expect(YY,bell_state('10')))