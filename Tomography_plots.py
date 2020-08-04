from qutip import*
import numpy as np
from matplotlib import pyplot as plt
state = (basis(2,0) + 1j*basis(2,1))/np.sqrt(2)
matrix_histogram_complex(ket2dm(state))
plt.show()


sI = np.array([[1,0],[0,1]])
sX = np.array([[0,1],[1,0]])
sY = np.array([[0,-1j],[1j,0]])
sZ = np.array([[1,0],[0,-1]])

theta = np.pi/2
sX2 = np.cos(theta/2)*sI - 1j*np.sin(theta/2)*sX
print(sX2)
CZ = np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,-1]])

II = np.kron(sI,sI)
IX = np.kron(sI,sX)
IY = np.kron(sI,sY)
IZ = np.kron(sI,sZ)
XI = np.kron(sX,sI)
XX = np.kron(sX,sX)
XY = np.kron(sX,sY)
XZ = np.kron(sX,sZ)
YI = np.kron(sY,sI)
YX = np.kron(sY,sX)
YY = np.kron(sY,sY)
YZ = np.kron(sY,sZ)
ZI = np.kron(sZ,sI)
ZX = np.kron(sZ,sX)
ZY = np.kron(sZ,sY)
ZZ = np.kron(sZ,sZ)
X2X2 = np.kron(sX2,sX2)
X2I = np.kron(sX2,sI)


iState = np.array([[1,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]])
state = X2X2.dot(iState).dot(X2X2.conjugate().transpose())
state = CZ.dot(state).dot(CZ.conjugate().transpose())
state = X2I.dot(state).dot(X2I.conjugate().transpose())
rho = Qobj(state)
matrix_histogram_complex(rho)

# state = tensor(basis(2,0), basis(2,0))
# state = tensor(sigmax(), qeye(2))*state
# rho = ket2dm(state)
# matrix_histogram_complex(rho)



