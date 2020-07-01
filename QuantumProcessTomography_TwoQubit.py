import numpy as np
from matplotlib import pyplot as plt
from qutip import*

### Single qubit quantum process tomography ###
sI = np.array([[1,0],[0,1]])
sX = np.array([[0,1],[1,0]])
sY = np.array([[0,-1j],[1j,0]])
sZ = np.array([[1,0],[0,-1]])

sX2p = np.cos(np.pi/4)*sI - 1j*np.sin(np.pi/4)*sX
sY2p = np.cos(np.pi/4)*sI - 1j*np.sin(np.pi/4)*sY
sX2m = np.cos(-np.pi/4)*sI - 1j*np.sin(-np.pi/4)*sX
sY2m = np.cos(-np.pi/4)*sI - 1j*np.sin(-np.pi/4)*sY

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
E = np.array([II,IX,IY,IZ,XI,XX,XY,XZ,YI,YX,YY,YZ,ZI,ZX,ZY,ZZ])

IY2p = np.kron(sI,sY2p)
IX2p = np.kron(sI,sX2p)
IX2m = np.kron(sI,sX2m)
XY2p = np.kron(sX,sY2p)
XX2m = np.kron(sX,sX2m)
Y2pI = np.kron(sY2p,sI)
Y2pX = np.kron(sY2p,sX)
Y2pY2p = np.kron(sY2p,sY2p)
Y2pX2m = np.kron(sY2p,sX2m)
X2mI = np.kron(sX2m,sI)
X2mX = np.kron(sX2m,sX)
X2mY2p = np.kron(sX2m,sY2p)
X2mX2m = np.kron(sX2m,sX2m)
#Define the gate
phi = np.pi/2
CZ = np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,-1]])
CNOT = np.array([[1,0,0,0], [0,1,0,0], [0,0,0,1], [0,0,1,0]])
ZZp = np.array([[1,0,0,0], [0,-1,0,0], [0,0,-1,0], [0,0,0,1]])
ZZ2p = np.array([[1,0,0,0], [0,1j,0,0], [0,0,1j,0], [0,0,0,1]])
gate = np.kron(sX2p,sX2p)
gate = ZZ2p

#Input state
rho0 = np.array([[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
rho1 = IX.dot(rho0).dot(np.conj(IX.transpose()))
rho2 = IY2p.dot(rho0).dot(np.conj(IY2p.transpose()))
rho3 = IX2m.dot(rho0).dot(np.conj(IX2m.transpose()))
rho4 = XI.dot(rho0).dot(np.conj(XI.transpose()))
rho5 = XX.dot(rho0).dot(np.conj(XX.transpose()))
rho6 = XY2p.dot(rho0).dot(np.conj(XY2p.transpose()))
rho7 = XX2m.dot(rho0).dot(np.conj(XX2m.transpose()))
rho8 = Y2pI.dot(rho0).dot(np.conj(Y2pI.transpose()))
rho9 = Y2pX.dot(rho0).dot(np.conj(Y2pX.transpose()))
rho10 = Y2pY2p.dot(rho0).dot(np.conj(Y2pY2p.transpose()))
rho11 = Y2pX2m.dot(rho0).dot(np.conj(Y2pX2m.transpose()))
rho12 = X2mI.dot(rho0).dot(np.conj(X2mI.transpose()))
rho13 = X2mX.dot(rho0).dot(np.conj(X2mX.transpose()))
rho14 = X2mY2p.dot(rho0).dot(np.conj(X2mY2p.transpose()))
rho15 = X2mX2m.dot(rho0).dot(np.conj(X2mX2m.transpose()))
rho_input = np.array([rho0, rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8, rho9, rho10, rho11, rho12, rho13, rho14, rho15])

#Process, output experimentally determined by quantum state tomography
rho_output = np.full_like(rho_input, 0)
for idx in range (len(rho_input)):
    rho_output[idx] = gate.dot(rho_input[idx]).dot(np.conj(gate.transpose()))

#Decoherence process

#Quantum state tomography of output state.
n = 2
d = 2**n
la = np.zeros(d**4, dtype = complex)
chi = np.zeros(d**4, dtype = complex)
beta = np.zeros((d**4,d**4), dtype = complex)
for j in range(16):
    for k in range(16):
        la[k+16*j] = np.trace(rho_output[j].dot(rho_input[k]))
        for m in range(16):
            for n in range(16):
                 beta[k+16*j, n+16*m] = np.trace(E[m].dot(rho_input[j]).dot(np.conj(E[n].transpose())).dot(rho_input[k]))

kappa = np.linalg.inv(beta)
chi_ideal = kappa.dot(la)
chi_ideal = np.reshape(chi_ideal, (16,16)).transpose()

op_label = [['$I$','$X$','$Y$','$Z$'] for i in range (2)]
chi_ideal_plot = Qobj(chi_ideal)
qpt_plot_combined(chi_ideal_plot, op_label)
plt.show()