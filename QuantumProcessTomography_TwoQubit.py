import numpy as np
from matplotlib import pyplot as plt
from qutip import*

### Single qubit quantum process tomography ###
sI = np.array([[1,0],[0,1]])
sX = np.array([[0,1],[1,0]])
sY = np.array([[0,-1j],[1j,0]])
sZ = np.array([[1,0],[0,-1]])

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

#Define the gate
phi = np.pi/2
sX2 = np.cos(phi/2)*sI - 1j*np.sin(phi/2)*sX
sY2 = np.cos(phi/2)*sI - 1j*np.sin(phi/2)*sY
sZ2 = np.cos(phi/2)*sI - 1j*np.sin(phi/2)*sZ
CZ = np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,-1]])
CNOT = np.array([[1,0,0,0], [0,1,0,0], [0,0,0,1], [0,0,1,0]])
gate = XX

#Input state
rho0 = np.array([[1,0],[0,0]])
rho1 = np.array([[0,0],[0,1]])
rho2 = np.array([[0.5,0.5],[0.5,0.5]])
rho3 = np.array([[0.5,-0.5j],[0.5j,0.5]])
rho_input = np.array([rho0, rho1, rho2, rho3])

#Process, output experimentally determined by quantum state tomography
rho_output = np.full_like(rho_input, 0)
for idx in range (len(rho_input)):
    rho_output[idx] = gate.dot(rho_input[idx]).dot(np.conj(gate.transpose()))

#Decoherence process

#Quantum state tomography of output state.
la = np.zeros(16, dtype = complex)
chi = np.zeros(16, dtype = complex)
beta = np.zeros((16,16), dtype = complex)
for j in range(4):
    for k in range(4):
        la[k+4*j] = np.trace(rho_output[j].dot(rho_input[k]))
        for m in range(4):
            for n in range(4):
                 beta[k+4*j, n+4*m] = np.trace(E[m].dot(rho_input[j]).dot(np.conj(E[n].transpose())).dot(rho_input[k]))

kappa = np.linalg.inv(beta)
chi = kappa.dot(la)
chi = np.reshape(chi, (4,4))
# for m in range(4):
#     for n in range(4):
#         for j in range(4):
#             for k in range(4):
#                 chi[n+4*m] = chi[n+4*m]+ kappa[m,n,j,k]*la[j,k]

op_label = [['$I$','$X$','$Y$','$Z$'] for i in range (1)]
chi = Qobj(chi)
qpt_plot_combined(chi, op_label)

plt.show()