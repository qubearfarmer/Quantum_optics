import numpy as np
from matplotlib import pyplot as plt
from qutip import*

### Single qubit quantum process tomography ###
sI = np.array([[1,0],[0,1]])
sX = np.array([[0,1],[1,0]])
sY = np.array([[0,-1j],[1j,0]])
sZ = np.array([[1,0],[0,-1]])
E = np.array([sI,sX,sY,sZ])

#Define the gate
phi = np.pi/2
sX2 = np.cos(phi/2)*sI - 1j*np.sin(phi/2)*sX
sY2 = np.cos(phi/2)*sI - 1j*np.sin(phi/2)*sY
sZ2 = np.cos(phi/2)*sI - 1j*np.sin(phi/2)*sZ
sX2m = np.cos(phi/2)*sI - 1j*np.sin(-phi/2)*sX
gate = sZ2
path = 'C:\\Users\\nguyen89\Google Drive\Research\Illustration\Thesis\Chapter 2\\processtomo_t2.pdf'

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
#Relaxation
gamma1 = 1
gamma2 = 1
t = 5
rho_output[0] = np.array([[np.exp(-gamma1*t), 0], [0, 1-np.exp(-gamma1*t)]])
rho_output[1] = np.array([[0,0],[0,1]])
rho_output[2] = np.array([[0.5*np.exp(-gamma1*t),0.5*np.exp(-gamma1*t/2)],[0.5*np.exp(-gamma1*t/2), 1-0.5*np.exp(-gamma1*t)]])
rho_output[3] = np.array([[0.5*np.exp(-gamma1*t),-0.5j*np.exp(-gamma1*t/2)],[0.5j*np.exp(-gamma1*t/2),1-0.5*np.exp(-gamma1*t)]])
#Dephasing
# rho_output[0] = np.array([[1,0],[0,0]])
# rho_output[1] = np.array([[0,0],[0,1]])
# rho_output[2] = np.array([[0.5,0.5*np.exp(-gamma2*t)],[0.5*np.exp(-gamma2*t),0.5]])
# rho_output[3] = np.array([[0.5,-0.5j*np.exp(-gamma2*t)],[0.5j*np.exp(-gamma2*t),0.5]])

#Quantum state tomography of output state.
n = 1
d =2**n
la = np.zeros(d**4, dtype = complex)
chi = np.zeros(d**4, dtype = complex)
beta = np.zeros((d**4,d**4), dtype = complex)
for j in range(4):
    for k in range(4):
        la[k+4*j] = np.trace(rho_output[j].dot(rho_input[k]))
        for m in range(4):
            for n in range(4):
                 beta[k+4*j, n+4*m] = np.trace(E[m].dot(rho_input[j]).dot(np.conj(E[n].transpose())).dot(rho_input[k]))

kappa = np.linalg.inv(beta)
chi = kappa.dot(la)
chi = np.reshape(chi, (4,4)).transpose()
# for m in range(4):
#     for n in range(4):
#         for j in range(4):
#             for k in range(4):
#                 chi[n+4*m] = chi[n+4*m]+ kappa[m,n,j,k]*la[j,k]

op_label = [["","$I$", "$X$", "$Y$", "$Z$"] for i in range (1)]

chi = Qobj(chi)
qpt_plot_combined(chi, op_label)
# plt.savefig(path, dpi=300)
plt.show()