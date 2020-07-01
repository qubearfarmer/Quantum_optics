# Project: Quantum tomography
# author: Long Nguyen
# Date: 12-1-2019
import numpy as np
from qutip import*
from matplotlib import pyplot as plt
from scipy.optimize import minimize

# x = np.array([[1,0,0,1], [0,0,0,0], [0,0,0,0], [1,0,0,1]])
def density_matrix(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15):
    t16 = 1 - t1 -t2 - t3
    tau = np.array([[t1, t4+1j*t5, t6+1j*t7, t8+1j*t9], [0, t2, t10+1j*t11, t12+1j*t13], [0, 0, t3, t14+1j*t15], [0,0,0,t16]])
    rho = np.conj(tau.transpose()).dot(tau)
    rho = rho/np.trace(rho)
    return rho

x = density_matrix(t1=0.5,t2=0,t3=0,t4=0,t5=0,t6=0,t7=0,t8=0.5,t9=0,t10=0,t11=0,t12=0,t13=0,t14=0,t15=0)
#random states
betaII = (np.random.randint(1,1000) + 1j*np.random.randint(1,1000))*1e-6
betaZI = (np.random.randint(1,1000) + 1j*np.random.randint(1,1000))*1e-6
betaIZ = (np.random.randint(1,1000) + 1j*np.random.randint(1,1000))*1e-6
betaZZ = (np.random.randint(1,1000) + 1j*np.random.randint(1,1000))*1e-6
#Labber
gate_sequence = np.array(['I-I','Xp-I','I-Xp',\
                          'X2p-I','X2p-X2p','X2p-Y2p','X2p-Xp',\
                          'Y2p-I','Y2p-X2p','Y2p-Y2p','Y2p-Xp',\
                          'I-X2p','Xp-X2p','I-Y2p','Xp-Y2p'])
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

M = np.zeros((len(gate_sequence), 4, 4), dtype = complex)
M[0,:,:] = betaII*II + betaZI*ZI + betaIZ*IZ + betaZZ*ZZ
M[1,:,:] = betaII*II - betaZI*ZI + betaIZ*IZ - betaZZ*ZZ
M[2,:,:] = betaII*II + betaZI*ZI - betaIZ*IZ - betaZZ*ZZ

M[3,:,:] = betaII*II + betaZI*YI + betaIZ*IZ + betaZZ*YZ
M[4,:,:] = betaII*II + betaZI*YI + betaIZ*IY + betaZZ*YY
M[5,:,:] = betaII*II + betaZI*YI - betaIZ*IX - betaZZ*YX
M[6,:,:] = betaII*II + betaZI*YI - betaIZ*IZ - betaZZ*YZ

M[7,:,:] = betaII*II - betaZI*XI + betaIZ*IZ - betaZZ*XZ
M[8,:,:] = betaII*II - betaZI*XI + betaIZ*IY - betaZZ*XY
M[9,:,:] = betaII*II - betaZI*XI - betaIZ*IX + betaZZ*XX
M[10,:,:] = betaII*II - betaZI*XI - betaIZ*IZ + betaZZ*XZ

M[11,:,:] = betaII*II + betaZI*ZI + betaIZ*IY + betaZZ*ZY
M[12,:,:] = betaII*II - betaZI*ZI + betaIZ*IY - betaZZ*ZY
M[13,:,:] = betaII*II + betaZI*ZI - betaIZ*IX - betaZZ*ZX
M[14,:,:] = betaII*II - betaZI*ZI - betaIZ*IX + betaZZ*ZX

#check
m = np.zeros(len(gate_sequence), dtype = complex) #measurement
epsilon = 0.1
rho = x*(1-epsilon) + epsilon*np.random.random((4, 4))
for idx in range (len(gate_sequence)):
    m[idx] = np.trace(M[idx,:,:].dot(rho))

matrix_histogram_complex(x)
matrix_histogram_complex(rho)

def likelihood(x):
    dist = 0
    for idx in range(len(gate_sequence)):
        dist = dist + abs((m[idx] - np.trace(M[idx, :, :].dot(density_matrix(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10], x[11],x[12],x[13],x[14])))))**2
    return dist
guess = np.ones(15)
# guess[0] = 1
res = minimize(likelihood, guess, method='powell',tol=1.e-10,
            options={'maxiter': 10000})
t = res.x
rho_reconstructed = Qobj(density_matrix(*t))
matrix_histogram_complex(rho_reconstructed)
plt.show()

