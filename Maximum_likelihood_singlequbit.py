from qutip import*
import numpy as np
from scipy.optimize import minimize
state = 0.5*np.array([[1,1], [1,1]])

def rho_single_qubit(t0,t1,t2):
    t = np.array([[t0, t1+1j*t2], [0, 1-t0]])
    rho = t.conjugate().transpose().dot(t) / np.trace(t.conjugate().transpose().dot(t))
    return rho

sI = np.array([[1,0],[0,1]])
sX = np.array([[0,1],[1,0]])
sY = np.array([[0,-1j],[1j,0]])
sZ = np.array([[1,0],[0,-1]])

M = np.zeros((3,2,2), dtype = complex)
#Define measurement operator
M[0,:,:] = sI + sZ
M[1,:,:] = sI + sY
M[2,:,:] = sI + sX
#Now we measure the state
m = np.zeros(3, dtype = complex) #measurement
epsilon = 0.1
rho = state*(1-epsilon) + epsilon*np.random.random((2, 2))
for idx in range (3):
    m[idx] = np.trace(M[idx,:,:].dot(rho))

def likelihood(x):
    dist = 0
    for gate in range(3):
        dist = dist + abs(m[gate] - np.trace(M[gate, :,:].dot(rho_single_qubit(x[0], x[1], x[2]))))**2
    return dist
guess = np.ones(3)
res = minimize(likelihood, guess, method='powell',tol=1.e-10,
            options={'maxiter': 10000})
rho = rho_single_qubit(*guess)
t = res.x
rho_reconstructed = Qobj(rho_single_qubit(*t))
matrix_histogram_complex(rho_reconstructed)