import numpy as np
import scipy.linalg

rho0 = np.array([[1,0],[0,0]])
rho1 = np.array([[0,0],[0,1]])

def fidelity(rho,rho_ideal):
    f = abs(np.trace(scipy.linalg.sqrtm(scipy.linalg.sqrtm(rho_ideal).dot(rho).dot(scipy.linalg.sqrtm(rho_ideal)))))
    return f

print (fidelity(rho0,rho0))