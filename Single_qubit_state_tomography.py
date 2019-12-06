# Project: Quantum tomography
# author: Long Nguyen
# Date: 12-1-2019
import numpy as np
from qutip import *
from matplotlib import pyplot as plt

# random states
rho1 = rand_dm(2)
rho2 = rand_dm(2)
# ground state
# rho1 = ket2dm(basis(2, 0))
# rho2 = ket2dm(basis(2, 0))
# superposition state
# rho1 = rx(np.pi / 2) * rho1 * (rx(np.pi / 2)).dag()
# rho2 = rx(np.pi/2)*rho2*(rx(np.pi/2)).dag()

##############single qubit tomography##############
#ideal measurement
rho1_measure = 0.5*(qeye(2)+expect(sigmax(),rho1)*sigmax()+expect(sigmay(),rho1)*sigmay()+expect(sigmaz(),rho1)*sigmaz())
#real measurement with voltages
#beta coeffcients are experimentally determined from either single shot readout or from Rabi calibration + temperature measurement
betaI = (np.random.randint(1,1000) + 1j*np.random.randint(1,1000))*1e-6
betaZ = (np.random.randint(1,1000) + 1j*np.random.randint(1,1000))*1e-6
#Gate sequence in Labber is I, X2p, Y2m
mI = 0.5*expect(betaI*qeye(2)-betaZ*sigmaz(), rho1)
mX2p = 0.5*expect(betaI*qeye(2)-betaZ*sigmay(), rho1) #measurement after doing X pi/2 pulse
mY2p = 0.5*expect(betaI*qeye(2)+betaZ*sigmax(), rho1) #measurement after doing Y pi/2 pulse

#smarter math
measurement_matrix = 0.5*np.array([[0, 0, -betaZ], [0, -betaZ, 0], [betaZ, 0, 0]])
# print (measurement_matrix)
avgX, avgY, avgZ = np.linalg.inv(measurement_matrix).dot(np.array([mI, mX2p, mY2p]).transpose()-0.5*betaI).transpose()
rho1_reconstructed = 0.5*(qeye(2) + avgX*sigmax() + avgY*sigmay() + avgZ*sigmaz())
matrix_histogram_complex(rho1)
matrix_histogram_complex(rho1_reconstructed)
plt.show()