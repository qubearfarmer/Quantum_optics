# Project: Quantum tomography
# author: Long Nguyen
# Date: 12-1-2019
import numpy as np
from qutip import*
from matplotlib import pyplot as plt

#random states
rho1 = rand_dm(2)
rho2 = rand_dm(2)
rho = tensor(rho1, rho2)
#Labber
gate_sequence = np.array(['I-I','Xp-I','I-Xp',\
                          'X2p-I','X2p-X2p','X2p-Y2p','X2p-Xp',\
                          'Y2p-I','Y2p-X2p','Y2p-Y2p','Y2p-Xp',\
                          'I-X2p','Xp-X2p','I-Y2p','Xp-Y2p'])
betaII = (np.random.randint(1,1000) + 1j*np.random.randint(1,1000))*1e-6
betaZI = (np.random.randint(1,1000) + 1j*np.random.randint(1,1000))*1e-6
betaIZ = (np.random.randint(1,1000) + 1j*np.random.randint(1,1000))*1e-6
betaZZ = (np.random.randint(1,1000) + 1j*np.random.randint(1,1000))*1e-6
M0 = 0.25*(betaII*tensor(qeye(2), qeye(2)) + betaIZ*tensor(qeye(2), sigmaz()) + betaZI*tensor(sigmaz(), qeye(2)) + betaZZ*tensor(sigmaz(), sigmaz()))
M1 = 0.25*(betaII*tensor(qeye(2), qeye(2)) + betaIZ*tensor(qeye(2), sigmaz()) - betaZI*tensor(sigmaz(), qeye(2)) - betaZZ*tensor(sigmaz(), sigmaz()))
M2 = 0.25*(betaII*tensor(qeye(2), qeye(2)) - betaIZ*tensor(qeye(2), sigmaz()) + betaZI*tensor(sigmaz(), qeye(2)) - betaZZ*tensor(sigmaz(), sigmaz()))

M3 = 0.25*(betaII*tensor(qeye(2), qeye(2)) + betaIZ*tensor(qeye(2), sigmaz()) + betaZI*tensor(sigmay(), qeye(2)) + betaZZ*tensor(sigmay(), sigmaz()))
M4 = 0.25*(betaII*tensor(qeye(2), qeye(2)) + betaIZ*tensor(qeye(2), sigmay()) + betaZI*tensor(sigmay(), qeye(2)) + betaZZ*tensor(sigmay(), sigmay()))
M5 = 0.25*(betaII*tensor(qeye(2), qeye(2)) - betaIZ*tensor(qeye(2), sigmax()) + betaZI*tensor(sigmay(), qeye(2)) - betaZZ*tensor(sigmay(), sigmax()))
M6 = 0.25*(betaII*tensor(qeye(2), qeye(2)) - betaIZ*tensor(qeye(2), sigmaz()) + betaZI*tensor(sigmay(), qeye(2)) - betaZZ*tensor(sigmay(), sigmaz()))

M7 = 0.25*(betaII*tensor(qeye(2), qeye(2)) + betaIZ*tensor(qeye(2), sigmaz()) - betaZI*tensor(sigmax(), qeye(2)) - betaZZ*tensor(sigmax(), sigmaz()))
M8 = 0.25*(betaII*tensor(qeye(2), qeye(2)) + betaIZ*tensor(qeye(2), sigmay()) - betaZI*tensor(sigmax(), qeye(2)) - betaZZ*tensor(sigmax(), sigmay()))
M9 = 0.25*(betaII*tensor(qeye(2), qeye(2)) - betaIZ*tensor(qeye(2), sigmax()) - betaZI*tensor(sigmax(), qeye(2)) + betaZZ*tensor(sigmax(), sigmax()))
M10 = 0.25*(betaII*tensor(qeye(2), qeye(2)) - betaIZ*tensor(qeye(2), sigmaz()) - betaZI*tensor(sigmax(), qeye(2)) + betaZZ*tensor(sigmax(), sigmaz()))

M11 = 0.25*(betaII*tensor(qeye(2), qeye(2)) + betaIZ*tensor(qeye(2), sigmay()) + betaZI*tensor(sigmaz(), qeye(2)) + betaZZ*tensor(sigmaz(), sigmay()))
M12 = 0.25*(betaII*tensor(qeye(2), qeye(2)) + betaIZ*tensor(qeye(2), sigmay()) - betaZI*tensor(sigmaz(), qeye(2)) - betaZZ*tensor(sigmaz(), sigmay()))
M13 = 0.25*(betaII*tensor(qeye(2), qeye(2)) - betaIZ*tensor(qeye(2), sigmax()) + betaZI*tensor(sigmaz(), qeye(2)) - betaZZ*tensor(sigmaz(), sigmax()))
M14 = 0.25*(betaII*tensor(qeye(2), qeye(2)) - betaIZ*tensor(qeye(2), sigmax()) - betaZI*tensor(sigmaz(), qeye(2)) + betaZZ*tensor(sigmaz(), sigmax()))

m0 = expect(M0,rho)
m1 = expect(M1,rho)
m2 = expect(M2,rho)
m3 = expect(M3,rho)
m4 = expect(M4,rho)
m5 = expect(M5,rho)
m6 = expect(M6,rho)
m7 = expect(M7,rho)
m8 = expect(M8,rho)
m9 = expect(M9,rho)
m10 = expect(M10,rho)
m11 = expect(M11,rho)
m12 = expect(M12,rho)
m13 = expect(M13,rho)
m14 = expect(M14,rho)


measurement_matrix = np.empty((len(gate_sequence), 15), dtype=complex)
gate_matrix = np.zeros((4, 4), dtype = complex)
gate_matrix[0,0] = betaII
for idx, gate in enumerate(gate_sequence):
    gate_matrix = np.zeros((4, 4), dtype=complex)
    gate1,gate2 = gate.split('-')
    if gate1 == 'I':
        row_index = 3
        sign1 = 1
    elif gate1[0:2] == 'Y2':
        row_index = 1
        if gate1[-1] == 'p':
            sign1 = -1
        elif gate1[-1] == 'm':
            sign1 = 1
    elif gate1[0:2] == 'X2':
        row_index = 2
        if gate1[-1] == 'm':
            sign1 = -1
        elif gate1[-1] == 'p':
            sign1 = 1
    elif gate1 in ['Xp', 'Xm', 'Yp', 'Ym']:
        row_index = 3
        sign1 = -1

    if gate2 == 'I':
        column_index = 3
        sign2 = 1
    elif gate2[0:2] == 'Y2':
        column_index = 1
        if gate2[-1] == 'p':
            sign2 = -1
        elif gate2[-1] == 'm':
            sign2 = 1
    elif gate2[0:2] == 'X2':
        column_index = 2
        if gate2[-1] == 'm':
            sign2 = -1
        elif gate2[-1] == 'p':
            sign2 = 1
    elif gate2 in ['Xp', 'Xm', 'Yp', 'Ym']:
        column_index = 3
        sign2 = -1

    gate_matrix[0, column_index] = sign2*betaIZ
    gate_matrix[row_index, 0] = sign1*betaZI
    gate_matrix[row_index, column_index] = sign1 * sign2 * betaZZ
    gate_matrix[0, 0] = betaII
    measurement_matrix[idx, 0:3] = 0.25*gate_matrix[0, 1:]
    measurement_matrix[idx, 3:7] = 0.25*gate_matrix[1, :]
    measurement_matrix[idx, 7:11] = 0.25*gate_matrix[2, :]
    measurement_matrix[idx, 11:15] = 0.25*gate_matrix[3, :]


    # print (gate1, gate2, row_index, column_index, sign1, sign2)
    # print (measurement_matrix[idx,:] - measurement_matrix_test[idx,:])

avgIX, avgIY, avgIZ, \
avgXI, avgXX, avgXY, avgXZ, \
avgYI, avgYX, avgYY, avgYZ, \
avgZI, avgZX, avgZY, avgZZ \
    = np.linalg.inv(measurement_matrix).dot(np.array([  
    m0, m1, m2, m3,
    m4, m5, m6, m7,
    m8, m9, m10, m11,
    m12, m13, m14]).transpose()-0.25*betaII).transpose()


rho_reconstructed = 0.25*(tensor(qeye(2), qeye(2)) + avgIX*tensor(qeye(2), sigmax()) + avgIY*tensor(qeye(2), sigmay()) + avgIZ*tensor(qeye(2), sigmaz()) +
                          avgXI*tensor(sigmax(), qeye(2)) + avgXX*tensor(sigmax(), sigmax()) + avgXY*tensor(sigmax(), sigmay()) + avgXZ*tensor(sigmax(), sigmaz()) +
                          avgYI*tensor(sigmay(), qeye(2)) + avgYX*tensor(sigmay(), sigmax()) + avgYY*tensor(sigmay(), sigmay()) + avgYZ*tensor(sigmay(), sigmaz()) +
                          avgZI*tensor(sigmaz(), qeye(2)) + avgZX*tensor(sigmaz(), sigmax()) + avgZY*tensor(sigmaz(), sigmay()) + avgZZ*tensor(sigmaz(), sigmaz()))

#plotting
# rho_target = ket2dm(bell_state(state='10'))
matrix_histogram_complex(rho)
# matrix_histogram_complex(rho_target)
# matrix_histogram_complex(rho_measure)
matrix_histogram_complex(rho_reconstructed)
# print (fidelity(rho_target, rho))
# print (fidelity(rho_target, rho_reconstructed))
# print ((rho_reconstructed.dag()*rho_reconstructed).tr())
print (expect(sigmaz(),basis(2,1)))
plt.show()
