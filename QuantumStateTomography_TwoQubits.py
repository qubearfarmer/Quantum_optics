# Project: Quantum tomography
# author: Long Nguyen
# Date: 12-1-2019
import numpy as np
from qutip import*
from matplotlib import pyplot as plt

x = np.array([[1,0,0,0], [0,0,0,0], [1,0,0,0], [1,0,0,0]])
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
M[0,:,:] = betaII*II - betaZI*ZI - betaIZ*IZ + betaZZ*ZZ
M[1,:,:] = betaII*II + betaZI*ZI - betaIZ*IZ - betaZZ*ZZ
M[2,:,:] = betaII*II - betaZI*ZI + betaIZ*IZ - betaZZ*ZZ

M[3,:,:] = betaII*II - betaZI*YI - betaIZ*IZ + betaZZ*YZ
M[4,:,:] = betaII*II - betaZI*YI - betaIZ*IY + betaZZ*YY
M[5,:,:] = betaII*II - betaZI*YI + betaIZ*IX - betaZZ*YX
M[6,:,:] = betaII*II - betaZI*YI + betaIZ*IZ - betaZZ*YZ

M[7,:,:] = betaII*II + betaZI*XI - betaIZ*IZ - betaZZ*XZ
M[8,:,:] = betaII*II + betaZI*XI - betaIZ*IY - betaZZ*XY
M[9,:,:] = betaII*II + betaZI*XI + betaIZ*IX + betaZZ*XX
M[10,:,:] = betaII*II + betaZI*XI + betaIZ*IZ + betaZZ*XZ

M[11,:,:] = betaII*II - betaZI*ZI - betaIZ*IY + betaZZ*ZY
M[12,:,:] = betaII*II + betaZI*ZI - betaIZ*IY - betaZZ*ZY
M[13,:,:] = betaII*II - betaZI*ZI + betaIZ*IX - betaZZ*ZX
M[14,:,:] = betaII*II + betaZI*ZI + betaIZ*IX + betaZZ*ZX

#check
m = np.zeros(len(gate_sequence), dtype = complex)
for idx in range (len(gate_sequence)):
    m[idx] = np.trace(M[idx,:,:].dot(x))

epsilon = 0.1
rho = x*(1-epsilon) + epsilon*np.random.random((4, 4))   #initial guess
matrix_histogram_complex(x)
matrix_histogram_complex(rho)
dist = 0
for idx in range (len(gate_sequence)):
    dist = dist + m[idx] - np.trace(M[idx,:,:].dot(rho))
plt.show()
# measurement_matrix = np.array([[],
#                                [],
#                                [],
#                                [],
#                                [],
#                                [],
#                                ])
# measurement_matrix = np.empty((len(gate_sequence), 15), dtype=complex)
# gate_matrix = np.zeros((4, 4), dtype = complex)
# gate_matrix[0,0] = betaII
# for idx, gate in enumerate(gate_sequence):
#     gate_matrix = np.zeros((4, 4), dtype=complex)
#     gate1,gate2 = gate.split('-')
#     if gate1 == 'I':
#         row_index = 3
#         sign1 = 1
#     elif gate1[0:2] == 'Y2':
#         row_index = 1
#         if gate1[-1] == 'p':
#             sign1 = -1
#         elif gate1[-1] == 'm':
#             sign1 = 1
#     elif gate1[0:2] == 'X2':
#         row_index = 2
#         if gate1[-1] == 'm':
#             sign1 = -1
#         elif gate1[-1] == 'p':
#             sign1 = 1
#     elif gate1 in ['Xp', 'Xm', 'Yp', 'Ym']:
#         row_index = 3
#         sign1 = -1
#
#     if gate2 == 'I':
#         column_index = 3
#         sign2 = 1
#     elif gate2[0:2] == 'Y2':
#         column_index = 1
#         if gate2[-1] == 'p':
#             sign2 = -1
#         elif gate2[-1] == 'm':
#             sign2 = 1
#     elif gate2[0:2] == 'X2':
#         column_index = 2
#         if gate2[-1] == 'm':
#             sign2 = -1
#         elif gate2[-1] == 'p':
#             sign2 = 1
#     elif gate2 in ['Xp', 'Xm', 'Yp', 'Ym']:
#         column_index = 3
#         sign2 = -1
#
#     gate_matrix[0, column_index] = sign2*betaIZ
#     gate_matrix[row_index, 0] = sign1*betaZI
#     gate_matrix[row_index, column_index] = sign1 * sign2 * betaZZ
#     gate_matrix[0, 0] = betaII
#     measurement_matrix[idx, 0:3] = 0.25*gate_matrix[0, 1:]
#     measurement_matrix[idx, 3:7] = 0.25*gate_matrix[1, :]
#     measurement_matrix[idx, 7:11] = 0.25*gate_matrix[2, :]
#     measurement_matrix[idx, 11:15] = 0.25*gate_matrix[3, :]
#
#
#     # print (gate1, gate2, row_index, column_index, sign1, sign2)
#     # print (measurement_matrix[idx,:] - measurement_matrix_test[idx,:])
#
# avgIX, avgIY, avgIZ, \
# avgXI, avgXX, avgXY, avgXZ, \
# avgYI, avgYX, avgYY, avgYZ, \
# avgZI, avgZX, avgZY, avgZZ \
#     = np.linalg.inv(measurement_matrix).dot(np.array([
#     m0, m1, m2, m3,
#     m4, m5, m6, m7,
#     m8, m9, m10, m11,
#     m12, m13, m14]).transpose()-0.25*betaII).transpose()
#
#
# rho_reconstructed = 0.25*(tensor(qeye(2), qeye(2)) + avgIX*tensor(qeye(2), sigmax()) + avgIY*tensor(qeye(2), sigmay()) + avgIZ*tensor(qeye(2), sigmaz()) +
#                           avgXI*tensor(sigmax(), qeye(2)) + avgXX*tensor(sigmax(), sigmax()) + avgXY*tensor(sigmax(), sigmay()) + avgXZ*tensor(sigmax(), sigmaz()) +
#                           avgYI*tensor(sigmay(), qeye(2)) + avgYX*tensor(sigmay(), sigmax()) + avgYY*tensor(sigmay(), sigmay()) + avgYZ*tensor(sigmay(), sigmaz()) +
#                           avgZI*tensor(sigmaz(), qeye(2)) + avgZX*tensor(sigmaz(), sigmax()) + avgZY*tensor(sigmaz(), sigmay()) + avgZZ*tensor(sigmaz(), sigmaz()))
#
# #plotting
# # rho_target = ket2dm(bell_state(state='10'))
# matrix_histogram_complex(rho)
# # matrix_histogram_complex(rho_target)
# # matrix_histogram_complex(rho_measure)
# matrix_histogram_complex(rho_reconstructed)
# # print (fidelity(rho_target, rho))
# # print (fidelity(rho_target, rho_reconstructed))
# # print ((rho_reconstructed.dag()*rho_reconstructed).tr())
# print (expect(sigmaz(),basis(2,1)))
# plt.show()
