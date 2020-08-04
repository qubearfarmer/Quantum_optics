from qutip import*
from qutip.qip import*
import sys
import numpy as np
from matplotlib import pyplot as plt
sys.path.append('C:\\Users\\nguyen89\Documents\Python Codes\ImageMagick')
# import imagemagick


# psi_0 = tensor(basis(2,1), basis(2,1))
# psi_1 = tensor(ry(phi=np.pi/2),qeye(2))*psi_0
# psi_2 = cnot()*psi_1
# matrix_histogram_complex(ket2dm(psi_2))
matrix_histogram_complex(ket2dm(bell_state('01')))
plt.show()

# state = (basis(2,0) -1j* basis(2,1))/np.sqrt(2)
# dm = ket2dm(state)
# print (dm)
# matrix_histogram_complex(dm)
# plt.show()