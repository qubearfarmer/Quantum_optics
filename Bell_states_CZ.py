from qutip import*
from qutip.qip import*
import sys
import numpy as np
from matplotlib import pyplot as plt
sys.path.append('C:\\Users\\nguyen89\Documents\Python Codes\ImageMagick')
# import imagemagick


psi_0 = tensor(basis(2,1), basis(2,1))
psi_1 = tensor(rx(phi=np.pi/2), rx(phi=np.pi/2))*psi_0
psi_2 = csign()*psi_1
psi_3 = tensor(rx(phi=np.pi/2), qeye(2))*psi_2
matrix_histogram_complex(ket2dm(psi_3))
matrix_histogram_complex(ket2dm(bell_state('01')))
plt.show()

