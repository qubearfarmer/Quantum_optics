from qutip import*
from qutip.qip import*
import sys
import numpy as np
sys.path.append('C:\\Users\\nguyen89\Documents\Python Codes\ImageMagick')
# import imagemagick


psi_0 = tensor(basis(2,0), basis(2,0))
psi_1 = tensor(ry(phi=np.pi/2), ry(phi=np.pi/2))*psi_0
psi_2 = csign()*psi_1
psi_3 = tensor(ry(phi=np.pi/2), qeye(2))*psi_2
# print(psi_3)
# matrix_histogram_complex((ket2dm(psi_3)))

