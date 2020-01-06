from qutip import*
from qutip.qip import*
import sys
import numpy as np
sys.path.append('C:\\Users\\nguyen89\Documents\Python Codes\ImageMagick')
# import imagemagick


psi_0 = tensor(basis(2,1), basis(2,1))
psi_1 = tensor(ry(phi=np.pi/2),qeye(2))*psi_0
psi_2 = cnot()*psi_1
print(psi_2)
