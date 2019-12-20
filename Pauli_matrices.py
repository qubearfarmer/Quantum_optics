from qutip import*
from matplotlib import pyplot as plt
import numpy as np

bSphere = Bloch()
up = basis(2,0)
down = basis(2,1)

#basic
bSphere.add_states(up)
bSphere.add_states(down)
sup_state = (up+down)
bSphere.add_states(sup_state)

theta = np.pi/4
phi = np.pi/3
# arState = np.cos(theta/2)*up + np.sin(theta/2)*down*np.exp(1.0j*phi)
arVec = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]
#Arbitrary states
# bSphere.add_states(up)
# bSphere.add_vectors(arVec)
bSphere.make_sphere()

#Pauli matrices
# bSphere = Bloch()
# bSphere.add_states({arState,sigmax()*arState})
# bSphere.make_sphere()
plt.show()