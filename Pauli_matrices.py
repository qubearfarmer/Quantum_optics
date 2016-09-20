from qutip import*
from matplotlib import pyplot as plt

bSphere = Bloch()
up = basis(2,0)
down = basis(2,1)


theta = pi/4
phi = pi/3
arState = np.cos(theta/2)*up + np.sin(theta/2)*down*exp(1.0j*phi)
arVec = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]
#Arbitrary states
bSphere.add_states({up,down, arState})
bSphere.add_vectors(arVec)
bSphere.make_sphere()

#Pauli matrices
bSphere = Bloch()
bSphere.add_states({arState,sigmax()*arState})
bSphere.make_sphere()
plt.show()