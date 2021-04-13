from qutip import*
from matplotlib import pyplot as plt
import numpy as np

#basic
# plt.figure(1)
bSphere = Bloch()
up = basis(2,0)
down = basis(2,1)

# bSphere.add_states(up)
# bSphere.add_states(down)
# #
sup_state1 = (up+down)/np.sqrt(2)
bSphere.add_states(sup_state1)
bSphere.vector_color = ['b','r']
bSphere.add_annotation((up+down)/np.sqrt(2), text=r"$\frac{|0\rangle + |1\rangle}{\sqrt{2}}$")
#
sup_state2 = (up-down)/np.sqrt(2)
bSphere.add_states(sup_state2)
bSphere.add_annotation((up-down)/np.sqrt(2), text=r"$\frac{|0\rangle - |1\rangle}{\sqrt{2}}$")
#
sup_state3 = (up+np.exp(1j*np.pi/2)*down)/np.sqrt(2)
bSphere.add_states(sup_state3)
bSphere.add_annotation((up+np.exp(1j*np.pi/2)*down)/np.sqrt(2.2), text=r"$\frac{|0\rangle + i|1\rangle}{\sqrt{2}}$")
#
sup_state4 = (up-np.exp(1j*np.pi/2)*down)/np.sqrt(2)
bSphere.add_states(sup_state4)
bSphere.add_annotation((up-np.exp(1j*np.pi/2)*down)/np.sqrt(2), text=r"$\frac{|0\rangle - i|1\rangle}{\sqrt{2}}$")
#
bSphere.set_label_convention('01')
bSphere.sphere_color = '#FFFFFF'
bSphere.make_sphere()
#
# path = 'C:\\Users\\nguyen89\Google Drive\Research\Illustration\Thesis\Chapter 2\BlochSphere1.pdf'
# plt.savefig(path, dpi=300)

##################################################################################################
# theta = np.pi/3
# phi = np.pi/3
# arVec = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]
# bSphere.add_vectors(arVec)
# bSphere.add_annotation([np.sin(theta)*np.cos(phi)*1.1, np.sin(theta)*np.sin(phi), np.cos(theta)*1.3], text=r"$|\psi\rangle$")
# bSphere.sphere_color = '#FFFFFF'
# bSphere.make_sphere()

#################################################################################################
# sup_state1 = (up+down)/np.sqrt(2)
# bSphere.add_states(sup_state1)
# sup_state2 = rz(phi=0.15)*sup_state1
# bSphere.add_states(sup_state2)
# sup_state2 = rz(phi=0.3)*sup_state1
# bSphere.add_states(sup_state2)
# sup_state2 = rz(phi=-0.15)*sup_state1
# bSphere.add_states(sup_state2)
# sup_state2 = rz(phi=-0.3)*sup_state1
# bSphere.add_states(sup_state2)
# bSphere.vector_color = ['g']
# bSphere.sphere_color = '#FFFFFF'
# bSphere.make_sphere()


path = 'C:\\Users\\nguyen89\Google Drive\Research\Illustration\Thesis\Chapter 2\BlochSphere4.pdf'
plt.savefig(path, dpi=300)

plt.show()