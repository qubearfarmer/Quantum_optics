from qutip import*
import numpy as np
from matplotlib import pyplot as plt

gates = [['C-NOT', cnot()],
         ['C-Z', csign()],
         ['I-I', tensor(qeye(2), qeye(2))],
         ['Xp-Xp', tensor(sigmax(), sigmax())],
         ['X2p-X2p', tensor(rx(np.pi/2), rx(np.pi/2))],
         ['SWAP', swap()],
         ['$i$SWAP', iswap()],
         ['$\sqrt{i\mathrm{SWAP}}$', sqrtiswap()],
         ['S-NOT', snot()],
         ['$\pi/2$ phase gate', phasegate(np.pi/2)],
         ['Toffoli', toffoli()],
         ['Fredkin', fredkin()],
         ['Xp', sigmax()],
         ['Z2p', rz(phi=np.pi/2)]]

def plt_qpt_gate(gate, figsize=(8,6)):

    name = gate[0]
    U_psi = gate[1]

    N = len(U_psi.dims[0])  # number of qubits

    # create a superoperator for the density matrix
    # transformation rho = U_psi * rho_0 * U_psi.dag()
    U_rho = spre(U_psi) * spost(U_psi.dag())

    # operator basis for the process tomography
    op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()] for i in range(N)]

    # labels for operator basis
    op_label = [['',"$I$", "$X$", "$Y$", "$Z$"] for i in range(N)]
    print(op_label)

    # calculate the chi matrix
    chi = qpt(U_rho, op_basis)

    # visualize the chi matrix
    fig, ax = qpt_plot_combined(chi, op_label, name, figsize=figsize)

    ax.set_title(name)

    return fig, ax

plt_qpt_gate(gates[4])

plt.show()