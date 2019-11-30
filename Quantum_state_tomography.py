import numpy as np
from qutip import*
from matplotlib import pyplot as plt

rho1 = rand_dm(2)
rho2 = rand_dm(2)

matrix_histogram_complex(rho1)
rho = tensor(rho1,rho2)
matrix_histogram_complex(rho)

rho1_measure = 0.5*(qeye(2)+expect(sigmax(),rho1)*sigmax()+expect(sigmay(),rho1)*sigmay()+expect(sigmaz(),rho1)*sigmaz())
rho_measure = 0.25*(tensor(qeye(2),qeye(2))
                   +expect(tensor(qeye(2),sigmax()),rho)*tensor(qeye(2),sigmax())
                   +expect(tensor(qeye(2),sigmay()),rho)*tensor(qeye(2),sigmay())
                   +expect(tensor(qeye(2),sigmaz()),rho)*tensor(qeye(2),sigmaz())
                   +expect(tensor(sigmax(),qeye(2)),rho)*tensor(sigmax(),qeye(2))
                   +expect(tensor(sigmax(),sigmax()),rho)*tensor(sigmax(),sigmax())
                   +expect(tensor(sigmax(),sigmay()),rho)*tensor(sigmax(),sigmay())
                   +expect(tensor(sigmax(),sigmaz()),rho)*tensor(sigmax(),sigmaz())
                   +expect(tensor(sigmay(),qeye(2)),rho)*tensor(sigmay(),qeye(2))
                   +expect(tensor(sigmay(),sigmax()),rho)*tensor(sigmay(),sigmax())
                   +expect(tensor(sigmay(),sigmay()),rho)*tensor(sigmay(),sigmay())
                   +expect(tensor(sigmay(),sigmaz()),rho)*tensor(sigmay(),sigmaz())
                   +expect(tensor(sigmaz(),qeye(2)),rho)*tensor(sigmaz(),qeye(2))
                   +expect(tensor(sigmaz(),sigmax()),rho)*tensor(sigmaz(),sigmax())
                   +expect(tensor(sigmaz(),sigmay()),rho)*tensor(sigmaz(),sigmay())
                   +expect(tensor(sigmaz(),sigmaz()),rho)*tensor(sigmaz(),sigmaz())
                   )
matrix_histogram_complex(rho1_measure)
matrix_histogram_complex(rho_measure)


plt.show()
