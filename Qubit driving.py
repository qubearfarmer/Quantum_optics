import numpy as np
from matplotlib import pyplot as plt
from qutip import*

def qubit_integrate(epsilon, delta, g1, g2, solver):
    H = epsilon / 2.0 * sigmaz() + delta / 2.0 * sigmax()

    # collapse operators
    c_ops = []

    if g1 > 0.0:
        c_ops.append(np.sqrt(g1) * sigmam())

    if g2 > 0.0:
        c_ops.append(np.sqrt(g2) * sigmaz())

    e_ops = [sigmax(), sigmay(), sigmaz()]

    if solver == "me":
        output = mesolve(H, psi0, tlist, c_ops, e_ops)
    elif solver == "es":
        output = essolve(H, psi0, tlist, c_ops, e_ops)
    elif solver == "mc":
        ntraj = 250
        output = mcsolve(H, psi0, tlist, ntraj, c_ops, [sigmax(), sigmay(), sigmaz()])
    else:
        raise ValueError("unknown solver")

    return output.expect[0], output.expect[1], output.expect[2]

epsilon = 0.0 * 2 * np.pi   # cavity frequency
delta   = 0.0 * 2 * np.pi   # atom frequency
g2 = 0.0
g1 = 1

# intial state
psi0 = basis(2,0)

tlist = np.linspace(0,5,200)

sx1, sy1, sz1 = qubit_integrate(epsilon, delta, g1, g2, "me")
plt.plot(tlist , np.real(sz1))

plt.show()