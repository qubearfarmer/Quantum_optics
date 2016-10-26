import matplotlib.pyplot as plt
import numpy as np
width = 100
x = np.linspace(0,4*width,4*width)
gaussian_arr = np.power(2, -1 * np.power(2*(x - 2*width)/width, 2))
plt.plot(x,gaussian_arr)
plt.show()