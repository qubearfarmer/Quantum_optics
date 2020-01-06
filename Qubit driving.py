# Project: Quantum tomography
# author: Long Nguyen
# Date: 12-1-2019
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack

# Number of samplepoints
# N = 600
# # sample spacing
# T = 1.0 / 800.0
# x = np.linspace(0.0, N*T, N)
# y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)
# yf = scipy.fftpack.fft(y)
# xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
#
# fig, ax = plt.subplots()
# ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))

t = np.linspace(0,10,1001)
x = np.cos(t)
y = np.cos(10*t)
# plt.plot(t,x)
# plt.plot(t,y)
plt.plot(t,x*y,linewidth = 2.0)
# plt.plot(t,x+y,linewidth = 2.0)
plt.show()