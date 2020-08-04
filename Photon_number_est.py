import numpy as np

power = -130 #dBm
power = 10**(power/10) * 1e-3
kappa = 5e6
energy = power/kappa
freq = 7.5e9
h = 6.626e-34
n = energy/(h*freq)
print (n)
print (power)