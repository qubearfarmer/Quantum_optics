import numpy as np
from matplotlib import pyplot as plt

def gaussian_shape(x, sigma, T_gate):
    width = sigma*T_gate
    x_0 = width * 2
    y = np.zeros_like(x)
    for idx, t in enumerate (x):
        if t<=T_gate:
            y[idx] = np.exp(-0.5*(t-x_0)**2/width**2)
    return y

def gaussian_flattop(x, sigma, T_edge, T_gate):
    width = sigma * T_edge
    T_flat = T_gate - T_edge
    t_0 = sigma ** -1 * 0.5 * width
    xi_x = np.zeros_like(x)
    for idx, t in enumerate(x):
        if t <= t_0:
            xi_x[idx] = np.exp(-0.5 * (t - t_0) ** 2 / width ** 2)
        elif (t > t_0) and (t <= T_gate - t_0):
            xi_x[idx] = 1
        elif (t > T_gate - t_0) and (t <= T_gate):
            xi_x[idx] = np.exp(-0.5 * (t - (t_0 + T_flat)) ** 2 / width ** 2)
        else:
            xi_x[idx] = 0
    return xi_x

def drag_gaussian(x, sigma, T_gate, beta):
    y = gaussian_shape(x, sigma, T_gate)
    y_drag = beta*np.gradient(y)
    for idx, t in enumerate(x):
        if t>=T_gate:
            y_drag[idx] = 0
    return y_drag

def drag_gaussian1(x, sigma, T_gate, beta):
    y = np.zeros_like(x)
    width = sigma * T_gate
    x_0 = width * 2
    for idx, t in enumerate(x):
        if t<=T_gate:
            y[idx] = -np.exp(-0.5 * (t - x_0) ** 2 / width ** 2) * (t -x_0)/width**2
    return y*beta

def drag_gaussian_flattop(x, sigma, T_edge, T_gate, beta):
    y = gaussian_flattop(x, sigma, T_edge, T_gate)
    y_drag = beta * np.gradient(y)
    for idx, t in enumerate(x):
        if t >= T_gate:
            y_drag[idx] = 0
    return y_drag

def drag_gaussian_flattop1(x, sigma, T_edge, T_gate, beta):
    width = sigma * T_edge
    T_flat = T_gate - T_edge
    t_0 = sigma ** -1 * 0.5 * width
    xi_x = np.zeros_like(x)
    xi_y = np.zeros_like(x)
    for idx, t in enumerate(x):
        if t <= t_0:
            xi_x[idx] = np.exp(-0.5 * (t - t_0) ** 2 / width ** 2)
            xi_y[idx] = -np.exp(-0.5 * (t - t_0) ** 2 / width ** 2) * (t - t_0) / width ** 2
        elif (t > t_0) and (t <= T_gate - t_0):
            xi_x[idx] = 1
        elif (t > T_gate - t_0) and (t <= T_gate):
            xi_x[idx] = np.exp(-0.5 * (t - (t_0 + T_flat)) ** 2 / width ** 2)
            xi_y[idx] = -np.exp(-0.5 * (t - (t_0 + T_flat)) ** 2 / width ** 2) * (t - (t_0 + T_flat)) / width ** 2
    return xi_y*beta
x = np.linspace(0,300,301)
# y = np.linspace (1000, 1500, 101)
# print (y[np.where(300>x>200)])
freq = 5
# plt.plot(x, gaussian_shape(x, 0.25, 80), linewidth = 2.0)
# plt.plot(x, gaussian_shape(x, 0.25, 80)*np.cos(2*np.pi*freq*x))
plt.plot(x, gaussian_flattop(x, 0.25, 80, 200))
# plt.plot(x, drag_gaussian(x, 0.25, 80, 1.9))
# plt.plot(x, drag_gaussian1(x, 0.25, 80, 1.9))
plt.plot(x, drag_gaussian_flattop(x, 0.25, 80, 200, 1.9))
# plt.plot(x, drag_gaussian_flattop1(x, 0.25, 80, 200, 1.9))
# plt.tick_params(labelsize = 18)
plt.show()