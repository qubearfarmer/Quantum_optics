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

def gaussian_haonan(x, width):
    sigma = width / np.sqrt(2 * np.pi)
    y = np.zeros_like(x)
    y_offset = np.exp(-(-width)**2/(2*sigma**2))
    rescale_factor = 1/(1-y_offset)
    for idx, t in enumerate(x):
        if t<=2*width:
            y[idx] = (np.exp(-(x[idx]-width)**2/(2*sigma**2))-y_offset)*rescale_factor
    return y

def gaussian_flattop_haonan(x, width, T_flat):
    sigma = width / np.sqrt(2 * np.pi)
    y = np.zeros_like(x)
    y_offset = np.exp(-(-width) ** 2 / (2 * sigma ** 2))
    rescale_factor = 1 / (1 - y_offset)
    for idx, t in enumerate(x):
        if t <= width:
            y[idx] = (np.exp(-(x[idx] - width) ** 2 / (2 * sigma ** 2)) - y_offset)*rescale_factor
        elif t>width and t<width +T_flat:
            y[idx] = 1
        elif t>=width +T_flat and t <= 2*width +T_flat:
            y[idx] = (np.exp(-(x[idx] - width-T_flat) ** 2 / (2 * sigma ** 2))-y_offset)*rescale_factor
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
    width = sigma * T_gate
    x_0 = width * 2
    y = np.zeros_like(x)
    for idx, t in enumerate(x):
        if t <= T_gate:
            y[idx] = -np.exp(-0.5 * (t - x_0) ** 2 / width ** 2)*(t-x_0)/width**2
    return y

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
    y_offset = np.exp(-(-width) ** 2 / (2 * sigma ** 2))
    rescale_factor = 1 / (1 - y_offset)
    for idx, t in enumerate(x):
        if t <= t_0:
            xi_x[idx] = (np.exp(-0.5 * (t - t_0) ** 2 / width ** 2)-y_offset)*rescale_factor
            xi_y[idx] = -rescale_factor*np.exp(-0.5 * (t - t_0) ** 2 / width ** 2) * (t - t_0) / width ** 2
        elif (t > t_0) and (t <= T_gate - t_0):
            xi_x[idx] = 1
        elif (t > T_gate - t_0) and (t <= T_gate):
            xi_x[idx] = (np.exp(-0.5 * (t - (t_0 + T_flat)) ** 2 / width ** 2)-y_offset)*rescale_factor
            xi_y[idx] = -rescale_factor*np.exp(-0.5 * (t - (t_0 + T_flat)) ** 2 / width ** 2) * (t - (t_0 + T_flat)) / width ** 2
    return xi_y*beta

def drag_gaussian_flattop_haonan(x, width, T_flat, alpha):
    xi_x = np.zeros_like(x)
    xi_y = np.zeros_like(x)
    amplitude = 1
    sigma = width / np.sqrt(2 * np.pi)
    y_offset =  np.exp(-(-width) ** 2 / (2 * sigma ** 2))
    rescale_factor = 1.0 / (1.0 - y_offset)
    for idx, t in enumerate (x):
        if t <= width:
            xi_x[idx] = (np.exp(-0.5 * (t - width) ** 2 / sigma ** 2) - y_offset) * rescale_factor
            xi_y[idx] = (-np.exp(-0.5 * (t - width) ** 2 / sigma ** 2) * (t - width) / sigma ** 2) * rescale_factor
        elif (t > width) and (t < width + T_flat):
            xi_x[idx] = 1
            xi_y[idx] = 0
        elif (t >= width + T_flat) and (t <= 2 * width + T_flat):
            xi_x[idx] = (np.exp(-0.5 * (t - width - T_flat) ** 2 / sigma ** 2) - y_offset) * rescale_factor
            xi_y[idx] = (-np.exp(-0.5 * (t - width - T_flat) ** 2 / sigma ** 2) * (
                        t - width - T_flat) / sigma ** 2) * rescale_factor
        else:
            xi_x[idx] = 0
            xi_y[idx] = 0
    xi_y = alpha * xi_y
    return xi_x, xi_y

x = np.linspace(0,200,201)
# y = np.linspace (1000, 1500, 101)
# print (y[np.where(300>x>200)])
freq = 5
# plt.plot(x, gaussian_shape(x, 0.25, 40), linewidth = 2.0)
# plt.plot(x, gaussian_shape(x, 0.25, 80)*np.cos(2*np.pi*freq*x))
# plt.plot(x, gaussian_haonan(x, 10), linewidth = 2.0)
# plt.plot(x, gaussian_flattop(x, 0.25, 40, 100))
# plt.plot(x, gaussian_flattop_haonan(x, 10, 50))
# plt.plot(x, drag_gaussian(x, 0.25, 40, 1)/3)
# plt.plot(x, np.gradient(gaussian_shape(x, 0.25, 40)))
# plt.plot(x, drag_gaussian(x, 0.25, 40, 1)/np.gradient(gaussian_shape(x, 0.25, 40)))
# plt.plot(x, gaussian_shape(x, 0.25, 80))
# plt.plot(x[1:], np.diff(gaussian_shape(x, 0.25, 80)))
# plt.plot(x, drag_gaussian1(x, 0.25, 80, 1.9))
# plt.plot(x, drag_gaussian_flattop(x, 0.25, 80, 200, 1.9))
# plt.plot(x, drag_gaussian_flattop1(x, 0.25, 40, 100, 1))
plt.plot(x, drag_gaussian_flattop_haonan(x, 20, 70, 1.9)[0], label = 'Long pulse')
plt.plot(x, drag_gaussian_flattop_haonan(x, 20, 70, 1.9)[1], label = 'Long DRAG')
plt.plot(x, 1.9*np.gradient(drag_gaussian_flattop_haonan(x, 20, 70, 1.9)[0]), label = 'Long gradient')
# plt.tick_params(labelsize = 18)

t = np.linspace(0, 200, 201)
num_t = len(t)
t_step = t[1] - t[0]
t0_ind = int(num_t / 2)
t0 = t[t0_ind]
width = 20
truncation_range = 2
plateau = 70
total_duration = truncation_range * width + plateau
std = width / np.sqrt(2 * np.pi)
values1 = np.exp(-(t - t0 - plateau / 2) ** 2 / (2 * std ** 2))
values1[t < (t0 + plateau / 2)] = 1
values2 = np.exp(-(t - t0 + plateau / 2) ** 2 / (2 * std ** 2))
values2[t > (t0 - plateau / 2)] = 1
values = values1 * values2
# plateau_ind = (t > (t0 - plateau / 2)) * (t < (t0 + plateau / 2))
# values[plateau_ind] = values.max()
# plt.plot(t, np.exp(-(t - t0 - plateau / 2) ** 2 / (2 * std ** 2)))
values[t < (t0 - total_duration / 2)] = 0
values[t > (t0 + total_duration / 2)] = 0
non_zero_value = values[values != 0]
values = values - non_zero_value.min()
values[t < (t0 - total_duration / 2)] = 0
values[t > (t0 + total_duration / 2)] = 0
pulse_envelope = values / values.max()
plt.plot(t,pulse_envelope, label = 'Haonan pulse')
plt.plot(t,1.9*np.gradient(pulse_envelope), label = 'Haonan gradient')
plt.legend()
plt.show()