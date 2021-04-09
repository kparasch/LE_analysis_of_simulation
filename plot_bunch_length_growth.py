import numpy as np
import matplotlib.pyplot as plt
import pickle
plt.style.use('kostas')

plt.close('all')
fig1 = plt.figure(1, figsize=[10, 9/1.6])
ax1 = fig1.add_subplot(111)
sets = []
#sets.append({'q':1.20, 'intensity': 1.20, 'mode': 'MBMQ', 'color': 'b',
#             'label': '$\\mathbf{1.20\cdot 10^{11}\,ppb}$'})
#sets.append({'q':1.20, 'intensity': 1.00, 'mode': 'MBMQ', 'color': 'g',
#             'label': '$\\mathbf{1.00\cdot 10^{11}\,ppb}$'})
#sets.append({'q':1.20, 'intensity': 0.80, 'mode': 'MBMQ', 'color': 'm',
#             'label': '$\\mathbf{0.80\cdot 10^{11}\,ppb}$'})
#sets.append({'q':1.20, 'intensity': 0.60, 'mode': 'MBMQ', 'color': 'r',
#             'label': '$\\mathbf{0.60\cdot 10^{11}\,ppb}$'})
#sets.append({'q':1.20, 'intensity': 0.00, 'mode': 'NOEC', 'color': 'k',
#             'label': '$\\mathbf{0.00\cdot 10^{11}\,ppb}$'})
sets.append({'q':1.20, 'intensity': 1.20, 'mode': 'MBMQ', 'color': 'b',
             'label': '$\\mathbf{1.20\cdot 10^{11}\,ppb}$'})
sets.append({'q':1.20, 'intensity': 1.20, 'mode': 'MBMBC', 'color': 'm',
             'label': '$\\mathbf{1.20\cdot 10^{11}\,ppb, +0.01\,e^-/m^3}$', 'MBc': 0.01})
sets.append({'q':1.20, 'intensity': 1.20, 'mode': 'MBMBC', 'color': 'g',
             'label': '$\\mathbf{1.20\cdot 10^{11}\,ppb, +0.10\,e^-/m^3}$', 'MBc': 0.10})
sets.append({'q':1.20, 'intensity': 0.00, 'mode': 'NOEC', 'color': 'k',
             'label': '$\\mathbf{0.00\cdot 10^{11}\,ppb}$'})

for thisset in sets:
    q=thisset['q']
    intensity = thisset['intensity']
    mode = thisset['mode']
    label = thisset['label']
    color = thisset['color']
    if 'MBc' in thisset.keys():
        MBc = thisset['MBc']
    else:
        MBc = 0.

    if mode == 'MBMQ' and intensity==1.20 and q==1.20:
        emit = np.array([1.5, 2, 2.5, 3.0, 3.5, 4., 4.5])
        bl_growth = np.array([-1.3e-4, -2.6e-4,-1.25e-3, -3.58e-3, -6.17e-3, -8.50e-3, -13.84e-3])
        bl_growth_err = np.array([5.3e-4, 5.33e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.2e-4])
    elif mode == 'MBMQ' and intensity == 1.00 and q == 1.20:
        emit = np.array([1.5, 2, 2.5, 3.0, 3.5, 4., 4.5])
        bl_growth = np.array([-4.4e-4, -0.2e-4, -1.08e-3, -3.03e-3, -5.42e-3, -8.89e-3, -15.08e-3])
        bl_growth_err = np.array([5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4])
    elif mode == 'MBMQ' and intensity == 0.80 and q == 1.20:
        emit = np.array([1.5, 2, 2.5, 3.0, 3.5, 4., 4.5])
        bl_growth = np.array([-5.6e-4, -15.6e-4, -1.07e-3, -3.10e-3, -4.51e-3, -10.84e-3, -15.51e-3])
        bl_growth_err = np.array([5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.2e-4, 5.2e-4])
    elif mode == 'MBMQ' and intensity == 0.60 and q == 1.20:
        emit = np.array([1.5, 2, 2.5, 3.0, 3.5, 4., 4.5])
        bl_growth = np.array([-3.5e-4, -7.6e-4, -2.00e-3, -3.21e-3, -6.77e-3, -12.06e-3, -21.68e-3])
        bl_growth_err = np.array([5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4])
    elif mode == 'MBMBC' and intensity == 1.20 and q == 1.20 and MBc == 0.01:
        emit = np.array([1.5, 2, 2.5, 3.0, 3.5, 4., 4.5])
        bl_growth = np.array([-8.3e-4, -0.7e-4, -1.19e-3, -1.68e-3, -4.24e-3, -5.54e-3, -11.68e-3])
        bl_growth_err = np.array([5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.2e-4])
    elif mode == 'MBMBC' and intensity == 1.20 and q == 1.20 and MBc == 0.10:
        emit = np.array([1.5, 2, 2.5, 3.0, 3.5, 4., 4.5])
        bl_growth = np.array([-3.4e-4, -0.7e-4, -1.20e-3, -1.93e-3, -4.29e-3, -6.44e-3, -11.71e-3])
        bl_growth_err = np.array([5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.2e-4])
    elif mode == 'MBMBC' and intensity == 1.20 and q == 1.20 and MBc == 0.50:
        emit = np.array([1.5, 2, 2.5, 3.0, 3.5, 4., 4.5])
        bl_growth = np.array([-5.1e-4, -1.45e-3, -2.94e-3, -5.27e-3, -11.75e-3, -19.24e-3, -33.49e-3])
        bl_growth_err = np.array([5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.2e-4, 5.2e-4, 5.0e-4])
    elif mode == 'NOEC' and intensity==0. and q==1.20:
        emit = np.array([1.5, 2, 2.5, 3.0, 3.5, 4.0, 4.5])
        bl_growth = np.array([-3.6e-4, -6.5e-4, -9.4e-4, -8.8e-4, -3.4e-3, -4.22e-3, -8.76e-3])
        bl_growth_err = np.array([5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4, 5.3e-4])
    ax1.errorbar(emit, 1.e3*bl_growth, 1.e3*bl_growth_err, c=color, fmt='.-', capsize=5, label=label)



ax1.set_xlabel('$\\mathbf{Emittance\ [\mu m]}$')
ax1.set_ylabel('$\\mathbf{\\frac{d\\sigma_{\\tau}}{dt}\ [mm/h]}$')
ax1.legend()
plt.show()