import numpy as np
import matplotlib.pyplot as plt

intensity = [0.4, 0.6, 0.8, 1.0, 1.2]
MBMQ_e1 = [0.036, 0.328, 0.112, 0.090, 0.110]
MBMQ_e1_err = [0.016, 0.017, 0.021, 0.024, 0.018]
MBMQ_e2 = [0.123, 0.519, 0.153, 0.069, 0.109]
MBMQ_e2_err = [0.024, 0.016, 0.017, 0.020, 0.010]

MB_e1 = [0.022, 0.235, 0.018, 0.006, 0.001]
MB_e1_err = [0.011, 0.019, 0.016, 0.019, 0.012]
MB_e2 = [0.014, 0.577, 0.099, -0.020, 0.001]
MB_e2_err = [0.022, 0.014, 0.019, 0.010, 0.0028]

NOEC_e1 = 0.001
NOEC_e1_err = 0.012
NOEC_e2 = 0.001
NOEC_e2_err = 0.015

plt.close('all')
fig1 = plt.figure(1, figsize=[10, 9 / 1.6 * 2])
ax1 = fig1.add_subplot(211)
ax1.fill_between(intensity, np.ones_like(intensity) * (NOEC_e1 - NOEC_e1_err),
                            np.ones_like(intensity) * (NOEC_e1 + NOEC_e1_err),
                 alpha=0.5, color='k', label='$\\mathbf{No\ e\\textbf{-}cloud}$')
ax1.errorbar(intensity, MB_e1, MB_e1_err, c='b', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB}$')
ax1.errorbar(intensity, MBMQ_e1, MBMQ_e1_err, c='r', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB,\ MQ}$')

ax1.legend(loc='upper right', frameon=True)
ax1.set_ylabel('$\\mathbf{\\frac{d\\epsilon_x}{dt}\ [\\mu m/h]}$')
ax1.set_xlabel('$\\mathbf{Bunch\ intensity\ [10^{11}\,ppb]}$')
ax1.set_xlim(0.35, 1.25)
ax1.set_ylim(-0.05, 0.6)

ax2 = fig1.add_subplot(212)
ax2.fill_between(intensity, np.ones_like(intensity) * (NOEC_e2 - NOEC_e2_err),
                 np.ones_like(intensity) * (NOEC_e2 + NOEC_e2_err),
                 alpha=0.5, color='k', label='$\\mathbf{No\ e\\textbf{-}cloud}$')
ax2.errorbar(intensity, MB_e2, MB_e2_err, c='b', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB}$')
ax2.errorbar(intensity, MBMQ_e2, MBMQ_e2_err, c='r', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB,\ MQ}$')

ax2.legend(loc='upper right', frameon=True)
ax2.set_ylabel('$\\mathbf{\\frac{d\\epsilon_y}{dt}\ [\\mu m/h]}$')
ax2.set_xlabel('$\\mathbf{Bunch\ intensity\ [10^{11}\,ppb]}$')
ax2.set_xlim(0.35, 1.25)
ax2.set_ylim(-0.05, 0.6)
plt.show()
