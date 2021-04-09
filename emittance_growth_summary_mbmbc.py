import numpy as np
import matplotlib.pyplot as plt

intensity = 1.2
MBc = [0., 0.01, 0.10, 0.30, 0.50]
MBMBC_e1 = [0.001, 0.004, -0.005, 0.167, 0.291]
MBMBC_e1_err = [0.012, 0.035, 0.015, 0.029, 0.022]
MBMBC_e2 = [0.001, 0.000, 0.051, 0.294, 0.941]
MBMBC_e2_err = [0.028, 0.018, 0.021, 0.021, 0.016]

MQMBc = [0., 0.10, 0.30]
MBMQMBC_e1 = [0.110, 0.122, 0.243]
MBMQMBC_e1_err = [0.018, 0.017, 0.015]
MBMQMBC_e2 = [0.109, 0.013, 0.401]
MBMQMBC_e2_err = [0.010, 0.016, 0.027]

NOEC_e1 = 0.001
NOEC_e1_err = 0.012
NOEC_e2 = 0.001
NOEC_e2_err = 0.015

plt.close('all')
fig1 = plt.figure(1, figsize=[10, 9 / 1.6 * 2])
ax1 = fig1.add_subplot(211)
ax1.fill_between(MBc, np.ones_like(MBc) * (NOEC_e1 - NOEC_e1_err),
                            np.ones_like(MBc) * (NOEC_e1 + NOEC_e1_err),
                 alpha=0.5, color='k', label='$\\mathbf{No\ e\\textbf{-}cloud}$')
ax1.errorbar(MBc, MBMBC_e1, MBMBC_e1_err, c='b', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB}$')
ax1.errorbar(MQMBc, MBMQMBC_e1, MBMQMBC_e1_err, c='r', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB,\ MQ}$')

ax1.legend(loc='upper right', frameon=True)
ax1.set_ylabel('$\\mathbf{\\frac{d\\epsilon_x}{dt}\ [\\mu m/h]}$')
ax1.set_xlabel('$\\mathbf{Additional\ e\\textbf{-}density\, [10^{12}\,e^-/m^3]}$')
ax1.set_xlim(-0.05, 0.65)
ax1.set_ylim(-0.05, 0.7)

ax2 = fig1.add_subplot(212)
ax2.fill_between(MBc, np.ones_like(MBc) * (NOEC_e2 - NOEC_e2_err),
                 np.ones_like(MBc) * (NOEC_e2 + NOEC_e2_err),
                 alpha=0.5, color='k', label='$\\mathbf{No\ e\\textbf{-}cloud}$')
ax2.errorbar(MBc, MBMBC_e2, MBMBC_e2_err, c='b', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB}$')
ax2.errorbar(MQMBc, MBMQMBC_e2, MBMQMBC_e2_err, c='r', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB,\ MQ}$')

ax2.legend(loc='upper right', frameon=True)
ax2.set_ylabel('$\\mathbf{\\frac{d\\epsilon_y}{dt}\ [\\mu m/h]}$')
ax2.set_xlabel('$\\mathbf{Additional\ e\\textbf{-}density\, [10^{12}\,e^-/m^3]}$')
ax2.set_xlim(-0.05, 0.65)
ax2.set_ylim(-0.05, 1.7)
plt.show()
