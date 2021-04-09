import numpy as np
import matplotlib.pyplot as plt
plt.style.use('kostas')

intensity = 1.2
MBc = [0., 0.01, 0.10, 0.30, 0.50]
MBMBC_lr = [0.038, 0.029, 0.029, 0.146, 2.072]
MBMBC_lr_err = [0.005, 0.003, 0.003, 0.018, 0.410]

MQMBc = [0., 0.10, 0.30]
MBMQMBC_lr = [0.060, 0.082, 0.366]
MBMQMBC_lr_err = [0.008, 0.012, 0.054]

intensity_MB = [0.40, 0.60, 0.80, 1.00, 1.20]
intensity_MBMQ = [0.40, 0.60, 0.80, 1.00, 1.20]
MBMQ_lr = [0.070, 1.622, 0.182, 0.102, 0.060]
MBMQ_lr_err = [0.008, 0.370, 0.036, 0.017, 0.008]

MB_lr = [0.078, 0.413, 0.043, 0.032, 0.038]
MB_lr_err = [0.009, 0.064, 0.004, 0.004, 0.005]


NOEC_lr_int = np.array([0.072, 0.048, 0.036, 0.029, 0.024])
NOEC_lr_err_int = np.array([0.010, 0.007, 0.005, 0.004, 0.003])
NOEC_lr = NOEC_lr_int[-1]
NOEC_lr_err = NOEC_lr_err_int[-1]

plt.close('all')
fig1 = plt.figure(1, figsize=[10, 9 / 1.6 * 1])
ax1 = fig1.add_subplot(111)
ax1.fill_between(MBc, np.ones_like(MBc) * (NOEC_lr - NOEC_lr_err),
                            np.ones_like(MBc) * (NOEC_lr + NOEC_lr_err),
                 alpha=0.5, color='k', label='$\\mathbf{No\ e\\textbf{-}cloud}$')
ax1.errorbar(MBc, MBMBC_lr, MBMBC_lr_err, c='b', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB}$')
ax1.errorbar(MQMBc, MBMQMBC_lr, MBMQMBC_lr_err, c='r', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB,\ MQ}$')

ax1.legend(loc='upper left', frameon=True)
ax1.set_ylabel('$\\mathbf{Loss\ rate\ [\\%/h]}$')
ax1.set_xlabel('$\\mathbf{Additional\ e\\textbf{-}density\, [10^{12}\,e^-/m^3]}$')
ax1.set_title('$\\mathbf{1.20\\cdot 10^{11}\,ppb}$')
ax1.set_xlim(-0.05, 0.65)
ax1.set_ylim(-0.05, 2.5)

fig2 = plt.figure(2, figsize=[10, 9 / 1.6 * 1])
ax2 = fig2.add_subplot(111)
ax2.fill_between(intensity_MBMQ, (NOEC_lr_int - NOEC_lr_err_int),
                 (NOEC_lr_int + NOEC_lr_err_int),
                 alpha=0.5, color='k', label='$\\mathbf{No\ e\\textbf{-}cloud}$')
ax2.errorbar(intensity_MB, MB_lr, MB_lr_err, c='b', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB}$')
ax2.errorbar(intensity_MBMQ, MBMQ_lr, MBMQ_lr_err, c='r', capsize=5, label='$\\mathbf{e\\textbf{-}cloud\ in\ MB,\ MQ}$')

ax2.legend(loc='upper right', frameon=True)
ax2.set_ylabel('$\\mathbf{Loss\ rate\ [\\%/h]}$')
ax2.set_xlabel('$\\mathbf{Bunch\ intensity\ [10^{11}\,p]}$')
ax2.set_title('$\\mathbf{}$')
ax2.set_xlim(0.35, 1.25)
ax2.set_ylim(-0.05, 2.5)

plt.show()
