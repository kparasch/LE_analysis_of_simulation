import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.style.use('kostas')

mode = 'MBMQMBC'
intensity = 1.20
pdf = PdfPages(f'emittance_growth_{mode}.pdf')
for MBc in [0., 0.10, 0.30]:
#for MBc in [0., 0.01, 0.10, 0.30, 0.50]:

    if MBc != 0.:
        fname = f'Emittance_DATA{mode}_{MBc:.2f}MBc_{intensity:.2f}e11ppb'
    else:
        fname = f'Emittance_DATAMBMQ_{intensity:.2f}e11ppb'

    emit_dict = pickle.load(open('data/' + fname + '.pkl','rb'))

    if mode == 'MBMBC' and MBc == 0.50:
        mask = emit_dict['time']*60. < 13.0
        emit_dict['time'] = emit_dict['time'][mask]
        emit_dict['e1'] = emit_dict['e1'][mask]
        emit_dict['e2'] = emit_dict['e2'][mask]
        emit_dict['e1_err'] = emit_dict['e1_err'][mask]
        emit_dict['e2_err'] = emit_dict['e2_err'][mask]

    emit_dict['e1'] *= 2.
    emit_dict['e2'] *= 2.
    emit_dict['e1_err'] *= 2.
    emit_dict['e2_err'] *= 2.

    plt.close('all')
    fig1 = plt.figure(1, figsize=[10,9/1.6*2])
    ax1 = fig1.add_subplot(211)
    ax1.fill_between(emit_dict['time']*60., emit_dict['e1'] - emit_dict['e1_err'], emit_dict['e1'] + emit_dict['e1_err'], alpha=0.5, color='b')
    ax1.plot(emit_dict['time']*60., emit_dict['e1'], 'b')
    p_e1, cov_e1 = np.polyfit(emit_dict['time'], emit_dict['e1'], 1, w=1. / emit_dict['e1_err'],
                              cov='unscaled')
    ax1.plot(emit_dict['time']*60., p_e1[0] * emit_dict['time'] + p_e1[1], 'r',
             label='$\\mathbf{\\frac{\\mathbf d \\epsilon_x}{\\mathbf dt} = '+f'{p_e1[0]:.3f} \\pm {np.sqrt(cov_e1[0,0]):.3f} '+'\ \\mu m/h}$')
    ax1.legend(loc='upper left', frameon=True)
    ax1.set_ylabel('$\\mathbf{\\epsilon_x\ [\\mu m]}$')
    ax1.set_xlabel('$\\mathbf{Time\ [min]}$')
    ax1.set_xlim(-0.5, 16.)

    ax2 = fig1.add_subplot(212)
    ax2.fill_between(emit_dict['time']*60., emit_dict['e2'] - emit_dict['e2_err'], emit_dict['e2'] + emit_dict['e2_err'], alpha=0.5, color='b')
    ax2.plot(emit_dict['time']*60., emit_dict['e2'], 'b')
    p_e2, cov_e2 = np.polyfit(emit_dict['time'], emit_dict['e2'], 1, w=1. / emit_dict['e2_err'],
                              cov='unscaled')
    ax2.plot(emit_dict['time']*60., p_e2[0] * emit_dict['time'] + p_e2[1], 'r',
             label='$\\mathbf{\\frac{\\mathbf d \\epsilon_y}{\\mathbf dt} = '+f'{p_e2[0]:.3f} \\pm {np.sqrt(cov_e2[0,0]):.3f} '+'\ \\mu m/h}$')
    ax2.legend(loc='upper left', frameon=True)
    ax2.set_ylabel('$\\mathbf{\\epsilon_y\ [\\mu m]}$')
    ax2.set_xlabel('$\\mathbf{Time\ [min]}$')
    ax2.set_xlim(-0.5, 16.)
    fig1.suptitle('$\\mathbf{'+ mode + f', \ 1.20e12ppb, \ +{MBc:.2f}e12\,e^-/m^3' +'}$')
    pdf.savefig(fig1)

pdf.close()
plt.show()