import numpy as np
import matplotlib.pyplot as plt
import pickle
import h5py
#import ScientificColourMaps6 as SCM6
import scipy.special
from scipy.constants import c
import scipy.optimize
import LE_analysis
import longitudinal_distribution

#plt.style.use('kostas')
#plt.close('all')

#mycmap = SCM6.bamako
mode = 'MBMQ'
#mycolors = mycmap(list(map(int,np.linspace(0, 256, ncurves))))
ptau_list = [0.00, 2.94, 4.13, 5.01, 5.73, 6.34, 6.87, 7.33, 7.73, 8.06]
#ptau = 8.06
intensity=1.20
#for MBc in [0.01, 0.10, 0.30, 0.50]:
MBc=0.
for intensity in [0.40, 0.60, 0.80, 1.00, 1.20]:
    print(f'Mode: {mode}, intensity = {intensity:.2f} e11ppb, MBc = {MBc:.2f}')
    # if mode == 'MBMBC' or mode == 'MBMQMBC':
    #     out_fname = f'data/Emittance_DATA{mode}_{MBc:.2f}MBc_{intensity:.2f}e11ppb.pkl'
    # else:
    #     out_fname = f'data/Emittance_DATA{mode}_{intensity:.2f}e11ppb.pkl'
    # print(out_fname)
    # for intensity in [0.40, 0.60, 0.80, 1.00, 1.20]:
    ##intensity=0.60
    plt.close('all')
    in_folder = '/eos/project/e/ecloud-simulations/kparasch/Norm_tracking_data/'
    #in_file = f'LE_LHC_450GeV_DATAseyMB_DATAseyMQ_{intensity:.2f}e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau:.2f}e-4ptau_norm.h5'
    #in_files = [f'LE_LHC_450GeV_DATAseyMB_{intensity:.2f}e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau:.2f}e-4ptau_norm.h5' for ptau in ptau_list]
    if mode == 'MB':
        in_files = [f'LE_LHC_450GeV_DATAseyMB_{intensity:.2f}e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau:.2f}e-4ptau_norm.h5' for ptau in ptau_list]
    if mode == 'MBMBC':
        in_files = [
            f'LE_LHC_450GeV_DATAseyMB_{MBc:.2f}MBc_{intensity:.2f}e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau:.2f}e-4ptau_norm.h5'
            for ptau in ptau_list]
    if mode == 'MBMQMBC':
        in_files = [f'LE_LHC_450GeV_DATAseyMB_DATAseyMQ_{MBc:.2f}MBc_{intensity:.2f}e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau:.2f}e-4ptau_norm.h5'
            for ptau in ptau_list]
    elif mode == 'MBMQ':
        in_files = [f'LE_LHC_450GeV_DATAseyMB_DATAseyMQ_{intensity:.2f}e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau:.2f}e-4ptau_norm.h5' for ptau in ptau_list]
    elif mode == 'NOEC':
        #intensity=0.00
        in_files = [f'LE_LHC_450GeV_0.00e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau:.2f}e-4ptau_norm.h5' for ptau in ptau_list]



    files = [h5py.File(in_folder + in_file, 'r') for in_file in in_files]

    e1_ref = 2.e-6/479.603977
    e2_ref = 2.e-6/479.603977
    q_ref = 1.2
    bunch_length_ref = 0.09

    fig1 = plt.figure(1, figsize=[10*2, 9/1.6*2])
    ax1 = fig1.add_subplot(221)
    ax2 = fig1.add_subplot(222)
    ax3 = fig1.add_subplot(223)
    ax4 = fig1.add_subplot(224)

    ptau = np.array(ptau_list)
    this_file = files[0]
    J1_0 = this_file['turn0']['J1'][()]
    J2_0 = this_file['turn0']['J2'][()]

    optics = this_file['optics']
    A = optics['rf_volt_V'][()]/(2.*np.pi*optics['rf_freq_Hz'][()]*optics['p0c_eV'][()]/c*optics['length'][()])
    B = 2.*np.pi*optics['rf_freq_Hz'][()]/c
    D = (optics['alfa'][()] - 1./optics['gamma0'][()]**2)/2.

    m = D/2./A*(ptau*1.e-4)**2
    K = scipy.special.ellipk(m)
    Jz = 4*np.sqrt(2*A)/(np.pi*B*np.sqrt(D)) * ( scipy.special.ellipe(m) - ( 1 - m ) * scipy.special.ellipk(m))

    loss_dict = LE_analysis.combine_losses(list(map(lambda x: in_folder + x, in_files)), Jz,
                                           e1=e1_ref, e2=e2_ref, bunch_length=bunch_length_ref, q=q_ref)
    loss_rate = np.mean(loss_dict['tot_L'][-5:] ** 2) ** 0.5
    loss_rate_err = np.mean(loss_dict['tot_Lerr'][-5:] ** 2) ** 0.5
    #e9p/h -> perc/h
    loss_rate /= intensity
    loss_rate_err /= intensity
    print(f'{loss_rate:.3f} +- {loss_rate_err:.3f}')
plt.show()
