import numpy as np
import matplotlib.pyplot as plt
import pickle
import h5py
#import ScientificColourMaps6 as SCM6
import scipy.special
from scipy.constants import c
import scipy.optimize
import longitudinal_distribution

#plt.style.use('kostas')
#plt.close('all')

#mycmap = SCM6.bamako
ncurves = 1001
#ncurves = 201
#mycolors = mycmap(list(map(int,np.linspace(0, 200, ncurves))))
mode = 'MBMQMBC'
#mycolors = mycmap(list(map(int,np.linspace(0, 256, ncurves))))
ptau_list = [0.00, 2.94, 4.13, 5.01, 5.73, 6.34, 6.87, 7.33, 7.73, 8.06]
#ptau = 8.06
intensity=1.20
for MBc in [0.10, 0.30]:
    if mode == 'MBMBC' or mode == 'MBMQMBC':
        out_fname = f'data/Emittance_DATA{mode}_{MBc:.2f}MBc_{intensity:.2f}e11ppb.pkl'
    else:
        out_fname = f'data/Emittance_DATA{mode}_{intensity:.2f}e11ppb.pkl'
    print(out_fname)
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
        intensity=0.00
        in_files = [f'LE_LHC_450GeV_0.00e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau:.2f}e-4ptau_norm.h5' for ptau in ptau_list]



    #this_file = h5py.File(in_folder + in_file, 'r')
    files = [h5py.File(in_folder + in_file, 'r') for in_file in in_files]

    e1_ref = 2.e-6/479.603977
    e2_ref = 2.e-6/479.603977

    bins = 50
    last_sigma = 6
    density=True
    Jrange = (0, last_sigma**2/2.)
    xrange = (-last_sigma, last_sigma)
    fig1 = plt.figure(1, figsize=[10*2, 9/1.6*2])
    ax1 = fig1.add_subplot(221)
    ax2 = fig1.add_subplot(222)
    ax3 = fig1.add_subplot(223)
    ax4 = fig1.add_subplot(224)

    rho = {}
    for ptau, this_file in zip(ptau_list, files):
        J1_0 = this_file['turn0']['J1'][()]
        J2_0 = this_file['turn0']['J2'][()]
        e1_init = 2.0 * (e1_ref/2.)
        e2_init = 2.0 * (e2_ref/2.)

        optics = this_file['optics']
        A = optics['rf_volt_V'][()]/(2.*np.pi*optics['rf_freq_Hz'][()]*optics['p0c_eV'][()]/c*optics['length'][()])
        B = 2.*np.pi*optics['rf_freq_Hz'][()]/c
        D = (optics['alfa'][()] - 1./optics['gamma0'][()]**2)/2.

        m = D/2./A*(ptau*1.e-4)**2
        K = scipy.special.ellipk(m)
        Jz = 4*np.sqrt(2*A)/(np.pi*B*np.sqrt(D)) * ( scipy.special.ellipe(m) - ( 1 - m ) * scipy.special.ellipk(m))
        #print(Jz, longitudinal_distribution.h_exp(Jz))
        rho[ptau] = np.exp(- (J1_0 / e1_init + J2_0 / e2_init) ) * longitudinal_distribution.h_exp(Jz)

    turn_list = []
    e1_list = []
    e2_list = []
    e1err_list = []
    e2err_list = []
    for ii, turn in enumerate(np.linspace(0, 10000000, ncurves)):
        for ptau, this_file in zip(ptau_list, files):
            this_turn = this_file[f'turn{int(turn):d}']
            alive = this_turn['at_turn'][()] == turn
            J1 = this_turn['J1'][()][alive] / e1_init
            J2 = this_turn['J2'][()][alive] / e2_init
            xn = this_turn['xn'][()][alive] / np.sqrt(e1_init)
            yn = this_turn['yn'][()][alive] / np.sqrt(e2_init)
            if ptau == 0.00:
                countsJ1, binsJ1 = np.histogram(J1, weights=rho[ptau][alive], bins=bins, range=Jrange)
                countsJ2, binsJ2 = np.histogram(J2, weights=rho[ptau][alive], bins=bins, range=Jrange)
                countsX, binsX = np.histogram(xn, weights=rho[ptau][alive], bins=bins, range=xrange)
                countsY, binsY = np.histogram(yn, weights=rho[ptau][alive], bins=bins, range=xrange)
            else:
                countsJ1_t, binsJ1 = np.histogram(J1, weights=rho[ptau][alive], bins=bins, range=Jrange)
                countsJ2_t, binsJ2 = np.histogram(J2, weights=rho[ptau][alive], bins=bins, range=Jrange)
                countsX_t, binsX   = np.histogram(xn, weights=rho[ptau][alive], bins=bins, range=xrange)
                countsY_t, binsY   = np.histogram(yn, weights=rho[ptau][alive], bins=bins, range=xrange)
                countsJ1 += countsJ1_t
                countsJ2 += countsJ2_t
                countsX += countsX_t
                countsY += countsY_t


        countsJ1 /= np.sum(countsJ1) * (binsJ1[1]-binsJ1[0])
        countsJ2 /= np.sum(countsJ2) * (binsJ2[1]-binsJ2[0])
        countsX /= np.sum(countsX) * (binsX[1]-binsX[0])
        countsY /= np.sum(countsY) * (binsY[1]-binsY[0])

        cbinsJ1 =0.5*(binsJ1[1:]+binsJ1[:-1])
        cbinsJ2 =0.5*(binsJ2[1:]+binsJ2[:-1])
        cbinsX =0.5*(binsX[1:]+binsX[:-1])
        cbinsY =0.5*(binsY[1:]+binsY[:-1])
        fit_range = cbinsJ1 < 3.**2/2.
        fitexp = lambda x,p0,p1: p0*np.exp(-x/p1)
        popt_e1, pcov_e1 = scipy.optimize.curve_fit(fitexp, cbinsJ1[fit_range], countsJ1[fit_range], p0=[countsJ1[0],1.],
                                              sigma=np.sqrt(countsJ1[fit_range]))
        popt_e2, pcov_e2 = scipy.optimize.curve_fit(fitexp, cbinsJ2[fit_range], countsJ2[fit_range], p0=[countsJ2[0],1.],
                                              sigma=np.sqrt(countsJ2[fit_range]))
        if ii%100 == 0:
            print(f'{ii}, e1 = {popt_e1[1]:.3f} +- {np.sqrt(pcov_e1[1,1]):.3f}, e2 = {popt_e2[1]:.3f} +- {np.sqrt(pcov_e2[1,1]):.3f}')
        turn_list.append(turn)
        e1_list.append(popt_e1[1])
        e2_list.append(popt_e2[1])
        e1err_list.append(np.sqrt(pcov_e1[1,1]))
        e2err_list.append(np.sqrt(pcov_e2[1,1]))

        ax1.plot(cbinsJ1, countsJ1, '.-' )
        ax2.plot(cbinsJ2, countsJ2, '.-' , label = f'{int(turn / 1e6):d}' + ' M')
        ax3.plot(cbinsX, countsX, '.-' )
        ax4.plot(cbinsY, countsY, '.-')

        # ax1.hist(J1, weights=rho[alive], bins=bins, range=Jrange, histtype='step', color=mycolors[ii], density=density)
        # ax2.hist(J2, weights=rho[alive], bins=bins, range=Jrange, histtype='step', color=mycolors[ii], density=density,
        #          label='$\\mathbf{'+f'{int(turn/1e6):d}'+'\ M}$')
        # ax3.hist(xn, weights=rho[alive], bins=bins, range=xrange, histtype='step', color=mycolors[ii], density=density)
        # ax4.hist(yn, weights=rho[alive], bins=bins, range=xrange, histtype='step', color=mycolors[ii], density=density)


    ax1.set_xlabel('J_1 [s^2/2]}$')
    ax2.set_xlabel('J_2 [s^2/2]}$')
    ax3.set_xlabel('x [s]}$')
    ax4.set_xlabel('y [s]}$')

    #ax1.set_yscale('log')
    #ax2.set_yscale('log')
    #ax3.set_yscale('log')
    #ax4.set_yscale('log')

    Jp = np.linspace(0, last_sigma**2/2., 1000)
    Jexp = np.exp(-Jp)
    xp = np.linspace(-last_sigma, last_sigma, 1000)
    xgaus = np.exp(-xp**2 / 2.) / np.sqrt(2. * np.pi)

    ax1.plot(Jp,Jexp, 'k')
    ax2.plot(Jp,Jexp, 'k')
    ax3.plot(xp, xgaus, 'k')
    ax4.plot(xp, xgaus, 'k')
    ax2.legend(loc='upper right', fontsize=20, bbox_to_anchor=(1.0, 1.00), bbox_transform=plt.gcf().transFigure,
               title='Turns')
    #ax2.legend()

    #fig1.suptitle('$\\mathbf{'+f'{intensity:.2f}'+'\\cdot 10^{11}\, ppb,\ p_{\\tau}=' + f'{ptau:.2f}' + '\\cdot 10^{-4},\ MB+MQ\ e\\textbf{-}clouds}$')

    turn2hour = 26658.883 / c / 3600.
    p_e1, cov_e1 = np.polyfit(np.array(turn_list)*turn2hour, np.array(e1_list), 1, w=1./np.array(e1err_list), cov='unscaled')
    p_e2, cov_e2 = np.polyfit(np.array(turn_list)*turn2hour, np.array(e2_list), 1, w=1./np.array(e2err_list), cov='unscaled')
    e1_growth = p_e1[0]*2.0
    e1_growth_err = np.sqrt(cov_e1[0,0])*2.0
    e2_growth = p_e2[0]*2.0
    e2_growth_err = np.sqrt(cov_e2[0,0])*2.0

    print(f'Horizontal emittance growth: {e1_growth:.3f} +- {e1_growth_err:.3f} um/h')
    print(f'Vertical emittance growth: {e2_growth:.3f} +- {e2_growth_err:.3f} um/h')

    out_dict = {'turns': turn_list, 'time': np.array(turn_list)*turn2hour,
                'e1': np.array(e1_list), 'e1_err': np.array(e1err_list),
                'e2': np.array(e2_list), 'e2_err': np.array(e2err_list),
                'e1_growth': e1_growth, 'e1_growth_err': e1_growth_err,
                'e2_growth': e2_growth, 'e2_growth_err': e2_growth_err
               }

    pickle.dump(out_dict, open(out_fname, 'wb'))
plt.show()
