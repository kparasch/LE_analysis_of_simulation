import numpy as np
import h5py
import scipy
from scipy.constants import c
import scipy.special
import scipy.optimize
import LE_analysis
import pickle

import matplotlib.pyplot as plt


def integrate(xx, yy, yyerr):
    M, N = yy.shape
    tot_L = np.zeros(N)
    tot_VarL = np.zeros(N)
    for ii in range(1, M):
        x0 = xx[ii - 1]
        x1 = xx[ii]
        y0 = yy[ii - 1, :]
        y1 = yy[ii, :]
        erry0 = yyerr[ii - 1, :]
        erry1 = yyerr[ii, :]
        tot_L += (0.5 * (x1 - x0)) * (y1 + y0)
        tot_VarL += (0.5 * (x1 - x0)) ** 2 * (erry1 ** 2 + erry0 ** 2)
    tot_errL = np.sqrt(tot_VarL)
    return tot_L, tot_errL

def rejection_sampling(f, n, xmin, xmax):
    distr = np.zeros(n)
    for ii in range(n):
        x = np.random.random() * (xmax - xmin) + xmin
        u = np.random.random()
        while u > f(x):
            x = np.random.random() * (xmax - xmin) + xmin
            u = np.random.random()
        distr[ii] = x
    return distr


mode='MBMQ'
intensity=1.20
for mode in ['MBMBC']:
    for MBc in [0.1, 0.5]:
        energy='450GeV'
        q=1.2
        rfvoltage=6.
        # emit = 2.0
        print(f'Mode: {mode}, MBc = {MBc:.2f}, Intensity = {intensity:.2f}e11ppb, q = {q:.2f}')
        emit_list = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        bl_growth_list=[]
        bl_growth_err_list=[]
        for emit in emit_list:
            e1 = emit*1.e-6 / 479.603977  # emittance geometric for 2 micrometers normalized at 450GeV
            e2 = emit*1.e-6 / 479.603977  # emittance geometric for 2 micrometers normalized at 450GeV
            bunch_length = 0.090  # rms bunch length 0.09m -> full (4sigma) bunch length 1.2ns

            in_folder = '/eos/project/e/ecloud-simulations/kparasch/Norm_tracking_data/'
            ptaus = ['0.00', '2.94', '4.13', '5.01', '5.73', '6.34', '6.87', '7.33', '7.73', '8.06']
            if mode == 'MB':
                fname_list = [
                    f'LE_LHC_450GeV_DATAseyMB_{intensity:.2f}e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau}e-4ptau_norm.h5'
                    for ptau in ptaus]
            elif mode == 'MBMQ':
                fname_list = [
                    f'LE_LHC_450GeV_DATAseyMB_DATAseyMQ_{intensity:.2f}e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau}e-4ptau_norm.h5'
                    for ptau in ptaus]
            elif mode == 'MBMBC':
                fname_list = [
                    f'LE_LHC_450GeV_DATAseyMB_{MBc:.2f}MBc_{intensity:.2f}e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau}e-4ptau_norm.h5'
                    for ptau in ptaus]
            elif mode == 'NOEC':
                intensity = 0.00
                fname_list = [f'LE_LHC_450GeV_0.00e11ppb_62.270qx_60.295qy_15qprime_40IMO_6VRF_{ptau}e-4ptau_norm.h5' for ptau in
                            ptaus]

            ptaus_f = [float(ptau) * 1.e-4 for ptau in ptaus]
            thisfile = h5py.File(in_folder + fname_list[0], 'r')
            optics = thisfile['optics']
            A = optics['rf_volt_V'][()] / (2. * np.pi * optics['rf_freq_Hz'][()] * optics['p0c_eV'][()] / c * optics['length'][()])
            B = 2. * np.pi * optics['rf_freq_Hz'][()] / c
            D = (optics['alfa'][()] - 1. / optics['gamma0'][()] ** 2) / 2.
            m = D / 2. / A * (np.array(ptaus_f)) ** 2
            K = scipy.special.ellipk(m)
            Jz_list = 4 * np.sqrt(2 * A) / (np.pi * B * np.sqrt(D)) * (scipy.special.ellipe(m) - (1 - m) * scipy.special.ellipk(m))

            L_list = []
            Lerr_list = []

            for fname in fname_list:
                thisfile = h5py.File(in_folder + fname, 'r')

                sigma1 = 0.
                sigma2 = thisfile['args/n_sigma'][()]

                t0 = thisfile['lost/turn_lost'][()]
                J1 = thisfile['turn0/J1'][()]
                J2 = thisfile['turn0/J2'][()]

                tc, L, Lerr = LE_analysis.get_loss_rate_one_set(t0, J1, J2, e1=e1, e2=e2, s1=sigma1, s2=sigma2, q=q)
                L_list.append(L)
                Lerr_list.append(Lerr)

            M = len(Jz_list) + 1
            N = len(tc)

            Jz_array = np.array(Jz_list)
            losses_dict = LE_analysis.combine_losses(list(map(lambda x: in_folder + x, fname_list)), Jz_list, e1=e1, e2=e2, bunch_length=bunch_length, q=q, rfvoltage=rfvoltage)
            npoints=25
            xx1, yy1, yyerr1 = LE_analysis.get_xx_yy_yyerr(L_list, Lerr_list, Jz_array, distribution='parabolic',
                                               bunch_length=bunch_length, rfvoltage=rfvoltage)

            Jz = losses_dict['Jz_exp']
            turn2sec = 26658.883 / c
            sec2min = 1. / 60
            min2hour = 1. / 60
            conv = turn2sec * sec2min * min2hour
            dt = tc[1] - tc[0]
            dJz = np.mean(np.diff(Jz))
            conv_factor = (100. / dt) / conv # * dJz


            Hdist = LE_analysis.h_exp(Jz, bunch_length=bunch_length)

            L_by_Jz = losses_dict['L_exp']
            Lerr_by_Jz = losses_dict['Lerr_exp']
            mask = Hdist > 0
            for ii in range(losses_dict['L_exp'].shape[1]):
                L_by_Jz[:, ii][mask] /= Hdist[mask]
                Lerr_by_Jz[:, ii][mask] /= Hdist[mask]


            I_by_Jz = np.empty([L_by_Jz.shape[0], L_by_Jz.shape[1] + 1])
            Ierr_by_Jz = np.empty([L_by_Jz.shape[0], L_by_Jz.shape[1] + 1])
            for nn in range(I_by_Jz.shape[0]):
                I_by_Jz[nn, :], Ierr_by_Jz[nn, :] = LE_analysis.intensity_from_losses(tc, L_by_Jz[nn, :] * conv_factor, Lerr_by_Jz[nn, :] * conv_factor)
            #    I1[:, nn], Ierr1[:, nn] = LE_analysis.intensity_from_losses(tc, L_list[nn] * conv_factor, Lerr_list[nn] * conv_factor)
            # #print(xx1)
            # print(I1, I1.shape, tc.shape)
            # print(xx1.shape, yy1.shape, yyerr1.shape)

            plt.close('all')
            fig1 = plt.figure(2, figsize=[10,9/1.6*2])
            ax1 = fig1.add_subplot(211)
            tstart = losses_dict['tc'][0]
            dt2 = (losses_dict['tc'][1] - tstart)/2.
            tend = losses_dict['tc'][-1]
            N_I = len(losses_dict['tot_I'])
            te = np.linspace(tstart - dt2, tend + dt2, N_I)
            ax1.plot(te, losses_dict['tot_I'], '.-')
            ax1.set_ylabel('$\\mathbf{Relative\ intensity}$')
            ax1.set_xlabel('$\\mathbf{Time\ [min]}$')

            ax2 = fig1.add_subplot(212)
            ax2.plot(losses_dict['tc'], losses_dict['tot_L'], 'k.-')
            ax2.plot(losses_dict['tc'], losses_dict['tot_L_para'], 'r.-')
            ax2.plot(losses_dict['tc'], losses_dict['tot_L_exp'], 'g.-')
            ax2.plot(losses_dict['tc'], losses_dict['tot_L_qexp'], 'b.-')
            ax2.set_ylabel('$\\mathbf{Loss\ rate\ [\%/h]}$')
            ax2.set_xlabel('$\\mathbf{Time\ [min]}$')

            fig2 = plt.figure(3, figsize=[10,9/1.6*2])
            ax3 = fig2.add_subplot(211)
            for ii in range(I_by_Jz.shape[1]):
                ax3.plot(Jz*1.e5, I_by_Jz[:, ii])
            ax3.set_xlabel('$\\mathbf{J_z\ [10^{-5}]}$')

            for ii in range(I_by_Jz.shape[1]):
                I_by_Jz[:, ii] *= Hdist
                Ierr_by_Jz[:, ii] *= Hdist

            tot_I, tot_Ierr = integrate(Jz, I_by_Jz, Ierr_by_Jz)
            ax4 = fig2.add_subplot(212)
            for ii in range(I_by_Jz.shape[1]):
                ax4.plot(Jz*1.e5, I_by_Jz[:, ii]/1.e5)
            ax4.set_xlabel('$\\mathbf{J_z\ [10^{-5}]}$')
            ax4.set_ylabel('$\\mathbf{I\ [10^{5}/J_z]}$')
            ax4.set_yscale('log')
            # ax4.plot(tot_I[:-1]/tot_I[0])
            ax1.plot(te[:-1], tot_I[:-1]/tot_I[0], 'r.-')


            action = lambda x: 4*np.sqrt(2*A)/(np.pi*B*np.sqrt(D)) * ( scipy.special.ellipe(x) - ( 1 - x ) * scipy.special.ellipk(x))
            action_sep = 4*np.sqrt(2*A)/(np.pi*B*np.sqrt(D)) *  scipy.special.ellipe(1)
            NN = 100000
            theta = np.random.uniform(size=NN)*2.*np.pi

            tau_rms = []
            tau_rms_err = []
            for jj in range(I_by_Jz.shape[1]):
                yy = I_by_Jz[:, jj]/I_by_Jz[0, jj]
                #np.interp(x, Jz, yy)
                action_expo = rejection_sampling(lambda x: np.interp(x, Jz, yy), NN, 0, Jz[-1])
                tau_expo = np.empty(NN)
                ptau_expo = np.empty(NN)
                for ii,J in enumerate(action_expo):
                    m = scipy.optimize.fsolve( lambda y: (action(y) - J), J/action_sep)
                    K = scipy.special.ellipk(m)
                    G = 2.*K/np.pi
                    sn, cn, dn, ph = scipy.special.ellipj(G*theta[ii],m)
                    tau_expo[ii] = 2./B*np.arcsin(np.sqrt(m)*sn)
                    ptau_expo[ii] = np.sqrt(2*A*m/D)*cn

                #print(f'{jj}, {np.std(tau_expo):.4f} +- {np.std(tau_expo)/np.sqrt(2.*NN):.4f}')
                tau_rms.append(np.std(tau_expo))
                tau_rms_err.append(np.std(tau_expo)/np.sqrt(2.*NN))


            tau_rms = np.array(tau_rms)
            tau_rms_err = np.array(tau_rms_err)
            par, cov = np.polyfit(np.array(te)/60., tau_rms, 1, w=1./tau_rms_err, cov='unscaled')
            bl_growth = par[0]
            bl_growth_err = np.sqrt(cov[0,0])
            plt.figure(5)
            plt.errorbar(te, tau_rms, tau_rms_err)
            plt.plot(te, par[0]*te/60 + par[1],'r')

            print(f'Emit: {emit}, Bunch length growth: {bl_growth:.5f} +- {bl_growth_err:.5f} m/h')
            bl_growth_list.append(bl_growth)
            bl_growth_err_list.append(bl_growth_err)
        pickle.dump({'emittance': np.array(emit_list),
                     'bl_growth': np.array(bl_growth_list),
                     'bl_growth_err': np.array(bl_growth_err_list)},
                    open(f'data/Bunch_length_DATA{mode}_{MBc:.2f}MBc_{intensity:.2f}e11ppb.pkl', 'wb')
                    #open(f'data/Bunch_length_DATA{mode}_{intensity:.2f}e11ppb.pkl', 'wb')
                    )




plt.show()