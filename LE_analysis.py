from scipy.constants import c
import numpy as np
import h5py

def combine_losses(fname_list, Jz_list, e1=None, e2=None, bunch_length=None, energy='450GeV', q=1., rfvoltage=6):
    if e1 is None:
        e1 = 2.e-6/479.603977 # emittance geometric for 2 micrometers normalized at 450GeV
    if e2 is None:
        e2 = 2.e-6/479.603977 # emittance geometric for 2 micrometers normalized at 450GeV
    if bunch_length is None:
        bunch_length = 0.090 # rms bunch length 0.09m -> full (4sigma) bunch length 1.2ns
    if energy == '450GeV':
        Jufp = 1.456e-4
        J0 = 8.382e-5
    elif energy == '6500GeV':
        Jufp = 4.663e-5
        J0 = 2.684e-5

    L_list = []
    Lerr_list = []

    for fname in fname_list:
        thisfile = h5py.File(fname, 'r')

        sigma1 = 0.
        sigma2 = thisfile['args/n_sigma'][()]

        t0 = thisfile['lost/turn_lost'][()]
        J1 = thisfile['turn0/J1'][()]
        J2 = thisfile['turn0/J2'][()]

        tc, L, Lerr = get_loss_rate_one_set(t0, J1, J2, e1=e1, e2=e2, s1=sigma1, s2=sigma2, q=q)
        L_list.append(L)
        Lerr_list.append(Lerr)


    M = len(Jz_list) + 1
    N = len(tc)

    Jz_array = np.array(Jz_list)

    xx1, yy1, yyerr1 = get_xx_yy_yyerr(L_list, Lerr_list, Jz_array, distribution='parabolic', energy='450GeV', bunch_length=bunch_length, rfvoltage=rfvoltage)
    xx2, yy2, yyerr2 = get_xx_yy_yyerr(L_list, Lerr_list, Jz_array, distribution='q-exponential', energy='450GeV', bunch_length=bunch_length, rfvoltage=rfvoltage)
    xx3, yy3, yyerr3 = get_xx_yy_yyerr(L_list, Lerr_list, Jz_array, distribution='exponential', energy='450GeV', bunch_length=bunch_length, rfvoltage=rfvoltage)
    
    dt = tc[1] - tc[0]
    tot_L1, tot_Lerr1 = integrate_with_linear_interpolation(xx1,yy1,yyerr1,dt)
    tot_L2, tot_Lerr2 = integrate_with_linear_interpolation(xx2,yy2,yyerr2,dt)
    tot_L3, tot_Lerr3 = integrate_with_linear_interpolation(xx3,yy3,yyerr3,dt)

    L_per_h = np.array([tot_L1, tot_L2, tot_L3])
    Lerr_per_h = np.array([tot_Lerr1, tot_Lerr2, tot_Lerr3])

    tot_L = np.mean(L_per_h, axis=0)
    syst_err = np.std(L_per_h, axis=0)
    stat_err = np.sqrt(np.mean(Lerr_per_h**2,axis=0))

    tot_Lerr = np.sqrt(syst_err**2 + stat_err**2)
    
    tot_I, tot_Ierr = intensity_from_losses(tc, tot_L, tot_Lerr)

    turn2sec = 26658.883 / c
    sec2min = 1./60.
    min2hour = 1./60.

    output_dict = { 'tc' : tc*turn2sec*sec2min,
                    'dt' : dt,
                    'Jz_para' : xx1,
                    'L_para' : yy1,
                    'Lerr_para' : yyerr1,
                    'tot_L_para' : tot_L1, 
                    'tot_Lerr_para' : tot_Lerr1, 
                    'Jz_qexp' : xx2,
                    'L_qexp' : yy2,
                    'Lerr_qexp' : yyerr2,
                    'tot_L_qexp' : tot_L2, 
                    'tot_Lerr_qexp' : tot_Lerr2, 
                    'Jz_exp' : xx3,
                    'L_exp' : yy3,
                    'Lerr_exp' : yyerr3,
                    'tot_L_exp' : tot_L3, 
                    'tot_Lerr_exp' : tot_Lerr3, 
                    'tot_L' : tot_L,
                    'tot_Lerr' : tot_Lerr,
                    'tot_I' : tot_I,
                    'tot_Ierr' : tot_Ierr
                  }
    return output_dict

def Jz_last_func(energy='450GeV', rfvoltage=6, distribution='parabolic', bunch_length=0.090):
    if energy != '450GeV':
        raise Exception(f'Error: Energy {energy} not implemented!')

    if distribution == 'parabolic':
        if rfvoltage == 8:
            J0_polynomial = np.polynomial.polynomial.Polynomial([7.387e-5, 1.333e-4, 3.967e-5,
                                                                 -4.751e-5, -1.643e-5, 5.120e-5,
                                                                 4.033e-5], domain=[2.723e-3, 1.374e-1])
        elif rfvoltage == 6:
            J0_polynomial = np.polynomial.polynomial.Polynomial([6.783e-5, 1.233e-4, 3.539e-5,
                                                                 -6.123e-5, -2.693e-5, 7.507e-5,
                                                                 6.144e-5], domain=[2.927e-3, 1.417e-1],
                                                                window=[-1.,1.])
        else:
            raise Exception(f'Error: RF voltage {rfvoltage:d} not implemented!')
        J0 = J0_polynomial(bunch_length)
        return J0
    else:
        if rfvoltage == 8:
            J_ufp = 1.456e-4
        elif rfvoltage == 6:
            J_ufp = 1.261e-4
        else:
            raise Exception(f'Error: RF voltage {rfvoltage:d} not implemented!')
        return J_ufp


def get_xx_yy_yyerr(L_list, Lerr_list, Jz_array, distribution='parabolic', energy='450GeV', bunch_length=0.090, rfvoltage=6):
    if distribution == 'parabolic':
        h_dist = h_parabolic(Jz_array, energy=energy, bunch_length=bunch_length, rfvoltage=rfvoltage)
    elif distribution == 'q-exponential':
        h_dist = h_qexp(Jz_array, energy=energy, bunch_length=bunch_length, rfvoltage=rfvoltage)
    elif distribution == 'exponential':
        h_dist = h_exp(Jz_array, energy=energy, bunch_length=bunch_length, rfvoltage=rfvoltage)

    Jz_last =  Jz_last_func(energy=energy, rfvoltage=rfvoltage, distribution=distribution, bunch_length=bunch_length)

    M = len(Jz_array) + 1
    N = len(L_list[0])

    xx = np.zeros([M])
    yy = np.zeros([M,N])
    yyerr = np.zeros([M,N])

    xx[:-1] = Jz_array
    xx[-1] = Jz_last
    yy[-1,:] = 0.
    yyerr[-1,:] = 0.
    for ii in range(M-1):
        yy[ii,:] = L_list[ii] * h_dist[ii]
        yyerr[ii,:] = Lerr_list[ii] * h_dist[ii]

    return xx, yy, yyerr

def integrate_with_linear_interpolation(xx,yy,yyerr,dt):
    M, N = yy.shape
    tot_L = np.zeros(N)
    tot_VarL = np.zeros(N)
    for ii in range(1, M):
        x0 = xx[ii-1]
        x1 = xx[ii]
        y0 = yy[ii-1,:]
        y1 = yy[ii,:]
        erry0 = yyerr[ii-1,:]
        erry1 = yyerr[ii,:]
        tot_L += ( 0.5 * (x1 - x0) ) * (y1 + y0)
        tot_VarL += ( 0.5 * (x1 - x0) )**2  * (erry1**2 + erry0**2)
    tot_errL = np.sqrt(tot_VarL)

    turn2sec = 26658.883 / c
    sec2min = 1./60
    min2hour = 1./60
    conv = turn2sec*sec2min*min2hour

    tot_L *= (100./dt)/conv
    tot_errL *= (100./dt)/conv

    return tot_L, tot_errL

def get_loss_rate_one_set(t0, J1, J2, e1=None, e2=None, s1=0., s2=5.5, q=1.):

    if e1 is None:
        e1 = 2.e-6/479.603977 # emittance geometric for 2 micrometers normalized at 450GeV
    if e2 is None:
        e2 = 2.e-6/479.603977 # emittance geometric for 2 micrometers normalized at 450GeV

    J1e = J1/e1
    J2e = J2/e2

    mask = np.logical_and(s1**2 / 2. < J1 + J2,
                          s2**2 / 2. > J1 + J2,
                         )

    tc, L, Lerr = losses(t0[mask], J1e[mask], J2e[mask], sigma1=s1, sigma2=s2, q=q)
    #tc = 0.5*(t_i[:-1]+t_i[1:])

    return tc, L, Lerr


def K(t0, t1, t2):
    return 1*np.logical_and(t0 >= t1, t0 < t2)

def losses(t0, J1, J2, sigma1=0., sigma2=5.5, q=1., n_points=25, last_turn=1e7):
    M = t0.shape[0]
    V_1   = 1/2 * np.pi**2 * sigma1**4 /16
    V_2   = 1/2 * np.pi**2 * sigma2**4 /16
    V = V_2 - V_1

    q_corr = (-0.239 * (q-1) + 1)**2
    J1 /= q_corr
    J2 /= q_corr
    g_cst =  4 / np.pi**2 * (2-q)**2 #/ (1 - np.exp(-sigma2**2/2) * (1 + sigma2**2/2))
    if q==1.:
        g_fct  = np.exp(-J1-J2)
    else:
        g_fct = np.power(1 - ( 1 - q) * (J1) ,1./(1-q)) * np.power(1 - ( 1 - q) * (J2) ,1./(1-q))
    t_i = np.linspace(0, last_turn, n_points + 1)
    losses = np.zeros(n_points,)
    losses_squared = np.zeros(n_points,)
    for i in range(n_points):
        losses[i] = (V / M * g_cst) * np.sum( np.multiply( K(t0, t_i[i], t_i[i+1]), g_fct))
        losses_squared[i] = (V**2 / M * g_cst**2) * np.sum( np.multiply( K(t0, t_i[i], t_i[i+1]), g_fct**2))

    variance_losses = (losses_squared - losses**2)/M
    losses_std = np.sqrt(variance_losses)
    tc = 0.5*(t_i[:-1]+t_i[1:])

    return tc, losses, losses_std


def intensity_from_losses(tc,L,Lerr):
    turn2sec = 26658.883/ c
    sec2min = 1./60
    min2hour = 1./60

    T = tc*turn2sec*sec2min*min2hour
    dt = T[1]-T[0]
    A = L*dt/100.
    Aerr = Lerr*dt/100.

    meanI = np.empty(len(A)+1)
    errI = np.empty(len(A)+1)
    VarA = Aerr**2
    meanI[0]=1
    errI[0]=0
    for i in range(len(A)):
        losses_at_t = np.sum(A[:i+1])
        meanI[i+1] = 1 - losses_at_t
        variance_losses_total = VarA
        variance_intensity_at_t = np.sum(variance_losses_total[:i+1])
        errI[i+1] = np.sqrt(variance_intensity_at_t)
    
    return meanI, errI

def h_parabolic(J, bunch_length=0.090, energy='450GeV', rfvoltage=6):
    if energy == '450GeV':
        if rfvoltage == 8:
            Jufp = 1.456e-4
            J0_polynomial = np.polynomial.polynomial.Polynomial([7.387e-5, 1.333e-4, 3.967e-5,
                                                                 -4.751e-5, -1.643e-5, 5.120e-5,
                                                                 4.033e-5], domain=[2.723e-3, 1.374e-1])
        elif rfvoltage == 6:
            Jufp = 1.261e-4
            J0_polynomial = np.polynomial.polynomial.Polynomial([6.783e-5, 1.233e-4, 3.539e-5,
                                                                 -6.123e-5, -2.693e-5, 7.507e-5,
                                                                 6.144e-5], domain=[2.927e-3, 1.417e-1],
                                                                window=[-1.,1.])
        J0 = J0_polynomial(bunch_length)
        #J0 = 8.382e-5
    elif energy == '6500GeV':
        Jufp = 4.663e-5
        J0 = 2.684e-5

    hj = 3./J0*(1.-J/J0)**2*(1-np.heaviside(J/J0-1,0))
    
    return hj

def h_qexp(J, bunch_length=0.090, energy='450GeV', rfvoltage=6):
    if energy == '450GeV':
        if rfvoltage == 8:
            Jufp = 1.456e-4
            #J0=2.694e-5
            J0_polynomial = np.polynomial.polynomial.Polynomial([2.180e-5, 3.808e-5, 1.146e-5,
                                                                 -4.729e-6, 5.966e-6, 1.134e-5,
                                                                 5.538e-6], domain=[2.702e-3, 1.312e-1],
                                                                window=[-1., 1.])
        elif rfvoltage == 6:
            Jufp = 1.261e-4
            #J0=2.694e-5
            J0_polynomial = np.polynomial.polynomial.Polynomial([1.979e-5, 3.523e-5, 1.176e-5,
                                                                 -8.317e-6, -2.184e-7, 1.701e-5,
                                                                 1.278e-5], domain=[2.872e-3, 1.350e-1],
                                                                window=[-1., 1.])
        J0 = J0_polynomial(bunch_length)
    elif energy == '6500GeV':
        Jufp = 4.663e-5
        J0=8.626e-6

    q=0.85
    qp = 1./(2-q)
    const = 1./(1 - np.power(1 - (1 - qp)*Jufp/J0, 1./(1-qp)))
    hj= const*(2-q)/J0 * np.power( 1 - (1 - q)*J/J0, 1./(1-q))
 
    return hj

def h_exp(J, bunch_length=0.090, energy='450GeV', rfvoltage=6):
    if energy == '450GeV':
        if rfvoltage == 8:
            Jufp = 1.456e-4
            #J0=2.039e-5
            J0_polynomial = np.polynomial.polynomial.Polynomial([1.543e-5, 2.703e-5, 8.535e-6,
                                                                 7.399e-7, 9.595e-6, 5.617e-6,
                                                                 -1.883e-7], domain=[2.659e-3, 1.260e-1],
                                                                window=[-1., 1.])
        elif rfvoltage == 6:
            Jufp = 1.261e-4
            #J0=2.039e-5
            J0_polynomial = np.polynomial.polynomial.Polynomial([1.438e-5, 2.532e-5, 9.490e-6,
                                                                 7.264e-7, 6.577e-6, 7.994e-6,
                                                                 3.636e-6], domain=[2.849e-3, 1.313e-1],
                                                                window=[-1., 1.])
        J0 = J0_polynomial(bunch_length)
    elif energy == '6500GeV':
        Jufp = 4.663e-5
        J0=6.528e-6

    const = 1./(1 - np.exp(-Jufp/J0))
    hj = const/J0 * np.exp(-J/J0)

    return hj
