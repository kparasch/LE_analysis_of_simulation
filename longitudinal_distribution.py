import numpy as np

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
