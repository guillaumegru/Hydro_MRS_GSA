"""
This python library is used to compute AMA senisitivity indices for the first two statistical moments.

@author: Guillaume Gru
"""

import numpy as np 

def compute_stat_moments(y):
    '''
    Computes the first 2 statistical moments of the distribution sampled in array y.

    Parameters
    ----------
    y : ndarray (float64)
        Distribution sample.

    Returns
    -------
    mm1: float64
        Empirical expected value.
    mm2: float64
        Empirical variance.
    '''
    
    mm1 = np.mean(y) # empirical mean
    mm2 = np.var(y) # empirical variance

    return mm1, mm2

def Ama_indices(nb_bins, x_par, y, mm1, mm2):

    '''
    Computes AMA sensitivity indices for the first 2 statistical moments.
    Ref: Dell'Oca et al., 2017 - Moment-based metrics for global sensitivity analysis of hydrological systems
         doi : https://doi.org/10.5194/hess-21-6219-2017
    
    Parameters
    ----------
    nb_bins: int
        Number of sub-intervals to consider for the paramter x variation range.
    x_par: ndarray of float64
        Array containing the samples of one model parameter.
        Size : n_simulations
    y: ndarray of float64
        Array containing the model outputs samples.
        Size : n_simulations

    mm1: float64
        y first statistical moment: expected value
    mm2: float64
        y first statistical moment: variance 
    
    Returns
    -------
    amae: float64
        AMAE index for the parameter x

    amav: float64
        AMAV index for the parameter x

    '''

    m1 = np.zeros(nb_bins, dtype=np.float64)
    m2 = np.zeros(nb_bins, dtype=np.float64)
    nb_val = np.zeros(nb_bins, dtype=int)

    MinP = np.min(x_par)
    MaxP = np.max(x_par)
    delta = (MaxP - MinP) / nb_bins * 1.000001

    for i in range(len(y)):
        k = int((x_par[i] - MinP) / delta)
        if k >= nb_bins:
            k = nb_bins - 1
        nb_val[k] += 1
        m1[k] += y[i]
        m2[k] += y[i] * y[i]

    for k in range(nb_bins):
        if nb_val[k] > 0:
            m1[k] /= nb_val[k]
            m2[k] = m2[k] / nb_val[k] - m1[k] * m1[k]
        else:
            m1[k] = 0.0
            m2[k] = 0.0

    amae = 0.0
    amav = 0.0
    for k in range(nb_bins):
        amae += abs(mm1 - m1[k])
        amav += abs(mm2 - m2[k]) / mm2

    amae = amae / nb_bins / mm1
    amav = amav / nb_bins

    return amae, amav


def ishigami(x1, x2, x3, a, b):
    '''
    Computes Ishigami function's output.
    x1, x2, x3 can be arrays.
    a and b are scalars.
    
    Parameters
    ----------
    x1: ndarray of float64
        first dimension.
    x2: ndarray of float64
        second dimension.
    x3: ndarray of float64
        third dimension.

    a: float64
        first parameter of Ishigami function
    b: float64
        second parameter of Ishigami function
    
    Returns
    -------
    y : ndarray of float64
        y(x1,x2,x3) = ishigami(x1,x2,x3)
    '''

    t1 = np.sin(2*np.pi * x1 - np.pi)
    t2 = a * np.sin(2*np.pi * x2 - np.pi)**2
    t3 = b * np.sin(2*np.pi*x1 - np.pi) * (2*np.pi * x3 - np.pi)**4 

    y = t1+t2+t3
    return(y)

    