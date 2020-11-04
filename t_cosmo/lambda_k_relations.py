import numpy as np
from scipy.special import gamma

n = 0.8
a_k_i = np.loadtxt('fit_values.txt').T

def get_lambda_0_k(k):
    def _get_lambda_k(lambda_0_0):
        a_i = a_k_i[k - 1]
        lam_minus_one_fifth = lambda_0_0**(-1./5.)
        lam_terms = np.array([lam_minus_one_fifth,
                              lam_minus_one_fifth**2,
                              lam_minus_one_fifth**3])
        prefactor = gamma(k + 10./(3. - n))/gamma(10./(3. - n))
        postfactor = np.dot(a_i, lam_terms)
        return prefactor * lambda_0_0 * postfactor
    return _get_lambda_k

def get_lambda_from_mass(source_frame_mass, lambda_0_0, M0=1.4):
    """Based on Yagi & Yunes (2018) expansion (Eq. 22)"""
    lambda_tot = lambda_0_0

    for k, lam_k in enumerate(map(get_lambda_0_k, (1, 2, 3))):
        r = lam_k(lambda_0_0) * (1 - source_frame_mass/M0)**(k + 1)
        lambda_tot += r
    return lambda_tot
