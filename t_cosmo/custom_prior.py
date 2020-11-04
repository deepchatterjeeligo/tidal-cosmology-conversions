import h5py

import numpy as np
import scipy as sp

from bilby import prior


class CustomMassPrior(prior.Prior):    
    def __init__(self, *args, **kwargs):
        """Prior using on masses.
        
        Notes
        -----
        First two arguments are `scipy.interpolate.interp1d`
        instances for the PDF and CDF of the mass posterior.
        samples.
        """
        interpolant_pdf, interpolant_cdf, *args = args

        kwargs['minimum'] = interpolant_pdf.x.min()
        kwargs['maximum'] = interpolant_pdf.x.max()

        self.interpolant_pdf = interpolant_pdf
        self.interpolant_cdf = interpolant_cdf

        super().__init__(*args, **kwargs)

    def prob(self, val):
        return self.interpolant_pdf(val)

    def ln_prob(self, val):
        return np.log(self.prob(val))

    def rescale(self, val):
        self.test_valid_for_rescaling(val)
        if not hasattr(val, 'size') or val.size==1:
            return self._rescale(val)
        else:
            return np.array(
                [self._rescale(v) for v in val]
            )
    def _rescale(self, val):
        try:
            res = sp.optimize.brentq(
                lambda z: self.interpolant_cdf(z) - val,
                self.minimum, self.maximum
            )
        except ValueError as e:  # pathological issue at the boundaries
            if e.__str__() != 'f(a) and f(b) must have different signs':
                raise
            res = np.random.uniform(self.minimum, self.maximum)
        return res


class ZSquaredPrior(prior.Prior):
    def prob(self, z):
        norm = (1./3.) * (self.maximum**3 - self.minimum**3)
        res = z**2 / norm
        if not hasattr(z, 'size') or z.size==1:
            return res if self.minimum <= z <= self.maximum else 0.
        res[z <= self.minimum] = 0.
        res[z >= self.maximum] = 0.
        return res

    def rescale(self, p):
        self.test_valid_for_rescaling(p)
        norm = self.maximum**3 - self.minimum**3
        res = np.cbrt(p * norm + self.minimum**3)
        return res
