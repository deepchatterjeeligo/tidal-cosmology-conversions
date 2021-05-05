import io
import os
from urllib import parse, request

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


# create a custom mass prior for GW170817 masses
_hdf5_url = 'https://dcc.ligo.org/public/0157/P1800370/005/GW170817_GWTC-1.hdf5'
try:
    fp = open(os.path.basename(parse.urlparse(_hdf5_url).path), "rb")
except FileNotFoundError:
    fp = io.BytesIO(request.urlopen(_hdf5_url).read())
_posterior_sample_file = h5py.File(fp, 'r')

_mass1_detector_frame = _posterior_sample_file[
    'IMRPhenomPv2NRT_lowSpin_posterior']['m1_detector_frame_Msun']
_mass2_detector_frame = _posterior_sample_file[
    'IMRPhenomPv2NRT_lowSpin_posterior']['m2_detector_frame_Msun']

# FIXME bins are hardcoded
_mm, _bb = np.histogram(_mass1_detector_frame, bins=20)
_freq_mass_1 = _mm/np.sum(_mm)
_bin_centers_freq_mass_1 = 0.5*(_bb[1:] + _bb[:-1])

_mm, _bb = np.histogram(_mass1_detector_frame, bins=100)
_cum_freq_mass_1 = np.cumsum(_mm)/np.sum(_mm)
_bin_centers_cum_freq_mass_1 = 0.5*(_bb[1:] + _bb[:-1])
# detector frame mass_1 and mass_2 interpolants
_m1_freq_interpolant = sp.interpolate.interp1d(
    _bin_centers_freq_mass_1, _freq_mass_1,
    fill_value=(0.,  0.), bounds_error=False
)
_m1_cum_freq_interpolant = sp.interpolate.interp1d(
    _bin_centers_cum_freq_mass_1, _cum_freq_mass_1,
    fill_value=(.0,  1.), bounds_error=False
)
_mm, _bb = np.histogram(_mass2_detector_frame, bins=20)
_freq_mass_2 = _mm/np.sum(_mm)
_bin_centers_freq_mass_2 = 0.5*(_bb[1:] + _bb[:-1])

_mm, _bb = np.histogram(_mass2_detector_frame, bins=100)
_cum_freq_mass_2 = np.cumsum(_mm)/np.sum(_mm)
_bin_centers_cum_freq_mass_2 = 0.5*(_bb[1:] + _bb[:-1])

_m2_freq_interpolant = sp.interpolate.interp1d(
    _bin_centers_freq_mass_2, _freq_mass_2,
    fill_value=(0.,  0.), bounds_error=False
)
_m2_cum_freq_interpolant = sp.interpolate.interp1d(
    _bin_centers_cum_freq_mass_2, _cum_freq_mass_2,
    fill_value=(0.,  1.), bounds_error=False
)

mass_1_prior_gw170817 = CustomMassPrior(
    _m1_freq_interpolant,
    _m1_cum_freq_interpolant
)
"""Primary mass prior using GW170817 posterior samples."""
mass_2_prior_gw170817 = CustomMassPrior(
    _m2_freq_interpolant,
    _m2_cum_freq_interpolant
)
"""Secondary mass prior using GW170817 posterior samples."""
