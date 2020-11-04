from astropy import constants as const, cosmology, units as u
import bilby
from gwpy.timeseries import TimeSeries
import numpy as np
import scipy as sp


def hubble_constant_to_lambda(parameters):
    """Convert Hubble constant and luminosity distance
    to component tidal deformabilities using lambda-m
    relations

    Parameters
    ----------
    parameters: dict
        parameter dictionary typically passed to a bilby
        conversion function

    Notes
    -----
    `parameters` is assumed to constain following keys:
    `h0`, `luminosity_distance`, `lambda_0_0`. FlatLambdaCDM
    is used with provided `h0` and `Om0=0.3` to obtain redshift.
    """
    converted_parameters = parameters.copy()
    added_keys = list()

    hubble_constant = parameters['h0'] * u.km / u.s / u.Mpc
    d_l = parameters['luminosity_distance'] * u.Mpc

    cosmo = cosmology.FlatLambdaCDM(H0=hubble_constant, Om0=0.3)

    try:
        REDSHIFT = cosmology.z_at_value(cosmo.luminosity_distance, d_l)
    except cosmology.func.CosmologyError:
        # Raised when unphysically small luminosity distance
        # is provided. Redshift is considered to be zero.
        REDSHIFT = 0.

    converted_parameters['z'] = REDSHIFT
    added_keys += ['z']

    added_keys += ['mass_1', 'mass_2']
    total_mass = bilby.gw.conversion.chirp_mass_and_mass_ratio_to_total_mass(
        parameters['chirp_mass'], parameters['mass_ratio']
    )

    converted_parameters['mass_1'] = total_mass / (1 + converted_parameters['mass_ratio'])
    converted_parameters['mass_2'] = total_mass - converted_parameters['mass_1']
    # FIXME: Hard code m_0 value
    M0 = 1.4

    mass_1_source = converted_parameters['mass_1'] / (1 + REDSHIFT)
    mass_2_source = converted_parameters['mass_2'] / (1 + REDSHIFT)

    lambda_0_0 = parameters['lambda_0_0']
    lambda_1 = get_lambda_from_mass(
        mass_1_source, lambda_0_0, M0=M0
    )
    lambda_2 = get_lambda_from_mass(
        mass_2_source, lambda_0_0, M0=M0
    )

    converted_parameters['lambda_1'] = 0 if lambda_1 < 0 else lambda_1
    converted_parameters['lambda_2'] = 0 if lambda_2 < 0 else lambda_2

    added_keys += ['lambda_1', 'lambda_2']
    return converted_parameters, added_keys
