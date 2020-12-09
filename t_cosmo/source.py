from astropy import (
    constants as const,
    cosmology,
    units as u
)
import bilby

from .lambda_k_relations import get_lambda_from_mass


def lambda_0_gw170817_lal_binary_neutron_star(
        frequency_array, chirp_mass, mass_ratio, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, lambda_0_0,
        **kwargs):
    """Use lambda-m relations for lambda_0_0 measurement of GW170817

    Parameters
    ----------
    parameters: dict
        parameter dictionary typically passed to a bilby
        conversion function

    Notes
    -----
    Redshift is fixed at a z = 0.0099, following
    `LIGO DCC-P1800115 <https://dcc.ligo.org/public/0152/P1800115/012/README.md>`_
    """
    REDSHIFT = 0.0099

    total_mass = bilby.gw.conversion.chirp_mass_and_mass_ratio_to_total_mass(
        chirp_mass, mass_ratio
    )

    mass_1 = total_mass / (1 + mass_ratio)
    mass_2 = total_mass - mass_1
    # FIXME: Hard code m_0 value
    M0 = 1.4

    mass_1_source = mass_1 / (1 + REDSHIFT)
    mass_2_source = mass_2 / (1 + REDSHIFT)

    lambda_1 = get_lambda_from_mass(
        mass_1_source, lambda_0_0, M0=M0
    )
    lambda_2 = get_lambda_from_mass(
        mass_2_source, lambda_0_0, M0=M0
    )
    # set lambda to 0 if negative
    lambda_1 = 0 if lambda_1 < 0 else lambda_1
    lambda_2 = 0 if lambda_2 < 0 else lambda_2

    return bilby.gw.source.lal_binary_neutron_star(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, lambda_1, lambda_2,
        **kwargs
    )


def lambda_0_h0_lal_binary_neutron_star(
        frequency_array, chirp_mass, mass_ratio, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, lambda_0_0, h0,
        **kwargs):
    """Convert Hubble constant and luminosity distance
    to component tidal deformabilities using lambda-m
    relations.

    Parameters
    ----------
    parameters: dict
        parameter dictionary typically passed to a bilby
        conversion function

    Notes
    -----
    FlatLambdaCDM is used with provided `h0` and `Om0=0.3` to obtain
    redshift.
    """
    hubble_constant = h0 * u.km / u.s / u.Mpc

    cosmo = cosmology.FlatLambdaCDM(H0=hubble_constant, Om0=0.3)

    try:
        REDSHIFT = cosmology.z_at_value(
            cosmo.luminosity_distance,
            luminosity_distance * u.Mpc    
        )
    except cosmology.func.CosmologyError:
        # Raised when unphysically small luminosity distance
        # is provided. Redshift is considered to be zero.
        REDSHIFT = 0.

    total_mass = bilby.gw.conversion.chirp_mass_and_mass_ratio_to_total_mass(
        chirp_mass, mass_ratio
    )

    mass_1 = total_mass / (1 + mass_ratio)
    mass_2 = total_mass - mass_1
    # FIXME: Hard code m_0 value
    M0 = 1.4

    mass_1_source = mass_1 / (1 + REDSHIFT)
    mass_2_source = mass_2 / (1 + REDSHIFT)

    lambda_1 = get_lambda_from_mass(
        mass_1_source, lambda_0_0, M0=M0
    )
    lambda_2 = get_lambda_from_mass(
        mass_2_source, lambda_0_0, M0=M0
    )
    # set lambda to 0 if negative
    lambda_1 = 0 if lambda_1 < 0 else lambda_1
    lambda_2 = 0 if lambda_2 < 0 else lambda_2

    return bilby.gw.source.lal_binary_neutron_star(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, lambda_1, lambda_2,
        **kwargs
    )


def lambda_0_z_lal_binary_neutron_star(
        frequency_array, chirp_mass, mass_ratio, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, lambda_0_0, z,
        **kwargs):
    """Source function to sample redshift directly"""
    total_mass = bilby.gw.conversion.chirp_mass_and_mass_ratio_to_total_mass(
        chirp_mass, mass_ratio
    )

    mass_1 = total_mass / (1 + mass_ratio)
    mass_2 = total_mass - mass_1
    # FIXME: Hard code m_0 value
    M0 = 1.4

    mass_1_source = mass_1 / (1 + z)
    mass_2_source = mass_2 / (1 + z)

    lambda_1 = get_lambda_from_mass(
        mass_1_source, lambda_0_0, M0=M0
    )
    lambda_2 = get_lambda_from_mass(
        mass_2_source, lambda_0_0, M0=M0
    )

    return bilby.gw.source.lal_binary_neutron_star(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, lambda_1, lambda_2,
        **kwargs
    )
