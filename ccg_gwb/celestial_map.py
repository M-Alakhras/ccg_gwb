# celestial_map.py
"""
functions to investigate pulsars distribution across the sky.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates.angles import Angle
from pint import pint_units


def _pulsars_sky_coordinates(psrs):
    ra_rad = np.array([Angle(psr._raj * pint_units["H:M:S"]).to_value(u.rad) for psr in psrs])
    dec_rad = np.array([Angle(psr._decj * pint_units["D:M:S"]).to_value(u.rad) for psr in psrs])
    sky_coordinates = SkyCoord(ra=ra_rad * u.rad, dec=dec_rad * u.rad, frame="gcrs")
    return sky_coordinates


def _pulsars_radec(psrs):
    sky_coordinates = _pulsars_sky_coordinates(psrs)
    ra = sky_coordinates.ra.wrap_at(180.0 * u.degree).rad
    dec = sky_coordinates.dec.rad
    return ra, dec


def _pulsars_separations(psrs):
    npsrs = len(psrs)
    sky_coordinates = _pulsars_sky_coordinates(psrs)
    sep = np.zeros((npsrs, npsrs))
    for i in range(npsrs):
        sep[i, :] = sky_coordinates[i].separation(sky_coordinates)
    return sep
