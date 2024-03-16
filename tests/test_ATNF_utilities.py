#!/usr/bin/env python

"""
test_ATNF_utilities
--------------------------------------

Tests for ATNF utility functions.
"""

import os

from ccg_gwb import CCG_CACHEDIR
from ccg_gwb.simulation.ATNF_utilities import ATNF_ephemeris, parse_ephem, query_ATNF
from tests.ccg_gwb_test_data import datadir


def test_ATNF_ephemeris():
    """Check querying ATNF."""

    CACHE_file = CCG_CACHEDIR + "/P0 < 0.002_ephems.pkl"
    query = query_ATNF(condition="P0 < 0.002")
    msg = "querying ATNF failed"
    assert len(query.get_pulsars()) > 0, msg

    ephems = ATNF_ephemeris(query, cache=False)
    msg = "could not read ephemeris"
    assert len(ephems) > 0, msg
    ephems = ATNF_ephemeris(query, cache=True)
    msg = "cache file not written"
    assert os.path.exists(CACHE_file), msg


def test_parse_ephem():
    """Check parse ephem function"""

    ephem_file = datadir + "/J0023+0923.ephem"
    with open(ephem_file, "r") as file:
        ephem = file.read()
    "\n".join(ephem)
    params = parse_ephem(ephem, quiet=True)
    PSR = [param.value for param in params if param.name == "PSR"][0]
    msg = "incorrect pulsar name"
    assert PSR == "J0023+0923", msg
