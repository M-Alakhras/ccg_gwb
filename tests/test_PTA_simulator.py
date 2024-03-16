#!/usr/bin/env python

"""
test_PTA_simulator
--------------------------------------

Tests for PTA simulation process.
"""

import glob
import os

from ccg_gwb.simulation.PTA_Simulator import PTA_Simulator, list_all_simulations, load_Simulator
from tests.ccg_gwb_test_data import datadir


def test_PTA_Simulator():
    Sim = PTA_Simulator("test_simulation")
    assert Sim.outdir == os.getcwd()
    Sim.ATNF = True
    assert Sim.ATNF is True
    Sim.ATNF_Condition = "P0 < 0.002"
    assert Sim.ATNF_Condition == "P0 < 0.002"
    Sim.outdir = datadir + "/test_simulation"
    assert Sim.outdir == datadir + "/test_simulation"
    assert Sim.TimingModel_Simulator.outdir == datadir + "/test_simulation/par"
    assert Sim.TOAs_Simulator.pardir == datadir + "/test_simulation/par"
    assert Sim.TOAs_Simulator.outdir == datadir + "/test_simulation/toas"
    Sim.signals = ["gw"]
    assert Sim.signals == ["wn", "gw"]
    Sim.signals = ["gw", "wn"]
    assert Sim.signals == ["wn", "gw"]
    Sim.signals = ["wn"]
    assert Sim.signals == ["wn"]
    Sim.start(restart=True)
    par_files = glob.glob(Sim.TimingModel_Simulator.outdir + "/*.par")
    msg = "Timing models not simulated"
    assert len(par_files) > 0, msg
    # toas_files = glob.glob(Sim.TOAs_Simulator.outdir+"/*.toas")
    list_all_simulations()
    Sim_loaded = load_Simulator("badname__simulation_&123$67!!")
    Sim_loaded = load_Simulator("test_simulation")
    msg = "Simulator not loaded"
    assert Sim_loaded.name == "test_simulation", msg

    # cleaning
    for parfile in par_files:
        os.remove(parfile)
    os.rmdir(Sim.TimingModel_Simulator.outdir)
    os.rmdir(Sim.TOAs_Simulator.outdir)
    os.rmdir(Sim.outdir)
