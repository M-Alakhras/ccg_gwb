# PTA_Simulator.py
"""
Define PTA Simulator class.
"""

import datetime
import glob
import os
import pickle

import numpy as np
from enterprise.pulsar import Pulsar
from pint.models import get_model
from tqdm.auto import tqdm

from ccg_gwb import CCG_CACHEDIR
from ccg_gwb.simulation.timing_model_simulation import TimingModel_Simulator
from ccg_gwb.simulation.toas_simulation import TOAs_Simulator


class PTA_Simulator(object):

    def __init__(
        self,
        name,
        npsrs=20,
        nrealizations=200,
        signals=None,
        outdir=None,
        ATNF=True,
        ATNF_Condition="P0 < 0.03",
        quiet=True,
    ):

        if outdir is None:
            outdir = os.getcwd()

        # setting Simulator parameters
        self._name = name
        self._npsrs = npsrs
        self._nrealizations = nrealizations
        self._signals = signals
        self.validate_signals()
        self._outdir = os.path.abspath(outdir)
        self._psrdir = self._outdir + "/psr"
        self._status = "init"
        self._current_signal = self.signals[0]
        self._begin = "Not started yet."
        self._finish = "Not finished yet."
        self.quiet = quiet
        self._psrs = None
        self._exclude = []

        self.TimingModel_Simulator = TimingModel_Simulator(
            ATNF=ATNF, ATNF_Condition=ATNF_Condition, outdir=self.outdir + "/par", quiet=self.quiet
        )
        self.TOAs_Simulator = TOAs_Simulator(pardir=self.outdir + "/par", outdir=self.outdir + "/toas")

        self.save()

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        print("Warning:: Simulation name can't be changed.")

    @property
    def npsrs(self):
        return self._npsrs

    @npsrs.setter
    def npsrs(self, value):
        self._npsrs = value

    @property
    def nrealizations(self):
        return self._nrealizations

    @nrealizations.setter
    def nrealizations(self, value):
        self._nrealizations = value

    @property
    def signals(self):
        return self._signals

    @signals.setter
    def signals(self, value):
        self._signals = value
        self.validate_signals()

    @property
    def outdir(self):
        return self._outdir

    @outdir.setter
    def outdir(self, value):
        full_path = os.path.abspath(value)
        self._outdir = full_path
        self._psrdir = full_path + "/psr"
        self.TimingModel_Simulator.outdir = full_path + "/par"
        self.TOAs_Simulator.pardir = full_path + "/par"
        self.TOAs_Simulator.outdir = full_path + "/toas"

    @property
    def psrdir(self):
        return self._psrdir

    @psrdir.setter
    def psrdir(self, value):
        print(f"Warning:: psr directory are set automatically to: {self._psrdir}")

    @property
    def ATNF(self):
        return self.TimingModel_Simulator.ATNF

    @ATNF.setter
    def ATNF(self, value):
        self.TimingModel_Simulator.ATNF = value

    @property
    def ATNF_Condition(self):
        return self.TimingModel_Simulator.ATNF_Condition

    @ATNF_Condition.setter
    def ATNF_Condition(self, value):
        self.TimingModel_Simulator.ATNF_Condition = value

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, value):
        print("Warning:: Simulator status is read-only information.")

    @property
    def current_signal(self):
        return self._current_signal

    @current_signal.setter
    def current_signal(self, value):
        print("Warning:: Simulator current signal is read-only information.")

    @property
    def psrs(self):
        return self._psrs

    @psrs.setter
    def psrs(self, value):
        print("Warning:: PTA pulsars are loaded automatically.")

    @psrs.deleter
    def psrs(self):
        del self._psrs
        self._psrs = None

    @property
    def exclude(self):
        return self._exclude

    @exclude.setter
    def exclude(self, value):
        self._exclude = value

    def start(self, restart=False):
        # initialization
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        if not os.path.exists(self.psrdir):
            os.mkdir(self.psrdir)
        if restart:
            self._begin = "Not started yet."
            self._finish = "Not finished yet."
            self._status = "init"

        if self._begin == "Not started yet.":
            self._begin = datetime.datetime.now()
        self.summary()
        print("Start simulating...")

        if self.status == "init":
            print("Simulating timing models...")
            self.TimingModel_Simulator.start()
            print(f'Timing models have been simulated and saved into: "{self.TimingModel_Simulator.outdir}"')
            self._status = "pars"
            self.save()

        if self.status == "pars":
            print("Simulating times of arrivals...")
            self.TOAs_Simulator.start()
            print(f'Times of arrivals have been simulated and saved into: "{self.TOAs_Simulator.outdir}"')
            self._status = "toas"
            self.save()

        if self.status == "toas":
            print("Creating pulsar objects...")
            self.create_pulsars()
            self._status = "psrs"
            self.save()

        if self.status == "psrs":
            # loop on signals
            current_index = self.signals.index(self.current_signal)
            for signal in self.signals[current_index:]:
                print(f"Simulating PTA realizations with signals: {self.current_signal}")
                for realization in tqdm(range(self.nrealizations)):
                    pass
                    # TODO:
                    # 5- define signals
                    # 4- create PTA
                    # 6- make realizations
                self._current_signal = signal
                self.save()

        # finalizing
        if self._finish == "Not finished yet.":
            self._finish = datetime.datetime.now()
        self._status = "finish"
        self.save()
        print("Simulating has been finished.")
        self.summary()

    def validate_signals(self):
        if self.signals is None:
            self._signals = ["wn"]
        if "wn" not in self.signals:
            self._signals = ["wn"] + self._signals
        if self.signals[0] != "wn":
            self._signals.remove("wn")
            self._signals = ["wn"] + self._signals

    def create_pulsars(self):
        toasfiles = glob.glob(self.TOAs_Simulator.outdir + "/*.toas")
        for toasfile in tqdm(toasfiles):
            PSR = os.path.basename(toasfile)[:-5]
            psrfile = self.psrdir + "/" + PSR + ".psr"
            if not os.path.exists(psrfile):
                parfile = self.TimingModel_Simulator.outdir + "/" + PSR + ".par"
                parfile_tdb = parfile.replace(".par", "_tdb.par")
                if os.path.exists(parfile_tdb):
                    parfile = parfile_tdb
                model = get_model(parfile)
                with open(toasfile, "rb") as file:
                    toas = pickle.load(file)
                psr = Pulsar(model, toas, planets=False)
                with open(psrfile, "wb") as file:
                    pickle.dump(psr, file)
        self.load_pulsars(all_pulsars=True)
        self._exclude = self._exclude + [
            psr2.name
            for psr_idx, psr1 in enumerate(self.psrs)
            for psr2 in self.psrs[psr_idx + 1 :]
            if np.all(psr1.pos == psr2.pos)
        ]
        self.load_pulsars()

    def load_pulsars(self, all_pulsars=False):
        psrfiles = glob.glob(self.psrdir + "/*.psr")
        if len(psrfiles) > 0:
            psrs = []
            for psrfile in psrfiles:
                PSR = os.path.basename(psrfile)[:-4]
                if PSR in self.exclude:
                    continue
                with open(psrfile, "rb") as file:
                    psrs.append(pickle.load(file))
                if len(psrs) == self.npsrs and not all_pulsars:
                    break
            self._psrs = psrs
        else:
            self._psrs = None

    def save(self):
        file_name = CCG_CACHEDIR + "/" + self.name + ".sim"
        del self.psrs
        with open(file_name, "wb") as file:
            pickle.dump(self, file)
        self.load_pulsars()

    def summary(self):
        msg = f"""======================================
Simulation: {self.name}
======================================
Simulation started at: {self._begin}
Simulation finished at: {self._finish}
Number of pulsars: {self.npsrs}
Signals to be simulated: {self.signals}
Number of realizations: {self.nrealizations}
Output folder: '{self.outdir}'
Simulation stage: {self.status}
Current_signal: {self.current_signal}
======================================
"""
        print(msg)


def list_all_simulations():
    simulator_files = glob.glob(CCG_CACHEDIR + "/*.sim")
    print("Available simulators:")
    for simulator_idx, simulator_file in enumerate(simulator_files):
        simulator = os.path.basename(simulator_file)[:-4]
        print(f"{simulator_idx}- {simulator}")


def load_Simulator(name):
    file_name = CCG_CACHEDIR + "/" + name + ".sim"
    if os.path.exists(file_name):
        with open(file_name, "rb") as file:
            Simulator = pickle.load(file)
            Simulator.load_pulsars()
        return Simulator
    else:
        print(f"Warning:: There is no simulator named {name}.")
        list_all_simulations()
