# PTA_Simulator.py
"""
Define PTA Simulator class.
"""

import os
import glob
import datetime
import pickle
from tqdm.auto import tqdm
from ccg_gwb import (CCG_ENV, CCG_PATH, CCG_CACHEDIR)
from ccg_gwb.simulation.timing_model_simulation import TimingModel_Simulator
from ccg_gwb.simulation.toas_simulation import TOAs_Simulator

        
class PTA_Simulator(object):

    def __init__(self, name, npsrs=20, nrealizations=200, signals=None, outdir=None):

        if outdir is None:
            outdir = os.getcwd()

        # setting Simulator parameters
        self._name = name
        self._npsrs = npsrs
        self._nrealizations = nrealizations
        self._signals = signals
        self.validate_signals()
        self._outdir = os.path.abspath(outdir)
        self._status = 'init'
        self._current_signal = self.signals[0]
        self._begin = 'Not started yet.'
        self._finish = 'Not finished yet.'

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        self.TimingModel_Simulator = TimingModel_Simulator(ATNF=True, ATNF_Condition='P0 < 0.03', outdir=self.outdir+"/par")
        self.TOAs_Simulator = TOAs_Simulator(pardir=self.outdir+"/par",outdir=self.outdir+"/toas")

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

    @property
    def outdir(self):
        return self._outdir

    @outdir.setter
    def outdir(self, value):
        self._outdir = os.path.abspath(value)

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

    def start(self, restart=False):
        # initialization
        if restart:
            self._begin = 'Not started yet.'
            self._finish = 'Not finished yet.'
            self._status = 'init'
        
        if self._begin == 'Not started yet.':
            self._begin = datetime.datetime.now()
        self.summary()
        print("Start simulating...")

        # TODO:
        if self.status == 'init':
            print("Simulating timing models...")
            self.TimingModel_Simulator.start()
            print(f'Timing models have been simulated and saved into: "{self.TimingModel_Simulator.outdir}"')
            self._status = 'pars'
            self.save()

        if self.status == 'pars':
            print("Simulating times of arrivals...")
            self.TOAs_Simulator.start()
            print(f'Times of arrivals have been simulated and saved into: "{self.TOAs_Simulator.outdir}"')
            self._status = 'toas'
            self.save()

        if self.status == 'toas':
            # 3- create pulsar objects
            self._status = 'psrs'
            self.save()

        if self.status == 'psrs':
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
        if self._finish == 'Not finished yet.':
            self._finish = datetime.datetime.now()
        self._status = 'finish'
        self.save()
        print("Simulating has been finished.")
        self.summary()

    def validate_signals(self):
        if self.signals is None:
            self._signals = ['wn']
        if 'wn' not in self.signals:
            self._signals = ['wn'] + self._signals
        if self.signals[0] != 'wn':
            self._signals.remove('wn')
            self._signals = ['wn'] + self._signals

    def save(self):
        file_name = CCG_CACHEDIR + "/" + self.name + ".sim"
        with open(file_name, 'wb') as file:
            pickle.dump(self, file)

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
    simulator_files = glob.glob(CCG_CACHEDIR +"/*.sim")
    print("Available simulators:")
    for simulator_idx, simulator_file in enumerate(simulator_files):
        simulator = os.path.basename(simulator_file)[:-4]
        print(f"{simulator_idx}- {simulator}")
    
    
def load_Simulator(name):
    file_name = CCG_CACHEDIR + "/" + name + ".sim"
    if os.path.exists(file_name):
        with open(file_name, 'rb') as file:
            Simulator = pickle.load(file)
        return Simulator
    else:
        print(f"Warning:: There is no simulator named {name}.")
        list_all_simulations()