# timing_model_simulation.py
"""
Simulate times of arrivals (TOAs).
"""

import os
import glob
from astropy import units as u
from astropy.time import Time

class TOAs_Simulator(object):

    def __init__(self, 
                 start_time='2000-01-01 00:00:00',
                 end_time='2020-01-01 00:00:00',
                 cadence=30,
                 fuzz=0.0,
                 radio_frequency=1400.0,
                 observatory='GBT',
                 add_noise=False,
                 add_correlated_noise=False,
                 wideband=False,
                 wideband_dm_error=0.0,
                 name_prefix='fake',
                 include_bipm=False,
                 include_gps=False,
                 multi_freqs_in_epoch=False,
                 flags=None,
                 subtract_mean=True,
                 pardir=None, outdir=None):

        if outdir is None:
            outdir = os.getcwd()

        if pardir is None:
            pardir = os.getcwd()

        parfiles = glob.glob(pardir+"/*.par")
        if len(parfiles) == 0:
            print(f'Warning: timing model files not found in: "{pardir}"')

        self._start_time = Time(start_time, format='iso')
        self._end_time = Time(end_time, format='iso')
        self._cadence = cadence * u.d
        self._fuzz = fuzz * u.d
        self._radio_frequency = radio_frequency * u.MHz
        self._observatory = observatory
        self._add_noise = add_noise
        self._add_correlated_noise = add_correlated_noise
        self._wideband = wideband
        self._wideband_dm_error = wideband_dm_error
        self._name_prefix = name_prefix
        self._include_bipm = include_bipm
        self._include_gps = include_gps
        self._multi_freqs_in_epoch = multi_freqs_in_epoch
        self._flags = flags
        self._subtract_mean = subtract_mean
        self._parfiles = parfiles
        self._nmodels = len(parfiles)
        self._outdir = os.path.abspath(outdir)

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

    @property
    def start_time(self):
        return self._start_time
    @start_time.setter
    def start_time(self, value):
        self._start_time = Time(value, format='iso')
    
    @property
    def end_time(self):
        return self._end_time
    @end_time.setter
    def end_time(self, value):
        self._end_time = Time(value, format='iso')
    
    @property
    def cadence(self):
        return self._cadence
    @cadence.setter
    def cadence(self, value):
        self._cadence = value * u.d
    
    @property
    def fuzz(self):
        return self._fuzz
    @fuzz.setter
    def fuzz(self, value):
        self._fuzz = value * u.d
    
    @property
    def radio_frequency(self):
        return self._radio_frequency
    @radio_frequency.setter
    def radio_frequency(self, value):
        self._radio_frequency = value * u.MHz
    
    @property
    def observatory(self):
        return self._observatory
    @observatory.setter
    def observatory(self, value):
        self._observatory = value
    
    @property
    def add_noise(self):
        return self._add_noise
    @add_noise.setter
    def add_noise(self, value):
        self._add_noise = value
    
    @property
    def add_correlated_noise(self):
        return self._add_correlated_noise
    @add_correlated_noise.setter
    def add_correlated_noise(self, value):
        self._add_correlated_noise = value
    
    @property
    def wideband(self):
        return self._wideband
    @wideband.setter
    def wideband(self, value):
        self._wideband = value
    
    @property
    def wideband_dm_error(self):
        return self._wideband_dm_error
    @wideband_dm_error.setter
    def wideband_dm_error(self, value):
        self._wideband_dm_error = value
    
    @property
    def name_prefix(self):
        return self._name_prefix
    @name_prefix.setter
    def name_prefix(self, value):
        self._name_prefix = value
    
    @property
    def include_bipm(self):
        return self._include_bipm
    @include_bipm.setter
    def include_bipm(self, value):
        self._include_bipm = value
    
    @property
    def include_gps(self):
        return self._include_gps
    @include_gps.setter
    def include_gps(self, value):
        self._include_gps = value
    
    @property
    def multi_freqs_in_epoch(self):
        return self._multi_freqs_in_epoch
    @multi_freqs_in_epoch.setter
    def multi_freqs_in_epoch(self, value):
        self._multi_freqs_in_epoch = value
    
    @property
    def flags(self):
        return self._flags
    @flags.setter
    def flags(self, value):
        self._flags = value
    
    @property
    def subtract_mean(self):
        return self._subtract_mean
    @subtract_mean.setter
    def subtract_mean(self, value):
        self._subtract_mean = value
    
    @property
    def parfiles(self):
        return self._parfiles
    @parfiles.setter
    def parfiles(self, value):
        print("Warning:: Searching for timing model files is done automaticaly.")

    @property
    def nmodels(self):
        return self._nmodels
    @nmodels.setter
    def nmodels(self, value):
        print("Warning:: Number of timing model files is automaticaly calculated.")

    @property
    def outdir(self):
        return self._outdir
    @outdir.setter
    def outdir(self, value):
        self._outdirt = os.path.abspath(value)

    def start(self):
        pass