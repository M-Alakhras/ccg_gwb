# ccg_gwb
"""Main module."""

import glob
import json
import os
import pickle

import numpy as np
from enterprise.signals.utils import set_residuals
from enterprise_extensions.frequentist import optimal_statistic as optstat
from tqdm.auto import tqdm

from ccg_gwb.celestial_map import _pulsars_radec, _pulsars_separations, _pulsars_sky_coordinates
from ccg_gwb.multi_pulsar_tools import (
    _calculate_cross_correlation,
    _calculate_noise_weighted_cross_correlation,
    _compute_os,
    pta_model,
)
from ccg_gwb.one_pulsar_tools import _calculate_power_spectrum, _estimate_white_noise
from ccg_gwb.plotter import (
    _plot_celestial_map,
    _plot_cross_correlation,
    _plot_noise_weighted_cross_correlation,
    _plot_ORF,
    _plot_power_spectrum,
    _plot_rho,
    _plot_separations_hist,
)


class Analyzer:

    def __init__(self, datadir=None):
        if datadir is None:
            datadir = os.getcwd()

        self._datadir = os.path.abspath(datadir)
        self._pardir = self.datadir + "/par"
        self._toasdir = self.datadir + "/toas"
        self._psrdir = self.datadir + "/psr"
        self._signalsdir = self.datadir + "/signals"
        self._signals = glob.glob(self.signalsdir + "/*")
        self._psrs = None
        self._res = None
        self._wn_params = None
        self._OS = None
        self._xi = None
        self._ORF = None
        self._rho = None
        self._sig = None
        self._OptStat = None
        self._OptStat_sig = None
        self._SNR = None

    @property
    def datadir(self):
        return self._datadir

    @datadir.setter
    def datadir(self, value):
        self._datadir = os.path.abspath(value)

    @property
    def pardir(self):
        return self._pardir

    @pardir.setter
    def pardir(self, value):
        self._pardir = os.path.abspath(value)

    @property
    def toasdir(self):
        return self._toasdir

    @toasdir.setter
    def toasdir(self, value):
        self._toasdir = os.path.abspath(value)

    @property
    def psrdir(self):
        return self._psrdir

    @psrdir.setter
    def psrdir(self, value):
        self._psrdir = os.path.abspath(value)

    @property
    def signalsdir(self):
        return self._signalsdir

    @signalsdir.setter
    def signalsdir(self, value):
        self._signalsdir = os.path.abspath(value)
        self._signals = glob.glob(self._signalsdir + "/*")

    @property
    def signals(self):
        return self._signals

    @signals.setter
    def signals(self, value):
        print("Warning:: Signals are loaded automatically.")

    @property
    def psrs(self):
        return self._psrs

    @psrs.setter
    def psrs(self, value):
        print("Warning:: Pulsars are loaded automatically.")

    @property
    def res(self):
        return self._res

    @res.setter
    def res(self, value):
        print("Warning:: use load_realization().")

    @property
    def wn_params(self):
        return self._wn_params

    @wn_params.setter
    def wn_params(self, value):
        print("Warning:: use load_realization().")

    @property
    def OS(self):
        return self._OS

    @OS.setter
    def OS(self, value):
        print("Warning:: use load_realization().")

    @property
    def xi(self):
        return self._xi

    @xi.setter
    def xi(self, value):
        print("Warning:: use compute_optimal_statistic().")

    @property
    def ORF(self):
        return self._ORF

    @ORF.setter
    def ORF(self, value):
        print("Warning:: use compute_optimal_statistic().")

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, value):
        print("Warning:: use compute_optimal_statistic().")

    @property
    def sig(self):
        return self._sig

    @sig.setter
    def sig(self, value):
        print("Warning:: use compute_optimal_statistic().")

    @property
    def OptStat(self):
        return self._OptStat

    @OptStat.setter
    def OptStat(self, value):
        print("Warning:: use compute_optimal_statistic().")

    @property
    def OptStat_sig(self):
        return self._OptStat_sig

    @OptStat_sig.setter
    def OptStat_sig(self, value):
        print("Warning:: use compute_optimal_statistic().")

    @property
    def SNR(self):
        return self._SNR

    @SNR.setter
    def SNR(self, value):
        print("Warning:: use compute_optimal_statistic().")

    def load_pulsars(self, all_pulsars=False):
        psrsfiles = glob.glob(self.psrdir + "/*.psrs")
        psrfiles = glob.glob(self.psrdir + "/*.psr")
        if len(psrsfiles) > 0 and not all_pulsars:
            print(f'loading pulsars from: "{psrsfiles[0]}"')
            with open(psrsfiles[0], "rb") as file:
                self._psrs = pickle.load(file)
        elif len(psrfiles) > 0:
            psrs = []
            for psrfile in psrfiles:
                with open(psrfile, "rb") as file:
                    psrs.append(pickle.load(file))
            self._psrs = psrs
        else:
            self._psrs = None

    def load_realization(self, signal="wn", realization=1):
        signaldir = self.datadir + "/signals/" + signal
        resfile = signaldir + f"/realization_{realization}.res"
        wn_params_file = self.datadir + f"/signals/wn/params_{realization}.json"
        if os.path.exists(resfile):
            with open(resfile, "rb") as file:
                self._res = pickle.load(file)
                set_residuals(self.psrs, self.res)
        else:
            print(f"Warning:: realization {resfile} not found.")
            self._res = None
            self._wn_params = None
            return
        if os.path.exists(wn_params_file):
            with open(wn_params_file, "r") as file:
                self._wn_params = json.load(file)
        else:
            print(f"Warning:: white noise dictionary {wn_params_file} not found.")
            self._wn_params = None
            return
        self.initialize_OS()

    def initialize_OS(self):
        pta = pta_model(self.psrs, noisedict=self.wn_params, gamma_common=13 / 3)
        self._OS = optstat.OptimalStatistic(psrs=self.psrs, pta=pta)

    def pulsars_sky_coordinates(self):
        return _pulsars_sky_coordinates(self.psrs)

    def pulsars_radec(self):
        return _pulsars_radec(self.psrs)

    def pulsars_separations(self):
        return _pulsars_separations(self.psrs)

    def calculate_power_spectrum(self, psr):
        return _calculate_power_spectrum(psr.toas, psr.residuals)

    def plot_celestial_map(self):
        _plot_celestial_map(self.psrs)

    def plot_separations_hist(self):
        _plot_separations_hist(self.psrs)

    def plot_power_spectrum(self, psr):
        _plot_power_spectrum(psr)

    def estimate_white_noise(self, psr):
        return _estimate_white_noise(psr)

    def calculate_cross_correlation(self):
        return _calculate_cross_correlation(self.psrs)

    def plot_cross_correlation(self):
        _plot_cross_correlation(self.psrs)

    def compute_optimal_statistic(self, loud=True):
        self._xi, self._ORF, self._rho, self._sig, self._OptStat, self._OptStat_sig = _compute_os(self.OS)
        self._SNR = self._OptStat / self._OptStat_sig
        if loud:
            print("A^2 = {} +- {}".format(self._OptStat, self._OptStat_sig))
            print("SNR = {}".format(self._SNR))

    def plot_ORF(self):
        if self.ORF is None:
            self.compute_optimal_statistic(loud=False)
        _plot_ORF(self.xi, self.ORF)

    def plot_rho(self):
        if self.rho is None:
            self.compute_optimal_statistic(loud=False)
        _plot_rho(self.xi, self.rho, self.sig)

    def calculate_noise_weighted_cross_correlation(self):
        if self.rho is None:
            self.compute_optimal_statistic(loud=False)
        return _calculate_noise_weighted_cross_correlation(self.psrs, self.rho)

    def plot_noise_weighted_cross_correlation(self):
        if self.rho is None:
            self.compute_optimal_statistic(loud=False)
        _plot_noise_weighted_cross_correlation(self.psrs, self.rho)

    def realization_MCMC(self, nsamples=1000):
        if self.res is None:
            print("Warning:: load realization first.")
            return
        chain = {"rho": [], "sig": [], "OptStat": [], "OptStat_sig": []}
        for step in tqdm(range(nsamples)):
            self.compute_optimal_statistic(loud=False)
            chain["rho"].append(self.rho)
            chain["sig"].append(self.sig)
            chain["OptStat"].append(self.OptStat)
            chain["OptStat_sig"].append(self.OptStat_sig)
        chain["rho"] = np.array(chain["rho"])
        chain["sig"] = np.array(chain["sig"])
        chain["OptStat"] = np.array(chain["OptStat"])
        chain["OptStat_sig"] = np.array(chain["OptStat_sig"])
        rho_avg = np.average(chain["rho"], weights=1 / chain["sig"] ** 2, axis=0)
        sig_avg = 1 / np.sum(1 / chain["sig"] ** 2, axis=0)
        OptStat_avg = np.average(chain["OptStat"], weights=1 / chain["OptStat_sig"] ** 2, axis=0)
        OptStat_sig_avg = 1 / np.sum(1 / chain["OptStat_sig"] ** 2, axis=0)
        SNR_avg = OptStat_avg / OptStat_sig_avg
        print("A^2 = {} +- {}".format(OptStat_avg, OptStat_sig_avg))
        print("SNR = {}".format(SNR_avg))
        _plot_rho(self.xi, rho_avg, sig_avg)
        _plot_noise_weighted_cross_correlation(self.psrs, rho_avg)
        return chain, [rho_avg, sig_avg, OptStat_avg, OptStat_sig_avg]
