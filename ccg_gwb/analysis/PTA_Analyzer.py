# PTA_Analyzer.py
"""Define Analyzer class."""

import glob
import json
import os
import pickle

import numpy as np
from enterprise.signals.utils import set_residuals
from enterprise_extensions.frequentist import optimal_statistic as optstat
from tqdm.auto import tqdm

from ccg_gwb.celestial_map import _pulsars_separations, _pulsars_sky_coordinates
from ccg_gwb.multi_pulsar_tools import (
    _binning_pulsars_cross_correlation,
    _calculate_cross_correlation,
    _calculate_noise_weighted_cross_correlation,
    _compute_os,
    pta_model,
)
from ccg_gwb.network_tools import (
    _calculate_network_measures,
    _create_multilayer_network,
    _create_network,
    _thresholding_network,
    all_correlation_network_measures,
    all_multilayer_network_measures,
)
from ccg_gwb.one_pulsar_tools import _estimate_white_noise, _lowpass_filter
from ccg_gwb.plotter import (
    _plot_celestial_map,
    _plot_cross_correlation,
    _plot_network,
    _plot_network_measures,
    _plot_noise_weighted_cross_correlation,
    _plot_power_spectrum,
    _plot_residuals,
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
        self.refresh()

    def refresh(self):
        self._res = None
        self._wn_params = None
        self._cor = None
        self._cor_mode = "valid"
        self._cor_filter = "None"
        self._cor_cutoff = 5e-9
        self._cor_nlags = 21
        self._OS = None
        self._xi = None
        self._ORF = None
        self._rho = None
        self._sig = None
        self._OptStat = None
        self._OptStat_sig = None
        self._SNR = None
        self._network = None

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
    def cor(self):
        return self._cor

    @cor.setter
    def cor(self, value):
        print("Warning:: use calculate_cross_correlation().")

    @property
    def cor_mode(self):
        return self._cor_mode

    @cor_mode.setter
    def cor_mode(self, value):
        print("Warning:: use calculate_cross_correlation().")

    @property
    def cor_filter(self):
        return self._cor_filter

    @cor_filter.setter
    def cor_filter(self, value):
        print("Warning:: use calculate_cross_correlation().")

    @property
    def cor_cutoff(self):
        return self._cor_cutoff

    @cor_cutoff.setter
    def cor_cutoff(self, value):
        print("Warning:: use calculate_cross_correlation().")

    @property
    def cor_nlags(self):
        return self._cor_nlags

    @cor_nlags.setter
    def cor_nlags(self, value):
        print("Warning:: use calculate_cross_correlation().")

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

    @property
    def network(self):
        return self._network

    @network.setter
    def network(self, value):
        print("Warning:: use create_network().")

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
        self.refresh()
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
        if "gw" in signal:
            log10_A = signal.split("_")[1]
            log10_A = -float(log10_A[6:])
        else:
            log10_A = None
        self.initialize_OS(log10_A=log10_A)

    def initialize_OS(self, log10_A=None):
        pta = pta_model(self.psrs, noisedict=self.wn_params, log10_A_common=log10_A, gamma_common=13 / 3)
        self._OS = optstat.OptimalStatistic(psrs=self.psrs, pta=pta, noisedict=self.wn_params)

    def pulsars_sky_coordinates(self):
        return _pulsars_sky_coordinates(self.psrs)

    def pulsars_separations(self):
        return _pulsars_separations(self.psrs)

    def plot_celestial_map(self):
        _plot_celestial_map(self.psrs)

    def plot_separations_hist(self):
        _plot_separations_hist(self.psrs)

    def plot_residuals(self, psr, filter_residuals="None", cutoff=5e-9, order=6):
        _plot_residuals(psr, filter_residuals=filter_residuals, cutoff=cutoff, order=order)

    def plot_power_spectrum(self, psr, filter_residuals="None", cutoff=5e-9, order=6):
        _plot_power_spectrum(psr, filter_residuals=filter_residuals, cutoff=cutoff, order=order)

    def estimate_white_noise(self, psr):
        return _estimate_white_noise(psr)

    def calculate_cross_correlation(self, mode="valid", filter_residuals="None", cutoff=5e-9, order=6, nlags=21):
        if self._cor is None or self._cor_mode != mode or self._cor_cutoff != cutoff:
            self._cor_mode = mode
            self._cor_filter = filter_residuals
            self._cor_cutoff = cutoff
            self._cor = _calculate_cross_correlation(
                self.psrs, mode=mode, filter_residuals=filter_residuals, cutoff=cutoff, order=order, nlags=nlags
            )

    def plot_cross_correlation(self):
        if self._cor is None:
            print("Warning:: use calculate_cross_correlation() first.")
        _plot_cross_correlation(self.psrs, self._cor)

    def calculate_optimal_statistic(self, loud=True):
        if self._SNR is None:
            self._xi, self._ORF, self._rho, self._sig, self._OptStat, self._OptStat_sig = _compute_os(
                self.OS, params=self.wn_params
            )
            self._SNR = self._OptStat / self._OptStat_sig
        if loud:
            print("A^2 = {} +- {}".format(self._OptStat, self._OptStat_sig))
            print("SNR = {}".format(self._SNR))

    def plot_noise_weighted_cross_correlation(self):
        if self.rho is None:
            print("Warning:: use calculate_optimal_statistic() first.")
            return
        rho_matrix = _calculate_noise_weighted_cross_correlation(self.psrs, self.rho)
        _plot_noise_weighted_cross_correlation(self.psrs, rho_matrix)

    def plot_rho(self, binning=False, bins=9):
        if self.rho is None:
            print("Warning:: use calculate_optimal_statistic() first.")
            return
        xi = self.xi
        rho = self.rho
        sig = self.sig
        if binning:
            xi, rho, sig = _binning_pulsars_cross_correlation(self.xi, self.rho, self.sig, bins=bins)
        _plot_rho(xi, rho, sig, self.OptStat)

    def create_network(self, method="correlation"):
        if method == "correlation":
            if self._cor is None:
                print("Warning:: run calculate_cross_calculation() first.")
                return
            psr_names = [psr.name for psr in self.psrs]
            self._network = _create_network(self._cor, psr_names)
        elif method == "multilayer":
            self._network = _create_multilayer_network(self.psrs)

    def plot_network(self):
        if self._network is None:
            print("Warning:: use create_network() first.")
            return
        _plot_network(self._network)

    def thresholding_network(self, th):
        return _thresholding_network(self._network, th)

    def plot_threshold_network(self, th):
        if self._network is None:
            print("Warning:: use create_network() first.")
            return
        G = self.thresholding_network(th)
        _plot_network(G)

    def calculate_network_measures(self):
        if self._network is None:
            print("Warning:: use create_network() first.")
        return _calculate_network_measures(self._network)

    def calculate_threshold_network_measures(self, th):
        if self._network is None:
            print("Warning:: use create_network() first.")
        G = self.thresholding_network(th)
        return _calculate_network_measures(G)

    def calculate_thresholds(self, n=21):
        nrealizations = len(glob.glob(self.signalsdir + "/wn/*.res"))
        th_min = 999999
        th_max = -999999
        print("calculating thresholds...")
        for realization in tqdm(range(nrealizations)):
            self.load_realization(signal="wn", realization=realization + 1)
            cor_matrix = _calculate_cross_correlation(
                self.psrs,
                mode=self.cor_mode,
                filter_residuals=self.cor_filter,
                cutoff=self.cor_cutoff,
                order=6,
                nlags=self.cor_nlags,
            )
            if th_min > np.min(cor_matrix):
                th_min = np.min(cor_matrix)
            if th_max < np.max(cor_matrix):
                th_max = np.max(cor_matrix)
        thresholds = np.linspace(th_min, th_max, n)
        print("Thresholds = ", thresholds)
        return thresholds

    def calculate_all_realizations_network_measures(
        self,
        method="correlation",
        thresholds=None,
        correlation_mode="valid",
        filter_residuals="None",
        cutoff=5e-9,
        order=6,
        nlags=21,
    ):
        if filter_residuals == "lowpass":
            res = [_lowpass_filter(psr, cutoff=cutoff, order=order) for psr in self.psrs]
            set_residuals(self.psrs, res)
        if method == "correlation" and correlation_mode in ["valid", "full", "partial"] and thresholds is None:
            thresholds = self.calculate_thresholds(n=21)
        measures = {}
        SNRs = {}
        for signal in self.signals:
            if "wn" not in signal:
                continue
            realizations = glob.glob(signal + "/*.res")
            nrealizations = len(realizations)
            signal_name = signal.split("/")[-1]
            print(signal_name)
            measures.update({signal_name: {}})
            SNRs.update({signal_name: np.zeros(nrealizations)})
            if method == "correlation":
                for measure in all_correlation_network_measures:
                    if correlation_mode in ["valid", "full", "partial"]:
                        measures[signal_name].update({measure: np.zeros((nrealizations, len(thresholds)))})
                    elif correlation_mode == "SSS":
                        measures[signal_name].update({measure: np.zeros(nrealizations)})
                    else:
                        print("Warning:: Unknown correlation mode.")
                        return
                for realization in tqdm(range(nrealizations)):
                    self.load_realization(signal=signal_name, realization=realization + 1)
                    self.calculate_cross_correlation(
                        mode=correlation_mode,
                        filter_residuals=filter_residuals,
                        cutoff=cutoff,
                        order=order,
                        nlags=nlags,
                    )
                    self.create_network()
                    self.calculate_optimal_statistic(loud=False)
                    SNRs[signal_name][realization] = self.SNR
                    if correlation_mode in ["valid", "full", "partial"]:
                        for ii, th in enumerate(thresholds):
                            network_measures = self.calculate_threshold_network_measures(th)
                            for measure in all_correlation_network_measures:
                                measures[signal_name][measure][realization, ii] = network_measures[measure]
                    elif correlation_mode == "SSS":
                        network_measures = self.calculate_network_measures()
                        for measure in all_correlation_network_measures:
                            measures[signal_name][measure][realization] = network_measures[measure]
            elif method == "multilayer":
                for measure in all_multilayer_network_measures:
                    measures[signal_name].update({measure: np.zeros(nrealizations)})
                for realization in tqdm(range(nrealizations)):
                    self.load_realization(signal=signal_name, realization=realization + 1)
                    self.create_network(method="multilayer")
                    self.calculate_optimal_statistic(loud=False)
                    SNRs[signal_name][realization] = self.SNR
                    network_measures = self.calculate_network_measures()
                    for measure in all_multilayer_network_measures:
                        measures[signal_name][measure][realization] = network_measures[measure]
        set_residuals(self.psrs, self.res)
        return thresholds, SNRs, measures

    def plot_all_realizations_network_measures(self, thresholds, SNRs, measures, nth_plot=None):
        if thresholds is None:
            mask = None
        else:
            if nth_plot is None:
                nth_plot = len(thresholds)
            mask = np.linspace(0, len(thresholds) - 1, nth_plot, dtype=int)
        _plot_network_measures(thresholds=thresholds, SNRs=SNRs, measures=measures, mask=mask)
