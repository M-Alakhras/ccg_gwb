# one_pulsar_tools.py
"""
functions to analyze one pulsar signal.
"""
import os
import shutil

import enterprise.signals.parameter as parameter
import numpy as np
from astropy import units as u
from enterprise.signals import gp_signals, signal_base, utils, white_signals
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
from scipy.signal import periodogram


def _calculate_sampling_frequency(toas):
    diff = np.diff(toas)
    fs = 1 / np.median(diff)
    return fs * u.Hz


def _calculate_power_spectrum(toas, res):
    fs = _calculate_sampling_frequency(toas)
    (f, S) = periodogram(res, fs=fs.value, scaling="spectrum")
    f = f * u.Hz
    S = S * u.s**2
    return f.to(u.Hz), S.to(u.us**2)


def _estimate_white_noise(psr):
    # Uniform prior on EFAC
    efac = parameter.Uniform(0.01, 10.0)
    log10_equad = parameter.Uniform(-8.5, -5)
    log10_ecorr = parameter.Uniform(-8.5, -5)
    # red noise parameters
    # Uniform in log10 Amplitude and in spectral index
    log10_A = parameter.Uniform(-20, -11)
    gamma = parameter.Uniform(0, 7)

    # Set up signals
    # white noise
    ef = white_signals.MeasurementNoise(efac=efac, log10_t2equad=log10_equad)
    ec = white_signals.EcorrKernelNoise(log10_ecorr=log10_ecorr)
    # red noise (powerlaw with 30 frequencies)
    pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
    rn = gp_signals.FourierBasisGP(spectrum=pl, components=30)
    # timing model
    # tm = gp_signals.TimingModel()

    # full model is sum of components
    model = ef + ec + rn  # + tm

    # initialize PTA
    pta = signal_base.PTA([model(psr)])
    x0 = np.hstack([p.sample() for p in pta.params])
    ndim = len(x0)
    cov = np.diag(np.ones(ndim))
    sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov, outDir=os.getcwd() + "/chains", resume=False)

    # sampler for N steps
    N = int(1e5)

    # SCAM = Single Component Adaptive Metropolis
    # AM = Adaptive Metropolis
    # DE = Differential Evolution
    # You can keep all these set at default values
    sampler.sample(
        x0,
        N,
        SCAMweight=30,
        AMweight=15,
        DEweight=50,
    )
    chain = np.loadtxt(os.getcwd() + "/chains/chain_1.txt")

    burn = int(0.1 * chain.shape[0])
    ind_efac = list(pta.param_names).index(psr.name + "_efac")
    ind_log10_t2equad = list(pta.param_names).index(psr.name + "_log10_t2equad")
    ind_log10_ecorr = list(pta.param_names).index(psr.name + "_log10_ecorr")
    efac = np.median(chain[burn:, ind_efac])
    log10_t2equad = np.median(chain[burn:, ind_log10_t2equad])
    log10_ecorr = np.median(chain[burn:, ind_log10_ecorr])
    shutil.rmtree(os.getcwd() + "/chains")

    return efac, log10_t2equad, log10_ecorr
