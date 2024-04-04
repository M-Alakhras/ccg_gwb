# multi_pulsar_tools.py
"""
functions to analyze multi pulsar signals.
"""

import numpy as np
import pandas as pd
import pingouin as pg
import scipy.linalg as sl
from enterprise.signals import parameter, selections, signal_base, utils, white_signals
from enterprise_extensions import model_utils
from enterprise_extensions.blocks import common_red_noise_block

from ccg_gwb.one_pulsar_tools import _lowpass_filter


def sss(res):
    ii = np.arange(0, len(res), dtype=int)
    g = np.random.randn(len(res))
    ip = ii + 1.0 * g
    ihat = np.argsort(ip)
    return res[ihat]


def _calculate_cross_correlation(psrs, mode="valid", filter_residuals="None", cutoff=5e-9, order=6, nlags=301):
    npsrs = len(psrs)
    correlation_matrix = np.zeros((npsrs, npsrs))
    psr_names = [psr.name for psr in psrs]
    if filter_residuals == "None":
        psr_res = {psr.name: psr.residuals for psr in psrs}
    elif filter_residuals == "lowpass":
        psr_res = {psr.name: _lowpass_filter(psr, cutoff=cutoff, order=order) for psr in psrs}
    df = pd.DataFrame(psr_res)
    for ii, psr1 in enumerate(psr_names):
        for jj, psr2 in enumerate(psr_names):
            if ii == jj:
                continue
            res1 = psr_res[psr1]
            res2 = psr_res[psr2]
            sig1 = np.correlate(res1, res1, mode="valid")[0]
            sig2 = np.correlate(res2, res2, mode="valid")[0]
            if mode == "valid":
                cor = np.correlate(res1, res2, mode="valid") / np.sqrt(sig1 * sig2)
                correlation_matrix[ii, jj] = cor[0]
            elif mode == "full":
                cor = np.correlate(res1, res2, mode="full") / np.sqrt(sig1 * sig2)
                cor_idx = np.arange((len(cor) - nlags) // 2, (len(cor) + nlags) // 2, dtype=int)
                correlation_matrix[ii, jj] = np.max(cor[cor_idx])
            elif mode == "SSS":
                cor = np.correlate(res1, res2, mode="full") / np.sqrt(sig1 * sig2)
                cor_idx = np.arange((len(cor) - nlags) // 2, (len(cor) + nlags) // 2, dtype=int)
                cor_sss = []
                for nsss in range(99):
                    res1_sss = sss(res1)
                    res2_sss = sss(res2)
                    sig1_sss = np.correlate(res1_sss, res1_sss, mode="valid")[0]
                    sig2_sss = np.correlate(res2_sss, res2_sss, mode="valid")[0]
                    cor_sss.append(np.correlate(res1_sss, res2_sss, mode="full") / np.sqrt(sig1_sss * sig2_sss))
                cor_sss = np.array(cor_sss)
                correlation_matrix[ii, jj] = np.any(np.all(cor[cor_idx] > cor_sss[:, cor_idx], axis=0))
            elif mode == "partial":
                r = pg.partial_corr(data=df, x=psr1, y=psr2)
                correlation_matrix[ii, jj] = r.r["pearson"]
            else:
                print(f"Warning:: Uknown correlation mode '{mode}'")
    return correlation_matrix


def _white_noise_block():
    backend = selections.Selection(selections.no_selection)

    # white noise parameters
    efac = parameter.Constant()
    equad = parameter.Constant()
    ecorr = parameter.Constant()

    # white noise signals
    efeq = white_signals.MeasurementNoise(efac=efac, log10_t2equad=equad, selection=backend, name=None)
    ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=backend, name=None)
    # combine signals
    s = efeq + ec
    return s


def pta_model(
    psrs, psd="powerlaw", noisedict=None, n_gwbfreqs=30, log10_A_common=None, gamma_common=None, logmin=-17, logmax=-9
):

    # find the maximum time span to set GW frequency sampling
    Tspan = model_utils.get_tspan(psrs)

    # common red noise block
    s = common_red_noise_block(
        psd="powerlaw",
        prior="log-uniform",
        Tspan=Tspan,
        components=n_gwbfreqs,
        log10_A_val=log10_A_common,
        gamma_val=gamma_common,
        orf="hd",
        name="gw",
        logmin=logmin,
        logmax=logmax,
        combine=False,
    )

    # adding white-noise, and acting on psr objects
    models = []
    for p in psrs:
        s3 = s + _white_noise_block()
        models.append(s3(p))

    # set up PTA
    pta = signal_base.PTA(models)

    # set white noise parameters
    if noisedict is None:
        print("No noise dictionary provided!...")
    else:
        noisedict = noisedict
        pta.set_default_params(noisedict)

    return pta


def get_phiinvs_diag(array, nrows, ncols):
    """Split a matrix into sub-matrices."""
    r, h = array.shape
    a = array.reshape(h // nrows, nrows, -1, ncols).swapaxes(1, 2).reshape(-1, nrows, ncols)
    b = np.diag(
        np.array(np.linspace(0, ((r / nrows) * (h / ncols)) - 1, int((r / nrows) * (h / ncols)), dtype=int)).reshape(
            int(r / nrows), int(h / ncols)
        )
    )
    return a[b]


def _compute_os(OS, params=None, loud=False):
    npsrs = len(OS.pta._signalcollections)

    if params is None:
        params = {name: par.sample() for name, par in zip(OS.pta.param_names, OS.pta.params)}
    else:
        # check to see that the params dictionary includes values
        # for all of the parameters in the model
        for p in OS.pta.param_names:
            if p not in params.keys():
                msg = "{0} is not included ".format(p)
                msg += "in the parameter dictionary. "
                msg += "Drawing a random value."
                if loud:
                    print(msg)

    TNrs = OS.get_TNr(params=params)
    TNTs = OS.get_TNT(params=params)
    FNrs = OS.get_FNr(params=params)
    FNFs = OS.get_FNF(params=params)
    FNTs = OS.get_FNT(params=params)

    r, h = TNTs[0].shape
    phiinvs = OS.pta.get_phiinv(params, logdet=False)
    phiinvs = get_phiinvs_diag(phiinvs, r, h)

    X, Z = [], []
    for TNr, TNT, FNr, FNF, FNT, phiinv in zip(TNrs, TNTs, FNrs, FNFs, FNTs, phiinvs):
        Sigma = TNT + (np.diag(phiinv) if phiinv.ndim == 1 else phiinv)
        try:
            cf = sl.cho_factor(Sigma)
            SigmaTNr = sl.cho_solve(cf, TNr)
            SigmaTNF = sl.cho_solve(cf, FNT.T)
        except np.linalg.LinAlgError:
            SigmaTNr = np.linalg.solve(Sigma, TNr)
            SigmaTNF = np.linalg.solve(Sigma, FNT.T)

        FNTSigmaTNr = np.dot(FNT, SigmaTNr)
        X.append(FNr - FNTSigmaTNr)
        Z.append(FNF - np.dot(FNT, SigmaTNF))

    rho, sig, ORF, xi = [], [], [], []
    for ii in range(npsrs):
        for jj in range(ii + 1, npsrs):
            if OS.gamma_common is None and "gw_gamma" in params.keys():
                phiIJ = utils.powerlaw(OS.freqs, log10_A=0, gamma=params["gw_gamma"])
            else:
                phiIJ = utils.powerlaw(OS.freqs, log10_A=0, gamma=OS.gamma_common)
            top = np.dot(X[ii], phiIJ * X[jj])
            bot = np.trace(np.dot(Z[ii] * phiIJ[None, :], Z[jj] * phiIJ[None, :]))
            # cross correlation and uncertainty
            rho.append(top / bot)
            sig.append(1 / np.sqrt(bot))
            # Overlap reduction function for PSRs ii, jj
            ORF.append(OS.orf(OS.psrlocs[ii], OS.psrlocs[jj]))
            # angular separation
            xi.append(np.arccos(np.dot(OS.psrlocs[ii], OS.psrlocs[jj])))

    rho = np.array(rho)
    sig = np.array(sig)
    ORF = np.array(ORF)
    xi = np.array(xi)
    OptStat = np.sum(rho * ORF / sig**2) / np.sum(ORF**2 / sig**2)
    OptStat_sig = 1 / np.sqrt(np.sum(ORF**2 / sig**2))

    return xi, ORF, rho, sig, OptStat, OptStat_sig


def _calculate_noise_weighted_cross_correlation(psrs, rho):
    npsrs = len(psrs)
    rho_matrix = np.zeros((npsrs, npsrs))
    idx = 0
    for ii in range(npsrs):
        for jj in range(ii + 1, npsrs):
            rho_matrix[ii, jj] = rho_matrix[jj, ii] = rho[idx]
            idx += 1
    return rho_matrix


def _binning_pulsars_cross_correlation(xi, rho, sig, bins=9):
    bins_edge = np.linspace(0, np.pi, bins + 1)
    bins_center = bins_edge[:-1] + (bins_edge[1] - bins_edge[0]) / 2
    bins_num = np.zeros(bins)
    rho_bins = np.zeros(bins)
    sig_bins = np.zeros(bins)
    for bb in range(bins):
        for ii, x in enumerate(xi):
            if x >= bins_edge[bb] and x < bins_edge[bb + 1]:
                bins_num[bb] += 1
                rho_bins[bb] += rho[ii] / sig[ii] ** 2
                sig_bins[bb] += 1 / sig[ii] ** 2
    mask = np.where(bins_num != 0)
    bins_center = bins_center[mask]
    bins = len(bins_center)
    rho_bins = rho_bins[mask].reshape(bins)
    sig_bins = sig_bins[mask].reshape(bins)
    rho_bins = rho_bins / sig_bins
    sig_bins = 1 / np.sqrt(sig_bins)
    return bins_center, rho_bins, sig_bins
