# plotter.py
"""
plotting tools.
"""

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from ccg_gwb.celestial_map import _pulsars_radec, _pulsars_separations
from ccg_gwb.network_tools import all_network_measures
from ccg_gwb.one_pulsar_tools import _calculate_power_spectrum, _lowpass_filter


def hd_func(xii):
    if xii == 0:
        return 0.5
    else:
        omc2 = (1 - np.cos(xii)) / 2
        return 1.5 * omc2 * np.log(omc2) - 0.25 * omc2 + 0.5


def _plot_celestial_map(psrs, figsize=(14, 4.2), color="green"):
    ra, dec = _pulsars_radec(psrs)
    plt.figure(figsize=figsize)
    plt.subplot(111, projection="mollweide")
    plt.title("Celestial Map")
    plt.plot(ra, dec, "*", markersize=10, alpha=0.5, color=color)
    plt.subplots_adjust(top=0.95, bottom=0.0)
    plt.grid(True, alpha=0.3)
    plt.show()


def _plot_separations_hist(psrs, bins=20, figsize=(8, 5), color="orange"):
    sep = _pulsars_separations(psrs)
    plt.figure(figsize=figsize)
    plt.hist(sep[np.triu_indices(len(psrs), k=1)], bins=bins, histtype="step", color=color, linewidth=2)
    plt.xlim([0, 180])
    plt.title("Angular Separations")
    plt.xlabel("Angular seperation [$\degree$]")
    plt.ylabel("Number of pulsar-pairs")
    plt.grid(alpha=0.3)
    plt.show()


def _plot_residuals(psr, filter_residuals="None", cutoff=2e-9, order=6):
    plt.plot(psr.toas, psr.residuals / 1e-6, label="Original data")
    if filter_residuals == "lowpass":
        res = _lowpass_filter(psr, cutoff=cutoff, order=order)
        plt.plot(psr.toas, res / 1e-6, label="Filtered data")
    plt.title("Pulsar Timing Residuals")
    plt.xlabel("Times of arrivals [s]")
    plt.ylabel(r"Timing residuals [$\mu s$]")
    plt.ylim([-10, 10])
    plt.grid(which="major", alpha=0.8)
    plt.grid(which="minor", alpha=0.2)
    plt.legend()
    plt.show()


def _plot_power_spectrum(psr, filter_residuals="None", cutoff=2e-9, order=6):
    f, S = _calculate_power_spectrum(psr.toas, psr.residuals)
    plt.loglog(f, S, label="Original data")
    if filter_residuals == "lowpass":
        res = _lowpass_filter(psr, cutoff=cutoff, order=order)
        f, S = _calculate_power_spectrum(psr.toas, res)
        plt.loglog(f, S, label="Filtered data")
    plt.title("Power Spectrum")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel(r"Timing residuals Power [$\mu s^2$]")
    plt.ylim([0.01 * np.min(S[1:].value), 100 * np.max(S[1:].value)])
    plt.grid(which="major", alpha=0.8)
    plt.grid(which="minor", alpha=0.2)
    plt.legend()
    plt.show()


def _plot_cross_correlation(psrs, cor):
    npsrs = len(psrs)
    psr_names = [psr.name for psr in psrs]
    plt.imshow(cor, norm=colors.CenteredNorm(), cmap="RdBu", origin="lower")
    plt.title("Cross Correlation")
    plt.xticks(ticks=np.arange(npsrs), labels=psr_names, rotation=90)
    plt.xticks(ticks=np.arange(-0.5, npsrs, 1), minor=True)
    plt.yticks(ticks=np.arange(npsrs), labels=psr_names)
    plt.yticks(ticks=np.arange(-0.5, npsrs, 1), minor=True)
    plt.grid(which="minor", color="w", linestyle="-", linewidth=2)
    plt.tick_params(which="minor", bottom=False, left=False)
    plt.colorbar()
    plt.show()


def _plot_noise_weighted_cross_correlation(psrs, rho_matrix):
    npsrs = len(psrs)
    psr_names = [psr.name for psr in psrs]
    plt.imshow(rho_matrix, norm=colors.CenteredNorm(), cmap="RdBu", origin="lower")
    plt.title("Noise Weighted Cross Correlation")
    plt.xticks(ticks=np.arange(npsrs), labels=psr_names, rotation=90)
    plt.xticks(ticks=np.arange(-0.5, npsrs, 1), minor=True)
    plt.yticks(ticks=np.arange(npsrs), labels=psr_names)
    plt.yticks(ticks=np.arange(-0.5, npsrs, 1), minor=True)
    plt.grid(which="minor", color="w", linestyle="-", linewidth=2)
    plt.tick_params(which="minor", bottom=False, left=False)
    plt.colorbar()
    plt.show()


def _plot_rho(xi, rho, sig, OptStat):
    xii = np.linspace(0, 180, 181)
    ORF = np.array([hd_func(x * np.pi / 180.0) for x in xii])

    isort = np.argsort(xi)
    plt.errorbar(
        xi[isort] * 180 / np.pi,
        rho[isort],
        yerr=sig[isort],
        ls="None",
        marker="o",
        ms=3,
        color="r",
        mec="k",
        ecolor="k",
        capsize=2,
        label=r"$\rho_{ab}",
    )
    plt.plot(xii, ORF * OptStat, ls="--", color="k", label=r"$\hat{A}^2 \times$ ORF")
    plt.title("Noise Weighted Cross Correlation")
    plt.xlabel(r"$\xi$ [$\degree$]")
    plt.ylabel(r"$\rho_{ab}$ [$\hat{A}^2 \times$ ORF]")
    plt.xlim([0.0, 180.0])
    plt.grid(alpha=0.3)
    plt.show()


def _plot_network(G):
    edges, weights = zip(*nx.get_edge_attributes(G, "weight").items())
    weights = tuple([(1 + abs(x)) ** 2 for x in weights])
    d = nx.degree(G)
    nodelist, node_sizes = zip(*dict(d).items())
    positions = nx.circular_layout(G)
    plt.figure(figsize=(8, 8))
    nx.draw_networkx_nodes(
        G,
        positions,
        node_color="#DA70D6",
        nodelist=nodelist,
        node_size=tuple([2 * 19**2 for x in node_sizes]),
        alpha=0.8,
    )
    nx.draw_networkx_labels(G, positions, font_size=8, font_family="sans-serif")
    edge_colour = plt.cm.PuRd
    nx.draw_networkx_edges(
        G,
        positions,
        edgelist=edges,
        style="solid",
        width=weights,
        edge_color=weights,
        edge_cmap=edge_colour,
        edge_vmin=min(weights),
        edge_vmax=max(weights),
    )
    plt.axis("off")
    plt.title("Correlation Network")
    plt.show()


def _plot_network_measures(thresholds=None, SNRs=None, measures=None, mask=None):
    SNRs_med = np.zeros(len(measures))
    for idx, signal in enumerate(measures.keys()):
        SNRs_med[idx] = np.median(SNRs[signal])
    isort = np.argsort(SNRs_med)
    SNRs_med = SNRs_med[isort]

    plt.figure(figsize=(12, 8))
    signals = np.array(list(measures.keys()))
    signals = signals[isort]

    for signal in signals:
        plt.hist(SNRs[signal], bins=20, histtype="step", linewidth=2, label=signal)
    plt.xscale("log")
    plt.title("Signal to Noise Ratio")
    plt.xlabel("SNR")
    plt.ylabel("Number of realizations")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()

    for measure in all_network_measures:
        if thresholds is not None:
            plt.figure(figsize=(12, 8))
            for signal in signals:
                measure_avg = np.nanmean(measures[signal][measure], axis=0)
                measure_std = np.nanstd(measures[signal][measure], axis=0)
                plt.errorbar(thresholds, measure_avg, yerr=measure_std, label=signal, marker="o", mec="k", capsize=2)
            plt.title(measure)
            plt.xlabel("threshold")
            plt.ylabel(measure)
            plt.legend()
            plt.show()

            measure_th_avg = np.zeros((len(thresholds), len(signals)))
            measure_th_std = np.zeros((len(thresholds), len(signals)))
            for idx, signal in enumerate(signals):
                measure_avg = np.nanmean(measures[signal][measure], axis=0)
                measure_std = np.nanstd(measures[signal][measure], axis=0)
                for th_idx, th in enumerate(thresholds):
                    measure_th_avg[th_idx, idx] = measure_avg[th_idx]
                    measure_th_std[th_idx, idx] = measure_std[th_idx]
            plt.figure(figsize=(12, 8))
            for th_idx in mask:
                plt.errorbar(
                    SNRs_med,
                    measure_th_avg[th_idx, :],
                    yerr=measure_th_std[th_idx, :],
                    label=str(thresholds[th_idx]),
                    marker="o",
                    mec="k",
                    capsize=2,
                )
            plt.title(measure)
            plt.xscale("log")
            plt.xlabel("SNR")
            plt.ylabel(measure)
            plt.legend()
            plt.show()
        else:
            measure_avg = np.zeros(len(signals))
            measure_std = np.zeros(len(signals))
            for idx, signal in enumerate(signals):
                measure_avg[idx] = np.nanmean(measures[signal][measure])
                measure_std[idx] = np.nanstd(measures[signal][measure])
            plt.figure(figsize=(12, 8))
            plt.errorbar(SNRs_med, measure_avg, yerr=measure_std, marker="o", mec="k", capsize=2)
            plt.title(measure)
            plt.xscale("log")
            plt.xlabel("SNR")
            plt.ylabel(measure)
            plt.show()
