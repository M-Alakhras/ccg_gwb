# plotter.py
"""
plotting tools.
"""

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

from ccg_gwb.celestial_map import _pulsars_radec, _pulsars_separations
from ccg_gwb.multi_pulsar_tools import _calculate_cross_correlation, _calculate_noise_weighted_cross_correlation
from ccg_gwb.one_pulsar_tools import _calculate_power_spectrum


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


def _plot_power_spectrum(psr):
    f, S = _calculate_power_spectrum(psr.toas, psr.residuals)
    plt.loglog(f, S)
    plt.title("Power Spectrum")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel(r"Timing residuals Power [$\mu s^2$]")
    plt.ylim([0.01 * np.min(S[1:].value), 100 * np.max(S[1:].value)])
    plt.grid(which="major", alpha=0.8)
    plt.grid(which="minor", alpha=0.2)
    plt.show()


def _plot_cross_correlation(psrs):
    npsrs = len(psrs)
    psrname = [psr.name for psr in psrs]
    correlation_matrix = _calculate_cross_correlation(psrs)
    plt.imshow(correlation_matrix, norm=colors.CenteredNorm(), cmap="RdBu", origin="lower")
    plt.title("Cross Correlation")
    plt.xticks(ticks=np.arange(npsrs), labels=psrname, rotation=90)
    plt.xticks(ticks=np.arange(-0.5, npsrs, 1), minor=True)
    plt.yticks(ticks=np.arange(npsrs), labels=psrname)
    plt.yticks(ticks=np.arange(-0.5, npsrs, 1), minor=True)
    plt.grid(which="minor", color="w", linestyle="-", linewidth=2)
    plt.tick_params(which="minor", bottom=False, left=False)
    plt.colorbar()
    plt.show()


def _plot_ORF(xi, ORF):
    plt.plot(xi * 180 / np.pi, ORF, ls="None", marker="o", ms=3, color="r", mec="k")
    plt.title("Overlap Reduction Function")
    plt.xlabel(r"$\xi$ [$\degree$]")
    plt.ylabel("ORF [dimensionless]")
    plt.xlim([0.0, 180.0])
    plt.ylim([-0.5, 1.0])
    plt.grid(alpha=0.3)
    plt.show()


def _plot_rho(xi, rho, sig):
    plt.errorbar(
        xi * 180 / np.pi, rho, yerr=sig, ls="None", marker="o", ms=3, color="r", mec="k", ecolor="k", capsize=2
    )
    plt.title("Noise Weighted Cross Correlation")
    plt.xlabel(r"$\xi$ [$\degree$]")
    plt.ylabel(r"$\rho_{ab}$ [$\hat{A}^2 \times$ ORF]")
    plt.xlim([0.0, 180.0])
    plt.grid(alpha=0.3)
    plt.show()


def _plot_noise_weighted_cross_correlation(psrs, rho):
    npsrs = len(psrs)
    psrname = [psr.name for psr in psrs]
    rho_matrix = _calculate_noise_weighted_cross_correlation(psrs, rho)
    plt.imshow(rho_matrix, norm=colors.CenteredNorm(), cmap="RdBu", origin="lower")
    plt.title("Noise Weighted Cross Correlation")
    plt.xticks(ticks=np.arange(npsrs), labels=psrname, rotation=90)
    plt.xticks(ticks=np.arange(-0.5, npsrs, 1), minor=True)
    plt.yticks(ticks=np.arange(npsrs), labels=psrname)
    plt.yticks(ticks=np.arange(-0.5, npsrs, 1), minor=True)
    plt.grid(which="minor", color="w", linestyle="-", linewidth=2)
    plt.tick_params(which="minor", bottom=False, left=False)
    plt.colorbar()
    plt.show()
