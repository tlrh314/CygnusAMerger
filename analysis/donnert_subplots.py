"""
File: donnert_subplots.py
Author: Timo L. R. Halbesma <timohalbesma@gmail.com>
Date created: Mon Apr 18, 2016 12:33 pm
Last modified: Thu Apr 21, 2016 01:46 PM

"""

import numpy
from matplotlib import pyplot
# pyplot.rcParams.update({"font.size": 18})

from amuse.units import units
from amuse.units import constants
from amuse.units.quantities import new_quantity
from amuse.units.quantities import VectorQuantity
from amuse.units import nbody_system
from amuse import datamodel

import amuse.plot as amuse_plot

def amuse_nth_root(quant, n):
    # Simply telling AMUSE quant**(1./3) breaks the units :-(
    return new_quantity((quant.number)**(1./n), (quant.unit ** (1./n)).to_simple_form())


def donnert_2014_figure1(analytical_P500=False):
    # Well, I will just eyeball the value of rho_0 in Figure 1 in Donnert (2014) and adopt this value, I suppose...
    # disturbed = InitialCluster(rho_0=(9e-27 | units.g/units.cm**3),
    #                            plot=False, do_print=False,
    #                            disturbed_cluster=True, cool_core_cluster=False)
    # coolcore = InitialCluster(rho_0=(3e-25 | units.g/units.cm**3),
    #                           plot=False, do_print=False,
    #                           disturbed_cluster=False, cool_core_cluster=True)
    # rickersarazin = InitialCluster(plot=False, do_print=False, disturbed_cluster=False, cool_core_cluster=False)

    r = VectorQuantity.arange(units.kpc(1), units.kpc(10000), units.parsec(100))

    #fig, ((ax1, ax2), (ax3, ax4)) = pyplot.subplots(2, 2, figsize=(20, 20), dpi=500)
    fig, ((ax1, ax2), (ax3, ax4)) = pyplot.subplots(2, 2)

    # Plot the density situation
    pyplot.sca(ax1)
    amuse_plot.hist(r, bins=int(numpy.sqrt(len(r))), label="Hanky")
    # amuse_plot.loglog(coolcore.r, coolcore.rho_gas.as_quantity_in(units.g / units.cm**3),
    #                   c='k', ls='dotted', label="Gas, cool core")
    # amuse_plot.loglog(disturbed.r, disturbed.rho_gas.as_quantity_in(units.g / units.cm**3),
    #                   c='k', ls='dashed', label="Gas, disturbed")
    # amuse_plot.loglog(disturbed.r, disturbed.rho_dm.as_quantity_in(units.g / units.cm**3),
    #                   c='k', ls='solid', label="Dark Matter")
    amuse_plot.ylabel(r'$\rho$')
    amuse_plot.xlabel(r'$r$')
    pyplot.legend(loc=3, frameon=False, fontsize=12)

    ax1.set_ylim(ymin=1e-28, ymax=1e-23)

    # pyplot.axvline(x=disturbed.r_200.value_in(units.kpc), lw=2, c='k')
    # pyplot.text(disturbed.r_200.value_in(units.kpc), 7e-24, r'$r_{200}$', fontsize=12)
    # pyplot.axvline(x=disturbed.r_500.value_in(units.kpc), lw=2, c='k')
    # pyplot.text(disturbed.r_500.value_in(units.kpc), 7e-24, r'$r_{500}$', fontsize=12)

    # # a_hernq
    # intersect = disturbed.dm_density(disturbed.a)
    # ymin, ymax = ax1.get_ylim()
    # ymin = (numpy.log(intersect.value_in(units.g/units.cm**3)/ymin)) / (numpy.log(ymax/ymin))
    # pyplot.axvline(x=disturbed.a.value_in(units.kpc), ymin=ymin, ymax=1, lw=2, c='k')
    # pyplot.text(disturbed.a.value_in(units.kpc), 5e-24, r'$a_{\rm Hernq}$', fontsize=12)

    # # r_core,dist
    # intersect = disturbed.gas_density(disturbed.r_c)
    # ymin, ymax = ax1.get_ylim()
    # ymin = (numpy.log(intersect.value_in(units.g/units.cm**3)/ymin)) / (numpy.log(ymax/ymin))
    # pyplot.axvline(x=disturbed.r_c.value_in(units.kpc), ymin=ymin, ymax=1, lw=2, ls='--', c='k')
    # pyplot.text(disturbed.r_c.value_in(units.kpc), 5e-24, r'$r_{\rm core,dist}$', fontsize=12)

    # # r_core,cc
    # intersect = coolcore.gas_density(coolcore.r_c)
    # ymin, ymax = ax1.get_ylim()
    # ymin = (numpy.log(intersect.value_in(units.g/units.cm**3)/ymin)) / (numpy.log(ymax/ymin))
    # pyplot.axvline(x=coolcore.r_c.value_in(units.kpc), ymin=ymin, ymax=1, lw=2, ls=':', c='k')
    # pyplot.text(coolcore.r_c.value_in(units.kpc), 5e-24, r'$r_{\rm core,cc}$', fontsize=12)

    # Plot the mass situation
    pyplot.sca(ax2)
    amuse_plot.hist(r, bins=int(numpy.sqrt(len(r))), label="Hanky")
    # amuse_plot.loglog(coolcore.r, coolcore.M_gas_below_r.as_quantity_in(units.MSun),
    #                   c='k', ls='dotted', label="Gas, cool core")
    # amuse_plot.loglog(disturbed.r, disturbed.M_gas_below_r.as_quantity_in(units.MSun),
    #                   c='k', ls='dashed', label="Gas, disturbed")
    # amuse_plot.loglog(disturbed.r, disturbed.M_dm_below_r.as_quantity_in(units.MSun),
    #                   c='k', ls='solid', label="Dark Matter")
    amuse_plot.ylabel(r'$M(<r)$')
    amuse_plot.xlabel(r'$r$')
    pyplot.legend(loc=8, frameon=False, fontsize=12)

    ax2.set_ylim(ymin=1e10, ymax=5e15)

    # pyplot.axvline(x=disturbed.r_200.value_in(units.kpc), lw=2, c='k')
    # pyplot.text(disturbed.r_200.value_in(units.kpc), 3e15, r'$r_{200}$', fontsize=12)
    # pyplot.axvline(x=disturbed.r_500.value_in(units.kpc), lw=2, c='k')
    # pyplot.text(disturbed.r_500.value_in(units.kpc), 3e15, r'$r_{500}$', fontsize=12)

    # # a_hernq
    # intersect = disturbed.dm_cummulative_mass(disturbed.a)
    # ymin, ymax = ax2.get_ylim()
    # ymin = (numpy.log(intersect.value_in(units.MSun)/ymin)) / (numpy.log(ymax/ymin))
    # pyplot.axvline(x=disturbed.a.value_in(units.kpc), ymin=ymin, ymax=1, lw=2, c='k')
    # pyplot.text(disturbed.a.value_in(units.kpc), 2e15, r'$a_{\rm Hernq}$', fontsize=12)

    # # r_core,dist
    # intersect = disturbed.dm_cummulative_mass(disturbed.r_c)
    # ymin, ymax = ax2.get_ylim()
    # ymin = (numpy.log(intersect.value_in(units.MSun)/ymin)) / (numpy.log(ymax/ymin))
    # pyplot.axvline(x=disturbed.r_c.value_in(units.kpc), ymin=ymin, ymax=1, lw=2, ls='--', c='k')
    # pyplot.text(disturbed.r_c.value_in(units.kpc), 2e15, r'$r_{\rm core,dist}$', fontsize=12)

    # # r_core,cc
    # intersect = coolcore.dm_cummulative_mass(coolcore.r_c)
    # ymin, ymax = ax2.get_ylim()
    # ymin = (numpy.log(intersect.value_in(units.MSun)/ymin)) / (numpy.log(ymax/ymin))
    # pyplot.axvline(x=coolcore.r_c.value_in(units.kpc), ymin=ymin, ymax=1, lw=2, ls=':', c='k')
    # pyplot.text(coolcore.r_c.value_in(units.kpc), 2e15, r'$r_{\rm core,cc}$', fontsize=12)

    # Plot the temperature situation
    pyplot.sca(ax3)
    amuse_plot.hist(r, bins=int(numpy.sqrt(len(r))), label="Hanky")
    # amuse_plot.loglog(disturbed.r, disturbed.T_r.as_quantity_in(units.K),
    #                   c='k', ls='solid', label="disturbed")
    # amuse_plot.loglog(disturbed.r, disturbed.T_r_dm.as_quantity_in(units.K),
    #                   c='k', ls='dashdot', label="from DM, disturbed")
    # amuse_plot.loglog(disturbed.r, disturbed.T_r_gas.as_quantity_in(units.K),
    #                   c='k', ls='dashed', label="from Gas, disturbed")
    # amuse_plot.loglog(coolcore.r, coolcore.T_r.as_quantity_in(units.K),
    #                   c='k', ls='dotted', label="cool core")
    amuse_plot.ylabel(r'$T$')
    amuse_plot.xlabel(r'$r$')
    pyplot.legend(loc=10, frameon=False, fontsize=12)

    ax3.set_ylim(ymin=6e6, ymax=3e8)

    # pyplot.axvline(x=disturbed.r_200.value_in(units.kpc), lw=2, c='k')
    # pyplot.text(disturbed.r_200.value_in(units.kpc), 2.6e8, r'$r_{200}$', fontsize=12)
    # pyplot.axvline(x=disturbed.r_500.value_in(units.kpc), lw=2, c='k')
    # pyplot.text(disturbed.r_500.value_in(units.kpc), 2.6e8, r'$r_{500}$', fontsize=12)

    # pyplot.axvline(x=disturbed.a.value_in(units.kpc), lw=2, c='k')
    # pyplot.text(disturbed.a.value_in(units.kpc), 2.4e8, r'$a_{\rm Hernq}$', fontsize=12)

    # # r_core,dist
    # intersect = disturbed.temperature(disturbed.r_c)
    # ymin, ymax = ax3.get_ylim()
    # ymin = (numpy.log(intersect.value_in(units.K)/ymin)) / (numpy.log(ymax/ymin))
    # pyplot.axvline(x=disturbed.r_c.value_in(units.kpc), ymin=ymin, ymax=1, lw=2, ls='--', c='k')
    # pyplot.text(disturbed.r_c.value_in(units.kpc), 2.4e8, r'$r_{\rm core,dist}$', fontsize=12)

    # # r_core,cc
    # intersect = coolcore.temperature(coolcore.r_c)
    # ymin, ymax = ax3.get_ylim()
    # ymin = (numpy.log(intersect.value_in(units.K)/ymin)) / (numpy.log(ymax/ymin))
    # pyplot.axvline(x=coolcore.r_c.value_in(units.kpc), ymin=ymin, ymax=1, lw=2, ls=':', c='k')
    # pyplot.text(coolcore.r_c.value_in(units.kpc), 2.4e8, r'$r_{\rm core,cc}$', fontsize=12)

    # TODO: pressure situation
    pyplot.sca(ax4)
    amuse_plot.hist(r, bins=int(numpy.sqrt(len(r))), label="Hanky")
    colours=[(255./255, 127./255, 0./255),
             (152./255, 78./255, 163./255),
             (77./255, 175./255, 74./255),
             (52./255, 126./255, 184./255),
             (228./255, 26./255, 28./255)]

    # print "Warning, the pressure plots make little sense because for each M_DM is the same but M_200 differs"
    # for i, M_200 in enumerate(VectorQuantity([3e15, 1.5e15, 1e15, 0.5e15, 0.1e15], units.MSun)):
    #     disturbed = InitialCluster(rho_0=(9e-27 | units.g/units.cm**3), M_200=M_200,
    #                                plot=False, do_print=False,
    #                                disturbed_cluster=True, cool_core_cluster=False)
    #     coolcore = InitialCluster(rho_0=(3e-25 | units.g/units.cm**3), M_200=M_200,
    #                               plot=False, do_print=False,
    #                               disturbed_cluster=False, cool_core_cluster=True)
    #     if analytical_P500:
    #         amuse_plot.loglog(disturbed.r/disturbed.r_500, (disturbed.P_gas/disturbed.find_p500_analytically()),
    #                           c=colours[i], ls='solid', label=r'{0} $M_\odot$'.format(disturbed.M_200.value_in(units.MSun)))
    #         amuse_plot.loglog(coolcore.r/coolcore.r_500, (coolcore.P_gas/coolcore.find_p500_analytically()),
    #                           c=colours[i], ls='dashed')
    #     else:
    #         amuse_plot.loglog(disturbed.r/disturbed.r_500, (disturbed.P_gas/disturbed.P_500),
    #                           c=colours[i], ls='solid', label=r'{0} $M_\odot$'.format(disturbed.M_200.value_in(units.MSun)))
    #         amuse_plot.loglog(coolcore.r/coolcore.r_500, (coolcore.P_gas/coolcore.P_500),
    #                           c=colours[i], ls='dashed')

    amuse_plot.ylabel(r'$P(r)/P_{500}$')
    amuse_plot.xlabel(r'$r/r_{500}$')
    legend = pyplot.legend(loc=3, frameon=False, fontsize=16)
    # Set the color situation
    # for colour, text in zip(colours, legend.get_texts()):
    #     text.set_color(colour)

    ax4.set_ylim(ymin=1e-2, ymax=1e3)
    ax4.set_xticks((0.01, 0.10, 1.00))
    ax4.set_xticklabels(("0.01", "0.10", "1.00"))
    ax4.set_xlim(xmin=0.01, xmax=2.0)
    ax4.minorticks_on()
    ax4.tick_params('both', length=10, width=2, which='major')
    ax4.tick_params('both', length=5, width=1, which='minor')
    ax4.text(0.015, 4, "disturbed", color="grey", fontsize=15)
    ax4.text(0.02, 110, "cool cores", color="grey", fontsize=15)
    ax4.text(0.2, 50, "Arnaud et al. 2010", color="black", fontsize=15)

    # Set the xticks situation
    for ax in [ax1, ax2, ax3]:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks((10, 100, 1000))
        ax.set_xticklabels(("10", "100", "1000"))
        ax.set_xlim(xmin=10, xmax=5000)
        ax.minorticks_on()
        ax.tick_params('both', length=10, width=2, which='major')
        ax.tick_params('both', length=5, width=1, which='minor')

    # pyplot.savefig("../img/Donnert2014_Figure1_by_TLRH.pdf", format="pdf", dpi=1000)
    pyplot.show()


if __name__ == '__main__':
    donnert_2014_figure1()
