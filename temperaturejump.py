""" Python script for temperature jump of simulation and observation """

import os
import argparse
from collections import OrderedDict
import pandas
import numpy
import scipy
from scipy import ndimage
from scipy.ndimage import filters
from astropy.io import fits
import matplotlib
from matplotlib import pyplot
import matplotlib.gridspec as gridspec

import peakutils
from deco import concurrent, synchronized

from ioparser import parse_gadget_parms
from ioparser import SimulationOutputParser
from cosmology import CosmologyCalculator
cc = CosmologyCalculator(z=0.0562)
arcsec2kpc = cc.kpc_DA
kB_kev_per_kelvin = 8.6173427909e-08
from macro import *
from plotsettings import PlotSettings
style = PlotSettings()


def find_centroids(sim, snapnr, method="DMhist", debug=False):
    """ Find the centroids of the simulated clusters.
        @param sim   : SimulationOutputParser instance
        @param snapnr: snapshot number to use
        @param method: string to choose method, options: "DMrho", "DMmedian"
        @param xlen  : when using DMmedian, also provide xlen
        @return      : centroids of both clusters in pixelvalues

    """

    if method == "DMhist":
        print "WARNING: DMhist method gives incorrect position if",
        print "barycenter is set incorrectly in P-Smac2!"
        if "xlen" not in dir(sim):
            print "ERROR: xlen and ylen not defined for simulation. Set by hand"
            return None
        sim.read_gadget(snapnr)
        if "pixelscale" not in dir(sim):
            sim.pixelscale = sim.parms.get("BoxSize")/sim.xlen
            print sim.pixelscale
        # NumericalCluster does not place particles in datamodel for two clusters
        # We are interested in the dark matter only anyway
        dmpos = sim.numerical.raw_data.pos[sim.numerical.raw_data.Ngas:
                                           sim.numerical.raw_data.N]
        xdm = dmpos[:,0]
        ydm = dmpos[:,1]
        zdm = dmpos[:,2]

        xhist, xbins = numpy.histogram(xdm, bins=sim.xlen)
        xcenters = xbins[:-1] + 0.5 * (xbins[1:] - xbins[:-1])
        yhist, ybins = numpy.histogram(ydm, bins=sim.ylen)
        ycenters = ybins[:-1] + 0.5 * (ybins[1:] - ybins[:-1])

        # Now find the peaks in the histogram of the dark matter density
        xpeaks = peakutils.indexes(xhist, thres=0.5)
        ypeaks = peakutils.indexes(yhist)

        found_peaks = False
        if len(xpeaks) == 2 and len(ypeaks) == 1:
            # Further optimize peakfinding by interpolating
            try:
                xpeaks = peakutils.interpolate(range(0, len(xhist)), xhist, ind=xpeaks)
                ypeaks = peakutils.interpolate(range(0, len(yhist)), yhist, ind=ypeaks)
            except RuntimeError as err:
                if "Optimal parameters not found: Number of calls to function has reached" in str(err):
                    print "peakutils.interpolate broke"
                else:
                    raise
            found_peaks = True
            CygA = xpeaks[0], ypeaks[0]
            CygB = xpeaks[1], ypeaks[0]
            distance = numpy.sqrt((CygA[0]-CygB[0])**2+(CygA[1]-CygB[1])**2)*sim.pixelscale
            print "Success: found 2 xpeaks, and 1 ypeak!"
            print "  CygA: (x, y) = {0}".format(CygA)
            print "  CygB: (x, y) = {0}".format(CygB)
            print "  distance     = {0:.2f}".format(distance)
            print
            return xpeaks[0], xpeaks[1], ypeaks[0]
        else:
            distance = numpy.nan

    # Find core separation by splitting the box in half, in x-direction to find
    # both haloes along the merger axis, and in y-direction to find the impact
    # parameter, should any be used. Then find the core separation by calculating
    # the median value of the dark matter particle positions, under the
    # assumption that the gas halo does not lag behind the DM halo.
    if method == "DMmedian":
        print "WARNING: DMmedian method gives incorrect position if",
        print "barycenter is set incorrectly in P-Smac2!"
        if "xlen" not in dir(sim):
            print "ERROR: xlen and ylen not defined for simulation. Set by hand"
            return None
        if debug: print "Reading dm particle positions..."
        sim.read_gadget(snapnr)
        if "pixelscale" not in dir(sim):
            sim.pixelscale = sim.parms.get("BoxSize")/sim.xlen
        # NumericalCluster does not place particles in datamodel for two clusters
        # We are interested in the dark matter only anyway
        dmpos = sim.numerical.raw_data.pos[sim.numerical.raw_data.Ngas:
                                           sim.numerical.raw_data.N]
        xdm = dmpos[:,0]
        ydm = dmpos[:,1]
        zdm = dmpos[:,2]

        boxhalf = sim.pixelscale*sim.xlen/2
        left = xdm[numpy.where(xdm < boxhalf)]
        right = xdm[numpy.where(xdm > boxhalf)]
        bottom = ydm[numpy.where(ydm < boxhalf)]
        top = ydm[numpy.where(ydm > boxhalf)]

        x1 = numpy.median(left)/sim.pixelscale
        x2 = numpy.median(right)/sim.pixelscale
        y = numpy.median(numpy.concatenate((bottom,top)))/sim.pixelscale
        x1 += (5052.37-4927.44)/sim.pixelscale * 2
        x2 += (5052.37-4927.44)/sim.pixelscale * 2

        return x1, x2, y

    # Find core separation by using a peakfinder to measure the distance between
    # the summed dark matter density output of P-Smac2, under the assumption
    # that the gas halo does not lag behind the DM halo.
    if method == "DMrho":
        sim.read_smac("dm-density_projection-z.fits.fz")

        xpix = numpy.arange(0, sim.xlen, 1, dtype=numpy.int)
        ypix = numpy.arange(0, sim.ylen, 1, dtype=numpy.int)

        # The sum of the snapshot in the fits cube of the darkmatter density
        # along the x-axis/y-axis also gives cluster centroids
        xsum = numpy.sum(sim.rhodmdata[snapnr], axis=0)
        ysum = numpy.sum(sim.rhodmdata[snapnr], axis=1)
        pyplot.figure()
        pyplot.plot(xpix, xsum, "ro", label="x")
        pyplot.plot(ysum, ypix, label="y")
        pyplot.legend()
        pyplot.show()


def create_panda(xlen, ylen, xc, yc, r, a1, a2):
    """ Create mask of spherical wedge with radius r between angles a1 and a2
            for a numpy array of shape xlen, ylen
        @param xlen, ylen: shape of the array to mask
        @param xy, yc    : coordinates to use as center for the wedges
        @param r         : radius of the spherical shell
        @param a1, a2    : angles, where horizontal = 0; counterclockwise;
                           values given must be in degrees
        @return          : numpy mask containing indices
                           obtaining the masked values works as array[mask]
    """

    # Thickness of the spherical wedge shell section is three pixels now
    dr = 1

    # Center mask shape of (xlen, ylen) around (xc, yc)
    y,x = numpy.ogrid[-yc:ylen-yc, -xc:xlen-xc]

    # Convert degrees to radians
    a1 *= numpy.pi/180
    a2 *= numpy.pi/180

    # Some trial and error angle magic. Be careful with arctan2 domain!
    a1 = a1%(2*numpy.pi)
    a2 = a2%(2*numpy.pi)

    angles = numpy.arctan2(y,x)
    angles += 2*numpy.pi
    angles = angles % (2*numpy.pi)

    if a2 < a1:
        anglemask = ((a1 <= angles) | (angles <= a2))
    else:
        anglemask = ((a1 <= angles) & (angles <= a2))

    return (x**2 + y**2 <= r**2) & (x**2 + y**2 >= (r-dr)**2) & anglemask


def add_suzaku():
    # Eyeballed from Sarazin, Finoguenov & Wik (2012) Fig. 4. Black points.
    r_arcmin_Suzaku = numpy.array([-2, -0.1, 2, 4, 5.9, 8, 10, 12.3])
    r_error_plus = numpy.array([1, 1, 1, 1, 1, 1, 1, 1.5])
    r_error_min = numpy.array([1, 1, 1, 1, 1, 1, 1, 1.5])
    T_kev_Suzaku = numpy.array([5, 5.5, 6.2, 7, 8.2, 8.5, 7.2, 6.4])
    T_error_plus = numpy.array([1, 0.3, 1, 2, 1.5, 1.6, 1.1, 0.7])
    T_error_min = numpy.array([1, 0.7, 1, 2, 1, 1.3, 1.3, 0.8])
    pyplot.errorbar(r_arcmin_Suzaku*60, T_kev_Suzaku,
        xerr=[r_error_min*60, r_error_plus*60],
        yerr=[T_error_min, T_error_plus], label="Suzaku",
        c="k", fmt="+", elinewidth=4, capsize=0)

    # Eyeballed from Sarazin, Finoguenov & Wik (2012) Fig. 4. Grey points.
    r_arcmin_Suzaku = numpy.array([0.1, 6.1, 11.1])
    r_error_plus = numpy.array([3.5, 1.8, 3])
    r_error_min = numpy.array([3.1, 2.8, 3.2])
    T_kev_Suzaku = numpy.array([5.6, 8.8, 7])
    T_error_plus = numpy.array([0.3, 0.7, 0.5])
    T_error_min = numpy.array([0.4, 0.7, 0.4])
    pyplot.errorbar(r_arcmin_Suzaku*60, T_kev_Suzaku,
        xerr=[r_error_min*60, r_error_plus*60],
        yerr=[T_error_min, T_error_plus], label="Suzaku (average)",
        c="grey", fmt="+", elinewidth=4, capsize=0)


def add_chandra_average():
    # Average profile is over all wedges but the merger wedge. Obtained from
    # /scratch/martyndv/cygnus/combined/spectral/maps/radial/sn100_new/plots_cygA
    # Last edit of Martijn's fitfile was August 20, 2016 at 21:12
    fitfile = "data/20160829/cygA_1T_fixnH_fitresults.dat"
    # Last edit of Martijn's sbfile was August 20, 2016 at 21:10
    sbfile = "data/20160829/cygA_sb_sn100.dat"

    fitresults = pandas.read_csv(fitfile, delim_whitespace=True, skiprows=1)
    sbresults = pandas.read_csv(sbfile, delim_whitespace=True, skiprows=18)

    average_radii = 0.5*(sbresults[u"Radius1"] +
                         sbresults[u"Radius2"] ) * arcsec2kpc
    average_binsize = (sbresults[u"Radius2"] -
                       sbresults[u"Radius1"]) * arcsec2kpc
    average_kT = fitresults[u"kT1"]
    average_kT_loerr = numpy.abs(fitresults[u"kT1_loerr"])
    average_kT_hierr = fitresults[u"kT1_hierr"]
    # Fix error values where fit did not converge
    fail_indices = average_kT == average_kT_loerr
    average_kT_loerr[fail_indices] = 0.05*average_kT
    average_kT_hierr[fail_indices] = 0.05*average_kT

    pyplot.errorbar(average_radii, average_kT,
                    yerr=[average_kT_loerr, average_kT_hierr],
                    xerr=average_binsize/2, marker="o", ls="", c="b", ms=4,
                    alpha=0.2, elinewidth=2, label="Observation: Average")


def add_chandra_sector(merger=True, hot=True, cold=True, wedges=False):
    # Hot, cold and merger wedges sector radial profiles, T-fit to Chandra
    # observation. Data analysis was done by Martijn de Vries. Files obtained
    # from /scratch/martyndv/cygnus/combined/spectral/maps/sector/
    # Last edit of Martijn's sectorfile was August 19, 2016 at 00:06
    sectorfile = "data/20160829/cygnus_sector_fitresults.dat"
    #  Last edit of Martijn's sbfile was June 18, 2016 at 11:50
    sbfile = "data/20160829/cygnus_sn100_sbprofile.dat"

    sectorfit = pandas.read_csv(sectorfile, delim_whitespace=True, skiprows=1)
    sbprofile = pandas.read_csv(sbfile, delimiter="|")

    sectorfit_merger = sectorfit[0:149]
    sbprofile_merger = sbprofile[0:149]
    sectorfit_hot = sectorfit[149:308]
    sbprofile_hot = sbprofile[149:308]
    sectorfit_cold = sectorfit[308:372]
    sbprofile_cold = sbprofile[308:372]

    merger_radii = 0.5*(sbprofile_merger[u" Inner radius (arcsec) "] +
                        sbprofile_merger[u" Outer radius (arcsec) "] ) * arcsec2kpc
    merger_binsize = (sbprofile_merger[u" Outer radius (arcsec) "] -
                      sbprofile_merger[u" Inner radius (arcsec) "]) * arcsec2kpc
    merger_kT = sectorfit_merger[u"kT1"]
    merger_kT_loerr = numpy.abs(sectorfit_merger[u"kT1_loerr"])
    merger_kT_hierr = sectorfit_merger[u"kT1_hierr"]
    hot_radii = 0.5*(sbprofile_hot[u" Inner radius (arcsec) "] +
                         sbprofile_hot[u" Outer radius (arcsec) "] ) * arcsec2kpc
    hot_binsize = (sbprofile_hot[u" Outer radius (arcsec) "] -
                       sbprofile_hot[u" Inner radius (arcsec) "]) * arcsec2kpc
    hot_kT = sectorfit_hot[u"kT1"]
    hot_kT_loerr = numpy.abs(sectorfit_hot[u"kT1_loerr"])
    hot_kT_hierr = sectorfit_hot[u"kT1_hierr"]
    cold_radii = 0.5*(sbprofile_cold[u" Inner radius (arcsec) "] +
                      sbprofile_cold[u" Outer radius (arcsec) "] ) * arcsec2kpc
    cold_binsize = (sbprofile_cold[u" Outer radius (arcsec) "] -
                    sbprofile_cold[u" Inner radius (arcsec) "]) * arcsec2kpc
    cold_kT = sectorfit_cold[u"kT1"]
    cold_kT_loerr = numpy.abs(sectorfit_cold[u"kT1_loerr"])
    cold_kT_hierr = sectorfit_cold[u"kT1_hierr"]
    if merger:
        pyplot.errorbar(merger_radii, merger_kT, xerr=merger_binsize/2,
                        yerr=[merger_kT_loerr, merger_kT_hierr],
                        marker="o", ls="", c="g", ms=4, alpha=0.2,
                        elinewidth=2, label="Chandra: Merger")
    if hot:
        pyplot.errorbar(hot_radii, hot_kT, xerr=hot_binsize/2,
                        yerr=[hot_kT_loerr, hot_kT_hierr],
                        marker="o", ls="", c="r", ms=6, alpha=0.2,
                        elinewidth=2, label="Chandra: Hot")
    if cold:
        pyplot.errorbar(cold_radii, cold_kT, xerr=cold_binsize/2,
                        yerr=[cold_kT_loerr, cold_kT_hierr],
                        marker="o", ls="", c="purple", ms=6, alpha=0.2,
                        elinewidth=2, label="Chandra: Cold")

    # Plot that shows the wedges and the extraction radii
    if wedges:
        fig = pyplot.gcf()
        wedges = [merger_radii, hot_radii, cold_radii, hot_radii]
        colors = ["g", "r", "purple", "r"]
        for lim in (50, 500, 2500):
            pyplot.figure(2, figsize=(12,12))
            for i, (th1, th2) in enumerate(zip([6, 96, 225, 315],
                                               [96, 225, 315, 366])):
                x = numpy.cos(numpy.linspace(2*numpy.pi/360*th1,
                                             2*numpy.pi/360*th2, 100))
                y = numpy.sin(numpy.linspace(2*numpy.pi/360*th1,
                                             2*numpy.pi/360*th2, 100))
                for r in wedges[i]:
                    pyplot.plot(r*x, r*y, c=colors[i])
            pyplot.xlim(-lim, lim)
            pyplot.ylim(-lim, lim)
            pyplot.savefig("out/chandra_adaptivebinning_regions_{0}.png"\
                .format(lim), dpi=300)
            pyplot.close(2)
        pyplot.scf(fig)


def add_ds9_wedge(sim, vel, i):
    # data files creates using ds9, region/panda centered on CygA,
    # MergerAxis has angles -45 to 45; NotMergerAxis has angles 45 to -45
    # Annuli in range(0, 200, 42), centered at (1004, 1024) in cubenumber 78

    timestamp = "data/ds9wedges/"+sim.timestamp
    MergerAxis = pandas.read_csv(timestamp+"_MergerAxis.dat", delimiter=" ",
        header=None, names=["Radius", "Temperature", "TemperatureErr"])

    pyplot.plot(MergerAxis["Radius"]*sim.pixelscale,
                MergerAxis["Temperature"]*kB_kev_per_kelvin, c="k",
                lw=4, label="Simulation {0:1.1f}: Merger".format(vel))

    NotMergerAxis = pandas.read_csv(timestamp+"_NotMergerAxis.dat", delimiter=" ",
        header=None, names=["Radius", "Temperature", "TemperatureErr"])
    pyplot.plot(NotMergerAxis["Radius"]*sim.pixelscale,
                NotMergerAxis["Temperature"]*kB_kev_per_kelvin, c="k",
                ls="dashed", lw=4, label="Simulation {0:1.1f}: Average".format(vel))


def add_own_wedge(sim, snapnr, i, wedges=False):
    """ Add own wedge using the create_panda method """

    # Quick 'n' Dirty find peaks in X-ray surface brightness
    # xsumleft = numpy.sum(snap[:,0:sim.xlen/2+200], axis=0)
    # xsumright = numpy.sum(snap[:,sim.xlen/2+200:], axis=0)
    # ysum = numpy.sum(snap, axis=1)
    # xpeakleft = peakutils.indexes(xsumleft)[0]
    # xpeakright = peakutils.indexes(xsumright)[0] + sim.xlen/2+200
    # ypeak = peakutils.indexes(ysum)[0]

    # x1, x2, y = find_centroids(sim, snapnr, method="DMmedian")
    peaks = {0: (985, 1217, 1024), 3: (1004, 1134, 1024), 8: (1000, 1138, 1024)}

    xpeak1, xpeak2, ypeak = peaks.get(i)

    radii = numpy.arange(1, 200, 200./42)
    quiescent_temperature = numpy.zeros(len(radii))
    quiescent_temperature_std = numpy.zeros(len(radii))
    merger_temperature = numpy.zeros(len(radii))
    merger_temperature_std = numpy.zeros(len(radii))

    fig = pyplot.gcf()
    if wedges:
        snap = sim.temdata[snapnr]
        fig2 = pyplot.figure()
        pyplot.imshow(snap*kB_kev_per_kelvin, vmin=0.1, vmax=9,
                      origin="lower", cmap="afmhot")
        pyplot.plot(xpeak1, ypeak, "wo", mfc=None, ms=5)
        pyplot.plot(xpeak2, ypeak, "wo", mfc=None, ms=5)
        pyplot.xlim(xpeak1-300, xpeak2+300)
        pyplot.ylim(ypeak-300, ypeak+300)

    for j, r in enumerate(radii):
        print r
        quiescent_mask = create_panda(sim.xlen, sim.ylen, xpeak1, ypeak, r, 45, -45)
        quiescent_temperature[j] = numpy.median(snap[quiescent_mask])
        quiescent_temperature_std[j] = numpy.std(snap[quiescent_mask])
        if wedges:
            y, x = numpy.where(quiescent_mask)
            pyplot.scatter(x, y, s=1, c="r", edgecolor="face", alpha=1)
        merger_mask = create_panda(sim.xlen, sim.ylen, xpeak1, ypeak, r, -45, 45)
        merger_temperature[j] = numpy.median(snap[merger_mask])
        merger_temperature_std[j] = numpy.std(snap[merger_mask])
        if wedges:
            y, x = numpy.where(merger_mask)
            pyplot.scatter(x, y, s=1, c="w", edgecolor="face", alpha=1)
            pyplot.xlabel("x [pixel]")
            pyplot.ylabel("y [pixel]")
            pyplot.gca().set_aspect("equal")
            pyplot.tight_layout()
            pyplot.savefig(sim.outdir+"extraction.png")
    pyplot.figure(fig.number)

    radii *= sim.pixelscale
    quiescent_temperature*=kB_kev_per_kelvin
    quiescent_temperature_std*=kB_kev_per_kelvin
    merger_temperature*=kB_kev_per_kelvin
    merger_temperature_std*=kB_kev_per_kelvin
    pyplot.errorbar(radii, quiescent_temperature, yerr=quiescent_temperature_std,
                    c="k", lw=4, label="quiescent")
    pyplot.errorbar(radii, merger_temperature, yerr=merger_temperature_std,
                    c="k", lw=4, label="merger")



def plot_temperature_profile(xscale="log"):
    #if "dt" not in dir(self):
    #    self.set_timebetsnap()

    timestamps = ["20160819T2322", "20160819T2357", "20160820T0032",
                  "20160820T0107", "20160820T0142", "20160820T0218",
                  "20160820T0252", "20160820T0328", "20160820T0403",
                  "20160820T0438", "20160820T0513"]
    bestsnaps = numpy.array([78, 77, 68, 60, 53, 47, 43, 39, 35, 33])
    ZeroEOrbitFrac = numpy.arange(0.0, 1.1, 0.1)
    MergerVelocity = numpy.array([])  # km/s


    for i, timestamp in enumerate(timestamps):
        if i not in [0, 3, 8]:
            continue
        print i

        pyplot.figure()

        # add_suzaku()
        add_chandra_average()
        add_chandra_sector(merger=True, hot=True, cold=True)

        sim = SimulationOutputParser("/Volumes/SURFlisa", timestamp)
        sim.read_smac("temperature-emission-weighted_projection-z.fits.fz")

        add_ds9_wedge(sim, ZeroEOrbitFrac[i], i)
        add_own_wedge(sim, bestsnaps[i], i, wedges=True)

        pyplot.axvline(72, ls="dashed", c="k")
        #pyplot.xlabel("Distance Along Merger Axis [kpc]")
        #pyplot.ylabel("Temperature [K]")
        pyplot.xlabel("Radius [kpc]")
        pyplot.ylabel("kT [keV]")
        pyplot.xlim(5, 1000)
        pyplot.ylim(2, 11)
        if xscale=="log":
            pyplot.xscale("log")
            pyplot.xticks([10, 100, 1000], [r"$10^1$", r"$10^2$", r"$10^3$"])
        pyplot.legend(loc="upper left", fontsize=12)
        pyplot.tight_layout()
        pyplot.savefig(sim.outdir+"temperature.png", dpi=300)
    pyplot.show()


def new_argument_parser():
    parser = argparse.ArgumentParser(
        description="Plot Temperature Jump")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-t", "--timestamp", dest="timestamp", nargs=1,
        help="string of the Simulation ID")

    return parser


if __name__ == "__main__":
    # arguments = new_argument_parser().parse_args()
    # sim = SimulationOutputParser("/Volumes/SURFlisa", arguments.timestamp[0])

    # sim.xlen = 2048
    # sim.ylen = 2048
    # find_centroids(sim, 0, method="DMrho")

    #sim.read_smac("temperature-emission-weighted_projection-z.fits.fz")
    plot_temperature_profile()
