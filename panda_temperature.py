from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 18})
import numpy
import pandas


arcsec2kpc = 1.091
pyplot.figure(figsize=(12,9))

# First, data: everything but merger wedge
# --------------------------------------------------------------------------- #
# Average profile is over all wedges but the merger wedge. Obtained from
# /scratch/martyndv/cygnus/combined/spectral/maps/radial/sn100_new/plots_cygA
# Last edit of Martijn's files was August 20, 2016 at 21:1{0,2}
fitfile = "../MartijnStruis/plot_cygA/cygA_1T_fixnH_fitresults.dat"
sbfile = "../MartijnStruis/plot_cygA/cygA_sb_sn100.dat"

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
average_kT_loerr[fail_indices] = 0.1*average_kT
average_kT_hierr[fail_indices] = 0.1*average_kT

pyplot.errorbar(average_radii, average_kT, yerr=[average_kT_loerr, average_kT_hierr],
                xerr=average_binsize/2, marker="o", ls="", c="b", ms=4, alpha=0.2,
                elinewidth=2, label="Observation: Average")

# second, data: merger wedge
# --------------------------------------------------------------------------- #
# Hot, cold and merger wedges sector radial profiles, T-fit to Chandra observation
# Data analysis was done by Martijn de Vries. Files obtained from
# /scratch/martyndv/cygnus/combined/spectral/maps/sector/
# Last edit of Martijn's files was August 19, 2016 at 00:06 (fit)
#                              and June 18, 2016 at 11:50 (sb)
sectorfile = "../MartijnStruis/cygnus_sector_fitresults.dat"
sbfile = "../MartijnStruis/cygnus_sn100_sbprofile.dat"

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

# Plot that shows the wedges and the extraction radii
show_pizza = False
if show_pizza:
    wedges = [merger_radii, hot_radii, cold_radii, hot_radii]
    colors = ["g", "r", "purple", "r"]
    pyplot.figure(figsize=(12,12))
    for i, (th1, th2) in enumerate(zip([6, 96, 225, 315], [96, 225, 315, 366])):
        x = numpy.cos(numpy.linspace(2*numpy.pi/360*th1, 2*numpy.pi/360*th2, 100))
        y = numpy.sin(numpy.linspace(2*numpy.pi/360*th1, 2*numpy.pi/360*th2, 100))
        for r in wedges[i]:
            pyplot.plot(r*x, r*y, c=colors[i])
    pyplot.xlim(-500, 500)
    pyplot.ylim(-500, 500)
    pyplot.show()

pyplot.errorbar(merger_radii, merger_kT, yerr=[merger_kT_loerr, merger_kT_hierr],
                xerr=merger_binsize/2, marker="o", ls="", c="g", ms=4, alpha=0.2,
                elinewidth=2, label="Observation: Merger Axis")
#pyplot.errorbar(hot_radii, hot_kT, yerr=[hot_kT_loerr, hot_kT_hierr],
#                xerr=hot_binsize/2, marker="o", ls="", c="r", ms=6,
#                elinewidth=2, label="hot")
#pyplot.errorbar(cold_radii, cold_kT, yerr=[cold_kT_loerr, cold_kT_hierr],
#                xerr=cold_binsize/2, marker="o", ls="", c="purple", ms=6,
#                elinewidth=2, label="cold")
pyplot.axvline(72, ls="dashed", c="k")

# Fourth, Sarazin Suzaku observation
# --------------------------------------------------------------------------- #
suzaku = True
if suzaku:
    # Eyeballed from Sarazin, Finoguenov & Wik (2012) Fig. 4. Black points.
    r_arcmin_Suzaku = numpy.array([-2, -0.1, 2, 4, 5.9, 8, 10, 12.3])
    r_error_plus = numpy.array([1, 1, 1, 1, 1, 1, 1, 1.5])
    r_error_min = numpy.array([1, 1, 1, 1, 1, 1, 1, 1.5])
    T_kev_Suzaku = numpy.array([5, 5.5, 6.2, 7, 8.2, 8.5, 7.2, 6.4])
    T_error_plus = numpy.array([1, 0.3, 1, 2, 1.5, 1.6, 1.1, 0.7])
    T_error_min = numpy.array([1, 0.7, 1, 2, 1, 1.3, 1.3, 0.8])
    pyplot.errorbar(r_arcmin_Suzaku*60, T_kev_Suzaku,
        xerr=[r_error_min*60, r_error_plus*60],
        yerr=[T_error_min, T_error_plus],
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
        yerr=[T_error_min, T_error_plus],
        c="grey", fmt="+", elinewidth=4, capsize=0)


# Finally, simulation: merger wedge and everything but merger wedge
# --------------------------------------------------------------------------- #
# data files creates using ds9, region/panda centered on CygA,
# MergerAxis has angles -45 to 45; NotMergerAxis has angles 45 to -45
# Annuli in range(0, 200, 42), centered at (1004, 1024) in cubenumber 78

snapshots = numpy.array([78, 77, 68, 60, 53, 47])
ZeroEOrbitFrac = numpy.arange(0.0, 0.8, 0.1)
MergerVelocity = numpy.array([])  # km/s
for vel, sim in zip(ZeroEOrbitFrac,
        ["20160819T2322", "20160819T2357", "20160820T0032", "20160820T0107",
            "20160820T0142", "20160820T0218"]):#, "20160820T0252", "20160820T0328"]):
    MergerAxis = pandas.read_csv(sim+"_MergerAxis.dat", delimiter=" ",
            header=None, names=["Radius", "Temperature", "TemperatureErr"])

    kB_kev_per_kelvin = 8.6173427909e-08
    pixelscale = 10105./2048

    pyplot.plot(MergerAxis["Radius"]*pixelscale,
                MergerAxis["Temperature"]*kB_kev_per_kelvin,
                lw=2, label="Simulation {0:1.1f}: Merger Axis".format(vel))

vel, sim = 0.4, "20160820T0142"
NotMergerAxis = pandas.read_csv(sim+"_NotMergerAxis.dat", delimiter=" ",
        header=None, names=["Radius", "Temperature", "TemperatureErr"])
pyplot.plot(NotMergerAxis["Radius"]*pixelscale,
            NotMergerAxis["Temperature"]*kB_kev_per_kelvin, ls="dashed",
            c="k", lw=2, label="Simulation {0:1.1f}: Average".format(vel))
#pyplot.xlabel("Distance Along Merger Axis [kpc]")
#pyplot.ylabel("Temperature [K]")
pyplot.xlabel("Radius [kpc]")
pyplot.ylabel("kT [keV]")
#pyplot.xscale("log")
pyplot.xlim(5, 1000)
pyplot.ylim(2, 11)
pyplot.xticks([10, 100, 1000], [r"$10^1$", r"$10^2$", r"$10^3$"])
#pyplot.legend(loc="upper left", fontsize=12)
pyplot.show()
# --------------------------------------------------------------------------- #
