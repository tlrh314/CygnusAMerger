import numpy
from matplotlib import pyplot
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const

from cluster import ObservedCluster
from cosmology import CosmologyCalculator
cc = CosmologyCalculator(z=0.0562)
arcsec2kpc = cc.kpc_DA
kB_kev_per_kelvin = 8.6173427909e-08

from plotsettings import PlotSettings
style = PlotSettings()

def t_cool(n_p, T_g):
    t_cool = 8.5e10*u.yr * (1e-3/u.cm**3)/n_p * numpy.sqrt(T_g/(1e8*u.K))
    return t_cool


def plot_cygA_cooling_curve(sbresults, neresults):

    average_radii = 0.5*(sbresults["Radius1"] +
                         sbresults["Radius2"] ) * arcsec2kpc
    average_binsize = (sbresults["Radius2"] -
                       sbresults["Radius1"]) * arcsec2kpc
    # ["Bin", "V", "kT", "fkT", "n", "fn", "P", "fP", "Yparm"]
    average_kT = neresults["kT"]
    average_fkT = neresults["fkT"]
    # Fix bad bins
    average_kT[0:2] = numpy.nan
    average_fkT[0:2] = numpy.nan
    average_kT[-3:] = numpy.nan
    average_fkT[-3:] = numpy.nan

    # Fix error values where fit did not converge
    # fail_indices = average_kT == average_kT_loerr
    # average_kT_loerr[fail_indices] = 0.05*average_kT[fail_indices]
    # average_kT_hierr[fail_indices] = 0.05*average_kT[fail_indices]

    average_n = neresults["n"]
    average_fn = neresults["fn"]
    average_n[0:2] = numpy.nan
    average_fn[0:2] = numpy.nan
    average_n[-3:] = numpy.nan
    average_fn[-3:] = numpy.nan

    Tcool = t_cool(average_n/u.cm**3, average_kT*u.keV/const.k_B.to(u.keV/u.K))
    fTcool = t_cool(average_fn/u.cm**3, average_fkT*u.keV/const.k_B.to(u.keV/u.K))

    pyplot.figure(figsize=(12, 9))
    # TODO: looks weird
    # pyplot.errorbar(average_radii, Tcool.value, yerr=fTcool.value,
    #                marker="o", ls="", c="b", ms=4)
    pyplot.plot(average_radii, Tcool.value)
    pyplot.xlabel("Radius [arcsec?]")
    pyplot.ylabel("Cooling Time [yr]")
    pyplot.xscale("log")
    pyplot.yscale("log")
    pyplot.tight_layout()

    pyplot.figure(figsize=(12, 9))
    pyplot.errorbar(average_radii, average_n, yerr=average_fn,
                    marker="o", ls="", c="b", ms=4)
    pyplot.xlabel("Radius [arcsec?]")
    pyplot.ylabel("Number Density [1/cm^3]")
    pyplot.xscale("log")
    pyplot.yscale("log")
    pyplot.tight_layout()

    pyplot.figure(figsize=(12, 9))
    pyplot.errorbar(average_radii, (average_kT/const.k_B.to(u.keV/u.K)).value,
                    yerr=(average_fkT/const.k_B.to(u.keV/u.K)).value,
                    marker="o", ls="", c="b", ms=4)
    pyplot.xlabel("Radius [arcsec?]")
    pyplot.ylabel("Temperature [K]")
    pyplot.xscale("log")
    pyplot.yscale("log")
    pyplot.tight_layout()
    pyplot.show()


if __name__ == "__main__":
    datadir = "data/20161108/"

    # AVERAGE / 'quiescent' profile

    # /scratch/martyndv/cygnus/combined/spectral/maps/radial/sn100/cygA_maps
    # Last edit: Oct 18 10:57
    # 252 bins (CygA). kT1, Z1, norm1, nH_Gal, red_chisq, counts, num_bins, model
    cygA_fitresults = ascii.read(datadir+"cygA_1T_fixnH_fitresults.dat")
    # /scratch/martyndv/cygnus/combined/spectral/maps/radial/sn100/cygB_maps
    # Last edit: Oct 18 12:02
    # 36 bins (CygNW)
    cygB_fitresults = ascii.read(datadir+"cygB_1T_fixnH_fitresults.dat")

    # /scratch/martyndv/cygnus/combined/spectral/maps/radial/sn100/cygA_plots
    # Last edit: Oct 18 09:27 (CygA), and Oct 18 11:37 (CygNW).
    # Edit by TLRH after copy:
        # header of datafile: i) removed spaces, ii) renamed Error to avoid double
    # 252 bins (CygA). Radius1, Radius2, SB, SBError, BGRD, BGRDError, AREA
    cygA_sbresults = ascii.read(datadir+"cygA_sb_sn100.dat")
    # 36 bins (CygNW)
    cygB_sbresults = ascii.read(datadir+"cygB_sb_sn100.dat")

    # /scratch/martyndv/cygnus/combined/spectral/maps/radial/pressure_sn100
    # Last edit: Nov  2 14:16 (CygA), and Nov  2 14:21 (CygNW).
    # 252 bins (CygA). Volume, Temperature, number density, Pressure, Compton-Y
    # Edit by TLRH after copy: removed | af beginning and end of each line
    # Header has units and spaces, so rename and override
    header = ["Bin", "V", "kT", "fkT", "n", "fn", "P", "fP", "Yparm"]
    cygA_neresults = ascii.read(datadir+"cygA_sn100_therm_profile.dat",
                                names=header, data_start=1)
    cygB_neresults = ascii.read(datadir+"cygB_sn100_therm_profile.dat",
                                names=header, data_start=1)

    plot_cygA_cooling_curve(cygA_sbresults, cygA_neresults)

    # TODO: find EXPTOT for the above profiles. Is this 1.03Msec?
    # TODO: check differences between 982 ksec, and latest data copied at 20161108
    # TODO: fit betamodel etc
    # TODO: what are the differences between cygA_1T_fixnH_fitresults and the
    #       other file cygA_sn100_therm_profile ? I
    # TODO: what is cygA_sn100_radial_therm_profile.dat? It seems older than
    #       cygA_sn100_therm_profile.dat
    # TODO: why is there no cygB_sn100_radial_therm_profile.dat

    # MERGER / HOT / COLD profiles
    # /scratch/martyndv/cygnus/combined/spectral/maps/sector/plots
    # Last edit: Oct 18 12:33
    # TODO: parse cygnus_sector_fitresults.dat

    # /scratch/martyndv/cygnus/combined/spectral/maps/sector/pressure/
    # Last edit: Oct 18 12:26
    # TODO: parse cygnus_sector_sn100_sbprofile.dat

    # /scratch/martyndv/cygnus/combined/spectral/maps/sector/pressure
    # Last edit:  Nov  2 14:35
    # TODO: parse cygnus_sector_therm_profile.dat

    # TODO: Difference between cygnus_sector_therm_profile.dat and
    #       cygnus_sector_fitresults?
