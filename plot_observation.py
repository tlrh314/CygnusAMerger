import matplotlib
matplotlib.use("Qt4Agg")
from matplotlib import pyplot
pyplot.rcParams.update({"font.size": 22})
pyplot.rcParams.update({"xtick.major.size": 16})
pyplot.rcParams.update({"xtick.minor.size": 8})
pyplot.rcParams.update({"ytick.major.size": 16})
pyplot.rcParams.update({"ytick.minor.size": 8})

from cluster import ObservedCluster


def plot_observation(cygA, cygB):
    poster_style = True
    if poster_style:
        pyplot.style.use(["dark_background"])
        # magenta, dark blue, orange, green, light blue (?)
        data_colour = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                       (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                       (41./255, 239./255, 239./255)]
        fit_colour = "white"
    else:
        data_colour = ["g", "r", "b"]
        fit_colour = "k"

    fig = pyplot.figure(figsize=(12, 9))

    arcsec2kpc = cygA.cc.kpc_DA
    pyplot.errorbar(cygA.radius/arcsec2kpc+cygA.binsize/2, cygA.number_density,
                    xerr=cygA.binsize/2, yerr=cygA.number_density_std,
                    marker="o", ms=5 if poster_style else 3,
                    elinewidth=3 if poster_style else 1,
                    ls="", c=data_colour[0], label="CygA")
    arcsec2kpc = cygB.cc.kpc_DA
    pyplot.errorbar(cygB.radius/arcsec2kpc+cygB.binsize/2, cygB.number_density,
                    xerr=cygB.binsize/2, yerr=cygB.number_density_std,
                    marker="o", ms=5 if poster_style else 3,
                    elinewidth=3 if poster_style else 1,
                    ls="", c=data_colour[2], label="CygB")

    pyplot.gca().set_xscale("log")
    pyplot.gca().set_yscale("log")
    # pyplot.tight_layout()
    pyplot.xlim(0.5, 2000)
    pyplot.ylabel(r"$n_e(r)$ [1/cm$^{-3}$]")
    pyplot.xlabel(r"$r$ [arcsec]")

    pyplot.legend()
    pyplot.savefig("out/raw_observed_clusters{0}.png"\
        .format("_dark" if poster_style else ""), dpi=300)
    pyplot.show()


if __name__ == "__main__":
    cygA_observed_900ksec = ObservedCluster("cygA")
    cygB_observed_900ksec = ObservedCluster("cygB")

    plot_observation(cygA_observed_900ksec, cygB_observed_900ksec)
