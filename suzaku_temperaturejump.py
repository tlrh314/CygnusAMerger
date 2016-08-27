import matplotlib
from matplotlib import pyplot
matplotlib.use("Qt4Agg", warn=False)
pyplot.rcParams.update({"font.size": 22})
pyplot.rcParams.update({"xtick.major.size": 8})
pyplot.rcParams.update({"xtick.minor.size": 4})
pyplot.rcParams.update({"ytick.major.size": 8})
pyplot.rcParams.update({"ytick.minor.size": 4})
pyplot.rcParams.update({"xtick.major.width": 2})
pyplot.rcParams.update({"xtick.minor.width": 2})
pyplot.rcParams.update({"ytick.major.width": 2})
pyplot.rcParams.update({"ytick.minor.width": 2})
pyplot.rcParams.update({"xtick.major.pad": 8})
pyplot.rcParams.update({"xtick.minor.pad": 8})
pyplot.rcParams.update({"ytick.major.pad": 8})
pyplot.rcParams.update({"ytick.minor.pad": 8})
pyplot.rcParams.update({"legend.loc": "best"})
pyplot.rcParams.update({"figure.autolayout": True})

def plot_sarazin_suzaku():
    # Eyeballed from Sarazin, Finoguenov & Wik (2012) Fig. 4. Black points.
    r_arcmin_Suzaku = [-2, -0.1, 2, 4, 5.9, 8, 10, 12.3]
    r_error_plus = [1, 1, 1, 1, 1, 1, 1, 1.5]
    r_error_min = [1, 1, 1, 1, 1, 1, 1, 1.5]
    T_kev_Suzaku = [5, 5.5, 6.2, 7, 8.2, 8.5, 7.2, 6.4]
    T_error_plus = [1, 0.3, 1, 2, 1.5, 1.6, 1.1, 0.7]
    T_error_min = [1, 0.7, 1, 2, 1, 1.3, 1.3, 0.8]
    pyplot.errorbar(r_arcmin_Suzaku, T_kev_Suzaku,
        xerr=[r_error_min, r_error_plus],
        yerr=[T_error_min, T_error_plus],
        c="k", fmt="+", elinewidth=4, capsize=0)

    # Eyeballed from Sarazin, Finoguenov & Wik (2012) Fig. 4. Grey points.
    r_arcmin_Suzaku = [0.1, 6.1, 11.1]
    r_error_plus = [3.5, 1.8, 3]
    r_error_min = [3.1, 2.8, 3.2]
    T_kev_Suzaku = [5.6, 8.8, 7]
    T_error_plus = [0.3, 0.7, 0.5]
    T_error_min = [0.4, 0.7, 0.4]
    pyplot.errorbar(r_arcmin_Suzaku, T_kev_Suzaku,
        xerr=[r_error_min, r_error_plus],
        yerr=[T_error_min, T_error_plus],
        c="grey", fmt="+", elinewidth=4, capsize=0)

if __name__ == "__main__":
    pyplot.figure(figsize=(12,9))

    plot_sarazin_suzaku()

    pyplot.text(0, 4, "center", fontsize=36, ha="center", va="top")
    pyplot.text(4.7, 4, "merger\nshock", fontsize=36, ha="left", va="bottom")
    pyplot.text(9, 4, "subcluster", fontsize=36, ha="left", va="bottom")

    pyplot.xlim(-3, 13.5)
    pyplot.xlabel("Distance Along Merger Axis (arcmin)")
    pyplot.xticks([0, 5, 10], ["0", "5", "10"])
    pyplot.minorticks_on()
    pyplot.ylim(0, 10.1)
    pyplot.ylabel(r"$T$ (keV)")

    pyplot.savefig("out/suzaku_temperaturejump.png", dpi=300)
    # pyplot.show()
