import matplotlib
from matplotlib import pyplot
matplotlib.use("Qt4Agg", warn=False)
#from matplotlib import rc
#matplotlib.rc("font", **{"family":"serif", "serif":["Computer Modern Roman"],
#                         "size":28, "weight":"bold"})


class PlotSettings(object):
    def __init__(self, poster_style=False):
        # magenta, dark blue, orange, green, light blue
        self.c = [(255./255, 64./255, 255./255), (0./255, 1./255, 178./255),
                  (255./255, 59./255, 29./255), (45./255, 131./255, 18./255),
                  (41./255, 239./255, 239./255)]

        # use the same colour for CygA and CygB observation
        self.cygA= "green"
        self.cygB = "blue"

        if not poster_style:
            self.fit = "black"
            self.params = {
                "figure.figsize": (12,9),
                "font.size": 28,
                "xtick.major.size": 8,
                "xtick.minor.size": 4,
                "ytick.major.size": 8,
                "ytick.minor.size": 4,
                "xtick.major.width": 2,
                "xtick.minor.width": 2,
                "ytick.major.width": 2,
                "ytick.minor.width": 2,
                "xtick.major.pad": 8,
                "xtick.minor.pad": 8,
                "ytick.major.pad": 8,
                "ytick.minor.pad": 8,
                "lines.linewidth": 1,
                "lines.markersize": 2,
                "axes.linewidth": 1,
                "legend.loc": "best",
                # "axes.labelsize": 12,
                # "text.fontsize": 12,
                # "legend.fontsize": 12,
                # "xtick.labelsize": 10,
                # "ytick.labelsize": 10,
                # "text.usetex": True,
                # "text.latex.preamble": r"\boldmath",
                # "axes.unicode_minus": True,
            }
        else:
            pyplot.style.use(["dark_background"])
            self.fit = "white"
            self.params = {
                "figure.figsize": (12,9),
                "font.size": 28,
                "font.weight": "bold",
                "xtick.major.size": 8,
                "xtick.minor.size": 4,
                "ytick.major.size": 8,
                "ytick.minor.size": 4,
                "xtick.major.width": 2,
                "xtick.minor.width": 2,
                "ytick.major.width": 2,
                "ytick.minor.width": 2,
                "xtick.major.pad": 8,
                "xtick.minor.pad": 8,
                "ytick.major.pad": 8,
                "ytick.minor.pad": 8,
                "lines.linewidth": 5,
                "lines.markersize": 10,
                "axes.linewidth": 5,
                "legend.loc": "best",
                # "axes.labelsize": 12,
                # "text.fontsize": 12,
                # "legend.fontsize": 12,
                # "xtick.labelsize": 10,
                # "ytick.labelsize": 10,
                # "text.usetex": True,
                # "text.latex.preamble": r"\boldmath",
                # "axes.unicode_minus": True,
            }
        matplotlib.rcParams.update(self.params)


if __name__ == "__main__":
    style = PlotSettings()
    print style.params
    import numpy
    x = numpy.arange(-10, 10, 0.1)
    pyplot.figure()
    pyplot.plot(x, x**2, label="x**2")
    pyplot.legend()
    pyplot.show()
