import os
import numpy
from matplotlib import pyplot
import matplotlib.gridspec as gridspec
# pyplot.rcParams.update({'font.size': 22})
# pyplot.rcParams.update({"text.usetex": True})
from astropy.io import fits

def make_video(analysis_dir):
    """ Make a video of four physical properties (options) """
    def plot_panels(chosen):
        """ Fill four panels with one snapshot, save png """
        projection = "_projection-z"

        headers = []
        data = []
        for option in chosen:
            with fits.open(analysis_dir+option+projection+".fits.fz") as f:
                headers.append(f[0].header)
                data.append(f[0].data)

        number_of_snapshots = headers[0]['NAXIS3']
        # number_of_pixels_x = headers[0]['NAXIS1']
        # number_of_pixels_y = headers[0]['NAXIS2']
        # See entire header, including comments starting with "/"
        # for line in repr(headers[0]).split("\n"):
        #    print line
        for n in range(number_of_snapshots):
            # Set up four-panel plot, stitched together
            pyplot.figure(figsize=(16,16))
            gs1 = gridspec.GridSpec(2, 2)
            gs1.update(wspace=0, hspace=0) # remove spacing between axes.

            for i in range(4):  # gridspec indexes start at 0
                ax = pyplot.subplot(gs1[i])
                pyplot.axis('on')
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.set_aspect('equal')

                # Plot every panel
                ax.imshow(data[i][n])

            pyplot.tight_layout()
            pyplot.savefig("out/snapshot_{0:03d}.png".format(n))
            pyplot.close()

    options = {0: "dm-density",
               1: "dm-annihilation",
               2: "physical-density",
               3: "pressure",
               4: "sz-compton-y",
               5: "sz-thermal-dt-over-t",
               6: "temperature-emission-weighted",
               7: "temperature-spectroscopic",
               8: "xray-surface-brightness",
               }

    chosen = (options[0], options[7], options[8], options[4])
    number_of_snapshots = plot_panels(chosen)
    os.system("./make_movie.sh")


if __name__ == "__main__":
    make_video(analysis_dir="../runs/20160526T1354/analysis/")
