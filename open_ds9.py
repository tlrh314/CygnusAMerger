#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import argparse
import logging
try:
    import pyds9
except ImportError:
    print("Program requires installation of ds9")
    exit(1)

try:
    import IPython
except ImportError:
    print("Program requires IPython installation")
    exit(1)
try:
    from astropy.io import fits as pyfits
except ImportError:
    try:
        import pyfits
    except ImportError:
        has_fits = False
    else:
        has_fits = True

else:
    has_fits = True


logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s|%(name)s|%(levelname)s|%(message)s")
logger = logging.getLogger("open_ds9")


def open_file(d, fname):
    try:
        d.set("file {}".format(fname))
    except ValueError as err:
        print str(err)
        print has_fits
        if "XPA$ERROR Unable to load fits" in str(err) and has_fits:
            logger.exception(
                "Cannot load fits file natively, using pyfits")
            d.set_np2arr(pyfits.getdata(fname).astype(float))
        else:
            raise


def main(args):
    # ds9_targets returns None if no ds9 window is opened
    # ds9_running = pyds9.ds9_targets()
    # if ds9_running:
    #     print "ds9 is already running"
    #     for instance in ds9_running:
    #         targetid = instance.replace("DS9:", "")
    #         d = pyds9.DS9(targetid)
    #         print "  instance =", targetid
    #         while True:
    #             close = raw_input("    Close [y/n]: ")
    #             if close.lower() in "yn":
    #                 if close == "y":
    #                     d.set("exit")
    #                 break
    #             else:
    #                 print "    Invalid character"
    # else:
    # List of commands: http://ds9.si.edu/doc/ref/command.html

    if args.chandra:
        rundir = "/media/SURFlisa/runs/"
        mosaic = "ChandraObservation/StruisMosaics/mosaic/cygnus_tot_flux.fits"
        args.filename = rundir+mosaic


    d = pyds9.DS9()
    time.sleep(0.5)
    open_file(d, args.filename)

    d.set("preserve pan yes")
    d.set("preserve scale yes")
    print d.id
    print d.target

    header = "DS9 instance in object `d`"
    if has_fits:
        f = pyfits.open(args.filename)
        header += "\npyfits instance in object `f`"

    if "temperature" in args.filename:
        d.set("scale linear")
        d.set("scale limits 2e6 2e8")
        d.set("cmap bb")
        time.sleep(0.5)
        a = d.get_arr2np()
        d.set("frame new")
        kB = 8.6173427909e-08  # keV/K
        d.set_np2arr(a*kB)
        d.set("scale limits 0.1 9")
    if "xray" in args.filename:
        d.set("scale log")
        d.set("scale limits 2e-16 0.02")
        d.set("cmap sls")
    if "density" in args.filename:
        d.set("scale log")
        d.set("scale limits 2e-16 0.02")
        d.set("cmap sls")
    if "velocity" in args.filename:
        d.set("scale log")
        d.set("scale limits 1.1e7 2.2e8")
        d.set("cmap bb")

    d.set("scale mode User")
    # d.set("scale scope global")  # takes forever
    # d.set("scale open")

    if "fits.fz" in args.filename:
        # Simulation
        # seems slow until al frames have been loaded at least once
        d.set("cube interval 0.125")
        d.set("cube play")
        #d.set("cube close")
        # d.set("cube 35")

    if "cygnus_tot_flux.fits" in args.filename:
        # Observation
        d.set("scale log")
        d.set("scale limits 7.0e-10 1.0e-6")
        d.set("cmap sls")
        d.set("smooth radius 9")
        d.set("smooth yes")
        d.set("zoom to fit")
        # contrast and bias
        d.set("cmap value 3.08889 0.30112")
        # default
        # d.set("cmap value 1 0.5")

    # reg = "/Users/timohalbesma/Desktop/TemperatureJump/Simulation/"+\
    #    "20160819T2357.reg"
    # d.set("regions load '{0}'".format(reg))

    # This takes forever
    # fits = d.get_arr2np()
    # print fits.shape
    # print fits.dtype

    #pyplot.figure(figsize=(12,9))
    #pyplot.imshow(fits[0], origin="lower")
    #pyplot.show()

    IPython.embed(banner1="", header=header)

    # Close the ds9 window if it is still open
    if pyds9.ds9_targets() and "DS9:{0} {1}".format(d.target, d.id) in pyds9.ds9_targets():
        d.set("exit")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("--zscale", action="store_true")
    parser.add_argument("--chandra", action="store_true")

    main(parser.parse_args())
