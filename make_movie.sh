cd out
# -y                    overwrite if outputfile exists
# -r 5                  framerate
# -i snapshot_%3d.png   name of input files /w 3 digits (e.g. 000, 001, 002)
ffmpeg -y -r 5 -i snapshot_%3d.png movie.mp4
