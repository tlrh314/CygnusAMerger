cd out
# -y                    overwrite if outputfile exists
# -r 5                  framerate
# -i snapshot_%3d.png   name of input files /w 3 digits (e.g. 000, 001, 002)

ffmpeg -y -r 5 -i snapshot_%3d.png movie.mp4

# Should work for Telegram.. but only works on MBP, not iOS
# ffmpeg -y -r 5 -i snapshot_%3d.png -acodec libfaac -b:a 128k -vcodec mpeg4 -b:v 1200k -flags +aic+mv4 movie.mp4 
# Other attempt... does screw up quality
# ffmpeg -y -r 5 -i snapshot_%3d.png -acodec aac -vcodec mpeg4 -s 2560x1280 -strict experimental movie.mp4 
