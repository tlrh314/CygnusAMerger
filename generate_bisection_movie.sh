#!/usr/bin/env bash

set -e 
FPS=2
nice -n 19 ffmpeg -y -r $FPS -i out/findmass_cygA_freebeta_900ksec_dark%3d.png \
    -profile:v high444 -level 4.1 -c:v libx264 -preset slow -crf 25 -s '2000:2000' \
    -an "out/findmass_cygA_freebeta_900ksec_dark.mp4"
nice -n 19 ffmpeg -y -r $FPS -i out/findmass_cygA_900ksec_dark%3d.png \
    -profile:v high444 -level 4.1 -c:v libx264 -preset slow -crf 25 -s '2000:2000' \
    -an "out/findmass_cygA_900ksec_dark.mp4"


nice -n 19 ffmpeg -y -r $FPS -i out/findmass_cygB_freebeta_900ksec_dark%3d.png \
    -profile:v high444 -level 4.1 -c:v libx264 -preset slow -crf 25 -s '2000:2000' \
    -an "out/findmass_cygB_freebeta_900ksec_dark.mp4"
nice -n 19 ffmpeg -y -r $FPS -i out/findmass_cygB_900ksec_dark%3d.png \
    -profile:v high444 -level 4.1 -c:v libx264 -preset slow -crf 25 -s '2000:2000' \
    -an "out/findmass_cygB_900ksec_dark.mp4"
