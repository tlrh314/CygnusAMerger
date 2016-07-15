#!/usr/bin/env bash

set -e

FPS=2
TIMESTAMP="20160707T0034"

# for i in {0..360..10}; do
#     leadingzeros="$(printf "%03d" ${i})"
#     mv "../runs/20160707T0034/out/rotation_xray_${i}_0_0.pdf" \
#         "../runs/20160707T0034/out/rotation_xray_${leadingzeros}_000_000.pdf"
#     mv "../runs/20160707T0034/out/rotation_xray_0_${i}_0.pdf" \
#         "../runs/20160707T0034/out/rotation_xray_000_${leadingzeros}_000.pdf"
#     mv "../runs/20160707T0034/out/rotation_xray_0_0_${i}.pdf" \
#         "../runs/20160707T0034/out/rotation_xray_000_000_${leadingzeros}.pdf"
# done

for i in {0..360..10}; do
    i="$(printf "%03d" ${i})"
    nice -n 19 ffmpeg -y -r $FPS -i ../runs/${TIMESTAMP}/out/rotation_xray_${i}_000_000.png \
        -profile:v high444 -level 4.1 -c:v libx264 -preset slow -crf 25 -s '2000:2000' \
        -an "../runs/${TIMESTAMP}/out/rotation_xray_x_000_000.mp4"

    nice -n 19 ffmpeg -y -r $FPS -i ../runs/${TIMESTAMP}/out/rotation_xray_000_${i}_000.png \
        -profile:v high444 -level 4.1 -c:v libx264 -preset slow -crf 25 -s '2000:2000' \
        -an "../runs/${TIMESTAMP}/out/rotation_xray_000_x_000.mp4"

    nice -n 19 ffmpeg -y -r $FPS -i ../runs/${TIMESTAMP}/out/rotation_xray_000_000_${i}.png \
        -profile:v high444 -level 4.1 -c:v libx264 -preset slow -crf 25 -s '2000:2000' \
        -an "../runs/${TIMESTAMP}/out/rotation_xray_000_000_x.mp4"
done
