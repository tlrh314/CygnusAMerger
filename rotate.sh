#!/usr/bin/env bash

set -e

FPS=3
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
    #touch ../runs/${TIMESTAMP}/out/rotation_xray_000_000_${i}.png  
    #touch ../runs/${TIMESTAMP}/out/rotation_xray_000_${i}_000.png  
    #touch ../runs/${TIMESTAMP}/out/rotation_xray_${i}_000_000.png  
    touch ../runs/${TIMESTAMP}/out/rotation_xray_180_${i}_000.png  
done

convert -delay 10 -loop 0  \
    ../runs/${TIMESTAMP}/out/rotation_xray_180_*_000.png  \
    ../runs/${TIMESTAMP}/out/rotation_xray_180_x_000.gif

exit 0

convert -delay 10 -loop 0  \
    ../runs/${TIMESTAMP}/out/rotation_xray_000_000_*.png  \
    ../runs/${TIMESTAMP}/out/rotation_xray_000_000_x.gif

convert -delay 10 -loop 0  \
    ../runs/${TIMESTAMP}/out/rotation_xray_000_*_000.png  \
    ../runs/${TIMESTAMP}/out/rotation_xray_000_x_000.gif

convert -delay 10 -loop 0  \
    ../runs/${TIMESTAMP}/out/rotation_xray_*_000_000.png  \
    ../runs/${TIMESTAMP}/out/rotation_xray_x_000_000.gif
