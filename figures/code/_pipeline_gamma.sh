#!/bin/sh
# Runs the 193M b010 plate extraction, then the figure scripts that need the gamma and b010 caches.
# usage: nohup sh _pipeline_gamma.sh > /mnt/Elements/DNS/analysis/asme2027_figs/gamma/pipeline.log 2>&1 &
cd /home/sk999/Documents/ASME_2027/figures/code
L=/mnt/Elements/DNS/analysis/asme2027_figs/gamma
python3 _bl_extract.py c193b > $L/bl_c193b_b010.log 2>&1 &
BL=$!
until [ -f $L/c193b_z010.npz ] && [ -f $L/c193b_x10.npz ] && [ -f $L/c20bw_z010.npz ]; do sleep 30; done
python3 fig_spots.py > $L/fig_spots.log 2>&1
wait $BL
for f in fig_wall_state fig_bl_params fig_bl_profiles fig_eta_x; do python3 $f.py > $L/$f.log 2>&1; done
echo pipeline done > $L/pipeline.done
