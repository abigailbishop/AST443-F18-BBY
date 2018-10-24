#! /bin/bash -u

for file in $(ls /astrolab/Fall_18/ebiermann/Lab1_archive/A1-cal_images/unsolved/*)
do
    solve-field --ra 300.1833 --dec 22.7108 --radius 0.5 --continue ${file}
done
