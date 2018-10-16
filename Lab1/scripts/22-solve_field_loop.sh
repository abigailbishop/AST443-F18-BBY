#! /bin/bash -u

for file in $(ls /astrolab/Fall_18/ebiermann/AST443-F18-BBY-data/Lab1/A1-cal_images/*)
do
    solve-field --ra 348.5 --dec 8.70 --radius 0.5 --continue ${file}
done
