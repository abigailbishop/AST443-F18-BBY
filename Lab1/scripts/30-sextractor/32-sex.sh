#! /bin/bash -u

for file in $(ls /astrolab/Fall_18/ebiermann/AST443-F18-BBY-data/Lab1/A1-cal_images)
do
    sex ${file} -c default.se -31-Lab1.cat ../cat/${file%%.*}.cat
done
