#! /bin/bash -u

for file in $(ls /astrolab/Fall_18/ebiermann/Lab1_archive/A2-solved_images/*)
do
    sex ${file} -c default.se -CATALOG_NAME ${file%%.*}.cat
done
