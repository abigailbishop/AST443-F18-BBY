#! /bin/bash -u

for file in $(ls /astrolab/Fall_18/ebiermann/AST443-F18-BBY-data/Lab1/A2-solved_images/*)
do
    sex ${file} -c default.se -CATALOG_NAME ${file%%.*}.cat
done
