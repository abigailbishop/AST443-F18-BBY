#! /bin/bash -u

for file in $(ls)
do
    solve-field --ra 348.5 --dec 8.70 --radius 0.5 ${file}
done
