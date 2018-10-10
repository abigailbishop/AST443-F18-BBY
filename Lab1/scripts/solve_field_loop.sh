#! /bin/bash -u

for file in $(ls)
do
    solve-field ${file}
done
