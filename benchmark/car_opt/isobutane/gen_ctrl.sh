#!/bin/bash

## Generate control files for cartesian optimization benchmark
## Example:
# name = "ethane_energy"
# path = "benchmark/examples/ethane.mol2"
# calc_type = "car"
# output = "./ethane_energy.log"

for mol2_file in ../examples/*.mol2
do
    mol=$(basename "$mol2_file" .mol2)
    mkdir -p $mol
    cd $mol
    mol2_file="../$mol2_file"
    cp $mol2_file .
    filename="./$(basename "$mol2_file" .mol2).ctrl"
    name="$(basename "$mol2_file" .mol2)"
    path="./$(basename "$mol2_file")"
    calc_type="car"
    output="./$(basename "$mol2_file" .mol2).log"
    
    echo "name = \"$name\"" > $filename
    echo "path = \"$path\"" >> $filename
    echo "calc_type = \"$calc_type\"" >> $filename
    echo "output = \"$output\"" >> $filename
    cd ..
done