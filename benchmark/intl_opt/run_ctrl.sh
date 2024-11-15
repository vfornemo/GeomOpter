#!/bin/bash

for mol_dir in $(find . -mindepth 1 -type d)
do
cd $mol_dir
for ctrl in ./*.ctrl
do
../../../target/release/GeomOpter $ctrl
done
cd ..
done
