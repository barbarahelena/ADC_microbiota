#!/bin/bash
./XGBeast.py -name amyloid -path amycsf -n 200 -param param_grid_tuned.json -x class
./XGBeast.py -name amyloid -path amycsf -n 200 -param param_grid_tuned.json -x class -permute

./XGBeast.py -name ptau -path ptau -n 200 -param param_grid_tuned.json -x class
./XGBeast.py -name ptau -path ptau -n 200 -param param_grid_tuned.json -x class -permute

./XGBeast.py -name mta_tot -path mta -n 200 -param param_grid_tuned.json -x class
./XGBeast.py -name mta_tot -path mta -n 200 -param param_grid_tuned.json -x class -permute

./XGBeast.py -name gca_tot -path gca -n 200 -param param_grid_tuned.json -x class
./XGBeast.py -name gca_tot -path gca -n 200 -param param_grid_tuned.json -x class -permute

./XGBeast.py -name cmb_tot -path cmb -n 200 -param param_grid_tuned.json -x class
./XGBeast.py -name cmb_tot -path cmb -n 200 -param param_grid_tuned.json -x class -permute

./XGBeast.py -name faz_tot -path faz -n 200 -param param_grid_tuned.json -x class
./XGBeast.py -name faz_tot -path faz -n 200 -param param_grid_tuned.json -x class -permute
