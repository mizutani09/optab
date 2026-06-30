#!/bin/bash

set -euo pipefail

# # check around 1000 K
# label="withline_11_11_lowtemp"
# tmin=2.9d0
# tmax=3.1d0
# lmax=11
# pmin=-5.5d0
# pmax=14.5d
# jmax=11

# # goal
# label="withline_11_11_test"
# tmin=1d0
# tmax=6.0d0
# lmax=11
# pmin=-8.0d0
# pmax=14.5d0
# jmax=11

# label="test_default"
# tmin=2.5d0
# tmax=7.5d0
# lmax=21
# pmin=-4d0
# pmax=9d0
# jmax=21

# # test around 100 K
label="withline_11_11_verylowtemp_lnK"
tmin=1.85d0
tmax=1.90d0
lmax=11
pmin=-8.0d0
pmax=14.5d0
jmax=11

cd $OPTAB/work/FastChem-lnk_interpolate_dev/input

# 'table_detail3'	label
# 2.5d0, 6.0d0, 11	tmin, tmax, lmax
# -8.0d0, 14.5d0, 11		pmin, pmax, jmax

sed -i "1s|.*|'${label}'\tlabel|" prep_FastChem.dat
sed -i "2s|.*|${tmin}\, ${tmax}\, ${lmax}\ttmin, tmax, lmax|" prep_FastChem.dat
sed -i "3s|.*|${pmin}\, ${pmax}\, ${jmax}\t\tpmin, pmax, jmax|" prep_FastChem.dat
./prep_FastChem

cd $FASTCHEM
./fastchem input/config.input_${label}

exit 0

# Convert the output from FastChem to the format required by optab.
cd $OPTAB/eos/FastChem
$OPTAB/eos/src/convert_FastChem $FASTCHEM/output/${label}.dat

# # For checking the output
# python ../python/eos.py ${label}.h5 mmw --syms=100

# Rebuild executable to avoid running stale src/a.out after merges.
make -C $OPTAB/src -j"$(nproc)"

cd $OPTAB/sample/
sed -i "5s|.*|export EOS='/home/kosuke/simulation/optab/eos/FastChem/${label}.h5'|" sample.sh
bash sample.sh
# if there is ${label} directory, remove it before moving the sample output to ${label}
if [ -d "${label}" ]; then
    rm -rf "${label}"
fi

mv sample ${label}

# For checking the output
# python python/opac.py table ross 150


cd $OPTAB/sample/python
python make_opacity_table_hdf5.py ../${label} --out opacity_${label}_athena.h5 --nlogt 128 --nlogrho 128 --crop full
# python plot_opacity_table_hdf5.py --mean ross --prefix test opacity_${label}_athena.h5
# python plot_opacity_table_hdf5.py --mean pla --prefix test opacity_${label}_athena.h5
python plot_opacity_table_hdf5_4panel.py opacity_${label}_athena.h5 --out optab_4panel_${label}.png
