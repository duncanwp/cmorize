#!/bin/bash
# the first argument should be the experiment name
seasons=( winter summer )
for season in "${seasons[@]}"
do
    members=( 0 1 2 3 4 5 )
    for m in "${members[@]}"
    do
        exp=abrupt_${season}_kick_co2_${m}
        rundir=/gws/nopw/j04/eo_shared_data_vol2/scratch/dwatsonparris/${exp}/
        mkdir ${rundir}/Raw
        mv ${rundir}/* ${rundir}/Raw/
        mkdir ${rundir}/Data
        /home/users/dwatsonparris/miniconda3/envs/cis_env3/bin/python cmorize_ham.py ${rundir}/Raw/${exp}'*_.nc' ${rundir}/Data/ECHAM6.3-HAM2.2_abrupt_zero_co2_${season}_${m} --pdrmip --pdrmip_format --monthly --time '20500101-20600101'
    done
done
