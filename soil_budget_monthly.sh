#!/bin/bash

module load daint-mc
module load CDO
module load NCO
module load NCL

SCRIPTDIR=/users/regenass/Soil_Scripts

if [ $# -lt 3 ]; then
    echo "Usage: sh soil_budget_monthly.sh <simulation> <startyear> <endyear>"
fi 


simulation=$1
startyear=$2
endyear=$3


cd /scratch/snx3000/regenass/cosmo5_validation/${simulation}/postprocessed/monmean

for year in $(seq $startyear $endyear); do
 
    cd $year


    for stream in 03 04 05; do
  
        cdo mergetime ${simulation}_monmean_out${stream}_${year}* monmean_out${stream}_${year}.nc 
   
    done

    stream=01
    cd /scratch/snx3000/regenass/cosmo5_validation/$simulation/output/out01/$year

    cdo mergetime lffd${year}010100.nc lffd${year}020100.nc lffd${year}030100.nc lffd${year}040100.nc\
 lffd${year}050100.nc lffd${year}060100.nc lffd${year}070100.nc lffd${year}080100.nc\
 lffd${year}090100.nc lffd${year}100100.nc lffd${year}110100.nc lffd${year}120100.nc\
 /scratch/snx3000/regenass/cosmo5_validation/${simulation}/postprocessed/monmean/${year}/monmean_out01_${year}.nc

    cd /scratch/snx3000/regenass/cosmo5_validation/${simulation}/postprocessed/monmean
done


for stream in 01 03 04 05; do
    cdo mergetime 19*/monmean_out${stream}_*.nc monmean_out${stream}.nc
done

cd $SCRIPTDIR

sed s/SIMULATION/\"${simulation}\"/g  plot_soil_model_50km_cclm_TEMPLATE.ncl > plot_soil_model_50km_cclm.ncl

ncl plot_soil_model_50km_cclm.ncl
