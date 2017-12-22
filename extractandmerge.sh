#!/bin/bash

module load daint-mc
module load CDO
module load NCO

#input arguments
simulation=$1
year=$2


for streamno in 01 03 04 05; do
   
    if [ $streamno -eq 01 ]; then
        vars="W_SO,W_SO_ICE,W_I,W_SNOW"
    elif [ $streamno -eq 03 ]; then
        vars="RAIN_CON,RAIN_GSP,TOT_PREC"
    elif [ $streamno -eq 04 ]; then
        vars="ALHFL_S,ALHFL_BS,ALHFL_PL,T_2M"
    elif [ $streamno -eq 05 ]; then
        vars="RUNOFF_S,RUNOFF_G,SNOW_MELT"
    else
        echo "Stream not valid."
    fi



    outstream=out$streamno
 
    #set directory
    directory=/scratch/snx3000/regenass/cosmo5_validation/$simulation/output/$outstream/$year
    cd $directory 


    # loop outputfiles
   for filename in lffd*.nc; do
     echo "Extracting soil variables from $filename"
     cdo selname,$vars $filename ${filename}_extracted.nc
   done

   echo "Mergetime"
   # split to 4 parts first to avoid opening of too many files.
   cdo mergetime lffd${year}01*_extracted.nc lffd${year}02*_extracted.nc lffd${year}03*_extracted.nc ${outstream}_1.nc
   cdo mergetime lffd${year}04*_extracted.nc lffd${year}05*_extracted.nc lffd${year}06*_extracted.nc ${outstream}_2.nc
   cdo mergetime lffd${year}07*_extracted.nc lffd${year}08*_extracted.nc lffd${year}09*_extracted.nc ${outstream}_3.nc
   cdo mergetime lffd${year}10*_extracted.nc lffd${year}11*_extracted.nc lffd${year}12*_extracted.nc ${outstream}_4.nc
   cdo mergetime ${outstream}_*.nc ${outstream}s.nc
   rm *extracted*
   rm ${outstream}_*

   if [ $streamno -eq 03 ]; then
       # sum up to 6hr values
       cdo timselmean,6,1 ${outstream}s.nc tmp1
       # reintroduce first timestep 
       cdo seltimestep,1 ${outstream}s.nc tmp2 
       cdo mergetime tmp1 tmp2 tmp3
       # cut out last (incomplete) timestep
       cdo seldate,${year}-01-01T00:00:00,${year}-12-31T18:00:00 tmp3 ${outstream}s_summed.nc
   fi
 
   vars=""
   echo "DONE Stream ${outstream}."

done
