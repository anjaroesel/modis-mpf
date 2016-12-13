#!/bin/bash


FILES=`ls  MOD09A1.A2001*.hdf`

PROG="modis_hdf2tif_polster_wgs84.sh"
COMMAND1="sh"



for f in $FILES; do
        echo $f
        FOLDERNAME="`echo ${f} | cut -d'.' -f2,3 | tr '.' '_'`"
        echo $FOLDERNAME
        mkdir $FOLDERNAME
        #$COMMAND $PROG $f
        #mv MOD_Grid*.tif $FOLDERNAME
        echo "now moving files"

        
done



#FILE="MOD09A1.A2008161.h14v01.005.2008172103719.hdf"
#NEWNAME="`echo ${i} | cut -d':' -f4,5 | tr ':' '_' | tr -s '|' '_'`"
