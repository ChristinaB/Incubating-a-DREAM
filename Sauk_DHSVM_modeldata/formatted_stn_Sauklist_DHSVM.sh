#!/bin/sh

# use awk to format input for DHSVM
 ## run command example "./formatted_stn_Sauklist_DHSVM.sh Sauk_94CentroidVICpoints_UTM.csv Sauk_94CentroidVICpoints_inputlist.txt"

# gmauger,cband Jan-Feb 2016
# cband, cbev Nov 2016

incsv=$1
outfl=$2

# clear output file, if exists already:
\rm -f $outfl

# skip header line, format output:
#	*** NOTE: assuming an 8-col csv file, with column order:
#		FID,OBJECTID_1,OBJECTID,LAT,LONG_,ClimateDEM,NORTHING,EASTING

awk -F, '(NR>1){stnnum=$1+1; \
		printf("Station Name %d = data_%.5f_%.5f\n",stnnum,$3,$4); \
		printf("North Coordinate %d = %.6f\n",stnnum,$6); \
		printf("East Coordinate %d = %.6f\n",stnnum,$7); \
		printf("Elevation %d = %.0f\n",stnnum,$5); \
		printf("Station File %d = ",stnnum); \
		printf("/civil/shared/ecohydrology/SkagitSauk/DHSVM-Glacier/DHSVM/inputLivneh2013_WRFbc_historic/bc_2_WRF_delta_1500to3000/data_%.5f_%.5f\n",$3,$4); \
		printf("\n");}' $incsv > $outfl

