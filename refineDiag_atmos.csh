#!/bin/tcsh
###################################################################
#
##        Refine Atmospheric Diagnostics for CMIP6
#
####################################################################

set source_dir = $1     # CODE_DIRECTORY
set refineDiagDir = $2  # OUTPUT DIRECTORY

set refineAtmosErrors = 0

foreach INFILE (`/bin/ls *.{atmos,aerosol}*cmip*.nc`)
  # replace cmip in output filename by refined
  set OUTFILE = `echo $INFILE | sed -e 's/_cmip/_refined/'`
  $source_dir/refine_Atmos_no_redux.py $INFILE -o $refineDiagDir/$OUTFILE -v
  if ($?) @ refineAtmosErrors++
end

foreach INFILE (`/bin/ls *.atmos_daily*cmip*.nc`)
  # replace cmip in output filename by refined
  set OUTFILE = `echo $INFILE | sed -e 's/_cmip/_refined/' -e 's/daily/month/'`
  $source_dir/refine_Atmos_redux_daily_to_monthly.py $INFILE -o $refineDiagDir/$OUTFILE -v
  if ($?) @ refineAtmosErrors++
end
