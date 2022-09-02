#!/usr/bin/env tcsh
###################################################################
#
#        Refine Atmospheric Diagnostics for CMIP6
#
###################################################################
#  required files -> output file
#  -----------------------------
#  *.atmos_month_cmip.tile?.nc  -> *.atmos_month_refined.tile?.nc
#  *.atmos_daily_cmip.tile?.nc  -> *.atmos_month_refined.tile?.nc
#  *.atmos_daily_cmip.tile?.nc  -> *.atmos_daily_refined.tile?.nc
#  *.atmos_tracer_cmip.tile?.nc -> *.atmos_tracer_refined.tile?.nc
#
#  require variables (set outside this script)
#  -----------------
#  set refineDiagDir = ???   # output directory
#                            # input directory = `cwd`
###################################################################

#set refineDiagScriptDir = `dirname $0`

set CODE_DIRECTORY = ./refine_scripts
git clone https://github.com/raphaeldussin/xrefineDiags.git $CODE_DIRECTORY

ln -s $CODE_DIRECTORY/*.csh .
ln -s $CODE_DIRECTORY/*.py .

# check that output directory is defined

if ($?refineDiagDir) then

# run script to create refined fields

   module load python

   ./refineDiag_atmos.csh `pwd` $refineDiagDir
   # or in python
   #./refineDiag_atmos.py -w `pwd` -r $refineDiagDir -s . --only_cmip
   set return_status = $?

else
   echo "ERROR: refineDiagDir is not defined - cannot run refineDiag scripts"
   set return_status = 1
endif

# set the status for frepp wrapper script to access

set status = $return_status

