#!/bin/bash
# This is an example runscript for NITROGEN/VALENCE.
# You might need to adjust the paths.
#
# This sets the appropriate environment variables and paths
# and then calls NITROGEN, using a VSVB input 
#
# Usage:
#    $ ./run_script.sh path_to_nitrogen_executable nitrogen_input VALENCE_input
#   For example: 
#    $ ./run_script.sh ~/nitrogen_v1.10dev/nitrogen/bin/nitrogen VSVB_STAT.job h2o.inp

export VALENCE_ROOT=/home/user/VALENCE/
# put the directories for libvalence.so and libsimint.so in the LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$VALENCE_ROOT/lib/:$VALENCE_ROOT/simint/lib64/:$LD_LIBRARY_PATH
# this variable needs to be set to the path of a script to run
# VALENCE in parallel during the NITROGEN calculation.
# if you're not running in parallel, it's not important
export VALENCE_SCRIPT=$VALENCE_ROOT/nitrogen/parallel_run_script.sh

if [ "$#" -ne 3 ]; then
    echo "This sets the appropriate environment variables and paths"
    echo "and then calls NITROGEN, using a VSVB input "
    echo "  Usage:"
    echo "   $ ./run_script.sh path_to_nitrogen_executable nitrogen_input VALENCE_input"
    echo "   For example: "
    echo "   $ ./run_script.sh ~/nitrogen_v1.10dev/nitrogen/bin/nitrogen VSVB_STAT.job h2o.inp"
    exit 1
fi

$1 $2 < $3
