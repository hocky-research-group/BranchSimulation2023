#!/bin/bash

tprfile=$1 # GROMACS tpr file
deffnm=$2 # output name (can specify the directory and output prefix)
steps=$3 # number of steps to run

# Need the tprfile, output directory/prefix, number of steps to run #
if [ -z "$tprfile" ];then
echo "Usage: $0 tprfile deffnm steps"
exit
fi

# Load GROMACS module #
module load gromacs/openmpi/intel/2020.4

gmxexe=gmx_mpi
$gmxexe mdrun -s $tprfile -deffnm $deffnm -nsteps $steps -ntomp 1
