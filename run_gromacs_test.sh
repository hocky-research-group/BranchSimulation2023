#!/bin/bash

tprfile=$1
deffnm=$2
steps=$3

# Need the tprfile, output directory/prefix, number of steps to run #
if [ -z "$tprfile" ];then
echo "Usage: $0 tprfile deffnm steps"
exit
fi

# Load GROMACS module #
module load gromacs-plumed/openmpi/intel/2022.3

gmxexe=gmx_mpi
mpirun -np 1 $gmxexe mdrun -s $tprfile -deffnm $deffnm -nsteps $steps -ntomp 1
