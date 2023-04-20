**An example command to run a MD simulation (after calling the GROMACS module) with PLUMED is the following**\
gmx_mpi mdrun -s $tprfile$ -deffnm $outname$ -nsteps $steps$ -ntomp 1 -plumed $plumedfile$

**In this case, the $tprfile$ refers to "mammalian_junction_nodaughter_ATP_c27_ionized_npt_5000ns_every50ps-g2020.4-UNCHECKED_new_restrained.tpr" which can be accessed through the link below**
https://drive.google.com/drive/folders/1pJPzwOkZ_i67u2cEODVmURXGOIU-LB6x

**$plumedfile$ refers to any one of the three input files provided**
1) plumed_steerMD_domaindist_nodaughter_k10000_60ns.dat
2) plumed_steerMD_domaindist_nodaughter_k10000_100ns.dat
3) plumed_steerMD_domaindist_nodaughter_k10000_150ns.dat
