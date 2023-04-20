**An example command to run a MD simulation (after calling the GROMACS module) with PLUMED is the following**\
gmx_mpi mdrun -s $tprfile$ -deffnm $outname$ -nsteps $steps$ -ntomp 1 -plumed $plumedfile$

**In this case, the $tprfile$ refers to "mammalian_junction_justcomplex_ADP_c27_ionized_shift_npt_5000ns_every50ps-g2020.4-MODIFIED.tpr" which can be accessed through the link below**
https://drive.google.com/drive/folders/1zkX30tUwU5p8Hfa7pNKjD4esiA7Pu89U

**$plumedfile$ refers to any one of the three input files provided**
1) mammalian_junction_justcomplex_ATP_k10000_60ns_plumed_domain_dist.dat
2) mammalian_junction_justcomplex_ATP_k10000_100ns_plumed_domain_dist.dat
3) mammalian_junction_justcomplex_ATP_k10000_150ns_plumed_domain_dist.dat
