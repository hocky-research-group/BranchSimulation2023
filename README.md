# Molecular Simulations Reveal Structural Changes Supporting A Multi-step Model Of Arp2/3 Activation (In prep.)
### Yuvraj Singh, Brad J. Nolen, Glen M. Hocky
**This repository sample trajectory (.xtc format) and topology (.psf format) files provided to peform preliminary analysis as well as the necessary simulation and analysis scripts. The repository consists of the following directories**

1) **7TPT**: Branch junction pdb
2) **BtArp23_splay**: Inactive Arp2/3 complex pdb
3) **full_junction**: Full branch junction sample trajectory and topology files
4) **just_complex**: Active Arp2/3 complex sample trajectory and topology files
5) **morph**: Snapshots of 31 cryo-em structures of Arp2/3 complex transitioning from the splayed configuration to shortpitch configuration
6) **no_daughter**: Active Arp2/3 complex (bound to mother filament only) sample trajectory and topology files
7) **splayed**: Inactive Arp2/3 complex sample trajectory and topology files

## Analysis and Simulation files

**Two Jupyter notebooks (.ipynb) are provided for the following functions**

$\bullet$ BranchSimulation_calculations.ipynb: To run relevant calculations and save the data as text files, that can be read

$\bullet$ BranchSimulation_plotting.ipynb: To visualize the output data from "BranchSimulation_calculations.ipynb"

**Due to the bulk of the simulation input files, not all files have been uploaded to this repository. Any additional files can be provided upon request. The GROMACS input files (.tpr) for the Branch junction, Active Arp2/3 complex (bound to mother filament only),  Active Arp2/3 complex, and Inactive Arp2/3 complex systems can be accessed through the following link**\
https://drive.google.com/drive/folders/10LAyCBMfzWyq-rSxwd-FKUzbqMsSdXj6 

**An example command to run a MD simulation (after calling the GROMACS module) is the following**\
gmx_mpi mdrun -s $tprfile$ -deffnm $outname$ -nsteps $steps$ -ntomp 1

**The above command can also be executed using the example bash script (run_gromacs_test.sh) provided**\
bash run_gromacs_test.sh $tprfile$ $outname$ $steps$

In addition, we have provided the following scripts to facilitate with this repository

$\bullet$ cgmap.py: Python script to map all-atom groups to coarse-grained quatities (https://software.rcc.uchicago.edu/mscg/docs/intro.html)

$\bullet$ mammalian_functions.py: All pertinent functions used in "BranchSimulation_calculations.ipynb"

$\bullet$ unzip_files.sh: To unzip any zip files

## Package rquirements
The following python packages should be installed before running any of the jupyter notebooks

$\bullet$ MDTraj https://mdtraj.org/1.9.4/installation.html 

$\bullet$ Seaborn https://seaborn.pydata.org/installing.html 
