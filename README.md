# Key conformation transitions during activation of Arp2/3 complex revealed by molecular dynamics simulations (In prep.)
### Yuvraj Singh, Brad J. Nolen, Glen M. Hocky
**This repository includes a Jupyter Notebook (Arp23_analysis.ipynb) to run relevant analysis on sample trajectories (.xtc formate) and topology (.psf format) files provided. The repository consists of the following directories**

1) **7TPT**; Branch junction pdb
2) **BtArp23_splay**; Inactive Arp2/3 complex pdb
3) **full_junction**; Full branch junction sample trajectory files and input files
4) **just_complex**; Active Arp2/3 complex sample trajectory files and input files
5) **morph**; Snapshots of 31 cryo-em structures of Arp2/3 complex transitioning from the splayed configuration to shortpitch configuration
6) **no_daughter**; Active Arp2/3 complex (bound to mother filament only) sample trajectory files and input files
7) **splayed**; Inactive Arp2/3 complex sample trajectory files and input files

**Further instructions for to run analysis are provided in Arp23_analysis.ipynb.**

**Due to the bulk of the simulation input files, not all files have been uploaded in the repository. Any additional files can always be provided upon requets. The GROMACS input files (.tpr) for the full junction, no daughter, complex alone, and splayed systems can be accessed through the following link**

https://drive.google.com/drive/folders/10LAyCBMfzWyq-rSxwd-FKUzbqMsSdXj6 

**An example command to run a MD simulation (after calling the GROMACS module) if the following**

mpirun -np 1 gmx_mpi mdrun -s $tprfile$ -deffnm $outprefix$ -nsteps $steps$ -ntomp 1

**The above command can also be executed using the example bash script provided**

bash run_gromacs_test.sh $tprfile$ $outprefix$ $steps$
