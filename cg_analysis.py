import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mdtraj as md
import cgmap as cgmap
import seaborn as sns

if len(sys.argv) != 3:
    print("Usage:", sys.argv[0])
    print("       topfile")
    print("       coordfile")
    exit(0)

topfile=sys.argv[1]
coordfile=sys.argv[2]

if coordfile=="":
    print("Working with Morph trajectory")
    if topfile.startswith("pombe"):
        print("Working with Pombe pdb")
        from pombe_functions import *
    else:
        print("Working with Mammalian morph")
        from mammalian_functions import *
elif coordfile.startswith("mammalian"):
    print("Working with Mammalian Arp2/3")
    from mammalian_functions import *
elif coordfile.startswith("pombe"):
    print("Working with Pombe Arp2/3")
    from pombe_functions import *
else:
    exit(0)

# Load the trajectory
if coordfile=="":
    trj = md.load(topfile)
    for a in trj.top.atoms: a.mass = a.element.mass 
    for a in trj.top.atoms: a.charge = 0
    prefix=os.path.splitext(topfile)[0]
else:
    trj = md.load(coordfile, top=topfile)
    for a in trj.top.atoms: a.mass = a.element.mass # computing the masses of each atom
    for a in trj.top.atoms: a.charge = 0 # setting all charges to neutral
    prefix=os.path.splitext(coordfile)[0]

print(prefix)

#print(get_angle_dist_arp23(trj))
np.savetxt("%s_cg_quantities.txt" % prefix,np.array(get_angle_dist_arp23(trj)).T,fmt="%.10f")

