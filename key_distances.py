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

# Twisting/Flattening
print("=========Twisting/Flattening=========")
print(other_distances(trj)[0])
np.savetxt("%s_Twisting_Flattening.txt" % prefix,np.array(other_distances(trj)).T,fmt="%.10f")
print("")

## Clamp twisting (OLD)
#print("=========Clamp twisting=========")
#print(clamp_twisting_OLD(trj))
#np.savetxt("%s_Clamp_twisting_OLD.txt" % prefix,np.array(clamp_twisting_OLD(trj)).T,fmt="%.10f")
#print("")

# Clamp twisting
if coordfile.startswith("mammalian"):
    print(clamp_twisting_OLD(trj))
    np.savetxt("%s_Clamp_twisting_OLD.txt" % prefix,np.array(clamp_twisting_OLD(trj)).T,fmt="%.10f")
    print(clamp_twisting_NEW(trj))
    np.savetxt("%s_Clamp_twisting_NEW.txt" % prefix,np.array(clamp_twisting_NEW(trj)).T,fmt="%.10f")
else:
    print(clamp_twisting_NEW(trj))
    np.savetxt("%s_Clamp_twisting.txt" % prefix,np.array(clamp_twisting_NEW(trj)).T,fmt="%.10f")

## Clamp twisting (NEW)
#print("=========Clamp twisting=========")
#print(clamp_twisting_NEW(trj))
#np.savetxt("%s_Clamp_twisting_NEW.txt" % prefix,np.array(clamp_twisting_NEW(trj)).T,fmt="%.10f")
#print("")
print("")
# Helix bending
print("=========Helix bending=========")
print(helix_bending(trj))
np.savetxt("%s_Helix_bending.txt" % prefix,np.array(helix_bending(trj)).T,fmt="%.10f")
print("")

# Movement of C-term groove
print("=========C-term groove=========")
print(c_term_groove(trj)[0])
np.savetxt("%s_C-term_groove.txt" % prefix,np.array(c_term_groove(trj)).T,fmt="%.10f")
print("")

# Opening of BEG, W-loop uncurling
print("=========BEG, W-loop uncurling=========")
print(BEG_Wloop(trj)[0])
np.savetxt("%s_BEG_Wloop.txt" % prefix,np.array(BEG_Wloop(trj)).T,fmt="%.10f")
print("")

# Cleft opening/closing
print(cleft_opening_closing(trj)[0])
np.savetxt("%s_Cleft_opening_closing.txt" % prefix,np.array(cleft_opening_closing(trj)).T,fmt="%.10f")







