#!/usr/bin/python
import numpy as np
import mdtraj as md
import cgmap as cgmap
import os

def load_trj(topfile,coordfile):
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
    return trj, prefix

# 1st function
selstring="backbone"
def mammalian_arp23_four_site(selstring,arp2_name="'ARP2'",arp3_name="'ARP3'"):
    arp2_string_list = ["%s and (residue 7 to 33 or residue 74 to 150 or residue 352 to 387) ",
                       "%s and (residue 34 to 38 or residue 55 to 73)",
                       "%s and (residue 151 to 185 or residue 277 to 326 or residue 339 to 351)",
                       "%s and (residue 186 to 265)",
                       "%s"]

    arp3_string_list = ["%s and (residue 6 to 32 or residue 78 to 153 or residue 375 to 408)",
                       "%s and (residue 33 to 37 or residue 60 to 77)",
                       "%s and (residue 154 to 196 or residue 295 to 344 or residue 362 to 374)",
                       "%s and (residue 197 to 282)",
                       "%s"
                       ]

    # Note the arp3 selection does not contain the dloop!
    monomer_list = [1]
    
    bead_string_list = [ bead_string%("segment_id %s and %s"% (arp2_name,selstring)) for bead_string in arp2_string_list ]
    selection_string_list = [s for i in monomer_list for s in bead_string_list ]
    label_list = ['C%i'%(j+1) for i in monomer_list for j in range(len(bead_string_list)) ]
    resSeq_list = [len(bead_string_list)*i+j+1 for i in range(len(monomer_list)) for j in range(len(bead_string_list)) ]
    segment_id_list = ["A%i"%i for i in monomer_list for s in bead_string_list ]
    chain_list = [i for i in monomer_list for s in bead_string_list ]
    arp2_results = (selection_string_list, label_list, resSeq_list, segment_id_list, chain_list)

    bead_string_list = [ bead_string%("segment_id %s and %s"%(arp3_name,selstring)) for bead_string in arp3_string_list ]
    selection_string_list = [s for i in monomer_list for s in bead_string_list ]
    label_list = ['C%i'%(j+1) for i in monomer_list for j in range(len(bead_string_list)) ]
    resSeq_list = [len(bead_string_list)*i+j+1 for i in range(len(monomer_list)) for j in range(len(bead_string_list)) ]
    segment_id_list = ["A%i"%i for i in monomer_list for s in bead_string_list ]
    chain_list = [i for i in monomer_list for s in bead_string_list ]
    arp3_results = (selection_string_list, label_list, resSeq_list, segment_id_list, chain_list)

    return arp2_results, arp3_results

def mammalian_arp23_two_site(selstring,arp2_name="'ARP2'",arp3_name="'ARP3'"):
    arp2_string_list = ["%s and (residue 151 to 185 or residue 277 to 326 or residue 339 to 351)",
                       "%s and (residue 186 to 265)"]

    arp3_string_list = ["%s and (residue 154 to 196 or residue 295 to 344 or residue 362 to 374)",
                       "%s and (residue 197 to 282)"]

    # Note the arp3 selection does not contain the dloop!
    monomer_list = [1]

    bead_string_list = [ bead_string%("segment_id %s and %s"% (arp2_name,selstring)) for bead_string in arp2_string_list ]
    selection_string_list = [s for i in monomer_list for s in bead_string_list ]
    label_list = ['C%i'%(j+1) for i in monomer_list for j in range(len(bead_string_list)) ]
    resSeq_list = [len(bead_string_list)*i+j+1 for i in range(len(monomer_list)) for j in range(len(bead_string_list)) ]
    segment_id_list = ["A%i"%i for i in monomer_list for s in bead_string_list ]
    chain_list = [i for i in monomer_list for s in bead_string_list ]
    arp2_results = (selection_string_list, label_list, resSeq_list, segment_id_list, chain_list)

    bead_string_list = [ bead_string%("segment_id %s and %s"% (arp3_name,selstring)) for bead_string in arp3_string_list ]
    selection_string_list = [s for i in monomer_list for s in bead_string_list ]
    label_list = ['C%i'%(j+1) for i in monomer_list for j in range(len(bead_string_list)) ]
    resSeq_list = [len(bead_string_list)*i+j+1 for i in range(len(monomer_list)) for j in range(len(bead_string_list)) ]
    segment_id_list = ["A%i"%i for i in monomer_list for s in bead_string_list ]
    chain_list = [i for i in monomer_list for s in bead_string_list ]
    arp3_results = (selection_string_list, label_list, resSeq_list, segment_id_list, chain_list)

    return arp2_results, arp3_results


def get_angle_dist_arp23(trj):
    selstring="backbone"
    arp23_results_list = mammalian_arp23_four_site(selstring)
    arp23_dist_dihed = []
    com_positions = []
    for subunit_info in arp23_results_list:
        sel_list,label_list, resSeq_list, segment_id_list, chain_list = subunit_info
        #print(subunit_info)
        trj_cg = cgmap.cg_by_selection(trj,sel_list, bead_label_list=label_list, resSeq_list=resSeq_list, segment_id_list=segment_id_list )
        dihedral_map = ( md.compute_dihedrals(trj_cg, ((1,0,2,3),),periodic=True) )
        dihed_trj = np.degrees(dihedral_map.mean(axis=1))
        dist_map = ( md.compute_distances(trj_cg, ((1,3),),periodic=True) )
        dist_trj = dist_map.mean(axis=1)
        arp23_dist_dihed.append((dist_trj,dihed_trj))

    arp2_dist = arp23_dist_dihed[0][0]
    arp2_dihed = arp23_dist_dihed[0][1]
    arp3_dist = arp23_dist_dihed[1][0]
    arp3_dihed = arp23_dist_dihed[1][1]

    return arp2_dist, arp2_dihed, arp3_dist, arp3_dihed

# Using only S3 and S4 for the COM distance
def get_com_dist_arp23(trj):
    selstring="backbone"
    arp23_results_list = mammalian_arp23_two_site(selstring)
    arp23_dist_dihed = []
    com_positions = []
    for subunit_info in arp23_results_list:
        #print(subunit_info)
        sel_list,label_list, resSeq_list, segment_id_list, chain_list = subunit_info
        trj_cg = cgmap.cg_by_selection(trj,sel_list, bead_label_list=label_list, resSeq_list=resSeq_list, segment_id_list=segment_id_list )
        com_positions.append(md.compute_center_of_mass(trj_cg))

    com_array = np.array(com_positions)
    dr = com_array[0]-com_array[1]
    arp2_arp3_dist = np.sqrt( (dr*dr).sum(axis=-1) ) 
    
    return arp2_arp3_dist

# 3rd function
def get_residue_distance(trj,segname, res1,res2,atom1="CA",atom2="CA"):
    atom1 = trj.top.select("segname '%s' and residue %i \
                               and name '%s'"%(segname,res1,atom1))[0]
    atom2 = trj.top.select("segname '%s' and residue %i \
                               and name %s"%(segname,res2,atom2))[0]
    distances = md.compute_distances(trj,((atom1,atom2),),periodic=True)
    return distances[:,0]

# 4th function
def get_residue_dist_angle_dihedral(trj,segname_list, residue_list, atom_name_list=None):
    atom_list = []
    for i in range(len(segname_list)):
        if atom_name_list is None:
            atom_name = "CA"
        else:
            atom_name = atom_name_list[i]
            
        atom_list.append(
            trj.top.select("segname '%s' and residue %i \
                               and name '%s'"%(segname_list[i],residue_list[i],atom_name))[0]
        )

    if len(segname_list) == 4:
        dihedrals = np.degrees(
            md.compute_dihedrals(trj,(atom_list,),periodic=True)
        )
        return dihedrals[:,0]
    elif len(segname_list) == 3:
        angles = np.degrees(
            md.compute_angles(trj,(atom_list,),periodic=True)
        )
        return angles[:,0]
    elif len(segname_list) == 2:
        distances = md.compute_distances(trj,(atom_list,),periodic=True)
        return distances[:,0]
    else:
        return None

def get_residue_distance(trj,segname, res1,res2,atom1='CA',atom2='CA'):
    atom1 = trj.top.select("segname '%s' and residue %i \
                               and name '%s'"%(segname,res1,atom1))[0]
    atom2 = trj.top.select("segname '%s' and residue %i \
                               and name '%s'"%(segname,res2,atom2))[0]
    distances = md.compute_distances(trj,((atom1,atom2),),periodic=True)
    return distances[:,0]

# Other distances
def other_distances(trj):
    d2_116_187 = get_residue_distance(trj, 'ARP2', 116, 187)
    d2_72_211 = get_residue_distance(trj,'ARP2',72,211)
    d3_119_198 = get_residue_distance(trj,'ARP3',119,198)
    d3_76_222 = get_residue_distance(trj,'ARP3',76,222)
    return d2_116_187, d2_72_211, d3_119_198, d3_76_222

# Clamp twisting
def clamp_twist(trj):
    bead_string_list = ["segment_id 'RPC4' and (residue 14 to 50 or residue 63 to 73 or residue 117 to 132) and backbone",
                       "segment_id 'RPC4' and (residue 3 to 13 or residue 51 to 62 or residue 74 to 77 or residue 135 to 142) and backbone",
                       "segment_id 'RPC2' and (residue 125 to 262) and backbone",
                       "segment_id 'ARP3' and (residue 6 to 32 or residue 78 to 153 or residue 375 to 408 or residue 33 to 37 or residue 60 to 77) and backbone"]

    # Note the arp3 selection does not contain the dloop!
    monomer_list = [1]

    selection_string_list = [s for i in monomer_list for s in bead_string_list ]
    label_list = ['C%i'%(j+1) for i in monomer_list for j in range(len(bead_string_list)) ]
    resSeq_list = [len(bead_string_list)*i+j+1 for i in range(len(monomer_list)) for j in range(len(bead_string_list)) ]
    segment_id_list = ["A%i"%i for i in monomer_list for s in bead_string_list ]
    chain_list = [i for i in monomer_list for s in bead_string_list ]
    clamp_results = (selection_string_list, label_list, resSeq_list, segment_id_list, chain_list)

    sel_list,label_list, resSeq_list, segment_id_list, chain_list = clamp_results
    trj_cg = cgmap.cg_by_selection(trj,sel_list, bead_label_list=label_list, resSeq_list=resSeq_list, segment_id_list=segment_id_list )
    dihedral_map = np.degrees(( md.compute_dihedrals(trj_cg, ((0,1,2,3),),periodic=True) ))
    
    return dihedral_map

# Helix bending
def helix_bending(trj):
    helix_bending = get_residue_dist_angle_dihedral(trj,("RPC4","RPC4","RPC4"), (166,143,130))
    return helix_bending

# Movement of C-term groove
def c_term_groove(trj):
    d2_152_387 = get_residue_distance(trj,'ARP2',152,387)
    d3_163_409 = get_residue_distance(trj,'ARP3',163,409)
    return d2_152_387, d3_163_409

# Opening of BEG, W-loop uncurling
def BEG_Wloop(trj):
    # BEG1
    d2_173_140 = get_residue_distance(trj,'ARP2',173,140)
    d2_173_371 = get_residue_distance(trj,'ARP2',173,371)

    # BEG2
    d3_184_143 = get_residue_distance(trj,'ARP3',184,143)
    d3_184_393 = get_residue_distance(trj,'ARP3',184,393)
    return d2_173_140, d2_173_371, d3_184_143, d3_184_393

# Cleft opening/closing
def cleft_opening_closing(trj):
    # Arp2
    d2_15_162 = get_residue_distance(trj,'ARP2',15,162)
    d2_16_161 = get_residue_distance(trj,'ARP2',16,161)
    d2_62_211 = get_residue_distance(trj,'ARP2',62,211)

    # Arp3
    d3_14_173 = get_residue_distance(trj,'ARP3',14,173)
    d3_15_172 = get_residue_distance(trj,'ARP3',15,172)
    d3_67_222 = get_residue_distance(trj,'ARP3',67,222)

    return d2_15_162, d2_16_161, d2_62_211, d3_14_173, d3_15_172, d3_67_222












