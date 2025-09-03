from Bio.PDB import *
from scipy.spatial import distance_matrix
from scipy.spatial.transform import Rotation as Rot

import os
import numpy as np
import pandas as pd
import copy
from pymol import cmd

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

okay_h_names = ['OAH','NH1','NH2','OH','CH2',]

BOND_CUT = 2
INTERACT = 5.0

VDW_CUT = 1.75
INTERACT_RING = 4.5

aa_inv = {
    "H": "HIS",
    "K": "LYS",
    "R": "ARG",
    "D": "ASP",
    "E": "GLU",
    "S": "SER",
    "T": "THR",
    "N": "ASN",
    "Q": "GLN",
    "A": "ALA",
    "V": "VAL",
    "L": "LEU",
    "I": "ILE",
    "M": "MET",
    "F": "PHE",
    "Y": "TYR",
    "W": "TRP",
    "P": "PRO",
    "G": "GLY",
    "C": "CYS",
    "X": "MSE",
    "Z": "CYD",
    "CSO": "CSO",
    'SNN': "SNN"
}
aa_inv.keys()
PROT_AA = list(aa_inv.values())
ope = copy.deepcopy(PROT_AA)
#allow secondary conformations (A + B)
for ii in ope:
    PROT_AA.append('A' + ii)
    PROT_AA.append('B' + ii)

def convert_to_pdb(PRED, WT):
    '''
    input:
        PRED: (str) file location of PREDICTED structure
        WT  : (str) file location of WILD TYPE structure
    output:
        temporary files TEMP_WT and TEMP_PRED in PDB format
    returns:
        none
    '''
    #load the files
    cmd.delete('all')

    cmd.load(PRED,'p0')
    cmd.save('./TEMP_PRED.pdb')
    cmd.delete('all')

    cmd.load(WT,'w0')
    cmd.save('./TEMP_WT.pdb')
    cmd.delete('all')

def get_wt_res(wt_file='./TEMP_WT.pdb',pred_file='./TEMP_PRED.pdb'):
    '''
    input:
        wt_file  : (str) file location of WT structure
        pred_file: (str) file location of Pred structure
    returns:
        cint_align: (list) all residues that interact within 10 Ang of the glycan
        residue_diff: (int) offset between wt numbering and pdb numbering
    '''
    parser=PDBParser()
    structure=parser.get_structure("prot", wt_file)
    pr,pc, res, coor = get_ligand_coor(structure) #protein-res, prot-coor, ligand res, ligand coor
    cint_align = find_interactChains(coor, pc,pc,pr,INTERACT = 10.0)
    #print(cint_align)

    structure=parser.get_structure("prot", pred_file)
    pr_,pc_, res_, coor_ = get_ligand_coor(structure)
    new_res, residue_diff = fix_num(pr_,pr)

    return cint_align, residue_diff

def align_pred(res,diff, wt_file='./TEMP_WT.pdb',pred_file='./TEMP_PRED.pdb'):
    '''
    input:
        res (list): residues to be aligned (WT)
        diff (int): residue alignment (for PRED)
        wt_file  : (str) file location of WT structure
        pred_file: (str) file location of Pred structure
    returns:
        cint_align: (list) all residues that interact within 10 Ang of the glycan
        residue_diff: (int) offset between wt numbering and pdb numbering
    '''
    myres = []
    for ii in res:
        if ii[0] not in myres:
            myres.append(ii[0])
    res = myres
    #print(diff)

    pred_res = ''
    xtal_res = ''

    #get the residues
    for jj in res:
        #print(jj)
        xtal_res += ' or resi ' + jj
        pred_res += ' or resi ' + str(int(jj) - diff)
    xtal_res = xtal_res[3:]
    pred_res = pred_res[3:]

    cmd.delete('all')
    cmd.load(pred_file,'pdb1')
    cmd.load(wt_file,'wt1')
    cmd.align('pdb1 and (' + pred_res + ')','wt1 and (' + xtal_res + ')')
    cmd.delete('wt1')
    cmd.save('./TEMP_PRED_ALIGN.pdb')

    return;

def obtain_ligands(WT_PROT='A',WT_CARB='B',PRED_PROT='A',PRED_CARB='B',wt_file='./TEMP_WT.pdb',pred_file='./TEMP_PRED_ALIGN.pdb'):
    '''
    input:
        WT_CARB : (list) all carb chains in WT struct
        PRED_CARB: (list) all carb chains in pred struct
        wt_file  : (str) file location of WT structure
        pred_file: (str) file location of Pred structure
    output:
        temporary files for carb chains
    returns:
        none
    '''

    cmd.delete('all')

    cmd.load(wt_file,'wt2')
    for ii in WT_PROT:
        #print('del wt ' + ii)
        cmd.remove('wt2 and chain ' + ii)
    cmd.save('TEMP_WT_LIG.pdb')

    cmd.delete('all')

    cmd.load(pred_file,'pdb2')
    #print(pred_file)
    if len(PRED_CARB) == 1:
        for ii in PRED_PROT:
            #print('del pred ' + ii)
            cmd.remove('pdb2 and chain ' + ii)
        cmd.save('TEMP_PRED_LIG_' + PRED_CARB[0] + '.pdb')
    else:
        striter = 10
        cmd.delete('all')
        for ii in PRED_CARB:
            striter += 1
            cmd.load(pred_file,'pdb' + str(striter) )
            for kk in PRED_PROT:
                cmd.remove('pdb' + str(striter) + ' and chain ' + kk)
            for jj in PRED_CARB:
                if ii != jj:
                    cmd.remove('pdb' + str(striter) + ' and chain ' + jj)
            cmd.save('TEMP_PRED_LIG_' + ii + '.pdb')
    return;

def get_ligand_coor(structure):
    '''
    input:
        structure (BioPDB structure)
    output:
        res_p: (list) protein residues
        coor_p: (list) protein coordinates
        res_c: (list) carb residues
        coor_c: (list) carb coordinates
    '''

    coor_c = []
    coor_p = []
    res_c = []
    res_p = []

    coor_x = []
    res_x = []

    models = structure.get_models()
    models = list(models)
    for m in range(len(models)):
        chains = list(models[m].get_chains())
        for c in range(len(chains)):
            residues = list(chains[c].get_residues())
            for r in range(len(residues)):
                res = residues[r].get_resname()
                if res == 'HOH':
                    continue;

                atoms = list(residues[r].get_atoms())

                for a in range(len(atoms)):
                    at = atoms[a]

                    if 'H' == at.element: continue;

                    if str(residues[r].get_resname()) in PROT_AA:
                        coor_p.append( at.get_coord() )
                        res_p.append( [ str(residues[r].id[1]).strip(), str(chains[c].id).strip(), str(residues[r].get_resname()), str(at.get_name()) ] )

                    else:
                        coor_c.append( at.get_coord() )
                        res_c.append( [ str(residues[r].id[1]).strip(), str(chains[c].id).strip(), str(residues[r].get_resname()), str(at.get_name()) ] )


    return res_p, coor_p, res_c, coor_c

class glycan():
    """
    Class object for a GLYCAN

    Args:
        coors (arr nx3): coordinates of heavy atoms
        atom_names (arr str): names of the atoms

    Variables:
        name, coor, atom_names
        adj_mat (nxn): One-hot of bonded atoms
        edges (nx?): array of arrays of the non-sparse edge connections
        ring_atom (arr nx[5,6]x1 or nx6x1): defines which atoms are in the ring
    """

    def __init__(self,coor,atom_names,BOND_CUTOFF=1.85):

        self.coor = coor
        #print(len(self.coor))
        self.atom_names = atom_names

        self.BOND_CUTOFF = BOND_CUTOFF


        #initialize empty variables
        self.adj_mat = []
        self.edges = []
        self.ring_atom = []
        self.com = []
        self.ring_atom_plus = []

        self.calc_adjacency()

        ope = []

        for jj in range(len(self.coor)):
            o = self.calc_ring(jj)
            ope.append(o)
            self.calc_adjacency()

        ring = []
        ring.append(ope[0])
        for jj in range(1,len(ope)):
            if type(ope[jj]) == bool:
                continue;

            skip = False
            for kk in range(len(ring)):
                if ring[kk][0] == ope[jj][0]:
                    skip = True
            if skip:
                continue;

            ring.append(ope[jj])

        #print(len(self.coor))
        ring_plus = []
        all_ring_atom = [];
        #print(ring)
        for ii in ring:
            for jj in ii:
                all_ring_atom.append(jj)

        for jj in range(len(ring)):
            ring_plus.append([])
            r = []
            for kk in ring[jj]:
                r.append(self.coor[kk])
            #print(self.coor[kk])
            #print(self.coor)
            d = distance_matrix(r,np.array(self.coor))
            d = d < self.BOND_CUTOFF
            #print(np.shape(d))
            d = np.sum(d,axis=0) >= 1
            #print(np.shape(d))
            for ll in range(len(d)):
                if d[ll] and ll not in ring_plus[jj]:
                    if ll not in all_ring_atom:
                        ring_plus[jj].append(ll)

            for uu in range(10):
                try:
                    r = []
                    for kk in ring_plus[jj]:
                        r.append(self.coor[kk])
                    #print(self.coor[kk])
                    #print(self.coor)
                    d = distance_matrix(r,np.array(self.coor))
                    d = d < self.BOND_CUTOFF
                    #print(np.shape(d))
                    d = np.sum(d,axis=0) >= 1
                    #print(np.shape(d))
                    for ll in range(len(d)):
                        if d[ll] and ll not in ring_plus[jj]:
                            if ll not in all_ring_atom:
                                ring_plus[jj].append(ll)
                except:
                    break;

        self.ring_atom_plus = ring_plus

        self.ring_atom = ring
        self.ring_atom_name, self.ring_com = self.get_ring_atom_name()


    def calc_adjacency(self):
        #get the adjacency matrix and edge list of the carb

        #calculate atom-atom distances and set cutoffs
        dm = distance_matrix(self.coor,self.coor)

        adj_mat = dm < self.BOND_CUTOFF;
        #no self interactions
        for i in range(len(adj_mat)):
            adj_mat[i,i] = 0

        #get the list of the adjacency matrix
        edge_list = [];
        for ii in range(len(adj_mat)):
            edge_list.append([])
            for jj in range(len(adj_mat)):
                if adj_mat[ii,jj]:
                    edge_list[ii].append(jj)

        #store local variables into class variables
        self.adj_mat = adj_mat
        self.edges = edge_list
        return

    #recursive algo to get cycle of the graph
    def visit(self,n,edge_list,visited,st):
        """
        Args:
            n - node we are searching from
            edge_list - adjacency of each node, is periodically
                modified to remove connection to parent coming from
            st - start node
        Returns:
            arr - array of the cycle found
        """

        if n == st and visited[st] == True:
            return [n]

        visited[n] = True
        r = False
        arr = []
        for e in edge_list[n]:
            try:
                edge_list[e].remove(n)
            except:
                continue;

            r = self.visit(e,edge_list,visited,st)

            if type(r) != bool:
                arr.append(n)
                for j in r:
                    arr.append(j)
        if arr == []:
            return False
        return arr

    def calc_ring(self,i):
        #gets the ring atoms, calls recursive visit function
        ring = self.visit(i,copy.deepcopy(self.edges),np.zeros(len(self.coor)),i)
        ind = 0
        while type(ring) == bool:
            ring = self.visit(ind,copy.deepcopy(self.edges),np.zeros(len(self.coor)),ind)
            ind += 1;
            if ind >= len(self.coor):
                break;

        self.ring_atom = np.unique(ring).astype(int)

        return self.ring_atom

    def get_ring_atom_name(self):
        #gets the ring_atom_names in PDB notation and the com of each ring
        r = []
        com = []
        for jj in self.ring_atom:
            r.append([])
            com.append(np.array([0.,0.,0.]))
            for kk in jj:
                r[-1].append(self.atom_names[kk])
                com[-1] += np.array(self.coor[kk])
            com[-1] /= len(r[-1])
        return r, np.array(com)

def find_interactChains(coor_c,coor_p,res_c,res_p,INTERACT=5.0):
    '''
        Finds all protein residues within INTERACT Ang of carbs - for ALIGNMENT
    input:
        coor_c : (list) all coordinates of carbohydrate ligand
        coor_p : (list) all coordinates of protein
        res_c  : (list) residue information of carbohydrate ligand
        res_p  : (list) residue information of protein
        INTERACT: (float) distance to determine interactions
    returns:
        chain_int: (list) all protein residues that interact with the carb
    '''
    #determine chain-chain interactions
    d = distance_matrix(coor_c,coor_p) < INTERACT
    a = np.array( np.where(d == 1) )
    a = np.array(a)

    chain_int = []
    for ii in range(a.shape[1]):
        res2 = res_p[ a[1,ii] ]
        chain_int.append(res2)

    return chain_int

def find_interactRingAtomRes(rcom,gly,pc,pr,INTERACT=6.0):
    '''
        Finds all protein residues within INTERACT Ang of carb residues - For FNAT
    input:
        rcom : (list) all coordinates of carbohydrate residues (rings)
        gly : (glycan instance) glycan information
        pc  : (list) protein coordiantes
        pr  : (list) protein residues
        INTERACT: (float) distance to determine interactions
    returns:
        cint: (list) all protein residues that interact with the individual residues of the carb
    '''
    #determine chain-chain interactions
    cint = []
    for jj in range(len(gly.ring_atom_plus)):
        at = []
        for ii in gly.ring_atom_plus[jj]:
            at.append(gly.coor[ii])
        try:
            d = distance_matrix(at,pc) < INTERACT
            cint.append([])
            for ii in range(len(at)):
                for jj in range(len(pr)):
                    if d[ii,jj]:
                        res2 = int(pr[jj][0])
                        if res2 not in cint[-1]:
                            cint[-1].append(res2)
        except:
            cint.append([])
    return cint

def fnat_full_lig(wt_r,pred_r):
    '''
        Calculates Fnat_full
    input:
        wt_r : (list) residues of experimental (WT) protein that interact with carb
        pred_r : (list) residues of predicted protein that interact with carb
    returns:
        fnat: (float) Fnat full
    '''
    y, y_hat = [], []

    for ii in wt_r:
        if ii[0] not in y:
            y.append(ii[0])
    for ii in pred_r:
        if ii[0] not in y_hat:
            y_hat.append(ii[0])
    y = np.sort(np.array(y).astype(int))
    y_hat = np.sort(np.array(y_hat).astype(int))

    try:
        a = np.max(y_hat)
    except:
        a = 200
    if np.max(y) > a:
        a = np.max(y)

    y_arr = np.zeros(a + 500)
    y_pred_arr = np.zeros(a + 500)

    y_arr[y] = 1
    y_pred_arr[y_hat] = 1

    fnat = np.sum(y_arr * y_pred_arr) / np.sum(y_arr)
    return fnat

def hungarian_fnat(cint,cint_,ba=0):
    #gets Fnat - fraction of natural contacts

    curr = cint
    curr_ = cint_

    #print('fnat_res:')
    #print(curr,'\n\n',curr_)

    f = np.zeros((len(curr),len(curr_)))
    n = np.zeros((len(curr),len(curr_)))

    for ii in range(len(curr_)):

        for jj in range(len(curr_[ii])):
            n[:,ii] += 1;

            for aa in range(len(curr)):
                for bb in range(len(curr[aa])):
                    if curr[aa][bb] + ba == curr_[ii][jj]:
                        f[aa,ii] += 1

    #print('F:')
    #print(f,'\n\nN:')
    #print(n)
    rolling_f = 0;
    rolling_n = 0;

    #print(np.shape(f))

    while True:
        a = np.argmax(f)
        #print(a)
        r = a // len(curr_)
        c = a % len(curr_)
        #print(r,c)

        rolling_f += f[r,c]
        rolling_n += n[r,c]


        skip = False
        curr_no = []

        f[:,c] = 0
        f[r,:] = 0
        n[:,c] = 0
        n[r,:] = 0

        if np.sum(f) < 1:
            break;

    #pick up the remaining ones
    while np.sum(n) > 0:

        a = np.argmax(n)
        r = a // len(curr_)
        c = a % len(curr_)

        rolling_n += np.sum(n[r,c])
        n[:,c] = 0
        n[r,:] = 0

    #print(f,n,f/n)
    #print(rolling_f, rolling_n)
    return rolling_f / rolling_n

def get_all_info(file):
    """
    input:
        file (str): file name string
    return:
        prot_res (arr, str): PDB names of protein CA atoms
        prot_coor (arr, float): coordinates of the CA atoms
        int_res (arr, str): PDB names of protein residues interacting with rings
        ring_atom_name (arr, str): PDB names of glycan ring ATOMS
        ring_com (arr, float): Center of Mass (COM) of all ring atoms
        gly_coor (arr,float): all atom coordinates of the glycan
        gly (class): raw glycan class for further analysis if needed
    """
    parser=PDBParser()
    structure=parser.get_structure("prot", file)
    pr,pc, res, coor = get_ligand_coor(structure)

    #get glycan_info
    at = []
    for ii in res:
        at.append(ii[-1])
    gly = glycan(coor,at)

    #get protein info
    prot_res, prot_coor = [], []
    for ii in range(len(pr)):
        if pr[ii][-1] == 'CA':
            prot_res.append(pr[ii])
            prot_coor.append(pc[ii])

    #get interact info
    int_ = find_interactChains(coor,pc,res,pr)
    cint = find_interactRingAtomRes(gly.ring_com,gly,pc,pr,INTERACT = 5.0)

    return pr, pc, int_, gly.ring_atom_name, gly.ring_com, gly.coor, gly, cint

def hungarian_lrms(gc,gc_,g,g_):
    #Deprecated version to calculate LRMS
    #Still calcuated ; however, is not the appropriate way.
    #code is provided for legacy purposes only

    big_no = []
    rms = []
    dm = distance_matrix(gc,gc_)
    iter = 0;

    while True:
        a = np.argmin(dm)
        r = a // len(gc_)
        c = a % len(gc_)
        skip = False
        curr_no = []
        if g.atom_names[r][0] == g_.atom_names[c][0]:
            rms.append(dm[r,c] ** 2)
            dm[:,c] = 1e10
            dm[r,:] = 1e10
        else:
            curr_no.append(a)
        dm[r,c] = 1e10

        if np.sum(dm < 1e9) < 1:
            break;

    return np.sqrt( np.sum(rms) / len(rms) )

def hungarian_rirms(gc,gc_):
    """
        Calculates the ring-ring RMS using a hungarian approach for ring identity
    input:
        gc (list): glycan coordinates of predicted structure
        gc_ (list): glycan coordinates of experimental structure
    returna:
        (float): Ring-ring RMS (rirms)
    """

    big_no = []
    rms = []
    dm = distance_matrix(gc,gc_)
    iter = 0;

    while True:
        a = np.argmin(dm)
        r = a // len(gc_)
        c = a % len(gc_)

        skip = False
        curr_no = []

        rms.append(dm[r,c] ** 2)
        dm[:,c] = 1e10
        dm[r,:] = 1e10
        dm[r,c] = 1e10

        if np.sum(dm < 1e9) < 1:
            break;

    return np.sqrt( np.sum(rms) / len(rms) )

def fix_num(pr,pr_):
    """
        fixes any possible numbering issues between the predicted and wt structs
    input:
        pr (list): protein residues of predicted struct
        pr_ (list): protein residues of experimentally solved struct
    return:
        pr (list): updated protein residues for predicted structure in alignment with exp struct
        best_adj (int): numeric alignment between pred and exp
    """
    #simplify to the simple CA only
    ca, ca_ = [], []
    for ii in pr:

        if 'CA' in ii[-1]:
            ca.append([int(ii[0]) , ii[2] ])
    for ii in pr_:
        if 'CA' in ii[-1]:
            ca_.append([int(ii[0]) , ii[2] ])

    best_corr = 0
    best_adj = 0

    for ii in range(-1000,1000):
        corr = 0
        new_ca = []
        for jj in ca:
            new_ca.append([jj[0]+ii,jj[1]])
        for jj in new_ca:
            if jj in ca_:
                corr += 1
        if corr > best_corr:
            best_corr = corr
            best_adj = ii
        if corr > len(ca_) - 10:
            best_corr = corr
            best_adj = ii
            break;

    my_d = []
    for ii in range(len(pr)):
        pr[ii][0] = str( int(pr[ii][0]) + best_adj )

    return pr, best_adj

def get_clash(pc,gc,vdw=1.85):
    """
        calculates number of clashes between carb and prot.
        Not reported in paper - because it was so low (except RFAA)
    input:
        pc (list): protein coordinates
        gc (list): glycan coordiantes
        vdw (float): Van Der Waals radius
    return:
        n_clash (int): number of clashes
    """
    dm = distance_matrix(gc,pc)
    dm = dm < vdw
    n_clash = np.sum(dm)
    return n_clash

def get_sc_lrms(ref_file, mol_file):
    """
        calculates LRMS the correct way using RDKIT
    input:
        ref_file (str): experimental ligand structure
        mol_file (str): predicted ligand structure
    return:
        (float): LRMS
    """
    # Load molecules
    ref_mol = Chem.MolFromPDBFile(ref_file) #returns mol obj
    mol_mol = Chem.MolFromPDBFile(mol_file,sanitize=False)
    mcs = rdFMCS.FindMCS([ref_mol, mol_mol])
    mcs_smarts = mcs.smartsString
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)
    # Get the atom indices of the MCS in both molecules
    ref_match = ref_mol.GetSubstructMatch(mcs_mol)# returns a tuple of integers mol.GetAtoms()
    mol_match = mol_mol.GetSubstructMatch(mcs_mol)
    mmap = list(zip(mol_match, ref_match))
    # Calculate the RMSD for the MCS atoms
    return Chem.rdMolAlign.CalcRMS(mol_mol, ref_mol,map=[mmap])#, map=[list(mol_match),[list(ref_match)]])


def calc_metrics(decoy,native,same_ligand=True,is_align=True,is_same_num=True):
    """
    input:
        decoy (str): file name string of predicted structure
        native (str): file name string of native structure
        same_ligand (bool): if the ligand used is longer than the native ligand then False
    return:
        d (float): Dice of the prediction
        rirms (float): Ring RMS
        lrms (float): Ligand RMS
        dockq (float): dockq score
        s (str): string of d,rirms,lrms,dockq for easy printing
    """

    #print('o')

    pr, pc, i, ran, rcom, gc, g, cint = get_all_info(decoy)
    pr_, pc_, i_, ran_, rcom_, gc_, g_, cint_ = get_all_info(native)

    #for readability
    for ii in range(len(cint)):
        cint[ii] = list(np.sort(cint[ii]))
    for ii in range(len(cint_)):
        cint_[ii] = list(np.sort(cint_[ii]))

    #nonredundant residues of binding pocket
    nrr = []
    #print(i_)
    for ii in i:
        if int(ii[0]) not in nrr:
            nrr.append(int(ii[0]))

    nrr.sort()

    ab_clash = get_clash(pc,gc,vdw=VDW_CUT) // 1
    aa_clash = (get_clash(gc,gc,vdw=1) - len(gc)) // 2

    ba = 0
    if is_same_num == False:
        pr, ba = fix_num(pr,pr_)

    o = ''
    o += decoy + ','
    for ii in nrr:
        o += str(ii) + '|'

    rres = []
    rres_ = []

    rirms = hungarian_rirms(rcom,rcom_)

    lrms = hungarian_lrms(gc,gc_,g,g_)

    f_res = hungarian_fnat(cint,cint_,ba=ba)
    f_full = fnat_full_lig(i_,i)

    return f_full,f_res,lrms,rirms,ab_clash,aa_clash
