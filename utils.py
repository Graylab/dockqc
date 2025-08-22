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

a1 = 1.5
a2 = 8.5

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
#print(PROT_AA)
#allow secondary conformations (A + B)
for ii in ope:
    PROT_AA.append('A' + ii)
    PROT_AA.append('B' + ii)

#print(PROT_AA)

#Determine which ones are glycosylated and which are non-covalent


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
    cmd.load(PRED,'pdb1')
    cmd.load(WT,'wt1')
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
            cmd.load(PRED,'pdb' + str(striter) )
            for ii in PRED_PROT:
                cmd.remove('pdb' + str(striter) + ' and chain ' + ii)
            for jj in PRED_CARB:
                if ii != jj:
                    cmd.remove('pdb' + str(striter) + ' and chain ' + ii)
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

        #for jj in range(len(ring)):
        #    print('a',ring[jj],'\n','p',ring_plus[jj])


        self.ring_atom = ring
        self.ring_atom_name, self.ring_com = self.get_ring_atom_name()

        #print(len(self.coor),'\n\n\n')
        #self.print_variables()

    def calc_adjacency(self):
        #get the adjacency matrix and edge list of the carb

        #calculate atom-atom distances and set cutoffs
        dm = distance_matrix(self.coor,self.coor)
        #print(dm)
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
        #print(edge_list)
        #print(n)

        if n == st and visited[st] == True:
            return [n]

        visited[n] = True
        r = False
        arr = []
        #print(n,edge_list[n],visited)

        for e in edge_list[n]:
            #if n in edge_list[e]:
            try:
                edge_list[e].remove(n)
            except:
                continue;
            #print('\t',e)

            r = self.visit(e,edge_list,visited,st)
            #print('\t\t',r)

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
            #print('ringboi',ind,len(self.coor))
            ring = self.visit(ind,copy.deepcopy(self.edges),np.zeros(len(self.coor)),ind)
            ind += 1;
            if ind >= len(self.coor):
                break;

        #print(ring)
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

def get_rings(file):
    """
    input:
        file (str): file name string
    output:
        ring_atom_name (arr, str): PDB names of glycan ring ATOMS
        ring_com (arr, float): Center of Mass (COM) of all ring atoms
        gly (class): raw glycan class for further analysis if needed
    """
    parser=PDBParser()
    structure=parser.get_structure("prot", file)
    _,_, res, coor = get_ligand_coor(structure)
    at = []
    for ii in res:
        at.append(ii[-1])
    #print(at)
    gly = glycan(coor,at)
    return gly.ring_atom_name, gly.ring_com, gly

def find_interactChains(coor_c,coor_p,res_c,res_p,INTERACT=5.0):
    #determine chain-chain interactions
    d = distance_matrix(coor_c,coor_p) < INTERACT
    a = np.array( np.where(d == 1) )
    a = np.array(a)
    #print(a)

    #chain_int = {}
    chain_int = []
    for ii in range(a.shape[1]):
        #res1 = res_c[ a[0,ii] ]
        res2 = res_p[ a[1,ii] ]
        chain_int.append(res2)

    #print(np.shape(chain_int))
    #print(np.unique( np.array(chain_int)[:,0] ))

    return chain_int

def find_interactRingRes(rcom,pc,pr,INTERACT=6.0):
    #determine chain-chain interactions
    d = distance_matrix(rcom,pc) < INTERACT

    cint = []
    for ii in range(len(rcom)):
        cint.append([])
        for jj in range(len(pr)):
            if d[ii,jj]:
                res2 = int(pr[jj][0])
                if res2 not in cint[-1]:
                    cint[-1].append(res2)
    return cint

def find_interactRingAtomRes(rcom,gly,pc,pr,INTERACT=6.0):
    #determine chain-chain interactions
    cint = []
    for jj in range(len(gly.ring_atom_plus)):

        #print(gly.ring_atom[jj])

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
    y, y_hat = [], []

    for ii in wt_r:
        if ii[0] not in y:
            y.append(ii[0])
            #print('_',ii)
    for ii in pred_r:
        if ii[0] not in y_hat:
            y_hat.append(ii[0])
            #print('\t',ii)
    y = np.sort(np.array(y).astype(int))
    y_hat = np.sort(np.array(y_hat).astype(int))

    #print(y_hat,'\n')
    #print(y)

    #print(y,'\n',y_hat)
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

    #print(np.sum( y_arr * y_pred_arr), np.sum(y_arr), np.sum(y_pred_arr) )

    d = np.sum(y_arr * y_pred_arr) / np.sum(y_arr)
    return d

def fnat(cint,cint_):
    #gets Fnat - fraction of natural contacts
    f ,n = 0,0
    for ii in range(len(cint_)):
        for jj in range(len(cint_[ii])):
            n += 1
            if ii < len(cint):
                if cint_[ii][jj] in cint[ii]:
                    f += 1
    #print(f,n,f/n)
    return f / n

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


def adjusted_fnat(cint,cint_):
    #gets Fnat - fraction of natural contacts
    f = 0;
    for ii in range(-1,len(cint)):
        new_c = []
        for jj in range(len(cint)):
            new_c.append(cint[(jj + ii) % len(cint)])
        o = fnat(new_c,cint_)
        if o > f:
            f = o

    return f

def rms(x,y):
    if len(x) != len(y):
        #print('lengths of ligands do not match')
        return -1
    r = 0
    for ii in range(len(x)):
        r += np.linalg.norm(x[ii] - y[ii]) ** 2
    r /= len(x)
    return np.sqrt(r)

def get_all_info(file):
    """
    input:
        file (str): file name string
    output:
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

    #print(coor)

    #get glycan_info
    at = []
    for ii in res:
        at.append(ii[-1])
    #print(at)
    gly = glycan(coor,at)

    #get protein info
    prot_res, prot_coor = [], []
    for ii in range(len(pr)):
        if pr[ii][-1] == 'CA':
            prot_res.append(pr[ii])
            prot_coor.append(pc[ii])

    #get interact info
    int_ = find_interactChains(coor,pc,res,pr)
    #print('ope:',int_)




    cint = find_interactRingAtomRes(gly.ring_com,gly,pc,pr,INTERACT = 5.0)

    #print(file, cint)

    return pr, pc, int_, gly.ring_atom_name, gly.ring_com, gly.coor, gly, cint

def adjusted_lrms(gc,gc_,g,g_):

    r = 1e+10
    big_no = []

    for jj in range(len(gc_)):
        dm = distance_matrix(gc,gc_)
        no = []
        ope = []
        curr_no = []
        u_gc_ = []
        for kk in range(len(gc_)):
            ii = (jj + kk) % len(gc_)
            #print(ii)

            a = np.argmin(dm[:,ii])

            #print(dm[:,ii])

            skip = False
            curr_no = []
            if g.atom_names[a][0] == g_.atom_names[ii][0]:
                no.append(a)
                ope.append([a,ii])
                dm[:,ii] = 1e10
                dm[a,:] = 1e10
                skip = True
                u_gc_.append(gc_[ii])
            else:
                curr_no.append(a)

            cnt = 0
            while a in no or a in curr_no:
                #print('\t',a)
                if skip:
                    break;

                dm[a,ii] = 1e10
                a = np.argmin(dm[:,ii])

                if g.atom_names[a][0] == g_.atom_names[ii][0]:
                    no.append(a)
                    ope.append([a,ii])
                    u_gc_.append(gc_[ii])
                    dm[:,ii] = 1e10
                    dm[a,:] = 1e10
                    skip = True
                    break;
                else:
                    curr_no.append(a)

                cnt += 1
                #escape if something is broken
                if cnt > len(gc) + 50:
                    no.append(-1)
                    ope.append([-1,-1])
                    break;

        cgc = []
        for kk in no:
            if kk != -1:
                cgc.append(gc[kk])
        cgc = np.array(cgc)
        #print(no)
        crms = rms(cgc,np.array(u_gc_))
        if crms < r:
            r = crms
            big_no = copy.deepcopy(no)

    o = []
    for ii in big_no:
        if ii != -1:
            o.append(ii)
        #print(jj,crms)
    return r, o

def hungarian_lrms(gc,gc_,g,g_):

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

    big_no = []
    rms = []
    dm = distance_matrix(gc,gc_)
    iter = 0;

    #print(dm)

    while True:
        a = np.argmin(dm)
        r = a // len(gc_)
        c = a % len(gc_)

        #print(r,c,round(dm[r,c]))

        skip = False
        curr_no = []

        rms.append(dm[r,c] ** 2)
        dm[:,c] = 1e10
        dm[r,:] = 1e10
        dm[r,c] = 1e10

        #print(dm,'\n\n')

        if np.sum(dm < 1e9) < 1:
            break;

    return np.sqrt( np.sum(rms) / len(rms) )

def adjusted_rrms(gc,gc_):

    r = 1e+10
    nrc = []

    for jj in range(len(gc_)):
        dm = distance_matrix(gc,gc_)
        no = []
        ope = []
        curr_no = []
        for kk in range(len(gc_)):
            ii = (jj + kk) % len(gc_)
            #print(ii)

            a = np.argmin(dm[:,ii])

            skip = False
            curr_no = []

            no.append(a)
            ope.append([a,ii])
            dm[:,ii] = 1e10
            dm[a,:] = 1e10


        cgc = []
        for kk in no:
            cgc.append(gc[kk])
        cgc = np.array(cgc)
        #print(no)
        crms = rms(cgc,gc_)
        if crms < r:
            r = crms
            nrc = cgc

        #print(jj,crms)
    return r, nrc

def fix_num(pr,pr_):
    #simplify to the simple CA only
    ca, ca_ = [], []
    for ii in pr:

        if 'CA' in ii[-1]:
            #print(ii)
            ca.append([int(ii[0]) , ii[2] ])
    for ii in pr_:

        if 'CA' in ii[-1]:
            #print('\t',ii)
            ca_.append([int(ii[0]) , ii[2] ])

    #print(ca[:10],'\n',ca_[:10])
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

    #print(best_adj)

    return pr, best_adj

def fnat_dice(cint,cint_):
    #print(cint)
    #print(cint_)
    c, c_ = [], []
    for ii in cint:
        for jj in ii:
            c.append(jj)
    for ii in cint_:
        for jj in ii:
            c_.append(jj)
    y = np.unique(np.array(c_).astype(int))
    y_hat = np.unique(np.array(c).astype(int))

    a = np.max(y_hat)
    if np.max(y) > a:
        a = np.max(y)

    y_arr = np.zeros(a + 500)
    y_pred_arr = np.zeros(a + 500)

    y_arr[y] = 1
    y_pred_arr[y_hat] = 1

    d = 2 * np.sum(y_arr * y_pred_arr) / (np.sum(y_arr) + np.sum(y_pred_arr))
    return d

def get_clash(pc,gc,vdw=1.85):
    dm = distance_matrix(gc,pc)
    dm = dm < vdw
    n_clash = np.sum(dm)
    return n_clash

def get_sc_lrms(ref_file, mol_file):
    # Load molecules
    ref_mol = Chem.MolFromPDBFile(ref_file) #returns mol obj
    #print("hi\n",ref_mol)
    mol_mol = Chem.MolFromPDBFile(mol_file,sanitize=False)
    #print("dddddd")
    mcs = rdFMCS.FindMCS([ref_mol, mol_mol])
    #print("GOTTA CATCHEM ALLL")
    mcs_smarts = mcs.smartsString
    #print(mcs_smarts)
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)
    # Get the atom indices of the MCS in both molecules
    ref_match = ref_mol.GetSubstructMatch(mcs_mol)# returns a tuple of integers mol.GetAtoms()
    mol_match = mol_mol.GetSubstructMatch(mcs_mol)
    mmap = list(zip(mol_match, ref_match))
    # Calculate the RMSD for the MCS atoms
    #print('oof')
    return Chem.rdMolAlign.CalcRMS(mol_mol, ref_mol,map=[mmap])#, map=[list(mol_match),[list(ref_match)]])
    #takes symmetry into account! maybe try without atomIds


def calc_metrics(decoy,native,same_ligand=True,is_align=True,is_same_num=True):
    """
    input:
        decoy (str): file name string of predicted structure
        native (str): file name string of native structure
        same_ligand (bool): if the ligand used is longer than the native ligand then False
    output:
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

    prms = get_prot_rmsd(pc,pr,pc_,pr_)

    rirms = hungarian_rirms(rcom,rcom_)

    lrms = hungarian_lrms(gc,gc_,g,g_)

    f_res = hungarian_fnat(cint,cint_,ba=ba)
    f_full = fnat_full_lig(i_,i)

    return f_full,f_res,lrms,rirms,ab_clash,aa_clash
