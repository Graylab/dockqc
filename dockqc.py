import os
import numpy as np
import pandas as pd

from Bio.PDB import *
from scipy.spatial import distance_matrix
from utils import *

import os
import numpy as np
import pandas as pd
import argparse

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

from pymol import cmd

PRED = "./ex/pred.cif"
WT = './ex/wt.pdb'
WT_PROT = ['B']
WT_CARB = ['H']

PRED_PROT = ['A']
PRED_CARB = ['B']

def dockqc(fnat_res, fnat_full, lrms, rirms, d1 = 4.0, d2=2.0, fnat_res_pref=0.5):

    #print(eps)
    lrm_scaled = 1 / (1 + (lrms / d1)**2)

    rrm_scaled = 1 / (1 + (rirms / d2)**2)

    fnat = fnat_res_pref * fnat_res + (1-fnat_res_pref) * fnat_full

    dqc = (fnat + rrm_scaled + lrm_scaled) / 3
    return dqc


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('-x', '--exp_file', default='./ex/wt.pdb', type=str, help="experimental TRUE structure")
    parser.add_argument('-p', "--pred_file", default='./ex/pred.pdb', type=str, help="predicted structure to compare to")
    parser.add_argument('-xp', '--exp_prot', default='A', type=str,  help="experimental protein chain IDs (comma seperated)")
    parser.add_argument('-xc', "--exp_carb", default='B', type=str, help="experimental carb chain IDs (comma seperated)")
    parser.add_argument('-pp', '--pred_prot', default='A', type=str, help="predicted protein chain IDs (comma seperated)")
    parser.add_argument('-pc', '--pred_carb', default='B', type=str, help="predicted carb chain IDs (comma seperated)")

    global Options
    Options = parser.parse_args()
    #parser.print_help()

    WT = Options.exp_file.split(',')[0]
    PRED = Options.pred_file.split(',')[0]

    WT_PROT = Options.exp_prot.split(',')
    WT_CARB = Options.exp_carb.split(',')
    PRED_PROT = Options.pred_prot.split(',')
    PRED_CARB = Options.pred_carb.split(',')

    #print('converting to pdbs')
    convert_to_pdb(PRED=PRED,WT=WT)

    #print('getting interacting residues of WT')
    wt_res, res_diff = get_wt_res()

    #print('aligning predicted to wild type')
    align_pred(wt_res, res_diff)

    #print('seperating ligands')
    obtain_ligands(WT_PROT=WT_PROT,WT_CARB=WT_CARB,PRED_PROT=PRED_PROT,PRED_CARB=PRED_CARB)

    #print('doing dockqc metrics')
    f_full,f_res,lrms,rirms,ab_clash,aa_clash = calc_metrics(
        decoy='./TEMP_PRED_ALIGN.pdb', native='./TEMP_WT.pdb',is_align=True,is_same_num=False)

    scrms = []
    for jj in PRED_CARB:
        scrms.append(get_sc_lrms('./TEMP_WT_LIG.pdb', './TEMP_PRED_LIG_' + jj + '.pdb'  ) )
        #print(scrms)
    lrms = np.mean(scrms)

    dqc = dockqc(f_res, f_full, lrms, rirms)

    #delete temp files
    ls = os.listdir('./')
    for ii in ls:
        if 'TEMP_' in ii and '.pdb' in ii:
            os.remove(ii)

    print('NATIVE:   \t' + WT)
    print('PREDICTED:\t' + PRED)
    print('----------')
    print('Fnat_full:\t',round(f_res,3))
    print('Fnat_res :\t',round(f_res,3))
    print('lrms :    \t',round(lrms,3))
    print('rirms:    \t',round(rirms,3))
    print('----------')
    print('DockQC:   \t',round(dqc,3))

    print('\n\nFin.')


if __name__ == '__main__':
    main()
