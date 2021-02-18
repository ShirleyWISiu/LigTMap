"""
    Module checks interactions between two molecules and
    creates interacion fingerprints.

"""
from __future__ import division
from six.moves import zip_longest
from itertools import chain
from collections import OrderedDict
import numpy as np
from scipy.sparse import csr_matrix, isspmatrix_csr
import oddt
from oddt.utils import is_openbabel_molecule
from oddt.interactions import (pi_stacking,
                               hbond_acceptor_donor, hbonds,
                               salt_bridge_plus_minus,
                               hydrophobic_contacts,
                               acceptor_metal,
                               close_contacts)
from oddt import toolkit
from oddt import shape
from oddt import fingerprints
from rdkit.Chem import Draw
from oddt import interactions
import sys
import os
#in_num = sys.argv[1]
#file_name = 'index/index_'+in_num
sdb = sys.argv[1]
 
__all__ = ['InteractionFingerprint',
           'SimpleInteractionFingerprint',
           'SPLIF',
           'similarity_SPLIF',
           'ECFP',
           'PLEC',
           'dice',
           'tanimoto',
           'close_contacts',
           'hbond_acceptor_donor',
           'hbonds',
           'halogenbond_acceptor_halogen',
           'halogenbonds',
           'pi_stacking',
           'salt_bridge_plus_minus',
           'salt_bridges',
           'hydrophobic_contacts',
           'pi_cation',
           'acceptor_metal',
           'pi_metal']


#os.makedirs(sdb)
#error_file = open(sdb+'/error.dat','a')
with open('/home/user/share/tid/TID_WWW_RUNS/Program/Index/IFP/'+sdb+'/error.dat','r') as index_file:
	for line in index_file:
		pid = line[0:4]
		crystal_protein = next(oddt.toolkit.readfile('pdb', '/home/user/share/tid/TID_WWW_RUNS/Program/docking_files/protein/'+pid+'_protein.pdb'))
		crystal_protein.protein = True
		crystal_ligand = next(oddt.toolkit.readfile('pdb', '/home/user/share/tid/TID_WWW_RUNS/Program/docking_files/ligand/'+pid+'_ligand.pdb'))

		IFP = fingerprints.InteractionFingerprint(crystal_ligand, crystal_protein)
		IFP.tofile(sdb+'/'+pid+'.bin')
