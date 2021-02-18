#!/usr/bin/python
#!/home/shaikh/Downloads/pychem-1.0
import sys
from rdkit import Chem
from rdkit.Chem.EState import Fingerprinter
from rdkit.Chem import Descriptors
from rdkit.Chem.rdmolops import RDKFingerprint
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
#from sklearn import cross_validation
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_score
from keras import regularizers
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import Ridge, LinearRegression
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.model_selection import GridSearchCV

from pychem.pychem import Chem
from pychem import pychem
from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from pychem import fingerprint

from rdkit.Chem.rdMolDescriptors import GetHashedAtomPairFingerprintAsBitVect,GetHashedTopologicalTorsionFingerprintAsBitVect
from rdkit.Chem.AtomPairs.Sheridan import GetBPFingerprint
from rdkit.Chem.EState.Fingerprinter import FingerprintMol
from rdkit.Avalon.pyAvalonTools import GetAvalonFP #GetAvalonCountFP  #int vector version
from rdkit.Chem.AllChem import  GetMorganFingerprintAsBitVect, GetErGFingerprint
from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray
import rdkit.DataStructs.cDataStructs
from sklearn.metrics import r2_score
from numpy.random import randn
from numpy.random import seed
from scipy.stats import pearsonr
from sklearn.model_selection import LeaveOneOut
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pickle
from keras.models import load_model
import csv 
import subprocess
import pickle
from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from pychem import pychem
from pychem import fingerprint 
import os
#db = ['Anticogulant','Beta_secretase','Bromodomain','Carbonic_Anhydrase','Cholera','Diabetes','Estrogen','HCV','HIV','Hpyroli','Hydrolase','Influenza','Isomerase','Peroxisome','Kinase','Ligase','Transferase','Tuberculosis']
#db = ["Anticogulant","Apoptosis","Beta_secretase","Bromodomain","Carbonic_Anhydrase","Cellcycle_inhibitor","Cholera","Diabetes","Estrogen","HCV","Heat_Shock_protein","HIV","Hpyroli","Hydrolase","Immune_system","Influenza","Isomerase","Kinase","Ligase","Lyase","Oxidoreductase","Peroxisome","Protein_binding","Signaling_protein","Transcription","Transferase","Transport_protein","Tuberculosis"]
db = ["Protein_binding","Signaling_protein","Transcription","Transferase","Transport_protein","Tuberculosis"]
for sdb in db:
	count_line = 0
	if os.path.isdir(sdb) == False:
		os.makedirs(sdb)
	with open('../../'+sdb+'.csv','r') as input_file:
		for line in input_file:
			if count_line != 0:
				line = line.rstrip()
				line = line.split(';')
				pdbid = line[0]
				print sdb,pdbid
				if True:
					inp = line[1]
					#factory = Gobbi_Pharm2D.factory 
					mol1=pychem.Chem.MolFromSmiles(inp)
					AllChem.EmbedMolecule(mol1,useRandomCoords=True)	
					#FP_inp=fingerprint.CalculateMACCSFingerprint(mol1)
					FP_inp = fingerprint.CalculateDaylightFingerprint(mol1)
					#print type(FP_inp)
					#print type(test)
					with open(sdb+'/'+pdbid+'.dat','r') as dl_file:
						data_daylight_bin = dl_file.readline().rstrip()
						data_daylight = DataStructs.CreateFromBitString(data_daylight_bin)

					Tanifing_Daylight = fingerprint.CalculateSimilarity(FP_inp[2],data_daylight,'Tanimoto')

					print Tanifing_Daylight
					break

					print 'SUCCESS'
			count_line += 1
	break
