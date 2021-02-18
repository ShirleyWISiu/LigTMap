from __future__ import division
import sys

def main_part():
	smile_inf_name = sys.argv[1]
	input_num = sys.argv[2]
	tanifing_str = '0.4'
	select_db = sys.argv[3]
	total_db = int(sys.argv[4])
	rootpa = sys.argv[5]
	output_path = "Output/Input_"+input_num+"/"+select_db
	summary_path = "Output/Input_"+input_num

	smile_inf = open(smile_inf_name,'r')
	smile = smile_inf.readline().rstrip()
	smile_inf.close()
	inp = smile
	tanifingcut = float(tanifing_str)
	output_list = []

	try:
		status_f = open('Output/status.txt','a')
                import rdkit
                import openbabel
                import pychem
                import pybel
                import oddt
		from rdkit import Chem
		from rdkit.Chem.EState import Fingerprinter
		from rdkit.Chem import Descriptors
		from rdkit.Chem.rdmolops import RDKFingerprint
		import pandas as pd
		import numpy as np
		from sklearn.preprocessing import StandardScaler

		from pychem.pychem import Chem
		from pychem import pychem
		from rdkit import Chem, DataStructs, RDConfig
		from rdkit.Chem import AllChem
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
		import csv 
		import subprocess
		import pickle
		import copy

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
		import os

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

	
		print "Step 1: Ligand Similarity Search"	
		#### PART1: Similarity <--> Tanifing + Tanipharm Score ####
	
		try:	
			# Generate the fingerprint
			#factory = Gobbi_Pharm2D.factory 
			mol1=pychem.Chem.MolFromSmiles(inp)
			AllChem.EmbedMolecule(mol1,useRandomCoords=True)	
			#ph1=Generate.Gen2DFingerprint(mol1,factory,dMat=Chem.Get3DDistanceMatrix(mol1))
			FP_inp_Morgan=fingerprint.CalculateMorganFingerprint(mol1,radius=2)
			FP_inp_MACCSF = fingerprint.CalculateMACCSFingerprint(mol1)
			FP_inp_Daylight = fingerprint.CalculateDaylightFingerprint(mol1)
		except:
			smile_err_file = open(summary_path+'/smile_err.dat','w')
			smile_err_file.close()

		id_ligname = {}
		id_proname = {}
		with open(rootpa+"Index/"+select_db+".csv",'r') as index_file:
			label_line = index_file.readline()
			for y in index_file:
				try:
					y = y.rstrip()
					data_pdb = y.split(';')[0]
					data_smile = y.split(';')[1]
					affinity = y.split(';')[2]
					id_proname[data_pdb] = y.split(';')[3]
					id_ligname[data_pdb] = y.split(';')[4]
					# Fing part
					with open(rootpa+'Index/fing/Morgan/'+select_db+'/'+data_pdb+'.bin','r') as data_morgan_file:
						data_morgan_bin = data_morgan_file.readline().rstrip()
						data_morgan = DataStructs.UIntSparseIntVect(data_morgan_bin)

					with open(rootpa+'Index/fing/MACCSF/New/'+select_db+'/'+data_pdb+'.dat','r') as data_maccsf_file:
						data_maccsf_bin = data_maccsf_file.readline().rstrip()
						#data_fing_MACCSF = DataStructs.ExplicitBitVect(data_fing_bin_MACCSF)
						data_maccsf = DataStructs.CreateFromBitString(data_maccsf_bin)

					# Calculate the Tani Score
					Tanifing_Morgan = fingerprint.CalculateSimilarity(FP_inp_Morgan[2],data_morgan,'Tanimoto')
					Tanifing_MACCSF = fingerprint.CalculateSimilarity(FP_inp_MACCSF[2],data_maccsf,'Tanimoto')
					
					if total_db > 1:
						# Daylight part
						with open(rootpa+'Index/fing/Daylight/'+select_db+'/'+data_pdb+'.dat','r') as data_daylight_file:
							data_daylight_bin = data_daylight_file.readline().rstrip()
							data_daylight = DataStructs.CreateFromBitString(data_daylight_bin)
						Tanifing_Daylight = fingerprint.CalculateSimilarity(FP_inp_Daylight[2],data_daylight,'Tanimoto')
						sum_score= (Tanifing_Morgan+Tanifing_MACCSF+Tanifing_Daylight)/3
					else:
						sum_score= (Tanifing_Morgan+Tanifing_MACCSF)/2

					if sum_score >= tanifingcut:
						row = [inp,data_smile,data_pdb,affinity,str(Tanifing_Morgan),str(Tanifing_MACCSF),str(sum_score)]
						output_list.append(row)				
				except:
					error_file = open(output_path+"/error.txt",'a')
					error_file.write(data_smile+'\n')
					error_file.close()
			
		#### PART1 COMPLETED ####

		print "Step 2: Docking"
		#### PART2: Docking ####

		# Prepare the input ligand
		#command = os.environ.get('OBABEL')+"/bin/obabel -:"+smile+" -opdb -O "+output_path+"/input.pdb --gen3d" #utako
                command = "obabel -:"+smile+" -opdb -O "+output_path+"/input.pdb --gen3d "
                subprocess.check_output(command.split())

		command = os.environ.get('MGLTools')+"/bin/pythonsh "+os.environ.get('MGLTools')+"/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l "+output_path+"/input.pdb -o "+output_path+"/input.pdbqt -A 'hydrogens' -U 'nphs_lps_waters'"
		subprocess.check_output(command.split())

		# Get the PDBID
		out_file=open(output_path+'/file_list','w')

		pid_list = []
		for ele in output_list:
			out_file.write(ele[2]+'\n')
			pid_list.append(ele[2])

		out_file.close()

		# Docking
		command = rootpa+"docking_files/dock.sh pso "+output_path+" 1 1 "+rootpa
		subprocess.check_call(command.split())   ##utako
                ##subprocess.check_output(command.split())   ##utako
                
		# Write the docking score to the output.csv
		score_file=open(output_path+'/DOCK_LOG/score_1.dat','r')
		score_list = score_file.readlines()
		score_file.close()

		out_file=open(output_path+'/output.csv','w')
		out_file.write("Input;Data;PDB;Affinity;Tani_Morgan;Tani_MACCSF;LigandScore;ILDScore\n")

		i = 0
		for line in output_list:
			out_line = ';'.join(line)
			idscore = score_list[i].rstrip()
			if idscore != '':
				out_line = out_line +';'+idscore+'\n'
			else:
				out_line = out_line + ';n.a.\n'
			out_file.write(out_line)
			i = i + 1

		out_file.close()

		# Create complex file
		for pdbid in pid_list:
				ligand_filename = output_path+'/DOCK_LOG/'+pdbid+'/'+pdbid+'_ligand_1.pdb'
				protein_filename = rootpa+'docking_files/protein/'+pdbid+'_protein.pdb'
				filenames = [protein_filename,ligand_filename]
				output_filename = output_path+'/Complex/complex_'+pdbid+'.pdb'
				with open(output_filename,'w') as outfile:
					for fname in filenames:
						if os.path.isfile(fname):
							with open(fname,'r') as infile:
								outfile.write(infile.read())
		#### PART2 COMPLETED ####

		print "Step 3: Activity Prediction"
		#### PART3: Activity prediction ####
		# Read the data, affinity prediction
		csv_name = output_path+'/output.csv'
		#data = pd.read_table(csv_name, sep=';', header= [0])  utako
		data = pd.read_csv(csv_name, sep=';', header= [0])

		# Add some new columns
		data['Mol'] = data['Input'].apply(Chem.MolFromSmiles)
		num_mols = len(data)

		def MorganFingerprint(mol):
		   	return FingerprintMol(mol)[0]

		# Scale X to unit variance and zero mean
		data['Fingerprint'] = data['Mol'].apply(MorganFingerprint)

		X = np.array(list(data['Fingerprint']))

		st = StandardScaler()
		X = np.array(list(data['Fingerprint']))
		Test = X

		def ExplicitBitVect_to_NumpyArray(bitvector):
			bitstring = bitvector.ToBitString()
			intmap = map(int, bitstring)
			return np.array(list(intmap))


		class fingerprint():
			def __init__(self, fp_fun, name):
				self.fp_fun = fp_fun
				self.name = name
				self.x = []
	
			def apply_fp(self, mols):
				for mol in mols:
					fp = self.fp_fun(mol)
					if isinstance(fp, tuple):
						fp = np.array(list(fp[0]))
					if isinstance(fp, rdkit.DataStructs.cDataStructs.ExplicitBitVect):
						fp = ExplicitBitVect_to_NumpyArray(fp)
					if isinstance(fp,rdkit.DataStructs.cDataStructs.IntSparseIntVect):
						fp = np.array(list(fp))

					self.x += [fp]

					if (str(type(self.x[0])) != "<class 'numpy.ndarray'>"):
						print("WARNING: type for ", self.name, "is ", type(self.x[0]))

		# Load the model
		with open(rootpa+"Model/"+select_db+".sav", "rb") as f:
			Model = pickle.load(f)

		predictions = Model.predict(Test)

		'''
		# Write the activity score to the summary.csv (sort)
		out_filename = summary_path+'/summary.csv'

			#if True:
		if predictions[0] >= 4:		
			if os.path.isfile(out_filename) != True:
				out_file = open(out_filename,'w')
				out_file.write("Model;PIC50")
				out_file.write("\n"+select_db+';'+str(predictions[0]))
				out_file.close()
			else:
				tmp_list = []
				with open(out_filename,'r') as tmp_file:
					label_line = tmp_file.readline()
					for line in tmp_file:
						line = line.rstrip()
						line = line.split(';')
						tmp_list.append([line[0],line[1]])
				tmp_list.append([select_db,str(predictions[0])])
				tmp_list.sort(key=lambda s: s[1],reverse=True)
				out_file = open(out_filename,'w')
				out_file.write(label_line)
				for ele in tmp_list:
					out_file.write(';'.join(ele)+'\n')
				out_file.close()
		'''
		#### PART3 COMPLETED ####		

		countline = -1
		with open(output_path+'/output.csv','r') as countlinefile:
			for line in countlinefile:
				countline += 1

		print "Step 4: Binding Similarity Search (Total "+str(countline)+" target pdb, please wait...)"
		#### PART4: Similarity <--> Binding fingerprint ####

		def tanimoto(a, b, sparse=False):
			if sparse:
				a = np.unique(a)
				b = np.unique(b)
				a_b = float(len(np.intersect1d(a, b, assume_unique=True)))
				denominator = len(a) + len(b) - a_b
				if denominator > 0:
					return a_b / denominator
			else:
				a = a.astype(bool)
				b = b.astype(bool)
				a_b = (a & b).sum().astype(float)
				denominator = a.sum() + b.sum() - a_b
				if denominator > 0:
					return a_b / denominator
			return 0.

		with open(output_path+'/output.csv','r') as first_result_table:
			ignore_label = first_result_table.readline()
			for line in first_result_table:
				line_list = line.rstrip().split(';')
				r_pdbid = line_list[2]
				r_ligandscore = line_list[6]
				r_dockscore = line_list[7]

				if (os.path.isfile(rootpa+'Index/IFP/'+select_db+'/'+r_pdbid+'.bin') == True and os.path.isfile(output_path+'/DOCK_LOG/'+r_pdbid+'/'+r_pdbid+'_ligand_1.pdb') == True):
					crystal_IFP = np.fromfile(rootpa+'Index/IFP/'+select_db+'/'+r_pdbid+'.bin',dtype=np.uint8)
					bind_protein = next(oddt.toolkit.readfile('pdb',rootpa+'docking_files/protein/'+r_pdbid+'_protein.pdb'))
					bind_protein.protein = True
					bind_ligand = next(oddt.toolkit.readfile('pdb',output_path+'/DOCK_LOG/'+r_pdbid+'/'+r_pdbid+'_ligand_1.pdb'))
					IFP = fingerprints.InteractionFingerprint(bind_ligand, bind_protein)
					bindscore = tanimoto(crystal_IFP,IFP)

					ligtmapscore = 0.7*float(r_ligandscore)+0.3*bindscore
  
 
					# Write the binding fingerprint score to the IFP_result.csv (sort)
					IFP_filename = summary_path+'/IFP_result.csv'

					if os.path.isfile(IFP_filename) != True:
						IFP_file = open(IFP_filename,'w')
						IFP_file.write("PDB;Class;TargetName;LigandName;LigandSimilarityScore;BindingSimilarityScore;LigTMapScore;PredictedAffinity;DockingScore\n")
						IFP_file.write(r_pdbid+';'+select_db+';'+id_proname[r_pdbid]+';'+id_ligname[r_pdbid]+';'+str(round(float(r_ligandscore),6))+';'+str(round(bindscore,6))+';'+str(round(ligtmapscore,6))+';'+str(round(predictions[0],6))+';'+str(round(float(r_dockscore),6))+'\n')
						IFP_file.close()
					else:
						tmp_list = []
						tmp_file = open(IFP_filename,'r')
						label_line = tmp_file.readline()
						for line in tmp_file.readlines():
							line = line.rstrip()
							line = line.split(';')
							tmptmp_list = []
							for i in range(9):
								tmptmp_list.append(line[i])
							tmp_list.append(tmptmp_list)
						tmp_file.close()
						tmp_list.append([r_pdbid,select_db,id_proname[r_pdbid],id_ligname[r_pdbid],str(round(float(r_ligandscore),6)),str(round(bindscore,6)),str(round(ligtmapscore,6)),str(round(predictions[0],6)),str(round(float(r_dockscore),6))])
						tmp_list.sort(key=lambda s: s[6],reverse=True)
						IFP_file = open(IFP_filename,'w')
						IFP_file.write(label_line)
						for ele in tmp_list:
							write_str = ';'.join(ele)
							IFP_file.write(write_str+'\n')
						IFP_file.close()
				else:
					IFP_file = open(IFP_filename,'a')
					IFP_file.write(r_pdbid+';'+select_db+';'+id_proname[r_pdbid]+';'+id_ligname[r_pdbid]+';'+r_ligandscore+';n.a.;'+r_ligandscore+';'+str(predictions[0])+';'+r_dockscore+'\n')
					IFP_file.close()
		#### PART4 COMPLETED ####

		status_f.write(input_num+':'+select_db+':Complete\n')
		status_f.close()
		print "DONE"
	except Exception as e:
 		status_f = open('Output/status.txt','a')
		status_f.write(input_num+':'+select_db+':Fail\n')
		status_f.close()
		print "NO RESULT"

main_part()
