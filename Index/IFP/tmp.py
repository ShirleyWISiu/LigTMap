ifp_list = []
with open('IFP_error.csv','r') as ifp_file:
	for line in ifp_file:
		line = line.rstrip()
		line = line.split(';')
		pdbid = line[0]
		ifp_list.append(pdbid)

with open('../fing/MACCSF_error.csv','r') as pham_file:
	for line in pham_file:
		line = line.rstrip().split(';')
		pdbid = line[0]
		if pdbid not in ifp_list:
			print pdbid	
