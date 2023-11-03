"""
Author: Sumit Tarafder
"""

import argparse
import os
import time    
import numpy as np
import shutil
import torch
from tqdm import tqdm

purines = ['A' , 'G']
pyrimidines = ['U', 'C']

def createAndCheckCG(args, idname, seqlen, fullmodelpath,coordfilename, atompernt):

	outputdir = args.Path
	
	decoyfilename = fullmodelpath
	outfilename = outputdir + idname + f"/CG{atompernt}.txt"
		
	infoarr = np.zeros((seqlen, atompernt) , dtype=np.int32)

	savedlineslist = [None]*(atompernt+3) # 3 extra for 3 backups

	ntcur = 0
	count = 0
	terflag = 1

	with open(decoyfilename, 'r', encoding='UTF-8') as idfile,\
		open(outfilename, 'w', encoding='UTF-8') as outfile:
		
		while (linef := idfile.readline().rstrip()):
			
			linef = linef
			tokenlist = linef.split()
			
			if tokenlist[0] == "TER":

				if idfile.readline():
					break
				missinglist = []
				for i in range(len(savedlineslist)-3):
					if savedlineslist[i] is None:
						missinglist.append(i)
				
				for j in range(len(missinglist)):
							
					i = missinglist[j]

					if i == 0:
						savedlineslist[i] =  savedlineslist[atompernt]
						if savedlineslist[i] is None:
							savedlineslist[i] = "0.000 0.000 0.000"
						infoarr[nt-1][i] = 1
					elif i == 1:
						savedlineslist[i] =  savedlineslist[atompernt+1]
						if savedlineslist[i] is None:
							savedlineslist[i] = "0.000 0.000 0.000"
						infoarr[nt-1][i] = 1
					elif i == 2:
						savedlineslist[i] =  savedlineslist[atompernt+2]
						if savedlineslist[i] is None:
							savedlineslist[i] = "0.000 0.000 0.000"
						infoarr[nt-1][i] = 1

				for i in range(len(savedlineslist)-3):
					if savedlineslist[i] is not None:
						outfile.write(savedlineslist[i]+"\n")
					
				terflag = 0
				
			if tokenlist[0] == "ATOM":
				newtokenlist = []
				
				for t in range(len(tokenlist)):
					if tokenlist[t] != '' and tokenlist[t] != '\t':
						newtokenlist.append(tokenlist[t])

				atom = newtokenlist[2]

				ntname = newtokenlist[3]
				
				if ntname in purines:
					nttype = "Purine"
				else:
					nttype = "Pyrimidine"

				chain = linef[21:22]
				
				if chain != ' ':
					nt = newtokenlist[5]
					nt = int("".join(nt.split()))
				else:
					nt = newtokenlist[4]
					nt = int("".join(nt.split()))

				if count == 0 or nt == ntcur:
					
					x = linef[29:38].rstrip()
					x = "".join(x.split())
					y = linef[38:46].rstrip()
					y = "".join(y.split())
					n = len(str(y).split(".")[1])
					
					if n == 2:
						y = linef[38:47].rstrip()
						y = "".join(y.split())
						z = linef[47:55].rstrip()
						z = "".join(z.split())
					else:
						z = linef[46:55].rstrip()
						z = "".join(z.split())

					newline = x + " " + y + " " + z
					
					if atom == "P" and atompernt > 0:
						infoarr[nt-1][0] = 1
						savedlineslist[0] = newline
					if atom == "C4'" and atompernt > 1:
						infoarr[nt-1][1] = 1
						savedlineslist[1] = newline
					if  atompernt > 2:
						if (atom == "N9" and nttype == "Purine") or (atom == "N1" and nttype == "Pyrimidine"):
							infoarr[nt-1][2] = 1
							savedlineslist[2] = newline
					if atom == "C5'" and atompernt > 3:
						infoarr[nt-1][3] = 1
						savedlineslist[3] = newline
					if atom == "C3'" and atompernt > 4:
						infoarr[nt-1][4] = 1
						savedlineslist[4] = newline
					if atom == "C1'" and atompernt > 5:
						infoarr[nt-1][5] = 1
						savedlineslist[5] = newline
					if atom == "O3'" and atompernt > 6:
						infoarr[nt-1][6] = 1
						savedlineslist[6] = newline
					
					if atom == "C5'":
						savedlineslist[atompernt] = newline
					if atom == "C3'":
						savedlineslist[atompernt+1] = newline
					if atom == "C1'":
						savedlineslist[atompernt+2] = newline
						
				if nt != ntcur:
					ntcur = nt
					
					if count > 0:
						missinglist = []
						for i in range(len(savedlineslist)-3):
							if savedlineslist[i] is None:
								missinglist.append(i)
						for j in range(len(missinglist)):
							
							i = missinglist[j]

							if i == 0:
								savedlineslist[i] =  savedlineslist[atompernt]
								
								if savedlineslist[i] is None:
									savedlineslist[i] = "0.000 0.000 0.000"
								infoarr[nt-2][i] = 1
							elif i == 1:
								savedlineslist[i] =  savedlineslist[atompernt+1]
								if savedlineslist[i] is None:
									savedlineslist[i] = "0.000 0.000 0.000"

								infoarr[nt-2][i] = 1
							elif i == 2:
								savedlineslist[i] =  savedlineslist[atompernt+2]

								if savedlineslist[i] is None:
									savedlineslist[i] = "0.000 0.000 0.000"

								infoarr[nt-2][i] = 1
								
						
						for i in range(len(savedlineslist)-3):
							if savedlineslist[i] is not None:
								outfile.write(savedlineslist[i]+"\n")

						savedlineslist = [None]*(atompernt+3) # 3 extra for 3 backups
							
					x = linef[29:38].rstrip()
					x = "".join(x.split())
					y = linef[38:46].rstrip()
					y = "".join(y.split())
					n = len(str(y).split(".")[1])
					
					if n == 2:
						y = linef[38:47].rstrip()
						y = "".join(y.split())
						z = linef[47:55].rstrip()
						z = "".join(z.split())
					else:
						z = linef[46:55].rstrip()
						z = "".join(z.split())

					newline = x + " " + y + " " + z

					if atom == "P" and atompernt > 0:
						infoarr[nt-1][0] = 1
						savedlineslist[0] = newline
					if atom == "C4'" and atompernt > 1:
						infoarr[nt-1][1] = 1
						savedlineslist[1] = newline
					if  atompernt > 2:
						if (atom == "N9" and nttype == "Purine") or (atom == "N1" and nttype == "Pyrimidine"):
							infoarr[nt-1][2] = 1
							savedlineslist[2] = newline
					if atom == "C5'" and atompernt > 3:
						infoarr[nt-1][3] = 1
						savedlineslist[3] = newline
					if atom == "C3'" and atompernt > 4:
						infoarr[nt-1][4] = 1
						savedlineslist[4] = newline
					if atom == "C1'" and atompernt > 5:
						infoarr[nt-1][5] = 1
						savedlineslist[5] = newline
					if atom == "O3'" and atompernt > 6:
						infoarr[nt-1][6] = 1
						savedlineslist[6] = newline
					
					if atom == "O5'":
						savedlineslist[atompernt] = newline
					if atom == "C3'":
						savedlineslist[atompernt+1] = newline
					if atom == "C1'":
						savedlineslist[atompernt+2] = newline

				count += 1
				#break

		if terflag == 1:
			
			missinglist = []
			for i in range(len(savedlineslist)-3):
				if savedlineslist[i] is None:
					missinglist.append(i)
			
			for j in range(len(missinglist)):
				
				i = missinglist[j]

				if i == 0:
					savedlineslist[i] =  savedlineslist[atompernt]
					if savedlineslist[i] is None:
						savedlineslist[i] = "0.000 0.000 0.000"
					infoarr[nt-1][i] = 1
				elif i == 1:
					savedlineslist[i] =  savedlineslist[atompernt+1]
					if savedlineslist[i] is None:
						savedlineslist[i] = "0.000 0.000 0.000"
					infoarr[nt-1][i] = 1
				elif i == 2:
					savedlineslist[i] =  savedlineslist[atompernt+2]
					if savedlineslist[i] is None:
						savedlineslist[i] = "0.000 0.000 0.000"
					infoarr[nt-1][i] = 1

			for i in range(len(savedlineslist)-3):
				if savedlineslist[i] is not None:
					outfile.write(savedlineslist[i]+"\n")
			
	count = np.count_nonzero(infoarr)
	
	if count == infoarr.shape[0] * infoarr.shape[1]:
		coordtxt = np.loadtxt(outfilename)
		np.save(coordfilename, coordtxt)

		return 1
	else:
		return 0

def renumber(path, idname, pdbfilename):
	
	newpdbfilename = path + idname + f"/temprenum.pdb"
		
	ntnumnew = 0
	ntprev = "-100"
	counter = 1
	terflag = 1

	with open(pdbfilename, 'r', encoding='UTF-8') as iddfile,open(newpdbfilename, 'w', encoding='UTF-8') as outfile:
		while (linef := iddfile.readline()):
			#linef = linef.rstrip()
			
			tokenlist = linef.split(" ")

			if tokenlist[0] == "ATOM":
			
				newtokenlist = linef.split()
				
				if newtokenlist[0] != "ATOM":
					outfile.write(linef)
				else:
					linef = linef.rstrip()
					atom = newtokenlist[2]
				
					nt = newtokenlist[3]
					chain = linef[21:22]
					
					if chain == " ":
						chain = 'A'
					
					ntnum = linef[23:28]
					
					if ntnum != ntprev:
						ntprev = ntnum
						ntnumnew += 1

					x = linef[29:38].rstrip()
					x = "".join(x.split())
					
					y = linef[38:46].rstrip()
					y = "".join(y.split())
					
					n = len(str(y).split(".")[1])
					
					if n == 2:
						y = linef[38:47].rstrip()
						y = "".join(y.split())
						z = linef[47:55].rstrip()
						z = "".join(z.split())
					else:
						z = linef[46:55].rstrip()
						z = "".join(z.split())

					newstring = "ATOM"
					newstring += '{:7d}'.format(counter)

					if len(atom) <= 3:
						newstring += "  "
						newstring += '{}'.format(atom)
					else:
						newstring += " "
						newstring += '{}'.format(atom)
						newstring += " "

					atomlen = len(atom)
					space = 6 - atomlen

					for i in range(space):
						newstring += " "

					newstring += '{}'.format(nt)

					newstring += '{:>2}'.format(chain)
					newstring += '{:>4}'.format(ntnumnew)

					newstring += '{:>12}'.format(x)
					newstring += '{:>8}'.format(y)
					newstring += '{:>8}'.format(z)

					rest = linef[54:]

					if rest == "":
						newstring += '{:>6}'.format("1.00")
						newstring += '{:>6}'.format("0.00")
						newstring += '{:>12}'.format(atom[0])
					else:
						newstring += rest

					outfile.write(newstring + "\n")
					counter += 1
							
							
				#break
	os.remove(pdbfilename)
	os.rename(newpdbfilename, pdbfilename)

def getFastaAndChain(pdbfilename):
	fastaseq = ""

	with open(pdbfilename, 'r', encoding='UTF-8') as idfile:
		count=0
		ntcur = 0

		while (linef := idfile.readline()):
			linef = linef.rstrip()
			tokenlist = linef.split(" ")

			newtokenlist = []
			
			for t in range(len(tokenlist)):
				if tokenlist[t] != '':
					newtokenlist.append(tokenlist[t])

			if tokenlist[0] == "TER":
				break

			if tokenlist[0] == "ATOM":
				
				count += 1

				chain = newtokenlist[4]
					
				nt = linef[22:28].rstrip()
				nt = int("".join(nt.split()))

				if nt != ntcur:
					c = linef[19:20].rstrip()
					c = "".join(c.split())
					ntcur = nt

					fastaseq += c
					
	return fastaseq,chain	
	
def getOneHot(fastaseq):
	
	seq = fastaseq.upper()

	s = []

	for i in range(len(seq)):
		if seq[i] == 'A':
			s.append(0)
		elif seq[i] == 'U':
			s.append(1)
		elif seq[i] == 'G':
			s.append(2)
		elif seq[i] == 'C':
			s.append(3)
		else:
			s.append(4)

	onehot = np.eye(5)[s]
	
	return onehot

def CreateGraph(X, K, graphfilename):
	X_C4p = X[:,1,:] #1 because C4'
	
	X = torch.from_numpy(X_C4p)

	r1 = X.unsqueeze(0)
	r2 = X.unsqueeze(1)
    
	pdist = (r1 - r2).norm(dim=2)
    
	_, Edges = torch.topk(-pdist, min(K, X_C4p.shape[0]), axis=1)

	torch.save(Edges, graphfilename)
	
	return Edges

def CreateNodeFeat(args,seq,X,idname):
	
	nodefeat_filepath = args.Path + idname + "/nodefeat.pt"

	onehot = getOneHot(seq)

	L = len(seq)
	relpos = np.array([i/L for i in range(L)])
	relpos = np.expand_dims(relpos,axis=1)
	
	node_feat = relpos #1
	node_feat = np.concatenate((relpos, onehot), axis = 1) #5
	
	node_feat = torch.as_tensor(node_feat)
	
	torch.save(node_feat, nodefeat_filepath)


def rbf(D):
	
	"""
	Credit: PIPPack repository (https://github.com/Kuhlman-Lab/PIPPack)
	"""

	num_rbf = 16
	# Distance radial basis function
	D_min, D_max, D_count = 0., 100., num_rbf
	D_mu = torch.linspace(D_min, D_max, D_count)
	D_mu = D_mu.view([1,1,1,-1])
	D_sigma = (D_max - D_min) / D_count
	D_expand = torch.unsqueeze(D, -1)
	RBF = torch.exp(-((D_expand - D_mu) / D_sigma)**2)
	
	return RBF

def gather_edges(edges, neighbor_idx):

	"""
	Credit: PIPPack repository (https://github.com/Kuhlman-Lab/PIPPack)
	"""

	# Features [B,N,N,C] at Neighbor indices [B,N,K] => Neighbor features [B,N,K,C]
	neighbors = neighbor_idx.unsqueeze(-1).expand(-1, -1, -1, edges.size(-1))
	edge_features = torch.gather(edges, 2, neighbors)
	return edge_features

def get_rbf(A, B, Edges):
	
	"""
	Credit: PIPPack repository (https://github.com/Kuhlman-Lab/PIPPack)
	"""

	D_A_B = torch.sqrt(torch.sum((A[:,:,None,:] - B[:,None,:,:])**2,-1) + 1e-6) #[B, L, L]
	
	D_A_B_neighbors = gather_edges(D_A_B[:,:,:,None], Edges)[:,:,:,0] #[B,L,K]
	
	RBF_A_B = rbf(D_A_B_neighbors)
	
	return RBF_A_B

def get_atomic_distances(X, Edges):
	Edges = torch.unsqueeze(Edges, dim=0)
	X = torch.as_tensor(X)
	X = torch.unsqueeze(X, dim=0)

	RBF_all = []

	for i in range(X.shape[-2]):
		for j in range(X.shape[-2]):
			RBF_all.append(get_rbf(X[..., i, :], X[..., j, :], Edges))

	RBF_all = torch.cat(tuple(RBF_all), dim=-1)
	
	return RBF_all.squeeze(0)

def getSeqSepRNA(Edges):
    
    Edges = Edges.numpy()

    L = Edges.shape[0]
    K = Edges.shape[1]

    F = 5 #Number of encodings

    ssarray = np.zeros([L*K , F])
    count = 0

    for i in range(L):
        for j in range(K):

            node1 = i
            node2 = Edges[i][j]
            
            seq_sep = abs(node1 - node2)
            
            if seq_sep == 0:
                onehotindex = 0
            elif seq_sep == 1:
                onehotindex = 1
            elif seq_sep <= 5:
                onehotindex = 2
            elif seq_sep <= 24:
                onehotindex = 3
            else:
                onehotindex = 4

            ssarray[count][onehotindex] = 1

            count += 1
    
    ssarray = np.reshape(ssarray, (L,K,F))
    ssarray = torch.as_tensor(ssarray)

    return ssarray

def CreateEdgeFeat(Edges, X, idname):
    edgefeat_filepath = args.Path + idname + "/edgefeat.pt"

    E = getSeqSepRNA(Edges)
    R = get_atomic_distances(X, Edges)

    E = torch.cat((E, R), -1)

    torch.save(E, edgefeat_filepath)

def genFeatID(idname, args):

	path = args.Path

	fullmodelpath = path + idname + "/Fullmodel.pdb"

	#=====================Step-1=============================
	
	if not os.path.exists(path + idname):
		os.makedirs(path + idname)

	shutil.copy(f"Input/{idname}.pdb", fullmodelpath)
	
	renumber(path, idname, fullmodelpath)

	fastaseq,_ = getFastaAndChain(fullmodelpath)
	L = len(fastaseq)

	#=====================Step-2=============================
	# Generate a Y = 3 atom CG file from it for now, P,C4', N
	coordfilename = path + idname + "/Coord.npy"
	
	Y = 3
	
	allatomflag = createAndCheckCG(args, idname, L, fullmodelpath, coordfilename, Y)

	if allatomflag == 0:
		return 0
	
	X_s = np.load(coordfilename)#(L*Y, 3)
	X = np.reshape(X_s, (L,Y,3)) #(L,Y,3)
	
	#=====================Step-3 =============================
	# Create graph based on C4' atoms, and Top K neighbors, where K = 20
	K = args.num_neighbours
	graphfilename = args.Path + idname + "/adj_mat.pt"

	Edges = CreateGraph(X, K, graphfilename)
	
	#====================Step-4================================
	# Create Node embeddings

	CreateNodeFeat(args, fastaseq, X_s, idname)

	#===================Step-5=================================
	# Create edge features
	
	CreateEdgeFeat(Edges, X, idname)

def generateFeatures(args):

	idlist = []
	
	inputfile = "Input/input.txt"

	with open(inputfile, 'r', encoding='UTF-8') as file:

		while (idname := file.readline()):
			idname = idname.rstrip()

			idname = idname.split(".")[0]

			idlist.append(idname)
			
		for i in tqdm(range(len(idlist))):

			idname = idlist[i]

			flag = genFeatID(idname, args)

			if flag == 0:
				print(f"Skipping this target {idname}, some required atoms are missing")
				continue
			
			

if __name__ == '__main__':

	path = "Feature/" # '/' must at end
	starttime = time.time()

	parser = argparse.ArgumentParser("Feature.py")
	parser.add_argument('--Path', type = str,default = path, help = "Path of folder\
						where individual feature folders are present, default='./Feature'")
	parser.add_argument('--num_neighbours', type=int, default = 20, help="Number of nearest neigbours to consider in graph construction, top K in KNN")

	args = parser.parse_args()
	
	generateFeatures(args)

	end = time.time()

	#print(f"Total time taken = {end - starttime} seconds.")