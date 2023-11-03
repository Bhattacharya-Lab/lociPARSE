"""
Author: Sumit Tarafder
"""

import torch
import numpy as np

class RNADataset():

	def __init__(self,partition):
		self.data = []
		self.RNAids = []
		
		#print(partition)
		
		idlist = "Input/" + partition + ".txt"
		feature_dir = "Feature/"
		
		idlistall = []

		with open(idlist,'r', encoding = 'UTF-8') as idfile:
		
			while (RNAid := idfile.readline()):
			
				RNAid = RNAid.rstrip()

				RNAid = RNAid.split(".")[0]

				idlistall.append(RNAid)
			
		for i in range(len(idlistall)):
			
			RNAid = idlistall[i]
				
			node_feature_path = feature_dir + RNAid + "/nodefeat.pt"
			edge_feature_path = feature_dir + RNAid + "/edgefeat.pt"
			coordpath = feature_dir + RNAid + "/Coord.npy"
			graph_path = feature_dir + RNAid + "/adj_mat.pt"
			
			node_feats = torch.load(node_feature_path)
			L = node_feats.size()[0]
			
			Xs = np.load(coordpath)
			Y = int(Xs.shape[0] / L)
			X = np.reshape(Xs, (L,Y,3)) #(L,Y,3)
			edge_feats = torch.load(edge_feature_path)
			
			graph = torch.load(graph_path)
					
			self.data.append((node_feats, edge_feats, X, graph))
			self.RNAids.append(RNAid)

				
	# support indexing such that dataset[i] can be used to get i-th sample
	def __getitem__(self, index):
		#return self.RNAgraphs[index], self.Edgelabels[index]
		return self.data[index], self.RNAids[index]

	# we can call len(dataset) to return the size
	def __len__(self):
		return len(self.data)
