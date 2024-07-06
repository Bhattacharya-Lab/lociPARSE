"""
Author: Sumit Tarafder
"""

import torch
import numpy as np
from tqdm import tqdm

class RNADataset():

	def __init__(self,partition):
		self.data = []
		
		idlist = "Input/Dataset/Train" + partition + ".txt"
		feature_dir = "Feature/"

		idlistall = []

		with open(idlist,'r', encoding = 'UTF-8') as idfile:
		
			while (RNAid := idfile.readline()):
			
				RNAid = RNAid.rstrip()

				idlistall.append(RNAid)

		for i in tqdm(range(len(idlistall))):
			
			RNAid = idlistall[i]
				
			node_feature_path = feature_dir + RNAid + "/nodefeat.pt"
			label_path = feature_dir + RNAid + "/labels.npy"
			edge_feature_path = feature_dir + RNAid + "/edgefeat.pt"
			coordpath = feature_dir + RNAid + "/Coord.npy"
			graph_path = feature_dir + RNAid + "/adj_mat.pt"
			
			graph = torch.load(graph_path)
				
			node_feats = torch.load(node_feature_path)
			L = node_feats.size()[0]

			Xs = np.load(coordpath)
			Y = int(Xs.shape[0] / L)
			X = np.reshape(Xs, (L,Y,3)) #(L,Y,3)
			
			labels = np.load(label_path)
			
			if len(labels.shape) != 2:
				labels = labels[:, np.newaxis]
				
			labels = torch.as_tensor(labels)
			
			edge_feats = torch.load(edge_feature_path)
			
			self.data.append((node_feats, edge_feats, X, graph, labels))
				
	# support indexing such that dataset[i] can be used to get i-th sample
	def __getitem__(self, index):
		return self.data[index]

	# we can call len(dataset) to return the size
	def __len__(self):
		return len(self.data)
