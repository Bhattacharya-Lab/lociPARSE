from .feature_generation import genFeat
from . import model as ipa
from .utils import *

import torch
from torch import nn
import os, sys
import numpy as np

class pMoL:
    def __init__(self, pMoL):
        self.value = float("{:.2f}".format(pMoL))

    def show(self):
        print(self.value)

class pNuL:
    def __init__(self, values):
        
        pNul = [float("{:.2f}".format(value)) for value in values]
        
        self.values = pNul

    def show(self, index=None):
        
        if not index is None: 
            if index <= 0 or index > len(self.values):
                sys.exit("Index out of range")
            else:
                print(self.values[index-1])
        else:
            print(self.values)

class ScoreResult:
    def __init__(self, pMoL_value, pNuL_values):
        self.pMoL = pMoL(pMoL_value)
        self.pNuL = pNuL(pNuL_values)

    def save(self, filename):
        with open(filename, "w") as wfile:
            wfile.write(str(self.pMoL.value) + "\n")
            
            for i in range(len(self.pNuL.values)):
                wfile.write(str(i+1) + str(f"\t{self.pNuL.values[i]}\n"))

class lociparse:
    def __init__(self, model_path=None):
        self.check_cuda = self.check_cuda_compatibility()
        self.device = self.get_device()
        self.model = self.load_model(model_path)

    def check_cuda_compatibility(self):
        
        if not torch.cuda.is_available():
            return 0
        
        # Get the CUDA version that PyTorch was compiled with
        torch_compiled_cuda_version = torch.version.cuda
        
        # Get the CUDA version available on the system
        system_cuda_version = torch.cuda.get_arch_list()

        if torch_compiled_cuda_version is None:
            return 0

        # Check if the system CUDA version supports the CUDA version PyTorch was compiled with
        major_compiled_version = int(torch_compiled_cuda_version.split('.')[0])
        major_system_versions = [int(v.split('_')[1]) for v in system_cuda_version]

        if not major_compiled_version in major_system_versions:
            # Print versions for debugging
            # print(f"PyTorch compiled with CUDA version: {torch_compiled_cuda_version}")
            # print(f"System CUDA capabilities: {system_cuda_version}")
        
            # print("The installed PyTorch version is not compatible with the current CUDA version on the machine. Running on CPU")
            return 0
       
        return 1
    
    def get_device(self):
        
        device_str="cuda:0"
        device = torch.device(device_str if self.check_cuda else "cpu")
        return device

    def load_model(self, model_path=None):
        if model_path is None:
            model_path = os.path.join(os.path.dirname(__file__), 'QAmodel_lociPARSE.pt')
         
        model = ipa.RNAQA(in_node_nf= 6, hidden_nf=128, out_node_nf=1, in_edge_nf=149, droprate= 0.1, device=self.device, n_layers=4)

        model.load_state_dict(torch.load(model_path, map_location=self.device))
        model = model.to(self.device)
        model.eval()
        
        return model

    def extract_features(self, pdb_file_path):
        
        V, E, G, X = genFeat(pdb_file_path)
        
        return V, E, G, X

    def score(self, pdb_file_path):
        
        V, E, G, X = self.extract_features(pdb_file_path)

        V = V.unsqueeze(0)
        E = E.unsqueeze(0)
        G = G.unsqueeze(0)
        X = X.unsqueeze(0)
        
        V = V.to(torch.float32).to(self.device)
        X = X.to(torch.float32).to(self.device)
        
        E = E.to(torch.float32).to(self.device)

        G = G.to(self.device)
        
        with torch.no_grad():
            Pred_V= self.model(V, E, G, X)
            
            pNuL = nn.Sigmoid()(Pred_V.squeeze())
            pNuL = pNuL.detach().cpu().numpy()

            pMoL = np.average(pNuL)

        return ScoreResult(pMoL, pNuL)

        