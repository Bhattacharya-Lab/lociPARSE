"""
Author: Sumit Tarafder
"""

import torch
import argparse
from torch.utils.data import DataLoader
from torch import nn
import os
import numpy as np
from tqdm import tqdm

from dataloader import RNADataset
import model as ipa

parser = argparse.ArgumentParser(description='RNA QA')
parser.add_argument('--gpunum', type=int, default=0, metavar='N',
                help='GPU node for inference')
parser.add_argument('--modelpath', type=str, default="Model/QAmodel_paper.pt", metavar='N',
                help='GPU node for inference')
                

args = parser.parse_args()

def check_cuda_compatibility():
        
    if not torch.cuda.is_available():
        return 0
    
    torch_compiled_cuda_version = torch.version.cuda
    
    system_cuda_version = torch.cuda.get_arch_list()

    if torch_compiled_cuda_version is None:
        return 0

    major_compiled_version = int(torch_compiled_cuda_version.split('.')[0])
    major_system_versions = [int(v.split('_')[1]) for v in system_cuda_version]

    if not major_compiled_version in major_system_versions:
        # Print versions for debugging
        # print(f"PyTorch compiled with CUDA version: {torch_compiled_cuda_version}")
        # print(f"System CUDA capabilities: {system_cuda_version}")
    
        print("The installed PyTorch version is not compatible with the current CUDA version on the machine. Running on CPU")
        return 0
    
    return 1

check_cuda = check_cuda_compatibility()
device_str="cuda:"+str(args.gpunum)
device = torch.device(device_str if check_cuda else "cpu")

print(f"Device = {device}")

def main():

    torch.manual_seed(1)
    torch.cuda.random.manual_seed(1)
    
    model = ipa.RNAQA(in_node_nf= 6, hidden_nf=128, out_node_nf=1, in_edge_nf=149, droprate= 0.1, device=device, n_layers=4)
    
    #===================================================
    restore_path = args.modelpath

    model.load_state_dict(torch.load(restore_path, map_location=device))
    print(f"Model loaded for inference = {restore_path}")

    model = model.to(device)

    dataset_test = RNADataset(partition='input')
    loader_test = DataLoader(dataset_test, batch_size=1, shuffle=False, drop_last=False)
    
    test(model,loader_test)
    
def test(model, loader_test):

    model.eval()

    for batch_idx, (data, RNAid) in tqdm(enumerate(loader_test)):
        
        idstr=""
        
        for item in RNAid:
            idstr += item
            break
        
        V, E, X, G = data

        V = V.to(torch.float32).to(device)
        X = X.to(torch.float32).to(device)
        
        E = E.to(torch.float32).to(device)

        G = G.to(device)
        
        with torch.no_grad():
            Pred_V= model(V, E, G, X)
            
            pNuL = nn.Sigmoid()(Pred_V.squeeze())
            pNuL = pNuL.detach().cpu().numpy()

            pMoL = np.average(pNuL)

        pred_path = f"Prediction/{idstr}"
            
        if not os.path.exists(pred_path):
            os.makedirs(pred_path)
        
        with open(pred_path + "/score.txt", "w") as wfile:
            wfile.write("{:.2f}".format(pMoL) + "\n")
            
            for i in range(pNuL.shape[0]):
                wfile.write(str(i+1) + "\t{:.2f}".format(pNuL[i]) + "\n")

if __name__ == "__main__":
    main()

