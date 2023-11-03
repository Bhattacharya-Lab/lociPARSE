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

args = parser.parse_args()

cuda_on = torch.cuda.is_available()
device_str="cuda:"+str(args.gpunum)

device = torch.device(device_str if cuda_on else "cpu")

print(device)

def main():

    torch.manual_seed(1)
    torch.cuda.random.manual_seed(1)
    
    model = ipa.RNAQA(in_node_nf= 6, hidden_nf=128, out_node_nf=1, in_edge_nf=149, droprate= 0.1, device=device, n_layers=4)
    
    #===================================================
    restore_path = "Model/QAmodel.pt"

    model.load_state_dict(torch.load(restore_path))
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

