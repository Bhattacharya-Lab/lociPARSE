import torch
import argparse
from torch.utils.data import DataLoader
from dataloader_train import RNADataset

import model as ipa
import os, random
from torch import nn, optim

import numpy as np
from tqdm import tqdm

parser = argparse.ArgumentParser(description='RNA QA')
parser.add_argument('--epochs', type=int, default=200, metavar='N',
                    help='number of epochs to train (default: 10)')
parser.add_argument('--lr', type=float, default=5e-4, metavar='N',
                    help='learning rate')
parser.add_argument('--nf', type=int, default=64, metavar='N',
                    help='hidden feature dimension')
parser.add_argument('--n_layers', type=int, default=4, metavar='N',
                    help='number of layers for the EGCL')
parser.add_argument('--weight_decay', type=float, default=1e-16, metavar='N',
                    help='timing experiment')
parser.add_argument('--droprate', type=float, default= 0.1, metavar='N',
                    help='Dropout rate for dropout of neurons')
parser.add_argument('--gpunum', type=int, default=0, metavar='N',
                    help='node index of GPU  used for pytorch DDP training')
parser.add_argument('--numneighbour', type=int, default=20, metavar='N',
                help='Num neigh in graph construction')
parser.add_argument('--SEED', type=int, default=0, metavar='N',
                help='Seed for training reproducibility')
                        
args = parser.parse_args()

device_str="cuda:"+str(args.gpunum)
device = torch.device(device_str if torch.cuda.is_available() else "cpu")

print(device)

save_dir="Model"

try:
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
except OSError:
    pass

save_path_valid = os.path.join(save_dir, 'QAmodel_retrained.pt')

print(f"Saving model in - {save_path_valid}")

l1loss = nn.L1Loss()

def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.random.manual_seed(seed)
    torch.manual_seed(seed)
    torch.cuda.random.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False  
    
def main():

    SEED = args.SEED
    
    seed_all(SEED)

    dataset_train = RNADataset(partition='train')
    loader_train = DataLoader(dataset_train, batch_size=1, drop_last=False)
    
    print(len(loader_train))

    dataset_val = RNADataset(partition='valid')
    loader_val = DataLoader(dataset_val, batch_size=1, drop_last=False)

    print(len(loader_val))
    
    #IPA
    
    model = ipa.RNAQA(in_node_nf= 6, hidden_nf=args.nf, out_node_nf=1, in_edge_nf=149, droprate= args.droprate, device=device, n_layers=args.n_layers)
    
    #===================================================
    
    model = model.to(device)
    #print(model)

    optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    best_val_loss = 1e22
    best_val_epoch = 0
    
    for epoch in tqdm(range(0, args.epochs)):
    
        train(model, optimizer, epoch, loader_train)
        
        val_loss = train(model, optimizer, epoch, loader_val, backprop=False)
            
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_val_epoch = epoch

            save_path_valid = f"{save_dir}/QAmodel_retrained.pt"

            torch.save(model.state_dict(), save_path_valid)
            print("Saving model*** Best Val Loss: %.5f \t Best epoch %d" % (best_val_loss, best_val_epoch))
            print(save_path_valid)
            
    print(f"Best validation loss = {best_val_loss}")
    print(f"Best validation epoch =  {best_val_epoch}")


def train(model, optimizer, epoch, loader_train, backprop=True):
    
    if backprop:
        model.train()
    else:
        model.eval()
    
    res = {'epoch': epoch, 'loss': 0, 'counter': 0}

    for batch_idx, data in enumerate(loader_train):

        optimizer.zero_grad()

        batch_size = 1

        V, E, X, G, labels = data

        V = V.to(torch.float32).to(device)
        X = X.to(torch.float32).to(device)
        
        labels = labels.squeeze()
        labels = labels.to(device)
        
        E = E.to(torch.float32).to(device)
        G = G.to(device)
        
        Pred_V = model(V, E, G, X)
        
        per_node_plDDT = Pred_V.squeeze()

        loss = l1loss(nn.Sigmoid()(per_node_plDDT),labels) #l1_loss

        if backprop:
            loss.backward()
            optimizer.step()
            
        res['loss'] += loss.item() * batch_size
        res['counter'] += batch_size
        
        
    lossval = float("{:.3f}".format(res['loss'] / res['counter']))

    if backprop:
        print(f"Epoch = {epoch}, Loss = {lossval}")

    return lossval

if __name__ == "__main__":
    
    main()

