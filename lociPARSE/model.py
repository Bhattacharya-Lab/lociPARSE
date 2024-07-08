"""
Credit: PIPPack repository (https://github.com/Kuhlman-Lab/PIPPack)
"""

import math, random
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from .utils import get_rigid_from_three_points

def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.random.manual_seed(seed)
    torch.manual_seed(seed)
    torch.cuda.random.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def get_bb_frames(P: torch.Tensor, C4p: torch.Tensor, N: torch.Tensor):
    return get_rigid_from_three_points(P, C4p, N)
 
def gather_nodes(nodes, neighbor_idx):
    # Features [...,N,C] at Neighbor indices [...,N,K] => [...,N,K,C]
    is_batched = neighbor_idx.dim() == 3
    n_feat_dims = nodes.dim() - (1 + is_batched)

    # Flatten and expand indices per batch [...,N,K] => [...,NK] => [...,NK,C]
    neighbors_flat = neighbor_idx.view((*neighbor_idx.shape[:-2], -1))
    for _ in range(n_feat_dims):
        neighbors_flat = neighbors_flat.unsqueeze(-1)
    neighbors_flat = neighbors_flat.expand(*([-1] * (1 + is_batched)), *nodes.shape[-n_feat_dims:])
    
    # Gather and re-pack
    neighbor_features = torch.gather(nodes, -n_feat_dims - 1, neighbors_flat)
    neighbor_features = neighbor_features.view(list(neighbor_idx.shape) + list(nodes.shape[-n_feat_dims:]))
    return neighbor_features


def cat_neighbors_nodes(h_nodes, h_neighbors, E_idx):
    h_nodes = gather_nodes(h_nodes, E_idx)
    h_nn = torch.cat([h_neighbors, h_nodes], -1)
    return h_nn


def get_act_fxn(act: str):
    if act == 'relu':
        return F.relu
    elif act == 'gelu':
        return F.gelu
    elif act == 'elu':
        return F.elu
    elif act == 'selu':
        return F.selu
    elif act == 'celu':
        return F.celu
    elif act == 'leaky_relu':
        return F.leaky_relu
    elif act == 'prelu':
        return F.prelu
    elif act == 'silu':
        return F.silu
    elif act == 'sigmoid':
        return nn.Sigmoid()


class MLP(nn.Module):
    def __init__(self, num_in, num_inter, num_out, num_layers, act='relu', bias=True):
        super().__init__()
        
        # Linear layers for MLP
        self.W_in = nn.Linear(num_in, num_inter, bias=bias)
        self.W_inter = nn.ModuleList([nn.Linear(num_inter, num_inter, bias=bias) for _ in range(num_layers - 2)])
        self.W_out = nn.Linear(num_inter, num_out, bias=bias)
        
        # Activation function
        self.act = get_act_fxn(act)
        
    def forward(self, X):
        
        # Embed inputs with input layer
        X = self.act(self.W_in(X))
        
        # Pass through intermediate layers
        for layer in self.W_inter:
            X = self.act(layer(X))
            
        # Get output from output layer
        X = self.W_out(X)
        
        return X


class IPA(nn.Module):
    def __init__(self, node_dim, edge_dim, hidden_dim=16, n_heads=4, n_query_points=8, n_value_points=4, position_scale=10.0, dropout=0.1, act='relu'):
        super().__init__()
        
        self.hidden_dim = hidden_dim
        self.n_heads = n_heads
        self.n_query_points = n_query_points
        self.n_value_points = n_value_points
        self.position_scale = position_scale
        
        # Linear layers for queries, keys, and values
        self.linear_q = nn.Linear(node_dim, hidden_dim * n_heads)
        self.linear_kv = nn.Linear(node_dim, 2 * hidden_dim * n_heads)
        
        self.linear_q_points = nn.Linear(node_dim, n_heads * n_query_points * 3)
        self.linear_kv_points = nn.Linear(node_dim, n_heads * (n_query_points + n_value_points) * 3)
        
        self.linear_b = nn.Linear(edge_dim, n_heads)
        
        self.head_weights = nn.Parameter(torch.zeros((n_heads)))
        with torch.no_grad():
            self.head_weights.fill_(0.541324854612918)
        
        out_dim = n_heads * (edge_dim + hidden_dim + n_value_points * 4)
        self.linear_out = nn.Linear(out_dim, node_dim)
        
        self.dropout = nn.ModuleList([nn.Dropout(dropout) for _ in range(2)])
        self.norm = nn.ModuleList([nn.LayerNorm(hidden_dim) for _ in range(2)])
        
        self.node_dense = MLP(node_dim, node_dim * 4, node_dim, num_layers=2, act=act)
        
    def _get_node_update(self, h_V, h_E, E_idx, bb_to_global, mask_attend=None):

        # Generate queries, keys, and values from nodes
        q = self.linear_q(h_V) # [*, N_res, H * C]
        q = q.view(q.shape[:-1] + (self.n_heads, -1)) # [*, N_res, H, C]
        
        kv = self.linear_kv(h_V) # [*, N_res, 2 * H * C]
        kv = kv.view(kv.shape[:-1] + (self.n_heads, -1)) # [*, N_res, H, 2 * C]
        k, v = torch.split(kv, self.hidden_dim, dim=-1) # 2 [*, N_res, H, C]
        
        # Generate query, key, and value points from nodes
        q_pts = self.linear_q_points(h_V) # [*, N_res, H * P_q * 3]
        q_pts = torch.split(q_pts, q_pts.shape[-1] // 3, dim=-1) # 3 [*, N_res, H * P_q]
        q_pts = torch.stack(q_pts, dim=-1) # [*, N_res, H * P_q, 3]
        q_pts = bb_to_global[..., None].apply(q_pts) # [*, N_res, H * P_q, 3]
        q_pts = q_pts.view(
            q_pts.shape[:-2] + (self.n_heads, self.n_query_points, 3)
        ) # [*, N_res, H, P_q, 3]
        
        kv_pts = self.linear_kv_points(h_V) # [*, N_res, H * (P_q + P_v) * 3]
        kv_pts = torch.split(kv_pts, kv_pts.shape[-1] // 3, dim=-1) # 3 [*, N_res, H * (P_q + P_v)]
        kv_pts = torch.stack(kv_pts, dim=-1) # [*, N_res, H * (P_q + P_v), 3]
        kv_pts = bb_to_global[..., None].apply(kv_pts) # [*, N_res, H * (P_q + P_v), 3]
        kv_pts = kv_pts.view(kv_pts.shape[:-2] + (self.n_heads, -1, 3))
        k_pts, v_pts = torch.split(
            kv_pts, [self.n_query_points, self.n_value_points], dim=-2
        )# [*, N_res, H, P_q, 3], [*, N_res, H, P_v, 3]
        
        # Compute attention bias
        b = self.linear_b(h_E) # [*, N_res, K, H]
        
        # Compute attention weight
        a = torch.einsum("...ihc,...ijhc->...ijh", q, gather_nodes(k, E_idx))
        a *= math.sqrt(1.0 / (3 * self.hidden_dim))
        a += math.sqrt(1.0 / 3) * b # [*, N_res, K, H]
        
        pt_att = q_pts.unsqueeze(-4) - gather_nodes(k_pts, E_idx) # [*, N_res, K, H, P_q, 3]
        pt_att = torch.sum(pt_att ** 2, dim=-1) # [*, N_res, K, H, P_q]
        
        head_weights = F.softplus(self.head_weights).view(
            *((1,) * len(pt_att.shape[:-2]) + (-1, 1))
        ) # [*, 1, 1, H, 1]
        pt_att = math.sqrt(1.0 / (3 * (self.n_query_points * 9.0 / 2))) * head_weights * pt_att # [*, N_res, K, H, P_q]
        pt_att = torch.sum(pt_att, dim=-1) * -0.5 # [*, N_res, K, H]
        
        if mask_attend is not None:
            att_mask = 1e5 * (mask_attend - 1)
        else:
            att_mask = torch.zeros_like(E_idx)
        
        a = a + pt_att + att_mask[..., None] # [*, N_res, K, H]
        a = F.softmax(a, dim=-2) # [*, N_res, K, H]
        
        # Compute update
        # [*, N_res, H, C_hidden]
        o = torch.einsum('...ijh,...ijhc->...ihc', a, gather_nodes(v, E_idx))
        o = o.view(*o.shape[:-2], -1)

        o_pt = torch.einsum("...ijh,...ijhpx->...ihpx", a, gather_nodes(v_pts, E_idx))
        o_pt = bb_to_global[..., None, None].invert_apply(o_pt) # [*, N_res, H, P_v, 3]
        o_pt_norm = torch.sqrt(torch.sum(o_pt ** 2, dim=-1) + 1e-8).view(*o_pt.shape[:-3], -1)
        o_pt = o_pt.reshape(*o_pt.shape[:-3], -1, 3)
        
        o_pair = torch.einsum("...ijh,...ijc->...ihc", a, h_E) # [*, N_res, H, C_z]
        o_pair = o_pair.view(*o_pair.shape[:-2], -1)
        
        # Compute node update
        s = self.linear_out(
            torch.cat(
                (o, *torch.unbind(o_pt, dim=-1), o_pt_norm, o_pair), dim=-1
            )
        )
        
        return s
        
    def forward(self, h_V, h_E, E_idx, bb_to_global, mask_V=None, mask_attend=None):
        s = self._get_node_update(h_V, h_E, E_idx, bb_to_global, mask_attend)
        h_V = self.norm[0](h_V + self.dropout[0](s))
        node_m = self.node_dense(h_V)
        h_V = self.norm[1](h_V + self.dropout[1](node_m))
        
        if mask_V is not None:
            h_V = h_V * mask_V[..., None]
            
        return h_V


class RNAQA(nn.Module):
    def __init__(self, in_node_nf, hidden_nf, out_node_nf, in_edge_nf, droprate= 0.1, device='cpu', n_layers=4, k_neighbors = 20, act = "relu", n_points = 3, position_scale = 10.0):
        
        super().__init__()
        self.hidden_nf = hidden_nf
        self.device = device
        self.n_layers = n_layers
        self.k_neighbors = k_neighbors
        self.droprate = droprate
        self.position_scale = position_scale
        
        self.embedding_out = nn.Linear(self.hidden_nf, out_node_nf)
        
        # Normalization and embedding
        self.node_embedding = nn.Linear(in_node_nf,  self.hidden_nf, bias=True)
        self.norm_nodes = nn.LayerNorm(self.hidden_nf)
        self.edge_embedding = nn.Linear(in_edge_nf, self.hidden_nf, bias=True)
        self.norm_edges = nn.LayerNorm(self.hidden_nf)

        self.mpnn_layers = nn.ModuleList([
                IPA(hidden_nf, hidden_nf, hidden_nf, dropout=droprate, act=act)
                for _ in range(self.n_layers)
            ])            
        
        # Addition Sumit Tarafder
        # For down projecting hidden dimension from 128 to 1 
        self.embedding_out = nn.Linear(self.hidden_nf, out_node_nf)
            
        self.to(self.device)

        seed_all(1000)

    def forward(self, V, E, E_idx, X):
        
        # Embed nodes
        V = self.node_embedding(V)
        V = self.norm_nodes(V)
        
        # Embed edges
        E = self.edge_embedding(E)
        E = self.norm_edges(E)

        # Get backbone global frames from N, CA, and C
        scaled_X = X / self.position_scale
        
        bb_to_global = get_bb_frames(scaled_X[..., 0, :], scaled_X[..., 1, :], scaled_X[..., 2, :])
        
        for layer in self.mpnn_layers:
            V = layer(V, E, E_idx, bb_to_global)
            
        pred_node_emb = self.embedding_out(V)
        
        return pred_node_emb
        
