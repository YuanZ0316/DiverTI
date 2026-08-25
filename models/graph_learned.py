import torch.nn.functional as F
import os
import torch.nn as nn
import torch.sparse
import pandas as pd
import numpy as np
from models.model_utils import normalize_adj_symm, df2tensor  
from torch_geometric.nn import GCNConv


class GraphLearned(torch.nn.Module):
    def __init__(self, nlayers, isize, neighbor, gamma, adj, dis, device, omega, di_dir, expr, percent, root_node=None):
        super().__init__()
        self.adj = adj.to(device)
        d_matrix = torch.tensor(dis, dtype=torch.float32, device=device)
        d_sorted, _ = d_matrix.sort()
        c1 = d_matrix > 0
        d_cut = torch.median(d_sorted[:, neighbor])  # on the base of k to remove long distance
        c2 = d_matrix <= d_cut
        adj_mask = torch.logical_and(c1, c2)  # 1-k neighbor, no self-loop
        if percent == 0:
            self.adj_mask = adj_mask
        else:
            self.adj_mask = self.init_adj_mask_di(di_dir, adj_mask, device=device, exprs=expr, 
                                                 percent=percent, root_node=root_node)
        print('The number of useful edges is {}'.format(self.adj_mask.sum()))
        d_matrix = torch.where(self.adj_mask, d_matrix, torch.inf) / d_cut
        self.s_d = 1 / torch.exp(gamma * torch.pow(d_matrix, 2))

        self.convs = nn.ModuleList()
        self.convs = torch.nn.ModuleList(
            [GCNConv(in_channels=isize, out_channels=isize)] +
            [GCNConv(in_channels=isize, out_channels=isize) for _ in range(nlayers - 1)])
        self.input_dim = isize
        self.omega = omega

    def calculate_diversity_increment(self, expr_i, expr_j, expr_r):

        def diversity_measure(X):
            X_flat = X.flatten()
            X_flat = X_flat[X_flat > 0]  
            if len(X_flat) == 0:
                return 0
            N = np.sum(X_flat)
            ans = np.sum([n * np.log(n) for n in X_flat])
            Dx = N * np.log(N) - ans
            return Dx


        D_ir = diversity_measure(np.concatenate([expr_i, expr_r]))
        D_jr = diversity_measure(np.concatenate([expr_j, expr_r]))
        D_ijr = diversity_measure(np.concatenate([expr_i, expr_j, expr_r]))

        di = D_ijr - (D_ir + D_jr) / 2.0
        di_normalized = np.log1p(di) if di > 0 else 0.0
        
        return di_normalized

    def generate_DI_from_adj(self, adj_mask, gene_exprs, root_node, save_path):

        import itertools

        n_cells = adj_mask.shape[0]
        cells = range(n_cells)

        root_neighbors = torch.where(adj_mask[root_node])[0].cpu().numpy()

        tri_pairs = []
        for i, j in itertools.combinations(root_neighbors, 2):
            if i != j:
                tri_pairs.append([root_node, i, j])
        
        print(f'Generating DI matrix with {len(tri_pairs)} triplets for root {root_node}')

        di_values = []
        for triple in tri_pairs:
            r, i, j = triple
            di_val = self.calculate_diversity_increment(
                gene_exprs[i], gene_exprs[j], gene_exprs[r]
            )
            di_values.append([i, j, di_val])

        di_df = pd.DataFrame(di_values, columns=['Cell1', 'Cell2', 'DI'])

        di_df.to_csv(save_path + f'DI_root_{root_node}.csv', index=False)
        return di_df

    def init_adj_mask_di(self, path, adj, device, exprs, percent=0.5, root_node=None):

        if root_node is None:
            root_node = 0  
        di_file = path + f'DI_root_{root_node}.csv'
        
        if not os.path.exists(di_file):
            if not os.path.exists(path):
                os.makedirs(path)
            print(f'Generating DI matrix for root {root_node}...')
            di_df = self.generate_DI_from_adj(
                adj_mask=adj, gene_exprs=exprs, root_node=root_node, save_path=path
            )
        else:
            di_df = pd.read_csv(di_file)
        
        group_DI = di_df.groupby(['Cell1', 'Cell2'])['DI'].mean()
        group_DI = group_DI.loc[group_DI > -np.inf]

        if len(group_DI) > 0:
            threshold = group_DI.quantile(percent)
            group_DI = group_DI.loc[group_DI > threshold]
        else:
            threshold = 0.1
            
        group_DI.fillna(np.mean(group_DI.values) if len(group_DI) > 0 else 0.1, inplace=True)
        
        if len(group_DI) > 0 and group_DI.min() < 0:
            group_DI = group_DI - group_DI.min()
        
        adj_di = df2tensor(group_DI.reset_index(), adj)
        adj_di = torch.Tensor(adj_di.toarray())
        adj_di = adj_di > 0
        adj_di = adj_di.to(device)
        
        print(f'DI-based adjacency mask has {adj_di.sum()} edges')
        return adj_di

    def internal_forward(self, h):
        for i, conv in enumerate(self.convs):
            h = conv(h, self.adj.t())
            if i != (len(self.convs) - 1):
                h = F.relu(h, inplace=True)
        return h

    def forward(self, features, eps=1e-8):
        h = self.internal_forward(features)
        h_norm = torch.linalg.vector_norm(h, ord=2, dim=1, keepdim=True)
        s1 = (h @ h.t()) / (h_norm @ h_norm.t() + eps)
        s1 = torch.where(s1 >= 0, s1, 0)  # symmetric, [0,1]

        s2 = self.s_d
        s = self.omega * s1 + (1 - self.omega) * s2
        s = torch.where(self.adj_mask, s, 0)
        s_norm = normalize_adj_symm(s)
        return s_norm, s


def generate_DI_matrices_for_all_roots(adj_mask, gene_exprs, save_path, potential_roots=None):

    if potential_roots is None:
        degrees = adj_mask.sum(dim=1)
        potential_roots = torch.where(degrees > torch.median(degrees))[0].cpu().numpy()
    
    di_matrices = {}
    for root in potential_roots:
        di_file = save_path + f'DI_root_{root}.csv'
        if not os.path.exists(di_file):
            model = GraphLearned(nlayers=2, isize=gene_exprs.shape[1], neighbor=30, 
                               gamma=1.0, adj=adj_mask, dis=adj_mask.float().numpy(),
                               device='cpu', omega=0.5, di_dir=save_path, 
                               expr=gene_exprs, percent=0.8, root_node=root)
            model.generate_DI_from_adj(adj_mask, gene_exprs, root, save_path)
    
    return di_matrices