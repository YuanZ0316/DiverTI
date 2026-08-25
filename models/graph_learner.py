import torch.nn as nn
import torch.nn.functional as F
import torch
from torch_geometric.nn import GCNConv
from .model_utils import normalize_adj_symm
import os  


class GraphLearner(nn.Module):
    def __init__(self, nlayers, isize, neighbor, gamma, adj, dis, device, omega, 
                 di_mode=False, di_dir=None, expr=None, percent=0.1, threshold=0.05):
        super().__init__()
        
        self.device = device
        self.omega = omega
        self.di_mode = di_mode
        self.neighbor = neighbor
        
        self.adj = adj.to(device)
        d_matrix = torch.tensor(dis, dtype=torch.float32, device=device)

        if di_mode:
            print("🚀 Initializing GraphLearner in DI mode...")
            self.di_matrix = self._load_di_matrix(di_dir, neighbor)
            self.adj_mask = self._create_di_adj_mask(self.di_matrix, neighbor)
            print(f"✅ DI matrix loaded with shape: {self.di_matrix.shape}")
        else:

            d_sorted, _ = d_matrix.sort()
            c1 = d_matrix > 0
            d_cut = torch.median(d_sorted[:, neighbor])
            c2 = d_matrix <= d_cut
            self.adj_mask = torch.logical_and(c1, c2)
        
        print('The number of useful edges is {}'.format(self.adj_mask.sum()))

        d_matrix = torch.where(self.adj_mask, d_matrix, torch.inf)
        d_cut = torch.median(d_matrix[d_matrix < torch.inf])
        d_matrix = d_matrix / d_cut
        self.s_d = 1 / torch.exp(gamma * torch.pow(d_matrix, 2))
        

        self.convs = nn.ModuleList()
        self.convs = torch.nn.ModuleList(
            [GCNConv(in_channels=isize, out_channels=isize)] +
            [GCNConv(in_channels=isize, out_channels=isize) for _ in range(nlayers - 1)])
        self.input_dim = isize

    def _load_di_matrix(self, di_dir, neighbor):
        di_path = os.path.join(di_dir, 'global_DI_matrix.csv')
        
        if os.path.exists(di_path):
            print(f"📥 Loading DI matrix from {di_path}")
            di_matrix = np.loadtxt(di_path, delimiter=',')
            di_matrix = torch.tensor(di_matrix, dtype=torch.float32, device=self.device)
            return di_matrix
        else:
            print(f"❌ DI matrix not found at {di_path}")
            raise FileNotFoundError(f"DI matrix not found at {di_path}")

    def _create_di_adj_mask(self, di_matrix, neighbor):
        n_nodes = di_matrix.shape[0]

        mask = torch.zeros_like(di_matrix, dtype=torch.bool)
        
        for i in range(n_nodes):
            _, top_indices = torch.topk(di_matrix[i], min(neighbor + 1, n_nodes))
            top_indices = top_indices[top_indices != i]
            if len(top_indices) > neighbor:
                top_indices = top_indices[:neighbor]
            
            mask[i, top_indices] = True
            mask[top_indices, i] = True  
        print(f"✅ DI-based adjacency mask created with {mask.sum()} edges")
        return mask

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
        s1 = torch.where(s1 >= 0, s1, 0)

        s2 = self.s_d  
        
        if self.di_mode:
            s_di = torch.where(self.adj_mask, self.di_matrix, 0)
            s = 0.7 * s_di + 0.2 * s1 + 0.1 * s2
        else:
            s = self.omega * s1 + (1 - self.omega) * s2
        
        s = torch.where(self.adj_mask, s, 0)
        s_norm = normalize_adj_symm(s)
        return s_norm, s