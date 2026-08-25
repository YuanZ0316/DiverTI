import networkx as nx
import pandas as pd
import numpy as np
import os
import math
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

class BinaryDITree:
    def __init__(self, adata, root, predict_key, threshold=0.1, seed=0, save_dir=None):
        self.seed = seed
        self.tree = nx.DiGraph()
        self.root = root
        self.save_dir = save_dir
        self.threshold = threshold  # DI阈值（小的值表示相似度高）
        self.predict_key = predict_key
        self.adata = adata
        self.debug = False
        self.connectivity_matrix = None  # 新增：存储最终的有向树邻接矩阵
    def read_express(self, cluster_node=None):
        if not isinstance(self.adata.X, np.ndarray):
            express = pd.DataFrame(self.adata.X.toarray(), index=self.adata.obs.index, 
                                 columns=self.adata.var.index, dtype=np.float64)
        else:
            express = pd.DataFrame(self.adata.X, index=self.adata.obs.index, 
                                 columns=self.adata.var.index, dtype=np.float64)
        
        express = pd.concat([express, self.adata.obs[self.predict_key]], axis=1)
        express = express.groupby(self.predict_key).mean()
        express[self.predict_key] = express.index
        express = express.drop(self.predict_key, axis=1)
        express.columns = express.columns.astype(str)
        express.index = express.index.astype(str)
        
        if cluster_node is not None:
            express = express.loc[cluster_node]
        else:
            cluster_node = list(set(express.index.values.tolist()))
        
        return express, cluster_node

    def D(self, X):
        N = 0
        ans = 0
        for x in X:
            if x == 0:  # 忽略log(0)
                continue
            else:
                n = x
                N = N + n
                ans += n * math.log(n)
        
        if N == 0:
            return 0
        
        Dx = N * math.log(N) - ans
        return Dx

    def calculate_di(self, node_i, node_j):
        try:
            expr_i = self.tree.nodes[node_i]['express']
            expr_j = self.tree.nodes[node_j]['express']

            if len(expr_i.shape) > 1:
                expr_i = expr_i.flatten()
            if len(expr_j.shape) > 1:
                expr_j = expr_j.flatten()

            expr_combined = expr_i + expr_j 

            D_i = self.D(expr_i)
            D_j = self.D(expr_j)
            D_combined = self.D(expr_combined)

            di = D_combined - D_i - D_j
            
            if self.debug:
                print(f" {node_i} vs {node_j}")
                print(f"  D_i={D_i:.2f}, D_j={D_j:.2f}, D_combined={D_combined:.2f}, DI={di:.6f}")
            return di
            
        except Exception as e:
            print(f"Error calculating DI for {node_i}, {node_j}: {e}")
            return float('inf')  
    def cal_global_di_matrix(self, nodes):
        if not os.path.exists(f'{self.save_dir}/trial{self.seed}/'):
            os.makedirs(f'{self.save_dir}/trial{self.seed}/')
        
        di_path = f'{self.save_dir}/trial{self.seed}/global_DI_matrix.csv'
        print(f" {di_path}")
        
        if os.path.exists(di_path):
            di_df = pd.read_csv(di_path, index_col=0, header=0)
            di_df.columns = di_df.columns.astype(str)
            di_df.index = di_df.index.astype(str)
        else:
            di_df = pd.DataFrame(np.full((len(nodes), len(nodes)), float('inf')), 
                               index=nodes, columns=nodes)
            
            for i, node_i in enumerate(nodes):
                for j, node_j in enumerate(nodes):
                    if i < j:  
                        di_value = self.calculate_di(node_i, node_j)
                        di_df.loc[node_i, node_j] = di_value
                        di_df.loc[node_j, node_i] = di_value
                if i % 5 == 0:  
                    print(f"{i+1}/{len(nodes)}")
            
            np.fill_diagonal(di_df.values, 0)
            
            di_df.to_csv(di_path, index=True, header=True)
        return di_df

    def build_tree_from_di(self, nodes):

        di_matrix = self.cal_global_di_matrix(nodes)

        G_full = nx.Graph()
        for node in nodes:
            G_full.add_node(node)

        for i, node_i in enumerate(nodes):
            for j, node_j in enumerate(nodes):
                if i < j:
                    weight = di_matrix.loc[node_i, node_j]
                    G_full.add_edge(node_i, node_j, weight=weight)
        
        print(f"{G_full.number_of_nodes()} , {G_full.number_of_edges()} ")

        mst = nx.minimum_spanning_tree(G_full, weight='weight')
        print(f" {mst.number_of_edges()} ")

        directed_tree = nx.bfs_tree(mst, self.root)

        self.tree.clear_edges()
        for edge in directed_tree.edges():
            self.tree.add_edge(edge[0], edge[1])
            if self.debug:
                di_value = di_matrix.loc[edge[0], edge[1]]
                print(f"{edge[0]} -> {edge[1]}, {di_value:.6f}")
        
        print(f"{self.tree.number_of_edges()} ")

    def binary_cut_tree(self):
        print("")

        for node in list(self.tree.nodes()):
            children = list(self.tree.successors(node))
            
            if len(children) > 2:
                print(f" {node} {len(children)} ")
                self._split_node(node, children)


    def _split_node(self, parent, children):

        if len(children) < 3:
            return

        di_values = {}
        for i, child_i in enumerate(children):
            for j, child_j in enumerate(children):
                if i < j:
                    if hasattr(self, 'global_di_matrix'):
                        di_val = self.global_di_matrix.loc[child_i, child_j]
                    else:
                        di_val = self.calculate_di(child_i, child_j)
                    di_values[(child_i, child_j)] = di_val
        
        groups = self._hierarchical_clustering(children, di_values)
        
        if len(groups) == 2:
            rep1 = self._find_representative(parent, groups[0])
            rep2 = self._find_representative(parent, groups[1])

            self._reorganize_connections(parent, rep1, rep2, groups[0], groups[1])

    def _hierarchical_clustering(self, nodes, di_values):
        if len(nodes) <= 2:
            return [nodes]
        
        clusters = [[node] for node in nodes]
        
        while len(clusters) > 2:
            min_di = float('inf')
            merge_i, merge_j = -1, -1
            
            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    avg_di = self._average_di_between_clusters(clusters[i], clusters[j], di_values)
                    if avg_di < min_di:
                        min_di = avg_di
                        merge_i, merge_j = i, j
            
            if merge_i != -1 and merge_j != -1:
                clusters[merge_i].extend(clusters[merge_j])
                clusters.pop(merge_j)
        
        return clusters

    def _average_di_between_clusters(self, cluster1, cluster2, di_values):
        total_di = 0
        count = 0
        
        for node1 in cluster1:
            for node2 in cluster2:
                key = (node1, node2) if node1 < node2 else (node2, node1)
                if key in di_values:
                    total_di += di_values[key]
                    count += 1
        
        return total_di / count if count > 0 else float('inf')

    def _find_representative(self, parent, nodes):
        if not nodes:
            return None

        similarities = {}
        for node in nodes:
            if hasattr(self, 'global_di_matrix'):
                di_val = self.global_di_matrix.loc[parent, node]
            else:
                di_val = self.calculate_di(parent, node)
            similarities[node] = di_val

        return min(similarities, key=similarities.get)

    def _reorganize_connections(self, parent, rep1, rep2, group1, group2):

        for child in list(self.tree.successors(parent)):
            self.tree.remove_edge(parent, child)

        self.tree.add_edge(parent, rep1)
        self.tree.add_edge(parent, rep2)

        for node in group1:
            if node != rep1:

                for pred in list(self.tree.predecessors(node)):
                    self.tree.remove_edge(pred, node)
                self.tree.add_edge(rep1, node)
        
        for node in group2:
            if node != rep2:
                for pred in list(self.tree.predecessors(node)):
                    self.tree.remove_edge(pred, node)
                self.tree.add_edge(rep2, node)

    def init_tree(self, cluster_list=None, debug=False):

        self.debug = debug
        express, nodes = self.read_express(cluster_list)
        print(f"{nodes}")
        
        self.tree.clear()
        self.tree.add_nodes_from(nodes)
        
        for node in nodes:
            self.tree.nodes[node]['express'] = express.loc[node].values

        self.tree.nodes[self.root]['is_root'] = True

        self.build_tree_from_di(nodes)

        self.binary_cut_tree()
        
        print(f" {self.tree.number_of_nodes()} , {self.tree.number_of_edges()} ")



    def generate_connectivity_matrix(self):

        if self.tree.number_of_edges() == 0:
            print("Warning: Directed tree not built; unable to generate adjacency matrix.")
            return None
        nodes = sorted([str(node) for node in self.tree.nodes()])
        conn_matrix = pd.DataFrame(
            0.0,
            index=nodes,
            columns=nodes,
            dtype=np.float64
        )
        for start_node, end_node in self.tree.edges():
            start_node_str = str(start_node)
            end_node_str = str(end_node)
            conn_matrix.loc[start_node_str, end_node_str] = 1.0  
            if self.debug:
                print(f"{start_node_str} → {end_node_str}")
        
        self.connectivity_matrix = conn_matrix
        self.adata.uns['cascat_directed_tree_conn'] = conn_matrix  
        print(f"adata.uns['cascat_directed_tree_conn']")
        print(conn_matrix)
        return conn_matrix

    def construct_tree(self):
        
        # 1. 初始化树（原有逻辑不变）
        self.init_tree()
        
        # 2. 确保所有边都有权重（原有逻辑不变）
        for edge in self.tree.edges():
            self.tree.edges[edge]['weight'] = 1
        
        # 3. 新增：自动生成绘图用的邻接矩阵
        self.generate_connectivity_matrix()
        
        return self.tree