import numpy as np
from collections import defaultdict

def norm_rowcol(matrix):
    # 按行求和
    row_norm=np.sum(matrix,axis=1).reshape(-1,1)
    # 行归一化
    matrix=matrix/row_norm
    # 按列求和
    col_norm=np.sum(matrix,axis=0)
    return matrix/col_norm

def sub_sample(graph,GAS, sampling_size,gene_size,gene_shape,cell_shape):
    cell_indexs=gene_shape+np.random.choice(np.arange(cell_shape),sampling_size,replace=False)
    sub_matrix=GAS[:,cell_indexs-gene_shape]
    gene_indexs=np.nonzero(np.sum(sub_matrix,axis=1))[0]

    sub_matrix=GAS[gene_indexs,:][:,cell_indexs-gene_shape]

    sub_matrix=norm_rowcol(sub_matrix)
    
    _indexs=np.argsort(np.sum(sub_matrix,axis=1))[::-1]
    gene_indexs=gene_indexs[_indexs]
    gene_indexs=gene_indexs[:gene_size]
    
    feature={
        'gene':graph.node_feature['gene'][gene_indexs,:],
        'cell':graph.node_feature['cell'][cell_indexs-gene_shape,:],
    }

    times={
        'gene': np.ones(gene_size),
        'cell':np.ones(sampling_size)
    }

    indxs={
        'gene':gene_indexs,
        'cell':cell_indexs-gene_shape
    }

    edge_list = defaultdict(  # target_type
        lambda: defaultdict(  # source_type
            lambda: defaultdict(  # relation_type
                lambda: []  # [target_id, source_id]
            )))

    for i in range(gene_size):
        edge_list['gene']['gene']['self'].append([i,i])

    for i in range(sampling_size):
        edge_list['cell']['cell']['self'].append([i,i])

    for i,cell_id in enumerate(cell_indexs):
        for j,gene_id in enumerate(gene_indexs):
            if gene_id in graph.edge_list['cell']['gene']['g_c'][cell_id]:
                edge_list['cell']['gene']['g_c'].append([i,j])
                edge_list['gene']['cell']['rev_g_c'].append([j,i])

    return feature, times, edge_list, indxs