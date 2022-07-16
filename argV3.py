import os
import random
import time
import argparse
import dill as pickle 

import numpy as np
import pandas as pd
import torch
import torch.utils.data as data
from torch import nn, optim
from torch.nn import functional as F

from reduction import reduction
from utils import debuginfoStr,loadGAS ,build_data, build_graph
from sub_sample import sub_sample
from pyHGT.model import GNN, GNN_from_raw

from warnings import filterwarnings
filterwarnings("ignore")

seed = 0
random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)

torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

print('cuda version: ',torch.version.cuda)

# list_in = r.list_in

# arguments
parser = argparse.ArgumentParser(description='Training GNN on gene cell graph')
parser.add_argument('--data_path', type=str)
parser.add_argument('--epoch', type=int, default=100)
# sampling times
parser.add_argument('--n_batch', type=int, default=25,
                    help='Number of batch (sampled graphs) for each epoch')

parser.add_argument('--cell_rate', type=float, default=0.9)
parser.add_argument('--gene_rate', type=float, default=0.3)

# Result
parser.add_argument('--data_name', type=str,
                    help='The name for dataset')
parser.add_argument('--result_dir', type=str,
                    help='The address for storing the models and optimization results.')
parser.add_argument('--reduction', type=str, default='AE',
                    help='the method for feature extraction, pca, raw, AE')
parser.add_argument('--in_dim', type=int, default=256,
                    help='Number of hidden dimension (AE)')
# GAE
parser.add_argument('--n_hid', type=int,
                    help='Number of hidden dimension')
parser.add_argument('--n_heads', type=int,
                    help='Number of attention head')
parser.add_argument('--n_layers', type=int, default=2,
                    help='Number of GNN layers')
parser.add_argument('--dropout', type=float, default=0,
                    help='Dropout ratio')
parser.add_argument('--lr', type=float,
                    help='learning rate')

parser.add_argument('--batch_size', type=int,
                    help='Number of output nodes for training')
parser.add_argument('--layer_type', type=str, default='hgt',
                    help='the layer type for GAE')
parser.add_argument('--loss', type=str, default='kl',
                    help='the loss for GAE')
parser.add_argument('--factor', type=float, default='0.5',
                    help='the attenuation factor')
parser.add_argument('--patience', type=int, default=5,
                    help='patience')
parser.add_argument('--rf', type=float, default='0.0',
                    help='the weights of regularization')
parser.add_argument('--cuda', type=int, default=0,
                    help='cuda 0 use GPU0 else cpu ')
parser.add_argument('--rep', type=str, default='T',
                    help='precision truncation')
parser.add_argument('--AEtype', type=int, default=1,
                    help='AEtype:1 embedding node autoencoder 2:HGT node autoencode')
parser.add_argument('--optimizer', type=str, default='adamw',
                    help='optimizer')

args = parser.parse_args()

# gene_cell=list_in['cell_gene']
# args.data_name=list_in['data_name']
# args.n_hid=int(list_in['n_hid'])
# args.n_heads=int(list_in['n_heads'])
# args.n_layers=int(list_in['n_layers'])
# args.lr=float(list_in['lr'])
# args.n_batch=int(list_in['n_batch'])
# args.epoch=int(list_in['epoch'])
# args.result_dir=list_in['result_dir']
# args.cuda=int(list_in['cuda'])
# args.data_type=list_in['data_type']
# args.cell_rate=float(list_in['cell_rate'])
# args.gene_rate=float(list_in['gene_rate'])
# args.gene_name=list_in['gene_name']
# args.cell_name=list_in['cell_name']

gene_cell=loadGAS(args.data_path)

print(f'GAS loaded!',gene_cell.shape)
# dataname,n_hid,n_heads,n_layers,lr,n_batch,batch_size
file0=f'epoch_{args.epoch}_n_hid_{args.n_hid}_nheads_{args.n_heads}_lr_{args.lr}n_batch{args.n_batch}'
print(f'\n{file0}')

gene_dir = args.result_dir+'/gene/'
cell_dir = args.result_dir+'/cell/'
model_dir = args.result_dir+'/model/'
att_dir = args.result_dir+'/att/'

start_time = time.time()
print('---0:00:00---scRNA starts.')

if args.cuda == 0:
    device = torch.device("cuda:" + "0")
    print("cuda>>>")
else:
    device = torch.device("cpu")
print(device)

encoded,encoded2= reduction(args.reduction,gene_cell,device)
debuginfoStr('Feature extraction finished')

graph=build_graph(gene_cell,encoded,encoded2)
debuginfoStr('Build Graph finished')

print("Start sampling!")
np.random.seed(seed)
jobs = []
cell_num=int((gene_cell.shape[1]*args.cell_rate)/args.n_batch)
gene_num=int((gene_cell.shape[0]*args.gene_rate)/args.n_batch)
print(f'cell_num: {cell_num}, gene_num: {gene_num}')
for _ in range(args.n_batch):
    p = sub_sample(graph,
                    gene_cell,
                    cell_num,
                    gene_num,
                    gene_cell.shape[0],
                    gene_cell.shape[1])
    jobs.append(p)
print("Sampling end!")
debuginfoStr('Cell Graph constructed and pruned')


if (args.reduction != 'raw'):
    gnn = GNN(conv_name=args.layer_type, in_dim=encoded.shape[1],
              n_hid=args.n_hid, n_heads=args.n_heads, n_layers=args.n_layers, dropout=args.dropout,
              num_types=2, num_relations=2, use_RTE=False).to(device)
else:
    gnn = GNN_from_raw(conv_name=args.layer_type, in_dim=[encoded.shape[1], encoded2.shape[1]],
                       n_hid=args.n_hid, n_heads=args.n_heads, n_layers=args.n_layers, dropout=args.dropout,
                       num_types=2, num_relations=2, use_RTE=False,
                       AEtype=args.AEtype).to(device)

# default: adamw
if args.optimizer == 'adamw':
    optimizer = torch.optim.AdamW(gnn.parameters(), lr=args.lr)
elif args.optimizer == 'adam':
    optimizer = torch.optim.Adam(gnn.parameters(), lr=args.lr)
elif args.optimizer == 'sgd':
    optimizer = torch.optim.SGD(gnn.parameters(), lr=args.lr)
elif args.optimizer == 'adagrad':
    optimizer = torch.optim.Adagrad(gnn.parameters(), lr=args.lr)

scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, 'min', factor=args.factor, patience=args.patience, verbose=True)

gnn.train()
for epoch in np.arange(args.epoch):
    L = 0
    for job in jobs:
        feature,time,edge_list,indxs = job
        node_dict = {}
        node_feature = []
        node_type = []
        node_time = []
        edge_index = []
        edge_type = []
        edge_time = []

        node_num = 0
        types = graph.get_types()   # ['gene','cell']
        for t in types:
            #print("t in types "+str(t)+"\n")
            node_dict[t] = [node_num, len(node_dict)]
            node_num += len(feature[t])
            if args.reduction == 'raw':
                node_feature.append([])
        # node_dict: {'gene':[0,0],'cell':[134,1]}
        for t in types:
            t_i = node_dict[t][1]
            #print("feature t:\n")
            #print("t_i="+str(t_i)+" t="+str(t)+"\n")
            # print(feature[t].shape)
            if args.reduction != 'raw':
                node_feature += list(feature[t])
            else:
                node_feature[t_i] = torch.tensor(
                    feature[t], dtype=torch.float32).to(device)

            node_time += list(time[t])
            node_type += [node_dict[t][1] for _ in range(len(feature[t]))]
        edge_dict = {e[2]: i for i, e in enumerate(graph.get_meta_graph())}
        edge_dict['self'] = len(edge_dict)
        # {'g_c': 0, 'rev_g_c': 1 ,'self': 2}
        for target_type in edge_list:
            for source_type in edge_list[target_type]:
                for relation_type in edge_list[target_type][source_type]:
                    for ii, (ti, si) in enumerate(edge_list[target_type][source_type][relation_type]):
                        tid, sid = ti + \
                            node_dict[target_type][0], si + \
                            node_dict[source_type][0]
                        edge_index += [[sid, tid]]
                        edge_type += [edge_dict[relation_type]]

                        # Our time ranges from 1900 - 2020, largest span is 120.
                        # edge_time += [node_time[tid] - node_time[sid] + 120]
                        edge_time += [120]

        if (args.reduction != 'raw'):
            node_feature = torch.stack(node_feature)
            node_feature = torch.tensor(node_feature, dtype=torch.float32)
            node_feature = node_feature.to(device)

        #node_feature = torch.trunc(node_feature*10000)/10000
        node_type = torch.LongTensor(node_type)
        edge_time = torch.LongTensor(edge_time)
        edge_index = torch.LongTensor(edge_index).t()
        edge_type = torch.LongTensor(edge_type)
        if (args.reduction == 'raw'):
            node_rep, node_decoded_embedding = gnn.forward(node_feature, 
                                                           node_type.to(device),
                                                           edge_time.to(device),
                                                           edge_index.to(device),
                                                           edge_type.to(device))
        else:
            node_rep = gnn.forward(node_feature, 
                                   node_type.to(device),
                                   edge_time.to(device),
                                   edge_index.to(device),
                                   edge_type.to(device))

        if args.rep == 'T':
            node_rep = torch.trunc(node_rep*10000000000)/10000000000
            if args.reduction == 'raw':
                for t in types:
                    t_i = node_dict[t][1]
                    # print("t_i="+str(t_i))
                    node_decoded_embedding[t_i] = torch.trunc(
                        node_decoded_embedding[t_i]*10000000000)/10000000000


        gene_matrix = node_rep[node_type == 0, ]
        cell_matrix = node_rep[node_type == 1, ]

        regularization_loss = 0
        for param in gnn.parameters():
            regularization_loss += torch.sum(torch.pow(param, 2))
        if (args.loss == "kl"):
            decoder = torch.mm(gene_matrix, cell_matrix.t())
            adj = gene_cell[indxs['gene'], ]
            adj = adj[:, indxs['cell']]
            adj = torch.tensor(adj, dtype=torch.float32).to(device)
            if args.reduction == 'raw':
                if epoch % 2 == 0:
                    loss = F.kl_div(decoder.softmax(
                        dim=-1).log(), adj.softmax(dim=-1), reduction='sum')+args.rf*regularization_loss
                else:
                    loss = nn.MSELoss()(
                        node_feature[0], node_decoded_embedding[0])+args.rf*regularization_loss
                    for t_i in range(1, len(types)):
                        loss += nn.MSELoss()(node_feature[t_i],
                                             node_decoded_embedding[t_i])
            else:
                loss = F.kl_div(decoder.softmax(dim=-1).log(),
                                adj.softmax(dim=-1), reduction='sum')

        if (args.loss == "cross"):
            # negative_sampling not defined
            print("negative_sampling not defined!")
            exit()
            pass

        L += loss.item()
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    scheduler.step(L/(int(gene_cell.shape[0])))
    print('Epoch :', epoch+1, '|', 'train_loss:%.12f' %
          (L/(int(gene_cell.shape[0]))/args.n_batch))


state = {'model': gnn.state_dict(), 'optimizer': scheduler.state_dict(),
         'epoch': epoch}
torch.save(state, model_dir+file0)
debuginfoStr('Graph Autoencoder training finished')

debuginfoStr('load training model')
state = torch.load(model_dir+file0, map_location=lambda storage, loc: storage)
device = torch.device("cpu")

if (args.reduction != 'raw'):
    gnn = GNN(conv_name=args.layer_type, in_dim=encoded.shape[1], n_hid=args.n_hid, n_heads=args.n_heads, n_layers=args.n_layers, dropout=args.dropout,
              num_types=2, num_relations=2, use_RTE=False).to(device)
else:
    gnn = GNN_from_raw(conv_name=args.layer_type, in_dim=[encoded.shape[1], encoded2.shape[1]], n_hid=args.n_hid, n_heads=args.n_heads, n_layers=args.n_layers, dropout=args.dropout,
                       num_types=2, num_relations=2, use_RTE=False,
                       AEtype=args.AEtype).to(device)


# model.eval()
if (gene_cell.shape[1]>10000):

    if (gene_cell.shape[0]>10000):
        ba = 500
    else:
        ba = gene_cell.shape[0]
else:
    if (gene_cell.shape[0]>10000):
        ba = 5000
    else:
        ba = gene_cell.shape[0]

gnn.load_state_dict(state['model'])
g_embedding = []
gene_name = []
cell_name = []
attention = []

with torch.no_grad():
    for i in range(0, gene_cell.shape[0], ba):
        adj = gene_cell[i:(i+ba), :]  
        x,node_type, edge_time, edge_index,edge_type=build_data(adj,encoded[i:(ba+i), :],encoded2)
        if args.reduction != 'raw':
            node_rep = gnn.forward((torch.cat((x['gene'], x['cell']), 0)).to(device), 
            node_type.to(device),edge_time.to(device),
            edge_index.to(device), edge_type.to(device))
        else:
            node_rep, _ = gnn.forward([x['gene'].to(device), x['cell'].to(device)], 
                                       node_type.to(device),edge_time.to(device), 
                                       edge_index.to(device), edge_type.to(device))

        gene_name = gene_name + list(np.array(edge_index[0]+i))
        cell_name = cell_name + list(np.array(edge_index[1]-adj.shape[0]))
        attention.append(gnn.att)
        gene_matrix = node_rep[node_type == 0, ]
        cell_matrix = node_rep[node_type == 1, ]
        g_embedding.append(gene_matrix)

if gene_cell.shape[0] % ba == 0:
    gene_matrix = np.vstack(g_embedding[0:int(gene_cell.shape[0]/ba)])
    attention = np.vstack(attention[0:int(gene_cell.shape[0]/ba)])
else:
    final_tensor = np.vstack(g_embedding[0:int(gene_cell.shape[0]/ba)])
    gene_matrix = np.concatenate((final_tensor, gene_matrix), 0)
    final_attention = np.vstack(attention[0:int(gene_cell.shape[0]/ba)])
    attention = np.concatenate((final_attention, gnn.att), 0)
cell_matrix = cell_matrix.detach().numpy()
np.savetxt(gene_dir+file0, gene_matrix, delimiter=' ')
np.savetxt(cell_dir+file0, cell_matrix, delimiter=' ')
'''
g = np.nonzero(gene_cell)[0]
c = np.nonzero(gene_cell)[1]+gene_cell.shape[0]
name1 = pd.DataFrame(
    gene_name[0:torch.tensor([g, c]).shape[1]], columns=['gene'])
name2 = pd.DataFrame(
    cell_name[0:torch.tensor([g, c]).shape[1]], columns=['cell'])
df = pd.DataFrame(attention)
df2 = pd.concat([name1, name2, df], axis=1)
attention = df2
df2.to_csv(att_dir+file0, sep=",", index=True)
import time
debuginfoStr(f'Finished! time: {(time.time()-start_time)/60} min')
'''
