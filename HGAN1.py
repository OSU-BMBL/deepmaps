# from memory_profiler import profile
import time
import resource
import datetime
import numpy as np
import random
import torch
from torch import nn, optim
from torch.nn import functional as F
import os
import sys
from pyHGT.data import *
from warnings import filterwarnings

filterwarnings("ignore")
from pyHGT.model import *
from torch_geometric.data import Data

seed = 0

random.seed(seed)
# np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
np.random.seed(seed)
os.environ["PYTHONHASHSEED"] = str(seed)

# torch.cuda.manual_seed(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
import torch.utils.data as data
import argparse

parser = argparse.ArgumentParser(description="Training GNN on gene cell graph")

parser.add_argument(
    "--result_dir",
    type=str,
    default="default.txt",
    help="The address for storing the models and optimization results.",
)
parser.add_argument(
    "--input_dir",
    type=str,
    default="default.txt",
    help="The address for storing the models and optimization results.",
)
parser.add_argument(
    "--in_dim", type=int, default=256, help="Number of hidden dimension"
)
parser.add_argument("--n_hid", type=int, default=64, help="Number of hidden dimension")
parser.add_argument("--n_heads", type=int, default=8, help="Number of attention head")
parser.add_argument("--n_layers", type=int, default=4, help="Number of GNN layers")
parser.add_argument("--dropout", type=float, default=0, help="Dropout ratio")
parser.add_argument(
    "--sample_depth", type=int, default=6, help="How many numbers to sample the graph"
)
parser.add_argument(
    "--sample_width",
    type=int,
    default=520,
    help="How many nodes to be sampled per layer per type",
)
parser.add_argument("--lr", type=float, default=1e-3, help="learning rate")
parser.add_argument(
    "--n_batch",
    type=int,
    default=64,
    help="Number of batch (sampled graphs) for each epoch",
)
parser.add_argument(
    "--batch_size", type=int, default=128, help="Number of output nodes for training"
)

args = parser.parse_args()


file0 = (
    "2w_2w_2_n_batch"
    + str(args.n_batch)
    + "_batch_size_"
    + str(args.batch_size)
    + "sample_depth_"
    + str(args.sample_depth)
    + "_nheads_"
    + str(args.n_heads)
    + "_nlayers_"
    + str(args.n_layers)
    + "_sample_width_"
    + str(args.sample_width)
    + "_lr_"
    + str(args.lr)
    + "_n_hid_"
    + str(args.n_hid)
)
print(file0)
gene_dir = args.result_dir + "/gene/"
cell_dir = args.result_dir + "/cell/"
model_dir = args.result_dir + "/model/"


def debuginfoStr(info):

    # print('---'+str(datetime.timedelta(seconds=int(time.time()-start_time)))+'---'+info)
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print("Mem consumption: " + str(mem))


def load_data(path, sep, col_name, row_name):
    f = open(path, "r")
    sourceInLine = f.readlines()
    f.close()
    gene_cell = []
    for line in sourceInLine:
        temp1 = line.strip("\n")
        temp2 = temp1.split(sep)
        gene_cell.append(temp2)
    if col_name == True:
        cell_name = gene_cell[0]
        del gene_cell[0]
    else:
        cell_name = ""

    if row_name == True:
        gene_cell = np.array(gene_cell)
        gene_name = gene_cell[:, 0]
        gene_cell = gene_cell[:, 1 : gene_cell.shape[1] + 1]
    else:
        gene_name = ""
    gene_cell = np.array(gene_cell)
    print(
        "The number of gene is {}, The number of cell is {}".format(
            gene_cell.shape[0], gene_cell.shape[1]
        )
    )
    return (gene_cell, gene_name, cell_name)


class EarlyStopping:
    def __init__(self, patience=7, verbose=False, delta=0):
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.delta = delta

    def __call__(self, val_loss, model):

        score = -val_loss

        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
        elif score < self.best_score + self.delta:
            self.counter += 1
            print(f"EarlyStopping counter: {self.counter} out of {self.patience}")
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
            self.counter = 0

    def save_checkpoint(self, val_loss, model):
        """
        Saves model when validation loss decrease.
        """
        if self.verbose:
            print(
                f"Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f}).  Saving model ..."
            )
        state = {
            "model": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "epoch": epoch,
        }
        torch.save(state, model_dir + file0)
        # torch.save(model.state_dict(), 'checkpoint.pt')
        # torch.save(model, 'finish_model.pkl')
        self.val_loss_min = val_loss


start_time = time.time()
print("---0:00:00---scRNA starts loading.")
debuginfoStr("scRNA has been successfully loaded")
gene_cell, gene_name, cell_name = load_data(
    args.input_dir, sep=" ", col_name=True, row_name=True
)
gene_cell = gene_cell.astype("float")
debuginfoStr("scRNA has been successfully loaded")
cuda = -1
if cuda == -1:
    device = torch.device("cuda:" + "0")
    print("cuda>>>")
else:
    device = torch.device("cpu")
print(device)


class AE(nn.Module):
    def __init__(self, dim):
        super(AE, self).__init__()
        self.dim = dim
        self.fc1 = nn.Linear(dim, 512)
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, 512)
        self.fc4 = nn.Linear(512, dim)

    def encode(self, x):
        h1 = F.relu(self.fc1(x))
        return F.relu(self.fc2(h1))
        return h1

    def decode(self, z):
        h3 = F.relu(self.fc3(z))
        return torch.relu(self.fc4(h3))
        # return torch.relu(self.fc4(z))

    def forward(self, x):
        z = self.encode(x.view(-1, self.dim))
        return self.decode(z), z


gene = torch.tensor(gene_cell, dtype=torch.float32).to(device)
if gene_cell.shape[0] < 5000:
    ba = gene_cell.shape[0]
else:
    ba = 5000
loader1 = data.DataLoader(gene, ba)

EPOCH_AE = 1500
model = AE(dim=gene.shape[1]).to(device)
optimizer = optim.Adam(model.parameters(), lr=1e-3)
loss_func = nn.MSELoss()
for epoch in range(EPOCH_AE):
    embedding1 = []
    for _, batch_x in enumerate(loader1):

        decoded, encoded = model(batch_x)
        # encoded1 , decoded1 = Coder2(cell)
        loss = loss_func(batch_x, decoded)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        embedding1.append(encoded)
    # print('Epoch :', epoch,'|','train_loss:%.4f'%loss.data)
if gene.shape[0] % ba != 0:
    torch.stack(embedding1[0 : int(gene.shape[0] / ba)])
    a = torch.stack(embedding1[0 : int(gene.shape[0] / ba)])
    a = a.view(ba * int(gene.shape[0] / ba), 256)
    encoded = torch.cat((a, encoded), 0)

else:
    encode = torch.stack(embedding1)
    encoded = encode.view(gene.shape[0], 256)

if gene_cell.shape[1] < 5000:
    ba = gene_cell.shape[1]
else:
    ba = 5000
cell = torch.tensor(np.transpose(gene_cell), dtype=torch.float32).to(device)
loader2 = data.DataLoader(cell, ba)
model2 = AE(dim=cell.shape[1]).to(device)
optimizer2 = optim.Adam(model2.parameters(), lr=1e-3)
for epoch in range(EPOCH_AE):
    embedding1 = []
    for _, batch_x in enumerate(loader2):
        decoded2, encoded2 = model2(batch_x)
        loss = loss_func(batch_x, decoded2)
        optimizer2.zero_grad()
        loss.backward()
        optimizer2.step()
        embedding1.append(encoded2)
    # print('Epoch :', epoch,'|','train_loss:%.4f'%loss.data)
if cell.shape[0] % ba != 0:
    torch.stack(embedding1[0 : int(cell.shape[0] / ba)])
    a = torch.stack(embedding1[0 : int(cell.shape[0] / ba)])
    a = a.view(ba * int(cell.shape[0] / ba), 256)
    encoded2 = torch.cat((a, encoded2), 0)
    # encode.shape
else:
    encode = torch.stack(embedding1)
    encoded2 = encode.view(cell.shape[0], 256)
debuginfoStr("Feature autoencoder training finished")
print(encoded)
print(encoded2)
gnn = GNN(
    conv_name="hgt",
    in_dim=args.in_dim,
    n_hid=args.n_hid,
    n_heads=args.n_heads,
    n_layers=args.n_layers,
    dropout=args.dropout,
    num_types=2,
    num_relations=2,
).to(device)
args_optimizer = "adamw"
if args_optimizer == "adamw":
    optimizer = torch.optim.AdamW(gnn.parameters(), lr=args.lr)
elif args_optimizer == "adam":
    optimizer = torch.optim.Adam(gnn.parameters(), lr=0.1)
elif args_optimizer == "sgd":
    optimizer = torch.optim.SGD(gnn.parameters(), lr=0.01, weight_decay=0.01)
elif args_optimizer == "adagrad":
    optimizer = torch.optim.Adagrad(gnn.parameters(), lr=0.001, weight_decay=0.1)

debuginfoStr("Start construct cell grpah")
target_nodes = np.arange(gene_cell.shape[1] + gene_cell.shape[0])
# gene cell
g = np.nonzero(gene_cell)[0]
c = np.nonzero(gene_cell)[1] + gene_cell.shape[0]
edge1 = list(g)
edge2 = list(c)
# print(len(edge1))
# print(len(edge2))
# node_feature = torch.cat( (encoded, encoded2), 0)
# node_feature = (node_feature-torch.mean(node_feature))/torch.std(node_feature)
edge_index = torch.tensor([edge1, edge2], dtype=torch.long)
# x={'gene': torch.tensor(node_feature[0:encoded.shape[0],:], dtype=torch.float),
#    'cell': torch.tensor(node_feature[encoded.shape[0]:(encoded2.shape[0]+encoded.shape[0]),:], dtype=torch.float),
# }
x = {
    "gene": torch.tensor(encoded, dtype=torch.float),
    "cell": torch.tensor(encoded2, dtype=torch.float),
}
# x = torch.tensor(np.random.randn(len(target_nodes),64), dtype=torch.float)
# h = Data(edge_index_dict=edge_index, x=x)
edge_index_dict = {("gene", "g_c", "cell"): torch.tensor([g, c], dtype=torch.long)}
edge_reltype = {("gene", "g_c", "cell"): torch.tensor([g, c]).shape[1]}
num_nodes_dict = {"gene": gene_cell.shape[0], "cell": gene_cell.shape[1]}
data = Data(
    edge_index_dict=edge_index_dict,
    edge_reltype=edge_reltype,
    num_nodes_dict=num_nodes_dict,
    x=x,
)
graph = Graph()
edg = graph.edge_list
edge_index_dict = data.edge_index_dict
for key in edge_index_dict:
    # print(key)
    edges = edge_index_dict[key]
    s_type, r_type, t_type = key[0], key[1], key[2]
    elist = edg[t_type][s_type][r_type]
    rlist = edg[s_type][t_type]["rev_" + r_type]
    for s_id, t_id in edges.t().tolist():
        year = 1
        elist[t_id][s_id] = year
        rlist[s_id][t_id] = year
edg = {}
deg = {key: np.zeros(data.num_nodes_dict[key]) for key in data.num_nodes_dict}

for k1 in graph.edge_list:
    if k1 not in edg:
        edg[k1] = {}
    for k2 in graph.edge_list[k1]:
        if k2 not in edg[k1]:
            edg[k1][k2] = {}
        for k3 in graph.edge_list[k1][k2]:
            if k3 not in edg[k1][k2]:
                edg[k1][k2][k3] = {}
            for num1, e1 in enumerate(graph.edge_list[k1][k2][k3]):
                if len(graph.edge_list[k1][k2][k3][e1]) == 0:
                    continue

                edg[k1][k2][k3][num1] = {}
                for num2, e2 in enumerate(graph.edge_list[k1][k2][k3][e1]):
                    edg[k1][k2][k3][num1][num2] = graph.edge_list[k1][k2][k3][e1][e2]
                deg[k1][num1] += len(edg[k1][k2][k3][num1])
            # print(k1, k2, k3, len(edg[k1][k2][k3]))

graph.node_feature["gene"] = data.x["gene"]
graph.node_feature["cell"] = data.x["cell"]
print(graph.node_feature["gene"].shape)
print(graph.node_feature["cell"].shape)
debuginfoStr("Cell Graph constructed and pruned")


def sub_sample(samp_nodes, sampled_depth=2, sampled_number=10):
    inp = {
        "gene": np.concatenate([samp_nodes, graph.years[samp_nodes]])
        .reshape(2, -1)
        .transpose()
    }
    """
        Sample Sub-Graph based on the connection of other nodes with currently sampled nodes
        We maintain budgets for each node type, indexed by <node_id, time>.
        Currently sampled nodes are stored in layer_data.
        After nodes are sampled, we construct the sampled adjacancy matrix.
    """
    layer_data = defaultdict(lambda: {})  # target_type  # {target_id: [ser, time]}

    budget = defaultdict(  # source_type
        lambda: defaultdict(lambda: [0.0, 0])  # source_id  # [sampled_score, time]
    )
    new_layer_adj = defaultdict(  # target_type
        lambda: defaultdict(  # source_type
            lambda: defaultdict(lambda: [])  # relation_type  # [target_id, source_id]
        )
    )
    """
        For each node being sampled, we find out all its neighborhood, 
        adding the degree count of these nodes in the budget.
        Note that there exist some nodes that have many neighborhoods
        (such as fields, venues), for those case, we only consider 
    """

    def add_budget(te, target_id, target_time, layer_data, budget):
        for source_type in te:
            tes = te[source_type]
            for relation_type in tes:
                if relation_type == "self" or target_id not in tes[relation_type]:
                    continue
                adl = tes[relation_type][target_id]
                if len(adl) < sampled_number:
                    sampled_ids = list(adl.keys())
                else:
                    sampled_ids = np.random.choice(
                        list(adl.keys()), sampled_number, replace=False
                    )
                for source_id in sampled_ids:
                    source_time = adl[source_id]
                    if source_time == None:
                        source_time = target_time
                    if source_id in layer_data[source_type]:
                        continue
                    budget[source_type][source_id][0] += 1.0 / len(sampled_ids)
                    budget[source_type][source_id][1] = source_time

    """
        First adding the sampled nodes then updating budget.
    """
    for _type in inp:
        for _id, _time in inp[_type]:
            layer_data[_type][_id] = [len(layer_data[_type]), _time]
    for _type in inp:
        te = graph.edge_list[_type]
        for _id, _time in inp[_type]:
            add_budget(te, _id, _time, layer_data, budget)
    """
        We recursively expand the sampled graph by sampled_depth.
        Each time we sample a fixed number of nodes for each budget,
        based on the accumulated degree.
    """
    for layer in range(sampled_depth):
        sts = list(budget.keys())
        for source_type in sts:
            te = graph.edge_list[source_type]
            keys = np.array(list(budget[source_type].keys()))
            if sampled_number > len(keys):
                """
                Directly sample all the nodes
                """
                sampled_ids = np.arange(len(keys))
            else:
                """
                Sample based on accumulated degree
                """
                score = np.array(list(budget[source_type].values()))[:, 0] ** 2
                score = score / np.sum(score)
                sampled_ids = np.random.choice(
                    len(score), sampled_number, p=score, replace=False
                )
            sampled_keys = keys[sampled_ids]
            """
                First adding the sampled nodes then updating budget.
            """
            for k in sampled_keys:
                layer_data[source_type][k] = [
                    len(layer_data[source_type]),
                    budget[source_type][k][1],
                ]
            for k in sampled_keys:

                add_budget(te, k, budget[source_type][k][1], layer_data, budget)
                budget[source_type].pop(k)

    """
        Prepare feature, time and adjacency matrix for the sampled graph
    """
    feature = {}
    times = {}
    indxs = {}
    texts = []
    for _type in layer_data:
        print(_type)
        if len(layer_data[_type]) == 0:
            continue
        idxs = np.array(list(layer_data[_type].keys()), dtype=np.int)
        print(idxs)
        tims = np.array(list(layer_data[_type].values()))[:, 1]
        if _type == "cell":
            idxs = idxs - gene_cell.shape[0]
        feature[_type] = graph.node_feature[_type][idxs]
        times[_type] = tims
        indxs[_type] = idxs

    edge_list = defaultdict(  # target_type
        lambda: defaultdict(  # source_type
            lambda: defaultdict(lambda: [])  # relation_type  # [target_id, source_id]
        )
    )
    for _type in layer_data:
        for _key in layer_data[_type]:
            _ser = layer_data[_type][_key][0]
            edge_list[_type][_type]["self"] += [[_ser, _ser]]
    """
        Reconstruct sampled adjacancy matrix by checking whether each
        link exist in the original graph
    """
    for target_type in graph.edge_list:
        te = graph.edge_list[target_type]
        tld = layer_data[target_type]
        for source_type in te:
            tes = te[source_type]
            sld = layer_data[source_type]
            for relation_type in tes:
                tesr = tes[relation_type]
                for target_key in tld:
                    if target_key not in tesr:
                        continue
                    target_ser = tld[target_key][0]
                    for source_key in tesr[target_key]:
                        """
                        Check whether each link (target_id, source_id) exist in original adjacancy matrix
                        """
                        if source_key in sld:
                            source_ser = sld[source_key][0]
                            edge_list[target_type][source_type][relation_type] += [
                                [target_ser, source_ser]
                            ]
    # print("feature",feature)

    return feature, times, edge_list, indxs, texts


graph.years = np.ones(len(target_nodes))

np.random.seed(seed)
jobs = []

for batch_id in np.arange(args.n_batch):
    # samp_nodes=np.random.choice(np.arange(gene_cell.shape[0]), batch_size, replace = False)

    p = sub_sample(
        np.random.choice(np.arange(gene_cell.shape[0]), args.batch_size, replace=False)
    )
    jobs.append(p)
debuginfoStr("Start Graph Autoencoder training")
patience = 10
early_stopping = EarlyStopping(patience, verbose=True)
# (node_feature-torch.mean(node_feature))/torch.std(node_feature)
for epoch in np.arange(500):
    L = 0
    for job in jobs:
        # print(job)
        feature = job[0]
        time = job[1]
        edge_list = job[2]
        indxs = job[3]
        node_dict = {}
        node_feature = []
        node_type = []
        node_time = []
        edge_index = []
        edge_type = []
        edge_time = []

        node_num = 0
        types = graph.get_types()
        for t in types:
            node_dict[t] = [node_num, len(node_dict)]
            node_num += len(feature[t])

        for t in types:
            node_feature += list(feature[t])
            node_time += list(time[t])
            node_type += [node_dict[t][1] for _ in range(len(feature[t]))]

        edge_dict = {e[2]: i for i, e in enumerate(graph.get_meta_graph())}
        edge_dict["self"] = len(edge_dict)

        for target_type in edge_list:
            for source_type in edge_list[target_type]:
                for relation_type in edge_list[target_type][source_type]:
                    for ii, (ti, si) in enumerate(
                        edge_list[target_type][source_type][relation_type]
                    ):
                        tid, sid = (
                            ti + node_dict[target_type][0],
                            si + node_dict[source_type][0],
                        )
                        edge_index += [[sid, tid]]
                        edge_type += [edge_dict[relation_type]]
                        """
                            Our time ranges from 1900 - 2020, largest span is 120.
                        """
                        edge_time += [node_time[tid] - node_time[sid] + 120]
        # node_feature = torch.stack(node_feature)
        # print(node_feature)

        # print(node_type)
        # print(edge_time)
        # print(edge_index)
        # print(edge_type)
        # print(idxs)
        node_feature = torch.stack(node_feature)
        node_type = torch.LongTensor(node_type)
        edge_time = torch.LongTensor(edge_time)
        edge_index = torch.LongTensor(edge_index).t()
        edge_type = torch.LongTensor(edge_type)

        # node_feature, node_type, edge_time, edge_index, edge_type, node_dict, edge_dict = \
        #    to_torch(feature, times, edge_list, graph)
        node_rep = gnn.forward(
            node_feature.to(device),
            node_type.to(device),
            edge_time.to(device),
            edge_index.to(device),
            edge_type.to(device),
        )
        gene_matrix = node_rep[
            node_type == 0,
        ]
        cell_matrix = node_rep[
            node_type == 1,
        ]
        # print(node_rep.shape)
        decoder = torch.mm(gene_matrix, cell_matrix.t())
        # loss_func = nn.MSELoss(size_average=True)
        adj = gene_cell[
            indxs["gene"],
        ]
        adj = adj[:, indxs["cell"]]
        # adj= torch.tensor(adj,dtype=torch.float32).to(device)
        # loss=loss_func(adj,decoder)
        adj = torch.tensor(adj, dtype=torch.float32).to(device)
        # loss=loss_func(adj,decoder)
        loss = F.kl_div(
            decoder.softmax(dim=-1).log(), adj.softmax(dim=-1), reduction="sum"
        )

        L += loss.item()
        loss.requires_grad_(True)
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(gnn.parameters(), 1)
        optimizer.step()
    early_stopping(L / (int(gene_cell.shape[0] / args.n_batch)), gnn)
    if early_stopping.early_stop:
        print("Early stopping")
        break
    # print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    print("{}, {}".format(epoch + 1, L / (int(gene_cell.shape[0] / args.n_batch))))
debuginfoStr("Graph Autoencoder training finished")
# state = {'model':gnn.state_dict(), 'optimizer':optimizer.state_dict(), 'epoch':epoch}
# torch.save(state,model_dir+file0)
debuginfoStr("load training model")
checkpoint = torch.load(model_dir + file0, map_location=lambda storage, loc: storage)
# state = torch.load(model_dir+file0, map_location=lambda storage, loc: storage)
device = torch.device("cpu")
model = GNN(
    conv_name="hgt",
    in_dim=args.in_dim,
    n_hid=args.n_hid,
    n_heads=args.n_heads,
    n_layers=args.n_layers,
    dropout=args.dropout,
    num_types=2,
    num_relations=2,
).to(device)
model.eval()
# if (gene_cell.shape[1]<000):
#    ba = 5000
# else:
ba = 2000
optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr)
node_feature = torch.cat((encoded, encoded2), 0)
# node_feature = (node_feature-torch.mean(node_feature))/torch.std(node_feature)
optimizer.load_state_dict(checkpoint["optimizer"])
model.load_state_dict(checkpoint["model"])
g_embedding = []
if gene_cell.shape[0] % ba == 0:
    upper = gene_cell.shape[0] / ba + 1
else:
    upper = gene_cell.shape[0] / ba + 2

with torch.no_grad():
    for i in range(1, int(upper)):
        if i != int(gene_cell.shape[0] / ba + 1):
            batch_gene = node_feature[
                (i - 1) * ba : ba * i,
            ]
            adj = gene_cell[(i - 1) * ba : ba * i, :]
        else:
            batch_gene = node_feature[(i - 1) * ba + 1 : gene_cell.shape[0] + 1, :]
            adj = gene_cell[(i - 1) * ba : ba * i, :]
        g = np.nonzero(adj)[0]
        c = np.nonzero(adj)[1] + adj.shape[0]
        edge1 = list(g)
        edge2 = list(c)
        edge_index = torch.tensor([edge1, edge2], dtype=torch.long)
        x = {
            "gene": torch.tensor(
                node_feature[(i - 1) * ba : ba * i, :], dtype=torch.float
            ),
            "cell": torch.tensor(encoded2, dtype=torch.float),
        }
        edge_index_dict = {
            ("gene", "g_c", "cell"): torch.tensor([g, c], dtype=torch.long)
        }
        edge_reltype = {("gene", "g_c", "cell"): torch.tensor([g, c]).shape[1]}
        num_nodes_dict = {"gene": adj.shape[0], "cell": gene_cell.shape[1]}
        data = Data(
            edge_index_dict=edge_index_dict,
            edge_reltype=edge_reltype,
            num_nodes_dict=num_nodes_dict,
            x=x,
        )
        a = np.nonzero(adj)[0]
        b = np.nonzero(adj)[1]
        node_type = list(np.zeros(adj.shape[0])) + list(np.ones(adj.shape[1]))
        node_type = torch.LongTensor(node_type)
        edge_index = data["edge_index_dict"][("gene", "g_c", "cell")]
        edge_type = list(np.zeros(len(edge_index[1])))
        edge_time = torch.LongTensor(list(np.zeros(len(edge_index[1]))))
        edge_type = torch.LongTensor(edge_type)
        node_rep = model.forward(
            (
                torch.cat(
                    (
                        batch_gene,
                        node_feature[
                            (gene_cell.shape[0]) : (node_feature.shape[0] + 1), :
                        ],
                    ),
                    0,
                )
            ).to(device),
            node_type.to(device),
            edge_time.to(device),
            edge_index.to(device),
            edge_type.to(device),
        )
        gene_matrix = node_rep[
            node_type == 0,
        ]
        cell_matrix = node_rep[
            node_type == 1,
        ]
        g_embedding.append(gene_matrix)
if gene_cell.shape[0] % ba == 0:
    gene_matrix = np.vstack(g_embedding[0 : int(gene_cell.shape[0] / ba)])
else:
    final_tensor = np.vstack(g_embedding[0 : int(gene_cell.shape[0] / ba)])
    gene_matrix = np.concatenate((final_tensor, gene_matrix), 0)
# final_tensor=np.vstack(g_embedding[0:int(gene_cell.shape[0]/ba)])
# print(final_tensor.shape)
# gene_matrix=np.concatenate((final_tensor,gene_matrix),0)
# print(gene_matrix.shape)
np.savetxt(gene_dir + file0, gene_matrix, delimiter=" ")
np.savetxt(cell_dir + file0, cell_matrix.detach().numpy(), delimiter=" ")
debuginfoStr(" finished")
