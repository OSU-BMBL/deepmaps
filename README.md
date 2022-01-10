# deepmaps

This is the repository for the manuscript: todo.

If you have any questions or feedback, please contact Qin Ma <qin.ma@osumc.edu>.

## Dev environment

```{bash}
python: 3.8.5
pytorch: 1.9.1
torch-geometric: 2.0.1
NVIDIA Driver Version: 450.102.04
CUDA Version: 11.0
GPU: 2x A100-PCIE-40GB
System: Red Hat Enterprise Linux release 8.3 (Ootpa)
```

## How to install

```{bash}
# Python 3.8.5
pip install -r requirements.txt

```

If there exists any problem in pytorch-genomic package install, please do as follows:

Check your torch version, python version and cuda version, download “torch_cluster.whl” , “torch_scatter.whl”, “torch_sparse.whl” and “torch_spline_conv.whl” from https://pytorch-geometric.com/whl/, then pip install \*.whl, and install other package by pip.

For example:
The torch version is 1.5.0, python version is 3.7, linux and cuda is 10.1:

1. Step1: click torch-1.5.0+cu101
2. Step2:

```
wget https://data.pyg.org/whl/torch-1.5.0%2Bcu101/torch_cluster-1.5.7-cp37-cp37m-linux_x86_64.whl
wget https://data.pyg.org/whl/torch-1.5.0%2Bcu101/torch_scatter-2.0.5-cp37-cp37m-linux_x86_64.whl
wget https://data.pyg.org/whl/torch-1.5.0%2Bcu101/torch_sparse-0.6.7-cp37-cp37m-linux_x86_64.whl
wget https://data.pyg.org/whl/torch-1.5.0%2Bcu101/torch_spline_conv-1.2.0-cp37-cp37m-linux_x86_64.whl
```

3. Step3:

```
pip install torch_cluster-1.5.7-cp37-cp37m-linux_x86_64.whl
pip install torch_scatter-2.0.5-cp37-cp37m-linux_x86_64.whl
pip install torch_sparse-0.6.7-cp37-cp37m-linux_x86_64.whl
pip install torch_spline_conv-1.2.0-cp37-cp37m-linux_x86_64.whl
Step4: test “import torch_geometric”
If lack other packages when you run the code, please ‘pip install package’ directly.

```
