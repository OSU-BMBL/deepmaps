# DeepMAPS

This is the repository for the manuscript: [Single-cell biological network inference using a heterogeneous graph transformer](https://www.nature.com/articles/s41467-023-36559-0).

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

## Preparations

### Example data

We used a single-cell multiome ATAC+Gene expression dataset from [10X Genomics](https://www.10xgenomics.com/resources/datasets/fresh-frozen-lymph-node-with-b-cell-lymphoma-14-k-sorted-nuclei-1-standard-2-0-0). The raw data is derived from 14,566 cells diagnosed with diffuse small lymphocytic lymphoma (DSLL) of the lymph node lymph.

- [RNA+ATAC count matrix (.h5) (118 MB)](https://bmblx.bmi.osumc.edu/downloadFiles/deepmaps/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5)
- [ATAC fragments (.tsv.gz) (2.7 GB)](https://bmblx.bmi.osumc.edu/downloadFiles/deepmaps/lymph_node_lymphoma_14k_atac_fragments.tsv.gz)
- [ATAC fragments index (.tbi) (1 MB)](https://bmblx.bmi.osumc.edu/downloadFiles/deepmaps/lymph_node_lymphoma_14k_atac_fragments.tsv.gz.tbi)
- [RNA velocity matrix (.csv.gz) (434 MB)](https://bmblx.bmi.osumc.edu/downloadFiles/deepmaps/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv.gz)

### Manual installation

- python: 3.8
- pytorch: 1.9.0
- cuda: 10.2
- torch_geometric: 2.0.3

```{bash}
conda create -n deepmaps_env python=3.8.5
conda activate deepmaps_env
conda install pytorch==1.9.0 torchvision==0.10.0 torchaudio==0.9.0 cudatoolkit=10.2 -c pytorch
conda install pyg -c pyg -c conda-forge
pip install kneed==0.7.0
pip install seaborn==0.11.1
pip install dill==0.3.3
```

### Docker

The DeepMAPS docker image and tutorial can be found here: https://github.com/OSU-BMBL/deepmaps/tree/master/docker

### Troubleshooting

If there exists any problem in pytorch-genomic package install, please do as follows:

Check your torch version, python version and cuda version, download “torch_cluster.whl” , “torch_scatter.whl”, “torch_sparse.whl” and “torch_spline_conv.whl” from https://pytorch-geometric.com/whl/, then pip install \*.whl, and install other package by pip.

Check your torch version, python version and cuda version,

First, download the following packages from https://pytorch-geometric.com/whl/

1. torch_cluster.whl
2. torch_scatter.whl
3. torch_sparse.whl
4. torch_spline_conv.whl

then go to the download directory and `pip install \*.whl`

For example:
If your torch version is 1.5.0, python version is 3.7, linux and cuda is 10.1:

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

```

4. Step4: test if packages are installed

```
python -c "import torch_geometric"
```

If lack other packages when you are running the code, please run `pip install [package NAME]` directly.
