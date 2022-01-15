# DeepMAPS-docker

## Python Base image

This base image contains all necessary dependencies for DeepMAPS. Including PyTorch (1.7.0, CUDA11, GPU version), PyTorch Geometric, Velocity, Lisa2, etc.

You can pull from Docker Hub
```{bash, eval=FALSE}
docker pull osubmbl/deepmaps-python-base
```
or build docker image yourself from `Dockerfile` 
```{bash, eval=FALSE}
docker build -f Python-base.Dockerfile -t osubmbl/deepmaps-python-base .
```
Test what packages are installed
```{bash, eval=FALSE}
docker run osubmbl/deepmaps-python-base
```
Start Jupyter notebook
```{bash, eval=FALSE}
docker run -p 8888:8888 --gpus=all --ipc=host osubmbl/deepmaps-python-base
docker run -p 8888:8888 --gpus=all --ipc=host osubmbl/deepmaps-python-base jupyter notebook --allow-root --ip 0.0.0.0
```

## R Base image

To build the docker image, enter project root directory first.

This base image contains all necessary for the package. Including plumber, Seurat, Signac, tidyverse, BioConductor suite (GenomicRanges, SingleCellExperiment, etc.)

You can pull from Docker Hub
```{bash, eval=FALSE}
docker pull osubmbl/deepmaps-r-base
```
or build docker image yourself from `Dockerfile` 
```{bash, eval=FALSE}
docker build -f R-base.Dockerfile -t osubmbl/deepmaps-r-base .
```
Test what packages are installed
```{bash, eval=FALSE}
docker run osubmbl/deepmaps-r-base
```

### Other images

```{bash, eval=FALSE}
# whoami
docker run -d -P -p 9000:80 --restart unless-stopped containous/whoami

```
