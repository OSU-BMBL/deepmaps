# DeepMAPS-docker

## Python Base image

This base image contains all necessary dependencies for DeepMAPS. Including PyTorch (1.7.0, CUDA11, GPU version), PyTorch Geometric, Velocity, Lisa2, etc.

```{bash, eval=FALSE}
# Build
docker build -f Python-base.Dockerfile -t wangcankun100/deepmaps-python-base .

# Test
docker run wangcankun100/deepmaps-python-base

# Start Jupyter notebook
docker run -p 8888:8888 --gpus=all --ipc=host wangcankun100/deepmaps-python-base
docker run -p 8888:8888 --gpus=all --ipc=host wangcankun100/deepmaps-python-base jupyter notebook --allow-root --ip 0.0.0.0

# Push
docker push wangcankun100/deepmaps-python-base
```

## R Base image

To build the docker image, enter project root directory first.

This base image contains all necessary for the package. Including plumber, Seurat, Signac, tidyverse, BioConductor suite (GenomicRanges, SingleCellExperiment, etc.)

```{bash, eval=FALSE}
# Build
docker build --progress=plain -f R-base.Dockerfile -t wangcankun100/deepmaps-r-base .

# Test what packages are installed
docker run wangcankun100/deepmaps-r-base

# Push
docker push wangcankun100/deepmaps-r-base
```

### Other images

```{bash, eval=FALSE}
# whoami
docker run -d -P -p 9000:80 --restart unless-stopped containous/whoami

```
