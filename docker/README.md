# DeepMAPS-docker

This base image contains all Python-related packages for the HGT model. Including PyTorch, PyTorch Geometric, Velocity, Lisa2, etc. This base image also contains all necessary for the R package. Including plumber, Seurat, Signac, tidyverse, Bioconductor suite (GenomicRanges, SingleCellExperiment, etc.)

## Requirements

To use this image you must have Docker Engine installed. Instructions for setting up Docker Engine are [available on the Docker website](https://docs.docker.com/engine/installation/).

## Prebuilt image

An prebuilt image (~16.2GB) is available on Docker Hub under the name [osubmbl/deepmaps-base](https://hub.docker.com/r/osubmbl/deepmaps-base).

For example, you can pull the image (recommended) using:

```bash
$ docker pull osubmbl/deepmaps-base
```

or you can build a docker image yourself from:

```bash
$ docker build -f Deepmaps-base.Dockerfile -t osubmbl/deepmaps-base .
```

## Usage

It is possible to run DeepMAPS programs inside a container using the shell command. We provide a guide to run on example data.

First, clone the DeepMAPS repository to create a local copy on your computer:

```bash
$ git clone https://github.com/OSU-BMBL/deepmaps.git
```

Note that a database (~9GB) needs to be downloaded from cistrome.org when you run the LISA. Occasionally, a user may not be able to connect to cistrome.org from their institutional server due to some security measure. To circumvent this, one can manually install the data required to run LISA.

First, use the command below to issue the download URL for the required dataset.

```bash
$ lisa download hg38 oneshot --url
http://cistrome.org/~alynch/data/lisa_data/hg38_1000_2.0.h5
```

Then on your machine, download LISA's required data from cistrome.org to `[your_lisa_path]`.

```bash
LISA_PATH=[your_lisa_path]
cd LISA_PATH
wget http://cistrome.org/~alynch/data/lisa_data/hg38_1000_2.0.h5
```

Starting with an interactive bash (Change `[your_deepmaps_repository_path]` to your cloned DeepMAPS directory):

```bash
$ WORK_DIR=[your_deepmaps_repository_path]
$ docker run --rm -it --init \
  -v $WORK_DIR:/deepmaps \
  -v $LISA_PATH:/lisa_data \
  --gpus=all \
  --ipc=host \
  --network=host \
  osubmbl/deepmaps-base bash
```

The last step is to install the data in the docker package and install it to LISA's package directory:
```bash
$ lisa install hg38 oneshot /lisa_data/hg38_1000_2.0.h5
```

The LISA site package folder should now contain a directory called data with the downloaded dataset inside:

```bash
$ ls /home/user/miniconda/lib/python3.8/site-packages/lisa/data
```

Then run the test script:

```bash
$ bash /deepmaps/docker/test.sh
```

The script downloads example data and runs DeepMAPS. DeepMAPS applies Louvain a graph-based model to cluster cells based on the cell feature reduction matrix which returns from the HGT model, you can see the visualization from `/deepmaps/plot.png`

Also, you can start with a jupyter notebook:

```bash
$ docker run -p 8888:8888 \
    -v $WORK_DIR:/deepmaps \
    -v $LISA_PATH:/lisa_data \
    --gpus=all \
    --network=host \
    --ipc=host  osubmbl/deepmaps-base \
    jupyter notebook --allow-root --ip 0.0.0.0
```
