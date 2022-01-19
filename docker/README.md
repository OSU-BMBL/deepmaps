# DeepMAPS-docker

This base image contains all Python-related packages for the HGT model. Including PyTorch, PyTorch Geometric, Velocity, Lisa2, etc. This base image also contains all necessary for the R package. Including plumber, Seurat, Signac, tidyverse, Bioconductor suite (GenomicRanges, SingleCellExperiment, etc.)

## Requirements

To use this image you must have Docker Engine installed. Instructions for setting up Docker Engine are [available on the Docker website](https://docs.docker.com/engine/installation/).

## Prebuilt image

An prebuilt image (~16.3GB) is available on Docker Hub under the name [osubmbl/deepmaps-base](https://hub.docker.com/r/osubmbl/deepmaps-base).

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

Note that [a database (~9GB)](https://github.com/liulab-dfci/lisa2/tree/master/docs) from [LISA](https://github.com/liulab-dfci/lisa2) and [JASPAR TFBS database](https://github.com/OSU-BMBL/deepmaps/tree/master/jaspar) are required when you run the DeepMAPS's gene regulatory network workflow. The database was not packed into Docker image to reduce the size. Here we will download the database manually from cistrome.org to `[your_lisa_path]`:

```bash
$ LISA_PATH=[your_lisa_path]
$ cd $LISA_PATH
$ wget http://cistrome.org/~alynch/data/lisa_data/hg38_1000_2.0.h5
$ wget https://bmblx.bmi.osumc.edu/downloadFiles/deepmaps/jaspar_hg38_500.qsave
```

Start with an interactive bash (Change `[your_deepmaps_repository_path]` to your cloned DeepMAPS directory). The database file will be linked to LISA's installation directory:

```bash
$ WORK_DIR=[your_deepmaps_repository_path]
$ docker run --rm -it --init \
  -v $WORK_DIR:/deepmaps \
  -v $LISA_PATH:/home/user/miniconda/lib/python3.8/site-packages/lisa/data \
  --gpus=all \
  --ipc=host \
  --network=host \
  osubmbl/deepmaps-base bash
```

## Test workflow

Run the following testing script if you would like to test [DeepMAPS scRNA+ATACseq workflow](https://github.com/OSU-BMBL/deepmaps/blob/master/scRNA_scATAC_analyses_tutorial.html):

```bash
$ bash /deepmaps/docker/test.sh > /deepmaps/test_output.txt
```

The script will:

1. Download an example dataset
2. Load and process example data
3. Run HGT model and clustering
4. Generate cell-type-specific regulons

You may check the visualization from `/deepmaps/plot.png` and output logs when the testing script complete.

Also, you can start with a jupyter notebook if you would prefer working on an interactive UI:

```bash
$ docker run -p 8888:8888 \
    -v $WORK_DIR:/deepmaps \
    -v $LISA_PATH:/home/user/miniconda/lib/python3.8/site-packages/lisa/data \
    --gpus=all \
    --network=host \
    --ipc=host  osubmbl/deepmaps-base \
    jupyter notebook --allow-root --ip 0.0.0.0
```
