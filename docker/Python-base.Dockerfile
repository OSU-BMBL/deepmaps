# runtime:2.3GB, devel: 5.27GB
FROM anibali/pytorch:1.7.0-cuda11.0
LABEL maintainer="Cankun Wang <cankun.wang@osumc.edu>"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH ${PATH}:/opt
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/opt/conda/lib/
#RUN pip install -U tensorflow tensorboard tensorflow-estimator tensorflow-gpu tensorflow-probability torch torchvision keras
# Copy a script that we will use to correct permissions after running certain commands

USER root

RUN apt-get -q update \
	&& apt-get install -y --no-install-recommends \
	wget \
	ca-certificates \
	sudo \
	make \
	build-essential \
	&& apt-get clean && rm -rf /var/lib/apt/lists/*

# Python
RUN conda update -y conda \
	&& conda install -y -c anaconda pip \
	#&& pip install -U pip \
	&& conda config --add channels conda-forge

RUN conda install numpy numpy_groupies scipy matplotlib pandas seaborn scikit-learn notebook dash plotly black bokeh h5py click jupyter jupyterlab pytables \
	&& pip install -U --no-cache-dir jupyterthemes jupyter_contrib_nbextensions python-igraph umap-learn numba Cython transformers pyreadr dill redis\
	&& jupyter notebook --generate-config \
	&& jupyter lab clean

# Pytorch Geometric
RUN pip install torch-scatter -f https://pytorch-geometric.com/whl/torch-1.7.0+cu110.html
RUN pip install torch-sparse -f https://pytorch-geometric.com/whl/torch-1.7.0+cu110.html
RUN pip install torch-cluster -f https://pytorch-geometric.com/whl/torch-1.7.0+cu110.html
RUN pip install torch-spline-conv -f https://pytorch-geometric.com/whl/torch-1.7.0+cu110.html
RUN pip install torch-geometric

# Velocity
RUN pip install -U --no-cache-dir \
	velocyto \
	scvelo

# LISA2
RUN conda install -c liulab-dfci lisa2

# socket-io must be v4.6
RUN pip install -U --no-cache-dir python-socketio==4.6.1

# Clean up
RUN conda clean --all -f -y

# Add Tini
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

EXPOSE 8888
ENTRYPOINT ["/tini", "--"]
CMD ["conda", "list"]
#CMD ["/bin/bash"]
