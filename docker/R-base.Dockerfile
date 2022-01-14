FROM satijalab/seurat:4.0.1
LABEL maintainer="Cankun Wang <cankun.wang@osumc.edu>"

WORKDIR /tmp

# Ubuntu dependency found at /rocker_scripts/install_tidyverse.sh
RUN NCPUS=${NCPUS:-1}
RUN set -e
RUN apt-get -q  update && apt-get -y --no-install-recommends install \
	libxml2-dev \
	libcairo2-dev \
	libgit2-dev \
	default-libmysqlclient-dev \
	libpq-dev \
	libsasl2-dev \
	libsqlite3-dev \
	libssh2-1-dev \
	unixodbc-dev \
	python \
	python3-pip

# Found more libraries need to be installed during my debugging
RUN apt-get -y --no-install-recommends install \
	libbz2-dev \
	liblzma-dev \
	libsodium-dev \
	libhiredis-dev

###############
# HTSlib 1.11.0#
# Folk from https://github.com/chrisamiller/docker-r-seurat/blob/master/Dockerfile
###############
#ENV HTSLIB_INSTALL_DIR=/opt/htslib

RUN wget --no-check-certificate https://github.com/samtools/htslib/archive/1.11.0.zip && \
	unzip 1.11.0.zip && \
	rm 1.11.0.zip && \
	cd /tmp/htslib-1.11.0 && \
	#./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
	make && \
	make install && \ 
	cp -R * /usr/lib/

# Install Bioconductor dependencies
RUN R -e 'BiocManager::install(c("GenomicAlignments", "ggbio", "biovizBase", "fgsea", "ComplexHeatmap", "karyoploteR", "MAGeCKFlute"))'

# Install Bioconductor databases
RUN R -e 'BiocManager::install(c("JASPAR2020", "GO.db", "EnsDb.Hsapiens.v86","EnsDb.Mmusculus.v79", "org.Hs.eg.db", "org.Mm.eg.db"))'

# Install CRAN dependencies
RUN install2.r --error --skipinstalled -r $CRAN \
	Polychrome \
	qs \
	plumber \
	vroom \
	lintr \
	gert \
	Signac \ 
	logger \
	tictoc \
	msigdbr \
	Gmisc \ 
	rematch \
	readr \
	openxlsx \
	readxl \
	sp \
	statmod \
	statmod\
	nloptr\
	minqa\
	lme4\
	rio\
	maptools\
	pbkrtest\
	carData\
	car\
	corrplot\
	broom\
	rstatix\
	polynom\
	ggsignif\
	ggsci\
	ggpubr \
	spatstat \
	rlist \
	redux \
	devtools \
	enrichR \
	NMF \
	ggalluvial  \
	svglite \
	expm \
	sna \
	gg.gap

RUN installGithub.r -d FALSE -u FALSE\
	immunogenomics/presto \ 
	liulab-dfci/MAESTRO \
	sqjin/CellChat 

# Socket-io 
RUN pip3 install python-socketio[client]==4.6.1

# Clean up installation
RUN rm -rf /tmp/* 
RUN rm -rf /var/lib/apt/lists/*

# Add Tini
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
#ENTRYPOINT ["/tini", "--"]
ENTRYPOINT ["/tini", "--"]
CMD ["Rscript", "-e", "installed.packages()"]

