#/bin/bash

cd /deepmaps/
mkdir cell
mkdir gene
mkdir model
mkdir att

cd /deepmaps/docker
mkdir data
cd data 

mkdir /deepmaps/lisa_output

if [ ! -e lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5 ] 
then 
	wget https://bmbl.bmi.osumc.edu/downloadFiles/deepmaps/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5
else 
    echo "skipping lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5"
fi

if [ ! -e lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv ] 
then 
	wget https://bmbl.bmi.osumc.edu/downloadFiles/deepmaps/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv.gz
    gunzip lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv.gz
else 
    echo "skipping lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv.gz"
fi

cd /deepmaps/
Rscript /deepmaps/docker/test.R
