cd /tmp/deepmaps/
mkdir cell
mkdir gene
mkdir model
mkdir att

cd /tmp/deepmaps/docker
mkdir data
cd data 

wget https://bmbl.bmi.osumc.edu/downloadFiles/deepmaps/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5

wget https://bmbl.bmi.osumc.edu/downloadFiles/deepmaps/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv.gz

gunzip lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv.gz

cd /tmp/deepmaps/
Rscript /tmp/deepmaps/docker/test.R
