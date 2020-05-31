### Install QIIME_1.9
## If 'conda create -n qiime1 qiime -c bioconda' throws an error, try:
# conda create -n qiime1 python=2.7 numpy scipy matplotlib pandas biom-format qcli pyqi -c bioconda
# conda activate qiime1
# pip install qiime

## Put the 'FastTree' and 'uclust' in your $PATH

# Run QIIME:
pick_closed_reference_otus.py -i $PWD/ALL_SeqPrep-merged_no-PhiX.fna -o $PWD/qiime_pick_closed_reference -a -O 6 -p $PWD/QIIME.param

# Summarize the biom table
biom summarize-table -i $PWD/qiime_pick_closed_reference/otu_table.biom -o $PWD/qiime_pick_closed_reference/otu_table.summary.txt

#beta diversity:
core_diversity_analyses.py -a -O 6 -i $PWD/qiime_pick_closed_reference/otu_table.biom -o $PWD/core_diversity -m $PWD/map.txt -e 393514 -c "Treatment,Smoker,Region,BCG,Sex,Age,Sequencing" -t '~/anaconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/trees/97_otus.tree'

#Jackknifed_Beta_Diversity_and_Hierarchical_Clustering:
jackknifed_beta_diversity.py -i $PWD/qiime_pick_closed_reference/otu_table.biom -o $PWD/Jackknifed_Beta_Diversity_and_Hierarchical_Clustering -e 393514 -m $PWD/qiime_pick_closed_reference/otu_table.biom/map.txt -t '~/anaconda3/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/trees/97_otus.tree'
