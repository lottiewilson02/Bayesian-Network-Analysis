# Bayesian-Network-Analysis
Scripts used for bayesian network analysis
Uses output data from featureCounts from RNA-seq pipeline (DESeq2 not required as is incorporated into these scripts) 
Produces GRNs using both Tabu and Hillclimbing algorithms. 
If featureCounts has been run on data through RNA-seq pipeline, then extract_data.R script is not needed. 

Script order - Extract_data -> countNormalisation -> bnlearn
