################################################################################
# RNA-seq Read Count Normalisation
################################################################################

# Takes a parquet file containing expression data from extract_data.R as input  
# and normalises each library using the DESeq2 shrinkage estimate model. Also
# facilitates the extraction of genes by user-supplied .csv 

################################################################################
# Environment Setup
################################################################################

setwd(file.path("C:", "Users", "c21010903", "OneDrive - Cardiff University", "PhD", "Gene network", "Glen_Regulatory_Networks-master", "src"))
getwd() #check that you are now in the correct working directory

libs = (c('BiocManager', 
          'DESeq2', 
          'dplyr',
          'arrow'))

if (requireNamespace(libs, 
                     quietly = TRUE)) {
  install.packages(libs)
}

invisible(lapply(libs, 
                 library, 
                 character.only = TRUE))

################################################################################
# User-defined variables
################################################################################

#input_file = file.choose() # this can be the .csv merged file from feature counts
input_file = "AB_GlenHTWA_merged.csv" #needs to be a string not file - so will read into Values not Data, this is fineee
gene_list = "AB_gene_list.csv" # list of desired genes to extract, also same as above, string not file

output_file = "AB_countnormalised" # name of the generated normalised output file

norm_input_file = "AB_countnormalised.parquet" # name of norm file if just for gene 
                                         # extraction

################################################################################
# Function to normalise reads
################################################################################

normalise = function(input_file, output_file) {
  
  norm_func = function(expre_data) {
    
    countData = as.data.frame(expre_data)
    
    rownames(countData) = countData$GeneID
  
    countData$GeneID = NULL
    
    colData = data.frame(
      condition = rep("unknown", ncol(countData)),
      row.names = colnames(countData))
    
    dds = DESeqDataSetFromMatrix(countData = countData, 
                                 colData = colData, 
                                 design = ~ 1) # Design ~ 1 means no contrast
    
    dds = estimateSizeFactors(dds)
    
    sizeFactors(dds)
    
    normalised_counts = counts(dds, 
                               normalized = TRUE)
    
    normalised_counts = as.data.frame(normalised_counts)
    
    normalised_counts$id = as.data.frame(expre_data)$GeneID
    
    normalised_counts = normalised_counts %>% select(id, everything())
    
    return(normalised_counts)
}
  
  if (exists("data_set") == TRUE) {
    normalised_counts = norm_func(data_set[["expre_data"]])
    }
  else {
    expre_data = read.csv(input_file)
    normalised_counts = norm_func(expre_data)
    }
  
  write.csv(as.data.frame(
    normalised_counts), paste0(output_file, ".csv"))
  
  return(normalised_counts)
}

################################################################################
# Function to extract genes of interest
################################################################################

gene_extract = function(gene_list) {
  
  gene_list = read.csv("AB_gene_list.csv")

  extract = function(normalised_counts, gene_list) {
    
    GRN_genes = as.data.frame(normalised_counts) %>%
      #tibble::rownames_to_column(., var = "id") %>%
      filter(id %in% gene_list$id)
    
}

  if (exists("normalised_counts") == TRUE) {
    extract(normalised_counts, gene_list)
    } 
  else {
    normalised_counts = read.csv(norm_input_file)
    extract(normalised_counts, gene_list)}
}
  

################################################################################
# Using functions
################################################################################

normalised_counts = normalise(input_file, output_file)

extracted = gene_extract(gene_list)

################################################################################
# Write out file
################################################################################

write_parquet(as.data.frame(
  extracted), paste0(gene_list, ".parquet"))

  
