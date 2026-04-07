################################################################################
# dataextract
################################################################################

# Takes a h5 file as input, extracts expression data and samples names, lets the
# user filter samples by tissue type (using sample titles) and produces a matrix.
# Writes the output as a parquet file.

################################################################################
# Environment Setup
################################################################################

setwd(dir = "./")

libs = (c('BiocManager', 
          'rhdf5', 
          'tidyr',
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
###############################################################################

out_file_name = c("demo") # name of the ouptut file

tissue = c("shoot") # key word or list of key words to extract tissue types

################################################################################
# Function to read in data
################################################################################

readin = function(h5file){

  contents <- h5ls(file = h5file) #contents of the h5 file
  
  samples = h5read(h5file, 
                   "meta/Sample_geo_accession") # sample IDs
  
  genes = h5read(h5file, 
                 "meta/genes") # gene list
  
  sample_title = h5read(h5file, 
                        "meta/Sample_title") # contains tissue info
  
  return(list(file_name = h5file,
              contents = contents, 
              samples = samples, 
              genes = genes, 
              sample_title = sample_title))
  }

################################################################################
# Function to extract expression data from h5file
################################################################################

tissuetype = function(h5file, 
                      tissues, 
                      sample_title, 
                      samples,
                      genes){
  
  matching_indices <- grep(paste(tissues, collapse = "|"), 
                        sample_title, 
                        ignore.case = TRUE)
  
  full_expre_data = t(h5read(h5file, 
                             "data/expression"))
  
  rownames(full_expre_data) = samples
  
  expre_data = t(h5read(h5file, 
                        "data/expression",
                        index=list(1:length(genes),
                                   matching_indices)))
  
  samp = as.character(h5read(h5file, 
                             "meta/Sample_geo_accession"))
  
  sample_locations = samples[matching_indices]
  
  colnames(expre_data) = genes
  
  rownames(expre_data) = samples[matching_indices]
  
  expre_data = t(expre_data)
  #expre_data = cbind(gene = genes, expre_data) # needed to keep ids for parquet file
  expre_data = bind_cols(gene = genes, expre_data) 

  warning(print(length(grep(paste(tissues, collapse = "|"), 
                            sample_title, 
                            ignore.case = TRUE))))
  
  return(list(expre_data = expre_data, 
              full_expre_data = full_expre_data, 
              sample_locations = sample_locations))
}

################################################################################
# Function to data check extraction
################################################################################

data_check = function(data_set) {
  
  lapply(data_set[["sample_locations"]], function(i) {
    
    output = data_set[["expre_data"]][,i] == data_set[["full_expre_data"]][i,]
    
    if (all(output)) {
    print("Data check successful")}
    else {
      warning(("Data check failed! Data extraction err"))
      data_set[["expre_data"]][,i] == data_set[["full_expre_data"]][i,]
    }
    }
  )
}

################################################################################
# Using functions
################################################################################

file = readin(h5file = "data/Arabidopsis_thaliana_genecount_v1.h5")

data_set = tissuetype(h5file       = file$file_name,
                      tissues      = tissue,                         # <-- change to
                      sample_title = file$sample_title,              # tissue type
                      samples      = file$samples,                   # of interest
                      genes        = file$genes)

check = data_check(data_set = data_set)

################################################################################
# Writing out file
################################################################################

write_parquet(as.data.frame(
  data_set[["expre_data"]]), paste0("data/", out_file_name, ".parquet"))
