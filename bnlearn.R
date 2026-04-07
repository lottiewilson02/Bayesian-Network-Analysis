################################################################################
# bnlearn
################################################################################

# Takes genes expression data extracted from normalised RNA-seq countdata and a 
# gene list as input, discretises the data, and performs hill-climbing or TABU 
# Bayesian network learning algorithms to said data. Outputs a BN model as csv.

################################################################################
# Environment Setup
################################################################################

setwd(file.path("C:", "Users", "c21010903", "OneDrive - Cardiff University", "PhD", "Gene network", "Glen_Regulatory_Networks-master", "src"))
getwd() #check that you are now in the correct working directory


url <- "https://cran.r-project.org/src/contrib/Archive/Rgraphviz/Rgraphviz_1.18.1.tar.gz" #Need to install Rgraphviz manually from archive, as has been discontinued
  pkgFile <- "Rgraphviz1.18.1.tar.gz"
  download.file(url = url, destfile = pkgFile)

libs = (c('BiocManager', 
          'bnlearn',
          'dplyr',
          'arrow',
          'Rgraphviz',
          'gRain'))


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

extracted_file = c("MJ_gene_extracted.parquet")

outfile = c("MJ_gene_list") #this may be the wrong file

gene_list = read.csv("MJ_gene_list.csv")

################################################################################
# Function to discretise data
################################################################################

discretise = function(gene_list) {
  
  discretise_data = function(extracted) {
    
    extracted_t = as.data.frame(t(extracted[,-1])) # df needs to be transposed
    colnames(extracted_t) = extracted$id           #transfer gene ids
    
    disc_expre_data = discretize(extracted_t, 
                                 method = 'hartemink', 
                                 breaks = 3,         # number of bins
                                 ibreaks = 10)
    
    colnames(disc_expre_data)[match(gene_list$id, 
                                    colnames(disc_expre_data))] = gene_list$X
    
    return(disc_expre_data)
  }
  
  if (exists("extracted")) {
    result = discretise_data(extracted)
  } else {
    extracted = read_parquet(paste0(extracted_file))
    result = discretise_data(extracted)
    }

  return(result)
}


################################################################################
# Generating Bayesian Networks
################################################################################

# Hill-climbing method (fast and low mem)

hillclimbing = function(outfile, discretize_data) {
  
  net_hc = hc(discretize_data, 
           restart = 10000,
           max.iter = 10000,
           maxp = 5, 
           whitelist = NULL,
           blacklist = NULL,
           score = NULL,
           start = NULL)
  
  fittedbn = bn.fit(net_hc, discretize_data)
  
  graphviz.chart(fittedbn, type = "barprob", layout = "neato", scale = c(0.75, 0.9), grid = TRUE)
  
  plot(graphviz.plot(net_hc, layout = "dot", shape = "circle", fontsize = 20))
  
  
  arcs_df = as.data.frame(arcs(net_hc))
  
  arcs_df$interaction <- "regulation" # user needs to look at the graph to infer
  
  colnames(arcs_df) = c("from", "to", "type")
  
  write.csv(arcs_df, 
            paste0(outfile,"_hc.csv"),
            row.names=FALSE, 
            quote=FALSE)
  
  return(net_hc)

}


# Tabu search method (greedy)

tabu_search = function(outfile, discretize_data) {

  net_tabu = tabu(discretize_data, 
                  max.iter = Inf,
                  maxp = 10,  
                  tabu = 50, 
                  max.tabu = 50,
                  whitelist = NULL,
                  blacklist = NULL,
                  score = NULL,
                  start = NULL)
  
  plot(graphviz.plot(net_tabu, layout = "dot", shape = "circle", fontsize = 10))
  
  fittedbn = bn.fit(net_tabu, discretize_data)
  
  graphviz.chart(fittedbn, type = "barprob", layout = "neato", scale = c(0.75, 0.9), grid = TRUE)
  
  arcs_df = as.data.frame(arcs(net_tabu))
  
  arcs_df$interaction <- "regulation" # user needs to look at the graph to infer
  
  colnames(arcs_df) = c("from", "to", "type")
  
  write.csv(arcs_df, 
            paste0(outfile,"_tabu.csv"),
            row.names=FALSE, 
            quote=FALSE)
  
  return(net_tabu)

}



################################################################################
# Using functions
################################################################################

discretize_data = discretise(gene_list)

hc_net = hillclimbing(outfile, discretize_data)

tabu_net = tabu_search(outfile, discretize_data)

