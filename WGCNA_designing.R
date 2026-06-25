#####################################################
#SET UP AND LOADING PACKAGES
#####################################################

install.packages(c("tidyverse", "magrittr", "WGCNA"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("impute", force = TRUE)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore", force = TRUE)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("maaslin3", force = TRUE)

library(tidyverse)
library(magrittr)
library(WGCNA)
library(impute)
library(preprocessCore)
library(maaslin3)
library(dplyr)
library(DESeq2)

#point to the correct working directory
setwd("C:/Users/c21010903/OneDrive - Cardiff University/PhD/Gene network/GlenHTWA_network_mapping/AB")
getwd() #check that you are now in the correct working directory

#input data
data <-as.matrix(read.csv(
  file = paste ("AB_GlenHTWA_featcount_merge_WGCNA.csv", sep = ""),
  header = TRUE, 
  row.names = 1,
  sep = ","))

#Convert data from matrix to dataframe if needed(allows addition of GeneId rowname)
#Also sort out columns

dframe <- as.data.frame(data)

dframe["GeneId"] <- row.names(data)   #add in GeneId rowname
dframe1 <- dframe[, c(13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)] #moves GeneId column to first
dframe1

row.names(dframe1) <- NULL      #Removes greyed-out gene names that were duplicated rownames

data[1:5,1:10]      #Look at first 5 rows and 10 columns

#Tidy data to visualise any outliers
col_sel = names(dframe1)[-1] #get all but first column name
mdata <- dframe1 %>%
  tidyr::pivot_longer(
    .,
    col = all_of(col_sel)
  ) %>%
  mutate(
    group = gsub("-.*","", name) %>% gsub("[.].*","", .) #Get the shorter treatment names
  )

  mdata$group = factor(mdata$group,
                       levels = c("HT", "WA", "GA", "GD"))    #Order the groups in the plot

#Plot groups to visualise outgroups (sample groups vs RNA Seq counts)
(
  p <- mdata %>%
    ggplot(., aes(x = name, y = value)) +
    geom_violin() +
    geom_point(alpha = 0.2) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)
    ) +
    labs(x = "Treatment Groups", y = "RNA Seq Counts") +
    facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")
)

#Normalise the counts using DESeq2
ds_input = as.matrix(dframe1[,-1])
row.names(ds_input) = dframe1$GeneId
ds_input[1:5, 1:10]

meta_df <- data.frame(sample = names(dframe1[-1])) %>%
  mutate(
    Type = gsub("-.*","", sample) %>% gsub("[.].*","", .)
  )
  
  dds <- DESeqDataSetFromMatrix(round(ds_input),
                                meta_df,
                                design = ~Type)
  dds <- DESeq(dds)
  
  vsd <- varianceStabilizingTransformation(dds)
  
  #library(genefilter) 
  
  wpn_vsd <-getVarianceStabilizedData(dds)
  rv_wpn <- rowVars(wpn_vsd)
  summary(rv_wpn)
  
  q75dframe1 <- quantile( rowVars(wpn_vsd), .75)
  #q95dframe1 <- quantile(rowVars(wpn_vsd), .95) #use this if need 95% quantile expression
  expr_normalised <- wpn_vsd[ rv_wpn > q75dframe1, ]

#Visualise normalisation
  
expr_normalised[1:5, 1:10]

expr_normalised_df <- data.frame(expr_normalised) %>%
  mutate(
    Gene_Id = row.names(expr_normalised)
  ) %>%
  pivot_longer(-Gene_Id)

expr_normalised_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalised and 75 quantile expression",
    x = "treatment",
    y = "normalised expression"
  )
  
#For WGCNA, rows need to be treatment, and column = gene probes - if this is not the case then need to transpose the data
#Transpose the data

input_matrix = t(expr_normalised)
input_matrix[1:5,1:10] #should now see rows = treatments and columns = gene proves, so can start WGCNA

allowWGCNAThreads() #Allow multi-threading (optional)

#Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  input_matrix,
  #blocksize = 30,
  powerVector = powers,
  verbose = 5
)

#Visualisation
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale Independence")
)

text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)

abline(h = 0.9, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean Connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

#Select a soft threshold (power) near the curve of the plot, (ie. maybe 8, 9 or 10 for my AB GlenHTWA data?) but experiment with other powers to see how it affects your results. Can then use this to create the network using blockwiseModules

picked_power = 9 #make sure to adjust this based on output soft threshold graphs each time
temp_cor <- cor
cor <- WGCNA::cor

network <- blockwiseModules(input_matrix,
                            #Adjacency Function
                            power = picked_power,
                            networkType = "signed",
                            
                            #Tree and block options
                            deepSplit = 2,
                            pamRespectsDendro = F,
                            #detectCutHeight = 0.75,
                            minModuleSize = 30,
                            maxBlockSize = 4000,
                            
                            #Module Adjustments
                            reassignThreshold = 0,
                            mergeCutHeight = 0.25,
                            
                            #TOM - Archive the run results in TOM file (saves time)
                            saveTOMS = T,
                            saveTOMFileBase = "WGCNA_TOM",
                            
                            #Output Options
                            numericLabels = T,
                            verbose = 3)

cor <- temp_cor #Return cor function to original namespace

##Visualise the modules
#Convert labels to colours for plotting
mergedColors = labels2colors(network$colors)
#Plot the dendrogram and the module colours beneath
plotDendroAndColors(
  network$dendrograms[[1]],
  mergedColors[network$blockGenes[[1]]],
  "Module Colours",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

network$colors[network$blockGenes[[1]]]
table(network$colorrs)

#Relate Module (cluster) Assignments to Treatment Groups
module_df <- data.frame(
  GeneId = names(network$colors),
  colors = labels2colors(network$colors)
)
module_df[1:5,]

write_delim(module_df,
            file = "AB_WGCNA_gene_modules.txt",
            delim = "\t")

#^^This writes a tab delim file listing the genes and their modules. However, need to work out which modules are associated with each trait (ie. HT, WA, GA, GD). WGCNA will calculate an Eigengene (hypothetical central gene) for each module, so it is easier to determine if modules 

#Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_matrix, mergedColors)$eigengenes

#Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME", "", .)

#Add treatment names
MEs0$treatment = row.names(MEs0)

#Tidy and plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Cluster-Transcriptome Relationships", y = "Modules", fill = "corr") #Trait/treatment is the transcriptome here
#The clusters are basically the groups of potential genes of interest, it doesn't matter what they are at this stage, as can draw them out later I believe
#Red clusters (colour groups) show high correlation with the corresponding transcriptome so pick these out for the examination of expression profiles (next section)

#Examine Expression Profiles
#Pick out a few modules (clusters) of interest here, from looking at heat map
modules_of_interest = c("pink", "red", "lightgreen")

#Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$GeneId

#Get normalised expression for those genes
expr_normalised[1:5, 1:10]

subexpr = expr_normalised[submod$GeneId,]

submod_df = as.data.frame(subexpr) %>%
  mutate(
    GeneId = row.names(.)
  ) %>%
  pivot_longer(-GeneId) %>%
  mutate(
    module = module_df[GeneId,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=GeneId)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "transcriptome",
       y = "normalised expression")

##Generate and export networks
genes_of_interest = module_df %>% 
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalised[genes_of_interest$GeneId,]
expr_of_interest[1:5,1:5]

# Only recalculate TOM for modules of interest (faster, although there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

#Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)

#Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "AB_GlenHTWA_WGCNA_EdgeList.tsv",
            delim = "\t")

#The edgelist.txt exported is the complete correlation network for the chosen clusters (colours). The network still needs to be subsetted down (by weight or minimal spanning) to identify hub genes. Steps forward include identifying hub genes and gene ontology enrichment. The edgelist can be loaded into igraph (R package) or Cytoscape (standalone Java Application) for further network analysis.


