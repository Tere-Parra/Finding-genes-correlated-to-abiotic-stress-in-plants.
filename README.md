# Differential network analysis to find gene clusters correlated to abiotic stress in plants.

## Introduction

Weighted Correlation Network Analysis (WGCNA) can be used for finding clusters (modules) of highly correlated genes, summarizing such groups using the module eigengene or intramodular hub genes, for relating modules to one another and external sample traits, and calculating module membership measures. It is a type of analysis based on correlation coefficients between different samples. The correlation coefficient is the closest relationship that presents two variables, this is measured in Pearson, Spearman, and Kendal (Chengcheng et al., 2022; Sahai et al., 2018). It is used to evaluate gene expression patterns, where those most similar to each other are grouped into a module. This relationship also allows associating characteristics, as genes concentrated in a module may be involved in the same biological process or signaling pathway.

Here we use WGCNA to find correlated genes in multiple abiotic stress in plants (Cleopatra Mandarin and Carrizo citrange). We test seven different stress treatments in citrus plants:

1) Water stress (WS)
   
2) High irradiance/light (HL)
   
3) Heat stress (HS)
   
4) High irradiance and heat stress combination (HL + HS)
   
5) Water stress and high irradiance combination (WS + HL)
   
6) Water stress and heat stress combination (WS + HS)
    
7) The combination of water stress, high light, and heat stress (WS +HL +HS).

8) Control (CT)

Therefore, it is reasonable to expect that genes that cluster together in co-expressed modules within the WGCNA analysis would have some similarity in their function or regulation, as they are similarly expressed under these stress conditions. However, due to the conditions, it is expected to find gene modules that overlap to some extent, but there are also significant differences in gene expression between conditions. Also, it is important to consider both gene modules that are common to various conditions and those that are specific to the combinations.

In summary, it is likely to find co-expressed genes in the different abiotic stress conditions, which would suggest similarities in the gene response to these conditions. However, accurate interpretation of the results will require further analysis and consideration of the underlying biology of the experimental conditions being studied.

 ## About the data 
The count data and metadata were taken from: Damián Balfagón, Zandalinas, S. I., dos, T., Claudete Santa‐Catarina, & Gómez‐Cadenas, A. (2022). Reduction of heat stress pressure and activation of photosystem II repairing system are crucial for citrus tolerance to multiple abiotic stress combinations. Physiologia Plantarum, 174(6). https://doi.org/10.1111/ppl.13809 and from the Gene Expression Omnibus Data (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203331) 

The data contains 48 samples with 6 biological replicates for each stress condition. You can find the GDC data here (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203331)

## Methods
First and foremost, we are going to use the following libraries:

''' R
library(WGCNA);

library(tidyverse);

library(CorLevelPlot);

library(gridExtra);

library(DESeq2);

library(apeglm);

library(ggplot2);  

library(dynamicTreeCut);  

library(flashClust);  

library(lattice);  

library(survival); 

library(Hmisc); 

library(vsn);

library(GEOquery);


# Allowing only two threads (consult WGCNA packages for more information)

allowWGCNAThreads(2);
ALLOW_WGCNA_THREADS = 2;  

# Initial variables 
options(stringsAsFactors = FALSE);    # indicates whether strings in a data frame should be treated as factor variables or as just plain strings
enableWGCNAThreads(2)

#############################################################################################################################################################
#####################                  DATA MODIFICATION                   ##################################################################################
#############################################################################################################################################################

data <- read.table("counts.txt", header=TRUE, sep="\t")

#Delete columns we will not use
colnames(data)

data <- filter( data, Type=="protein_coding")
#12 658 coding genes

data <- select(data, -Feature_GID, -Feature_TID, -Type, -Gene_Symbol,
               -Gene_Synonym, -Protein_ID, -Product)

data <- data[, 1:49]



#Now we have the count data

########### Obtain meta-data ####################

#Load GEO Data Query metadata 

query <-getGEO("GSE203331", GSEMatrix = TRUE)

#Obtain matrix metadata
metadata <-pData(phenoData(query[[1]]))

#Let's Change samples names
colnames(metadata)

#1) Select samples name, geo_accession and treatment information
metadata <- select(metadata, c("title", "treatment:ch1"))


# 2) change samples name for data sample names

#assign data samples names to nuevos_nombres

nuevos_nombres <- as.vector(colnames(data)[2:49])


# To the metadata$title names 
metadata$title <- nuevos_nombres

#Apply the change
metadata$title <- colnames(data)[2:49]


#Test the correct names position 
print(nuevos_nombres)

names <- as.data.frame(nuevos_nombres)

phenodata <- data.frame(
  treatment = c("Control", "Control", "Control", "Control", "Control", "Control", 
                "High light and heat stress combination", "High light and heat stress combination", 
                "High light and heat stress combination", "High light and heat stress combination", 
                "High light", "High light", "High light", "High light", "High light", "High light", 
                "High light and heat stress combination", "High light and heat stress combination", 
                "Heat stress", "Heat stress", "Heat stress", "Heat stress", "Heat stress", "Heat stress", 
                "Water stress and high light combination", "Water stress and high light combination", 
                "Water stress and high light combination", "Water stress and high light combination", 
                "Water stress and high light combination", 
                "Water stress and heat stress combination", "Water stress and heat stress combination", 
                "Water stress and heat stress combination", "Water stress and heat stress combination", 
                "Water stress and heat stress combination", "Water stress and heat stress combination", 
                "Water stress", "Water stress", "Water stress", "Water stress", "Water stress", "Water stress", 
                "Water stress and high light combination", 
                "Water stress, high light and heat stress combination", 
                "Water stress, high light and heat stress combination", 
                "Water stress, high light and heat stress combination", 
                "Water stress, high light and heat stress combination", 
                "Water stress, high light and heat stress combination", 
                "Water stress, high light and heat stress combination")
  
  
)

# some samples name dont match very well 
metadata <- cbind(names, metadata, phenodata)
colnames(metadata)

#Now we have the same samples names
new_data <- select(metadata, "nuevos_nombres", "treatment")

#Now let's prepare the binary metadata 


traits <- new_data %>% 
  mutate(treatment = ifelse(grepl('Control', treatment), 1, 0))  
  

new_data$treatment <- factor(new_data$treatment , 
                                              levels = c("Control",
                                                        "High light and heat stress combination",
                                                        "High light",
                                                        "Heat stress",
                                                        "Water stress and high light combination",
                                                        "Water stress and heat stress combination",
                                                        "Water stress",
                                                        "Water stress, high light and heat stress combination"))

severity.out_1 <- binarizeCategoricalColumns(new_data$treatment,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)

#Join both data frames
PhenoData <- cbind( traits, severity.out_1)

colnames(PhenoData)

Pheno <- rownames_to_column(PhenoData)
PhenoData <- Pheno[,2:9]


#PhenoData <- PhenoData %>%
           #column_to_rownames(var='nuevos_nombres')

save(data, PhenoData,file="matrix_data.RData")    
load("matrix_data.RData")




''' 




‌
