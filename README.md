# Differential network analysis to find gene clusters correlated to abiotic stress in plants.

## Introduction
Weighted Correlation Network Analysis (WGCNA) can be used for finding clusters (modules) of highly correlated genes, summarizing such groups using the module eigengene or intramodular hub genes, for relating modules to one another and external sample traits, and calculating module membership measures. It is a type of analysis based on correlating coefficients between different samples. The correlation coefficient is the closest relationship that presents two variables, this is measured in Pearson, Spearman, and Kendal (Chengcheng et al., 2022; Sahai et al., 2018). It is used to evaluate gene expression patterns, where those most like each other are grouped into a module. This relationship also allows associating characteristics, as genes concentrated in a module may be involved in the same biological process or signaling pathway.
Here we use WGCNA to find correlated genes in multiple abiotic stress in plants (Cleopatra Mandarin and Carrizo citrange). We test seven different stress treatments in citrus plants:

1) Water stress (WS)
   
2) High irradiance/light (HL)
   
3) Heat stress (HS)
   
4) High irradiance and heat stress combination (HL + HS)
   
5) Water stress and high irradiance combination (WS + HL)
   
6) Water stress and heat stress combination (WS + HS)
    
7) The combination of water stress, high light, and heat stress (WS +HL +HS).

8) Control (CT)

Therefore, it is reasonable to expect that genes that cluster together in co-expressed modules within the WGCNA analysis would have some similarity in their function or regulation, as they are similarly expressed under these stress conditions. However, due to the conditions, it is expected to find gene modules that overlap to some extent, but there are also significant differences in gene expression between conditions. Also, it is important to consider both gene modules that are common to various conditions and those that are specific to the combinations.
In summary, it is likely to find co-expressed genes in the different abiotic stress conditions, which would suggest similarities in the gene response to these conditions. However, accurate interpretation of the results will require further analysis and consideration of the underlying biology of the experimental conditions being studied.

 ## About the data 
The count data and metadata were taken from: Damián Balfagón, Zandalinas, S. I., dos, T., Claudete Santa‐Catarina, & Gómez‐Cadenas, A. (2022). Reduction of heat stress pressure and activation of photosystem II repairing system are crucial for citrus tolerance to multiple abiotic stress combinations. Physiologia Plantarum, 174(6). https://doi.org/10.1111/ppl.13809 and from the Gene Expression Omnibus Data (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203331) 
The data contains 48 samples with 6 biological replicates for each stress condition. You can find the GDC data here (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203331)

## Methods
The WGCA package has three different ways to perform network construction and module discovery. 
1) An automatic detection of the construction of the network and its modules. This way is convenient for those people who seek to reach a result with minimal effort or with limited computational tools. 
2) A step-by-step construction of the network and module detection. Using this option it is possible to manage and choose our own parameters for the generation of the network and modules. It is also characteristic of not-so-large volumes of expression data. 
3) An automatic way to build block networks and module detection. Special for extremely large gene expression data. It is important to emphasize that this step requires the user to have available multiple cores and a graphical interface with R.  
In this tutorial we build on step 2) to build the network.  In this tutorial, I will present the necessary steps before starting with network construction and module discovery. 
To do this, we will need the packages listed below, to also allow the management of certain cores available on our PC. If you are working from your laptop I recommend you first check the configuration of how many cores your processor handles. In general, it is recommended to use 4 cores, however, I will use only 2. Another step before starting is to use the function (stringsAsFactors = FALSE) in R, which is convenient for handling data frames. R, automatically indicates that matrices should be treated as factorial variables, which can cause your data not to load properly when building or importing it. Therefore, if you do not implement this step your columns containing characters (text) become factors. 
After the construction of the network, we will identify highly related modules and assign them to the phenotypic traits of abiotic stress. We will test whether the associations between gene expression and traits are significant and perform a more in-depth analysis with those modules that are statistically significant.  

 ''' R
#Libraries we are going to use
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
library(Pheatmap);

# Allowing only two threads (consult WGCNA packages for more information)
allowWGCNAThreads(2);
ALLOW_WGCNA_THREADS = 2;  
# Initial variables 
Options (stringsAsFactors = FALSE);
enableWGCNAThreads(2)

'''

Subsequently, we will modify our data in such a way that the names of the samples of the counting data are the same as in the metadata. Because our data coming from the GEO Database is not fixed to start the analysis. We will have to modify them and then analyze their gene expression. 

''' R
#############################################################################################################################################################
##### DATA MODIFICATION                   ##################################################################################
#############################################################################################################################################################
#Import count data (you can download from the GEO page)
data <- read.table("counts.txt", header=TRUE, sep="\t")

#Delete columns we will not use
colnames(data)

#We are going to use only the genes protein coding.  
data <- filter( data, Type=="protein_coding")
#12 658 coding genes

#Besides, we need only the samples names and the Genes 
data <- select(data, -Feature_GID, -Feature_TID, -Type, -Gene_Symbol,
               -Gene_Synonym, -Protein_ID, -Product)

Data <- Data[, 1:49]
#Now we have the count data !


##### Obtain meta-data ####################

#Load GEO Data Query metadata 
query <-getGEO("GSE203331", GSEMatrix = TRUE)

#Obtain matrix metadata
metadata <-pData(phenoData(query[[1]]))

#Let's Change samples names
colnames(metadata)

#1) Select samples name, geo_accession and treatment information.
metadata <- select(metadata, c("title", "treatment:CH1"))

#2) Change samples name for data sample names

#assign data samples names to nuevos_nombres

nuevos_nombres <- as.vector(colnames(data)[2:49])


# To the metadata$title names 
metadata$title <- nuevos_nombres

#Apply the change
metadata$title <- colnames(data)[2:49]



#Test the position of the correct name 
print(nuevos_nombres)
names <- as.data.frame(nuevos_nombres)
#Complement with the different treatment data
phenodata <- data.frame(
  treatment = c("Control", "Control", "Control", "Control", "Control", 
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
                "Water stress", "Water stress", "Water stress", "Water stress", "Water stress", 
                "Water stress and high light combination", 
                "Water stress, high light and heat stress combination", 
                "Water stress, high light and heat stress combination", 
                "Water stress, high light and heat stress combination", 
                "Water stress, high light and heat stress combination", 
                "Water stress, high light and heat stress combination", 
                "Water stress, high light and heat stress combination")
  
  
)

# Join and apply the changes
metadata <- cbind(names, metadata, phenodata)
colnames(metadata)

#Now we have the same samples names
new_data <- select(metadata, "nuevos_nombres", "treatment")
'''

Now that we have modified the count and metadata data to work with WGCA, it is also necessary to modify the metadata and convert it to a binary array with the metadata.  This matrix will allow WGCNA to identify which metadata comes from which sample. For example, not all samples are heat stress, with the matrix with binary data it will be marked with 1 if that sample is heat stress and 0, otherwise. Therefore, it is expected to have a column with the samples and other columns with each of the stress conditions in the plant. 

''' R
#Now let's prepare the binary metadata, first we are going to binarize control 
traits <- new_data%>% 
  mutate(treatment = ifelse(grepl('Control', treatment), 1, 0))  
  
#And the rest of the metadata 
new_data$treatment <- factor(new_data$treatment , 
                                              levels = c("Control",
                                                        "High light and heat stress combination",
                                                        "High light",
                                                        "Heat stress",
                                                        "Water stress and high light combination",
                                                        "Water stress and heat stress combination",
                                                        "Water stress",
                                                        "Water stress, high light and heat stress combination"))
#We use binarizecategoricalColumns from WGCNA 
severity.out_1 <- binarizeCategoricalColumns(new_data$treatment,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)

#Join both data frames
PhenoData <- cbind(traits, severity.out_1)

colnames(PhenoData)

#Verify the number of columns 
Pheno <- rownames_to_column(PhenoData)
PhenoData <- Pheno[,2:9]

#Save R session
save(data, PhenoData,file="matrix_data. RData")    
load("matrix_data. RData")

 '''
Now we have the metadata modified for WGCNA, in the same way we have to make the recommendations of WGCNA with the counting data. According to the library, it is necessary to do a variance stabilization with the Deseq2 package before starting with the creation of the network and identification of modules.  This step is known as gene expression analysis and will serve to verify that the samples are not too dispersed and can harm the analysis, also to eliminate genes with very low counts. Specifically, WGCNA suggests a cleansing of the count and sample data prior to variance stabilization. 

''' R
##############################################################
####### Gene Expression Analysis ############

## WGCNA requires the application of variation from Deseq2 library into the data counts
    
#Delete Na from count matriz
data <- na.omit(data)

data <- data %>%
    column_to_rownames(var='Entrez_Gene_ID')

# Create dds (see Deseq2 tutorial for more details)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = PhenoData,
                              design = ~1) # not spcifying model

## remove all genes with counts < 10 in more than 90% of samples (48*0.90=43.2)
## suggested by WGCNA on RNAseq FAQ
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html

dds90 <- dds[rowSums(counts(dds) >= 10) > 43 ]
nrow(dds90) # from 12, 657 now we have -> 8,316 coding genes.

# Perform variance stabilization (You can find more info on WGCNA FQ)
dds_norm <- vst(dds90)

# Get normalized counts (transform norm data to WGCNA entry)
norm.counts <- assay(dds_norm) %>% 
  t()


save (dds_norm, norm.counts, file="expression_analysis. RData")
load("expression_analysis. RData")
##############################################################
###### Graphic visualization ############

#Let's visualize the dds_norm and the sample distribution
#Data frame transformation
rld_mat_wt <- Assay(dds_norm)

#Correlation values
vsd_cor_wt <- cor(rld_mat_wt)

###### Heatmap ####################

library(pheatmap)
# Prepare the data
data_for_heatmap <- as.matrix(vsd_cor_wt)

# Convert to character data
annotation_row <- as.character(new_data$treatment)

# Add word spaces
annotation_row_with_spaces <- paste(" ", annotation_row, " ")

# Graphic heatmap 
pheatmap(data_for_heatmap,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main="Treatment sample distribution",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         show_rownames = FALSE,
         show_colnames = TRUE,
         row_names_side = "left",
         annotation_colors = "black",
         annotation_names_row = FALSE,
         labels_row = annotation_row_with_spaces,
         fontsize_row = 8, # Adjust the font size of row labels
         fontsize_col = 12, # Adjust the font size of column labels
         angle_col = "45") # Set the angle of column labels to 45 degree

#########################################################################3
## Hclust graphic is another way to verify the sample distribution 
htree <- hclust(dist(t(rld_mat_wt)), method = "average")
plot(htree) 

#We can see the application of the variance over the samples with a boxplot 
boxplot(rld_mat_wt,outline=FALSE, main="After Stabilization", show_colnames = TRUE,
        xaxt="n")
'''

With the heatmap and the hclust we can visualize that the biological samples coincide and are not so far from each other. Therefore, we can continue with the analysis without removing any of the samples. 

Next, we will perform the first WGCA step for the construction of the network. In this step we will build a weighted gene network, for this we will have to choose a threshold potential that corresponds to the similarity of co-expression. The pickSoftThreshold function helps to perform network topology analysis and choose the appropriate similarity.  This first step is vital and we must rely on the best statistical interpretation obtained. 
In total we will have the following indices: 

$fitIndices
   Power SFT.R.sq slope truncated. R.sq mean.k.  median.k.max.k.

Where Powers is the choice of how to visualize the data, in this case we select that we want the results in 30 clusters. So it will show us from 1 to 30 the statistics. It is also the potency that is being evaluated. This is used to transform the gene expression matrix into a similarity matrix that will be used to build the co-expression network.

Sft.R.sq measures variation in the similarity matrix, so it can explain the relationships between genes. A high value indicates a better power adjustment and therefore a more precise network.

Slope: Presents the slope of the curve. It is used to see how Sft.R.sq increases as the power changes. A high value suggests that Sft.R.sq increases rapidly, i.e. it has an effective potency. It is also important to determine when our network stabilizes. 

Truncated.R.sa: Determines whether the power is high enough to capture relevant co-expression relationships

Mean.k (k mean): Average degree of nodes in the coexpression network built with the specified power. Therefore, it is the connection number that a gene has in the network. The higher, the denser a network would be.

Median.k: Median degree of nodes in the coexpression network constructed with the specified power. Like the middle k, the median k is the number of connections of a gene in the network (number of genes at that node)

Mx.k: maximum degree of network nodes. How connected is the gene to the network? 






‌
