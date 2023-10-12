
##########################################################################################################
#                                    LIBRARIES
##########################################################################################################
library(WGCNA);              
library(tidyverse);         
library(dynamicTreeCut);    
library(flashClust);        
library(lattice);           
library(survival);         
library(Hmisc);             
library(here);              
library(GEOquery);
library(DESeq2);
library(pheatmap);

#Work splace        
setwd("C:/Users/Teresita Parra/Pictures/Abiotic_Stress")    


allowWGCNAThreads(2);
ALLOW_WGCNA_THREADS = 2;                
# Initial variables 
options(stringsAsFactors = FALSE);    
enableWGCNAThreads(2);

##########################################################################################################
#                   LOAD DATA AND META-DATA; PREPROCESS AND FORMAT DATA
##########################################################################################################

#Count Data
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


print(metadata$`treatment:ch1`)

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
PhenoData <- Pheno[,2:10]


#PhenoData <- PhenoData %>%
           #column_to_rownames(var='nuevos_nombres')

save(data, PhenoData,file="matrix_data.RData")    
load("matrix_data.RData")


##############################################################
############### Gene Expression Anaylsis ############
###############################################################

## WGCNA requires the aplication of variantion from Deseq2 library into 
#the data counts
colnames(data)

         
#delete Na from count matriz

data <- na.omit(data)

data <- data %>%
    column_to_rownames(var='Entrez_Gene_ID')

# create dds
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = PhenoData,
                              design = ~1) # not spcifying model


## remove all genes with counts < 10 in more than 90% of samples (48*0.90=43.2)
## suggested by WGCNA on RNAseq FAQ
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html


dds90 <- dds[rowSums(counts(dds) >= 10) > 43 ]
nrow(dds90) # from 12, 657 now we have -> 8,316 coding genes


# perform variance stabilization (You can find more info en WGCNA FQ)
dds_norm <- vst(dds90)


# get normalized counts (transform norm data to WGCNA entry)
norm.counts <- assay(dds_norm) %>% 
  t()


save(dds_norm, norm.counts, file="expression_analysis.RData")


load("expression_analysis.RData")
#Let's visualize the dds_norm an the sample distribution

#Data frame transformation
rld_mat_wt <- assay(dds_norm)

#Correlation values
vsd_cor_wt <- cor(rld_mat_wt)
############### Heatmap ####################

library(pheatmap)
# Preparar los datos
data_for_heatmap <- as.matrix(vsd_cor_wt)

# Convert tissueType a un vector de caracter
annotation_row <- as.character(new_data$treatment)

# Añadir espacios entre palabras 
annotation_row_with_spaces <- paste(" ", annotation_row, " ")

# Graficar el heatmap usando la librearia pheatmap 
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
         fontsize_row = 8,     # Adjust the font size of row labels
         fontsize_col = 12,    # Adjust the font size of column labels
         angle_col = "45")       # Set the angle of column labels to 45 degree

#########################################################################3
## We can appreciate biological replicates from CT together and CAR together

htree <- hclust(dist(t(rld_mat_wt)), method = "average")
plot(htree) 

#We can see the application of the variance over the samples with a boxplot too

boxplot(rld_mat_wt,outline=FALSE, main="After Stabilization", show_colnames = TRUE,
        xaxt="n")


#####################################################################
#             TOPOLOGY FOR SOFT-THRESHOLDING
######################################################################


 
powers = c(c(1:10), seq(from =10, to=30, by=1)); 

sft = pickSoftThreshold(norm.counts,
                        dataIsExpr = TRUE,
                        powerVector = powers,
                        corFnc = cor,                             
                        verbose = 5,                              
                        networkType = "signed");                  


write.table(sft, "results_SFT_corr.txt")

# Plot the results
sizeGrWindow(5, 5)
par(mfrow = c(1,2));      
cex1 = 1 

# Scale-free topology fit index as a function of the soft-threshold power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",
     ylab="Scale Free Topology",
     type="n", main = paste("Scale independence for \n Cleopatra"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
 abline(h=0.70,col="red");  
 

# Mean connectivity as a function of the soft-threshold power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=160,col="red");


softPower = 28; 

##########################################################################
#   Generating adjacency and TOM 
###################################################

#calculate the adjacency matrix
adjacency = adjacency(norm.counts,type = "signed",   
                     power = softPower);
dim(adjacency)
head(adjacency)


TOM = TOMsimilarityFromExpr(adjacency,                         
                          TOMType = "signed", 
                          power = softPower);
TOM[1:5,1:5]
dim(TOM)

dissTOM = 1-TOM

save(TOM,dissTOM, softPower,adjacency, file="DissTOM.RData")

load("DissTOM.RData")

##########################################
##    Generate Modules 
##########################################


geneTree = hclust(as.dist(dissTOM), method="average")

## Plot the results
sizeGrWindow(9, 12)
#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="Gene clusters", 
     main="Gene clustering on TOM−based dissimilarity", 
     cex.main = 1,
     cex.axis = 1,
     cex=0.3)   

########################################################################################
##    Module identification using dynamic tree cut 
########################################################################################

## Set the minimum module size
minModuleSize = 20;
diff(geneTree$height)


dynamicMods = cutreeDynamic(dendro= geneTree, 
                            distM = dissTOM,   
                            deepSplit=2,   
                            pamRespectsDendro= FALSE,     
                            minClusterSize = minModuleSize)


## when cutHeight not given, for method=="tree" it defaults to 0.99, for method=="hybrid" it defaults to 99% of the range between the 5th percentile and the maximum of the joining heights on the dendrogram

table(dynamicMods)

write.table(table(dynamicMods), file = "results_dynamicMods.txt", sep = ",", quote = FALSE, row.names = F)

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)      
sort(table(dynamicColors), decreasing = TRUE)   
#24 diff colors 

write.table(sort(table(dynamicColors), decreasing = TRUE), file = "results_dynamicColors.txt", sep = ",", quote = FALSE, row.names = F)

# Plot the dendrogram and colors 
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.05,
                    addGuide = TRUE, guideHang = 0.1,
                    cex.main = 1,
                    main = "Gene dendrogram and Module Colors")

save(dynamicMods, dynamicColors, geneTree, file="DynamicMods.RData")

load("DynamicMods.RData")


###########################################################################
# Merging of modules whose expression profiles are very similar
###############################################################

# Calculate eigengenes 
MEList= moduleEigengenes(norm.counts, colors= dynamicColors,
                         excludeGrey = TRUE,
                         softPower = softPower)

MEs = MEList$eigengenes
MEs[1:5,]


MEDiss = 1-cor(MEs)          

METree = hclust(as.dist(MEDiss),method= "average")

save(MEs, MEDiss, METree, file= "Module_Identification.RData")

load("Module_Identification.RData")

# Clustering of module eigengenes
sizeGrWindow(6,14)

## plots tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes", 
     cex.main = 1.2,
     xlab= "",      
     ylab = "",
     cex = 0.5,    
     sub= "");     

abline(h=0.50,col="red"); 

#################################################################
##   merging modules. 
###################################################################


MEDissThres = 0.00

abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(norm.counts, dynamicColors, 
                          cutHeight= MEDissThres, 
                          verbose =3)
merge$dendro
merge$oldDendro
merge$newMEs
dim(merge$newMEs)
merge$cutHeight

## The merged module colors
mergedColors = merge$colors

## Eigengenes of the new merged modules
mergedMEs = merge$newMEs
mergedMEs$MEblack[1:5]

length(mergedMEs)   
sort(table(mergedColors), decreasing = TRUE)
write.table(sort(table(mergedColors), decreasing = TRUE), file = "mergedColors.txt", sep = ",", quote = FALSE, row.names = F)



## plot dendrogram with module colors below it
plotDendroAndColors(geneTree, main= "Clustering of MEs ", 
                    cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, 
                    cex.main = 1,
                    hang=0.03, addGuide= TRUE, 
                    guideHang=0.05)


colorOrder = c("grey", standardColors(50));   # standardColors(50)
moduleLabels = match(mergedColors, colorOrder)-1     # moduleColors
MEs = mergedMEs


save(MEs, moduleLabels, moduleColors, file= "MergedMods.RData")

load("MergedMods.RData.RData")


##########################################################################
##Calculation TOM plot 
########################################################################

sort(table(mergedColors), decreasing = TRUE)


diag(dissTOM) = NA;

TOMplot(dissTOM,
        geneTree,
        as.character(mergedColors[mergedColors]),
        ColorsLeft = mergedColors,
        terrainColors = TRUE,
        main = "TOM graphic");

#########################################################################
##       Correlate traits
################################################################################

# Define the number of genes and samples
nGenes = ncol(norm.counts)
nSamples = nrow(norm.counts)

## Recalculate MEs with color labels.   
MEs0 = moduleEigengenes(norm.counts, moduleColors)$eigengenes  
MEs = orderMEs(MEs0)

# Link the external values to the correlation matrix

PhenoData <- PhenoData %>%
             column_to_rownames(var="nuevos_nombres")
moduleTraitCor2 = cor(MEs, PhenoData, use= "p")

# Calculates Student asymptotic p-value for given correlations.
moduleTraitPvalue2 = corPvalueStudent(moduleTraitCor2, nSamples)   

textMatrix2= paste(signif(moduleTraitCor2, 2))

dim(textMatrix2) = dim(moduleTraitCor2)


sizeGrWindow(9, 12)
par(mar= c(3.5, 10, 2, 1))   


labeledHeatmap(Matrix= moduleTraitCor2,
               xLabels= names(PhenoData),
               yLabels = names(MEs),
               yLabelsPosition = "left",
               yColorWidth = 0.1,
               cex.lab.y = 0.8,
               cex.lab.x = 0.7,
               ySymbols= names(MEs),
               colors= blueWhiteRed(50),
               textMatrix= (textMatrix2),
               setStdMargins= FALSE,
               cex.text= 0.6,
               zlim= c(-1,1),
               main= paste("Module-trait for Plant Abiotic Stress"))


 

