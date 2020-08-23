##############################################################################################################
### LOAD PACKAGES AND UPLOAD DATA
##############################################################################################################

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 

install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") ) 
source("http://bioconductor.org/biocLite.R") 
biocLite("impute")
install.packages("WGCNA")

# Load the package
library(dynamicTreeCut)
library(cluster)
library(flashClust)
library(Hmisc)
library(reshape)
library(foreach)
library(doParallel)
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


#Read in the transcriptomics dataset
transcriptomics_data = read.csv("RNA_COUNTS_LOW_REMOVED_VARIANCE_STABILIZING_TRANSFORMED.csv", head=TRUE, row.names=1);
t_data <- as.data.frame(t(transcriptomics_data), head=TRUE)

#pick a threshold
powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(t_data,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


softPower = 8;

SubGeneNames=rownames(t_data)
#calclute the adjacency matrix
adj= adjacency(t_data,type = "unsigned", power = softPower);

#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(t_data,networkType = "unsigned", TOMType = "unsigned", power = softPower);

colnames(TOM) =rownames(TOM)
dissTOM=1-TOM

### MODULE DETECTION
#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average");

#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="",cex=0.3);

# Set the minimum module size
minModuleSize = 20;

# Module identification using dynamic tree cut

dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

#Plot the module assignment under the dendrogram; note: The grey color is reserved for unassigned genes
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#discard the unassigned genes, and focus on the rest
restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(t_data[,restGenes], power = softPower)

colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

#set the diagonal of the dissimilarity to NA 
diag(diss1) = NA;

#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
sizeGrWindow(7,7)
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))


#EXTRACT MODULES OF GENES IN TEXT FILES
module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}


#LOOK AT EXPRESSION PATTERNS OF GENES AS THEY ARE CLUSTERED
module.order <- unlist(tapply(1:ncol(t_data),as.factor(dynamicColors),I))
m<-t(t(t_data[,module.order])/apply(t_data[,module.order],2,max))
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])

#QUANTIFY MODULE SIMILARITY BY EIGENGENE CORRELATION (eigengenes=module representatives)
MEList = moduleEigengenes(t_data, colors = dynamicColors)
MEs = MEList$eigengenes
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))