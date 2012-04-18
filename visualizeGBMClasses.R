# visualizeGBMClasses.R

# Erich S. Huang
# Sage Bionetworks
# Seattle, Washington
# erich.huang@sagebase.org

# NOTES
# Code to visualize GBM molecular classes as defined by the class discovery
# object.


### LOAD IN REQUIRED LIBRARIES
require(synapseClient)
require(HDclassif)
require(mclust)
require(ggplot2)
require(corpcor)
require(survival)

## LOAD DATA FROM PREVIOUS STEPS
theseData <- loadEntity("275017")
coxRes <- theseData$objects$coxRes
gbmMat <- theseData$objects$gbmMat
tmpSurv <- theseData$objects$tmpSurv

theseResults <- loadEntity("275029")
bigClust <- theseResults$objects$bigClust
subsetMat <- theseResults$objects$subsetMat



### TAKING ALL THE DATA AND LOOKING AT THE PC SPACE
fullSVD <- fast.svd(subsetMat)
fullClPlot <- clPairs(fullSVD$v[ , 1:4], bigClust$class)

### A FUNCTION TO QUICKLY PRODUCE VISUALS OF CLASSES
visSubspace <- function(K){
  qMat <- bigClust$Q[[K]]
  kTx <- sort(abs(qMat[ , 1]), decreasing = T)
  kTopTx <- kTx[1:150]
  kSubMat <- subsetMat[names(kTopTx), ]
  kSVD <- fast.svd(kSubMat)
  kSvdDF <- as.data.frame(kSVD$v)
  colnames(kSvdDF) <- paste("PrinComp", 1:4, sep = "")
  pc12Plot <- ggplot(kSvdDF, aes(PrinComp1, PrinComp2)) +
    geom_point(aes(colour = factor(bigClust$class),
                   shape = factor(bigClust$class),
                   size = 20)) +
                     scale_size(guide = "none") +
                     scale_shape(guide = "none")
  clPairPlot <- clPairs(kSVD$v[ , 1:4], bigClust$class)
  return(list("pcDataFrame" = kSvdDF, 
              "ggplotObj" = pc12Plot, 
              "clPairPlot" = clPairPlot))
}

classOneObjects <- visSubspace(1)
classTwoObjects <- visSubspace(2)
classThreeObjects <- visSubspace(3)

### GENERATE KAPLAN MEIER PLOT OF GBM CLASSES
gbmSurvFit <- survfit(tmpSurv ~ bigClust$class)
gbmStrata <- bigClust$class
ggkm(gbmSurvFit, 
     ystratalabs = (c("ClassOne", "ClassTwo", "ClassThree")), 
     timeby = 365,
     main = "GBM K-M Plot By Class")

###
# Alternative visualization (density plots)
fullSvdDF <- as.data.frame(t(rbind(fullSVD$v[ , 1], bigClust$class)))
colnames(fullSvdDF) <- c("PrinComp1", "MolecularClass")
densityPlotONE <- ggplot(fullSvdDF, 
                        aes(x = PrinComp1,
                            fill = factor(MolecularClass))) + 
                              geom_density(alpha = 0.3) +
                              opts(title = "Density Plot GBM Classes on PC1") +
                              ylab("Density") +
                              xlab("Molecular Class") +
                              opts(plot.title = theme_text(size = 14))

fullSvdDF <- as.data.frame(t(rbind(fullSVD$v[ , 2], bigClust$class)))
colnames(fullSvdDF) <- c("PrinComp2", "MolecularClass")
densityPlotTWO <- ggplot(fullSvdDF, 
                      aes(x = PrinComp2,
                          fill = factor(MolecularClass))) + 
                            geom_density(alpha = 0.3) +
                            opts(title = "Density Plot GBM Classes on PC2") +
                            ylab("Density") +
                            xlab("Molecular Class") +
                            opts(plot.title = theme_text(size = 14))

fullSvdDF <- as.data.frame(t(rbind(fullSVD$v[ , 4], bigClust$class)))
colnames(fullSvdDF) <- c("PrinComp4", "MolecularClass")
densityPlotTHREE <- ggplot(fullSvdDF, 
                      aes(x = PrinComp4,
                          fill = factor(MolecularClass))) + 
                            geom_density(alpha = 0.3) +
                            opts(title = "Density Plot GBM Classes on PC4") +
                            ylab("Density") +
                            xlab("Molecular Class") +
                            opts(plot.title = theme_text(size = 14))

