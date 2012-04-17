# visualizeGBMClasses.R

# Erich S. Huang
# Sage Bionetworks
# Seattle, Washington
# erich.huang@sagebase.org

# NOTES
# Code to visualize GBM molecular classes as defined by the class discovery
# object.

### LOAD IN REQUIRED LIBRARIES
require(HDclassif)
require(mclust)
require(ggplot2)
require(corpcor)

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

