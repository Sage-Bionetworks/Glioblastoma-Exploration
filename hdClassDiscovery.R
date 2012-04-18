# hdClassDiscovery.R

# Erich S. Huang
# Sage Bionetworks
# Seattle, Washington
# erich.huang@sagebase.org

# NOTES
# Here we take the data subsetted to transcripts significantly associated with 
# survival, perform p-value adjustment on these, subset the transcripts

### LOAD IN OTHER REQUIRED LIBRARIES
require(synapseClient)
require(HDclassif)
require(mclust)
require(Biobase)
require(ggplot2)
require(corpcor)
require(survival)

## LOAD DATA FROM PREVIOUS STEPS
theseData <- loadEntity("275017")
coxRes <- theseData$objects$coxRes
gbmMat <- theseData$objects$gbmMat
tmpSurv <- theseData$objects$tmpSurv


## SELECT THE TRANSCRIPTS WITH SIGNIFICANT PVALUES AFTER BH ADJUSTMENT
adjPvals <- p.adjust(coxRes[2, ], method = "BH")
topInd <- grep("TRUE", adjPvals <= 0.05)

## SUBSET THE DATA
subsetMat <- gbmMat[topInd, ]

### PERFORM HIGH DIMENSIONAL GAUSSIAN MIXTURE MODELING
bigClust <- hddc(t(subsetMat),
                 K = 1:10)



## STORE OUTPUT FROM HDDC IN SYNAPSE
myResults <- createEntity(Data(list(name = "results of high dimensional classification algorithm",
                                 parentId = "275012")))
myResults <- addObject(myResults, bigClust)
myResults <- addObjects(myResults, subsetMat)
myResults <- storeEntity(myResults)
