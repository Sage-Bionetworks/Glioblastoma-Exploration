# hdClassDiscovery.R

# Erich S. Huang
# Sage Bionetworks
# Seattle, Washington
# erich.huang@sagebase.org

# NOTES
# Here we take the data subsetted to transcripts significantly associated with 
# survival, perform p-value adjustment on these, subset the transcripts

### LOAD IN OTHER REQUIRED LIBRARIES
print("Loading required packages")
require(synapseClient)
require(HDclassif)
require(mclust)
require(Biobase)
require(ggplot2)
require(corpcor)
require(survival)

## LOAD DATA FROM PREVIOUS STEPS
print("Loading intermediate data from Synapse")
theseData <- loadEntity("299114")
coxRes <- theseData$objects$coxRes
gbmMat <- theseData$objects$gbmMat
tmpSurv <- theseData$objects$tmpSurv


## SELECT THE TRANSCRIPTS WITH SIGNIFICANT PVALUES AFTER BH ADJUSTMENT
print("Subsetting the data")
adjPvals <- p.adjust(coxRes[2, ], method = "BH")
topInd <- grep("TRUE", adjPvals <= 0.05)

## SUBSET THE DATA
subsetMat <- gbmMat[topInd, ]

### PERFORM HIGH DIMENSIONAL GAUSSIAN MIXTURE MODELING
print("Running high dimensional Gaussian mixture modeling")
bigClust <- hddc(t(subsetMat),
                 K = 1:10)


return(list("bigClust" = bigClust,
            "subsetMat" = subsetMat))


## STORE OUTPUT FROM HDDC IN SYNAPSE
# myResults <- createEntity(Data(list(name = "results of high dimensional classification algorithm",
#                                  parentId = "275012")))
# myResults <- addObject(myResults, bigClust)
# myResults <- addObjects(myResults, subsetMat)
# myResults <- storeEntity(myResults)
