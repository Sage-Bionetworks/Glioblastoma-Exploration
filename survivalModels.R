#########################################################################
## PULL IN GBM DATA FROM THE COHERENT DATA LIBRARY (CDL) WITHIN SYNAPSE
#####
##   - AFFYMETRIX U133A EXPRESSION DATA
##   - CLINICAL DATA
#####
##   - EXPLORE ASSOCIATION OF GENE LEVEL EXPRESSION WITH SURVIVAL
##   - CREATE P-VALUE HISTOGRAMS (RAW AND ADJUSTED)
##   - STORE 'SIGNIFICANT' GENES FOR FURTHER CLASSIFICATION
#########################################################################

require(synapseClient)

## PULL IN GBM DATA BY SOURCING "populate data" CODE ENTITY FROM SYNAPSE
loadGBM <- loadEntity("275016")


#####
## ASSOCIATION OF EXPRESSION WITH SURVIVAL
#####
plot(survfit(tmpSurv ~ 1))

coxRes <- apply(gbmMat, 1, function(x){
  summary(coxph(tmpSurv ~ x))$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
})

## P-VALUE HISTOGRAMS
pvalHist <- qplot(coxRes[2, ], geom = "histogram") + 
  opts(title = "GBM Transcripts and Survival Unadjusted p-values")
adjPvals <- p.adjust(coxRes[2, ], method = "BH")
adjPvalHist <- qplot(adjPvals, geom = "histogram") + 
  opts(title = "GBM Transcripts and Survival B-H Adjusted p-values")

## VOLCANO PLOT
vplotDF <- as.data.frame(t(rbind(log2(coxRes[1, ]), -1*log10(coxRes[2, ]))))
colnames(vplotDF) <- (c("Column1", "Column2"))

volcanoPlot <- ggplot(vplotDF, aes(Column1, Column2)) + geom_point() +
  opts(title = "Volcano Plot GBM Transcripts") +
  ylab("- log10 P-values") +
  xlab("Transcripts") +
  opts(plot.title = theme_text(size = 14))


## STORE DATA IN SYNAPSE
myData <- createEntity(Data(list(name = "results of association of expression with survival",
                                 parentId = "275012")))
myData <- addObject(myData, gbmMat)
myData <- addObject(myData, tmpSurv)
myData <- addObject(myData, coxRes)
myData <- storeEntity(myData)

