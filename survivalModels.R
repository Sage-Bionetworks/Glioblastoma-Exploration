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

survivalModels <- function(){

print("Loading required packages")
require(synapseClient)
require(ggplot2)
require(survival)

## PULL IN GBM DATA BY SOURCING "intermediate data a" ENTITY FROM SYNAPSE
print("Loading data from Synapse")
dataReturn <- loadEntity(299091)
gbmPat <- dataReturn$objects$gbmPat
gbmClin <- dataReturn$objects$gbmClin
gbmMat <- dataReturn$objects$gbmMat
tmpSurv <- dataReturn$objects$tmpSurv
 
#####
## ASSOCIATION OF EXPRESSION WITH SURVIVAL
#####
print("Building Cox models and assessing transcript significance")
plot(survfit(tmpSurv ~ 1))

coxRes <- apply(gbmMat, 1, function(x){
  summary(coxph(tmpSurv ~ x))$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
})

## P-VALUE HISTOGRAMS
print("Generating p-value histograms")
pvalHist <- qplot(coxRes[2, ], geom = "histogram") + 
  opts(title = "GBM Transcripts and Survival Unadjusted p-values")
adjPvals <- p.adjust(coxRes[2, ], method = "BH")
adjPvalHist <- qplot(adjPvals, geom = "histogram") + 
  opts(title = "GBM Transcripts and Survival B-H Adjusted p-values")

## VOLCANO PLOT
print("Generating volcano plot")
vplotDF <- as.data.frame(t(rbind(log2(coxRes[1, ]), -1*log10(coxRes[2, ]))))
colnames(vplotDF) <- (c("Column1", "Column2"))

volcanoPlot <- ggplot(vplotDF, aes(Column1, Column2)) + geom_point() +
  opts(title = "Volcano Plot GBM Transcripts") +
  ylab("- log10 P-values") +
  xlab("Transcripts") +
  opts(plot.title = theme_text(size = 14))

print("Generating output list")
print("To inspect returned objects, look at 'yourReturn$objects'")
return(list("gbmMat" = gbmMat,
            "tmpSurv" = tmpSurv,
            "coxRes" = coxRes,
            "pvalHist" = pvalHist,
            "adjPvalHist" = adjPvalHist,
            "volcanoPlot" = volcanoPlot))
## These objects are saved to Synapse Entity 299114
}

