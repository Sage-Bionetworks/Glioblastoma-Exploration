# hdClassDiscovery.R

# Erich S. Huang
# Sage Bionetworks
# Seattle, Washington
# erich.huang@sagebase.org

# NOTES
# Here we take the data subsetted to transcripts significantly associated with 
# survival, perform p-value adjustment on these, subset the transcripts

### LOAD IN REQUIRED LIBRARIES
require(HDclassif)
require(mclust)
require(survival)
require(Biobase)
require(ggplot2)
require(corpcor)

### ASSOCIATION OF EXPRESSION WITH SURVIVAL
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

## SELECT THE TRANSCRIPTS WITH SIGNIFICANT PVALUES AFTER BH ADJUSTMENT
topInd <- grep("TRUE", adjPvals <= 0.05)

## SUBSET THE DATA
subsetMat <- gbmMat[topInd, ]

### PERFORM HIGH DIMENSIONAL GAUSSIAN MIXTURE MODELING
bigClust <- hddc(t(subsetMat),
                 K = 1:10)

