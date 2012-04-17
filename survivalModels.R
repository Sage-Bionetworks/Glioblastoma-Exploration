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

## PULL IN GBM DATA
## THIS WILL BE A loadEntity(12345) OF CODE ENTITY STORING getGBMdata.R
source("./getGBMdata.R")


#####
## ASSOCIATION OF EXPRESSION WITH SURVIVAL
#####
plot(survfit(tmpSurv ~ 1))

coxRes <- apply(gbmMat, 1, function(x){
  summary(coxph(tmpSurv ~ x))$coefficients[1, c("exp(coef)", "Pr(>|z|)")]
})

## P-VALUE HISTOGRAMS
hist(coxRes[2, ], main = "p-value histogram",
                  xlab = "p-value of gene expression on overall survival")
hist(p.adjust(coxRes[2, ], method="BH"), main = "adjusted p-value histogram",
                                         xlab = "adjusted p-value of gene expression on overall survival")



