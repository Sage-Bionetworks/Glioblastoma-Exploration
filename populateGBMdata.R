#########################################################################
## PULL IN GBM DATA FROM THE COHERENT DATA LIBRARY (CDL) WITHIN SYNAPSE
#####
##   - AFFYMETRIX U133A EXPRESSION DATA
##   - CLINICAL DATA
#########################################################################

populateGBMdata <- function(){

require(synapseClient)
require(Biobase)
require(survival)


#####
## EXPRESSION DATA
#####
ge <- loadEntity("274865")
gbmExprSet <- ge$objects$coherentEset
gbmMat <- exprs(gbmExprSet)

## GET THE ROW NAMES FROM ANOTHER LAYER
ne <- loadEntity("274544")
tmpNames <- featureNames(ne$objects$eset)

## SUBSET TO THOSE WITHOUT NAS
id <- which(colSums(is.na(gbmMat)) != nrow(gbmMat))
gbmMat <- gbmMat <- gbmMat[, id]
rownames(gbmMat) <- tmpNames

## GET RID OF NON-TUMOR SAMPLES
theseTissues <- sapply(strsplit(colnames(gbmMat), "-", fixed=T), function(x){
  substr(x[4], 1, 2)
})
gbmMat <- gbmMat[ , theseTissues == "01" ]

## LOOK AT PATIENT LEVEL AND REMOVE DUPS
thesePats <- sapply(strsplit(colnames(gbmMat), "-", fixed=T), function(x){
  paste(x[1:3], collapse="-")
})
idd <- duplicated(thesePats)
gbmMat <- gbmMat[ , !idd]
thesePats <- thesePats[!idd]


#####
## CLINICAL DATA
#####
gc <- loadEntity("274426")
gbmClin <- gc$objects$clinAll
gbmPat <- gbmClin$clinical_patient_public_gbm
rownames(gbmPat) <- gbmPat$bcr_patient_barcode
gbmPat <- gbmPat[ thesePats, ]

gbmPat$vital_status[ gbmPat$vital_status == "[Not Available]" ] <- NA
gbmPat$survTime <- as.numeric(gbmPat$days_to_death)
gbmPat$surv <- ifelse( gbmPat$vital_status=="DECEASED", 1, 0)
gbmPat$survTime[ which(gbmPat$vital_status == "LIVING") ] <- as.numeric(gbmPat$days_to_last_followup)[ which(gbmPat$vital_status == "LIVING") ]
tmpSurv <- Surv(gbmPat$survTime, gbmPat$surv)



rm(list=setdiff(ls(), c("gbmPat", "gbmClin", "gbmMat", "tmpSurv")))

return(list("gbmPat" = gbmPat,
            "gbmClin" = gbmClin,
            "gbmMat" = gbmMat,
            "tmpSurv" = tmpSurv))
}


