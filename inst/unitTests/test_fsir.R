## Testing some stuff related to the familial standardized incidence ratio
## (FSIR) from Kerber (1995).
data(minnbreast)
mbsub <- minnbreast
##mbsub <- minnbreast[minnbreast$famid %in% 4:19, ]
mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(mbped) <- FamAgg:::.PEDCN
fad <- FAData(pedigree=mbped)
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id
trait(fad) <- tcancer


notrun_test_fsir <- function(){
    do.plot <- FALSE
    affected <- fad$affected
    ped <- pedigree(fad)
    kin <- kinship(fad)
    nas <- is.na(affected)
    affected <- affected[!nas]
    ped <- ped[!nas,]
    kin <- kin[!nas, !nas]
    nas <- is.na(ped$sex)
    affected <- affected[!nas]
    ped <- ped[!nas,]
    kin <- kin[!nas, !nas]
    ## get the endage for the individuals:
    endages <- minnbreast$endage
    names(endages) <- as.character(minnbreast$id)
    endages <- endages[as.character(ped$id)]
    ## Remove all those that have an NA endage
    nas <- is.na(endages)
    endages <- endages[!nas]
    kin <- kin[!nas, !nas]
    ped <- ped[!nas, ]
    affected <- affected[!nas]
    stratMat <- FamAgg:::factor2matrix(ped$sex)
    stratMatAge <- stratMat * endages
    ## lambda would be the incidence ratio in the population.
    ## using incidence ratio from http://www.cancerresearchuk.org/health-professional/cancer-statistics/statistics-by-cancer-type/breast-cancer/incidence-invasive#heading-Zero
    ## these are the incidence rates in UK!
    FemaleBreastC <- 155.3/100000
    MaleBreastC <- 1.1/100000
    ## and for prostate cancer in UK
    ## http://www.cancerresearchuk.org/health-professional/cancer-statistics/statistics-by-cancer-type/prostate-cancer/incidence#heading-Zero
    FemaleProstateC <- 0
    MaleProstateC <- 134.3/100000
    ## the same number for the US:
    ## http://seer.cancer.gov/statfacts/html/prost.html
    ## these are new cases per year and 100000 males.
    ## FemaleProstateCUS <- 0
    ## MaleProstateCUS <- 137.9/100000

    ## Consider:
    ## shouldn't be timeInStrata the actual time at risk for each individual? This
    ## might then be interpreted for the actual example that the strata is not just
    ## 1 and 0 for male/female, but these numbers multiplied with the endage for each
    ## participant. That way, e.g. male 1 is endage years at risk.

    lambda <- c(M=(1.1+134.4)/100000, F=155.3/100000)
    ## what would be the rates estimated from the pedigree data set?
    Counts <- table(ped$sex[ped$affected > 0])
    ## relate that to the number of M/F.
    lambdaInternal <- Counts / as.numeric(table(ped$sex))

    ## That would be the FSIR with strata being male/female
    FSIRs <- FamAgg:::fsir(affected=affected, kin=kin,
                           lambda=lambda, timeInStrata=stratMat)
    ## the same with strata being the actual time at risk in the sex strata:
    FSIRsAge <- FamAgg:::fsir(affected=affected, kin=kin, lambda=lambda,
                              timeInStrata=stratMatAge)
    if(do.plot){
        par(mfrow=c(1, 2))
        plot(density(FSIRs, na.rm=TRUE))
        abline(v=1, col="grey")
        plot(density(FSIRsAge, na.rm=TRUE))
        abline(v=1, col="grey")
    }
    ## we get unexpectedly large FSIRs if we don't consider the age/time at risk.
    sum(FSIRs < 20 & FSIRs > 0, na.rm=TRUE)
    ## looks better for those that are
    sum(FSIRsAge > 1, na.rm=TRUE)
    sum(FSIRsAge > 2, na.rm=TRUE)

    if(do.plot){
        ## check: boxplot for FSIRs of affected and of unaffected.
        boxplot(split(FSIRsAge, f=affected), ylab="FSIR", xlab="affected")
    }

    ## Now doing it with the method.
    ## First constract time at risk strata matrix
    stratMat <- FamAgg:::factor2matrix(fad$sex)
    stratMat <- stratMat * mbsub$endage
    allFsirs <- FamAgg:::FSIR(fad, trait=trait(fad), lambda=lambda,
                              timeInStrata=stratMat)
    intAllFsirs <- FamAgg:::FSIR(fad, lambda=lambdaInternal,
                                 timeInStrata=stratMat)
    ## compare the calculated FSIRs.
    if(do.plot){
        par(mfrow=c(1, 2))
        plot(density(allFsirs, na.rm=TRUE), main="External lambda")
        abline(v=1, col="grey")
        plot(density(intAllFsirs, na.rm=TRUE), main="Internal lambda")
        abline(v=1, col="grey")
        ## Something is wrong with the internal lambda...
    }

    ## Check the family for some of the guys with the highest risk...
    id <- names(sort(allFsirs, decreasing=TRUE))[1]
    fam <- family(fad, id=id)
    if(do.plot){
        plotPed(fad, family=fam[1, "family"], label1=allFsirs[as.character(fam$id)])
    }
    ## OK, so male individuals have a higher FSIR, and parents of a single affected
    ## child that are not in kinship with other individuals in the family.

    id <- names(sort(allFsirs, decreasing=TRUE))[3]
    fam <- family(fad, id=id)
    if(do.plot){
        plotPed(fad, family=fam[1, "family"], label1=allFsirs[as.character(fam$id)], cex=0.5)
    }

    ## Check family 432
    fam <- family(fad, family="432")
    if(do.plot){
        plotPed(fad, family=fam[1, "family"], label1=allFsirs[as.character(fam$id)], cex=0.5,
                only.phenotyped=TRUE)
    }

    ## Aggregate the FSIR per family.
    fsirPerFamMat <- aggregate(allFsirs, by=list(fad$family), mean, na.rm=TRUE)
    fsirPerFam <- fsirPerFamMat[, "x"]
    names(fsirPerFam) <- fsirPerFamMat[, 1]

    head(sort(fsirPerFam, decreasing=TRUE))
    if(do.plot){
        plot(density(fsirPerFam), main="mean FSIR per family", xlab="mean FSIR")
        abline(v=1, col="grey")
    }

    ## plot for some of the top families.

    ## Compare FSIR to FR


    ## binning by time interval.
    endages <- mbsub$endage
    ## stratify in < 40 > 40
}


## calculate person-time at risk:
test_slice_age <- function(){
    ## given: ages, slice that into time span.
    Ages <- c(13, 65, 45, 35, 19, 34, 49, 45, 40, 39, 17)
    AgeTable <- sliceAge(Ages)
    checkEquals(max(AgeTable[, 1]), 40)
}



