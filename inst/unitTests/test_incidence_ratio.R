## incidence ratio test as defined by Kerber 1995
data(minnbreast)
mbsub <- minnbreast[minnbreast$famid %in% 4:19, ]
PedDf <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
colnames(PedDf) <- FamAgg:::.PEDCN
PedDf$sex <- FamAgg:::sanitizeSex(PedDf$sex)
## generate the FAData.
fad <- FAData(pedigree=PedDf)
## specify the trait.
tcancer <- mbsub$cancer
names(tcancer) <- mbsub$id

## the time at risk...
tar <- mbsub$endage


test_incidence_ratio <- function(){

    fr <- FamAgg:::.FR(ped=pedigree(fad), kin=kinship(fad), trait=tcancer,
                       timeAtRisk=tar)
    ## the same but per family
    frFam <- FamAgg:::.FR(ped=pedigree(fad), kin=kinship(fad), trait=tcancer,
                          timeAtRisk=tar, perFamilyTest=TRUE)
    ## So, it's not the same if we run this within family or for the whole pedigree...
    ## frVals <- unlist(fr, use.names=FALSE)
    ## frFamVals <- unlist(frFam, use.names=FALSE)
    ## checkEquals(frVals[!is.na(frVals)], frFamVals[!is.na(frFamVals)])
    ## checkEquals(unlist(fr, use.names=FALSE), unlist(frFam, use.names=FALSE))
    ## fr[[1]][1:10]
    ## frFam[["4"]]
}

test_estimate_time_at_risk <- function(){
    sdates <- c("2012-04-17", "2014-05-29", "1999-12-31", "2002-10-10")
    edates <- c("2015-09-15", "2015-09-15", "2005-09-15", "2015-09-15")
    idates <- c(NA, NA, "2007-07-13", "2013-12-23")
    ddates <- c(NA, NA, NA, "2014-03-14")

    checkException(estimateTimeAtRisk(sdates))
    ests <- estimateTimeAtRisk(startDate=sdates, endDate=edates)
    ests2 <- estimateTimeAtRisk(startDate=sdates, endDate=edates,
                                incidenceDate=idates)
    ests3 <- estimateTimeAtRisk(startDate=sdates, endDate=edates,
                                incidenceDate=idates, deathDate=ddates)
}

notrun_test_fr_simulation <- function(){
    fr <- FamAgg:::.FR(ped=pedigree(fad), kin=kinship(fad), trait=tcancer,
                       timeAtRisk=tar)
    ## the same using simulations.
    system.time(
        frSim <- FamAgg:::.FRSimulation(pedigree(fad), kinship(fad), tcancer,
                                        tar, nsim=1000)
    )
}


