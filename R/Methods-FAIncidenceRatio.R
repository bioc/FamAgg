## This is the implementation of the FI (familial incidence ratio) and the FSIR
## (familial standardized incidence ratio) from Kerber (1995)

## calculate the time at risk
## startDate: can be either the birth date or the start date for the study.
## endDate: the end date of the study.
## incidenceDate: the date of the incidence.
## deathDate: the date of death.
## the time at risk will be the minimum (earliest) date of (endDate, incidenceDate, deathDate)
## minus the startDate
## need at least endData, in that case endDate is assumed to be "all in one"
## allowNegative: if TRUE return also negative differences. If FALSE, negative values are replaced by 0
## Note: we will subtract half of incidence.subtract for each affected individual.
estimateTimeAtRisk <- function(startDate=NULL, startDateFormat="%Y-%m-%d",
                               endDate=NULL, endDateFormat="%Y-%m-%d",
                               incidenceDate=NULL, incidenceDateFormat="%Y-%m-%d",
                               deathDate=NULL, deathDateFormat="%Y-%m-%d",
                               allowNegative=FALSE, affected=NULL,
                               incidenceSubtract=0.5){
    if(is.null(startDate))
        stop("'startDate' has to be specified (either the start of the study, or the birth date).")
    if(!is(startDate, "Date")){
        message("Transforming startDate into Date format...", appendLF=FALSE)
        startDate <- as.Date(startDate, format=startDateFormat)
        message("OK")
    }
    ## check if we've got at least one of the following:
    if(is.null(endDate)){
        stop("'endDate' has to be specified!")
    }else{
        if(is.null(incidenceDate) & is.null(deathDate))
            message("'endDate' is assumed to specify, for each individual, either the end year of the study, the individual's death date or the individual's incidence date.")
        if(length(endDate) != length(startDate))
            stop("'startDate' and 'endDate' have to have the same length!")
    }
    message("Transforming endDate into Date format...", appendLF=FALSE)
    endDate <- as.Date(endDate, format=endDateFormat)
    message("OK")
    if(!is.null(affected)){
        if(length(affected) != length(startDate)){
            affected <- NULL
            warning("length of affected does not match length of startData; dropping.")
        }else{
            affected <- as.logical(affected)
        }
    }
    ## if we've got an incidenceDate, process that.
    if(!is.null(incidenceDate)){
        message("Transforming incidenceDate into Date format...", appendLF=FALSE)
        incidenceDate <- as.Date(incidenceDate, format=incidenceDateFormat)
        message("OK")
        if(length(incidenceDate) != length(startDate))
            stop("'incidenceDate', 'startDate' and 'endDate' all have to have the same length!")
        ## now get the earlier time point: incidence or end:
        DF <- data.frame(endDate, incidenceDate)
        endDate <- as.Date(apply(DF, MARGIN=1, min, na.rm=TRUE))
        affected <- rep(FALSE, length(incidenceDate))
        affected[!is.na(incidenceDate)] <- TRUE
    }
    if(!is.null(deathDate)){
        message("Transforming deathDate into Date format...", appendLF=FALSE)
        deathDate <- as.Date(deathDate, format=deathDateFormat)
        message("OK")
        if(length(deathDate) != length(startDate))
            stop("'deathDate', 'startDate' and 'endDate' all have to have the same length!")
        ## now get the earlier time point: incidence or end:
        DF <- data.frame(endDate, deathDate)
        endDate <- as.Date(apply(DF, MARGIN=1, min, na.rm=TRUE))
    }
    ## OK, now we're set.
    Diff <- endDate - startDate
    if(!is.null(affected)){
        Diff[affected] <- Diff[affected] - incidenceSubtract
    }
    Diff <- as.numeric(Diff)
    if(!allowNegative)
        Diff[which(Diff < 0)] <- 0
    return(Diff)
}


## takes a numeric vector of ages as input and slices it into
sliceAge <- function(x, slices=c(0, 40, Inf)){
    SliceDiff <- diff(slices)
    Res <- lapply(x, function(z){
        tmp <- z - slices
        tmp <- tmp[1:(length(tmp)-1)]
        bm <- which(tmp > SliceDiff)
        tmp[bm] <- SliceDiff[bm]
        tmp[which(tmp < 0)] <- 0
        return(tmp)
    })
    Res <- do.call(rbind, Res)
    CN <- paste0(slices[1:(length(slices)-1)], ", ", slices[2:length(slices)])
    CN <- paste0("(", CN, "]")
    colnames(Res) <- CN
    return(Res)
}



##****************************************************************************
##
##  familial incidence rate as described in Kerber
##
##
##****************************************************************************
setMethod("familialIncidenceRate", "FAData",
          function(object, trait=NULL, timeAtRisk=NULL, perFamilyTest=FALSE, prune=TRUE){
              return(.FR(ped=pedigree(object), kin=kinship(object), trait=trait,
                         timeAtRisk=timeAtRisk, perFamilyTest=perFamilyTest, prune=prune))
          })


##****************************************************************************
##
##  familial incidence rate along with Monte Carlo simulations to
##  estimate the significance of the ratio.
##  Note: we skipped here the perFamilyTest option as that makes a) no sense,
##  and b) would make the simulations eventually problematic (sampling from
##  small sample sets).
##
##
##****************************************************************************
.FRSimulation <- function(ped=NULL, kin=NULL, trait=NULL, timeAtRisk=NULL,
                          prune=TRUE, nsim=50000, ...){
    ## Input argument checking.
    if(is.null(ped))
        stop("No pedigree submitted!")
    if(is.null(kin))
        stop("No kinship matrix submitted!")
    if(is.null(trait))
        stop("No trait information submitted!")
    if(length(trait) != nrow(ped))
        stop("Argument 'trait' has to have the same length then there are rows (individuals) in the pedigree 'ped'!")
    ped <- cbind(ped, AFF=trait)
    allIds <- as.character(ped$id)
    if(is.null(timeAtRisk)){
        stop("Argument 'timeAtRisk' has to be submitted!")
    }
    if(length(timeAtRisk) != nrow(ped))
        stop("Argument 'timeAtRisk' has to have the same length than there are rows (individuals) in the pedigree 'ped'!")
    timeAtRisk <- as.numeric(timeAtRisk)
    ped <- cbind(ped, TAR=timeAtRisk)
    ## Subset the data set to individuals with valid values for trait, timeAtRisk.
    nas <- is.na(ped$AFF)
    if(any(nas))
        warning(paste0("Removed ", sum(nas), " individuals from the pedigree ",
                       "as they have missing values in the trait."))
    ped <- ped[!nas, , drop=FALSE]
    nas <- is.na(ped$TAR)
    if(any(nas))
        warning(paste0("Removed ", sum(nas), " individuals from the pedigree ",
                       "as they have missing values for the time at risk."))
    ped <- ped[!nas, , drop=FALSE]
    if(prune)
        ped <- doPrunePed(ped)    ## CHECK: Do I really have to do this?
    ## Calculate the kinship for the observed.
    idx <- match(as.character(ped$id), colnames(kin))
    if(any(is.na(idx)))
        stop("Some of the individuals in the pedigree are missing in the kinship matrix!")
    kin <- as.matrix(kin[idx, idx])
    ## Remove self-self kinship and set all affected > 1 to 1.
    diag(kin) <- 0
    Denomi <- ped$TAR %*% kin
    affected <- rep(0, nrow(ped))
    affected[ped$AFF > 0] <- 1
    nAff <- sum(affected)
    sampleFrom <- 1:length(affected)
    obsFr <- as.vector((affected %*% kin) / Denomi)
    names(obsFr) <- ped$id
    ## Simulations: generate random sets of affected individuals, same size then
    ## the observed individuals.
    SimFr <- lapply(1:nsim, function(z){
        simIdx <- sample(sampleFrom, nAff)
        affected[] <- 0
        affected[simIdx] <- 1
        return(as.vector((affected %*% kin) / Denomi))
        ## Would be less memory demanding to just return whether the expected
        ## is larger than the observed
    })
    ## Performing the stuff on the large matrix.
    ## Want to compare the observed FR against the simulated FOR EACH INDIVIDUAL.
    ## Why: that way we keep the kinship the same and just calculate an expexted
    ## FR for random occurence of cases.
    SimFr <- do.call(cbind, SimFr)
    PVals <- rowSums(SimFr >= obsFr)/nsim
    denses <- apply(SimFr, MARGIN=1, density)
    ## Result should match all ids in the (original) pedigree.
    fr <- rep(NA, length(allIds))
    names(fr) <- allIds
    idx <- match(names(obsFr), allIds)
    fr[idx] <- obsFr
    ## p-values
    pv <- rep(NA, length(allIds))
    names(pv) <- allIds
    pv[idx] <- PVals
    ## density
    dens <- vector("list", length(allIds))
    names(dens) <- allIds
    dens[idx] <- denses
    return(list(FR=fr, pvalueFR=pv, expDensity=dens))
}

## prune: remove un-connected individuals from the pedigree AFTER removing individuals without
##        valid value in trait or a missing value in timeAtRisk.
.FR <- function(ped=NULL, kin=NULL, trait=NULL, timeAtRisk=NULL, perFamilyTest=FALSE,
                prune=TRUE, ...){
    ## subset the data to all the individuals for which we do have the required data.
    ## * trait not NA
    ## * timeAtRist not NA
    ## * did not get removed after pruning for un-connected individuals.
    if(is.null(ped))
        stop("No pedigree submitted!")
    if(is.null(kin))
        stop("No kinship matrix submitted!")
    if(is.null(trait))
        stop("No trait information submitted!")
    if(length(trait) != nrow(ped))
        stop("Argument 'trait' has to have the same length then there are rows (individuals) in the pedigree 'ped'!")
    ped <- cbind(ped, AFF=trait)
    if(is.null(timeAtRisk)){
        stop("Argument 'timeAtRisk' has to be submitted!")
    }
    if(length(timeAtRisk) != nrow(ped))
        stop("Argument 'timeAtRisk' has to have the same length than there are rows (individuals) in the pedigree 'ped'!")
    timeAtRisk <- as.numeric(timeAtRisk)
    ped <- cbind(ped, TAR=timeAtRisk)
    if(perFamilyTest){
        entityIs <- "family"
    }else{
        entityIs <- "pedigree"
        ped[, "family"] <- 1
    }
    ## OK, now we can go on: lapply on the splitted ped by family
    pedL <- split(ped, ped$family)
    Res <- lapply(pedL, function(z){
        allIds <- as.character(z$id)
        ## remove NAs:
        ## NA in affected/trait
        nas <- is.na(z$AFF)
        if(any(nas)){
            z <- z[!nas, , drop=FALSE]
            warning(paste0("Removed ", sum(nas), " individuals from ",
                           entityIs, " ", z[1, "family"],
                           " as they have missing values in the trait."))
        }
        ## NA in time at risk:
        nas <- is.na(z$TAR)
        if(any(nas)){
            z <- z[!nas, , drop=FALSE]
            warning(paste0("Removed ", sum(nas), " individuals from ",
                           entityIs, " ", z[1, "family"],
                           " as they have missing values for the time at risk."))
        }
        ## OK, now starting to reduce the data.
        if(prune)
            z <- doPrunePed(z)
        ## Subset the kinship matrix and ensure correct ordering.
        idx <- match(as.character(z$id), colnames(kin))
        if(any(is.na(idx)))
            stop("Some of the individuals in the pedigree are missing in the kinship matrix!")
        kinSub <- as.matrix(kin[idx, idx])
        ## remove self-self kinship.
        diag(kinSub) <- 0
        affected <- rep(0, nrow(z))
        affected[z$AFF > 0] <- 1
        ## The magic formula (3) in Kerber 1995 as an R matrix multiplication.
        FR <- as.vector((affected %*% kinSub)/(z$TAR %*% kinSub))
        names(FR) <- z$id
        retFR <- rep(NA, length(allIds))
        names(retFR) <- allIds
        retFR[names(FR)] <- FR
        return(retFR)
    })
    return(Res)
}


## The tricky thing is labda and timeInStrata. Both have to have the same unit.
setMethod("FSIR", "FAData", function(object, trait=NULL, lambda=NULL, timeInStrata=NULL){
    ## First we need to do a lot of checking and testing.
    if(is.null(trait)){
        ## check internal trait...
        if(length(object@.trait) == 0)
            stop("trait is missing!")
        trait <- trait(object)
    }
    ## Check lambda
    if(is.null(lambda))
        stop("lambda missing!")
    if(is.null(names(lambda)))
        stop("lambda has to be a named vector with the names corresponding to the strata names!")
    ## Check timeInStrata
    if(is.null(timeInStrata))
        stop("timeInStrata missing!")
    if(!is.matrix(timeInStrata))
        stop("timeInStrata has to be a matrix!")
    if(nrow(timeInStrata) != length(object$id))
        stop("timeInStrata has to have the same number of rows as there are individuals in the pedigree!")
    if(length(lambda)!=ncol(timeInStrata))
        stop("length of lambda has to match the number of columns of timeInStrata!")
    if(!all(colnames(timeInStrata) %in% names(lambda)))
        stop("Names of lambda does not match the colnames of timeInStrata!")
    lambda <- lambda[colnames(timeInStrata)]
    ## Done.
    trait(object) <- trait
    trait <- trait(object)   # that way we ensure that we have the same ordering.
    kin <- kinship(object)
    ## Just to be on the save side... ensure that the ordering of id/Trait matches the kin
    kin <- kin[names(trait), names(trait)]
    diag(kin) <- 0
    ## Now start subsetting the data:
    ## * NA in trait
    nas <- is.na(trait)
    if(any(nas)){
        warning(paste0("Excluding ", sum(nas),
                       " individuals because of a missing value in the trait."))
        trait <- trait[!nas]
        kin <- kin[!nas, !nas]
        timeInStrata <- timeInStrata[!nas, , drop=FALSE]
    }
    ## * NA in timeInStrata.
    nas <- apply(timeInStrata, MARGIN=1, function(z){
        any(is.na(z))
    })
    if(any(nas)){
        warning(paste0("Excluding ", sum(nas),
                       " individuals because of a missing value in timeInStrata."))
        trait <- trait[!nas]
        kin <- kin[!nas, !nas]
        timeInStrata <- timeInStrata[!nas, , drop=FALSE]
    }
    ## * Not related, i.e. individuals with a kinship sum of 0
    nas <- colSums(kin) == 0
    if(any(nas)){
        warning(paste0("Excluding ", sum(nas),
                       " individuals because they do not share kinship with any individual in the pedigree."))
        trait <- trait[!nas]
        kin <- kin[!nas, !nas]
        timeInStrata <- timeInStrata[!nas, , drop=FALSE]
    }
    ## Well done. Now let's do the test:
    fsirs <- fsir(affected=trait, kin=kin, lambda=lambda, timeInStrata=timeInStrata)
    ## Prepare the results vector.
    allIds <- object$id
    allFsirs <- rep(NA, length(allIds))
    names(allFsirs) <- allIds
    allFsirs[names(trait)] <- as.numeric(fsirs)
    return(allFsirs)
})


##************************************************
##
## FSIR Kerber 1995
## the familial standardized incidence ratio.
## What do we need to calculate it:
## * affected vector: can get the cases from that.
## * the kinship matrix: setting diag(kin) <- 0, so we exclude self-self kinship.
## * strata: some sort of stratification of cases.
## * time in strata: related to the strata: the time each individual spend in strata.
## * the population incidence rate per stratum: that's really tricky! how do we get the
##   population incidence
##   if we have only the cohort available... sampling? that will not work, if, only with
##   replacement... but still
##
## Side note: what if we use the incidence ratio from the whole pedigree? the FSIR is anyway
##            just calculated based
##            on individuals from the same family (got a kinship value of 0 for other
##            individuals). Such lambda would
##            represent the expected number of affected individuals in the family given the
##            data from the whole pedigree.
##************************************************
## Note: this function assumes that all elements are in the correct order!
## Arguments:
## * affected: 0, 1 coding, length corresponds to the phenotyped; length n.
## * kin: nxn kinship matrix. ordering has to match individuals in affected.
## * lambda: population incidence rate. numeric vector, length has to match the number of strata.
## * timeInStrata: ncol=length(lambda), nrow=n. has to have the same unit than lambda.
fsir <- function(affected, kin, lambda, timeInStrata){
    if(missing(affected))
        stop("affected missing!")
    if(missing(kin))
        stop("kinship matrix kin missing!")
    if(missing(lambda))
        stop("lambda missing!")
    if(missing(timeInStrata))
        stop("timeInStrata missing!")
    if(length(unique(c(length(affected), nrow(kin), ncol(kin), nrow(timeInStrata)))) != 1)
        stop("Length of affected has to match the number of rows and cols of kin and the number of rows of timeInStrata!")

    ## match names of lambda with colnames of timeInStrata:
    if(is.null(names(lambda)))
        stop("lambda has to be a named numeric vector with the names corresponding to the colnames of timeInStrata!")
    if(is.null(colnames(timeInStrata)))
        stop("timeInStrata has to be a matrix with column names corresponding to the names of argument lambda!")
    if(length(lambda) != ncol(timeInStrata))
        stop("length of lambda does not match number of columns of timeInStrata!")
    if(!all(names(lambda) %in% colnames(timeInStrata)))
        stop("names of lambda do not match with colnames of timeInStrata!")
    lambda <- lambda[colnames(timeInStrata)]

    ## setting the diagonal of kin to 0 to avoid self-self kinship
    diag(kin) <- 0
    ## calculate nominator of Kerber 1995 formula (4): the observed cases
    ## this results in a vector length affected, with 0 if the individual is not related to any
    ## affected.
    nomi <- as.numeric(affected %*% kin)
    ## funny thing is the denominator that lists the expected cases.
    ## lambda should be a vector of length k (k=number of strata), timeInStrat a matrix with k columns, n rows
    ## with n being the number of individuals and kin a nxn matrix. lambda and timeInStrat should be in the same unit,
    ## i.e. either proportion, or proportion per person-year.
    denomi <- lambda %*% t(timeInStrata) %*% kin
    if(nrow(denomi) > 1)
        stop("Did not get the expected single-row matrix!")
    denomi <- as.numeric(denomi)
    res <- nomi/denomi
    names(res) <- names(affected)
    return(res)
}

## convert a factor to a matrix, levels
factor2matrix <- function(x){
    if(!is.factor(x))
        x <- factor(x)
    Ls <- lapply(levels(x), function(z){
        return(as.numeric(x == z))
    })
    Mat <- do.call(cbind, Ls)
    colnames(Mat) <- levels(x)
    if(!is.null(names(x)))
        rownames(Mat) <- names(x)
    return(Mat)
}


##
.FSIR <- function(ped=NULL, kin=NULL, trait=NULL, timeAtRisk=NULL, prune=TRUE, strata=NULL, ...){
    ## subset the data to all the individuals for which we do have the required data.
    ## * trait not NA
    ## * timeAtRist not NA
    ## * strata not NA
    ## * did not get removed after pruning for un-connected individuals.
    if(is.null(ped))
        stop("No pedigree submitted!")
    if(is.null(kin))
        stop("No kinship matrix submitted!")
    if(is.null(trait))
        stop("No trait information submitted!")
    if(length(trait) != nrow(ped))
        stop("Argument 'trait' has to have the same length then there are rows (individuals) in the pedigree 'ped'!")
    ped <- cbind(ped, affected=trait)
    ## check strata
    haveStrata <- FALSE
    if(!is.null(strata)){
        if(length(strata)!=nrow(ped))
            stop("Argument 'strata' has to have the same length than there are rows (individuals) in the pedigree 'ped'!")
        ## add strata.
        ped <- cbind(ped, STRATA=strata)
        haveStrata <- TRUE
    }
    if(is.null(timeAtRisk)){
        stop("Argument 'timeAtRisk' has to be submitted!")
    }
    if(length(timeAtRisk) != nrow(ped))
        stop("Argument 'timeAtRisk' has to have the same length than there are rows (individuals) in the pedigree 'ped'!")
    timeAtRisk <- as.numeric(timeAtRisk)
    ped <- cbind(ped, TAR=timeAtRisk)
    ## OK, now starting to reduce the data.
    if(prune)
        ped <- doPrunePed(ped)
    ## NA in affected/trait
    nas <- is.na(ped$affected)
    if(any(nas)){
        ped <- ped[!nas, , drop=FALSE]
        warning(paste0("Removed ", sum(nas), " individuals with missing values in the trait."))
    }
    ## NA in time at risk:
    nas <- is.na(ped$TAR)
    if(any(nas)){
        ped <- ped[!nas, , drop=FALSE]
        warning(paste0("Removed ", sum(nas), " individuals with missing values for the time at risk."))
    }
    ## NA in strata
    if(haveStrata){
        nas <- is.na(ped$STRATA)
        if(any(nas)){
            ped <- ped[!nas, , drop=FALSE]
            warning(paste0("Removed ", sum(nas), " individuals with missing values in strata."))
        }
    }
    affIds <- as.character(ped[which(ped$affected > 0), "id"])
    if(length(affIds) < 2){
        warning("Can not perform the test: less than 2 affected individuals in trait.")
    }
    ## Subset the kinship matrix and ensure correct ordering.
    idx <- match(as.character(ped$id), colnames(kin))
    if(any(is.na(idx)))
        stop("Some of the individuals in the pedigree are missing in the kinship matrix!")
    kin <- as.matrix(kin[idx, idx])
    ## remove self-self kinship.
    diag(kin) <- 0
    ## The magic formula (3) in Kerber 1995:
    affected <- rep(0, nrow(ped))
    affected[ped$affected > 0] <- 1
    FR <- as.vector((affected %*% kin)/(ped$TAR %*% kin))
    names(FR) <- ped$id
    return(FR)
}


