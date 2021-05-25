##************************************************
##
##       FAKinSumResults
##       Most of the code here is from Chris X. Weichenberger
##
##************************************************
setMethod("show", "FAKinSumResults", function(object){
    callNextMethod()
    cat(paste0("Result info:\n"))
    cat(paste0(" * Number of rows of the result data.frame: ",
               nrow(result(object)),".\n"))
    cat(paste0(" * Number of simulations: ", object@nsim, ".\n"))
})
## calling the trait replacement method from FAResult and in addition
## reset the simulation result.
setReplaceMethod("trait", "FAKinSumResults", function(object, value){
    object <- callNextMethod()
    ## reset the result
    object@sim <- list()
    object@nsim <- 0
    object@traitname <- character()
    return(object)
})

####
## the analysis method.
setMethod("runSimulation", "FAKinSumResults",
          function(object, nsim=50000, strata=NULL, ...){
              if(length(trait(object)) == 0)
                  stop("No trait information available!")
              ind <- TRUE
              kin <- kinship(object)
              ## check strata
              ped <- pedigree(object)
              affectedIds <- as.character(ped[which(ped$affected > 0), "id"])
              phenotypedIds <-as.character(ped[!is.na(ped$affected), "id"])
              message("Cleaning data set (got in total ", nrow(ped),
                      " individuals):")
              ## Dummy for removing unaffected:
              message(" * not phenotyped individuals...")
              if(length(phenotypedIds) != nrow(ped)){
                  message(" ", nrow(ped)-length(phenotypedIds), " removed")
              }else{
                  message(" none present.")
              }
              if(!is.null(strata)){
                  if(length(strata)!=nrow(ped))
                      stop("Argument 'strata' has to have the same length than",
                           " there are individuals in the pedigree!")
                  ## add strata.
                  ped <- cbind(ped, STRATA=strata)
                  ## Subset the affectedIds and phenotypedIds to those for
                  ## which strata is not NA.
                  ped <- ped[!is.na(ped$STRATA) & !is.na(ped$affected), ,
                             drop=FALSE]
                  strata <- as.character(ped$STRATA)
                  names(strata) <- as.character(ped$id)
                  newAffIds <- as.character(ped[which(ped$affected > 0), "id"])
                  newPhenoIds <- as.character(ped[, "id"])
                  message(paste0(" * unaffected individuals without valid ",
                                 "strata values..."), appendLF=FALSE)
                  if(length(newPhenoIds)!=length(phenotypedIds)){
                      message(" ", length(phenotypedIds)-length(newPhenoIds),
                              " removed.")
                      phenotypedIds <- newPhenoIds
                  }else{
                      message(" none present.")
                  }
                  message(paste0(" * affected individuals without valid strata",
                                 " values..."), appendLF=FALSE)
                  if(length(newAffIds)!=length(affectedIds)){
                      message(" ", length(affectedIds)-length(newAffIds),
                              " removed.")
                      affectedIds <- newAffIds
                  }else{
                      message(" none present.")
                  }
              }
              message("Done")
              ## now forcing the kinship matrix to match ordering of
              ## phenotypedIds and eventually strata.
              kin <- kin[phenotypedIds, phenotypedIds]
              ## removing self-self kinship:
              diag(kin) <- NA
              object@nsim <- nsim
              ## break...
              if(length(affectedIds) <= 1){
                  warning(paste0("Trait ",
                                 object@traitname,
                                 ": only one or none affected individual!"))
                  return(object)
              }
              if(ind){
                  res <- significant.individuals.althyp(kin, affectedIds,
                                                        phenotypedIds,
                                                        nr.sim=nsim,
                                                        strata=strata,
                                                        ...)
              }else{
                  ## res <- significan.population(kin, affectedIds, phenotypedIds,
                  ##                              nr.sim=nsim)
              }
              object@sim <- res
              return(object)
          })

## results; get the results table.
setMethod("result", "FAKinSumResults",
          function(object, method="BH", cutoff=0.05, rmKinship=0){
    method <- match.arg(method, p.adjust.methods)
    TraitName <- object@traitname
    if(length(TraitName)==0)
        TraitName <- NA
    if(length(object@sim) == 0){
        ## generate a dummy matrix...
        MyRes <- data.frame(
            trait_name = TraitName,
            total_phenotyped = length(phenotypedIndividuals(object)),
            total_affected = length(affectedIndividuals(object)),
            affected_id = NA,
            family = NA,
            affected = NA,
            ksgrp = NA,
            kinship_sum = NA,
            freq = NA,
            pvalue = NA,
            padj = NA,
            check.names = FALSE, stringsAsFactors = FALSE)
        warning("No simulation data available! This means that either no ",
                "simulation was run yet (using the kinshipClusterTest ",
                "function or runSimulation method) or that the simulation ",
                "returned no result (i.e. too few affected individuals in ",
                "the trait).")
        return(MyRes)
    }
    resList <- object@sim
    affIds <- as.character(names(resList$sumKinship))
    ped <- pedigree(object)
    ped <- ped[ped$id %in% affIds, c("id", "family")]
    pedL <- split(ped, f=ped$id)
    fams <- unlist(lapply(pedL, function(z){
        return(paste(unique(z[, "family"]), collapse=", "))
    }))
    fams <- fams[affIds]
    res <- data.frame(
        trait_name=TraitName,
        total_phenotyped=length(phenotypedIndividuals(object)),
        total_affected=length(affectedIndividuals(object)),
        affected_id=affIds,
        family=fams,
        affected=rep(length(resList$affected), length(resList$sumKinship)),
        ksgrp=NA,
        kinship_sum=resList$sumKinship,
        freq=resList$frequencyKinship,
        pvalue=resList$pvalueKinship,
        padj=p.adjust(resList$pvalueKinship, method=method),
        stringsAsFactors=FALSE)
    rownames(res) <- as.character(res$affected_id)
    res <- res[order(res$pvalue, -res$kinship_sum), , drop=FALSE]
    padj <- res$padj
    names(padj) <- rownames(res)
    res$ksgrp <- define.result.groups.ks(kinship(object), padj, th=cutoff,
                                         rmKinship=rmKinship)
    return(res)
})

##' Define groups of related affected individuals for KS test.
##' @param kin Kinship matrix.
##' @param padj Named vector with adjusted P values of all affected individuals.
##' @param th  Threshold of `padj` for inclusion of affected individuals.
##' @param rmKinship Ignore all pairs with kinship <= rmKinship.
##' @return A named vector of group numbers or NAs, the names correspond to the
##' IDs of the affected individuals. An NA indicates that the affected
##' individual was above the threshold. The named vector corresponds to the
##' entries in vector `padj`.
define.result.groups.ks <- function(kin, padj, th=0.05, rmKinship=0)
{
    stopifnot("Vector padj must have associated names."=!is.null(names(padj)))
    allIDs <- names(padj)
    thIDs <- allIDs[padj<th]
    grpNr  <- 1
    grp    <- rep(NA, length(allIDs))
    names(grp) <- allIDs
    while( length(thIDs)>0 ) {
        kinIDs <- intersect(
            doShareKinship(kin=kin, id=thIDs[1], rmKinship=rmKinship), allIDs)
        ## Don't overwrite previously assigned group numbers.
        grp[is.na(grp) & names(grp) %in% kinIDs] <- grpNr
        grpNr <- grpNr + 1
        thIDs <- setdiff(thIDs, kinIDs)
    }
    grp
}

## null hypothesis. H0: sum of kinship values of affected
## with all other affected is random. Test: is the sum of kinship values of
## affected with all other affected larger or equal that we would expect by
## chance (randomly selected individuals being affected)? The result is
## essentially identical to the significant.individuals2 function.
## NOTE: this function is robust against cases in which the observed kinship
## sum was never sampled in the simulation! The significant.individuals.2
## function would return a NA for the p-value and the frequency, while the
## p-value could be non-NA if larger kinship sums were sampled in the simulation!
## strata is supposed to be ordered the same way than pool
## tableSimVals: if TRUE a table of the simulated kinship values is returned.
significant.individuals.althyp = function(ks, affected, pool, nr.sim,
                                          strata=NULL, tableSimVals = FALSE){
    if( length(affected)<=1 ) {
        res = data.frame(id=NA, kc=NA, freq=NA, pval=NA);
        return(res);
    }
    if(!all(affected %in% colnames(ks))){
        stop("Got affected ids for which I don't have a kinship value! ",
             "This shouldn't happen!")
    }
    obs <- colSums(ks[affected, affected], na.rm=TRUE)
    names(obs) <- affected

    ## We compute small chunks of simulated data in order to keep the memory
    ## footprint low. This involves strategies to compute incremental tables,
    ## densities and histograms.
    nr.sim.todo <- nr.sim
    estimates <- matrix(0, nrow = length(obs), ncol = 2,
                        dimnames = list(names(obs), c("freq", "pval")))
    expDensity <- expHist <- expTable <- NULL
    simCount <- 0
    nr.inc.sim <- 1000
    while( nr.sim.todo>0 ) {
        ## Last simulation loop: run only the leftover number of steps.
        if( nr.sim.todo<nr.inc.sim )
            nr.inc.sim <- nr.sim.todo
        ## return the actual values from the simulation.
        if(is.null(strata)){
            simVals <- bg.ks.distr(ks, affected, pool, nr.inc.sim)
        }else{
            simVals <- bg.ks.distr.strat(ks, affected, pool, nr.inc.sim, strata)
        }
        ## Perform all the incremental steps.
        simCount <- simCount + length(simVals)
        estimates <- estimates + do.call(rbind,lapply(obs, function(z){
            oFreq <- sum(simVals == z)
            oPval <- sum(simVals >= z)
            return(c(freq = if( is.na(oFreq) ) 0 else oFreq,
                     pval = if( is.na(oPval) ) 0 else oPval))
        }))
        expDensity <- inc.density(simVals, expDensity)
        expHist <- inc.hist(simVals, expHist)
        if( tableSimVals )
            expTable <- inc.table(simVals, expTable)
        nr.sim.todo <- nr.sim.todo - nr.inc.sim
    }

    if( expHist$nreplace>0 )
        warning("Histogram approximation (overflow): replaced ",
                expHist$nreplace, " value(s) with highest possible value, ",
                expHist$breaks[length(expHist$breaks)])
    sumKin <- obs
    names(sumKin) <- affected
    freqKin <- estimates[, "freq"]*100/simCount
    names(freqKin) <- affected
    pKin <- estimates[, "pval"]/simCount
    names(pKin) <- affected
    res <- list(sumKinship=sumKin,
                pvalueKinship=pKin,
                frequencyKinship=freqKin,
                expDensity=expDensity,
                expHist=expHist,
                affected=affected)
    if (tableSimVals) {
        res <- c(res, list(tableSimVals = expTable))
    }

    return(res)
}


## Am just sampling using the index 1:nrow of kinship matrix.
bg.ks.distr = function(ks, affIds, pool, nr.sim){
    ## Assuming that kinship and pool are ordered the same way!
    n <- length(affIds)
    colIdx <- 1:ncol(ks)
    bg <- lapply(1:nr.sim, function(x){
        affidx <- sample(colIdx, n, replace=FALSE)
        return(colSums(ks[affidx, affidx], na.rm=TRUE))
    })
    return(unlist(bg, use.names=FALSE))
}

## same as bg.ks.distr, but for stratified sampling.
## Important is that kinship, pool and strata are in the same ordering!
bg.ks.distr.strat = function(ks, affIds, pool, nr.sim, strata=NULL){
    ksColnames <- colnames(ks)
    names(strata) <- as.character(pool)
    affStrataCounts <- table(strata[affIds])
    affStrataCounts <- affStrataCounts[affStrataCounts > 0]
    bg <- lapply(1:nr.sim, function(x){
        ## subsetting by index is faster...
        affidx <- stratsample(strata, counts=affStrataCounts)
        return(colSums(ks[affidx, affidx], na.rm=TRUE))
    })
    return(unlist(bg, use.names=FALSE))
}



#############
## plotting method...
## plotPed representing the results from the kinship clustering test.
## plotPed for FAKinSumResults does only support plotting for id.
setMethod("plotPed", "FAKinSumResults",
          function(object, id=NULL,
                   family=NULL, filename=NULL,
                   device="plot", only.phenotyped=FALSE, ...){
                       if(!is.null(family))
                           stop("Generating a pedigree for a family is not ",
                                "supported for FAKinshipResult. See help for ",
                                "more information.")
                       callNextMethod(object=object, id=id, family=family,
                                      filename=filename,
                                      device=device, proband.id=id,
                                      only.phenotyped=only.phenotyped, ...)
                       ## alternatively use highlight.ids
                   })


setMethod("[", "FAKinSumResults", function(x, i, j, ..., drop){
    stop("Subsetting of a FAKinSumResults object is not supported!")
})

## This will be a crazy funky method to plot the simulation results.
setMethod("plotRes", "FAKinSumResults",
          function(object, id=NULL,
                   family=NULL, addLegend=TRUE, type="density", ...){
              type <- match.arg(type, c("density", "hist"))
              if(length(object@sim) == 0)
                  stop("No analysis performed yet!")
              if(is.null(id))
                  stop("The id of the affected individual for whom the ",
                       "simulation results should be displayed has to be ",
                       "specified!")
              id <- as.character(id)
              ## result on family is not allowed.
              if(!is.null(family) | length(id) > 1)
                  stop("plotRes for FAKinSumResults does only support ",
                       "specifying a single id!")
              ## check if the id is an affected for which we tested...
              if(!any(names(object@sim$pvalueKinship) == id))
                  stop("No simulation result is available for individual ", id,
                       ". id should be one of result(object)$affected_id.")
              ## OK, now we can really start doing things...
              kinSum <- object@sim$sumKinship[id]
              kinPval <- object@sim$pvalueKinship[id]
              fam <- family(object, id=id)[1, "family"]
              par(xpd=FALSE)
              if(type == "density"){
                  toplot <- object@sim$expDensity
                  XL <- range(c(range(toplot$x), kinSum))
                  plot(toplot, main=paste0("Affected: ", id, ", family: ", fam),
                       xlab="Kinship sum", type="h",
                       lwd=3, col="lightgrey", xlim=XL)
                  points(toplot, col="grey", type="l", lwd=2)
              }
              if(type == "hist"){
                  toplot <- object@sim$expHist
                  XL <- range(c(range(toplot$mids), kinSum))
                  plot(toplot, main=paste0("Affected: ", id, ", family: ", fam),
                       xlab="Kinship sum",
                       col="lightgrey", border="grey", xlim=XL)
              }
              Blue <- "#377EB8"
              abline(v=kinSum, col=Blue)
              if(addLegend){
                  legend("topright", legend=c(paste0("kinship sum: ",
                                                     format(kinSum, digits=2)),
                                              paste0("p-value     : ",
                                                     format(kinPval, digits=3))
                                              ))
              }
          })


