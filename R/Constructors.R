## functions related to FunSim.

##********************************************************
##
##   FAData
##   that's the constructor for the FunAggrData class
##   * pedigree: the pedigree, can be a character string or a
##     data.frame
##   * age: the age, can
##   * header: for the internal read.table function that loads data.
##   * sep:
##
##********************************************************
FAData <- function(pedigree, age, trait, traitName, header=FALSE, sep="\t", id.col="id",
                   family.col="family", father.col="father", mother.col="mother",
                   sex.col="sex"){
    if(missing(pedigree)){
        stop("Argument pedigree has to be specified!")
    }
    FSD <- new("FAData")
######################
    ## processing pedigree
    if(is.character(pedigree)){
        ## means, we're loading the data.frame ourselfs.
        message(paste("Reading pedigree information from file ", pedigree, "..."),
                appendLF=FALSE)
        pedigree <- read.table(pedigree, header=header, sep=sep)
        ## if I've got a header, check and subset to columns...
        if(header){
            ## check if I have all required column names...
            if(!all(c( id.col, family.col, father.col, mother.col, sex.col) %in%
                    colnames(pedigree)))
                stop(paste0("One or more required column names not found! Expect: ",
                            family.col, ", ", id.col, ", ", father.col, ", ", mother.col, ", ",
                            sex.col, "!"))
            pedigree <- pedigree[, c(family.col, id.col, father.col, mother.col, sex.col)]
            colnames(pedigree) <- .PEDCN
        }
        message("OK\n")
    }
    if(is.matrix(pedigree)){
        ## transform to data.frame
        pedigree <- as.data.frame(pedigree, stringsAsFactors=FALSE)
    }
    ##    if(!(class(pedigree) %in% c("data.frame", "pedigree", "pedigreeList"))){
    if(!(is.data.frame(pedigree) | is(pedigree, "pedigree") | is(pedigree, "pedigreeList"))){
        stop(paste0("Argument pedigree has to be a data.frame, a character string specifying the",
                    " file containing the pedigree, a pedigree object or a pedigreeList object!"))
    }else{
        if(is.data.frame(pedigree)){
            CN <- .PEDCN
            if(!noColNames(pedigree)){
                ## check if we have all column names required.
                if(sum(CN %in% colnames(pedigree)) != length(CN)){
                    stop("pedigree is expected to have column names ", paste(CN, collapse=", "))
                }
                pedigree <- pedigree[, CN]
            }else{
                if(ncol(pedigree)!=length(CN))
                    stop("pedigree is expected to have ", length(CN), " columns but I got",
                         ncol(pedigree), "!")
                warning("Forcing column names ", paste(CN, collapse=";"), " on pedigree.")
                colnames(pedigree) <- CN
            }
        }
    }
    pedigree(FSD) <- pedigree
    ## processing age
    if(!missing(age)){
        if(is.character(age)){
            message(paste("Reading age information from file ", age, "..."), appendLF=FALSE)
            FN <- age
            age <- read.table(age, header=header, sep=sep, as.is=TRUE)
            message("OK\n")
            ## assuming first column is the ID, second the age.
            warning("Assuming first column in ", FN, " contains IDs and second age.")
            ageNames <- age[, 1]
            age <- age[ , 2]
            names(age) <- ageNames
        }
        if(!is.numeric(age)){
            stop("Argument age has to be a named numeric vector or a character string specifying the file containing the age!")
        }
        age(FSD) <- age
    }
    validateFAData(FSD)
    if(!missing(traitName))
        FSD@traitname <- traitName
    if(!missing(trait)){
        if(length(trait) == nrow(pedigree)){
            if(is.null(names(trait)))
                names(trait) <- as.character(pedigree$id)
        }
        trait(FSD) <- trait
    }
    return(FSD)
}


