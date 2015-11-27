\name{getAll}
\Rdversion{1.1}
\alias{getAll}
\alias{getAll,data.frame-method}
\alias{getAll,pedigree-method}
\alias{getAll,pedigreeList-method}
\alias{getAll,FAData-method}
\alias{getExternalMatched}
\alias{getExternalMatched,data.frame-method}
\alias{getExternalMatched,pedigree-method}
\alias{getExternalMatched,pedigreeList-method}
\alias{getExternalMatched,FAData-method}
\alias{getGenerationMatched}
\alias{getGenerationMatched,data.frame-method}
\alias{getGenerationMatched,pedigree-method}
\alias{getGenerationMatched,pedigreeList-method}
\alias{getGenerationMatched,FAData-method}
\alias{getGenerationSexMatched}
\alias{getGenerationSexMatched,data.frame-method}
\alias{getGenerationSexMatched,pedigree-method}
\alias{getGenerationSexMatched,pedigreeList-method}
\alias{getGenerationSexMatched,FAData-method}
\alias{getSexMatched}
\alias{getSexMatched,data.frame-method}
\alias{getSexMatched,pedigree-method}
\alias{getSexMatched,pedigreeList-method}
\alias{getSexMatched,FAData-method}

\title{Define sets of control individuals for one or more given individuals}

\description{
  These functions allow to define sets of (eventually matched) control
  individuals for one or more given individuals using pedigree data. By
  default, controls from the same family than the specified
  individual(s) are returned.
  They are used e.g. in the \code{\link{genealogicalIndexTest}} method.
}
\usage{

\S4method{getAll}{data.frame}(object, id=NULL, ...)

\S4method{getExternalMatched}{data.frame}(object, id=NULL, match.using, ...)

\S4method{getGenerationMatched}{data.frame}(object, id=NULL,
                                           include.anc=0, include.off=0, ...)

\S4method{getGenerationSexMatched}{data.frame}(object, id=NULL,
                                               include.anc=0, include.off=0, ...)

\S4method{getSexMatched}{data.frame}(object, id=NULL, ...)

}
\arguments{

  (in alphabetic order)

  \item{id}{
    A character vector with ids of individuals for whom matched
    control individuals (from the same family) should be defined.
  }

  \item{include.anc}{
    For \code{getGenerationMatched} and \code{getGenerationSexMatched}:
    number of ancestor generations of individuals \code{id} in which
    control individuals should be defined, in addition to the actual
    generation of the individuals \code{id}.
    By default the functions define control individuals only in the
    actual generation of the individuals in \code{id}.
  }

  \item{include.off}{
    For \code{getGenerationMatched} and \code{getGenerationSexMatched}:
    number of offspring generations of individuals \code{id} in which
    control individuals should be defined, in addition to the actual
    generation of the individuals \code{id}.
    By default the functions define control individuals only in the
    actual generation of the individuals in \code{id}.
  }

  \item{match.using}{
    For \code{getExternalMatched}: a named vector for matching of
    controls. The names of the vector have to correspond to the ids of
    the individuals in the pedigree (but can have a different ordering).
  }

  \item{object}{
    Either a \code{data.frame}, \code{pedigree} or a \code{pedigreeList}
    object specifying the pedigree. If a \code{data.frame} is submitted,
    the columns \code{id}, \code{family}, \code{father}, \code{mother} and
    \code{sex} are required.
  }

  \item{...}{
    Additional arguments; not used at present.
  }

}
\section{Get matched controls}{

  \describe{

    \item{getAll}{
      Simply returns the ids of all individuals in the family
      (i.e. individuals with the same family id in the pedigree) of the
      specified ids in \code{id}.
      Returns a list with the ids (of type \code{character}) of the
      controls. The names of the list correspond to the family id.
    }

    \item{getExternalMatched}{
      Returns the ids of matched individuals from the same family for
      the specified ids in \code{id}.
      The \code{match.using} vector is for the matching, i.e. the
      function first extracts the values for the individuals in
      \code{id} from \code{match.using} and returns the ids of all
      individuals whose value in \code{match.using} matches the value of
      the individuals specified by \code{id}. Individuals with a missing
      value in \code{match.using} are excluded.
      Returns a list with the ids (of type \code{character}) of the
      controls. The names of the list correspond to the family id.
    }

    \item{getGenerationMatched}{
      Returns the ids of individuals matched by the generation of the
      individuals in \code{id}. The function returns the ids of all
      individuals from the same generation(s) than the individuals in
      \code{id}. The arguments \code{include.anc} and \code{include.off}
      can be used to increase the range of generation from which
      individuals are selected.
      Returns a list with the ids (of type \code{character}) of the
      controls. The names of the list correspond to the family id.
    }

    \item{getGenerationSexMatched}{
      Same as \code{getGenerationMatched}, but matches in addition
      individuals by sex (see \code{getSexMatched}).
      Returns a list with the ids (of type \code{character}) of the
      controls. The names of the list correspond to the family id.
    }

    \item{getSexMatched}{
      Returns ids of individuals from the same family matching the sex
      of the individuals in \code{id}. Individuals with sex being
      \code{NA} are excluded. The result is essentially identical to the
      \code{getAll} if \code{id} contains ids of male and female
      individuals.
      Returns a list with the ids (of type \code{character}) of the
      controls. The names of the list correspond to the family id.
    }

  }
}
\value{
  Refer to the method and function description above for detailed
  information on the returned result object.
}
\author{
  Johannes Rainer.
}
\seealso{
  \code{\link[kinship2]{pedigree}},
  \code{\linkS4class{FAData}},
  \code{\linkS4class{FAProbResults}},
  \code{\linkS4class{FAKinGroupResults}},
  \code{\linkS4class{FAKinSumResults}},
  \code{\linkS4class{FAGenIndexResults}},
  \code{\link{genealogicalIndexTest}}
}
\examples{
##########################
##
##  Defining a small pedigree
##
## load the Minnesota Breast Cancer record and subset to the
## first families.
data(minnbreast)
mbsub <- minnbreast[minnbreast$famid==4 | minnbreast$famid==5, ]
mbped <- mbsub[, c("famid", "id", "fatherid", "motherid", "sex")]
## renaming column names
colnames(mbped) <- c("family", "id", "father", "mother", "sex")

## Get ids of individuals from the same generation than "4"
getGenerationMatched(mbped, id="4")

## Get the ids of gemale individuals from the same generation than
## "4" and "22"
getGenerationSexMatched(mbped, id=c("4", "22"))

## Generate a FAData
fad <- FAData(mbped)

## Plot the pedigree so we can evaluate the results of the function
plotPed(fad, family="4")

## We're using getExternalMatched on the sex of the individuals
## in the pedigree, thus the results are identical to the
## getSexMatched function.
## Extracting the sex using $sex returns a named vector just
## as we need for getExternalMatched.
head(fad$sex)
getExternalMatched(fad, id="4", match.using=fad$sex)

getSexMatched(fad, id="4")

}
\keyword{classes}

