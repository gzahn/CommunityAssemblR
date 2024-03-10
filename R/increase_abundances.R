#' Manipulates taxon or sample abundances in an arbitrary manner.
#'
#' Manipulates a community matrix with samples as rows and taxa as columns.
#'
#' @param dat Numeric matrix. The starting community matrix that will be edited.
#' @param prop Positive numeric vector of length 1, between 0,1. The proportion of taxa affected by the arbitrary abundance changes.
#' @param increase.scale Positive numeric vector of length 1 and >=0 and <=1. The proportionate amount to increase a given taxon or sample.
#' @param margin Character vector of length 1. Must be either 'taxa' or 'samples'. Partial matching not allowed. Determines whether the gradient increases or decreases taxa abundances. Default = 'negative'.
#' @param transplant.only Logical vector of length 1. If TRUE, the gradient will only affect new taxa that persisted after transplantation. This only makes sense if performing this function on a matrix that is the result of a previous transplantation step, such as a matrix provided by transplant_w_facilitation(). If margin = "samples" then this is ignored. Default = FALSE.
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'. Random taxon or random sample abundance values will be increased according to parameters selected.
#'
#' @examples
#' comm <- build_even_community(n.taxa = 100, n.samples = 44, n.reads = 3000, taxa.sd = 30)
#' increase_abundance(comm,prop = .5, increase.scale = .95, margin = "taxa")
#'
#' @export


increase_abundances <- function(dat, # a matrix, samples are rows and columns are taxa
                                prop = .1, # needs to be between 0,1
                                increase.scale = .9,
                                margin = "taxa", # needs to be "taxa" or "samples"
                                transplant.only = FALSE){
  # tests ####
  stopifnot("matrix" %in% class(dat))
  stopifnot(class(prop) %in% c("numeric","integer") & prop >= 0 & prop <= 1)
  stopifnot(class(increase.scale) %in% c("numeric","integer") & increase.scale >= 0)
  stopifnot(margin %in% c("taxa","samples"))
  # convert some proportion of taxa to be highly abundant

  if(margin == "taxa"){

    newtaxa.names <- colnames(dat)[grepl("^newtaxon_",x=colnames(dat))]

    if(transplant.only){
      n.abund.taxa <- round(prop * length(newtaxa.names))
      abund.taxa <- sample(newtaxa.names,n.rare.taxa)
    }
    if(!transplant.only){
      n.abund.taxa <- round(ncol(dat) * prop)
      abund.taxa <- sample(colnames(dat),n.abund.taxa,replace = FALSE)
    }
    # setup ####
    increase_counts <- function(x){round(x[,abund.taxa] * (1+increase.scale))}

    dat[,abund.taxa] <- increase_counts(dat)
    return(dat)
  }

  if(margin == "samples"){
    # setup ####
    n.abund.samples <- round(nrow(dat) * prop)
    abund.samples <- sample(rownames(dat),n.abund.samples,replace = FALSE)
    increase_counts <- function(x){round(dat[abund.samples,] * (1+increase.scale))}

    dat[abund.samples,] <- increase_counts(dat)
    return(dat)
  }
}

