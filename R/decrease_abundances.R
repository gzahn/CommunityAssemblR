#' Manipulates taxon or sample abundances in an arbitrary manner.
#'
#' Manipulates a community matrix with samples as rows and taxa as columns.
#'
#' @param dat Numeric matrix. The starting community matrix that will be edited.
#' @param prop Positive numeric vector of length 1, between 0,1. The proportion of taxa affected by the arbitrary abundance changes.
#' @param reduction.scale Positive numeric vector of length 1 and >=0 and <=1. The proportionate amount to decrease a given taxon or sample.
#' @param margin Character vector of length 1. Must be either 'taxa' or 'samples'. Partial matching not allowed. Determines whether the gradient increases or decreases taxa abundances. Default = 'negative'.
#' @param transplant.only Logical vector of length 1. If TRUE, the gradient will only affect new taxa that persisted after transplantation. This only makes sense if performing this function on a matrix that is the result of a previous transplantation step, such as a matrix provided by transplant_w_facilitation(). If margin = "samples" then this is ignored. Default = FALSE.
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'. Random taxon or random sample abundance values will be decreased according to parameters selected.
#'
#' @examples
#' comm <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#' decrease_abundance(comm,prop = .5,reduction.scale = .95,margin = "samples")
#'
#' @export

decrease_abundances <- function(dat, # a matrix, samples are rows and columns are taxa
                             prop = .1, # needs to be between 0,1
                             reduction.scale = .9, # needs to be between 0,1
                             transplant.only=FALSE, # if TRUE, only newly introduced taxa from a transplantation will be affected
                             margin = "taxa" # needs to be "taxa" or "samples"
                             ){
  # tests ####
  stopifnot("matrix" %in% class(dat))
  stopifnot(class(prop) %in% c("numeric","integer") & prop >= 0 & prop <= 1)
  stopifnot(class(reduction.scale) %in% c("numeric","integer") & reduction.scale >= 0 & reduction.scale <= 1)
  stopifnot(margin %in% c("taxa","samples"))
  # convert some proportion of taxa to be rarely seen

  if(margin == "taxa"){
    # setup ####

    newtaxa.names <- colnames(dat)[grepl("^newtaxon_",x=colnames(dat))]

    if(transplant.only){
      n.rare.taxa <- round(prop * length(newtaxa.names))
      rare.taxa <- sample(newtaxa.names,n.rare.taxa)
    }
    if(!transplant.only){
      n.rare.taxa <- round(ncol(dat) * prop)
      rare.taxa <- sample(colnames(dat),n.rare.taxa,replace = FALSE)
    }


    reduce_counts <- function(x){round(x[,rare.taxa] * (1-reduction.scale))}

    dat[,rare.taxa] <- reduce_counts(dat)
    return(dat)
  }

  if(margin == "samples"){
    # setup ####
    n.rare.samples <- round(nrow(dat) * prop)
    rare.samples <- sample(rownames(x),n.rare.samples,replace = FALSE)
    reduce_counts <- function(x){round(x[rare.samples,] * (1-reduction.scale))}

    dat[rare.samples,] <- reduce_counts(dat)
    return(dat)
  }
}



