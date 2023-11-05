#' Build even community matrix
#'
#' Builds a community matrix with samples as rows and taxa as columns.
#' Taxa abundances are drawn from a normal distribution for each sample (this is currently the only distribution supported).
#'
#'
#'
#' @param n.taxa Positive numeric vector of length 1. The number of taxa (columns) to generate. Default = 5.
#' @param n.samples Positive numeric vector of length 1. The number of samples (rows) to generate. Default = 5.
#' @param n.reads Positive numeric vector of length 1. The number of observations within each sample. This is analogous to per-sample sequencing depth.
#' @param taxa.dist Currently, only the 'normal' distribution is accepted.
#' @param taxa.sd The standard deviation for the randomization function that generates taxa abundances. Higher values generate more zeros for each taxon.
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'
#'
#' @examples
#' comm <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#'
#' @export

build_even_community <-
  function(n.taxa = 5, # must be positive integer (or coerceable)
           n.samples = 5, # must be positive integer (or coerceable)
           n.reads = 5000, # must be positive integer (or coerceable)
           taxa.dist = "normal", # type of species distribution with a sample
           taxa.sd = 0 # standard deviation for distributions
           ){
    # tests
    stopifnot(class(n.taxa) %in% c("numeric","double","integer") & n.taxa > 0)
    stopifnot(class(n.samples) %in% c("numeric","double","integer") & n.samples > 0)
    stopifnot(class(n.reads) %in% c("numeric","double","integer") & n.reads > 0)
    stopifnot(taxa.dist %in% c("normal") & class(taxa.dist) == "character")
    stopifnot(class(taxa.sd) %in% c("numeric","logical","integer"))

    # setup
    # n.readspersample = n.reads / n.samples


    # internal functions ####

    ## for normal species abundances ####
    unif_vect <- function(N, M, sd = 1, pos.only = TRUE) {
      vec <- rnorm(N, M/N, sd)
      if (abs(sum(vec)) < 0.01) vec <- vec + 1
      vec <- round(vec / sum(vec) * M)
      deviation <- M - sum(vec)
      for (. in seq_len(abs(deviation))) {
        vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
      }
      if (pos.only) while (any(vec < 0)) {
        negs <- vec < 0
        pos  <- vec > 0
        vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
        vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
      }
      vec
    }


    # make "otu table" ####
    sample.abundances <- list()

    community <-
      if(taxa.dist == "normal"){
        ## for uniform species abundances ####
        for(i in 1:n.samples){
          sample.abundances[[i]] <- unif_vect(N=n.taxa,M=n.reads,sd = taxa.sd)
        }

        mat <- matrix(unlist(sample.abundances),nrow = n.samples,byrow = TRUE)
        # mat <- t(mat)

        row.names(mat) <- paste0("sample_",1:n.samples)
        colnames(mat) <- paste0("taxon_",1:n.taxa)
        mat
      }

    # return community table
    return(community)
}

