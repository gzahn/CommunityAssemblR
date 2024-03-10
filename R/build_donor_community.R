#' Build a community matrix to simulate community transplantation
#'
#' Builds a community matrix with samples as rows and taxa as columns.
#' The number of samples will match the starting (resident) community provided.
#' Taxon abundances will be on a similar scale to the provided resident community.
#'
#' @param resident.comm Numeric matrix. The starting community matrix that this one will be based on. Typically, the result of build_even_community() or a similar starting community.
#' @param n.transplant.taxa Positive numeric vector of length 1. The number of taxa in the donor community. This represents the microbiome community that will be transplanted into the resident community.
#' @param overlap Positive numeric vector of length 1 and >=0 and <=1. The proportion (between 0,1) of donor taxa that will also be present in the recipient community.
#' @param prop.dominant The proportion of taxa that are dominant in the donor community; percent relative abundance. Default = 0; suggesting a cultured community, not a wild one.
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'. Non-overlapping taxa will have names beginning with "newtaxon_".
#'
#' @examples
#' comm <- build_even_community(n.taxa = 3000, n.samples = 100, n.reads = 500000, prop.dominant = .02)
#' donor <- build_donor_community(resident.comm = comm, n.transplant.taxa = 30,overlap = .75)
#'
#' @export

build_donor_community <- function(resident.comm,n.transplant.taxa,overlap=0,prop.dominant=0){

  # set up parameters for build_even_community
  n.taxa <- round(n.transplant.taxa)
  stopifnot(class(n.transplant.taxa) %in% c("numeric","integer") & n.transplant.taxa > 0)
  stopifnot("matrix" %in% class(resident.comm))
  n.samples <- nrow(resident.comm)
  stopifnot(class(overlap) %in% c("numeric","integer") & overlap >= 0 & overlap <= 1)
  n.reads <- round(sum(resident.comm) / ncol(resident.comm),0) * n.taxa
  # prop.dominant <- .2 # the default assumes, essentially, a 'cultured' community, not a wild one. Could change in future
  x <- resident.comm
  x[x > 0] <- 1
  prop.absent <- mean((nrow(x)-colSums(x)) / nrow(x),na.rm=TRUE)
  message("Assuming that samples are rows / taxa are columns...")
  message(paste0("Building new 'donor' community with ",overlap*100,"% overlap in community membership."))


  # build basic community of 'transplant' microbes
  transplant.comm <-
    build_even_community(n.taxa,n.samples,n.reads,prop.dominant = prop.dominant,prop.absent = prop.absent)


  if(overlap == 0){
    # make ALL distinct taxon names
    colnames(transplant.comm) <- sub(pattern = "taxon_",replacement = "newtaxon_",colnames(transplant.comm))
  }


  if(overlap > 0){
    # make SOME distinct taxon names
    n.to.rename <- round(n.taxa) * (1-overlap)
    taxa.to.rename <- sample(colnames(transplant.comm),n.to.rename,replace = FALSE)
    taxa.to.rename
    all.names <- colnames(transplant.comm)
    new.names <- sub(pattern = "taxon_",replacement = "newtaxon_",colnames(transplant.comm[,taxa.to.rename]))
    all.names[all.names %in% taxa.to.rename] <- new.names
    colnames(transplant.comm) <- all.names
  }
  return(transplant.comm)
}

