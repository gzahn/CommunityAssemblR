#' Simulates abiotic filtration of taxa along a gradient.
#'
#' Manipulates a community matrix with samples as rows and taxa as columns.
#' Simulates an environmental gradient and adjusts a proportion of taxon abundances along it.
#' Multiple groups can be provided so that the gradient is replicated, simulating experimental replication.
#' Where the gradient value is 0, no change is made to taxon abundances.
#'
#' @param dat Numeric matrix. The starting community matrix that will be edited.
#' @param prop Positive numeric vector of length 1, between 0,1. The proportion of taxa affected by the abiotic gradient.
#' @param gradient.min Positive numeric vector of length 1 and >=0 and <=1. The minimum starting value of the abiotic gradient.
#' @param gradient.max Positive numeric vector of length 1 and >=0 and <=1. The minimum ending value of the abiotic gradient.
#' @param gradient.strength Positive numeric vector of length 1. A multiplier for abundance transformations along the gradient. Values between 0 and 1 will reduce the environmental effect.
#' @param groups Positive whole number vector of length 1. Indicates how many gradient simulations to divide samples into. If groups is not a factor of the number of samples in the matrix, the remaining samples will receive a partial gradient adjustment.
#' @param association Character vector of length 1. Must be either 'positive' or 'negative'. Partial matching not allowed. Determines whether the gradient increases or decreases taxa abundances. Default = 'negative'.
#' @param transplant.only Logical vector of length 1. If TRUE, the gradient will only affect new taxa that persisted after transplantation. This only makes sense if performing this function on a matrix that is the result of a previous transplantation step, such as a matrix provided by transplant_w_facilitation(). Default = FALSE.
#' @param show.gradient Logical vector of length 1. If TRUE, a simple plot of the gradient strength along samples in the matrix will be returned as a side effect.
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'. Taxon abundance values will be reduced or increased along a gradient according to parameters selected.
#'
#' @examples
#' comm <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#' donor <- build_donor_community(resident.comm = comm, n.transplant.taxa = 30,overlap = .25)
#' transplanted <- transplant_w_niche_occupation(comm, donor, niche.size.donor = 2, niche.size.resident = 3, process.type = "permission")
#' filter_taxa_along_gradient(transplanted, prop=.5, groups=2, transplant.only=TRUE)
#'
#' @export


filter_taxa_along_gradient <- function(dat, # community matrix with rows=samples & cols=taxa (raw abundance counts)
                                       prop=.1,  # proportion of taxa affected by this gradient
                                       gradient.min=0, # minimum gradient value (strength)
                                       gradient.max=1, # maximum gradient value (strength)
                                       gradient.strength=1, # multiplier for effect of gradient on taxa abundances
                                       groups=1, # N to split samples into (3*groups must be < nrow(dat)/groups)
                                       # remainder samples not within main group size will begin a new group
                                       association='negative', ## does the gradient increase or decrease abundance of affected taxa?
                                       transplant.only=FALSE, # if TRUE, only newly introduced taxa from a transplantation will be affected
                                       show.gradient=FALSE){ # if TRUE, a plot showing fitted gradient levels will be shown

  stopifnot("matrix" %in% class(dat))
  stopifnot(class(prop) %in% c("numeric","integer") & prop >= 0 & prop <= 1)
  stopifnot(class(gradient.min) %in% c("numeric","integer") & gradient.min >= 0 & gradient.min <= 1)
  stopifnot(class(gradient.max) %in% c("numeric","integer") & gradient.max >= 0 & gradient.max <= 1)
  stopifnot(class(groups) %in% c("numeric","integer") & nrow(dat)/groups >= 3)
  stopifnot(association %in% c("positive","negative"))
  stopifnot(class(gradient.strength) %in% c("numeric","integer") & gradient.strength >= 0)

  # get vector(s) for assigning gradient values to
  v <- 1:nrow(dat)
  main.group.length <- length(v) %/% groups
  remainder.group.length <- length(v) %% groups

  # build gradient sequence
  grad <- seq(gradient.min,gradient.max,length.out=main.group.length)
  main.grad <- rep(grad,(groups))
  remainder.grad <- grad[1:remainder.group.length]
  if(remainder.group.length == 0){
    gradient.values <- c(main.grad)
  }
  if(remainder.group.length > 0){
    gradient.values <- c(main.grad,remainder.grad)
    warning("Number of samples is not a multiple of the number of groups. Extra partial gradient group has been added.")
  }


  # add a touch of random noise
  gradient.values <- gradient.values + rnorm(length(gradient.values),0,.1)
  gradient.values <- sapply(gradient.values,FUN = function(x){ifelse(x < 0,0,x)})
  gradient.values <- sapply(gradient.values,FUN = function(x){ifelse(x > 1,1,x)})

  if(show.gradient){
    plot(gradient.values,main = "fitted gradient")
  }

  # select taxa to reduce
  newtaxa.names <- colnames(dat)[grepl("^newtaxon_",x=colnames(dat))]

  if(transplant.only){
    n.taxa <- round(prop * length(newtaxa.names))
    affected.taxa <- sample(newtaxa.names,n.taxa)
  }
  if(!transplant.only){
    n.taxa <- round(prop * ncol(dat))
    affected.taxa <- sample(colnames(dat),n.taxa)
  }

  # get matrix of just affected taxa and transform
  if(association == 'positive'){
    dat[,affected.taxa] <- round(dat[,affected.taxa] * (1+gradient.values) * gradient.strength)
  }

  if(association == 'negative'){
    dat[,affected.taxa] <- round(dat[,affected.taxa] * (1-(gradient.values^(1/gradient.strength))))
  }

  return(dat)

}
