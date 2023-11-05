#' Force increased similarity between random groups of samples in terms of taxa counts
#'
#' Manipulates a community matrix with samples as rows and taxa as columns.
#'
#' @param dat Numeric matrix. The starting community matrix that will be edited.
#' @param n.links Positive numeric vector of length 1. The number of groupings to force. Must be >= 1 and < sample.number/group.size
#' @param group.size Positive numeric vector of length 1. Determines how many samples to make similar in each grouping. Must be >= 2 and <= sample.number/n.links
#' @param link.strength Numeric vector of length 1. how strong should the similarity function be. 0=no change; 1=make identical.
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'. Random sample abundance values will be altered to show increased similarity according to parameters selected.
#'
#' @examples
#' comm <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#' link_sample_communities(comm, n.links = 3, group.size = 4, link.strength = 0.25)
#'
#' @export


link_sample_communities <- function(dat,
                                    n.links = 1,
                                    group.size = 3,
                                    link.strength = .9){
  sample.number = nrow(dat)
  group.size <- round(group.size)

  stopifnot(n.links >= 1 & class(n.links) == "numeric" & n.links <= sample.number/group.size)
  stopifnot(group.size >= 2 & class(group.size) == "numeric" & group.size <= sample.number/n.links)
  stopifnot(link.strength <= 1 & link.strength >= 0 & class(link.strength) %in% c("numeric","integer"))

  for(i in 1:n.links){
  # pick random samples to link up
  grouped.samples <- sample(row.names(dat),group.size)
  first.sample <- grouped.samples[1]
  other.samples <- grouped.samples[-1]

  # get difference values between first and other samples
  # use link.strength to make them match first sample
  for(i in seq_along(other.samples)){
    sample.diff <- dat[first.sample,] - dat[other.samples[i],]
    dat[other.samples[i],] <- dat[other.samples[i],] + round(sample.diff * link.strength)
    }
  message(paste0("Samples: ",paste(other.samples,collapse = ", ")," were made to resemble Sample: ",first.sample))
  }
  return(dat)
}

