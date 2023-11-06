#' Simulates the actions of niche overlap and differential fitness to alter a taxa abundance table.
#'
#' Manipulates a community matrix with samples as rows and taxa as columns. Taxa are assigned random ranges within a number of niche distributions. Taxa are also assigned random "fitness" values. Each possible pairing of taxa are matched and competition is simulated based on any niche overlap and differential fitness. Taxon abundances are altered accordingly for both taxa.
#'
#' @param dat Numeric matrix. The starting community matrix that will be edited.
#' @param n.niches Positive whole numeric vector of length 1. The number of niches to simulate. All taxa will be assigned to one of these niches. Default=ncol(dat)/2.
#' @param niche.distr Character vector of length 1. Must currently be in c("normal"). The distribution from which to draw niche overlap values. Default='normal'.
#' @param niche.separation Numeric vector of length 1. Indicates the standard deviation for the initial seeds that generate niches. Larger values lead to less niche overlap. Default=n.niches^2.
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'. Abundance values are updated to reflect all competitive outcomes.
#'
#' @examples
#' comm <- build_even_community(n.taxa = 100, n.samples = 44, n.reads = 3000, taxa.sd = 30)
#' sim <- simulate_competition(comm, n.niches=15, niche.separation=5)
#' colSums(comm);colSums(sim)
#'
#' @export

simulate_competition <- function(dat,n.niches=ncol(dat)/2,niche.distr="normal",niche.separation=n.niches^2){

  #setup
  n.niches <- ifelse(round(n.niches) < 1, 1,round(n.niches))

  #tests
  stopifnot("matrix" %in% class(dat))
  stopifnot(class(n.niches) %in% c("numeric","integer") & n.niches > 0)
  stopifnot(niche.distr %in% c("normal"))
  stopifnot(class(niche.separation) %in% c("numeric","integer") & niche.separation > 0)


  # vector of niche 'names'
  niche.names <- paste0("niche_",1:n.niches)

  # select random parameters for niche distributions
  # THESE UNIFORM DISTRIBUTIONS DON'T WORK WELL YET!!!!!!!!!!!!!!!!!!
  unif.dist.mins <- -abs(rnorm(n.niches,sd=niche.separation))
  unif.dist.maxs <- abs(rnorm(n.niches,sd=niche.separation))

  # for normal dist
  # pick N evenly spaced values as means
  norm.dist.means <- seq(niche.separation,n.niches*(0.1+niche.separation),length.out=n.niches)
  norm.dist.sds <- rep(n.niches,n.niches)

  # add niche names to niche parameters
  names(norm.dist.means) <- niche.names
  names(norm.dist.sds) <- niche.names
  names(unif.dist.mins) <- niche.names
  names(unif.dist.maxs) <- niche.names


  # assign each taxon to a niche
  niche.assignments <- sample(niche.names,length(colnames(dat)),replace = TRUE)
  names(niche.assignments) <- colnames(dat)

  # build out niche ranges for each taxon

  # for normally distributed niches
  if(niche.distr == "normal"){
    # for-loop to get each taxon's range in the "niche landscape"
    ranges <- sapply(colnames(dat), function(x) NULL)
    for(i in colnames(dat)){
      # get two values from that distribution
      x <-
        rnorm(2,
              mean = norm.dist.means[niche.assignments[i]],
              sd = norm.dist.sds[niche.assignments[i]])
      # arrange them smallest to largest
      x <- x[order(x)]
      # assign to the list, iteratively
      ranges[[i]] <- x
    }
  }

  # for uniformly distributed niches
  if(niche.distr == "uniform"){
    # for-loop to get each taxon's range in the "niche landscape"
    ranges <- sapply(colnames(dat), function(x) NULL)
    for(i in colnames(dat)){
      # get two values from that distribution
      x <-
        runif(2,
              min = unif.dist.mins[niche.assignments[i]],
              max = unif.dist.maxs[niche.assignments[i]])
      # arrange them smallest to largest
      x <- x[order(x)]
      # assign to the list, iteratively
      ranges[[i]] <- x
    }
  }

  # function to calcluate range overlap...if it exists
  calculate_range_overlap <- function(start1, end1, start2, end2) {
    # Calculate the length of the overlapping region
    overlap_length <- min(end1, end2) - max(start1, start2)

    # Check for non-overlapping ranges
    if (overlap_length < 0) {
      overlap_length <- 0
    }

    return(overlap_length)
  }

  # for-loop to find niche overlap for all possible pairings
  pairing <- c()
  for(i in names(ranges)){
    x <- ranges[[i]]
    for(j in names(ranges)[names(ranges) != i]){
      y <- ranges[[j]]

      pairing.name <- paste(i,j,sep = "|")
      z <- calculate_range_overlap(x[1],x[2],y[1],y[2])
      pairing[pairing.name] <- z
    }
  }

  # compile matches in data frame
  match1 <- unlist(strsplit(names(pairing),split = "\\|"))[seq(1,(2*length(names(pairing))),by=2)]
  match2 <- unlist(strsplit(names(pairing),split = "\\|"))[seq(2,(2*length(names(pairing))),by=2)]
  matches <- data.frame(overlap=pairing,pairing=names(pairing),match1,match2)

  # assign fitness values to all taxa
  taxa.fitnesses <- runif(length(colnames(dat)))
  names(taxa.fitnesses) <- colnames(dat)
  matches$fit.1 <- taxa.fitnesses[matches$match1]
  matches$fit.2 <- taxa.fitnesses[matches$match2]
  matches$winner <- ifelse(matches$fit.1 > matches$fit.2,matches$match1,matches$match2)
  matches$loser <- ifelse(matches$fit.1 < matches$fit.2,matches$match1,matches$match2)

  # scale niche overlap to percentage
  scale01 <- function(x){
    if(max(x) == 0){return(x)}
    if(max(x) > 0){return((x - min(x)) / (max(x) - min(x)))}
  }
  matches$overlap <- scale01(matches$overlap)

  # Now...use this lookup table to build successive pairings out of the abundance matrix and transform the values
  for(i in 1:nrow(matches)){
    mat_pair <- dat[,c(matches$match1[i],matches$match2[i])]
    mat_pair[,matches$winner[i]] <- ifelse(mat_pair[,matches$winner[i]] == 0, 0,
                                           mat_pair[,matches$winner[i]] + (mat_pair[,matches$loser[i]] * matches$overlap[i]))
    mat_pair[,matches$loser[i]] <- mat_pair[,matches$loser[i]] - (mat_pair[,matches$loser[i]] * matches$overlap[i])
    dat[,c(matches$match1[i],matches$match2[i])] <- round(mat_pair)
  }

  return(dat)
}
