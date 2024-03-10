#' Build even community matrix
#'
#' Builds a community matrix with samples as rows and taxa as columns.
#' Taxa abundances are distributed in decreasing value after the dominant taxa are apportioned their shares.
#'
#' @param n.taxa Positive numeric vector of length 1. The number of taxa (columns) to generate. Default = 500.
#' @param n.samples Positive numeric vector of length 1. The number of samples (rows) to generate. Default = 100.
#' @param n.reads Positive numeric vector of length 1. The number of reads (observations) for the entire study. This is analogous to overall sequencing depth. Default = 5,000,000
#' @param prop.dominant The proportion of taxa that are dominant; represented as > 1 percent relative abundance. Default = .01
#' @param prop.absent The rough proportion of taxa that are missing from a given sample. Default = 0.2
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'
#'
#' @examples
#' comm <- build_even_community(n.taxa = 3000, n.samples = 100, n.reads = 500000, prop.dominant = .02)
#'
#' @export


build_even_community <-
  function(n.taxa = 3000, # must be positive integer (or coerceable)
           n.samples = 100, # must be positive integer (or coerceable)
           n.reads = 5000000, # must be positive integer (or coerceable)
           prop.dominant = .01, # must be postive value between 0:1
           prop.absent = .2 # must be positive value between 0:1
  ){

    stopifnot(class(n.taxa) %in% c("numeric","integer") & n.taxa > 0)
    n.taxa <- as.integer(n.taxa)

    stopifnot(class(n.samples) %in% c("numeric","integer") & n.samples > 0)
    n.samples <- as.integer(n.samples)

    stopifnot(class(n.reads) %in% c("numeric","integer") & n.reads > 0)
    n.reads <- as.integer(n.reads)

    stopifnot(class(prop.dominant) == "numeric" & prop.dominant >= 0 & prop.dominant < 1)
    stopifnot(class(prop.absent) == "numeric" & prop.absent >= 0 & prop.absent < 1)

    # build empty matrix
    mat <- matrix(0,nrow = n.samples, ncol = n.taxa)

    # create taxa and sample names
    taxa_names <- paste0("taxon_",1:ncol(mat))
    colnames(mat) <- taxa_names
    sample_names <- paste0("sample_",1:nrow(mat))
    rownames(mat) <- sample_names

    # select dominant taxa (order from left-to-right in descending abundance)
    n.dom <- round(n.taxa * prop.dominant,0)

    # add relabund values for dominant taxa
    dom.relabunds <- sort(runif(n.dom,min = .01,max=.1),decreasing = TRUE)
    # find remainders
    remain.relabund <- 1-sum(dom.relabunds)
    remain.taxa <- (n.taxa - n.dom)

    if(remain.relabund <= .5){
      dom.relabunds <- dom.relabunds * .1
      remain.relabund <- 1-sum(dom.relabunds)
      remain.taxa <- (n.taxa - n.dom)
    }

    if(remain.relabund <= 0){
      n.dom <- 1
      # add relabund values for dominant taxa
      dom.relabunds <- sort(runif(n.dom,min = .01,max=.1),decreasing = TRUE)
      # find remainders
      remain.relabund <- 1-sum(dom.relabunds)
      remain.taxa <- (n.taxa - n.dom)
      warning("Too many dominant taxa for a matrix this size. Setting number of dominant taxa to 1.")
    }

    # split remaining relabund into n.taxa decreasing values
    split_into_decreasing_values <- function(num, n) {
      if (n <= 1) {
        stop("Number of splits must be greater than 1")
      }
      if (num <= 0) {
        stop("Number to split must be greater than 0")
      }
      # Calculate the splitting ratio
      ratio <- num / (1:n)
      # Make sure all splits are non-negative
      if (any(ratio <= 0)) {
        stop("Unable to split the number into non-negative values")
      }
      # Adjust the ratio to ensure sum equals num
      ratio <- ratio / sum(ratio) * num
      # Output the decreasing values
      decreasing_values <- sort(ratio, decreasing = TRUE)
      return(decreasing_values)
    }
    remaining_relabunds <- split_into_decreasing_values(num=remain.relabund, n = remain.taxa)

    # build vector of relative abundance values
    total_relabunds <- c(dom.relabunds,remaining_relabunds)
    read_counts <- round(n.reads * total_relabunds)
    read_counts[read_counts < 1] <- 1


    # calculate and add values (with zeros)
    split_randomly_with_zero_proportion <- function(num, n, zero_proportion) {
      zero_proportion <- abs(prop.absent + rnorm(1,mean = 0,sd=.1))
      if (n <= 0) {
        stop("Total number of splits must be greater than 0")
      }
      if (num <= 0) {
        stop("Number to split must be greater than 0")
      }
      if (zero_proportion < 0 || zero_proportion > 1) {
        stop("Zero proportion must be between 0 and 1")
      }
      # Calculate the number of zeroes
      num_zeroes <- round(n * zero_proportion)
      # Generate random values
      random_values <- c(rep(0, num_zeroes), runif(n - num_zeroes, min = 0, max = num))
      # Adjust values to sum up to num
      random_values <- random_values / sum(random_values) * num
      # Shuffle the values
      random_values <- sample(random_values)
      random_values <- round(random_values,0)
      return(random_values)
    }
    # fill matrix
    for(taxon in 1:ncol(mat)){
      mat[,taxon] <- split_randomly_with_zero_proportion(num = read_counts[taxon], n = n.samples,zero_proportion = prop.absent)
    }

    return(mat)

  }
