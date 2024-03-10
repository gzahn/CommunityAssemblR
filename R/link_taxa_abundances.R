#' Simulate co-occurrence relationships between taxa
#'
#' Manipulates a community matrix with samples as rows and taxa as columns.
#'
#' @param dat Numeric matrix. The starting community matrix that will be edited.
#' @param n.taxa Positive numeric vector of length 1. The number of taxa to link together. If n.taxa=1, nothing happens. If n.taxa=2, two taxa will be linked. If n.taxa > 2, a hub network will be created (you must select relationship='hub'). Default = 2.
#' @param relationship Character vector of length 1. Must be one of c("positive","negative","hub"). If "positive", then taxon 2 will be altered to increase proportionally with taxon 1. If "negative", the opposite will occurr. If "hub", approximately half of the taxa will have a positive and half a negative relationship with a random hub taxon. Default = "positive".
#' @param link.scale Numeric vector of length 1. how strong should the similarity function be. Multipliciative. 0=no change.
#' @param n.links Positive whole numeric vector of length 1. The number of taxa that will be linked with other taxa. Default = 20 percent of all samples.
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'. Random taxon abundance values will be altered to simulate co-occurrence according to parameters selected.
#'
#' @examples
#' comm <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#' link_taxa_abundances(dat = comm,n.taxa = 10,relationship = "hub",link.scale = .5,n.links = 2)
#'
#' @export

link_taxa_abundances <- function(dat, # community matrix with samples as columns and rows as taxa (assumes row and col names)
                                 n.taxa = 2, # if greater than 2, do the hub taxa thing
                                 relationship = "positive", # positive or negative or hub
                                 link.scale = 0.2, # must be positive numeric; multiplicative amount for transformations
                                 n.links = round(ncol(dat) * .2) / 2){ # default is 20% of taxa are lined with some other taxon
# setup and tests
stopifnot("matrix" %in% class(dat))
stopifnot("numeric" %in% class(n.taxa) & n.taxa > 0)
stopifnot(relationship %in% c("positive","negative","hub"))
stopifnot(link.scale >= 0 & class(link.scale) %in% c("numeric","integer"))

n.taxa = round(n.taxa)
if(n.taxa > 2 & relationship != "hub"){
  stop("If more than two taxa are to be linked, you must select relationship = 'hub'")
}

if(n.taxa < 3 & relationship == "hub"){
  stop("If relationship = 'hub' n.linked must be greater than 2")
}

# three different functions (all these update dat in-place)

# do this n.links times
for(i in 1:n.links){
# positive
if(relationship == "positive"){
  # two linked taxa (positive link)
  #       pick two random taxa
  linked.taxa <- sample(colnames(dat),n.taxa,replace = FALSE)
  # where one is zero, make both zero
  taxa1_zeros <- unname(which(dat[,linked.taxa[1]] == 0))
  taxa2_zeros <- unname(which(dat[,linked.taxa[2]] == 0))
  dat[,linked.taxa[1]][taxa2_zeros] <- 0
  dat[,linked.taxa[2]][taxa1_zeros] <- 0
  # raise their abundances by an equal factor
  dat[,linked.taxa] <- dat[,linked.taxa] * round((1 + link.scale))
  }
}

for(i in 1:n.links){
# negative
if(relationship == "negative"){
  # two linked taxa (negative link)
  #       pick two random taxa
  linked.taxa <- sample(colnames(dat),n.taxa,replace = FALSE)
  # where one is zero, make the other one more abundant
  # find the more-abundant taxon
  greater.abund.taxa <- ifelse(sum(dat[,linked.taxa[1]]) > sum(dat[,linked.taxa[2]]),1,2)
  lesser.abund.taxa <- ifelse(sum(dat[,linked.taxa[1]]) < sum(dat[,linked.taxa[2]]),1,2)
  greater.abund.zeros <- unname(which(dat[,linked.taxa[greater.abund.taxa]] == 0))
  lesser.abund.zeros <- unname(which(dat[,linked.taxa[lesser.abund.taxa]] == 0))
  # where the less abundant taxon is zero, make the more abundant taxon even more abundant
  dat[lesser.abund.zeros,linked.taxa][,greater.abund.taxa] <-
    round(dat[lesser.abund.zeros,linked.taxa][,greater.abund.taxa] * (1 + link.scale))
  # where both taxa have observations, reduce the abundance of the less abundant taxon even more
  greater.abund.observations <- unname(which(dat[,linked.taxa[greater.abund.taxa]] != 0))
  lesser.abund.observations <- unname(which(dat[,linked.taxa[lesser.abund.taxa]] == 0))
  samples.both.positive <-
    which(dat[,linked.taxa][,greater.abund.taxa] > 0 & dat[,linked.taxa][,lesser.abund.taxa] > 0)
  # raise abundance of more abuntant taxon
  dat[,linked.taxa][,greater.abund.taxa][samples.both.positive] <-
  round(dat[,linked.taxa][,greater.abund.taxa][samples.both.positive] * (1 + link.scale))
  # reduce abundance of less abundant taxon
  dat[,linked.taxa][,lesser.abund.taxa][samples.both.positive] <-
    round(dat[,linked.taxa][,lesser.abund.taxa][samples.both.positive] * (1 - link.scale))
  }
}


# hub
for(i in 1:n.links){
if(relationship == "hub"){
  #       how many links in the hub?
  #       pick N taxa (n.taxa = 11 for testing)
  linked.taxa <- sample(colnames(dat),n.taxa,replace = FALSE)
  # divide roughly in half and select hub candidate
  first.half <- linked.taxa[1:round(length(linked.taxa) / 2)]
  second.half <- linked.taxa[(1 + round(length(linked.taxa)) / 2) : length(linked.taxa)]
  hub.taxon <- first.half[1]
  first.half <- first.half[-1]
  # where is hub taxon present?
  hub.present <- unname(which(dat[,hub.taxon] > 0))
  # raise first half where hub is present
  dat[hub.present,first.half] <- round(dat[hub.present,first.half] * (1 + link.scale))
  # lower second half when hub is present
  dat[hub.present,second.half] <- round(dat[hub.present,second.half] * (1 - link.scale))
  # raise hub taxon
  dat[,hub.taxon] <- round(dat[,hub.taxon] * (1 + link.scale))
  }
}

# in case of any negative values, replace with 0
dat[dat < 0] <- 0

# end of function
return(dat)
}


