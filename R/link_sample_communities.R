
# Force similarity between groups of samples in terms of taxa counts
# n.links = number of groupings to force. Must be >= 1 and < sample.number/group.size
# group.size = how many samples to make similar in each grouping. Must be >= 2 and <= sample.number/n.links
# link.strength = how strong should the similarity function be. 0 = no change; 1=make identical.

link_sample_communities <- function(dat, n.links = 1, group.size = 3, link.strength = .9){
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

