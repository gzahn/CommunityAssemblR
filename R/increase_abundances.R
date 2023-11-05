
increase_abundances <- function(dat = NULL, # a matrix, samples are rows and columns are taxa
                                     prop = .1, # needs to be between 0,1
                                     increase.scale = .9,
                                     margin = "taxa" # needs to be "taxa" or "samples"
){
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


# example ####
# source("./R/comm_even.R")
# 
# # build a default taxa table with even distribution using comm_even
# x <- comm_even(n.taxa = 22,n.samples = 33,n.reads = 30455, taxa.sd = 3000) 

# rowSums(reduce_abundance(x,prop = .5,reduction.scale = .95,margin = "samples"))




