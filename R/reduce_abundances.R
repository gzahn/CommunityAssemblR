
reduce_abundances <- function(dat = NULL, # a matrix, samples are rows and columns are taxa
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


# example ####
# source("./R/comm_even.R")
# 
# # build a default taxa table with even distribution using comm_even
# x <- comm_even(n.taxa = 22,n.samples = 33,n.reads = 30455, taxa.sd = 3000) 

# rowSums(reduce_abundance(x,prop = .5,reduction.scale = .95,margin = "samples"))




