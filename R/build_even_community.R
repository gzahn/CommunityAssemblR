# MAKE A TAXA ABUNDANCE TABLE
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
    stopifnot(class(taxa.sd) %in% c("numeric","logical"))

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

# build_even_community(n.taxa = 22,n.samples = 33,n.reads = 30455, taxa.sd = 3000) 
