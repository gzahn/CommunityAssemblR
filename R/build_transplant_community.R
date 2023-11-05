# build a donor community based on a recipient community, matching sample numbers, etc.
# provide resident/recipient community object, and build a transplant/donor community matrix with N taxa


# samples are rows, taxa are cols

build_donor_community <- function(resident.comm,n.transplant.taxa,overlap=0){
  
  # set up parameters for build_even_community
  n.taxa <- round(n.transplant.taxa)
  stopifnot(class(n.transplant.taxa) %in% c("numeric","integer") & n.transplant.taxa > 0)
  stopifnot("matrix" %in% class(resident.comm))
  n.samples <- nrow(resident.comm)
  stopifnot(class(overlap) %in% c("numeric","integer") & overlap >= 0 & overlap <= 1)
  n.reads <- round(mean(rowSums(resident.comm)))
  taxa.sd <- round(sd(colSums(resident.comm)))
  taxa.dist = "normal"
  message("Assuming that samples are rows / taxa are columns...")
  message(paste0("Building new 'donor' community with ",overlap*100,"% overlap in community membership."))
  
  # build basic community of 'transplant' microbes
  transplant.comm <- 
    build_even_community(n.taxa,n.samples,n.reads,taxa.dist,taxa.sd)
  
  
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

# example
# build_donor_community(resident.comm = x,n.transplant.taxa = 20,overlap = .2)
