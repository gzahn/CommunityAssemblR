# takes 2 communities (recipient/resident and donor)
# donor community should be pretty even, and same number of samples as recipient

# variables

transplant_w_niche_occupation <- function(n.niches = 1, # number of rounds of niche selection
                                          niche.shape = 'normal',
                                          niche.size.resident = 3, # number of taxa in resident community that occupy the niche
                                          niche.size.donor = 1, # number of taxa in donor community (new taxa only) that occupy the niche
                                          abundance.threshold = .05, # relative abundance threshold for resident taxa in niche to trigger effect
                                          # this will be mean of distribution...if niche.size.resident > 1, thresholds for additional taxa
                                          # will be based on the distribution of niche.shape
                                          process.type = "preemption"){ # 'preemption' or 'permission'
  
  stopifnot("matrix" %in% class(donor))
  stopifnot("matrix" %in% class(recipient))
  
  if(!niche.shape %in% c("normal")){
    stop("Currently, only 'normal' distribution accepted for niche.shape")
  }
  if(abundance.threshold <= 0 | abundance.threshold > 1){
    stop("abundance.threshold must be between 0 and 1")
  }
  
  stopifnot(class(niche.size.donor) == "numeric")
  niche.size.donor <- round(niche.size.donor)
  stopifnot(niche.size.donor > 0)
  
  stopifnot(class(niche.size.resident) == "numeric")
  niche.size.resident <- round(niche.size.resident)
  stopifnot(niche.size.resident > 0)
  
  stopifnot(class(n.niches) == "numeric")
  n.niches <- round(n.niches)
  stopifnot(n.niches > 0)
  
  if(!process.type %in% c("preemption","permission")){
    stop("process.type must be either 'preemtion' or 'permission'.")
  }
  
  # rescale donor community to same range as recipient community
  # scaled against mean from recipient
  donor <- round(apply(donor,2,scale01) * (round(median(colSums(recipient)))))
  
  # where taxa overlap, just take sum of both samples
  shared_donor <- donor[,colnames(donor) %in% colnames(recipient)]
  shared_recipient <- recipient[,colnames(recipient) %in% colnames(shared_donor)]
  new_shared_taxa <- shared_donor + shared_recipient
  
  recipient[,colnames(new_shared_taxa)] <- new_shared_taxa
  
  # find novel taxa from donor community
  novel_taxa <- colnames(donor)[!colnames(donor) %in% colnames(recipient)]
  # rearrange in order because I'm OCD
  novel_taxa <- novel_taxa[order(as.numeric(sub("newtaxon_","",x = novel_taxa)))]
  novel.taxa.matrix <- donor[,novel_taxa]
  
  
  ### Niche occupation process, repeat n.niches times
  
  for(X in 1:n.niches){
    # Select resident abundance thresholds for N taxa, based on niche.shape distribution
    abund.thresholds <- abs(rnorm(n=niche.size.resident,mean = abundance.threshold,sd = .05))
    
    # select recipient taxa to act as preemtion or permission gates
    resident.niche.taxa <- sample(colnames(recipient),niche.size.resident)
    
    # select donor taxa that will be influenced by resident niche occupation
    donor.niche.taxa <- sample(novel_taxa,niche.size.donor)
    
    # convert to relative abundance
    recipient_ra <- t(apply(recipient,1,function(x){x/sum(x)}))
    
    # subset to just the resident niche taxa
    ra <- recipient_ra[,resident.niche.taxa]
    
    # compare each column with respective abundance threshold, give T/F
    tests <- matrix(nrow = nrow(ra),ncol = niche.size.resident)
    for(i in 1:niche.size.resident){
      tests[,i] <- ra[,i] > abund.thresholds[i]
    }
    
    # if any one of the rows meets threshold, then effect can take place (either prevention or allowance)
    met.threshold <- as.logical(rowSums(tests))
    
    # pre-emption: if threshold is met in resident community, it prevents donor taxa from occupying
    if(process.type == "preemption"){
      novel.taxa.matrix[,donor.niche.taxa] <- (novel.taxa.matrix[,donor.niche.taxa] * !met.threshold)
    }
    
    # permission: if threshold is met in resident community, it allows donor taxa to occupy
    if(process.type == "permission"){
      novel.taxa.matrix[,donor.niche.taxa] <- (novel.taxa.matrix[,donor.niche.taxa] * met.threshold)
    }
    
  }
  
  final_community <- cbind(novel.taxa.matrix,recipient)  
  
  return(final_community)
  
}

