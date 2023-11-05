# simulate transplantation with recipient antagonism

# give the some recipient community members an antagonistic relationship with some donor community members

# takes 2 communities (recipient/resident and donor)
# donor community should be pretty even, and same number of samples as recipient

# antag.ubiq (what proportion of your donor taxa have antagonists in the resident community?)
# antag.strength (what scale of antagonism is simulated?)

scale01 <- function(x){
  if(max(x) == 0){return(x)}
  if(max(x) > 0){return((x - min(x)) / (max(x) - min(x)))}
}

transplant_w_antagonism <- function(recipient,donor,antag.ubiq=.5,antag.strength=.1,antag.abundant=TRUE){
  
  stopifnot("matrix" %in% class(donor))
  stopifnot("matrix" %in% class(recipient))
  stopifnot(class(antag.ubiq) %in% c("numeric","integer") & antag.ubiq >= 0 & antag.ubiq <= 1)
  stopifnot(class(antag.strength) %in% c("numeric","integer") & antag.strength >= 0)
  
  
  
  # rescale donor community to same range as recipient community
  # scaled against mean from recipient
  donor <- round(apply(donor,2,scale01) * (round(median(colSums(recipient) / nrow(recipient)))))
  
  
  # where taxa overlap, just take sum of both samples
  shared_donor <- donor[,colnames(donor) %in% colnames(recipient)]
  shared_recipient <- recipient[,colnames(recipient) %in% colnames(shared_donor)]
  new_shared_taxa <- shared_donor + shared_recipient
  
  recipient[,colnames(new_shared_taxa)] <- new_shared_taxa
  
  # where new taxa arriving, reduce the new taxa, scaled by antag.strength
  novel_taxa <- colnames(donor)[!colnames(donor) %in% colnames(recipient)]
  # rearrange in order because I'm OCD
  novel_taxa <- novel_taxa[order(as.numeric(sub("newtaxon_","",x = novel_taxa)))]
  
  # select antagonists from resident community
  num.antagonists <- round(length(novel_taxa)*antag.ubiq)
  if(!antag.abundant){ # choose random facilitator taxa
    antag.taxa <- sample(colnames(recipient),num.antagonists)  
  }
  if(antag.abundant){ # choose most abundant facilitator taxa
    antag.taxa <- names(colSums(recipient)[order(colSums(recipient),decreasing = TRUE)][1:num.antagonists])
  }
  antagonists <- recipient[,antag.taxa]
  
  # assign them to a paired novel taxon from the donor community
  antagonized <- donor[,sample(novel_taxa,num.antagonists)]
  
  # reduce 'antagonized' taxa by multiplier based on antagonists' scaled abundances
  antag.multiplier <- 1 - (t(apply(recipient,1,function(x){x/sum(x)}))[,antag.taxa] * antag.strength)
  antagonized <- round(antagonized * antag.multiplier)
  antagonized[antagonized < 0] <- 0 # take care of any potential negative values
  newdonor <- donor[,novel_taxa]
  newdonor[,colnames(antagonized)] <- antagonized
  
  
  # combine reduced new taxa with augmented recipient community
  
  final_community <- cbind(newdonor,recipient)
  
  return(final_community)
  
}

# # examples
# 
# 
# # make a resident community matrix (even)
# even <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30) 
# # add some hub taxa relationships
# recipient <- link_taxa_abundances(even,n.taxa = 20, relationship = 'hub',link.scale = 3)
# # build a donor community based on even resident community
# donor <- build_donor_community(resident.comm = even,n.transplant.taxa = 10,overlap = .3)
# # perform transplantation with resident antagonism
# transplant_w_antagonism(recipient, donor, antag.ubiq = .5, antag.strength = .5)

