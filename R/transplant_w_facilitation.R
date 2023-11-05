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

transplant_w_facilitation <- function(recipient,donor,facil.ubiq=.5,facil.strength=1,facil.abundant=TRUE){
  
  stopifnot("matrix" %in% class(donor))
  stopifnot("matrix" %in% class(recipient))
  stopifnot(class(facil.ubiq) %in% c("numeric","integer") & facil.ubiq >= 0 & facil.ubiq <= 1)
  stopifnot(class(facil.strength) %in% c("numeric","integer") & facil.strength >= 0)
  
  # rescale donor community to same range as recipient community
  # scaled against mean from recipient
  donor <- round(apply(donor,2,scale01) * (round(median(colSums(recipient)))))
  
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
  num.facilitators <- round(length(novel_taxa)*facil.ubiq)
  if(!facil.abundant){ # choose random facilitator taxa
    facil.taxa <- sample(colnames(recipient),num.facilitators)  
  }
  if(facil.abundant){ # choose most abundant facilitator taxa
    facil.taxa <- names(colSums(recipient)[order(colSums(recipient),decreasing = TRUE)][1:num.facilitators])
  }
  facilitators <- recipient[,facil.taxa]
  
  # assign them to a paired novel taxon from the donor community
  facilitated <- donor[,sample(novel_taxa,num.facilitators)]
  
  # increase 'facilitated' taxa by multiplier based on antagonists' scaled abundances
  facil.multiplier <- 1 + (t(apply(recipient,1,function(x){x/sum(x)}))[,facil.taxa] * facil.strength)
  facilitated <- round(facilitated * facil.multiplier)
  
  newdonor <- donor[,novel_taxa]
  newdonor[,colnames(facilitated)] <- facilitated
  
  
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

