library(CommunityAssemblR)

# build a basic even resident community
comm <- CommunityAssemblR::build_even_community(n.taxa = 350,n.samples = 100,n.reads = 10000,taxa.sd = 50)
# build donor communities on same relative scale as the resident community (same number of samples)
donor <- build_donor_community(resident.comm = comm,n.transplant.taxa = 50,overlap = 0)

# do a simple transplantation with known ecological process
# in this case, there are 8 resident taxa that occupy a niche that enables survival of 3 donor taxa
# when the resident taxa surpass 0.01 relative abundance in a given sample, the selected donor taxa
# are allowed to persist. If that threshold is not met, the donor taxa abundance drops to 0 in that sample
trans <- transplant_w_niche_occupation(recipient = comm,
                                       donor = donor,
                                       n.niches = 5,
                                       niche.size.resident = 8,
                                       niche.size.donor = 3,
                                       abundance.threshold = 0.01,
                                       process.type = "permission",
                                       print.niche.taxa = TRUE)

trans

# pull out the donor taxa from the final community
after <- trans[,grep("newtaxon",colnames(trans))]

# convert both to relative abundance so they can be compared at same scale
before <- t(apply(donor,1,function(x){x/sum(x)}))
after <- t(apply(after,1,function(x){x/sum(x)}))

# 20% of the donor taxa had a random recipient taxon that increased their abundance proportionally
heatmap(before,Rowv = NA,Colv = NA)
heatmap(after,Rowv = NA,Colv = NA)

# Can we build a model that can detect which recipient
# taxa were occupying niches that 'permitted' donor taxa to persist?

# these are the recipient taxa that facilitated some of the donor taxa
# it's possible that some of them overlap between 'rounds' of the transplant_ function
# because some taxa may occupy niche space that affects multiple other taxa in different ways
tmp_selected_niche.taxa
tmp_selected_niche.taxa$resident %>% as.data.frame()
tmp_selected_niche.taxa$donor %>% as.data.frame()

# looking closely at how round_1 of niche selection played out...

# pull out taxa names from round_1
r.taxa <- tmp_selected_niche.taxa$resident["round_1"] %>% unlist %>% unname
d.taxa <- tmp_selected_niche.taxa$donor["round_1"] %>% unlist %>% unname

# subset main resident/recipient matrix (as relative abundance matrix)
comm_relabund <- t(apply(comm,1,function(x){x/sum(x)}))
r.matrix <- comm_relabund[,r.taxa]

# in which samples did at least one of the niche-occupying taxa cross the abundance threshold?
t(apply(r.matrix,1,function(x){ifelse(x >= 0.01,TRUE,FALSE)}))
permissive.samples <- which(rowSums(t(apply(r.matrix,1,function(x){ifelse(x >= 0.01,TRUE,FALSE)}))) > 0)
# in only these samples should the 3 donor taxa have any abundance at all
# (with some gaussian variation, and disregarding the other 2 rounds of selection that were performed, of course)

heatmap(before[,d.taxa],Rowv = NA,Colv = NA)
heatmap(after[,d.taxa],Rowv = NA,Colv = NA)

