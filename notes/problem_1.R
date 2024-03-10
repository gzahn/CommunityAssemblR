library(CommunityAssemblR)

# build a basic even resident community
comm <- CommunityAssemblR::build_even_community(n.taxa = 350,n.samples = 100,n.reads = 10000,taxa.sd = 50)
# build donor communities on same relative scale as the resident community (same number of samples)
donor <- build_donor_community(resident.comm = comm,n.transplant.taxa = 50,overlap = 0)

# do a simple transplantation with known ecological process
# in this case, there are resident taxa that facilitate 20% of the donor taxa
trans <- transplant_w_facilitation(recipient = comm,
                          donor = donor,
                          facil.ubiq = .2,
                          facil.strength = 1,
                          facil.abundant = TRUE,
                          print.facil.taxa = TRUE)
trans

# pull out the donor taxa from the final community
after <- trans[,grep("newtaxon",colnames(trans))]

# convert both to relative abundance so they can be compared at same scale
before <- t(apply(donor,1,function(x){x/sum(x)}))
after <- t(apply(after,1,function(x){x/sum(x)}))

# 20% of the donor taxa had a random recipient taxon that increased their abundance proportionally
heatmap(before,Rowv = NA,Colv = NA)
heatmap(after,Rowv = NA,Colv = NA)

# Can we build a model that can detect which recipient taxa were facilitators?
tmp_facil.taxa # these are the recipient taxa that facilitated some of the donor taxa

# The main model idea would be something like:
after ~ comm # final composition of donor taxa as a function of starting composition of resident taxa

# I think it makes sense to transform/scale values in both matrices to reflect compositionality before modeling

# This is a simple test case with only one ecological process going on (facilitation) and no interaction
# between taxa in that process. A starting point to get into the headspace.
