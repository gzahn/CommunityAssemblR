library(CommunityAssemblR)


# build a basic even resident community
comm <- CommunityAssemblR::build_even_community(n.taxa = 350,n.samples = 100,n.reads = 10000,taxa.sd = 50)
# build donor communities on same relative scale as the resident community (same number of samples)
donor <- build_donor_community(resident.comm = comm,n.transplant.taxa = 50,overlap = 0)

# do a simple transplantation with known ecological process
# in this case, there are 8 resident taxa that occupy a niche that enables survival of 3 donor taxa
# when the resident taxa surpass 0.1 relative abundance in a given sample, the selected donor taxa
# are allowed to persist. If that threshold is not met, the donor taxa abundance drops to 0 in that sample
trans <- transplant_w_niche_occupation(recipient = comm,
                                       donor = donor,
                                       n.niches = 5,
                                       niche.size.resident = 8,
                                       niche.size.donor = 3,
                                       abundance.threshold = 0.1,
                                       process.type = "permission",
                                       print.niche.taxa = TRUE)

# retrieve the donor taxa so we can run another ecological process on them
updated <- retrieve_donor_taxa(trans)

# run a different process (another transplantation round)
trans <- transplant_w_niche_occupation(recipient = updated$new.recipients,
                              donor = updated$new.donors,
                              n.niches = 1,
                              niche.size.resident = 8,
                              niche.size.donor = 3,
                              abundance.threshold = 0.1,
                              process.type = 'preemption',
                              print.niche.taxa = TRUE)

# so now, we have a transplated community that was subject to 'niche permission' AND 'niche preemption'

# add some stochasiticity
updated <- retrieve_donor_taxa(trans)
trans <- transplant_w_stochasticity(updated$new.recipients,updated$new.donors)

# add an environmental gradient that negatively affects 25% of the donor taxa
updated <- retrieve_donor_taxa(trans)
newdonor <- filter_taxa_along_gradient(updated$new.donors,prop=.25,gradient.strength = 8,groups = 3)

final <- transplant_w_stochasticity(updated$new.recipients,newdonor,stochasticity = 0)


# now we can pull out our taxa of interest that have been subjected to those 3 processes
new <- final[,grepl("newtaxon",colnames(final))]

heatmap(donor, Rowv = NA)
heatmap(new, Rowv = NA)

# we lost some donor taxa completely
new[,colSums(new) == 0]
