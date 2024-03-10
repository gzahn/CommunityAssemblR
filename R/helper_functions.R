
# plot igraph, showing hub taxa articulation points as larger dots, sized by authority score
plot_hubs <- function(graph,bigpoint=10,littlepoint=3){
  am.coord <- layout.auto(graph)
  art <- articulation_points(graph)
  plot(graph,
       layout = am.coord,
       vertex.size=(scale01(abs(igraph::authority_score(graph)$vector)) * 10)+3,
       vertex.label=NA)
}


# # custom scaling function (0,1) that can handle vectors of all zeros
# scale01 <- function(x){
#   if(max(x) == 0){return(x)} # if all zeros, do nothing
#   if(max(x) > 0){return((x - min(x)) / (max(x) - min(x)))} # if some positive, scale as normal
#   if(sum(x) < 0){x <- abs(x) # not sure about this
#                  scaled <- (x - min(x)) / (max(x) - min(x))
#                  return(-scaled)} # if all negative, use absolute values, then replace negative sign
# }








######## WHEN you do the relative abundance transformation matters #########
# Since RA transformation is taking place in the final community with ALL the other taxa present,
# that can throw off the ability to compare before and after
# maybe make an if() that changes where to do the transformation?


# # function to gather up donor microbiome community relative abundances before and after transplanting
# check_transplant_success <- function(donor,transplant,newtaxa.only=FALSE,tidy=TRUE,keep.all.taxa=TRUE){
#
#
#   if(keep.all.taxa == TRUE){
#     donor_ra <- t(apply(donor,1,function(x){x/sum(x)}))
#     recipient_ra <- t(apply(transplant,1,function(x){x/sum(x)}))
#   }
#
#   if(keep.all.taxa == FALSE){
#     donor_ra <- t(apply(donor,1,function(x){x/sum(x)}))
#     recipient_ra <- t(apply(transplant[,colnames(donor_ra)],1,function(x){x/sum(x)}))
#   }
#
#
#
#   # if newtaxa.only=TRUE, subset to just the "newtaxa"
#   if(newtaxa.only == TRUE){
#     initial <- donor_ra[,grep(colnames(donor_ra),pattern="newtaxon")]
#     final <- recipient_ra[,colnames(initial)]
#   }
#
#   # in newtaxa.only=FALSE, then just subset to the taxa found in the donor community
#   if(newtaxa.only == FALSE){
#     initial <- donor_ra
#     final <- recipient_ra[,colnames(initial)]
#   }
#
#   # convert to data.frames
#   initial <- as.data.frame(initial)
#   final <- as.data.frame(final)
#
#   # add columns indicating timepoint
#   initial$timepoint <- "Initial"
#   final$timepoint <- "Final"
#
#   # join data frames
#   initial <- initial[,names(initial)[order(names(initial))]]
#   final <- final[,names(initial)[order(names(initial))]]
#   df <- rbind(initial,final)
#
#   df$sample <- rep(row.names(initial),2)
#
#   if(tidy){
#     df <- reshape(df,grep("tax",names(df),value = TRUE),direction = 'long',
#                   v.names = "realative_abundance",times = grep("tax",names(df),value = TRUE))
#     df$id <- NULL
#     row.names(df) <- NULL
#     names(df) <- c("timepoint","sample","taxon","relative_abundance")
#   }
#
#   return(df)
# }

#

# library(tidyverse)
# comm <- build_even_community(n.taxa = 30, n.samples = 10,taxa.sd = 300)
# donor <- build_donor_community(comm,n.transplant.taxa = 10, overlap = .5)
# transplant <- transplant_w_facilitation(comm,donor,facil.strength = 10,facil.ubiq = 1)
# transplant <- transplant_w_niche_occupation(comm,donor,n.niches = 1,niche.size.donor = 2,process.type = "permission") %>%
#   simulate_competition(niche.separation = 4)
# x <- check_transplant_success(donor,transplant,tidy = TRUE,newtaxa.only = FALSE,keep.all.taxa = FALSE)
#
# x %>%
#   group_by(timepoint,sample,taxon) %>%
#   summarize(total=sum(relative_abundance)) %>%
#   ggplot(aes(x=sample,y=taxon,fill=total)) +
#   geom_tile() +
#   facet_wrap(~factor(timepoint,levels = c("Initial","Final"))) +
#   theme(axis.text.x = element_blank())
# #
#
# #
#
#
