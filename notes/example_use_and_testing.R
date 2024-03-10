library(tidyverse)
library(igraph); packageVersion("igraph")
library(SpiecEasi); packageVersion("SpiecEasi")

# load CommunityAssemblR functions
source("./R/build_even_community.R")
source("./R/link_taxa_abundances.R")
source("./R/reduce_abundances.R")
source("./R/increase_abundances.R")
source("./R/build_transplant_community.R")
source("./R/transplant_w_antagonism.R")
source("./R/transplant_w_facilitation.R")
source("./R/transplant_w_niche_occupation.R")
source("./R/helper_functions.R")

# build an even community table
x <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
y <- link_taxa_abundances(dat = x,n.taxa = 50,relationship = "hub",link.scale = 3,n.links = 5)
z <- link_taxa_abundances(dat = x,n.taxa = 2,relationship = "positive",link.scale = 3,n.links = 10) %>% 
  link_taxa_abundances(n.taxa = 50,relationship = "hub",link.scale = 3,n.links = 5)


# make SpiecEasi networks
se.params <- list(rep.num=20, ncores=(parallel::detectCores()-1))

# altered with hub taxa
se.mb.x <- SpiecEasi::spiec.easi(data = x,
                                 method='mb',
                                 sel.criterion = "bstars",
                                 pulsar.params=se.params)

se.mb.y <- SpiecEasi::spiec.easi(data = y,
                                    method='mb',
                                    sel.criterion = "bstars",
                                    pulsar.params=se.params)

se.mb.z <- SpiecEasi::spiec.easi(data = z,
                                 method='mb',
                                 sel.criterion = "bstars",
                                 pulsar.params=se.params)

# refit to igraph objects
se.igraph.x <- adj2igraph(getRefit(se.mb.x))
se.igraph.y <- adj2igraph(getRefit(se.mb.y))
se.igraph.z <- adj2igraph(getRefit(se.mb.z))

# plot graphs
am.coord.x <- layout.auto(se.igraph.x)
plot(se.igraph.x, layout=am.coord.x, vertex.label=NA)

am.coord.y <- layout.auto(se.igraph.y)
plot(se.igraph.y, layout=am.coord.y, vertex.label=NA)

am.coord.z <- layout.auto(se.igraph.z)
plot(se.igraph.z, layout=am.coord.z)

# inspect graph properties
igraph::degree.distribution(se.igraph.x) %>% plot
igraph::degree.distribution(se.igraph.y) %>% plot
hubx <- igraph::hub_score(se.igraph.x)
huby <- igraph::hub_score(se.igraph.y)
hubz <- igraph::hub_score(se.igraph.z)
hubx$value
huby$value
hubz$value
# find articulation points
artpoints.x <- igraph::articulation.points(se.igraph.x)
artpoints.y <- igraph::articulation.points(se.igraph.y)
artpoints.z <- igraph::articulation.points(se.igraph.z)
y[,artpoints.y]

igraph::assortativity_degree(se.igraph.z)

# re-plot, but with 'hub taxa' identified
# plot_igraph_w_hubs <- function(graph,artpoints,big.size=10,small.size=3){
#   plot(graph,)
# }
plot(se.igraph.x, 
     layout=am.coord.x,
     vertex.size=ifelse(1:ncol(z) %in% artpoints.x,10,3), 
     vertex.label=NA)

plot(se.igraph.y, 
     layout=am.coord.y,
     vertex.size=ifelse(1:ncol(z) %in% artpoints.y,10,3), 
     vertex.label=NA)

plot(se.igraph.z, 
     layout=am.coord.z,
     vertex.size=ifelse(1:ncol(z) %in% artpoints.z,10,3), 
     vertex.label=NA)


# build N artificial communities with some set function parameters
# inspect outcomes
# store outcomes to compare with how parameters affect them

# network property values to calculate
hubscores <- list()
assortativity <- list()
betweenness <- list()
avg_path_length <- list()


for(links in 1:5){
  
hub_list <- list()
for(i in 1:10){
  hub_list[[i]] <- link_taxa_abundances(dat = even,n.taxa = 10,relationship = "hub",link.scale = 2*as.numeric(i),n.links = 2*as.numeric(i)) 
}

hub_se_list <- list()
hub_se_graph <- list()
for(i in seq_along(hub_list)){
  hub_se_list[[i]] <- SpiecEasi::spiec.easi(data = hub_list[[i]],
                        method='mb',
                        sel.criterion = "bstars",
                        pulsar.params=se.params)
  hub_se_graph[[i]] <- adj2igraph(getRefit(hub_se_list[[i]]))
}

for(i in 1:length(hub_se_graph)){
  print(i)
  print(plot_hubs(hub_se_graph[[i]]))
}
mygraph=hub_se_graph[[1]]
plot_hubs(mygraph)


# export network property results
hubscores[[links]] <- 
  map(hub_se_graph,igraph::hub_score) %>% map_dbl("value")
assortativity[[links]] <- 
  map(hub_se_graph,igraph::assortativity_degree)
betweenness[[links]] <- 
  map(hub_se_graph,igraph::betweenness)
avg_path_length[[links]] <- 
  map(hub_se_graph,igraph::average.path.length)


}

map(hub_se_graph,plot_hubs)


# plot hub scores
df <- hubscores %>% reduce(cbind) %>% as.data.frame
names(df) <- paste0("link_",1:5)
df %>% 
  pivot_longer(everything(),names_to = "num_links",
               values_to = "network_hub_score",names_prefix = "link_") %>% 
  ggplot(aes(x=network_hub_score,fill=num_links)) +
  ggridges::geom_density_ridges(aes(y=as.numeric(num_links))) +
  theme_minimal() +
  scale_fill_viridis_d(end=.8) +
  labs(y="Num. links in\nlink_taxa_abundances()",
       x="Network hub score",
       fill="Num. links")


# plot assortativity
data.frame(
  n_links = rep(1:5,each=10),
  assortativity = assortativity %>% unlist
) %>% 
  ggplot(aes(x=assortativity)) +
  ggridges::geom_density_ridges(aes(y=factor(n_links))) +
  theme_minimal() +
  labs(y="Num. links in\nlink_taxa_abundances()",
       x="Assortativity",
       fill="Num. links")

# plot betweenness
s <- summary(unlist(betweenness))

data.frame(
  betweenness = unlist(betweenness),
  n_links = rep(1:5,each=1000)
) %>% 
  filter(betweenness > 0) %>% 
  ggplot(aes(x=betweenness)) +
  ggridges::geom_density_ridges(aes(y=factor(n_links))) +
  theme_minimal() +
  labs(y="Num. links in\nlink_taxa_abundances()",
       x="Betweenness",
       fill="Num. links")

# plot average path length
data.frame(
  n_links = rep(1:5,each=10),
  avg_path_length = avg_path_length %>% unlist
) %>% 
  ggplot(aes(x=avg_path_length)) +
  ggridges::geom_density_ridges(aes(y=factor(n_links))) +
  theme_minimal() +
  labs(y="Num. links in\nlink_taxa_abundances()",
       x="Avg. path length",
       fill="Num. links")



#










# # make a resident community matrix (even)
even <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30) 
# # add some hub taxa relationships
recipient <- link_taxa_abundances(even,n.taxa = 35, relationship = 'hub',link.scale = 3,n.links = 3) %>% 
  increase_abundances(prop = .2,increase.scale = 3,margin = "taxa")
# # build a donor community based on even resident community
donor <- build_donor_community(resident.comm = even, n.transplant.taxa = 30,overlap = .75)
# # perform transplantation with resident antagonism
final <- transplant_w_antagonism(recipient, donor, antag.ubiq = .5, antag.strength = 10, antag.abundant = TRUE)
final2 <- transplant_w_facilitation(recipient, donor, facil.ubiq = .5, facil.strength = 10,facil.abundant = TRUE)
final3 <- transplant_w_niche_occupation(niche.size.donor = 2, niche.size.resident = 3, process.type = "permission")

x <- final3[,grep(pattern = "newtaxon_",colnames(final3))]
y <- donor[,colnames(donor) %in% colnames(x)] 

heatmap(x,Rowv = NA,Colv = NA);heatmap(y,Rowv = NA,Colv = NA)
rowSums(x)
rowSums(y)

a <- 
x %>% 
  as.data.frame() %>% 
  mutate(community="After",sample=row.names(x)) %>%
  pivot_longer(starts_with("newtaxon")) 
b <- 
y %>% 
  as.data.frame() %>% 
  mutate(community = "Before",sample=row.names(y)) %>% 
  pivot_longer(starts_with("newtaxon"))
full_join(a,b) %>% 
    ggplot(aes(x=sample,y=name,fill=value)) +
    geom_tile() +
  facet_wrap(~community)



colnames(donor)
# compare original donor newtaxa to final donor newtaxa
### Make this a function ###


# transform to relative abundance
final_ra <- t(apply(final,1,function(x){x/sum(x)})) 
final2_ra <- t(apply(final2,1,function(x){x/sum(x)}))
donor_ra <- t(apply(donor,1,function(x){x/sum(x)}))


rowSums(donor_ra)
rowSums(final2_ra)
final_transplants <- final_ra[,grepl("newtaxon_",colnames(final_ra))]
final2_transplants <- final2_ra[,grepl("newtaxon_",colnames(final2_ra))]
rowSums(final2_transplants)
starting_transplants <- donor_ra[,grepl("newtaxon_",colnames(donor_ra))]
rowSums(starting_transplants)
matrix_diff <- (final_transplants - starting_transplants)
matrix2_diff <- (final2_transplants - starting_transplants)
final_transplants - final2_transplants

# apply(matrix_diff,2,scale01)
as.data.frame(matrix_diff) %>% 
  mutate(sample=row.names(.)) %>% 
  pivot_longer(starts_with("newtaxon_")) %>% 
  ggplot(aes(x=sample,y=name,fill=value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_blank())

heatmap(
  matrix_diff,Colv = NA  
)

as.data.frame(matrix_diff) %>% 
  mutate(sample=row.names(.)) %>% 
  pivot_longer(starts_with("newtaxon_")) %>% 
  ggplot(aes(x=value)) +
  ggridges::geom_density_ridges(aes(y=name)) +
  theme_minimal()

matrix_diff %>% sum
matrix2_diff %>% sum

as.data.frame(matrix2_diff) %>% 
  mutate(sample=row.names(.)) %>% 
  pivot_longer(starts_with("newtaxon_")) %>% 
  ggplot(aes(x=value)) +
  ggridges::geom_density_ridges(aes(y=name)) +
  theme_minimal()
### NOTE: SINCE WE'RE LOOKING AT RELABUND HERE, IF THE RECIPIENT COMMUNITY HAS A LOT MORE TAXA THAN THE DONOR,
###       THEN THE RELABUND OF THE NEW TAXA WILL BE LOWER BY DEFAULT IN THE FINAL COMMUNITY
###       HOW CAN WE DEAL WITH THIS????!
###       In other words, neither antagonism or facilitation have much effect on relative abundances,
###       compared to the overwhelming effect of recipient community richness
###
###       Can we scale against number of taxa, perhaps, before calculating relative abundance?
