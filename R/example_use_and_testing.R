library(tidyverse)
library(igraph); packageVersion("igraph")
library(SpiecEasi); packageVersion("SpiecEasi")

# load CommunityAssemblR functions
source("./R/build_even_community.R")
source("./R/link_taxa_abundances.R")
source("./R/reduce_taxa_abundances.R")
View(y)
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
  hub_list[[i]] <- link_taxa_abundances(dat = x,n.taxa = 50,relationship = "hub",link.scale = 3,n.links = links) 
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
assortativity
# plot betweenness
betweenness
# plot average path length
avg_path_length