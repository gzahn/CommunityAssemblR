library(tidyverse)
library(igraph); packageVersion("igraph")
library(SpiecEasi); packageVersion("SpiecEasi")

x <- comm_even(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
y <- link_taxa_abundances(dat = x,n.taxa = 50,relationship = "hub",link.scale = 3,n.links = 5)
z <- link_taxa_abundances(dat = x,n.taxa = 2,relationship = "positive",link.scale = 3,n.links = 10) %>% 
  link_taxa_abundances(n.taxa = 50,relationship = "hub",link.scale = 3,n.links = 5)




se.params <- list(rep.num=20, ncores=(parallel::detectCores()-1))

View(x)
View(y)
y # altered with hub taxa
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

se.igraph.x <- adj2igraph(getRefit(se.mb.x))
se.igraph.y <- adj2igraph(getRefit(se.mb.y))
se.igraph.z <- adj2igraph(getRefit(se.mb.z))


am.coord.x <- layout.auto(se.igraph.x)
plot(se.igraph.x, layout=am.coord.x, vertex.label=NA)

am.coord.y <- layout.auto(se.igraph.y)
plot(se.igraph.y, layout=am.coord.y, vertex.label=NA)

am.coord.z <- layout.auto(se.igraph.z)
plot(se.igraph.z, layout=am.coord.z)

igraph::degree.distribution(se.igraph.x) %>% plot
igraph::degree.distribution(se.igraph.y) %>% plot
hubx <- igraph::hub_score(se.igraph.x)
huby <- igraph::hub_score(se.igraph.y)
hubz <- igraph::hub_score(se.igraph.z)
hubx$value
huby$value
hubz$value
artpoints.x <- igraph::articulation.points(se.igraph.x)
artpoints.y <- igraph::articulation.points(se.igraph.y)
artpoints.z <- igraph::articulation.points(se.igraph.z)
y[,artpoints.y]

igraph::assortativity_degree(se.igraph.z)


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

