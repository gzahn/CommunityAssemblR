
# plot igraph, showing hub taxa articulation points as larger dots, sized by authority score
plot_hubs <- function(graph,bigpoint=10,littlepoint=3){
  am.coord <- layout.auto(graph)
  art <- articulation_points(graph)
  plot(graph,
       layout = am.coord,
       vertex.size=(scale01(abs(igraph::authority_score(graph)$vector)) * 10)+3,
         #ifelse(1:length(graph) %in% artpoints.z, bigpoint, littlepoint),
       vertex.label=NA)
}


# custom scaling function (0,1) that can handle vectors of all zeros
scale01 <- function(x){
  if(max(x) == 0){return(x)} # if all zeros, do nothing
  if(max(x) > 0){return((x - min(x)) / (max(x) - min(x)))} # if some positive, scale as normal
  if(sum(x) < 0){x <- abs(x) # not sure about this
                 scaled <- (x - min(x)) / (max(x) - min(x))
                 return(-scaled)} # if all negative, use absolute values, then replace negative sign
}
