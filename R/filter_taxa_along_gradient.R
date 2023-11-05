# reduce or increase taxa abundances along "abiotic gradient"


# variables
dat=even 
prop=.1
gradient.min=0 
gradient.max=1
groups=1 


filter_taxa_along_gradient <- function(dat, # community matrix with rows=samples & cols=taxa (raw abundance counts)
                                       prop=.1,  # proportion of taxa affected by this gradient
                                       gradient.min=0, # minimum gradient value (strength)
                                       gradient.max=1, # maximum gradient value (strength)
                                       gradient.strength=1, # multiplier for effect of gradient on taxa abundances
                                       groups=1, # N to split samples into (3*groups must be < nrow(dat)/groups)
                                       # remainder samples not within main group size will begin a new group
                                       association='negative', ## does the gradient increase or decrease abundance of affected taxa?
                                       transplant.only=FALSE, # if TRUE, only newly introduced taxa from a transplantation will be affected
                                       show.gradient=FALSE){ # if TRUE, a plot showing fitted gradient levels will be shown
  
  stopifnot("matrix" %in% class(dat))
  stopifnot(class(prop) %in% c("numeric","integer") & prop >= 0 & prop <= 1)
  stopifnot(class(gradient.min) %in% c("numeric","integer") & gradient.min >= 0 & gradient.min <= 1)
  stopifnot(class(gradient.max) %in% c("numeric","integer") & gradient.max >= 0 & gradient.max <= 1)
  stopifnot(class(groups) %in% c("numeric","integer") & nrow(dat)/groups >= 3)
  stopifnot(association %in% c("positive","negative"))
  stopifnot(class(gradient.strength) %in% c("numeric","integer") & gradient.strength >= 0)
  
  # get vector(s) for assigning gradient values to
  v <- 1:nrow(dat)
  main.group.length <- length(v) %/% groups
  remainder.group.length <- length(v) %% groups
  
  # build gradient sequence
  grad <- seq(gradient.min,gradient.max,length.out=main.group.length)
  main.grad <- rep(grad,(groups))
  remainder.grad <- grad[1:remainder.group.length]
  if(remainder.group.length == 0){
    gradient.values <- c(main.grad)
  }
  if(remainder.group.length > 0){
    gradient.values <- c(main.grad,remainder.grad)
    warning("Number of samples is not a multiple of the number of groups. Extra partial gradient group has been added.")
  }
  
  
  # add a touch of random noise
  gradient.values <- gradient.values + rnorm(length(gradient.values),0,.1)
  gradient.values <- sapply(gradient.values,FUN = function(x){ifelse(x < 0,0,x)})
  gradient.values <- sapply(gradient.values,FUN = function(x){ifelse(x > 1,1,x)})
  
  if(show.gradient){
    plot(gradient.values,main = "fitted gradient")
  }
    
  # select taxa to reduce
  newtaxa.names <- colnames(dat)[grepl("^newtaxon_",x=colnames(dat))]
  
  if(transplant.only){
    n.taxa <- round(prop * length(newtaxa.names))
    affected.taxa <- sample(newtaxa.names,n.taxa)
  }
  if(!transplant.only){
    n.taxa <- round(prop * ncol(dat))
    affected.taxa <- sample(colnames(dat),n.taxa)  
  }
  
  # get matrix of just affected taxa and transform
  if(association == 'positive'){
    dat[,affected.taxa] <- round(dat[,affected.taxa] * (1+gradient.values) * gradient.strength)
  }
  
  if(association == 'negative'){
    dat[,affected.taxa] <- round(dat[,affected.taxa] * (1-(gradient.values^(1/gradient.strength))))
  }
  
  return(dat)

}
