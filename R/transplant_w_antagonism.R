#' Simulate community transplantation with recipient antagonism
#'
#' Manipulates a community matrix with samples as rows and taxa as columns. Donor taxa will be reduced correspondingly with the abundance of random antagonists present in the recipient community.
#'
#' @param recipient Numeric matrix. The starting community matrix that represents the resident community in N samples.
#' @param donor Numeric matrix. The donor community matrix that represents the donor community. Typically, the result of build_donor_community().
#' @param antag.ubiq Positive numeric vector of length 1, between 0 and 1. The proportion of your donor taxa that have antagonistic taxa in the recipient community.
#' @param antag.strength Positive numeric vector of length 1, between 0 and 1. Determines the scale of antagonism. i.e., how much reduction to do to donor taxa, proportionate to antagonist relative abundance.
#' @param antag.abundant Logical. If TRUE, recipient antagonists will be selected from the most abundant taxa. If FALSE, they will be selected randomly from all taxa.
#' @param print.antag.taxa Logical. If TRUE, the random recipient taxa selected as antagonists will be printed to the console, and also exported to an object named "tmp_facil.taxa" in the Global environment (yeah, I know that's sloppy).
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'.
#'
#' @examples
#' recipient <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#' donor <- build_donor_community(resident.comm = even,n.transplant.taxa = 10,overlap = .3)
#' transplant_w_antagonism(recipient, donor, antag.ubiq = .5, antag.strength = .5)
#'
#' @export

transplant_w_antagonism <- function(recipient,donor,
                                    antag.ubiq=.5,
                                    antag.strength=.1,
                                    antag.abundant=TRUE,
                                    print.antag.taxa=FALSE){

  scale01 <- function(x){
    if(max(x) == 0){return(x)}
    if(max(x) > 0){return((x - min(x)) / (max(x) - min(x)))}
  }

  stopifnot("matrix" %in% class(donor))
  stopifnot("matrix" %in% class(recipient))
  stopifnot(class(antag.ubiq) %in% c("numeric","integer") & antag.ubiq >= 0 & antag.ubiq <= 1)
  stopifnot(class(antag.strength) %in% c("numeric","integer") & antag.strength >= 0)



  # rescale donor community to same range as recipient community
  # scaled against mean from recipient
  donor <- round(apply(donor,2,scale01) * (round(median(colSums(recipient) / nrow(recipient)))))


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
  num.antagonists <- round(length(novel_taxa)*antag.ubiq)
  if(!antag.abundant){ # choose random facilitator taxa
    antag.taxa <- sample(colnames(recipient),num.antagonists)
  }
  if(antag.abundant){ # choose most abundant facilitator taxa
    antag.taxa <- names(colSums(recipient)[order(colSums(recipient),decreasing = TRUE)][1:num.antagonists])
  }
  antagonists <- recipient[,antag.taxa]

  if(print.antag.taxa == TRUE){
    print(antag.taxa)
    assign("tmp_antag.taxa", value = antag.taxa, envir = .GlobalEnv)
  }

  # assign them to a paired novel taxon from the donor community
  antagonized <- donor[,sample(novel_taxa,num.antagonists)]

  # reduce 'antagonized' taxa by multiplier based on antagonists' scaled abundances
  antag.multiplier <- 1 - (t(apply(recipient,1,function(x){x/sum(x)}))[,antag.taxa] * antag.strength)
  antagonized <- round(antagonized * antag.multiplier)
  antagonized[antagonized < 0] <- 0 # take care of any potential negative values
  newdonor <- donor[,novel_taxa]
  newdonor[,colnames(antagonized)] <- antagonized


  # combine reduced new taxa with augmented recipient community

  final_community <- cbind(newdonor,recipient)

  return(final_community)

}


