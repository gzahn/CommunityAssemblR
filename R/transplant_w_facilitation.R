#' Simulate community transplantation with recipient facilitation
#'
#' Manipulates a community matrix with samples as rows and taxa as columns. Donor taxa will be increased correspondingly with the abundance of random facilitative taxa present in the recipient community.
#'
#' @param recipient Numeric matrix. The starting community matrix that represents the resident community in N samples.
#' @param donor Numeric matrix. The donor community matrix that represents the donor community. Typically, the result of build_donor_community().
#' @param facil.ubiq Positive numeric vector of length 1, between 0 and 1. The proportion of your donor taxa that have facilitative taxa in the recipient community.
#' @param facil.strength Positive numeric vector of length 1, between 0 and 1. Determines the scale of facilitation. i.e., how much increase to do to donor taxa, proportionate to facilitator relative abundance.
#' @param facil.abundant Logical. If TRUE, recipient facilitators will be selected from the most abundant taxa. If FALSE, they will be selected randomly from all taxa.
#' @param print.facil.taxa Logical. If TRUE, the random recipient taxa selected as facilitators will be printed to the console, and also exported to an object named "tmp_facil.taxa" in the Global environment (yeah, I know that's sloppy).
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'.
#'
#' @examples
#' recipient <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#' donor <- build_donor_community(resident.comm = even,n.transplant.taxa = 10,overlap = .3)
#' transplant_w_facilitation(recipient, donor, facil.ubiq = .5, facil.strength = .5)
#'
#' @export


transplant_w_facilitation <- function(recipient,donor,
                                      facil.ubiq=.5,
                                      facil.strength=1,
                                      facil.abundant=TRUE,
                                      print.facil.taxa=FALSE){

  scale01 <- function(x){
    if(max(x) == 0){return(x)}
    if(max(x) > 0){return((x - min(x)) / (max(x) - min(x)))}
  }

  stopifnot("matrix" %in% class(donor))
  stopifnot("matrix" %in% class(recipient))
  stopifnot(class(facil.ubiq) %in% c("numeric","integer") & facil.ubiq >= 0 & facil.ubiq <= 1)
  stopifnot(class(facil.strength) %in% c("numeric","integer") & facil.strength >= 0)

  # rescale donor community to same range as recipient community
  # scaled against mean from recipient
  donor <- round(apply(donor,2,scale01) * (round(median(colSums(recipient)))))

  # where taxa overlap, just take sum of both samples
  shared_donor <- donor[,colnames(donor) %in% colnames(recipient)]
  shared_recipient <- recipient[,colnames(recipient) %in% colnames(shared_donor)]
  new_shared_taxa <- shared_donor + shared_recipient

  recipient[,colnames(new_shared_taxa)] <- new_shared_taxa

  # where new taxa arriving, reduce the new taxa, scaled by facil.strength
  novel_taxa <- colnames(donor)[!colnames(donor) %in% colnames(recipient)]
  # rearrange in order because I'm OCD
  novel_taxa <- novel_taxa[order(as.numeric(sub("newtaxon_","",x = novel_taxa)))]

  # select antagonists from resident community
  num.facilitators <- round(length(novel_taxa)*facil.ubiq)
  if(!facil.abundant){ # choose random facilitator taxa
    facil.taxa <- sample(colnames(recipient),num.facilitators)
  }
  if(facil.abundant){ # choose most abundant facilitator taxa
    facil.taxa <- names(colSums(recipient)[order(colSums(recipient),decreasing = TRUE)][1:num.facilitators])
  }
  facilitators <- recipient[,facil.taxa]

  if(print.facil.taxa == TRUE){
    print(facil.taxa)
    assign("tmp_facil.taxa", value = facil.taxa, envir = .GlobalEnv)
  }

  # assign them to a paired novel taxon from the donor community
  facilitated <- donor[,sample(novel_taxa,num.facilitators)]

  # increase 'facilitated' taxa by multiplier based on antagonists' scaled abundances
  facil.multiplier <- 1 + (t(apply(recipient,1,function(x){x/sum(x)}))[,facil.taxa] * facil.strength)
  facilitated <- round(facilitated * facil.multiplier)

  newdonor <- donor[,novel_taxa]
  newdonor[,colnames(facilitated)] <- facilitated


  # combine reduced new taxa with augmented recipient community

  final_community <- cbind(newdonor,recipient)

  return(final_community)

}
