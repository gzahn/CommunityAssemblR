#' Simulate community transplantation with stochastisicity in donor taxa abundances
#'
#' Manipulates a community matrix with samples as rows and taxa as columns. Donor taxa that overlap with recipient taxa will be joined additively. Novel taxa will be added with random noise.
#'
#' @param recipient Numeric matrix. The starting community matrix that represents the resident community in N samples.
#' @param donor Numeric matrix. The donor community matrix that represents the donor community. Typically, the result of build_donor_community().
#' @param stochasticity Positive numeric vector of length 1, between 0 and 1. The level of random noise to add to novel taxa abundances. If set to 0, no random noise will be added. Default=.05
#'
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'.
#'
#' @examples
#' recipient <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#' donor <- build_donor_community(resident.comm = recipient,n.transplant.taxa = 10,overlap = .3)
#' transplant_w_stochasticity(recipient, donor, stochasticity=0.05)
#'
#' @export

transplant_w_stochasticity <- function(recipient,donor,stochasticity=.05){

  scale01 <- function(x){
    if(max(x) == 0){return(x)}
    if(max(x) > 0){return((x - min(x)) / (max(x) - min(x)))}
  }

  stopifnot("matrix" %in% class(donor))
  stopifnot("matrix" %in% class(recipient))
  stopifnot(class(stochasticity) %in% c("numeric","integer") & stochasticity >= 0 & stochasticity <= 1)

  # rescale donor community to same range as recipient community
  # scaled against mean from recipient
  donor <- round(apply(donor,2,scale01) * (round(median(colSums(recipient) / nrow(recipient)))))


  # where taxa overlap, just take sum of both samples
  shared_donor <- donor[,colnames(donor) %in% colnames(recipient)]
  shared_recipient <- recipient[,colnames(recipient) %in% colnames(shared_donor)]
  new_shared_taxa <- shared_donor + shared_recipient

  recipient[,colnames(new_shared_taxa)] <- new_shared_taxa

  # where new taxa
  novel_taxa <- colnames(donor)[!colnames(donor) %in% colnames(recipient)]
  # rearrange in order because I'm OCD
  novel_taxa <- novel_taxa[order(as.numeric(sub("newtaxon_","",x = novel_taxa)))]

  # add random noise to novel taxa
  novel_taxa_dat <- donor[,novel_taxa]
  if(stochasticity > 0){
    multiplier <- matrix(rnorm(length(novel_taxa_dat),sd = 1),ncol=ncol(novel_taxa_dat),byrow = FALSE) * stochasticity
    new.values <- round(novel_taxa_dat * (1+multiplier))
    new.values <- (ifelse(new.values < 0, 0, new.values))
  }
  if(stochasticity == 0 ){
    new.values <- novel_taxa_dat
  }

  # combine reduced new taxa with augmented recipient community
  final_community <- cbind(new.values,recipient)

  return(final_community)
}


