#' Simulate community transplantation with recipient antagonism
#'
#' Manipulates a community matrix with samples as rows and taxa as columns. Donor taxa abundances will be increased or decreased correspondingly with the niche occupancy of taxa in the recipient community.
#'
#' @param recipient Numeric matrix. The starting community matrix that represents the resident community in N samples.
#' @param donor Numeric matrix. The donor community matrix that represents the donor community. Typically, the result of build_donor_community().
#' @param n.niches Positive whole numeric vector of length 1. Number of rounds of niche selection to perform.
#' @param niche.shape Character vector of length 1. Currently only "normal" is supported. The niche distribution shape to draw from for niche occupation values.
#' @param niche.size.resident Positive whole numeric vector of length 1. The number of taxa in resident community that occupy the niche. Default = 3.
#' @param niche.size.donor Positive whole numeric vector of length 1. The number of taxa in donor community that occupy the niche. Default = 1.
#' @param abundance.threshold Positive numeric vector of length 1, between 0 and 1. The relative abundance threshold for resident taxa in niche to trigger effect. This will be mean of distribution...if niche.size.resident > 1, thresholds for additional taxa will be based on the distribution of niche.shape. Default = 0.05.
#' @param process.type Character vector of length 1. Must currently be in c("preemption","permission"). Under "preemption", when the niche occupancy threshold is reached in a given sample, the resident taxon will begin to decrease the abundance of the donor taxon. Under "permission", when that threshold is reached, the resident taxon will enable the presence of the donor taxon (where the threshold is not reached, the donor taxon abundance will drop to 0). Default = "preemption".
#' @param print.niche.taxa Logical. If TRUE, the random recipient and donor taxa selected to occupy niches will be printed to the console, and also exported to a list object named "tmp_selected_niche.taxa", in the Global environment (yeah, I know that's sloppy).
#' @return Taxon abundance matrix with samples as rows and taxa as columns. class='matrix'.
#'
#' @examples
#' recipient <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#' donor <- build_donor_community(resident.comm = even,n.transplant.taxa = 10,overlap = .3)
#' transplant_w_niche_occupation(recipient, donor, n.niches = 1, niche.shape = 'normal',
#' niche.size.resident = 3, niche.size.donor = 1, abundance.threshold = .05,
#' process.type = "preemption")
#'
#' @export

transplant_w_niche_occupation <- function(recipient,
                                          donor,
                                          n.niches = 1, # number of rounds of niche selection
                                          niche.shape = 'normal',
                                          niche.size.resident = 3, # number of taxa in resident community that occupy the niche
                                          niche.size.donor = 1, # number of taxa in donor community (new taxa only) that occupy the niche
                                          abundance.threshold = .05, # relative abundance threshold for resident taxa in niche to trigger effect
                                          # this will be mean of distribution...if niche.size.resident > 1, thresholds for additional taxa
                                          # will be based on the distribution of niche.shape
                                          process.type = "preemption", # 'preemption' or 'permission'
                                          print.niche.taxa = FALSE){

  stopifnot("matrix" %in% class(donor))
  stopifnot("matrix" %in% class(recipient))

  if(!niche.shape %in% c("normal")){
    stop("Currently, only 'normal' distribution accepted for niche.shape")
  }
  if(abundance.threshold <= 0 | abundance.threshold > 1){
    stop("abundance.threshold must be between 0 and 1")
  }

  stopifnot(class(niche.size.donor) == "numeric")
  niche.size.donor <- round(niche.size.donor)
  stopifnot(niche.size.donor > 0)

  stopifnot(class(niche.size.resident) == "numeric")
  niche.size.resident <- round(niche.size.resident)
  stopifnot(niche.size.resident > 0)

  stopifnot(class(n.niches) == "numeric")
  n.niches <- round(n.niches)
  ifelse(n.niches == 0,1,n.niches)
  stopifnot(n.niches > 0)

  if(!process.type %in% c("preemption","permission")){
    stop("process.type must be either 'preemtion' or 'permission'.")
  }

  # rescale donor community to same range as recipient community
  # scaled against mean from recipient
  donor <- round(apply(donor,2,scale01) * (round(median(colSums(recipient)))))

  # where taxa overlap, just take sum of both samples
  shared_donor <- donor[,colnames(donor) %in% colnames(recipient)]
  shared_recipient <- recipient[,colnames(recipient) %in% colnames(shared_donor)]
  new_shared_taxa <- shared_donor + shared_recipient

  recipient[,colnames(new_shared_taxa)] <- new_shared_taxa

  # find novel taxa from donor community
  novel_taxa <- colnames(donor)[!colnames(donor) %in% colnames(recipient)]
  # rearrange in order because I'm OCD
  novel_taxa <- novel_taxa[order(as.numeric(sub("newtaxon_","",x = novel_taxa)))]
  novel.taxa.matrix <- donor[,novel_taxa]


  ### Niche occupation process, repeat n.niches times

  # fill this list with taxa names, if print.niche.taxa==TRUE
  if(print.niche.taxa == TRUE){
    selected.taxa <- list()
  }

  for(X in 1:n.niches){
    # Select resident abundance thresholds for N taxa, based on niche.shape distribution
    abund.thresholds <- abs(rnorm(n=niche.size.resident,mean = abundance.threshold,sd = .05))

    # select recipient taxa to act as preemtion or permission gates
    resident.niche.taxa <- sample(colnames(recipient),niche.size.resident)

    # select donor taxa that will be influenced by resident niche occupation
    donor.niche.taxa <- sample(novel_taxa,niche.size.donor)

    # if print.niche.taxa, store them in nested list
    if(print.niche.taxa == TRUE){
      selected.taxa[["resident"]][[paste0("round_",X)]] <- resident.niche.taxa
      selected.taxa[["donor"]][[paste0("round_",X)]] <- donor.niche.taxa
      print(resident.niche.taxa); print(donor.niche.taxa)
    }

    # convert to relative abundance
    recipient_ra <- t(apply(recipient,1,function(x){x/sum(x)}))

    # subset to just the resident niche taxa
    ra <- recipient_ra[,resident.niche.taxa]

    # compare each column with respective abundance threshold, give T/F
    tests <- matrix(nrow = nrow(ra),ncol = niche.size.resident)
    for(i in 1:niche.size.resident){
      tests[,i] <- ra[,i] > abund.thresholds[i]
    }

    # if any one of the rows meets threshold, then effect can take place (either prevention or allowance)
    met.threshold <- as.logical(rowSums(tests))

    # pre-emption: if threshold is met in resident community, it prevents donor taxa from occupying
    if(process.type == "preemption"){
      novel.taxa.matrix[,donor.niche.taxa] <- (novel.taxa.matrix[,donor.niche.taxa] * !met.threshold)
    }

    # permission: if threshold is met in resident community, it allows donor taxa to occupy
    if(process.type == "permission"){
      novel.taxa.matrix[,donor.niche.taxa] <- (novel.taxa.matrix[,donor.niche.taxa] * met.threshold)
    }

  }

  final_community <- cbind(novel.taxa.matrix,recipient)

  if(print.niche.taxa == TRUE){
    assign("tmp_selected_niche.taxa",selected.taxa,envir = .GlobalEnv)
  }

  return(final_community)

}

