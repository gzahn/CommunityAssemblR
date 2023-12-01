#' Retrieve the new donor taxa from a simulated transplantation experiment
#'
#' Pulls out the transplanted donor taxa after a simulated transplantation process. If you perform this on a matrix that was produced with a transplant_with_() function, it will remove them from the 'transplanted' matrix and return them as-is. This allows you to perform successive transplant_() experiments, thereby simulating multiple ecological processes.
#'
#' @param transplanted.community Numeric matrix. The community matrix that is the result of a simulated transplantation process.
#' @param rm.empty.taxa Should taxa with no observations be removed? Default = FALSE.
#' @return List of two matrices. The 'new.donors' element contains the new donor community, and the 'new.recipients' element contains the updated resident community.
#'
#' @details
#' This returns a list of two new matrices which represent the updated taxon abundances for the donor and recipient. The purpose of this is so that you can then pass these to another transplantation simulation or other simulated ecological process. Thus, you can simulate a range of processes into one 'final' transplantation experiment.
#'
#' @examples
#' recipient <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#' donor <- build_donor_community(resident.comm = even,n.transplant.taxa = 10,overlap = .3)
#' trans <- transplant_w_niche_occupation(recipient, donor, n.niches = 1, niche.shape = 'normal',
#' niche.size.resident = 3, niche.size.donor = 1, abundance.threshold = .05,
#' process.type = "preemption")
#' new_communities <- retrieve_donor_taxa(trans)
#'
#' transplant_w_stochasticity(recipient = new_communities$new.recipients,recipient = new_communities$new.donor)
#'
#' @export

retrieve_donor_taxa <- function(transplanted.community, rm.empty.taxa = FALSE){

  stopifnot("matrix" %in% class(transplanted.community))
  if(!any(grepl("newtaxon",colnames(trans)))){
    stop("There are no newtaxa in your community matrix. Was this matrix the output from a simulated transplantation?")
  }

  x <- transplanted.community[,grepl("newtaxon",colnames(transplanted.community))]

  new.donors <- x

  if(rm.empty.taxa){
  new.donors <- x[,colSums(transplanted.community[,grepl("newtaxon",colnames(transplanted.community))]) > 0]
  }

  x <- transplanted.community[,!grepl("newtaxon",colnames(transplanted.community))]
  new.recipients <- x

  return(list(new.donors = new.donors,
       new.recipients = new.recipients))

}
