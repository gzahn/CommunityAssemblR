#' Build a community matrix to simulate community transplantation
#'
#' Convenience funtion that gathers up donor microbiome community relative abundances before and after transplanting
#' Returns a data frame that can be easily used to examine the success of microbial community establishement
#'
#' @param donor Numeric matrix. The donor community, typically created with make_donor_community()
#' @param transplant Numeric matrix. The final community that is the result of a transplant_w_() function.
#' @param newtaxa.only Logical. If newtaxa.only=TRUE, only novel taxa that persisted after transplantation will be retained.
#' @param tidy Logical. If tidy=TRUE, the function will return a 'tidy' data set, as opposed to a 'wide' one. The tidy version has columns named c("timepoint","sample","taxon","relative_abundance"). Default=TRUE.
#' @param keep.all.taxa Logical. If keep.all.taxa=TRUE, all overlapping taxa between the donor and recipient communities will be used to calculate relative abundance values, not just the taxa that are eventually kept. Default=FALSE.
#'
#'
#' @return Either a 'wide' or 'long' data frame with relative abundance values from the original donor community and the associated values for those taxa after transplantation.
#'
#' @examples
#' comm <- build_even_community(n.taxa = 100,n.samples = 44,n.reads = 3000, taxa.sd = 30)
#' donor <- build_donor_community(resident.comm = comm, n.transplant.taxa = 30,overlap = .75)
#' transplant <- transplant_w_niche_occupation(comm,donor)
#' check_transplant_success(donor,transplant,newtaxa.only=TRUE)
#'
#' @export

check_transplant_success <- function(donor,transplant,newtaxa.only=FALSE,tidy=TRUE,keep.all.taxa=TRUE){

  if(keep.all.taxa == TRUE){
    donor_ra <- t(apply(donor,1,function(x){x/sum(x)}))
    recipient_ra <- t(apply(transplant,1,function(x){x/sum(x)}))
  }

  if(keep.all.taxa == FALSE){
    donor_ra <- t(apply(donor,1,function(x){x/sum(x)}))
    recipient_ra <- t(apply(transplant[,colnames(donor_ra)],1,function(x){x/sum(x)}))
  }

  # if newtaxa.only=TRUE, subset to just the "newtaxa"
  if(newtaxa.only == TRUE){
    initial <- donor_ra[,grep(colnames(donor_ra),pattern="newtaxon")]
    final <- recipient_ra[,colnames(initial)]
  }

  # in newtaxa.only=FALSE, then just subset to the taxa found in the donor community
  if(newtaxa.only == FALSE){
    initial <- donor_ra
    final <- recipient_ra[,colnames(initial)]
  }

  # convert to data.frames
  initial <- as.data.frame(initial)
  final <- as.data.frame(final)

  # add columns indicating timepoint
  initial$timepoint <- "Initial"
  final$timepoint <- "Final"

  # join data frames
  initial <- initial[,names(initial)[order(names(initial))]]
  final <- final[,names(initial)[order(names(initial))]]
  df <- rbind(initial,final)

  df$sample <- rep(row.names(initial),2)

  if(tidy){
    df <- reshape(df,grep("tax",names(df),value = TRUE),direction = 'long',
                  v.names = "realative_abundance",times = grep("tax",names(df),value = TRUE))
    df$id <- NULL
    row.names(df) <- NULL
    names(df) <- c("timepoint","sample","taxon","relative_abundance")
  }

  return(df)
}
