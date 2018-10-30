
#' @title Get posterior inclusion probabilities per feature
#' @name getPIPs
#' @description Function to obtain estimated posterior inclusion probabilities per feature from a fitted grpRR model.
#' @param object fitted grpRR model as obtained from  \code{\link{grpRR}}
#' @return 1-Column matrix of estimated posterior inclusion probabilities.
#' @importFrom methods is
#' @export
#' @examples
#' # create data
#' dat <- makeExampleData()
#' fit <- grpRR(dat$X, dat$y, dat$annot)
#' getPIPs(fit)

getPIPs <- function(object){
  # sanity check
  if(!is(object, "grpRR")){
    stop("object needs to be a grpRR object.")
  }
  if(!object$Options$spikeslab) {
    stop("object needs to be a sparse grpRR object. Use spikeslab = TRUE in grpRR.")
  }
  pips <- object$EW_s
  if(!is.null( object$Options$featurenames)){
    rownames(pips) <- object$Options$featurenames
  }

  return(pips)
}

