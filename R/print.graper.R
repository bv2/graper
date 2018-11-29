#' @title Print a graper object
#' @name print.graper
#' @description Function to print a fitted graper model.
#' @param x fitted graper model as obtained from  \code{\link{graper}}
#' @param ... additional print arguments
#' @return Print output.
#' @importFrom methods is
#' @export
#' @examples
#' # create data
#' dat <- makeExampleData()
#' fit <- graper(dat$X, dat$y, dat$annot)
#' print(fit)

print.graper <- function(x, ...){
  # sanity check
  if(!is(x, "graper")) {
    stop("object needs to be a graper object.")
  }

  cat(ifelse(x$Options$spikeslab, "Sparse", "Dense"), "graper object for a",
      ifelse(x$Options$family == "gaussian", "linear", "logistic"),
      "regression model with", nrow(x$EW_beta),
      "predictors in", length(unique(x$annot)), "groups.\n",
      "Group-wise shrinkage:\n",
      paste(unique(as.character(x$annot)), collapse="\t"), "\n",
      paste(round(x$EW_gamma, 2), collapse="\t"), "\n")

  if(x$Options$spikeslab) {
  cat("Group-wise sparsity (1 = dense, 0 = sparse):\n")
  cat(paste(unique(as.character(x$annot)), collapse="\t"), "\n")
  cat(paste(round(x$EW_pi, 2), collapse="\t"))
  }

}
