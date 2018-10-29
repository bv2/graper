#' @title Print a grpRR object
#' @name print.grpRR
#' @description Function to print a fitted grpRR model.
#' @param x fitted grpRR model as obtained from  \code{\link{grpRR}}
#' @param ... additional print arguments
#' @return Print output.
#' @export
#' @examples
#' # create data
#' dat <- makeExampleData()
#' fit <- grpRR(dat$X, dat$y, dat$annot)
#' print(fit)

print.grpRR <- function(x, ...){
  # sanity check
  if(class(x) != "grpRR") {
    stop("object needs to be a grpRR object.")
  }

  cat(paste(ifelse(x$Options$spikeslab, "Sparse", "Dense"), " grpRR object for a ", ifelse(x$Options$family == "gaussian", "linear", "logistic"),
              " regression model with ", nrow(x$EW_beta), " predictors in ", length(unique(x$annot))," groups.\n",
              "Group-wise shrinkage:\n", paste(unique(as.character(x$annot)), collapse="\t"), "\n",
              paste(round(x$EW_gamma,2), collapse="\t"), sep=""), "\n")
  if(x$Options$spikeslab) {
  cat(paste("Group-wise sparsity (1=dense, 0=sparse):\n"))
  cat(paste(unique(as.character(x$annot)), collapse="\t"), "\n")
  cat(paste(round(x$EW_pi,2), collapse="\t"))
  }
  
}
