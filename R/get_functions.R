
#' @title Get posterior inclusion probabilities per feature
#' @name getPIPs
#' @description Function to obtain estimated posterior inclusion
#'  probabilities per feature from a fitted graper model.
#' @param object fitted graper model as obtained from  \code{\link{graper}}
#' @return 1-Column matrix of estimated posterior inclusion probabilities.
#' @importFrom methods is
#' @export
#' @examples
#' # create data
#' dat <- makeExampleData()
#' fit <- graper(dat$X, dat$y, dat$annot)
#' getPIPs(fit)

getPIPs <- function(object){
    # sanity check
    if(!is(object, "graper")){
        stop("object needs to be a graper object.")
    }
    if(!object$Options$spikeslab) {
        stop("object needs to be a sparse graper object.
            Use spikeslab = TRUE in graper.")
    }

    pips <- object$EW_s
    if(!is.null(object$Options$featurenames)){
        rownames(pips) <- object$Options$featurenames
    }
    return(pips)
}

