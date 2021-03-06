#' @title Plot posterior distributions
#' @name plotPosterior
#' @description Function to plot the posterior of the model parameters
#'  obtained by graper from the variational inference framework.
#' @param fit fit as produced by \code{\link{graper}}
#' @param param2plot which parameter of the graper model to
#' plot (gamma, beta, tau or s)
#' @param beta0 true beta (if known)
#' @param gamma0 true gamma (if known)
#' @param pi0 true pi (if known)
#' @param tau0 true tau (if known)
#' @param s0 true s (if known)
#' @param jmax maximal number of components per
#' group to plot (for beta and s)
#' @param range plotting range (x-axis)
#' @export
#' @importFrom stats dbinom dnorm dgamma
#' @importFrom methods is
#' @import ggplot2
#' @return a ggplot object
#' @examples
#' # create data
#' dat <- makeExampleData()
#' # fit the graper model
#' fit <- graper(dat$X, dat$y, dat$annot)
#' # plot posterior distribution of the gamma parameter
#' plotPosterior(fit, param2plot="gamma")

plotPosterior <- function(fit, param2plot, beta0 = NULL,
                            gamma0 = NULL, tau0 = NULL, pi0=NULL,
                            s0=NULL, jmax=2, range=NULL) {

    # sanity check
    if(!is(fit, "graper")) {
        stop("fit needs to be a graper object.")
    }

    if (param2plot == "beta") {
        .plotBeta(fit=fit, beta0=beta0, range=range, jmax=jmax)
    } else if (param2plot=="s") {
        .plotS(fit=fit, s0=s0, range=range, jmax=jmax)
    } else if (param2plot=="pi") {
        .plotPi(fit=fit, pi0=pi0, range=range)
    } else if (param2plot=="gamma") {
        .plotGamma(fit, gamma0=gamma0, range=range)
    } else if (param2plot=="tau") {
        .plotTau(fit, tau0=tau0, range=range)
    } else {
        stop("param2plot needs to be beta, s, pi, gamma or tau.")
    }
}

#' @title Plot evidence lower bound
#' @name plotELBO
#' @description Function to plot the evidence lower bound (ELBO)
#' over iterations to monitor the convergence of the algorithm.
#' @param fit fit as produced by \code{\link{graper}}
#' @import ggplot2
#' @importFrom methods is
#' @export
#' @return a ggplot object
#' @examples
#' dat <- makeExampleData()
#' fit <- graper(dat$X, dat$y, dat$annot)
#' plotELBO(fit)

plotELBO <- function(fit){
    # sanity check
    if(!is(fit, "graper")) {
        stop("fit needs to be a graper object.")
    }
    iteration <- ELBO <- NULL # avoid notes in check
    if(is.null(fit$ELB_trace)) {
        stop("ELBO was not computed for this fit.")
    }
    df <- data.frame(iteration = seq_along(fit$ELB_trace),
                    ELBO = fit$ELB_trace)
    ggplot(df, aes(x=iteration, y=ELBO)) + geom_line()
}

#' @title Plot group-wise penalties
#' @name plotGroupPenalties
#' @description Function to plot the group-wise penalty factors
#'  (gamma) and sparsity levels.
#' @param fit fit as produced by \code{\link{graper}}
#' @import ggplot2
#' @importFrom methods is
#' @importFrom cowplot plot_grid
#' @export
#' @return a ggplot object
#' @examples
#' dat <- makeExampleData()
#' fit <- graper(dat$X, dat$y, dat$annot)
#' plotGroupPenalties(fit)

plotGroupPenalties <- function(fit){
    # sanity check
    if(!is(fit, "graper")) {
        stop("fit needs to be a graper object.")
    }
    group <- gamma <- pi <- NULL  # avoid notes in check
    df <- data.frame(group = unique(fit$annot),
                    gamma = fit$EW_gamma)
    gg1 <- ggplot(df, aes(x=group, y=gamma, fill=group)) +
        geom_bar(stat="identity") +
        ggtitle("Group-wise penalty factor") +
        theme_bw()

    if(fit$Options$spikeslab) {
        df$pi = fit$EW_pi
        gg2 <- ggplot(df, aes(x=group, y=pi, fill=group)) +
            geom_bar(stat="identity") +
            ggtitle("Group-wise sparsity level") +
            ylim(c(0,1)) +
            theme_bw()
        cowplot::plot_grid(gg1, gg2)
    } else {
        print(gg1)
    }
}


.plotGamma <- function(fit, gamma0, range){
        # avoid notes on global varibale binding in check
        gamma <- true_gamma <- density <- k <- mean_val <- NULL

        # extract fitted parameters
        gr <- lapply(seq_along(fit$alpha_gamma), function(k) {
            mean_val <- fit$EW_gamma[k]
            if(is.null(range)) {
                x <- seq(0, max(gamma0[k], 5 * mean_val), , length.out=1000)
            } else {
                x <- seq(range[1], range[2], length.out = 1000)
            }
            data.frame(gamma = x, k = k, mean_val = mean_val)
        })
        gr <- do.call(rbind, gr)
        gr$density <- stats::dgamma(gr$gamma, fit$alpha_gamma[gr$k],
                                    fit$beta_gamma[gr$k])
        if(!is.null(gamma0)) {
            gr$true_gamma <- gamma0[gr$k]
        }

        # plot
        gg <- ggplot(gr, aes(x=gamma, y=density)) +
            geom_line() + facet_wrap(~k, scales="free") +
            geom_vline(aes(xintercept=mean_val, col="mean"), linetype="dashed")
        if(!is.null(gamma0)) {
            gg <- gg + geom_vline(aes(xintercept=true_gamma, col="true_gamma"),
                                    linetype="dashed")
        }
        gg <- gg + scale_color_manual(name="statistics",
                            values=c(true_gamma = "blue", mean = "red"))
        print(gg)
}


.plotS <- function(fit, s0, range, jmax){
        if(!fit$Options$spikeslab) {
            stop("fit needs to be a sparse graper object.
                Use spikeslab = TRUE in graper.")
        }
        s <- true_s <- density <- mean_val <- psi <-  NULL

        # extract fitted parameters
        message("Only plotting the first ", jmax, " s per group.")
        js <- Reduce(c,lapply(unique(fit$annot),
                        function(gr) {
                            which(fit$annot==gr)[seq_len(min(jmax,
                                                    sum(fit$annot==gr)))]
                        }))
        gr <- lapply(js, function(j) {
            psi <- fit$EW_s[j]
            mean_val <- psi
            group <- fit$annot[j]
            x <- c(0,1)
            data.frame(s = x, j = j, psi = psi,
                mean_val = mean_val, group = group)
        })
        gr <- do.call(rbind, gr)
        gr$density <- gr$psi ^ gr$s * (1 - gr$psi) ^ (1 - gr$s)
        if(!is.null(s0)) {
            gr$true_s <- s0[gr$j]
        }

        # plot
        gg <- ggplot(gr, aes(xend=s, x=s, yend=0, y=density)) +
            geom_segment() +geom_point() + facet_wrap(~j, ncol=jmax) +
            geom_vline(aes(xintercept=mean_val, col="mean"), alpha=0.5,
                        linetype="dashed")
        if(!is.null(s0)) {
            gg <- gg + geom_vline(aes(xintercept=true_s, col="true_s"),
                                    linetype="dashed", alpha=0.5)
        }
        gg <- gg + scale_color_manual(name = "statistics",
                            values = c(true_s = "blue", mean = "red")) +
                theme_bw() + theme(strip.background = element_blank(),
                                    strip.text.x = element_blank())
        print(gg)
}


.plotPi <- function(fit, pi0, range){
        if(!fit$Options$spikeslab) {
            stop("fit needs to be a sparse graper object.
                Use spikeslab = TRUE in graper.")
        }
        # avoid notes on global varibale binding in check
        pi <- true_pi <- density <- mean_val <- k <- NULL

        # extract fitted parameters
        gr <- lapply(seq_along(fit$EW_pi), function(k) {
            mean_val <- fit$EW_pi[k]
            x <- seq(0, 1, length.out=1000)
            data.frame(pi = x, k = k, mean_val = mean_val)
        })
        gr <- do.call(rbind, gr)
        gr$density <- stats::dbeta(gr$pi, fit$alpha_pi[gr$k], fit$beta_pi[gr$k])
        if(!is.null(pi0)) {
            gr$true_pi <- pi0[gr$k]
        }

        # plot
        gg <- ggplot(gr, aes(x=pi, y=density)) + geom_line() +
            facet_wrap(~ k, scales="free") +
            geom_vline(aes(xintercept=mean_val, col="mean"), linetype="dashed")
        if(!is.null(pi0)) {
            gg <- gg + geom_vline(aes(xintercept=true_pi, col="true_pi"),
                                    linetype="dashed")
        }
        gg <- gg +  scale_color_manual(name = "statistics",
                                values = c(true_pi = "blue", mean = "red"))
        print(gg)
}

.plotBeta <- function(fit, beta0, range, jmax){
        # avoid notes on global varibale binding in check
        beta <- true_beta <- density <- mean_val <- mu_slab <- NULL
        va_slab <- psi <- j <- NULL

        # extract fitted parameters
        message("Only plotting the first ", jmax, " coefficients per group.")
        js <- Reduce(c,lapply(unique(fit$annot),
            function(gr) {
                which(fit$annot == gr)[seq_len(min(jmax, sum(fit$annot == gr)))]
            }))
        gr <- lapply(js, function(j) {
            va_slab <- fit$sigma2_tildebeta_s1[j]
            mu_slab <- fit$EW_tildebeta_s1[j]
            mean_val <- fit$EW_beta[j]
            psi <- fit$EW_s[j]
            group <- fit$annot[j]
            if(is.null(range)) {
                x <- c(0, seq(min(0, mu_slab - abs(mu_slab) / 2, beta0[j]),
                        max(mu_slab + abs(mu_slab) / 2, beta0[j], 0),
                        length.out=1000))
            } else {
                x <- c(0, seq(range[1], range[2], length.out=1000))
            }
            data.frame(beta = x, j = j, va_slab = va_slab,
                mean_val = mean_val, mu_slab = mu_slab,
                psi = psi, group = group)
            })
        gr <- do.call(rbind, gr)
        gr$density <- gr$psi * stats::dnorm(gr$beta, gr$mu_slab,
                                sqrt(gr$va_slab)) + (1 - gr$psi)*(gr$beta == 0)
        if(!is.null(beta0)) {
            gr$true_beta <- beta0[gr$j]
        }
        # plot
        gg <- ggplot(gr, aes(x=beta, y=density)) + geom_line() +
            facet_wrap(~j, scales="free", ncol=jmax) +
            geom_vline(aes(xintercept=mean_val, col="mean"),
                    alpha=0.5, linetype="dashed")
        if(!is.null(beta0)) {
            gg <- gg + geom_vline(aes(xintercept=true_beta, col="true_beta"),
                                linetype="dashed", alpha=0.5)
        }
        gg <- gg + scale_color_manual(name="statistics",
                                values=c(true_beta = "blue", mean = "red")) +
            theme(strip.background=element_blank(),
                strip.text.x=element_blank())
        print(gg)
}


.plotTau <- function(fit, tau0, range){
        if(fit$Options$family == "binomial") {
            stop("fit is from a logistic regression model. No tau to plot.")
        }
        # avoid notes on global varibale binding in check
        tau <- true_tau <- density <- mean_val <- NULL

        # extract parameters
        mean_val <- fit$EW_tau
        if(is.null(range)) {
            x <- seq(0, max(tau0,5 * mean_val), , length.out=1000)
        } else {
            x <- seq(range[1], range[2], length.out=1000)
        }
        df <- data.frame(tau = x,
                        density = stats::dgamma(x, shape=fit$alpha_tau,
                                                rate=fit$beta_tau),
                        mean_val = mean_val)
        if(!is.null(tau0)) {
            df$true_tau <- tau0
        }

        # plot
        gg <- ggplot(df, aes(x=tau, y=density)) + geom_line() +
            geom_vline(aes(xintercept=mean_val, col="mean"), linetype="dashed")
        if(!is.null(tau0)) {
            gg <- gg + geom_vline(aes(xintercept=true_tau, col="true_tau"),
                                linetype="dashed")
        }
        gg <- gg + scale_color_manual(name="statistics",
                            values=c(true_tau = "blue", mean = "red"))
        print(gg)
}
