# plotting_functions

#' @title plotPosterior 
#' @name plotPosterior
#' @description plot the posterior of the model parameters
#' @param fit fit as produced by fit_grpRR
#' @param param2plot which parametet to plot (gamma, beta, tau or s)
#' @param beta0 true beta (if known)
#' @param gamma0 true gamma (if known)
#' @param tau0 true tau (if known)
#' @param s0 true s (if known)
#' @param jmax maximal number of betas per group to plot
#' @import 
#' @export
# ---------------------------
plotPosterior <- function(fit, param2plot, beta0 = NULL, gamma0 = NULL, tau0 = NULL, pi0=NULL, s0=NULL, jmax=2, range=NULL) {

    if (param2plot=="beta") {
        message("Only plotting the first ", jmax, " coefficients per group.")
        js <- sapply(unique(fit$annot), function(gr) which(fit$annot==gr)[1:min(jmax, sum(fit$annot==gr))])
        gr <- lapply(js, function(j) {
            va_slab <- fit$sigma2_tildebeta_s1[j]
            mu_slab <- fit$EW_tildebeta_s1[j]
            mean = fit$EW_beta[j]
            psi <- fit$EW_s[j]
            group <- fit$annot[j]
            if(is.null(range)) x <- c(0,seq(min(0, mu_slab - 200*va_slab, beta0[j]), max(mu_slab + 200*va_slab, beta0[j]), length.out=1000))
            else x <- c(0,seq(range[1], range[2], length.out = 1000))
            data.frame(beta=x, j=j, va_slab=va_slab, mean=mean, mu_slab=mu_slab, psi=psi, group=group)
            }) %>% bind_rows()
        gr %<>% mutate(density = psi* dnorm(beta, mu_slab, sqrt(va_slab)) + (1-psi)*(beta==0))
        if(!is.null(beta0)) gr %<>% mutate(true_beta = beta0[j])
        gg <- ggplot(gr, aes(x=beta, y=density)) + geom_line() + facet_wrap(~j, scales="free", ncol=jmax) +
         geom_vline(aes(xintercept=mean, col="mean"),alpha=0.5, linetype="dashed")
        if(!is.null(beta0)) gg <- gg + geom_vline(aes(xintercept=true_beta, col="true_beta"), linetype="dashed",alpha=0.5) 
        gg <- gg+ scale_color_manual(name = "statistics", values = c(true_beta = "blue", mean = "red"))
        print(gg)
    }

        if (param2plot=="s") {
        message("Only plotting the first ", jmax, " s per group.")
        js <- sapply(unique(fit$annot), function(gr) which(fit$annot==gr)[1:min(jmax, sum(fit$annot==gr))])
        gr <- lapply(js, function(j) {
            psi <- fit$EW_s[j]
            mean <- psi
            group <- fit$annot[j]
            x <- c(0,1)
            data.frame(s=x, j=j, psi=psi,mean=mean, group=group)
            }) %>% bind_rows()
        gr %<>% mutate(density = psi^s *(1-psi)^(1-s))
        if(!is.null(pi0)) gr %<>% mutate(true_s = pi0[j])
        gg <- ggplot(gr, aes(xend=s, x=s, yend=0, y=density)) + geom_segment() +geom_point() + facet_wrap(~j, ncol=jmax) +
         geom_vline(aes(xintercept=mean, col="mean"),alpha=0.5, linetype="dashed")
        if(!is.null(pi0)) gg <- gg + geom_vline(aes(xintercept=true_pi, col="true_pi"), linetype="dashed",alpha=0.5) 
        gg <- gg+ scale_color_manual(name = "statistics", values = c(true_pi = "blue", mean = "red"))
        print(gg)
    }

        if (param2plot=="pi") {
        gr <- lapply(1:length(fit$EW_pi), function(k) {
            mea <- fit$EW_pi[k]
            x <- seq(0, 1, length.out=1000)
            data.frame(pi=x, k=k)
            }) %>% bind_rows()
        gr %<>% mutate(density = dbeta(pi, fit$alpha_pi[k], fit$beta_pi[k]),
                        mean = fit$EW_pi[k])
        if(!is.null(pi0)) gr %<>% mutate(true_pi = pi0[k])
        gg <- ggplot(gr, aes(x=pi, y=density)) + geom_line() + facet_wrap(~k, scales="free") +
         geom_vline(aes(xintercept=mean, col="mean"), linetype="dashed")
        if(!is.null(pi0)) gg <- gg + geom_vline(aes(xintercept=true_pi, col="true_pi"), linetype="dashed") 
        gg <- gg+ scale_color_manual(name = "statistics", values = c(true_pi = "blue", mean = "red"))
        print(gg)
    }
    if (param2plot=="gamma") {
        gr <- lapply(1:length(fit$alpha_gamma), function(k) {
            va <- fit$alpha_gamma[k]/fit$beta_gamma[k]^2
            mea <- fit$EW_gamma[k]
            if(is.null(range))  x <- seq(0, max(gamma0[k],5 * mea), , length.out=1000)
            else x <- seq(range[1], range[2], length.out = 1000)
            data.frame(gamma=x, k=k)
            }) %>% bind_rows()
        gr %<>% mutate(density = dgamma(gamma, fit$alpha_gamma[k], fit$beta_gamma[k]),
                        mean = fit$EW_gamma[k])
        if(!is.null(gamma0)) gr %<>% mutate(true_gamma = gamma0[k])
        gg <- ggplot(gr, aes(x=gamma, y=density)) + geom_line() + facet_wrap(~k, scales="free") +
         geom_vline(aes(xintercept=mean, col="mean"), linetype="dashed")
        if(!is.null(gamma0)) gg <- gg + geom_vline(aes(xintercept=true_gamma, col="true_gamma"), linetype="dashed") 
        gg <- gg+ scale_color_manual(name = "statistics", values = c(true_gamma = "blue", mean = "red"))
        print(gg)
    }

    if (param2plot=="tau") {
        va <- fit$alpha_tau/fit$alpha_tau^2
        mea <- fit$EW_tau
        if(is.null(range))  x <- seq(0, max(tau0,5 * mea), , length.out=1000)
        else x <- seq(range[1], range[2], length.out = 1000)
        df <- data.frame(tau = x, 
                         density = dgamma(x, shape = fit$alpha_tau, rate = fit$beta_tau),
                         mean = fit$EW_tau)
        if(!is.null(tau0)) df$true_tau <- tau0
        gg <- ggplot(df, aes(x=tau, y=density)) + geom_line() + geom_vline(aes(xintercept=mean, col="mean"), linetype="dashed")
        if(!is.null(tau0)) gg <- gg + geom_vline(aes(xintercept=true_tau, col="true_tau"), linetype="dashed") 
        gg <- gg+ scale_color_manual(name = "statistics", values = c(true_tau = "blue", mean = "red"))
        print(gg)
    }
}


#' @title plotELBO 
#' @name plotELBO
#' @description plot the ELBO over iteration to monitor convergence of the algorithm 
#' @param fit fit as produced by fit_grpRR
#' @export
plotELBO <- function(fit){
    df <- data.frame(iteration=1:length(fit$ELB_trace),
                     ELBO = fit$ELB_trace)
    ggplot(df, aes(x=iteration, y=ELBO)) +geom_line()
}

