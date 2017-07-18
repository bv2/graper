# plotting_functions
#
# This file contains functions used for plotting diagnositcs for the grpRR model fit
#         - plotPosterior: Plot posterior ditribution on parameters beta, gamma and tau
#         - plotVBFit: Plot convergence of parameters (parameter values vs. iteration)
#

#' @import ggplot2
#' @import gridExtra
#' @import reshape2
#' @export

# ---------------------------
plotPosterior<-function(result, beta0=NULL, gamma0=NULL, tau0=NULL, whichParam=c("tau", "beta", "gamma")){

  if("beta" %in% whichParam){
    par(mfrow=c(3,4))
    if( !is.null(dim(result$Sigma_beta))) {
      if(ncol(result$Sigma_beta)==1 |nrow(result$Sigma_beta)==1) result$Sigma_beta<-as.numeric(result$Sigma_beta)
    }
    if(is.null(dim(result$Sigma_beta))) result$Sigma_beta<-diag(result$Sigma_beta) #transfrom to diagnoal matrix in fully fatorized case
    for(i in 1:length(result$EW_beta)) {
      sx<-seq(min(-10, beta0[i]),max(10, beta0[i]))
      plot(x=sx, dnorm(sx,result$EW_beta[i], sd=sqrt(result$Sigma_beta[i,i])),'l', main=paste("beta", i))
      abline(v=result$EW_beta[i], col="blue")
      if(!is.null(beta0)) abline(v=beta0[i], col="red")
    }
  }

  if("gamma" %in% whichParam){
    par(mfrow=c(3,4))
    for(k in 1:length(result$alpha_gamma)) {
      sx<-seq(min(0, gamma0[k]),max(100,gamma0[k]))
      plot(sx, dgamma(sx,shape=result$alpha_gamma[k], rate=result$beta_gamma[k]),'l', main=paste("gamma", k))
      abline(v=result$EW_gamma[k], col="blue")
      if(!is.null(gamma0)) abline(v=gamma0[k], col="red")

    }
  }

  if("tau" %in% whichParam){
    par(mfrow=c(3,4))
    sx<-seq(min(0, tau0),max(100, tau0))
    plot(x=sx, dgamma(sx,shape=result$alpha_tau, rate=result$beta_tau),'l', main="tau")
    abline(v=result$EW_tau[k], col="blue")
    if(!is.null(tau0)) abline(v=tau0, col="red")
  }
}



# ---------------------------
plotVBFit<-function(result, whichParam=c("ELB","tau", "beta", "gamma")){
  #parameter convergence
  if("ELB" %in% whichParam) gg_ELB<-ggplot(result$trace_param, aes(x=iter, y=ELB))+geom_line()+ggtitle("Evidence lower bound")
  if("tau" %in% whichParam) gg_tau<-ggplot(result$trace_param, aes(x=iter, y=tau))+geom_line()+ggtitle("Noise precision")
  if("ELB" %in% whichParam &"ELB" %in% whichParam) grid.arrange(gg_ELB, gg_tau, nrow=2) else if("ELB" %in% whichParam) print(gg_ELB) else if("tau" %in% whichParam) print(gg_tau)

  if("beta" %in% whichParam )  {
    tracebeta<-melt(result$trace_param, id.vars = "iter", measure.vars = colnames(result$trace_param)[grep("beta",colnames(result$trace_param))])
    gg_beta<-ggplot(tracebeta, aes(x=iter, y=value))+geom_line()+facet_wrap(~variable)+ggtitle("Linear model coefficients")
    print(gg_beta)
  }

  if("gamma" %in% whichParam )  {
    tracegamma<-melt(result$trace_param, id.vars = "iter", measure.vars = colnames(result$trace_param)[grep("gamma",colnames(result$trace_param))])
    gg_gamma<-ggplot(tracegamma, aes(x=iter, y=value))+geom_line(size=1.5)+facet_wrap(~variable, ncol=4)+ggtitle("Prior precision per group")+
      theme(text = element_text(size=20),axis.text.x = element_text(angle=60, hjust=1))+xlab("iteration")+ylab("")
    print(gg_gamma)
  }
}

