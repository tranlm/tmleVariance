###############################################################################
# Description: Add comment
#
# Author: Linh Tran <tranlm@berkeley.edu>
# Date: Nov 17, 2015
###############################################################################

#' @export
psi.hat <- function(time.pt, d, i, aiptw=FALSE, tau=.001) {

	QstarFromE <- function(e, off, X, is.deterministic) {
		Qstar <- plogis(off + e*X)
		Qstar[is.deterministic] = rep(1, sum(is.deterministic))
		return(Qstar)
	}
	CalcScore <- function(e, Qstar.kplus1, Q, h.g.ratio, uncensored, intervention.match, is.deterministic) {
		Qstar <- QstarFromE(e, qlogis(Bound(Q, c(tau, 1-tau))), h.g.ratio, is.deterministic)
		ICtemp <- calcIC(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match)
		return(sum(ICtemp)^2)
	}

	data = d[i,]
	# small bug fix
	if (time.pt==1) {
		colnames(data)[which(colnames(data)=="Qbar.a_0")] = "Qbar.a_0.Y.1"
		colnames(data)[which(colnames(data)=="Qbar.a_1")] = "Qbar.a_1.Y.1"
		colnames(data)[which(colnames(data)=="cum.g.a_0")] = "cum.g.a_0.1"
		colnames(data)[which(colnames(data)=="cum.g.a_1")] = "cum.g.a_1.1"
	}

	## Weights ##
	Q.kplus1.a_0 = Q.kplus1.a_1 = data$Q.kplus1
	for(k in time.pt:1) {
		Qstar.data = data.frame(Q.kplus1.a_0=Q.kplus1.a_0, Q.kplus1.a_1=Q.kplus1.a_1, A=data[,paste0("A.",k-1)], a_1=data[,"A.0"], S1=1, off.a_0=qlogis(Bound(data[,paste0("Qbar.a_0.Y.",k)], c(tau,1-tau))), weight.a_0=1/data[,paste0("cum.g.a_0.",k)], off.a_1=qlogis(Bound(data[,paste0("Qbar.a_1.Y.",k)], c(tau,1-tau))), weight.a_1=1/data[,paste0("cum.g.a_1.",k)])
		Qstar.data$H.a_0 = ifelse(Qstar.data$A==0, 1, 0)
		Qstar.data$H.a_1 = ifelse(Qstar.data$a_1==1, 1, 0)
		Qstar.a_0 = suppressWarnings(glm(Q.kplus1.a_0 ~ -1 + offset(off.a_0) + H.a_0, data=Qstar.data, weights=weight.a_0, family="quasibinomial", subset=(data[,paste0("Qbar.a_0.Y.",k)]!=1)))
		Qstar.a_1 = suppressWarnings(glm(Q.kplus1.a_1 ~ -1 + offset(off.a_1) + H.a_1, data=Qstar.data, weights=weight.a_1, family="quasibinomial", subset=(data[,paste0("Qbar.a_1.Y.",k)]!=1)))
		score.a_0 = CalcScore(Qstar.a_0$coefficients[[1]], Q.kplus1.a_0, data[,paste0("Qbar.a_0.Y.",k)], rep(1,nrow(data)), rep(TRUE,nrow(Qstar.data)), Qstar.data$IndTx.a_0==1, data[,paste0("Qbar.a_0.Y.",k)]==1)
		if (score.a_0>0.001) {
			m <- nlminb(start=0, objective=CalcScore, Qstar.kplus1=Q.kplus1.a_0, Q=data[,paste0("Qbar.a_0.Y.",k)], h.g.ratio=rep(1,nrow(data)), uncensored=rep(TRUE,nrow(Qstar.data)), intervention.match=Qstar.data$IndTx.a_0==1, is.deterministic=data[,paste0("Qbar.a_0.Y.",k)]==1, control=list(abs.tol=0.0001^2, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
			if (m$objective> 0.0001^2) warning("Minimizer failed.")
			Qstar.a_0$coefficients[[1]] <- m$par
		}
		score.a_1 = CalcScore(Qstar.a_1$coefficients[[1]], Q.kplus1.a_1, data[,paste0("Qbar.a_1.Y.",k)], rep(1,nrow(data)), rep(TRUE,nrow(Qstar.data)), Qstar.data$IndTx.a_1==1, data[,paste0("Qbar.a_1.Y.",k)]==1)
		if (score.a_1>0.001) {
			m <- nlminb(start=0, objective=CalcScore, Qstar.kplus1=Q.kplus1.a_1, Q=data[,paste0("Qbar.a_1.Y.",k)], h.g.ratio=rep(1,nrow(data)), uncensored=rep(TRUE,nrow(Qstar.data)), intervention.match=Qstar.data$IndTx.a_1==1, is.deterministic=data[,paste0("Qbar.a_1.Y.",k)]==1, control=list(abs.tol=0.0001^2, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
			if (m$objective> 0.0001^2) warning("Minimizer failed.")
			Qstar.a_1$coefficients[[1]] <- m$par
		}
		if (score.a_0>0.001 | score.a_1>0.001) return(-Inf)
		Q.kplus1.a_0 = plogis(Qstar.data$off.a_0 + Qstar.a_0$coef[[1]])
		Q.kplus1.a_1 = plogis(Qstar.data$off.a_1 + Qstar.a_1$coef[[1]])
		Q.kplus1.a_0[data[,paste0("Qbar.a_0.Y.",k)]==1] = Q.kplus1.a_1[data[,paste0("Qbar.a_1.Y.",k)]==1] = 1
	}
	out.wt = c(tmle.wt.psi.a_0=mean(Q.kplus1.a_0), tmle.wt.psi.a_1=mean(Q.kplus1.a_1))

	## Covariate ##
	Q.kplus1.a_0 = Q.kplus1.a_1 = data$Q.kplus1
	for(k in time.pt:1) {
		Qstar.data = data.frame(Q.kplus1.a_0=Q.kplus1.a_0, Q.kplus1.a_1=Q.kplus1.a_1, A=data[,paste0("A.",k-1)], a_1=data[,"A.0"], S1=1, off.a_0=qlogis(Bound(data[,paste0("Qbar.a_0.Y.",k)], c(tau,1-tau))), weight.a_0=1/data[,paste0("cum.g.a_0.",k)], off.a_1=qlogis(Bound(data[,paste0("Qbar.a_1.Y.",k)], c(tau,1-tau))), weight.a_1=1/data[,paste0("cum.g.a_1.",k)])
		Qstar.data$IndTx.a_0 = ifelse(Qstar.data$A==0,1,0)
		Qstar.data$IndTx.a_1 = ifelse(Qstar.data$a_1==1,1,0)
		Qstar.data$H.a_0 = Qstar.data$IndTx.a_0*Qstar.data$weight.a_0
		Qstar.data$H.a_1 = Qstar.data$IndTx.a_1*Qstar.data$weight.a_1
		Qstar.a_0 = suppressWarnings(glm(Q.kplus1.a_0 ~ -1 + offset(off.a_0) + H.a_0, data=Qstar.data[data[,paste0("Qbar.a_0.Y.",k)]!=1,], family="quasibinomial"))
		Qstar.a_1 = suppressWarnings(glm(Q.kplus1.a_1 ~ -1 + offset(off.a_1) + H.a_1, data=Qstar.data[data[,paste0("Qbar.a_1.Y.",k)]!=1,], family="quasibinomial"))
		score.a_0 = CalcScore(Qstar.a_0$coefficients[[1]], Q.kplus1.a_0, data[,paste0("Qbar.a_0.Y.",k)], Qstar.data$weight.a_0, rep(TRUE,nrow(Qstar.data)), Qstar.data$IndTx.a_0==1, data[,paste0("Qbar.a_0.Y.",k)]==1)
		if (score.a_0>0.001) {
			m <- nlminb(start=0, objective=CalcScore, Qstar.kplus1=Q.kplus1.a_0, Q=data[,paste0("Qbar.a_0.Y.",k)], h.g.ratio=Qstar.data$weight.a_0, uncensored=rep(TRUE,nrow(Qstar.data)), intervention.match=Qstar.data$IndTx.a_0==1, is.deterministic=data[,paste0("Qbar.a_0.Y.",k)]==1, control=list(abs.tol=0.0001^2, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
			if (m$objective> 0.0001^2) warning("Minimizer failed.")
			Qstar.a_0$coefficients[[1]] <- m$par
		}
		score.a_1 = CalcScore(Qstar.a_1$coefficients[[1]], Q.kplus1.a_1, data[,paste0("Qbar.a_1.Y.",k)], Qstar.data$weight.a_1, rep(TRUE,nrow(Qstar.data)), Qstar.data$IndTx.a_1==1, data[,paste0("Qbar.a_1.Y.",k)]==1)
		if (score.a_1>0.001) {
			m <- nlminb(start=0, objective=CalcScore, Qstar.kplus1=Q.kplus1.a_1, Q=data[,paste0("Qbar.a_1.Y.",k)], h.g.ratio=Qstar.data$weight.a_1, uncensored=rep(TRUE,nrow(Qstar.data)), intervention.match=Qstar.data$IndTx.a_1==1, is.deterministic=data[,paste0("Qbar.a_1.Y.",k)]==1, control=list(abs.tol=0.0001^2, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
			if (m$objective> 0.0001^2) warning("Minimizer failed.")
			Qstar.a_1$coefficients[[1]] <- m$par
		}
		Q.kplus1.a_0 = QstarFromE(Qstar.a_0$coef[[1]], Qstar.data$off.a_0, Qstar.data$weight.a_0, data[,paste0("Qbar.a_0.Y.",k)]==1)
		Q.kplus1.a_1 = QstarFromE(Qstar.a_1$coef[[1]], Qstar.data$off.a_1, Qstar.data$weight.a_1, data[,paste0("Qbar.a_1.Y.",k)]==1)
	}
	out.cov = c(tmle.cov.psi.a_0=mean(Q.kplus1.a_0), tmle.cov.psi.a_1=mean(Q.kplus1.a_1))

	out = c(out.wt, out.cov)
	if (aiptw & time.pt==1) {
		aiptw.psi_0 = mean(ifelse(data[,"A.0"]==0,1,0)/data[,"cum.g.a_0.1"]*(data[,"Q.kplus1"]-data[,"Qbar.a_0.Y.1"]) + data[,"Qbar.a_0.Y.1"])
		aiptw.psi_1 = mean(ifelse(data[,"A.0"]==1,1,0)/data[,"cum.g.a_1.1"]*(data[,"Q.kplus1"]-data[,"Qbar.a_1.Y.1"]) + data[,"Qbar.a_1.Y.1"])
		out = c(out, aiptw.psi_0=aiptw.psi_0, aiptw.psi_1=aiptw.psi_1)
	} else {
		out = c(out, aiptw.psi_0=NA, aiptw.psi_1=NA)
	}
	return(out)

}

