###############################################################################
# Description: Add comment
#
# Author: Linh Tran <tranlm@berkeley.edu>
# Date: Jan 11, 2016
###############################################################################


#' @export
appliedResults = function(tmp, gmatrix.abar0, gmatrix.abar1, SL.library, time.pt, bootIndicies) {
	
	####################
	## ARGS for LTMLE ##
	####################
	Lnodes = Cnodes = gform.uadj = Qform.uadj = Qform = gform = cum.g.name = NULL
	for(i in 1:time.pt){
		## NODES ##
		Lnodes = c(Lnodes, paste0(vars.time,".", i))
		Cnodes = c(Cnodes, paste0("transfer.",i-1), paste0("eofu.",i-1))
		cum.g.name = c(cum.g.name, paste0("enroll.",i-1), paste0("transfer.",i-1), paste0("eofu.",i-1))
		
		## UNADJUSTED MODELS ##
		gform.uadj = c(gform.uadj, paste("enroll.", i-1, " ~ 1", sep=""))
		gform.uadj = c(gform.uadj, paste("transfer.", i-1, " ~ 1", sep=""))
		gform.uadj = c(gform.uadj, paste("eofu.", i-1, " ~ 1", sep=""))
		Qform.uadj = c(Qform.uadj, paste("Q.kplus1 ~ 1 ", sep=""))
		
		gform = c(gform, paste("enroll.", i-1, " ~ ", paste(c(vars.base, paste0(c(vars.time),".",i-1)), collapse=" + "), sep=""))
		gform = c(gform, paste("transfer.", i-1, " ~ ", paste(c(vars.base, paste0(c(vars.time, 'enroll'),".",i-1)), collapse=" + "), sep=""))
		gform = c(gform, paste("eofu.", i-1, " ~ ", paste(c(vars.base, paste0(c(vars.time),".",i-1)), collapse=" + "), sep=""))
		Qform = c(Qform, paste("Q.kplus1 ~ ", paste(c(vars.base, paste0(c(vars.time, "enroll"),".",i-1)), collapse=" + ")))
		names(Qform.uadj)[length(Qform.uadj)] = names(Qform)[length(Qform)] = paste("arvadhere.", i, sep="")
	}
	Ynodes = grep("^dead.ltfu.", names(tmp)) 
	Anodes = grep("^enroll.", names(tmp))
#	Qform = NULL	

	##################
	## INITIAL Qbar ##
	##################
	ltmle.ate = suppressMessages(ltmle(data=tmp, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE, Qform=Qform, gform=gform, abar=list(rep(0,time.pt), rep(1,time.pt)), gbounds=c(0.001,1), deterministic.g.function=MaintainTreatment, estimate.time=FALSE, stratify=FALSE, SL.library=SL.library, IC.variance.only=FALSE))
	saveRDS(ltmle.ate, file=sprintf('./inst/results/applied-ltmle-t%d.RDS', time.pt))
	gcomp.fit_0 = ltmle(data=tmp, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE, Qform=Qform, gform=gmatrix.abar0[,1:length(c(Anodes, Cnodes))], abar=rep(0,time.pt), gbounds=c(0.001,1), SL.library=SL.library, estimate.time=FALSE, stratify=FALSE, IC.variance.only=TRUE, gcomp=TRUE)
	gcomp.fit_1 = ltmle(data=tmp, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE, Qform=Qform, gform=gmatrix.abar1[,1:length(c(Anodes, Cnodes))], abar=rep(1,time.pt), gbounds=c(0.001,1), SL.library=SL.library, estimate.time=FALSE, stratify=FALSE, IC.variance.only=TRUE, gcomp=TRUE)
	Qbar.matrix.a_0 = gcomp.fit_0$fit$Qbar.matrix
	Qbar.matrix.a_1 = gcomp.fit_1$fit$Qbar.matrix
	Qbar.matrix.a_0[is.na(Qbar.matrix.a_0)] = Qbar.matrix.a_1[is.na(Qbar.matrix.a_1)] = 1
	Qbar.matrix = list(Qbar.matrix.a_0=Qbar.matrix.a_0, Qbar.matrix.a_1=Qbar.matrix.a_1)
	saveRDS(Qbar.matrix, file=sprintf('./inst/results/applied-Qbar-t%d.RDS', time.pt))

	
	###############
	## BOOTSTRAP ##
	###############
	psi.hat <- function(time.pt, d, i) {
		
		data = d[i,]
		
		## Weights ##
		Q.kplus1.a_0 = Q.kplus1.a_1 = data$Q.kplus1
		for(k in time.pt:1) {
			Qstar.data = data.frame(Q.kplus1.a_0=Q.kplus1.a_0, Q.kplus1.a_1=Q.kplus1.a_1, A=data[,paste0("enroll.",k-1)], C=data[,paste0("transfer.",k-1)]=="uncensored" & data[,paste0("eofu.",k-1)]=="uncensored", a_1=data[,"enroll.0"], S1=1, off.a_0=qlogis(Bound(data[,paste0("Qbar.a_0.dead.ltfu.",k)], c(.0001,.9999))), weight.a_0=1/data[,paste0("cum.g.a_0.",k)], off.a_1=qlogis(Bound(data[,paste0("Qbar.a_1.dead.ltfu.",k)], c(.0001,.9999))), weight.a_1=1/data[,paste0("cum.g.a_1.",k)])
			Qstar.a_0 = glm(Q.kplus1.a_0 ~ -1 + offset(off.a_0) + S1, data=Qstar.data, weights=scale(weight.a_0, center=FALSE), family="quasibinomial", subset=(C & A==0 & data[,paste0("Qbar.a_0.dead.ltfu.",k)]!=1))
			Qstar.a_1 = glm(Q.kplus1.a_1 ~ -1 + offset(off.a_1) + S1, data=Qstar.data, weights=scale(weight.a_1, center=FALSE), family="quasibinomial", subset=(C & a_1==1 & data[,paste0("Qbar.a_1.dead.ltfu.",k)]!=1))
			Q.kplus1.a_0 = plogis(Qstar.data$off.a_0 + Qstar.a_0$coef[[1]])
			Q.kplus1.a_1 = plogis(Qstar.data$off.a_1 + Qstar.a_1$coef[[1]])
			Q.kplus1.a_0[data[,paste0("Qbar.a_0.dead.ltfu.",k)]==1] = Q.kplus1.a_1[data[,paste0("Qbar.a_1.dead.ltfu.",k)]==1] = 1
		}
		out.wt = c(tmle.wt.psi.a_0=mean(Q.kplus1.a_0), tmle.wt.psi.a_1=mean(Q.kplus1.a_1))
		
		## Covariate ##
		Q.kplus1.a_0 = Q.kplus1.a_1 = data$Q.kplus1
		for(k in time.pt:1) {
			Qstar.data = data.frame(Q.kplus1.a_0=Q.kplus1.a_0, Q.kplus1.a_1=Q.kplus1.a_1, A=data[,paste0("enroll.",k-1)], a_1=data[,"enroll.0"], S1=1, off.a_0=qlogis(Bound(data[,paste0("Qbar.a_0.dead.ltfu.",k)], c(.0001,.9999))), weight.a_0=1/data[,paste0("cum.g.a_0.",k)], off.a_1=qlogis(Bound(data[,paste0("Qbar.a_1.dead.ltfu.",k)], c(.0001,.9999))), weight.a_1=1/data[,paste0("cum.g.a_1.",k)])
			Qstar.data$H.a_1 = ifelse(Qstar.data$a_1==1,1,0)*Qstar.data$weight.a_1
			Qstar.data$H.a_0 = ifelse(Qstar.data$A==0,1,0)*Qstar.data$weight.a_0
			Qstar.a_0 = glm(Q.kplus1.a_0 ~ -1 + offset(off.a_0) + H.a_0, data=Qstar.data, family="quasibinomial", subset=(data[,paste0("Qbar.a_0.dead.ltfu.",k)]!=1))
			Qstar.a_1 = glm(Q.kplus1.a_1 ~ -1 + offset(off.a_1) + H.a_1, data=Qstar.data, family="quasibinomial", subset=(data[,paste0("Qbar.a_1.dead.ltfu.",k)]!=1))
			Q.kplus1.a_0 = plogis(Qstar.data$off.a_0 + Qstar.a_0$coef[[1]]*Qstar.data$weight.a_0)
			Q.kplus1.a_1 = plogis(Qstar.data$off.a_1 + Qstar.a_1$coef[[1]]*Qstar.data$weight.a_1)
			Q.kplus1.a_0[data[,paste0("Qbar.a_0.dead.ltfu.",k)]==1] = Q.kplus1.a_1[data[,paste0("Qbar.a_1.dead.ltfu.",k)]==1] = 1
		}
		out.cov = c(tmle.cov.psi.a_0=mean(Q.kplus1.a_0), tmle.cov.psi.a_1=mean(Q.kplus1.a_1))
		
		out = c(out.wt, out.cov)
		return(out)
		
	}
	boot.data = cbind(Q.kplus1=tmp[[Ynodes[length(Ynodes)]]], tmp[,c(Anodes),drop=F], tmp[,c(Cnodes),drop=F], cum.g.a_0=ltmle.ate$cum.g[,seq(1, 3*time.pt) %% 3 == 0, 1, drop=FALSE], Qbar.a_0=Qbar.matrix.a_0, cum.g.a_1=ltmle.ate$cum.g[,seq(1, 3*time.pt) %% 3 == 0, 1, drop=FALSE], Qbar.a_1=Qbar.matrix.a_1)
	if(time.pt==1) names(boot.data)[c(5:8)] = c("cum.g.a_0.1", "Qbar.a_0.dead.ltfu.1", "cum.g.a_1.1", "Qbar.a_1.dead.ltfu.1")
	bootEsts = do.call("rbind", lapply(1:ncol(bootIndicies), function(i) psi.hat(time.pt, boot.data, bootIndicies[,i])))
	#bootVar = apply(bootEsts, 2, var)*n
	se.boot.wt = sd(bootEsts[,'tmle.wt.psi.a_1'] - bootEsts[,'tmle.wt.psi.a_0'])
	se.boot.cov = sd(bootEsts[,'tmle.cov.psi.a_1'] - bootEsts[,'tmle.cov.psi.a_0'])
	
	#############
	## RESULTS ##
	#############
	fitSummary = summary(ltmle.ate)
	output = list(
			psi.tmle = fitSummary$effect.measures$ATE$estimate,
			se.ic = fitSummary$measures.IC$ATE$std.dev,
			se.robust = fitSummary$measures.variance.estimate$ATE$std.dev,
			se.boot.wt = se.boot.wt,
			se.boot.cov = se.boot.cov,
			p.ic = 2*pnorm(-abs(fitSummary$effect.measures$ATE$estimate / fitSummary$measures.IC$ATE$std.dev)),
			p.robust = 2*pnorm(-abs(fitSummary$effect.measures$ATE$estimate / fitSummary$measures.variance.estimate$ATE$std.dev)),
			p.boot.wt = 2*pnorm(-abs(fitSummary$effect.measures$ATE$estimate / se.boot.wt)),
			p.boot.cov = 2*pnorm(-abs(fitSummary$effect.measures$ATE$estimate / se.boot.cov))
	)
	saveRDS(output, sprintf('./inst/results/applied-output-t%d.RDS', time.pt))
	return(output)
	
}


