###############################################################################
# Description: Runs a single simulation and returns results
#
# Author: Linh Tran <tranlm@berkeley.edu>
# Date: Dec 25, 2016
###############################################################################

RunSimulation = function(n, time.pt, beta_pos, beta_ate, gform, Qform, delta_0, delta_b, bootIndicies, seed) {
	
	if (!missing(seed)) {
		set.seed(seed)
	}
	
	# P_n #
	O = generateData(n=n, time.pt=time.pt, beta_pos=beta_pos, beta_ate=beta_ate)$O
	
	# g(A|Pa) #
	#nb. I pooled, while the package stratifies by time
	gform.t = NULL
	for(i in 1:time.pt) {
		## TRUE MODELS ##
		gform.t = c(gform.t, paste("A.", i-1, " ~ W1 + W2", paste(" + L1.", i-1, sep=""), paste(" + L2.", i-1, sep=""), paste(" + L1.", i-1, ":L2.", i-1, sep=""), sep=""))
	}
	if (time.pt>1) {
		O.long = reshape(O, direction="long", varying=list(grep("^L1.", names(O)), grep("^L2.", names(O)), grep("^A.", names(O)), grep("^Y.", names(O))), v.names=c("L1", "L2", "A", "Y"))
		O.long = O.long[order(O.long$id, O.long$time),]
		nodes = list(A=which(names(O) %in% Anodes), C=NULL, L=which(names(O) %in% Lnodes), Y=which(names(O) %in% Ynodes), AC=which(names(O) %in% Anodes), LY=which(names(O) %in% Ynodes))
		prevA = prevY = rep(NA,nrow(O.long))
		prevA[2:nrow(O.long)] = O.long$A[1:(nrow(O.long)-1)]; prevA[O.long$time==1] = 0 
		prevY[2:nrow(O.long)] = O.long$Y[1:(nrow(O.long)-1)]; prevY[O.long$time==1] = 0
		g.data = subset(O.long, prevA==0 & prevY==0)
		g.fit.true = glm(as.formula(gform), data=g.data, family="binomial")
		O.long$predA_1 = predict(g.fit.true, newdata=O.long, type="response")
		O.long$predA_1[prevY==1] = ifelse(prevA[prevY==1]==1,1,0)
		gmatrix.a_0 = reshape(subset(O.long, select=c("id", "time", "predA_1")), direction="wide", idvar="id")
		gmatrix.a_0 = as.matrix(gmatrix.a_0[,grep("^predA_1", colnames(gmatrix.a_0))])
		gmatrix.a_0[is.na(gmatrix.a_0)] = 0
		colnames(gmatrix.a_0) = paste0("A.", c(1:time.pt)-1)
		gmatrix.a_1 = gmatrix.a_0; gmatrix.a_1[,c(2:time.pt)] = 1
	} else {
		g.fit.true = glm(as.formula(gform), data=O, family="binomial")
		gmatrix.a_1 = matrix(g.fit.true$fitted.values, ncol=1, dimnames=list (NULL, "A.1"))
		gmatrix.a_0 = 1 - gmatrix.a_1; colnames(gmatrix.a_0) = "A.0 "
	}

	# Q(Y|Pa) #
	#nb. For convenience, just used the ltmle package
	#    https://github.com/joshuaschwab/ltmle/tree/0982a33e34b6455ce7ad470facdaf45d012f9a9f
	#    I had to manually update the package to output the Qstar arrays
	Qmatrix.a_0 = ltmle(data=O, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE, Qform=Qform, gform=gmatrix.a_1, abar=rep(0,time.pt), gbounds=c(delta_0,1), SL.library=NULL, estimate.time=FALSE, stratify=FALSE, gcomp=TRUE)$Qstar.array
	Qmatrix.a_1 = ltmle(data=O, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE, Qform=Qform, gform=gmatrix.a_1, abar=rep(1,time.pt), gbounds=c(delta_0,1), SL.library=NULL, estimate.time=FALSE, stratify=FALSE, gcomp=TRUE)$Qstar.array
	if (!is.matrix(Qmatrix.a_0)) {
		Qmatrix.a_0 = matrix(Qmatrix.a_0, ncol=time.pt)
		Qmatrix.a_1 = matrix(Qmatrix.a_1, ncol=time.pt)
	}
	colnames(Qmatrix.a_0) = colnames(Qmatrix.a_1) = Ynodes
	
	# EIF variance estimate #
	ltmle.ic.fit_0 = ltmle(data=O, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE, Qform=Qform, gform=gmatrix.a_1, abar=rep(0,time.pt), gbounds=c(delta_0,1), SL.library=NULL, estimate.time=FALSE, stratify=FALSE, variance.method="ic")
	ltmle.ic.fit_1 = ltmle(data=O, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE, Qform=Qform, gform=gmatrix.a_1, abar=rep(1,time.pt), gbounds=c(delta_0,1), SL.library=NULL, estimate.time=FALSE, stratify=FALSE, variance.method="ic")

	# Robust variance estimates #
	ltmle.fit_0 = ltmle(data=O, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE, Qform=Qform, gform=gform.t, abar=rep(0,time.pt), gbounds=c(delta_0,1), SL.library=NULL, estimate.time=FALSE, stratify=FALSE, variance.method="tmle")
	ltmle.fit_1 = ltmle(data=O, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE, Qform=Qform, gform=gform.t, abar=rep(1,time.pt), gbounds=c(delta_0,1), SL.library=NULL, estimate.time=FALSE, stratify=FALSE, variance.method="tmle")

	# Bootstrap variance estimate #
	# Update (20220530: Newer version of R doesn't automatically prepend names
	if (time.pt==1){
		colnames(Qmatrix.a_0) = paste0("Qbar.a_0.", colnames(Qmatrix.a_0))
		colnames(Qmatrix.a_1) = paste0("Qbar.a_1.", colnames(Qmatrix.a_1))		
	}
	boot.data = cbind(Q.kplus1=O[,paste0("Y.", time.pt)], O[,Anodes,drop=F], cum.g.a_0=apply(ltmle.ic.fit_0$cum.g.unbounded, 2, function(x) pmax(x, delta_b)),
	                  Qbar.a_0=Qmatrix.a_0, cum.g.a_1=apply(ltmle.ic.fit_1$cum.g.unbounded, 2, function(x) pmax(x, delta_b)), Qbar.a_1=data.frame(Qmatrix.a_1))
	bootEsts = do.call("rbind", lapply(1:ncol(bootIndicies), function(i) psi.hat(time.pt, boot.data, bootIndicies[,i], aiptw=TRUE)))
	bootVar = apply(bootEsts, 2, var)
	
	# aiptw #
	aiptw.fit_0 = aiptw(data=O, Ynodes=Ynodes, Anodes=Anodes, Cnodes=NULL, abar=rep(0,time.pt), cum.g=ltmle.ic.fit_0$cum.g, Qform=Qmatrix.a_0, SL.library=NULL, stratify=FALSE)
	aiptw.fit_1 = aiptw(data=O, Ynodes=Ynodes, Anodes=Anodes, Cnodes=NULL, abar=rep(1,time.pt), cum.g=ltmle.ic.fit_1$cum.g, Qform=Qmatrix.a_1, SL.library=NULL, stratify=FALSE)

	est.a_0 = c(psi.aiptw.a_0=aiptw.fit_0$estimate, 
			psi.tmle.a_0=ltmle.ic.fit_0$estimates[["tmle"]],
			variance.aiptw.ic.a_0=var(aiptw.fit_0$IC)/n, 
			variance.bootAIPTW.a_0=bootVar[["aiptw.psi_0"]],
			variance.tmle.ic.a_0=var(ltmle.ic.fit_0$IC$tmle)/n,
			variance.robust.a_0=ltmle.fit_0$variance.estimate/n, 
			variance.bootTMLE.wt.a_0=bootVar[["tmle.wt.psi.a_0"]], 
			variance.bootTMLE.cov.a_0=bootVar[["tmle.cov.psi.a_0"]])
	est.a_1 = c(psi.aiptw.a_1=aiptw.fit_1$estimate,
			psi.tmle.a_1=ltmle.ic.fit_1$estimates[["tmle"]], 
			variance.aiptw.ic.a_1=var(aiptw.fit_1$IC)/n,
			variance.bootAIPTW.a_1=bootVar[["aiptw.psi_1"]], 
			variance.tmle.ic.a_1=var(ltmle.ic.fit_1$IC$tmle)/n,
			variance.robust.a_1=ltmle.fit_1$variance.estimate/n, 
			variance.bootTMLE.wt.a_1=bootVar[["tmle.wt.psi.a_1"]], 
			variance.bootTMLE.cov.a_1=bootVar[["tmle.cov.psi.a_1"]])
	
	out = rbind(est.a_0, est.a_1)
	colnames(out) = c("psi.aiptw", "psi.tmle", "sigma2.aiptw.ic", "sigma2.aiptw.boot", "sigma2.tmle.ic", "sigma2.tmle.robust", "sigma2.tmle.boot_wt", "sigma2.tmle.boot_cov")
	return(out)
}


