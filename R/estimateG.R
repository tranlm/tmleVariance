###############################################################################
# Description: Add comment
#
# Author: Linh Tran <tranlm@berkeley.edu>
# Date: Jun 17, 2015
###############################################################################


#' @export 
estimateG = function(data, vars.time, gform) {

	## TRANSPOSES ##
	data$id = rownames(data)
	long.data = reshape(data, direction="long", idvar="id", varying=lapply(vars.time, function(x) grep(paste("^",x,sep=""), names(data))), v.names=vars.time)
	long.data = reshape::rename(long.data, c(time="rank"))
	long.data$rank = long.data$rank - 1
	long.data$id = as.numeric(long.data$id)
	long.data = long.data[order(long.data$id, long.data$rank),]
	rownames(long.data) = NULL
	prevY = rep(NA, nrow(long.data))
	prevY[2:length(prevY)] = long.data$Y[1:(nrow(long.data)-1)]
	prevY[long.data$rank==0] = 0
	data = subset(long.data, prevY==0)
	
	## SUBSETS ##
	prevA = rep(NA, nrow(data))
	prevA[2:length(prevA)] = data$A[1:(nrow(data)-1)]
	prevA[data$rank==0] = 0
	A1.data = subset(data, prevA==0)
	
	# Treatments of interest
	A1.data.psi = data

	###########
	## g-fit ##
	###########
	A1.glm = glm(as.formula(gform), data=A1.data, family="binomial")
	gA1.pred = predict(A1.glm, newdata=A1.data.psi, type="response")
	
	
	##############
	## g-matrix ##
	##############
	
	## abar_0 ##
	# n.b. We want P(not enroll), but ltmle package takes P(enroll)
	gA1.abar0 = gA1.pred
	
	## abar_1 ##
	gA1.abar1 = gA1.pred
	gA1.abar1[A1.data.psi$rank>0] = 1
	
	## Combines ##
	gA1 = cbind(A1.data.psi[,c("id", "rank")], gA1.abar0, gA1.abar1)
	merge.data = list(subset(data, select=c('id','rank','Y')), gA1)
	gALL = Reduce(function(x, y) merge(x, y, all=T, by=c("id", "rank")), merge.data, accumulate=F)
	gALL = gALL[order(gALL$id, gALL$rank),]
		
	## Reshapes to wide ##
	gmatrix.abar0 = gALL %>%
			subset(select=c("id", "rank", "gA1.abar0")) %>%
			rename(c(gA1.abar0="A1")) %>%
			reshape(direction="wide", idvar="id", timevar="rank")
	gmatrix.abar1 = gALL %>%
			subset(select=c("id", "rank", "gA1.abar1")) %>%
			rename(c(gA1.abar1="A1")) %>%
			reshape(direction="wide", idvar="id", timevar="rank") 
	rownames(gmatrix.abar0) = gmatrix.abar0$id; rownames(gmatrix.abar1) = gmatrix.abar1$id
	gmatrix.abar0$id = gmatrix.abar1$id = NULL
	gmatrix.abar0 = as.matrix(gmatrix.abar0); gmatrix.abar1 = as.matrix(gmatrix.abar1)
	
	## OUTPUT ##
	out = list(abar0=gmatrix.abar0, abar1=gmatrix.abar1, fit=A1.glm)

	return(out)
	
}

