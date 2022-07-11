###############################################################################
# Description: Add comment
#
# Author: Linh Tran <tranlm@berkeley.edu>
# Date: Aug 27, 2015
###############################################################################

#' @export 
SL.gbm.LT = function (Y, X, newX, family, obsWeights, gbm.trees = 500, interaction.depth = 1, cv.folds = 2, ...) {
	require(gbm)
	gbm.model = as.formula(paste("Y ~", paste(colnames(X), collapse = "+")))
	if (family$family == "binomial" & !all(Y %in% c(0,1))) {
		fit.gbm = gbm(formula = gbm.model, data = X, distribution = "gaussian", 
				n.trees = gbm.trees, interaction.depth = interaction.depth, 
				cv.folds = cv.folds, keep.data = TRUE, weights = obsWeights, 
				verbose = FALSE)
		best.iter = gbm.perf(fit.gbm, method = "cv", plot.it = FALSE)
		pred = bound(predict(fit.gbm, newdata = newX, best.iter, type = "response"), c(0,1))
		fit = list(object = fit.gbm, n.trees = best.iter)
	}
	else if (family$family == "binomial") {
		newY = Y
		fit.gbm = gbm(formula = gbm.model, data = X, distribution = "bernoulli", 
				n.trees = gbm.trees, interaction.depth = interaction.depth, 
				cv.folds = cv.folds, keep.data = TRUE, verbose = FALSE, 
				weights = obsWeights)
		best.iter = gbm.perf(fit.gbm, method = "cv", plot.it = FALSE)
		pred = predict(fit.gbm, newdata = newX, best.iter, type = "response")
		fit = list(object = fit.gbm, n.trees = best.iter)
	}
	out = list(pred = pred, fit = fit)
	class(out$fit) = c("SL.gbm")
	return(out)
}

#' @export
SL.svm.LT = function (Y, X, newX, family, type.reg = "eps-regression", type.class = "C-classification", nu = 0.5, gamma = 0.1, scale=FALSE,...) {
	library("e1071")
	if (family$family == "binomial" & !all(Y %in% c(0,1))) {
		fit.svm = try(svm(y = Y, x = X, nu = nu, type = type.reg, fitted = FALSE, gamma=gamma, scale=scale), silent=TRUE)
		if(inherits(fit.svm, "try-error")) {
			pred = rep(mean(Y), nrow(newX))
			fit = list(object=fit.svm[1])
		} else {
			pred = Bound(predict(fit.svm, newdata = newX), c(0,1))
			fit = list(object = fit.svm)
		}
	}
	else if (family$family == "binomial") {
		newY = as.factor(Y)
		fit.svm = try(svm(y = newY, x = X, nu = nu, type = type.class, fitted = FALSE, gamma=gamma, probability = TRUE, scale=scale), silent=TRUE)
		if(inherits(fit.svm, "try-error")) {
			pred = rep(mean(Y), nrow(newX))
			fit = list(object=fit.svm[1])
		} else {
			pred = try(attr(predict(fit.svm, newdata = newX, probability = TRUE), "prob")[, "1"])
			if(inherits(pred, "try-error")) {
				pred = rep(mean(Y), nrow(newX))
			}
			fit = list(object = fit.svm)
		}
	}
	out = list(pred = pred, fit = fit)
	class(out$fit) = c("SL.svm")
	return(out)
}

#' @export
SL.polymars.LT = function(Y, X, newX, family, obsWeights, cv=2, seed=9999, ...){
	library("polspline")
	if (family$family == "binomial" & !all(Y %in% c(0,1))) {
		fit.mars = try(polymars(Y, X, weights = obsWeights), silent=TRUE)
		if(inherits(fit.mars, "try-error")) {
			pred = rep(mean(Y), nrow(newX))
			fit = list(object=fit.mars[1])
		} else {
			pred = Bound(predict(fit.mars, x = newX), c(0,1))
			fit = list(object = fit.mars)
		}
	}
	else if (family$family == "binomial") {
		newY = Y
		fit.mars = try(polyclass(newY, X, cv = cv, weight = obsWeights, seed=seed), silent=TRUE)
		if(inherits(fit.mars, "try-error")) {
			pred = rep(mean(Y), nrow(newX))
			fit = list(object=fit.mars[1])
		} else {
			pred = ppolyclass(cov = newX, fit = fit.mars)[, 2]
			fit = list(fit = fit.mars)
		}
	}
	out = list(pred = pred, fit = fit)
	class(out$fit) = c("SL.polymars")
	return(out)
}

#' @export
SL.nnet.LT = function (Y, X, newX, family, obsWeights, size = 2, maxit = 1000, ...) {
	library("nnet")
	if (family$family == "binomial" & !all(Y %in% c(0,1))) {
		fit.nnet = try(nnet(x = X, y = Y, size = size, trace = FALSE, maxit = maxit, linout = TRUE, weights = obsWeights), silent=TRUE)
		if(inherits(fit.nnet, "try-error")) {
			pred = rep(mean(Y), nrow(newX))
			fit = list(object=fit.nnet[1])
		} else {
			pred = Bound(predict(fit.nnet, newdata = newX, type = "raw"), c(0,1))
			fit = list(object = fit.nnet)
		}
	}
	else if (family$family == "binomial") {
		newY = Y
		fit.nnet = try(nnet(x = X, y = newY, size = size, trace = FALSE, maxit = maxit, linout = FALSE, weights = obsWeights), silent=TRUE)
		if(inherits(fit.nnet, "try-error")) {
			pred = rep(mean(Y), nrow(newX))
			fit = list(object=fit.nnet[1])
		} else {
			pred = predict(fit.nnet, newdata = newX, type = "raw")
			fit = list(object = fit.nnet)
		}
	}
	out = list(pred = pred, fit = fit)
	class(out$fit) = c("SL.nnet")
	return(out)
}

#' @export
SL.lasso.LT = function(Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 4, nlambda = 100, useMin = TRUE, ...) {
	library("glmnet")
	# X must be a matrix, should we use model.matrix or as.matrix
	if(!is.matrix(X)) {
		X = model.matrix(~ -1 + ., X)
		newX = model.matrix(~ -1 + ., newX)
	}
	# now use CV to find lambda
	Y.matrix = cbind(1-Y,Y)
	fitCV = try(cv.glmnet(x = X, y = Y.matrix, weights = obsWeights, lambda = NULL, type.measure = 'deviance', nfolds = nfolds, family = family$family, alpha = alpha, nlambda = nlambda), silent=TRUE)
	if(inherits(fitCV, "try-error")) {
		pred = rep(mean(Y), nrow(newX))
		fit = list(object=fitCV[1])
	} else {
		# two options for lambda, fitCV$lambda.min and fitCV$lambda.1se
		pred = predict(fitCV$glmnet.fit, newx = newX, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), type = 'response')
		fit = list(object = fitCV, useMin = useMin)
	}
	class(fit) = c('SL.glmnet')
	out = list(pred = pred, fit = fit)
	return(out)
}

#' @export
SL.ridge.LT = function (Y, X, newX, family, obsWeights, id, alpha = 0, nfolds = 4, nlambda = 100, useMin = TRUE, ...) {
	library("glmnet")
	# X must be a matrix, should we use model.matrix or as.matrix
	if(!is.matrix(X)) {
		X = model.matrix(~ -1 + ., X)
		newX = model.matrix(~ -1 + ., newX)
	}
	# now use CV to find lambda
	Y.matrix = cbind(1-Y,Y)
	fitCV <- try(cv.glmnet(x = X, y = Y.matrix, weights = obsWeights, lambda = NULL, type.measure = "deviance", nfolds = nfolds, family = family$family, alpha = alpha, nlambda = nlambda), silent=TRUE)
	if(inherits(fitCV, "try-error")) {
		pred = rep(mean(Y), nrow(newX))
		fit = list(object=fitCV[1])
	} else {
		# two options for lambda, fitCV$lambda.min and fitCV$lambda.1se
		pred <- predict(fitCV$glmnet.fit, newx = newX, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), type = "response")
		fit <- list(object = fitCV, useMin = useMin)
	}
	class(fit) <- c("SL.glmnet")
	out <- list(pred = pred, fit = fit)
	return(out)
}

#' @export
SL.glmnet.LT = function (Y, X, newX, family, obsWeights, id, alpha = 0.5, nfolds = 4, nlambda = 100, useMin = TRUE, ...) {
	library("glmnet")
	# X must be a matrix, should we use model.matrix or as.matrix
	if(!is.matrix(X)) {
		X = model.matrix(~ -1 + ., X)
		newX = model.matrix(~ -1 + ., newX)
	}
	# now use CV to find lambda
	Y.matrix = cbind(1-Y,Y)
	fitCV <- try(cv.glmnet(x = X, y = Y.matrix, weights = obsWeights, lambda = NULL, type.measure = "deviance", nfolds = nfolds, family = family$family, alpha = alpha, nlambda = nlambda), silent=TRUE)
	if(inherits(fitCV, "try-error")) {
		pred = rep(mean(Y), nrow(newX))
		fit = list(object=fitCV[1])
	} else {
		# two options for lambda, fitCV$lambda.min and fitCV$lambda.1se
		pred <- predict(fitCV$glmnet.fit, newx = newX, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), type = "response")
		fit <- list(object = fitCV, useMin = useMin)
	}
	class(fit) <- c("SL.glmnet")
	out <- list(pred = pred, fit = fit)
	return(out)
}

#' @export 
method.NNloglik.LT = function() {
	out = list(
			require = NULL,
			computeCoef = function(Z, Y, libraryNames, verbose, obsWeights, control, ...) {
				# compute cvRisk
				cvRisk = apply(Z, 2, function(x) { -mean(obsWeights * (Y*log(x) + (1-Y)*log(1-x))) } )
				names(cvRisk) = libraryNames
				# compute coef
				.NNloglik = function(x, y, wt, start = rep(0, ncol(x))) {
					# adapted from MASS pg 445
					fmin = function(beta, X, y, w) {
						p = plogis(crossprod(t(X), beta))
						-sum(2 * w * (y*log(p) + (1-y)*log(1-p)))
					}
					gmin = function(beta, X, y, w) {
						eta = X %*% beta
						p = plogis(eta)
						-2 * t(w * dlogis(eta) * (y/p + -1*(1-y)/(1-p))) %*% X
					}
					fit = optim(start, fmin, gmin, X = x, y = y, w = wt, method = "L-BFGS-B", lower = 0, control=list(parscale= rep(1/length(start), length(start))), ...)
					invisible(fit)
				}
				tempZ = trimLogit(Z, trim = control$trimLogit)
				fit.nnloglik = .NNloglik(x = tempZ, y = Y, wt = obsWeights)
				if(verbose) {
					message(paste("Non-Negative log-likelihood convergence: ", fit.nnloglik$convergence == 0))
				}
				initCoef = fit.nnloglik$par
				initCoef[initCoef < 0] = 0.0
				initCoef[is.na(initCoef)] = 0.0
				# normalize so sum(coef) = 1 if possible
				if(sum(initCoef) > 0) {
					coef = initCoef/sum(initCoef)
				} else {
					warning("All algorithms have zero weight", call. = FALSE)
					coef = initCoef
				}
				out = list(cvRisk = cvRisk, coef = coef)
				return(out)
			},
			computePred = function(predY, coef, control, ...) {
				out = plogis(crossprod(t(trimLogit(predY, trim = control$trimLogit)), coef))
				return(out)
			}
	)
	invisible(out)
}
