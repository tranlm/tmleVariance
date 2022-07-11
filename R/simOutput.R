###############################################################################
# Description: Outputs results into excel file
#
# Author: Linh Tran <tranlm@berkeley.edu>
# Date: Oct 7, 2015
###############################################################################


#' @export
simOutput = function(file, results, gmat, time.pt) {
	
	parms = loadWorkbook(file=file)
	sheets = getSheets(parms)
	colStyle1 = list(
			'1'=CellStyle(parms) + DataFormat("#,##0.000") + Alignment(h="ALIGN_CENTER"),
			'2'=CellStyle(parms) + DataFormat("#,##0.000") + Alignment(h="ALIGN_CENTER"),
			'3'=CellStyle(parms) + DataFormat("#,##0.000") + Alignment(h="ALIGN_CENTER"),
			'4'=CellStyle(parms) + DataFormat("#,##0.000") + Alignment(h="ALIGN_CENTER"),
			'5'=CellStyle(parms) + DataFormat("#,##0.000") + Alignment(h="ALIGN_CENTER"),
			'6'=CellStyle(parms) + DataFormat("#,##0.000") + Alignment(h="ALIGN_CENTER"),
			'7'=CellStyle(parms) + DataFormat("#,##0.000") + Alignment(h="ALIGN_CENTER")
	)
	
	###############
	## g(A=a|Pa) ##
	###############
	addDataFrame(x=rbind(mean(gmat$gmat0<0.001), cbind(summary(gmat$gmat0))), sheet=sheets[["Type1Error"]], row.names = FALSE, col.names = FALSE, startColumn = 2+time.pt, startRow=12, colStyle=list('1'=CellStyle(parms) + DataFormat("#,##0.0000") + Alignment(h="ALIGN_CENTER")))
	addDataFrame(x=rbind(mean(gmat$gmat1<0.001), cbind(summary(gmat$gmat1))), sheet=sheets[["Type1Error"]], row.names = FALSE, col.names = FALSE, startColumn = 2+time.pt, startRow=12+11, colStyle=list('1'=CellStyle(parms) + DataFormat("#,##0.0000") + Alignment(h="ALIGN_CENTER")))

	#############
	## RESULTS ##
	#############
	psi_0 = 0
	
	estimates = do.call("rbind", mclapply(results, function(x) return(c(unadj=x$unadj[["estimate.(Intercept)"]], iptw=x$iptw[["estimate.(Intercept)"]], ltmle=x$ltmle.package[["estimate.(Intercept)"]]))))
	psi.est = paste0(formatC(apply(estimates, 2, mean), digits=3, format="f"), " (", formatC(apply(estimates, 2, sd), digits=3, format="f"), ")")
	addDataFrame(x=psi.est, sheet=sheets[["Type1Error"]], row.names = FALSE, col.names = FALSE, startColumn = 2+time.pt, startRow=37, colStyle=colStyle1)	
	
	iptw = do.call("rbind", mclapply(results, function(x) return(c(se=x$iptw[["std.dev"]], coverage=ifelse(x$iptw[["CI1"]] <= psi_0 & psi_0 <= x$iptw[["CI2"]], 1, 0), sig=x$iptw[["pvalue.(Intercept)"]]<=.05))))
	iptw.se = c(sd(iptw[,"se"]), sqrt(var(iptw[,"coverage"] - mean(iptw[,"coverage"]))/nrow(iptw)), sqrt(var(iptw[,"sig"] - mean(iptw[,"sig"]))/nrow(iptw)))
	iptw.est = paste0(formatC(apply(iptw, 2, mean), digits=3, format="f"), " (", formatC(iptw.se, digits=3, format="f"), ")")
	addDataFrame(x=iptw.est, sheet=sheets[["Type1Error"]], row.names = FALSE, col.names = FALSE, startColumn = 2+time.pt, startRow=41, colStyle=colStyle1)
	
	ltmle.ic = do.call("rbind", mclapply(results, function(x) return(c(se=x$ltmle.ic[["std.dev"]], coverage=ifelse(x$ltmle.ic[["CI1"]] <= psi_0 & psi_0 <= x$ltmle.ic[["CI2"]], 1, 0), sig=x$ltmle.ic[["pvalue.(Intercept)"]]<=.05))))
	ltmle.ic.se = c(sd(ltmle.ic[,"se"]), sqrt(var(ltmle.ic[,"coverage"] - mean(ltmle.ic[,"coverage"]))/nrow(ltmle.ic)), sqrt(var(ltmle.ic[,"sig"] - mean(ltmle.ic[,"sig"]))/nrow(ltmle.ic)))
	ltmle.ic.est = paste0(formatC(apply(ltmle.ic, 2, mean), digits=3, format="f"), " (", formatC(ltmle.ic.se, digits=3, format="f"), ")")
	addDataFrame(x=ltmle.ic.est, sheet=sheets[["Type1Error"]], row.names = FALSE, col.names = FALSE, startColumn = 2+time.pt, startRow=45, colStyle=colStyle1)
	
	ltmle.robust = do.call("rbind", mclapply(results, function(x) return(c(se=x$ltmle.robust[["std.dev"]], coverage=ifelse(x$ltmle.robust[["CI1"]] <= psi_0 & psi_0 <= x$ltmle.robust[["CI2"]], 1, 0), sig=x$ltmle.robust[["pvalue.(Intercept)"]]<=.05))))
	ltmle.robust.se = c(sd(ltmle.robust[,"se"]), sqrt(var(ltmle.robust[,"coverage"] - mean(ltmle.robust[,"coverage"]))/nrow(ltmle.robust)), sqrt(var(ltmle.robust[,"sig"] - mean(ltmle.robust[,"sig"]))/nrow(ltmle.robust)))
	ltmle.robust.est = paste0(formatC(apply(ltmle.robust, 2, mean), digits=3, format="f"), " (", formatC(ltmle.robust.se, digits=3, format="f"), ")")
	addDataFrame(x=ltmle.robust.est, sheet=sheets[["Type1Error"]], row.names = FALSE, col.names = FALSE, startColumn = 2+time.pt, startRow=49, colStyle=colStyle1)
	
	ltmle.package = do.call("rbind", mclapply(results, function(x) return(c(se=x$ltmle.package[["std.dev"]], coverage=ifelse(x$ltmle.package[["CI1"]] <= psi_0 & psi_0 <= x$ltmle.package[["CI2"]], 1, 0), sig=x$ltmle.package[["pvalue.(Intercept)"]]<=.05))))
	ltmle.package.se = c(sd(ltmle.package[,"se"]), sqrt(var(ltmle.package[,"coverage"] - mean(ltmle.package[,"coverage"]))/nrow(ltmle.package)), sqrt(var(ltmle.package[,"sig"] - mean(ltmle.package[,"sig"]))/nrow(ltmle.package)))
	ltmle.package.est = paste0(formatC(apply(ltmle.package, 2, mean), digits=3, format="f"), " (", formatC(ltmle.package.se, digits=3, format="f"), ")")
	addDataFrame(x=ltmle.package.est, sheet=sheets[["Type1Error"]], row.names = FALSE, col.names = FALSE, startColumn = 2+time.pt, startRow=53, colStyle=colStyle1)
	
	## EXTRA INFO ##
	icMax = do.call("rbind", mclapply(results, function(x) return(c(icMax=x$ltmle.package[["std.dev"]]==x$ltmle.ic[["std.dev"]]))))
	addDataFrame(x=mean(icMax), sheet=sheets[["Type1Error"]], row.names = FALSE, col.names = FALSE, startColumn = 2+time.pt, startRow=64, colStyle=list('1'=CellStyle(parms) + DataFormat("#,##0.000") + Alignment(h="ALIGN_CENTER", wrapText=FALSE)))

	## SAVES / CLOSES ##
	saveWorkbook(parms, file=file)
	
	invisible(NULL)
	
}
