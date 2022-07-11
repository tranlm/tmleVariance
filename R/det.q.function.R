###############################################################################
# Description: Add comment
#
# Author: Linh Tran <tranlm@berkeley.edu>
# Date: Nov 17, 2015
###############################################################################


#' @export 
det.q.function <- function(data, current.node, nodes, called.from.estimate.g) {
	if (!any(nodes$Y < current.node)) return(NULL)
	prev.Y <- data[, nodes$Y[nodes$Y < current.node], drop=F]
	prev.Y[is.na(prev.Y)] <- 0
	is.deterministic <- apply(prev.Y == 1, 1, any)
	Q.value <- data[is.deterministic, max(nodes$Y)] #this is 0 before scaling but may be nonzero after scaling
	return(list(is.deterministic=is.deterministic, Q.value=Q.value))   
}

