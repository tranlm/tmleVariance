###############################################################################
# Description: Add comment
#
# Author: Linh Tran <tranlm@berkeley.edu>
# Date: Nov 17, 2015
###############################################################################


#' @export
Bound = function(x, bounds) {
	stopifnot(length(bounds) == 2 && !any(is.na(bounds)))
	x[x < min(bounds)] <- min(bounds)
	x[x > max(bounds)] <- max(bounds)
	return(x)
}

