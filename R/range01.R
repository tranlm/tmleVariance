###############################################################################
# Description: Add comment
#
# Author: Linh Tran <tranlm@berkeley.edu>
# Date: Nov 17, 2015
###############################################################################


#' @export
range01 = function(x, ...) (x - min(x, ...)) / (max(x, ...) - min(x, ...))

#' @export
range01.inverse = function(y, x, ...) y*(max(x, ...) - min(x, ...)) + min(x, ...) 

