
#' @name calcQuality
#' @aliases calcQuality,GeneAberExpr-method
#' @rdname calcQuality-methods
#' @docType methods
#' @description Calculate the sqared difference between two p value models for the HMM analysis.
#' @param model the model matrix obtained from a GetTestModel() run.
#' @title description of function calcQuality
#' @export 
calcQuality <- function(x){
	unlist( lapply( seq( 1,nrow(x),2), function(id){sum((model[id,] - model[id+1,])^2)}) )
}


#' @name calcCancerCount
#' @aliases calcCancerCount,GeneAberExpr-method
#' @rdname calcCancerCount-methods
#' @docType methods
#' @description Calculate difference between very likely (>0.99) and very unlikely
#'  (<0.01) p value of the IdentifyStates result matrix.
#' @param res the result matrix obtained from a IdentifyStates() run.
#' @title description of function calcCancerCount
#' @export 
calcCancerCount <- function(res) {
	apply( res,2, function(x) {length(which(x > .99))- length(which(x < 0.01))} )
}