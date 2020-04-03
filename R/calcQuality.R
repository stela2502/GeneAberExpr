#' GeneAberExpr: A package to identify chromosomal aberationd based on single cell expression data.
#'
#' In short the package used the spase matrix most single cell data is stored in.
#' The main parts are implemented in c++ for speed.
#' 
#' Short way to run the analysis:
#' 
#' myData # is a dgCMatrix containing ordered single cell expression from one chromosome
#' cancer # the id's of a sample of interesting cells
#' normal # the ids of a comare group of cells
#' min=min(myData@x)
#' if( min > 0 ){ min = 0 } ## no '0's in this object...
#' max = quantile(myData@x, .99)
#' model = GetTestModel( myData, seq( min, max,(max - min)/9 ), cancer, normal )
#' quality = unlist( lapply( seq( 1,nrow(model),2), calcQuality ) )
#' res = IdentifyStates( chrSums[ sort(order(quality,decreasing=T)[1:800]),] , seq( min, max,(max - min)/9 ) , cancer, normal )
#' cancerCount = apply( res,2, calcCancerCount )
#' 
#'
#' @docType package
#' @name GeneAberExpr
#' @useDynLib GeneAberExpr, .registration=TRUE
#' @importFrom Rcpp sourceCpp

NULL;

#' @name calcQuality
#' @aliases calcQuality,GeneAberExpr-method
#' @rdname calcQuality-methods
#' @docType methods
#' @description Calculate the sqared difference between two p value models for the HMM analysis.
#' @param model the model matrix obtained from a GetTestModel() run.
#' @title description of function calcQuality
#' @export 
calcQuality <- function(model){
	unlist( lapply( seq( 1,nrow(model),2), function(id){sum((model[id,] - model[id+1,])^2)}) )
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