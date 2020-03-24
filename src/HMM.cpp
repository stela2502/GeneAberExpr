#include "../inst/include/GeneAberExpr.h"

#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

using namespace Rcpp;


// here I need to allow R to run a HMM analysis.
// I need the forwards/backward and all other HMM algorithms in place here and work on sparse matrices.
// should be quite straight forward based on my perl hmm class.
// all probability estimations are in the <ProbList> class.
// get a functional one with ProbList( <sparseMatrix>, idsA, idsB, starts )
// the starts should range from min to max of the matrix.


//' @title IdentifyStates tries to infere cell states from a sparse Matrix and a list of example cells.
//' @aliases IdentifyStates,GeneAberExpr-method
//' @rdname IdentifyStates
//' @description Mased on the example data a model for each genomic area is developed
//' Based on these per area models the probability to be in either of these states is calculated.
//' @param X the sparse matrix (tests are applied to columns!)
//' @param range the starting points in the size distribution of the data (10!!)
//' @param interest row IDs for the group of interest
//' @param background row IDS for the background
//' @return a matrix with the HMM internals (at the moment)
//' @export
// [[Rcpp::export]]
NumericMatrix IdentifyStates ( Eigen::MappedSparseMatrix<double> data, std::vector<double> range,  std::vector<int> interest,
		std::vector<int> background ){

	ProbList* model = new ProbList(); // 10 by default!
	Rcout << "estimating 'interest' probability functions " << std::endl;
	std::string a ="interest";
	model->estimate( data, interest, range, a );

	Rcout << "estimating 'background' probability functions " << std::endl;
	a = "Background";
	model->estimate( data, background, range, a );

	model->transmissionProb = new TransmissionProb( 1 );
	//transmissionProb
	std::vector<double> start = {0.5, 0.5};
	model->transmissionProb->setStart ( start );
	model->transmissionProb->setEnd   ( start );
	// from Interest to Interest and to Background
	start = {0.9, 0.00001};
	model->transmissionProb->setState ( 0, start);
	// from Background to Interest and to Background
	start = {0.00001, 0.9};
	model->transmissionProb->setState ( 1, start);

	NumericMatrix ret;
	MarcowChain* mc;
	mc = new MarcowChain ( data.outerSize(), 2);
	std::vector<double>D ( data.outerSize() );

	for ( int c_ = 0; c_< data.outerSize(); ++c_ ) {
		std::fill(D.begin(), D.end(), 0.0);
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, c_); it; ++it){
			D[it.row()] =  it.value();
		}
		ret = mc->run( D, model);
		return (ret);
	}

	return (model->as_Matrix());
}

//' @title GetTestModel runs sthe estimation only and returns a matrix
//' @aliases GetTestModel,GeneAberExpr-method
//' @rdname GetTestModel
//' @description Mased on the example data a model for each genomic area is developed
//' Based on these per area models the probability to be in either of these states is calculated.
//' @param X the sparse matrix (tests are applied to columns!)
//' @param range the starting points in the size distribution of the data (10!!)
//' @param interest row IDs for the group of interest
//' @param background row IDS for the background
//' @return a matrix with the HMM internals (at the moment)
//' @export
// [[Rcpp::export]]
NumericMatrix GetTestModel ( Eigen::MappedSparseMatrix<double> data, std::vector<double> range,  std::vector<int> interest,
		std::vector<int> background ){

	ProbList* model = new ProbList(); // 10 by default!
	Rcout << "estimating 'interest' probability functions " << std::endl;
	std::string a ="interest";
	model->estimate( data, interest, range, a );

	//model->list.at(1)->print();

	//model->print();

	Rcout << "estimating 'background' probability functions " << std::endl;
	a = "Background";
	model->estimate( data, background, range, a );

	//model->print();

	model->transmissionProb = new TransmissionProb( 1 );
	//transmissionProb
	std::vector<double> start = {0.5, 0.5};
	model->transmissionProb->setStart ( start );
	model->transmissionProb->setEnd   ( start );
	// from Interest to Interest and to Background
	start = {0.9, 0.00001};
	model->transmissionProb->setState ( 0, start);
	// from Background to Interest and to Background
	start = {0.00001, 0.9};
	model->transmissionProb->setState ( 1, start);

	NumericMatrix ret;
	ret = model->as_Matrix();
	return ( ret ) ;
}