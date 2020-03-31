#include "../inst/include/GeneAberExpr.h"

#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

using namespace Rcpp;

typedef Eigen::SparseMatrix<double>::InnerIterator SpMatIt; 

// here I need to allow R to run a HMM analysis.
// I need the forwards/backward and all other HMM algorithms in place here and work on sparse matrices.
// should be quite straight forward based on my perl hmm class.
// all probability estimations are in the <ProbList> class.
// get a functional one with ProbList( <sparseMatrix>, idsA, idsB, starts )
// the starts should range from min to max of the matrix.


//' @title IdentifyStates tries to infere cell states from a sparse Matrix and a list of example cells.
//' @aliases IdentifyStates,GeneAberExpr-method
//' @rdname IdentifyStates
//' @description Based on the example data a model for each genomic area is developed.
//' Based on these per area models the probability to be in either of these states is calculated.
//' @param X the sparse matrix (tests are applied to columns!)
//' @param range the starting points in the size distribution of the data (10!!)
//' @param interest row IDs for the group of interest
//' @param background row IDS for the background
//' @param phony print status messages (default false)
//' @return a matrix with the probabilites for each region to be in 'interest' state.
//' @export
// [[Rcpp::export]]
NumericMatrix IdentifyStates ( Eigen::SparseMatrix<double> data, std::vector<double> range,  std::vector<int> interest,
		std::vector<int> background, bool phony ){

	ProbList* model = new ProbList(); // 10 by default!
	Rcout << "estimating 'interest' probability functions " << std::endl;
	std::string a ="interest";
	model->estimate( data, interest, range, a, phony );

	Rcout << "estimating 'background' probability functions " << std::endl;
	a = "Background";
	model->estimate( data, background, range, a, phony );

	model->transmissionProb = new TransmissionProb( 2 );
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

	NumericMatrix tmp;
	MarcowChain* mc;
	Rcout << "Create MarcowChain object " << std::endl;
	mc = new MarcowChain ( data.innerSize(), 2);
	std::vector<double>D ( data.innerSize() );

	int modFac = data.outerSize()/100;
	Progress p( round(data.outerSize() / modFac), true);

	NumericMatrix ret( data.innerSize(), data.outerSize());
	
	for ( int c_ = 0;  c_ < data.outerSize(); ++c_ ) {
		std::fill(D.begin(), D.end(), 0.0);
		//for (SpMatIt it(data, i); it; ++it){
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, c_); it; ++it){
			D[it.row()] =  it.value();
		}
		tmp = mc->run( D, model, false); // returns only the p values not log!
		for ( int i = 0; i < data.innerSize(); i++){
			ret(i, c_) = tmp(i,0);
		}
		if ( c_ % modFac == 0) {
			p.increment();
		}
		
	}

	return ( ret );
}


//' @title IdentifyStatesTest Test version of 'IdentifyStates' returing both state probabilites for the first cell.
//' @aliases IdentifyStatesTest,GeneAberExpr-method
//' @rdname IdentifyStatesTest
//' @description Based on the example data a model for each genomic area is developed.
//' Based on these per area models the probability to be in either of these states is calculated.
//' @param X the sparse matrix (tests are applied to columns!)
//' @param range the starting points in the size distribution of the data (10!!)
//' @param interest row IDs for the group of interest
//' @param background row IDS for the background
//' @param phony print status messages (default false)
//' @return a matrix with the probabilites for each region to be in 'interest' state.
//' @export
// [[Rcpp::export]]
NumericMatrix IdentifyStatesTest ( Eigen::SparseMatrix<double> data, std::vector<double> range,  std::vector<int> interest,
		std::vector<int> background, bool phony ){

	ProbList* model = new ProbList(); // 10 by default!
	Rcout << "estimating 'interest' probability functions " << std::endl;
	std::string a ="interest";
	model->estimate( data, interest, range, a, phony );

	Rcout << "estimating 'background' probability functions " << std::endl;
	a = "Background";
	model->estimate( data, background, range, a , phony);

	model->transmissionProb = new TransmissionProb( 2 );
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

	//NumericMatrix ret( data.innerSize(), data.outerSize() );
	NumericMatrix tmp;
	MarcowChain* mc;
	Rcout << "Create MarcowChain object " << std::endl;
	mc = new MarcowChain ( data.innerSize(), 2);
	std::vector<double>D ( data.innerSize() );


	std::fill(D.begin(), D.end(), 0.0);
	//for (SpMatIt it(data, i); it; ++it){
	for (Eigen::SparseMatrix<double>::InnerIterator it(data, 0); it; ++it){
		D[it.row()] =  it.value();
	}
	tmp = mc->run( D, model, false); // returns only the p values not log!
	return(tmp);
	

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
//' @param phony print status messages (default false)
//' @return a matrix with the HMM internals (at the moment)
//' @export
// [[Rcpp::export]]
NumericMatrix GetTestModel ( Eigen::MappedSparseMatrix<double> data, std::vector<double> range,  std::vector<int> interest,
		std::vector<int> background, bool phony){


	ProbList* model = new ProbList(); // 10 by default!
	if ( phony ){
		Rcout << "estimating 'interest' probability functions " << std::endl;
	}
	std::string a ="interest";
	model->estimate( data, interest, range, a, phony );

	//model->list.at(1)->print();

	//model->print();
	if ( phony ){
		Rcout << "estimating 'background' probability functions " << std::endl;
	}
	a = "Background";
	model->estimate( data, background, range, a, phony );

	//model->print();

	model->transmissionProb = new TransmissionProb( 2 );
	//transmissionProb
	std::vector<double> start = {0.5, 0.5};
	//Rcout << "add starts vector of size " << start.size()<<  std::endl;
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


//' @title IceCreamTest processes the Ice cream HMM example
//' @aliases IceCreamTest,GeneAberExpr-method
//' @rdname IceCreamTest
//' @description the publication 'An Interactive Spreadsheet for Teaching the Forward-Backward Algorithm'
//' is a beatuful example of how to calculate an HMM. This function re-runs the excel data.
//' @return a matrix with the HMM internals (at the moment)
//' @export
// [[Rcpp::export]]
NumericMatrix IceCreamTest ( ){

	ProbList* model = new ProbList(); // 10 by default!
	//Rcout << "preparing probability functions" << std::endl;
	model->prepareIceCream();
	NumericMatrix ret;
	//Rcout << "preparing probability functions 2" << std::endl;

	//ret = model->as_Matrix();
	//return (ret);

	std::vector<double> ice { 2.0, 3.0, 3.0, 2.0, 3.0, 2.0, 3.0, 2.0, 2.0, 3.0, 1.0, 3.0, 3.0,
		1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 3.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 2.0, 
		3.0, 2.0,2.0 };

	//Rcout << "transition prob H->H .8 " << exp(model->transmissionProb->getTransmissionProb( 1, 1 ))<< std::endl;
	//Rcout << "transition prob H->C .1 " << exp(model->transmissionProb->getTransmissionProb( 1, 0 ))<< std::endl;
	//Rcout << "transition prob H->H .8 " << exp(model->transmissionProb->getTransmissionProb( 0, 0 ))<< std::endl;
	//Rcout << "transition prob C->H .1 " << exp(model->transmissionProb->getTransmissionProb( 0, 1 ))<< std::endl;



	MarcowChain* mc;
	//Rcout << "Create MarcovChain object " << std::endl;
	mc = new MarcowChain ( 33, 2);

	//Rcout << "Process iteration 1" << mc->chainLength << std::endl;

	ret = mc->run( ice, model, true ); // test returning all internal values

	return ( ret ) ;
}


//' @title goodGenes processes the Ice cream HMM example
//' @aliases goodGenes,GeneAberExpr-method
//' @rdname goodGenes
//' @description check in the sparse matrix if at least 40% of the row values are less than max
//' @return bool vector of rows with less than 40% over max columns.
//' @export
// [[Rcpp::export]]
std::vector<bool> goodGenes (Eigen::SparseMatrix<double> data, double max ){

	data = data.transpose();

	std::vector<bool> ret( data.outerSize() );
	//Rcout << "dims outer: " << data.outerSize() << " ; inner: " << data.innerSize() << std::endl;

	int out;
	for ( int c_ = 0;  c_ < data.outerSize(); ++c_ ) {
		out = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, c_); it; ++it){
			if (it.value() > max ){
				out ++;
			}
		}
		ret[c_] = static_cast<double>(out) / data.innerSize() < 0.4;	
	}

	return ( ret ) ;
}