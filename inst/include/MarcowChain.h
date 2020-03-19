#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

class MarcowChain{
public:
	std::vector<double> forwardResultsA;
	std::vector<double> forwardResultsB;
	std::vector<double> backwardResultsA;
	std::vector<double> backwardResultsB;
	std::vector<double> a_HA;
	std::vector<double> a_HB;
	std::vector<double> C;
	

	NumericMatrix run ( Eigen::SparseMatrix<double> data, ProbList model ){
		/*$
		here we need to asigne each column a probability for A and for B

		marcovChain->CalculateForwardProbability();
		$marcovChain->CalculateBackwardProbability();
		$marcovChain->CalculateTotalProbabilityFromStartToEnd();
		*/


	}

	void CalculateForwardProbability( std::vector<double> chain ) {

	}

	void CalculateBackwardProbability( std::vector<double> chain ) {

	}

	void CalculateTotalProbabilityFromStartToEnd( std::vector<double> chain ) {

	}

	double Add2Logs (double logA, double logB ) {
		double ret =  logA + log( 1 + exp( logB - logA ) );
		return (ret);
	}

}