
#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

using namespace Rcpp;

class ProbList {
public:

	std::vector<ProbEntry*> list;
	TransmissionProb* transmissionProb;
	//a vector of chars with max 10 entries...
	std::vector<std::string> states;

	std::vector<double>  starts;

	ProbList ();

	//method to fill the obejct with the ice cream data
	void prepareIceCream ();

	//std::vector<double> getStarts( );

	NumericMatrix as_Matrix ( );

	void estimate( Eigen::SparseMatrix<double> data, std::vector<double> range, std::string name, bool phony );

	void estimate( Eigen::SparseMatrix<double> data, std::vector<int> cols,	std::vector<double> range, std::string name, bool phony );

	double Prob_4_value ( int i, double val, int state);

	double Prob_4_value ( double val, int state);

	void print ();

};

