
#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

using namespace Rcpp;

class ProbEntry{
public:
	
	static const int s = 10;
	std::vector<double> * start;
	std::vector<double> prob;

	ProbEntry ();


	ProbEntry (std::vector<double> starts );

	void prepareIceCream();

	int pos4val ( double val, int state );

	void estimate ( std::vector<double> values, int state );

	// this return value will be logged in the ProbList
	double Prob_4_value ( double val, int state) ;
	void print( );
};