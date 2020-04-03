#ifndef PROBLIST_H
#define PROBLIST_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

using namespace Rcpp;

#include "ProbEntry.h"
#include "TransmissionProb.h"


class ProbList {
public:


	TransmissionProb* transmissionProb;
	std::vector<ProbEntry*> list;
	//a vector of chars with max 10 entries...
	std::vector<std::string> states;

	std::vector<double> * starts;

	ProbList () {};

	ProbList (int marcowLength, int states, std::vector<double> * s );

	void setup(int marcowLength, int states, std::vector<double> * s );

	//method to fill the obejct with the ice cream data
	void prepareIceCream (bool phony);

	//std::vector<double> getStarts( );

	NumericMatrix as_Matrix ( );

	void estimate( Eigen::SparseMatrix<double> data, std::vector<double> range, std::string name, bool phony );

	void estimate( Eigen::SparseMatrix<double> data, std::vector<int> cols,	std::vector<double> range, std::string name, bool phony );

	double Prob_4_value ( int i, double val, int state);

	double Prob_4_value ( double val, int state);

	void print ();

};

#endif