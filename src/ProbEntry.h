#ifndef PROBENTRY_H
#define PROBENTRY_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

using namespace Rcpp;

class ProbEntry{
public:
	
	static const int s = 10;
	std::vector<double> start;
	std::vector<double> prob;
	int states;

	ProbEntry ();

	ProbEntry (std::vector<double>* starts, int states );

	void prepareIceCream(bool phony);

	int pos4val ( double val, int state );

	void estimate ( std::vector<double> values, int state );

	// this return value will be logged in the ProbList
	double Prob_4_value_e ( double val, int state) ;
	void print( );
};

#endif