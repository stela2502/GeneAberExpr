#ifndef MARKOWCHAIN_H
#define MARKOWCHAIN_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

using namespace Rcpp;

#include "ProbList.h"


class MarcowChain{
public:
	std::vector<double> forwardResults;
	std::vector<double> backwardResults;
	std::vector<double> a_H;
	std::vector<double> p_H;
	std::vector<double> C;
	int chainLength;

	MarcowChain();

	MarcowChain( int chainLength, int hiddenStates );

	void setup( int chainLength, int hiddenStates );

	int at( int id, int state );

	NumericMatrix results( ProbList* model );

	NumericMatrix to_matrix( ProbList* model );

	NumericMatrix run ( std::vector<double> chain, ProbList* model, bool test );

	void CalculateForwardProbability( std::vector<double> chain, ProbList* model  );

	void CalculateBackwardProbability( std::vector<double> chain, ProbList* model  ) ;

	void CalculateTotalProbabilityFromStartToEnd( std::vector<double> chain,  ProbList* model  );

	double Add2Logs (double logA, double logB ) ;

};

#endif