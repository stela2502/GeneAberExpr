#ifndef TRANSMISSIONPROB_H
#define TRANSMISSIONPROB_H

#include <Rcpp.h>

using namespace Rcpp;

class TransmissionProb{
public:

	std::vector<double> startProbability;
	std::vector<double> endProbability;
	std::vector<double> transitionProb;
	int size = 0;

	TransmissionProb();
	TransmissionProb( int states);
	void setup( int states);
	void setStart ( std::vector<double> starts );
	void setEnd ( std::vector<double> ends );
	void setState ( int from, std::vector<double> pVals );
	double getTransmissionProb( int from, int to);
	double getStartProb( int state );
	double getEndProb( int state ) ;
	void print( std::vector<std::string> states ) ;
	void prepareIceCream(bool phony);
	int position( int from, int to );
};

#endif