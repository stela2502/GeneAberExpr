#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

using namespace Rcpp;

class TransmissionProb{
public:

	TransmissionProb();
	TransmissionProb( int states);
	void setStart ( std::vector<double> starts );
	void setEnd ( std::vector<double> ends );
	void setState ( int from, std::vector<double> pVals );
	double getTransmissionProb( int from, int to);
	double getStartProb( int state );
	double getEndProb( int state ) ;
	void print() ;
	void prepareIceCream();


private:
	std::vector<double> startProbability;
	std::vector<double> endProbability;
	std::vector<double> transitionProb;
	int size;
	int position( int from, int to );
};
