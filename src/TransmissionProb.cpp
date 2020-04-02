#include "../inst/include/GeneAberExpr.h"


TransmissionProb::TransmissionProb() {} 
TransmissionProb::TransmissionProb( int states) {
	this->size = states;
	this->startProbability.resize(states);
	this->endProbability.resize(states);
	this->transitionProb.resize(states* states);
	for ( int i = 0; i < states; i++) {
		this->startProbability.push_back(0.0);
		this->endProbability.push_back(0.0);
	}
	for ( int i = 0; i < states * states; i++) {
		this->transitionProb.push_back(0.0);
	}
}

void TransmissionProb::setStart ( std::vector<double> starts ) {
	if ( this->size != starts.size() ){
		 Rcout << "TransmissionProb not initialized correctly " << this->size << " != " << starts.size()<< std::endl;
		::Rf_error( "<- ERROR" );
	}
	this->startProbability = starts;
}
void TransmissionProb::setEnd ( std::vector<double> ends ) {
	this->endProbability = ends;
}
void TransmissionProb::setState ( int from, std::vector<double> pVals ) {
	for ( int to = 0; to < pVals.size(); to ++) {
		this->transitionProb[position( from, to)] = pVals[to];
	}
}
double TransmissionProb::getTransmissionProb( int from, int to){
	return( log(this->transitionProb[position(from, to)]));
}
double TransmissionProb::getStartProb( int state ) {
	return ( log( this->startProbability[state]));
}
double TransmissionProb::getEndProb( int state ) {
	return ( log(this->endProbability[state]));
}
void TransmissionProb::print() {
	Rcout << "A c++ TransmissionProb object containing:"  << std::endl;
	//Rcout << "Start probabilities: "  << std::ostream_iterator<double>(
	// this->startProbability, " ") << std::endl;
	//for ( int i =0; i < this->transitionProb.size() ){
	//	std::ostream_iterator<double>( this->transitionProb[i], " ") << std::endl;
	//}
	//Rcout << "End probabilities: "  << std::ostream_iterator<double>( 
//		this->endProbability, " ") << std::endl;
}
void TransmissionProb::prepareIceCream(){
	// 0 == C and 1 == H
	this->startProbability[0] = 0.5;
	this->startProbability[1] = 0.5;
	this->endProbability[0] = 0.1;
	this->endProbability[1] = 0.1;
	std::fill(this->transitionProb.begin(),  this->transitionProb.end(), 0.0);
	this->transitionProb[ position( 0,0) ] = 0.8;
	this->transitionProb[ position( 1,1) ] = 0.8;
	this->transitionProb[ position( 0,1) ] = 0.1;
	this->transitionProb[ position( 1,0) ] = 0.1;
}


int TransmissionProb::position( int from, int to ) {
	int states;
	states = this->startProbability.size();
	return ( from + to * states);
}

