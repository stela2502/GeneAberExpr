#include "TransmissionProb.h"

TransmissionProb::TransmissionProb() {} 

TransmissionProb::TransmissionProb( int states) {
	setup( states );
}

void TransmissionProb::setup(int states) {
	if ( this->size > 0 ) {
		Rcpp::stop("TransmissionProb must not be set up twice!" );
	}
	//Rcout << "TransmissionProb setup() 1"<< std::endl;
	this->size = states;
	this->startProbability.resize( size);
	this->endProbability.resize( size);
	this->transitionProb.resize( size * size );

	//Rcout << "TransmissionProb setup() 2"<< std::endl;
	this->size = states;
	double null = 0.0;
	std::fill(this->startProbability.begin(), this->startProbability.end(), null);
	std::fill(this->endProbability.begin(), this->endProbability.end(), null);
	std::fill(this->transitionProb.begin(), this->transitionProb.end(), null);
	/*Rcout << "TransmissionProb setup() 3"<< std::endl;
	std::vector<std::string> tmp(states);
	std::fill(tmp.begin(), tmp.end(), "unknown" );
	print(tmp);*/
}
void TransmissionProb::setStart ( std::vector<double> starts ) {
	if ( this->size != starts.size() ){
		Rcout << "TransmissionProb not initialized correctly " << this->size << " != " << starts.size()<< std::endl;
		Rcpp::stop( "<- ERROR" );
	}
	for ( int i = 0; i < starts.size(); i ++ ){
		this->startProbability.at(i) = starts.at(i);
	}
}
void TransmissionProb::setEnd ( std::vector<double> ends ) {
	for ( int i = 0; i < ends.size(); i ++ ){
		this->endProbability.at(i) = ends.at(i);
	}
}
void TransmissionProb::setState ( int from, std::vector<double> pVals ) {
	for ( int to = 0; to < pVals.size(); to ++) {
		this->transitionProb.at(position( from, to)) = pVals.at(to);
	}
}
double TransmissionProb::getTransmissionProb( int from, int to){
	return( log(this->transitionProb.at(position(from, to)) ));
}
double TransmissionProb::getStartProb( int state ) {
	return ( log( this->startProbability.at(state) ));
}
double TransmissionProb::getEndProb( int state ) {
	return ( log(this->endProbability.at(state) ));
}
void TransmissionProb::print( std::vector<std::string> states ) {
	Rcout << "A c++ TransmissionProb object containing " << this->size << " states:" << std::endl;
	Rcout << "Start probabilities: " ;
	for ( int i = 0; i <this->size; i ++ ) {
		Rcout <<  this->startProbability.at(i) << " ";
	}
	Rcout << std::endl << "End probabilities: "  ;
	for ( int i = 0; i <this->size; i ++ ) {
		Rcout <<  this->endProbability.at(i) << " ";
	}
	Rcout << std::endl << "And transition probabilities: " << std::endl ;
	Rcout << "state ";
	states.resize( this->size, ".NA." );

	for ( int i = 0; i <this->size; i ++ ) {
		Rcout << states[i] << " ";
	}
	Rcout << std::endl;

	for ( int i = 0; i <this->size; i ++ ) {
		Rcout << states.at(i) << " ";
		
		for ( int a = 0; a <this->size; a ++ ) {
			Rcout << this->transitionProb.at(position(i, a)) << " ";
		}
		Rcout << std::endl;
	}
}
void TransmissionProb::prepareIceCream( bool phony){
	// 0 == C and 1 == H
	// setup( 2 ); // is already set up!
	//Rcout << "TransmissionProb::prepareIceCream 1" <<std::endl;
	this->startProbability.at(0) = 0.5;
	this->startProbability.at(1) = 0.5;
	this->endProbability.at(0) = 0.1;
	this->endProbability.at(1) = 0.1;
	//Rcout << "TransmissionProb::prepareIceCream 2" <<std::endl;

	this->transitionProb.at( position( 0,0) ) = 0.8;
	this->transitionProb.at( position( 1,1) ) = 0.8;
	this->transitionProb.at( position( 0,1) ) = 0.1;
	this->transitionProb.at( position( 1,0) ) = 0.1;

	//Rcout << "TransmissionProb::prepareIceCream 3" <<std::endl;
}


int TransmissionProb::position( int from, int to ) {
	return ( from + to * this->startProbability.size());
}

