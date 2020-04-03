#include "ProbEntry.h"


ProbEntry::ProbEntry () {
	Rcpp::stop("do never init like that!" );
}

ProbEntry::ProbEntry (std::vector<double> * starts, int states ) {
	//Rcout << "ProbEntry::ProbEntry( std::vector<double>* , int): [1]" << std::endl;
	if ( starts->size() != this->s ) {
		Rcpp::stop("the start sizes need to be exactly 10!" );
	};
	//Rcout << "ProbEntry::ProbEntry( std::vector<double>* , int): [2]" << std::endl;
	this->start.resize(this->s);
	for ( int i = 0; i < this->s; i ++ ){
		this->start.at(i) = starts->at(i);
	}
	//Rcout << "ProbEntry::ProbEntry( std::vector<double>* , int): [3]" << std::endl;
	this->states = states;
	this->prob.resize( this->s * states );
	double null = 0.0;
	std::fill(this->prob.begin(), this->prob.end(), null);
	//Rcout << "ProbEntry::ProbEntry( std::vector<double>* , int): [4]" << std::endl;
	//print();
}

void ProbEntry::prepareIceCream( bool phony ){
	// start has already been set up!		
	if ( this->start.size() != 10 ){
		Rcpp::stop("ProbEntry::PrepareIceCream() - starts are empty WHY!? ");
	}
	
	this->prob.at(pos4val(1,0)) = 0.7;
	this->prob.at(pos4val(2,0)) = 0.2;
	this->prob.at(pos4val(3,0)) = 0.1;
	this->prob.at(pos4val(1,1)) = 0.1;
	this->prob.at(pos4val(2,1)) = 0.2;
	this->prob.at(pos4val(3,1)) = 0.7;

}

int ProbEntry::pos4val ( double val, int state ){
	if (  this->start.size() != 10 ){
		// first time usage...
		Rcout << "I currently have "<<this->start.size()<<
		" start positions - not the correct number of 10!" <<std::endl;
		int max = this->start.size();
		if ( max > 10 ) {
			max = 10;
		}
		for ( int i = 0; i < max; i ++ ) {
			Rcout << this->start.at(i)<<" ";
		}
		Rcout <<std::endl;
		Rcpp::stop("object not initialized correctly - starts are empty / too many!? " + 
			std::to_string( this->start.size() )  );
	}
	int ret = -1;
	if ( val < this->start.at(1) ){
		ret = 0 + this->s * state;
	}else {
		for (int i =1; i< this->s-1; i++){
			if ( this->start.at(i) <= val and this->start.at(i+1) > val ){
				ret = i + this->s * state;
			}
		}
	}
	if ( ret == -1){
		ret = this->s - 1 + this->s * state;
	}
	//Rcout << "ProbEntry::pos4val(" << val <<", " << state << ") = "<< ret << std::endl;
	return ret;
}

void ProbEntry::estimate ( std::vector<double> values, int state ) {

	//Rcout << "filling ProbEntry with size "<<  this->s << " and values size "<< values.size() << " for starte " << state <<std::endl;
	
	if ( this->prob.size() < state * this->s ){
		Rcpp::stop("ProbEntry::estimate - object is not initialized correctly!");
	}
	//this->print();
	std::vector<double> tmp (this->s);
	std::vector<double> prob (this->s);

	std::fill(tmp.begin(), tmp.end(), 0.0 );
	double sum = 0.0;
	int i;
	// fill an initial table
	//Rcout << "processing initial values "  << values.size() << std::endl;
	for( i = 0; i < values.size(); i ++) {
		//Rcout << "for loop "  << i << " value " << values[i]<< " pos4val = "<< pos4val(values[i]) << std::endl;
		tmp[ pos4val(values[i], 0) ] ++;

	}

	// smooth this
	//Rcout << "smoothing table "  << tmp.size() << std::endl;
	double min = values.size();
	sum = values.size();
	min= tmp[this->s-1]; // educated guess here -really
	for( i = 0; i < this->s -1; i ++) {
		if ( min > tmp[i] ){
			min = tmp[i];
		}
	}
	//Rcout << "rescaled and smoothed the table"  << std::endl;
	for ( i = 0; i < this->s; i++ ){
		//Rcout << "(prob[i] - min)" << (prob[i] - min) << " / sum " <<  sum << "+ 1e-9" << std::endl;
		tmp[i] = ((tmp[i] - min) / sum ) + 1e-9 ; //we can not have 0 values later on!
	}


	//Rcout << "add the probabilites at "  << state << std::endl;
	for ( int i = 0; i < this->s; i++ ){
		this->prob.at(pos4val( this->start.at(i), state) ) = tmp[i];
	}
	//Rcout << "correct size?" << this->prob.size() <<  std::endl;
	//print();
	//Rcout << "survived!" << std::endl;
}

// this return value will be logged in the ProbList
double ProbEntry::Prob_4_value_e ( double val, int state) {
	return(this->prob.at( pos4val( val, state ) ) );
}

void ProbEntry::print( ) {
	
	Rcout << "A c++ ProbEntry object containing " << (this->prob.size() / this->s) << 
	" state(s) and " <<  this->start.size() << " start points ("<< &this->start <<")" <<std::endl;
	
	Rcout << "start " << ": ";
	for ( int a=0; a < this->s; a ++) {
		Rcout << this->start.at(a) << " ";
	}
	Rcout << std::endl;

	for ( int state=0; state < this->states; state ++) {
		Rcout << "state "<< state << ": ";
		for ( int a =0; a< this->s; a++ ) {
			Rcout << this->prob.at(pos4val( this->start.at(a), state)) << " ";
		}
		Rcout << std::endl;
	}
}
