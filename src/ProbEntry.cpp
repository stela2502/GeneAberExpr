#include "../inst/include/GeneAberExpr.h"


ProbEntry::ProbEntry () {
	::Rf_error("do never init like that!" );
}

ProbEntry::ProbEntry (std::vector<double> starts ) {
	if ( starts.size() != this->s ) {
		::Rf_error("the start sizes need to be exactly 10!" );
	};
	this->start = &starts;
	this->prob.reserve( this->s *2 ); // no one would use this without 2 states
	for ( int i = 0; i < this->s ; i++) {
		this->prob.push_back(0.0);
		this->prob.push_back(0.0);// no one would use this without 2 states
	}
}

void ProbEntry::prepareIceCream( ){
	// start has already been set!		

	this->prob.reserve( 2 * this->s );
	for ( int i = 0; i < 2 * this->s; i++){
		this->prob.push_back(0.0);
	}
	this->prob[ pos4val(1,0) ] = 0.7;
	this->prob[ pos4val(2,0) ] = 0.2;
	this->prob[ pos4val(3,0) ] = 0.1;
	this->prob[ pos4val(1,1) ] = 0.1;
	this->prob[ pos4val(2,1) ] = 0.2;
	this->prob[ pos4val(3,1) ] = 0.7;

}

int ProbEntry::pos4val ( double val, int state ){
	if (  this->start->size() != 10 ){
		// first time usage...
		Rcout << "I currently have "<<this->start->size()<<
		" start positions - not the correct number of 10!" <<std::endl;
		for ( int i = 0; i < this->start->size(); i ++ ) {
			Rcout << this->start->at(i)<<" ";
		}
		Rcout <<std::endl;
		::Rf_error("object not initialized correctly - starts are empty / too many!? ");
	}

	if ( val < this->start->at(1) ){
		return 0 + this->s * state;
	}
	for (int i =1; i< this->s-1; i++){
		if ( this->start->at(i) <= val and this->start->at(i+1) > val ){
			return i + this->s * state;
		}
	}
	return this->s - 1 + this->s * state;
}

void ProbEntry::estimate ( std::vector<double> values, int state ) {

	//Rcout << "filling ProbEntry with size "<<  this->s << " and values size "<< values.size() << " for starte " << state <<std::endl;
	
	if ( this->prob.size() < state * this->s ){
		this->prob.reserve( state * this->s );
		for ( int i = this->prob.size(); i <state * this->s; i ++ ){
			this->prob.push_back(0.0);
		}
	}
	//this->print();
	std::vector<double> tmp (this->s);
	std::vector<double> prob (this->s);

	std::fill(tmp.begin(), tmp.end(), 1e-9 );
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
	sum = tmp[0]+ static_cast<double>(tmp[this->s-1]);
	min= tmp[this->s-1]; // educated guess here -really
	for( i = 1; i < this->s -1; i ++) {
		//the smoothing might be a horrible idea...
		//tmp[i] = (tmp[i] + tmp[i-1] + tmp[i+1]) / 3;
		sum += tmp[i];
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
		this->prob.at(pos4val( this->start->at(i), state) ) = tmp[i];
	}
	//Rcout << "correct size?" << this->prob.size() <<  std::endl;
	//print();
	//Rcout << "survived!" << std::endl;
}

// this return value will be logged in the ProbList
double ProbEntry::Prob_4_value ( double val, int state) {
	return(this->prob[ pos4val( val, state ) ]);
}

void ProbEntry::print( ) {
	
	Rcout << "A c++ ProbEntry object containing " << this->prob.size() / this->s << 
	" state(s) and " <<  this->start->size() << " start points" <<std::endl;
	
	Rcout << "start " << ": ";
	for ( int a=0; a < this->s; a ++) {
		Rcout << this->start->at(a) << " ";
	}
	Rcout << std::endl;

	for ( int state=0; state < this->prob.size() / this->s; state ++) {
		Rcout << "state "<< state << ": ";
		for ( int a =0; a< this->s; a++ ) {
			Rcout << this->prob[pos4val( this->start->at(a), state)] << " ";
		}
		Rcout << std::endl;
	}
}
