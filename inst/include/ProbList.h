#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>


using namespace Rcpp;

typedef std::pair<int, double> paired;

class ProbEntry{
public:
	
	static const int s = 10;
	std::vector<double> start;
	std::vector<double> prob;

	ProbEntry () {
		this->start.reserve(10);
		this->prob.reserve(0);
		::Rf_error("do never init like that!" );
	};

	ProbEntry (std::vector<double> starts ) {
		if ( starts.size() != this->s ) {
			::Rf_error("the start sizes need to be exactly 10!" );
		}

		this->prob.reserve( this->s *2 ); // no one would use this without 2 states
		this->start.reserve( this->s );
		for ( int i = 0; i < this->s ; i++) {
			this->start.push_back( starts[i] );
			this->prob.push_back(0.0);
			this->prob.push_back(0.0);// no one would use this without 2 states
		}
	};

	void prepareIceCream(){
		// start has already been set!		
		this->start.reserve(10);
		for ( int i = 1; i < 11; i++ ){
			//Rcout << "ading ice start "<< i << " at id " << i-1 <<std::endl;
			this->start[i-1] = i;
		}

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
	int pos4val ( double val, int state ){
		if (  this->start.size() != 10 ){
			// first time usage...
			Rcout << "I currently have "<<this->start.size()<<
				" start positions - not the correct number of 10!" <<std::endl;
			for ( int i = 0; i < this->start.size(); i ++ ) {
				Rcout << this->start[i]<<" ";
			}
			Rcout <<std::endl;
			::Rf_error("object not initialized correctly - starts are empty / too many!? ");
		}

		for (int i =0; i< this->s; i++){
			if ( this->start[i] <= val & this->s == i +1 ){
				this->lastVal = val;
				this->lastID = i+ this->s * state;
				return i + this->s * state;
			}
			if ( this->start[i] <= val & this->start[i+1] > val ){
				this->lastVal = val;
				this->lastID = i+ this->s * state;
				return i + this->s * state;
			}
		}
	};

	void estimate ( std::vector<double> values, int state ) {

		//Rcout << "filling ProbEntry with size "<<  this->s << 
		//" and values size "<< values.size() << " for starte " << state <<std::endl;
		
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
		sum = tmp[0]+ static_cast<double>(tmp[this->s]);
		min= tmp[this->s]; // educated guess here -really
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
			tmp[i] = (tmp[i] - min) / sum  ; //we can not have 0 values later on!
		}


		//Rcout << "add the probabilites at "  << state << std::endl;
		for ( int i = 0; i < this->s; i++ ){
			this->prob.at(pos4val( this->start[i], state) ) = tmp[i];
		}
		//Rcout << "correct size?" << this->prob.size() <<  std::endl;
		//print();
		//Rcout << "survived!" << std::endl;
	};
	// this return value will be logged in the ProbList
	double Prob_4_value ( double val, int state) {
		return(this->prob[ pos4val( val, state ) ]);
	};
	void print( ) {
		
		Rcout << "A c++ ProbEntry object containing " << this->prob.size() / this->s << 
			" state(s) and " <<  this->start.size() << " start points" <<std::endl;
		
		Rcout << "start " << ": ";
		for ( int a=0; a < this->s; a ++) {
			Rcout << this->start[a] << " ";
		}
		Rcout << std::endl;

		for ( int state=0; state < this->prob.size() / this->s; state ++) {
			Rcout << "state "<< state << ": ";
			for ( int a =0; a< this->s; a++ ) {
				Rcout << this->prob[pos4val( this->start[a], state)] << " ";
			}
			Rcout << std::endl;
		}
	};

private:
	//these two get set in the fixOrder call:
	int lastID;
	double lastVal;

	void fixOrder() {
		this->lastID= 0;
		this->lastVal = (this->start[0] + this->start[1])/2;
	};

};


class TransmissionProb{
public:

	TransmissionProb() {};
	TransmissionProb( int states) {
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
	};
	void setStart ( std::vector<double> starts ) {
		if ( this->size != starts.size() ){
			 Rcout << "TransmissionProb not initialized correctly " << this->size << " != " << starts.size()<< std::endl;
			::Rf_error( "<- ERROR" );
		}
		this->startProbability = starts;
	};
	void setEnd ( std::vector<double> ends ) {
		this->endProbability = ends;
	};
	void setState ( int from, std::vector<double> pVals ) {
		for ( int to = 0; to < pVals.size(); to ++) {
			this->transitionProb[position( from, to)] = pVals[to];
		}
	};
	double getTransmissionProb( int from, int to){
		return( log(this->transitionProb[position(from, to)]));
	};
	double getStartProb( int state ) {
		return ( log( this->startProbability[state]));
	};
	double getEndProb( int state ) {
		return ( log(this->endProbability[state]));
	};
	void print() {
		Rcout << "A c++ TransmissionProb object containing:"  << std::endl;
		//Rcout << "Start probabilities: "  << std::ostream_iterator<double>(
		// this->startProbability, " ") << std::endl;
		//for ( int i =0; i < this->transitionProb.size() ){
		//	std::ostream_iterator<double>( this->transitionProb[i], " ") << std::endl;
		//}
		//Rcout << "End probabilities: "  << std::ostream_iterator<double>( 
	//		this->endProbability, " ") << std::endl;
	};
	void prepareIceCream(){
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


private:
	std::vector<double> startProbability;
	std::vector<double> endProbability;
	std::vector<double> transitionProb;
	int size;
	int position( int from, int to ) {
		int states;
		states = this->startProbability.size();
		return ( from + to * states);
	}
};



class ProbList {
public:

	std::vector<ProbEntry*> list;
	TransmissionProb* transmissionProb;

	//a vector of chars with max 10 entries...

	std::vector<std::string> states;


	ProbList () {};

	//method to fill the obejct with the ice cream data
	void prepareIceCream () {
		this->states.reserve(2);
		this->states.push_back("Hot");
		this->states.push_back("Cold");
		this->list.reserve( 33 );
		std::vector<double> range {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0} ;
		for ( int i =0; i< 33; i++) {
			this->list.push_back(new ProbEntry( range ));
			//Rcout << "ProbEntry add to object "<<i << " size: " << this->list.size() << std::endl;
			this->list[ i ]-> prepareIceCream();
		}
		//Rcout << "ProbEntry locations -> "; 
		//for ( int i = 0; i < 33; i++){
		//	Rcout << this->list[i] << " ";
		//}
		//Rcout << std::endl;
		
		this->transmissionProb = new TransmissionProb(2);
		this->transmissionProb->prepareIceCream();
		//print();
	}


	NumericMatrix as_Matrix ( ) {

		//Rcout << "creating out object of size" << 
		//this->list.size() * this->states.size() << " * 10" << std::endl;
		NumericMatrix ret ( this->list.size() * this->states.size(), 10 );
		std::vector<double> range = this->list.at(0)->start;
		int r_ = 0;
		for ( int i = 0; i < this->list.size(); i++){
			for ( int state = 0; state < this->states.size(); state++) {
				//Rcout << "result line " << r_ << " of " << this->list.size() * this->states.size()<< std::endl;
				for ( int id = 0; id < 10; id ++){
					//Rcout << "   result line "<< r_ << " and id "<< id <<std::endl;
					ret(r_, id) = this->list[i]->prob[ this->list[i]->pos4val( range[id], state)];
				}
				r_ ++;
			}
		}
		return (ret);
	};

	void estimate( Eigen::SparseMatrix<double> data, std::vector<double> range, 
		std::string name, bool phony ){
		this->states.reserve(this->states.size() +1 );
		this->states.push_back( name );

		if (this->list.size() < data.innerSize() ){
			// initialize the whole model
			this->list.reserve( data.innerSize() );
			for ( int i =this->list.size(); i < data.innerSize(); i++){
				this->list.push_back( new ProbEntry(range) ) ;
			}
		}
		
		std::vector<double> D (data.outerSize());
		if ( phony ) {
			Rcout << "estimating probabilities for state " << name << std::endl;
		}
		//Progress p(data.innerSize(), true);

		data = data.transpose();

		for ( int c_=0; c_ < data.outerSize(); ++c_ ){
			std::fill(D.begin(), D.end(), 0.0);
			for (Eigen::SparseMatrix<double>::InnerIterator it(data, c_); it; ++it){
				D[it.row()] =  it.value();
			}
			this->list[c_]->estimate( D, this->states.size()-1 );
		}

		data = data.transpose();
	};


	void estimate( Eigen::SparseMatrix<double> data, std::vector<int> cols,
		std::vector<double> range, std::string name, bool phony ){
		
		if ( this->states.size() == 0 ){
			this->list.reserve( data.innerSize() );
			//std::generate(this->list.begin(), this->list.end() , &ProbEntry(Range) );
			for( int i = 0; i < data.innerSize(); i ++){
				this->list.push_back( new ProbEntry( range ) ) ;
			}
		}
		//this->print();
		this->states.push_back( name );

		std::vector<double> D (data.outerSize());
		std::vector<double> A ( cols.size() );

		if ( phony ) {
			Rcout << "estimating probabilities for state " << name << " ("<< this->states.size() << ")" << std::endl;
		}
		/*
		if ( this->states.size() > 1 ){
			this->list.at(0)->print();
		}
		*/

		Progress p(data.innerSize(), phony);
		data = data.transpose();

		//Rcout << "just before the for loop " << std::endl;
		//this->print();

		for ( int c_=0; c_ < data.outerSize(); ++c_ ){
			std::fill(D.begin(), D.end(), 0.0);
			for (Eigen::SparseMatrix<double>::InnerIterator it(data, c_); it; ++it){
				D[it.row()] =  it.value();
			}
			for ( int i =0; i < cols.size(); i++){
				A[i] = D[cols[i]]; 
			}
			this->list.at(c_)->estimate( A, this->states.size()-1);
			try {
				//Rcout << "At position " << c_ << std::endl;
				//this->list.at(c_)->print();
				this->list.at(c_)->estimate( A, this->states.size()-1);
				//this->list.at(c_)->print();
  			}
  			catch (const std::out_of_range& oor) {
  				Rcout << "At position " << c_ << " you tried to use an uninitialized object" << std::endl;
  				::Rf_error("Creating a ProbEntry object?! NOT HERE!") ;
  				//Rcout << "object creation fincished" << std::endl;
 			}
 			//Rcout << "results prob functions length " << this->list[c_]->prob.size() << std::endl;
			p.increment();
		}

		//Rcout << "results prob functions length " << this->list[0]->prob.size() << std::endl;
		//Rcout << "state " << name << " finished" << std::endl;
	};
	double Prob_4_value ( int i, double val, int state) {
		return (log(list[i]->Prob_4_value( val, state )));
	};
	double Prob_4_value ( double val, int state) {
		return (log (list[0]->Prob_4_value( val, state )) );
	};
	void print (){
		Rcout << "A c++ ProbList object containing"  << std::endl;
		Rcout << "Transition probabitlities for states:"  << std::endl;
		for ( int i = 0; i < this->states.size(); i++ ) {
			Rcout << "'" <<this->states[i] <<  "' ";
		}
		Rcout << std::endl;
		this->transmissionProb->print();
		Rcout << "prob distributions (n="<< this->list.size()<<  "):" << std::endl;		
		for ( int i = 0; i < 5; i++) {
			try {
			Rcout << "Prop_state "<< i <<  ":" << "(" << this->list.at(i) << ")" <<std::endl;		
			this->list.at(i)->print();
			}
			catch (const std::out_of_range& oor) {
				break;
			}
		}
		for ( int i = 0; i < 5; i++) {
			Rcout << "      ."  << std::endl;		
		}
		for ( int i = this->list.size()-5; i < this->list.size(); i++) {
			try {
			Rcout << "Prop_state"<< i <<  ":" << "(" << this->list.at(i) << ")" << std::endl;		
			this->list.at(i)->print();
			}
			catch (const std::out_of_range& oor) {
				break;
			}
		}
	};

};


















