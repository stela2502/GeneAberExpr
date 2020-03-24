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


	int pos4val ( double val, int state ){
		if (  this->start.size() != 10 ){
			// first time usage...
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
		
		this->prob.reserve( state * this->s );
		for ( int i = this->prob.size(); i <state * this->s; i ++ ){
			this->prob.push_back(0.0);
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
		sum = tmp[0]+ static_cast<double>(tmp[this->s]);
		min= tmp[this->s]; // educated guess here -really
		for( i = 1; i < this->s -1; i ++) {
			tmp[i] = (tmp[i] + tmp[i-1] + tmp[i+1]) / 3;
			sum += tmp[i];
			if ( min > tmp[i] ){
				min = tmp[i];
			}
		}
		//Rcout << "rescaled and smoothed the table"  << std::endl;
		for ( i = 0; i < this->s; i++ ){
			//Rcout << "(prob[i] - min)" << (prob[i] - min) << " / sum " <<  sum << "+ 1e-9" << std::endl;
			tmp[i] = (tmp[i] - min) / sum + 1e-9 ; //we can not have 0 values later on!
		}


		//Rcout << "add the probabilites at "  << state << std::endl;
		for ( int i = 0; i < this->s; i++ ){
			//try{
			this->prob.at(pos4val( this->start[i], state) ) = tmp[i];
  			/*
  			}
  			catch (const std::out_of_range& oor) {
	  			this->prob.push_back(tmp[i]);
	  			Rcout << "prob.push_back(tmp[i]) "  << tmp[i] << " at i " << 
	  			i << "( size="<<this->prob.size() << "; position=" << 
	  			pos4val( this->start[i], state) <<")" <<std::endl;
			}
			*/
		}
		//Rcout << "correct size?" << this->prob.size() <<  std::endl;
		//print();
		//Rcout << "survived!" << std::endl;
	};
	// this return value will be logged in the ProbList
	double Prob_4_value ( double val, int state) {
		if ( this->lastVal == val ){
			return ( this->prob[ this->lastID ] );
		}
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
	std::vector<double> tmp;
	//these two get set in the fixOrder call:
	int lastID;
	double lastVal;

	void fixOrder() {
		/*
		std::vector<paired> pairs;
		pairs.reserve(this->s);
		int i;

		for ( i = 0; i < this->s; i++ ){
			double sum = 0.0;
			for ( int state =0; state< this->prob.size(); state++){
				sum += this->prob[state][i];
			}
			pairs.push_back(std::make_pair(i, sum) );
		}		

		std::sort(pairs.begin(), pairs.end(), [](const std::pair<int,int> &left, const std::pair<int,int> &right) {
		    return left.second > right.second;
			}
		);

		Rcout << "new order: ";

		for ( i = 0; i < this->s; i++ ){
			try {
				this->order.at(i) = pairs[i].first;
  			}
  			catch (const std::out_of_range& oor) {
  				this->order.push_back(pairs[i].first);
			}
			Rcout << this->order[i]<< " " ;

		}
		Rcout << std::endl;
		*/
		this->lastID= 0;
		this->lastVal = (this->start[0] + this->start[1])/2;
	};

};


class TransmissionProb{
public:

	TransmissionProb() {};
	TransmissionProb( int states) {
		this->startProbability.resize(states);
		this->endProbability.resize(states);
		this->transitionProb.resize(states* states);
		std::fill(this->startProbability.begin(), this->startProbability.end(), 0.0 );
		std::fill(this->endProbability.begin(), this->endProbability.end(), 0.0 );
		std::fill(this->transitionProb.begin(), this->transitionProb.end(), 0.0 );
	};
	void setStart ( std::vector<double> starts ) {
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
	}

private:
	std::vector<double> startProbability;
	std::vector<double> endProbability;
	std::vector<double> transitionProb;
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


	ProbList () {
	};

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
					ret(r_, id) = this->list[i]->prob[ this->list[i]->pos4val( range[id], state)];
				}
				r_ ++;
			}
		}
		//Rcout << "finished" << std::endl;
		return (ret);
	};

	void estimate( Eigen::SparseMatrix<double> data, std::vector<double> range, 
		std::string name ){
		this->states.push_back( name );

		this->list.reserve( data.innerSize() );
		std::vector<double> D (data.outerSize());

		Rcout << "estimating probabilities for state " << name << std::endl;
		//Progress p(data.innerSize(), true);

		data = data.transpose();

		for ( int c_=0; c_ < data.outerSize(); ++c_ ){
			std::fill(D.begin(), D.end(), 0.0);
			for (Eigen::SparseMatrix<double>::InnerIterator it(data, c_); it; ++it){
				D[it.row()] =  it.value();
			}
			ProbEntry tmp = ProbEntry( range );
			this->list.push_back( &tmp ) ;
			this->list[c_]->estimate( D, this->states.size()-1 );
			
		}

		data = data.transpose();
	};


	void estimate( Eigen::SparseMatrix<double> data, std::vector<int> cols,
		std::vector<double> range, std::string name ){
		
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

		Rcout << "estimating probabilities for state " << name << " ("<< this->states.size() << ")" << std::endl;
		
		/*
		if ( this->states.size() > 1 ){
			this->list.at(0)->print();
		}
		*/

		Progress p(data.innerSize(), true);
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
			//this->list.at(c_)->estimate( A, this->states.size()-1);
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


















