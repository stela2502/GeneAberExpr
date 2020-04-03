
#include "ProbList.h"


ProbList::ProbList (int marcowLength, int states, std::vector<double> * s ){
	setup(marcowLength, states, s );
}

void ProbList::setup( int marcowLength, int states, std::vector<double> * s ){
	if ( this->list.size() > 0 ){
		Rcpp::stop("ProbList must not be set up twice!" );
	}
	this->starts = s;
	this->states.reserve(states);
	this->list.reserve(marcowLength);

	//Rcout << "ProbList::ProbList( int, int): [end]" << std::endl;
	this->transmissionProb = new TransmissionProb( states );

	//Rcout << "ProbList::ProbList( int, int): [-1]" << std::endl;
	for ( int i = 0; i < marcowLength; i ++ ) {
		//Rcout << "ProbList::ProbList( int, int): ["<<i<<"]" << std::endl;
		ProbEntry * tmp = new ProbEntry( s, states);
		this->list.push_back(tmp);
	}
}


void ProbList::prepareIceCream (bool phony) {

	if ( phony) {
		Rcout << "ProbList::prepareIceCream:" << std::endl;
	}
	std::vector<double> range(10);
	std::iota(range.begin(), range.end(), 1.0);

	setup( 33, 2, &range );

	if ( phony) {
		Rcout << "Setup finished:" << std::endl;
	}
	this->states.push_back ("Hot");
	this->states.push_back ("Cold");

	
	this->transmissionProb->prepareIceCream(phony);
	
	if ( phony) {	
		Rcout << "Initialize the ProbEntries list :" << std::endl;
	}
	//this->starts = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
	for ( int i =0; i< 33; i++) {
		this->list.at(i)-> prepareIceCream(phony);
		//Rcout << "ProbEntry add to object "<<i << "          ID: "<< this->list.at(i) << std::endl;
		//this->list.at(i)->print();
	}
	if ( phony) {
		Rcout << "ProbList list ptm's: ";
		for ( int i =0; i < 33; i++) {
			Rcout << this->list.at(i) << " ";
		}
		Rcout << "ProbList::prepareIceCream <model>:" << std::endl;
		print();
	}
}


//NumericMatrix as_Matrix ( ) {
NumericMatrix ProbList::as_Matrix ( ) {

		//Rcout << "creating out object of size" << 
		//this->list.size() * this->states.size() << " * 10" << std::endl;
	NumericMatrix ret ( this->list.size() * this->states.size(), this->list[0]->s);

	int r_ = 0;
	
	CharacterVector cn ( this->list[0]->s );
	CharacterVector rn ( this->list.size() * this->states.size() );

	for ( int i = 0; i < this->list.size(); i++){
		for ( int state = 0; state < this->states.size(); state++) {
				//Rcout << "result line " << r_ << " of " << this->list.size() * this->states.size()<< std::endl;
			for ( int id = 0; id < this->list[0]->s; id ++){
					//Rcout << "   result line "<< r_ << " and id "<< id <<std::endl;
				ret(r_, id) = Prob_4_value(i, this->starts->at(id), state);
			}
			rn(r_) = std::to_string(i) + this->states[state];
			r_ ++;
		}
	}
	for ( int id = 0; id < this->list[0]->s; id++){
		cn(id) = std::to_string(this->starts->at(id));
	}
	colnames( ret ) = cn;
	rownames( ret ) = rn;
	return (ret);
}

//void estimate( Eigen::SparseMatrix<double> data, std::vector<double> range, 
//	std::string name, bool phony ){
void ProbList::estimate( Eigen::SparseMatrix<double> data, std::vector<double> range, 
	std::string name, bool phony ){

	try {this->states.push_back( name );}
	catch( const std::exception& e ){
		Rcpp::stop("The ProbeList has not been initialized to contain that much information!");
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
}


void ProbList::estimate( Eigen::SparseMatrix<double> data, std::vector<int> cols,
	std::vector<double> range, std::string name, bool phony ){

	try {this->states.push_back( name );}
	catch( const std::exception& e ){
		Rcpp::stop("The ProbeList has not been initialized to contain that much information!");
	}
	

	std::vector<double> D (data.outerSize());
	std::vector<double> A ( cols.size() );

	if ( phony ) {
		Rcout << "estimating probabilities for state " << name << " ("<< this->states.size() << ")" << std::endl;
	}

	Progress p(data.innerSize(), phony);
	data = data.transpose();


	for ( int c_=0; c_ < data.outerSize(); ++c_ ){
		std::fill(D.begin(), D.end(), 0.0);
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, c_); it; ++it){
			D[it.row()] =  it.value();
		}
		for ( int i =0; i < cols.size(); i++){
			A[i] = D[cols[i]]; 
		}
		try {
			this->list.at(c_)->estimate( A, this->states.size()-1);
		}
		catch (const std::out_of_range& oor) {
			Rcpp::stop("At position " +std::to_string(c_)+" you tried to use an uninitialized object" ) ;
		}
	
		p.increment();
	}
}

double ProbList::Prob_4_value ( int i, double val, int state) {
	try{
		return (log(this->list.at(i)->Prob_4_value_e( val, state )));
	}
	catch( const std::exception& e ){
		Rcpp::stop("exception: list.at(i="+ std::to_string(i) + ") Prob_4_value_e: " + e.what() );
	}
}
double ProbList::Prob_4_value ( double val, int state) {
	return (log (this->list.at(0)->Prob_4_value_e( val, state )) );
}
void ProbList::print (){
	Rcout << "A c++ ProbList object containing"  << std::endl;
	Rcout << "Transition probabitlities for states:"  << std::endl;
	for ( int i = 0; i < this->states.size(); i++ ) {
		Rcout << "'" <<this->states[i] <<  "' ";
	}
	Rcout << std::endl;
	this->transmissionProb->print( this->states );
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
}