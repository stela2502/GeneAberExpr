
#include "../inst/include/GeneAberExpr.h"


void ProbList::prepareIceCream () {
	this->states.reserve(2);
	this->states.push_back("Hot");
	this->states.push_back("Cold");
	this->list.reserve( 33 );
	std::vector<double>  range({ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0});
	this->starts = range;
	//std::iota(range.begin(), range.end(), 1);
	//this->starts = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
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
				ret(r_, id) = this->list[i]->prob[
				 this->list[i]->pos4val( this->starts[id], state)];
			}
			rn(r_) = std::to_string(i) + this->states[state];
			r_ ++;
		}
	}
	for ( int id = 0; id < this->list[0]->s; id++){
		cn(id) = std::to_string(this->starts[id]);
	}
	colnames( ret ) = cn;
	rownames( ret ) = rn;
	return (ret);
}

//void estimate( Eigen::SparseMatrix<double> data, std::vector<double> range, 
//	std::string name, bool phony ){
void ProbList::estimate( Eigen::SparseMatrix<double> data, std::vector<double> range, 
	std::string name, bool phony ){

	this->states.reserve(this->states.size() +1 );
	this->states.push_back( name );
	

	if (this->list.size() < data.innerSize() ){
		this->starts = range;
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
}


void ProbList::estimate( Eigen::SparseMatrix<double> data, std::vector<int> cols,
	std::vector<double> range, std::string name, bool phony ){

	if ( this->states.size() == 0 ){
		this->starts = range;
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
}
double ProbList::Prob_4_value ( int i, double val, int state) {
	return (log(list[i]->Prob_4_value( val, state )));
}
double ProbList::Prob_4_value ( double val, int state) {
	return (log (list[0]->Prob_4_value( val, state )) );
}
void ProbList::print (){
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
}