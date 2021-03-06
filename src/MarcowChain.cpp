#include "MarcowChain.h"

MarcowChain::MarcowChain() {
	Rcpp::stop("Use the other constructor with chainLength and states..");
}

MarcowChain::MarcowChain( int chainLength, int hiddenStates ) {
	setup( chainLength, hiddenStates);
	//Rcout << "prepare finished: " << this->chainLength << std::endl;
}

void MarcowChain::setup( int chainLength, int hiddenStates ){
	int total = chainLength*hiddenStates;
	//Rcout << "prepare Object with long and short vectors: " << total << " ;"<< 
	//  chainLength << std::endl;
	this->chainLength = chainLength;
	this->forwardResults.resize( total );
	this->backwardResults.resize( total );
	this->a_H.resize( total );
	this->p_H.resize( total );
	this->C.resize( chainLength );
	double null = 0.0;
	std::fill(this->forwardResults.begin(), this->forwardResults.end(), null);
	std::fill(this->backwardResults.begin(), this->backwardResults.end(), null);
	std::fill(this->a_H.begin(), this->a_H.end(), null);
	std::fill(this->p_H.begin(), this->p_H.end(), null);
	std::fill(this->C.begin(), this->C.end(), null);
}

int MarcowChain::at( int id, int state ) {
	//Rcout << "at("<< id<<", "<<state<<") = " << state <<" * "<< 
	 //  this->chainLength<< " + " << id << std::endl;
	return state * this->chainLength + id;
}

NumericMatrix MarcowChain::run ( std::vector<double> chain, ProbList* model, bool test ){
	
	//here we need to asigne each column a probability for A and for B
	if ( this->chainLength != chain.size() ) {
		Rcout << "setup was not correct?! " << this->chainLength << std::endl;
		setup(chain.size(), model->states.size() );
		//::Rf_error("Not initialized correctly!") ;
	}

	Rcout  << std::endl << "CalculateForwardProbability " << std::endl;
	CalculateForwardProbability(chain, model);
	Rcout << "CalculateBackwardProbability " << std::endl;
	CalculateBackwardProbability(chain, model);
	Rcout << "CalculateTotalProbabilityFromStartToEnd " << std::endl;
	CalculateTotalProbabilityFromStartToEnd(chain, model);
	NumericMatrix ret;
	if ( test ) {
		ret = to_matrix( model );
	}else {
		ret = results( model );
	}
	
	return ( ret );
}

NumericMatrix MarcowChain::results( ProbList* model ){

	int columns = model->states.size();

	//Rcout << "prepare results " << std::endl;
	NumericMatrix ret ( this->chainLength, columns );
	std::fill( ret.begin(), ret.end(), NumericVector::get_na() ) ;
	int id;
	int c_;
	for ( int r_ = 0; r_ < this->chainLength; r_++){
		c_ = 0;
		for( int i = 0; i < model->states.size() ; i++ )  {
			id = at( r_, i );
			//ret(r_, c_++) = this->forwardResults.at(id);
			//ret(r_, c_++) = this->backwardResults.at(id);
			//ret(r_, c_++) = this->a_H.at(id);
			ret(r_, c_++) = exp(this->p_H.at(id));
		}
		//ret(r_, c_++) = this->C.at(r_);
	}

	Rcpp::StringVector StringV2(model->states.size());
	StringV2 = model->states;
	colnames( ret ) = StringV2;
	//Rcout << "finished!" << std::endl;
	//colnames(ret) = CharacterVector::create (cn);
	return ( ret );

}

NumericMatrix MarcowChain::to_matrix( ProbList* model ){

	int states = model->states.size();
	int columns = model->states.size() * 4 + 1;

	//Rcout << "prepare results " << std::endl;
	NumericMatrix ret ( this->chainLength, columns );
	std::fill( ret.begin(), ret.end(), NumericVector::get_na() ) ;
	int id;
	int c_;
	for ( int r_ = 0; r_ < this->chainLength; r_++){
		c_ = 0;
		for( int i = 0; i < model->states.size() ; i++ )  {
			id = at( r_, i );
			ret(r_, c_++) = this->forwardResults.at(id);
			ret(r_, c_++) = this->backwardResults.at(id);
			ret(r_, c_++) = this->a_H.at(id);
			ret(r_, c_++) = this->p_H.at(id);
		}
		ret(r_, c_++) = this->C.at(r_);
	}

	//Rcout << "finished!" << std::endl;
	//colnames(ret) = CharacterVector::create (cn);
	return ( ret );

}

void MarcowChain::CalculateForwardProbability( std::vector<double> chain, ProbList* model  ) {
	if ( this->forwardResults.size() < chain.size() * model->states.size() ){
		setup( chain.size(), model->states.size() );
	}
	
	int state;
	int from;

	for ( state = 0;  state < model->states.size(); state ++  ){
		// the point 0 needs to be set as log(0) is undefined!
		/*Rcout << "forward state start 0:"<<state<<" tying to add to id "<< 
			at(0,state) << " : " << model->transmissionProb->getStartProb(state) << " + "<< 
			model->Prob_4_value(0, chain[0], state ) << std::endl;*/
		this->forwardResults.at(at(0,state)) =  
			model->transmissionProb->getStartProb(state) + 
				model->Prob_4_value(0, chain[0], state );
	}
	
	for ( int i =1 ;i < chain.size(); i ++ ) {
		//this forward C = (Last forward Prob C * p(C|C) + last forward Prob H * p(C|H) )
		//                 * Prob_4_value( i, chain[i], 'C' )

		for ( state = 0; state < model->states.size(); state ++ ) {
			/*Rcout << "forward state "<< state << " = (" <<
			exp(model->transmissionProb->getTransmissionProb( 0, state )) <<
			" * " << exp(this->forwardResults.at(at(i-1, 0))) << " + " << 
			exp(model->transmissionProb->getTransmissionProb( 1, state )) <<
			" * " << exp(this->forwardResults.at(at(i-1, 1))) << ") * "<<
			exp( model->Prob_4_value( i, chain[i], state ));*/

			// the first calculation without add needs to be kept separate (log(0)!)
			this->forwardResults.at(at(i, state)) = 
				model->transmissionProb->getTransmissionProb( 0, state )  +
				this->forwardResults.at(at(i-1, 0));
			// calculate the rest for to->state!
			for ( from = 1; from < model->states.size(); from ++ ) {
				this->forwardResults.at(at(i,state)) = Add2Logs(
					this->forwardResults.at(at(i,state)) , 
					model->transmissionProb->getTransmissionProb( from, state )  
					+ 
					this->forwardResults.at(at(i-1, from))
				);
			}
			this->forwardResults.at(at(i,state)) +=  
				model->Prob_4_value( i, chain[i], state );
			//Rcout << "forwardResults at () "<< i<< ", "<< state<<"="<< exp(this->forwardResults.at(at(i,state)))<<std::endl;
		} 
	}
}

void MarcowChain::CalculateBackwardProbability( std::vector<double> chain, ProbList* model  ) {
	
	//this->backwardResults.resize( chain.size() * model->states.size() );
	int from;
	int state;

	for ( state = 0; state <  model->states.size(); state ++) {
		this->backwardResults.at(at(chain.size()-1, state)) = 
			model->transmissionProb->getEndProb(state);
	}

	for ( int i = chain.size()-2; i > -1; i -- ) {
		// =C$14*E29*INDEX(C$11:C$13;$B29;1)+C$15*F29*INDEX(D$11:D$13;$B29;1)
		// this backward C = ( 
		//        p(C|C) * last backward Prob C * Prob_4_value( i+1, chain[i+1], 'C' ) 
		//          +
		//        p(H|C) * last backward Prob H * Prob_4_value( i+1, chain[i+1], 'H' ))
		for ( state = 0; state <  model->states.size(); state ++) {
			/*Rcout << "backward state "<< state << " = (" <<
			exp(model->transmissionProb->getTransmissionProb( 0, state )) <<
			" * " << exp(this->backwardResults.at(at(i+1, 0))) << " * " << 
			exp(model->Prob_4_value( i+1, chain[i+1], 0 )) <<
			" + " << exp(model->transmissionProb->getTransmissionProb( 0, 1 ))<<
			exp(this->backwardResults.at(at(i+1, 1))) << " * "<<
			exp( model->Prob_4_value( i+1, chain[i], 1 ));*/

			this->backwardResults.at(at(i,state)) = 
				model->transmissionProb->getTransmissionProb( 0, state )  +
				this->backwardResults.at(at(i+1, 0)) + model->Prob_4_value( i+1, chain[i+1], 0 );
			for ( from = 1; from < model->states.size(); from ++ ) {
				this->backwardResults.at(at(i,state)) = Add2Logs(
					this->backwardResults.at(at(i, state)) , 
					model->transmissionProb->getTransmissionProb( from, state )  +
				this->backwardResults.at(at(i+1, from)) + model->Prob_4_value( i+1, chain[i+1], from )
				);
			}
			//Rcout << " = "<< exp(this->backwardResults.at(at(i,state)))<<std::endl;
		}
	}

}

void MarcowChain::CalculateTotalProbabilityFromStartToEnd( std::vector<double> chain,  ProbList* model  ) {
	//this->a_H.resize( chain.size() * model->states.size() );
	//this->p_H.resize( chain.size() * model->states.size() );
	//this->C.resize( chain.size());
	int state;
	for ( int i =0; i < chain.size(); i++ ){
		int id;
		
		for ( state = 0; state < model->states.size(); state ++) {
			id =at(i,state);
			this->a_H[id] = this->forwardResults.at(id) + 
			  this->backwardResults.at(id);
		}
		this->C[i] = this->a_H.at(at(i,0));
		for (  state = 1; state < model->states.size(); state ++) {
			this->C.at(i) = Add2Logs(this->C[i], this->a_H.at(at(i, state)) );
		}
		//scale results
		for ( state = 0; state < model->states.size(); state ++) {
			this->p_H.at(at(i, state)) = this->a_H[at(i, state)] - this->C.at(i);
		}
	}
}

double MarcowChain::Add2Logs (double logA, double logB ) {
	double ret =  logA + log( 1 + exp( logB - logA ) );
	return (ret);
}