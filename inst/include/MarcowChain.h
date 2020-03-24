#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>
#include <iterator>

using namespace Rcpp;

class MarcowChain{
public:
	std::vector<double> forwardResults;
	std::vector<double> backwardResults;
	std::vector<double> a_H;
	std::vector<double> p_H;
	std::vector<double> C;
	int chainLength;

	MarcowChain() {	};

	MarcowChain( int chainLength, int hiddenStates ) {
		int total = chainLength*hiddenStates;
		this->forwardResults.reserve( total );
		this->backwardResults.reserve( total );
		this->a_H.reserve( total );
		this->p_H.reserve( total );
		this->C.reserve( chainLength );
		this->chainLength = chainLength;
	};

	int at( int id, int state ) {
		return state * this->chainLength + id;
	};

	NumericMatrix run ( std::vector<double> chain, ProbList* model ){
		
		//here we need to asigne each column a probability for A and for B

		CalculateForwardProbability(chain, model);
		CalculateBackwardProbability(chain, model);
		CalculateTotalProbabilityFromStartToEnd(chain, model);


		int states = model->states.size();
		int columns = model->states.size() * states + 1;

		/*
		std::vector<char &> cn(columns);

		for( int i = 0; i < states ; i++ )  {
			std::string charV = model->states[i] + " forward";
			cn[i*states] = charV;
			std::string charV = model->states[i] + " backward";
			cn[i*states+1] = charV;
			std::string charV = model->states[i] + " a_H";
			cn[i*states+2] = charV;
			std::string charV = model->states[i] + " p_H";
			cn[i*states+3] = charV;
		}
		cn[columns] = "C";
		*/

		NumericMatrix ret ( this->chainLength, columns );
		int id;
		for ( int r_ = 0; r_ < this->chainLength; r_++){
			for( int i = 0; i < model->states.size() ; i++ )  {
				id = at( r_, i );
				ret(r_, i) = this->forwardResults[id];
				ret(r_, i*states+1) = this->backwardResults[id];
				ret(r_, i*states+2) = this->a_H[id];
				ret(r_, i*states+3) = this->p_H[id];
			}
			ret(r_, columns) = this->C[r_];
		}


		//colnames(ret) = CharacterVector::create (cn);
		return ( ret );

	};

	double _forwardAt(int i, double val, int state, ProbList* model   ) {
		double res = 1;
		// sum (p from other to state * Prob_4_Value(state))
		double probStateI = model->Prob_4_value( i, val, state );
		if ( i == 1){
			res = model->transmissionProb->getStartProb(state) + probStateI; // multiply
		}else {
			for( int i = 0; i < model->states.size(); i++ ){
				// multiply probabilities from all elements to this state
				res = Add2Logs( res, model->transmissionProb->getTransmissionProb( i, state) + probStateI);
			}
		}
		return (res);
	};

	void CalculateForwardProbability( std::vector<double> chain, ProbList* model  ) {
		this->forwardResults.reserve( chain.size() * model->states.size() );
		
		int state;
		int from;
		double probStateI;
		for ( state = 0;  state < model->states.size(); state ++  ){
			// the point 0 needs to be set as log(0) is undefined!
			this->forwardResults[at(0,state)] =  
				model->transmissionProb->getStartProb(state) + 
					model->Prob_4_value(0, chain[0], state );
		}
		
		for ( int i =1 ;i < chain.size(); i ++ ) {
			for ( state = 0; state < model->states.size(); state ++ ) {
				// the forst calculation without add needs to be kept separate (log(0)!)
				this->forwardResults[at(i, state)] = 
					model->transmissionProb->getTransmissionProb( 0, state )  +
					this->forwardResults[0, state];
				// calculate the rest!
				for ( from = 1; from < model->states.size(); from ++ ) {
					this->forwardResults[at(i,state)] = Add2Logs(
						this->forwardResults[at(i,state)] , 
						model->transmissionProb->getTransmissionProb( from, state )  
						+ 
						this->forwardResults[at(i-1, from)]
					);
				}
				this->forwardResults[at(i,state)] +=  
					model->Prob_4_value( i, chain[i], state );
			}
		}	
	};

	void CalculateBackwardProbability( std::vector<double> chain, ProbList* model  ) {
		
		this->backwardResults.reserve( chain.size() * model->states.size() );
		int i;
		int from;
		int state;

		for ( state = 0; state <  model->states.size(); state ++) {
			this->backwardResults[at(chain.size()-1, state)] = 
				model->transmissionProb->getEndProb(state);
		}

		for ( int i = chain.size()-2; i > -1; i -- ) {
			for ( state = 0; state <  model->states.size(); state ++) {
				this->backwardResults[at(i,state)] = 
					model->transmissionProb->getTransmissionProb( 0, state )  +
					this->backwardResults[at(i+1,0)];
				for ( from = 1; from < model->states.size(); from ++ ) {
					this->backwardResults[at(i,state)] = Add2Logs(
						this->backwardResults[at(i, state)] , 
						model->transmissionProb->getTransmissionProb( from, state )  
						+ 
						this->backwardResults[at(i+1,from)]
					);
				}
				this->backwardResults[at(i,state)] +=  
					model->Prob_4_value( i, chain[i], state );
			}
		}
 
	};

	void CalculateTotalProbabilityFromStartToEnd( std::vector<double> chain,  ProbList* model  ) {
		this->a_H.reserve( chain.size() * model->states.size() );
		this->p_H.reserve( chain.size() * model->states.size() );
		this->C.reserve( chain.size());
		int i;
		int state;
		int a;

		for ( int i =0; i < chain.size(); i++ ){
			int id;
			
			for ( state = 0; state < model->states.size(); state ++) {
				id =at(i,state);
				this->a_H[id] = this->forwardResults[id] + 
				  this->backwardResults[id];
			}
			this->C[i] = this->a_H[at(i,0)];
			for (  state = 1; state < model->states.size(); state ++) {
				this->C[i] = Add2Logs(this->C[i], this->a_H[at(i, state)]);
			}
			//scale results
			for ( state = 0; state < model->states.size(); state ++) {
				this->p_H[at(i, state)] -= this->C[i];
			}
		}
	};

	double Add2Logs (double logA, double logB ) {
		double ret =  logA + log( 1 + exp( logB - logA ) );
		return (ret);
	};

};