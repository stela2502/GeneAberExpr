#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

class MarcowChain{
public:
	std::vector<std::vector<double>> forwardResults;
	std::vector<std::vector<double>> backwardResults;
	std::vector<std::vector<double>> a_H;
	std::vector<double> C;
	

	NumericMatrix run ( Eigen::SparseMatrix<double> data, ProbList model ){
		/*$
		here we need to asigne each column a probability for A and for B

		marcovChain->CalculateForwardProbability();
		$marcovChain->CalculateBackwardProbability();
		$marcovChain->CalculateTotalProbabilityFromStartToEnd();
		*/

	}

	double _forwardAt(int i, double val, int state, ProbList model   ) {
		double res = 1;
		// sum (p from other to state * Prob_4_Value(state))
		double probStateI = model.Prob_4_value( i, val, state );
		if ( i == 1){
			res = model.transmissionProb.startProbability(state) + probStateI; // multiply
		}else {
			for( int i = 0; i < model.states.size(); i++ ){
				// multiply probabilities from all elements to this state
				res = Add2Logs( res, model.transmissionProb.getTransmissionProb( i, state) + probStateI);
			}
		}
		return (res);
	}

	void CalculateForwardProbability( std::vector<double> chain, ProbList model  ) {
		this->forwardResults.reserve( chain.size() );
		for ( in i =0; i < chain.size(); i++) {
			this->forwardResults[i].reserve(model.states.size());
		}
		int state;
		int from;
		double probStateI;
		for ( state = 0;  state < model.states.size(); state ++  ){
			this->forwardResults[0][state] =  
				model.StartProbability(state) + model.ProbList.Prob_4_value(0, chain[i], state );
		}
		
		for ( int i =1 ;i < chain.size(); i ++ ) {

			std::fill(this->forwardResults[i].begin(), this->forwardResults[i].end(), 1.0);
			for ( state = 0; state < model.states.seize(); state ++ ) {
				for ( from = 0; from < model.states.seize(); from ++ ) {
					this->forwardResults[i][state] = Add2Logs(
						this->forwardResults[i][state] , 
						model.transmissionProb.getTransmissionProb( from, state )  
						+ 
						this->forwardResults[i-1][from]
					);
				}
				this->forwardResults[i][state] +=  model.Prob_4_value( i, chain[i], state );
			}
		}	
	}

	void CalculateBackwardProbability( std::vector<double> chain, ProbList model  ) {
		
		this->backwardResults.reserve( chain.size() );
		int i;
		int state;
		for ( i = 0; i < chain.size(); i++) {
			this->backwardResults[i].reserve(model.states.size());
		}
		for ( state = 0; state <  model.states.size(); state ++) {
			this->backwardResults[state][ chain.size()-1 ] = model.getEndProb(state);
		}

		for ( int i = chain.size()-2; i > -1; i -- ) {

			std::fill(this->backwardResults[i].begin(), this->backwardResults[i].end(), 1.0);
			for ( state = 0; state <  model.states.size(); state ++) {
				for ( from = 0; from < model.states.seize(); from ++ ) {
					this->backwardResults[i][state] = Add2Logs(
						this->backwardResults[i][state] , 
						model.transmissionProb.getTransmissionProb( from, state )  
						+ 
						this->backwardResults[i+1][from]
					);
				}
				this->backwardResults[i][state] +=  model.Prob_4_value( i, chain[i], state );
			}
		}
 
	}

	void CalculateTotalProbabilityFromStartToEnd( std::vector<double> chain ) {

	}

	double Add2Logs (double logA, double logB ) {
		double ret =  logA + log( 1 + exp( logB - logA ) );
		return (ret);
	}

}