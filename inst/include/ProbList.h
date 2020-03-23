#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <math.h>

class ProbEntry{
public:
	std::vector<std::vector<double>> prob;
	static int size = 10;
	std::vector<double> start (size);
	std::vector<int>    order (size); // in which order should the start be evaluated?

	typedef std::pair<int, double> paired;

	ProbEntry () { };

	ProbEntry (std::vector<double> starts ) {
		if ( starts.size() != this->size ) {
			::Rf_error("the start sizes need to be exactly 10!" );
		}
		this->start = starts;
		for ( int i = 0; i <starts.size() +1; i++ ){
			order[i] = i;
		}
	};

	int pos4val ( double val ){
		for (int i =0; i< this->order.size(); i++){
			if ( this->start[i] <= val & this->order.size() == i+1 ){
				this->lastVal = val;
				this->lastID = i;
				return i;
			}
			if ( this->start[i] <= val & this->start[i+1] > val ){
				this->lastVal = val;
				this->lastID = i;
				return i;
			}
		}
	};

	std::vector<double> estimate ( std::vector<double> values, int state ) {
		std::vector<double> tmp (this->size);
		std::vector<double> prob (this->size);

		std::fill(tmp.begin(), tmp.end(), 0.0 );
		std::fill(prob.begin(), prob.end(), 0.0 );
		<double> sum = 0.0;
		int i;
		// fill an initial table
		for( i = 0; i < values.size(); i ++) {
			tmp[ pos4val(values[i]) ] ++;
		}
		sum = std::accumulate(tmp.begin(), tmp.end(), 0.0);
		//calculate raw prob
		for( i = 0; i < tmp.size(); i ++) {
			prob[i] = tmp[i] / sum;
		}

		// smooth this
		double min = 100;
		sum = 0.0;
		for( i = 1; i < tmp.size() -1; i ++) {
			prob[i] = (prob[i] + prob[i-1] + prob[i+1]) / 3;
			sum += prob[i];
			if ( prob[i] < min){
				min = prob[i];
			}
		}

		for ( i = 0; i < tmp.size(); i++ ){
			tmp[i] = log((prob[i] - min) / sum + 1e-9) ; //we can not have 0 values later on!
		}
		if ( this->prob.size() < state ){
			this->prob.resize(state);
		}
		this->prob[state] = tmp;
		fixOrder();
	}
	// this return value will be logged in the ProbList
	double Prob_4_value ( double val, int state) {
		if ( this->lastVal == val ){
			return ( this->proA[ this->lastID ] );
		}
		return(this->pro[state][ pos4val( val ) ]);
	};
	void print( ) {
		
		Rcout << "A c++ ProbEntry object containing:"  << std::endl;

	}

private:
	std::vector<double> tmp;
	//these two get set in the fixOrder call:
	int lastID;
	double lastVal;

	bool cmp_second(const paired & left, const paired & right) {
    	return left.second < right.second;
	};


	void fixOrder() {
		std::vector<paired> pairs;
		pairs.reserve(this->size);
		for ( i = 0; i < this->size; i++ ){
			double sum = 0.0;
			for ( int state =0; state< this->prob.size(); state++){
				sum += this->prob[state][i];
			}
			pairs.push_back(std::make_pair(i, sum) );
		}		
		std::sort(pairs.begin(), pairs.end(), cmp_second<paired>);

		for ( i = 0; i < tmp.size(); i++ ){
			this->order[i] = pairs[i].first;
		}
		this->lastID= 0;
		this->lastVal = (this->start[0] + this->start[1])/2;
	}

};

class ProbList {
public:

	std::vector<ProbEntry> list;
	TransmissionProb transmissionProb;

	std::vector<string> states;

	//a vector of chars with max 10 entries...
	std::vector<std::string> states;

	ProbList () {};


	void estimate( Eigen::SparseMatrix<double> data, std::vector<double> range, std::string name ){
		this->states.push_back( name );

		this->list.reserve( data.innerSize() );
		std::vector<double> D (data.outerSize());

		Rcout << "estimating probabilities for state " << state << std::endl;
		Progress p(data.innerSize(), true);

		data = data.transpose();

		for ( int c_=0; c_ < data.outerSize(); ++c_ ){
			std::fill(D.begin(), D.end(), 0.0);
			for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
				D[it.row()] =  it.value();
			}
			list[i]=ProbEntry( range ) ;
			this->list[c_].estimate( D, this->states.size() );
			p.increment();
		}
	}

	void estimate( Eigen::SparseMatrix<double> data, std::vector<int> cols,
		std::vector<double> range, std::string name ){
		this->states.push_back( name );

		this->list.reserve( data.innerSize() );
		std::vector<double> D (data.outerSize());
		std::vector<double> A ( col.size() );

		Rcout << "estimating probabilities for state " << state << std::endl;
		Progress p(data.innerSize(), true);
		data = data.transpose();

		for ( int c_=0; c_ < data.outerSize(); ++c_ ){
			std::fill(D.begin(), D.end(), 0.0);
			for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
				D[it.row()] =  it.value();
			}
			for ( int i =0; i < col.size(); i++){
				A[i] = D[col[i]]; 
			}
			list[i]=ProbEntry( range ) ;
			this->list[c_].estimate( A, this->states.size() );
			p.increment();
		}
	}
	double Prob_4_value ( int i, double val, int state) {
		return (log(list[i].Prob_4_value( val, state )));
	};
	double Prob_4_value ( double val, int state) {
		return (log (list[0].Prob_4_value( val, state )) );
	};
	void print (){
		Rcout << "A c++ ProbList object containing"  << std::endl;
		Rcout << "Transition probabitlities for states:"  << std::endl;
		Rcout << std::ostream_iterator<double>( this->states, " ") << std::endl;
		this->transmissionProb.print();
		Rcout << "prob distributions (n="<< this->list.size()<<  "):" << std::endl;		
		for ( int i = 0; i < 5; i++) {
			this->list[i].print(i);
		}
		for ( int i = 0; i < 5; i++) {
			Rcout << "   ." << std::endl;
		}
		for ( int i = this->list.size()-5; i < this->list.size(); i++) {
			this->list[i].print(i);
		}
	}
};

class TransmissionProb{
public:

	TransmissionProb() {};
	TransmissionProb( int states) {
		this.startProbability.resize(states);
		this.endProbability.resize(states);
		this.transitionProb.resize(states);
		for ( int i = 0; i < states; i++){
			this.transitionProb[i].resize(size);
		}
	};
	void setStart ( std::vector<double> starts ) {
		this->startProbability = starts;
	};
	void setEnd ( std::vector<double> ends ) {
		this->endProbability = ends;
	};
	void setState ( int from, std::vector<double> pVals ) {
		this.transitionProb[from] = pVals;
	};
	double getTransmissionProb( int from, int to){
		return( log(this.transitionProb[from][to]));
	};
	double getStartProb( int state ) {
		return ( log( this.startProbability[state]));
	};
	double getEndProb( int state ) {
		return ( log(this.endProbability[state]));
	};
	void print() {
		Rcout << "A c++ TransmissionProb object containing:"  << std::endl;
		Rcout << "Start probabilities: "  << std::ostream_iterator<double>( this->startProbability, " ") << std::endl;
		for ( int i =0; i < this->transitionProb.size() ){
			std::ostream_iterator<double>( this->transitionProb[i], " ") << std::endl;
		}
		Rcout << "End probabilities: "  << std::ostream_iterator<double>( this->endProbability, " ") << std::endl;
	}

private:
	std::vector<double> startProbability;
	std::vector<double> endProbability;
	std::vector<std::vector<double>> transitionProb;
}

















