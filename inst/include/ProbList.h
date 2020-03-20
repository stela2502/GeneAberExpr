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

	pos4val ( double val ){
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

	double Prob_4_value ( double val, int state) {
		if ( this->lastVal == val ){
			return ( this->proA[ this->lastID ] );
		}
		this->pro[state][ pos4val( val ) ];
	};


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
	double startProbability;
	double endProbability;
	double probability_for_change_to_state_A;
	double probability_for_change_to_state_B;
	//a vector of chars with max 10 entries...
	std::vector<std::string> states;

	ProbList () {};

	void estimate( Eigen::SparseMatrix<double> data, std::vector<double> range, std::string name ){
		this->states.push_back( name );

		this->list.reserve( data.innerSize() );
		Rcout << "estimating probabilities for state " << state << std::endl;
		Progress p(data.innerSize(), true);
		data = data.transpose();

		for ( int c_=0; c_ < data.outerSize(); ++c_ ){
			std::fill(D.begin(), D.end(), 0.0);
			for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
				D[it.row()] =  it.value();
			}
			list[i]=ProbEntry( range ) ;
			this->list[c_].estimate( D, this->states.size()-1 );
			p.increment();
		}
	}

	double Prob_4_value ( int i, double val, int state) {
		return (list[i].Prob_4_value( val, state ));
	};
	double Prob_4_value ( double val, int state) {
		return (list[0].Prob_4_value( val, state ));
	};
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
		this.transitionProb[from][to];
	};
	double getStartProb( int state ) {
		this.startProbability[state];
	};
	double getEndProb( int state ) {
		this.endProbability[state];
	};

private:
	std::vector<double> startProbability;
	std::vector<double> endProbability;
	std::vector<std::vector<double>> transitionProb;
}

















