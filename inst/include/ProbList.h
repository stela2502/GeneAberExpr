#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>

class ProbEntry{
	public:

	std::vector<double> probA (10);
	std::vector<double> probB (10);
	std::vector<double> start (10);
	std::vector<int>    order (10); // in which order should the start be evaluated?

	typedef std::pair<int, double> paired;

	ProbEntry () { };

	ProbEntry (std::vector<double> starts ) {
		if ( starts.size() != 10 ) {
			::Rf_error("the start sizes need to be exactly 10!" );
		}
		this->start = starts;
		//this->proA.reserve(starts.size() );
		//this->proB.reserve(starts.size() );
		//this->tmp.reserve(starts.size() );
		for ( int i = 0; i <starts.size() +1; i++ ){
			if (this->probA.size() < i ){
				this->probA.push_back( 0.0 );
				this->probB.push_back( 0.0 );
				this->tmp.push_back( 0.0 );
				this->order.push_back(i);
			}else {
				probA[i] = 0.0;
				probB[i] = 0.0;
				order[i] = i;
			}
		}
	};

	pos4val ( double val ){
		for (int i =0; i< this->order.size(); i++){
			if ( this->start[i] <= val & this->order.size() == i+1 ){
				return i;
			}
			if ( this->start[i] <= val & this->start[i+1] > val ){
				return i;
			}
		}
	};

	void estimateA ( std::vector<double> values ) {
		this->probA = estimate( values );
		fixOrder();
	};

	void estimateB ( std::vector<double> values ) {
		this->probB = estimate( values );
		fixOrder();
	};

	double Prob_B_4_value ( double val) {
		this->proB[ pos4val( val ) ];
	};
	double Prob_A_4_value( double val) {
		this->proA[ pos4val( val ) ];
	};


private:
	std::vector<double> tmp;


	bool cmp_second(const paired & left, const paired & right) {
    	return left.second < right.second;
	};

	std::vector<double> estimate ( std::vector<double> values ) {
		std::vector<double> tmp (this->proB.size());
		std::vector<double> prob (this->proB.size());

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
			tmp[i] = (prob[i] - min) / sum + 1e-9 ; //we can not have 0 values later on!
		}
		return (tmp);
	}

	void fixOrder() {
		std::vector<paired> pairs;
		pairs.reserve(this->probA.size());
		for ( i = 0; i < this->probA.size(); i++ ){
			pairs.push_back(std::make_pair(i, this->probA[i] + this->probB[i] ) );
		}		
		std::sort(pairs.begin(), pairs.end(), cmp_second<paired>);

		for ( i = 0; i < tmp.size(); i++ ){
			this->order[i] = pairs[i].first;
		}
	}

};

class ProbList {
public:

	std::vector<ProbEntry> list;

	ProbList () {};

	ProbList (Eigen::SparseMatrix<double> data,  
		std::vector<int>A; std::vector<int>B, std::vector<double> range ) 
	{
		if ( range.size() != 10 ) {
			::Rf_error("the ranges length needs to be exactly 10!" );
		}
		//here I need to initialize the list to data.innerSize();
		this->list.reserve( data.innerSize() );

		Rcout << "estimating probabilities" << std::endl;
		Progress p(data.innerSize(), true);

		data = data.transpose();
		std::vector<double> D(X.innerSize());
		std::vector<double>valA (A.size());
		std::vector<double>valB (B.size());
		for ( int c_=0; c_ < X.outerSize(); ++c_ ){
			std::fill(D.begin(), D.end(), 0.0);
			for (Eigen::SparseMatrix<double>::InnerIterator it(X, c_); it; ++it){
				D[it.row()] =  it.value();
			}
			for( int i = 0; i < A.size(); i ++){
				valA[i] = D[A[i]];
			}
			for( int i = 0; i < A.size(); i ++){
				valB[i] = D[B[i]];
			}
			list[i]=ProbEntry( range ) ;
			list[i].estimateA( valA );
			list[i].estimateB( valB );
			p.increment();
		}
	};
	double Prob_B_4_value ( int i, double val) {
		return (list[i].Prob_B_4_value( val ));
	};
	double Prob_A_4_value ( int i, double val) {
		return (list[i].Prob_A_4_value( val ));
	};
};


















