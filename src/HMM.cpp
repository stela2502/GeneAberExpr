#include "GeneAberExpr.h"

#include <Rcpp.h>
#include <RcppEigen.h>
#include <progress.hpp>

// here I need to allow R to run a HMM analysis.
// I need the forwards/backward and all other HMM algorithms in place here and work on sparse matrices.
// should be quite straight forward based on my perl hmm class.
// all probability estimations are in the <ProbList> class.
// get a functional one with ProbList( <sparseMatrix>, idsA, idsB, starts )
// the starts should range from min to max of the matrix.
