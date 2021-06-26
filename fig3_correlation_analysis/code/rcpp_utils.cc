#include <algorithm>
#include <vector>

#include <Rcpp.h>

//' Worker function to complement knn_pseudobulk function.
//' Selects a set of non-overlapping pseudobulks given the distance matrix
//' @param mat Logical symmetric matrix with TRUE for any pairs that overlap
//' @return Logical vector with TRUE for a set of non-overlapping pseudobulks
// [[Rcpp::export]]
Rcpp::LogicalVector filter_pseudobulk_helper_cpp(const Rcpp::LogicalMatrix &mat) {
    if (mat.nrow() != mat.ncol())
        Rcpp::stop("Mat must have same number of rows and columns");
    for (int i = 0; i < mat.nrow(); i++) {
        if(mat(i,i))
            Rcpp::stop("Mat must have FALSE along diagonal");
    }

    Rcpp::LogicalVector keep(mat.nrow(), true);
    Rcpp::IntegerVector overlap_counts = Rcpp::colSums(mat);
    
    bool done;
    size_t discard;
    int max_overlaps;
    while(true) {
        // Find the pseudobulk with the most overlaps
        max_overlaps = 0;
        for (int i = 0; i < mat.nrow(); i++) {
            if (keep[i] && overlap_counts[i] > max_overlaps) {
                discard = i;
                max_overlaps = overlap_counts[i];
            }
        }
        if (max_overlaps == 0) break;
        
        // Discard the pseudobulk and update overlap_counts
        keep[discard] = false;
        for (int i = 0; i < mat.nrow(); i++) {
            overlap_counts[i] -= mat(i, discard);
        }
    }
    return keep;
}



typedef Rcpp::NumericVector::iterator DoubleIterator;
typedef Rcpp::NumericVector::const_iterator const_DoubleIterator;

//' Computes ranks on a sparse vector x, equivalent to base R rank(x, ties.method="average")
//' Shifts ranks in x such that the implicit zero entries will have rank 0.
//' @param in_start input start iterator
//' @param in_end input end iterator
//' @param out_start output start iterator
//' @param out_end output end iterator
//' @param ord Working space integer vector for sorting ranks
//' @param len length of sparse vector x (number of zeros is assumed to be len-x.size())
void fractional_rank(const_DoubleIterator in_start, const_DoubleIterator in_end,
                     DoubleIterator out_start, DoubleIterator out_end,
                     std::vector<int> &ord, size_t len) {
    // 1. Get sorted order of indices
    // 2. Iterate through sorted indices, and copy over ranks
    // Notes: 
    //  - zero rank will be nneg + (1+nzero)/2
    //  - negative values will shift rank down by zero_rank
    //  - positive values will shift rank down by nzero_explicit, then up by
    //    zero_rank
    int ties, nneg, nzero_implicit, nzero_explicit;
    size_t sz = in_end - in_start;
    
    if (ord.size() < sz) ord.resize(sz);
    for (size_t i = 0; i < sz; i++) ord[i] = i;
    
    nzero_implicit = len - sz;
    nzero_explicit = 0;
    nneg = 0;
    for (size_t i = 0; i < sz; i++) {
        nzero_explicit += in_start[i] == 0;
        nneg += in_start[i] < 0;
    }
    
    std::sort(ord.begin(), ord.begin() + sz, [&](int a, int b){ 
        return in_start[a] < in_start[b];
    });
    
    double zero_rank = nneg + (1+nzero_explicit+nzero_implicit)/2.0;
    size_t i = 0;
    // Loop over negatives
    for (; i < sz && in_start[ord[i]] < 0; i += ties) {
        ties = 1;
        while (i + ties < sz && in_start[ord[i]] == in_start[ord[i+ties]]) ties++;
        for (size_t j = 0; j < ties; j++) {
            out_start[ord[i+j]] = i + 1 + (ties - 1)/2.0 - zero_rank;
        }
    }
    // Loop over explicit zeros
    for (; i < sz && in_start[ord[i]] == 0; i += ties) {
        ties = 1;
        while (i + ties < sz && in_start[ord[i]] == in_start[ord[i+ties]]) ties++;
        for (size_t j = 0; j < ties; j++) {
            out_start[ord[i+j]] = 0.0;
        }
    }
    // Loop over positives
    for (; i < sz; i += ties) {
        ties = 1;
        while (i + ties < sz && in_start[ord[i]] == in_start[ord[i+ties]]) ties++;
        for (size_t j = 0; j < ties; j++) {
            out_start[ord[i+j]] = i + 1 + (ties - 1)/2.0 + nzero_implicit - zero_rank;
        }
    }
}

//' Worker function to calculate columwise ranks on a sparse matrix. Ties
//' are computed as with rank(x, ties.method="average")
//' @param mat dgCMatrix
//' @return dgCMatrix with column-wise ranks for the entries of mat, shifted
//'    down such that rank(0) = 0
// [[Rcpp::export]]
Rcpp::S4 column_rank(const Rcpp::S4 &mat) {
    if (!mat.inherits("dgCMatrix")) Rcpp::stop("mat must be a dgCMatrix");
    Rcpp::S4 ret = clone(mat);
    
    Rcpp::NumericVector x = mat.slot("x");
    Rcpp::IntegerVector p = mat.slot("p");
    Rcpp::NumericVector ret_x = ret.slot("x");

    size_t cols = p.size() - 1;
    size_t rows = Rcpp::as<Rcpp::NumericVector>(mat.slot("Dim"))[0];
    std::vector<int> ord(cols);

    for (int i = 0; i < cols; i++) {
        fractional_rank(
            x.begin()+p[i], x.begin()+p[i+1], 
            ret_x.begin()+p[i], ret_x.begin()+p[i+1], 
            ord, rows
        );        
    }

    return ret;
}

/*** 
# Comprehensive testing example, finishing with showing how
# these ranked matrices can be used for calculation of spearman correlation

# Generate a dense matrix with 50% zero, and 50% numbers in range -0.5, 0.5
set.seed(1252)
nrow <- 6
ncol <- 100
dense <- matrix(sample(c(rep(0,10),runif(10)-0.5), nrow*ncol, replace=TRUE), ncol=ncol)
sparse <- as(dense, "dgCMatrix")

sparse_explicit_zeros <- sparse
sparse_explicit_zeros@x[rbernoulli(length(sparse_explicit_zeros@x),p=0.2)] <- 0
dense_explicit_zeros <- as.matrix(sparse_explicit_zeros)

dense_rank <- apply(dense, 2, function(x) rank(x, ties.method="average")-sum(x<0)-(1+sum(x==0))/2)
dense_explicit_zeros_rank <- apply(dense_explicit_zeros, 2, function(x) rank(x, ties.method="average")-sum(x<0)-(1+sum(x==0))/2)

sparse_rank <- column_rank(sparse)
sparse_explicit_zeros_rank <- column_rank(sparse_explicit_zeros)

stopifnot(all.equal(
  as(dense_rank, "dgCMatrix"),
  sparse_rank
))

stopifnot(all.equal(
  as(dense_explicit_zeros_rank, "dgCMatrix"),
  sparse_explicit_zeros_rank
))

stopifnot(all.equal(
  cor(dense, method="spearman"),
  cor(as.matrix(sparse_rank), method="pearson")
))
*/





// Overall goal: Provide a subroutine useful for finding correlations for a subset
// of pairwise combinations of matrix columns. In this case, the easiest subroutine
// I can think of calculates a selected subset of entries from the matrix multiplication.
// Then it should be reasonably efficient to make the adjustments based on per-column
// means and standard deviations.

// To handle possible combinations of sparse & dense matrices, I've implemented
// sparse-sparse, sparse-dense, and dense-dense ops.



//' Compute partial entries of t(mat1) %*% mat2, for 2 sparse matrices
//' @param mat1 dgCMatrix 1
//' @param mat2 dgCMatrix 2
//' @param idx1 Desired rows of output entries (columns of mat1). 0-based indexes
//' @param idx2 Desired columns of output entries (columns of mat2). 0-based indexes
//' @return Values of mat1[,idx1] %*% mat2[,idx2] (same length as idx1 and idx2)
// [[Rcpp::export]]
Rcpp::NumericVector partialDot_sparse_sparse_cpp(
        const Rcpp::S4 &mat1, const Rcpp::S4 &mat2,
        const Rcpp::IntegerVector &idx1, const Rcpp::IntegerVector &idx2) {
    if (idx1.size() != idx2.size()) Rcpp::stop("idx1 and idx2 must have same length");
    size_t sz = idx1.size();

    // EXTRACT DATA INTO C++
    Rcpp::NumericVector x1 = mat1.slot("x");
    Rcpp::IntegerVector i1 = mat1.slot("i");
    Rcpp::IntegerVector p1 = mat1.slot("p");
    size_t rows1 = Rcpp::as<Rcpp::NumericVector>(mat1.slot("Dim"))[0];

    Rcpp::NumericVector x2 = mat2.slot("x");
    Rcpp::IntegerVector i2 = mat2.slot("i");
    Rcpp::IntegerVector p2 = mat2.slot("p");
    size_t rows2 = Rcpp::as<Rcpp::NumericVector>(mat2.slot("Dim"))[0];

    if (rows1 != rows2) Rcpp::stop("mat1 and mat2 must have same number of rows");

    // CONFIRM COLUMNS ARE ORDERED
    for (size_t i = 0; i < p1.size()-1; i++) {
        for (size_t j = p1[i]; j < p1[i+1]-1; j++) {
            if (i1[j] >= i1[j+1]) Rcpp::stop("mat1 must have each column sorted by row index");
        }
    }

    for (size_t i = 0; i < p2.size()-1; i++) {
        for (size_t j = p2[i]; j < p2[i+1]-1; j++) {
            if (i2[j] >= i2[j+1]) Rcpp::stop("mat2 must have each column sorted by row index");
        }
    }

    // COMPUTE ACTUAL DOT PRODUCTS
    Rcpp::NumericVector ret(sz);
    double accumulator;
    size_t entry1, entry2;
    for (size_t i = 0; i < sz; i++) {
        // Compute single dot product. Walk through the corresponding columns of
        // mat1 and mat2 in order by row, adding up the value whenever their
        // columns overlap
        accumulator = 0.0;
        entry1 = p1[idx1[i]];
        entry2 = p2[idx2[i]];
        while (entry1 < p1[idx1[i] + 1] && entry2 < p2[idx2[i] + 1]) {
            if (i1[entry1] == i2[entry2]) {
                accumulator += x1[entry1] * x2[entry2];
                entry1++;
                entry2++;
            } else if(i1[entry1] < i2[entry2]) {
                entry1++;
            } else {
                entry2++;
            }
        }
        ret[i] = accumulator;
    }
    return ret;
}

//' Compute partial entries of t(mat1) %*% mat2, for mat1 sparse, mat2 dense
//' @param mat1 dgCMatrix 1
//' @param mat2 dense matrix 2
//' @param idx1 Desired rows of output entries (columns of mat1). 0-based indexes
//' @param idx2 Desired columns of output entries (columns of mat2). 0-based indexes
//' @return Values of mat1[,idx1] %*% mat2[,idx2] (same length as idx1 and idx2)
// [[Rcpp::export]]
Rcpp::NumericVector partialDot_sparse_dense_cpp(
        const Rcpp::S4 &mat1, const Rcpp::NumericMatrix &mat2,
        const Rcpp::IntegerVector &idx1, const Rcpp::IntegerVector &idx2) {
    if (idx1.size() != idx2.size()) Rcpp::stop("idx1 and idx2 must have same length");
    size_t sz = idx1.size();

    // EXTRACT DATA INTO C++
    Rcpp::NumericVector x1 = mat1.slot("x");
    Rcpp::IntegerVector i1 = mat1.slot("i");
    Rcpp::IntegerVector p1 = mat1.slot("p");
    size_t rows1 = Rcpp::as<Rcpp::NumericVector>(mat1.slot("Dim"))[0];

    if (rows1 != mat2.nrow()) Rcpp::stop("mat1 and mat2 must have same number of rows");

    // COMPUTE ACTUAL DOT PRODUCTS
    Rcpp::NumericVector ret(sz);
    double accumulator;
    for (size_t i = 0; i < sz; i++) {
        accumulator = 0.0;
        for (size_t j = p1[idx1[i]]; j < p1[idx1[i]+1]; j++) {
            accumulator += x1[j] * mat2(i1[j], idx2[i]);
        }
        ret[i] = accumulator;
    }
    return ret;
}

//' Compute partial entries of t(mat1) %*% mat2, for mat1 dense, mat2 dense
//' @param mat1 dense matrix 1
//' @param mat2 dense matrix 2
//' @param idx1 Desired rows of output entries (columns of mat1). 0-based indexes
//' @param idx2 Desired columns of output entries (columns of mat2). 0-based indexes
//' @return Values of mat1[,idx1] %*% mat2[,idx2] (same length as idx1 and idx2)
// [[Rcpp::export]]
Rcpp::NumericVector partialDot_dense_dense_cpp(
        const Rcpp::NumericMatrix &mat1, const Rcpp::NumericMatrix &mat2,
        const Rcpp::IntegerVector &idx1, const Rcpp::IntegerVector &idx2) {
    if (idx1.size() != idx2.size()) Rcpp::stop("idx1 and idx2 must have same length");
    size_t sz = idx1.size();

    if (mat1.nrow() != mat2.nrow()) Rcpp::stop("mat1 and mat2 must have same number of rows");
    size_t nrow = mat1.nrow();

    // COMPUTE ACTUAL DOT PRODUCTS
    Rcpp::NumericVector ret(sz);
    double accumulator;
    for (size_t i = 0; i < sz; i++) {
        accumulator = 0.0;
        for (size_t j = 0; j < nrow; j++) {
            accumulator += mat1(j, idx1[i]) * mat2(j, idx2[i]);
        }
        ret[i] = accumulator;
    }
    return ret;
}

/*** 
# Comprehensive testing example for normal functioning of C++ code
library(Matrix)
# Generate a dense matrix with 50% zero, and 50% numbers in range -0.5, 0.5
set.seed(1252)
nrow <- 6
ncol <- 100
comparison_count <- 100
dense1 <- matrix(sample(c(rep(0,10),runif(10)-0.5), nrow*ncol, replace=TRUE), ncol=ncol)
dense2 <- matrix(sample(c(rep(0,10),runif(10)-0.5), nrow*ncol, replace=TRUE), ncol=ncol)
sparse1 <- as(dense1, "dgCMatrix")
sparse2 <- as(dense2, "dgCMatrix")

comparisons <- matrix(sample.int(ncol, comparison_count*2, replace=TRUE), ncol=2) - 1

res1 <- partialDot_dense_dense(dense1, dense2, comparisons[,1], comparisons[,2])
res2 <- partialDot_sparse_dense(sparse1, dense2, comparisons[,1], comparisons[,2])
res3 <- partialDot_sparse_sparse(sparse1, sparse2, comparisons[,1], comparisons[,2])

res_official <- (t(dense1) %*% dense2)[comparisons + 1]
stopifnot(all.equal(res1, res_official))
stopifnot(all.equal(res2, res_official))
stopifnot(all.equal(res3, res_official))
*/