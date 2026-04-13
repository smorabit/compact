#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Helper function: Calculate cosine similarity natively on sparse vectors
double sparse_cosine(const arma::sp_mat& v1, const arma::sp_mat& v2) {
    // Calculate dot product
    double dot_prod = arma::as_scalar(v1.t() * v2);
    
    // Calculate vector magnitudes (Euclidean norm)
    double norm1 = arma::norm(v1);
    double norm2 = arma::norm(v2);
    
    // Prevent division by zero
    if (norm1 == 0.0 || norm2 == 0.0) return 0.0;
    
    return dot_prod / (norm1 * norm2);
}

// [[Rcpp::export]]
arma::sp_mat SparseColDeltaCor(const arma::sp_mat& e, 
                               const arma::sp_mat& d, 
                               const arma::sp_mat& knn) {
    
    int n_cells = e.n_cols;
    
    // Initialize an empty sparse matrix for the output
    arma::sp_mat out_cor(n_cells, n_cells);
    
    // Iterate ONLY over the non-zero edges in the KNN graph
    for(arma::sp_mat::const_iterator it = knn.begin(); it != knn.end(); ++it) {
        
        // it.row() is the target cell, it.col() is the source cell
        int r = it.row();
        int c = it.col();
        
        // Extract the sparse column vectors
        arma::sp_mat e_col = e.col(c);
        arma::sp_mat d_col = d.col(r);
        
        // Compute cosine similarity and assign to the output matrix
        out_cor(r, c) = sparse_cosine(e_col, d_col);
    }
    
    return out_cor;
}