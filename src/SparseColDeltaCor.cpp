#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat SparseColDeltaCor(const arma::sp_mat& e, 
                               const arma::sp_mat& d, 
                               const arma::sp_mat& knn) {
    
    int n_cells = e.n_cols;
    arma::sp_mat out_cor(n_cells, n_cells);
    
    // Iterate ONLY over the edges in the KNN graph
    for(arma::sp_mat::const_iterator it = knn.begin(); it != knn.end(); ++it) {
        
        int i = it.row(); // Source cell
        int j = it.col(); // Target cell
        
        // 1. Direct initialization: cast the sparse subviews directly to dense vectors
        // (This is much cleaner and bypasses the conv_to template error)
        arma::vec e_source(e.col(i));
        arma::vec e_target(e.col(j));
        arma::vec d_source(d.col(i));
        
        // 2. Calculate the Displacement Vector (e_j - e_i)
        arma::vec displacement = e_target - e_source;
        
        // 3. Calculate Pearson Correlation
        double pearson_cor = arma::as_scalar(arma::cor(displacement, d_source));
        
        // 4. Assign to the sparse output matrix
        out_cor(i, j) = pearson_cor;
    }
    
    return out_cor;
}