#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat SparseColDeltaCor(const arma::sp_mat& e, 
                               const arma::sp_mat& d, 
                               const arma::sp_mat& knn) {
    
    int n_cells = e.n_cols;
    
    // 1. Transpose the KNN mask so we iterate efficiently: 
    // Row becomes Target (j), Col becomes Source (i)
    arma::sp_mat knn_t = knn.t();
    
    arma::sp_mat out_cor(n_cells, n_cells);
    
    // Iterate over the transposed graph
    for(arma::sp_mat::const_iterator it = knn_t.begin(); it != knn_t.end(); ++it) {
        
        int j = it.row(); // Target cell
        int i = it.col(); // Source cell
        
        arma::vec e_source(e.col(i));
        arma::vec e_target(e.col(j));
        arma::vec d_source(d.col(i));
        
        // Calculate Displacement Vector (e_j - e_i)
        arma::vec displacement = e_target - e_source;
        
        // Calculate Pearson Correlation
        double pearson_cor = arma::as_scalar(arma::cor(displacement, d_source));
        
        // Velocyto Parity: Row = Target, Col = Source
        out_cor(j, i) = pearson_cor;
    }
    
    return out_cor;
}