#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat SparseEmbArrows(const arma::mat& emb, const arma::sp_mat& tp, double arrowScale=1.0, int nthreads=1) {
    
    int n_cells = emb.n_rows; 
    
    // Output matrix: 2 rows (x,y dims) by N cols (to match Velocyto's original structure)
    arma::mat dm(2, n_cells, arma::fill::zeros);
    
    #pragma omp parallel for num_threads(nthreads)
    for(int i = 0; i < n_cells; ++i) {
        
        // In Armadillo, sparse matrices are compressed by column. 
        // We can instantly fetch just the neighbors for cell i.
        arma::sp_mat::const_col_iterator it = tp.begin_col(i);
        arma::sp_mat::const_col_iterator it_end = tp.end_col(i);
        
        // 1. Count the neighbors (K) to establish the "null" baseline probability
        int K = 0;
        for(arma::sp_mat::const_col_iterator temp_it = it; temp_it != it_end; ++temp_it) {
            K++;
        }
        
        if (K == 0) continue; // If no neighbors, the cell doesn't move.
        
        double p_null = 1.0 / K;
        double ds_x = 0.0;
        double ds_y = 0.0;
        
        double emb_i_x = emb(i, 0);
        double emb_i_y = emb(i, 1);
        
        // 2. Calculate spatial displacement ONLY for actual neighbors
        for(; it != it_end; ++it) {
            int j = it.row(); // Target cell
            double p_actual = (*it); // Actual transition probability
            
            // Calculate coordinate difference
            double dx = emb(j, 0) - emb_i_x;
            double dy = emb(j, 1) - emb_i_y;
            
            // Calculate distance
            double dist = std::sqrt(dx*dx + dy*dy);
            
            if (dist > 0.0) { 
                // Normalize the distance vector and apply user scaling
                double norm_x = (dx / dist) * arrowScale;
                double norm_y = (dy / dist) * arrowScale;
                
                // Subtract the null background (diffusion) 
                // from the actual probability to find the true biological direction
                ds_x += norm_x * (p_actual - p_null);
                ds_y += norm_y * (p_actual - p_null);
            }
        }
        
        // Store the final net vector for cell i
        dm(0, i) = ds_x;
        dm(1, i) = ds_y;
    }
    
    return dm;
}