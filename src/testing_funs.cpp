#include <RcppEigen.h>
#include <vector> //std::vector
#include <string>
// #include <math> //for NAN
#include <iostream>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

//' Deprecated, slower old version of calc_overlap. Only included for testing.
//' @name calc_overlap_deprecated
//' @param cc_adj the cell-cell graph adjacency matrix.
//' @param cg_adj The cell-gene graph adjacency matrix.
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> calc_overlap_deprecated(Eigen::SparseMatrix<int> cc_adj,
                                                    Eigen::SparseMatrix<int> cg_adj) {
    
    // initialize vector to store triplets
    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> trp;
    Eigen::SparseMatrix<int> overlap_mat_all = cc_adj * cg_adj;
    Eigen::SparseMatrix<int> cc_tadj = cc_adj.transpose();
    
    // calcualte the rowSums of matrix cc_adj which is also the number of neighbourhoods of each cell
    std::vector<double> cell_nn_nums;
    for (int i=0; i < cc_tadj.outerSize(); i++){
        int k = 0;
        
        for (Eigen::SparseMatrix<int>::InnerIterator it(cc_tadj, i); it; ++it){  // Iterate over rows
            k += 1;
        }
        
        cell_nn_nums.push_back(k);
    }
    
    // loop over genes (its faster in column major matrix)
    for (int i=0; i < cg_adj.outerSize(); i++){
        
        // only preserve the edges which are shown in cg_adj matrix
        for (Eigen::SparseMatrix<int>::InnerIterator it(cg_adj, i); it; ++it){  // Iterate over rows
            
            double value = overlap_mat_all.coeffRef(it.row(), i)/cell_nn_nums[it.row()];
            
            trp.push_back(Trip(it.row(),
                               i, 
                               value));
        }
    }
    
    // build Matrix from triplets
    Eigen::SparseMatrix<double> overlap(cg_adj.rows(), cg_adj.cols());
    overlap.setFromTriplets(trp.begin(), trp.end());
    
    return overlap;
}


