#include <RcppEigen.h>
#include <vector> //std::vector
#include <string>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]


//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> calc_overlap(Eigen::SparseMatrix<int> cc_adj, 
                                         Eigen::SparseMatrix<int> cg_adj) {

  // store the results in sparse matrix
  Eigen::SparseMatrix<double> overlap(cg_adj.rows(), cg_adj.cols());
  
  // loop over genes (its faster in column major matrix)
  for (int i=0; i < cg_adj.outerSize(); i++){
    
    std::vector<int> cell_idxs;
    
    // get the idxs of all the cells that have edge to gene
    for (Eigen::SparseMatrix<int>::InnerIterator it(cg_adj, i); it; ++it){  // Iterate over rows
      cell_idxs.push_back(it.row());
    
    }
    
    // for each cell that has an edge, get the nearest neighbours
    for(int k = 0; k < cell_idxs.size(); k++){

      std::vector<int> cell_NNs;
      
      // loop over cell-cell Adj. matrix to get the neighbors.
      // is there a way so we dont have to loop over the whole matrix?
      for(int l=0; l < cc_adj.outerSize(); l++){
        for (Eigen::SparseMatrix<int>::InnerIterator it(cc_adj, l); it; ++it){
          
          if(it.row() == cell_idxs[k]){
            cell_NNs.push_back(it.col());
          } else {
            continue;
          }
          
        }
      }
    
    
      std::vector<int> n_overlap; //size of overlap
        
      std::sort(cell_NNs.begin(), cell_NNs.end());
      std::sort(cell_idxs.begin(), cell_idxs.end());
      
      // get the intersection of cells that have an edge to the gene
      // and that are neighbours of the cell with index cell_idxs[k]
      std::set_intersection(cell_NNs.begin(), cell_NNs.end(),
                            cell_idxs.begin(), cell_idxs.end(),
                            std::back_inserter(n_overlap));
      
      
      double perc_ov;
      double ov_size = n_overlap.size();
      double nn_size = cell_NNs.size();
      perc_ov = ov_size / nn_size; //percent of overlap
      
      // update entry in overlap matrix
      overlap.coeffRef(cell_idxs[k], i) = perc_ov;
  
    }
  }
  
  // std::cout << Eigen::MatrixXd(overlap) << std::endl;
  
  overlap.makeCompressed(); //removes an error with R Matrix
  return overlap;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// 
// library(Matrix)
// library(tictoc)
// cc_adj <- data.frame("c1" = c(1, 1, 0, 1, 0, 0),
//                      "c2" = c(1, 1, 0, 0, 1, 0),
//                      "c3" = c(0, 0, 1, 1, 0, 1),
//                      "c4" = c(0, 0, 1, 1, 0, 1),
//                      "c5" = c(0, 1, 0, 0, 1, 0),
//                      "c6" = c(0, 0, 1, 1, 0, 1))
// cc_adj <- Matrix(as.matrix(cc_adj))
// rownames(cc_adj) <- paste0("c", 1:6)
// 
// cg_adj <- data.frame(g1 = c(1, 0, 0, 0, 0, 0),
//                      g2 = c(1, 1, 0, 0, 0, 0),
//                      g3 = c(0, 1, 1, 0, 0, 0),
//                      g4 = c(0, 0, 1, 1, 0, 0),
//                      g5 = c(0, 0, 0, 1, 1, 1),
//                      g6 = c(0, 0, 0, 0, 1, 1),
//                      g7 = c(0, 0, 0, 1, 1, 1))
// 
// cg_adj <- Matrix(as.matrix(cg_adj))
// rownames(cg_adj) <- paste0("c", 1:6)
// 
// tic()
// overlap_mat <- determine_overlap(cg_adj,
//                                  cc_adj)
// toc()
// tic()
// test <- calc_overlap(cc_adj, cg_adj)
// toc()
// 
// cg<-as.matrix(cg_adj)
// cc <- as.matrix(cc_adj)
// microbenchmark::microbenchmark(
//   a<-determine_overlap(cg,
//                     cc),
//   b<-calc_overlap(cc_adj, cg_adj)
// )
// */
