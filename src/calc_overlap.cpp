
#include <RcppEigen.h>
#include <vector> //std::vector
#include <string>
// #include <math> //for NAN
#include <iostream>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

//' WIP replacement for `determine_overlap` function
//' @name calc_overlap
//' @description
//' c++ implementation for calculating the cell-neighour-overlapping among k nearest cell neighbours of each gene.
//' @param cc_adj sparse matrix (dgCMatrix), cell-cell adjacency matrix
//' @param cg_adj sparse matrix (dgCMatrix), cell-gene adjacency matrix
//' @export
// [[Rcpp::export]]

Eigen::SparseMatrix<double> calc_overlap_idx(Eigen::MatrixXd<int> cc_idx,
                                         Eigen::MatrixXd<int> cg_idx,
                                         double overlap) {


  //convert the kNN index matrices to sparse kNN adjacency matrices
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  int k = cg_idx.cols();
  tripletList.reserve(cg_idx.rows() * cg_idx.cols());
  for(int j=0; j< k; ++j){
    for(int i=0; i<cg_idx.rows(); ++i) {
      tripletList.push_back(T(i, cg_idx(i, j) - 1, 1));
    }
  }
  Eigen::SparseMatrix<int> cg_adj(cg_idx.rows(), cg_idx.rows());
  cg_adj.setFromTriplets(tripletList.begin(), tripletList.end());

  std::vector<T> tripletList;
  int k = cc_idx.cols();
  tripletList.reserve(cc_idx.rows() * cc_idx.cols());
  for(int j=0; j< k; ++j){
    for(int i=0; i<cc_idx.rows(); ++i) {
      tripletList.push_back(T(i, cc_idx(i, j) - 1, 1));
    }
  }
  Eigen::SparseMatrix<int> cc_adj(cc_idx.rows(), cc_idx.rows());
  cg_adj.setFromTriplets(tripletList.begin(), tripletList.end());

  // calculating the overlapping by producting two matrices
  Eigen::SparseMatrix<int> overlap_mat_all = cc_adj * cg_adj;
  // Eigen::SparseMatrix<int> cc_tadj = cc_adj.transpose();
  free(cc_adj)


  // calcualte the rowSums of matrix cc_adj which is also the number of neighbourhoods of each cell
  // the number of neighbourhoods of each cell is the number of columns in the matrix cc_idx, it is the same for all the cells
  std::vector<double> cell_nn_nums;

  for (int i = 0; i < k; k++){
    cell_nn_nums.push_back(k)
  }


  // The codes to output a dense MatrixXd which stores the indeces of cell-gene connections
  Eigen::MatrixXd<double> overlap_idx(cg_adj.rows(), k) = NA_INTEGER;
  for (int i=0; i < cg_adj.outerSize(); i++){
    int tmp = 0;
    // only preserve the edges which are shown in cg_adj matrix
    for (Eigen::SparseMatrix<int>::InnerIterator it(cg_adj, i); it; ++it){  // Iterate over rows

      double value = overlap_mat_all.coeffRef(it.row(), i)/cell_nn_nums[it.row()];

      if (value >= overlap){
        // overlap[it.row(), i] = value;
        overlap_idx[i, tmp] = it.row();
        tmp++;
      }

    }
    }

    return overlap_idx;



}


Eigen::SparseMatrix<double> calc_overlap(Eigen::SparseMatrix<int> cc_adj,
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


// Eigen::SparseMatrix<double> calc_overlap(Eigen::SparseMatrix<int> cc_adj,
//                                          Eigen::SparseMatrix<int> cg_adj) {
//
//   // initialize vector to store triplets
//   typedef Eigen::Triplet<double> Trip;
//   std::vector<Trip> trp;
//   Eigen::SparseMatrix<int> overlap_mat_all = cc_adj * cg_adj;
//   Eigen::SparseMatrix<int> cc_tadj = cc_adj.transpose();
//
//   // calcualte the rowSums of matrix cc_adj which is also the number of neighbourhoods of each cell
//   std::vector<double> cell_nn_nums;
//   for (int i=0; i < cc_tadj.outerSize(); i++){
//     int k = 0;
//
//     for (Eigen::SparseMatrix<int>::InnerIterator it(cc_tadj, i); it; ++it){  // Iterate over rows
//       k += 1;
//     }
//
//     cell_nn_nums.push_back(k);
//   }
//
//   // loop over genes (its faster in column major matrix)
//   for (int i=0; i < cg_adj.outerSize(); i++){
//
//     // only preserve the edges which are shown in cg_adj matrix
//     for (Eigen::SparseMatrix<int>::InnerIterator it(cg_adj, i); it; ++it){  // Iterate over rows
//
//       double value = overlap_mat_all.coeffRef(it.row(), i)/cell_nn_nums[it.row()];
//
//       trp.push_back(Trip(it.row(),
//                          i,
//                          value));
//     }
//     }
//
//   // build Matrix from triplets
//   Eigen::SparseMatrix<double> overlap(cg_adj.rows(), cg_adj.cols());
//   overlap.setFromTriplets(trp.begin(), trp.end());
//
//   return overlap;
// }

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
