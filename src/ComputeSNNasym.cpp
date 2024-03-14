
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <string>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// ' Calculates SNN from adjacency matrix with uneven number of neighbours per
// ' row.
// '
// ' @description
// ' Calculates SNN from adjacency matrix with uneven number of neighbours per
// ' row.
// '
// ' @param SNN A sparse matrix (adjacency matrix)
// ' @param prune numeric. Below which Jaccard similarity edges should be
// ' removed.
// ' @param  mode The type of neighboring vertices to use for calculating
// similarity '  scores(Jaccard Index). Three options: "out", "in" and "all": '
// * "out": Select neighbouring vertices by out-going edges; ' * "in": Selecting
// neighbouring vertices by in-coming edges; ' * "all": Selecting neigbouring
// vertices by both in-coming and out-going edges. ' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double>
ComputeSNNasym(Eigen::Map<Eigen::SparseMatrix<int>> SNN, double prune,
               String mode) {
  // row Sums (of edges)
  Eigen::VectorXi k_i(SNN.rows());
  // Contains the edges two nodes have in common
  // in absolute numbers
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> res_dense;

  // Only outgoing edges considered
  if (mode == "out") {
    k_i = SNN * Eigen::VectorXi::Ones(SNN.cols());
    res_dense = SNN * (SNN.transpose());

    // Only ingoing edges considered
  } else if (mode == "in") {

    k_i = SNN.transpose() * Eigen::VectorXi::Ones(SNN.rows());
    res_dense = SNN.transpose() * SNN;

    // In- and outgoing edges considered
  } else if (mode == "all") {

    // Sums up in- and outgoing edges
    Eigen::SparseMatrix<int> sym(SNN.rows(), SNN.cols());

    // it is proven that initializing this transposed matrix here is better than
    // initializing it earlier.
    Eigen::SparseMatrix<int> transposed = SNN.transpose();
    sym = SNN + transposed;

    for (int i = 0; i < sym.outerSize(); ++i) {

      for (Eigen::SparseMatrix<int>::InnerIterator it(sym, i); it; ++it) {

        it.valueRef() = 1;
      }
    }

    // the product of two matrices might not be
    // sparse anymore. Storing the results into a
    // sparse matrix might raise up bad_alloc error
    res_dense = sym * (sym.transpose());

    k_i = sym * Eigen::VectorXi::Ones(sym.cols());
  }

  // Number of non-zero elements in each column.
  std::vector<int> nzs(res_dense.cols());

  typedef Eigen::Triplet<double> Trip;
  std::vector<Trip> trp;
  double overlapping;
  double a;
  double b;

  for (int i = 0; i < res_dense.cols(); ++i) { // number of columns ?

    int cnt = 0;
    for (int j = 0; j < res_dense.rows(); ++j) { // Iterate over rows

      a = k_i(j);
      b = k_i(i);

      if ((a + (b - res_dense(j, i))) != 0) {
        overlapping = res_dense(j, i) / (a + (b - res_dense(j, i)));
      } else {
        overlapping = 0;
      }
      if (overlapping >= prune) {
        trp.push_back(Trip(j, i, overlapping));
        cnt++;
      }
    }

    nzs[i] = cnt;
  }

  Eigen::SparseMatrix<double> res(res_dense.rows(), res_dense.cols());
  res.reserve(nzs);
  res.setFromTriplets(trp.begin(), trp.end());

  return res;
}

// Eigen::SparseMatrix<double> ComputeSNNasym_idx(Eigen::MatrixXd<int> cc_idx,
//                                            Eigen::MatrixXd<int> cg_idx,
//                                            Eigen::MatrixXd<int> gg_idx,
//                                            Eigen::MatrixXd<int> gc_idx = 0,
//                                            double prune,
//                                            String mode) {
//
//   // convert four index matrices to a cell-gene sparse SNN matrix
//   typedef Eigen::Triplet<double> T;
//   std::vector<T> tripletList;
//   int ngene = gg_idx.rows();
//   int ncell = cc_idx.rows();
//   int k = ngene + ncell;
//   int gc_idx_dim = 0;
//
//   if(gc_idx == 0){
//     gc_idx_dim = cg_idx.rows() * cg_idx.cols();
//   }else{
//     gc_idx_dim = gc_idx.rows() * gc_idx.cols();
//   }
//
//   tripletList.reserve( cc_idx.rows()*cc_idx.cols() +
//   cg_idx.rows()*cg_idx.cols() + gg_idx.rows()*gg_idx.cols() + gc_idx_dim );
//
//   // write cc-knn graph into triplets, the SNN graph are composed of four
//   sub-subgraphs
//   // with the subgraphs in the order of cc-knn, cg-knn, gc-knn, gg-knn from
//   left to right, from up to down for (int j=0; j< cc_idx.cols(); ++j){
//     for(int i=0; i< cc_idx.rows(); ++i) {
//       tripletList.push_back(T(i, cc_idx(i, j) - 1, 1));
//     }
//   }
//   // write cg-knn graph into triplets
//   for (int j = 0; j < cg_idx.cols(); ++j){
//     for (int i = 0; i < cg_idx.rows(); ++i){
//       tripletList.push_back(T(i, cg_idx(i, j) + ncell -1, 1));
//     }
//   }
//   // write gc-knn graph into tripletList
//   if (gc_idx == 0){
//     for (int j = 0; j < cg_idx.cols(); ++j){
//       for (int i = 0; i < cg_idx.rows(); ++i){
//         tripletList.push_back(T(cg_idx(i,j) + ncell - 1, i, 1));
//       }
//     }
//   } else{
//     for (int j = 0; j < gc_idx.cols(); ++j){
//       for (int i = 0; i < gc_idx.rows(); ++i){
//         tripletList.push_back(T(i + ncell, gc_idx(i, j) -1, 1));
//       }
//     }
//   }
//   // write gg-knn graph into the sparse cell-gene tripletList
//   for (int j = 0; j < gg_idx.cols(); ++j){
//     for (int i = 0; i < gg_idx.rows(); ++i){
//       tripletList.push_back(T(i + ncell, gg_idx(i, j) + ncell -1, 1));
//     }
//   }
//   // convert tripletList to sparse matrix
//   Eigen::SparseMatrix<int> SNN( k*k, k*k );
//   cg_adj.setFromTriplets(tripletList.begin(), tripletList.end());
//
//   Eigen::VectorXd k_i(SNN.rows());
//   // Eigen::VectorXd k_j(SNN.cols());
//
//   if (mode == "out"){
//     k_i = SNN * Eigen::VectorXd::Ones(SNN.cols());
//     SNN = SNN * (SNN.transpose());
//   }else if (mode == "in"){
//     k_i = SNN.transpose() * Eigen::VectorXd::Ones(SNN.rows());
//     SNN =  (SNN.transpose()) * SNN;
//   }else if (mode == "all"){
//     Eigen::SparseMatrix<double> sym;
//     sym = SNN + (Eigen::SparseMatrix<double> (SNN.transpose()));
//
//     for (int i=0; i < sym.outerSize(); ++i){
//       for (Eigen::SparseMatrix<double>::InnerIterator it(sym, i); it; ++it){
//         int ki = it.row();
//         int kj = it.col();
//
//         // if (SNN(ki, kj) + SNN(kj,ji) >0){
//         it.valueRef() = 1;
//         // }
//       }
//     }
//
//     SNN = sym * (sym.transpose());
//     k_i = sym * Eigen::VectorXd::Ones(sym.cols());
//     // std::cout << "all" << std::endl
//     }
//
//   for (int i=0; i < SNN.outerSize(); ++i){  //number of columns ?
//     for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){ //
//     Iterate over rows
//
//       int ki = it.row();
//       int kj = it.col();
//
//       double a;
//       a = k_i(ki);
//
//       double b;
//       b = k_i(kj);
//
//       it.valueRef() = it.value()/(a + (b - it.value()));
//       if(it.value() < prune){
//         it.valueRef() = 0;
//       }
//     }
//   }
//   SNN.prune(0.0); // actually remove pruned values
//   return SNN;
// }

// Eigen::SparseMatrix<double> ComputeSNNasym_copy(Eigen::SparseMatrix<double>
// SNN,
//                                            double prune,
//                                            String mode) {
//
//   Eigen::VectorXd k_i(SNN.rows());
//   // Eigen::VectorXd k_j(SNN.cols());
//
//   if (mode == "out"){
//     k_i = SNN * Eigen::VectorXd::Ones(SNN.cols());
//     SNN = SNN * (SNN.transpose());
//   }else if (mode == "in"){
//     k_i = SNN.transpose() * Eigen::VectorXd::Ones(SNN.rows());
//     SNN =  (SNN.transpose()) * SNN;
//   }else if (mode == "all"){
//     Eigen::SparseMatrix<double> sym;
//     sym = SNN + (Eigen::SparseMatrix<double> (SNN.transpose()));
//
//     for (int i=0; i < sym.outerSize(); ++i){
//       for (Eigen::SparseMatrix<double>::InnerIterator it(sym, i); it; ++it){
//         int ki = it.row();
//         int kj = it.col();
//
//         // if (SNN(ki, kj) + SNN(kj,ji) >0){
//         it.valueRef() = 1;
//         // }
//       }
//     }
//
//     SNN = sym * (sym.transpose());
//     k_i = sym * Eigen::VectorXd::Ones(sym.cols());
//     // std::cout << "all" << std::endl
//     }
//   for (int i=0; i < SNN.outerSize(); ++i){  //number of columns ?
//     for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){ //
//     Iterate over rows
//
//       int ki = it.row();
//       int kj = it.col();
//
//       double a;
//       a = k_i(ki);
//
//       double b;
//       b = k_i(kj);
//
//       it.valueRef() = it.value()/(a + (b - it.value()));
//       if(it.value() < prune){
//         it.valueRef() = 0;
//       }
//     }
//   }
//   SNN.prune(0.0); // actually remove pruned values
//   return SNN;
// }
