
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
ComputeSNNasym(Eigen::Map<Eigen::SparseMatrix<int>> SNN,
               double prune,
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
