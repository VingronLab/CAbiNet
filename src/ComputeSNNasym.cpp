
#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

//' Calculates SNN from adjacency matrix with uneven number of neighbours per
//' row.
//' 
//' @description
//' Calculates SNN from adjacency matrix with uneven number of neighbours per
//' row.
//'
//' @param SNN A sparse matrix (adjacency matrix)
//' @param prune numeric. Below which Jaccard similarity edges should be 
//' removed.
//' 
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> ComputeSNNasym(Eigen::SparseMatrix<double> SNN, double prune) {

  // int k_i;
  // int k_j;

  // k_i = SNN.colwise().sum();
  // k_j = SNN.rowwise().sum();
  Eigen::VectorXd k_i(SNN.rows());
  Eigen::VectorXd k_j(SNN.cols());

  k_i = SNN * Eigen::VectorXd::Ones(SNN.cols());
  // k_j = Eigen::VectorXd::Ones(SNN.rows()).transpose() * SNN;

  // std::cout << "Row sums:\n" << k_i << std::endl;
  // std::cout << "col_sums:\n" << k_j << std::endl;

  SNN = SNN * (SNN.transpose());
  for (int i=0; i < SNN.outerSize(); ++i){  //number of columns ?
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){  // Iterate over rows

      int ki = it.row();
      int kj = it.col();

      double a;
      a = k_i(ki);

      double b;
      b = k_i(kj);

      // std::cout << "row: " << ki << std::endl;
      // std::cout << "column: " << kj << std::endl;
      // std::cout << "row sum: " << a << std::endl;
      // std::cout << "col sum: " << b << std::endl;
      // std::cout << "intersection: " << it.value() << std::endl;
      // std::cout << "\n" << std::endl;


      it.valueRef() = it.value()/(a + (b - it.value()));
      if(it.value() < prune){
        it.valueRef() = 0;
      }
    }
  }
  SNN.prune(0.0); // actually remove pruned values
  return SNN;
}
