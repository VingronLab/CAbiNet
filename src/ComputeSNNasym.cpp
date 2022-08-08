
#include <RcppEigen.h>
#include <string>
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
//' @param  mode The type of neighboring vertices to use for calculating similarity
//'  scores(Jaccard Index). Three options: "out", "in" and "all":
//' * "out": Select neighbouring vertices by out-going edges;
//' * "in": Selecting neighbouring vertices by in-coming edges;
//' * "all": Selecting neigbouring vertices by both in-coming and out-going edges.
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> ComputeSNNasym(Eigen::SparseMatrix<double> SNN, 
                                           double prune,
                                           String mode) {

  Eigen::VectorXd k_i(SNN.rows());
  // Eigen::VectorXd k_j(SNN.cols());
  
  if (mode == "out"){
    k_i = SNN * Eigen::VectorXd::Ones(SNN.cols());
    SNN = SNN * (SNN.transpose());
  }else if (mode == "in"){
    k_i = SNN.transpose() * Eigen::VectorXd::Ones(SNN.rows());
    SNN =  (SNN.transpose()) * SNN;
  }else if (mode == "all"){
    Eigen::SparseMatrix<double> sym;
    sym = SNN + (Eigen::SparseMatrix<double> (SNN.transpose()));
    
    for (int i=0; i < sym.outerSize(); ++i){
      for (Eigen::SparseMatrix<double>::InnerIterator it(sym, i); it; ++it){
        int ki = it.row();
        int kj = it.col();
        
        // if (SNN(ki, kj) + SNN(kj,ji) >0){
        it.valueRef() = 1;
        // }
      }
    }
    
    SNN = sym * (sym.transpose());
    k_i = sym * Eigen::VectorXd::Ones(sym.cols());
    // std::cout << "all" << std::endl
    }
  for (int i=0; i < SNN.outerSize(); ++i){  //number of columns ?
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){  // Iterate over rows

      int ki = it.row();
      int kj = it.col();

      double a;
      a = k_i(ki);

      double b;
      b = k_i(kj);

      it.valueRef() = it.value()/(a + (b - it.value()));
      if(it.value() < prune){
        it.valueRef() = 0;
      }
    }
  }
  SNN.prune(0.0); // actually remove pruned values
  return SNN;
}

