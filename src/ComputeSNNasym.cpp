
#include <RcppEigen.h>
#include <string>
#include <Eigen/Dense>
using namespace Rcpp;

// NOTE: Set sparse matrix to float?

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
// ' @param  mode The type of neighboring vertices to use for calculating similarity
// '  scores(Jaccard Index). Three options: "out", "in" and "all":
// ' * "out": Select neighbouring vertices by out-going edges;
// ' * "in": Selecting neighbouring vertices by in-coming edges;
// ' * "all": Selecting neigbouring vertices by both in-coming and out-going edges.
// ' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> ComputeSNNasym_sparse(Eigen::Map<Eigen::SparseMatrix<int>> SNN,
                                                  double prune,
                                                  String mode) {
    Eigen::VectorXi k_i(SNN.rows());
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>  res_dense;

    if (mode == "out"){
        k_i = SNN * Eigen::VectorXi::Ones(SNN.cols());
        res_dense = SNN * (SNN.transpose());
    }
    else if (mode == "in"){

        k_i = SNN.transpose() * Eigen::VectorXi::Ones(SNN.rows());
        res_dense =  SNN.transpose() * SNN;
    }else if (mode == "all"){

        Eigen::SparseMatrix<int> sym(SNN.rows(), SNN.cols());

        Eigen::SparseMatrix<int> transposed = SNN.transpose();
        // it is proved that initializing this transposed matrix here is better than initializing it earlier.
        sym = SNN + transposed;

        for (int i=0; i < sym.outerSize(); ++i){

            for (Eigen::SparseMatrix<int>::InnerIterator it(sym, i); it; ++it){

                it.valueRef() = 1;
            }
        }

        res_dense = sym * (sym.transpose()); // the product of two matrices might not be sparse anymore. Storing the results into a sparse matrix might raise up bad_alloc error
        k_i = sym * Eigen::VectorXi::Ones(sym.cols());
    }

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> trp;
    double overlapping;
    double a;
    double b;
    double val;

    for (int i = 0; i < res_dense.cols(); ++i){  //number of columns ?

        for (int j = 0; j < res_dense.rows(); ++j){  // Iterate over rows

            a = (double)k_i(j);
            b = (double)k_i(i);
            val = (double)res_dense(j, i);

            if ((a + (b - val)) != 0 ){
                overlapping = val/(a + (b - val));
            }else{
                overlapping = 0.0;
            }
            if(overlapping >= prune){
                trp.push_back(Trip(j,
                               i,
                               overlapping));;
            }
        }
    }

    Eigen::SparseMatrix<double> res(res_dense.rows(), res_dense.cols());
    res.setFromTriplets(trp.begin(), trp.end());

    return res;
}

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
Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>
ComputeSNNasym_dense(Eigen::Map<Eigen::SparseMatrix<int>> SNN,
                     float prune,
                     String mode) {

  Eigen::VectorXi k_i(SNN.rows());
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> res_dense;

  if (mode == "out") {
    k_i = SNN * Eigen::VectorXi::Ones(SNN.cols());
    res_dense = SNN * (SNN.transpose());
  } else if (mode == "in") {

    k_i = SNN.transpose() * Eigen::VectorXi::Ones(SNN.rows());
    res_dense = SNN.transpose() * SNN;
  } else if (mode == "all") {

    Eigen::SparseMatrix<int> sym(SNN.rows(), SNN.cols());

    Eigen::SparseMatrix<int> transposed = SNN.transpose();

    // it is proven that initializing this transposed matrix here is better than
    // initializing it earlier.
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

  Eigen::MatrixXf res_dense_flt = res_dense.cast<float>();

  for (int i = 0; i < res_dense_flt.rows(); ++i) {   // number of rows
    for (int j = 0; j < res_dense_flt.cols(); ++j) { // Iterate over cols

      float a;
      a = (float)k_i(i);

      float b;
      b = (float)k_i(j);

      float val = res_dense_flt(i, j);
      float ov = val / (a + (b - val));

      if (ov < prune) {
        res_dense_flt.coeffRef(i, j) = 0.0f;
      } else {
        res_dense_flt.coeffRef(i, j) = ov;
      }
    }
  }
  return res_dense_flt;
}
