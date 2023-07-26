#include <Rcpp.h>
#include <faiss.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

NumericMatrix faiss_test(){

    int d = 64;                            // dimension
    int nb = 100000;                       // database size
    int nq = 10000;                        // nb of queries
    float *xb = new float[d * nb];
    float *xq = new float[d * nq];
    for(int i = 0; i < nb; i++) {
        for(int j = 0; j < d; j++) xb[d * i + j] = drand48();
        xb[d * i] += i / 1000.;
    }
    for(int i = 0; i < nq; i++) {
        for(int j = 0; j < d; j++) xq[d * i + j] = drand48();
        xq[d * i] += i / 1000.;
    }
    
    faiss::IndexFlatL2 index(d);           // call constructor
    printf("is_trained = %s\n", index.is_trained ? "true" : "false");
    index.add(nb, xb);                     // add vectors to the index
    printf("ntotal = %ld\n", index.ntotal);
    
    int k = 4;
    {       // sanity check: search 5 first vectors of xb
        idx_t *I = new idx_t[k * 5];
        float *D = new float[k * 5];
        index.search(5, xb, k, D, I);
        printf("I=\n");
        for(int i = 0; i < 5; i++) {
            for(int j = 0; j < k; j++) printf("%5ld ", I[i * k + j]);
            printf("\n");
        }
        ...
        delete [] I;
        delete [] D;
    }
    {       // search xq
        idx_t *I = new idx_t[k * nq];
        float *D = new float[k * nq];
        index.search(nq, xq, k, D, I);
             ...
    }
    
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
 timesTwo(42)
*/
