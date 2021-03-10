// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, Rcpp)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <iostream>

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

//' @title G_fun_cpp
//' @description Supporting function for gradient computation in \code{\link{SSDL}[CDL]}. 
//' @param x is a vector of length \eqn{2m}
//' @param y is a vector of length \eqn{2m}
//' @param W is a frequency matrix \eqn{m \times s}
//' @return a vector of length \eqn{s}
//' @keywords internal
//' @export
// [[Rcpp::export]]
arma::rowvec G_fun_cpp(arma::vec x, arma::vec y, arma::mat W){
  int m = W.n_rows;
  arma::vec sk_y(m);
  for(int i = 0; i<m; i++){
    sk_y[i] = x[i]*y[m+i] - x[m+i]*y[i];
  }
  return (sk_y.t()*W);
}

struct G_fun_worker : public Worker {
  arma::mat& input;
  arma::vec& y;
  arma::mat& W;
  arma::mat& output;
  //rowvec G_fun_cpp_1(vec input, vec y, mat W);
  G_fun_worker(arma::mat& input, arma::vec& y, arma::mat& W, arma::mat& output) :
    input(input), y(y), W(W), output(output) {}
  void operator()(size_t begin, size_t end) {
    for(size_t i = begin; i< end; i++) output.col(i) = G_fun_cpp(input.col(i), y, W).t();
  }
};

//'@title Supporting function for the gradient computation.
//'@description The gradient computation using a parallel application of the function \code{\link{G_fun_cpp}} 
//'to a vector \code{y}, which is defined as a difference between a data sketch and 
//'the decomposition sketch, \emph{i.e.} \eqn{SK - SK(D\cdot A)}.
//'@seealso \code{\link{Gradient_D_cpp_parallel}}
//'@param S is a matrix
//'@param y is a vector
//'@param W is a frequency matrix
//'@return a gradient stored in the matrix \code{S}.
//'@keywords internal
//'@export
//[[Rcpp::export]]
arma::mat gradient(arma::mat S, arma::vec y,arma::mat W){
  int K = S.n_cols;
  int d = W.n_cols;
  arma::mat output(d, K);
  G_fun_worker worker(S, y, W, output);
  parallelFor(0, K, worker);
  return output;
}
//' @title Rowsums_cpp
//'@description Implementation of a matrix rows' summation and a sketch vector subtraction
//'@details The sums over the rows of the matrix \eqn{P \in R^{2m \time n}} is a sketch of the decomposition, 
//'\emph{i.e.} \eqn{SK(D\cdot A)}. It is defined as \eqn{P = [cos(W\cdot D \cdot A);sin(W\cdot D \cdot A)]}, where \eqn{W} is a frequency matrix, \eqn{D} is 
//'a dictionary and \eqn{A} is a code matrix. This function computes a vector of the difference between a 
//'data sketch and the decomposition sketch: \eqn{SK - SK(D\cdot A)}. 
//'@seealso \code{\link{Gradient_D_cpp_parallel}} 
//'@param P is a matrix
//'@param SK is a data sketch
//'@return a vector of length \eqn{2m}
//'@keywords internal
//'@export
// [[Rcpp::export]]
arma::vec Rowsums_cpp(arma::mat P, arma::vec SK) {
  int m = P.n_rows;
  double N = (double) P.n_cols;
  arma::vec SK_DA(m);
  for(int i = 0; i< m;i++){
    SK_DA(i) = sum(P.row(i));
  }
  SK_DA = SK_DA/N;
  SK_DA -= SK;
  return SK_DA;
}

struct Rowsums_worker : public Worker {
  arma::mat& input;
  arma::vec& SK;
  double N = (double)input.n_cols;
  arma::vec& output;

  Rowsums_worker(arma::mat& input,arma::vec& SK,double N, arma::vec& output) :
    input(input),SK(SK),N(N), output(output) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++) output(i) = sum(input.row(i))/N - SK(i);
  }
};
//'@title Rowsums_cpp_parallel
//'@description Parallel implementation of a matrix rows' summation and a sketch vector subtraction
//'@param P is a matrix
//'@param SK is a data sketch
//'@return a vector of matrix row sums minus a data sketch
//'@seealso \code{\link[SSDL]{Rowsums_cpp}}
//'@keywords internal
//'@export
// [[Rcpp::export]]
arma::vec Rowsums_cpp_parallel(arma::mat P, arma::vec SK) {
  int m = P.n_rows;
  double N = (double)P.n_cols;
  arma::vec SK_DA(m);
  Rowsums_worker A(P,SK, N, SK_DA);
  parallelFor(0, m, A);
  return SK_DA;
}

//'@title Sparse matrix product
//'@description Implementation of the product a matrix by a sparse matrix
//'
//'@param D is a matrix
//'@param Z is a sparse matrix
//'@return a product matrix
//'@keywords internal
//'@export
// [[Rcpp::export]]
arma::mat Sparse_prod(arma::mat D, arma::mat Z) {
  int n = D.n_rows;
  int N = Z.n_cols;
  arma::mat R(n, N);
  for(int i =0;i<N; i++){
    
    arma::vec x = Z.col(i);
    arma::uvec ids = find(x > 0); // Find indices
    
    R.col(i)=D.cols(ids)* x.elem(ids);
  }
  return R;
}
//' @title Sparse product by rows
//'@description Implementation sparse matrix product by rows
//'
//'@param D is a matrix
//'@param Z is a matrix
//'@param CosSin indicates whether to compute cos and sin of the resulting matrix
//'@return a matrix
//'@keywords internal
//'@export
// [[Rcpp::export]]
arma::mat Sparse_prod_row(arma::mat D, arma::mat Z,  bool CosSin = false) {
  int n = D.n_rows;
  int K = Z.n_rows;
  arma::mat R(n, K);

    arma::mat P(2*n, K);
    arma::mat R_real(n, K);
    arma::mat R_im(n,K);

  for(int i =0;i<K; i++){
    
    arma::rowvec x = Z.row(i);
    arma::uvec ids = find(x > 0); // Find indices
    
    R.col(i)= D.cols(ids)*x.elem(ids);
    if(CosSin){
      R_real.col(i) = cos(R.col(i));
      R_im.col(i) = sin(R.col(i));
    }
  }
  P = join_cols(R_real, R_im);
  if(CosSin){
    return P;
  }else{
    return R;
  }
  
}
struct Sparse_prod_worker : public Worker {
  arma::mat& D;
  arma::mat& Z;
  arma::mat& output;

  Sparse_prod_worker(arma::mat& D, arma::mat& Z, arma::mat& output) :
    D(D),Z(Z), output(output) {

  }

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      arma::vec x = Z.col(i);
      arma::uvec ids = find(x > 0); // Find indices
      
      output.col(i)= D.cols(ids)*x.elem(ids);
    }
  }
};
//'@title Sparse matrix product
//'@description Parallel computation of a product of a matrix by a sparse matrix
//'@param D is a matrix
//'@param Z is a matrix
//'@return a product matrix
//'@keywords internal
//'@export
// [[Rcpp::export]]
arma::mat Sparse_prod_parallel(arma::mat D, arma::mat Z) {
  int m = D.n_rows;
  int N = Z.n_cols;
  arma::mat output(m, N);
  Sparse_prod_worker Prod(D, Z, output);
  parallelFor(0, N, Prod);
  return output;
}
//'@title Gradient_D_cpp_parallel
//'@description Parallel computation of the gradient with respect to a dictionary matrix and
//' the objective function computation.  
//' @details The objective function is given as \eqn{\|SK - SK(D\cdot A)\|^2}, where \eqn{SK} is 
//' a data sketch, \eqn{A = \{\alpha_1, \dots, \alpha_n\}} is a code matrix and \eqn{SK(D\cdot A)} 
//' denotes a decomposition sketch, which is defined as: 
//' \eqn{SK(D\cdot A) = \frac{1}{n}\left[\sum_{i=1}^n \cos(W\cdot D \cdot \alpha_i), 
//' \sum_{i=1}^n \sin(W\cdot D \cdot \alpha_i)\right]}.
//' The gradient of this objective function with respect to a dictionary element \eqn{d_l \in R^{s}} is given as:
//' \eqn{- 2 \left( \nabla_{d_l} SK(D\cdot A) \right)^{\top} \cdot r},
//' where \eqn{r = SK - SK(D \cdot A)}, \eqn{\frac{\delta}{\delta d_l} SK^j(D\cdot A) = 1i \cdot \left( \frac{1}{n} 
//' \sum_{i = 1}^n  A_{li}\cdot \prod_{k=1}^K SK^j(A_{ki}\cdot d_k) \right)\cdot w_j^{\top}},
//' and \eqn{SK^j(\cdot )} is the \eqn{j^{th}} coordinate of the sketch vector.
//'@param D is a dictionary \eqn{s \times K}.
//'@param A is a code matrix \eqn{K \times n}.
//'@param W is a frequency matrix \eqn{m \times s} with frequency vectors in matrix rows.
//'@param SK is a data sketch. It is a \eqn{2m}-dimensional vector.
//'@param ComputeGrad indicates whether to compute the gradient or only the objective function value
//'@return a list
//'\itemize{
//'   \item \code{grad} is a computed gradient
//'   \item \code{ObjFun} is a objective function value
//'   \item \code{diff} is a vector of the difference between the data sketch and the decomposition sketch
//'} 
//'@examples
//' RcppParallel::setThreadOptions(numThreads = 2)
//' X = matrix(abs(rnorm(n = 1000)), ncol = 100, nrow = 10)
//' X_fbm = bigstatsr::as_FBM(X)$save()
//' W = chickn::GenerateFrequencies(Data = X_fbm, m = 64, N0 = ncol(X_fbm),
//'                                 ncores = 1, niter= 3, nblocks = 2, sigma_start = 0.001)$W
//' SK= chickn::Sketch(X_fbm, W)
//' D = X_fbm[, sample(ncol(X_fbm), 10)]
//' A = sapply(sample(ncol(X_fbm), 5), function(i){
//'     as.vector(glmnet::glmnet(x = D, y = X_fbm[,i],
//'               lambda = 0, intercept = FALSE, lower.limits = 0)$beta)})
//' G = Gradient_D_cpp_parallel(D, A, W, SK)$grad                                    
//'@export
//[[Rcpp::export]]
List Gradient_D_cpp_parallel(arma::mat D, arma::mat A, arma::mat W, arma::vec SK, 
                             bool ComputeGrad = true) {
  int m = W.n_rows;
  int N = A.n_cols;
  int K = D.n_cols;
  int d = D.n_rows;
 // clock_t tStart = clock();
  arma::mat prodWD = W*D;
  arma::mat prodWDA = Sparse_prod_parallel(prodWD,A);

  //printf("Prod WDA: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  //tStart = clock();
  arma::mat P = join_cols(cos(prodWDA), sin(prodWDA));
  // NumericMatrix P2(2*m,N);
  // for(int i=0;i<N;i++){
  //   for(int j=0; j<m;j++){
  //     P2(j,i)= cos(prodWDA(j,i));
  //     P2(j+m, i) =  sin(prodWDA(j,i));
  //   }
  //   
  // }
 
  //printf("P cos sin: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);


   arma::vec diff(2*m);
  // tStart = clock();
  // diff = Rowsums_cpp_parallel(P,SK);
  // printf("Row sum: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  
  arma::vec diff2(2*m);
  NumericMatrix P2 = wrap(P);
  //tStart = clock();
  diff2 = rowSums(P2);
  diff = diff2/double(N) - SK;
  //printf("Row sum: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

  List Output;
  if (ComputeGrad){//CHECK
    arma::mat grad(d, K);
    
   // tStart = clock();
    arma::mat S = Sparse_prod_parallel(P,A.t());
    //printf("Sparse prod P AT: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

    //tStart = clock(); 
    grad = gradient(S, diff, W);
    //printf("Gradient: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    Output["grad"] = grad*2./(double)N;
 
    
  }
  Output["ObjFun"] = as_scalar(diff.t()*diff);
  Output["diff"] = diff;
  return Output;
}

//'@title Gradient D
//'@description Gradient with respect to dictionary using sparse product
//'@param D is a dictionary \eqn{s \times K}.
//'@param A is a code matrix \eqn{K \times n}.
//'@param W is a frequency matrix \eqn{m \times s} with frequency vectors in matrix rows.
//'@param SK is a data sketch. It is a \eqn{2m}-dimensional vector.
//'@return a gradient matrix
//'@keywords internal
//'@export
// [[Rcpp::export]]
arma::mat Gradient_D_cpp(arma::mat D, arma::mat A, arma::mat W, arma::vec SK) { // m, z, W, K, N
  int m = W.n_rows;
  int N = A.n_cols;
  int K = D.n_cols;
  int d = D.n_rows;
  arma::mat prodWD = W*D;
  arma::mat prodWDA = Sparse_prod(prodWD,A);//prodWD*A;
  arma::mat P = join_cols(cos(prodWDA), sin(prodWDA));
  arma::vec SK_DA(2*m);
  for(int i = 0; i< 2*m;i++){
    SK_DA(i) = sum(P.row(i));
  }
  
  // NumericMatrix P2 = wrap(P);
  // SK_DA = rowSums(P2);
  // 
  SK_DA = SK_DA/(double)N;
  arma::vec diff = SK_DA - SK;
  arma::mat S = Sparse_prod(P,A.t());//P*A.t();
  arma::mat grad(d, K);

  for(int i = 0; i< K; i++){
    grad.col(i) = G_fun_cpp(S.col(i), diff, W).t();
  }
  return grad*2./(double)N;
}
//'@title COMP objective function  
//'@description Computation of the objective function from the Compressive Orthogonal Matching Pursuit algorithm.
//'@details The objective function of the Compressive Orthogonal Matching Pursuit is defined as:
//'\eqn{-\frac{SK(d)\cdot r}{\|SK(d)\|}}, where \eqn{SK(d)} denotes a 
//'sketch of the dictionary element \code{d} and \eqn{r} is the residue vector, which is defined as the difference between
//'the data sketch \eqn{SK} and the weighted sum of the dictionary elements' sketches, \emph{i.e.}
//'\eqn{SK - \sum_{i=1}^K \beta_i \cdot SK(d_i)}.
//'This function is involved in \code{\link{COMP_initialization}} routine. 
//'@inheritParams Gradient_D_cpp_parallel
//'@param d is a dictionary element 
//'@param residue is a residue vector.
//'@return an objective function value
//'@examples
//' X = matrix(abs(rnorm(n = 1000)), ncol = 100, nrow = 10)
//' X_fbm = bigstatsr::as_FBM(X)$save()
//' W = chickn::GenerateFrequencies(Data = X_fbm, m = 64, N0 = ncol(X_fbm),
//'                                 ncores = 1, niter= 3, nblocks = 2, sigma_start = 0.001)$W
//' SK= chickn::Sketch(X_fbm, W)
//' D = X_fbm[, sample(ncol(X_fbm), 10)]
//' weights = sample(10, 10)/10
//' SK_D = rbind(cos(W%*%D), sin(W%*%D))
//' d = D[,1]
//' r = SK - SK_D%*%weights
//' OF = ObjFun_COMP_cpp(d, W, r)
//' @export
//'@seealso \code{\link{COMP_initialization}}, \code{\link{Gradient_COMP_cpp}}
// [[Rcpp::export]]
double ObjFun_COMP_cpp(arma::vec d, arma::mat W, arma::vec residue){
  int m = W.n_rows;
  arma::vec prodWd = W*d;
  arma::vec SK_d(2*m);
  SK_d.subvec(0, m-1) = cos(prodWd);
  SK_d.subvec(m, 2*m-1) = sin(prodWd);
  double norm_SK_d = as_scalar(sqrt(SK_d.t()*SK_d));
  SK_d = SK_d/norm_SK_d;
  return as_scalar(-residue.t()*SK_d);
}

//' @title COMP Gradient
//' @description The gradient of the objective function from the Compressive Orthogonal Matching Pursuit
//' with respect to a dictionary element. 
//' @details \code{Gradient_COMP_cpp} computes the gradient of the objective function 
//' \eqn{OF(d) = -\frac{SK(d)\cdot r}{\|SK(d)\|}}, where \eqn{SK(d)} denotes a 
//' sketch of the dictionary element \code{d} and \eqn{r} is the residue vector. 
//' The gradient is given as \eqn{\nabla_d OF(d) = \frac{-G(SK(d), y, W)}{\|SK(d)\|}}, where 
//' a vector \eqn{y  = r-\left(r^{\top} \cdot SK(d)\right)\cdot SK(d)} and 
//' a function \eqn{G(x, y, W)} is
//' given as: 
//' \eqn{G(x,y, W) = \left(x[1:m]\odot y[m+1:2m] - x[m+1]\odot y[1:m]\right)^{\top}\cdot W},
//' where \eqn{\odot} denotes an element-wise vector multiplication.   
//' @inheritParams ObjFun_COMP_cpp
//' @return a gradient vector
//' @examples
//' X = matrix(abs(rnorm(n = 1000)), ncol = 100, nrow = 10)
//' X_fbm = bigstatsr::as_FBM(X)$save()
//' W = chickn::GenerateFrequencies(Data = X_fbm, m = 64, N0 = ncol(X_fbm),
//'                                 ncores = 1, niter= 3, nblocks = 2, sigma_start = 0.001)$W
//' SK= chickn::Sketch(X_fbm, W)
//' D = X_fbm[, sample(ncol(X_fbm), 10)]
//' weights = sample(10, 10)/10
//' SK_D = rbind(cos(W%*%D), sin(W%*%D))
//' d = D[,1]
//' r = SK - SK_D%*%weights
//' Grad = Gradient_COMP_cpp(d, W, r)
//' @export
//' @seealso \code{\link{ObjFun_COMP_cpp}}, \code{\link{COMP_initialization}}
// [[Rcpp::export]]
arma::rowvec Gradient_COMP_cpp(arma::vec d, arma::mat W, arma::vec residue){
  int m = W.n_rows;
  arma::vec prodWd = W*d;
  arma::vec SK_d(2*m);
  SK_d.subvec(0, m-1) = cos(prodWd);
  SK_d.subvec(m, 2*m-1) = sin(prodWd);
  //vec SK_d =  join_cols(cos(prodWd), sin(prodWd));
  double norm_SK_d = as_scalar(sqrt(SK_d.t()*SK_d));
  SK_d = SK_d/norm_SK_d;
  double val = as_scalar(-residue.t()*SK_d);
  arma::vec y = val*SK_d + residue;
  arma::rowvec G  = G_fun_cpp(SK_d, y, W);
  return -G/norm_SK_d;
}
