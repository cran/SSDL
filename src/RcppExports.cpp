// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// G_fun_cpp
arma::rowvec G_fun_cpp(arma::vec x, arma::vec y, arma::mat W);
RcppExport SEXP _SSDL_G_fun_cpp(SEXP xSEXP, SEXP ySEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(G_fun_cpp(x, y, W));
    return rcpp_result_gen;
END_RCPP
}
// gradient
arma::mat gradient(arma::mat S, arma::vec y, arma::mat W);
RcppExport SEXP _SSDL_gradient(SEXP SSEXP, SEXP ySEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(gradient(S, y, W));
    return rcpp_result_gen;
END_RCPP
}
// Rowsums_cpp
arma::vec Rowsums_cpp(arma::mat P, arma::vec SK);
RcppExport SEXP _SSDL_Rowsums_cpp(SEXP PSEXP, SEXP SKSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type SK(SKSEXP);
    rcpp_result_gen = Rcpp::wrap(Rowsums_cpp(P, SK));
    return rcpp_result_gen;
END_RCPP
}
// Rowsums_cpp_parallel
arma::vec Rowsums_cpp_parallel(arma::mat P, arma::vec SK);
RcppExport SEXP _SSDL_Rowsums_cpp_parallel(SEXP PSEXP, SEXP SKSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type SK(SKSEXP);
    rcpp_result_gen = Rcpp::wrap(Rowsums_cpp_parallel(P, SK));
    return rcpp_result_gen;
END_RCPP
}
// Sparse_prod
arma::mat Sparse_prod(arma::mat D, arma::mat Z);
RcppExport SEXP _SSDL_Sparse_prod(SEXP DSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(Sparse_prod(D, Z));
    return rcpp_result_gen;
END_RCPP
}
// Sparse_prod_row
arma::mat Sparse_prod_row(arma::mat D, arma::mat Z, bool CosSin);
RcppExport SEXP _SSDL_Sparse_prod_row(SEXP DSEXP, SEXP ZSEXP, SEXP CosSinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< bool >::type CosSin(CosSinSEXP);
    rcpp_result_gen = Rcpp::wrap(Sparse_prod_row(D, Z, CosSin));
    return rcpp_result_gen;
END_RCPP
}
// Sparse_prod_parallel
arma::mat Sparse_prod_parallel(arma::mat D, arma::mat Z);
RcppExport SEXP _SSDL_Sparse_prod_parallel(SEXP DSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(Sparse_prod_parallel(D, Z));
    return rcpp_result_gen;
END_RCPP
}
// Gradient_D_cpp_parallel
List Gradient_D_cpp_parallel(arma::mat D, arma::mat A, arma::mat W, arma::vec SK, bool ComputeGrad);
RcppExport SEXP _SSDL_Gradient_D_cpp_parallel(SEXP DSEXP, SEXP ASEXP, SEXP WSEXP, SEXP SKSEXP, SEXP ComputeGradSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type SK(SKSEXP);
    Rcpp::traits::input_parameter< bool >::type ComputeGrad(ComputeGradSEXP);
    rcpp_result_gen = Rcpp::wrap(Gradient_D_cpp_parallel(D, A, W, SK, ComputeGrad));
    return rcpp_result_gen;
END_RCPP
}
// Gradient_D_cpp
arma::mat Gradient_D_cpp(arma::mat D, arma::mat A, arma::mat W, arma::vec SK);
RcppExport SEXP _SSDL_Gradient_D_cpp(SEXP DSEXP, SEXP ASEXP, SEXP WSEXP, SEXP SKSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type SK(SKSEXP);
    rcpp_result_gen = Rcpp::wrap(Gradient_D_cpp(D, A, W, SK));
    return rcpp_result_gen;
END_RCPP
}
// ObjFun_COMP_cpp
double ObjFun_COMP_cpp(arma::vec d, arma::mat W, arma::vec residue);
RcppExport SEXP _SSDL_ObjFun_COMP_cpp(SEXP dSEXP, SEXP WSEXP, SEXP residueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type residue(residueSEXP);
    rcpp_result_gen = Rcpp::wrap(ObjFun_COMP_cpp(d, W, residue));
    return rcpp_result_gen;
END_RCPP
}
// Gradient_COMP_cpp
arma::rowvec Gradient_COMP_cpp(arma::vec d, arma::mat W, arma::vec residue);
RcppExport SEXP _SSDL_Gradient_COMP_cpp(SEXP dSEXP, SEXP WSEXP, SEXP residueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type residue(residueSEXP);
    rcpp_result_gen = Rcpp::wrap(Gradient_COMP_cpp(d, W, residue));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SSDL_G_fun_cpp", (DL_FUNC) &_SSDL_G_fun_cpp, 3},
    {"_SSDL_gradient", (DL_FUNC) &_SSDL_gradient, 3},
    {"_SSDL_Rowsums_cpp", (DL_FUNC) &_SSDL_Rowsums_cpp, 2},
    {"_SSDL_Rowsums_cpp_parallel", (DL_FUNC) &_SSDL_Rowsums_cpp_parallel, 2},
    {"_SSDL_Sparse_prod", (DL_FUNC) &_SSDL_Sparse_prod, 2},
    {"_SSDL_Sparse_prod_row", (DL_FUNC) &_SSDL_Sparse_prod_row, 3},
    {"_SSDL_Sparse_prod_parallel", (DL_FUNC) &_SSDL_Sparse_prod_parallel, 2},
    {"_SSDL_Gradient_D_cpp_parallel", (DL_FUNC) &_SSDL_Gradient_D_cpp_parallel, 5},
    {"_SSDL_Gradient_D_cpp", (DL_FUNC) &_SSDL_Gradient_D_cpp, 4},
    {"_SSDL_ObjFun_COMP_cpp", (DL_FUNC) &_SSDL_ObjFun_COMP_cpp, 3},
    {"_SSDL_Gradient_COMP_cpp", (DL_FUNC) &_SSDL_Gradient_COMP_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SSDL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
