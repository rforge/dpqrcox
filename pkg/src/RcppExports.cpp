
#include <Rcpp.h>

using namespace Rcpp;

NumericVector dpqr_cox(NumericVector);
RcppExport SEXP dpqrcox_dpqr_cox(SEXP Sxa) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::NumericVector xa(Sxa);
    __result = Rcpp::wrap(dpqr_cox(xa));
    return __result;
END_RCPP
}
