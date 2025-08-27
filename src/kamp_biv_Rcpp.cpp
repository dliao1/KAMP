#include <RcppArmadillo.h>
#include <cmath>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List kamp_biv_Rcpp(const arma::mat& Wr,
                           const arma::ivec& marks, 
                           const int mark1_code,
                           const int mark2_code,
                           const double areaW) {

  const uword n = Wr.n_rows;

  const double R0 = accu(Wr); //sum(Wr), total weights
  const double R1 = accu(Wr % Wr); //squared weights
  const vec rowsum = sum(Wr, 1); // returns n x 1, each entry = total weight of all neighbors for point i
  const double R2 = accu(rowsum % rowsum) - R1; // square of row sums
  const double R3 = R0*R0 - 2.0*R1 - 4.0*R2;

  uvec idx1 = find(marks == mark1_code);
  uvec idx2 = find(marks == mark2_code);

  const double m1 = (double) (idx1.n_elem);
  const double m2 = (double) (idx2.n_elem);

  if (m1 == 0 || m2 == 0) {
    stop("Selected marks have zero count.");
  }

  const double Knum = accu(Wr(idx1, idx2));

  const double K     = areaW * Knum / (m1 * m2);
  const double nd    = (double) (n);
  const double mu_K  = areaW * (R0 / (nd * (nd - 1.0)));

  const double f1 = (m1 * m2) / (nd * (nd - 1.0));
  const double f2 = f1 * (m1 + m2 - 2.0) / (nd - 2.0);
  const double f3 = f1 * (m1 - 1.0) * (m2 - 1.0) / ((nd - 2.0) * (nd - 3.0));

  const double var_K = areaW*areaW * (R1*f1 + R2*f2 + R3*f3) / (m1*m1*m2*m2) - mu_K*mu_K;

  double z = NA_REAL, p = NA_REAL;
  if (var_K > 0) {
    z = (K - mu_K) / std::sqrt(var_K);
    p = pnorm(-z, 0.0, 1.0, /*lower_tail=*/ 1, /*log_p=*/ 0);
  }

  return List::create(
    _["K"]      = K,
    _["mu_K"]   = mu_K,
    _["var_K"]  = var_K,
    _["z"]      = z,
    _["p"]      = p
  );
}
