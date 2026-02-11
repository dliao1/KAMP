#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List kamp_pair_sums_trans_biv(IntegerVector i,
                                    IntegerVector j,
                                    NumericVector d,
                                    NumericVector w,
                                    NumericVector rvals,
                                    LogicalVector is_mark1,
                                    LogicalVector is_mark2, // adding 2nd mark vec for bivariate
                                    int npts) {

  int n_r = rvals.size();
  int npairs = d.size();

  NumericVector R0(n_r);
  NumericVector R1(n_r);
  NumericVector R2(n_r);
  NumericVector Ksum(n_r);
  NumericVector deg(npts, 0.0);

  double curr_R0 = 0; // R0 = sum(Wr)
  double curr_R1 = 0; // R1 = sum(Wr^2)
  double curr_R2 = 0;
  double curr_Ksum = 0;
  int pair_index = 0;

  for (int r_index = 0; r_index < n_r; r_index++) {
    double current_r = rvals[r_index];

    while (pair_index < npairs && d[pair_index] <= current_r) {
      double wk = w[pair_index];
      int ii = i[pair_index] - 1;
      int jj = j[pair_index] - 1;

      if (ii < 0 || ii >= npts || jj < 0 || jj >= npts)
        stop("Index out of bounds");

      curr_R0 += wk;
      curr_R1 += wk * wk;

      double prev_deg_i = deg[ii];
      deg[ii] += wk;
      curr_R2 += (deg[ii] * deg[ii] - prev_deg_i * prev_deg_i);

      // change for bivariate
      if (is_mark1[ii] && is_mark2[jj]) {
        curr_Ksum += wk;
      }

      pair_index++;
    }

    R0[r_index] = curr_R0;
    R1[r_index] = curr_R1;
    R2[r_index] = curr_R2 - curr_R1;
    Ksum[r_index] = curr_Ksum;
  }

  return List::create(
    _["R0"] = R0,
    _["R1"] = R1,
    _["R2"] = R2,
    _["Ksum"] = Ksum
  );
}
