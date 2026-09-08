#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List kamp_pair_sums_bord(IntegerVector i,
                         IntegerVector j,
                         NumericVector d,
                         NumericVector rvals,
                         LogicalVector is_mark1,
                         NumericVector bdist,
                         int npts) {

  // i, j, d are pairs sorted by distance, bdist is bdist.points() per point.
  //
  // border weight is just 1, but eligibility (bdist[i] >= r) depends on r,
  // and shrinks as r grows -- a pair counted at small r can drop out at a
  // bigger r. trans/iso can just accumulate as r sweeps up since their
  // weights don't change, but border can't, so we redo the sums per r here.

  int n_r = rvals.size();
  int npairs = d.size();

  NumericVector R0(n_r);
  NumericVector R1(n_r);
  NumericVector R2(n_r);
  NumericVector Ksum(n_r);

  std::vector<double> deg(npts);

  for (int r_index = 0; r_index < n_r; r_index++) {
    double current_r = rvals[r_index];

    std::fill(deg.begin(), deg.end(), 0.0);
    double curr_R0 = 0.0;
    double curr_Ksum = 0.0;

    for (int p = 0; p < npairs; p++) {
      if (d[p] > current_r) break; // pairs sorted by distance; none further qualify

      int ii = i[p] - 1; // converts R index to C index
      int jj = j[p] - 1;

      if (ii < 0 || ii >= npts || jj < 0 || jj >= npts) {
        stop("Index out of bounds");
      }

      if (bdist[ii] < current_r) continue; // reference point ii not eligible at this r

      curr_R0 += 1.0;   // weight is always 1 for border correction
      deg[ii] += 1.0;   // rowSums(Wr) for point ii

      if (is_mark1[ii] && is_mark1[jj]) {
        curr_Ksum += 1.0;
      }
    }

    double sum_deg2 = 0.0;
    for (int k = 0; k < npts; k++) {
      sum_deg2 += deg[k] * deg[k];
    }

    R0[r_index] = curr_R0;
    R1[r_index] = curr_R0; // weight is 1, so w^2 sums equal w sums
    R2[r_index] = sum_deg2 - curr_R0; // sum(rowSums(Wr)^2) - R1
    Ksum[r_index] = curr_Ksum;
  }

  return List::create(
    _["R0"] = R0,
    _["R1"] = R1,
    _["R2"] = R2,
    _["Ksum"] = Ksum
  );
}
