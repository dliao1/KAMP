#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List kamp_pair_sums_trans(IntegerVector i,
                          IntegerVector j,
                          NumericVector d,
                          NumericVector w,
                          NumericVector rvals,
                          LogicalVector is_mark1,
                          int npts) {

  // i =
  // j =
  // d = distance for each pair
  // w = weight for each pair
  // rvals = vector of radius values at which to compute sums
  // is_mark1 = logical vector indicating whether each point has mark 1
  // npts = total number of points

  int n_r = rvals.size(); // num radii
  int npairs = d.size(); // num pairs

  NumericVector R0(n_r);
  NumericVector R1(n_r);
  NumericVector R2(n_r);
  NumericVector Ksum(n_r);
  NumericVector deg(npts, 0.0); // sum of all weights of pairs involving current point

  double curr_R0 = 0; // R0 = sum(Wr)
  double curr_R1 = 0; // R1 = sum(Wr^2)
  double curr_R2 = 0; // R2 = sum(rowSums(Wr)^2))
  double curr_Ksum = 0;
  int pair_index = 0;

  for (int r_index = 0; r_index < n_r; r_index++) { // loop over all radii - where main speedup compared to R comes from
    double current_r = rvals[r_index];

    while (pair_index < npairs && d[pair_index] <= current_r) { // while curr pair num < total num of pairs and distance at that index <= current r
      double wk = w[pair_index]; // weight at current pair
      int ii = i[pair_index] - 1; // converts R index to C index
      int jj = j[pair_index] - 1;

      if (ii < 0 || ii >= npts || jj < 0 || jj >= npts) {
        stop("Index out of bounds");
      }

      curr_R0 += wk; // sums all R0
      curr_R1 += wk * wk;

      // building sum(rowSums(Wr)^2))
      // ** NOTE THAT Wr IS ASYMMETRIC SO ONLY ONE ROW SUM IS UPDATED PER PAIR **

      // Proof:
      // deg[ii] = a
      // R2_old = a^2 + sum(deg(other_pts)^2)
      // R2_new = (a + wk)^2 + sum(deg(curr_pt)^2)
      // R2_new - R2_old = (a + wk)^2 - a^2
      // so R2_new = R2_old + difference between new and old squared degree of current point ii

      double prev_deg_i = deg[ii]; // stores previous degree of point ii
      deg[ii] += wk; //updates sum of weights for current point  ii

      curr_R2 += (deg[ii]*deg[ii] - prev_deg_i*prev_deg_i);

      // Filters which pairs contribute to K for the current radius
      if (is_mark1[ii] && is_mark1[jj]) {
        curr_Ksum += wk;
      }

      pair_index++; //not resetting pair index, just adding new ones as r increases
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
