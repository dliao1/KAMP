#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

static double log_choose(double a, double b) {
  if (b < 0 || b > a) return R_NegInf;
  return R::lgammafn(a + 1) - R::lgammafn(b + 1) - R::lgammafn(a - b + 1);
}

// [[Rcpp::export]]
List gamp_sums(IntegerVector nbr,
               NumericVector d,
               IntegerVector start,
               NumericVector rvals,
               int npts,
               int m1,
               int m2,
               bool univariate) {

  // Close pairs sorted by (i, d): point p's neighbors are
  // start[p] .. start[p + 1] - 1, nearest first
  // "mark2" below means mark1 for univariate G

  int n_r = rvals.size();
  double N = npts;

  // number of mark2 points among the others, given p (t1) or p and q (t2) are mark1
  int t1 = univariate ? m1 - 1 : m2;
  int t2 = univariate ? m1 - 2 : m2;

  // B1[k] = P(none of k points is mark2 | p is mark1)
  std::vector<double> B1(npts + 1, 0.0);
  double lden1 = log_choose(N - 1, t1);
  for (int k = 0; k <= npts; k++) {
    B1[k] = std::exp(log_choose(N - 1 - k, t1) - lden1);
  }

  // B2[k] = P(none of k points is mark2 | p and q are mark1)
  // cf = P(p and q are both mark1)
  bool pairs = m1 >= 2;
  std::vector<double> B2(2 * npts + 1, 0.0);
  double cf = 0.0;
  if (pairs) {
    double lden2 = log_choose(N - 2, t2);
    for (int k = 0; k <= 2 * npts; k++) {
      B2[k] = std::exp(log_choose(N - 2 - k, t2) - lden2);
    }
    cf = ((double) m1 / N) * ((double) (m1 - 1) / (N - 1));
  }

  NumericVector mu_G(n_r);
  NumericVector var_G(n_r);

  std::vector<int> deg(npts);
  std::vector<int> cnt(npts, 0);     // neighbors shared with p
  std::vector<char> isnb(npts, 0);   // is q a neighbor of p
  std::vector<int> touched;

  for (int r_index = 0; r_index < n_r; r_index++) {
    double current_r = rvals[r_index];

    int maxdeg = 0;
    for (int p = 0; p < npts; p++) {
      deg[p] = std::upper_bound(d.begin() + start[p], d.begin() + start[p + 1], current_r) -
               (d.begin() + start[p]);
      if (deg[p] > maxdeg) maxdeg = deg[p];
    }

    // mu_p = P(p is mark1 and has a mark2 neighbor within r)
    double sum_mu = 0.0, sum_mu2 = 0.0;
    for (int p = 0; p < npts; p++) {
      double mu_p = ((double) m1 / N) * (1.0 - B1[deg[p]]);
      sum_mu += mu_p;
      sum_mu2 += mu_p * mu_p;
    }
    mu_G[r_index] = sum_mu / m1;

    // Var(G) = (1/m1^2) * [sum mu_p (1 - mu_p) + 2 * sum_{p<q} (pi_pq - mu_p mu_q)]
    // pi_pq = P(p and q are mark1 and both have a mark2 neighbor)
    // First treat every pair as far apart (no shared neighbors), then correct
    // the pairs that are neighbors or share neighbors

    double sum_pi = 0.0;
    if (pairs) {
      std::vector<double> hist(maxdeg + 1, 0.0);
      double sumB2 = 0.0;
      for (int p = 0; p < npts; p++) {
        hist[deg[p]] += 1.0;
        sumB2 += B2[deg[p]];
      }
      double pairB2 = 0.0; // sum_{p<q} B2[n_p + n_q]
      for (int a = 0; a <= maxdeg; a++) {
        if (hist[a] == 0) continue;
        pairB2 += hist[a] * (hist[a] - 1) / 2.0 * B2[2 * a];
        for (int b = a + 1; b <= maxdeg; b++) {
          if (hist[b] == 0) continue;
          pairB2 += hist[a] * hist[b] * B2[a + b];
        }
      }
      sum_pi = cf * (N * (N - 1) / 2.0 - (N - 1) * sumB2 + pairB2);

      for (int p = 0; p < npts; p++) {
        int dp = deg[p];
        for (int k = start[p]; k < start[p] + dp; k++) isnb[nbr[k] - 1] = 1;

        // count neighbors shared by p and q (p -> u -> q)
        for (int k = start[p]; k < start[p] + dp; k++) {
          int u = nbr[k] - 1;
          for (int k2 = start[u]; k2 < start[u] + deg[u]; k2++) {
            int q = nbr[k2] - 1;
            if (q == p) continue;
            if (cnt[q] == 0) touched.push_back(q);
            cnt[q]++;
          }
        }

        // neighbor pairs (univariate: p and q are each other's mark1 neighbor, so pi = cf)
        for (int k = start[p]; k < start[p] + dp; k++) {
          int q = nbr[k] - 1;
          if (q <= p) continue;
          int dq = deg[q];
          double pi_pq = univariate ? cf :
            cf * (1.0 - B2[dp - 1] - B2[dq - 1] + B2[dp + dq - 2 - cnt[q]]);
          sum_pi += pi_pq - cf * (1.0 - B2[dp] - B2[dq] + B2[dp + dq]);
        }

        // non-neighbor pairs that share neighbors
        for (int q : touched) {
          if (q > p && !isnb[q]) {
            int dq = deg[q];
            sum_pi += cf * (B2[dp + dq - cnt[q]] - B2[dp + dq]);
          }
          cnt[q] = 0;
        }
        touched.clear();
        for (int k = start[p]; k < start[p] + dp; k++) isnb[nbr[k] - 1] = 0;
      }
    }

    double sum_mumu = (sum_mu * sum_mu - sum_mu2) / 2.0; // sum_{p<q} mu_p mu_q
    double v = ((sum_mu - sum_mu2) + 2.0 * (sum_pi - sum_mumu)) / ((double) m1 * m1);
    var_G[r_index] = v > 0 ? v : 0.0;
  }

  return List::create(
    _["mu_G"] = mu_G,
    _["var_G"] = var_G
  );
}
