#include <RcppArmadillo.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
#include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace R;
using namespace Rcpp;
using namespace std;

uword get_k(const int& i, const int& ii, const double& I) {
     // Precompute reusable terms
     double term1 = I * (I - 1.0) / 2.0;  // Total number of pairs
     double term2 = (I - i) * (I - i - 1.0) / 2.0;  // Remaining pairs after row i
     return static_cast<uword>(term1 - term2 + ii - i - 1.0);
}

void tunning(double& del2, int& n_tun, const double& mix_rate, const int& b) {
     // Declare constants within the function
     const double target_rate = 0.35;
     const double tolerance = 0.05;
     const double max_del2 = 10.0;
     const double min_del2 = 1e-10;
     
     if (b % n_tun == 0) {
          double diff = mix_rate - target_rate;
          if (std::abs(diff) > tolerance) {
               del2 *= std::exp(0.1 * diff);
               del2 = std::clamp(del2, min_del2, max_del2);
               n_tun = 100;
          } else {
               n_tun += 100;
          }
     }
}

// [[Rcpp::export]]
double loglik(const double& I, const double& zeta, const mat& U, const vec& Y) {
     double out = 0.0;
     uword k;
     double distance, prob;
     
     // Precompute number of pairs
     for (uword i = 0; i < I - 1; ++i) {
          const rowvec& row_i = U.row(i); // Cache row i for efficiency
          for (uword ii = i + 1; ii < I; ++ii) {
               k = get_k(i, ii, I);
               distance = norm(row_i - U.row(ii));
               prob = R::pnorm(zeta - distance, 0.0, 1.0, true, false);
               out += R::dbinom(Y[k], 1, prob, 1);
          }
     }
     
     return out;
}

double lfcd_U(const rowvec& x, const uword& i, const double& I, const double& K, const double& sigsq, 
              const double& zeta, const mat& U, const vec& Y) {
     double out = -pow(norm(x), 2.0) / (2.0 * sigsq);
     double distance, prob;
     uword k;
     
     // Precompute distances and log-likelihood contributions for i < I - 1
     if (i < I - 1) {
          for (uword ii = i + 1; ii < I; ++ii) {
               distance = norm(x - U.row(ii));
               prob = R::pnorm(zeta - distance, 0.0, 1.0, true, false); // Use the CDF of the standard normal
               k = get_k(i, ii, I);
               out += R::dbinom(Y[k], 1, prob, 1);
          }
     }
     
     // Precompute distances and log-likelihood contributions for i > 0
     if (i > 0) {
          for (uword ii = 0; ii < i; ++ii) {
               distance = norm(x - U.row(ii));
               prob = R::pnorm(zeta - distance, 0.0, 1.0, true, false); // Use the CDF of the standard normal
               k = get_k(ii, i, I);
               out += R::dbinom(Y[k], 1, prob, 1);
          }
     }
     
     return out;
}

// [[Rcpp::export]]
List sample_U(const double& b, int n_tun_U, double del2_U, int n_U, const int& n_burn, 
              const double& I, const double& K, const double& sigsq, const double& zeta, mat U, const vec& Y) {
     double sqrt_del2_U = sqrt(del2_U);  // Precompute square root of del2_U
     rowvec u_c(K), u_p(K);             // Predeclare row vectors to avoid reallocation
     double current_lfcd, proposed_lfcd, accept_prob;
     
     // Metropolis step
     for (uword i = 0; i < I; ++i) {
          u_c = U.row(i);                // Cache current row
          u_p = u_c + sqrt_del2_U * randn<rowvec>(K);  // Generate proposed row
          
          // Compute log-likelihood contributions
          current_lfcd = lfcd_U(u_c, i, I, K, sigsq, zeta, U, Y);
          proposed_lfcd = lfcd_U(u_p, i, I, K, sigsq, zeta, U, Y);
          
          // Compute acceptance probability
          accept_prob = exp(proposed_lfcd - current_lfcd);
          
          // Accept or reject the proposal
          if (R::runif(0, 1) < accept_prob) {
               U.row(i) = u_p;  // Update row if accepted
               n_U++;           // Increment accepted updates
          }
     }
     
     // Adjust tuning parameter during burn-in
     if (b < n_burn) {
          double mix_rate = static_cast<double>(n_U) / (b * I);
          tunning(del2_U, n_tun_U, mix_rate, b);
     }
     
     // Return updated values
     return List::create(Named("U") = U,
                         Named("del2_U") = del2_U,
                         Named("n_U") = n_U,
                         Named("n_tun_U") = n_tun_U);
}

double lfcd_zeta(const double& x, const double& I, const double& omesq, const mat& U, const vec& Y) {
     double out = -pow(x, 2.0) / (2.0 * omesq);  // Prior term
     double distance, prob;
     uword k;
     
     // Loop over all unique pairs of nodes
     for (uword i = 0; i < I - 1; ++i) {
          const rowvec& row_i = U.row(i);  // Cache row i
          for (uword ii = i + 1; ii < I; ++ii) {
               distance = norm(row_i - U.row(ii));
               prob = R::pnorm(x - distance, 0.0, 1.0, true, false); // Use the CDF of the standard normal
               k = get_k(i, ii, I);
               out += R::dbinom(Y[k], 1, prob, 1);
          }
     }
     
     return out;
}


// [[Rcpp::export]]
List sample_zeta(const double& b, int n_tun_zeta, double del2_zeta, int n_zeta, const int& n_burn, 
                 const double& I, const double& omesq, double zeta, const mat& U, const vec& Y) {
     double sqrt_del2_zeta = sqrt(del2_zeta);  // Precompute square root of del2_zeta
     double zeta_p = R::rnorm(zeta, sqrt_del2_zeta);  // Generate proposal
     double current_lfcd = lfcd_zeta(zeta, I, omesq, U, Y);  // Current likelihood
     double proposed_lfcd = lfcd_zeta(zeta_p, I, omesq, U, Y);  // Proposed likelihood
     double accept_prob = exp(proposed_lfcd - current_lfcd);  // Acceptance probability
     
     // Metropolis step: Accept or reject the proposal
     if (R::runif(0, 1) < accept_prob) {
          zeta = zeta_p;  // Accept proposal
          n_zeta++;       // Increment accepted updates
     }
     
     // Adjust tuning parameter during burn-in
     if (b < n_burn) {
          double mix_rate = static_cast<double>(n_zeta) / b;
          tunning(del2_zeta, n_tun_zeta, mix_rate, b);
     }
     
     // Return updated values
     return List::create(Named("zeta") = zeta,
                         Named("del2_zeta") = del2_zeta,
                         Named("n_zeta") = n_zeta,
                         Named("n_tun_zeta") = n_tun_zeta);
}

// [[Rcpp::export]]
double sample_sigsq(const double& I, const double& K, const double& a_sig, const double& b_sig, const mat& U) {
     // Precompute reused values
     double shape = a_sig + 0.5 * K * I;
     double scale = 1.0 / (b_sig + 0.5 * accu(pow(U, 2))); // Accumulate squared elements efficiently
     return 1.0 / R::rgamma(shape, scale);
}

// [[Rcpp::export]]
double sample_omesq(const double& a_ome, const double& b_ome, const double& zeta) {
     // Precompute reused values
     double shape = a_ome + 0.5;
     double scale = 1.0 / (b_ome + 0.5 * pow(zeta, 2)); // Simplify scaling
     return 1.0 / R::rgamma(shape, scale);
}

// [[Rcpp::export]]
vec sample_Y(const double& I, const double& zeta, const mat& U, const vec& na_indices, vec Yna) {
     uword k;
     double distance, prob;
     
     // Loop over unique pairs of indices
     for (uword i = 0; i < I - 1; ++i) {
          const rowvec& row_i = U.row(i); // Cache row i
          for (uword ii = i + 1; ii < I; ++ii) {
               k = get_k(i, ii, I);
               if (static_cast<bool>(na_indices[k])) {
                    distance = norm(row_i - U.row(ii));
                    prob = R::pnorm(zeta - distance, 0.0, 1.0, true, false);
                    Yna[k] = R::rbinom(1, prob);
               }
          }
     }
     
     return Yna;
}

// [[Rcpp::export]]
List WAIC(const double& I, const double& K, const double& B, const vec& Y, 
          const vec& zeta_chain, const mat& U_chain) {
     
     uword i, ii, m, b;
     double tmp, sum_a_ib_sq, sum_a_ib;
     double zeta;
     rowvec u_i(K), u_ii(K);
     double lppd = 0.0, slp = 0.0, pWAIC2 = 0.0;
     
     for (i = 0; i < I - 1; i++) {
          for (ii = i + 1; ii < I; ii++) {
               m = get_k(i, ii, I);
               tmp = 0.0;
               sum_a_ib_sq = 0.0;
               sum_a_ib = 0.0;
               
               for (b = 0; b < B; b++) {
                    zeta = zeta_chain[b];
                    
                    for (uword k = 0; k < K; k++) {
                         u_i[k] = U_chain.at(b, k * I + i);
                         u_ii[k] = U_chain.at(b, k * I + ii);
                    }
                    
                    double a_ib = R::dbinom(Y[m], 1, 
                                            R::pnorm(zeta - norm(u_i - u_ii), 0.0, 1.0, true, false), 1);
                    
                    // WAIC computations
                    tmp += exp(a_ib) / B;
                    slp += a_ib / B;
                    sum_a_ib_sq += pow(a_ib, 2);
                    sum_a_ib += a_ib;
               }
               
               lppd += log(tmp);
               pWAIC2 += (sum_a_ib_sq - B * pow(sum_a_ib / B, 2)) / (B - 1.0);
          }
     }
     
     double pWAIC1 = 2.0 * (lppd - slp);
     double waic1 = -2.0 * (lppd - pWAIC1);
     double waic2 = -2.0 * (lppd - pWAIC2);
     
     return List::create(Named("lppd") = lppd,
                         Named("pWAIC1") = pWAIC1,
                         Named("pWAIC2") = pWAIC2,
                         Named("waic1") = waic1,
                         Named("waic2") = waic2);
}

// [[Rcpp::export]]
rowvec interaction_probs0(const double& I, const double& K, const double& B, const vec& zeta_chain, const mat& U_chain) {
     rowvec out(0.5 * I * (I - 1.0), fill::zeros);
     double zeta, distance;
     rowvec u_i(K), u_ii(K);
     
     // Loop through chains
     for (uword b = 0; b < B; ++b) {
          zeta = zeta_chain[b];
          for (uword i = 0; i < I - 1; ++i) {
               for (uword ii = i + 1; ii < I; ++ii) {
                    u_i = U_chain.row(b).cols(K * i, K * (i + 1) - 1);
                    u_ii = U_chain.row(b).cols(K * ii, K * (ii + 1) - 1);
                    
                    // Compute interaction probability using CDF of the standard normal
                    distance = norm(u_i - u_ii);
                    out[get_k(i, ii, I)] += R::pnorm(zeta - distance, 0.0, 1.0, true, false) / B;
               }
          }
     }
     
     return out;
}

// [[Rcpp::export]]
mat simulate_data(const double& I, const double& zeta, const mat& U) {
     mat Y(I, I, fill::zeros);
     
     for (uword i = 0; i < I - 1; ++i) {
          const rowvec& u_i = U.row(i);  // Cache row i for efficiency
          for (uword ii = i + 1; ii < I; ++ii) {
               Y(i, ii) = R::rbinom(1, R::pnorm(zeta - norm(u_i - U.row(ii)), 0.0, 1.0, true, false));
               Y(ii, i) = Y(i, ii);  // Ensure symmetry
          }
     }
     
     return Y;
}
