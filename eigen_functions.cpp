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

uword get_k(const int& i, const int& ii, const double& I) 
{
     // Precompute reusable terms
     double term1 = I * (I - 1.0) / 2.0;  // Total number of pairs
     double term2 = (I - i) * (I - i - 1.0) / 2.0;  // Remaining pairs after row i
     return static_cast<uword>(term1 - term2 + ii - i - 1.0);
}

void tunning(double& del2, int& n_tun, const double& mix_rate, const int& b) 
{
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
double loglik(const double& I, const double& zeta, const mat& U, 
              const rowvec& Lambda, const vec& Y) 
{
     double out = 0.0;
     double eta, prob;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               eta = zeta + accu(U.row(i)%U.row(ii)%Lambda);
               prob = R::pnorm(eta, 0.0, 1.0, 1, 0);
               out += R::dbinom(Y[get_k(i, ii, I)], 1, prob, 1);
          }
     }
     
     return out;
}

// [[Rcpp::export]]
double lfcd_U(const rowvec& x, const uword& i, const double& I, 
              const double& K, const double& sigsq, const double& zeta, 
              const mat& U, const rowvec& Lambda, const vec& Y) 
{
     double out = -std::pow(arma::norm(x), 2.0) / (2.0 * sigsq);
     double eta, prob;
     
     // Loop over upper-triangle pairs where i < ii
     if (i < I - 1) {
          for (uword ii = i + 1; ii < I; ii++) {
               eta = zeta + accu(x%U.row(ii)%Lambda);
               prob = R::pnorm(eta, 0.0, 1.0, 1, 0); // Standard normal CDF
               out += R::dbinom(Y[get_k(i, ii, I)], 1, prob, 1);
          }
     }
     
     // Loop over lower-triangle pairs where ii < i
     if (i > 0) {
          for (uword ii = 0; ii < i; ii++) {
               eta = zeta + accu(x%U.row(ii)%Lambda);
               prob = R::pnorm(eta, 0.0, 1.0, 1, 0); // Standard normal CDF
               out += R::dbinom(Y[get_k(ii, i, I)], 1, prob, 1);
          }
     }
     
     return out;
}

// [[Rcpp::export]]
List sample_U(const double& b, int n_tun_U, double del2_U, int n_U, 
              const int& n_burn, const double& I, const double& K, 
              const double& sigsq, const double& zeta, mat U, 
              const rowvec& Lambda, const vec& Y) 
{
     // Metropolis-Hastings step
     rowvec u_c(K), u_p(K);
     double lfcd_c, lfcd_p, acceptance_prob;
     
     for (uword i = 0; i < I; i++) {
          u_c = U.row(i);
          u_p = u_c + std::sqrt(del2_U) * arma::randn<rowvec>(K);
          
          // Compute log-likelihoods only once per proposal
          lfcd_c = lfcd_U(u_c, i, I, K, sigsq, zeta, U, Lambda, Y);
          lfcd_p = lfcd_U(u_p, i, I, K, sigsq, zeta, U, Lambda, Y);
          
          acceptance_prob = std::exp(lfcd_p - lfcd_c);
          
          if (R::runif(0, 1) < acceptance_prob) {
               U.row(i) = u_p;
               n_U++;
          }
     }
     
     // Tuning step (only during burn-in)
     if (b < n_burn) {
          double mix_rate = static_cast<double>(n_U) / (b * I);
          tunning(del2_U, n_tun_U, mix_rate, b);
     }
     
     return List::create(Named("U")       = U,
                         Named("del2_U")  = del2_U,
                         Named("n_U")     = n_U,
                         Named("n_tun_U") = n_tun_U);
}

// [[Rcpp::export]]
double lfcd_Lambda(const double& x, const uword& k, const double& I, 
                   const double& kapsq, const double& zeta, const mat& U, 
                   rowvec Lambda, const vec& Y) 
{
     double out = -std::pow(x, 2.0) / (2.0 * kapsq);
     
     // Update Lambda[k] with the new proposal
     Lambda[k] = x;
     
     double eta, prob;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               eta = zeta + accu(U.row(i)%U.row(ii)%Lambda);
               prob = R::pnorm(eta, 0.0, 1.0, 1, 0); // Standard normal CDF
               out += R::dbinom(Y[get_k(i, ii, I)], 1, prob, 1);
          }
     }
     
     return out;
}

// [[Rcpp::export]]
List sample_Lambda(const double& b, int n_tun_Lambda, double del2_Lambda, int n_Lambda, 
                   const int& n_burn, const double& I, const double& K, 
                   const double& kapsq, const double& zeta, const mat& U, 
                   rowvec Lambda, const vec& Y) 
{
     double lambda_c, lambda_p, lfcd_c, lfcd_p, acceptance_ratio;
     
     for (int k = 0; k < K; ++k) {
          lambda_c = Lambda[k];
          lambda_p = R::rnorm(lambda_c, std::sqrt(del2_Lambda));
          
          // Compute log-likelihoods only once per proposal
          lfcd_c = lfcd_Lambda(lambda_c, k, I, kapsq, zeta, U, Lambda, Y);
          lfcd_p = lfcd_Lambda(lambda_p, k, I, kapsq, zeta, U, Lambda, Y);
          
          acceptance_ratio = std::exp(lfcd_p - lfcd_c);
          
          if (R::runif(0, 1) < acceptance_ratio) {
               Lambda[k] = lambda_p;
               n_Lambda++;
          }
     }
     
     // Adaptive tuning during burn-in
     if (b < n_burn) {
          double mix_rate = static_cast<double>(n_Lambda) / (b * K);
          tunning(del2_Lambda, n_tun_Lambda, mix_rate, b);
     }
     
     return List::create(Named("Lambda")       = Lambda,
                         Named("del2_Lambda")  = del2_Lambda,
                         Named("n_Lambda")     = n_Lambda,
                         Named("n_tun_Lambda") = n_tun_Lambda);
}

// [[Rcpp::export]]
double lfcd_zeta(const double& x, const double& I, const double& omesq, 
                 const mat& U, const rowvec& Lambda, const vec& Y) 
{
     double out = -std::pow(x, 2.0) / (2.0 * omesq);
     double eta, prob;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               eta = x + accu(U.row(i)%U.row(ii)%Lambda);
               prob = R::pnorm(eta, 0.0, 1.0, 1, 0); // Standard normal CDF
               out += R::dbinom(Y[get_k(i, ii, I)], 1, prob, 1);
          }
     }
     
     return out;
}

// [[Rcpp::export]]
List sample_zeta(const double& b, int n_tun_zeta, double del2_zeta, int n_zeta, 
                 const int& n_burn, const double& I, const double& omesq, 
                 double zeta, const mat& U, const rowvec& Lambda, 
                 const vec& Y) 
{
     // Metropolis-Hastings step
     double zeta_p = R::rnorm(zeta, std::sqrt(del2_zeta));
     
     // Compute log-likelihoods only once per proposal
     double lfcd_c = lfcd_zeta(zeta, I, omesq, U, Lambda, Y);
     double lfcd_p = lfcd_zeta(zeta_p, I, omesq, U, Lambda, Y);
     
     double acceptance_ratio = std::exp(lfcd_p - lfcd_c);
     
     if (R::runif(0, 1) < acceptance_ratio) {
          zeta = zeta_p;
          n_zeta++;
     }
     
     // Adaptive tuning during burn-in
     if (b < n_burn) {
          double mix_rate = static_cast<double>(n_zeta) / b;
          tunning(del2_zeta, n_tun_zeta, mix_rate, b);
     }
     
     return List::create(Named("zeta")       = zeta,
                         Named("del2_zeta")  = del2_zeta,
                         Named("n_zeta")     = n_zeta,
                         Named("n_tun_zeta") = n_tun_zeta);
}

// [[Rcpp::export]]
double sample_sigsq(const double& I, const double& K, const double& a_sig, 
                    const double& b_sig, const mat& U) 
{
     double shape = a_sig + 0.5 * K * I;
     double rate  = b_sig + 0.5 * arma::accu(arma::square(U)); 
     return 1.0 / R::rgamma(shape, 1.0 / rate);
}

// [[Rcpp::export]]
double sample_kapsq(const double& K, const double& a_kap, const double& b_kap, 
                    const rowvec& Lambda) 
{
     double shape = a_kap + 0.5 * K;
     double rate  = b_kap + 0.5 * arma::accu(arma::square(Lambda));
     return 1.0 / R::rgamma(shape, 1.0 / rate);
}

// [[Rcpp::export]]
double sample_omesq(const double& a_ome, const double& b_ome, const double& zeta) 
{
     double shape = a_ome + 0.5;
     double rate  = b_ome + 0.5 * std::pow(zeta, 2.0);
     return 1.0 / R::rgamma(shape, 1.0 / rate);
}

// [[Rcpp::export]]
vec sample_Y(const double& I, const double& zeta, const mat& U, 
                   const rowvec& Lambda, const vec& na_indices, 
                   vec Yna) 
{
     uword k;
     double eta, prob;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               k = get_k(i, ii, I);
               if (na_indices[k]) {
                    eta  = zeta + accu(U.row(i)%U.row(ii)%Lambda);
                    prob = R::pnorm(eta, 0.0, 1.0, 1, 0); // Standard normal CDF
                    Yna[k] = R::rbinom(1, prob);
               }
          }
     }
     
     return Yna;
}

// [[Rcpp::export]]
List WAIC(const double& I, const double& K, const double& B, const vec& Y, 
          const vec& zeta_chain, const mat& U_chain, const mat& Lambda_chain) 
{
     double lppd = 0.0, slp = 0.0, pWAIC2 = 0.0;
     rowvec u_i(K), u_ii(K), lam(K);
     double tmp, a_ib, a_ib_ssq, a_ib_sum, eta, prob;
     uword m;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               m = get_k(i, ii, I);
               tmp = 0.0;
               a_ib_ssq = 0.0;
               a_ib_sum = 0.0;
               
               for (uword b = 0; b < B; ++b) {
                    double zeta = zeta_chain[b];
                    
                    // Efficiently extract U and Lambda values for sample `b`
                    for (uword k = 0; k < K; ++k) {
                         u_i[k]  = U_chain.at(b, k * I + i);
                         u_ii[k] = U_chain.at(b, k * I + ii);
                         lam[k]  = Lambda_chain.at(b, k);
                    }
                    
                    // Compute eta and probability using the standard normal CDF
                    eta  = zeta + arma::dot(u_i % u_ii, lam);
                    prob = R::pnorm(eta, 0.0, 1.0, 1, 0);
                    a_ib = R::dbinom(Y[m], 1, prob, 1);
                    
                    // Compute WAIC components
                    tmp += std::exp(a_ib) / B;
                    slp += a_ib / B;
                    a_ib_ssq += std::pow(a_ib, 2);
                    a_ib_sum += a_ib;
               }
               
               lppd   += std::log(tmp);
               pWAIC2 += (a_ib_ssq - B * std::pow(a_ib_sum / B, 2)) / (B - 1.0);
          }
     }
     
     // Compute WAIC measures
     double pWAIC1 =  2.0 * lppd - 2.0 * slp;
     double waic1  = -2.0 * lppd + 2.0 * pWAIC1;
     double waic2  = -2.0 * lppd + 2.0 * pWAIC2;
     
     return List::create(Named("lppd")   = lppd,
                         Named("pWAIC1") = pWAIC1,
                         Named("pWAIC2") = pWAIC2,
                         Named("waic1")  = waic1,
                         Named("waic2")  = waic2);
}

// [[Rcpp::export]]
mat simulate_data(const double& I, const double& zeta, 
                  const mat& U, const rowvec& Lambda) 
{
     mat Y(I, I, arma::fill::zeros);
     double eta, prob;
     
     for (uword i = 0; i < I - 1; i++) {
          for (uword ii = i + 1; ii < I; ii++) {
               eta  = zeta + accu(U.row(i)%U.row(ii)%Lambda);
               prob = R::pnorm(eta, 0.0, 1.0, 1, 0); // Standard normal CDF
               Y(i, ii) = R::rbinom(1, prob);
               Y(ii, i) = Y(i, ii);  // Ensure symmetry
          }
     }
     
     return Y;
}
