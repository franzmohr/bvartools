// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026 Franz X. Mohr

#include "kalman_durbin_koopman_2002.h"
#include "bayests/arma.h"

arma::mat kalman_durbin_koopman_2002(arma::mat &y, arma::mat &z,
                                     arma::mat sigma_u, arma::mat sigma_v,
                                     arma::mat B,
                                     arma::vec &a_init, arma::mat &P_init) {
  
  arma::uword k = y.n_rows;
  int t = y.n_cols;
  arma::uword nvars = z.n_cols;
  
  arma::mat sigma_u_temp = arma::zeros<arma::mat>(k * t, k);
  if (sigma_u.n_rows == k){
    for (int i = 0; i < t; i++){
      sigma_u_temp.rows(i * k, (i + 1) * k - 1) = sigma_u;
    }
    sigma_u = sigma_u_temp;
  }
  
  arma::mat sigma_v_temp = arma::zeros<arma::mat>(nvars * t, nvars);
  if (sigma_v.n_rows == nvars){
    for (int i = 0; i < t; i++){
      sigma_v_temp.rows(i * nvars, (i + 1) * nvars - 1) = sigma_v;
    }
    sigma_v = sigma_v_temp;
  }
  
  arma::mat B_temp = arma::zeros<arma::mat>(nvars * t, nvars);
  if (B.n_rows == nvars){
    for (int i = 0; i < t; i++){
      B_temp.rows(i * nvars, (i + 1) * nvars - 1) = B;
    }
    B = B_temp;
  }
  
  arma::mat yplus = y * 0;
  arma::mat aplus = arma::zeros<arma::mat>(nvars, t + 1);
  
  arma::mat U;
  arma::vec s;
  
  // Step 1
  arma::eig_sym(s, U, P_init);
  arma::mat A = U * arma::diagmat(sqrt(s)) * arma::trans(U);
  aplus.col(0) = A * arma::randn<arma::vec>(nvars); // cf. Jarocinski (2015)
  
  int p1, p2, pA1, pA2;
  for (int i = 0; i < t; i++){
    p1 = i * k;
    p2 = (i + 1) * k - 1;
    pA1 = i * nvars;
    pA2 = (i + 1) * nvars - 1;
    
    arma::eig_sym(s, U, sigma_u.rows(p1, p2));
    A = U * arma::diagmat(sqrt(s)) * arma::trans(U);
    yplus.col(i) = z.rows(p1, p2) * aplus.col(i) + A * arma::randn<arma::vec>(k);
    
    arma::eig_sym(s, U, sigma_v.rows(pA1, pA2));
    A = U * arma::diagmat(sqrt(s)) * arma::trans(U);
    aplus.col(i + 1) = B.rows(pA1, pA2) * aplus.col(i) + A * arma::randn<arma::vec>(nvars);
  }
  
  // Kalman filtering
  arma::mat ystar = y - yplus;
  arma::mat a = arma::zeros<arma::mat>(nvars, t + 1);
  a.col(0) = a_init; // cf. DK (2002, p.606)
  arma::mat P = P_init;
  arma::mat v = y * 0;
  arma::mat Fi = arma::zeros<arma::mat>(k * t, k);
  arma::mat K = arma::zeros<arma::mat>(nvars * t, k);
  arma::mat L = arma::zeros<arma::mat>(nvars * t, nvars);
  for (int i = 0; i < t ; i++){
    p1 = i * k;
    p2 = (i + 1) * k - 1;
    pA1 = i * nvars;
    pA2 = (i + 1) * nvars - 1;
    v.col(i) = ystar.col(i) - z.rows(p1, p2) * a.col(i);
    Fi.rows(p1, p2) = arma::inv(z.rows(p1, p2) * P * arma::trans(z.rows(p1, p2)) + sigma_u.rows(p1, p2));
    K.rows(pA1, pA2) = B.rows(pA1, pA2) * P * arma::trans(z.rows(p1, p2)) * Fi.rows(p1, p2);
    L.rows(pA1, pA2) = B.rows(pA1, pA2) - K.rows(pA1, pA2) * z.rows(p1, p2);
    a.col(i + 1) = B.rows(pA1, pA2) * a.col(i) + K.rows(pA1, pA2) * v.col(i);
    P = B.rows(pA1, pA2) * P * arma::trans(L.rows(pA1, pA2)) + sigma_v.rows(pA1, pA2);
  }
  
  // Backward smoothing
  arma::mat r = arma::zeros<arma::mat>(nvars, t);
  for (int i = (t - 1); i > 0; i--){
    r.col(i - 1) = arma::trans(z.rows(i * k, (i + 1) * k - 1)) * Fi.rows(i * k, (i + 1) * k - 1) * v.col(i) + arma::trans(L.rows(i * nvars, (i + 1) * nvars - 1)) * r.col(i);
  }
  arma::vec r0 = arma::trans(z.rows(0, k - 1)) * Fi.rows(0, k - 1) * v.col(0) + arma::trans(L.rows(0, nvars - 1)) * r.col(0);
  
  a.col(0) = a_init + P_init * r0;
  for (int i = 0; i < t; i++){
    a.col(i + 1) = B.rows(i * nvars, (i + 1) * nvars - 1) * a.col(i) + sigma_v.rows(i * nvars, (i + 1) * nvars - 1) * r.col(i);
  }
  
  // Step 3
  return a + aplus;
}
