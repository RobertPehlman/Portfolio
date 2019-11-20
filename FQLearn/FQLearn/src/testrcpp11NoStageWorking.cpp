#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include<iostream>
#include<armadillo>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <vector>
using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]

arma::vec newProjectionNoStage(arma::vec input, int p, int L){
  mat gamma0(input);
  gamma0 = gamma0.rows(0,(p * L - 1));
  gamma0.reshape(p, L); // MAY NEED TO ALTER THIS
  
  
  mat U0;
  vec s0;
  mat V0;
  svd(U0,s0,V0,gamma0);
  V0 = V0.cols(0,p-1);
  
  mat ortho0 = U0 * V0.t();
  //ortho0 = ortho0.t();
  ortho0.reshape(L * p,1);
  
  
  mat alphas(input);
  alphas = alphas.rows((L*p), ((L+4)*p-1));
  
  mat variances(input);
  variances = variances.rows(((L+4)*p), ((L+5)*p + 2));
  int varlength = variances.n_elem;
  bool status1 = all(vectorise(variances) >= 0.001);
  if(!status1){
    for(int i = 0; i < varlength; i++){
      variances(i,0) = std::max(variances(i,0), 0.05);
    }
  }
  
  variances((p + 2), 0) = std::max(variances((p + 2),0), 0.2);

  mat both = join_vert(ortho0,alphas);
  
  both = join_vert(both,variances);
  vec proj = vectorise(both);
  return(proj);
}


// [[Rcpp::export]]

double newLikelihoodNoStage(arma::vec input, arma::uvec indices, double rho, int p, int L, arma::vec m, std::vector<arma::mat> B_T0, std::vector<arma::vec> data, std::vector<arma::vec> indexVec){
  
  int length = indices.n_rows;
  int indexcounter = 0;
  double likelihood = 0;
  double det = 0;
  double kernel = 0;
  
  vec variances(input);
  
  mat gamma0(input);
  gamma0 = gamma0.rows(0,(p * L - 1));
  gamma0.reshape(p,L); // MAY NEED TO ALTER THIS

  vec varpars0 = variances.rows((p*(L+4)), (p*(L+5) - 1));
  mat sigma0 = square(diagmat(varpars0));
  
  while(indexcounter < length){
    int i = indices(indexcounter);
    indexcounter++;

    mat tempB_T0(L,m(i));
    tempB_T0.zeros();
    
    for(int j = 0; j < L; j++){
      for(int l = 0; l < m(i); l++){
        tempB_T0 = B_T0[i];
      }
    }
    
    mat tau(m(i),m(i));
    tau.eye();
    tau = tau * input((p*(L+5) + 1)) * input((p*(L+5) + 1)); // tau matrix
    mat phi0 = gamma0*B_T0[i];
    mat B1 = phi0.t();
    mat Omega = B1*sigma0*B1.t() + tau;
 
    mat cholO(m(i),m(i));
    
    bool success = false;
    int countFails = 0;
    while(success == false)
    {
      success = chol(cholO, Omega);
      
      if(success == false)
      {
        countFails++;
        Omega += eye(Omega.n_rows,Omega.n_rows) * 5 * 1e-6;
      }
      if(countFails > 5){
        return(1000000);
      }
    }
    
    mat rooti = trans(inv(trimatu(cholO)));
    
    det = 0;
    
    for(int j=0; j<m(i); j++){
      det+=log(rooti(j,j));
    }
    
    kernel = 0;
    vec z = rooti * (data[i]);
    kernel = sum(z%z);
    likelihood += .5 * kernel - det + m(i) / 2 * log(2*M_PI);
  }
  
  mat absvals = abs(gamma0);
  double l1 = accu(absvals);
  likelihood += rho * l1;
  
  return(likelihood);
}


// [[Rcpp::export]]

arma::vec newGradientNoStage(arma::vec input, arma::uvec indices, double rho, int p, int L, arma::vec m, std::vector<arma::mat> B_T0, std::vector<arma::vec> data, std::vector<arma::vec> indexVec){
  
  int length = indices.n_rows;
  int indexcounter = 0;

  vec variances(input);
  
  mat gamma0(input);
  gamma0 = gamma0.rows(0,(p * L - 1));
  gamma0.reshape(p,L); // MAY NEED TO ALTER THIS
  //  gamma0 = gamma0.t();
  
  
  vec varpars0 = variances.rows((p*(L+4)), (p*(L+5) - 1));
  mat lambda0 = square(diagmat(varpars0));
  
  mat gamma0Gradient(L,p);
  gamma0Gradient.zeros();
  vec lambda0Gradient(p);
  lambda0Gradient.zeros();

  double tauGradient = 0;
  
  while(indexcounter < length){
    int i = indices(indexcounter);
    indexcounter++;

    mat tempB_T0(L,m(i));
    tempB_T0.zeros();
    tempB_T0 = B_T0[i];

    
    mat tauMat(m(i),m(i));
    tauMat.eye();
    double tau = input((p*(L+5) + 1));
    tauMat = tauMat * tau * tau; // tau matrix

    mat O_T = B_T0[i];
    
    mat phi0 = (gamma0*B_T0[i]).t();
    
    vec W_A = (data[i]);
    mat B1 = phi0;
    mat Omega = B1*lambda0*B1.t() + tauMat;
    

    mat invOmega(m(i),m(i));
    bool success = false;
    int countFails = 0;
    while(success == false){
      success = inv_sympd(invOmega, Omega);
      
      if(success == false)
      {
        Omega += eye(Omega.n_rows,Omega.n_rows) * 1e-5;
        countFails++;
      }
      
    }
    
    
    
    mat dOmega = .5 *(invOmega - invOmega * W_A * W_A.t() * invOmega);
    mat phi0Lambda = phi0 * lambda0;
    
    
    
    gamma0Gradient -= 2 * (O_T * dOmega * phi0Lambda);
    

    lambda0Gradient -= diagvec(B1.t() * dOmega * B1 * 2 * diagmat(varpars0));
    
    tauGradient -= trace(dOmega) * 2 * tau;
  
  }
  

  gamma0Gradient = gamma0Gradient.t();
  mat l1gradient0 = sign(gamma0);
  gamma0Gradient = gamma0Gradient - rho * l1gradient0;
  mat RiemannianGradient0 = gamma0Gradient - .5 * gamma0 * (gamma0.t() * gamma0Gradient + gamma0Gradient.t() * gamma0);
  RiemannianGradient0.reshape((L * p),1);
  
  
  
  
  mat lambda0GradientMat(lambda0Gradient);
  
  
  mat tauGradientMat(1,1);
  tauGradientMat(0,0) = tauGradient;
  
  // Placeholder matrices to avoid changing indices
  
  mat alphaGradientMat((4*p),1);
  alphaGradientMat.zeros();
  mat nuGradientMat(1,1);
  nuGradientMat.zeros();
  mat hGradientMat(1,1);
  hGradientMat.zeros();
  
  mat fullGradientMat = join_vert(RiemannianGradient0, alphaGradientMat);
  fullGradientMat = join_vert(fullGradientMat, lambda0GradientMat);
  fullGradientMat = join_vert(fullGradientMat, nuGradientMat);
  fullGradientMat = join_vert(fullGradientMat, tauGradientMat);
  fullGradientMat = join_vert(fullGradientMat, hGradientMat);

  vec fullGradient = -1 * vectorise(fullGradientMat);
  return (fullGradient);
}



// [[Rcpp::export]]

arma::vec newDescentFullOnlyNoStage(arma::mat fulldata, arma::mat Obasis, arma::vec input, double rho, int p, int iter2, arma::vec stepsizes, double tol){
  int T = Obasis.n_rows;
  int L = Obasis.n_cols;
  //int p = 5;
  int n = fulldata.n_cols;
  Obasis = Obasis.t();
  
  vector<vec> data (n);
  vector<vec> indices (n);
  vector<mat> B_T0 (n);
  vec m(n);
  m.zeros();
  
  for(int j = 0; j < n; j++){
    for(int i = 0; i<T; i++){
      if(!std::isnan(fulldata(i,j)))
        m(j)=m(j)+1;
    }
  }
  
  int indexcounter=0;
  
  for(int j = 0; j < n; j++){
    vec temp(m(j));
    temp.zeros();
    vec tempind(m(j));
    tempind.zeros();
    for(int i = 0; i<T; i++){
      if(!std::isnan((double)fulldata(i,j))){
        temp(indexcounter) = fulldata(i,j);
        tempind(indexcounter) = i;
        indexcounter++;
      }
    }
    data[j] = temp;
    indices[j] = tempind;
    indexcounter = 0;
  }
  
  indexcounter = 0;
  
  for(int k = 0; k < n; k++){
    
    mat tempB_T0(L,m(k));
    tempB_T0.zeros();

    
    for(int j = 0; j < L; j++){
      for(int i = 0; i < m(k); i++){
        tempB_T0(j,i) = Obasis(j,indices[k](i));
        
      }
      indexcounter = 0;
    }
    B_T0[k] = tempB_T0;
  }
  
  uvec allindices = linspace<uvec>(0,(n-1), n);
  
  double rho0 = rho;
 
  vec both = input;
  
  both = newProjectionNoStage(both, p, L);
  
  
  int n_stepsizes = stepsizes.n_rows;
  vec values ((n_stepsizes * n_stepsizes));
  
  double fOld = newLikelihoodNoStage(both, allindices, rho0, p, L, m, B_T0, data, indices);
  cout << "First Likelihood =" << fOld << "\n";
  for(int i = 0; i < iter2; i++){
    values.ones();
    values = values * 1000000;
    
    vec grad = newGradientNoStage(both, allindices, rho0, p, L, m, B_T0, data, indices);
    grad = grad / n;
    
    for(int j = 0; j < n_stepsizes-3; j++){
      int upperk = std::min(j+10, n_stepsizes);
      int lowerk = std::min(j+1, n_stepsizes);
      for(int k = lowerk; k < upperk; k++){
        mat gradGamma(grad);
        gradGamma = gradGamma.rows(0,(p*L - 1));
        gradGamma = gradGamma * stepsizes(j);
        mat gradNonOrthos(grad);
        gradNonOrthos = gradNonOrthos.rows((p*L), ((L+5)*p + 2));
        gradNonOrthos = gradNonOrthos * stepsizes(k);
        mat temp = join_vert(gradGamma,gradNonOrthos);
        vec increment = vectorise(temp);
        values((j*n_stepsizes + k)) = newLikelihoodNoStage(newProjectionNoStage(both - increment, p, L), allindices, rho0, p, L, m, B_T0, data, indices);  //CHANGE BACK TO sampL !!!
      }
    }
    
    int minpos1 = 0;
    int minpos2 = 0;
    double min = values(0);
    
    for(int j = 0; j < n_stepsizes; j++){
      for(int k = 0; k < n_stepsizes; k++){
        if(values((j*n_stepsizes + k)) < min){
          min = values((j*n_stepsizes + k));
          minpos1 = j;
          minpos2 = k;
        }
      }
    }
    
    mat gradGamma(grad);
    gradGamma = gradGamma.rows(0,(p*L - 1));
    gradGamma = gradGamma * stepsizes(minpos1);
    mat gradNonOrthos(grad);
    gradNonOrthos = gradNonOrthos.rows((p*L), ((L+5)*p + 2));
    gradNonOrthos = gradNonOrthos * stepsizes(minpos2);
    mat matincrement = join_vert(gradGamma,gradNonOrthos);
    vec increment = vectorise(matincrement);
    
    vec temp = newProjectionNoStage(both - increment, p, L);
    both = temp;
    values.zeros();
    
    if((i % 5) == 0){

      double fNew = newLikelihoodNoStage(both, allindices, rho0, p, L, m, B_T0, data, indices);
      cout << fNew << "\n";
      if(abs(fNew - fOld) < tol){
        break;
      }
      fOld = fNew;
    }
    
    
  }
  vec estAndLikeVec(both);
  estAndLikeVec.insert_rows(((L+5)*p + 3),1);
  estAndLikeVec(((L+5)*p + 3)) = newLikelihoodNoStage(both, allindices, rho0, p, L, m, B_T0, data, indices);
  return (estAndLikeVec);
  return (both);
}
