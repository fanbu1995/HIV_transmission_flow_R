// Generalized Spatial Fusion Model Framework for 
// Joint Analysis of Point and Areal Data, Appendix B

// Please cite: Wang, C., Puhan, M.A., and Furrer. R, (2017) Generalized Spatial Fusion Model Framework for Joint Analysis of Point and Areal Data, Spatial Statistics.

// Contains the Stan model used to fit our proposed fusion model
// Authors: Craig Wang, Milo Puhan, Reinhard Furrer (reinhard.furrer@math.uzh.ch).	
 
// The nngp_w_lpdf function is modified based on the original NNGP implementation shared on Stan discussion forum by Lu Zhang

functions{
  real nngp_w_lpdf(vector w_r, real sigmasq, real phi, matrix neardist, matrix neardistM,
                   int[,] nearind, int n, int M){
    vector[n] V;
    vector[n] Uw;
    vector[n] w;
    int dim;
    int h;
    real out;
    Uw = w_r;
    
    for (i in 2:n) {
      matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M] temp_neardistM;
      matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M] L;
      vector[ i < (M + 1)? (i - 1) : M] u;
      vector[ i < (M + 1)? (i - 1) : M] v;
      row_vector[i < (M + 1)? (i - 1) : M] v2;
      dim = (i < (M + 1))? (i-1) : M;
      
      // get exp(-phi * neardistM)
      if(dim == 1){temp_neardistM[1, 1] = 1;}
      else{
        h = 0;
        for (j in 1:(dim - 1)){
          for (k in (j + 1):dim){
            h = h + 1;
            temp_neardistM[j, k] = exp(- neardistM[(i - 1), h] / phi);
            temp_neardistM[k, j] = temp_neardistM[j, k];
          }
        }
        for(j in 1:dim){
          temp_neardistM[j, j] = 1;
        }
      }
      L = cholesky_decompose(temp_neardistM);
      u = to_vector(exp(- neardist[(i - 1), 1:dim] / phi));
      v = mdivide_left_tri_low(L, u);
      V[i] = 1 - (v' * v);
      v2 = mdivide_right_tri_low(v', L);
      Uw[i] = Uw[i] - v2 * w_r[nearind[(i - 1), 1:dim]];
    }
    V[1] = 1;
    out = sum(- 0.5 * log(V) - 0.5 / sigmasq * ((Uw .* Uw) ./ V)) - 0.5 * n * log(sigmasq);
    return out;
  }
  }

data {
  int n;
  int a;
  int M;
  int aL;
  vector[n] Y;
  int pn; // number of coefficient for point
  int pa; // number of coefficient for area
  matrix[n, pn] Xn; // design matrix for point
  matrix[a, pa] Xa; // design matrix for area
  int nearind[(n+aL) - 1, M];
  matrix[(n+aL) - 1, M] neardist;
  matrix[(n+aL) - 1, (M * (M - 1) / 2)] neardistM;
  matrix[n, n+aL] A0; // aggregation matrix for point
  matrix[a, n+aL] A1; // aggregation matrix for areal
  int Q[a]; // response at areas
  
}

parameters{
  vector[pn] beta; // coefficients
  vector[pa] alpha;
  real<lower = 0> sigma_sq;
  real<lower = 0> tau_sq;
  real<lower = 0> phi;
  vector[n+aL] w;
}

transformed parameters{
  vector[n] A0w;
  vector[a] A1w;  
  A0w = w[1:n];
  A1w = A1*w;
}

model{
  beta ~ normal(0, 5);
  alpha ~ normal(0, 5);
  sigma_sq ~ inv_gamma(2, 1);
  tau_sq ~ inv_gamma(2, 1);
  phi ~ normal(300,100);
  w ~ nngp_w(sigma_sq, phi, neardist, neardistM, nearind, n+aL, M);
  Y ~ normal(Xn * beta + A0w, sqrt(tau_sq));
  Q ~ poisson_log(Xa * alpha + A1w);
}

