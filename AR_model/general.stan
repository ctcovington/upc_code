data {
  int<lower=0> n;  // number of time points
  vector[n] y;     // data
  real phi_center; // center parameter for phi
  real phi_mult; // multiplicative transformation factor for phi
  real sigma_center; // center parameter for sigma
  real sigma_mult; // multiplicative transformation factor for sigma
}

parameters {
    array[2] real<lower=0,upper=1> u;
}

transformed parameters {
    real phi = (2*phi_center - phi_mult)/2 + phi_mult*u[1]; 
    real sigma = (2*sigma_center - sigma_mult)/2 + sigma_mult*u[2]; 
}

model {
    // set priors
    phi ~ uniform((2*phi_center - phi_mult)/2 + phi_mult*0, 
                  (2*phi_center - phi_mult)/2 + phi_mult*1);
    sigma ~ uniform((2*sigma_center - sigma_mult)/2 + sigma_mult*0, 
                  (2*sigma_center - sigma_mult)/2 + sigma_mult*1);
  
    // establish AR(1) inferential model
    for (i in 2:n) {
        y[i] ~ normal(phi * y[i-1], sigma);
    }
}
