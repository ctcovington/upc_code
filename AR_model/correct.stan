data {
  int<lower=0> n;  // number of time points
  vector[n] y;     // data
  real s;
}

parameters {
  real <lower=-0.5,upper=0.5> phi;
  real <lower=1,upper=2> sigma;
}

model {
    // set priors
    phi ~ normal(0, 0.4) T[-0.5, 0.5]; # truncated normal
    sigma ~ normal(1.5, 0.4) T[1, 2]; # truncated normal
  
    // establish AR(1) inferential model
    for (i in 2:n) {
        y[i] ~ normal(phi * y[i-1], sigma);
    }
}
