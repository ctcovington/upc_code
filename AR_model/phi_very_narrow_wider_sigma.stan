data {
  int<lower=0> n;  // number of time points
  vector[n] y;     // data
  real s;
}

parameters {
    array[2] real<lower=0,upper=1> u;
}

transformed parameters {
  real phi = u[1]/10 - 0.05;
  real sigma = 3*u[2]+1;
}

model {
    // set priors
    phi ~ uniform(-0.05, 0.05);
    sigma ~ uniform(1, 4);
  
    // establish AR(1) inferential model
    for (i in 2:n) {
        y[i] ~ normal(phi * y[i-1], sigma);
    }
}
