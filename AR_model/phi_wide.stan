data {
  int<lower=0> n;  // number of time points
  vector[n] y;     // data
}

parameters {
    array[2] real<lower=0,upper=1> u;
}

transformed parameters {
  real phi = u[1]*4/5 - 0.4;
  real sigma = u[2]+1;
}

model {
    // set priors
    phi ~ uniform(-0.4, 0.4);
    sigma ~ uniform(1, 2);
  
    // establish AR(1) inferential model
    for (i in 2:n) {
        y[i] ~ normal(phi * y[i-1], sigma);
    }
}
