data {
  int<lower=0> n;  // number of time points
  vector[n] y;     // data
}

parameters {
    array[2] real<lower=0,upper=1> u;
}

transformed parameters {
  real phi = u[1]-0.5;  
  real sigma = 3*u[2];  
}

model {
    // set priors
    phi ~ uniform(-0.5, 0.5);
    sigma ~ uniform(0, 3);
  
    // establish AR(1) inferential model
    for (i in 2:n) {
        y[i] ~ normal(phi * y[i-1], sigma);
    }
}
