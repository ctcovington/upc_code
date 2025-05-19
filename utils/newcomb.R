#       i think it is based on https://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall07/reading/NG.pdf 
generate_normalgamma_params <- function(n,mu,kappa,alpha,beta) {
    lambda_samples <- rep(0, n)
    mu_samples <- rep(0, n)
    for (i in 1:n) {
        lambda_samples[i] <- rgamma(1, shape = alpha, rate = beta)
        mu_samples[i] <- rnorm(1, mean = mu, sd = (kappa*lambda_samples[i])^(-1/2))
    }
    param_matrix <- matrix( c(mu_samples, lambda_samples), ncol = 2 )
    colnames(param_matrix) <- c('mu', 'lambda')
    return( data.table(param_matrix) )
}