sim_normal <- function(n, mu, sigma) {
    return( rnorm(n, mu, sigma) )
}

sim_bernoulli <- function(n, p, p_repeat) {
    # generate dependent bernoulli
    y <- rep(0, n)
    y[1] <- sample(c(0,1), prob = c(1-p,p), size = 1)
    for (i in 2:n) {
        u <- runif(1, 0, 1)
        # with probability "p_repeat", simply repeat last observation
        # otherwise, draw from Bernoulli(p)
        if (u <= p_repeat) {
            y[i] <- y[i-1]
        } else {
            y[i] <- sample(c(0,1), prob = c(1-p,p), size = 1)
        }
    }
    return(y)
}

create_smoking_outcome <- function(smoking, betas, alphas) {
    p_i <- inv.logit( betas[1]*smoking$wave + betas[2]*smoking$parsmk + betas[3]*smoking$wave*smoking$parsmk + alphas[smoking$newid] )
    smkreg <- rbinom(n = nrow(smoking), size = 1, prob = p_i)
    return(smkreg)
}

create_random_smkreg <- function(smoking) {
    n_new_id <- length(unique(smoking$newid))

    # generate latent variables 
    p <- rbeta(1, 1, 1) # probability that each individual is in group Z=1
    Z <- rbinom(n_new_id, 1, p) # assign each newid a latent class

    # draw from priors on random effect parameter distributions
    mu <- rnorm(2, 0, 5)
    tau <- rgamma(2, shape = 1, rate = 1/5)

    # draw coefficients from priors
    betas <- rnorm(3, 0, 5) # coefficients on wave/parsmk
    alphas <- rep(0, n_new_id) # random effect coefficients
    alphas[Z == 0] <- rnorm( sum(Z==0), mu[1], 1/tau[1] )
    alphas[Z == 1] <- rnorm( sum(Z==1), mu[2], 1/tau[2] )

    # create simulated smoking data
    parsmk <- create_smoking_outcome(smoking, betas, alphas)   

    return(parsmk)
}