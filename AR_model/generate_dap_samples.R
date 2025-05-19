# generate AR(1) data
gen_AR1 <- function(n, phi_1, sigma) {
    Y <- rep(0, n)
    Y[1] <- rnorm(1, mean = 0, sd = sigma)
    for (i in 2:n) {
        Y[i] <- rnorm(1, mean = phi_1 * Y[i-1], sd = sigma)
    }
    return(Y)
}

# generate AR(1) data with heteroskedasticity
gen_AR1_heteroskedastic <- function(n, phi_1, sigma) {
    Y <- rep(0, n)
    Y[1] <- rnorm(1, mean = 0, sd = sigma)
    for (i in 2:n) {
        Y[i] <- rnorm(1, mean = phi_1 * Y[i-1], sd = ( (n+i)/(2*n) )*sigma)
    }
    return(Y)
}

gen_AR1_heteroskedastic_v1 <- function(n, phi_1, sigma, s) {
    Y <- rep(0, n)
    Y[1] <- rnorm(1, mean = 0, sd = sigma)
    for (i in 2:n) {
        # c_i <- 1 + s*(2*i-n-1)/(2*n)
        c_i <- 1 + s*(2*i-n-1)/(n)
        Y[i] <- rnorm(1, mean = phi_1 * Y[i-1], sd = sqrt(c_i)*sigma)
    }
    return(Y)
}

gen_AR1_heteroskedastic_v2 <- function(n, phi_1, sigma, s) {
    Y <- rep(0, n)
    Y[1] <- rnorm(1, mean = 0, sd = sigma)
    for (i in 2:n) {
        # c_i <- 1 + s*(2*i-n-1)/(2*n)
        c_i <- 1 + s*(2*i-n-1)/(n)
        Y[i] <- rnorm(1, mean = phi_1 * Y[i-1], sd = c_i*sigma)
    }
    return(Y)
}

# generate AR(2) data
gen_AR2 <- function(n, phi_1, sigma, s) {
    Y <- rep(0, n)
    Y[1] <- rnorm(1, mean = 0, sd = sigma)
    Y[2] <- rnorm(1, mean = phi_1*Y[1], sd = sigma)
    for (i in 3:n) {
        Y[i] <- rnorm(1, mean = phi_1 * Y[i-1] + s*sign(phi_1)*(abs(phi_1))**(1/4) * Y[i-2], sd = sigma)
    }
    return(Y)
}

gen_Cauchy <- function(n, phi_1, sigma, s) {
    Y <- rep(0, n)
    Y[1] <- rcauchy(1, location = 0, scale = sigma)
    Y[2] <- rcauchy(1, location = phi_1*Y[1], scale = sigma)
    for (i in 3:n) {
        Y[i] <- rcauchy(1, location = phi_1 * Y[i-1] + s*sign(phi_1)*(abs(phi_1))**(1/4) * Y[i-2], 
                           scale = sigma)
    }
    return(Y)
}

generate_dap_samples <- function(model.stan, n, DGP, n_samples, s) {

    phi_1 <- rtruncnorm(n = 1, a = -0.5, b = 0.5, mean = 0, sd = 0.4)
    sigma <- rtruncnorm(n = 1, a = 1, b = 2, mean = 1.5, sd = 0.4)


    # generate data
    if (DGP == 'AR1') {
        y <- gen_AR1(n, phi_1, sigma)
    } else if (DGP == 'AR1_heteroskedastic_V1') {
        y <- gen_AR1_heteroskedastic_v1(n, phi_1, sigma, s)
    } else if (DGP == 'AR1_heteroskedastic_V2') {
        y <- gen_AR1_heteroskedastic_v2(n, phi_1, sigma, s)
    } else if (DGP == 'AR2') {
        y <- gen_AR2(n, phi_1, sigma, s)
    } else if (DGP == 'Cauchy') {
        y <- gen_Cauchy(n, phi_1, sigma, s)
    }

    stan_data <- list(n = n, y = y, s = s)

    # fit model in stan
    n_burnin <- 1000
    n_chains <- 1

    # fit model and draw from posterior to estimate DAP
    model.fit <- sampling(model.stan, data = stan_data, iter = n_burnin+n_samples, 
                            warmup = n_burnin, chains = n_chains, 
                            verbose = FALSE, refresh = 0)

    dap_samples <- data.table(as.data.frame(model.fit))
    dap_samples[, lp__ := NULL]
    colnames(dap_samples) <- c('phi', 'sigma')

    return(list('dap_samples' = dap_samples, 'y' = y))
}