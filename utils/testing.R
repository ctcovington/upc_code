# chisqtest <- function(x, p = rep(1, length(x)) / length(x)) {
#     n <- sum(x)
#     k <- length(x)
#     Tee <- sum((x - n * p)^2 / (n * p))
#     pvalue <- pchisq(Tee, df = k - 1, lower.tail = FALSE)
#     return(pvalue)
# }

pacman::p_load('dplyr', 'nortest', 'goftest')

unif_chisq_test <- function(x, n_bins) {
    binned_x <- cut(x, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE, labels = FALSE)
    counts <- table( factor(binned_x, levels = 1:n_bins) )
    p <- as.numeric( chisq.test(counts, p = rep(1/n_bins, n_bins))$p.value )
    return(p)
}

uniformity_test <- function(mat) {
    # take in matrix of u-values from "n_runs" runs (# runs will generally be # of posterior samples)
    # and generate p-value for test of uniformity for each run
    n_runs <- nrow(mat)
    p_values <- rep(0, n_runs)
    for (i in 1:n_runs) {
        # p_values[i] <- unif_chisq_test(mat[i,], 30)
        ad_test <- goftest::ad.test( qnorm(mat[i,], 0, 1), null = 'pnorm', mean = 0, sd = 1 )
        p_values[i] <- ad_test$p.value
        # resid_ks <- ks.test(mat[i,], punif)
        # p_values[i] <- resid_ks$p.value
    }
    return( pmax(0, pmin(1, p_values) ) )
}

# NOTE: this version combines all three tests
# resid_independence_test <- function(resid_mat) {
#     nc <- ncol(resid_mat)
#     indep_p_values <- rep(0, nrow(resid_mat))
#     for (i in 1:length(indep_p_values)) {
#         dcov_p <- dcov.test(resid_mat[i, 1:(nc-1)], resid_mat[i, 2:nc], R = 2000)$p
#         hoeffd_p <- hoeffd(resid_mat[i, 1:(nc-1)], resid_mat[i, 2:nc])$P[1,2]
#         dw_p <- dwtest(lm(resid_mat[i,] ~ 1))$p
#         indep_p_values[i] <- p_merger(c(dcov_p, hoeffd_p, dw_p), M = 10)
#     }
#     return(indep_p_values)
# }

build_hoeffD_null <- function(n, n_sims) {
    # this builds null for Hoeffding test of independence between two uniforms
    null <- rep(0, n_sims)
    for (i in 1:n_sims) {
        print(glue::glue('Building Hoeffding null: sim {i} of {n_sims}'))
        null[i] <- hoeffd(runif(n), runif(n))$D[1,2]
    }
    return(null)
}

build_lag_1_hoeffD_null <- function(n, n_sims) {
    # this builds null for Hoeffding test of independence between two uniforms
    null <- rep(0, n_sims)
    for (i in 1:n_sims) {
        y <- runif(n)
        null[i] <- hoeffd(y[1:(n-1)], y[2:n])$D[1,2]
    }
    return(null)
}

build_lag_2_hoeffD_null <- function(n, n_sims) {
    # this builds null for Hoeffding test of independence between two uniforms
    null <- rep(0, n_sims)
    for (i in 1:n_sims) {
        y <- runif(n)
        null[i] <- hoeffd(y[1:(n-2)], y[3:n])$D[1,2]
    }
    return(null)
}

build_hoeffD_index_null <- function(n, n_sims) {
    # this builds null for Hoeffding test of independence between uniform and indices
    null <- rep(0, n_sims)
    for (i in 1:n_sims) {
        print(glue::glue('Building Hoeffding index test null: sim {i} of {n_sims}'))
        null[i] <- hoeffd(1:n, runif(n) )$D[1,2]
    }
    return(null)
}

### TODO: this is version that would work for empirical null hoeffding
resid_independence_test <- function(resid_mat, null) {
    nc <- ncol(resid_mat)
    n <- length(null)
    indep_p_values <- rep(0, nrow(resid_mat))
    for (i in 1:length(indep_p_values)) {
        hoeff <- hoeffd(resid_mat[i, 1:(nc-1)], resid_mat[i, 2:nc])$D[1,2]
        pre_p <- mean(hoeff < null)
        indep_p_values[i] <- max( 1 / (n+1), min(pre_p, n / (n+1)))
    }
    return(indep_p_values)
}

resid_second_order_independence_test <- function(resid_mat, null) {
    nc <- ncol(resid_mat)
    n <- length(null)
    indep_p_values <- rep(0, nrow(resid_mat))
    for (i in 1:length(indep_p_values)) {
        hoeff <- hoeffd(resid_mat[i, 1:(nc-2)], resid_mat[i, 3:nc])$D[1,2]
        pre_p <- mean(hoeff < null)
        indep_p_values[i] <- max( 1 / (n+1), min(pre_p, n / (n+1)))
    }
    return(indep_p_values)
}

resid_index_independence_test <- function(resid_mat, null) {
    n <- ncol(resid_mat)
    ps <- rep(0, nrow(resid_mat))
    n_null <- length(null)
    for (i in 1:nrow(resid_mat)) {
        hoeff <- hoeffd(1:n, resid_mat[i, ])$D[1,2]
        pre_p <- mean(hoeff < null)
        ps[i] <- max( 1 / (n_null+1), min(pre_p, n_null / (n_null+1)))
    }
    return(ps)
}

# NOTE: this version uses only one test
# resid_independence_test <- function(resid_mat) {
#     nc <- ncol(resid_mat)
#     indep_p_values <- rep(0, nrow(resid_mat))
#     for (i in 1:length(indep_p_values)) {
#         # indep_p_values[i] <- hoeffd(resid_mat[i, 1:(nc-1)], resid_mat[i, 2:nc])$P[1,2]
#         dw_p <- dwtest(lm(resid_mat[i,] ~ 1))$p
#         indep_p_values[i] <- dw_p
#     }
#     return(indep_p_values)
# }

u_to_p_converter <- function(u_values) {
    # converts U values to p-values for test of uniformity
    return( pmax(0, pmin(1, 2*pmin(u_values, 1-u_values))) )
}

grid_harmonic_fn <- function(x, K) {
    inv_vec <- 1/seq(1, K)
    l <- sum(inv_vec)
    num <- K * (l*x <= 1)
    denom <- ceiling( K * l * x )
    return(num/denom)
}

# harmonic_mean_p_function <- function(p_vec) {
#     K <- length(p_vec)
#     return( K / (sum(1/p_vec)) )
# }

exchangeable_harmonic_mean_p_merger <- function(p_values) {
    # F_{EH} from middle of page 20 in https://arxiv.org/pdf/2404.03484
    K <- length(p_values)
    T_K <- log(K) + log(log(K)) + 1
    merged_p_value <- Inf
    for (l in 1:K) {
        # print(glue::glue('l={l} of {K}'))
        sorted_p_l <- sort(p_values[1:l])
        for (m in 1:l) {
            candidate_p_value <- (l*T_K/m + 1) * harmonic.mean(sorted_p_l[1:m])
            merged_p_value <- min(merged_p_value, candidate_p_value)
        }
    }
    return( min(1, merged_p_value) )
}

# exchangeable_harmonic_mean_p_merger <- function(p_values) {
#     # F_{EH} from middle of page 20 in https://arxiv.org/pdf/2404.03484
#     K <- length(p_values)
#     T_K <- log(K) + log(log(K)) + 1
#     merged_p_value <- Inf
    
#     harmonic_means <- numeric(K)
#     cumulative_sum <- 0
    
#     for (i in 1:K) {
#         cumulative_sum <- cumulative_sum + 1 / sort(p_values)[i]
#         harmonic_means[i] <- i / cumulative_sum
#     }
    
#     for (l in 1:K) {
#         for (m in 1:l) {
#             candidate_p_value <- (l * T_K / m + 1) * harmonic_means[m]
#             merged_p_value <- min(merged_p_value, candidate_p_value)
#         }
#     }
    
#     return(merged_p_value)
# }


p_merger <- function(p_values, M) {
    # approximate grid harmonic merging function up to error 2^{-M}
    # valid under arbitrary dependence
    # Algorithm 1 in Vovk, Wang, & Wang (2021)
    K <- length(p_values)
    L <- 0
    R <- 1
    for (m in 1:M) {
        epsilon <- (L+R)/2
        grid_harmonic_eval <- sapply( p_values/epsilon, grid_harmonic_fn, K )
        if (mean(grid_harmonic_eval) >= 1) {
            R <- epsilon
        } else {
            L <- epsilon
        }
    }
    return(R)
}

exchangeable_p_merger <- function(p_values, M) {
    # approximate exchangeable grid harmonic merging function up to error 2^{-M}
    # valid under exchangeability
    # Algorithm 1 in Appendix E of Gasparin, Wang, & Ramdas (2024)
    K <- length(p_values)
    L <- 0
    R <- 1
    for (m in 1:M) {
        epsilon <- (L+R)/2
        grid_harmonic_eval <- sapply( p_values/epsilon, grid_harmonic_fn, K )
        cumulative_mean <- cummean(grid_harmonic_eval)
        if (max(cumulative_mean) >= 1) {
            R <- epsilon
        } else {
            L <- epsilon
        }
    }
    return(R)
}

Cauchy_p_merger <- function(p_values) {
    return( 1-pcauchy(mean(tan(pi*(0.5-p_values)))) )
    # return(ACAT::ACAT(p_values))
}