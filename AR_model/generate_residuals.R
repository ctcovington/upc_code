# resid_u_generation <- function(res_dt, y_mat) {
#     n_sims <- nrow(y_mat)
#     n <- ncol(y_mat)
#     u_resid_mat <- matrix(0, n_sims, n)
#     for (i in 1:n_sims) {
#         y <- as.numeric(y_mat[i,])
#         phi <- as.numeric(res_dt[i, 'phi'])
#         sigma <- as.numeric(res_dt[i, 'sigma'])

#         # calculate residuals and map through implied normal cdf
#         resids <- c( y[1], y[2:n] - phi*y[1:(n-1)] )
#         u_resids <- pnorm(resids, mean = 0, sd = sigma)
#         u_resid_mat[i, ] <- u_resids
#     }
#     return(u_resid_mat)
# }

resid_u_generation <- function(res_dt, y) {
    n <- length(y)
    n_samples <- nrow(res_dt)
    u_resid_mat <- matrix(0, n_samples, n)
    for (i in 1:n_samples) {
        phi <- as.numeric(res_dt[i, 'phi'])
        sigma <- as.numeric(res_dt[i, 'sigma'])

        # calculate residuals and map through implied normal cdf
        resids <- c( y[1], y[2:n] - phi*y[1:(n-1)] )
        u_resid_mat[i,] <- pnorm(resids, mean = 0, sd = sigma)
    }

    return(u_resid_mat)
}