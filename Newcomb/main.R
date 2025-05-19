library(pacman)
pacman::p_load('rstan', 'data.table', 'ggplot2', 'forcats', 'energy', 
               'gridExtra', 'EnvStats', 'stringr', 'Hmisc', 'lmtest',
               'latex2exp', 'glue', 'scales')
               
source('../utils/generate_plots.R')
source('../utils/testing.R')
source('../utils/sim_data.R')
source('../utils/newcomb.R')

# parent directories
data_dir <- file.path('..', '..', 'data', 'Newcomb')
results_dir <- file.path('..', '..', 'results', 'Newcomb')

dgps <- c( 
            'sim', 
            'BDA'
        )

priors <- c(
            'empirical_prior', 
            'weakly_informative_prior', 
            'very_bad_prior'
            )

# plot histograms of actual data
real_data <- fread( file.path(data_dir, 'newcomb.txt') )
real_data <- data.table( 1000*(real_data$V1 - 24.8) ) # transform data as done in BDA

sim_data <- sim_normal(n = nrow(real_data), mu = mean(real_data$V1), sigma = sd(real_data$V1))
sim_data <- data.table(sim_data)
colnames(sim_data) <- 'V1'

real_data_hist <- ggplot(real_data, aes(x = V1)) + 
                    geom_histogram(position = 'identity', colour = 'black', fill = 'blue') +
                    labs(x = '(relative) speed of light (in nanoseconds)') + 
                    theme(
                        axis.title.x = element_text(size = 36),
                        axis.text.x = element_text(size = 36),
                        # axis.ticks.x = element_blank(),
                        axis.title.y = element_text(size = 36),
                        axis.text.y = element_text(size = 36),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )
ggsave(file.path(results_dir, 'real_data_histogram.png'), real_data_hist, width = 12, height = 8, dpi = 300)

sim_data_hist <- ggplot(sim_data, aes(x = V1)) + 
                    geom_histogram(position = 'identity', colour = 'black', fill = 'blue') +
                    labs(x = '(relative) speed of light (in nanoseconds)') + 
                    theme(
                        axis.title.x = element_text(size = 36),
                        axis.text.x = element_text(size = 36),
                        # axis.ticks.x = element_blank(),
                        axis.title.y = element_text(size = 36),
                        axis.text.y = element_text(size = 36),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )
ggsave(file.path(results_dir, 'sim_data_histogram.png'), sim_data_hist, width = 12, height = 8, dpi = 300)

n_posterior_samples <- 500000

for (dgp in dgps) {
    for (prior in priors) {
        print( glue::glue('dgp: {dgp}, prior: {prior}') )
        # # load model from stan file
        # model.stan <- stan_model(file = paste0(prior, '.stan'))

        # # set up directory structure
        prior_results_dir <- file.path(results_dir, dgp, prior)
        dir.create(prior_results_dir, showWarnings = FALSE, recursive = TRUE)

        # total_mu_p_val_path <- file.path(prior_results_dir, 'total_mu_p_val_dt.csv')
        # total_sigma_p_val_path <- file.path(prior_results_dir, 'total_sigma_p_val_dt.csv')
        # total_data_p_val_path <- file.path(prior_results_dir, 'total_data_p_val_dt.csv')
        # single_sample_p_val_path <- file.path(prior_results_dir, 'single_sample_p_val_dt.csv')
        u_resid_mat_path <- file.path(prior_results_dir, 'u_resid_mat.csv')
        # y_path <- file.path(prior_results_dir, 'y_mat.csv')
        # p_mu_density_mat_path <- file.path(prior_results_dir, 'p_mu_density_mat.csv')
        # p_sigma_density_mat_path <- file.path(prior_results_dir, 'p_sigma_density_mat.csv')

        # read Newcomb data
        newcomb <- fread( file.path(data_dir, 'newcomb.txt') )
        newcomb <- data.table( 1000*(newcomb$V1 - 24.8) ) # transform data as done in BDA
        if (dgp == 'sim') {
            newcomb <- sim_normal(n = nrow(newcomb), mu = mean(newcomb$V1), sigma = sd(newcomb$V1))
            newcomb <- data.table(newcomb)
            colnames(newcomb) <- 'V1'
        }
        y <- newcomb$V1 # save data simply as y
        n <- length(y)

        # set prior parameters for NormalGamma model
        # TODO: change these depending on the prior
        mu_0 <- case_when(
                                prior == 'empirical_prior' ~ mean(y),
                                prior == 'very_bad_prior' ~ 179,
                                prior == 'weakly_informative_prior' ~ 0,
                                TRUE ~ NA_real_
                                )
        kappa_0 <- case_when(
                                prior == 'empirical_prior' ~ n,
                                prior == 'very_bad_prior' ~ n,
                                prior == 'weakly_informative_prior' ~ 0.1,
                                TRUE ~ NA_real_
                                )

        alpha_0 <- case_when(
                                prior == 'empirical_prior' ~ n/2,
                                prior == 'very_bad_prior' ~ n/2,
                                prior == 'weakly_informative_prior' ~ 2,
                                TRUE ~ NA_real_
                                )

        beta_0 <- case_when(
                                prior == 'empirical_prior' ~ sd(y)^2*alpha_0,
                                prior == 'very_bad_prior' ~ 42^2*kappa_0*(alpha_0-1),
                                prior == 'weakly_informative_prior' ~ 300,
                                TRUE ~ NA_real_
                                )
        
        
        # mu_0 <- 0

        # if (X,T) ~ NormalGamma, s.t. T = precision, we've been saying that:
        # X | T ~ N(\mu, 1/sqrt(T))
        # weakly informative: \mu ~ N(0, 100^2), 1/sqrt(T) ~ Unif(0,100)
        #                     E X = 0,
        #                     E T \approx 8.62
        # data-dependent: \mu ~ N( \bar{X}, \hat{\sigma}^2 ), 1/sqrt(T) ~ Unif(0, 2 \hat{\sigma}))
        #                     E X = \bar{X}
        #                     E T \approx 205
        # poorly chosen: \mu ~ N( 179, 42^2 ), 1/sqrt(T) ~ Unif(0, 100)
        #                     E X = 179
        #                     E T \approx 8.62



        # let's set kappa_0 = 1, alpha_0 = 1 for weakly informative and poorly chosen, 
        # set kappa_0 = n, alpha_0 = n/2 for data-dependent
        # then set beta_0 = E[sd]^2 * alpha_0 where E[sd] = 50 for weakly informative and poorly chosen,
        #                                           E[sd] = sd(y) for data-dependent
        # we know E X = mu_0
        #         E T = beta_0^(-1)
        #
        #        Var X = 
        #        Var T = beta_0^(-2)



        # kappa_0 <- 1
        # alpha_0 <- 1
        # beta_0 <- 1



        #### test marginals for prior ####
        prior_params <- generate_normalgamma_params(n_posterior_samples, mu_0, kappa_0, alpha_0, beta_0)
        prior_params$sigma <- sqrt(1/prior_params$lambda)
        prior_mu_hist <- ggplot(prior_params, aes(x = pmax(-100, pmin(100, mu)))) +
                                geom_histogram(bins = 100, colour = 'black', fill = 'blue')
        prior_sigma_hist <- ggplot(prior_params, aes(x = pmin(100, sigma))) + 
                                geom_histogram(bins = 100, colour = 'black', fill = 'blue')
        ##################################
        
        observed_min <- min(y)

        # calculate analytic posterior
        mu_n <- (kappa_0*mu_0 + sum(y)) / (kappa_0 + n)
        kappa_n <- kappa_0 + n
        alpha_n <- alpha_0 + n/2
        beta_n <- beta_0 + (1/2) * (kappa_0*mu_0^2 - kappa_n*mu_n^2 + sum(y^2))

        # sample from posterior (params and predictive) and generate residual u-values
        posterior_samples <- generate_normalgamma_params(n_posterior_samples, mu_n, kappa_n, alpha_n, beta_n)
        posterior_samples[, sigma := lambda^(-1/2)]
        u_resid_mat <- matrix(0, n_posterior_samples, n)
        y_pred_mat <- matrix(0, n_posterior_samples, n)
        for (j in 1:n_posterior_samples) {
            mu_sample <- as.numeric(posterior_samples[j, 'mu'])
            sigma_sample <- as.numeric(posterior_samples[j, 'sigma'])

            y_pred_mat[j, ] <- rnorm(n, mu_sample, sigma_sample)
            

            u_resid_mat[j,] <- pnorm( y, mean = mu_sample, sd = sigma_sample )            
        }
        # if (i == 1) {
        fwrite(u_resid_mat, u_resid_mat_path)
        # }

        # convert posterior samples to u-values
        # TODO: use prior marginals from https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf
        # prior_lambda_samples <- rgamma(nrow(u_samples), shape = alpha_0, rate = beta_0)^(-1/2)
        u_samples <- copy(posterior_samples)
        u_samples[, lambda := pgamma(lambda, shape = alpha_0, rate = beta_0)]
        u_samples[, mu := pnorm(mu, mu_0, kappa_0^(-1/2) * posterior_samples[, sigma]) ] 
        # u_samples[, mu := pnorm(mu, mu_n, kappa_0 * prior_lambda_samples) ]
        u_samples[, sigma := lambda]

        # generate p-values for each posterior sample
        resid_p_values <- uniformity_test(u_resid_mat) # p-values for resid uniformity testing
        mu_p_values <- u_to_p_converter(u_samples$mu) # convert parameter u values to p values
        sigma_p_values <- u_to_p_converter(u_samples$sigma) # convert parameter u values to p values

        # going back for now to empirical density estimates, bc I haven't worked out analytic for p-values yet
        min_posterior_predictive_p <- mean( apply(y_pred_mat, 1, function(row) {
            min(row) <= observed_min
        }) )


    ####################################################################################################
    ## below is old stuff 
    #####################################################################################################


        # get minimum predicted value from each posterior
        y_pred_mins <- apply(y_pred_mat, 1, function(row) min(row) )
        y_pred_min_dt <- data.table(y_pred_mins)

        # get actual observed minimum value
        # observed_min <- min(newcomb$V1)

        # # generate residual u values
        # u_resid_mat <- matrix(0, nrow(parameter_us), stan_data$n)
        # for (i in 1:nrow(u_resid_mat)) {
        #     u_resid_mat[i,] <- pnorm( newcomb$V1 - mu_samples[i], mean = 0, sd = sigma_samples[i] )
        # }

        # plot histogram of min values 
        min_value_histogram(y_pred_min_dt, observed_min, prior_results_dir)
        min_value_density(y_pred_min_dt, observed_min, prior_results_dir)

        # generate residual trajectory plot
        # NOTE: remove because we don't use it and it scales poorly 
        # residual_trajectory_plot(u_resid_mat, prior_results_dir)

        # generate histograms of parameter u values
        parameter_us <- data.table('u1' = u_samples$mu, 'u2' = u_samples$sigma)
        parameter_names <- c('mu', 'sigma')
        LaTeX_parameter_names <- c('$\\mu$', '$\\sigma$')
        parameter_hist_plots(parameter_us, parameter_names, LaTeX_parameter_names, prior_results_dir)
        parameter_density_plots(parameter_us, parameter_names, LaTeX_parameter_names, prior_results_dir)

        ### generate uniformity test results ###
        # # generate p-values for each posterior sample
        # resid_p_values <- uniformity_test(u_resid_mat) # p-values for resid uniformity testing
        # mu_p_values <- u_to_p_converter(parameter_us$u1) # convert parameter u values to p values
        # sigma_p_values <- u_to_p_converter(parameter_us$u2) # convert parameter u values to p values

        # combine p-values using Cauchy method
        resid_p_value <- Cauchy_p_merger(resid_p_values)
        mu_p_value <- Cauchy_p_merger(mu_p_values)
        sigma_p_value <- Cauchy_p_merger(sigma_p_values)

        # generate plots of data u values
        u_resid_mat_sub <- u_resid_mat[sample(1:nrow(u_resid_mat), size = 1000, replace = FALSE), ]
        p_val_string <- signif(resid_p_value, 2)
        one_posterior_sample_data_u_density(u_resid_mat, prior_results_dir)
        one_posterior_sample_data_u_histogram(u_resid_mat, prior_results_dir)
        one_posterior_sample_residual_ecdf_plot(u_resid_mat, prior_results_dir)
        one_posterior_sample_residual_tilted_ecdf_plot(u_resid_mat, prior_results_dir)
        residual_ecdf_plot(u_resid_mat_sub, 50, p_val_string, prior_results_dir)
        residual_tilted_ecdf_plot(u_resid_mat_sub, 50, p_val_string, prior_results_dir, 0.35)

        # plot posterior p-values
        p_val_dt <- data.table('mu' = mu_p_values, 'sigma' = sigma_p_values, 'data' = resid_p_values)
        Newcomb_p_value_histograms(p_val_dt, prior_results_dir, plot_min_posterior_predictive = FALSE)
        Newcomb_p_value_densities(p_val_dt, prior_results_dir, plot_min_posterior_predictive = FALSE)
        # TODO: do approximated versions

        # # approximate p_value densities
        # mu_0 <- case_when(
        #                         prior == 'empirical_prior' ~ mean(newcomb$V1),
        #                         prior == 'very_bad_prior' ~ 179,
        #                         prior %in% c('weakly_informative_prior') ~ 0,
        #                         TRUE ~ NA_real_
        #                         )
        # sigma_0 <- case_when(
        #                         prior == 'empirical_prior' ~ sd(newcomb$V1),
        #                         prior == 'very_bad_prior' ~ 42,
        #                         prior %in% c('weakly_informative_prior') ~ 100,
        #                         TRUE ~ NA_real_
        #                         )
        # max_sigma <- case_when(
        #                         prior == 'empirical_prior' ~ 2*sd(newcomb$V1),
        #                         prior %in% c('weakly_informative_prior', 'very_bad_prior') ~ 100,
        #                         TRUE ~ NA_real_
        #                       )
        # y <- newcomb$V1 

        # # TODO: work in progress -- p_sigma doesn't necessarily seem correct (gives 0 for many priors -- 
        # #       might be mathematically correct but experiencing underflow)
        # # density_dt_list <- approximated_Newcomb_p_value_densities(y, mu_samples, sigma_samples, mu_0, sigma_0, max_sigma, prior_results_dir, plot_boolean = TRUE)
        # # p_mu_density_dt <- density_dt_list[[1]]
        # # p_sigma_density_dt <- density_dt_list[[2]]

        test_file <- file.path(prior_results_dir, 'test_output.txt')
        output <- glue::glue('residual uniformity test: p = {resid_p_value}\nmu uniformity test: p = {mu_p_value}\nsigma uniformity test: p = {sigma_p_value}')
        writeLines(output, test_file)
    }
}

################################################################################################################################

### get distribution of posterior predictive p-value for min (as is done in BDA) under the correct model ###
set.seed(1)

dgp <- 'sim'
prior <- 'correct_prior'

# set up directory structure
prior_results_dir <- file.path(results_dir, dgp, prior)
dir.create(prior_results_dir, showWarnings = FALSE, recursive = TRUE)

# load model from stan file
# model.stan <- stan_model(file = paste0(prior, '.stan'))

# load data 
real_data <- fread( file.path(data_dir, 'newcomb.txt') )
real_data <- data.table( 1000*(real_data$V1 - 24.8) ) # transform data as done in BDA

# establish correct parameters for mu/sigma
n <- nrow(real_data)

n_sims <- 10^5
# n_sims <- 20
resid_ps <- rep(0, n_sims)
mu_ps <- rep(0, n_sims)
sigma_ps <- rep(0, n_sims)
single_sample_resid_ps <- rep(0, n_sims)
single_sample_mu_ps <- rep(0, n_sims)
single_sample_sigma_ps <- rep(0, n_sims)
min_posterior_predictive_ps <- rep(0, n_sims)
# y_diff_posterior_predictive_ps <- rep(0, n_sims)

p_val_dir <- file.path(prior_results_dir, glue::glue('n_sims={n_sims}'))
dir.create(p_val_dir, showWarnings = FALSE, recursive = TRUE)
p_val_path <- file.path(p_val_dir, 'p_val_dt.csv')
total_mu_p_val_path <- file.path(p_val_dir, 'total_mu_p_val_dt.csv')
total_sigma_p_val_path <- file.path(p_val_dir, 'total_sigma_p_val_dt.csv')
total_data_p_val_path <- file.path(p_val_dir, 'total_data_p_val_dt.csv')
single_sample_p_val_path <- file.path(p_val_dir, 'single_sample_p_val_dt.csv')
u_resid_mat_path <- file.path(p_val_dir, 'u_resid_mat.csv')
y_path <- file.path(p_val_dir, 'y_mat.csv')
p_mu_density_mat_path <- file.path(p_val_dir, 'p_mu_density_mat.csv')
p_sigma_density_mat_path <- file.path(p_val_dir, 'p_sigma_density_mat.csv')
p_data_density_mat_path <- file.path(p_val_dir, 'p_data_density_mat.csv')

# set desired number of posterior samples
lower_n_posterior_samples <- 50
higher_n_posterior_samples <- 100000

# store posterior p values 
total_resid_ps <- matrix(0, n_sims, lower_n_posterior_samples)
total_mu_ps <- matrix(0, n_sims, lower_n_posterior_samples)
total_sigma_ps <- matrix(0, n_sims, lower_n_posterior_samples)
y_mat <- matrix(0, n_sims, n)

n_posterior_samples_plotting <- 20
samples <- sample(1:nrow(total_mu_ps), size = n_posterior_samples_plotting, replace = FALSE)
colors <- hue_pal()(n_posterior_samples_plotting) 

# initialize objects for average density calculations

# repeating 0 bc I couldn't figure out another good way for the empirical density to enforce 0 at the ends and capture densities in the tails
# NOTE: this is necessary only for empirical density calculations -- analytic do not need the extra 0
# x_vals <- c(0, seq(0, 1, by = 0.01)) 
x_vals <- seq(0, 1, by = 0.01)

p_mu_density_mat <- matrix(0, n_sims, length(x_vals))
p_sigma_density_mat <- matrix(0, n_sims, length(x_vals))
p_data_density_mat <- matrix(0, n_sims, length(x_vals)+1)

if (file.exists(p_val_path)) {
    p_val_dt <- fread( p_val_path )
    single_sample_p_val_dt <- fread( single_sample_p_val_path )
    total_mu_ps <- fread( total_mu_p_val_path )
    total_sigma_ps <- fread( total_sigma_p_val_path )
    total_resid_ps <- fread( total_data_p_val_path )
    u_resid_mat <- fread( u_resid_mat_path )
    y_mat <- fread( y_path )
    p_mu_density_mat <- read.table( p_mu_density_mat_path )
    p_sigma_density_mat <- read.table( p_sigma_density_mat_path )
    p_data_density_mat <- read.table( p_data_density_mat_path )
} else {
    for (i in 1:n_sims) {
        print( glue::glue('sim {i} of {n_sims}') )

        ####################################################################
        ### TODO: trying here to derive analytic posteriors for mu/sigma ###
        ####################################################################
        
        # set prior parameters for NormalGamma model
        mu_0 <- 0
        kappa_0 <- 1
        alpha_0 <- 1
        beta_0 <- 1

        # sample data using prior parameters
        dgp_params <- generate_normalgamma_params(1, mu_0, kappa_0, alpha_0, beta_0)
        y <- rnorm( n, dgp_params[[1]], dgp_params[[2]]^(-1/2) )

        observed_min <- min(y)

        # calculate analytic posterior
        mu_n <- (kappa_0*mu_0 + sum(y)) / (kappa_0 + n)
        kappa_n <- kappa_0 + n
        alpha_n <- alpha_0 + n/2
        beta_n <- beta_0 + (1/2) * (kappa_0*mu_0^2 - kappa_n*mu_n^2 + sum(y^2))

        # sample from posterior (params and predictive) and generate residual u-values
        if (i %in% samples) {
            n_posterior_samples <- higher_n_posterior_samples # bump up posterior simulations for sim where I generate plots
        } else {
            n_posterior_samples <- lower_n_posterior_samples
        }
        posterior_samples <- generate_normalgamma_params(n_posterior_samples, mu_n, kappa_n, alpha_n, beta_n)
        posterior_samples[, sigma := lambda^(-1/2)]
        u_resid_mat <- matrix(0, n_posterior_samples, n)
        y_pred_mat <- matrix(0, n_posterior_samples, n)
        for (j in 1:n_posterior_samples) {
            mu_sample <- as.numeric(posterior_samples[j, 'mu'])
            sigma_sample <- as.numeric(posterior_samples[j, 'sigma'])

            y_pred_mat[j, ] <- rnorm(n, mu_sample, sigma_sample)
            

            u_resid_mat[j,] <- pnorm( y, mean = mu_sample, sd = sigma_sample )            
        }

        if (i == 10) {
            fwrite(u_resid_mat, u_resid_mat_path)
        }

        # convert posterior samples to u-values
        # TODO: use prior marginals from https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf
        # prior_lambda_samples <- rgamma(nrow(u_samples), shape = alpha_0, rate = beta_0)^(-1/2)
        u_samples <- copy(posterior_samples)
        u_samples[, lambda := pgamma(lambda, shape = alpha_0, rate = beta_0)]
        u_samples[, mu := pnorm(mu, mu_0, kappa_0^(-1/2) * posterior_samples[, sigma]) ] # TODO: this version uses posterior samples from lambda
        u_samples[, sigma := lambda]

        # generate p-values for each posterior sample
        resid_p_values <- uniformity_test(u_resid_mat) # p-values for resid uniformity testing
        mu_p_values <- u_to_p_converter(u_samples$mu) # convert parameter u values to p values
        sigma_p_values <- u_to_p_converter(u_samples$sigma) # convert parameter u values to p values

        if (i %in% samples) {
            y_dt <- data.table('y' = y)
            y_ecdf <- ggplot(y_dt, aes(x = y)) + 
                        stat_ecdf(position = 'identity', linewidth = 1, alpha = 1, colour = 'blue', show.legend = FALSE)

            u_resid_dt <- data.table( t(u_resid_mat) )
            u_resid_dt_long <- melt(u_resid_dt)
            u_ecdf_samples <- sample(u_resid_dt_long$variable, size = 1000, replace = FALSE)
            u_ecdf <- ggplot(u_resid_dt_long[variable %in% u_ecdf_samples], aes(x = value, group = variable)) + 
                        stat_ecdf(position = 'identity', linewidth = 1, alpha = 0.01, colour = 'blue', show.legend = FALSE) +
                        geom_abline(intercept = 0, slope = 1, colour = 'black') + 
                        coord_cartesian(xlim = c(0,1), ylim = c(0,1))
            
            u_hist <- ggplot(u_resid_dt_long, aes(x = value)) + 
                        geom_histogram(breaks = seq(0, 1, 0.005, colour = 'black', fill = 'blue'))

            unif_p_dt <- data.table('p_unif' = resid_p_values)
            unif_p_hist <- ggplot(unif_p_dt, aes(x = p_unif)) + 
                                geom_histogram(breaks = seq(0, 1, 0.01), colour = 'black', fill = 'blue', show.legend = FALSE)
            
            centered_posterior_mu <- posterior_samples$mu - dgp_params$mu
            centered_posterior_sigma <- posterior_samples$sigma - dgp_params$lambda^(-1/2)

            centered_dt <- data.table('mu_posterior' = posterior_samples$mu,
                                      'sigma_posterior' = posterior_samples$sigma,
                                      'mu_posterior_error' = centered_posterior_mu,
                                      'sigma_posterior_error' = centered_posterior_sigma
                                     )
            
            mu_hist <- ggplot(centered_dt, aes(x = mu_posterior)) + 
                                geom_histogram(bins = 100, colour = 'black', fill = 'blue')
            sigma_hist <- ggplot(centered_dt, aes(x = sigma_posterior)) + 
                                geom_histogram(bins = 100, colour = 'black', fill = 'blue')
            mu_error_hist <- ggplot(centered_dt, aes(x = mu_posterior_error)) + 
                                geom_histogram(bins = 100, colour = 'black', fill = 'blue')
            sigma_error_hist <- ggplot(centered_dt, aes(x = sigma_posterior_error)) + 
                                geom_histogram(bins = 100, colour = 'black', fill = 'blue')

            ggsave(file.path(p_val_dir, glue::glue('y_ecdf_{i}.png')), y_ecdf, width = 12, height = 8, dpi = 300)
            ggsave(file.path(p_val_dir, glue::glue('u_ecdf_{i}.png')), u_ecdf, width = 12, height = 8, dpi = 300)
            ggsave(file.path(p_val_dir, glue::glue('u_hist_{i}.png')), u_hist, width = 12, height = 8, dpi = 300)
            ggsave(file.path(p_val_dir, glue::glue('data_unif_p_hist_{i}.png')), unif_p_hist, width = 12, height = 8, dpi = 300)
            ggsave(file.path(p_val_dir, glue::glue('mu_hist_{i}.png')), mu_hist, width = 12, height = 8, dpi = 300)
            ggsave(file.path(p_val_dir, glue::glue('sigma_hist_{i}.png')), sigma_hist, width = 12, height = 8, dpi = 300)
            ggsave(file.path(p_val_dir, glue::glue('mu_error_hist_{i}.png')), mu_error_hist, width = 12, height = 8, dpi = 300)
            ggsave(file.path(p_val_dir, glue::glue('sigma_error_hist_{i}.png')), sigma_error_hist, width = 12, height = 8, dpi = 300)
        }

        # going back for now to empirical density estimates, bc I haven't worked out analytic for p-values yet
        min_posterior_predictive_ps[i] <- mean( apply(y_pred_mat, 1, function(row) {
            min(row) <= observed_min
        }) )

        # store p-values in aggregate table
        total_resid_ps[i, 1:lower_n_posterior_samples] <- resid_p_values[1:lower_n_posterior_samples]
        total_mu_ps[i, 1:lower_n_posterior_samples] <- mu_p_values[1:lower_n_posterior_samples]
        total_sigma_ps[i, 1:lower_n_posterior_samples] <- sigma_p_values[1:lower_n_posterior_samples]

        # combine p-values using Cauchy method
        single_sample_resid_ps[i] <- resid_p_values[1]
        single_sample_mu_ps[i] <- mu_p_values[1]
        single_sample_sigma_ps[i] <- sigma_p_values[1]
        resid_ps[i] <- Cauchy_p_merger(resid_p_values)
        mu_ps[i] <- Cauchy_p_merger(mu_p_values)
        sigma_ps[i] <- Cauchy_p_merger(sigma_p_values)

        # get individual posterior p_mu/p_sigma density estimates 
        # NOTE: this "works" but isn't analytic -- trying analytic below
        # mu_p_value_dt <- data.table('mu' = mu_p_values)
        # mu_dens_data <- ggplot_build(
        #             ggplot(mu_p_value_dt, aes(x = mu)) +
        #                 stat_bin(aes(y = after_stat(density)),
        #                          breaks = seq(0, 1, by = 0.01)
        #                         )
        #             )$data[[1]]
        # mu_dens_data <- data.table(mu_dens_data[, c('xmin', 'density')])
        # setnames(mu_dens_data, c('x', 'density'))
        # mu_dens_data <- rbindlist( list( data.table('x' = 0, 'density' = 0), mu_dens_data, data.table('x' = 1, 'density' = 0) ) )
        # p_mu_density_mat[i, ] <- mu_dens_data$density

        # ### TEST: try analytic p_mu density calculation ###
        # TODO: this seems wrong but maybe isn't
        # p_mu_density <- sapply(seq(0, 1, by = 0.01), FUN = Newcomb_NormalGamma_f_p_mu_given_y, y = y, lambda_samples = prior_lambda_samples, mu_0 = mu_0, kappa_0 = kappa_0)
        p_mu_density <- sapply(seq(0, 1, by = 0.01), FUN = Newcomb_NormalGamma_f_p_mu_given_y, y = y, lambda_samples = posterior_samples$lambda, mu_0 = mu_0, kappa_0 = kappa_0)
        p_mu_density[1] <- 0
        p_mu_density_mat[i, ] <- p_mu_density
        # mu_density_mat[i, ] <- mu_density


        # ##################################
        # ## test p_mu | lambda,y density ##
        # ##################################
        # lambda_posterior <- lambda_samples[1]
        # mu_posterior_dist_given_lambda_y <- rnorm(100000, mean = mu_n, sd = (kappa_n*lambda)^(-1/2))
        # U_mu_posterior_dist_given_lambda_y <- pnorm(mu_posterior_dist_given_lambda_y, mean = mu_0, sd = lambda_posterior^(-1/2))
        # p_mu_posterior_dist_given_lambda_y <- 2*pmin(U_mu_posterior_dist_given_lambda_y, 1-U_mu_posterior_dist_given_lambda_y)
        
        # # NOTE: I believe analytic density for mu | lambda,y is correct
        # mu_posterior_analytic_density <- sapply(seq(-2.5, 2.5, by = 0.01), FUN = Newcomb_NormalGamma_f_mu_given_lambda_y, y = y, lambda = lambda_posterior, mu_0 = mu_0, kappa_0 = kappa_0)
        
        # # this is wrong
        # p_mu_posterior_analytic_density <- sapply(seq(0, 1, by = 0.01), FUN = Newcomb_NormalGamma_f_p_mu_given_lambda_y, y = y, lambda = lambda_posterior, mu_0 = mu_0, kappa_0 = kappa_0)

        # NOTE: this "works" but isn't analytic -- trying analytic below
        # sigma_p_value_dt <- data.table('sigma' = sigma_p_values)
        # sigma_dens_data <- ggplot_build(
        #             ggplot(sigma_p_value_dt, aes(x = sigma)) +
        #                 stat_bin(aes(y = after_stat(density)),
        #                          breaks = seq(0, 1, by = 0.01)
        #                         )
        #             )$data[[1]]
        # sigma_dens_data <- data.table(sigma_dens_data[, c('xmin', 'density')])
        # setnames(sigma_dens_data, c('x', 'density'))
        # sigma_dens_data <- rbindlist( list( data.table('x' = 0, 'density' = 0), sigma_dens_data, data.table('x' = 1, 'density' = 0) ) )
        # p_sigma_density_mat[i, ] <- sigma_dens_data$density

        # TODO: this may or may not work -- testing now
        p_sigma_density <- sapply(seq(0, 1, by = 0.01), FUN = Newcomb_NormalGamma_f_p_lambda_given_y, alpha_I = alpha_n, beta_I = beta_n)
        p_sigma_density[1] <- 0
        p_sigma_density_mat[i, ] <- p_sigma_density

        # get p_data density
        data_p_value_dt <- data.table('data' = resid_p_values)
        data_dens_data <- ggplot_build(
                    ggplot(data_p_value_dt, aes(x = data)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
        data_dens_data <- data.table(data_dens_data[, c('xmin', 'density')])
        setnames(data_dens_data, c('x', 'density'))
        data_dens_data <- rbindlist( list( data.table('x' = 0, 'density' = 0), data_dens_data, data.table('x' = 1, 'density' = 0) ) )
        p_data_density_mat[i, ] <- data_dens_data$density

        # get posterior predictive values
        # NOTE: converting from completely empirical estimate to partially analytic estimate
        

        # min_posterior_predictive_ps[i] <- approx_min_posterior_predictive_p(posterior_samples$mu, posterior_samples$sigma, length(y), observed_min)

        # # TODO: finish this! want to have running average over both of these (and maybe p_T as well -- need to revisit)
        # plot_boolean <- ifelse(i == 1, TRUE, FALSE)
        # density_dt_list <- approximated_Newcomb_p_value_densities(newcomb$V1, posterior_samples$mu, posterior_samples$sigma, 0, 1, 1, p_val_dir, plot_boolean)
        # p_mu_density_dt <- density_dt_list[[1]]
        # p_sigma_density_dt <- density_dt_list[[2]]
        # p_mu_density_mat[i, ] <- p_mu_density_dt[, density]
        # p_sigma_density_mat[i, ] <- p_sigma_density_dt[, density]
    }

    # get avg density for mu/sigma
    # p_mu_avg_density <- p_mu_avg_density / n_sims
    # p_sigma_avg_density <- p_sigma_avg_density / n_sims
    # p_mu_avg_density_dt <- data.table('p_mu' = x_vals, 'density' = p_mu_avg_density)
    # p_sigma_avg_density_dt <- data.table('p_sigma' = x_vals, 'density' = p_sigma_avg_density)
    # p_mu_avg_density_dt[, standardized_density := density / max(density)]
    # p_sigma_avg_density_dt[, standardized_density := density / max(density)]

    # save results dts
    p_val_dt <- data.table('mu' = mu_ps, 'sigma' = sigma_ps, 
                           'data' = resid_ps, 'min_posterior_predictive' = min_posterior_predictive_ps
                           )

    single_sample_p_val_dt <- data.table('mu' = single_sample_mu_ps, 'sigma' = single_sample_sigma_ps, 
                                         'data' = single_sample_resid_ps, 'min_posterior_predictive' = min_posterior_predictive_ps
                                        )

    fwrite( p_val_dt, p_val_path )
    fwrite( total_resid_ps, total_data_p_val_path )
    fwrite( total_mu_ps, total_mu_p_val_path )
    fwrite( total_sigma_ps, total_sigma_p_val_path )
    fwrite( single_sample_p_val_dt, single_sample_p_val_path )
    fwrite( y_mat, y_path )
    write.table( p_mu_density_mat, p_mu_density_mat_path, row.names = FALSE, col.names = FALSE )
    write.table( p_sigma_density_mat, p_sigma_density_mat_path, row.names = FALSE, col.names = FALSE )
    write.table( p_data_density_mat, p_data_density_mat_path, row.names = FALSE, col.names = FALSE )
}

Newcomb_p_value_histograms(single_sample_p_val_dt, p_val_dir, plot_min_posterior_predictive = TRUE)
Newcomb_p_value_densities(single_sample_p_val_dt, p_val_dir, plot_min_posterior_predictive = TRUE)

# ### plot p-values ###

# p_val_dt_long <- melt(p_val_dt)
# single_sample_p_val_dt_long <- melt(single_sample_p_val_dt)

# var_names <- unique(p_val_dt_long$variable)
# presentation_var_names <- c(
#                              TeX('$p_{\\mu}$'),
#                              TeX('$p_{\\sigma}$'),
#                              expression('p'['data']),
#                              expression('p'['T'])
#                            )

# for (i in 1:length(var_names)) {
#     var_name <- var_names[i]
#     presentation_var_name <- presentation_var_names[i]

#     temp_dt <- p_val_dt_long[variable == var_name]
#     dens_data <- ggplot_build(
#                     ggplot(temp_dt, aes(x = value)) +
#                         stat_bin(aes(y = after_stat(density)),
#                                  breaks = seq(0, 1, by = 0.01)
#                                 # bins = 30
#                                 )
#                     )$data[[1]]
#     dens_data <- dens_data[, c('xmin', 'density')]
#     # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
#     dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
#                                  dens_data,
#                                  data.frame(x = 1, density = 0))
#                           )

#     indiv_p_val_plot <- ggplot(dens_data, aes(x = x, y = density)) + 
#                 geom_line(color = 'blue', size = 1) + 
#                 coord_cartesian(xlim = c(0, 1)) +
#                 labs(x = presentation_var_name, y = 'Density') + 
#                 theme(
#                     axis.title.x = element_text(size = rel(3.0)),
#                     axis.text.x = element_text(size = rel(3.0)),
#                     axis.title.y = element_text(size = rel(3.0)),
#                     axis.text.y = element_text(size = rel(3.0)),
#                     panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                     panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                     panel.background = element_rect(fill = "white") # White background
#                 )
#     ggsave( file.path(p_val_dir, glue::glue('{var_name}_p_value_density.png')), indiv_p_val_plot, width = 12, height = 8, dpi = 300 )

#     indiv_p_val_plot <- ggplot(p_val_dt_long[variable == var_name], aes(x = value)) + 
#                             geom_histogram(position = 'identity', binwidth = 0.01, boundary = 0, colour = 'black', fill = 'blue') +  
#                             coord_cartesian(xlim = c(0,1)) +
#                             labs(x = presentation_var_name) +
#                             theme(
#                                     axis.title.x = element_text(size = rel(3.0)),
#                                     axis.text.x = element_blank(),
#                                     axis.ticks.x = element_blank(),
#                                     axis.title.y = element_blank(),
#                                     axis.text.y = element_blank(),
#                                     axis.ticks.y = element_blank(),
#                                     panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                                       panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                                       panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                                       panel.background = element_rect(fill = "white") # White background
#                                     ) 
#     ggsave( file.path(p_val_dir, glue::glue('{var_name}_p_value_histogram.png')), indiv_p_val_plot, width = 12, height = 8, dpi = 300 )

#     indiv_p_val_plot <- ggplot(single_sample_p_val_dt_long[variable == var_name], aes(x = value)) + 
#                             stat_bin(
#                                         aes(y = after_stat(density)), position = 'identity',
#                                         # bins = 30,
#                                         breaks = seq(0, 1, by = 0.01),
#                                         # binwidth = 0.01, boundary = 0,
#                                          geom = "line",
#                                         color = "blue",
#                                         size = 1
#                                     ) +
#                             coord_cartesian(xlim = c(0,1)) +
#                             labs(x = presentation_var_name, y = 'Density') +
#                             theme(
#                                     axis.title.x = element_text(size = rel(3.0)),
#                                     axis.text.x = element_blank(),
#                                     axis.ticks.x = element_blank(),
#                                     axis.title.y = element_blank(),
#                                     axis.text.y = element_blank(),
#                                     axis.ticks.y = element_blank(),
#                                     panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                                       panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                                       panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                                       panel.background = element_rect(fill = "white") # White background
#                                     ) 
#     ggsave( file.path(p_val_dir, glue::glue('single_sample_{var_name}_p_value_density.png')), indiv_p_val_plot, width = 12, height = 8, dpi = 300 )

#     indiv_p_val_plot <- ggplot(single_sample_p_val_dt_long[variable == var_name], aes(x = value)) + 
#                             geom_histogram(position = 'identity', binwidth = 0.01, boundary = 0, colour = 'black', fill = 'blue') +
#                             coord_cartesian(xlim = c(0,1)) +
#                             labs(x = presentation_var_name) +
#                             theme(
#                                     axis.title.x = element_text(size = rel(3.0)),
#                                     axis.text.x = element_blank(),
#                                     axis.ticks.x = element_blank(),
#                                     axis.title.y = element_blank(),
#                                     axis.text.y = element_blank(),
#                                     axis.ticks.y = element_blank(),
#                                     panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                                       panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                                       panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                                       panel.background = element_rect(fill = "white") # White background
#                                     ) 
#     ggsave( file.path(p_val_dir, glue::glue('single_sample_{var_name}_p_value_histogram.png')), indiv_p_val_plot, width = 12, height = 8, dpi = 300 )
# }

# p_val_plot <- ggplot(p_val_dt_long, aes(x = value)) + 
#                 geom_histogram(position = 'identity', binwidth = 0.01, boundary = 0, colour = 'black', fill = 'blue') +  
#                 coord_cartesian(xlim = c(0,1)) +
#                 labs(x = 'p-value') +
#                 theme(
#                         axis.title.x = element_text(size = rel(3.0)),
#                         axis.text.x = element_text(size = rel(3.0)),,
#                         # axis.ticks.x = element_blank(),
#                         axis.title.y = element_blank(),
#                         axis.text.y = element_blank(),
#                         axis.ticks.y = element_blank(),
#                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                         panel.background = element_rect(fill = "white") # White background
#                         ) + 
#                 facet_wrap(~variable, scales = 'free')
# ggsave(file.path(p_val_dir, 'full_p_value_plot_histogram.png'), p_val_plot, width = 12, height = 8, dpi = 300)

# p_val_plot <- ggplot(p_val_dt_long, aes(x = value)) + 
#                 stat_bin(
#                                         aes(
#                                             y = after_stat(density / max(density))),
#                                         position = 'identity',
#                                         # bins = 30,
#                                         # binwidth = 0.01, boundary = 0,
#                                         breaks = seq(0, 1, by = 0.01),
#                                         geom = "line",
#                                         color = "blue",
#                                         size = 1
#                                     ) +
#                 coord_cartesian(xlim = c(0,1)) +
#                 labs(x = 'p-value') +
#                 theme(
#                         axis.title.x = element_text(size = rel(3.0)),
#                         axis.text.x = element_text(size = rel(3.0)),,
#                         # axis.ticks.x = element_blank(),
#                         axis.title.y = element_blank(),
#                         axis.text.y = element_blank(),
#                         axis.ticks.y = element_blank(),
#                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                         panel.background = element_rect(fill = "white") # White background
#                         ) + 
#                 facet_wrap(~variable, scales = 'free')
# ggsave(file.path(p_val_dir, 'full_p_value_plot_density.png'), p_val_plot, width = 12, height = 8, dpi = 300)

# BDA_p_val_plot <- ggplot(p_val_dt_long[variable %in% c('min_posterior_predictive')], aes(x = value)) + 
#                 geom_histogram(position = 'identity', binwidth = 0.01, boundary = 0, colour = 'black', fill = 'blue') + 
#                 coord_cartesian(xlim = c(0,1)) +
#                 labs(x = 'p-value') +
#                 theme(
#                         axis.title.x = element_text(size = rel(3.0)),
#                         axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         axis.title.y = element_blank(),
#                         axis.text.y = element_blank(),
#                         axis.ticks.y = element_blank(),
#                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                         panel.background = element_rect(fill = "white") # White background
#                         )
# ggsave(file.path(p_val_dir, 'BDA_p_value_plot_histogram.png'), BDA_p_val_plot, width = 12, height = 8, dpi = 300)

# BDA_p_val_plot <- ggplot(p_val_dt_long[variable %in% c('min_posterior_predictive')], aes(x = value)) + 
#                     stat_bin(
#                                         aes(y = after_stat(density / max(density))), position = 'identity',
#                                         # bins = 30,
#                                         # binwidth = 0.01, boundary = 0,
#                                         breaks = seq(0, 1, by = 0.01),
#                                         geom = "line",
#                                         color = "blue",
#                                         size = 1
#                                     ) +
#                 coord_cartesian(xlim = c(0,1)) +
#                 labs(x = 'p-value') +
#                 theme(
#                         axis.title.x = element_text(size = rel(3.0)),
#                         axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         axis.title.y = element_blank(),
#                         axis.text.y = element_blank(),
#                         axis.ticks.y = element_blank(),
#                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                         panel.background = element_rect(fill = "white") # White background
#                         )
# ggsave(file.path(p_val_dir, 'BDA_p_value_plot_density.png'), BDA_p_val_plot, width = 12, height = 8, dpi = 300)

# p_val_ecdf <- ggplot(p_val_dt_long, aes(x = value)) + 
#                 stat_ecdf(position = 'identity', linewidth = 1, alpha = 1, color = 'blue') +
#                 geom_abline(intercept = 0, slope = 1, colour = 'black') +
#                 coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
#                 labs(x = 'p-value') +
#                 facet_wrap(~variable) + 
#                 theme(
#                         axis.title.x = element_text(size = rel(3.0)),
#                         axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         axis.title.y = element_blank(),
#                         axis.text.y = element_blank(),
#                         axis.ticks.y = element_blank(),
#                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                         panel.background = element_rect(fill = "white") # White background
#                         )
# ggsave(file.path(p_val_dir, 'p_value_ecdf.png'), p_val_ecdf, width = 12, height = 8, dpi = 300)

# BDA_p_val_ecdf <- ggplot(p_val_dt_long[variable %in% c('min_posterior_predictive')], aes(x = value)) + 
#                 stat_ecdf(position = 'identity', linewidth = 1, alpha = 1, color = 'blue') +
#                 geom_abline(intercept = 0, slope = 1, colour = 'black') +
#                 coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
#                 labs(x = 'p-value') +
#                 theme(
#                         axis.title.x = element_text(size = rel(3.0)),
#                         axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         axis.title.y = element_blank(),
#                         axis.text.y = element_blank(),
#                         axis.ticks.y = element_blank(),
#                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                         panel.background = element_rect(fill = "white") # White background
#                         )
# ggsave(file.path(p_val_dir, 'BDA_p_value_ecdf.png'), BDA_p_val_ecdf, width = 12, height = 8, dpi = 300)

# single_sample_p_val_plot <- ggplot(single_sample_p_val_dt_long, aes(x = value)) + 
#                 stat_bin(
#                                         aes(y = after_stat(density / max(density))), position = 'identity',
#                                         binwidth = 0.01, boundary = 0,
#                                         geom = "line",
#                                         color = "blue",
#                                         size = 1
#                                     ) +
#                 labs(x = 'p-value') +
#                 theme(
#                         axis.title.x = element_text(size = rel(3.0)),
#                         axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         axis.title.y = element_blank(),
#                         axis.text.y = element_blank(),
#                         axis.ticks.y = element_blank(),
#                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                         panel.background = element_rect(fill = "white") # White background
#                         ) + 
#                 facet_wrap(~variable, scales = 'free')
# ggsave(file.path(p_val_dir, 'full_single_sample_p_value_plot_density.png'), single_sample_p_val_plot, width = 12, height = 8, dpi = 300)

# single_sample_p_val_plot <- ggplot(single_sample_p_val_dt_long, aes(x = value)) + 
#                 geom_histogram(position = 'identity', binwidth = 0.01, boundary = 0, colour = 'black', fill = 'blue') +
#                 labs(x = 'p-value') +
#                 theme(
#                         axis.title.x = element_text(size = rel(3.0)),
#                         axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         axis.title.y = element_blank(),
#                         axis.text.y = element_blank(),
#                         axis.ticks.y = element_blank(),
#                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                         panel.background = element_rect(fill = "white") # White background
#                         ) + 
#                 facet_wrap(~variable, scales = 'free')
# ggsave(file.path(p_val_dir, 'full_single_sample_p_value_plot_histogram.png'), single_sample_p_val_plot, width = 12, height = 8, dpi = 300)

# single_sample_p_val_ecdf <- ggplot(single_sample_p_val_dt_long, aes(x = value)) + 
#                 stat_ecdf(position = 'identity', linewidth = 1, alpha = 1, color = 'blue') +
#                 geom_abline(intercept = 0, slope = 1, colour = 'black') +
#                 coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
#                 labs(x = 'p-value') +
#                 facet_wrap(~variable) + 
#                 theme(
#                         axis.title.x = element_text(size = rel(3.0)),
#                         axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank(),
#                         axis.title.y = element_blank(),
#                         axis.text.y = element_blank(),
#                         axis.ticks.y = element_blank(),
#                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                         panel.background = element_rect(fill = "white") # White background
#                         )
# ggsave(file.path(p_val_dir, 'single_sample_p_value_ecdf.png'), single_sample_p_val_ecdf, width = 12, height = 8, dpi = 300)

# histograms/densities across individual posteriors
# n_posterior_samples_plotting <- 50
# samples <- sample(1:nrow(total_mu_ps), size = n_posterior_samples_plotting, replace = FALSE)
# colors <- hue_pal()(n_posterior_samples_plotting) 

### mu ###
p_mu_analytic_dt <- data.table(p_mu_density_mat)
colnames(p_mu_analytic_dt) <- as.character(as.numeric(str_replace_all(colnames(p_mu_analytic_dt), 'V', '')))
p_mu_analytic_dt[, sim_sample := as.character( 1:nrow(p_mu_analytic_dt) )]
p_mu_analytic_dt_long <- melt( p_mu_analytic_dt, variable.name = 'p_mu', value.name = 'density', id.vars = 'sim_sample' )
# p_mu_analytic_dt_long[, p_mu := pmax(0, as.numeric(p_mu)/100 - 0.02)] # empirical only 
p_mu_analytic_dt_long[, p_mu := as.numeric(p_mu)/100 - 0.01]
p_mu_analytic_dt_long[, standardized_density := density / max(density), by = sim_sample]

avg_densities_analytic <- data.table(colMeans(p_mu_density_mat))
avg_densities_analytic[, density := V1]
# avg_densities_analytic[, standardized_density := V1 / max(V1)]
# avg_densities_analytic[, p_mu := c(0, seq(0, 1, 0.01))]
avg_densities_analytic[, p_mu := seq(0, 1, 0.01)]
avg_densities_analytic <- avg_densities_analytic[! (p_mu %in% c(0,1))] # remove points at 0,1

mu_indiv_posterior_density_plot <-  ggplot(p_mu_analytic_dt_long[sim_sample %in% samples], aes(x = p_mu, y = standardized_density, colour = sim_sample, group = sim_sample)) + 
                                        geom_line(position = 'identity', alpha = 0.5, size = 1) +
                                        scale_color_manual(values = colors) +
                                        geom_line(data = avg_densities_analytic, 
                                                  aes(x = p_mu, y = density, group = 1), position = 'identity', color = 'black', size = 1) +
                                        labs(x = TeX('$p_{\\mu}$'), y = 'standardized density') +
                                        theme(
                                                        legend.position = 'none',
                                                        axis.title.x = element_text(size = rel(3.0)),
                                                        axis.text.x = element_text(size = rel(3.0)),
                                                        # axis.ticks.x = element_text(size = rel(3.0)),
                                                        axis.title.y = element_text(size = rel(3.0)),
                                                        axis.text.y = element_text(size = rel(3.0)),
                                                        # axis.ticks.y = element_blank(),
                                                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                                        panel.background = element_rect(fill = "white") # White background
                                                        ) 
                                                        
ggsave(file.path(p_val_dir, 'mu_individual_posterior_density_plot.png'), mu_indiv_posterior_density_plot, width = 12, height = 8, dpi = 300)

### sigma ###
p_sigma_analytic_dt <- data.table(p_sigma_density_mat)
colnames(p_sigma_analytic_dt) <- as.character(as.numeric(str_replace_all(colnames(p_sigma_analytic_dt), 'V', '')))
p_sigma_analytic_dt[, sim_sample := as.character( 1:nrow(p_sigma_analytic_dt) )]
p_sigma_analytic_dt_long <- melt( p_sigma_analytic_dt, variable.name = 'p_sigma', value.name = 'density', id.vars = 'sim_sample' )
# p_sigma_analytic_dt_long[, p_sigma := pmax(0, as.numeric(p_sigma)/100 - 0.02)]
p_sigma_analytic_dt_long[, p_sigma := as.numeric(p_sigma)/100 - 0.01]
p_sigma_analytic_dt_long[, standardized_density := density / max(density), by = sim_sample]

avg_densities_analytic <- data.table(colMeans(p_sigma_density_mat))
avg_densities_analytic[, density := V1]
# avg_densities_analytic[, standardized_density := V1 / max(V1)]
# avg_densities_analytic[, p_sigma := c(0, seq(0, 1, 0.01))]
avg_densities_analytic[, p_sigma := seq(0, 1, 0.01)]
avg_densities_analytic <- avg_densities_analytic[! (p_sigma %in% c(0,1) )] # remove points at 0,1

sigma_indiv_posterior_density_plot <-  ggplot(p_sigma_analytic_dt_long[sim_sample %in% samples], aes(x = p_sigma, y = standardized_density, colour = sim_sample, group = sim_sample)) + 
                                        geom_line(position = 'identity', alpha = 0.5, size = 1) +
                                        scale_color_manual(values = colors) +
                                        geom_line(data = avg_densities_analytic, 
                                                  aes(x = p_sigma, y = density, group = 1), position = 'identity', color = 'black', size = 1) +
                                        labs(x = TeX('$p_{\\sigma}$'), y = 'standardized density') +
                                        theme(
                                                        legend.position = 'none',
                                                        axis.title.x = element_text(size = rel(3.0)),
                                                        axis.text.x = element_text(size = rel(3.0)),
                                                        # axis.ticks.x = element_text(size = rel(3.0)),
                                                        axis.title.y = element_text(size = rel(3.0)),
                                                        axis.text.y = element_text(size = rel(3.0)),
                                                        # axis.ticks.y = element_blank(),
                                                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                                        panel.background = element_rect(fill = "white") # White background
                                                        ) 
                                                        
ggsave(file.path(p_val_dir, 'sigma_individual_posterior_density_plot.png'), sigma_indiv_posterior_density_plot, width = 12, height = 8, dpi = 300)


# sigma_p_dt <- data.table( t(total_sigma_ps) )
# sigma_p_dt_long <- melt(sigma_p_dt, variable.name = 'sim_sample')
# sigma_p_dt_long[, sim_sample := str_replace_all(sim_sample, 'V', '')]

# # get averaged densities
# all_densities_plot <- ggplot(sigma_p_dt_long, aes(x = value, group = sim_sample)) + 
#                         stat_bin(aes(y = after_stat(density)), binwidth = 0.01, boundary = 0)
# all_densities_data <- data.table( ggplot_build(all_densities_plot)$data[[1]] )
# # avg_densities <- all_densities_data[, mean(density), by = x]
# # avg_densities[, standardized_density := V1 / max(V1)]
# avg_densities <- data.table(colMeans(p_sigma_density_mat))
# avg_densities[, standardized_density := V1 / max(V1)]
# avg_densities[, p_sigma := seq(0, 1, 0.01)]

# sigma_indiv_posterior_density_plot <- ggplot(sigma_p_dt_long[sim_sample %in% samples,], aes(x = value, colour = sim_sample, group = sim_sample)) + 
#                                         stat_bin(aes(y = after_stat(ndensity)), position = 'identity',
#                                                                         binwidth = 0.01, boundary = 0,
#                                                                         alpha = 0.5,
#                                                                         geom = "line",
#                                                                         size = 1
#                                                                     ) +
#                                         scale_color_manual(values = colors) +
#                                         geom_line(data = avg_densities, 
#                                                   aes(x = p_sigma, y = standardized_density, group = 1), position = 'identity', color = 'black', size = 1) +
#                                         labs(x = TeX('$p_{\\sigma}$'), y = 'Standardized Density') +
#                                                 theme(
#                                                         legend.position = 'none',
#                                                         axis.title.x = element_text(size = rel(3.0)),
#                                                         axis.text.x = element_text(size = rel(3.0)),
#                                                         # axis.ticks.x = element_text(size = rel(3.0)),
#                                                         axis.title.y = element_text(size = rel(3.0)),
#                                                         axis.text.y = element_text(size = rel(3.0)),
#                                                         # axis.ticks.y = element_blank(),
#                                                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                                                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                                                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                                                         panel.background = element_rect(fill = "white") # White background
#                                                         ) 
# ggsave(file.path(p_val_dir, 'sigma_individual_posterior_density_plot.png'), sigma_indiv_posterior_density_plot, width = 12, height = 8, dpi = 300)

# # data
# data_p_dt <- data.table( t(total_resid_ps) )
# data_p_dt <- data_p_dt[, sample(1:ncol(data_p_dt), 250000, replace = FALSE), with = FALSE]
# data_p_dt_long <- melt(data_p_dt, variable.name = 'sim_sample')
# data_p_dt_long[, sim_sample := str_replace_all(sim_sample, 'V', '')]
# samples <- sample(unique(data_p_dt_long$sim_sample), size = n_posterior_samples_plotting, replace = FALSE)


# # get averaged densities
# all_densities_plot <- ggplot(data_p_dt_long, aes(x = value, group = sim_sample)) + 
#                         stat_bin(aes(y = after_stat(density)), binwidth = 0.02, boundary = 0)
# all_densities_data <- data.table( ggplot_build(all_densities_plot)$data[[1]] )
# avg_densities <- all_densities_data[, mean(density), by = x]
# avg_densities[, density := V1]
# # avg_densities[, standardized_density := V1 / max(V1)]
# avg_densities[, xmin := x - 0.01]
# # avg_densities <- avg_densities[, c('xmin', 'standardized_density')]
# # avg_densities <- rbindlist( list(
# #                             data.table('xmin' = 0, 'standardized_density' = 0),
# #                             avg_densities,
# #                             data.table('xmin' = 1, 'standardized_density' = 0) 
# #                           )
# # )

# data_indiv_posterior_density_plot <- ggplot(data_p_dt_long[sim_sample %in% samples,], aes(x = value, colour = sim_sample, group = sim_sample)) + 
#                                         stat_bin(aes(y = after_stat(ndensity)), position = 'identity',
#                                                                         binwidth = 0.01, boundary = 0,
#                                                                         alpha = 0.50,
#                                                                         geom = "line",
#                                                                         # color = "blue",
#                                                                         size = 1
#                                                                     ) +
#                                         scale_color_manual(values = colors) +
#                                         geom_line(data = avg_densities, 
#                                                   aes(x = xmin, y = density, group = 1), position = 'identity', color = 'black', size = 1) +
#                                         labs(x = expression('p'['data,unif']), y = 'Standardized Density') +
#                                                 theme(
#                                                         legend.position = 'none',
#                                                         axis.title.x = element_text(size = rel(3.0)),
#                                                         axis.text.x = element_text(size = rel(3.0)),
#                                                         # axis.ticks.x = element_text(size = rel(3.0)),
#                                                         axis.title.y = element_text(size = rel(3.0)),
#                                                         axis.text.y = element_text(size = rel(3.0)),
#                                                         # axis.ticks.y = element_blank(),
#                                                         panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                                                         panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                                                         panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                                                         panel.background = element_rect(fill = "white") # White background
#                                                         ) 
# ggsave(file.path(p_val_dir, 'data_individual_posterior_density_plot.png'), data_indiv_posterior_density_plot, width = 12, height = 8, dpi = 300)



##################################################

p_data_analytic_dt <- data.table(p_data_density_mat)
colnames(p_data_analytic_dt) <- as.character(as.numeric(str_replace_all(colnames(p_data_analytic_dt), 'V', '')))
p_data_analytic_dt[, sim_sample := as.character( 1:nrow(p_data_analytic_dt) )]
p_data_analytic_dt_long <- melt( p_data_analytic_dt, variable.name = 'p_data', value.name = 'density', id.vars = 'sim_sample' )
p_data_analytic_dt_long[, p_data := pmax(0, as.numeric(p_data)/100 - 0.02)] # empirical only 
# p_data_analytic_dt_long[, p_data := as.numeric(p_data)/100 - 0.01]
p_data_analytic_dt_long[, standardized_density := density / max(density), by = sim_sample]

avg_densities_analytic <- data.table(colMeans(p_data_density_mat))
avg_densities_analytic[, density := V1]
# avg_densities_analytic[, standardized_density := V1 / max(V1)]
avg_densities_analytic[, p_data := c(0, seq(0, 1, 0.01))]
# avg_densities_analytic[, p_data := seq(0, 1, 0.01)]
avg_densities_analytic <- avg_densities_analytic[! (p_data %in% c(0,1))] # remove points at 0,1

data_indiv_posterior_density_plot <-  ggplot(p_data_analytic_dt_long[sim_sample %in% samples], aes(x = p_data, y = standardized_density, colour = sim_sample, group = sim_sample)) + 
                                        geom_line(position = 'identity', alpha = 0.5, size = 1) +
                                        scale_color_manual(values = colors) +
                                        geom_line(data = avg_densities_analytic, 
                                                  aes(x = p_data, y = density, group = 1), position = 'identity', color = 'black', size = 1) +
                                        labs(x = expression('p'['data,unif']), y = 'standardized density') +
                                        theme(
                                                        legend.position = 'none',
                                                        axis.title.x = element_text(size = rel(3.0)),
                                                        axis.text.x = element_text(size = rel(3.0)),
                                                        # axis.ticks.x = element_text(size = rel(3.0)),
                                                        axis.title.y = element_text(size = rel(3.0)),
                                                        axis.text.y = element_text(size = rel(3.0)),
                                                        # axis.ticks.y = element_blank(),
                                                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                                        panel.background = element_rect(fill = "white") # White background
                                                        ) 
                                                        
ggsave(file.path(p_val_dir, 'data_individual_posterior_density_plot.png'), data_indiv_posterior_density_plot, width = 12, height = 8, dpi = 300)


##################################################












# get results from correct model
u_resid_mat <- as.matrix(u_resid_mat)
resid_p_values <- uniformity_test(u_resid_mat)
resid_p_value <- Cauchy_p_merger(resid_p_values)
p_val_string <- signif(resid_p_value, 2)
one_posterior_sample_data_u_histogram(u_resid_mat, p_val_dir)
one_posterior_sample_data_u_density(u_resid_mat, p_val_dir)
one_posterior_sample_residual_ecdf_plot(u_resid_mat, p_val_dir)
one_posterior_sample_residual_tilted_ecdf_plot(u_resid_mat, p_val_dir)
residual_ecdf_plot(u_resid_mat, 50, p_val_string, p_val_dir)
residual_tilted_ecdf_plot(u_resid_mat, 50, p_val_string, p_val_dir, 0.35)