set.seed(1)

setwd('/Users/christiancovington/research/uniform_parameterization_checks/code/bernoulli')

library(pacman)
pacman::p_load('rstan', 'data.table', 'ggplot2', 'forcats', 'energy', 
               'gridExtra', 'EnvStats', 'stringr', 'Hmisc', 'lmtest',
               'latex2exp', 'glue', 'scales', 'invgamma')

source('../utils/generate_plots.R')
source('../utils/testing.R')
source('../utils/sim_data.R')

build_hoeffD_null <- function(n, n_sims) {
    null <- rep(0, n_sims)
    for (i in 1:n_sims) {
        # if (i %% 10000 == 0) {
        #     print(glue::glue('run {i} of {n_sims}'))
        # }
        x <- runif(n)
        null[i] <- hoeffd(x[1:(n-1)], x[2:n])$D[1,2]
    }
    return(null)
}

# parent directories
data_dir <- file.path('..', '..', 'data', 'bernoulli')
results_dir <- file.path('..', '..', 'results', 'bernoulli')

# dgps <- c('sim', 'BDA')

dgps <- c('sim')
priors <- c('uniform_prior', 
            'uninformative_prior', 
            'very_bad_prior')

null_dist_66 <- build_hoeffD_null(66, 10^6)
# null_dist_100 <- build_hoeffD_null(100, 10^6)

for (dgp in dgps) {
    # read bernoulli data
    if (dgp == 'sim') {
        if (file.exists(file.path(data_dir, 'sim.txt'))) {
            bernoulli <- fread( file.path(data_dir, 'sim.txt') )
        } else {
            bernoulli <- sim_bernoulli(n = 100, p = 0.5, p_repeat = 0.8)
            bernoulli <- data.table(bernoulli)
            colnames(bernoulli) <- 'V1'
            fwrite( bernoulli, file.path(data_dir, 'sim.txt') )
        }
    } else {
        bernoulli <- fread( file.path(data_dir, 'bernoulli.txt') )
    }

    # record number of times vector switches values
    true_switches <- sum(abs(diff(bernoulli$V1)))
    
    for (prior in priors) {
        print(glue('dgp: {dgp}, prior: {prior}'))
        # load model from stan file
        model.stan <- stan_model(file = paste0(prior, '.stan'))

        # set up directory structure
        prior_results_dir <- file.path(results_dir, dgp, prior)
        dir.create(prior_results_dir, showWarnings = FALSE, recursive = TRUE)
    
        stan_data <- list(n = nrow(bernoulli),
                            y = bernoulli$V1)
        
        # fit model in stan
        n_burnin <- 1000
        n_posterior_samples <- 10^6
        n_chains <- 1

        dap_samples_file <- file.path(prior_results_dir, 'dap_samples.csv')
        if (file.exists(dap_samples_file)) {
            dap_samples <- readRDS(dap_samples_file)
        } else {
            # fit model
            model.fit <- sampling(model.stan, data = stan_data, iter = n_burnin+n_posterior_samples, 
                                    warmup = n_burnin, chains = n_chains, 
                                    verbose = FALSE, refresh = 0)
            
            # extract posterior estimates 
            dap_samples <- extract(model.fit)
            saveRDS(dap_samples, dap_samples_file)
        }
        
        theta_samples <- dap_samples$theta
        parameter_us <- data.table(dap_samples$u)
        y_pred_samples <- dap_samples$y_pred # get posterior predictive values
        names(parameter_us) <- c('u1')

        # get number of switches for each posterior predictive vector
        y_pred_switches <- apply(y_pred_samples, 1, function(row) sum(abs(diff(row))) )
        y_pred_switch_dt <- data.table('y_pred_switches' = y_pred_switches)
        bernoulli_switch_histogram(y_pred_switch_dt, true_switches, prior_results_dir)

        # generate residual u values
        u_resid_mat_file <- file.path(prior_results_dir, 'u_resid_mat.csv')
        if (file.exists(u_resid_mat_file)) {
            u_resid_mat <- as.matrix( fread(u_resid_mat_file) )
        } else {
            u_resid_mat <- matrix(0, nrow(parameter_us), stan_data$n)
            for (i in 1:nrow(u_resid_mat)) {
                for (j in 1:ncol(u_resid_mat)) {
                    if (bernoulli[j, 'V1'] == 1) {
                        u_resid_mat[i,j] <- runif(1,1-theta_samples[i],1)
                    } else {
                        u_resid_mat[i,j] <- runif(1,0,1-theta_samples[i])
                    }
                }
            }
            fwrite(u_resid_mat, u_resid_mat_file)
        }

        # generate residual trajectory plot
        # residual_trajectory_plot(u_resid_mat, prior_results_dir)

        # generate histograms of parameter u values
        parameter_names <- c('theta')
        LaTeX_parameter_names <- c('$\\theta$')
        parameter_hist_plots(parameter_us, parameter_names, LaTeX_parameter_names, prior_results_dir)

        ### generate test results ###
        p_val_dt_file <- file.path(prior_results_dir, 'p_val_dt.csv')
        if (file.exists(p_val_dt_file)) {
            p_val_dt <- fread(p_val_dt_file)
        } else {
            # generate uniformity p-values for each posterior sample
            resid_p_values <- uniformity_test(u_resid_mat) # p-values for resid uniformity testing
            theta_p_values <- u_to_p_converter(parameter_us$u1) # convert parameter u values to p values

            # test for independence of residual p-values
            indep_p_values <- resid_independence_test(u_resid_mat, null_dist_66)

            # plot posterior p-values
            p_val_dt <- data.table('theta' = theta_p_values, 
                                'data' = resid_p_values,
                                'indep' = indep_p_values
                                )
        }
        
        if (prior == 'uniform_prior') {
            beta_a <- 1
            beta_b <- 1
        } else if (prior == 'uninformative_prior') {
            beta_a <- 1/2
            beta_b <- 1/2
        } else if (prior == 'very_bad_prior') {
            beta_a <- 1
            beta_b <- 50
        }

        bernoulli_p_value_densities_single_posterior(p_val_dt, bernoulli$V1, beta_a, beta_b, prior_results_dir, plot_switch = FALSE)

        # combine p-values using Cauchy method
        resid_p_value <- Cauchy_p_merger(resid_p_values)
        theta_p_value <- Cauchy_p_merger(theta_p_values)
        indep_p_value <- Cauchy_p_merger(indep_p_values)

        # generate plots of data u values
        u_resid_mat_sub <- u_resid_mat[sample(1:nrow(u_resid_mat), size = 1000, replace = FALSE), ]
        p_val_string <- signif(resid_p_value, 2)
        one_posterior_sample_data_u_density(u_resid_mat, prior_results_dir)
        one_posterior_sample_data_u_histogram(u_resid_mat, prior_results_dir)
        one_posterior_sample_residual_ecdf_plot(u_resid_mat, prior_results_dir)
        one_posterior_sample_residual_tilted_ecdf_plot(u_resid_mat, prior_results_dir)
        residual_ecdf_plot(u_resid_mat_sub, 50, p_val_string, prior_results_dir)
        residual_tilted_ecdf_plot(u_resid_mat_sub, 50, p_val_string, prior_results_dir, 0.35)

        test_file <- file.path(prior_results_dir, 'test_output.txt')
        output <- glue::glue('residual uniformity test: p = {resid_p_value}\ntheta uniformity test: p = {theta_p_value}\nresidual independence test: p = {indep_p_value}')
        writeLines(output, test_file)
    }
}


###########################################################
## simulate from correct model, observe frequentist properties
## over large number of simulations 
###########################################################

set.seed(1)

n <- 100
n_sims <- 10^5
n_sims_string <- glue::glue('n_sims={n_sims}')
dgp <- 'sim'
prior <- 'correct_prior'
null_dist_100 <- build_hoeffD_null(100, 10^6)

# # load model from stan file
# model.stan <- stan_model(file = paste0(prior, '.stan'))

# set up directory structure
prior_results_dir <- file.path(results_dir, dgp, prior, n_sims_string)
dir.create(prior_results_dir, showWarnings = FALSE, recursive = TRUE)
p_val_path <- file.path(prior_results_dir, 'p_val_dt.csv')
total_theta_p_val_path <- file.path(prior_results_dir, 'total_theta_p_val_dt.csv')
total_unif_p_val_path <- file.path(prior_results_dir, 'total_unif_p_val_dt.csv')
total_indep_p_val_path <- file.path(prior_results_dir, 'total_indep_p_val_dt.csv')
p_theta_density_mat_path <- file.path(prior_results_dir, 'p_theta_density_mat.csv')
p_unif_density_mat_path <- file.path(prior_results_dir, 'p_unif_density_mat.csv')
p_indep_density_mat_path <- file.path(prior_results_dir, 'p_indep_density_mat.csv')

theta_0_values <- rep(0, n_sims)
y_values <- matrix(0, n_sims, n)
theta_p_values <- rep(0, n_sims)
switch_p_values <- rep(0, n_sims)
unif_switch_p_values <- rep(0, n_sims)
data_unif_p_values <- rep(0, n_sims)
data_indep_p_values <- rep(0, n_sims)

# store posterior p values 
n_posterior_sims <- 50
total_theta_ps <- matrix(0, n_sims, n_posterior_sims)
total_unif_ps <- matrix(0, n_sims, n_posterior_sims)
total_indep_ps <- matrix(0, n_sims, n_posterior_sims)
y_mat <- matrix(0, n_sims, n)

x_vals <- seq(0, 1, by = 0.01)
p_theta_density_mat <- matrix(0, n_sims, length(x_vals)+1)
p_unif_density_mat <- matrix(0, n_sims, length(x_vals)+1)
p_indep_density_mat <- matrix(0, n_sims, length(x_vals)+1)

n_posterior_samples_plotting <- 20
samples <- sample(1:nrow(total_theta_ps), size = n_posterior_samples_plotting, replace = FALSE)
colors <- hue_pal()(n_posterior_samples_plotting) 

if (file.exists(p_val_path)) {
    p_val_dt <- fread(p_val_path)
    total_theta_ps <- fread(total_theta_p_val_path)
    total_unif_ps <- fread(total_unif_p_val_path)
    total_indep_ps <- fread(total_indep_p_val_path)
    p_theta_density_mat <- read.table(p_theta_density_mat_path)
    p_unif_density_mat <- read.table(p_unif_density_mat_path)
    p_indep_density_mat <- read.table(p_indep_density_mat_path)
} else {
    for (i in 1:n_sims) {
        print(glue('sim {i} of {n_sims}'))

        # simulate data from correct model (Beta(1,1) (aka Unif(0,1)) prior on theta)
        a_0 <- 1
        b_0 <- 1
        theta_0 <- rbeta(1, a_0, b_0)
        theta_0_values[i] <- theta_0

        y <- rbinom(n, 1, theta_0)
        y_values[i, ] <- y

        T_y <- sum(diff(y) != 0)

        # get analytic posterior and calculate posterior predictive T_y
        if (i %in% samples) {
            n_posterior_sims <- 10^6 # bump up posterior simulations for sim where I generate plots
            null_dist_100 <- build_hoeffD_null(100, 10^6) # update null distribution every time we plot a new curve
        } else {
            n_posterior_sims <- 50
            # null_dist_100 <- build_hoeffD_null(100, 10^5)
        }
        theta <- rbeta(n_posterior_sims, a_0 + sum(y), b_0 + n - sum(y) ) # posterior draws of theta
        u_theta <- theta # prior is uniform
        pp_T_y <- rep(0, n_posterior_sims) # simulate data using each posterior draw of theta and calculate T
        for (k in 1:n_posterior_sims) {
            pp_T_y[k] <- sum(diff( rbinom(n, 1, theta[k]) ) != 0)
        }

        # generate residual u values
        u_resid_mat <- matrix(0, length(theta), length(y))
        for (j in 1:nrow(u_resid_mat)) {
            for (k in 1:ncol(u_resid_mat)) {
                if (y[k] == 1) {
                    u_resid_mat[j,k] <- runif(1,1-theta[j],1)
                } else {
                    u_resid_mat[j,k] <- runif(1,0,1-theta[j])
                }
            }
        }

        # get p-values 
        theta_single_p_vals <- 2*pmin(u_theta, 1-u_theta)
        theta_p_values[i] <- theta_single_p_vals[1]
        # theta_p_values[i] <- 1-u_theta[1]
        switch_p_values[i] <- mean(pp_T_y > T_y)

        resid_p_values <- uniformity_test(u_resid_mat) # p-values for resid uniformity testing
        # test for independence of residual p-values
        indep_p_values <- resid_independence_test(u_resid_mat, null_dist_100)

        data_unif_p_values[i] <- resid_p_values[1]
        data_indep_p_values[i] <- indep_p_values[1]

        single_p_val_dt <- data.table(
                           'theta' = theta_single_p_vals,
                           'data' = resid_p_values,
                           'indep' = indep_p_values
                           )
        total_theta_ps[i,] <- 1-u_theta[1:ncol(total_theta_ps)]
        total_unif_ps[i,] <- resid_p_values[1:ncol(total_unif_ps)]
        total_indep_ps[i,] <- indep_p_values[1:ncol(total_indep_ps)]

        # perform density calculations for each set of p-values
        # TODO: could write function to do this one analytically 
        #       NOTE: bernoulli_f_p_theta_given_y is wrong and, I think, quite hard to calculate analytically
        #             trying empirical version instead
        # p_theta_density <- sapply(seq(0, 1, by = 0.01), FUN = bernoulli_f_p_theta_given_y, y = y)
        # p_theta_density_mat[i, ] <- p_theta_density

        theta_p_value_dt <- data.table('theta' = theta_single_p_vals)
        theta_dens_data <- ggplot_build(
                    ggplot(theta_p_value_dt, aes(x = theta)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
        theta_dens_data <- data.table(theta_dens_data[, c('xmin', 'density')])
        setnames(theta_dens_data, c('x', 'density'))
        theta_dens_data <- rbindlist( list( data.table('x' = 0, 'density' = 0), theta_dens_data, data.table('x' = 1, 'density' = 0) ) )
        p_theta_density_mat[i, ] <- theta_dens_data$density

        unif_p_value_dt <- data.table('unif' = resid_p_values)
        unif_dens_data <- ggplot_build(
                    ggplot(unif_p_value_dt, aes(x = unif)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
        unif_dens_data <- data.table(unif_dens_data[, c('xmin', 'density')])
        setnames(unif_dens_data, c('x', 'density'))
        unif_dens_data <- rbindlist( list( data.table('x' = 0, 'density' = 0), unif_dens_data, data.table('x' = 1, 'density' = 0) ) )
        p_unif_density_mat[i, ] <- unif_dens_data$density

        indep_p_value_dt <- data.table('indep' = indep_p_values)
        indep_dens_data <- ggplot_build(
                    ggplot(indep_p_value_dt, aes(x = indep)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
        indep_dens_data <- data.table(indep_dens_data[, c('xmin', 'density')])
        setnames(indep_dens_data, c('x', 'density'))
        indep_dens_data <- rbindlist( list( data.table('x' = 0, 'density' = 0), indep_dens_data, data.table('x' = 1, 'density' = 0) ) )
        p_indep_density_mat[i, ] <- indep_dens_data$density

        if (i == samples[1]) {
            bernoulli_p_value_densities_single_posterior(single_p_val_dt, y, a_0, b_0, prior_results_dir, plot_switch = FALSE)
        }
    }

    p_val_dt <- data.table('theta_0' = theta_0_values,
                           'switch' = switch_p_values, 
                           'theta' = theta_p_values,
                           'data' = data_unif_p_values,
                           'indep' = data_indep_p_values
                           )
    fwrite(p_val_dt, p_val_path)
    fwrite(total_theta_ps, total_theta_p_val_path)
    fwrite(total_unif_ps, total_unif_p_val_path)
    fwrite(total_indep_ps, total_indep_p_val_path)
    write.table(p_theta_density_mat, p_theta_density_mat_path, col.names = FALSE, row.names = FALSE)
    write.table(p_unif_density_mat, p_unif_density_mat_path, col.names = FALSE, row.names = FALSE)
    write.table(p_indep_density_mat, p_indep_density_mat_path, col.names = FALSE, row.names = FALSE)
}

# plot distribution of p-values
bernoulli_p_value_densities_sims(p_val_dt, prior_results_dir, plot_switch = TRUE)


###### histograms/densities across individual posteriors ########

### theta ###
p_theta_analytic_dt <- data.table(p_theta_density_mat)
colnames(p_theta_analytic_dt) <- as.character(as.numeric(str_replace_all(colnames(p_theta_analytic_dt), 'V', '')))
p_theta_analytic_dt[, sim_sample := as.character( 1:nrow(p_theta_analytic_dt) )]
p_theta_analytic_dt_long <- melt( p_theta_analytic_dt, variable.name = 'p_theta', value.name = 'density', id.vars = 'sim_sample' )
# p_theta_analytic_dt_long[, p_theta := pmax(0, as.numeric(p_theta)/100 - 0.02)] # empirical only 
p_theta_analytic_dt_long[, p_theta := pmax(0, as.numeric(p_theta)/100 - 0.02)]
p_theta_analytic_dt_long[, standardized_density := density / max(density), by = sim_sample]
p_theta_analytic_dt_long <- p_theta_analytic_dt_long[! (p_theta %in% c(0,1))] # TODO: trying this


avg_densities_analytic <- data.table(colMeans(p_theta_density_mat))
avg_densities_analytic[, density := V1]
# avg_densities_analytic[, p_theta := c(0, seq(0, 1, 0.01))]
avg_densities_analytic[, p_theta := c(0, seq(0, 1, 0.01))]
avg_densities_analytic <- avg_densities_analytic[! (p_theta %in% c(0,1))] # remove points at 0,1

theta_indiv_posterior_density_plot <-  ggplot(p_theta_analytic_dt_long[sim_sample %in% samples], aes(x = p_theta, y = standardized_density, colour = sim_sample, group = sim_sample)) + 
                                        geom_line(position = 'identity', alpha = 0.5, size = 1) +
                                        scale_color_manual(values = colors) +
                                        geom_line(data = avg_densities_analytic, 
                                                  aes(x = p_theta, y = density, group = 1), position = 'identity', color = 'black', size = 1) +
                                        coord_cartesian(xlim = c(0, 1), ylim = c(0, NA)) +
                                        labs(x = TeX('$p_{\\theta}$'), y = 'standardized density') +
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
                                                        
ggsave(file.path(prior_results_dir, 'theta_individual_posterior_density_plot.png'), theta_indiv_posterior_density_plot, width = 12, height = 8, dpi = 300)

### unif ###
p_unif_analytic_dt <- data.table(p_unif_density_mat)
colnames(p_unif_analytic_dt) <- as.character(as.numeric(str_replace_all(colnames(p_unif_analytic_dt), 'V', '')))
p_unif_analytic_dt[, sim_sample := as.character( 1:nrow(p_unif_analytic_dt) )]
p_unif_analytic_dt_long <- melt( p_unif_analytic_dt, variable.name = 'p_unif', value.name = 'density', id.vars = 'sim_sample' )
p_unif_analytic_dt_long[, p_unif := pmax(0, as.numeric(p_unif)/100 - 0.02)] # empirical only 
# p_unif_analytic_dt_long[, p_unif := as.numeric(p_unif)/100 - 0.01]
p_unif_analytic_dt_long[, standardized_density := density / max(density), by = sim_sample]
p_unif_analytic_dt_long <- p_unif_analytic_dt_long[! (p_unif %in% c(0,1))] # TODO: trying this

avg_densities_analytic <- data.table(colMeans(p_unif_density_mat))
avg_densities_analytic[, density := V1]
# avg_densities_analytic[, standardized_density := V1 / max(V1)]
avg_densities_analytic[, p_unif := c(0, seq(0, 1, 0.01))]
# avg_densities_analytic[, p_unif := seq(0, 1, 0.01)]
avg_densities_analytic <- avg_densities_analytic[! (p_unif %in% c(0,1))] # remove points at 0,1

unif_indiv_posterior_density_plot <-  ggplot(p_unif_analytic_dt_long[sim_sample %in% samples], aes(x = p_unif, y = standardized_density, colour = sim_sample, group = sim_sample)) + 
                                        geom_line(position = 'identity', alpha = 0.5, size = 1) +
                                        scale_color_manual(values = colors) +
                                        geom_line(data = avg_densities_analytic, 
                                                  aes(x = p_unif, y = density, group = 1), position = 'identity', color = 'black', size = 1) +
                                        coord_cartesian(xlim = c(0, 1), ylim = c(0, NA)) +
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
                                                        
ggsave(file.path(prior_results_dir, 'unif_individual_posterior_density_plot.png'), unif_indiv_posterior_density_plot, width = 12, height = 8, dpi = 300)

### indep ###
p_indep_analytic_dt <- data.table(p_indep_density_mat)
colnames(p_indep_analytic_dt) <- as.character(as.numeric(str_replace_all(colnames(p_indep_analytic_dt), 'V', '')))
p_indep_analytic_dt[, sim_sample := as.character( 1:nrow(p_indep_analytic_dt) )]
p_indep_analytic_dt_long <- melt( p_indep_analytic_dt, variable.name = 'p_indep', value.name = 'density', id.vars = 'sim_sample' )
p_indep_analytic_dt_long[, p_indep := pmax(0, as.numeric(p_indep)/100 - 0.02)] # empirical only 
# p_indep_analytic_dt_long[, p_indep := as.numeric(p_indep)/100 - 0.01]
p_indep_analytic_dt_long[, standardized_density := density / max(density), by = sim_sample]
p_indep_analytic_dt_long <- p_indep_analytic_dt_long[! (p_indep %in% c(0,1))] # TODO: trying this

avg_densities_analytic <- data.table(colMeans(p_indep_density_mat))
avg_densities_analytic[, density := V1]
# avg_densities_analytic[, standardized_density := V1 / max(V1)]
avg_densities_analytic[, p_indep := c(0, seq(0, 1, 0.01))]
# avg_densities_analytic[, p_indep := seq(0, 1, 0.01)]
avg_densities_analytic <- avg_densities_analytic[! (p_indep %in% c(0,1))] # remove points at 0,1

indep_indiv_posterior_density_plot <-  ggplot(p_indep_analytic_dt_long[sim_sample %in% samples], aes(x = p_indep, y = standardized_density, colour = sim_sample, group = sim_sample)) + 
                                        geom_line(position = 'identity', alpha = 0.5, size = 1) +
                                        scale_color_manual(values = colors) +
                                        geom_line(data = avg_densities_analytic, 
                                                  aes(x = p_indep, y = density, group = 1), position = 'identity', color = 'black', size = 1) +
                                        coord_cartesian(xlim = c(0, 1), ylim = c(0, NA)) +
                                        labs(x = expression('p'['data,indep']), y = 'standardized density') +
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
                                                        
ggsave(file.path(prior_results_dir, 'indep_individual_posterior_density_plot.png'), indep_indiv_posterior_density_plot, width = 12, height = 8, dpi = 300)
