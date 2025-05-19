library(pacman)
pacman::p_load('rstan', 'data.table', 'ggplot2', 'forcats', 
'energy', 
'gridExtra', 'EnvStats', 'stringr', 'latex2exp',
'Hmisc', 'readr', 'truncnorm', 'stringr', 'purrr', 'patchwork')

setwd('/Users/christiancovington/research/uniform_parameterization_checks/code/AR_model')

source('generate_dap_samples.R')
source('generate_residuals.R')
source('plot_creation.R')
source('../utils/generate_plots.R')
source('../utils/testing.R')
# source('power_plot_creation.R')

# parent directories
data_dir <- file.path('..', '..', 'data', 'AR')
results_dir <- file.path('..', '..', 'results', 'AR')

set.seed(1)

# hyperparameters
n_sims <- 100000
n_samples <- 1
n <- 500

qualitative_priors <- c('correct', 
                        'phi_very_narrow',
                        'sigma_very_narrow'
                        )

# establish sets of priors whose behavior we are interested in
# qualitative_priors <- c('very_narrow', 'narrow', 'correct', 'wide', 'very_wide',
#                         'phi_very_narrow', 'phi_narrow', 'phi_wide', 'phi_very_wide',
#                         'sigma_very_narrow', 'sigma_narrow', 'sigma_wide', 'sigma_very_wide'
#                     )



# phi_centers <- seq(-0.3, 0.3, by = 0.15) # truth is 0
# phi_mults <- c(0.1, 0.5, 0.5, 0.75, 1) # truth is 0.5
# sigma_centers <- seq(1, 2, by = 0.5) # truth is 1.5
# sigma_mults <- c(0.1, 0.5, 1, 1.5, 2) # truth is 1

# # generate all quantitative priors
# quantitative_priors <- data.table(expand.grid(phi_centers, phi_mults, 
#                                               sigma_centers, sigma_mults)
#                                  )
# # keep only the priors for which either phi or sigma is correct or either both centers or both mults are correct
# # CC NOTE: this is for plotting purposes because it's difficult to 
# #          show more than two dimensions
# quantitative_priors <- quantitative_priors[(Var1 == 0 & Var2 == 0.5) | (Var3 == 1.25 & Var4 == 1) | 
#                                            (Var1 == 0 & Var3 == 1.5) | (Var2 == 0.5 & Var4 == 1)
#                                           ]

hoeffd_null_lag_1 <- build_lag_1_hoeffD_null(n, 10000)
hoeffd_null_lag_2 <- build_lag_2_hoeffD_null(n, 10000)
hoeffd_index_null <- build_hoeffD_index_null(n, 10000)

#############################################
############# null testing ##################
#############################################
# big_n_sims <- n_sims*10
# lag_1_ps <- rep(0, big_n_sims)
# lag_2_ps <- rep(0, big_n_sims)
# for (i in 1:big_n_sims) {
#     # print(glue::glue('run {i} of {big_n_sims}'))
#     y <- runif(n)
#     lag_1_D <- hoeffd(y[1:(n-1)], y[2:n])$D[1,2]
#     lag_2_D <- hoeffd(y[1:(n-2)], y[3:n])$D[1,2]

#     lag_1_ps[i] <- mean(hoeffd_null_lag_1 < lag_1_D)
#     lag_2_ps[i] <- mean(hoeffd_null_lag_2 < lag_2_D)
# }
#################################################
############# end null testing ##################
#################################################

DGPs <- c('AR1', 'AR1_heteroskedastic_V1', 'AR1_heteroskedastic_V2', 'AR2')

# run model over qualitative set of priors
for (prior in qualitative_priors) {
    model.stan <- stan_model(file = paste0(prior, '.stan'))

    if (prior == 'correct') {
        relevant_DGPs <- c('AR1', 'AR1_heteroskedastic_V1', 'AR1_heteroskedastic_V2', 'AR2')
    } else if (prior == 'phi_very_narrow') {
        relevant_DGPs <- c('AR1')
    } else if (prior == 'sigma_very_narrow') {
        relevant_DGPs <- c('AR1')
    }

    for (DGP in relevant_DGPs) {
        print(glue::glue('generating data/plots for {prior} prior and {DGP} DGP'))
        # directory setup
        model_results_dir <- file.path(results_dir, glue::glue('n = {n}'), glue::glue('n_sims = {n_sims}'), glue::glue('n_posterior = {n_samples}'), prior, DGP)
        dir.create(model_results_dir, showWarnings = FALSE, recursive = TRUE)

        #######################
        ### sample from DAP ###
        #######################
        dap_samples_file <- file.path(model_results_dir, 'dap_samples.rds')
        y_file <- file.path(model_results_dir, 'y.rds')
        if (file.exists(dap_samples_file)) {
            # res_dt <- fread(dap_samples_file)
            # y <- fread(y_file)
            dap_samples_list <- read_rds(dap_samples_file)
            y_list <- read_rds(y_file)
        } else {
            dap_samples_list <- vector('list', length = n_sims)
            y_list <- vector('list', length = n_sims)
            for (i in 1:n_sims) {
                print(glue::glue('sim {i} of {n_sims}'))
                res <- generate_dap_samples(model.stan, n, DGP, n_samples, s = 1)
                
                # map samples back through prior
                dap_samples <- res['dap_samples']$dap_samples
                if (prior == 'correct') {
                    dap_samples[, u1 := ptruncnorm(phi, a = -0.5, b = 0.5, mean = 0, sd = 0.4)]
                    dap_samples[, u2 := ptruncnorm(sigma, a = 1, b = 2, mean = 1.5, sd = 0.4)]
                } else if (prior == 'phi_very_narrow') {
                    dap_samples[, u1 := ptruncnorm(phi, a = -0.5, b = 0.5, mean = 0, sd = 0.1)]
                    dap_samples[, u2 := ptruncnorm(sigma, a = 1, b = 2, mean = 1.5, sd = 0.4)]
                } else if (prior == 'sigma_very_narrow') {
                    dap_samples[, u1 := ptruncnorm(phi, a = -0.5, b = 0.5, mean = 0, sd = 0.4)]
                    dap_samples[, u2 := ptruncnorm(sigma, a = 1, b = 2, mean = 1.5, sd = 0.1)]
                }

                dap_samples_list[[i]] <- dap_samples
                y_list[[i]] <- res['y']$y
            }
            write_rds(dap_samples_list, dap_samples_file)
            write_rds(y_list, y_file)
        }

        ####################################
        ### perform tests on DAP samples ###
        ####################################
        p_val_file <- file.path(model_results_dir, 'p_val_dt.csv')

        if (file.exists(p_val_file)) {
            p_val_dt <- fread(p_val_file)
        } else {
            p_val_dt_list <- vector('list', length = n_sims)

            for (i in 1:n_sims) {
                print(glue::glue('perform tests for sim {i} of {n_sims}'))
                res_dt <- dap_samples_list[[i]]
                y <- y_list[[i]]

                # test uniformity of parameter u values 
                phi_p_values <- u_to_p_converter(res_dt$u1)
                sigma_p_values <- u_to_p_converter(res_dt$u2)

                phi_p_value <- Cauchy_p_merger(phi_p_values)
                sigma_p_value <- Cauchy_p_merger(sigma_p_values)

                # generate data-u values and perform uniformity test
                u_resid_mat <- resid_u_generation(res_dt, y)
                data_unif_p_values <- uniformity_test(u_resid_mat)
                data_unif_p_value <- Cauchy_p_merger(data_unif_p_values)
                data_u_unif_p_val_string <- signif(data_unif_p_value, 2)

                # test first/second order independence of data u-values
                first_order_indep_p_values <- resid_independence_test(u_resid_mat, hoeffd_null_lag_1)
                second_order_indep_p_values <- resid_second_order_independence_test(u_resid_mat, hoeffd_null_lag_2)
                first_order_indep_p_value <- Cauchy_p_merger(first_order_indep_p_values)
                second_order_indep_p_value <- Cauchy_p_merger(second_order_indep_p_values)

                # test independence of data u values and index
                index_indep_p_values <- resid_index_independence_test(u_resid_mat, hoeffd_index_null)
                index_indep_p_value <- Cauchy_p_merger(index_indep_p_values)

                # save test results
                p_val_dt_list[[i]] <- data.table('phi_unif' = phi_p_value,
                                                'sigma_unif' = sigma_p_value,
                                                'data_unif' = data_unif_p_value,
                                                'index_indep' = index_indep_p_value,
                                                'lag1_indep' = first_order_indep_p_value,
                                                'lag2_indep' = second_order_indep_p_value
                                                )
                # p_val_dt_list[[i]] <- data.table('phi_unif' = phi_p_values[1],
                #                                 'sigma_unif' = sigma_p_values[1],
                #                                 'data_unif' = data_unif_p_values[1],
                #                                 'index_indep' = index_indep_p_values[1],
                #                                 'lag1_indep' = first_order_indep_p_values[1],
                #                                 'lag2_indep' = second_order_indep_p_values[1],
                #                                 'Cauchy_phi_unif' = phi_p_value,
                #                                 'Cauchy_sigma_unif' = sigma_p_value,
                #                                 'Cauchy_data_unif' = data_unif_p_value,
                #                                 'Cauchy_index_indep' = index_indep_p_value,
                #                                 'Cauchy_lag1_indep' = first_order_indep_p_value,
                #                                 'Cauchy_lag2_indep' = second_order_indep_p_value
                #                                 )
                # produce plots for "single run" of interest
                if (i == 1) {
                    # plot data u-value ecdfs
                    residual_ecdf_plot(u_resid_mat, min(nrow(u_resid_mat), 50), data_u_unif_p_val_string, model_results_dir)
                    residual_tilted_ecdf_plot(u_resid_mat, min(nrow(u_resid_mat), 50), data_u_unif_p_val_string, model_results_dir, 0.35)
                } 
            }

            # combine results into single table 
            p_val_dt <- rbindlist(p_val_dt_list)
            fwrite(p_val_dt, p_val_file)
        }


        # plot distribution of p-values for each test 
        test_names <- c(
                        TeX('$p_{\\phi}$'),
                        TeX('$p_{\\sigma}$'),
                        quote('p'['data,unif']),
                        quote('p'['data,index']),
                        quote('p'['data,lag1']),
                        quote('p'['data,lag2'])
                       )
        file_names <- c('phi_unif', 'sigma_unif', 'data_unif',
                        'index_indep', 'data_lag1_indep', 'data_lag2_indep'
                        )
        AR_p_val_dist(p_val_dt, test_names, file_names, model_results_dir)
    }
}

#############################################################
## make omnibus plot for 5 different scenarios we consider ##
#############################################################

omnibus_p_val_dt_list <- vector('list', length = 5)
priors  <- c('correct', 'phi_very_narrow', 'sigma_very_narrow', 'correct', 'correct')
DGPs <- c('AR1', 'AR1', 'AR1', 'AR1_heteroskedastic_V1', 'AR2')
scenarios <- c(glue::glue('scenario {1:5}'))
plotmath_scenarios <- c( quote('scenario 1'), quote('scenario 2'), quote('scenario 3'), quote('scenario 4'), quote('scenario 5') )
for (i in 1:5) {
    model_results_dir <- file.path(results_dir, glue::glue('n = {n}'), glue::glue('n_sims = {n_sims}'), glue::glue('n_posterior = {n_samples}'), priors[i], DGPs[i])
    p_val_dt <- fread( file.path(model_results_dir, 'p_val_dt.csv') )
    p_val_dt[, scenario := scenarios[i]]
    omnibus_p_val_dt_list[[i]] <- p_val_dt
}
omnibus_p_val_dt <- rbindlist(omnibus_p_val_dt_list)
omnibus_p_val_dt_long <- melt(omnibus_p_val_dt, id.vars = 'scenario', variable.name = 'test', value.name = 'p-value')

omnibus_p_val_dt[, scenario := factor(scenario, labels = plotmath_scenarios)]
omnibus_p_val_dt_long[, scenario := factor(scenario, labels = plotmath_scenarios)]
omnibus_p_val_dt_long[, test := factor(test, labels = test_names)]

plot_list <- vector('list', 30)
y_label_expression <- quote('CDF of ')
i <- 1
for (t_var in test_names) {
    for (s_var in scenarios) {
        if (i %in% 26:30) {
            x_lab <- glue::glue('{s_var}: p-value')
        } else {
            x_lab <- ''
        }

        if (i %% 5 == 1) {
            y_lab <- substitute(a * b, list(a = y_label_expression, b = t_var))
        } else {
            y_lab <- ''
        }
        t_num <- which(test_names == t_var)
        sub_dt <- omnibus_p_val_dt[scenario == s_var, t_num, with = FALSE]
        colnames(sub_dt) <- c('p-value')

        plot_list[[i]] <- ggplot(sub_dt, aes(x = `p-value`)) + 
                            stat_ecdf(position = 'identity', linewidth = 1, alpha = 1, colour = 'blue', show.legend = FALSE) + 
                            geom_abline(intercept = 0, slope = 1, colour = 'black') + 
                            coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
                            labs(x = x_lab, y = y_lab) +
                            theme(
                                legend.position = 'none',
                                strip.text = element_text(size = 12),
                                axis.title.x = element_text(size = 12),
                                axis.text.x = element_blank(),
                                axis.title.y = element_text(size = 12),
                                axis.text.y = element_blank(),
                                panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                panel.background = element_rect(fill = "white") # White background
                                )
        i <- i + 1
    }
}

omnibus_ecdf_plot <- wrap_plots(plot_list, nrow = 6) # each row is a test, each column is a scenario
ggsave( file.path(results_dir, glue::glue('n = {n}'), glue::glue('n_sims = {n_sims}'), glue::glue('n_posterior = {n_samples}'), 'omnibus_ecdf.png'), height = 12, width = 10, dpi = 300 )





##############
## NOTE: code below makes first version, but I can't figure out how to give each 
##       plot its own label 
##############
# omnibus_ecdf_plot <- ggplot(omnibus_p_val_dt_long, aes(x = `p-value`)) + 
#                         stat_ecdf(position = 'identity', linewidth = 1, alpha = 1, colour = 'blue', show.legend = FALSE) + 
#                         geom_abline(intercept = 0, slope = 1, colour = 'black') + 
#                         coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
#                         labs(x = 'p-value', y = 'empirical cdf') +
#                         theme(
#                             legend.position = 'none',
#                             strip.text = element_text(size = 12),
#                             axis.title.x = element_text(size = 12),
#                             axis.text.x = element_blank(),
#                             axis.title.y = element_text(size = 12),
#                             axis.text.y = element_blank(),
#                             panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
#                             panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
#                             panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
#                             panel.background = element_rect(fill = "white") # White background
#                             ) +
#                         facet_grid(test ~ scenario, labeller = label_parsed, scales = 'free')

# ggsave( file.path(results_dir, glue::glue('n = {n}'), glue::glue('n_sims = {n_sims}'), glue::glue('n_posterior = {n_samples}'), 'omnibus_ecdf.png'), height = 12, width = 12, dpi = 300 )

###################################################
### calculate power over different alternatives ###
###################################################
set.seed(1)

n <- 200
n_sims <- 10^4
n_posteriors <- c(1)

hoeffd_null_lag_1 <- build_lag_1_hoeffD_null(n, 10000)
hoeffd_null_lag_2 <- build_lag_2_hoeffD_null(n, 10000)
hoeffd_index_null <- build_hoeffD_index_null(n, 10000)

alternatives <- c('phi_additive_narrow_power',
                  'sigma_additive_narrow_power',
                  'phi_multiplicative_narrow_power',
                  'sigma_multiplicative_narrow_power',
                  'heteroskedastic_V1_power',
                  'heteroskedastic_V2_power',
                  'AR_2_power'
                 )

s_vec <- seq(0, 1, 0.05)

# generate p-values for each combination of alternative/s
for (alternative in alternatives) {
    for (s in s_vec) {
        for (n_posterior in n_posteriors) {
            # change s from 0/1 to 0.01/0.99 if s in {0,1} would cause prior to be point mass
            if (alternative %in% c('phi_multiplicative_narrow_power', 'sigma_multiplicative_narrow_power')) {
                if (s == 0) {
                    s <- 0.01
                } else if (s == 1) {
                    s <- 0.99
                }
            }

            print(glue::glue('alternative: {alternative}, s = {s}'))
            # directory setup
            model_results_dir <- file.path(results_dir, 'power_calculations', glue::glue('n = {n}'), glue::glue('n_sims = {n_sims}'), glue::glue('n_posterior = {n_posterior}'), alternative)
            dir.create(model_results_dir, showWarnings = FALSE, recursive = TRUE)

            # set prior/model
            if (alternative %in% c('phi_additive_narrow_power', 'sigma_additive_narrow_power',
                                'phi_multiplicative_narrow_power', 'sigma_multiplicative_narrow_power')) {
                model.stan <- stan_model(file = paste0(alternative, '.stan'))
            } else {
                model.stan <- stan_model(file = 'correct.stan')
            }

            # set DGP
            if (alternative %in% c('phi_additive_narrow_power', 'sigma_additive_narrow_power',
                                'phi_multiplicative_narrow_power', 'sigma_multiplicative_narrow_power')) {
                DGP <- 'AR1'
            } else if (alternative == 'heteroskedastic_V1_power') {
                DGP <- 'AR1_heteroskedastic_V1'
            } else if (alternative == 'heteroskedastic_V2_power') {
                DGP <- 'AR1_heteroskedastic_V2'
            } else if (alternative == 'AR_2_power') {
                DGP <- 'AR2'
            } else {
                except('Alternative must be correctly specified')
            }


            #######################
            ### sample from DAP ###
            #######################
            dap_samples_file <- file.path(model_results_dir, glue::glue('dap_samples_s_{s}.rds'))
            y_file <- file.path(model_results_dir, glue::glue('y_s_{s}.rds'))
            if (file.exists(dap_samples_file)) {
                dap_samples_list <- read_rds(dap_samples_file)
                y_list <- read_rds(y_file)
            } else {
                dap_samples_list <- vector('list', length = n_sims)
                y_list <- vector('list', length = n_sims)
                for (i in 1:n_sims) {
                    if (i %% 100 == 0) {
                        print(glue::glue('sim {i} of {n_sims}'))
                    }
                    res <- generate_dap_samples(model.stan, n, DGP, n_samples, s)

                    # if alternative is misspecified prior, manually map samples back through prior
                    # NOTE: this is done manually (rather than automatically in Stan) because there is, afaik, 
                    #       no truncated normal cdf/quantile function in Stan
                    # TODO: make sure that these match whatever parametrization we end up using in the 
                    #       Stan files
                    dap_samples <- res['dap_samples']$dap_samples
                    if (alternative == 'phi_additive_narrow_power') {
                        dap_samples[, u1 := ptruncnorm(phi, a = -0.5, b = 0.5, mean = 0-3*s, sd = 0.4)]
                        dap_samples[, u2 := ptruncnorm(sigma, a = 1, b = 2, mean = 1.5-3*s, sd = 0.4)]
                    } else if (alternative == 'sigma_additive_narrow_power') {
                        dap_samples[, u1 := ptruncnorm(phi, a = -0.5, b = 0.5, mean = 0-3*s, sd = 0.4)]
                        dap_samples[, u2 := ptruncnorm(sigma, a = 1, b = 2, mean = 1.5-3*s, sd = 0.4)]
                    } else if (alternative == 'phi_multiplicative_narrow_power') {
                        dap_samples[, u1 := ptruncnorm(phi, a = -0.5, b = 0.5, mean = 0, sd = (1-s)*0.4)]
                        dap_samples[, u2 := ptruncnorm(sigma, a = 1, b = 2, mean = 1.5, sd = 0.4)]
                    } else if (alternative == 'sigma_multiplicative_narrow_power') {
                        dap_samples[, u1 := ptruncnorm(phi, a = -0.5, b = 0.5, mean = 0, sd = 0.4)]
                        dap_samples[, u2 := ptruncnorm(sigma, a = 1, b = 2, mean = 1.5, sd = (1-s)*0.4)]
                    }

                    dap_samples_list[[i]] <- dap_samples
                    y_list[[i]] <- res['y']$y
                }
                write_rds(dap_samples_list, dap_samples_file)
                write_rds(y_list, y_file)
            }

            ####################################
            ### perform tests on DAP samples ###
            ####################################
            p_val_file <- file.path(model_results_dir, glue::glue('p_val_dt_s_{s}.csv'))

            if (file.exists(p_val_file)) {
                p_val_dt <- fread(p_val_file) # NOTE: this does nothing at the moment
            } else {
                p_vals <- rep(0, n_sims)

                for (i in 1:n_sims) {
                    # print(glue::glue('perform tests for sim {i} of {n_sims}'))
                    res_dt <- dap_samples_list[[i]]
                    y <- y_list[[i]]

                    if (alternative %in% c('phi_additive_narrow_power', 'phi_multiplicative_narrow_power')) { 
                        p_vals[i] <- Cauchy_p_merger(u_to_p_converter(res_dt$u1))
                    } else if (alternative %in% c('sigma_additive_narrow_power', 'sigma_multiplicative_narrow_power')) {
                        p_vals[i] <- Cauchy_p_merger(u_to_p_converter(res_dt$u2))
                    } else if (alternative %in% c('heteroskedastic_V1_power', 'heteroskedastic_V2_power')) {
                        u_resid_mat <- resid_u_generation(res_dt, y)
                        p_vals[i] <- Cauchy_p_merger(resid_index_independence_test(u_resid_mat, hoeffd_index_null))
                    } else if (alternative == 'AR_2_power') {
                        u_resid_mat <- resid_u_generation(res_dt, y)
                        p_vals[i] <- Cauchy_p_merger(resid_second_order_independence_test(u_resid_mat, hoeffd_null_lag_2))
                    }
                }
                write.csv(p_vals, p_val_file)
            }
        }
    }
}

###########################################
### create combined p-value dt and plot ###
###########################################

test_names <- c(expression(Pr(p[phi] <= 0.05)),
                expression(Pr(p[phi] <= 0.05)),
                expression(Pr(p[sigma] <= 0.05)),
                expression(Pr(p[sigma] <= 0.05)),
                expression(Pr(p['data, index'] <= 0.05)),
                expression(Pr(p['data, index'] <= 0.05)),
                expression(Pr(p['data, lag2'] <= 0.05))
                )

for (n_posterior in n_posteriors) {
    for (i in seq_along(alternatives)) {
        alternative <- alternatives[i]
        test_name <- test_names[i]

        
        # load p-values
        model_results_dir <- file.path(results_dir, 'power_calculations', glue::glue('n = {n}'), glue::glue('n_sims = {n_sims}'), glue::glue('n_posterior = {n_posterior}'), alternative)
        p_val_dt_list <- vector('list', length = length(s_vec))

        # get power for each value of s
        for (j in seq_along(s_vec)) {
            s <- s_vec[j]
            if (alternative %in% c('phi_multiplicative_narrow_power', 'sigma_multiplicative_narrow_power')) {
                if (s == 0) {
                    s <- 0.01
                } else if (s == 1) {
                    s <- 0.99
                }
            }
            p_vals <- fread( file.path(model_results_dir, glue::glue('p_val_dt_s_{s}.csv')) )
            p_val_dt_list[[j]] <- data.table('s' = s, 'p_val' = p_vals$x)
        }
        p_val_dt <- rbindlist(p_val_dt_list)
        p_val_dt[, power := mean(p_val <= 0.05), by = s]
        p_val_dt <- unique(p_val_dt[, c('s', 'power')])

        # plot power as function of s 
        plot <- ggplot(p_val_dt, aes(x = s, y = power)) + 
                    geom_line(colour = 'blue') + 
                    labs(x = 's', y = test_name) + 
                    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
                    theme(
                        axis.title.x = element_text(size = 24),
                        axis.text.x = element_text(size = 24),
                        axis.title.y = element_text(size = 24),
                        axis.text.y = element_text(size = 24),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.5), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                    )
        ggsave(file.path(model_results_dir, 'power_plot.png'), plot, width = 12, height = 8, dpi = 300)
    }

}
