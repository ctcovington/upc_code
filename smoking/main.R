library(pacman)
pacman::p_load('rstan', 'data.table', 'ggplot2', 'forcats', 'energy', 
               'gridExtra', 'EnvStats', 'stringr', 'Hmisc', 'lmtest', 'dplyr',
               'parallel', 'rjags', 'boot', 'psych', 'broom', 'latex2exp', 'glue', 'MetBrewer')
               
source('../utils/generate_plots.R')
source('../utils/testing.R')
source('../utils/sim_data.R')

setwd('~/research/uniform_parameterization_checks/code/smoking')

# parent directories
data_dir <- file.path('..', '..', 'data', 'smoking')
results_dir <- file.path('..', '..', 'results', 'smoking')
dir.create(results_dir, showWarnings = FALSE)

# perform entire experiment for both simulated and real data

# data_sources <- c('sim_data')
data_sources <- c('real_data')

wave_colors <- met.brewer('Juarez')
covar_colors <- met.brewer('Egypt')

for (data_source in data_sources) {
    # load data
    if (data_source == 'real_data') {
        smoking_no_scale <- fread(file.path(data_dir, 'smoking.txt'))
        smoking_no_scale[, wave := paste0('wave', wave)]
        smoking_no_scale[, sex := ifelse(`sex(1=F)` == 1, 'female', 'male')]
        smoking_no_scale[, parsmk := ifelse(parsmk == 1, 'parsmk = 1', 'parsmk = 0')]

        smoking <- fread(file.path(data_dir, 'smoking.txt'))
    } else if (data_source == 'sim_data') {
        smoking <- fread(file.path(data_dir, 'smoking.txt'))
        smoking[, smkgreg := create_random_smkreg(smoking)]
    }

    # Calculate marginal means for smkreg by wave and sex
    data_sex_distribution <- smoking %>%
        group_by(wave, `sex(1=F)`) %>%
        summarise(
            p_smkreg = mean(smkreg),
            .groups = 'drop'
        ) %>%
        mutate(variable = 'sex')

    # Calculate marginal means for smkreg by wave and parental smoking
    data_parsmk_distribution <- smoking %>%
        group_by(wave, parsmk) %>%
        summarise(
            p_smkreg = mean(smkreg),
            .groups = 'drop'
        ) %>%
        mutate(variable = 'parsmk') 

    # Combine the two datasets for plotting
    data_distribution_dt <- bind_rows(data_sex_distribution, data_parsmk_distribution)
    data_distribution_dt <- data_distribution_dt %>%
                                mutate(
                                    variable_value = ifelse(!is.na(`sex(1=F)`), `sex(1=F)`, parsmk)
                                )
    data_distribution_dt <- data.table(data_distribution_dt)
    data_distribution_dt[, 'grouping' := '']
    data_distribution_dt[variable == 'sex' & variable_value == 0, grouping := 'male']
    data_distribution_dt[variable == 'sex' & variable_value == 1, grouping := 'female']
    data_distribution_dt[variable == 'parsmk' & variable_value == 0, grouping := 'parsmk = 0']
    data_distribution_dt[variable == 'parsmk' & variable_value == 1, grouping := 'parsmk = 1']

    # Plot distribution of average smkreg by covariates/wave
    data_distribution_plot <- ggplot(data_distribution_dt, aes(x = wave, 
                                                            y = p_smkreg, 
                                                            colour = grouping
                                                            )) + 
        geom_line(linewidth = 2) + 
        scale_colour_manual(values = covar_colors) +
        scale_x_continuous(breaks = 1:6) +
        labs(
            x = 'wave',
            y = 'proportion smoking regularly',
            colour = ''
        ) + 
        theme(
                        legend.position = 'top',
                        legend.key.width = unit(3, 'line'),
                        legend.title = element_text(size = 24),
                        legend.text = element_text(size = 24),
                        axis.title.x = element_text(size = 24),
                        axis.text.x = element_text(size = 24),
                        # axis.ticks.x = element_text(size = 24),
                        axis.title.y = element_text(size = 24),
                        axis.text.y = element_text(size = 24),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )
    ggsave(file.path(results_dir, 'data_distribution.png'), data_distribution_plot, width = 12, height = 8, dpi = 300)
    
    # scale covariates to be mean 0, sd 1
    smoking[, wave := scale(wave)]
    smoking[, parsmk := scale(parsmk)]
    smoking[, `sex(1=F)` := scale(`sex(1=F)`)]

    parsmk_0_map <- smoking[, min(parsmk)]
    sex_0_map <- smoking[, min(`sex(1=F)`)]
    parsmk_1_map <- smoking[, max(parsmk)]
    sex_1_map <- smoking[, max(`sex(1=F)`)]

    # get covariate vectors matched w individual (not observations)
    parsmk_vec <- as.numeric(unlist(smoking[, .SD[1], newid][, 'parsmk']))
    sex_vec <- as.numeric(unlist(smoking[, .SD[1], newid][, 'sex(1=F)']))

    parsmk_vec_named <- parsmk_vec 
    parsmk_vec_named[ parsmk_vec_named == parsmk_0_map ] <- 'parsmk = 0'
    parsmk_vec_named[ parsmk_vec_named == parsmk_1_map ] <- 'parsmk = 1'

    sex_vec_named <- sex_vec 
    sex_vec_named[ sex_vec_named == sex_0_map ] <- 'male'
    sex_vec_named[ sex_vec_named == sex_1_map ] <- 'female'

    covariate_map_dt <- data.table('newid' = unique(smoking$newid),
                                'sex' = sex_vec,
                                'parsmk' = parsmk_vec
                                )

    fwrite( covariate_map_dt, file.path(results_dir, 'covariate_map.csv') )
    # rm(covariate_map_dt)

    if (data_source == 'real_data') {
        models <- c(
                # 'model_fixed_intercept',
                # 'model_intercept_only',
                'model_0', 
                'model_1',
                'model_2'
                )
    } else if (data_source == 'sim_data') {
        models <- c('model_0', 
                'model_1',
                'model_1.5',
                'model_2')
    }

    for (model in models) {
        print(glue::glue('running procedure for {model}'))

        # set up directory structure
        model_results_dir <- file.path(results_dir, data_source, model)
        dir.create(model_results_dir, showWarnings = FALSE, recursive = TRUE)

        # fit model
        stan_data <- list(n = nrow(smoking),
                            m = length(unique(smoking$newid)),
                            newid = smoking$newid,
                            wave = smoking$wave,
                            sex = smoking$`sex(1=F)`,
                            parsmk = smoking$parsmk,
                            smkreg = smoking$smkreg
                            )

        model.fit <- jags.model(file = glue::glue('{model}_JAGS.txt'),
                                data = stan_data,
                                n.chains = 1
                            )
        
        dap_samples_file <- file.path(model_results_dir, 'dap_samples.rds')

        if (file.exists(dap_samples_file)) {
            dap_samples <- readRDS(dap_samples_file)
        } else {
            # perform burn-in
            # update(model.fit, n.iter = 10000)
            update(model.fit, n.iter = 5000)

            # get posterior samples 
            dap_samples <- coda.samples(model.fit, n.iter = 1000*100,
                                            variable.names = c('mu', 'alpha', 'Z', 'g', 'q', 'p', 'mu', 'tau'),
                                            thin = 100
                                            )[[1]]
            saveRDS(dap_samples, dap_samples_file)
        }

        # separate posterior samples
        # posterior_beta <- dap_samples[, startsWith(colnames(dap_samples), 'beta')]
        posterior_mu <- dap_samples[, startsWith(colnames(dap_samples), 'mu')]
        posterior_alpha <- dap_samples[, startsWith(colnames(dap_samples), 'alpha')]
        posterior_p <- dap_samples[, startsWith(colnames(dap_samples), 'p')]
        if (model %in% c('model_2')) {
            posterior_Z <- dap_samples[, startsWith(colnames(dap_samples), 'Z')]
            # posterior_g <- dap_samples[, 'g']
            # posterior_q <- dap_samples[, 'q']
        }

        # map posterior samples back to u based on priors
        # beta_us <- pnorm(posterior_beta, mean = 0, sd = 5)
        if (model != 'model_2') {
            alpha_us <- matrix(0, nrow(posterior_alpha), ncol(posterior_alpha))
            for (i in 1:nrow(alpha_us)) {
               for (j in 1:ncol(alpha_us)) {
                    alpha_us[i,j] <- pnorm(posterior_alpha[i,j], 
                                           mean = posterior_mu[i] ,
                                           sd = 5
                                           )
               }
            }
        } else if (model == 'model_2') {
            # posterior_mu <- dap_samples[, startsWith(colnames(dap_samples), 'mu')]
            posterior_tau <- dap_samples[, startsWith(colnames(dap_samples), 'tau')]

            alpha_us <- matrix(0, nrow(posterior_alpha), ncol(posterior_alpha))
            for (i in 1:nrow(alpha_us)) {
               for (j in 1:ncol(alpha_us)) {
                    alpha_us[i,j] <- pnorm(posterior_alpha[i,j], 
                                           mean = posterior_mu[i, posterior_Z[i,j]+1] ,
                                           sd = 1 / sqrt(posterior_tau[i, posterior_Z[i,j]+1]) 
                                           )
               }
            }
        }

        n_mats <- 20
        resid_mat_file <- file.path(model_results_dir, 'data_u_object.rds')
        if (file.exists(resid_mat_file)) {
            u_resid_mats <- readRDS(resid_mat_file)
        } else {
            u_resid_mats <- list(type = 'vector', length = n_mats)
            for (i in 1:n_mats) {
                print(glue::glue('generating data u matrix {i} of {n_mats}'))
                u_resid_mat <- matrix(0, length(dap_samples[,'alpha[1]']), stan_data$n)
                for (j in 1:ncol(u_resid_mat)) {
                    if (smoking[j, 'smkreg'] == 1) {
                        u_resid_mat[,j] <- mapply(runif, n = 1, min = 1-posterior_p[,j], max = 1)
                    } else {
                        u_resid_mat[,j] <- mapply(runif, n = 1, min = 0, max = 1-posterior_p[,j])
                    }
                }
                u_resid_mats[[i]] <- u_resid_mat
            }
            saveRDS(u_resid_mats, resid_mat_file)
        }
        u_resid_mat <- u_resid_mats[[1]]
        # u_resid_mat <- Reduce('+', u_resid_mats) / length(u_resid_mats) # get avg u value over n_mats runs

        # fwrite( u_resid_mat, file.path(model_results_dir, 'resid_us.csv') )

        # generate residual trajectory plot
        # residual_trajectory_plot(u_resid_mat, model_results_dir)

        # generate residual scatter over wave
        scatter_dt <- data.table('newid' = smoking$newid,
                                'sex' = smoking$`sex(1=F)`,
                                'parsmk' = smoking$parsmk,
                                'wave' = smoking$wave,
                                'resid' = colMeans(u_resid_mat)
                                )
        scatter_dt[, sex := ifelse(sex == sex_1_map, 'female', 'male')]
        scatter_dt[, parsmk := ifelse(parsmk == parsmk_1_map, '1', '0')]

        scatter_plot_sex <- ggplot(scatter_dt, aes(x = as.factor(wave), y = resid, 
                            group = as.factor(newid), color = as.factor(sex)
                            )) +
            geom_line(alpha = 0.15) +
            geom_point(alpha = 0.15) +
            labs(x = 'wave', y = 'data u-value') + 
            scale_color_manual(name = 'sex', values = c('male' = 'blue', 'female' = 'red')) +
            theme(
                        axis.title.x = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        # axis.ticks.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )

        ggsave(file.path(model_results_dir, 'resid_wave_trends_sex.png'), scatter_plot_sex, width = 12, height = 8, dpi = 300)
        rm(scatter_plot_sex)

        scatter_plot_parsmk <- ggplot(scatter_dt, aes(x = as.factor(wave), y = resid, 
                            group = as.factor(newid), color = as.factor(parsmk)
                            )) +
            geom_line(alpha = 0.15) +
            geom_point(alpha = 0.15) +
            labs(x = 'wave', y = 'data u-value') + 
            scale_color_manual(name = 'parsmk', values = c('0' = 'blue', '1' = 'red')) + 
            theme(
                        axis.title.x = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        # axis.ticks.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )

        ggsave(file.path(model_results_dir, 'resid_wave_trends_parsmk.png'), scatter_plot_parsmk, width = 12, height = 8, dpi = 300)
        rm(scatter_plot_parsmk)

        # combine p-values
        print(glue::glue('get residual uniformity p-value'))
        M <- 10
        resid_p_values <- sapply(u_resid_mats, uniformity_test) # p-values for resid uniformity testing
        resid_p_value <- Cauchy_p_merger(apply(resid_p_values, MARGIN = 2, FUN = Cauchy_p_merger))
        exchangeable_resid_p_value <- exchangeable_p_merger(apply(resid_p_values, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)

        test_file <- file.path(model_results_dir, 'test_output.txt')

        ##############################
        ### TESTING !!!!!!!!!!!!!!!!!
        ##############################

        # plot alphas w/o averaging over posterior samples 
        alpha_dt <- melt( data.table(alpha_us), variable.name = 'newid', value.name = 'alpha_u' )
        alpha_dt[, newid := str_remove_all(newid, 'alpha\\[')]
        alpha_dt[, newid := str_remove_all(newid, '\\]')]
        alpha_dt[, newid := as.numeric(newid)]

        raw_alpha_dt <- melt( data.table(posterior_alpha), variable.name = 'newid', value.name = 'raw_alpha' )
        raw_alpha_dt[, newid := str_remove_all(newid, 'alpha\\[')]
        raw_alpha_dt[, newid := str_remove_all(newid, '\\]')]
        raw_alpha_dt[, newid := as.numeric(newid)]

        if (model %in% c('model_2', 'model_2', 'model_3')) {
            Z_dt <- melt( data.table(posterior_Z), variable.name = 'newid', value.name = 'Z' )
            Z_dt[, newid := str_remove_all(newid, 'Z\\[')]
            Z_dt[, newid := str_remove_all(newid, '\\]')]
            Z_dt[, newid := as.numeric(newid)]
        }

        indiv_hist_dt <- unique(smoking[, c('newid', 'sex(1=F)', 'parsmk')])
        colnames(indiv_hist_dt) <- c('newid', 'sex', 'parsmk')

        # indiv_hist_dt <- merge(indiv_vars_dt, alpha_dt, by = 'newid')
        
        indiv_hist_dt[, sex := ifelse(sex == sex_1_map, 'female', 'male')]
        indiv_hist_dt[, parsmk := ifelse(sex == parsmk_1_map, '1', '0')]
        indiv_hist_dt <- indiv_hist_dt[rep(seq_len(nrow(indiv_hist_dt)), each = nrow(posterior_alpha)), ]
        indiv_hist_dt[, sample := rep(1:nrow(posterior_alpha), times = ncol(posterior_alpha))]

        indiv_hist_dt[, alpha_u := alpha_dt$alpha_u]
        indiv_hist_dt[, raw_alpha := raw_alpha_dt$raw_alpha]

        if (model %in% c('model_2', 'model_2', 'model_3')) {
            indiv_hist_dt[, Z := as.factor(Z_dt$Z)]
        }

        true_smoker_status_dt <- smoking[, mean(smkreg), by = newid]

        indiv_hist_dt[, true_smoker_status := 'sometimes']
        indiv_hist_dt[newid %in% true_smoker_status_dt[V1 == 0, newid], true_smoker_status := 'never']
        indiv_hist_dt[newid %in% true_smoker_status_dt[V1 == 1, newid], true_smoker_status := 'always']

        indiv_hist_dt[, binary_true_smoker_status := ifelse(true_smoker_status == 'never', 'never', 'sometimes/always')]

        # write table of Z/binary_true_smoker_status to file
        if (model %in% c('model_2', 'model_2', 'model_3')) {
            table_file <- file.path(model_results_dir, 'Z_true_smoker_status_table.txt')
            Z_true_smoker_status_table <- table(indiv_hist_dt[, c('Z', 'true_smoker_status')])
            write.table(Z_true_smoker_status_table, table_file)
        
            # create table where each individual is assigned their most common Z
            table_2_file <- table_file <- file.path(model_results_dir, 'modal_Z_true_smoker_status_parsmk_table.txt')
            indiv_hist_dt[, modal_Z := (mean(Z == '1') >= 0.5), by = newid]
            unique_hist_dt <- unique(indiv_hist_dt[, c('newid', 'parsmk', 'true_smoker_status', 'modal_Z')])
            unique_hist_dt[, modal_Z := as.numeric(modal_Z)]
            modal_Z_true_smoker_status_parsmk_table <- table(unique_hist_dt[, c('true_smoker_status', 'parsmk', 'modal_Z')])
            write.table(modal_Z_true_smoker_status_parsmk_table, table_2_file)
        }

        # plot marginal hist/pdf/ecdf of alpha
        raw_alpha_min <- min(indiv_hist_dt$raw_alpha)
        raw_alpha_max <- max(indiv_hist_dt$raw_alpha)
        raw_alpha_diff <- raw_alpha_max - raw_alpha_min
        raw_alpha_hist_plot <- ggplot(indiv_hist_dt, aes(x = raw_alpha)) + 
                                geom_histogram(aes(y = after_stat(density)), breaks = seq(-20, 10, by = raw_alpha_diff/100), alpha = 1, fill = 'blue', colour = 'black') + 
                                labs(x = TeX('$\\alpha$'), y = 'density') +
                                theme(
                                    axis.title.x = element_text(size = 36),
                                    axis.text.x = element_text(size = 36),
                                    # axis.ticks.x = element_text(size = 36),
                                    axis.title.y = element_text(size = 36),
                                    axis.text.y = element_text(size = 36),
                                    # axis.ticks.y = element_blank(),
                                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                    panel.background = element_rect(fill = "white") # White background
                                    )
        alpha_avg_density_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u)) + 
                                geom_freqpoly(aes(y = after_stat(density)), breaks = seq(0, 1, by = 0.01), colour = 'blue')
        alpha_hist_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u)) + 
                                geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, by = 0.01), fill = 'blue', colour = 'black') + 
                                labs(x = TeX('$U_{\\alpha}$'), y = 'density') +
                                theme(
                        axis.title.x = element_text(size = 36),
                        axis.text.x = element_text(size = 36),
                        # axis.ticks.x = element_text(size = 36),
                        axis.title.y = element_text(size = 36),
                        axis.text.y = element_text(size = 36),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )
        alpha_pdf_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, group = sample)) +
                                        geom_density(linewidth = 0.1, alpha = 0.05, color = 'red')
        alpha_ecdf_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, group = sample)) +
                                        stat_ecdf(position = 'identity', linewidth = 1, alpha = 0.3, color = 'red') +
                                        geom_abline(intercept = 0, slope = 1, colour = 'black') +
                                        coord_cartesian(xlim = c(0,1), ylim = c(0,1))
        
        # indiv_hist_dt[, tilted_ecdf := ecdf(alpha_u)(alpha_u) - alpha_u, by = sample]
        # indiv_hist_dt <- indiv_hist_dt[order(sample, tilted_ecdf)]
        # alpha_tilted_ecdf_plot <- ggplot(indiv_hist_dt, aes(x = value, y = tilted_ecdf, group = sample)) + 
        #                                 geom_line(alpha = 0.3, colour = 'blue') + 
        #                                 geom_abline(intercept = 0, slope = 0, colour = 'black') + 
        #                                 annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
        #                                         x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 8) +
        #                                 coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
        #                                 labs(x = TeX('$U_{\\alpha}$'), y = 'Tilted Empirical CDF') +  
        #                                 guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        #                                 theme(axis.title.x = element_text(size = 12),
        #                                     axis.text.x = element_text(size = 12),,
        #                                     axis.title.y = element_text(size = 12),
        #                                     axis.text.y = element_text(size = 12),
        #                                     panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
        #                                     panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
        #                                     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
        #                                     panel.background = element_rect(fill = "white"), # White background
        #                                     plot.margin = margin(c(1,1,1,1), 'cm')
        #                                     )
        

        # alpha_tilted_ecdf_plot
        
        ggsave(file.path(model_results_dir, 'raw_alpha_hist.png'), raw_alpha_hist_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_avg_density.png'), alpha_avg_density_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_hist.png'), alpha_hist_plot, width = 16, height = 6, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_pdf.png'), alpha_pdf_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_ecdf.png'), alpha_ecdf_plot, width = 12, height = 8, dpi = 300)
        # ggsave(file.path(model_results_dir, 'alpha_tilted_ecdf.png'), alpha_tilted_ecdf_plot, width = 12, height = 8, dpi = 300)

        rm(raw_alpha_hist_plot, alpha_hist_plot, alpha_pdf_plot, alpha_ecdf_plot)

        # plot marginal hist/pdf/ecdf of alpha by true smoker status
        ss_colors <- met.brewer('Kandinsky', 3)

        ss_dt <- copy(indiv_hist_dt) 
        ss_dt[true_smoker_status == 'never', true_smoker_status := 'never smoker']
        ss_dt[true_smoker_status == 'sometimes', true_smoker_status := 'sometimes smoker']
        ss_dt[true_smoker_status == 'always', true_smoker_status := 'always smoker']
        ss_dt[, true_smoker_status := factor(true_smoker_status, levels = c('never smoker', 'sometimes smoker', 'always smoker'))]

        raw_alpha_ss_hist_plot <- ggplot(ss_dt, aes(x = raw_alpha, fill = true_smoker_status)) + 
                                geom_histogram(aes(y = after_stat(density)), breaks = seq(-20, 10, by = raw_alpha_diff/100), alpha = 1, colour = 'black') + 
                                scale_fill_manual(values = ss_colors) + 
                                guides(fill = 'none') +
                                labs(x = TeX('$\\alpha$'), y = 'density') +
                                theme(
                                    strip.text = element_text(size = 18),
                                    axis.title.x = element_text(size = 24),
                                    axis.text.x = element_text(size = 18),
                                    # axis.ticks.x = element_text(size = 24),
                                    axis.title.y = element_text(size = 24),
                                    axis.text.y = element_text(size = 24),
                                    # axis.ticks.y = element_blank(),
                                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                    panel.background = element_rect(fill = "white") # White background
                                    ) +
                                facet_grid(rows=vars(true_smoker_status))

        alpha_ss_hist_plot <- ggplot(ss_dt, aes(x = alpha_u, fill = true_smoker_status)) + 
                                geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, by = 0.01), alpha = 1, colour = 'black') + 
                                # guides(fill = guide_legend(title = 'Smoking Status')) + 
                                scale_fill_manual(values = ss_colors) + 
                                guides(fill = 'none') +
                                labs(x = TeX('$U_{\\alpha}$'), y = 'density') +
                                theme(
                                    strip.text = element_text(size = 24),
                                    axis.title.x = element_text(size = 24),
                                    axis.text.x = element_text(size = 18),
                                    # axis.ticks.x = element_text(size = 24),
                                    axis.title.y = element_text(size = 24),
                                    axis.text.y = element_text(size = 24),
                                    # axis.ticks.y = element_blank(),
                                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                    panel.background = element_rect(fill = "white") # White background
                                    ) + 
                                    facet_grid(~true_smoker_status)

        alpha_ss_pdf_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, colour = true_smoker_status, group = sample)) +
                                        geom_density(linewidth = 0.1, alpha = 0.05) + 
                                        facet_grid(~true_smoker_status)
        alpha_ss_ecdf_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, colour = true_smoker_status, group = sample)) +
                                        stat_ecdf(position = 'identity', linewidth = 1, alpha = 0.3) +
                                        geom_abline(intercept = 0, slope = 1, colour = 'black') +
                                        coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
                                        facet_grid(~true_smoker_status)
        
        ggsave(file.path(model_results_dir, 'raw_alpha_by_smoker_status_hist.png'), raw_alpha_ss_hist_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_by_smoker_status_hist.png'), alpha_ss_hist_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_by_smoker_status_pdf.png'), alpha_ss_pdf_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_by_smoker_status_ecdf.png'), alpha_ss_ecdf_plot, width = 12, height = 8, dpi = 300)

        rm(raw_alpha_ss_hist_plot, alpha_ss_hist_plot, alpha_ss_pdf_plot, alpha_ss_ecdf_plot)

        if (model %in% c('model_2', 'model_2', 'model_3')) {
            Z_colors <- met.brewer('Lakota', 2)
            indiv_hist_dt_copy <- copy(indiv_hist_dt)
            indiv_hist_dt_copy[Z == 0, Z_string := 'Z = 0']
            indiv_hist_dt_copy[Z == 1, Z_string := 'Z = 1']
            alpha_Z_hist_plot <- ggplot(indiv_hist_dt_copy, aes(x = alpha_u, fill = Z_string)) + 
                                geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, by = 0.01), colour = 'black') + 
                                scale_fill_manual(values = Z_colors) + 
                                guides(fill = 'none') +
                                labs(x = TeX('$U_{\\alpha}$'), y = 'density') +
                                theme(
                                    strip.text = element_text(size = 24),
                                    axis.title.x = element_text(size = 24),
                                    axis.text.x = element_text(size = 18),
                                    # axis.ticks.x = element_text(size = 24),
                                    axis.title.y = element_text(size = 24),
                                    axis.text.y = element_text(size = 24),
                                    # axis.ticks.y = element_blank(),
                                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                    panel.background = element_rect(fill = "white") # White background
                                    ) + 
                                facet_grid(~Z_string)

            alpha_Z_pdf_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, colour = Z, group = sample)) +
                                        geom_density(linewidth = 0.1, alpha = 0.05) +
                                        scale_colour_manual(name = 'Z', values = c('0' = 'blue', '1' = 'red')) +
                                        # geom_abline(intercept = 0, slope = 1, colour = 'black') +
                                        # coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
                                        facet_grid(~Z)
            
            alpha_Z_ecdf_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, colour = Z, group = sample)) +
                                        stat_ecdf(position = 'identity', linewidth = 1, alpha = 0.3) +
                                        scale_colour_manual(name = 'Z', values = c('0' = 'blue', '1' = 'red')) +
                                        geom_abline(intercept = 0, slope = 1, colour = 'black') +
                                        coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
                                        facet_grid(~Z)
            
            # raw_alpha_Z_hist_plot <- ggplot(indiv_hist_dt, aes(x = raw_alpha, fill = Z)) + 
            #                     geom_histogram(aes(y = after_stat(density)), breaks = seq(raw_alpha_min, raw_alpha_max, by = raw_alpha_diff/100), alpha = 0.3, colour = 'black') + 
            #                     scale_fill_manual(name = 'Z', values = c('0' = 'blue', '1' = 'red')) + 
            #                     facet_grid(~Z)
            raw_alpha_Z_hist_plot <- ggplot(indiv_hist_dt_copy, aes(x = raw_alpha, fill = Z_string)) + 
                                geom_histogram(aes(y = after_stat(density)), breaks = seq(-20, 10, by = raw_alpha_diff/100), alpha = 1, colour = 'black') + 
                                scale_fill_manual(values = Z_colors) + 
                                guides(fill = 'none') +
                                labs(x = TeX('$\\alpha$'), y = 'density') +
                                theme(
                                    strip.text = element_text(size = 24),
                                    axis.title.x = element_text(size = 24),
                                    axis.text.x = element_text(size = 18),
                                    # axis.ticks.x = element_text(size = 24),
                                    axis.title.y = element_text(size = 24),
                                    axis.text.y = element_text(size = 24),
                                    # axis.ticks.y = element_blank(),
                                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                    panel.background = element_rect(fill = "white") # White background
                                    ) + 
                                facet_grid(rows=vars(Z_string))


            ggsave(file.path(model_results_dir, 'raw_alpha_by_Z_hist.png'), raw_alpha_Z_hist_plot, width = 12, height = 8, dpi = 300)
            ggsave(file.path(model_results_dir, 'alpha_by_Z_hist.png'), alpha_Z_hist_plot, width = 12, height = 8, dpi = 300)
            ggsave(file.path(model_results_dir, 'alpha_by_Z_pdf.png'), alpha_Z_pdf_plot, width = 12, height = 8, dpi = 300)
            ggsave(file.path(model_results_dir, 'alpha_by_Z_ecdf.png'), alpha_Z_ecdf_plot, width = 12, height = 8, dpi = 300)

            rm(raw_alpha_Z_hist_plot, alpha_Z_hist_plot, alpha_Z_pdf_plot, alpha_Z_ecdf_plot)
        }
        
        alpha_sex_hist_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, fill = sex)) + 
                                geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, by = 0.01), alpha = 0.3, colour = 'black') + 
                                scale_fill_manual(name = 'sex', values = c('male' = 'blue', 'female' = 'red')) + 
                                facet_grid(~sex)

        alpha_sex_pdf_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, colour = sex, group = sample)) +
                                        geom_density(linewidth = 0.1, alpha = 0.05) +
                                        scale_colour_manual(name = 'sex', values = c('male' = 'blue', 'female' = 'red')) +
                                        facet_grid(~sex)

        alpha_sex_ecdf_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, colour = sex, group = sample)) +
                                    stat_ecdf(position = 'identity', linewidth = 1, alpha = 0.3) +
                                    scale_colour_manual(name = 'sex', values = c('male' = 'blue', 'female' = 'red')) +
                                    geom_abline(intercept = 0, slope = 1, colour = 'black') +
                                    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
                                    facet_grid(~sex)

        ggsave(file.path(model_results_dir, 'alpha_by_sex_hist.png'), alpha_sex_hist_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_by_sex_pdf.png'), alpha_sex_pdf_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_by_sex_ecdf.png'), alpha_sex_ecdf_plot, width = 12, height = 8, dpi = 300)

        rm(alpha_sex_hist_plot, alpha_sex_pdf_plot, alpha_sex_ecdf_plot)

        # indiv_hist_dt$parsmk <- ifelse(as.numeric(indiv_hist_dt$parsmk) == min(as.numeric(indiv_hist_dt$parsmk)), '0', '1')
        alpha_parsmk_hist_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, fill = parsmk)) + 
                                geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 1, by = 0.01), alpha = 0.3, colour = 'black') + 
                                scale_fill_manual(name = 'parsmk', values = c('0' = 'blue', '1' = 'red')) + 
                                facet_grid(~parsmk)
        
        alpha_parsmk_pdf_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, colour = parsmk, group = sample)) +
                                        geom_density(linewidth = 0.1, alpha = 0.05) +
                                        scale_color_manual(name = 'parsmk', values = c('0' = 'blue', '1' = 'red')) +
                                        facet_grid(~parsmk)

        alpha_parsmk_ecdf_plot <- ggplot(indiv_hist_dt, aes(x = alpha_u, colour = parsmk, group = sample)) +
                                    stat_ecdf(position = 'identity', linewidth = 1, alpha = 0.3) +
                                    scale_color_manual(name = 'parsmk', values = c('0' = 'blue', '1' = 'red')) +
                                    geom_abline(intercept = 0, slope = 1, colour = 'black') +
                                    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
                                    facet_grid(~parsmk)

        ggsave(file.path(model_results_dir, 'alpha_by_parsmk_hist.png'), alpha_parsmk_hist_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_by_parsmk_pdf.png'), alpha_parsmk_pdf_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, 'alpha_by_parsmk_ecdf.png'), alpha_parsmk_ecdf_plot, width = 12, height = 8, dpi = 300)

        rm(indiv_hist_dt, alpha_parsmk_hist_plot, alpha_parsmk_pdf_plot, alpha_parsmk_ecdf_plot)

        ###############################
        ### END TESTING 
        ###############################

        # test independence of covariates and random effects
        print(glue::glue('test independence of covariates and random effects'))
        parsmk_alpha_indep_file <- file.path(model_results_dir, 'parsmk_alpha_indep.rds')
        sex_alpha_indep_file <- file.path(model_results_dir, 'sex_alpha_indep.rds')

        if ( file.exists(parsmk_alpha_indep_file) ) {
           parsmk_alpha_ps <- readRDS(parsmk_alpha_indep_file)
           sex_alpha_ps <- readRDS(sex_alpha_indep_file) 
        } else {
            parsmk_alpha_ps <- rep(0, nrow(alpha_us))
            sex_alpha_ps <- rep(0, nrow(alpha_us))
            for (i in 1:nrow(alpha_us)) {
                parsmk_alpha_ps[i] <- wilcox.test(as.numeric(alpha_us[i,])[parsmk_vec == parsmk_0_map],
                                                as.numeric(alpha_us[i,])[parsmk_vec == parsmk_1_map]
                                                )$p.value # testing against parsmk_vec bc alpha is defined at the level of the individual
                sex_alpha_ps[i] <- wilcox.test(as.numeric(alpha_us[i,])[sex_vec == sex_0_map],
                                            as.numeric(alpha_us[i,])[sex_vec == sex_1_map]
                                            )$p.value
            }
            saveRDS(parsmk_alpha_ps, parsmk_alpha_indep_file)
            saveRDS(sex_alpha_ps, sex_alpha_indep_file)
        }

        print(glue::glue('test independence of covariates and residuals'))
        parsmk_resid_indep_file <- file.path(model_results_dir, 'parsmk_resid_indep.rds')
        sex_resid_indep_file <- file.path(model_results_dir, 'sex_resid_indep.rds')

        n_mats <- length(u_resid_mats)

        if ( file.exists(parsmk_resid_indep_file) ) {
           parsmk_resid_ps <- readRDS(parsmk_resid_indep_file)
           sex_resid_ps <- readRDS(sex_resid_indep_file) 
        } else {
            parsmk_resid_ps <- matrix(0, nrow(alpha_us), n_mats)
            sex_resid_ps <- matrix(0, nrow(alpha_us), n_mats)
            for (j in 1:n_mats) {
                urm <- u_resid_mats[[j]]
                print(glue::glue('run {j} of {n_mats}'))
                for (i in 1:nrow(urm)) {
                    parsmk_resid_ps[i,j] <- wilcox.test(as.numeric(urm[i,])[smoking$parsmk == parsmk_0_map],
                                                as.numeric(urm[i,])[smoking$parsmk == parsmk_1_map]
                                            )$p.value # testing against smoking$parsmk bc the residuals are defined at the level of the observation
                    sex_resid_ps[i,j] <- wilcox.test(as.numeric(urm[i,])[smoking$sex == sex_0_map],
                                                as.numeric(urm[i,])[smoking$sex == sex_1_map]
                                                )$p.value
                }
            }
            saveRDS(parsmk_resid_ps, parsmk_resid_indep_file)
            saveRDS(sex_resid_ps, sex_resid_indep_file)
        }

        parsmk_alpha_p <- Cauchy_p_merger(parsmk_alpha_ps)
        sex_alpha_p <- Cauchy_p_merger(sex_alpha_ps)
        parsmk_resid_p <- Cauchy_p_merger(apply(parsmk_resid_ps, MARGIN = 2, FUN = Cauchy_p_merger))
        sex_resid_p <- Cauchy_p_merger(apply(sex_resid_ps, MARGIN = 2, FUN = Cauchy_p_merger))

        exchangeable_parsmk_alpha_p <- exchangeable_p_merger(parsmk_alpha_ps, 10)
        exchangeable_sex_alpha_p <- exchangeable_p_merger(sex_alpha_ps, 10)
        exchangeable_parsmk_resid_p <- exchangeable_p_merger(apply(parsmk_resid_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)
        exchangeable_sex_resid_p <- exchangeable_p_merger(apply(sex_resid_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)

        rm(parsmk_alpha_ps, sex_alpha_ps, parsmk_resid_ps, sex_resid_ps)

        ### test independence of wave and residuals ###
        print(glue::glue('test independence of wave and residuals'))
        wave_resid_indep_file <- file.path(model_results_dir, 'wave_resid_indep.rds')

        # NOTE: hoeffd takes a very long time, so we'll evaluate a random 100 posterior samples

        sample_indices <- sample(1:nrow(u_resid_mat), size = 120, replace = FALSE)
        # sample_indices <- 1:nrow(u_resid_mat)

        if ( file.exists(wave_resid_indep_file) ) {
            wave_resid_ps <- readRDS(wave_resid_indep_file)
        } else {
            # wave_resid_ps <- matrix(0, nrow(u_resid_mat), n_mats)
            wave_resid_ps <- matrix(0, length(sample_indices), n_mats)
            for (j in 1:n_mats) {
                print(glue::glue('run {j} of {n_mats}'))
                urm <- u_resid_mats[[j]]
                for (i in 1:length(sample_indices)) {
                    # print(glue::glue('sample {i} of {nrow(u_resid_mat)}'))
                    # wave_resid_ps[i,j] <-  kruskal.test(urm[i,] ~ smoking$wave)$p.value
                    wave_resid_ps[i,j] <- hoeffd(urm[sample_indices[i], ], smoking$wave)$P[1,2]
                }
            }
            saveRDS(wave_resid_ps, wave_resid_indep_file)
        }

        wave_resid_p <- Cauchy_p_merger(apply(wave_resid_ps, MARGIN = 2, FUN = Cauchy_p_merger))
        exchangeable_wave_resid_p <- exchangeable_p_merger(apply(wave_resid_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)

        rm(wave_resid_ps)

        ### test independence of wave*sex and residuals ###
        print(glue::glue('test independence of wave*sex and residuals'))
        wave_x_sex_resid_indep_file <- file.path(model_results_dir, 'wave_x_sex_resid_indep.rds')

        if ( file.exists(wave_x_sex_resid_indep_file) ) {
            wave_x_sex_resid_ps <- readRDS(wave_x_sex_resid_indep_file)
        } else {
            # wave_x_sex_resid_ps <- matrix(0, nrow(u_resid_mat), n_mats)
            wave_x_sex_resid_ps <- matrix(0, length(sample_indices), n_mats)
            for (j in 1:n_mats) {
                print(glue::glue('run {j} of {n_mats}'))
                urm <- u_resid_mats[[j]]
                for (i in 1:length(sample_indices)) {
                    wave_x_sex_resid_ps[i,j] <- hoeffd(urm[sample_indices[i], ], smoking$wave*smoking$sex)$P[1,2]
                }
            }
            saveRDS(wave_x_sex_resid_ps, wave_x_sex_resid_indep_file)
        }

        wave_x_sex_resid_p <- Cauchy_p_merger(apply(wave_x_sex_resid_ps, MARGIN = 2, FUN = Cauchy_p_merger))
        exchangeable_wave_x_sex_resid_p <- exchangeable_p_merger(apply(wave_x_sex_resid_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)

        rm(wave_x_sex_resid_ps)

        ### test independence of wave*parsmk and residuals ###
        print(glue::glue('test independence of wave*parsmk and residuals'))
        wave_x_parsmk_resid_indep_file <- file.path(model_results_dir, 'wave_x_parsmk_resid_indep.rds')

        if ( file.exists(wave_x_parsmk_resid_indep_file) ) {
            wave_x_parsmk_resid_ps <- readRDS(wave_x_parsmk_resid_indep_file)
        } else {
            # wave_x_parsmk_resid_ps <- matrix(0, nrow(u_resid_mat), n_mats)
            wave_x_parsmk_resid_ps <- matrix(0, length(sample_indices), n_mats)
            for (j in 1:n_mats) {
                print(glue::glue('run {j} of {n_mats}'))
                urm <- u_resid_mats[[j]]
                for (i in 1:length(sample_indices)) {
                    wave_x_parsmk_resid_ps[i,j] <- hoeffd(urm[sample_indices[i], ], smoking$wave*smoking$parsmk)$P[1,2]
                }
            }
            saveRDS(wave_x_parsmk_resid_ps, wave_x_parsmk_resid_indep_file)
        }

        wave_x_parsmk_resid_p <- Cauchy_p_merger(apply(wave_x_parsmk_resid_ps, MARGIN = 2, FUN = Cauchy_p_merger))
        exchangeable_wave_x_parsmk_resid_p <- exchangeable_p_merger(apply(wave_x_parsmk_resid_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)

        rm(wave_x_parsmk_resid_ps)

        # ### subset to female ### 
        # print('...female only')
        # female_wave_resid_indep_file <- file.path(model_results_dir, 'wave_resid_indep_female.rds')

        # if ( file.exists(female_wave_resid_indep_file) ) {
        #     female_wave_resid_ps <- readRDS(female_wave_resid_indep_file)
        # } else {
        #     # female_wave_resid_ps <- matrix(0, nrow(u_resid_mat), n_mats)
        #     female_wave_resid_ps <- matrix(0, length(sample_indices), n_mats)
        #     for (j in 1:n_mats) {
        #         print(glue::glue('run {j} of {n_mats}'))
        #         urm <- u_resid_mats[[j]]
        #         for (i in 1:length(sample_indices)) {
        #             # wave_resid_ps[i,j] <- kruskal.test(urm[i, (smoking_sex_1)] ~ female_wave)$p.value
        #             female_wave_resid_ps[i,j] <- hoeffd(urm[sample_indices[i], (smoking_sex_1)], female_wave)$P[1,2]
        #         }
        #     }
        #     saveRDS(female_wave_resid_ps, female_wave_resid_indep_file)
        # }

        # female_wave_resid_p <- Cauchy_p_merger(apply(female_wave_resid_ps, MARGIN = 2, FUN = Cauchy_p_merger))
        # exchangeable_female_wave_resid_p <- exchangeable_p_merger(apply(female_wave_resid_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)

        # rm(female_wave_resid_ps)

        # ### subset to male ### 
        # print('...male only')
        # male_wave_resid_indep_file <- file.path(model_results_dir, 'wave_resid_indep_male.rds')

        # if ( file.exists(male_wave_resid_indep_file) ) {
        #     male_wave_resid_ps <- readRDS(male_wave_resid_indep_file)
        # } else {
        #     # male_wave_resid_ps <- matrix(0, nrow(u_resid_mat), n_mats)
        #     male_wave_resid_ps <- matrix(0, length(sample_indices), n_mats)
        #     for (j in 1:n_mats) {
        #         print(glue::glue('run {j} of {n_mats}'))
        #         urm <- u_resid_mats[[j]]
        #         for (i in 1:length(sample_indices)) {
        #             male_wave_resid_ps[i,j] <- hoeffd(urm[sample_indices[i], (smoking_sex_0)], male_wave)$P[1,2]
        #         }
        #     }
        #     saveRDS(male_wave_resid_ps, male_wave_resid_indep_file)
        # }

        # male_wave_resid_p <- Cauchy_p_merger(apply(male_wave_resid_ps, MARGIN = 2, FUN = Cauchy_p_merger))
        # exchangeable_male_wave_resid_p <- exchangeable_p_merger(apply(male_wave_resid_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)

        # rm(male_wave_resid_ps)

        ### test for interaction between wave*covariate ###
        print('test for interactions between wave and covariates')
        sex_wave_marginal_file <- file.path(model_results_dir, 'sex_wave_marginal_file.rds')
        sex_wave_interaction_file <- file.path(model_results_dir, 'sex_wave_interaction_file.rds')
        if ( file.exists(sex_wave_interaction_file) ) {
            sex_wave_marginal_ps <- readRDS(sex_wave_marginal_file)
            sex_wave_interaction_ps <- readRDS(sex_wave_interaction_file)
        } else {
            sex_wave_marginal_ps <- matrix(0, nrow(u_resid_mat), n_mats)
            sex_wave_interaction_ps <- matrix(0, nrow(u_resid_mat), n_mats)
            for (j in 1:n_mats) {
                print(glue::glue('run {j} of {n_mats}'))
                urm <- u_resid_mats[[j]]
                for (i in 1:nrow(u_resid_mat)) {
                    l <- lm(logit(urm[i,]) ~ smoking$wave + smoking$sex + smoking$sex*smoking$wave)
                    sex_wave_marginal_ps[i, j] <- as.numeric(tidy(l)[3,'p.value'])
                    sex_wave_interaction_ps[i, j] <- as.numeric(tidy(l)[4,'p.value'])
                }
            }
            saveRDS(sex_wave_marginal_ps, sex_wave_marginal_file)
            saveRDS(sex_wave_interaction_ps, sex_wave_interaction_file)
        }

        sex_wave_interaction_p <- Cauchy_p_merger(apply(sex_wave_interaction_ps, MARGIN = 2, FUN = Cauchy_p_merger))
        exchangeable_sex_wave_interaction_p <- exchangeable_p_merger(apply(sex_wave_interaction_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)
        sex_wave_marginal_p <- Cauchy_p_merger(apply(sex_wave_marginal_ps, MARGIN = 2, FUN = Cauchy_p_merger))
        exchangeable_sex_wave_marginal_p <- exchangeable_p_merger(apply(sex_wave_marginal_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)

        parsmk_wave_marginal_file <- file.path(model_results_dir, 'parsmk_wave_marginal_file.rds')
        parsmk_wave_interaction_file <- file.path(model_results_dir, 'parsmk_wave_interaction_file.rds')
        if ( file.exists(parsmk_wave_interaction_file) ) {
            parsmk_wave_marginal_ps <- readRDS(parsmk_wave_marginal_file)
            parsmk_wave_interaction_ps <- readRDS(parsmk_wave_interaction_file)
        } else {
            parsmk_wave_marginal_ps <- matrix(0, nrow(u_resid_mat), n_mats)
            parsmk_wave_interaction_ps <- matrix(0, nrow(u_resid_mat), n_mats)
            for (j in 1:n_mats) {
                print(glue::glue('run {j} of {n_mats}'))
                urm <- u_resid_mats[[j]]
                for (i in 1:nrow(u_resid_mat)) {
                    l <- lm(logit(urm[i,]) ~ smoking$wave + smoking$parsmk + smoking$parsmk*smoking$wave)
                    parsmk_wave_marginal_ps[i, j] <- as.numeric(tidy(l)[3,'p.value'])
                    parsmk_wave_interaction_ps[i, j] <- as.numeric(tidy(l)[4,'p.value'])
                }
            }
            saveRDS(parsmk_wave_marginal_ps, parsmk_wave_marginal_file)
            saveRDS(parsmk_wave_interaction_ps, parsmk_wave_interaction_file)
        }

        parsmk_wave_interaction_p <- Cauchy_p_merger(apply(parsmk_wave_interaction_ps, MARGIN = 2, FUN = Cauchy_p_merger))
        exchangeable_parsmk_wave_interaction_p <- exchangeable_p_merger(apply(parsmk_wave_interaction_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)
        
        parsmk_wave_marginal_p <- Cauchy_p_merger(apply(parsmk_wave_marginal_ps, MARGIN = 2, FUN = Cauchy_p_merger))
        exchangeable_parsmk_wave_marginal_p <- exchangeable_p_merger(apply(parsmk_wave_marginal_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10), M = 10)

        # test independence of covar and residuals in each wave j
        print(glue::glue('test independence of covariates and residuals in each wave'))

        parsmk_resid_wave_indep_file <- file.path(model_results_dir, 'parsmk_resid_wave_indep.rds')
        sex_resid_wave_indep_file <- file.path(model_results_dir, 'sex_resid_wave_indep.rds')
        exchangeable_parsmk_resid_wave_indep_file <- file.path(model_results_dir, 'exchangeable_parsmk_resid_wave_indep.rds')
        exchangeable_sex_resid_wave_indep_file <- file.path(model_results_dir, 'exchangeable_sex_resid_wave_indep.rds')

        if (file.exists(parsmk_resid_wave_indep_file)) {
            sex_wave_p_mat <- readRDS(sex_resid_wave_indep_file)
            parsmk_wave_p_mat <- readRDS(parsmk_resid_wave_indep_file)
            exchangeable_sex_wave_p_mat <- readRDS(exchangeable_sex_resid_wave_indep_file)
            exchangeable_parsmk_wave_p_mat <- readRDS(exchangeable_parsmk_resid_wave_indep_file)
        } else {
            sex_wave_p_mat <- matrix(0, n_mats, 6)
            parsmk_wave_p_mat <- matrix(0, n_mats, 6)
            exchangeable_sex_wave_p_mat <- matrix(0, n_mats, 6)
            exchangeable_parsmk_wave_p_mat <- matrix(0, n_mats, 6)
            
            for (k in 1:n_mats) {
                print(glue::glue('run {k} of {n_mats}'))
                urm <- u_resid_mats[[k]]

                sex_wave_ps <- matrix(0, nrow(urm), 6)
                parsmk_wave_ps <- matrix(0, nrow(urm), 6)
                for ( j in 1:6 ) {
                    for (i in 1:nrow(urm)) {
                        wave_map <- unique(smoking$wave)[j]
                        sex_wave_ps[i,j] <- wilcox.test(as.numeric(urm[i,])[smoking$sex == sex_0_map & smoking$wave == wave_map],
                                                        as.numeric(urm[i,])[smoking$sex == sex_1_map & smoking$wave == wave_map]
                                                    )$p.value
                        parsmk_wave_ps[i,j] <- wilcox.test(as.numeric(urm[i,])[smoking$parsmk == parsmk_0_map & smoking$wave == wave_map],
                                                        as.numeric(urm[i,])[smoking$parsmk == parsmk_1_map & smoking$wave == wave_map]
                                                    )$p.value
                    }
                }

                sex_wave_p_mat[k,] <- apply(sex_wave_ps, MARGIN = 2, FUN = Cauchy_p_merger)
                parsmk_wave_p_mat[k,] <- apply(parsmk_wave_ps, MARGIN = 2, FUN = Cauchy_p_merger)
                exchangeable_sex_wave_p_mat[k,] <- apply(sex_wave_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10)
                exchangeable_parsmk_wave_p_mat[k,] <- apply(parsmk_wave_ps, MARGIN = 2, FUN = exchangeable_p_merger, M = 10)
            }
            saveRDS(sex_wave_p_mat, sex_resid_wave_indep_file)
            saveRDS(parsmk_wave_p_mat, parsmk_resid_wave_indep_file)
            saveRDS(exchangeable_sex_wave_p_mat, exchangeable_sex_resid_wave_indep_file)
            saveRDS(exchangeable_parsmk_wave_p_mat, exchangeable_parsmk_resid_wave_indep_file)

            rm(sex_wave_ps, parsmk_wave_ps)
        }

        sex_wave_p <- apply(sex_wave_p_mat, MARGIN = 2, FUN = Cauchy_p_merger)
        parsmk_wave_p <- apply(parsmk_wave_p_mat, MARGIN = 2, FUN = Cauchy_p_merger)
        exchangeable_sex_wave_p <- apply(exchangeable_sex_wave_p_mat, MARGIN = 2, FUN = exchangeable_p_merger, M = 10)
        exchangeable_parsmk_wave_p <- apply(exchangeable_parsmk_wave_p_mat, MARGIN = 2, FUN = exchangeable_p_merger, M = 10)

        combined_sex_wave_p <- paste(sex_wave_p, collapse = ', ')
        combined_parsmk_wave_p <- paste(parsmk_wave_p, collapse = ', ')
        exchangeable_combined_sex_wave_p <- paste(exchangeable_sex_wave_p, collapse = ', ')
        exchangeable_combined_parsmk_wave_p <- paste(exchangeable_parsmk_wave_p, collapse = ', ')

        # test uniformity of alpha_us 
        alpha_uniformity_file <- file.path(model_results_dir, 'alpha_uniformity.rds')

        if (file.exists(alpha_uniformity_file)) {
            alpha_uniformity_ps <- readRDS(alpha_uniformity_file)
        } else {
            alpha_uniformity_ps <- uniformity_test(alpha_us)
            saveRDS(alpha_uniformity_ps, alpha_uniformity_file)
        }
        
        alpha_uniformity_p <- Cauchy_p_merger(alpha_uniformity_ps)
        exchangeable_alpha_uniformity_p <- exchangeable_p_merger(alpha_uniformity_ps, 10)

        # create trace plot of alpha_uniformity ps 
        alpha_uniformity_dt <- data.table('posterior_sample' = 1:length(alpha_uniformity_ps), 
                                          'alpha_uniformity_p' = alpha_uniformity_ps
                                         )
        alpha_uniformity_p_traceplot <- ggplot(alpha_uniformity_dt, aes(x = posterior_sample, y = alpha_uniformity_p)) + 
                                            geom_line()
        ggsave( file.path(model_results_dir, 'alpha_uniformity_p_traceplot.png'), alpha_uniformity_p_traceplot, width = 12, height = 8, dpi = 300 )

        rm(alpha_uniformity_ps, alpha_uniformity_dt, alpha_uniformity_p_traceplot)

        ### plot alpha us vs data us ###
        # get mean data u value for each individual
        u_data_dt <- data.table('newid' = smoking$newid,
                                'data' = colMeans(u_resid_mat)
                                )
        u_data_indiv_dt <- u_data_dt[, mean(data), by = newid]

        alpha_data_mean_u_dt <- data.table('alpha' = colMeans(alpha_us),
                                        'data' = u_data_indiv_dt$V1
                                        )
        if (model %in% c('model_2', 'model_2', 'model_3')) {
            alpha_data_mean_u_dt[, Z := colMeans(posterior_Z)]
        }

        alpha_data_mean_u_plot <- ggplot(alpha_data_mean_u_dt, aes(x = alpha, y = data)) + 
                                        geom_point() + 
                                        labs(x = 'avg alpha u', y = 'avg data u')
        
        if (model %in% c('model_2', 'model_2', 'model_3')) {
            Z_alpha_mean_u_plot <- ggplot(alpha_data_mean_u_dt, aes(x = Z, y = alpha)) + 
                                        geom_point() + 
                                        labs(x = 'avg Z', y = 'avg alpha u')
            Z_data_mean_u_plot <- ggplot(alpha_data_mean_u_dt, aes(x = Z, y = data)) + 
                                        geom_point() + 
                                        labs(x = 'avg Z', y = 'avg data u')
        }
        ggsave(file.path(model_results_dir, 'mean_data_u_by_alpha.png'), alpha_data_mean_u_plot, width = 12, height = 8, dpi = 300)
        
        if (model %in% c('model_2', 'model_2', 'model_3')) {
            ggsave(file.path(model_results_dir, 'mean_alpha_u_by_Z.png'), Z_alpha_mean_u_plot, width = 12, height = 8, dpi = 300)
            ggsave(file.path(model_results_dir, 'mean_data_u_by_Z.png'), Z_data_mean_u_plot, width = 12, height = 8, dpi = 300)

            rm(Z_alpha_mean_u_plot, Z_data_mean_u_plot)
        }

        rm(u_data_dt, u_data_indiv_dt, alpha_data_mean_u_dt, alpha_data_mean_u_plot)

        ### test independence of alpha u and mean data u within each individual ###

        # generate combined data set 
        u_resid_mat_long <- melt(data.table(u_resid_mat), variable.name = 'observation', value.name = 'data_u')
        u_resid_mat_long[, observation := str_remove(observation, 'V')]
        u_resid_mat_long[, posterior_index := rep(1:nrow(u_resid_mat), times = stan_data$n)]
        u_resid_mat_long[, newid := rep(smoking$newid, each = nrow(u_resid_mat))]
        u_resid_mat_long[, mean_data_u := mean(data_u), by = c('newid', 'posterior_index')]
        
        u_resid_mat_indiv <- unique(u_resid_mat_long[, c('posterior_index', 'newid', 'mean_data_u')])

        u_alpha_mat <- alpha_us
        u_alpha_mat_long <- melt(data.table(alpha_us), variable.name = 'newid', value.name = 'alpha_u')
        u_alpha_mat_long[, newid := str_remove(newid, 'V')]

        combined_resid_mat_indiv <- u_resid_mat_indiv
        combined_resid_mat_indiv[, alpha_u := u_alpha_mat_long$alpha_u]

        # test for independence within each posterior (use Hoeffding for this)
        resid_alpha_indep_p_values <- combined_resid_mat_indiv[, hoeffd(mean_data_u, alpha_u)$P[1,2], by = posterior_index]$V1
        resid_alpha_indep_p <- Cauchy_p_merger(resid_alpha_indep_p_values)
        exchangeable_resid_alpha_indep_p <- exchangeable_p_merger(resid_alpha_indep_p_values, 10)

        
        # add covariate info to u_resid_mat_long
        u_resid_mat_long[, wave := rep(smoking$wave, each = nrow(u_resid_mat_long)/stan_data$n )]
        u_resid_mat_long[, parsmk := rep(smoking$parsmk, each = nrow(u_resid_mat_long)/stan_data$n )]
        u_resid_mat_long[, sex := rep(smoking$sex, each = nrow(u_resid_mat_long)/stan_data$n )]
        u_resid_mat_long[, smkreg := rep(smoking$smkreg, each = nrow(u_resid_mat_long)/stan_data$n )]

        # ##################################
        # ### testing new visualizations ###
        # ##################################
        # sampled_indices <- sample( unique(u_resid_mat_long$observation), min( 2000, length(unique(u_resid_mat_long$observation)) ), replace = FALSE )
        # plot_mat_long <- u_resid_mat_long[ observation %in% sampled_indices ]
        # option_4_plot <- ggplot(plot_mat_long, aes(x = data_u, colour = as.factor(wave))) + 
        #     stat_ecdf(position = 'identity', linewidth = 0.5, alpha = 0.3) +
        #     geom_abline(intercept = 0, slope = 1, colour = 'black') +
        #     coord_cartesian(xlim = c(0,1), ylim = c(0,1))
        
        # resid_quantile_dt <- u_resid_mat_long[posterior_index == 1] %>%
        #     group_by(wave) %>%
        #     summarise(
        #         q1 = quantile(data_u, 0.25),
        #         q2 = quantile(data_u, 0.5),
        #         q3 = quantile(data_u, 0.75),
        #         .groups = 'drop'
        #     )
        # resid_quantile_dt_long <- melt(data.table(resid_quantile_dt), id.vars = c('wave'))
        
        # # option_2_plot <- ggplot(resid_quantile_dt_long, aes(x = as.factor(wave), y = value, colour = variable)) + 
        # #     geom_point()
        
        # resid_quantile_mean_dt <- u_resid_mat_long %>%
        #     group_by(wave, posterior_index) %>%
        #     summarise(
        #         q1 = quantile(data_u, 0.25),
        #         q2 = quantile(data_u, 0.5),
        #         q3 = quantile(data_u, 0.75),
        #         .groups = 'drop'
        #     ) %>%
        #     group_by(wave) %>%
        #     summarise(
        #         mean_q1 = mean(q1, na.rm = TRUE),
        #         mean_q2 = mean(q2, na.rm = TRUE),
        #         mean_q3 = mean(q3, na.rm = TRUE),
        #         .groups = 'drop'
        #     )
        # resid_quantile_mean_dt_long <- melt(data.table(resid_quantile_mean_dt), id.vars = c('wave'))
        # # option_3_plot <- ggplot(resid_quantile_mean_dt_long, aes(x = as.factor(wave), y = value, colour = variable)) + 
        # #     geom_point()

        ######################################
        ### end testing new visualizations ###
        ######################################

        # plot mean/se by wave/covar
        resid_wave_dt <- u_resid_mat_long %>%
            group_by(wave) %>%
            summarise(
                mean_resid = mean(data_u),
                se = sd(data_u) / sqrt(n()),
                .groups = 'drop'
            ) %>%
            mutate(
                lower_ci = mean_resid - qt(0.975, df=n()-1) * se,
                upper_ci = mean_resid + qt(0.975, df=n()-1) * se
            )

            wave_mean_resid_plot <- ggplot(resid_wave_dt, aes(x = as.factor(wave), y = mean_resid)) +
                                        geom_point(size = 3, position = position_dodge(width = 0.3)) +
                                        geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, 
                                                        position = position_dodge(width = 0.3)) +
                                        scale_x_discrete(labels = c('1', '2', '3', '4', '5', '6')) +
                                        labs(
                                            x = 'wave',
                                            y = 'mean data u-value'
                                        ) +
                                        theme(
                        axis.title.x = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        # axis.ticks.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )
            ggsave(file.path(model_results_dir, 'resid_trend_mean_se_wave.png'), wave_mean_resid_plot, width = 12, height = 8, dpi = 300)
        
        resid_parsmk_dt <- u_resid_mat_long %>%
            group_by(wave, parsmk) %>%
            summarise(
                mean_resid = mean(data_u),
                se = sd(data_u) / sqrt(n()),
                .groups = 'drop'
            ) %>%
            mutate(
                lower_ci = mean_resid - qt(0.975, df=n()-1) * se,
                upper_ci = mean_resid + qt(0.975, df=n()-1) * se
            )

            resid_parsmk_dt$parsmk <- ifelse(resid_parsmk_dt$parsmk == min(resid_parsmk_dt$parsmk), '0', '1')
            parsmk_mean_resid_plot <- ggplot(resid_parsmk_dt, aes(x = as.factor(wave), y = mean_resid, colour = as.factor(parsmk))) +
                                    geom_point(size = 3, position = position_dodge(width = 0.3)) +
                                    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, 
                                                    position = position_dodge(width = 0.3)) +
                                    labs(
                                        x = 'wave',
                                        y = 'mean data u-value'
                                    ) +
                                    scale_x_discrete(labels = c('1', '2', '3', '4', '5', '6')) +
                                    scale_color_manual(name = 'parsmk', values = c('0' = 'blue', '1' = 'red')) +
                                    theme(
                        axis.title.x = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        # axis.ticks.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )
            ggsave(file.path(model_results_dir, 'resid_trend_mean_se_parsmk.png'), parsmk_mean_resid_plot, width = 12, height = 8, dpi = 300)

            resid_sex_dt <- u_resid_mat_long %>%
                group_by(wave, sex) %>%
                summarise(
                    mean_resid = mean(data_u),
                    se = sd(data_u) / sqrt(n()),
                    .groups = 'drop'
                ) %>%
                mutate(
                    lower_ci = mean_resid - qt(0.975, df=n()-1) * se,
                    upper_ci = mean_resid + qt(0.975, df=n()-1) * se
                )

            resid_sex_dt$sex <- ifelse(resid_sex_dt$sex == min(resid_sex_dt$sex), '0', '1')
            sex_mean_resid_plot <- ggplot(resid_sex_dt, aes(x = as.factor(wave), y = mean_resid, colour = as.factor(sex))) +
                                    geom_point(size = 3, position = position_dodge(width = 0.3)) +
                                    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, 
                                                    position = position_dodge(width = 0.3)) +
                                    labs(
                                        x = 'wave',
                                        y = 'mean data u-value'
                                    ) +
                                    scale_x_discrete(labels = c('1', '2', '3', '4', '5', '6')) +
                                    scale_color_manual(name = 'sex', values = c('0' = 'blue', '1' = 'red')) +
                                    theme(
                        axis.title.x = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        # axis.ticks.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )
            ggsave(file.path(model_results_dir, 'resid_trend_mean_se_sex.png'), sex_mean_resid_plot, width = 12, height = 8, dpi = 300)

        # match u_resid with Z
        if (model == 'model_2') {
            # get long version of posterior Z 
            posterior_Z_long <- melt(data.table(posterior_Z), variable.name = 'newid', value.name = 'Z')
            posterior_Z_long[, newid := str_remove(newid, 'Z\\[')]
            posterior_Z_long[, newid := as.numeric(str_remove(newid, '\\]'))]
            posterior_Z_long[, posterior_index := rep(1:nrow(posterior_Z), times = ncol(posterior_Z))]

            # merge long Z and long resid
            resid_Z_long <- merge(u_resid_mat_long, posterior_Z_long, on = c('newid', 'posterior_index'))

            resid_Z_dt <- resid_Z_long %>%
                group_by(wave, Z) %>%
                summarise(
                    mean_resid = mean(data_u),
                    se = sd(data_u) / sqrt(n()),
                    .groups = 'drop'
                ) %>%
                mutate(
                    lower_ci = mean_resid - qt(0.975, df=n()-1) * se,
                    upper_ci = mean_resid + qt(0.975, df=n()-1) * se
                )

            Z_mean_resid_plot <- ggplot(resid_Z_dt, aes(x = as.factor(wave), y = mean_resid, colour = as.factor(Z))) +
                                    geom_point(size = 3, position = position_dodge(width = 0.3)) +
                                    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, 
                                                    position = position_dodge(width = 0.3)) +
                                    labs(
                                        x = 'wave',
                                        y = 'mean data u-value'
                                    ) +
                                    scale_x_discrete(labels = c('1', '2', '3', '4', '5', '6')) +
                                    scale_color_manual(name = 'Z', values = c('0' = 'blue', '1' = 'red')) +
                                    theme(
                        axis.title.x = element_text(size = 12),
                        axis.text.x = element_text(size = 12),
                        # axis.ticks.x = element_text(size = 12),
                        axis.title.y = element_text(size = 12),
                        axis.text.y = element_text(size = 12),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )
            ggsave(file.path(model_results_dir, 'resid_trend_mean_se_Z.png'), Z_mean_resid_plot, width = 12, height = 8, dpi = 300)   
        }

        # generate plots of residual ecdfs
        subset_size <- min(50, nrow(u_resid_mat))
        residual_ecdf_plot(u_resid_mat, subset_size, signif(resid_p_value, 2), model_results_dir)
        residual_tilted_ecdf_plot(u_resid_mat, subset_size, signif(resid_p_value, 2), model_results_dir, 0.05)
        if (model == 'model_0') {
            conditional_residual_tilted_ecdf_plot_vertical(u_resid_mat, subset_size, smoking_no_scale$wave, 'wave', model_results_dir, 0.1, wave_colors)
        } else {
            conditional_residual_tilted_ecdf_plot_vertical(u_resid_mat, subset_size, smoking_no_scale$wave, 'wave', model_results_dir, 0.1, wave_colors)
        }

        conditional_residual_tilted_ecdf_plot(u_resid_mat, subset_size, smoking_no_scale$sex, 'sex', model_results_dir, 0.05, covar_colors[1:2])
        conditional_residual_tilted_ecdf_plot(u_resid_mat, subset_size, smoking_no_scale$parsmk, 'parsmk', model_results_dir, 0.05, covar_colors[3:4])
        double_conditional_residual_tilted_ecdf_plot(u_resid_mat, subset_size, 
                                                     smoking_no_scale$parsmk, 'parsmk',
                                                     smoking_no_scale$wave, 'wave',
                                                     model_results_dir, 0.1, wave_colors
        )


        # generate plots of alpha ecdfs
        alpha_ecdf_plot_fn(u_alpha_mat, subset_size, signif(alpha_uniformity_p, 2), model_results_dir)
        alpha_tilted_ecdf_plot(u_alpha_mat, subset_size, signif(alpha_uniformity_p, 2), model_results_dir, 0.05)
        conditional_alpha_tilted_ecdf_plot(u_alpha_mat, subset_size, sex_vec_named, 'sex', model_results_dir, 0.05, covar_colors[1:2])
        conditional_alpha_tilted_ecdf_plot(u_alpha_mat, subset_size, parsmk_vec_named, 'parsmk', model_results_dir, 0.1, covar_colors[3:4])
        # if (model == 'model_0') {
        #     conditional_alpha_tilted_ecdf_plot(u_alpha_mat, subset_size, parsmk_vec_named, 'parsmk', model_results_dir, 0.1, covar_colors[3:4])
        # } else if (model == 'model_1') {
        #     conditional_alpha_tilted_ecdf_plot(u_alpha_mat, subset_size, parsmk_vec_named, 'parsmk', model_results_dir, 0.1, covar_colors[3:4])
        # } else {
        #     conditional_alpha_tilted_ecdf_plot(u_alpha_mat, subset_size, parsmk_vec_named, 'parsmk', model_results_dir, 0.05, covar_colors[3:4])
        # }
        

        # choose and note actual tests being looked at
        if (model %in% c('model_fixed_intercept', 'model_intercept_only')) {
            ps_of_interest <- c(wave_resid_p, sex_resid_p, parsmk_resid_p, sex_alpha_p, parsmk_alpha_p) 
            exchangeable_ps_of_interest <- c(exchangeable_wave_resid_p, exchangeable_sex_resid_p, exchangeable_parsmk_resid_p, exchangeable_sex_alpha_p, exchangeable_parsmk_alpha_p) 
            tests <- '(1-3): {wave, sex, parsmk} vs data u independence, (4-5): {sex, parsmk} vs alpha u independence'
        } else if (model == 'model_0') {
            ps_of_interest <- c(wave_resid_p, sex_resid_p, parsmk_resid_p, sex_alpha_p, parsmk_alpha_p) 
            exchangeable_ps_of_interest <- c(exchangeable_wave_resid_p, exchangeable_sex_resid_p, exchangeable_parsmk_resid_p, exchangeable_sex_alpha_p, exchangeable_parsmk_alpha_p) 
            tests <- '(1-3): {wave, sex, parsmk} vs data u independence, (4-5): {sex, parsmk} vs alpha u independence'
        } else if (model == 'model_1') {
            ps_of_interest <- c(alpha_uniformity_p, wave_x_parsmk_resid_p)
            exchangeable_ps_of_interest <- c(exchangeable_alpha_uniformity_p, exchangeable_wave_x_parsmk_resid_p)
            tests <- '(1): alpha u uniformity, (2): wave*parsmk vs data u independence'
        } else if (model == 'model_2') {
            ps_of_interest <- c(alpha_uniformity_p)
            exchangeable_ps_of_interest <- c(exchangeable_alpha_uniformity_p)
            tests <- 'alpha u uniformity'
        }
        adjusted_ps_of_interest <- p.adjust(ps_of_interest, method = 'holm')
        exchangeable_adjusted_ps_of_interest <- p.adjust(exchangeable_ps_of_interest, method = 'holm')

        combined_adjusted_ps_of_interest <- paste(adjusted_ps_of_interest, collapse = ', ')
        exchangeable_combined_adjusted_ps_of_interest <- paste(exchangeable_adjusted_ps_of_interest, collapse = ', ')

        output <- glue::glue('
                            Data Source: {data_source}\n
                            Model: {model}\n\n
                            Results Under Cauchy Assumptions\n
                            ----------------------------------\n
                            residual uniformity test: p = {resid_p_value}\n
                            alpha uniformity test: p = {alpha_uniformity_p}\n
                            parsmk/alpha independence test: p = {parsmk_alpha_p}\n
                            sex/alpha independence test: p = {sex_alpha_p}\n
                            parsmk/resid independence test: p = {parsmk_resid_p}\n
                            sex/resid independence test: p = {sex_resid_p}\n
                            wave/resid independence test: p = {wave_resid_p}\n
                            sex*wave/resid independence test: p = {wave_x_sex_resid_p}\n
                            parsmk*wave/resid independence test: p = {wave_x_parsmk_resid_p}\n
                            sex marginal test (linear regression): p = {sex_wave_marginal_p}\n
                            parsmk marginal test (linear regression): p = {parsmk_wave_marginal_p}\n
                            sex/wave interaction test (linear regression): p = {sex_wave_interaction_p}\n
                            parsmk/wave interaction test (linear regression): p = {parsmk_wave_interaction_p}\n
                            sex/resid by wave independence test: p = {combined_sex_wave_p}\n
                            parsmk/resid by wave independence test: p = {combined_parsmk_wave_p}\n
                            alpha/mean_resid independence test: p = {resid_alpha_indep_p}\n
                            tests of interest: {tests}: p = {combined_adjusted_ps_of_interest}\n\n
                            Results Under Exchangeability\n
                            -----------------------------\n
                            residual uniformity test: p = {exchangeable_resid_p_value}\n
                            alpha uniformity test: p = {exchangeable_alpha_uniformity_p}\n
                            parsmk/alpha independence test: p = {exchangeable_parsmk_alpha_p}\n
                            sex/alpha independence test: p = {exchangeable_sex_alpha_p}\n
                            parsmk/resid independence test: p = {exchangeable_parsmk_resid_p}\n
                            sex/resid independence test: p = {exchangeable_sex_resid_p}\n
                            wave/resid independence test: p = {exchangeable_wave_resid_p}\n
                            sex*wave/resid independence test: p = {exchangeable_wave_x_sex_resid_p}\n
                            parsmk*wave/resid independence test: p = {exchangeable_wave_x_parsmk_resid_p}\n
                            sex marginal test (linear regression): p = {exchangeable_sex_wave_marginal_p}\n
                            parsmk marginal test (linear regression): p = {exchangeable_parsmk_wave_marginal_p}\n
                            sex/wave interaction test (linear regression): p = {exchangeable_sex_wave_interaction_p}\n
                            parsmk/wave interaction test (linear regression): p = {exchangeable_parsmk_wave_interaction_p}\n
                            sex/resid by wave independence test: p = {exchangeable_combined_sex_wave_p}\n
                            parsmk/resid by wave independence test: p = {exchangeable_combined_parsmk_wave_p}\n
                            alpha/mean_resid independence test: p = {exchangeable_resid_alpha_indep_p}\n
                            tests of interest: {tests}: p = {exchangeable_combined_adjusted_ps_of_interest}
                            '
                            )
        writeLines(output, test_file)

        rm(u_resid_mat_long, u_resid_mat_indiv, combined_resid_mat_indiv, resid_alpha_indep_p_values)
    }
}

