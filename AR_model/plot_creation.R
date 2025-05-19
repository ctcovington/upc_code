plot_parameter_us <- function(res_dt, model_output_dir) {
    u_res_dt <- res_dt[, c('u1', 'u2')] # get u values for each parameter
    colnames(u_res_dt) <- c('phi', 'sigma') # rename u columns with proper parameter names
    plot_dt <- melt(u_res_dt)
    plot <- ggplot(plot_dt, aes(x = value)) + 
                stat_ecdf(position = 'identity', linewidth = 1, colour = 'red') + 
                geom_abline(intercept = 0, slope = 1, colour = 'black') + 
                coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
                facet_wrap(~variable, scale = 'free') + 
                labs(x = 'u')

    ggsave(file.path(model_output_dir, 'parameter_u_plots.png'), plot)

    # save individual plots to make combined plot
    plot_list <- lapply(sort(unique(plot_dt$variable)), function(i) {
                         ggplot(plot_dt[variable == i], aes(x = value)) + 
                                stat_ecdf(position = 'identity', linewidth = 1, colour = 'red', show.legend = FALSE) + 
                                geom_abline(intercept = 0, slope = 1, colour = 'black') + 
                                coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
                                ggtitle(paste0(i, ': u-values')) + 
                                theme(plot.title = element_text(size = 10))
                            }
                       )

    # look at proportion of time analyst would "reject" hypothesis of uniformity at varying levels 
    alpha_levels <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2)
    non_u_cols <- c('phi', 'sigma')
    rejection_rate_matrix <- matrix(0, nrow = length(alpha_levels), ncol = length(non_u_cols))
    for (i in 1:length(alpha_levels)) {
        for (j in 1:length(non_u_cols)) {
            alpha <- alpha_levels[i]
            colname <- non_u_cols[j]
            col <- u_res_dt[, ..colname]
            rejection_rate_matrix[i,j] <- mean( (col <= alpha/2) | (col >= 1-alpha/2) )
        }
    }

    rejection_rate_matrix <- cbind(alpha_levels, rejection_rate_matrix)
    rejection_rate_dt <- data.table(rejection_rate_matrix)
    colnames(rejection_rate_dt) <- c('alpha', non_u_cols)
    fwrite(rejection_rate_dt, file.path(model_output_dir, 'u_value_rejection_rates.csv'))

    rejection_rate_plot_dt <- melt(rejection_rate_dt, id.vars = 'alpha', variable.name = 'parameter', value.name = 'rejection_rate')
    rejection_rate_plot_dt[, parameter := fct_rev(parameter)]

    rejection_rate_plot <- ggplot(rejection_rate_plot_dt, aes(x = as.factor(alpha), y = parameter, fill = rejection_rate)) + 
                                geom_tile(color = 'black') + 
                                geom_text(aes(label = rejection_rate), color = 'black') +
                                scale_fill_gradient(name = 'Rejection rate', limits = c(0, 1), low = 'white', high = 'red') +
                                labs(
                                    # title = 'Rejection rate at various alpha levels',
                                     x = 'alpha',
                                     y = 'parameter')
    ggsave(file.path(model_output_dir, 'u_value_rejection_rate_plot.png'), rejection_rate_plot, width = 12, height = 8, dpi = 300)

    return(plot_list)
}

create_residual_plots <- function(u_resid_mat, y_mat, model_output_dir) {
    n_sims <- nrow(u_resid_mat)
    n <- ncol(u_resid_mat)

    resid_ks_ps <- rep(0, n_sims)
    indep_i_resid_ps <- rep(0, n_sims)
    # serial_corr_resid_ps <- rep(0, n_sims)
    lag_1_resid_ps <- rep(0, n_sims)
    lag_2_resid_ps <- rep(0, n_sims)
    # lag_3_resid_ps <- rep(0, n_sims)
    # lag_0_y_ps <- rep(0, n_sims)
    # lag_1_y_ps <- rep(0, n_sims)
    # lag_2_y_ps <- rep(0, n_sims)
    # lag_3_y_ps <- rep(0, n_sims)
    # lag_50_resid_ps <- rep(0, n_sims)
    # lag_50_y_ps <- rep(0, n_sims)

    if ( file.exists(file.path(model_output_dir, 'resid_pvals.csv')) ) {
        res_dt <- fread( file.path(model_output_dir, 'resid_pvals.csv') )
    } else {
        for (i in 1:n_sims) {
            # print( glue::glue('run {i} of {n_sims}') )
            resids <- u_resid_mat[i,]
            y <- as.numeric(y_mat[i,])

            # ks-test for residual us
            resid_ks <- ks.test(resids, punif)
            resid_ks_ps[i] <- resid_ks$p.value

            # p-values of tests of independence of index and residuals 
            indep_i_resid_ps[i] <- dcov.test(1:n, resids, R = 500)$p

            # p-values of tests of independence between lagged residuals
            lag_1_resid_ps[i] <- dcov.test(resids[1:(n-1)], resids[2:n], R = 100)$p
            lag_2_resid_ps[i] <- dcov.test(resids[1:(n-2)], resids[3:n], R = 100)$p
            # lag_3_resid_ps[i] <- dcov.test(resids[1:(n-3)], resids[4:n], R = 100)$p
            # lag_50_resid_ps[i] <- dcov.test(resids[1:(n-50)], resids[51:n], R = 100)$p

            # p-values of tests of independence between residuals and lagged y
            # lag_0_y_ps[i] <- dcov.test(y[1:n], resids[1:n], R = 100)$p
            # lag_1_y_ps[i] <- dcov.test(y[1:(n-1)], resids[2:n], R = 100)$p
            # lag_2_y_ps[i] <- dcov.test(y[1:(n-2)], resids[3:n], R = 100)$p
            # lag_3_y_ps[i] <- dcov.test(y[1:(n-3)], resids[4:n], R = 100)$p
            # lag_50_y_ps[i] <- dcov.test(y[1:(n-50)], resids[51:n], R = 100)$p

            # p-value of test of serial correlation in the residuals
            # TODO: look for a better test
            # try something like dcov.test(resids[1:498], resids[3:500], R = 100)$p
            # serial_corr_resid_ps[i] <- serialCorrelationTest(resids)$p.value

            # p-value of test of independence between 
        }
        res_dt <- data.table('KS' = resid_ks_ps, 
                            'index' = indep_i_resid_ps, 
                            'lag_1_resid' = lag_1_resid_ps,
                            'lag_2_resid' = lag_2_resid_ps
                            #  'lag_3_resid' = lag_3_resid_ps,
                            #  'lag_50_resid' = lag_50_resid_ps,
                            #  'lag_0_y' = lag_0_y_ps,
                            #  'lag_1_y' = lag_1_y_ps,
                            #  'lag_2_y' = lag_2_y_ps
                            #  'lag_3_y' = lag_3_y_ps
                            #  'lag_50_y' = lag_50_y_ps
                            #  'serial_corr' = serial_corr_resid_ps
                            )
    }
    res_dt <- res_dt[, c('KS', 'index', 'lag_1_resid', 'lag_2_resid')] # ensure only variables we want are needed
    plot_dt <- melt(res_dt, variable.name = 'test')

    # plot p-values for various tests
    resid_p_val_plot <- ggplot(plot_dt, aes(x = value, colour = test)) + 
                                    stat_ecdf(position = 'identity', linewidth = 1) + 
                                    geom_abline(intercept = 0, slope = 1, colour = 'black') + 
                                    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) + 
                                    scale_colour_manual(values = c('blue', 'red', 
                                                                '#2d8700', 
                                                                # '#905e00', 
                                                                '#53117c' 
                                                                # 'skyblue', 
                                                                # 'green', 
                                                                # 'orange', 
                                                                # 'purple')
                                                                )) +
                                    labs(
                                        # title = paste0('p-values from residual tests'),
                                         x = 'p-value')
    ggsave(file.path(model_output_dir, 'resid_tests.png'), resid_p_val_plot)
    fwrite(res_dt, file.path(model_output_dir, 'resid_pvals.csv'))

    # look at proportion of time analyst would "reject" hypothesis of uniformity at varying levels 
    alpha_levels <- c(0.01, 0.05, 0.1)
    # p_cols <- c('KS', 'indep', 'serial_corr')
    p_cols <- colnames(res_dt)
    rejection_rate_matrix <- matrix(0, nrow = length(alpha_levels), ncol = length(p_cols))
    for (i in 1:length(alpha_levels)) {
        for (j in 1:length(p_cols)) {
            alpha <- alpha_levels[i]
            colname <- p_cols[j]
            col <- res_dt[, ..colname]
            rejection_rate_matrix[i,j] <- mean( (col <= alpha/2) | (col >= 1-alpha/2) )
        }
    }

    rejection_rate_matrix <- cbind(alpha_levels, rejection_rate_matrix)
    rejection_rate_dt <- data.table(rejection_rate_matrix)
    colnames(rejection_rate_dt) <- c('alpha', p_cols)
    fwrite(rejection_rate_dt, file.path(model_output_dir, 'resid_test_rejection_rates.csv'))

    rejection_rate_plot_dt <- melt(rejection_rate_dt, id.vars = 'alpha', variable.name = 'test', value.name = 'rejection_rate')
    rejection_rate_plot_dt[, test := fct_rev(test)]

    rejection_rate_plot <- ggplot(rejection_rate_plot_dt, aes(x = as.factor(alpha), y = test, fill = rejection_rate)) + 
                                geom_tile(color = 'black') + 
                                geom_text(aes(label = rejection_rate), color = 'black') +
                                scale_fill_gradient(name = 'rejection rate', limits = c(0, 1), low = 'white', high = 'red') +
                                labs(
                                    # title = 'Rejection rate at various alpha levels',
                                     x = 'alpha',
                                     y = 'test')
    ggsave(file.path(model_output_dir, 'resid_test_rejection_rate_plot.png'), rejection_rate_plot)

    # create trajectory plots of residuals over index
    u_resid_dt <- data.table(u_resid_mat)
    u_resid_dt$sim <- 1:nrow(u_resid_dt)
    trajectory_dt <- melt(u_resid_dt, id.vars = 'sim')
    trajectory_dt$variable <- as.numeric( str_replace_all(trajectory_dt$variable, 'V', '') )
    names(trajectory_dt) <- c('sim', 'index', 'u')

    trajectory_plot <- ggplot(trajectory_dt, aes(x = index, y = u, group = as.factor(sim))) + 
                            stat_smooth(geom = 'line', alpha = 0.05)
    ggsave(file.path(model_output_dir, 'resid_trajectory_plot.png'), trajectory_plot, width = 12, height = 8, dpi = 300)

}

AR_p_val_dist <- function(p_val_dt, test_names, file_names, model_results_dir) {
    cols <- colnames(p_val_dt)
    for (i in 1:length(test_names)) {
        # make histogram
        hist_plot <- ggplot(p_val_dt, aes(x = p_val_dt[[cols[i]]])) + 
                    geom_histogram(breaks = seq(0, 1, 0.01), position = 'identity', colour = 'black', fill = 'blue') +
                    labs(x = test_names[i], y = 'count') +
                    theme(
                        axis.title.x = element_text(size = rel(3.0)),
                        axis.text.x = element_text(size = rel(3.0)),
                        axis.title.y = element_text(size = rel(3.0)),
                        axis.text.y = element_text(size = rel(3.0)),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )
        
        # make density plot
        dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = p_val_dt[[cols[i]]])) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
        dens_data <- dens_data[, c('x', 'density')]

        # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
        # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
        #                             dens_data,
        #                             data.frame(x = 1, density = 0))
        #                     )
        
        dens_plot <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', linewidth = 1) + 
                coord_cartesian(xlim = c(0, 1), ylim = c(0, NA)) +
                labs(x = test_names[i], y = 'density') + 
                theme(
                    axis.title.x = element_text(size = rel(3.0)),
                    axis.text.x = element_text(size = rel(3.0)),
                    axis.title.y = element_text(size = rel(3.0)),
                    axis.text.y = element_text(size = rel(3.0)),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
        
        ecdf_plot <- ggplot(p_val_dt, aes(x = p_val_dt[[cols[i]]])) + 
                    stat_ecdf(position = 'identity', linewidth = 1, alpha = 1, colour = 'blue', show.legend = FALSE) + 
                    geom_abline(intercept = 0, slope = 1, colour = 'black') + 
                    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
                    labs(x = test_names[i], y = 'empirical cdf') +
                    theme(
                        axis.title.x = element_text(size = rel(3.0)),
                        axis.text.x = element_text(size = rel(3.0)),
                        axis.title.y = element_text(size = rel(3.0)),
                        axis.text.y = element_text(size = rel(3.0)),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )

        ggsave(file.path(model_results_dir, glue::glue('{file_names[i]}_p_val_histogram.png')), hist_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, glue::glue('{file_names[i]}_p_val_density.png')), dens_plot, width = 12, height = 8, dpi = 300)
        ggsave(file.path(model_results_dir, glue::glue('{file_names[i]}_p_val_ecdf.png')), ecdf_plot, width = 12, height = 8, dpi = 300)
    }
}
