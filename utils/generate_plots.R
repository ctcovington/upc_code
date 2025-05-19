residual_trajectory_plot <- function(u_resid_mat, model_results_dir) {
    # create trajectory plots of residuals over index
    u_resid_dt <- data.table(u_resid_mat)
    u_resid_dt$sample <- 1:nrow(u_resid_dt)
    trajectory_dt <- melt(u_resid_dt, id.vars = 'sample')
    trajectory_dt$variable <- as.numeric( str_replace_all(trajectory_dt$variable, 'V', '') )
    names(trajectory_dt) <- c('sample', 'index', 'u')

    trajectory_plot <- ggplot(trajectory_dt, aes(x = index, y = u, group = as.factor(sample))) + 
                            # stat_smooth(geom = 'line', alpha = 0.05)
                            geom_point(alpha = 0.2) + 
                            labs(x = 'index', y = 'data u-value') + 
                            theme(axis.title.x = element_text(size = 36),
                                      axis.text.x = element_text(size = 36),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 36),
                                      axis.text.y = element_text(size = 36),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                      panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                      panel.background = element_rect(fill = "white") # White background
                                    #   axis.ticks.y = element_text(size = 36)
                                     )
                            # geom_vline(xintercept = 6, color = 'red', linetype = 'dashed') + 
                            # geom_vline(xintercept = 10, color = 'red', linetype = 'dashed') # place red lines at known outliers
    ggsave(file.path(model_results_dir, 'resid_trajectory_plot.png'), trajectory_plot, width = 12, height = 8, dpi = 300)
}

residual_ecdf_plot <- function(u_resid_mat, subset_size, p_value, model_results_dir) {
    u_resid_dt <- data.table(u_resid_mat)
    u_resid_dt$sample <- 1:nrow(u_resid_dt)
    if (subset_size > 1) {
        # get random sample
        u_resid_dt <- u_resid_dt[sample(.N, subset_size)]
    } else if (subset_size == 1) {
        # take first sample
        u_resid_dt <- u_resid_dt[1,]
    }
    u_resid_melt <- melt(u_resid_dt, id.vars = 'sample')[, c('sample', 'value')]

    # plot empirical cdf of residuals across DAP samples
    resid_ecdf_plot <- ggplot(u_resid_melt, aes(x = value, group = sample)) + 
                                stat_ecdf(position = 'identity', linewidth = 1, alpha = 0.1, colour = 'blue', show.legend = FALSE) + 
                                geom_abline(intercept = 0, slope = 1, colour = 'black') + 
                                coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
                                annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
                                         x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 13) +
                                labs(x = expression('U'['data']), y = 'empirical cdf') +  
                                theme(axis.title.x = element_text(size = 36),
                                    #   axis.text.x = element_text(size = 36),
                                      axis.text.x = element_text(size = 36),
                                    #   axis.ticks.x = element_blank(),
                                    #   axis.text.x = element_blank(),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 36),
                                      axis.text.y = element_text(size = 36),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                      panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                      panel.background = element_rect(fill = "white") # White background
                                     )
                                

    ggsave(file.path(model_results_dir, 'residual_ecdf_plot.png'), resid_ecdf_plot, width = 12, height = 8, dpi = 300)
}

alpha_ecdf_plot_fn <- function(u_alpha_mat, subset_size, p_value, model_results_dir) {
    u_alpha_dt <- data.table(u_alpha_mat)
    u_alpha_dt$sample <- 1:nrow(u_alpha_dt)
    if (subset_size > 1) {
        # get random sample
        u_alpha_dt <- u_alpha_dt[sample(.N, subset_size)]
    } else if (subset_size == 1) {
        # take first sample
        u_alpha_dt <- u_alpha_dt[1,]
    }
    u_alpha_melt <- melt(u_alpha_dt, id.vars = 'sample')[, c('sample', 'value')]

    # plot empirical cdf of alphas across DAP samples
    alpha_ecdf_plot <- ggplot(u_alpha_melt, aes(x = value, group = sample)) + 
                                stat_ecdf(position = 'identity', linewidth = 1, alpha = 0.1, colour = 'blue', show.legend = FALSE) + 
                                geom_abline(intercept = 0, slope = 1, colour = 'black') + 
                                coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
                                annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
                                         x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 13) +
                                labs(x = TeX('$U_{\\alpha}$'), y = 'empirical cdf') +  
                                theme(axis.title.x = element_text(size = 36),
                                    #   axis.text.x = element_text(size = 36),
                                      axis.text.x = element_text(size = 36),
                                    #   axis.ticks.x = element_blank(),
                                    #   axis.text.x = element_blank(),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 36),
                                      axis.text.y = element_text(size = 36),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                      panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                      panel.background = element_rect(fill = "white") # White background
                                     )
                                

    ggsave(file.path(model_results_dir, 'alpha_ecdf_plot.png'), alpha_ecdf_plot, width = 12, height = 8, dpi = 300)
}

one_posterior_sample_residual_ecdf_plot <- function(u_resid_mat, model_results_dir) {
    # u_resid_dt <- data.table(u_resid_mat[1,])
    # u_resid_dt$sample <- 1:nrow(u_resid_dt)
    # u_resid_melt <- melt(u_resid_dt, id.vars = 'sample')[, c('sample', 'value')]
    # u_resid_melt <- u_resid_melt[sample == min(u_resid_dt$sample),]

    u_resid_dt <- data.table(u_resid_mat[1,])
    # u_resid_dt$sample <- 1:nrow(u_resid_dt)
    u_resid_dt$sample <- 1
    u_resid_melt <- melt(u_resid_dt, id.vars = 'sample')
    u_resid_melt[, newid := str_remove(variable, 'V')]

    # plot empirical cdf of residuals across DAP samples
    resid_ecdf_plot <- ggplot(u_resid_melt, aes(x = value, group = sample)) + 
                                stat_ecdf(position = 'identity', linewidth = 2, alpha = 1, colour = 'blue', show.legend = FALSE) + 
                                geom_abline(intercept = 0, slope = 1, colour = 'black') + 
                                coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
                                labs(x = expression('U'['data']), y = 'empirical cdf') +  
                                theme(axis.title.x = element_text(size = 36),
                                    #   axis.text.x = element_text(size = 36),
                                      axis.text.x = element_text(size = 36),,
                                      axis.ticks.x = element_blank(),
                                    #   axis.text.x = element_blank(),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 36),
                                      axis.text.y = element_text(size = 36),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                      panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                      panel.background = element_rect(fill = "white") # White background
                                     )                           

    ggsave(file.path(model_results_dir, 'one_posterior_sample_residual_ecdf_plot.png'), resid_ecdf_plot, width = 12, height = 8, dpi = 300)
}

# TODO: write this to "tilt" the ecdf plot so that it's plotting F(x) - x
residual_tilted_ecdf_plot <- function(u_resid_mat, subset_size, p_value, model_results_dir, abs_y_max) {
    u_resid_dt <- data.table(u_resid_mat)
    u_resid_dt$sample <- 1:nrow(u_resid_dt)
    if (subset_size > 1) {
        # get random sample
        u_resid_dt <- u_resid_dt[sample(.N, subset_size)]
    } else if (subset_size == 1) {
        # take first sample
        u_resid_dt <- u_resid_dt[1,]
    }
    
    u_resid_melt <- melt(u_resid_dt, id.vars = 'sample')
    u_resid_melt[, newid := str_remove(variable, 'V')]

    u_resid_melt[, tilted_ecdf := ecdf(value)(value) - value, by = sample]
    u_resid_melt <- u_resid_melt[order(sample, tilted_ecdf)]
    u_resid_melt[, order_within_sample := seq_len(.N), by = sample]

    resid_tilted_ecdf_plot <- ggplot(u_resid_melt, aes(x = value, y = tilted_ecdf, group = sample)) + 
                                geom_line(alpha = 0.3, colour = 'blue') + 
                                geom_abline(intercept = 0, slope = 0, colour = 'black') + 
                                annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
                                         x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 13) +
                                coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
                                labs(x = expression('U'['data']), y = 'tilted empirical cdf') +  
                                guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                theme(axis.title.x = element_text(size = 36),
                                      axis.text.x = element_text(size = 36),,
                                      axis.title.y = element_text(size = 36),
                                      axis.text.y = element_text(size = 36),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                      panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                      panel.background = element_rect(fill = "white"), # White background
                                      plot.margin = margin(c(1,1,1,1), 'cm')
                                     )
    
    ggsave(file.path(model_results_dir, 'residual_tilted_ecdf_plot.png'), resid_tilted_ecdf_plot, width = 12, height = 8, dpi = 300)
}

alpha_tilted_ecdf_plot <- function(u_alpha_mat, subset_size, p_value, model_results_dir, abs_y_max) {
    u_alpha_dt <- data.table(u_alpha_mat)
    u_alpha_dt$sample <- 1:nrow(u_alpha_dt)
    if (subset_size > 1) {
        # get random sample
        u_alpha_dt <- u_alpha_dt[sample(.N, subset_size)]
    } else if (subset_size == 1) {
        # take first sample
        u_alpha_dt <- u_alpha_dt[1,]
    }
    
    u_alpha_melt <- melt(u_alpha_dt, id.vars = 'sample')
    u_alpha_melt[, newid := gsub("[^0-9.-]", "", variable)]

    u_alpha_melt[, tilted_ecdf := ecdf(value)(value) - value, by = sample]
    u_alpha_melt <- u_alpha_melt[order(sample, tilted_ecdf)]
    u_alpha_melt[, order_within_sample := seq_len(.N), by = sample]

    alpha_tilted_ecdf_plot <- ggplot(u_alpha_melt, aes(x = value, y = tilted_ecdf, group = sample)) + 
                                geom_line(alpha = 0.3, colour = 'blue') + 
                                geom_abline(intercept = 0, slope = 0, colour = 'black') + 
                                annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
                                         x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 13) +
                                coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
                                labs(x = TeX('$U_{\\alpha}$'), y = 'tilted empirical cdf') +  
                                guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                theme(axis.title.x = element_text(size = 36),
                                      axis.text.x = element_text(size = 36),,
                                      axis.title.y = element_text(size = 36),
                                      axis.text.y = element_text(size = 36),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                      panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                      panel.background = element_rect(fill = "white"), # White background
                                      plot.margin = margin(c(1,1,1,1), 'cm')
                                     )
    
    ggsave(file.path(model_results_dir, 'alpha_tilted_ecdf_plot.png'), alpha_tilted_ecdf_plot, width = 16, height = 6, dpi = 300)
}

one_posterior_sample_residual_tilted_ecdf_plot <- function(u_resid_mat, model_results_dir) {
    u_resid_dt <- data.table(u_resid_mat[1,])
    # u_resid_dt$sample <- 1:nrow(u_resid_dt)
    u_resid_dt$sample <- 1
    u_resid_melt <- melt(u_resid_dt, id.vars = 'sample')
    u_resid_melt[, newid := str_remove(variable, 'V')]
    # u_resid_melt <- u_resid_melt[sample == min(u_resid_dt$sample),]

    u_resid_melt[, tilted_ecdf := ecdf(value)(value) - value, by = sample]
    u_resid_melt <- u_resid_melt[order(sample, tilted_ecdf)]
    u_resid_melt[, order_within_sample := seq_len(.N), by = sample]

    resid_tilted_ecdf_plot <- ggplot(u_resid_melt, aes(x = value, y = tilted_ecdf, group = sample)) + 
                                geom_line(linewidth = 2, alpha = 1, colour = 'blue') + 
                                geom_abline(intercept = 0, slope = 0, colour = 'black') + 
                                coord_cartesian(xlim = c(0,1), ylim = c(-1,1)) +
                                labs(x = expression('U'['data']), y = 'tilted empirical cdf') +  
                                guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                theme(axis.title.x = element_text(size = 36),
                                      axis.text.x = element_text(size = 36),
                                    #   axis.ticks.x = element_blank(),
                                    #   axis.text.x = element_blank(),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 36),
                                      axis.text.y = element_text(size = 36),
                                    #   axis.ticks.y = element_blank(),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                      panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                      panel.background = element_rect(fill = "white"), # White background
                                      plot.margin = margin(c(1,1,1,1), 'cm')
                                     )
    
    ggsave(file.path(model_results_dir, 'one_posterior_sample_residual_tilted_ecdf_plot.png'), resid_tilted_ecdf_plot, width = 12, height = 8, dpi = 300)
}

conditional_residual_tilted_ecdf_plot <- function(u_resid_mat, subset_size, conditional_group, conditional_group_name, model_results_dir, abs_y_max, colors) {
    u_resid_dt <- data.table(u_resid_mat)
    u_resid_dt$sample <- 1:nrow(u_resid_dt)
    if (subset_size > 1) {
        # get random sample
        u_resid_dt <- u_resid_dt[sample(.N, subset_size)]
    } else if (subset_size == 1) {
        # take first sample
        u_resid_dt <- u_resid_dt[1,]
    }
    u_resid_melt <- melt(u_resid_dt, id.vars = 'sample')
    u_resid_melt[, newid := str_remove(variable, 'V')]

    # u_resid_melt[, tilted_ecdf := ecdf(value)(value) - value, by = sample]
    # u_resid_melt <- u_resid_melt[order(sample, tilted_ecdf)]
    # u_resid_melt[, order_within_sample := seq_len(.N), by = sample]
    u_resid_melt[, newid := as.numeric(newid)]

    conditional_group_dt <- data.table('newid' = 1:max( as.numeric(u_resid_melt$newid) ),
                                       'conditional_group' = as.factor(conditional_group)
                                      )
    
    u_resid_melt <- merge(u_resid_melt, conditional_group_dt, by = 'newid')
    u_resid_melt[, tilted_ecdf := ecdf(value)(value) - value, by = c('sample', 'conditional_group')]
    u_resid_melt <- u_resid_melt[order(sample, tilted_ecdf)]
    u_resid_melt[, order_within_sample := seq_len(.N), by = sample]

    resid_tilted_ecdf_plot <- ggplot(u_resid_melt, aes(x = value, 
                                                       y = tilted_ecdf, 
                                                       colour = conditional_group,
                                                       group = sample)) + 
                                geom_line(alpha = 0.3, size = 1) + 
                                scale_color_manual(values = colors) +
                                coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
                                geom_abline(intercept = 0, slope = 0, colour = 'black') + 
                                # annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
                                #          x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 8) +
                                labs(x = expression('U'['data']), 
                                     y = 'tilted empirical cdf',
                                     colour = conditional_group_name) +  
                                # facet_wrap(~conditional_group, nrow = 2) +
                                guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                theme(
                                    #   legend.position = 'none',
                                      axis.title.x = element_text(size = 24),
                                      axis.text.x = element_text(size = 24),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 24),
                                      axis.text.y = element_text(size = 24),
                                      axis.ticks.y = element_blank(),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                        panel.background = element_rect(fill = "white"), # White background
                                        plot.margin = margin(c(1,1,1,1), 'cm')
                                     ) 
    resid_tilted_ecdf_plot_facet <- ggplot(u_resid_melt, aes(x = value, 
                                                       y = tilted_ecdf, 
                                                       colour = conditional_group,
                                                       group = sample)) + 
                                geom_line(alpha = 0.3, size = 1) + 
                                scale_color_manual(values = colors) +
                                coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
                                geom_abline(intercept = 0, slope = 0, colour = 'black') + 
                                # annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
                                #          x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 8) +
                                labs(x = expression('U'['data']), 
                                     y = 'tilted empirical cdf',
                                     colour = conditional_group_name) +  
                                facet_wrap(~conditional_group, nrow = 2) +
                                guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                theme(
                                      legend.position = 'none',
                                      strip.text = element_text(size = 24),
                                      axis.title.x = element_text(size = 24),
                                      axis.text.x = element_text(size = 24),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 24),
                                      axis.text.y = element_text(size = 24),
                                      axis.ticks.y = element_blank(),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                        panel.background = element_rect(fill = "white"), # White background
                                        plot.margin = margin(c(1,1,1,1), 'cm')
                                     ) 
    
    ggsave(file.path(model_results_dir, glue('residual_tilted_ecdf_plot_by_{conditional_group_name}.png')), resid_tilted_ecdf_plot, width = 8, height = 6, dpi = 300)
    ggsave(file.path(model_results_dir, glue('residual_tilted_ecdf_plot_by_{conditional_group_name}_faceted.png')), resid_tilted_ecdf_plot_facet, width = 8, height = 6, dpi = 300)
}

conditional_residual_tilted_ecdf_plot_vertical <- function(u_resid_mat, subset_size, conditional_group, conditional_group_name, model_results_dir, abs_y_max, colors) {
    u_resid_dt <- data.table(u_resid_mat)
    u_resid_dt$sample <- 1:nrow(u_resid_dt)
    if (subset_size > 1) {
        # get random sample
        u_resid_dt <- u_resid_dt[sample(.N, subset_size)]
    } else if (subset_size == 1) {
        # take first sample
        u_resid_dt <- u_resid_dt[1,]
    }
    u_resid_melt <- melt(u_resid_dt, id.vars = 'sample')
    u_resid_melt[, newid := str_remove(variable, 'V')]
    
    # u_resid_melt[, tilted_ecdf := ecdf(value)(value) - value, by = sample]
    # u_resid_melt <- u_resid_melt[order(sample, tilted_ecdf)]
    # u_resid_melt[, order_within_sample := seq_len(.N), by = sample]
    u_resid_melt[, newid := as.numeric(newid)]
    
    conditional_group_dt <- data.table('newid' = 1:max( as.numeric(u_resid_melt$newid) ),
                                       'conditional_group' = as.factor(conditional_group)
    )
    
    u_resid_melt <- merge(u_resid_melt, conditional_group_dt, by = 'newid')
    u_resid_melt[, tilted_ecdf := ecdf(value)(value) - value, by = c('sample', 'conditional_group')]
    u_resid_melt <- u_resid_melt[order(sample, tilted_ecdf)]
    u_resid_melt[, order_within_sample := seq_len(.N), by = sample]
    
    resid_tilted_ecdf_plot <- ggplot(u_resid_melt, aes(x = value, 
                                                       y = tilted_ecdf, 
                                                       colour = conditional_group,
                                                       group = sample)) + 
        geom_line(alpha = 0.3, size = 1) + 
        scale_color_manual(values = colors) +
        coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
        geom_abline(intercept = 0, slope = 0, colour = 'black') + 
        # annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
        #          x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 8) +
        labs(x = expression('U'['data']), 
             y = 'tilted empirical cdf',
             colour = conditional_group_name) +  
        # facet_wrap(~conditional_group, nrow = 2) +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(
            #   legend.position = 'none',
            axis.title.x = element_text(size = 24),
            axis.text.x = element_text(size = 24),
            #   axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 24),
            axis.text.y = element_text(size = 24),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
            panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
            panel.background = element_rect(fill = "white"), # White background
            plot.margin = margin(c(1,1,1,1), 'cm')
        ) 
    resid_tilted_ecdf_plot_facet <- ggplot(u_resid_melt, aes(x = value, 
                                                             y = tilted_ecdf, 
                                                             colour = conditional_group,
                                                             group = sample)) + 
        geom_line(alpha = 0.3, size = 1) + 
        scale_color_manual(values = colors) +
        coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
        geom_abline(intercept = 0, slope = 0, colour = 'black') + 
        # annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
        #          x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 8) +
        labs(x = expression('U'['data']), 
             y = 'tilted empirical cdf',
             colour = conditional_group_name) +  
        facet_wrap(~conditional_group, nrow = length(unique(u_resid_melt$conditional_group))) +
        guides(colour = guide_legend(override.aes = list(alpha = 1))) +
        theme(
            legend.position = 'none',
            strip.text = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            axis.text.x = element_text(size = 18),
            #   axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
            panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
            panel.background = element_rect(fill = "white"), # White background
            plot.margin = margin(c(1,1,1,1), 'cm')
        ) 
    
    ggsave(file.path(model_results_dir, glue('residual_tilted_ecdf_plot_by_{conditional_group_name}.png')), resid_tilted_ecdf_plot, width = 12, height = 8, dpi = 300)
    ggsave(file.path(model_results_dir, glue('residual_tilted_ecdf_plot_by_{conditional_group_name}_faceted_vertical.png')), resid_tilted_ecdf_plot_facet, width = 8, height = 12, dpi = 300)
}


conditional_alpha_tilted_ecdf_plot <- function(u_alpha_mat, subset_size, conditional_group, conditional_group_name, model_results_dir, abs_y_max, colors) {
    u_alpha_dt <- data.table(u_alpha_mat)
    u_alpha_dt$sample <- 1:nrow(u_alpha_dt)
    if (subset_size > 1) {
        # get random sample
        u_alpha_dt <- u_alpha_dt[sample(.N, subset_size)]
    } else if (subset_size == 1) {
        # take first sample
        u_alpha_dt <- u_alpha_dt[1,]
    }
    u_alpha_melt <- melt(u_alpha_dt, id.vars = 'sample')
    u_alpha_melt[, newid := gsub("[^0-9.-]", "", variable)]

    # u_alpha_melt[, tilted_ecdf := ecdf(value)(value) - value, by = sample]
    # u_alpha_melt <- u_alpha_melt[order(sample, tilted_ecdf)]
    # u_alpha_melt[, order_within_sample := seq_len(.N), by = sample]
    u_alpha_melt[, newid := as.numeric(newid)]

    conditional_group_dt <- data.table('newid' = 1:max( as.numeric(u_alpha_melt$newid) ),
                                       'conditional_group' = as.factor(conditional_group)
                                      )
    
    u_alpha_melt <- merge(u_alpha_melt, conditional_group_dt, by = 'newid')
    u_alpha_melt[, tilted_ecdf := ecdf(value)(value) - value, by = c('sample', 'conditional_group')]
    u_alpha_melt <- u_alpha_melt[order(sample, tilted_ecdf)]
    u_alpha_melt[, order_within_sample := seq_len(.N), by = sample]

    alpha_tilted_ecdf_plot <- ggplot(u_alpha_melt, aes(x = value, 
                                                       y = tilted_ecdf, 
                                                       colour = conditional_group,
                                                       group = sample)) + 
                                geom_line(alpha = 0.3, size = 1) + 
                                coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
                                geom_abline(intercept = 0, slope = 0, colour = 'black') + 
                                # annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
                                #          x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 8) +
                                labs(x = TeX('$U_{\\alpha}$'), 
                                     y = 'tilted empirical cdf',
                                     colour = conditional_group_name) +  
                                # facet_wrap(~conditional_group, nrow = 2) +
                                scale_color_manual(values = colors) +
                                guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                theme(
                                    #   legend.position = 'none',
                                      axis.title.x = element_text(size = 24),
                                      axis.text.x = element_text(size = 24),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 24),
                                      axis.text.y = element_text(size = 24),
                                      axis.ticks.y = element_blank(),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                        panel.background = element_rect(fill = "white"), # White background
                                        plot.margin = margin(c(1,1,1,1), 'cm')
                                     ) 
    alpha_tilted_ecdf_plot_facet <- ggplot(u_alpha_melt, aes(x = value, 
                                                       y = tilted_ecdf, 
                                                       colour = conditional_group,
                                                       group = sample)) + 
                                geom_line(alpha = 0.3, size = 1) + 
                                coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
                                geom_abline(intercept = 0, slope = 0, colour = 'black') + 
                                scale_color_manual(values = colors) +
                                # annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
                                #          x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 8) +
                                labs(x = TeX('$U_{\\alpha}$'), 
                                     y = 'tilted empirical cdf',
                                     colour = conditional_group_name) +  
                                facet_wrap(~conditional_group, nrow = 2) +
                                guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                theme(
                                      legend.position = 'none',
                                      strip.text = element_text(size = 24),
                                      axis.title.x = element_text(size = 24),
                                      axis.text.x = element_text(size = 24),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 24),
                                      axis.text.y = element_text(size = 24),
                                      axis.ticks.y = element_blank(),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                        panel.background = element_rect(fill = "white"), # White background
                                        plot.margin = margin(c(1,1,1,1), 'cm')
                                     ) 
    
    ggsave(file.path(model_results_dir, glue('alpha_tilted_ecdf_plot_by_{conditional_group_name}.png')), alpha_tilted_ecdf_plot, width = 8, height = 6, dpi = 300)
    ggsave(file.path(model_results_dir, glue('alpha_tilted_ecdf_plot_by_{conditional_group_name}_faceted.png')), alpha_tilted_ecdf_plot_facet, width = 8, height = 6, dpi = 300)
}

double_conditional_residual_tilted_ecdf_plot <- function(u_resid_mat, subset_size, 
                                                         conditional_group_1, conditional_group_name_1,
                                                         conditional_group_2, conditional_group_name_2, 
                                                         model_results_dir, abs_y_max, colors) {
    u_resid_dt <- data.table(u_resid_mat)
    u_resid_dt$sample <- 1:nrow(u_resid_dt)
    if (subset_size > 1) {
        # get random sample
        u_resid_dt <- u_resid_dt[sample(.N, subset_size)]
    } else if (subset_size == 1) {
        # take first sample
        u_resid_dt <- u_resid_dt[1,]
    }
    u_resid_melt <- melt(u_resid_dt, id.vars = 'sample')
    u_resid_melt[, newid := str_remove(variable, 'V')]

    # u_resid_melt[, tilted_ecdf := ecdf(value)(value) - value, by = sample]
    # u_resid_melt <- u_resid_melt[order(sample, tilted_ecdf)]
    # u_resid_melt[, order_within_sample := seq_len(.N), by = sample]
    u_resid_melt[, newid := as.numeric(newid)]

    conditional_group_dt <- data.table('newid' = 1:max( as.numeric(u_resid_melt$newid) ),
                                       'conditional_group_1' = as.factor(conditional_group_1),
                                       'conditional_group_2' = as.factor(conditional_group_2)
                                      )
    
    u_resid_melt <- merge(u_resid_melt, conditional_group_dt, by = 'newid')
    u_resid_melt[, tilted_ecdf := ecdf(value)(value) - value, by = c('sample', 'conditional_group_1', 'conditional_group_2')]
    u_resid_melt <- u_resid_melt[order(sample, tilted_ecdf)]
    u_resid_melt[, order_within_sample := seq_len(.N), by = sample]

    resid_tilted_ecdf_plot <- ggplot(u_resid_melt, aes(x = value, 
                                                       y = tilted_ecdf, 
                                                       colour = conditional_group_2,
                                                       group = sample)) + 
                                geom_line(alpha = 0.3, size = 1) + 
                                scale_color_manual(values = colors) +
                                coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
                                geom_abline(intercept = 0, slope = 0, colour = 'black') + 
                                # annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
                                #          x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 8) +
                                labs(x = expression('U'['data']), 
                                     y = 'tilted empirical cdf',
                                     colour = conditional_group_name_2) +  
                                # facet_wrap(~conditional_group, nrow = 2) +
                                guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                theme(
                                    #   legend.position = 'none',
                                      axis.title.x = element_text(size = 18),
                                      axis.text.x = element_text(size = 18),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 18),
                                      axis.text.y = element_text(size = 18),
                                      axis.ticks.y = element_blank(),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                        panel.background = element_rect(fill = "white"), # White background
                                        plot.margin = margin(c(1,1,1,1), 'cm')
                                     )

    resid_tilted_ecdf_plot_facet <- ggplot(u_resid_melt, aes(x = value, 
                                                       y = tilted_ecdf, 
                                                       colour = conditional_group_2,
                                                       group = sample)) + 
                                geom_line(alpha = 0.3, size = 1) + 
                                scale_color_manual(values = colors) +
                                coord_cartesian(xlim = c(0,1), ylim = c(-abs_y_max,abs_y_max)) +
                                geom_abline(intercept = 0, slope = 0, colour = 'black') + 
                                # annotate(geom = 'text', label = glue::glue('p-value: {p_value}'), 
                                #          x = Inf, y = -Inf, hjust = 1, vjust = -0.5, size = 8) +
                                labs(x = expression('U'['data']), 
                                     y = 'tilted empirical cdf',
                                     colour = conditional_group_name_2) +  
                                facet_grid(conditional_group_2 ~ conditional_group_1) +
                                guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                                theme(
                                      legend.position = 'none',
                                      strip.text = element_text(size = 18),
                                    #   axis.title.x = element_text(size = 18),
                                    #   axis.text.x = element_text(size = 18),
                                    axis.title.x = element_blank(),
                                      axis.text.x = element_blank(),
                                    #   axis.ticks.x = element_blank(),
                                      axis.title.y = element_text(size = 18),
                                      axis.text.y = element_text(size = 18),
                                      axis.ticks.y = element_blank(),
                                      panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                        panel.background = element_rect(fill = "white"), # White background
                                        plot.margin = margin(c(1,1,1,1), 'cm')
                                     ) 
    
    ggsave(file.path(model_results_dir, glue('residual_tilted_ecdf_plot_by_{conditional_group_name_1}_x_{conditional_group_name_2}.png')), resid_tilted_ecdf_plot, width = 12, height = 8, dpi = 300)
    ggsave(file.path(model_results_dir, glue('residual_tilted_ecdf_plot_by_{conditional_group_name_1}_x_{conditional_group_name_2}_faceted.png')), resid_tilted_ecdf_plot_facet, width = 8, height = 12, dpi = 300)
}

parameter_hist_plots <- function(parameter_us, parameter_names, LaTeX_parameter_names, model_results_dir) {
    # plot histograms of posterior samples of parameters
    plot1 <- ggplot(parameter_us, aes(x = u1)) + 
                geom_histogram(position = 'identity', bins = 100, color = 'black', fill = 'blue') + 
                # geom_density(colour = 'blue') + 
                coord_cartesian(xlim = c(0,1)) +
                labs(x = TeX(glue('{LaTeX_parameter_names[1]} U-values')), y = 'count') + 
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
    ggsave( file.path(model_results_dir, glue::glue('posterior_samples_{parameter_names[1]}_histogram.png')), plot1, width = 12, height = 8, dpi = 300 )

    if (length(parameter_names) > 1) {
        plot2 <- ggplot(parameter_us, aes(x = u2)) + 
                    geom_histogram(position = 'identity', bins = 100, color = 'black', fill = 'blue') + 
                    # geom_density(colour = 'blue') + 
                    coord_cartesian(xlim = c(0,1)) +
                    labs(x = TeX(glue('{LaTeX_parameter_names[2]} U-values')), y = 'count') + 
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
        ggsave( file.path(model_results_dir, glue::glue('posterior_samples_{parameter_names[2]}_histogram.png')), plot2, width = 12, height = 8, dpi = 300 )
    }
}

parameter_density_plots <- function(parameter_us, parameter_names, LaTeX_parameter_names, model_results_dir) {
    # plot histograms of posterior samples of parameters
    dens_data <- ggplot_build(
                    ggplot(parameter_us, aes(x = u1)) +
                        stat_bin(aes(y = after_stat(density)),
                        # after_stat(density / max(density))),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    plot1 <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', linewidth = 1) + 
                # coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
                coord_cartesian(xlim = c(0, 1)) +
                labs(x = TeX(glue('{LaTeX_parameter_names[1]} U-values')), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave( file.path(model_results_dir, glue::glue('posterior_samples_{parameter_names[1]}_density.png')), plot1, width = 12, height = 8, dpi = 300 )

    # plot1 <- ggplot(parameter_us, aes(x = u1)) + 
    #             stat_bin(aes(y = after_stat(density / max(density))), position = 'identity',
    #                                                                     binwidth = 0.01, boundary = 0,
    #                                                                     alpha = 1,
    #                                                                     geom = "line",
    #                                                                     color = "blue",
    #                                                                     size = 2
    #                                                                 ) +
    #             coord_cartesian(xlim = c(0,1)) +
    #             labs(x = TeX(glue('{LaTeX_parameter_names[1]} U-values')), y = 'density') + 
    #             theme(
    #                 axis.title.x = element_text(size = 36),
    #                 axis.text.x = element_text(size = 36),
    #                 # axis.ticks.x = element_text(size = 36),
    #                 axis.title.y = element_text(size = 36),
    #                 axis.text.y = element_text(size = 36),
    #                 # axis.ticks.y = element_blank(),
    #                 panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
    #                     panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
    #                     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
    #                     panel.background = element_rect(fill = "white") # White background
    #                 )
    # ggsave( file.path(model_results_dir, glue::glue('posterior_samples_{parameter_names[1]}_density.png')), plot1, width = 12, height = 8, dpi = 300 )

    if (length(parameter_names) > 1) {
        dens_data <- ggplot_build(
                    ggplot(parameter_us, aes(x = u2)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
        dens_data <- dens_data[, c('x', 'density')]
        # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
        dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
                                    dens_data,
                                    data.frame(x = 1, density = 0))
                            )
        plot2 <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', linewidth = 1) + 
                coord_cartesian(xlim = c(0, 1)) +
                labs(x = TeX(glue('{LaTeX_parameter_names[1]} U-values')), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )

        ggsave( file.path(model_results_dir, glue::glue('posterior_samples_{parameter_names[2]}_density.png')), plot2, width = 12, height = 8, dpi = 300 )
        # plot2 <- ggplot(parameter_us, aes(x = u2)) + 
        #             stat_bin(aes(y = after_stat(density / max(density))), position = 'identity',
        #                                                                 binwidth = 0.05, boundary = 0,
        #                                                                 alpha = 1,
        #                                                                 geom = "line",
        #                                                                 color = "blue",
        #                                                                 size = 2
        #                                                             ) + 
        #             coord_cartesian(xlim = c(0,1)) +
        #             labs(x = TeX(glue('{LaTeX_parameter_names[2]} U-values')), y = 'count') + 
        #             theme(
        #                 axis.title.x = element_text(size = 36),
        #                 axis.text.x = element_text(size = 36),
        #                 # axis.ticks.x = element_text(size = 36),
        #                 axis.title.y = element_text(size = 36),
        #                 axis.text.y = element_text(size = 36),
        #                 axis.ticks.y = element_blank(),
        #                 panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
        #                 panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
        #                 panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
        #                 panel.background = element_rect(fill = "white") # White background
        #                 )
        # ggsave( file.path(model_results_dir, glue::glue('posterior_samples_{parameter_names[2]}_density.png')), plot2, width = 12, height = 8, dpi = 300 )
    }
}

min_value_histogram <- function(y_pred_min_dt, observed_min, prior_results_dir) {
    # Calculate the maximum y-value from the histogram
    hist_data <- ggplot_build(
    ggplot(y_pred_min_dt, aes(x = y_pred_mins)) +
        geom_histogram(position = 'identity', bins = 30, color = 'black', fill = 'white')
    )$data[[1]]

    max_y <- max(hist_data$count)  # Extract maximum count value

    plot <- ggplot(y_pred_min_dt, aes(x = y_pred_mins)) + 
                geom_histogram(position = 'identity', bins = 30, color = 'black', fill = 'blue') + 
                # geom_density(colour = 'blue') + 
                # labs(x = expression('T = min{Y'[1]*'...Y'[I]*'}'), y = 'count') +
                labs(x = 'T', y = 'count') +
                coord_cartesian(xlim = c(1.1*observed_min, 1.1*max(y_pred_min_dt$y_pred_mins))) +
                # coord_cartesian(ylim = c(0,max_y)) +
                theme(
                        axis.title.x = element_text(size = 36),
                        axis.text.x = element_text(size = 36),
                        # axis.ticks.x = element_blank(),
                        # axis.text.x = element_text(size = 36),
                        # axis.ticks.x = element_text(size = 36),
                        axis.title.y = element_text(size = 36),
                        axis.text.y = element_text(size = 36),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        ) + 
                geom_vline(aes(xintercept = observed_min), color = 'red', linetype = 'dashed', size = 2) 
                # + # Add red vertical line
                    # annotate(
                    #     'text', 
                    #     x = observed_min, 
                    #     y = max_y, 
                    #     label = "value of T on the observed data", 
                    #     color = 'red', 
                    #     vjust = 0.1,
                    #     hjust = -0.01,
                    #     size = 10
                    # )
        ggsave( file.path(prior_results_dir, 'posterior_min_value_histogram.png'), plot, width = 12, height = 8, dpi = 300 )
}

min_value_density <- function(y_pred_min_dt, observed_min, prior_results_dir) {
    dens_data <- ggplot_build(
                    ggplot(y_pred_min_dt, aes(x = y_pred_mins)) +
                        stat_bin(aes(y = after_stat(density / max(density))),
                                #  breaks = seq(0, 1, by = 0.01)
                                # binwidth = 0.5
                                bins = 30
                                )
                    )$data[[1]]
        dens_data <- dens_data[, c('x', 'ndensity')]
        # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
        # dens_data <- rbindlist( list(data.frame(x = 0, ndensity = 0), 
        #                             dens_data,
        #                             data.frame(x = 1, ndensity = 0))
        #                     )
    plot <- ggplot(dens_data, aes(x = x, y = ndensity)) + 
                geom_line(color = 'blue', linewidth = 1) + 
                labs(x = expression('T = min{Y'[1]*'...Y'[I]*'}'), y = 'density') +
                coord_cartesian(xlim = c(1.1*observed_min, 1.1*max(y_pred_min_dt$y_pred_mins))) +
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                ) + 
                geom_vline(aes(xintercept = observed_min), color = 'red', linetype = 'dashed', size = 2) + # Add red vertical line
                    annotate(
                        'text', 
                        x = observed_min, 
                        y = 1, 
                        label = "Minimum in Newcomb's data", 
                        color = 'red', 
                        vjust = 0,
                        hjust = -0.1,
                        size = 8
                    )
    ggsave( file.path(prior_results_dir, 'posterior_min_value_density.png'), plot, width = 12, height = 8, dpi = 300 )
}

approx_min_posterior_predictive_p <- function(mu_samples, sigma_samples, I, observed_min) {
    S <- length(mu_samples)
    approx_S <- rep(0, S)
    for (s in 1:S) {
        approx_S[s] <- 1 - ( 1 - pnorm(observed_min, mu_samples[s], sigma_samples[s]) )^I
    }
    approx <- mean(approx_S)

    return(approx)
}

Newcomb_p_value_histograms <- function(p_val_dt, dir, plot_min_posterior_predictive) {
    mu_p_hist <- ggplot(p_val_dt, aes(x = mu)) + 
                    geom_histogram(position = 'identity', binwidth = 0.01, boundary = 0, colour = 'black', fill = 'blue') +
                    labs(x = TeX('$p_{\\mu}$'), y = 'count') + 
                    coord_cartesian(xlim = c(0,1)) +
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
    ggsave(file.path(dir, 'mu_p_histogram.png'), mu_p_hist, width = 12, height = 8, dpi = 300)

    sigma_p_hist <- ggplot(p_val_dt, aes(x = sigma)) + 
                geom_histogram(position = 'identity', binwidth = 0.01, boundary = 0, colour = 'black', fill = 'blue') +
                labs(x = TeX('$p_{\\sigma}$'), y = 'count') + 
                coord_cartesian(xlim = c(0,1)) +
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
    ggsave(file.path(dir, 'sigma_p_histogram.png'), sigma_p_hist, width = 12, height = 8, dpi = 300)

    data_p_hist <- ggplot(p_val_dt, aes(x = data)) + 
                geom_histogram(position = 'identity', binwidth = 0.01, boundary = 0, colour = 'black', fill = 'blue') +
                labs(x = expression('p'['data,unif']), y = 'count') + 
                coord_cartesian(xlim = c(0,1)) +
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
    ggsave(file.path(dir, 'data_p_histogram.png'), data_p_hist, width = 12, height = 8, dpi = 300)

    if (plot_min_posterior_predictive == TRUE) {
        min_pp_p_hist <- ggplot(p_val_dt, aes(x = min_posterior_predictive)) + 
                geom_histogram(position = 'identity', binwidth = 0.02, boundary = 0, colour = 'black', fill = 'blue') +
                labs(x = expression('p'['T']), y = 'count') + 
                coord_cartesian(xlim = c(0,1)) +
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
        ggsave(file.path(dir, 'min_posterior_predictive_p_histogram.png'), min_pp_p_hist, width = 12, height = 8, dpi = 300)
    }
}

Newcomb_p_value_densities <- function(p_val_dt, dir, plot_min_posterior_predictive) {
    ### mu ###
    dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = mu)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    mu_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', size = 1) + 
                coord_cartesian(xlim = c(0, 1), ylim = c(0, NA)) +
                labs(x = TeX('$p_{\\mu}$'), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'mu_p_density.png'), mu_p_hist, width = 12, height = 8, dpi = 300)

    ### sigma ###
    dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = sigma)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    sigma_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', size = 1) + 
                coord_cartesian(xlim = c(0, 1), ylim = c(0, NA)) +
                labs(x = TeX('$p_{\\sigma}$'), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'sigma_p_density.png'), sigma_p_hist, width = 12, height = 8, dpi = 300)

    ### data ###
    dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = data)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    data_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', size = 1) + 
                coord_cartesian(xlim = c(0, 1), ylim = c(0, NA)) +
                labs(x = expression('p'['data,unif']), y = 'density') +
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'data_p_density.png'), data_p_hist, width = 12, height = 8, dpi = 300)

    ### min posterior predictive ###
    if (plot_min_posterior_predictive == TRUE) {
        dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = min_posterior_predictive)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.02)
                                )
                    )$data[[1]]
        dens_data <- dens_data[, c('x', 'density')]

        # # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
        # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
        #                             dens_data,
        #                             data.frame(x = 1, density = 0))
        #                     )
        
        min_pp_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                    geom_line(color = 'blue', size = 1) + 
                    coord_cartesian(xlim = c(0, 1)) +
                    labs(x = expression('p'['T']), y = 'density') +
                    theme(
                        axis.title.x = element_text(size = 36),
                        axis.text.x = element_text(size = 36),
                        axis.title.y = element_text(size = 36),
                        axis.text.y = element_text(size = 36),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                    )
        ggsave(file.path(dir, 'min_posterior_predictive_p_density.png'), min_pp_p_hist, width = 12, height = 8, dpi = 300)
    }

    # mu_p_hist <- ggplot(p_val_dt, aes(x = mu)) + 
    #                 stat_bin(aes(y = after_stat(density / max(density))), position = 'identity',
    #                                                                     binwidth = 0.05, boundary = 0,
    #                                                                     alpha = 1,
    #                                                                     geom = "line",
    #                                                                     color = "blue",
    #                                                                     size = 2
    #                                                                 ) +
    #                 labs(x = TeX('$p_{\\mu}$'), y = 'density') + 
    #                 coord_cartesian(xlim = c(0,1)) +
    #                 theme(
    #                     axis.title.x = element_text(size = 36),
    #                     axis.text.x = element_text(size = 36),
    #                     # axis.ticks.x = element_text(size = 36),
    #                     axis.title.y = element_text(size = 36),
    #                     axis.text.y = element_text(size = 36),
    #                     # axis.ticks.y = element_blank(),
    #                     panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
    #                     panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
    #                     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
    #                     panel.background = element_rect(fill = "white") # White background
    #                     )
    # ggsave(file.path(dir, 'mu_p_density.png'), mu_p_hist, width = 12, height = 8, dpi = 300)

    # sigma_p_hist <- ggplot(p_val_dt, aes(x = sigma)) + 
    #             # geom_histogram(position = 'identity', binwidth = 0.05, boundary = 0, colour = 'black', fill = 'white') +
    #             # geom_density(colour = 'blue') +
    #             stat_bin(aes(y = after_stat(density / max(density))), position = 'identity',
    #                                                                     binwidth = 0.05, boundary = 0,
    #                                                                     alpha = 1,
    #                                                                     geom = "line",
    #                                                                     color = "blue",
    #                                                                     size = 2
    #                                                                 ) +
    #             labs(x = TeX('$p_{\\sigma}$'), y = 'density') + 
    #             coord_cartesian(xlim = c(0,1)) +
    #             theme(
    #                 axis.title.x = element_text(size = 36),
    #                 axis.text.x = element_text(size = 36),
    #                 # axis.ticks.x = element_text(size = 36),
    #                 axis.title.y = element_text(size = 36),
    #                 axis.text.y = element_text(size = 36),
    #                 # axis.ticks.y = element_blank(),
    #                 panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
    #                 panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
    #                 panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
    #                 panel.background = element_rect(fill = "white") # White background
    #                 )
    # ggsave(file.path(dir, 'sigma_p_density.png'), sigma_p_hist, width = 12, height = 8, dpi = 300)

    # data_p_hist <- ggplot(p_val_dt, aes(x = data)) + 
    #             # geom_histogram(position = 'identity', binwidth = 0.05, boundary = 0, colour = 'black', fill = 'white') +
    #             # geom_density(colour = 'blue') + 
    #             stat_bin(aes(y = after_stat(density / max(density))), position = 'identity',
    #                                                                     binwidth = 0.05, boundary = 0,
    #                                                                     alpha = 1,
    #                                                                     geom = "line",
    #                                                                     color = "blue",
    #                                                                     size = 2
    #                                                                 ) +
    #             labs(x = expression('p'['data,unif']), y = 'density') + 
    #             coord_cartesian(xlim = c(0,1)) +
    #             theme(
    #                 axis.title.x = element_text(size = 36),
    #                 axis.text.x = element_text(size = 36),
    #                 # axis.ticks.x = element_text(size = 36),
    #                 axis.title.y = element_text(size = 36),
    #                 axis.text.y = element_text(size = 36),
    #                 # axis.ticks.y = element_blank(),
    #                 panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
    #                 panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
    #                 panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
    #                 panel.background = element_rect(fill = "white") # White background
    #                 )
    # ggsave(file.path(dir, 'data_p_density.png'), data_p_hist, width = 12, height = 8, dpi = 300)
}

one_posterior_sample_data_u_density <- function(u_resid_mat, dir) {
    u_dt <- data.table(u_resid_mat[1, ])
    dens_data <- ggplot_build(
                    ggplot(u_dt, aes(x = V1)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                             dens_data,
    #                             data.frame(x = 1, density = 0))
    #                     )
    
    u_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', size = 1) + 
                coord_cartesian(xlim = c(0, 1)) +
                labs(x = expression('U'['data']), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'one_posterior_sample_data_u_density.png'), u_hist, width = 12, height = 8, dpi = 300)
    

}

one_posterior_sample_data_u_histogram <- function(u_resid_mat, dir) {
    u_dt <- data.table(u_resid_mat[1, ])
    
    u_hist <- ggplot(u_dt, aes(x = V1)) + 
                geom_histogram(position = 'identity', 
                               color = 'black', fill = 'blue') + 
                coord_cartesian(xlim = c(0, 1)) +
                labs(x = expression('U'['data']), y = 'count') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'one_posterior_sample_data_u_histogram.png'), u_hist, width = 12, height = 8, dpi = 300)
    

}

Newcomb_f_mu_given_sigma_y <- function(mu, sigma, mu_0, sigma_0, y) {
    I <- length(y)
    tau <- sigma^2
    tau_0 <- sigma_0^2

    theta_2 <- 1 / ( tau_0 + I*tau ) # precision
    theta_1 <- (mu_0*tau_0 + tau*sum(y)) * theta_2 # mean
    # theta_1 <- ( sigma_0^(-2) * mu_0 + sigma^(-2) * sum(y) ) / ( sigma_0^(-2) + I*sigma^(-2) )
    # theta_2 <- 1 / sqrt( 1 / ( sigma_0^(-2) + I*sigma^(-2) ) )
    approx <- dnorm( mu, theta_1, 1/sqrt(theta_2) )
    return(approx)
}

Newcomb_f_p_mu_given_sigma_y <- function(p_mu, y, sigma, mu_0, sigma_0) {
    term_1 <- Newcomb_f_mu_given_sigma_y(mu_0 + sigma_0 * qnorm(p_mu/2), sigma, mu_0, sigma_0, y) * 
              (sigma_0/2) / dnorm( qnorm(p_mu/2) )
    term_2 <- Newcomb_f_mu_given_sigma_y(mu_0 + sigma_0 * qnorm(1 - p_mu/2), sigma, mu_0, sigma_0, y) * 
              (sigma_0/2) / dnorm( qnorm(1 - p_mu/2) )
    return(term_1 + term_2)
}

Newcomb_f_p_mu_given_y <- function(p_mu, y, sigma_samples, mu_0, sigma_0) {
    # p_mu: p-value for test of uniformity of mu at which we want to evaluate the density
    # y: length I vector of outcomes 
    # sigma_samples: length T vector of posterior samples of sigma
    # mu_0: prior mean for mu
    # sigma_0: prior sd for mu

    T <- length(sigma_samples)
    approx <- 0
    for (t in 1:T) {
        approx <- approx + Newcomb_f_p_mu_given_sigma_y(p_mu, y, sigma_samples[t], mu_0, sigma_0)
    }

    return(approx/T)
}

Newcomb_f_sigma_given_mu_y <- function(sigma, mu, y, alpha_0, beta_0) {
    tau <- sigma^(-2)
    n <- length(y)
    density <- dgamma(tau, alpha_0 + n/2, beta_0 + 1/2 * sum( (y-mu)^2 ))
    return(density)

    # if (sigma <= 0 || sigma > max_sigma) {
    #     return(0)
    # } else {
    #     I <- length(y)
    #     alpha <- (I - 1) / 2
    #     beta <- 0.5 * sum((y - mu)^2)
        
    #     # Logarithmic numerator
    #     log_num <- -I * log(sigma) - beta / sigma^2
        
    #     # Logarithmic denominator
    #     log_denom <- log(1/2) + lgamma(alpha) - alpha * log(beta) + 
    #                  pinvgamma(max_sigma^2, shape = alpha, scale = beta, log.p = TRUE)
        
    #     # Return density
    #     return(exp(log_num - log_denom))
    # }
}

Newcomb_f_p_sigma_given_mu_y <- function(p_sigma, mu, y, alpha_0, beta_0) {
    # NOTE: this assumes the Normal-Gamma model with alpha_0 = beta_0 = 1
    if (alpha_0 != 1 || beta_0 != 1) {
        return(NA)
    } else {
        t1_sigma <- qgamma(p_sigma/2, shape = alpha_0, rate = beta_0)
        t2_sigma <- qgamma(1-p_sigma/2, shape = alpha_0, rate = beta_0)

        term_1 <- log(1 / (2-p_sigma)) + log( Newcomb_f_sigma_given_mu_y(t1_sigma, mu, y, alpha_0, beta_0) )
        term_2 <- log(1 / p_sigma) + log( Newcomb_f_sigma_given_mu_y(t2_sigma, mu, y, alpha_0, beta_0) )
    
        log_sum <- log(exp(term_1) + exp(term_2))  # Avoid direct addition to prevent overflow
        return(exp(log_sum))
    }

    # half_max_sigma <- max_sigma / 2
    
    # Calculate the two terms

    # term_1 <- log(half_max_sigma) + 
    #           log(Newcomb_f_sigma_given_mu_y(half_max_sigma * p_sigma, mu, y, max_sigma))
    # term_2 <- log(half_max_sigma) + 
    #           log(Newcomb_f_sigma_given_mu_y(max_sigma - half_max_sigma * p_sigma, mu, y, max_sigma))
    
    # Combine terms in a numerically stable way
    log_sum <- log(exp(term_1) + exp(term_2))  # Avoid direct addition to prevent overflow
    return(exp(log_sum))
}

Newcomb_f_p_sigma_given_y <- function(p_sigma, mu_samples, y, alpha_0, beta_0) {
    T <- length(mu_samples)
    log_approx <- -Inf  # Initialize in log-space
    
    # Accumulate log probabilities
    for (t in 1:T) {
        log_approx <- log(exp(log_approx) + 
                          Newcomb_f_p_sigma_given_mu_y(p_sigma, mu_samples[t], y, alpha_0, beta_0))
    }
    
    # Normalize by the number of samples
    log_approx <- log_approx - log(T)
    return(exp(log_approx))  # Convert back to probability
}

inv.erf <- function (x) qnorm((1 + x)/2)/sqrt(2)

Newcomb_NormalGamma_f_mu_given_lambda_y <- function(mu, y, lambda, mu_0, kappa_0) {
    # NOTE: seems correct, but not for tiny lambda?
    n <- length(y)
    T <- length(lambda)

    mu_n <- (kappa_0*mu_0 + sum(y)) / (kappa_0 + n)
    kappa_n <- kappa_0 + n

    approx <- 0

    for (i in 1:T) {
        # print(glue::glue('run {i} or {T}'))
        lamb <- lambda[i]
        density <- dnorm(mu, mean = mu_n, sd = (kappa_n * lamb)^(-1/2) )
        approx <- approx + density
    }

    return(approx / T)
}

Newcomb_NormalGamma_f_p_mu_given_lambda_y <- function(p_mu, y, lambda, mu_0, kappa_0) {
    # NOTE/TODO: incorrect

    t1_a <- mu_0 + (kappa_0 * lambda)^(-1/2) * sqrt(2) * inv.erf(p_mu - 1)
    t1_b <- abs( (kappa_0 * lambda)^(-1/2) * sqrt(pi/2) * exp(inv.erf(p_mu - 1)^2) )

    t2_a <- mu_0 + (kappa_0 * lambda)^(-1/2) * sqrt(2) * inv.erf(1 - p_mu)
    t2_b <- abs( -(kappa_0 * lambda)^(-1/2) * sqrt(pi/2) * exp(inv.erf(1 - p_mu)^2) )

    density <- ( Newcomb_NormalGamma_f_mu_given_lambda_y(t1_a, y, lambda, mu_0, kappa_0)*t1_b + 
                 Newcomb_NormalGamma_f_mu_given_lambda_y(t2_a, y, lambda, mu_0, kappa_0)*t2_b
               )

    return(density)
}

Newcomb_NormalGamma_f_p_mu_given_y <- function(p_mu, y, lambda_samples, mu_0, kappa_0) {
    # p_mu: p-value for test of uniformity of mu at which we want to evaluate the density
    # y: length I vector of outcomes 
    # lambda_samples: length T vector of posterior samples of lambda (precision)
    # mu_0, kappa_0: first two prior parameters for Normal-Gamma model 

    T <- length(lambda_samples)
    approx <- 0
    for (t in 1:T) {
        approx <- approx + Newcomb_NormalGamma_f_p_mu_given_lambda_y(p_mu, y, lambda_samples[t], mu_0, kappa_0)
    }

    return(approx/T)
}

Newcomb_NormalGamma_f_p_lambda_given_y <- function(p_lambda, alpha_I, beta_I) {
    t1 <- dgamma( -log(1 - p_lambda/2), shape = alpha_I, rate = beta_I )
    t2 <- dgamma( -log(p_lambda/2), shape = alpha_I, rate = beta_I )

    return( t1*(1 / (2-p_lambda)) + t2*(1/p_lambda) )
}

approximated_Newcomb_p_value_densities <- function(y, mu_samples, sigma_samples, mu_0, sigma_0, max_sigma, dir, plot_boolean) {
    # p_mu
    x_vals <- seq(0, 1, by = 0.01)
    p_mu_density <- sapply(x_vals, FUN = Newcomb_f_p_mu_given_y, y = y, sigma_samples = sigma_samples, mu_0 = mu_0, sigma_0 = sigma_0)
    p_mu_density[1] <- 0; p_mu_density[length(x_vals)] <- 0
    p_mu_dt <- data.table('p_mu' = x_vals, 'density' = p_mu_density)
    if (plot_boolean == TRUE) {
        p_mu_density_plot <- ggplot(p_mu_dt, aes(x = p_mu, y = density)) + 
                            geom_line(color = 'blue', size = 1) + 
                                    coord_cartesian(xlim = c(0, 1)) +
                                    labs(x = TeX('$p_{\\mu}$'), y = 'density') + 
                                    theme(
                                        axis.title.x = element_text(size = 36),
                                        axis.text.x = element_text(size = 36),
                                        axis.title.y = element_text(size = 36),
                                        axis.text.y = element_text(size = 36),
                                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                        panel.background = element_rect(fill = "white") # White background
                                    )
        ggsave(file.path(dir, 'single_simulation_mu_p_density_analytic.png'), p_mu_density_plot, width = 12, height = 8, dpi = 300)
    }

    # p_sigma
    x_vals <- seq(0, 1, by = 0.01)
    p_sigma_density <- sapply(x_vals, FUN = Newcomb_f_p_sigma_given_y, y = y, mu_samples = mu_samples, max_sigma = max_sigma)
    p_sigma_density[1] <- 0; p_sigma_density[length(x_vals)] <- 0
    p_sigma_dt <- data.table('p_sigma' = x_vals, 'density' = p_sigma_density)
    if (plot_boolean == TRUE) {
        p_sigma_density_plot <- ggplot(p_sigma_dt, aes(x = p_sigma, y = density)) + 
                            geom_line(color = 'blue', size = 1) + 
                                    coord_cartesian(xlim = c(0, 1)) +
                                    labs(x = TeX('$p_{\\sigma}$'), y = 'density') + 
                                    theme(
                                        axis.title.x = element_text(size = 36),
                                        axis.text.x = element_text(size = 36),
                                        axis.title.y = element_text(size = 36),
                                        axis.text.y = element_text(size = 36),
                                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                                        panel.background = element_rect(fill = "white") # White background
                                    )
    ggsave(file.path(dir, 'single_simulation_sigma_p_density_analytic.png'), p_sigma_density_plot, width = 12, height = 8, dpi = 300)
    }

    return( list(p_mu_dt, p_sigma_dt) )
}

bernoulli_switch_histogram <- function(y_pred_switch_dt, true_switches, prior_results_dir) {
    # Calculate the maximum y-value from the histogram
    hist_data <- ggplot_build(
    ggplot(y_pred_switch_dt, aes(x = y_pred_switches)) +
        geom_histogram(position = 'identity', bins = 30, breaks = seq(true_switches, 1.1*max(y_pred_switch_dt$y_pred_switches), length.out = 29), color = 'black', fill = 'white')
    )$data[[1]]

    max_y <- max(hist_data$count)  # Extract maximum count value

    plot <- ggplot(y_pred_switch_dt, aes(x = y_pred_switches)) + 
                geom_histogram(position = 'identity', color = 'black', fill = 'blue', 
                               bins = 30,
                               breaks = seq(true_switches, 1.1*max(y_pred_switch_dt$y_pred_switches), length.out = 29)
                               ) + 
                # labs(x = expression(T == sum(1 * (Y[i] != Y[i+1]), i==1, I-1)), y = 'count') +
                labs(x = 'T', y = 'count') +
                coord_cartesian(ylim = c(0, 1.08*max_y)) + 
                # coord_cartesian(xlim = c(true_switches, 1.1*max(y_pred_switch_dt$y_pred_switches)),
                                # ylim = c(0, 1.1*max_y)) +
                theme(
                        axis.title.x = element_text(size = 24),
                        axis.text.x = element_text(size = 24),
                        # axis.ticks.x = element_blank(),
                        # axis.text.x = element_text(size = 24),
                        # axis.ticks.x = element_text(size = 24),
                        axis.title.y = element_text(size = 24),
                        axis.text.y = element_text(size = 24),
                        # axis.ticks.y = element_blank(),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        ) + 
                geom_vline(xintercept = true_switches, color = 'red', linetype = 'dashed', size = 2) 
                # + # Add red vertical line
                    # annotate(
                    #     'text', 
                    #     x = true_switches, 
                    #     y = 1.05*max_y, 
                    #     label = "value of T on the observed data", 
                    #     color = 'red', 
                    #     vjust = 0.1,
                    #     hjust = -0.01,
                    #     size = 10
                    # )
        ggsave( file.path(prior_results_dir, 'posterior_switch_histogram.png'), plot, width = 12, height = 8, dpi = 300 )
}

bernoulli_p_value_densities_single_posterior <- function(p_val_dt, y, beta_a, beta_b, dir, plot_switch) {
    if (plot_switch == TRUE) {
        ### switch p-value density ###
        dens_data <- ggplot_build(
                        ggplot(p_val_dt, aes(x = switch)) +
                            stat_bin(aes(y = after_stat(density)),
                                    breaks = seq(0, 1, by = 0.01)
                                    )
                        )$data[[1]]
        dens_data <- dens_data[, c('x', 'density')]

        # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
        # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
        #                              dens_data,
        #                              data.frame(x = 1, density = 0))
        #                       )
        
        switch_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                    geom_line(color = 'blue', size = 1) + 
                    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.1*max(dens_data$density))) +
                    labs(x = expression('p'['T']), y = 'density') + 
                    theme(
                        axis.title.x = element_text(size = 24),
                        axis.text.x = element_text(size = 24),
                        axis.title.y = element_text(size = 24),
                        axis.text.y = element_text(size = 24),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                    )
        ggsave(file.path(dir, 'switch_p_density.png'), switch_p_hist, width = 12, height = 8, dpi = 300)
    }

    # ### uniform switch p-value density (this is just a check -- it should always be roughly uniform) ###
    # dens_data <- ggplot_build(
    #                 ggplot(p_val_dt, aes(x = unif_switch)) +
    #                     stat_bin(aes(y = after_stat(density)),
    #                              breaks = seq(0, 1, by = 0.01)
    #                             )
    #                 )$data[[1]]
    # dens_data <- dens_data[, c('x', 'density')]
    # # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    # unif_switch_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
    #             geom_line(color = 'blue', size = 1) + 
    #             coord_cartesian(xlim = c(0, 1)) +
    #             labs(x = expression('p'['T']), y = 'density') + 
    #             theme(
    #                 axis.title.x = element_text(size = 36),
    #                 axis.text.x = element_text(size = 36),
    #                 axis.title.y = element_text(size = 36),
    #                 axis.text.y = element_text(size = 36),
    #                 panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
    #                 panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
    #                 panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
    #                 panel.background = element_rect(fill = "white") # White background
    #             )
    # ggsave(file.path(dir, 'unif_switch_p_density.png'), unif_switch_p_hist, width = 12, height = 8, dpi = 300)

    ### theta p-value density ###
    # p_theta_density <- sapply(seq(0, 1, by = 0.01), FUN = bernoulli_f_p_theta_given_y, y = y, beta_a = beta_a, beta_b = beta_b)
    # dens_data <- data.table('x' = seq(0, 1, 0.01), 'density' = p_theta_density)
    
    dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = theta)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    theta_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', size = 1) + 
                scale_x_continuous(limits = c(0,1)) +
                scale_y_continuous(limits = c(0, NA)) +
                # scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
                labs(x = TeX('$p_{\\theta}$'), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'theta_p_density_single_posterior.png'), theta_p_hist, width = 12, height = 8, dpi = 300)

    ### data uniformity p-value density ###
    dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = data)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    unif_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', size = 1) + 
                scale_x_continuous(limits = c(0,1)) +
                scale_y_continuous(limits = c(0, NA)) +
                # scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
                labs(x = expression('p'['data,unif']), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'data_unif_p_density_single_posterior.png'), unif_p_hist, width = 12, height = 8, dpi = 300)

    ### data independence p-value density ###
    dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = indep)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    independence_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', size = 1) + 
                scale_x_continuous(limits = c(0,1)) +
                scale_y_continuous(limits = c(0, NA)) +
                labs(x = expression('p'['data,indep']), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'data_independence_p_density_single_posterior.png'), independence_p_hist, width = 12, height = 8, dpi = 300)
}

bernoulli_p_value_densities_sims <- function(p_val_dt, dir, plot_switch) {
    if (plot_switch == TRUE) {
        ### switch p-value density ###
        dens_data <- ggplot_build(
                        ggplot(p_val_dt, aes(x = switch)) +
                            stat_bin(aes(y = after_stat(density)),
                                    breaks = seq(0, 1, by = 0.02)
                                    )
                        )$data[[1]]
        dens_data <- dens_data[, c('x', 'density')]

        # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
        # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
        #                              dens_data,
        #                              data.frame(x = 1, density = 0))
        #                       )
        
        switch_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                    geom_line(color = 'blue', size = 1) + 
                    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.1*max(dens_data$density))) +
                    labs(x = expression('p'['T']), y = 'density') + 
                    theme(
                        axis.title.x = element_text(size = 24),
                        axis.text.x = element_text(size = 24),
                        axis.title.y = element_text(size = 24),
                        axis.text.y = element_text(size = 24),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                    )
        ggsave(file.path(dir, 'switch_p_density.png'), switch_p_hist, width = 12, height = 8, dpi = 300)
    }

    ### theta p-value density ###
    # p_theta_density <- sapply(seq(0, 1, by = 0.01), FUN = bernoulli_f_p_theta_given_y, y = y)
    # dens_data <- data.table('x' = seq(0, 1, 0.01), 'density' = p_theta_density)

    dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = theta)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    theta_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', size = 1) + 
                scale_x_continuous(limits = c(0,1)) +
                scale_y_continuous(limits = c(0, NA)) +
                # scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
                labs(x = TeX('$p_{\\theta}$'), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'theta_p_density.png'), theta_p_hist, width = 12, height = 8, dpi = 300)

    ### data uniformity p-value density ###
    dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = data)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    data_uniformity_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', size = 1) + 
                scale_x_continuous(limits = c(0,1)) +
                scale_y_continuous(limits = c(0, NA)) +
                # scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
                labs(x = expression('p'['data,unif']), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'data_uniformity_p_density.png'), data_uniformity_p_hist, width = 12, height = 8, dpi = 300)

    ### data independence p-value density ###
    dens_data <- ggplot_build(
                    ggplot(p_val_dt, aes(x = indep)) +
                        stat_bin(aes(y = after_stat(density)),
                                 breaks = seq(0, 1, by = 0.01)
                                )
                    )$data[[1]]
    dens_data <- dens_data[, c('x', 'density')]

    # Add a row to ensure the curve starts at (0, 0) and ends at (1, 0)
    # dens_data <- rbindlist( list(data.frame(x = 0, density = 0), 
    #                              dens_data,
    #                              data.frame(x = 1, density = 0))
    #                       )
    
    data_independence_p_hist <- ggplot(dens_data, aes(x = x, y = density)) + 
                geom_line(color = 'blue', size = 1) + 
                scale_x_continuous(limits = c(0,1)) +
                scale_y_continuous(limits = c(0, NA)) +
                # scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
                labs(x = expression('p'['data,unif']), y = 'density') + 
                theme(
                    axis.title.x = element_text(size = 36),
                    axis.text.x = element_text(size = 36),
                    axis.title.y = element_text(size = 36),
                    axis.text.y = element_text(size = 36),
                    panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                    panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                    panel.background = element_rect(fill = "white") # White background
                )
    ggsave(file.path(dir, 'data_independence_p_density.png'), data_independence_p_hist, width = 12, height = 8, dpi = 300)
}

bernoulli_f_p_theta_given_y <- function(p_theta, y, beta_a, beta_b) {
    s <- sum(y)
    n <- length(y)
    t1 <- 1/2 * dbeta(p_theta/2, beta_a+s, beta_b+n-s)
    t2 <- 1/2 * dbeta(1-p_theta/2, beta_a+s, beta_b+n-s)

    return(t1+t2)
}