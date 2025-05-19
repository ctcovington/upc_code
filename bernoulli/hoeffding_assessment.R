library(pacman)
p_load('Hmisc', 'data.table', 'ggplot2')

results_dir <- file.path('..', '..', 'results', 'bernoulli', 'hoeffding_assessment')
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

###########################################
### build null distribution for given n ###
###########################################
build_hoeffD_null <- function(n, n_sims) {
    null <- rep(0, n_sims)
    for (i in 1:n_sims) {
        null[i] <- hoeffd(runif(n), runif(n))$D[1,2]
    }
    return(null)
}

n_sims_null <- 100000
n_sims <- 10000

ns <- c(50, 100, 250, 500)
dt_list <- vector('list', length = length(ns))
iter <- 1

for (n in ns) {
    print(glue::glue('runs for n = {n}'))
    null <- build_hoeffD_null(n, n_sims_null)
    theoretical_res <- rep(0, n_sims)
    empirical_res <- rep(0, n_sims)
    for (i in 1:n_sims) {
        x <- runif(n)
        y <- runif(n)
        hoeff <- hoeffd(x,y)
        theoretical_res[i] <- hoeff$P[1,2]
        empirical_res[i] <- mean(hoeff$D[1,2] < null)
    }

    dt_list[[iter]] <- data.table('n' = n, 'Hmisc null' = theoretical_res, 'empirical null' = empirical_res)
    iter <- iter + 1
}

n_string_levels = paste0('n = ', ns)
res_dt <- rbindlist(dt_list)
res_dt_long <- melt(res_dt, id.vars = 'n', variable.name = 'null', value.name = 'value')
res_dt_long[, n_string := factor(paste0('n = ', n), levels = n_string_levels)]

plot <- ggplot(res_dt_long, aes(x = value)) + 
            geom_histogram(breaks = seq(0, 1, 0.01), fill = 'blue', colour = 'black', alpha = 1) +
            labs(x = 'p-value', y = 'count') +
            facet_grid(null ~ n_string) + 
            theme(
                        strip.text = element_text(size = 24),
                        axis.title.x = element_text(size = 24),
                        axis.text.x = element_text(size = 12),
                        axis.title.y = element_text(size = 24),
                        axis.text.y = element_text(size = 24),
                        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Light gray major gridlines
                        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.25), # Light gray minor gridlines
                        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Thin black border
                        panel.background = element_rect(fill = "white") # White background
                        )
ggsave( file.path(results_dir, 'hoeffding_null_comparison.png'), plot, height = 8, width = 16, dpi = 300 )
