################################ Load libraries ################################
library(tidyverse)
library(Cairo)
library(gridExtra)

############################### Define functions ###############################
expected_divergence <- function(ne, migr) {
  
  # Calculate expected level of divergence given the effective population size
  # (ne) and the migration rate (migr).
  #
  # Parameters:
  # ne    : haploid effective population size; divide by two for diploid ne
  # migr  : migration rate
  
  res <- 1 / (1 + 4 * ne * migr)
  return(res)
}

plot_sims <- function(df, exp_div) {
  
  # Plot error bars showing mean and 95% confidence interval of divergence
  # estimates (on y-axis) versus sample size (on x-axis).
  #
  # Parameters:
  # df : data frame with columns containing mean, and 0.025 and 0.975 quantiles
  #     of divergence estimates.
  
  # define plot theme
  plot_theme <- theme(panel.background = element_rect(fill = "white", 
                                                      colour = "grey50"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.justification = c(1,1),
                      legend.position = c(1,1.05),
                      legend.title = element_blank(),
                      legend.key.size = unit(1, 'lines'),
                      legend.spacing.x = unit(0.5, 'lines'),
                      legend.key = element_blank(),
                      legend.background = element_blank())
  
  # create plot with ggplot2
  g <- ggplot(df, aes(x = factor(sample_size), y = mean)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.5) +
    geom_hline(yintercept = exp_div, linetype = 2) +
    xlab('Sample size') +
    ylim(0, 0.10) +
    plot_theme
  
  # return plot
  return(g)
}

##################################### Main ##################################### 

# set working dir to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read data
df <- readr::read_csv("sim_r1000/simulation_table.csv")

# summarize results; Hudson's estimator (Kst; allele counts)
hudson.fst <- df %>% 
  group_by(sample_size) %>%
    summarise('snps' = mean(n_snps),
              'mean' = mean(hudson_fst),
              'std' = sd(hudson_fst),
              'q025' = quantile(hudson_fst, probs = 0.025),
              'q975' = quantile(hudson_fst, probs = 0.975))

# summarize results; Weir & Cockerham's theta (Fst; allele frequencies)
wc.fst <- df %>% 
  group_by(sample_size) %>%
  summarise('snps' = mean(n_snps),
            'mean' = mean(wc_fst),
            'std' = sd(wc_fst),
            'q025' = quantile(wc_fst, probs = 0.025),
            'q975' = quantile(wc_fst, probs = 0.975))

# export summarized results
readr::write_csv(hudson.fst, "simulation_results_fst_huds.csv")
readr::write_csv(wc.fst, "simulation_results_fst_wc.csv")

# plot results
exp.div <- expected_divergence(ne = 20000, migr = 0.00025)
p.huds <- plot_sims(hudson.fst, exp.div) + 
  ggtitle('A') + 
  ylab(expression(italic("F")[ST]))
p.wc <- plot_sims(wc.fst, exp.div) +
  ggtitle('B') +
  ylab(expression(italic(theta)))

# export plots
plot_width <- 4.15
plot_height <- 5.15
CairoPDF('simulation_results.pdf', width = plot_width, height = plot_height)
grid.arrange(p.huds, p.wc)
dev.off()
