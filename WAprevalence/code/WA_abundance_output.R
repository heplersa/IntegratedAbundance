# EXAMINE WA ABUNDANCE MODEL MCMC CONVERGENCE & GENERATE FIGURES/TABLES #
# BRIAN N. WHITE #
# 2025-01-19 #

# LOAD R PACKAGES
library(tidyverse) # data manipulation and visualization
library(reshape2) # melt
library(cowplot) # combine plots
library(nimble) # Bayesian inference for multi-level models
library(stableGR) # new GR statistic
library(ggthemes) # theme_map
library(sf) # import shape files
library(tigris) # pull shape files from internet
library(flextable) # make pretty tables

# IMPORT PRE-PROCESSED DATA USED TO FIT MODEL. 
load("WAprevalence/data/data_for_analysis.Rda")

# IMPORT MCMC OUTPUT FROM MODEL
load("WAprevalence/output/MCMC_no_covariates.Rda")

# IMPORT SHAPE FILES FOR WA COUNTIES
load("WAprevalence/data/shape_county_WA.Rda")

# EXAMINE MCMC CONVERGENCE
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(1:234, 20), ", 1]"), ISB = F, filename = "pmp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(1:234, 20), ", 2]"), ISB = F, filename = "death", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("N[", sample(1:234, 20), "]"), ISB = F, filename = "N", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("lambda[", sample(1:234, 20), "]"), ISB = F, filename = "lambda", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:234, 20), ", 1]"), ISB = F, filename = "f_pmp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:234, 20), ", 2]"), ISB = F, filename = "f_death", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("v[", sample(1:234, 20), "]"), ISB = F, filename = "v", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("u[", sample(1:234, 20), "]"), ISB = F, filename = "u", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "beta", filename = "beta", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "beta.mu", filename = "beta.mu", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "mu", filename = "mu", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau", filename = "tau", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.f", filename = "tau.f", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.u", filename = "tau.u", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("eps[", sample(1:234, 20), ", 1]"), ISB = F, filename = "eps.pmp", "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("eps[", sample(1:234, 20), ", 2]"), ISB = F, filename = "eps.death", "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "cov.eps", filename = "cov.eps", wd = "WAprevalence/output/diagnostics")

# EXTRACT POSTERIOR MEANS, 95% CrI (QUANTILES), SD AND NEW GR DIAGNOSTIC STAT
results <- list(colMeans(samples),
                apply(samples,2,
                      quantile,probs=c(.025,.975)),
                apply(samples,2,sd), 
                apply(samples, 2, function(x) stable.GR(x, multivariate = F)$psrf))

# specify indices of parameters of interest
pmp_lwr <- which(names(results[[1]])=="pi[1, 1]")
pmp_upr <- which(names(results[[1]])=="pi[234, 1]")
death_lwr <- which(names(results[[1]])=="pi[1, 2]")
death_upr <- which(names(results[[1]])=="pi[234, 2]")
N_lwr <- which(names(results[[1]]) == "N[1]")
N_upr <- which(names(results[[1]]) == "N[234]")
lambda_lwr <- which(names(results[[1]]) == "lambda[1]")
lambda_upr <- which(names(results[[1]]) == "lambda[234]")
beta_lwr <- which(names(results[[1]]) == "beta[1, 1]")
beta_upr <- which(names(results[[1]]) == "beta[6, 2]")
mu_lwr <- which(names(results[[1]]) == "mu[1]")
mu_upr <- which(names(results[[1]]) == "mu[6]")

# create tidy data sets of estimates merged with corresponding spatio-temporal data
results_to_tibble <- function(results, par) {
  
  par_lwr <- get(paste(par, "_lwr", sep = ""), envir = .GlobalEnv)
  par_upr <- get(paste(par, "_upr", sep = ""), envir = .GlobalEnv)
  
  tibble(par = names(results[[1]][par_lwr:par_upr]),
         county = yfit$county,
         year = yfit$year,
         pop = yfit$pop,
         mean = results[[1]][par_lwr:par_upr],
         lwr95 = results[[2]][1, par_lwr:par_upr],
         upr95 = results[[2]][2, par_lwr:par_upr],
         sd = results[[3]][par_lwr:par_upr],
         gr = results[[4]][par_lwr:par_upr],
         pmp_obs_rate = (yfit$pmp/pop),
         death_obs_rate = (yfit$death/pop)
         
  )
  
}

pmp_results <- results_to_tibble(results, "pmp")
death_results <- results_to_tibble(results, "death")
lambda_results <- results_to_tibble(results, "lambda") %>% 
                    mutate(CrI = case_when(
                                          lwr95 > 1  ~ "95% CrI > 1",
                                          upr95 < 1  ~ "95% CrI < 1",
                                          .default = "95% CrI contains 1"),
                           CrI = fct_relevel(CrI, c("95% CrI < 1", "95% CrI contains 1", "95% CrI > 1"))
                           )
N_results <- results_to_tibble(results, "N") %>% 
                mutate(mean_prev = mean/pop,
                       lwr95_prev = lwr95/pop,
                       upr95_prev = upr95/pop)

# CREATE TABLES OF POSTERIOR MEANS W/ 95% CrI FOR SELECT PARAMETERS
pmp_results_csv <-  pmp_results %>%
                      select(county,
                             year,
                             mean,
                             lwr95,
                             upr95)

death_results_csv <-  death_results %>%
                        select(county,
                               year,
                               mean,
                               lwr95,
                               upr95)

N_results_csv <-  N_results %>%
                    select(county,
                           year,
                           mean,
                           lwr95,
                           upr95)

N_prev_results_csv <- N_results %>%
                        select(county,
                               year,
                               mean_prev,
                               lwr95_prev,
                               upr95_prev)

write.csv(pmp_results_csv,
          file = "WAprevalence/output/tables/pmp_results.csv",
          row.names = F)

write.csv(death_results_csv,
          file = "WAprevalence/output/tables/death_results.csv",
          row.names = F)

write.csv(N_results_csv,
          file = "WAprevalence/output/tables/N_results.csv",
          row.names = F)

write.csv(N_prev_results_csv,
          file = "WAprevalence/output/tables/N_prev_results.csv",
          row.names = F)

# EXAMINE GR STATISTICS FOR SELECT PARAMETERS
results[[4]][N_lwr:N_upr] %>% mean
results[[4]][N_lwr:N_upr] %>% median
results[[4]][N_lwr:N_upr] %>% sd
results[[4]][lambda_lwr:lambda_upr] %>% mean
results[[4]][lambda_lwr:lambda_upr] %>% median
results[[4]][lambda_lwr:lambda_upr] %>% sd

# CREATE CHOROPLETH MAPS

# function to create choropleth maps
create_choropleth_map <- function(data, value, colorbar_type = NULL, colorbar_title = NULL) {
  
  p <- shape_county_WA %>%
          mutate(NAME = tolower(NAME)) %>%
          rename(county = NAME) %>%
          left_join(data, by = c("county")) %>%
          ggplot() +
          geom_sf(aes(fill = {{value}})) 
  
  # is color bar monotonic or diverging?
  p <- {if(colorbar_type == "monotonic"){
    
    p + scale_fill_gradient(low = "white",
                            high = "red",
                            guide = guide_colorbar(barheight = 12))
    
  } else if(colorbar_type == "diverging") {
    
    p + scale_fill_gradient2(low = "blue",
                             mid = "white",
                             high = "red",
                             midpoint = 1,
                             guide = guide_colorbar(barheight = 12))
    
  } else {
    
    p + scale_fill_manual(values = c("95% CrI < 1" = "blue", 
                                     "95% CrI contains 1" = "white",
                                     "95% CrI > 1" = "red"))
    
    
    
  }
    
  } 
  
  # finish formatting
  p <- p + 
        labs(fill = colorbar_title) +
        theme_map() +
        theme(legend.position = "right") +
        facet_wrap(~year) +
        theme(strip.background = element_rect(fill = "white", color = NA),
              strip.text = element_text(color = "black",
                                        size = 12, 
                                        hjust = 0),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12)
        )
  
}

# generate and save maps

  # observed maps
  pmp_obs_rate_map <-  create_choropleth_map(data = pmp_results, value = pmp_obs_rate, colorbar_type = "monotonic")
  death_obs_rate_map <-  create_choropleth_map(data = death_results, value = death_obs_rate, colorbar_type = "monotonic")
  
  # model maps
  pmp_map <- create_choropleth_map(data = pmp_results, value = mean, colorbar_type = "monotonic")
  death_map <- create_choropleth_map(data = death_results, value = mean, colorbar_type = "monotonic")
  lambda_map <- create_choropleth_map(data = lambda_results, value = mean, colorbar_type = "diverging")
  lambda_CrI_map <- create_choropleth_map(data = lambda_results, value = CrI, colorbar_type = "other") + theme(legend.position = "right")
  N_map <- create_choropleth_map(data = N_results, value = mean_prev, colorbar_type = "monotonic")

  
ggsave(filename = "pmp_obs_rate.png", 
       plot = pmp_obs_rate_map, 
       path = "WAprevalence/output/maps", 
       bg = "White",
       dpi = "retina",
       height = 3,
       width = 10)

ggsave(filename = "death_obs_rate.png", 
       plot = death_obs_rate_map, 
       path = "WAprevalence/output/maps", 
       bg = "White",
       dpi = "retina",
       height = 3,
       width = 10) 

ggsave(filename = "pmp.png", 
       plot = pmp_map, 
       path = "WAprevalence/output/maps", 
       bg = "White",
       dpi = "retina",
       height = 3,
       width = 10)

ggsave(filename = "death.png", 
       plot = death_map, 
       path = "WAprevalence/output/maps", 
       bg = "White",
       dpi = "retina",
       height = 3,
       width = 10)

ggsave(filename = "lambda.png", 
       plot = lambda_map, 
       path = "WAprevalence/output/maps", 
       bg = "White",
       dpi = "retina",
       height = 3,
       width = 10)

ggsave(filename = "lambda_CrI.png", 
       plot = lambda_CrI_map, 
       path = "WAprevalence/output/maps", 
       bg = "White",
       dpi = "retina",
       height = 3,
       width = 10)

ggsave(filename = "N.png", 
       plot = N_map, 
       path = "WAprevalence/output/maps", 
       bg = "White",
       dpi = "retina",
       height = 3,
       width = 10)

# EXAMINE ESTIMATED TIME-VARYING INTERCEPTS FOR EACH OUTCOME

# note that this is on the log-odds scale; no covariates so this shows the estimated state-wide prevalence of each outcome among PWMO on the logit scale
tibble(pred_beta = results[[1]][beta_lwr:beta_upr],
       lwr95 = results[[2]][1, beta_lwr:beta_upr],
       upr95 = results[[2]][2, beta_lwr:beta_upr],
       year = rep(2017:2022, 2),
       outcome = rep(c("Buprenorphine prescription",
                       "Death due to opioid misuse"), each = 6)
) %>%
  ggplot(aes(x = year, y = pred_beta, fill = outcome)) +
  geom_point(aes(color = outcome)) +
  geom_errorbar(aes(ymin = lwr95, ymax = upr95, color = outcome),
                width = 0.05) +
  geom_line(linetype = "dashed", aes(color = outcome)) +
  geom_ribbon(aes(ymin = lwr95, ymax = upr95), alpha = 0.1) +
  #scale_y_continuous(sec.axis = sec_axis(trans=~.*1, 
  #                                       name="",
  #                                       breaks = c(-6 ,-4, -2 , 0), 
  #                                       labels = round(c(exp(-6), exp(-4), exp(-2), exp(0)), 2))
  #) +
  theme_classic() +
  labs(color = "Outcome",
       fill = "Outcome",
       x = "Year",
       y = "Estimated prevalence (log-odds)")

ggsave(filename = "beta_log_odds.png",
       path = "WAprevalence/output",
       dpi = "retina",
       width = 13,
       height = 10,
       units = "cm"
)

# EXAMINE ESTIMATED MULTI-YEAR AVERAGE STATE-WIDE RATE OF PWMO VS DATA
pred_mu_aggr <- samples[,mu_lwr:mu_upr] %>%
                    as_tibble() %>%
                    mutate(`2017-2018`= (`mu[1]` + `mu[2]`)/2,
                           `2018-2019`= (`mu[2]` + `mu[3]`)/2,
                           `2019-2020`= (`mu[3]` + `mu[4]`)/2,
                           `2020-2021`= (`mu[4]` + `mu[5]`)/2,
                           `2021-2022`= (`mu[5]` + `mu[6]`)/2) %>%
                    select(`2017-2018`,
                           `2018-2019`,
                           `2019-2020`,
                           `2020-2021`,
                           `2021-2022`) %>%
                    summarise_all(median) %>%
                    unlist()

CrI_aggr <- samples[,mu_lwr:mu_upr] %>%
              as_tibble() %>%
              mutate(`2017-2018`= (`mu[1]` + `mu[2]`)/2,
                     `2018-2019`= (`mu[2]` + `mu[3]`)/2,
                     `2019-2020`= (`mu[3]` + `mu[4]`)/2,
                     `2020-2021`= (`mu[4]` + `mu[5]`)/2,
                     `2021-2022`= (`mu[5]` + `mu[6]`)/2) %>%
              select(`2017-2018`,
                     `2018-2019`,
                     `2019-2020`,
                     `2020-2021`,
                     `2021-2022`) %>%
              summarise_all(quantile, probs = c(.025, .975))

lwr95_aggr <- CrI_aggr %>% slice(1) %>% unlist()
upr95_aggr <- CrI_aggr %>% slice(2) %>% unlist()

tibble(pred_mu_aggr = pred_mu_aggr,
       lwr95_aggr = lwr95_aggr,
       upr95_aggr = upr95_aggr,
       year = 2017:2021
       ) %>%
  ggplot() +
  geom_point(aes(x = year, y = S),
             data = tibble(year = c(2016:2018, 2021),
                           S = S[1:4])) +
  geom_errorbar(aes(x = year, y = S, ymin = S - 1.96*S.se, ymax = S + 1.96*S.se),
                width = 0.05,
                data = tibble(year = c(2016:2018, 2021),
                              S = S[1:4],
                              S.se = S.se[1:4])) +
  geom_line(aes(x = year, y = pred_mu_aggr), color = "blue") +
  geom_ribbon(aes(x = year, y = pred_mu_aggr, ymin = lwr95_aggr, ymax = upr95_aggr),
              color = "light blue",
              fill = "light blue",
              alpha = 0.5) +
  labs(x = "Year",
       y = "Statewide Prevalence of PWMO") +
  theme_classic()

ggsave("2_yr_mu_trend.png",
       device="png",
       path="WAprevalence/output",
       width = 12,
       height = 10,
       units = "cm")
