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
library(biscale) # create biscale plots
library(cowplot) # draw_plot
library(spatialEco) # crossCorrelation
library(spdep)
library(readxl)

# IMPORT PRE-PROCESSED DATA USED TO FIT MODEL. 
load("WAprevalence/data/data_for_analysis.Rda")

# IMPORT MCMC OUTPUT FROM MODEL
load("WAprevalence/output/mcmc/MCMC_no_covariates_2025_12_16.Rda")

# IMPORT SHAPE FILES FOR WA COUNTIES
load("WAprevalence/data/shape_county_WA.Rda")

# EXAMINE MCMC CONVERGENCE
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(1:234, 20), ", 1]"), ISB = F, filename = "pmp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(1:234, 20), ", 2]"), ISB = F, filename = "death", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(79:234, 20), ", 3]"), ISB = F, filename = "ed", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(1:234, 20), ", 4]"), ISB = F, filename = "hosp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("N[", sample(1:234, 20), "]"), ISB = F, filename = "N", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("lambda[", sample(1:234, 20), "]"), ISB = F, filename = "lambda", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:234, 20), ", 1]"), ISB = F, filename = "f_pmp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:234, 20), ", 2]"), ISB = F, filename = "f_death", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(79:234, 20), ", 3]"), ISB = F, filename = "f_ed", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:234, 20), ", 4]"), ISB = F, filename = "f_hosp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("v[", sample(1:234, 20), "]"), ISB = F, filename = "v", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("u[", sample(1:234, 20), "]"), ISB = F, filename = "u", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("beta[", 1:7, ", 1]"), ISB = F, filename = "beta_pmp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("beta[", 1:7, ", 2]"), ISB = F, filename = "beta_death", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("beta[", 3:7, ", 3]"), ISB = F, filename = "beta_ed", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("beta[", 1:7, ", 4]"), ISB = F, filename = "beta_hosp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "beta.mu", filename = "beta.mu", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "mu", filename = "mu", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.f", filename = "tau.f", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.u", filename = "tau.u", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("eps[", sample(1:234, 20), ", 1]"), ISB = F, filename = "eps_pmp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("eps[", sample(1:234, 20), ", 2]"), ISB = F, filename = "eps_death", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("eps[", sample(1:234, 20), ", 3]"), ISB = F, filename = "eps_ed", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("eps[", sample(1:234, 20), ", 4]"), ISB = F, filename = "eps_hosp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "cov.eps", filename = "cov.eps", wd = "WAprevalence/output/diagnostics")

# EXTRACT POSTERIOR MEANS, 95% CrI (QUANTILES), SD AND NEW GR DIAGNOSTIC STAT

# remove un-sampled MCMC parameters for ED outcome years 2017-2018 as these years missing for this outcome and thus not modeled; 39 counties x 2 years = 78
results <- list(colMeans(samples, na.rm = T),
                  apply(samples, 2,
                        quantile, probs=c(.025,.975), na.rm = T),
                  apply(samples, 2, sd, na.rm = T), 
                  apply(samples, 2, function(x) stable.GR(x, multivariate = F)$psrf))

# specify indices of parameters of interest
pmp_lwr <- which(names(results[[1]])=="pi[1, 1]")
pmp_upr <- which(names(results[[1]])=="pi[273, 1]")
death_lwr <- which(names(results[[1]])=="pi[1, 2]")
death_upr <- which(names(results[[1]])=="pi[273, 2]")
ed_lwr <- which(names(results[[1]])=="pi[1, 3]")
ed_upr <- which(names(results[[1]])=="pi[273, 3]")
hosp_lwr <- which(names(results[[1]])=="pi[1, 4]")
hosp_upr <- which(names(results[[1]])=="pi[273, 4]")
N_lwr <- which(names(results[[1]]) == "N[1]")
N_upr <- which(names(results[[1]]) == "N[273]")
lambda_lwr <- which(names(results[[1]]) == "lambda[1]")
lambda_upr <- which(names(results[[1]]) == "lambda[273]")
beta_lwr <- which(names(results[[1]]) == "beta[1, 1]")
beta_upr <- which(names(results[[1]]) == "beta[7, 4]")
mu_lwr <- which(names(results[[1]]) == "mu[1]")
mu_upr <- which(names(results[[1]]) == "mu[7]")

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
         death_obs_rate = (yfit$death/pop),
         ed_obs_rate = (yfit$ed/pop),
         hosp_obs_rate = (yfit$hosp/pop)
         
  )
  
}

pmp_results <- results_to_tibble(results, "pmp")
death_results <- results_to_tibble(results, "death")
ed_results <- results_to_tibble(results, "ed") %>% mutate(across(c("mean", "lwr95", "upr95", "gr"),  function(x) case_when(year %in% 2017:2018 ~ NA, .default = x)))
hosp_results <- results_to_tibble(results, "hosp")

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

ed_results_csv <-  ed_results %>%
  select(county,
         year,
         mean,
         lwr95,
         upr95)

hosp_results_csv <-  hosp_results %>%
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

write.csv(ed_results_csv,
          file = "WAprevalence/output/tables/ed_results.csv",
          row.names = F)

write.csv(hosp_results_csv,
          file = "WAprevalence/output/tables/hosp_results.csv",
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
        facet_wrap(~year, nrow = 2, ncol = 4) +
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
  ed_obs_rate_map <-  create_choropleth_map(data = ed_results, value = ed_obs_rate, colorbar_type = "monotonic")
  hosp_obs_rate_map <-  create_choropleth_map(data = hosp_results, value = hosp_obs_rate, colorbar_type = "monotonic")
  
  # model maps
  pmp_map <- create_choropleth_map(data = pmp_results, value = mean, colorbar_type = "monotonic")
  death_map <- create_choropleth_map(data = death_results, value = mean, colorbar_type = "monotonic")
  ed_map <- create_choropleth_map(data = ed_results, value = mean, colorbar_type = "monotonic")
  hosp_map <- create_choropleth_map(data = hosp_results, value = mean, colorbar_type = "monotonic")
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

ggsave(filename = "ed_obs_rate.png", 
       plot = ed_obs_rate_map, 
       path = "WAprevalence/output/maps", 
       bg = "White",
       dpi = "retina",
       height = 3,
       width = 10) 

ggsave(filename = "hosp_obs_rate.png", 
       plot = hosp_obs_rate_map, 
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

ggsave(filename = "ed.png", 
       plot = ed_map, 
       path = "WAprevalence/output/maps", 
       bg = "White",
       dpi = "retina",
       height = 3,
       width = 10)

ggsave(filename = "hosp.png", 
       plot = hosp_map, 
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

  # compute posterior distribution of statewide outcome prevalence
  
    # binomial outcomes
    post_binom_outcome_prev <- apply(samples[, beta_lwr:(beta_lwr+13)], 2, ilogit)
  
    # poisson outcomes
    post_pois_outcome_prev <- apply(samples[, (beta_lwr+14):(beta_upr)], 2, exp)
    
  # compute posterior statistics of interest: mean, 95% CrI
    
    # binomial outcomes 
    post_binom_outcome_prev <- list(colMeans(post_binom_outcome_prev),
                                    apply(post_binom_outcome_prev,2,
                                          quantile,probs=c(.025,.975)
                                          )
                         )
    
    # poisson outcomes
    post_pois_outcome_prev <- list(colMeans(post_pois_outcome_prev),
                                    apply(post_pois_outcome_prev,2,
                                          quantile,probs=c(.025,.975)
                                    )
    )
  
  # remove missing years for ED outcome
  post_pois_outcome_prev[[1]][1:2] <- NA 
  post_pois_outcome_prev[[2]][1, 1:2] <- NA
  post_pois_outcome_prev[[2]][2, 1:2] <- NA 
  
  # re-combine posterior means and credible intervals
  post_outcome_prev <- list(c(post_pois_outcome_prev[[1]], post_binom_outcome_prev[[1]]),
                            cbind(post_pois_outcome_prev[[2]], post_binom_outcome_prev[[2]]))
  
  # specify outcomes name and size
  outcomes <- c("ED visit due to opioid misuse",
                "Hospitalization due to opioid misuse",
                "Buprenorphine prescription",
                "Death due to opioid misuse"
  )
  
  K <- length(outcomes)
       
  # estimated prevalence among PWMO (binomial outcomes) & rate per-person among PWMO (poisson)
  # y-axis ticks mapped ot underlying prevalence or rate value
  tibble(pred_beta = post_outcome_prev[[1]],
           lwr95 = post_outcome_prev[[2]][1, ],
           upr95 = post_outcome_prev[[2]][2, ],
           year = rep(2017:2023, K),
           outcome = rep(outcomes, each = 7)
    ) %>%
      mutate(across(c(pred_beta, lwr95, upr95), log)) %>%
      ggplot(aes(x = year, y = pred_beta, fill = outcome)) +
      geom_point(aes(color = outcome)) +
      geom_errorbar(aes(ymin = lwr95, ymax = upr95, color = outcome),
                    width = 0.05) +
      geom_line(linetype = "dashed", aes(color = outcome)) +
      geom_ribbon(aes(ymin = lwr95, ymax = upr95), alpha = 0.1) +
      scale_y_continuous(
                         breaks = c(-10:0),
                         labels = round(c(0, exp(-9:0)), 3)
                        # sec.axis = sec_axis(trans=~.*1, 
                        #                     name="",
                        #                     breaks = c(-6:0), 
                        #                     labels = round(c(exp(-6), exp(-5), exp(-4), exp(-3), exp(-2), exp(-1), exp(0)), 2))
      ) +
      theme_classic() +
      labs(color = "Outcome",
           fill = "Outcome",
           x = "Year",
           y = "Estimated prevalence")
  
  ggsave(filename = "beta_log.png",
         path = "WAprevalence/output",
         dpi = "retina",
         width = 13,
         height = 10,
         units = "cm"
  )

# EXAMINE ESTIMATED MULTI-YEAR AVERAGE STATE-WIDE RATE OF PWMO VS DATA
  
  # compute posterior median and 95% CrI for 2-year state-wide rate of PWMO
  pred_mu_aggr <- samples[,mu_lwr:mu_upr] %>%
                      as_tibble() %>%
                      mutate(`2017-2018`= (`mu[1]` + `mu[2]`)/2,
                             `2018-2019`= (`mu[2]` + `mu[3]`)/2,
                             `2019-2020`= (`mu[3]` + `mu[4]`)/2,
                             `2020-2021`= (`mu[4]` + `mu[5]`)/2,
                             `2021-2022`= (`mu[5]` + `mu[6]`)/2,
                             `2022-2023`= (`mu[6]` + `mu[7]`)/2) %>%
                      select(`2017-2018`,
                             `2018-2019`,
                             `2019-2020`,
                             `2020-2021`,
                             `2021-2022`,
                             `2022-2023`) %>%
                      summarise_all(median) %>%
                      unlist()
  
  CrI_aggr <- samples[,mu_lwr:mu_upr] %>%
                as_tibble() %>%
                mutate(`2017-2018`= (`mu[1]` + `mu[2]`)/2,
                       `2018-2019`= (`mu[2]` + `mu[3]`)/2,
                       `2019-2020`= (`mu[3]` + `mu[4]`)/2,
                       `2020-2021`= (`mu[4]` + `mu[5]`)/2,
                       `2021-2022`= (`mu[5]` + `mu[6]`)/2,
                       `2022-2023`= (`mu[6]` + `mu[7]`)/2) %>%
                select(`2017-2018`,
                       `2018-2019`,
                       `2019-2020`,
                       `2020-2021`,
                       `2021-2022`,
                       `2022-2023`)  %>%
                summarise_all(quantile, probs = c(.025, .975))

  lwr95_aggr <- CrI_aggr %>% slice(1) %>% unlist()
  upr95_aggr <- CrI_aggr %>% slice(2) %>% unlist()
  
  # plot model estimate against NSDUH
  tibble(pred_mu_aggr = pred_mu_aggr,
         lwr95_aggr = lwr95_aggr,
         upr95_aggr = upr95_aggr,
         year = 2017:2022
  ) %>%
    ggplot() +
    geom_point(aes(x = year, y = S, color = "NSDUH Data"),
               data = tibble(year = c(2016:2018, 2021),
                             S = S[1:4])) +
    geom_errorbar(aes(x = year, y = S, ymin = S - 1.96*S.se, ymax = S + 1.96*S.se),
                  width = 0.05,
                  data = tibble(year = c(2016:2018, 2021),
                                S = S[1:4],
                                S.se = S.se[1:4])) +
    geom_line(aes(x = year, y = pred_mu_aggr, color = "Model")) +
    geom_ribbon(aes(x = year, y = pred_mu_aggr, ymin = lwr95_aggr, ymax = upr95_aggr),
                color = "light blue",
                fill = "light blue",
                alpha = 0.5) +
    labs(x = "Year",
         y = "Statewide Prevalence of PWMO") +
    scale_color_manual(name = "", values = c("NSDUH Data" = "black", "Model" = "blue")) +
    theme_classic()
  
  ggsave("2_yr_mu_trend.png",
         device="png",
         path="WAprevalence/output",
         width = 12,
         height = 10,
         units = "cm")

# CREATE BISCALE PLOT #

  # prepare data
  prev_pmp_data <-  N_results %>%
    select(county,
           year,
           mean_prev) %>%
    left_join(pmp_results,
              by = c("county", "year")) %>%
    select(county,
           year,
           mean_prev,
           mean) %>%
    rename(prev_est = mean_prev,
           pmp_est = mean)
  
  # for a given year, combine model estimates w/ spatial info & apply bi_class  
  biscale_data_year <- function(year) {
    
     shape_county_WA %>%
            mutate(NAME = tolower(NAME)) %>%
            rename(county = NAME) %>%
            left_join(prev_pmp_data[prev_pmp_data$year==year,],
                      by = c("county")) %>%
            bi_class(x = prev_est, 
                     y = pmp_est,
                     style = "quantile",
                     dim = 4)
    
    
    
  }
  
  # create biscale class variable for each year
  biscale_data_2017 <- biscale_data_year(2017)
  biscale_data_2018 <- biscale_data_year(2018)
  biscale_data_2019 <- biscale_data_year(2019)
  biscale_data_2020 <- biscale_data_year(2020)
  biscale_data_2021 <- biscale_data_year(2021)
  biscale_data_2022 <- biscale_data_year(2022)
  biscale_data_2023 <- biscale_data_year(2023)

  # stack data
  biscale_data <- biscale_data_2017 %>%
                    bind_rows(biscale_data_2018,
                              biscale_data_2019,
                              biscale_data_2020,
                              biscale_data_2021,
                              biscale_data_2022,
                              biscale_data_2023)

  # create biscale plot
  biscale_legend <- bi_legend(pal = "GrPink2",
                              dim = 4,
                              xlab = "Prevalence",
                              ylab = "Buprenorphine",
                              size = 5)
  
  biscale_map <- biscale_data %>%
                    ggplot() +
                    geom_sf(aes(fill = bi_class), 
                            color = "white",
                            size = 0.1, 
                            show.legend = F) +
                      bi_scale_fill(pal = "GrPink2", dim = 4) +
                      facet_wrap(~year, nrow = 2, ncol = 4) +
                      theme_map() +
                      theme(strip.background = element_rect(fill = "white", color = NA),
                            strip.text = element_text(color = "black",
                                                      size = 7.5, 
                                                      hjust = 0),
                            legend.text = element_text(size = 12),
                            legend.title = element_text(size = 12)
                      )
   
   ggdraw() +
     draw_plot(biscale_map, 0, 0, 1, 1) +
     draw_plot(biscale_legend, 0.77, .26, 0.2, 0.2) 
     #draw_plot(biscale_legend, 0.4, .8, 0.2, 0.2) 
   
   ggsave("biplot_4dim.png",
          device="png",
          path="WAprevalence/output/maps",
          width = 12,
          height = 10,
          units = "cm",
          bg = "white")
   
# COMPUTE SPATIAL CROSS CORRELATION USING LOCAL MORAN'S I #

   # compute adjacency matrix for NC counties
   WA_map <- shape_county_WA[order(shape_county_WA$COUNTYFP),] #convert to data frame
   
   nbmat <- poly2nb(WA_map)
   
   n <- length(shape_county_WA$NAME) # should be 39
   
   A <- matrix(0,n,n)
   
   for(i in 1:n){
     
     A[i,unlist(nbmat[[i]])]=1
     
   }
   
  # compute Moran's I for a given year of the data
  crossCorrelation_year <- function(year) {
    
    crossCorrelation_data <- prev_pmp_data[prev_pmp_data$year == year,]
    prev_est <- crossCorrelation_data %>% pull(prev_est)
    pmp_est <- crossCorrelation_data %>% pull(pmp_est)
    
    crossCorrelation(prev_est,
                     pmp_est,
                     w = A,
                     dist.function = "none")
    
  }
   
   # save crossCorrelation values
   crossCorrelation_2017 <- crossCorrelation_year(2017)
   crossCorrelation_2018 <- crossCorrelation_year(2018)
   crossCorrelation_2019 <- crossCorrelation_year(2019)
   crossCorrelation_2020 <- crossCorrelation_year(2020)
   crossCorrelation_2021 <- crossCorrelation_year(2021)
   crossCorrelation_2022 <- crossCorrelation_year(2022)
   crossCorrelation_2023 <- crossCorrelation_year(2023)
 
  # extract clusters and put into data frame
   cluster_data <- function(data, year) {
     
     prev_pmp_data_year <- prev_pmp_data[prev_pmp_data$year==year,] %>%
       mutate(cluster = data$clusters,
              bi_class = case_when(
                cluster == "Low.Low" ~ "1-1",
                cluster == "Low.High" ~ "1-2",
                cluster == "High.Low" ~ "2-1",
                cluster == "High.High" ~ "2-2"
               )
              )
     
     shape_county_WA %>%
       mutate(NAME = tolower(NAME)) %>%
       rename(county = NAME) %>%
       left_join(prev_pmp_data_year,
                 by = c("county"))
   
   }
   
   # extract tidy cluster data for each year
   cluster_data_2017 <- cluster_data(crossCorrelation_2017, 2017)
   cluster_data_2018 <- cluster_data(crossCorrelation_2018, 2018)
   cluster_data_2019 <- cluster_data(crossCorrelation_2019, 2019)
   cluster_data_2020 <- cluster_data(crossCorrelation_2020, 2020)
   cluster_data_2021 <- cluster_data(crossCorrelation_2021, 2021)
   cluster_data_2022 <- cluster_data(crossCorrelation_2022, 2022)
   cluster_data_2023 <- cluster_data(crossCorrelation_2023, 2023)
   
   # stack data
   cluster_data <- cluster_data_2017 %>%
                      bind_rows(cluster_data_2018,
                                cluster_data_2019,
                                cluster_data_2020,
                                cluster_data_2021,
                                cluster_data_2022,
                                cluster_data_2023)
   
   # create biscale plot using cluster from crossCorrelation
   cluster_legend <- bi_legend(pal = "GrPink",
                               dim = 2,
                               xlab = "Prevalence",
                               ylab = "Buprenorphine",
                               size = 5)
   
   cluster_map <- cluster_data %>%
     ggplot() +
     geom_sf(aes(fill = bi_class), 
             color = "white",
             size = 0.1, 
             show.legend = F) +
     bi_scale_fill(pal = "GrPink", dim = 2) +
     facet_wrap(~year, nrow=2, ncol =4) +
     theme_map() +
     theme(strip.background = element_rect(fill = "white", color = NA),
           strip.text = element_text(color = "black",
                                     size = 12, 
                                     hjust = 0),
           legend.text = element_text(size = 12),
           legend.title = element_text(size = 12)
     )
   
   ggdraw() +
     draw_plot(cluster_map, 0, 0, 1, 1) +
     draw_plot(cluster_legend, 0.4, .8, 0.2, 0.2) 
   
   ggsave("crossCorr_biplot_2dim.png",
          device="png",
          path="WAprevalence/output/maps",
          width = 12,
          height = 10,
          units = "cm",
          bg = "white")
   
# EXAMINE MEDICAID DATA ON MOUD FOR BUPRENORPHINE ESTIMATE QUALITY CONTROL #
# POPULATION: MEDICAID ENROLLEES WITH OUD IN WA FROM 2018-2024 #
  
  # import raw data
  medicaid_moud <- read_excel("WAprevalence/data/MAT_ResultsHCA_202409_ewa_smsuppressed.xlsx")
  
  # simplify data for QC purpose
  medicaid_moud_clean <- medicaid_moud %>%
    mutate(buprenorphine_naloxone = as.numeric(`Treatment of Medicaid enrollees with OUD with buprenorphine-naloxone`),
           buprenorphine = as.numeric(`Treatment of Medicaid enrollees with OUD with buprenorphine`)
    ) %>%
    separate_wider_delim(clndr_qtr, delim = "Q", names = c("year", "quarter")) %>%
    filter(!(county %in% c("ACH WIDE", "STATE WIDE"))) %>%
    select(year, quarter, county, buprenorphine_naloxone, buprenorphine)

  # visualize data 
  
    # merge in county spatial info and generate map of medicaid data
    map_medicaid_data <- function(outcome) {
      
      shape_county_WA %>%
      rename(county = NAME) %>%
      left_join(medicaid_moud_clean, by = c("county")) %>%
      ggplot() +
      geom_sf(aes(fill = {{outcome}})) +
      scale_fill_gradient(low = "white",
                          high = "red",
                          guide = guide_colorbar(barheight = 12)) +
      labs(fill = NULL) +
      theme_map() +
      theme(legend.position = "right") +
      facet_wrap(~year + quarter, nrow = 9, ncol = 4) +
      theme(strip.background = element_rect(fill = "white", color = NA),
            strip.text = element_text(color = "black",
                                      size = 12, 
                                      hjust = 0),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12)
      )
      
    }
    
    # generate maps
    buprenorphine_naloxone_map <- map_medicaid_data(buprenorphine_naloxone)
    buprenorphine_map <- map_medicaid_data(buprenorphine)
    
    # save maps
    ggsave("medicaid_moud_buprenorphine_naloxone_county.png",
           buprenorphine_naloxone_map,
           device="png",
           path="WAprevalence/output/maps",
           width = 15,
           height = 20,
           units = "cm",
           bg = "white")
    
    ggsave("medicaid_moud_buprenorphine_county.png",
           buprenorphine_map,
           device="png",
           path="WAprevalence/output/maps",
           width = 15,
           height = 20,
           units = "cm",
           bg = "white")
    
    # examine county-level trends via facet-wrap
    county_trend_medicaid_data <- function(outcome) {
      
      medicaid_moud_clean %>%
        mutate(across(c(year, quarter), as.numeric),
               date = ymd(paste(year, (quarter - 1)*3 + 1, 1, sep ="-"))) %>%
        ggplot(aes(x = date, y = {{outcome}})) +
        geom_line(alpha = 0.5) +
        facet_wrap(~county) +  
        theme_classic()
      
    }
    
    # generate county time trends
    county_trend_buprenorphine_naloxone <- county_trend_medicaid_data(buprenorphine_naloxone)
    county_trend_buprenorphine <- county_trend_medicaid_data(buprenorphine)
    
    # save county time trends
    ggsave("medicaid_moud_buprenorphine_naloxone_county_trend.png",
           county_trend_buprenorphine_naloxone,
           device="png",
           path="WAprevalence/output/maps",
           width = 30,
           height = 20,
           units = "cm",
           bg = "white")
    
    ggsave("medicaid_moud_buprenorphine_county_trend.png",
           county_trend_buprenorphine,
           device="png",
           path="WAprevalence/output/maps",
           width = 30,
           height = 20,
           units = "cm",
           bg = "white")
    
    # examine state-level trends
    state_trend_medicaid_data <- function(outcome) {
      
      medicaid_moud %>%
        mutate(buprenorphine_naloxone = as.numeric(`Treatment of Medicaid enrollees with OUD with buprenorphine-naloxone`),
               buprenorphine = as.numeric(`Treatment of Medicaid enrollees with OUD with buprenorphine`)
        ) %>%
        filter(county == "STATE WIDE") %>%
        separate_wider_delim(clndr_qtr, delim = "Q", names = c("year", "quarter")) %>%
        mutate(across(c(year, quarter), as.numeric),
               date = ymd(paste(year, (quarter - 1)*3 + 1, 1, sep ="-"))) %>%
        ggplot(aes(x = date, y = {{outcome}})) +
        geom_line() +
        theme_classic()
      
    }
    
    # generate state time trends
    state_trend_buprenorphine_naloxone <- state_trend_medicaid_data(buprenorphine_naloxone)
    state_trend_buprenorphine <- state_trend_medicaid_data(buprenorphine)
    
    # save state time trends
    ggsave("medicaid_moud_buprenorphine_naloxone_state_trend.png",
           state_trend_buprenorphine_naloxone,
           device="png",
           path="WAprevalence/output/maps",
           units = "cm",
           bg = "white")
    
    ggsave("medicaid_moud_buprenorphine_state_trend.png",
           state_trend_buprenorphine,
           device="png",
           path="WAprevalence/output/maps",
           units = "cm",
           bg = "white")
