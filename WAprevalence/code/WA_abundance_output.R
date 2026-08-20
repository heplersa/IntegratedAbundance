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

# IMPORT PRE-PROCESSED DATA USED TO FIT MODEL. 
load("WAprevalence/data/data_for_analysis.Rda")

# IMPORT MCMC OUTPUT FROM MODEL
load("WAprevalence/output/mcmc/MCMC_N_pois_3_chains_2_mill_each_2026_08_14.Rda")

# IMPORT SHAPE FILES FOR WA COUNTIES
load("WAprevalence/data/shape_county_WA.Rda")

# EXAMINE MCMC CONVERGENCE
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(1:273, 20), ", 1]"), ISB = F, filename = "pmp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(1:273, 20), ", 2]"), ISB = F, filename = "death", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(79:273, 20), ", 3]"), ISB = F, filename = "ed", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(1:273, 20), ", 4]"), ISB = F, filename = "hosp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("N[", sample(1:273, 20), "]"), ISB = F, filename = "N", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("lambda[", sample(1:273, 20), "]"), ISB = F, filename = "lambda", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:273, 20), ", 1]"), ISB = F, filename = "f_pmp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:273, 20), ", 2]"), ISB = F, filename = "f_death", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(79:273, 20), ", 3]"), ISB = F, filename = "f_ed", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:273, 20), ", 4]"), ISB = F, filename = "f_hosp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("v[", sample(1:273, 20), "]"), ISB = F, filename = "v", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("u[", sample(1:273, 20), "]"), ISB = F, filename = "u", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("beta[", 1:7, ", 1]"), ISB = F, filename = "beta_pmp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("beta[", 1:7, ", 2]"), ISB = F, filename = "beta_death", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("beta[", 3:7, ", 3]"), ISB = F, filename = "beta_ed", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("beta[", 1:7, ", 4]"), ISB = F, filename = "beta_hosp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "beta.mu", filename = "beta.mu", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "mu", filename = "mu", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.f", filename = "tau.f", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.u", filename = "tau.u", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("eps[", sample(1:273, 20), ", 1]"), ISB = F, filename = "eps_pmp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("eps[", sample(1:273, 20), ", 2]"), ISB = F, filename = "eps_death", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("eps[", sample(1:273, 20), ", 3]"), ISB = F, filename = "eps_ed", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("eps[", sample(1:273, 20), ", 4]"), ISB = F, filename = "eps_hosp", wd = "WAprevalence/output/diagnostics")
MCMCvis::MCMCtrace(samples, params = "prec.eps", filename = "prec.eps", wd = "WAprevalence/output/diagnostics")

# COMBINE CHAINS FOR POSTERIOR INFERENCE
samples <- rbind(samples[[1]], samples[[2]], samples[[3]])

# EXTRACT POSTERIOR MEANS, 95% CrI (QUANTILES), SD
results <- list(colMeans(samples, na.rm = T),
                  apply(samples, 2,
                        quantile, probs=c(.025,.975), na.rm = T),
                  apply(samples, 2, sd, na.rm = T))

# save posterior summaries for downstream use without the full MCMC samples
save(results, file = "WAprevalence/output/results.Rda")

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
         pmp_obs_rate = (yfit$pmp/pop),
         death_obs_rate = (yfit$death/pop),
         ed_obs_rate = (yfit$ed/pop),
         hosp_obs_rate = (yfit$hosp/pop)
         
  )
  
}

pmp_results <- results_to_tibble(results, "pmp")
death_results <- results_to_tibble(results, "death")
ed_results <- results_to_tibble(results, "ed")  %>% mutate(across(c("mean", "lwr95", "upr95"),  function(x) case_when(year %in% 2017:2018 ~ NA, .default = x)))
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

# CREATE STATE-WIDE PREVALENCE TABLE FOR MANUSCRIPT

  # prevalence column is the posterior for mu; PWUO-HR column sums the county-level estimates
  statewide_prevalence_table <- N_results %>%
                                  group_by(year) %>%
                                  summarise(N_est = sum(mean),
                                            pop = sum(pop)) %>%
                                  mutate(prev = 100*results[[1]][mu_lwr:mu_upr],
                                         lwr95 = 100*results[[2]][1, mu_lwr:mu_upr],
                                         upr95 = 100*results[[2]][2, mu_lwr:mu_upr]) %>%
                                  transmute(Year = year,
                                            `Prevalence % (95% CrI)` = sprintf("%.2f (%.2f, %.2f)", prev, lwr95, upr95),
                                            `Estimated PWUO-HR` = format(round(N_est), big.mark = ","),
                                            `Population 12+` = format(pop, big.mark = ","))

  write.csv(statewide_prevalence_table,
            file = "WAprevalence/output/tables/statewide_prevalence_table.csv",
            row.names = F)

  save_as_docx(autofit(flextable(statewide_prevalence_table)),
               path = "WAprevalence/output/tables/statewide_prevalence_table.docx")

# EXAMINE MODIFIED GELMAN-RUBIN (GR) STATISTICS FOR SELECT PARAMETERS

  # exclude 2017-2018 for ED outcome as these were not modeled and therefore not sampled
  unsampled_params <- c(which(colnames(samples)=="pi[1, 3]"):which(colnames(samples)=="pi[78, 3]"),
                        which(colnames(samples)=="f[1, 3]"):which(colnames(samples)=="f[78, 3]"),
                        which(colnames(samples)=="beta[1, 3]"):which(colnames(samples)=="beta[2, 3]"))
  
  # compute modified GR stat for each sampled parameter
  gr_stats <- apply(samples[ ,-unsampled_params], 2, function(x) stable.GR(x, multivariate = F)$psrf)
  
  # examine distribution of GR stats across counties and years for parameters of interest
  gr_stat_summary <- function(par) {
    
    par_lwr <- get(paste(par, "_lwr", sep = ""), envir = .GlobalEnv)
    par_upr <- get(paste(par, "_upr", sep = ""), envir = .GlobalEnv)
    
    c(mean(gr_stats[par_lwr:par_upr]), median(gr_stats[par_lwr:par_upr]), sd(gr_stats[par_lwr:par_upr]))
    
    
  }
  
  
  gr_stat_summary("pmp")
  gr_stat_summary("death")
  gr_stat_summary("ed")
  gr_stat_summary("hosp")
  gr_stat_summary("N")
  gr_stat_summary("lambda")
  gr_stat_summary("beta")
  gr_stat_summary("mu")

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

  # observed maps; rates displayed per 100,000 population aged 12+
  pmp_obs_rate_map <-  create_choropleth_map(data = pmp_results %>% mutate(pmp_obs_rate = 1e5*pmp_obs_rate), value = pmp_obs_rate, colorbar_type = "monotonic", colorbar_title = "Rate per\n100,000")
  death_obs_rate_map <-  create_choropleth_map(data = death_results %>% mutate(death_obs_rate = 1e5*death_obs_rate), value = death_obs_rate, colorbar_type = "monotonic", colorbar_title = "Rate per\n100,000")
  ed_obs_rate_map <-  create_choropleth_map(data = ed_results %>% mutate(ed_obs_rate = 1e5*ed_obs_rate), value = ed_obs_rate, colorbar_type = "monotonic", colorbar_title = "Rate per\n100,000")
  hosp_obs_rate_map <-  create_choropleth_map(data = hosp_results %>% mutate(hosp_obs_rate = 1e5*hosp_obs_rate), value = hosp_obs_rate, colorbar_type = "monotonic", colorbar_title = "Rate per\n100,000")
  
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

# CREATE COMPARISON MAP OF OBSERVED DEATH RATE FOR POP VS FOR PWMO IN 2023 #

  # compute county-level death rate for pop and for PWMO
  death_rate_data <- N_results %>%
    select(county,
           year,
           pop,
           mean) %>%
    rename(N_est = mean) %>%
    mutate(death = yfit$death,
           death_rate_pop = death/pop,
           death_rate_PWMO = death/N_est)

  # compute the statewide observed death rate for pop and for PWMO in 2023
  state_wide_death_rate <- death_rate_data %>%
    filter(year == 2023) %>%
    summarise(pop_rate = sum(death)/sum(pop),
              PWMO_rate = sum(death)/sum(N_est))

  write.csv(death_rate_data,
            file = "WAprevalence/output/tables/death_rate_data.csv",
            row.names = F)

  write.csv(state_wide_death_rate,
            file = "WAprevalence/output/tables/state_wide_death_rate.csv",
            row.names = F)

# EXAMINE ESTIMATED TIME-VARYING INTERCEPTS FOR EACH OUTCOME

  # compute posterior distribution of statewide outcome prevalence
  
    # binomial outcomes
    post_binom_outcome_prev <- apply(samples[, beta_lwr:(beta_lwr+13)], 2, ilogit)
  
    # poisson outcomes; exclude 2017-2018 for ED visit outcome
    post_pois_outcome_prev <- apply(samples[, (beta_lwr+16):(beta_upr)], 2, exp)
    
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
                                          quantile,probs=c(.025,.975),
                                    )
    )
  
  # re-combine posterior means and credible intervals; NA for 2017-2018 ED visit outcome
  post_outcome_prev <- list(c(post_binom_outcome_prev[[1]], c(NA, NA), post_pois_outcome_prev[[1]]),
                            cbind(post_binom_outcome_prev[[2]], matrix(c(NA, NA, NA, NA), nrow=2, ncol=2), post_pois_outcome_prev[[2]]))
  
  # specify outcomes name and size
  outcomes <- c("Buprenorphine receipt",
                "Overdose death",
                "ED visit",
                "Hospitalization"
  )

  # colorblind-safe palette (Okabe-Ito); overlapping death/hosp series get the most distinct pair
  outcome_cols <- c("Buprenorphine receipt" = "#009E73",
                    "Overdose death" = "#0072B2",
                    "ED visit" = "#CC79A7",
                    "Hospitalization" = "#D55E00")
  
  K <- length(outcomes)
       
  # estimated prevalence among PWMO (binomial outcomes) & rate per-person among PWMO (poisson)
  # y-axis ticks mapped ot underlying prevalence or rate value
  outcome_trend_plot <- tibble(pred_beta = post_outcome_prev[[1]],
         lwr95 = post_outcome_prev[[2]][1, ],
         upr95 = post_outcome_prev[[2]][2, ],
         year = rep(2017:2023, K),
         outcome = rep(outcomes, each = 7)
    ) %>%
      mutate(across(c(pred_beta, lwr95, upr95), log)) %>%
      ggplot(aes(x = year, y = pred_beta, fill = outcome)) +
      geom_point(aes(color = outcome, shape = outcome)) +
      geom_errorbar(aes(ymin = lwr95, ymax = upr95, color = outcome),
                    width = 0.05) +
      geom_line(linetype = "dashed", aes(color = outcome)) +
      geom_ribbon(aes(ymin = lwr95, ymax = upr95), alpha = 0.1) +
      scale_color_manual(values = outcome_cols) +
      scale_fill_manual(values = outcome_cols) +
      scale_y_continuous(
                         breaks = c(-10:0),
                         labels = round(c(0, exp(-9:0)), 3)
                        # sec.axis = sec_axis(trans=~.*1, 
                        #                     name="",
                        #                     breaks = c(-6:0), 
                        #                     labels = round(c(exp(-6), exp(-5), exp(-4), exp(-3), exp(-2), exp(-1), exp(0)), 2))
      ) +
      theme_bw(base_size = 11) +
      theme(panel.grid = element_blank(),
            axis.text = element_text(color = "black", size = 8.5)) +
      labs(color = "Outcome",
           fill = "Outcome",
           shape = "Outcome",
           x = "Year",
           y = "Statewide rate among PWUO-HR")
  
  ggsave(filename = "beta_log.png",
         plot = outcome_trend_plot,
         path = "WAprevalence/output",
         dpi = "retina",
         width = 13,
         height = 10,
         units = "cm"
  )

# EXAMINE ESTIMATED MULTI-YEAR AVERAGE STATE-WIDE RATE OF PWMO VS DATA
  
  # the 2016-2017 and 2023-2024 survey windows need model-implied prevalence for 2016 and 2024,
  # neither of which is a model node (the study period runs t = 1 (2017) through t = T = 7 (2023)).
  # mu_t is a deterministic function of beta.mu, so evaluate the trend at t = 0 and t = 8 on the
  # same posterior draws. this keeps the comparison joint and lets the fitted 2-year trend span
  # every observed NSDUH estimate rather than starting and stopping inside them
  mu_2016 <- plogis(samples[, "beta.mu[1]"])
  mu_2024 <- plogis(samples[, "beta.mu[1]"] + 8*samples[, "beta.mu[2]"])

  # compute posterior median and 95% CrI for 2-year state-wide rate of PWMO
  pred_mu_aggr <- samples[,mu_lwr:mu_upr] %>%
                      as_tibble() %>%
                      mutate(`2016-2017`= (mu_2016 + `mu[1]`)/2,
                             `2017-2018`= (`mu[1]` + `mu[2]`)/2,
                             `2018-2019`= (`mu[2]` + `mu[3]`)/2,
                             `2019-2020`= (`mu[3]` + `mu[4]`)/2,
                             `2020-2021`= (`mu[4]` + `mu[5]`)/2,
                             `2021-2022`= (`mu[5]` + `mu[6]`)/2,
                             `2022-2023`= (`mu[6]` + `mu[7]`)/2,
                             `2023-2024`= (`mu[7]` + mu_2024)/2) %>%
                      select(`2016-2017`,
                             `2017-2018`,
                             `2018-2019`,
                             `2019-2020`,
                             `2020-2021`,
                             `2021-2022`,
                             `2022-2023`,
                             `2023-2024`) %>%
                      summarise_all(median) %>%
                      unlist()

  CrI_aggr <- samples[,mu_lwr:mu_upr] %>%
                as_tibble() %>%
                mutate(`2016-2017`= (mu_2016 + `mu[1]`)/2,
                       `2017-2018`= (`mu[1]` + `mu[2]`)/2,
                       `2018-2019`= (`mu[2]` + `mu[3]`)/2,
                       `2019-2020`= (`mu[3]` + `mu[4]`)/2,
                       `2020-2021`= (`mu[4]` + `mu[5]`)/2,
                       `2021-2022`= (`mu[5]` + `mu[6]`)/2,
                       `2022-2023`= (`mu[6]` + `mu[7]`)/2,
                       `2023-2024`= (`mu[7]` + mu_2024)/2) %>%
                select(`2016-2017`,
                       `2017-2018`,
                       `2018-2019`,
                       `2019-2020`,
                       `2020-2021`,
                       `2021-2022`,
                       `2022-2023`,
                       `2023-2024`)  %>%
                summarise_all(quantile, probs = c(.025, .975))

  lwr95_aggr <- CrI_aggr %>% slice(1) %>% unlist()
  upr95_aggr <- CrI_aggr %>% slice(2) %>% unlist()
  
  # plot model estimate against NSDUH
  mu_trend_plot <- tibble(pred_mu_aggr = pred_mu_aggr,
         lwr95_aggr = lwr95_aggr,
         upr95_aggr = upr95_aggr,
         year = 2016:2023
  ) %>%
    ggplot() +
    geom_point(aes(x = year, y = S, color = "NSDUH Data"),
               data = tibble(year = c(2016:2018, 2021:2023),
                             S = S)) +
    geom_errorbar(aes(x = year, y = S, ymin = S - 1.96*S.se, ymax = S + 1.96*S.se),
                  width = 0.05,
                  data = tibble(year = c(2016:2018, 2021:2023),
                                S = S,
                                S.se = S.se)) +
    geom_line(aes(x = year, y = pred_mu_aggr, color = "Model")) +
    geom_ribbon(aes(x = year, y = pred_mu_aggr, ymin = lwr95_aggr, ymax = upr95_aggr),
                color = "light blue",
                fill = "light blue",
                alpha = 0.5) +
    labs(x = "Year",
         y = "Statewide prevalence of PWUO-HR") +
    scale_color_manual(name = "", values = c("NSDUH Data" = "black", "Model" = "blue")) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black", size = 8.5))
  
  ggsave("2_yr_mu_trend.png",
         mu_trend_plot,
         device="png",
         path="WAprevalence/output",
         width = 12,
         height = 10,
         units = "cm")

# CREATE PAIRED STATE-LEVEL FIGURE FOR MANUSCRIPT #

  # survey fit alongside the state-wide outcome rates
  state_pair_plot <- plot_grid(mu_trend_plot +
                                 theme(legend.position = c(0.76, 0.88),
                                       legend.background = element_blank()),
                               outcome_trend_plot,
                               labels = c("A)", "B)"),
                               label_size = 11,
                               rel_widths = c(0.44, 0.56))

  ggsave("state_level_pair.png",
         state_pair_plot,
         path = "WAprevalence/output",
         dpi = "retina",
         width = 26,
         height = 10,
         units = "cm",
         bg = "white")

# CREATE COUNTY TREND FIGURE FOR MANUSCRIPT #

  # stack prevalence and the four outcome rates; ED visits are only estimated from 2019
  trend_data <- bind_rows(N_results %>% transmute(county, year, est = mean_prev, panel = "A) Prevalence of higher-risk opioid use"),
                          pmp_results %>% transmute(county, year, est = mean, panel = "B) Buprenorphine rate among PWUO-HR"),
                          death_results %>% transmute(county, year, est = mean, panel = "C) Overdose death rate among PWUO-HR"),
                          ed_results %>% transmute(county, year, est = mean, panel = "D) ED visit rate among PWUO-HR"),
                          hosp_results %>% transmute(county, year, est = mean, panel = "E) Hospitalization rate among PWUO-HR"))

  # label the highest and lowest county in each panel's final year
  trend_labels <- trend_data %>%
                    filter(!is.na(est)) %>%
                    group_by(panel) %>%
                    filter(year == max(year)) %>%
                    filter(est == max(est) | est == min(est)) %>%
                    ungroup()

  # draw the labeled counties in black over the remaining counties in grey
  trend_highlight <- trend_data %>%
                       semi_join(trend_labels, by = c("panel", "county"))

  county_trend_plot <- trend_data %>%
                        ggplot(aes(x = year, y = est, group = county)) +
                        geom_line(linewidth = 0.3, color = "grey75") +
                        geom_line(data = trend_highlight, linewidth = 0.5) +
                        geom_text(aes(label = str_to_title(county)),
                                  data = trend_labels,
                                  hjust = -0.08,
                                  size = 2.3) +
                        facet_wrap(~panel, scales = "free_y", ncol = 3) +
                        coord_cartesian(clip = "off") + # county labels sit outside the panel
                        scale_x_continuous(breaks = seq(2017, 2023, 2),
                                           expand = expansion(mult = c(0.03, 0.32))) +
                        labs(x = "Year", y = NULL) +
                        theme_bw(base_size = 11) +
                        theme(panel.grid = element_blank(),
                              strip.background = element_blank(),
                              strip.text = element_text(color = "black", size = 8.5, hjust = 0),
                              axis.text = element_text(color = "black", size = 8.5),
                              panel.spacing = unit(1.2, "lines"))

  ggsave("trend_spaghetti.png",
         county_trend_plot,
         device = "png",
         path = "WAprevalence/output/maps",
         width = 11,
         height = 6.8,
         bg = "white")

# CREATE RATE DENOMINATOR COMPARISON FIGURE FOR MANUSCRIPT #

  # each outcome rate expressed per 100,000 population and per 100 PWUO-HR
  denominator_data <- bind_rows(pmp_results %>% mutate(outcome = "Buprenorphine"),
                                death_results %>% mutate(outcome = "Overdose death"),
                                ed_results %>% mutate(outcome = "ED visits"),
                                hosp_results %>% mutate(outcome = "Hospitalizations")) %>%
                        left_join(N_results %>% select(county, year, N_est = mean),
                                  by = c("county", "year")) %>%
                        filter(!is.na(mean)) %>%
                        transmute(county,
                                  year,
                                  outcome,
                                  `per 100,000 population 12+` = 10^5*mean*N_est/pop,
                                  `per 100 PWUO-HR` = 100*mean) %>%
                        pivot_longer(starts_with("per"),
                                     names_to = "denominator",
                                     values_to = "rate")

  # panel order: one outcome per row, population denominator on the left
  panel_order <- expand_grid(outcome = c("Buprenorphine", "Overdose death", "ED visits", "Hospitalizations"),
                             denominator = c("per 100,000 population 12+", "per 100 PWUO-HR")) %>%
                   mutate(panel = paste0(LETTERS[1:8], ") ", outcome, " ", denominator))

  denominator_plot <- denominator_data %>%
                        left_join(panel_order, by = c("outcome", "denominator")) %>%
                        mutate(panel = factor(panel, levels = panel_order$panel)) %>%
                        ggplot(aes(x = factor(year), y = rate)) +
                        geom_boxplot(linewidth = 0.3, outlier.size = 0.5, fill = "grey92") +
                        facet_wrap(~panel, scales = "free_y", ncol = 2) +
                        labs(x = "Year", y = NULL) +
                        theme_bw(base_size = 10) +
                        theme(panel.grid = element_blank(),
                              strip.background = element_blank(),
                              strip.text = element_text(color = "black", size = 8.5, hjust = 0),
                              axis.text = element_text(color = "black", size = 7.5),
                              panel.spacing = unit(0.8, "lines"))

  ggsave("rate_denominator_boxplots.png",
         denominator_plot,
         device = "png",
         path = "WAprevalence/output/maps",
         width = 8,
         height = 9.5,
         bg = "white")

# CREATE BISCALE PLOTS #

  # pair prevalence estimates w/ estimates for a given outcome
  prev_outcome_data <- function(outcome_results) {

    N_results %>%
      select(county,
             year,
             mean_prev) %>%
      left_join(outcome_results,
                by = c("county", "year")) %>%
      select(county,
             year,
             mean_prev,
             mean) %>%
      rename(prev_est = mean_prev,
             outcome_est = mean)

  }

  # for a given year, combine model estimates w/ spatial info & apply bi_class
  biscale_data_year <- function(data, year) {

     data_year <- shape_county_WA %>%
            mutate(NAME = tolower(NAME)) %>%
            rename(county = NAME) %>%
            left_join(data[data$year==year,],
                      by = c("county"))

     # years w/ no outcome estimates (ED visits in 2017-2018) are mapped as missing
     if(all(is.na(data_year$outcome_est))) {

       data_year %>% mutate(bi_class = NA)

     } else {

       data_year %>%
            bi_class(x = prev_est,
                     y = outcome_est,
                     style = "quantile",
                     dim = 4)

     }

  }

  # create and save biscale plot of prevalence vs a given outcome
  create_biscale_plot <- function(outcome_results, outcome_label, filename) {

    data <- prev_outcome_data(outcome_results)

    # create biscale class variable for each year and stack
    biscale_data <- biscale_data_year(data, 2017) %>%
                      bind_rows(biscale_data_year(data, 2018),
                                biscale_data_year(data, 2019),
                                biscale_data_year(data, 2020),
                                biscale_data_year(data, 2021),
                                biscale_data_year(data, 2022),
                                biscale_data_year(data, 2023))

    # create biscale plot
    biscale_legend <- bi_legend(pal = "GrPink2",
                                dim = 4,
                                xlab = "Prevalence",
                                ylab = outcome_label,
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

     # legend occupies the empty facet slot
     ggdraw() +
       draw_plot(biscale_map, 0, 0, 1, 1) +
       draw_plot(biscale_legend, 0.77, 0, 0.22, 0.528)

     ggsave(filename,
            device="png",
            path="WAprevalence/output/maps",
            width = 12,
            height = 5,
            units = "cm",
            bg = "white")

  }

  # create biscale plot for each outcome
  create_biscale_plot(pmp_results, "Buprenorphine", "biplot_4dim.png")
  create_biscale_plot(death_results, "Deaths", "biplot_4dim_death.png")
  create_biscale_plot(ed_results, "ED Visits", "biplot_4dim_ed.png")
  create_biscale_plot(hosp_results, "Hospitalizations", "biplot_4dim_hosp.png")

# COMPUTE POSTERIOR PROBABILITY OF GAP MEMBERSHIP BY RANK THRESHOLD (PREVALENCE VS BUPRENORPHINE) #

  # P(county is among the r highest-prevalence and r lowest-buprenorphine counties in a year),
  # classified within each posterior draw and averaged over draws; r = 10 matches the biscale quartiles
  gap_prob_year <- function(t) {

    prev_draws <- sweep(samples[, paste0("N[", (t - 1)*39 + 1:39, "]")], 2,
                        yfit$pop[yfit$year == 2016 + t], "/")
    bup_draws <- samples[, paste0("pi[", (t - 1)*39 + 1:39, ", 1]")]

    # rank each draw once; all thresholds reuse the same ranks
    prev_ranks <- t(apply(prev_draws, 1, rank))
    bup_ranks <- t(apply(bup_draws, 1, rank))

    map_dfr(1:19, function(r) {
      tibble(county = yfit$county[yfit$year == 2016 + t],
             year = 2016 + t,
             r = r,
             p_gap = colMeans(prev_ranks >= 40 - r & bup_ranks <= r))
    })

  }

  gap_probs <- map_dfr(1:7, gap_prob_year)

  write.csv(gap_probs,
            file = "WAprevalence/output/tables/gap_probability_by_threshold.csv",
            row.names = F)

  # print 2023 at the quartile-equivalent threshold for the manuscript
  gap_probs %>%
    filter(year == 2023, r == 10) %>%
    arrange(desc(p_gap)) %>%
    print(n = 10)

# CREATE THRESHOLD CURVES FIGURE FOR SUPPLEMENT #

  # highlight the seven counties with the highest quartile-threshold probabilities in 2023
  gap_2023 <- gap_probs %>% filter(year == 2023)

  gap_top <- gap_2023 %>%
               filter(r == 10) %>%
               slice_max(p_gap, n = 7) %>%
               pull(county)

  # cascade the end labels downward to avoid overlap; leader segments tie each label to its line
  gap_labels <- gap_2023 %>%
                  filter(r == 19, county %in% gap_top) %>%
                  arrange(desc(p_gap)) %>%
                  mutate(lab_y = accumulate(p_gap, function(prev, y) min(y, prev - 0.035)))

  gap_curve_plot <- ggplot(gap_2023, aes(x = r, y = p_gap, group = county)) +
                      geom_vline(xintercept = c(4, 10, 13, 19),
                                 linetype = "dashed",
                                 color = "grey60",
                                 linewidth = 0.3) +
                      geom_line(data = filter(gap_2023, !county %in% gap_top),
                                color = "grey75",
                                linewidth = 0.3) +
                      geom_line(data = filter(gap_2023, county %in% gap_top),
                                linewidth = 0.5) +
                      geom_segment(data = gap_labels,
                                   aes(x = 19, xend = 19.55, y = p_gap, yend = lab_y),
                                   color = "grey50",
                                   linewidth = 0.2) +
                      geom_text(data = gap_labels,
                                aes(x = 19.65, y = lab_y, label = str_to_title(county)),
                                hjust = 0,
                                size = 2.4) +
                      coord_cartesian(clip = "off") +
                      scale_x_continuous(breaks = c(1, 4, 10, 13, 19),
                                         expand = expansion(mult = c(0.02, 0.22)),
                                         sec.axis = dup_axis(breaks = c(4, 8, 10, 13, 19),
                                                             labels = c("deciles", "quintiles", "quartiles", "tertiles", "halves"),
                                                             name = "Equivalent quantile classification")) +
                      scale_y_continuous(limits = c(0, 1)) +
                      labs(x = "Classification threshold r (number of counties)",
                           y = "P(among r highest-prevalence and r lowest-coverage counties)") +
                      theme_bw(base_size = 11) +
                      theme(panel.grid = element_blank(),
                            axis.text = element_text(color = "black", size = 8.5),
                            axis.text.x.top = element_text(size = 8, face = "italic"),
                            axis.title.x.top = element_text(size = 8.5, face = "italic"))

  ggsave("gap_prob_curves_2023.png",
         gap_curve_plot,
         device = "png",
         path = "WAprevalence/output/maps",
         width = 7.5,
         height = 5.5,
         bg = "white")

# CREATE PREVALENCE VS COVERAGE SCATTER #

  # posterior mean prevalence against buprenorphine rate among PWUO-HR, by year;
  # the seven counties with the highest quartile-threshold gap probabilities are highlighted
  scatter_data <- N_results %>%
                    transmute(county, year, prev_pct = 100*mean_prev) %>%
                    left_join(pmp_results %>% transmute(county, year, pmp_pct = 100*mean),
                              by = c("county", "year")) %>%
                    mutate(hl = county %in% gap_top)

  # 2023 labels cascade down the panel's empty right side; leader segments tie them to points
  scatter_labels <- scatter_data %>%
                      filter(year == 2023, hl) %>%
                      arrange(desc(pmp_pct)) %>%
                      mutate(lab_y = accumulate(pmp_pct, function(prev, y) min(y, prev - 4)))

  prev_coverage_plot <- ggplot(scatter_data, aes(x = prev_pct, y = pmp_pct)) +
                          geom_point(data = filter(scatter_data, !hl), color = "grey60", size = 1.0) +
                          geom_point(data = filter(scatter_data, hl), size = 1.4) +
                          geom_segment(data = scatter_labels,
                                       aes(x = prev_pct, xend = 5.1, y = pmp_pct, yend = lab_y),
                                       color = "grey60", linewidth = 0.2) +
                          geom_text(data = scatter_labels,
                                    aes(x = 5.2, y = lab_y, label = str_to_title(county)),
                                    hjust = 0, size = 2.3) +
                          facet_wrap(~year, ncol = 4) +
                          labs(x = "Estimated prevalence of higher-risk opioid use (%)",
                               y = "Estimated buprenorphine receipt rate among PWUO-HR (%)") +
                          theme_bw(base_size = 12) +
                          theme(panel.grid = element_blank(),
                                strip.background = element_blank(),
                                strip.text = element_text(color = "black", size = 9.5, hjust = 0),
                                axis.text = element_text(color = "black", size = 9),
                                panel.spacing = unit(0.9, "lines"))

  ggsave("prev_coverage_scatter.png",
         prev_coverage_plot,
         device = "png",
         path = "WAprevalence/output/maps",
         width = 9,
         height = 5.2,
         bg = "white")

  # print 2023 coverage for the highest-prevalence counties and the county median for the manuscript
  scatter_data %>%
    filter(year == 2023) %>%
    mutate(prev_rank = rank(-prev_pct)) %>%
    filter(prev_rank <= 4) %>%
    arrange(desc(prev_pct)) %>%
    select(county, prev_pct, pmp_pct) %>%
    print()

  scatter_data %>%
    filter(year == 2023) %>%
    summarise(median_pmp_pct = median(pmp_pct)) %>%
    print()

# COMPUTE ESTIMATED NUMBER OF PWUO-HR NOT RECEIVING BUPRENORPHINE #

  # posterior of N minus the observed count of persons receiving buprenorphine;
  # the binomial likelihood guarantees N >= the observed count in every draw
  untreated_draws <- sweep(samples[, paste0("N[", 1:273, "]")], 2, yfit$pmp, "-")

  # same quantity as a proportion of the PWUO-HR population, per draw
  untreated_prop_draws <- untreated_draws/samples[, paste0("N[", 1:273, "]")]

  untreated_csv <- tibble(county = yfit$county,
                          year = yfit$year,
                          mean = colMeans(untreated_draws),
                          lwr95 = apply(untreated_draws, 2, quantile, probs = .025),
                          upr95 = apply(untreated_draws, 2, quantile, probs = .975),
                          prop_mean = colMeans(untreated_prop_draws),
                          prop_lwr95 = apply(untreated_prop_draws, 2, quantile, probs = .025),
                          prop_upr95 = apply(untreated_prop_draws, 2, quantile, probs = .975))

  write.csv(untreated_csv,
            file = "WAprevalence/output/tables/pwuohr_without_bup.csv",
            row.names = F)

  # statewide totals and proportions summed within each draw so the CrIs reflect joint uncertainty
  untreated_state_csv <- map_dfr(1:7, function(t) {

    state_draws <- rowSums(untreated_draws[, (t - 1)*39 + 1:39])
    N_state <- rowSums(samples[, paste0("N[", (t - 1)*39 + 1:39, "]")])

    tibble(year = 2016 + t,
           mean = mean(state_draws),
           lwr95 = quantile(state_draws, .025),
           upr95 = quantile(state_draws, .975),
           prop_mean = mean(state_draws/N_state),
           prop_lwr95 = quantile(state_draws/N_state, .025),
           prop_upr95 = quantile(state_draws/N_state, .975))

  })

  write.csv(untreated_state_csv,
            file = "WAprevalence/output/tables/pwuohr_without_bup_statewide.csv",
            row.names = F)

