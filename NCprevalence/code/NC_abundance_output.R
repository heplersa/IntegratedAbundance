## EXAMINE NC ABUNDANCE MODEL MCMC CONVERGENCE & GENERATE PAPER FIGURES/TABLES ##
## 2024-04-05 ##

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
setwd("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_abundance_model/code")
load("../data/data_for_analysis.Rda") # data_for_analysis.Rda file generated via NC_abundance_data.R file.

# IMPORT MCMC OUTPUT FROM MODEL
load("../output/model/full/runs/NC_abundance_output_full_logit_S_2024_02_12.Rda")

# EXAMINE MCMC CONVERGENCE
MCMCvis::MCMCtrace(samples, params = paste0("pi[", sample(1:600, 20), ", 3]"), ISB = F, filename = "B", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("N[", sample(1:600, 20), "]"), ISB = F, filename = "N", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("lambda[", sample(1:600, 20), "]"), ISB = F, filename = "lambda", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("W_ilogit[", sample(1:1000, 20), ", 1]"), ISB = F, filename = "W_ilogit_unemployment", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("W_ilogit[", sample(1:1000, 20), ", 2]"), ISB = F, filename = "W_ilogit_poverty", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("W_ilogit[", sample(1:1000, 20), ", 3]"), ISB = F, filename = "W_ilogit_bachelorsormore", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("g[", sample(1:600, 20), ", 1]"), ISB = F, filename = "g_unemployment", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("g[", sample(1:600, 20), ", 2]"), ISB = F, filename = "g_poverty", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("g[", sample(1:600, 20), ", 3]"), ISB = F, filename = "g_bachelorsormore", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:600, 20), ", 1]"), ISB = F, filename = "f_D", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:600, 20), ", 2]"), ISB = F, filename = "f_T", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("f[", sample(1:600, 20), ", 3]"), ISB = F, filename = "f_B", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("v[", sample(1:600, 20), "]"), ISB = F, filename = "v", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = paste0("u[", sample(1:600, 20), "]"), ISB = F, filename = "u", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "gamma", filename = "gamma", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "beta.D", filename = "beta.D", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "beta.T", filename = "beta.T", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "beta.B", filename = "beta.B", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "beta", filename = "beta", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "beta.mu", filename = "beta.mu", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "mu", filename = "mu", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "mu.W", filename = "mu.W", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "mu.W0", filename = "mu.W0", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau", filename = "tau", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.g", filename = "tau.g", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.f", filename = "tau.f", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.W", filename = "tau.W", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.W0", filename = "tau.W0", wd = "../output/model/full/diagnostics")
MCMCvis::MCMCtrace(samples, params = "tau.u", filename = "tau.u", wd = "../output/model/full/diagnostics")

# EXTRACT POSTERIOR MEANS, 95% CrI (QUANTILES), SD AND NEW GR DIAGNOSTIC STAT
results <- list(colMeans(samples),
                apply(samples,2,
                      quantile,probs=c(.025,.975)),
                apply(samples,2,sd), 
                apply(samples, 2, function(x) stable.GR(x, multivariate = F)$psrf))

    # specify indices of parameters of interest
    D_lwr <- which(names(results[[1]])=="pi[1, 1]")
    D_upr <- which(names(results[[1]])=="pi[600, 1]")
    T_lwr <- which(names(results[[1]])=="pi[1, 2]")
    T_upr <- which(names(results[[1]])=="pi[600, 2]")
    B_lwr <- which(names(results[[1]])=="pi[1, 3]")
    B_upr <- which(names(results[[1]])=="pi[600, 3]")
    N_lwr <- which(names(results[[1]]) == "N[1]")
    N_upr <- which(names(results[[1]]) == "N[600]")
    lambda_lwr <- which(names(results[[1]]) == "lambda[1]")
    lambda_upr <- which(names(results[[1]]) == "lambda[600]")
    beta_lwr <- which(names(results[[1]]) == "beta[1, 1]")
    beta_upr <- which(names(results[[1]]) == "beta[6, 3]")
    mu_lwr <- which(names(results[[1]]) == "mu[1]")
    mu_upr <- which(names(results[[1]]) == "mu[6]")
    beta_D_lwr <- which(names(results[[1]]) == "beta.D[1]")
    beta_D_upr <- which(names(results[[1]]) == "beta.D[2]")
    beta_T_lwr <- which(names(results[[1]]) == "beta.T[1]")
    beta_T_upr <- which(names(results[[1]]) == "beta.T[2]")
    beta_B_lwr <- which(names(results[[1]]) == "beta.B[1]")
    beta_B_upr <- which(names(results[[1]]) == "beta.B[2]")
    gamma_lwr <- which(names(results[[1]]) == "gamma[1]")
    gamma_upr <- which(names(results[[1]]) == "gamma[3]")
    W_ilogit_unemployment_lwr <- which(names(results[[1]])=="W_ilogit[401, 1]") # get estimates from 2016-2021, W starts in 2012
    W_ilogit_unemployment_upr <- which(names(results[[1]])=="W_ilogit[1000, 1]")
    W_ilogit_poverty_lwr <- which(names(results[[1]])=="W_ilogit[401, 2]")
    W_ilogit_poverty_upr <- which(names(results[[1]])=="W_ilogit[1000, 2]")
    W_ilogit_bachelorsormore_lwr <- which(names(results[[1]])=="W_ilogit[401, 3]")
    W_ilogit_bachelorsormore_upr <- which(names(results[[1]])=="W_ilogit[1000, 3]")
    
    
    # create tidy data sets of estimates merged with corresponding spatio-temporal data
    results_to_tibble <- function(results, par) {
      
                                    par_lwr <- get(paste(par, "_lwr", sep = ""), envir = .GlobalEnv)
                                    par_upr <- get(paste(par, "_upr", sep = ""), envir = .GlobalEnv)
      
                                    tibble(par = names(results[[1]][par_lwr:par_upr]),
                                           geoid = yfit$Geoid,
                                           county = yfit$Place,
                                           year = yfit$Year,
                                           pop = yfit$pop,
                                           mean = results[[1]][par_lwr:par_upr],
                                           lwr95 = results[[2]][1, par_lwr:par_upr],
                                           upr95 = results[[2]][2, par_lwr:par_upr],
                                           sd = results[[3]][par_lwr:par_upr]
                                           )
          
    }
    
    D_results <- results_to_tibble(results, "D")
    T_results <- results_to_tibble(results, "T")
    B_results <- results_to_tibble(results, "B")
    lambda_results <- results_to_tibble(results, "lambda") %>% mutate(CrI = case_when(
                                                                                            lwr95 > 1  ~ "95% CrI > 1",
                                                                                            upr95 < 1  ~ "95% CrI < 1",
                                                                                            .default = "95% CrI contains 1"),
                                                                      CrI = fct_relevel(CrI, c("95% CrI < 1", "95% CrI contains 1", "95% CrI > 1"))
    )
    N_results <- results_to_tibble(results, "N") %>% mutate(prev = mean/pop)
    W_ilogit_unemployment_results <- results_to_tibble(results, "W_ilogit_unemployment")
    W_ilogit_poverty_results <- results_to_tibble(results, "W_ilogit_poverty")
    W_ilogit_bachelorsormore_results <- results_to_tibble(results, "W_ilogit_bachelorsormore")
    
# CREATE CHOROPLETH MAPS
    
    # load shape files for NC counties
    shape_county <- st_read("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/data/Shapefiles/2020 County/cb_2020_us_county_500k.shp", stringsAsFactors = F)
    shape_county_NC <- shape_county %>% filter(STUSPS == 'NC') 
    # shape_county_NC <- counties(state = "NC") # if you need NC county shape files, can pull from internet using tigris package
    
    # function to create choropleth maps
    create_choropleth_map <- function(data, value, colorbar_type = NULL, colorbar_title = NULL) {
      
                                          p <- shape_county_NC %>%
                                                  mutate(geoid = as.numeric(GEOID)) %>%
                                                  left_join(data, by = c("geoid")) %>%
                                                  ggplot() +
                                                  geom_sf(aes(fill = {{value}})) 
                                          
                                          # is color bar monotonic or diverging?
                                          p <- {if(colorbar_type == "monotonic"){
                                            
                                                      p + scale_fill_gradient(low = "white",
                                                                              high = "blue",
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
                                          p <- p + labs(fill = colorbar_title) +
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
    D_map <- create_choropleth_map(data = D_results, value = mean, colorbar_type = "monotonic")
    T_map <- create_choropleth_map(data = T_results, value = mean, colorbar_type = "monotonic")
    B_map <- create_choropleth_map(data = B_results, value = mean, colorbar_type = "monotonic")
    lambda_map <- create_choropleth_map(data = lambda_results, value = mean, colorbar_type = "diverging")
    lambda_CrI_map <- create_choropleth_map(data = lambda_results, value = CrI, colorbar_type = "other") + theme(legend.position = "right")
    N_map <- create_choropleth_map(data = N_results, value = prev, colorbar_type = "monotonic")
    W_ilogit_unemployment_map <- create_choropleth_map(data = W_ilogit_unemployment_results, value = mean, colorbar_type = "monotonic")
    W_ilogit_poverty_map <- create_choropleth_map(data = W_ilogit_poverty_results, value = mean, colorbar_type = "monotonic")
    W_ilogit_bachelorsormore_map <- create_choropleth_map(data = W_ilogit_bachelorsormore_results, value = mean, colorbar_type = "monotonic")
    
    ggsave(filename = "D.png", 
           plot = D_map, 
           path = "../output/model/full/maps", 
           bg = "White",
           dpi = "retina",
           height = 3,
           width = 10)
    
    ggsave(filename = "T.png", 
           plot = T_map, 
           path = "../output/model/full/maps", 
           bg = "White",
           dpi = "retina",
           height = 3,
           width = 10)
    
    ggsave(filename = "B.png", 
           plot = B_map, 
           path = "../output/model/full/maps", 
           bg = "White",
           dpi = "retina",
           height = 3,
           width = 10)
    
    ggsave(filename = "lambda.png", 
           plot = lambda_map, 
           path = "../output/model/full/maps", 
           bg = "White",
           dpi = "retina",
           height = 3,
           width = 10)
    
    ggsave(filename = "lambda_CrI.png", 
           plot = lambda_CrI_map, 
           path = "../output/model/full/maps", 
           bg = "White",
           dpi = "retina",
           height = 3,
           width = 10)
    
    ggsave(filename = "N.png", 
           plot = N_map, 
           path = "../output/model/full/maps", 
           bg = "White",
           dpi = "retina",
           height = 3,
           width = 10)
    
    ggsave(filename = "W_ilogit_unemployment.png", 
           plot = W_ilogit_unemployment_map, 
           path = "../output/model/full/maps", 
           bg = "White",
           dpi = "retina",
           height = 3,
           width = 10)
    
    ggsave(filename = "W_ilogit_poverty.png", 
           plot = W_ilogit_poverty_map, 
           path = "../output/model/full/maps", 
           bg = "White",
           dpi = "retina",
           height = 3,
           width = 10)
    
    ggsave(filename = "W_ilogit_bachelorsormore.png", 
           plot = W_ilogit_bachelorsormore_map, 
           path = "../output/model/full/maps", 
           bg = "White",
           dpi = "retina",
           height = 3,
           width = 10)
    
# CREATE COMPARISON MAP OF OBSERVED DEATH RATE FOR POP VS FOR PWMO IN 2021
    
      # compute county-level death rate for pop and for PWMO
      yfit <- yfit %>% mutate(N = results[[1]][N_lwr:N_upr],
                              D_rate_pop = `Illicit opioid overdose deaths`/pop,
                              D_rate_PWMO = `Illicit opioid overdose deaths`/N)
    
      # compute the statewide observed death rate for pop and for PWMO in 2021
      state_wide_death_rate <- yfit %>%
                                  filter(Year == 2021) %>%
                                  summarise(pop_rate = sum(`Illicit opioid overdose deaths`)/sum(pop),
                                            PWMO_rate = sum(`Illicit opioid overdose deaths`)/sum(N)
                                  )
      
      # create choropleth maps
      
        # pop rate
         D_pop_map <-  shape_county_NC %>%
                            mutate(Geoid = as.numeric(GEOID)) %>%
                            left_join(yfit, by = c("Geoid")) %>%
                            filter(Year == 2021) %>%
                            ggplot() +
                            geom_sf(aes(fill = D_rate_pop)) +
                            scale_fill_gradient2(low = "blue",
                                                 mid = "white",
                                                 high = "red",
                                                 midpoint = state_wide_death_rate$pop_rate,
                                                 guide = guide_colorbar(barheight = 8)) +
                            labs(fill = NULL,
                                 title = "Overdose deaths per population") +
                            theme_map() +
                            theme(legend.position = "right") +
                            theme(strip.background = element_rect(fill = "white", color = NA),
                                  strip.text = element_text(color = "black",
                                                            size = 12, 
                                                            hjust = 0),
                                  legend.text = element_text(size = 12),
                                  legend.title = element_text(size = 12),
                                  plot.title = element_text(hjust = 0.7)
                            )
                
        # PWMO rate
        D_PWMO_map <- shape_county_NC %>%
                          mutate(Geoid = as.numeric(GEOID)) %>%
                          left_join(yfit, by = c("Geoid")) %>%
                          filter(Year == 2021) %>%
                          ggplot() +
                          geom_sf(aes(fill = D_rate_PWMO)) +
                          scale_fill_gradient2(low = "blue",
                                               mid = "white",
                                               high = "red",
                                               midpoint = state_wide_death_rate$PWMO_rate,
                                               guide = guide_colorbar(barheight = 8)) +
                          labs(fill = NULL,
                               title = "Overdose deaths per PWMO") +
                          theme_map() +
                          theme(legend.position = "right") +
                          theme(strip.background = element_rect(fill = "white", color = NA),
                                strip.text = element_text(color = "black",
                                                          size = 12, 
                                                          hjust = 0),
                                legend.text = element_text(size = 12),
                                legend.title = element_text(size = 12),
                                plot.title = element_text(hjust = 0.7)
                          )
          
      
        # make it so that scientific notation is disabled, undo this via option(scipen=0)
        options(scipen = 999)
        plot_grid(D_pop_map, D_PWMO_map,
                  align = "hv",
                  nrow = 1,
                  ncol = 2)
        
        ggsave(filename = "2021.png", 
               path = "../output/model/full/maps", 
               bg = "White",
               dpi = "retina",
               height = 3,
               width = 11)
    
# EXAMINE BOX-PLOTS AND TABLES OF OUTCOME RATES ACROSS COUNTIES FOR EACH YEAR

    # boxplot
    tibble("Unintentional overdose death" = results[[1]][D_lwr:D_upr],
           "Treatment program use" = results[[1]][T_lwr:T_upr],
           "Buprenorphine reception" = results[[1]][B_lwr:B_upr],
           year = as.factor(rep(2016:2021, each = 100))
    ) %>%
      select("Unintentional overdose death", year) %>%
      filter(year %in% c(2016, 2021)) %>%
      melt(id.vars = c("year"), value.name = "mean", variable.name = "outcome") %>%
      group_by(outcome, year) %>%
      summarise(median = median(mean)) %>%
      flextable() %>%
      save_as_docx(path = "../output/model/full/outcome_medians.docx")
    
     D_boxplot <- tibble(D_mean = results[[1]][D_lwr:D_upr],
                         T_mean = results[[1]][T_lwr:T_upr],
                         B_mean = results[[1]][B_lwr:B_upr],
                         year = as.factor(rep(2016:2021, each = 100))
                  ) %>%
                    ggplot(aes(x=year, y=D_mean)) +
                    geom_boxplot(outlier.fill = NA,
                                 outlier.shape = 21) +
                    labs(x = "Year",
                         y = "Rate",
                         title = "Overdose death") +
                    theme(plot.title = element_text(hjust = 0.5, size = 11))
     
     T_boxplot <- tibble(D_mean = results[[1]][D_lwr:D_upr],
                         T_mean = results[[1]][T_lwr:T_upr],
                         B_mean = results[[1]][B_lwr:B_upr],
                         year = as.factor(rep(2016:2021, each = 100))
                   ) %>%
                     ggplot(aes(x=year, y=T_mean)) +
                     geom_boxplot(outlier.fill = NA,
                                  outlier.shape = 21) +
                     labs(x = "Year",
                          y = NULL,
                          title = "Served by treatment program") +
                    theme(plot.title = element_text(hjust = 0.5, size = 11))
     
     B_boxplot <- tibble(D_mean = results[[1]][D_lwr:D_upr],
                         T_mean = results[[1]][T_lwr:T_upr],
                         B_mean = results[[1]][B_lwr:B_upr],
                         year = as.factor(rep(2016:2021, each = 100))
                   ) %>%
                     ggplot(aes(x=year, y=B_mean)) +
                     geom_boxplot(outlier.fill = NA,
                                  outlier.shape = 21) +
                     labs(x = "Year",
                          y = NULL,
                          title = "Received buprenorphine") +
                     theme(plot.title = element_text(hjust = 0.5, size = 11))
     
     plot_grid(D_boxplot, T_boxplot, B_boxplot,
               align = "hv",
               nrow = 1,
               ncol = 3)

# EXAMINE HOW MANY YEARS/COUNTIES WERE CENSORES (NA) FOR TREATMENT PROGRAM USE
yfit %>%
  group_by(Year) %>%
  summarise(no_na = sum(is.na(`People served by treatment programs`)),
            prop_na = no_na/n(),
            n = n()) %>%
  flextable() %>%
  save_as_docx(path = "../output/exploratory/missing_treatment_table.docx")

# EXAMINE ESTIMATED YEARLY STATE-WIDE RATE OF PWMO VS DATA
tibble(pred_mu = results[[1]][mu_lwr:mu_upr],
       lwr95 = results[[2]][1, mu_lwr:mu_upr],
       upr95 = results[[2]][2, mu_lwr:mu_upr],
       year = 2016:2021) %>%
      ggplot() +
      geom_point(aes(x = year, y = S),
                 data = tibble(year = c(2016:2019, 2018), # 2018 corrponds to the 2015-2018 4-year estimate
                               S = S)) +
      geom_errorbar(aes(x = year, y = S, ymin = S - 1.96*S.se, ymax = S + 1.96*S.se),
                    width = 0.05,
                    data = tibble(year = c(2016:2019, 2018), # 2018 corrponds to the 2015-2018 4-year estimate
                                  S = S,
                                  S.se = S.se)) +
      geom_line(aes(x = year, y = pred_mu), color = "blue") +
      geom_ribbon(aes(x = year, y = pred_mu, ymin = lwr95, ymax = upr95),
                  color = "light blue",
                  fill = "light blue",
                  alpha = 0.5) +
      labs(x = "Year",
           y = "Statewide Prevalence of PWMO") +
      theme_classic()


# EXAMINE ESTIMATED TIME-VARYING INTERCEPTS FOR EACH OUTCOME
    
    # log-odds scale
    tibble(pred_beta = results[[1]][beta_lwr:beta_upr],
           lwr95 = results[[2]][1, beta_lwr:beta_upr],
           upr95 = results[[2]][2, beta_lwr:beta_upr],
           year = rep(2016:2021, 3),
           outcome = rep(c("Illicit opioid overdose deaths",
                           "Patients served by treatment program",
                           "Patients receiving buprenorphine"), each = 6)
        ) %>%
           ggplot(aes(x = year, y = pred_beta, fill = outcome)) +
           geom_point(aes(color = outcome)) +
           geom_errorbar(aes(ymin = lwr95, ymax = upr95, color = outcome),
                        width = 0.05) +
           geom_line(linetype = "dashed", aes(color = outcome)) +
           geom_ribbon(aes(ymin = lwr95, ymax = upr95), alpha = 0.1) +
           theme_classic() +
           #theme(legend.position = c(.15, .4)) +
           labs(color = "Outcome",
                fill = "Outcome",
                x = "Year",
                y = "Estimated intercept (log-odds)")
    
    # odds scale
    tibble(pred_beta = colMeans(exp(samples[,beta_lwr:beta_upr])),
           lwr95 = apply(exp(samples[,beta_lwr:beta_upr]), 2, quantile, probs = 0.025),
           upr95 = apply(exp(samples[,beta_lwr:beta_upr]), 2, quantile, probs = 0.975),
           year = rep(2016:2021, 3),
           outcome = rep(c("Illicit opioid overdose deaths",
                           "Patients served by treatment program",
                           "Patients receiving buprenorphine"), each = 6)
    ) %>%
           ggplot(aes(x = year, y = pred_beta, fill = outcome)) +
           geom_point(aes(color = outcome)) +
           geom_errorbar(aes(ymin = lwr95, ymax = upr95, color = outcome),
                         width = 0.05) +
           geom_line(linetype = "dashed", aes(color = outcome)) +
           geom_ribbon(aes(ymin = lwr95, ymax = upr95), alpha = 0.1) +
           theme_classic() +
           #theme(legend.position = c(.15, .4)) +
           labs(color = "Outcome",
                fill = "Outcome",
                x = "Year",
                y = "Estimated intercept (odds)")

# CREATE TABLE OF FIXED EFFECT ESTIMATES
    
    # outcome level fixed effects
    tibble(odds = colMeans(exp(samples[,beta_B_lwr:(beta_B_lwr + 5)])),
           lwr95 = apply(exp(samples[,beta_B_lwr:(beta_B_lwr + 5)]), 2, quantile, probs=c(.025)),
           upr95 = apply(exp(samples[,beta_B_lwr:(beta_B_lwr + 5)]), 2, quantile, probs=c(.975)),
           covariate = c("MUA", "HPSA", "HIDTA", "MSA", "MUA", "HPSA"),
           outcome = rep(c("Buprenorphine reception", "Illicit opioid overdose death", "Treatment program use"), each = 2))%>%
      mutate(`Odds (95% CrI)` = paste(round(odds, 2), " (", round(lwr95, 2), ",", round(upr95, 2), ")", sep ="")) %>%
      select(outcome, covariate, `Odds (95% CrI)`) %>%
      rename(Outcome = outcome,
             Covariate = covariate) %>%
      flextable() %>%
      theme_box() %>%
      save_as_docx(path = "../output/model/full/beta_fixed.docx")

    # process level fixed effects (i.e. fixed effects for the model for lambda)
    tibble(odds = colMeans(exp(samples[,gamma_lwr:gamma_upr])),
           lwr95 = apply(exp(samples[,gamma_lwr:gamma_upr]), 2, quantile, probs=c(.025)),
           upr95 = apply(exp(samples[,gamma_lwr:gamma_upr]), 2, quantile, probs=c(.975)),
           covariate = c("Unemployment rate", "Poverty rate", "Educational attainment rate (Bachelor's or more)"))%>%
      mutate(`Relative Risk Ratio (95% CrI)` = paste(round(odds, 2), " (", round(lwr95, 2), ",", round(upr95, 2), ")", sep ="")) %>%
      select(covariate, `Relative Risk Ratio (95% CrI)`) %>%
      rename(Covariate = covariate) %>%
      flextable() %>%
      theme_box() %>%
      save_as_docx(path = "../output/model/full/gamma_fixed.docx")
    
# EXAMINE POSTERIOR STATE PREVELANCE
    
  # extract MCMC draws for N in 2021
  N_2021 <- samples[, (N_upr-99):N_upr] %>% as_tibble()
  
  # compute posterior distribution for sum of N in 2021
  N_sum_post <- apply(N_2021, 1, sum) %>% as_tibble()
            
  # examine the posterior distribution
  N_sum_post %>%
    ggplot(aes(x=value)) +
    geom_histogram(color = "black",
                   fill = "light blue",
                   alpha = 0.5) +
    geom_vline(xintercept = mean(N_sum_post$value), 
               linetype = "solid",
               color = "black") +
    theme_classic()
  
  # extract posterior mean
  N_sum_post_mean <- mean(N_sum_post$value)
  N_sum_post_CrI <- quantile(N_sum_post$value, prob = c(.025, .975))
  
  # extract NC pop in 2021
  NC_2021_pop <- yfit %>%
                    filter(Year == 2021) %>%
                    summarise(NC_pop = sum(pop)) %>%
                    pull()
  
  # compute posterior mean of prevelance of PWMO in NC in 2021
  N_sum_post_mean/NC_2021_pop*100
