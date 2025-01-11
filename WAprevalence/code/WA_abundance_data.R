# PRE-PROCESS WASHINGTON STATE DATA FOR USE IN INTEGRATED ABUNDANCE MODEL #
# BRIAN N. WHITE #
# 2025-01-10 #

# IMPORT R PACKAGES
library(tidyverse) # data manipulation and visualization
library(janitor) # clean_names
library(spdep) # poly2nb
library(tigris) # pull WA shape files from US Census
library(tidycensus) # pull pop data from US Census

# SPATIAL INFORMATION FOR WA

# pull and save county shape files
# source: tigris package which accesses https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html
  
  # un-comment out these two lines if you need to re-pull the shape file
  # shape_county_WA <- counties(state = "WA", year = 2024)
  # save(shape_county_WA, file = "WAprevalence/data/shape_county_WA.Rda")
  load("WAprevalence/data/shape_county_WA.Rda")
  
  # compute adjacency matrix for NC counties
  WA_map <- shape_county_WA[order(shape_county_WA$COUNTYFP),] #convert to data frame
  
  nbmat <- poly2nb(WA_map)
  
  n <- length(shape_county_WA$NAME) # should be 39
  
  A <- matrix(0,n,n)
  
  for(i in 1:n){
    
    A[i,unlist(nbmat[[i]])]=1
    
  }
  
  rownames(A) <- WA_map$NAME
  colnames(A) <- WA_map$NAME
  
  num <- colSums(A)
  
  adj <- NULL
  
  for(j in 1:nrow(A)){
    
    adj<-c(adj,which(A[j,]==1))
    
  }
  
  adj <- as.vector(adj)

# PREPARE OUTCOME VARIABLES 

# load county level population estimates from the WA state office of financial management (OFM); available from 2010-2022
# source: https://ofm.wa.gov/washington-data-research/population-demographics/population-estimates/small-area-estimates-program

  # import raw data
  WA_county_pop_raw <- read.csv("WAprevalence/data/population/saep_county20.csv")
  # process raw data
  WA_county_pop_processed <- WA_county_pop_raw %>%
                                slice(1:39) %>%
                                select(1, 4:13, 15:19) %>%
                                rename(county = `County.Name`) %>%
                                pivot_longer(cols = 2:16,
                                             names_to = c("year"),
                                             values_to = c("pop")) %>%
                                mutate(pop = str_remove_all(pop, ","),
                                       year = str_remove(year, "Estimated.Total.Population."),
                                       year = str_remove(year, "OFM.Adjusted.Total.Population.")) %>%
                                mutate(across(c(year, pop),
                                              as.numeric)) %>%
                                mutate(county = str_to_lower(county)) %>%
                                filter(year %in% 2017:2022)

# outcome variables: 3 outcomes (pmp, death, OUD), 6 years (2017-2022), 39 counties
# source: confidential data; pulled by Dave Kline Jan 2025
  
  # import raw data
  outcomes_raw <- read.csv("WAprevalence/data/outcomes/Single_year_crc_file.csv")
  # process raw data
  outcomes_processed <- outcomes_raw %>%
                          group_by(year, county) %>%
                          summarise(pmp = sum(pmp),
                                    death = sum(death),
                                    OUD = sum(OUD)) %>%
                          mutate(county = str_to_lower(county)) %>%
                          filter(county != "unknown") %>%
                          left_join(WA_county_pop_processed,
                                    by = c("county", "year"))
  
  # check that there is no missing data
  apply(outcomes_processed, 2, function(x) sum(is.na(x)))
  # check that all counties-years are present; 6 years x 39 counties = 234 rows
  nrow(outcomes_processed) == 6*39
  
  # rename for use in Bayesian model
  yfit <- outcomes_processed

# PREPARE STATE-LEVEL SURVEY DATA FOR PREVALENCE OF OPIOID MISUSE WITHIN LAST YEAR

# load state-level survey data from NSDUH using the SAMHSA datatools web application.
# source: https://datatools.samhsa.gov/

  # pull raw data  
  nsduh_2016_2017_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2016_2017.csv")
  nsduh_2017_2018_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2017_2018.csv")
  nsduh_2018_2019_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2018_2019.csv")
  nsduh_2021_2022_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2021_2022.csv")
  
  # process raw data
  process_nsduh <- function(nsduh_raw, year) {
    
          nsduh_processed <- nsduh_raw %>%
                                clean_names() %>%
                                filter(state_fips_code_numeric == "53 - Washington" &
                                       rc_opioids_past_year_misuse %in% c("1 - Yes" , "1 - Misused in the past year")
                                       ) %>%
                                mutate(year = year) %>%
                                       rename(prev = total,
                                              prev_se = total_se) %>%
                                       select(year,
                                              prev,
                                              prev_se)
          
          return(nsduh_processed)
          
  }
  
  nsduh_2016_2017_processed <- process_nsduh(nsduh_2016_2017_raw, "2016-2017")
  nsduh_2017_2018_processed <- process_nsduh(nsduh_2017_2018_raw, "2017-2018")
  nsduh_2018_2019_processed <- process_nsduh(nsduh_2018_2019_raw, "2018-2019")
  nsduh_2021_2022_processed <- process_nsduh(nsduh_2021_2022_raw, "2021-2022")
  
  nsduh_processed <- nsduh_2016_2017_processed %>%
                        bind_rows(nsduh_2017_2018_processed,
                                  nsduh_2018_2019_processed,
                                  nsduh_2021_2022_processed)
  
  # rename for Bayesian model
  S <- nsduh_processed$prev
  S.se <- nsduh_processed$prev_se
  
# compute logit of state-level survey prevalence
  
  # logit transformation & SE of logit transformation via delta method
  logit <- function(x) log(x/(1-x)) # logit transformation 
  logit_se <- function(prop, prop_se) prop_se/(prop*(1-prop)) # SE of logit transformed variable via delta-approximation
  
  logit_S <- logit(S)
  logit_S.se <- logit_se(S, S.se)

# compute intermediate values required to compute mean in normal model for the state-level prevalence
# see section 3.1.2 of model paper for reference

  # modeling data from 2017 - 2022
  T0 <- 2017 

  # compute slope (ell.rate) for the mean in the normal model
  ell.lb <- c(2016,2017,2018,2021)
  ell.ub <- c(2017,2018,2019,2022)
  ell.lb <- ell.lb - T0 + 1
  ell.ub <- ell.ub - T0 + 1
  ell.rate <- (ell.ub^02+ell.ub-ell.lb^2+ell.lb)/(2*(ell.ub-ell.lb+1))

# SAVE PREPARED DATA FOR USE IN NIMBLE MODEL
save(adj, num,
     yfit,
     S, S.se, logit_S, logit_S.se, ell.rate,
     file = "WAprevalence/data/data_for_analysis.Rda")
