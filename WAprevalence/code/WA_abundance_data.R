# PRE-PROCESS WASHINGTON STATE DATA FOR USE IN INTEGRATED ABUNDANCE MODEL #
# BRIAN N. WHITE #
# 2025-01-17 #

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
                                filter(year %in% 2017:2023)

# compute county-level population aged 12+ to match NSDUH survey universe
# approach: compute 12+ age share from ACS 5-year estimates, apply to OFM true yearly totals
# source: US Census Bureau ACS via tidycensus; table B01001 (Sex by Age)

  # ACS variables needed to compute 12+ age share
  # verify codes with: load_variables(2023, "acs5") %>% filter(str_detect(name, "B01001_0"))
  age_vars <- c(total = "B01001_001",
                male_under5 = "B01001_003",
                male_5to9 = "B01001_004",
                male_10to14 = "B01001_005",
                female_under5 = "B01001_027",
                female_5to9 = "B01001_028",
                female_10to14 = "B01001_029")

  pull_acs_age_share <- function(year) {

          get_acs(geography = "county",
                  state = "WA",
                  variables = age_vars,
                  year = year,
                  survey = "acs5",
                  output = "wide") %>%
            mutate(year = year) %>%
            rename(county = NAME) %>%
            mutate(county = str_remove(county, " County, Washington"),
                   county = str_to_lower(county)) %>%
            # compute 12+ share; subtract under-5, 5-9, and 2/5 of 10-14
            # (assumes uniform age distribution within 10-14 bin for ages 10-11)
            mutate(pop_under12 = (male_under5E + male_5to9E + (2/5)*male_10to14E) +
                                 (female_under5E + female_5to9E + (2/5)*female_10to14E),
                   share_12plus = (totalE - pop_under12) / totalE) %>%
            select(year,
                   county,
                   share_12plus)

  }

  acs_age_share <- map_dfr(2017:2023, pull_acs_age_share)

  # apply ACS 12+ share to OFM yearly totals; overwrite pop to match NSDUH 12+ universe
  WA_county_pop_processed <- WA_county_pop_processed %>%
                                left_join(acs_age_share,
                                          by = c("year", "county")) %>%
                                mutate(pop = round(pop * share_12plus)) %>%
                                select(county, year, pop)

# outcome variables: 3 outcomes (pmp, death, OUD), 6 years (2017-2022), 39 counties
# source: confidential data; pulled by Dave Kline Jan 2025
  
  # import raw data
  
    # hospitalization and ED visit outcomes
    any_opioid_overdose_ED <- read.csv("WAprevalence/data/outcomes/Overdose_Downloadable_ED.csv")
    any_opioid_overdose_hospitalization <- read.csv("WAprevalence/data/outcomes/Overdose_Downloadable_Hospitalizations.csv")
  
    # summer 2025 update; expand study period to 2017-2023 (so add 2023)
    outcomes_raw <- read.csv("WAprevalence/data/outcomes/final_county_data_single_year_dedupe_with_unknowns.csv")
  
  # process raw data
    
    # clean and combine hospitalization and ED visit outcomes
    any_opioid_overdose_ED_clean <- any_opioid_overdose_ED %>%
      filter(Time.Breakdown %in% 2017:2023 &
             Drug.Category == "Any Opioid" &
             Demographic.Category == "Overall",
             Geography == "County",
             !(Location %in% c("Benton-Franklin", "Chelan-Douglas", "Northeast Tri County"))) %>%
      rename(year = Time.Breakdown,
             county = Location,
             ed = ED.Visit.Count) %>%
      mutate(ed = if_else(ed == "*", NA, as.numeric(ed)),
             county = tolower(str_replace_all(county, " County", ""))) %>%
      select(year,
             county,
             ed)
    
    any_opioid_overdose_hospitalization_clean <- any_opioid_overdose_hospitalization %>%
      filter(Year %in% 2017:2023 &
             Time.Aggregation == "1 year rolling counts" &
             Drug.Category == "Any Opioid" &
             Demographic.Category == "Overall" &
             Geography == "County" &
             !(Location %in% c("Benton-Franklin", "Chelan-Douglas", "Northeast Tri County", "Unassigned Region", "Unassigned County"))) %>%
      rename(year = Year,
             county = Location,
             hosp = Hospitalization.Count) %>%
      mutate(hosp = if_else(hosp == "*", NA, as.numeric(hosp)),
             county = tolower(str_replace_all(county, " County", ""))) %>%
      select(year,
             county,
             hosp)

    # extract marginal county by year counts for pmp_oud and death_oud using latest data; merge in hospitalization and ED visit outcomes
    outcomes_processed <- outcomes_raw %>%
            filter(final_county != "Unknown") %>%
            group_by(final_county, year) %>% 
            summarise(pmp = sum(pmp_oud),
                      death = sum(death_oud)) %>%
            rename(county = final_county) %>%
            mutate(county = tolower(county)) %>%
            left_join(WA_county_pop_processed,
                      by = c("year", "county")) %>%
            left_join(any_opioid_overdose_ED_clean,
                      by = c("year", "county")) %>%
            right_join(any_opioid_overdose_hospitalization_clean,
                      by = c("year", "county")) %>%
            filter(year >= 2017) %>%
            arrange(year, county)
  
  # check that there is no missing data in marginal outcomes; for ed and hospitalization counts <5 are censored; account for this in the model
  apply(outcomes_processed[, c("pmp", "death", "ed", "hosp")], 2, function(x) sum(is.na(x))) == c(0, 0, 0, 0)
  # check that all counties-years are present; 6 years x 39 counties = 234 rows
  nrow(outcomes_processed) == 7*39

  # rename for use in Bayesian model
  yfit <- outcomes_processed

# PREPARE STATE-LEVEL SURVEY DATA FOR PREVALENCE OF OPIOID MISUSE WITHIN LAST YEAR

# load state-level survey data from NSDUH using the SAMHSA datatools web application.
# source: https://datatools.samhsa.gov/

  # pull raw data; no restricted use 2022-2023 data available as of early August 2025 (i.e. with the state level 2-year data)
  nsduh_2016_2017_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2016_2017.csv")
  nsduh_2017_2018_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2017_2018.csv")
  nsduh_2018_2019_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2018_2019.csv")
  nsduh_2021_2022_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2021_2022.csv")
  nsduh_2022_2023_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2022_2023.csv")
  
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
  nsduh_2022_2023_processed <- process_nsduh(nsduh_2022_2023_raw, "2022-2023")
  
  nsduh_processed <- nsduh_2016_2017_processed %>%
                        bind_rows(nsduh_2017_2018_processed,
                                  nsduh_2018_2019_processed,
                                  nsduh_2021_2022_processed,
                                  nsduh_2022_2023_processed)
  
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

  # modeling data from 2017 - 2023
  T0 <- 2017 

  # compute linear (ell.rate) and quadratic (ell.rate2) factors for the mean in the normal model
  ell.lb <- c(2016,2017,2018,2021, 2022)
  ell.ub <- c(2017,2018,2019,2022, 2023)
  ell.lb <- ell.lb - T0 + 1
  ell.ub <- ell.ub - T0 + 1
  ell.rate <- (ell.ub + ell.lb)/2
  ell.rate2 <- (ell.ub*(ell.ub + 1)*(2*ell.ub + 1) - ell.lb*(ell.lb - 1)*(2*ell.lb - 1))/(6*(ell.ub - ell.lb + 1))

# SAVE PREPARED DATA FOR USE IN NIMBLE MODEL
save(adj, num,
     yfit,
     S, S.se, logit_S, logit_S.se, ell.rate, ell.rate2,
     file = "WAprevalence/data/data_for_analysis.Rda")