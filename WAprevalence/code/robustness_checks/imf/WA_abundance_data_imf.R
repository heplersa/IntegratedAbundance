# PRE-PROCESS WASHINGTON STATE DATA FOR USE IN INTEGRATED ABUNDANCE MODEL #
# IMF SENSITIVITY ANALYSIS #
# BRIAN N. WHITE #
# 2025-01-17 #
#
# standalone copy of WA_abundance_data.R, kept separate so the sensitivity analysis is
# fully contained. The ONLY substantive difference from the primary analysis is the NSDUH
# input: the 2022-2023 and 2023-2024 survey estimates use OPIIMFNMYR (opioid misuse including
# illicitly manufactured fentanyl) in place of OPINMYR. OPIIMFNMYR is not available for the
# earlier survey periods, so 2016-2017 through 2021-2022 are unchanged. The model itself,
# the outcome data, and the population denominators are identical.
# Keep in sync with WA_abundance_data.R if that file changes.

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
  
  # compute adjacency matrix for WA counties
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

  # load county-level population estimates aged 12+ from Census Population Estimates Program
  # source: US Census Bureau PEP via tidycensus
  # approach: subtract under-12 pop from total, assuming uniformity within 10-14 age bin for ages 10-11

    # pull vintage 2019 PEP data; DATE codes 10-12 = July 1 estimates for 2017-2019
    # source: https://www.census.gov/data/developers/data-sets/popest-popproj/popest/popest-vars/2019.html
    pep_2017_2019 <- get_estimates(geography = "county",
                                   product = "characteristics",
                                   breakdown = "AGEGROUP",
                                   breakdown_labels = TRUE,
                                   state = "WA",
                                   time_series = TRUE,
                                   vintage = 2019) %>%
      filter(DATE %in% 10:12) %>%
      mutate(year = DATE - 10 + 2017) %>%
      select(-DATE)

    # pull vintage 2023 PEP data; year column gives 2020-2023 directly
    pep_2020_2023 <- get_estimates(geography = "county",
                                   product = "characteristics",
                                   breakdown = "AGEGROUP",
                                   breakdown_labels = TRUE,
                                   state = "WA",
                                   time_series = TRUE,
                                   vintage = 2023)

    # combine vintages, compute 12+ population, clean county names
    WA_county_pop_processed <- bind_rows(pep_2017_2019, pep_2020_2023) %>%
      filter(AGEGROUP %in% c("All ages", "Age 0 to 4 years",
                              "Age 5 to 9 years", "Age 10 to 14 years")) %>%
      mutate(county = str_remove(NAME, " County, Washington"),
             county = str_to_lower(county)) %>%
      select(county, year, AGEGROUP, value) %>%
      pivot_wider(names_from = AGEGROUP, values_from = value) %>%
      mutate(pop = round(`All ages` - `Age 0 to 4 years` - `Age 5 to 9 years` -
                           (2/5) * `Age 10 to 14 years`)) %>%
      select(county, year, pop)

# outcome variables: 4 outcomes (pmp, death, ED visit, hospitalization), 7 years (2017-2023), 39 counties
# sources: all four outcomes are published on the public WA DOH overdose dashboard, which suppresses counts under 10
# (https://doh.wa.gov/data-and-statistical-reports/washington-tracking-network-wtn/opioids/overdose-dashboard)
# ed and hosp counts here are downloads from that dashboard (pulled by BNW)
# pmp and death counts here come from a direct DOH data transfer with unsuppressed counts, pulled by Dave Kline Jan 2025
  
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
  
  # check that there is no missing data in marginal outcomes; for ed and hospitalization, counts under 10 are suppressed at source (raw value "*") and are censored in [1,9]; account for this in the model
  apply(outcomes_processed[, c("pmp", "death", "ed", "hosp")], 2, function(x) sum(is.na(x))) == c(0, 0, 0, 0)
  # check that all counties-years are present; 6 years x 39 counties = 234 rows
  nrow(outcomes_processed) == 7*39

  # rename for use in Bayesian model
  yfit <- outcomes_processed

# PREPARE STATE-LEVEL SURVEY DATA FOR PREVALENCE OF OPIOID MISUSE WITHIN LAST YEAR

# load state-level survey data from NSDUH using the SAMHSA datatools web application.
# source: https://datatools.samhsa.gov/

  # pull raw data; 2016-2017 through 2021-2022 are past year opioid misuse (OPINMYR),
  # 2022-2023 and 2023-2024 are the IMF-inclusive variant (OPIIMFNMYR), which SAMHSA does not
  # publish for the earlier survey periods
  # no two-year datasets span 2020 (NSDUH trend break), hence the gap between 2018-2019 and 2021-2022
  nsduh_2016_2017_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2016_2017.csv")
  nsduh_2017_2018_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2017_2018.csv")
  nsduh_2018_2019_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2018_2019.csv")
  nsduh_2021_2022_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2021_2022.csv")
  nsduh_2022_2023_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2022_2023_OPIIMFNMYR.csv")
  nsduh_2023_2024_raw <- read.csv("WAprevalence/data/nsduh/nsduh_2yr_2023_2024_OPIIMFNMYR.csv")

  # process raw data
  # the SAMHSA export is not consistent across survey periods: the misuse variable is column 1
  # in the newer pulls and column 2 in the older ones, its column name varies with the variable
  # label, and its "misused" level is labelled "1 - Yes", "1 - Misused in the past year", or
  # "1 - Misused within the past year". Identify the column positionally and the level by its "1" code.
  process_nsduh <- function(nsduh_raw, year) {

          nsduh_clean <- clean_names(nsduh_raw)

          misuse_var <- setdiff(names(nsduh_clean)[1:2], "state_fips_code_numeric")

          nsduh_processed <- nsduh_clean %>%
                                filter(state_fips_code_numeric == "53 - Washington" &
                                       str_detect(.data[[misuse_var]], "^1 - ")
                                       ) %>%
                                mutate(year = year) %>%
                                       rename(prev = total,
                                              prev_se = total_se) %>%
                                       select(year,
                                              prev,
                                              prev_se)

          # each export should yield exactly one WA row for the "misused" level
          stopifnot(length(misuse_var) == 1, nrow(nsduh_processed) == 1)

          return(nsduh_processed)

  }

  nsduh_2016_2017_processed <- process_nsduh(nsduh_2016_2017_raw, "2016-2017")
  nsduh_2017_2018_processed <- process_nsduh(nsduh_2017_2018_raw, "2017-2018")
  nsduh_2018_2019_processed <- process_nsduh(nsduh_2018_2019_raw, "2018-2019")
  nsduh_2021_2022_processed <- process_nsduh(nsduh_2021_2022_raw, "2021-2022")
  nsduh_2022_2023_processed <- process_nsduh(nsduh_2022_2023_raw, "2022-2023")
  nsduh_2023_2024_processed <- process_nsduh(nsduh_2023_2024_raw, "2023-2024")

  nsduh_processed <- nsduh_2016_2017_processed %>%
                        bind_rows(nsduh_2017_2018_processed,
                                  nsduh_2018_2019_processed,
                                  nsduh_2021_2022_processed,
                                  nsduh_2022_2023_processed,
                                  nsduh_2023_2024_processed)
  
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

  # compute survey period midpoints (ell.rate) at which the linear trend is evaluated.
  # the 2023-2024 survey straddles the end of the study period; its midpoint (7.5) is still
  # a valid evaluation point for the linear predictor, since mu is a function of beta.mu
  # rather than a node that must exist for every year covered by a survey
  ell.lb <- c(2016, 2017, 2018, 2021, 2022, 2023)
  ell.ub <- c(2017, 2018, 2019, 2022, 2023, 2024)
  ell.lb <- ell.lb - T0 + 1
  ell.ub <- ell.ub - T0 + 1
  ell.rate <- (ell.ub + ell.lb)/2

# SAVE PREPARED DATA FOR USE IN NIMBLE MODEL
save(adj, num,
     yfit,
     S, S.se, logit_S, logit_S.se, ell.rate,
     file = "WAprevalence/data/data_for_analysis_imf.Rda")