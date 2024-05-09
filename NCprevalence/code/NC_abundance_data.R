# IMPORT R PACKAGES
library(readxl) # import excel files
library(writexl) # write excel files
library(tidyverse) # data manipulation and visualization
library(stringr) # tidy string columns
library(reshape2) # melt
library(maptools) # map2SpatialPolygons
library(spdep) # nb2WB
library(rgdal) # import shape files

setwd("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_Abundance_model")

# SPATIAL INFORMATION FOR NC

# read shape files into sf files for use in R
shape_county <- st_read("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/data/Shapefiles/2020 County/cb_2020_us_county_500k.shp", stringsAsFactors = F)
shape_county_NC <- shape_county %>% filter(STUSPS == 'NC')

# compute adjacency matrix for NC counties
NC_map <- readOGR("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/data/Shapefiles/2020 County/cb_2020_us_county_500k.shp")

NC_map <- NC_map[NC_map$STATEFP=='37',]

NC_map <- NC_map[order(NC_map$COUNTYFP),] #convert to data frame

nbmat <- poly2nb(NC_map)

n <- 100

A <- matrix(0,n,n)

for(i in 1:n){
  
  A[i,unlist(nbmat[[i]])]=1
  
}

rownames(A) <- NC_map$NAME
colnames(A) <- NC_map$NAME

num <- colSums(A)

adj <- NULL

for(j in 1:nrow(A)){
  
  adj<-c(adj,which(A[j,]==1))
  
}

adj <- as.vector(adj)

# PREPARE OUTCOME VARIABLES

# outcome variables
OAP <- read.csv("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_abundance_model/data/OAP_2023_05_18.csv")

# what years are variables observed on
OAP %>% group_by(Measure) %>% summarise(range(Year))

pop <- OAP %>%
        filter(Place.Type == "Counties" & Measure == "Incarcerated individuals") %>%
        select(Place, Geoid, Year, Rate.Denom) %>%
        mutate(Rate.Denom = str_replace_all(Rate.Denom, ",", "")) %>%
        mutate(Rate.Denom = as.numeric(Rate.Denom)) %>%
        rename(pop = Rate.Denom)

OAP_outcomes <- OAP %>%
  filter(Measure %in% c("Patients receiving buprenorphine", "Illicit opioid overdose deaths", "Drug overdose ED visits", "People served by treatment programs") &
           Place.Type == "Counties") %>%
  select(Measure, Place, Geoid, Year, rate_label, Rate.Denom, rate_mult, Value.Count, Value.Rate) %>%
  mutate_at(.vars = c("Rate.Denom", "rate_mult", "Value.Count", "Value.Rate"), 
            str_replace_all, 
            pattern = ",",
            replacement = ""
  ) %>%
  mutate_at(.vars = c("Rate.Denom", "rate_mult", "Value.Count", "Value.Rate"), as.numeric)


# 2016-2021 (patients receiving buprenorphine) not-unique therefore Poisson
OAP_outcomes %>%
  filter(Measure == "Patients receiving buprenorphine") %>%
  group_by(Year) %>%
  summarise(n = n())

# 1999-2021 (illicit opioid overdose deaths) unique therefore Binomial
OAP_outcomes %>%
  filter(Measure == "Illicit opioid overdose deaths") %>%
  group_by(Year) %>%
  summarise(n = n())

# 2016-2022 (drug overdose ED visits) # not unique therefore Poisson
OAP_outcomes %>%
  filter(Measure == "Drug overdose ED visits") %>%
  group_by(Year) %>%
  summarise(n = n())

# 2013-2021 (people served by treatment programs): some missing counties, unique therefore Binomial
OAP_outcomes %>%
  filter(Measure == "People served by treatment programs") %>%
  group_by(Year) %>%
  summarise(n = n())

# temporal intersection is 2016-2021, spatial intersection is no less than 95 out of 100 counties across 2016-2021

# investigate which counties are missing in each year for treatment variable
missing_counties <- list()
treatment_years <- 2013:2021

for(i in 1:length(treatment_years)) {
  
  missing_counties[[i]] <- NC_map$NAME %in% OAP[OAP$Measure == "People served by treatment programs"  & OAP$Place.Type == "Counties" &  OAP$Year == treatment_years[i], 2]
  
  names(missing_counties[[i]]) <- NC_map$NAME
  
  
}

names(missing_counties) <- treatment_years
missing_counties

missing_counties <- lapply(missing_counties, function(x) which(x == FALSE)) %>% 
  lapply(names) %>% 
  unlist(use.names = F) %>%
  unique()

# create outcome matrix with missing outcomes for a given year/county given value 0 by pivot_wider; removed `Drug overdose ED visits` as too broad a category
yfit <- OAP_outcomes %>%
              select(Place, Year, Measure, Value.Count) %>%
              pivot_wider(names_from = Measure,
                          values_from = Value.Count,
                          values_fill = NA) %>%
              filter(Year %in% 2016:2021) %>%
              mutate(pop = pop[pop$Year %in% 2016:2021, 4],
                     Geoid = rep(OAP_outcomes[1:100, 3], 6)
              )%>%
              select(`Illicit opioid overdose deaths`, `People served by treatment programs`, `Patients receiving buprenorphine`, Place, Year, pop, Geoid) %>% 
              mutate(D_rate = `Illicit opioid overdose deaths`/pop*10^5,
                     T_rate = `People served by treatment programs`/pop*10^5,
                     B_rate = `Patients receiving buprenorphine`/pop*10^5,
              )

# IMPORT AND CLEAN STATE-LEVEL SURVEY DATA

# state-level survey data (all states)
two_year_state_survey_2015_16 <- read.csv("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_abundance_model/data/state_level_survey/2yr_opioid_misuse_2015-16.csv")
two_year_state_survey_2016_17 <- read.csv("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_abundance_model/data/state_level_survey/2yr_opioid_misuse_2016-17.csv")
two_year_state_survey_2017_18 <- read.csv("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_abundance_model/data/state_level_survey/2yr_opioid_misuse_2017-18.csv")
two_year_state_survey_2018_19 <- read.csv("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_abundance_model/data/state_level_survey/2yr_opioid_misuse_2018-19.csv")
four_year_state_survey_2015_18 <- read.csv("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_abundance_model/data/state_level_survey/4yr_opioid_misuse_2015-18.csv")

# extract NC data and combine different 2 year estimates
two_year_state_survey_2015_16 <- two_year_state_survey_2015_16 %>%
  mutate(years = rep("2015-2016", nrow(two_year_state_survey_2015_16))) %>%
  filter(STATE.FIPS.CODE..NUMERIC. == "37 - North Carolina")

two_year_state_survey_2016_17 <- two_year_state_survey_2016_17 %>%
  mutate(years = rep("2016-2017", nrow(two_year_state_survey_2016_17))) %>%
  filter(STATE.FIPS.CODE..NUMERIC. == "37 - North Carolina")

two_year_state_survey_2017_18 <- two_year_state_survey_2017_18 %>%
  mutate(years = rep("2017-2018", nrow(two_year_state_survey_2017_18))) %>%
  filter(STATE.FIPS.CODE..NUMERIC. == "37 - North Carolina")

two_year_state_survey_2018_19 <- two_year_state_survey_2018_19 %>%
  mutate(years = rep("2018-2019", nrow(two_year_state_survey_2018_19))) %>%
  filter(STATE.FIPS.CODE..NUMERIC. == "37 - North Carolina")

four_year_state_survey_2015_18 <- four_year_state_survey_2015_18 %>%
  mutate(years = rep("2015-2018", nrow(four_year_state_survey_2015_18))) %>%
  filter(STATE.FIPS.CODE..NUMERIC. == "37 - North Carolina")


# extract 4-year estimate
two_year_state_survey <- rbind(two_year_state_survey_2015_16,
                               two_year_state_survey_2016_17,
                               two_year_state_survey_2017_18,
                               two_year_state_survey_2018_19)

two_year_state_survey <- two_year_state_survey %>%
  filter(RC.OPIOIDS...PAST.YEAR.MISUSE == "1 - Yes" | RC.OPIOIDS...PAST.YEAR.MISUSE == "1 - Misused in the past year") %>%
  select(STATE.FIPS.CODE..NUMERIC., years,  Column.., Column...SE)

four_year_state_survey <- four_year_state_survey_2015_18 %>%
  filter(RC.OPIOIDS...PAST.YEAR.MISUSE == "1 - Yes" | RC.OPIOIDS...PAST.YEAR.MISUSE == "1 - Misused in the past year") %>%
  select(STATE.FIPS.CODE..NUMERIC., years,  Column.., Column...SE)

T0 <- 2016 # modeling data from 2016 - 2021

ell.lb <- c(2015,2016,2017,2018,2015)
ell.ub <- c(2016,2017,2018,2019,2018)
ell.lb <- ell.lb - T0 + 1
ell.ub <- ell.ub - T0 + 1
ell.rate <- (ell.ub^02+ell.ub-ell.lb^2+ell.lb)/(2*(ell.ub-ell.lb+1))
S <- c(two_year_state_survey$Column.., four_year_state_survey$Column.. )
S.se <- c(two_year_state_survey$Column...SE, four_year_state_survey$Column...SE)

logit <- function(x) log(x/(1-x)) # logit transformation 
logit_se <- function(prop, prop_se) prop_se/(prop*(1-prop)) # SE of logit transformed variable via delta-approximation

logit_S <- logit(S)
logit_S.se <- logit_se(S, S.se)

# IMPORT AND CLEAN DATA-LEVEL DESIGN MATRICES

# sterile-syringe program (SSP)
#SSP <- read_csv("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_abundance_model/data/SSP.csv")
#SSP <- SSP %>% slice(1:600) %>% select(GEOID, county, year, SSP, `Hepatitis C`)

# HRSA (medically under-served areas)
MUA <- read_excel("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_abundance_model/data/HRSA/MUA_DET.xlsx")
# HRSA (health professional shortage areas (primary care))
HPSA <- read_excel("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline/NC_abundance_model/data/HRSA/BCD_HPSA_FCT_DET_PC.xlsx")

MUA <- MUA %>%
  filter(`State FIPS Code` == "37" & `Medically Underserved Area/Population (MUA/P) Component Geographic Type Code` == "SCTY") %>%
  select(`Designation Date`,
         `MUA/P Status Code`,
         `Break in Designation`,
         `County Equivalent Name`) %>%
  rename(Place = `County Equivalent Name`,
         MUA_status =`MUA/P Status Code`,
         MUA_designation_date = `Designation Date`)

HPSA <- HPSA %>%
  filter(`Common State FIPS Code` == "37" & 
           `Designation Type` == "Geographic HPSA" &
           `HPSA Component Type Description` == "Single County") %>%
  # filter(`County Equivalent Name` %in%  c("Camden", "Hyde", "Onslow", "Person")) %>%
  select(`County Equivalent Name`,
         `HPSA Status Code`, 
         `HPSA Status`,
         `HPSA Designation Date`,
         `HPSA Designation Last Update Date`,
         `Withdrawn Date`) %>%
  filter(`HPSA Status` %in% c("Designated", "Proposed For Withdrawal")) %>%
  rename(Place = `County Equivalent Name`,
         HPSA_status = `HPSA Status Code`)



X.T_data <- data.frame(Place = yfit$Place, Year = yfit$Year)

X.T_data <- X.T_data %>%
  left_join(MUA, by = c("Place")) %>%
  left_join(HPSA, by = c("Place")) %>%
  mutate(MUA = if_else(MUA_status == "D", 1, 0, missing = 0),
         HPSA =  if_else(HPSA_status %in% c("D", "P"), 1, 0, missing = 0)
  ) %>%
  select(Place, Year, MUA, HPSA)


X.T <- model.matrix(~-1 + MUA + HPSA, data = X.T_data)
X.B <- X.T


# HIDTA (DEA: https://www.nhac.org/news/HIDTA_Counties.htm)
HIDTA_counties <- c("Alamance",
                    "Buncombe",
                    "Durham",
                    "Gaston",
                    "Guilford",
                    "Henderson",
                    "Johnston",
                    "McDowell",
                    "Mecklenburg",
                    "Randolph",
                    "Rockingham",
                    "Union",
                    "Wake",
                    "Wayne",
                    "Wilson")

non_MSA_counties <- c("Sampson",
                      "Duplin",
                      "Columbus",
                      "Macon",
                      "Bladen",
                      "Cherokee",
                      "Ashe",
                      "Montgomery",
                      "Hertford",
                      "Caswell",
                      "Martin",
                      "Greene",
                      "Polk",
                      "Warren",
                      "Bertie",
                      "Yancey",
                      "Avery",
                      "Mitchell",
                      "Chowan",
                      "Washington",
                      "Clay",
                      "Alleghany",
                      "Graham",
                      "Hyde",
                      "Tyrrell")

MSA_counties <- setdiff(unique(yfit$Place), non_MSA_counties)

X.D_data <- data.frame(Place = yfit$Place, Year = yfit$Year)
X.D_data <- X.D_data %>%
  mutate(HIDTA = if_else(Place %in% HIDTA_counties, 1, 0),
         MSA = if_else(Place %in% MSA_counties, 1, 0)
  )

X.D <- model.matrix(~-1 + HIDTA + MSA, data = X.D_data)

# IMPORT AND CLEAN AMERICAN COMMUNITY SURVEY (ACS) DATA

clean_ACS <- function(data, county) {
  
  if(county == T) {
    
    results <- data %>%
      filter(SumLev == "050") 
    
  }
  
  else { 
    
    results <- data %>%
      filter(SumLev == "040", AreaName == "North Carolina", Geocomp == "00")
    
  }
  
  results <- results %>%
    select(years, 
           geoid, 
           CivLabForce, CivLabForce_moe, 
           EmployedCLF, EmployedCLF_moe,
           TotPop, TotPop_moe, 
           Poor, Poor_moe, 
           Bachelorsormore, Bachelorsormore_moe,
           Over25, Over25_moe) %>%
    rename(Year = years,
           poor = Poor,
           poor_moe = Poor_moe) %>%
    mutate_at(.vars = c("CivLabForce", "CivLabForce_moe",
                        "EmployedCLF", "EmployedCLF_moe",
                        "TotPop", "TotPop_moe",
                        "poor", "poor_moe",
                        "Bachelorsormore", "Bachelorsormore_moe",
                        "Over25", "Over25_moe"),
              str_replace_all,
              pattern = ",",
              replacement = "") %>%
    mutate_at(.vars = c("CivLabForce", "CivLabForce_moe",
                        "EmployedCLF", "EmployedCLF_moe",
                        "TotPop", "TotPop_moe",
                        "poor", "poor_moe",
                        "Bachelorsormore", "Bachelorsormore_moe",
                        "Over25", "Over25_moe"),
              as.numeric)
  
  return(results)
  
}

# import ACS data 
for(i in 2012:2021) {
  data <- read_csv(paste("./data/NC/NC_5year_", i, ".csv", sep = ""))
  assign(paste("NC_5year_", i, sep = ""), data)
}

for(i in c(2012:2019, 2021)){
  data <- read_csv(paste("./data/NC/NC_1year_", i, ".csv", sep = ""))
  assign(paste("NC_1year_", i, sep = ""), data)
}

# change sumlev and areaname column names to SumLev to be consistent with other year's naming conventions
NC_1year_2021 <- NC_1year_2021 %>% rename(SumLev = sumlev, AreaName = areaname, Geocomp = geocomp)
NC_5year_2020 <- NC_5year_2020 %>% rename(SumLev = sumlev, AreaName = areaname, Geocomp = geocomp)
NC_5year_2021 <- NC_5year_2021 %>% rename(SumLev = sumlev, AreaName = areaname, Geocomp = geocomp)

# change Years column in NC_5_x to years to match the NC_1_x column naming and its use in the clean_ACS function
dataset_names <- paste0("NC_5year_", 2012:2019)

# iterate over each dataset name and rename the column
purrr::map(dataset_names, function(ds_name) {
  # Get the dataset from the name
  dataset <- get(ds_name)
  
  # Rename the Years column to years
  renamed_dataset <- dataset %>%
    rename(years = Years,)
  
  # Assign the renamed dataset back to the global environment
  assign(ds_name, renamed_dataset, envir = .GlobalEnv)
})

# remove quotations from AreaName in NC_1year_2016
NC_1year_2016 <- NC_1year_2016 %>% mutate(AreaName = str_replace_all(AreaName, '"', ""))

# clean 1-year county data
years <- setdiff(2012:2021, 2020)  # This generates the sequence 2012 to 2021 but skips 2020.

original_dataset_names <- paste0("NC_1year_", years)

# 2. Generate the new dataset names
new_dataset_names <- paste0("county_1year_", years)

# 3. Apply the function to each dataset
for (i in seq_along(original_dataset_names)) {
  
  # Get the dataset from the global environment
  dataset <- get(original_dataset_names[i], envir = .GlobalEnv)
  
  # Apply the function
  cleaned_dataset <- clean_ACS(dataset, county = T)
  
  # Assign the cleaned dataset to the new name
  assign(new_dataset_names[i], cleaned_dataset, envir = .GlobalEnv)
}

# clean 1-year state data
original_dataset_names <- paste0("NC_1year_", years)

# 2. Generate the new dataset names
new_dataset_names <- paste0("state_1year_", years)

# 3. Apply the function to each dataset
for (i in seq_along(original_dataset_names)) {
  
  # Get the dataset from the global environment
  dataset <- get(original_dataset_names[i], envir = .GlobalEnv)
  
  # Apply the function
  cleaned_dataset <- clean_ACS(dataset, county = F)
  
  # Assign the cleaned dataset to the new name
  assign(new_dataset_names[i], cleaned_dataset, envir = .GlobalEnv)
}

# clean 5-year county data
years <- 2012:2021

original_dataset_names <- paste0("NC_5year_", years)

# 2. Generate the new dataset names
new_dataset_names <- paste0("county_5year_", years)

# 3. Apply the function to each dataset
for (i in seq_along(original_dataset_names)) {
  
  # Get the dataset from the global environment
  dataset <- get(original_dataset_names[i], envir = .GlobalEnv)
  
  # Apply the function
  cleaned_dataset <- clean_ACS(dataset, county = T)
  
  # Assign the cleaned dataset to the new name
  assign(new_dataset_names[i], cleaned_dataset, envir = .GlobalEnv)
}

# combine years
countydata1 <- county_1year_2012 %>% bind_rows(county_1year_2013,
                                               county_1year_2014,
                                               county_1year_2015,
                                               county_1year_2016,
                                               county_1year_2017,
                                               county_1year_2018,
                                               county_1year_2019,
                                               county_1year_2021)

countydata5 <- county_5year_2012 %>% mutate(Year = rep(2012, 100)) %>%
  bind_rows(county_5year_2013,
            county_5year_2014,
            county_5year_2015,
            county_5year_2016,
            county_5year_2017,
            county_5year_2018,
            county_5year_2019,
            county_5year_2020,
            county_5year_2021)

statedata1 <- state_1year_2012 %>% bind_rows(state_1year_2013,
                                             state_1year_2014,
                                             state_1year_2015,
                                             state_1year_2016,
                                             state_1year_2017,
                                             state_1year_2018,
                                             state_1year_2019,
                                             state_1year_2021)

# add counties & years with missing data to the data frames above

# extract NC county geoids
county_geoids <- paste0("05000US", unique(yfit$Geoid), sep = "")

# create empty data frame of years and counties to merge with years/counties we have data for
test <- data.frame(geoid = rep(county_geoids, 10), Year = rep(2012:2021, each = 100))

countydata1 <- test %>% left_join(countydata1, by = c("geoid", "Year"))
statedata1 <- statedata1 %>% add_row(Year = 2020, geoid = "04000US37", .before = 9)

# compute rates for the ACS outcomes and their corresponding standard errors

# compute standard error of rate where numerator=x ,denominator=y
rate_se <- function(x, y, x_moe, y_moe) (1/y)*sqrt(x_moe^2 - (x^2/y^2)*y_moe^2)/1.645

clean_ACS2 <- function(data) { data %>% mutate(UnemployedCLF = CivLabForce - EmployedCLF,
                                               UnemploymentRate = UnemployedCLF/CivLabForce,
                                               UnemploymentRate_se = rate_se(UnemployedCLF, CivLabForce, EmployedCLF_moe, CivLabForce_moe),
                                               PovertyRate = poor/TotPop,
                                               PovertyRate_se = rate_se(poor, TotPop, poor_moe, TotPop_moe),
                                               BachelorsormoreRate = Bachelorsormore/Over25,
                                               BachelorsormoreRate_se = rate_se(Bachelorsormore, Over25, Bachelorsormore_moe, Over25_moe),
                                               logit_UnemploymentRate = logit(UnemploymentRate),
                                               logit_UnemploymentRate_se = logit_se(UnemploymentRate, UnemploymentRate_se),
                                               logit_PovertyRate = logit(PovertyRate),
                                               logit_PovertyRate_se = logit_se(PovertyRate, PovertyRate_se),
                                               logit_BachelorsormoreRate = logit(BachelorsormoreRate),
                                               logit_BachelorsormoreRate_se = logit_se(BachelorsormoreRate, BachelorsormoreRate_se)) %>%
    mutate_at(.vars = c("UnemploymentRate", "UnemploymentRate_se",
                        "PovertyRate", "PovertyRate_se",
                        "BachelorsormoreRate", "BachelorsormoreRate_se"),
              function(x) x*100) %>% # convert proportions to percentages, note that the logit transformed variables remain proportions
    select(geoid, Year, 
           UnemploymentRate, UnemploymentRate_se,
           PovertyRate, PovertyRate_se,
           BachelorsormoreRate, BachelorsormoreRate_se,
           logit_UnemploymentRate, logit_UnemploymentRate_se,
           logit_PovertyRate, logit_PovertyRate_se,
           logit_BachelorsormoreRate, logit_BachelorsormoreRate_se,
    )
}

countydata1 <- clean_ACS2(countydata1)
countydata5 <- clean_ACS2(countydata5)
statedata1 <- clean_ACS2(statedata1)

# which rows of W1obs and W1.state have data
ony_ind <- which(is.na(countydata1[, 3]) == F)
ony_ind_state <- which(is.na(statedata1[, 3]) == F)

# check that missing years and counties are consistent across the included ACS variables
unemployment_ony_ind <- which(is.na(countydata1[, 3]) == F) 
poverty_ony_ind <- which(is.na(countydata1[, 5]) == F)
bachelorsormore_ony_ind <- which(is.na(countydata1[, 7]) == F)

all(unemployment_ony_ind == poverty_ony_ind) 
all(unemployment_ony_ind == bachelorsormore_ony_ind) 
all(bachelorsormore_ony_ind == poverty_ony_ind) 

unemployment_se_ony_ind <- which(is.na(countydata1[, 4]) == F) 
poverty_se_ony_ind <- which(is.na(countydata1[, 6]) == F)
bachelorsormore_se_ony_ind <- which(is.na(countydata1[, 8]) == F)

all(unemployment_se_ony_ind == poverty_se_ony_ind) 
all(unemployment_se_ony_ind == bachelorsormore_se_ony_ind) 
all(bachelorsormore_se_ony_ind == poverty_se_ony_ind) 

# create data matrices for logit of ACS variables and their standard errors
logit_W1.obs <- countydata1 %>% select(logit_UnemploymentRate, logit_PovertyRate, logit_BachelorsormoreRate) %>% drop_na
logit_W1.se <- countydata1 %>% select(logit_UnemploymentRate_se, logit_PovertyRate_se, logit_BachelorsormoreRate_se) %>% drop_na

logit_W5.obs <- countydata5 %>% filter(Year >= 2016) %>% select(logit_UnemploymentRate, logit_PovertyRate, logit_BachelorsormoreRate)
logit_W5.se <- countydata5 %>% filter(Year >= 2016) %>% select(logit_UnemploymentRate_se, logit_PovertyRate_se, logit_BachelorsormoreRate_se) 

logit_W1.state <- statedata1 %>% select(logit_UnemploymentRate, logit_PovertyRate, logit_BachelorsormoreRate) %>% drop_na
logit_W1.state.se <- statedata1 %>% select(logit_UnemploymentRate_se,logit_PovertyRate_se, logit_BachelorsormoreRate_se) %>% drop_na

# create initial value matrix for logit_W
logit_W_init <- countydata5 %>% select(logit_UnemploymentRate, logit_PovertyRate, logit_BachelorsormoreRate)

# SAVE PREPARED DATA FOR USE IN NIMBLE MODEL
save(adj, num,
     yfit,
     X.D, X.T, X.B,
     S, S.se, logit_S, logit_S.se, ell.rate,
     logit_W1.obs, logit_W5.obs, logit_W1.state,
     logit_W1.se, logit_W5.se, logit_W1.state.se,
     ony_ind, ony_ind_state, logit_W_init,
     file = "./data/data_for_analysis.Rda")
