# WAprevalence

Analysis code and data for "Estimating Higher-Risk Opioid Use with a Bayesian Spatiotemporal Model to Inform MOUD Planning in Washington State": a Bayesian spatiotemporal integrated abundance model estimating the annual number and prevalence of people who use opioids in ways that increase health risks (PWUO-HR) in Washington State's 39 counties, 2017-2023.

## Code

Run in order:

1. `WA_abundance_data.R` builds the analysis data set: county outcome counts, NSDUH state-level survey estimates, and Census population denominators (pulled via tidycensus).
2. `WA_abundance_model.R` specifies the nimble model and runs the MCMC; `WA_abundance_model_parallel.R` runs the three chains in parallel.
3. `WA_abundance_output.R` checks convergence and produces the posterior summaries, figures, and tables reported in the manuscript and supplement.

## Data

- `nsduh/` NSDUH two-year state-level estimates of past-year opioid misuse (six survey periods, 2016-2024); files with an `_OPIIMFNMYR` suffix are the IMF-inclusive variant used only by the sensitivity analysis
- `population/` WA OFM small area population estimates
- `shape_county_WA.Rda` county shapefile
- County-level outcome counts (buprenorphine receipt, overdose deaths, ED visits, hospitalizations) are from the WA DOH overdose dashboard (https://doh.wa.gov/data-and-statistical-reports/washington-tracking-network-wtn/opioids/overdose-dashboard) and are not stored in this repository.

## Output

`WA_abundance_output.R` writes MCMC diagnostics, figures, and tables to `output/`.
