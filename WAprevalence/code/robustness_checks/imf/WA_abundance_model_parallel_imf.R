# PARALLEL 3-CHAIN RUN of WA_abundance_model_imf.R
# each worker builds its own compiled model and runs one chain. launch from the repo root.
# needs >= 3 cores and ~3x the RAM of a single run.

library(parallel)
library(coda)

wd <- getwd()

run_chain <- function(chain, wd){
  setwd(wd)
  BUILD_ONLY <- TRUE                                             # source build-only (skip the script's own run)
  source("WAprevalence/code/robustness_checks/imf/WA_abundance_model_imf.R", local = TRUE) # -> compiled_mcmc, inits_list, MCS
  runMCMC(compiled_mcmc,
          inits   = inits_list[[chain]],
          niter   = MCS,
          nburnin = MCS/2,
          thin    = 50,
          setSeed = chain,
          samplesAsCodaMCMC = TRUE,
          progressBar = FALSE)
}

cl <- makeCluster(3)
chains <- parLapply(cl, 1:3, run_chain, wd = wd)
stopCluster(cl)

samples <- do.call(mcmc.list, chains)

dir.create("WAprevalence/output/mcmc", recursive = TRUE, showWarnings = FALSE)
save(samples, file = "WAprevalence/output/mcmc/MCMC_N_pois_3_chains_2_mill_each_imf_2026_08_14.Rda")
