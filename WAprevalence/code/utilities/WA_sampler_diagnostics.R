# throwaway convergence check for the 3-chain run; run after the model (uses `samples`)
# if starting fresh: load("WAprevalence/output/mcmc/<run>.Rda") first
library(coda)

load("WAprevalence/output/mcmc/test.Rda")

T   <- length(grep("^mu\\[", varnames(samples)))
beta_cols <- grep("^beta\\[", varnames(samples), value = TRUE)
beta_cols <- beta_cols[apply(as.matrix(samples[, beta_cols]), 2, var) > 0]  # drop the inert ED early-year cells
beta_cols <- beta_cols[-which(is.na(beta_cols))]

# effective sample size, pooled across chains
print(round(effectiveSize(samples[, beta_cols])))

# R-hat (needs the multi-chain mcmc.list)
if(nchain(samples) > 1) print(round(gelman.diag(samples[, beta_cols], multivariate = FALSE)$psrf, 3))

# trace + density, one color per chain
plot(samples[, beta_cols])

