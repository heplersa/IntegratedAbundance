# throwaway convergence check for the 3-chain run; run after the model (uses `samples`), then delete
# if starting fresh: load("WAprevalence/output/mcmc/<run>.Rda") first
library(coda)

T   <- length(grep("^mu\\[", varnames(samples)))
beta_cols <- grep("^beta\\[", varnames(samples), value = TRUE)
beta_cols <- beta_cols[apply(as.matrix(samples[, beta_cols]), 2, var) > 0]  # drop the inert ED early-year cells
key <- c("beta.mu[1]", "beta.mu[2]", beta_cols, paste0("mu[", 1:T, "]"))

# effective sample size, pooled across chains
print(round(effectiveSize(samples[, key])))

# R-hat (needs the multi-chain mcmc.list)
if(nchain(samples) > 1) print(round(gelman.diag(samples[, key], multivariate = FALSE)$psrf, 3))

# trace + density, one color per chain
plot(samples[, key])
