# throwaway sampler diagnostics; run after WA_abundance_model.R, then delete
library(coda)

sm <- as.matrix(samples)
T <- length(grep("^mu\\[", colnames(sm)))    # years
R <- length(grep("^N\\[", colnames(sm))) / T # counties

# across-time correlation of beta, within each outcome
cat("\n== across-time corr of beta (per outcome) ==\n")
for(j in c(1,2,4)){ cat("outcome", j, "\n"); print(round(cor(sm[, paste0("beta[", 1:T, ", ", j, "]")]), 2)) }
cat("outcome 3\n"); print(round(cor(sm[, paste0("beta[", 3:T, ", 3]")]), 2))

# across-outcome correlation of beta, within each year
cat("\n== across-outcome corr of beta (per year) ==\n")
for(t in 1:T){
  cols <- if(t >= 3) 1:4 else c(1,2,4)
  cat("year", t, "\n"); print(round(cor(sm[, paste0("beta[", t, ", ", cols, "]")]), 2))
}

# beta vs statewide mu and total N (detection/abundance ridge)
cat("\n== cor(beta, mu) and cor(beta, total N), by year ==\n")
for(t in 1:T){
  Ntot <- rowSums(sm[, paste0("N[", (t-1)*R + 1:R, "]")])
  cols <- if(t >= 3) 1:4 else c(1,2,4)
  for(j in cols){
    b <- sm[, paste0("beta[", t, ", ", j, "]")]
    cat(sprintf("year %d outcome %d:  mu %5.2f   N %5.2f\n",
                t, j, cor(b, sm[, paste0("mu[", t, "]")]), cor(b, Ntot)))
  }
}

# effective sample size
cat("\n== ESS: beta ==\n");    print(round(effectiveSize(samples[, grep("^beta\\[", colnames(sm))])))
cat("\n== ESS: beta.mu ==\n"); print(round(effectiveSize(samples[, c("beta.mu[1]", "beta.mu[2]")])))
cat("\n== ESS: mu ==\n");      print(round(effectiveSize(samples[, grep("^mu\\[", colnames(sm))])))
