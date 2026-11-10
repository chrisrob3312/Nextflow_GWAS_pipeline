#!/usr/bin/env Rscript
# ============================================================================
# test_tractor_model.R: the score test in bin/tractor_model.R, without GENESIS
# ============================================================================
#   Rscript tests/test_tractor_model.R          (needs Matrix and survival)
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(survival) })
here <- local({ f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(f)) dirname(normalizePath(f[1])) else "tests" })
source(file.path(here, "..", "bin", "tractor_model.R"))
pass <- 0; fail <- 0
check <- function(name, ok) { if (isTRUE(ok)) { pass <<- pass + 1; cat("ok  ", name, "\n") } else { fail <<- fail + 1; cat("FAIL", name, "\n") } }

# --- minor allele counts --------------------------------------------------------------------
check("mac_minor: ALT 30 on 32 haplotypes is a minor allele count of 2", mac_minor(30, 32) == 2)
check("mac_minor: ALT 3 on 32 haplotypes is 3", mac_minor(3, 32) == 3)

set.seed(7)
n <- 200
time <- rexp(n, 0.1); event <- rbinom(n, 1, 0.7); age <- runif(n, 1, 18)
ops0 <- cox_ops(time, event, rep(0, n), NULL)                 # null without covariates: eta = 0

# Two ancestries, AFR and EUR (reference EUR); everyone below has two EUR haplotypes unless said
A <- c("AFR", "EUR"); zero <- rep(0, n); two <- rep(2, n)
tv <- function(dose_eur, hap_eur = two, dose_afr = zero, mac_min = 10, ops = ops0)
    test_variant_tractor(list(AFR = 2 - hap_eur), list(AFR = dose_afr, EUR = dose_eur), list(AFR = 2 - hap_eur, EUR = hap_eur),
                         ops, A, "EUR", mac_min = mac_min)
dose30 <- c(rep(2, 15), rep(0, n - 15))                         # ALT 30
r <- tv(dose30)
check("ALT 30 on 400 EUR haplotypes: MAC_EUR 30, tested", r$mac[["EUR"]] == 30 && r$status == "ok" && !is.na(r$marginal$EUR$beta))
hap16 <- c(rep(2, 16), rep(0, n - 16))                         # only 16 people carry EUR haplotypes: 32 of them
r <- tv(dose30, hap_eur = hap16)
check("ALT 30 on 16 people's 32 EUR haplotypes: MAC_EUR 2, not tested", r$mac[["EUR"]] == 2 && is.na(r$marginal$EUR$beta))
check("... no background left to test: status mac (MAC_STRATUM is 30)", r$status == "mac" && r$mac_stratum == 30)
r <- tv(c(1, 1, rep(0, n - 2)))
check("stratum-wide minor allele count 2 < 10: status mac, MAC_STRATUM 2", r$status == "mac" && r$mac_stratum == 2)

# --- the Cox score test equals survival's (one testable background, constant local ancestry) ---
g <- rbinom(n, 2, 0.3)
r <- tv(g, mac_min = 1)
sc <- unname(summary(coxph(Surv(time, event) ~ g, ties = "breslow"))$sctest[1])
check(sprintf("Cox score statistic %.6f equals coxph's score test %.6f", r$joint$stat, sc), abs(r$joint$stat - sc) < 1e-6)
check("one testable background: P_JOINT equals P_EUR (1 df)", r$joint$df == 1 && abs(r$joint$p - r$marginal$EUR$p) < 1e-12)

# with a covariate: the efficient score equals coxph's score test for g at the covariate-only fit
fit0 <- coxph(Surv(time, event) ~ age, ties = "breslow")
ops1 <- cox_ops(time, event, as.numeric(fit0$linear.predictors), matrix(age, ncol = 1))
r <- tv(g, mac_min = 1, ops = ops1)
f1 <- coxph(Surv(time, event) ~ age + g, ties = "breslow", init = c(coef(fit0), 0), control = coxph.control(iter.max = 0))
check(sprintf("with a covariate: score statistic %.6f equals coxph's %.6f", r$joint$stat, f1$score), abs(r$joint$stat - f1$score) < 1e-6)

# --- k = 2: joint, per-ancestry, heterogeneity -----------------------------------------------
la <- rbinom(n, 2, 0.5)                                         # AFR haplotypes per person
dA <- rbinom(n, la, 0.4); dE <- rbinom(n, 2 - la, 0.2)
r <- test_variant_tractor(list(AFR = la), list(AFR = dA, EUR = dE), list(AFR = la, EUR = 2 - la), ops0, c("AFR", "EUR"), "EUR", mac_min = 5)
check("k=2: joint test on 2 df, both backgrounds estimated, heterogeneity on 1 df",
      r$status == "ok" && r$joint$df == 2 && !is.na(r$marginal$AFR$beta) && !is.na(r$marginal$EUR$beta) && r$het$df == 1)
row <- result_row(list(CHR = "22", POS = 1L, ID = "22:1:A:G", REF = "A", ALT = "G"), r, c("AFR", "EUR"), "LAT1")
check("result_row: STRATUM first, then CHR POS ID REF ALT N MAC_STRATUM", identical(names(row)[1:8], c("STRATUM", "CHR", "POS", "ID", "REF", "ALT", "N", "MAC_STRATUM")))

cat(sprintf("\npassed %d, failed %d\n", pass, fail))
quit(status = as.integer(fail > 0))
