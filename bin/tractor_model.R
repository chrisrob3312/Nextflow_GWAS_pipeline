# ============================================================================
# tractor_model.R: the Tractor-GENESIS score test, apart from the I/O around it
# ============================================================================
# Sourced by tractor_genesis_adapter.R (from the same directory); also by tests/test_tractor_model.R,
# which runs it without GENESIS. Needs Matrix (and survival for the Cox null, in the adapter).
#   genesis_ops(nullmod, n)          P y and D'PD from a GENESIS null model (binary / quantitative)
#   cox_ops(time, event, eta, Xc)    the same from a Cox null: martingale residuals, Breslow information
#   mac_minor(alt, n_hap)            minor allele count
#   test_variant_tractor(...)        the conditional score test of one variant
#   result_row(...)                  one output row
#   slice_chunk / test_chunk         a block of variants: its matrices, and its rows + status
# ============================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a

# --- helpers for GENESIS null-model geometry ---------------------------------
# GENESIS stores cholSigmaInv = L with SigmaInv = L L'; CX = L'X; CXCXI = CX (CX'CX)^-1
.as_L <- function(nm, n) {
    L <- nm$cholSigmaInv
    if (is.null(L)) stop("GENESIS null model has no cholSigmaInv; update GENESIS (>= 2.16)")
    if (is.matrix(L) || inherits(L, "Matrix")) return(L)
    if (length(L) == n) return(Matrix::Diagonal(n, x = as.numeric(L)))
    Matrix::Diagonal(n, x = rep(as.numeric(L)[1], n))
}
.crossL <- function(L, M) as.matrix(Matrix::crossprod(L, M))   # L'M
.multL  <- function(L, v) as.numeric(L %*% v)                    # L v

genesis_ops <- function(nm, n) {
    L <- .as_L(nm, n)
    CX <- as.matrix(nm$CX); CXCXI <- as.matrix(nm$CXCXI)
    PY <- nm$fit$resid.PY %||% nm$resid
    if (is.null(PY)) {
        Y  <- nm$fit$workingY %||% nm$workingY
        Yt <- .crossL(L, Y)
        PY <- .multL(L, Yt - CXCXI %*% crossprod(CX, Yt))
    }
    list(
        type = "genesis",
        PY = as.numeric(PY),
        # projected cross-product  D'PD  for a design matrix D (n x p)
        proj_crossprod = function(D) {
            Dt <- .crossL(L, D)
            crossprod(Dt) - crossprod(Dt, CXCXI) %*% crossprod(CX, Dt)
        }
    )
}

cox_ops <- function(time, event, eta, Xc) {
    w <- exp(eta)
    n <- length(time)
    # Breslow baseline hazard and martingale residuals under the null
    ord_asc <- order(time, -event)
    t_a <- time[ord_asc]; w_a <- w[ord_asc]; e_a <- event[ord_asc]
    S0_a <- rev(cumsum(rev(w_a)))                      # sum_{j: t_j >= t_i} w_j
    S0_a <- S0_a[match(t_a, t_a)]                      # ties share one risk set
    dL   <- ifelse(e_a == 1, 1 / S0_a, 0)
    Lam  <- ave(cumsum(dL), t_a, FUN = max)            # Lambda0(t_i) incl. all events at t_i
    mart <- numeric(n); mart[ord_asc] <- e_a - Lam * w_a

    # Risk-set machinery in DESCENDING time: risk set of t_k = positions 1..last(tie group)
    ord_d <- order(-time, event)
    t_d <- time[ord_d]; w_d <- w[ord_d]; e_d <- event[ord_d]
    last_idx <- ave(seq_along(t_d), t_d, FUN = max)
    ev_pos <- which(e_d == 1)
    Xc_d <- if (!is.null(Xc)) Xc[ord_d, , drop = FALSE] else NULL

    # Observed information (Breslow) for full design A = [D | Xc], then project on Xc
    info <- function(A_d) {
        p <- ncol(A_d)
        S0 <- cumsum(w_d)[last_idx]
        S1 <- apply(w_d * A_d, 2, cumsum)[last_idx, , drop = FALSE]
        I  <- matrix(0, p, p)
        for (a in 1:p) for (b in a:p) {
            S2 <- cumsum(w_d * A_d[, a] * A_d[, b])[last_idx]
            v  <- sum((S2 / S0 - (S1[, a] / S0) * (S1[, b] / S0))[ev_pos])
            I[a, b] <- v; I[b, a] <- v
        }
        I
    }
    q <- if (is.null(Xc_d)) 0 else ncol(Xc_d)
    list(
        type = "cox",
        PY = mart,
        proj_crossprod = function(D) {
            A <- if (q > 0) cbind(D[ord_d, , drop = FALSE], Xc_d) else D[ord_d, , drop = FALSE]
            I <- info(A)
            p <- ncol(D)
            if (q == 0) return(I)
            iD <- 1:p; iX <- (p + 1):(p + q)
            I[iD, iD] - I[iD, iX] %*% solve(I[iX, iX] + diag(1e-8, q)) %*% I[iX, iD]
        }
    )
}


# ---------------------------------------------------------------------------
# Minor allele counts (CP-17)
# ---------------------------------------------------------------------------
# On an ancestry background a, the alleles are the haplotypes of that ancestry: sum of the
# hapcount (0, 1 or 2 per person). ALT_a = sum of the dosage on those haplotypes. The minor allele
# count is min(ALT_a, haplotypes_a - ALT_a): an allele at ALT frequency 0.94 on a background has
# a small MINOR count even though its ALT count is large. Stratum-wide: min(ALT, 2n - ALT).
mac_minor <- function(alt, n_hap) pmax(0, pmin(alt, n_hap - alt))

empty_result <- function(ancestries, n_dose, mac, mac_stratum, status) list(
    joint = list(stat = NA, df = n_dose, p = NA),
    marginal = setNames(lapply(ancestries, function(a) list(ancestry = a, beta = NA, se = NA, z = NA, p = NA)), ancestries),
    het = list(Q = NA, df = n_dose - 1, p = NA, I2 = NA),
    n_eff = NA, mac = mac, mac_stratum = mac_stratum,
    status = status          # "mac": below the MAC gate; "error": the test failed
)

# Conditional score test (identical for binary / quantitative / survival):
#   D = [LA_1..LA_{k-1} | Dose_1..Dose_k]
#   U = D' P y ,  M = D' P D      (P, Py supplied by null_ops)
#   Dose block conditioned on the LA block (Schur complement):
#     U_c = U_G - M_GZ M_ZZ^-1 U_Z ,   V_c = M_GG - M_GZ M_ZZ^-1 M_ZG
#   Joint (k-df):    U_c' V_c^-1 U_c  ~ chi2_k       (H0: all dose effects 0)
#   Per-ancestry:    beta = V_c^-1 U_c (one-step estimate), se = sqrt(diag V_c^-1)
#   Heterogeneity:   contrasts C beta (beta_j - beta_1),  Q = (Cb)'(C V_c^-1 C')^-1 (Cb) ~ chi2_{k-1}
#                    - uses the full covariance, so correlated dose terms are handled
# MAC gate (CP-17): the variant is tested only if its stratum-wide minor allele count reaches
# mac_min; within it, an ancestry's dose term only if that background's minor allele count does
# (and the dose varies).
test_variant_tractor <- function(
    la_vec,      # named list: LA counts per non-reference ancestry
    dose_vec,    # named list: dosages per ancestry (all k)
    hap_vec,     # named list: hapcounts per ancestry (all k), for the per-background MAC
    null_ops,    # list(PY, proj_crossprod) from the null model
    ancestries,  # all ancestry names
    ref_anc,     # reference ancestry
    mac_min = 10
) {
    n <- length(dose_vec[[1]])
    non_ref <- setdiff(ancestries, ref_anc)
    n_dose <- length(ancestries)

    # Build design: LA terms then Dose terms; mean-impute sporadic missingness
    D <- matrix(0, nrow = n, ncol = length(non_ref) + n_dose)
    colnames(D) <- c(paste0("LA_", non_ref), paste0("Dose_", ancestries))
    for (anc in non_ref)    D[, paste0("LA_", anc)]   <- as.numeric(la_vec[[anc]])
    for (anc in ancestries) D[, paste0("Dose_", anc)] <- as.numeric(dose_vec[[anc]])
    n_missing <- sum(!complete.cases(D))
    if (n_missing > 0) {
        for (j in seq_len(ncol(D))) { miss <- is.na(D[, j]); if (any(miss)) D[miss, j] <- mean(D[, j], na.rm = TRUE) }
    }

    # Minor allele counts: per background, and over the whole stratum
    alt <- vapply(ancestries, function(a) sum(D[, paste0("Dose_", a)]), 0)
    n_hap <- vapply(ancestries, function(a) sum(as.numeric(hap_vec[[a]]), na.rm = TRUE), 0)
    mac <- setNames(mac_minor(alt, n_hap), ancestries)
    mac_stratum <- mac_minor(sum(alt), 2 * n)
    if (mac_stratum < mac_min) return(empty_result(ancestries, n_dose, mac, mac_stratum, "mac"))
    testable <- ancestries[mac >= mac_min & apply(D[, paste0("Dose_", ancestries), drop = FALSE], 2, var) > 1e-10]
    if (length(testable) == 0) return(empty_result(ancestries, n_dose, mac, mac_stratum, "mac"))

    # Drop LA columns with no variation (e.g. no AFR tracts in this chunk)
    la_cols <- paste0("LA_", non_ref)
    la_cols <- la_cols[apply(D[, la_cols, drop = FALSE], 2, var) > 1e-10]
    g_cols  <- paste0("Dose_", testable)
    Dk <- D[, c(la_cols, g_cols), drop = FALSE]

    res <- tryCatch({
        U <- as.numeric(crossprod(Dk, null_ops$PY))
        M <- as.matrix(null_ops$proj_crossprod(Dk))
        iZ <- seq_along(la_cols); iG <- length(la_cols) + seq_along(g_cols)

        if (length(iZ) > 0) {
            MZZi <- solve(M[iZ, iZ, drop = FALSE] + diag(1e-8, length(iZ)))
            U_c <- U[iG] - M[iG, iZ, drop = FALSE] %*% MZZi %*% U[iZ]
            V_c <- M[iG, iG, drop = FALSE] - M[iG, iZ, drop = FALSE] %*% MZZi %*% M[iZ, iG, drop = FALSE]
        } else { U_c <- U[iG]; V_c <- M[iG, iG, drop = FALSE] }
        V_c <- (V_c + t(V_c)) / 2
        V_inv <- solve(V_c + diag(1e-8, nrow(V_c)))

        # Joint test
        stat_joint <- as.numeric(t(U_c) %*% V_inv %*% U_c)
        df_joint <- length(iG)
        p_joint <- pchisq(stat_joint, df_joint, lower.tail = FALSE)

        # Per-ancestry one-step estimates
        beta <- as.numeric(V_inv %*% U_c); se <- sqrt(diag(V_inv))
        z <- beta / se; p <- 2 * pnorm(-abs(z))
        marginal <- setNames(lapply(ancestries, function(a) list(ancestry = a, beta = NA, se = NA, z = NA, p = NA)), ancestries)
        for (i in seq_along(testable)) {
            marginal[[testable[i]]] <- list(ancestry = testable[i], beta = beta[i], se = se[i], z = z[i], p = p[i])
        }

        # Heterogeneity across ancestry backgrounds (covariance-aware)
        if (length(testable) >= 2) {
            k <- length(testable)
            C <- cbind(-1, diag(k - 1))               # beta_j - beta_1
            Cb <- C %*% beta
            Q <- as.numeric(t(Cb) %*% solve(C %*% V_inv %*% t(C) + diag(1e-10, k - 1)) %*% Cb)
            df_het <- k - 1
            p_het <- pchisq(Q, df_het, lower.tail = FALSE)
            I2 <- max(0, (Q - df_het) / Q * 100)
        } else { Q <- NA; df_het <- n_dose - 1; p_het <- NA; I2 <- NA }

        list(joint = list(stat = stat_joint, df = df_joint, p = p_joint),
             marginal = marginal,
             het = list(Q = Q, df = df_het, p = p_het, I2 = I2),
             n_eff = n - n_missing, mac = mac, mac_stratum = mac_stratum, status = "ok")
    }, error = function(e) empty_result(ancestries, n_dose, mac, mac_stratum, "error"))

    res
}

# result_row(var_info, res, ancestries, stratum): one output row (a list), STRATUM first
result_row <- function(var_info, res, ancestries, stratum = NULL) {
    vi <- as.list(var_info)
    row <- list()
    if (!is.null(stratum)) row$STRATUM <- stratum
    row$CHR <- as.character(vi$CHR)
    row$POS <- as.integer(vi$POS)
    row$ID  <- if (!is.null(vi$ID)) as.character(vi$ID) else paste0(vi$CHR, ":", vi$POS)
    row$REF <- if (!is.null(vi$REF)) as.character(vi$REF) else NA_character_
    row$ALT <- if (!is.null(vi$ALT)) as.character(vi$ALT) else NA_character_
    row$N <- as.integer(res$n_eff)
    row$MAC_STRATUM <- as.numeric(res$mac_stratum)
    row$STAT_JOINT <- as.numeric(res$joint$stat)
    row$DF_JOINT <- as.integer(res$joint$df)
    row$P_JOINT <- as.numeric(res$joint$p)
    for (anc in ancestries) {
        row[[paste0("MAC_", anc)]]  <- as.numeric(res$mac[[anc]])
        row[[paste0("BETA_", anc)]] <- as.numeric(res$marginal[[anc]]$beta)
        row[[paste0("SE_", anc)]]   <- as.numeric(res$marginal[[anc]]$se)
        row[[paste0("Z_", anc)]]    <- as.numeric(res$marginal[[anc]]$z)
        row[[paste0("P_", anc)]]    <- as.numeric(res$marginal[[anc]]$p)
    }
    row$Q_HET <- as.numeric(res$het$Q)
    row$DF_HET <- as.integer(res$het$df)
    row$P_HET <- as.numeric(res$het$p)
    row$I2 <- as.numeric(res$het$I2)
    row
}

# slice_chunk(td, idx): the rows idx of the Tractor matrices (a chunk the worker tests)
slice_chunk <- function(td, idx) list(
    variant_info = td$variant_info[idx],
    hapcount = lapply(td$hapcount, function(x) x[idx, , drop = FALSE]),
    dosage   = lapply(td$dosage,   function(x) x[idx, , drop = FALSE]))

# test_chunk(ch, ...): every variant of a chunk -> list(dt = its rows, status = per-variant status)
test_chunk <- function(ch, null_ops, ancestries, ref_anc, mac_min, stratum = NULL) {
    non_ref <- setdiff(ancestries, ref_anc)
    nv <- nrow(ch$variant_info)
    rows <- vector("list", nv); status <- character(nv)
    for (v in seq_len(nv)) {
        dose_vec <- lapply(setNames(ancestries, ancestries), function(a) ch$dosage[[a]][v, ])
        hap_vec  <- lapply(setNames(ancestries, ancestries), function(a) ch$hapcount[[a]][v, ])
        res <- test_variant_tractor(hap_vec[non_ref], dose_vec, hap_vec, null_ops, ancestries, ref_anc, mac_min)
        status[v] <- res$status
        rows[[v]] <- result_row(ch$variant_info[v], res, ancestries, stratum)
    }
    list(dt = data.table::rbindlist(rows), status = status)
}
