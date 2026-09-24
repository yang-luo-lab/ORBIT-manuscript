# ============================================================================
#  ORBIT — simulation figures
#
#  Fig A (Fig 2) : row 1  ORBIT_Rank Z-scatter, 4 scenarios
#          row 2  ORBIT FPR vs rho | K | NA rate, + power
#  Fig B (Sup Fig 3): row 1  RankProd FPR vs rho | K
#          row 2  RankProd FPR vs NA rate, na.rm = TRUE | na.rm = FALSE
#          row 3  -log10 P, ORBIT vs RankProd on discordant data | power
#                 (RankProd solid, ORBIT dashed; fully observed data)
#
#  Colour = signal proportion (no signal = grey); linetype = scenario
#  RankProd: one-class RankProducts, two-sided p = min(2*min(p_up, p_down), 1)
#    na.rm = TRUE  : NA replaced by the gene-wise median (NAreplace), all
#                    features then ranked in every column
#    na.rm = FALSE : NA ignored; per-column ranks among non-NA features; P from
#                    rankprodbounds with nrep = ncol for every feature
# ============================================================================

# ============================================================================
# 0. Setup — locate this script, read/write everything in ../Data
#    Works with Rscript, source(), and RStudio "Source"/"Run".
# ============================================================================
script_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  f <- gsub("~+~", " ", f, fixed = TRUE)
  if (length(f)) return(dirname(normalizePath(f[1])))                          # Rscript
  for (e in rev(sys.frames()))
    if (!is.null(e$ofile)) return(dirname(normalizePath(e$ofile)))             # source()
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path                             # RStudio
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()                                                                        # fallback
})
setwd(script_dir)

# repo layout:  <repo>/Code/<this script>   and   <repo>/Data/  (CSV caches + PDFs)
OUT_DIR <- normalizePath(file.path(script_dir, "..", "Data/Benchmark and simulation"),
                         mustWork = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("Script dir: ", getwd(), "  |  outputs -> ", OUT_DIR)

# dependencies (installed only if missing)
for (pkg in c("ggplot2", "dplyr", "patchwork"))
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("RankProd", quietly = TRUE)) BiocManager::install("RankProd", ask = FALSE)
if (!requireNamespace("ORBIT", quietly = TRUE))
  stop("Package 'ORBIT' not installed. Install it first, e.g. remotes::install_github('<user>/ORBIT').")

library(ORBIT)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RankProd)

# ============================================================================
# 1. Config
# ============================================================================
N_FEAT  <- 10000
N_REP   <- 20
ALPHA   <- 0.05
CACHE   <- TRUE      # reuse CSV if present
VERBOSE <- TRUE      # per-run progress

SCATTER_ORACLE_COR <- TRUE   # scatter panels use oracle ORBIT_cor

RHO_GRID <- seq(0, 0.9, by = 0.1)
K_GRID   <- 2:20
NA_GRID  <- seq(0, 0.9, by = 0.1)

# FPR scenarios: all lack concordant cross-omics signal => every feature is
# null under ORBIT H0, so FPR = mean(P < alpha)
SCENARIOS        <- c("null", "individual", "opposite")
FPR_SIGNAL_PROPS <- c(0.10, 0.20, 0.50)  # individual / opposite scenarios
FPR_MU           <- 2                    # multiplicative scale on signal features
FPR_K_FIXED      <- 2
FPR_RHO_FIXED    <- 0
FPR_SCEN_GRID    <- rbind(
  data.frame(scenario = "null", signal_prop = 0, stringsAsFactors = FALSE),
  expand.grid(scenario = c("individual", "opposite"),
              signal_prop = FPR_SIGNAL_PROPS, stringsAsFactors = FALSE)
)

POWER_PROPS      <- c(0, 0.05, 0.10, 0.20, 0.50)   # 0 = no-signal baseline
POWER_MU         <- 3
POWER_RHO_SIGNAL <- 0.9

# signal proportion -> label / colour (shared by FPR and power)
prop_label  <- function(p) ifelse(p == 0, "No signal", paste0(p * 100, "%"))
PROP_LEVELS <- prop_label(POWER_PROPS)
PROP_PAL    <- setNames(c("grey50", "#4C72B0", "#55A868", "#DD8452", "#C44E52"),
                        PROP_LEVELS)

scen_labels <- c(null = "No signal", individual = "Independent signals",
                 opposite = "Opposite signals")
scen_lty    <- c(null = "solid", individual = "dashed", opposite = "dotted")

RP_MODE_LABEL <- c(`TRUE`  = "na.rm = TRUE (median imputation)",
                   `FALSE` = "na.rm = FALSE (NA ignored)")

# styling
BASE_SIZE   <- 20
BASE_FAMILY <- "Arial"
SIG_COLOR   <- "#F77F00"     # ORBIT-significant points (Fig A scatter)
COL_NULL_PT <- "grey55"      # Fig B panel E
COL_OPP_PT  <- "#C44E52"

theme_panel <- function() {
  theme_minimal(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.line        = element_blank(),
      plot.title       = element_text(hjust = 0.5, family = BASE_FAMILY,
                                      size = BASE_SIZE + 9, face = "plain",
                                      color = "grey15", margin = margin(b = 4)),
      legend.position  = "top",
      legend.title     = element_text(size = BASE_SIZE + 5, face = "bold",
                                      family = BASE_FAMILY),
      legend.text      = element_text(size = BASE_SIZE + 5, family = BASE_FAMILY),
      axis.title       = element_text(size = BASE_SIZE + 9, family = BASE_FAMILY,
                                      color = "grey25"),
      axis.text        = element_text(size = BASE_SIZE, family = BASE_FAMILY),
      aspect.ratio     = 1
    )
}

# ============================================================================
# 2. Progress / caching helpers
# ============================================================================
fmt_hms <- function(sec) {
  if (!is.finite(sec)) return("--:--")
  sec <- max(0, round(sec))
  sprintf("%02d:%02d:%02d", sec %/% 3600, (sec %% 3600) %/% 60, sec %% 60)
}

progress_init <- function(total, label) {
  env <- new.env(parent = emptyenv())
  env$total <- total; env$i <- 0L
  env$label <- label; env$t0 <- Sys.time()
  if (VERBOSE) {
    message("")
    message("=== ", label, " : ", total, " runs ", format(env$t0, "%H:%M:%S"), " ===")
  }
  env
}

progress_step <- function(env, detail = "") {
  env$i <- env$i + 1L
  if (!VERBOSE) return(invisible(NULL))
  elapsed <- as.numeric(difftime(Sys.time(), env$t0, units = "secs"))
  per     <- elapsed / env$i
  eta     <- per * (env$total - env$i)
  frac    <- env$i / env$total
  nbar    <- 24L
  bar     <- paste0(strrep("=", round(frac * nbar)),
                    strrep("-", nbar - round(frac * nbar)))
  cat(sprintf("\r  [%s] %4d/%-4d %5.1f%% | %-34s | %5.2fs/run | el %s | eta %s",
              bar, env$i, env$total, frac * 100, detail, per,
              fmt_hms(elapsed), fmt_hms(eta)),
      file = stderr())
  utils::flush.console()
  if (env$i == env$total) {
    cat("\n", file = stderr())
    message("  done in ", fmt_hms(elapsed), "  (", sprintf("%.2f", per), " s/run)")
  }
  invisible(NULL)
}

# warn if a single run exceeds hb seconds
run_timed <- function(expr, detail, hb = 120) {
  t   <- Sys.time()
  res <- force(expr)
  el  <- as.numeric(difftime(Sys.time(), t, units = "secs"))
  if (VERBOSE && el > hb) message("\n  [slow] ", detail, " took ", fmt_hms(el), " s")
  res
}

run_grid <- function(jobs, label, FUN) {
  pb  <- progress_init(nrow(jobs), label)
  out <- vector("list", nrow(jobs))
  for (i in seq_len(nrow(jobs))) {
    detail <- paste(sprintf("%s=%s", names(jobs), format(jobs[i, ], trim = TRUE)),
                    collapse = " ")
    out[[i]] <- run_timed(FUN(i), detail)
    progress_step(pb, detail)
  }
  do.call(rbind, out)
}

run_or_load <- function(path, expr) {
  path <- file.path(OUT_DIR, path)
  if (CACHE && file.exists(path)) {
    message("[cache] ", path, " exists — skipping simulation")
    return(read.csv(path, stringsAsFactors = FALSE))
  }
  res <- force(expr)
  write.csv(res, path, row.names = FALSE)
  message("[write] ", path, "  (", nrow(res), " rows)")
  res
}

# ============================================================================
# 3. Utilities
# ============================================================================
build_exch <- function(rho, K) { S <- matrix(rho, K, K); diag(S) <- 1; S }

rmvn_exch <- function(nr, K, rho) {
  if (nr == 0) return(matrix(numeric(0), 0, K))
  matrix(rnorm(nr * K), nr, K) %*% chol(build_exch(rho, K))
}

# balanced +/- signs across K omics (exactly K/2 up when K even)
balanced_signs <- function(K) {
  n_up <- if (K %% 2 == 0) K / 2L else sample(c(K %/% 2L, K %/% 2L + 1L), 1)
  sample(rep(c(1, -1), c(n_up, K - n_up)))
}

# (p, sign) per omic -> signed Z matrix (features x omics)
omics_to_Z <- function(omics_list) {
  feats <- omics_list[[1]]$feature
  Z <- vapply(omics_list, function(df) {
    i <- match(feats, df$feature)
    df$sign[i] * qnorm(df$stat[i] / 2, lower.tail = FALSE)
  }, numeric(length(feats)))
  rownames(Z) <- feats
  Z
}

# Z matrix -> list of data.frames(feature, sign, stat = two-sided p)
Z_to_omics <- function(Z) {
  n <- nrow(Z); K <- ncol(Z)
  P        <- pmax(2 * pnorm(-abs(Z)), .Machine$double.xmin)
  sign_mat <- ifelse(Z >= 0, 1, -1)
  feat <- paste0("feature", seq_len(n)); omic <- paste0("omic", seq_len(K))
  setNames(lapply(seq_len(K), function(j)
    data.frame(feature = feat, sign = sign_mat[, j], stat = P[, j])), omic)
}

agg_ci <- function(df, x_var, y_var, extra_groups = character(0)) {
  df %>%
    group_by(across(all_of(c(x_var, extra_groups)))) %>%
    summarise(
      mean  = mean(.data[[y_var]], na.rm = TRUE),
      sd    = sd(.data[[y_var]],   na.rm = TRUE),
      n     = sum(!is.na(.data[[y_var]])),
      se    = sd / sqrt(n),
      ci_lo = pmax(mean - qt(0.975, df = n - 1) * se, 0),
      ci_hi = mean + qt(0.975, df = n - 1) * se,
      .groups = "drop"
    )
}

# extract a ggplot legend as a patchwork element
get_leg <- function(p, pos = c("bottom", "right")) {
  pos <- match.arg(pos)
  g <- ggplotGrob(
    p + theme(legend.position      = pos,
              legend.direction     = if (pos == "bottom") "horizontal" else "vertical",
              legend.justification = if (pos == "bottom") "center" else c(0, 0.5),
              legend.box.margin    = margin(0, 0, 0, 0),
              legend.margin        = margin(0, 0, 0, 0))
  )
  nms <- vapply(g$grobs, function(x) x$name, character(1))
  ix  <- which(nms == paste0("guide-box-", pos))
  if (!length(ix)) ix <- which(nms == "guide-box")
  if (!length(ix)) ix <- grep("^guide-box", nms)
  ix <- ix[vapply(g$grobs[ix], function(x) !inherits(x, "zeroGrob"), logical(1))]
  if (!length(ix)) return(plot_spacer())
  wrap_elements(full = g$grobs[[ix[1]]]) + theme(plot.margin = margin(0, 0, 0, 0))
}

# legend drawn inside the panel (ggplot2 >= 3.5 and older)
legend_inside <- function(p, xy = c(0.02, 0.02), just = c(0, 0)) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    p + theme(legend.position = "inside", legend.position.inside = xy,
              legend.justification.inside = just)
  } else {
    p + theme(legend.position = xy, legend.justification = just)
  }
}

tag_it <- function(p, tag) p + labs(tag = tag) +
  theme(plot.tag = element_text(size = BASE_SIZE + 10, family = BASE_FAMILY),
        plot.tag.position = c(0.01, 0.98))

# ============================================================================
# 4. Simulators
# ============================================================================

# --- K = 2, concordant / discordant signal (Fig A "Opposite Direction") -----
simulate_different_signal <- function(n = 1000,
                                      pi_null = 0.8, pi_concordant = 0.1,
                                      pi_discordant = 0.1, mu_signal = 2,
                                      rho_null = 0, rho_signal = 0, seed = NULL) {
  K <- 2L
  if (abs(sum(c(pi_null, pi_concordant, pi_discordant)) - 1) > 1e-8)
    stop("pi's must sum to 1.")
  if (!is.null(seed)) set.seed(seed)
  n <- as.integer(n)
  L_null   <- chol(build_exch(rho_null,   K))
  L_signal <- chol(build_exch(rho_signal, K))
  group_probs <- c(null = pi_null,
                   up_up = pi_concordant / 2, down_down = pi_concordant / 2,
                   up_down = pi_discordant / 2, down_up = pi_discordant / 2)
  group <- sample(names(group_probs), size = n, replace = TRUE, prob = group_probs)
  is_signal <- group != "null"
  sign_lookup <- rbind(null = c(0, 0), up_up = c(1, 1), down_down = c(-1, -1),
                       up_down = c(1, -1), down_up = c(-1, 1))
  sign_per_omic <- sign_lookup[group, , drop = FALSE]
  Z <- matrix(NA_real_, n, K)
  idx_null <- which(!is_signal); idx_sig <- which(is_signal)
  if (length(idx_null))
    Z[idx_null, ] <- matrix(rnorm(length(idx_null) * K), length(idx_null), K) %*% L_null
  if (length(idx_sig)) {
    Zs <- matrix(rnorm(length(idx_sig) * K), length(idx_sig), K) %*% L_signal
    Z[idx_sig, ] <- sign_per_omic[idx_sig, , drop = FALSE] * (Zs + mu_signal)
  }
  omics_list <- Z_to_omics(Z)
  names(is_signal) <- omics_list[[1]]$feature
  list(omics_list = omics_list, is_signal = is_signal)
}

# --- FPR scenarios ----------------------------------------------------------
#   null       : pure noise
#   individual : each signal feature scaled in one random omic only
#   opposite   : signal in all omics, signs forced balanced
simulate_fpr_scenario <- function(scenario, K = 2, rho = 0, n = N_FEAT,
                                  signal_prop = 0.2, mu = FPR_MU,
                                  na_rate = 0, seed = NULL) {
  scenario <- match.arg(scenario, SCENARIOS)
  if (!is.null(seed)) set.seed(seed)
  
  Z <- rmvn_exch(n, K, rho)
  feat <- paste0("feature", seq_len(n)); omic <- paste0("omic", seq_len(K))
  sig_mat <- matrix(0L, n, K, dimnames = list(feat, omic))
  n_sig <- round(signal_prop * n)
  
  if (scenario == "individual" && n_sig > 0) {
    idx <- seq_len(n_sig)
    j   <- sample.int(K, n_sig, replace = TRUE)
    Z[cbind(idx, j)] <- Z[cbind(idx, j)] * mu
    sig_mat[cbind(idx, j)] <- 1L
  } else if (scenario == "opposite" && n_sig > 0) {
    idx  <- seq_len(n_sig)
    Smat <- t(vapply(seq_len(n_sig), function(i) balanced_signs(K), numeric(K)))
    for (jj in seq_len(K)) Z[idx, jj] <- abs(Z[idx, jj]) * mu * Smat[, jj]
    sig_mat[idx, ] <- 1L
  }
  
  omics_list <- Z_to_omics(Z)
  
  if (na_rate > 0) {
    n_na <- round(na_rate * n)
    if (n_na > 0) {
      na_idx <- sample(seq_len(n), n_na)
      omics_list[[1]]$stat[na_idx] <- NA
      omics_list[[1]]$sign[na_idx] <- NA
    }
  }
  list(omics_list = omics_list, sig_mat = sig_mat)
}

# --- power data (shared by ORBIT and RankProd) ------------------------------
simulate_power_data <- function(K, signal_prop, seed) {
  set.seed(seed)
  ORBIT_simulate_signal(
    n = N_FEAT, K = K,
    pi_null = 1 - signal_prop, pi_up = signal_prop / 2, pi_down = signal_prop / 2,
    mu_signal = POWER_MU, rho_null = 0, rho_signal = POWER_RHO_SIGNAL
  )
}

# ============================================================================
# 5. Analyzers
# ============================================================================

# --- ORBIT ------------------------------------------------------------------
fpr_one <- function(scenario, K, rho, na_rate, seed, signal_prop = 0) {
  sim <- simulate_fpr_scenario(scenario, K = K, rho = rho, signal_prop = signal_prop,
                               na_rate = na_rate, seed = seed)
  direction <- setNames(rep(1, K), names(sim$omics_list))
  
  omics_clean <- lapply(sim$omics_list, function(df)
    df[!is.na(df$stat) & !is.na(df$sign), , drop = FALSE])
  
  # oracle significance flags for ORBIT_cor
  omics_for_cor <- Map(function(df, nm) {
    df$significance <- sim$sig_mat[df$feature, nm]; df
  }, omics_clean, names(omics_clean))
  
  rho_est <- suppressMessages(ORBIT_cor(omics_for_cor, direction, sig_mode = "column"))$rho
  res     <- suppressMessages(ORBIT_Rank(omics_clean, direction, rho = rho_est, seed = seed))
  
  data.frame(scenario = scenario, K = K, rho = rho, na_rate = na_rate,
             signal_prop = signal_prop, rho_est = rho_est,
             fpr = mean(res$P < ALPHA, na.rm = TRUE))
}

# --- RankProd: signed Z -> one-class rank product -> two-sided p ------------
# na_rm = TRUE  : median imputation (NAreplace) then rank all features per column
# na_rm = FALSE : NA ignored; per-column ranks among non-NA; nrep = ncol in the null
rp_two_sided_p <- function(Z, seed, na_rm = TRUE) {
  utils::capture.output(                       # RankProd cat()s an NA notice
    rp <- suppressMessages(suppressWarnings(
      RankProducts(Z, cl = rep(1L, ncol(Z)), logged = TRUE, na.rm = na_rm,
                   gene.names = rownames(Z), plot = FALSE, rand = seed,
                   MinNumOfValidPairs = 1)     # keep feature if >= 1 non-NA omic
    )))
  p2 <- pmin(2 * pmin(rp$pval[, 1], rp$pval[, 2]), 1)
  names(p2) <- rownames(Z)
  p2
}

fpr_one_rp <- function(scenario, K, rho, na_rate, seed, signal_prop = 0, na_rm = TRUE) {
  sim <- simulate_fpr_scenario(scenario, K = K, rho = rho, signal_prop = signal_prop,
                               na_rate = na_rate, seed = seed)
  p2  <- rp_two_sided_p(omics_to_Z(sim$omics_list), seed, na_rm = na_rm)
  data.frame(scenario = scenario, K = K, rho = rho, na_rate = na_rate,
             signal_prop = signal_prop, na_rm = na_rm,
             fpr = mean(p2 < ALPHA, na.rm = TRUE))
}

# --- power (BH-adjusted detection, shared summary) --------------------------
power_summary <- function(K, signal_prop, p, truth) {
  if (anyNA(truth)) stop("Features returned by the method are absent from is_signal.")
  det <- p.adjust(p, method = "BH") < ALPHA
  data.frame(K = K, signal_prop = signal_prop,
             power = sum(det & truth, na.rm = TRUE) / max(sum(truth), 1),
             fdp   = if (sum(det, na.rm = TRUE) > 0)
               sum(det & !truth, na.rm = TRUE) / sum(det, na.rm = TRUE) else 0)
}

power_one <- function(K, signal_prop, seed) {
  sim <- simulate_power_data(K, signal_prop, seed)
  direction <- setNames(rep(1, K), names(sim$omics_list))
  res <- ORBIT_Rank(sim$omics_list, direction, rho = 0, seed = seed)
  power_summary(K, signal_prop, res$P, sim$is_signal[as.character(res$Feature)])
}

power_one_rp <- function(K, signal_prop, seed) {
  sim <- simulate_power_data(K, signal_prop, seed)
  p2  <- rp_two_sided_p(omics_to_Z(sim$omics_list), seed)   # no NA: na.rm irrelevant
  power_summary(K, signal_prop, p2, sim$is_signal[names(p2)])
}

# ============================================================================
# 6. Simulations (cached to CSV)
# ============================================================================

# FPR grid: one parameter varies (`var`), the others fixed; identical seeds for
# ORBIT and RankProd
fpr_jobs <- function(var, grid) {
  g <- setNames(list(grid, seq_len(N_REP)), c(var, "rep"))
  merge(do.call(expand.grid, g), FPR_SCEN_GRID, by = NULL)
}

run_fpr <- function(file, label, var, grid, seed_base, FUN, ...) {
  run_or_load(file, {
    jobs <- fpr_jobs(var, grid)
    run_grid(jobs, label, function(i) {
      args <- list(scenario = jobs$scenario[i], K = FPR_K_FIXED, rho = FPR_RHO_FIXED,
                   na_rate = 0, seed = seed_base + jobs$rep[i],
                   signal_prop = jobs$signal_prop[i], ...)
      args[[var]] <- jobs[[var]][i]
      r <- do.call(FUN, args)
      r$rep <- jobs$rep[i]; r
    })
  })
}

run_power <- function(file, label, FUN) {
  run_or_load(file, {
    jobs <- expand.grid(K = K_GRID, signal_prop = POWER_PROPS, rep = seq_len(N_REP))
    run_grid(jobs, label, function(i) {
      r <- tryCatch(
        FUN(jobs$K[i], jobs$signal_prop[i], seed = 4000L + jobs$rep[i]),
        error = function(e) {
          message("\n  [FAIL] K=", jobs$K[i], " prop=", jobs$signal_prop[i],
                  " rep=", jobs$rep[i], " : ", conditionMessage(e))
          data.frame(K = jobs$K[i], signal_prop = jobs$signal_prop[i],
                     power = NA_real_, fdp = NA_real_)
        })
      r$rep <- jobs$rep[i]; r
    })
  })
}

# --- ORBIT FPR --------------------------------------------------------------
df_rho <- run_fpr("FPR_scenario_rho_v3.csv", "FPR 1/3  varying rho",
                  "rho", RHO_GRID, 1000L, fpr_one)
df_K   <- run_fpr("FPR_scenario_K_v3.csv",   "FPR 2/3  varying K",
                  "K", K_GRID, 2000L, fpr_one)
df_NA  <- run_fpr("FPR_scenario_NA_v3.csv",  "FPR 3/3  varying NA rate",
                  "na_rate", NA_GRID, 3000L, fpr_one)

# --- RankProd FPR (rho / K grids contain no NA, so na.rm is irrelevant there;
#     the NA grid is run once per na.rm mode) --------------------------------
df_rho_rp     <- run_fpr("FPR_RP_rho_v1.csv", "RankProd FPR 1/4  varying rho",
                         "rho", RHO_GRID, 1000L, fpr_one_rp)
df_K_rp       <- run_fpr("FPR_RP_K_v1.csv",   "RankProd FPR 2/4  varying K",
                         "K", K_GRID, 2000L, fpr_one_rp)
df_NA_rp      <- run_fpr("FPR_RP_NA_v1.csv",
                         "RankProd FPR 3/4  varying NA rate (na.rm = TRUE)",
                         "na_rate", NA_GRID, 3000L, fpr_one_rp, na_rm = TRUE)
df_NA_rp_skip <- run_fpr("FPR_RP_NA_narmFALSE_v1.csv",
                         "RankProd FPR 4/4  varying NA rate (na.rm = FALSE)",
                         "na_rate", NA_GRID, 3000L, fpr_one_rp, na_rm = FALSE)

# --- power ------------------------------------------------------------------
df_power    <- run_power("Power_K_signalprop_v2.csv",
                         "Power  varying K x signal prop", power_one)
df_power_rp <- run_power("Power_RP_K_signalprop_v1.csv",
                         "RankProd Power  varying K x signal prop", power_one_rp)

# ============================================================================
# 7. Aggregation + FPR / power panels
# ============================================================================
prep_fpr <- function(df, x_var) {
  df %>%
    agg_ci(x_var, "fpr", extra_groups = c("scenario", "signal_prop")) %>%
    mutate(scenario   = factor(scenario, levels = SCENARIOS),
           prop_label = factor(prop_label(signal_prop), levels = PROP_LEVELS))
}

prep_power <- function(df) {
  df %>%
    agg_ci("K", "power", extra_groups = "signal_prop") %>%
    mutate(signal_label = factor(prop_label(signal_prop), levels = PROP_LEVELS))
}

sum_rho <- prep_fpr(df_rho, "rho")
sum_K   <- prep_fpr(df_K,   "K")
sum_NA  <- prep_fpr(df_NA,  "na_rate")

sum_rho_rp     <- prep_fpr(df_rho_rp,     "rho")
sum_K_rp       <- prep_fpr(df_K_rp,       "K")
sum_NA_rp      <- prep_fpr(df_NA_rp,      "na_rate")
sum_NA_rp_skip <- prep_fpr(df_NA_rp_skip, "na_rate")

sum_power    <- prep_power(df_power)
sum_power_rp <- prep_power(df_power_rp)

# y-axis: ORBIT and RankProd scaled separately (RankProd FPR can be much higher)
fpr_axis <- function(...) {
  y_max <- max(c(unlist(lapply(list(...), `[[`, "ci_hi")), 0.10), na.rm = TRUE)
  y_max <- ceiling(y_max * 40) / 40
  brks  <- pretty(c(0, y_max), n = 5)
  list(max = y_max, breaks = brks[brks <= y_max])
}
ax_orbit <- fpr_axis(sum_rho, sum_K, sum_NA)
ax_rp    <- fpr_axis(sum_rho_rp, sum_K_rp, sum_NA_rp, sum_NA_rp_skip)

make_fpr_plot <- function(d, x_var, x_label, x_breaks, show_y_title = TRUE,
                          ax = ax_orbit) {
  fpr_prop_breaks <- c("No signal", prop_label(FPR_SIGNAL_PROPS))
  ggplot(d, aes(x = .data[[x_var]], y = mean,
                color = prop_label, fill = prop_label, linetype = scenario,
                group = interaction(scenario, prop_label))) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.10, color = NA,
                show.legend = FALSE) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = ALPHA, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    scale_color_manual(values = PROP_PAL, breaks = fpr_prop_breaks,
                       name = "Signal Proportion") +
    scale_fill_manual(values = PROP_PAL, breaks = fpr_prop_breaks,
                      name = "Signal Proportion") +
    scale_linetype_manual(values = scen_lty, labels = scen_labels,
                          breaks = SCENARIOS, name = NULL, drop = FALSE) +
    guides(color    = guide_legend(order = 1,
                                   override.aes = list(linewidth = 2, linetype = "solid")),
           fill     = "none",
           linetype = guide_legend(order = 2,
                                   override.aes = list(color = "black", linewidth = 1))) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = ax$breaks, limits = c(0, ax$max),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = x_label, y = if (show_y_title) "False Positive Rate" else NULL) +
    theme_panel()
}

make_power_plot <- function(sum_df) {
  d   <- filter(sum_df, signal_prop > 0)                
  lev <- PROP_LEVELS[PROP_LEVELS != "No signal"]  
  ggplot(d, aes(x = K, y = mean, color = signal_label,
                fill = signal_label, group = signal_label)) +
    geom_ribbon(aes(ymin = pmax(ci_lo, 0), ymax = pmin(ci_hi, 1)),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 1.5) +
    scale_x_continuous(breaks = seq(2, 20, by = 2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_color_manual(values = PROP_PAL, limits = lev, name = "Signal proportion") +
    scale_fill_manual(values = PROP_PAL, limits = lev, name = "Signal proportion") +
    guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) +
    labs(x = "Number of Omics (K)", y = "Power") +
    theme_panel() +
    theme(legend.title = element_text(family = BASE_FAMILY, face = "bold", size = 20))
}

# RankProd solid, ORBIT dashed
make_power_plot_both <- function(sum_orbit, sum_rp) {
  d <- bind_rows(mutate(sum_orbit, method = "ORBIT"),
                 mutate(sum_rp,    method = "RankProd")) %>%
    mutate(method = factor(method, levels = c("RankProd", "ORBIT")))
  
  ggplot(d, aes(x = K, y = mean, color = signal_label, fill = signal_label,
                linetype = method, group = interaction(method, signal_label))) +
    geom_ribbon(aes(ymin = pmax(ci_lo, 0), ymax = pmin(ci_hi, 1)),
                alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.6) +
    scale_x_continuous(breaks = seq(2, 20, by = 2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_color_manual(values = PROP_PAL, limits = PROP_LEVELS, name = "Signal proportion") +
    scale_fill_manual(values = PROP_PAL, limits = PROP_LEVELS, name = "Signal proportion") +
    scale_linetype_manual(values = c(RankProd = "solid", ORBIT = "dashed"), name = "Method") +
    guides(color    = guide_legend(order = 1, ncol = 1,
                                   override.aes = list(linewidth = 2, linetype = "solid")),
           fill     = "none",
           linetype = guide_legend(order = 2, ncol = 1,
                                   override.aes = list(color = "black", linewidth = 1.3))) +
    labs(x = "Number of Omics (K)", y = "Power") +
    theme_panel() +
    theme(legend.position  = "right",
          legend.direction = "vertical",
          legend.box       = "vertical",
          legend.title     = element_text(family = BASE_FAMILY, face = "bold", size = 22),
          legend.text      = element_text(family = BASE_FAMILY, size = 20),
          legend.key.width = unit(1.6, "cm"),
          legend.spacing.y = unit(0.6, "cm"),
          axis.title       = element_text(size = BASE_SIZE + 12, family = BASE_FAMILY,
                                          color = "grey25"),
          axis.text        = element_text(size = BASE_SIZE + 3, family = BASE_FAMILY))
}

# ORBIT panels (Fig A)
pr <- make_fpr_plot(sum_rho, "rho",     "Inter-Omics Correlation", RHO_GRID, TRUE)
pk <- make_fpr_plot(sum_K,   "K",       "Number of Omics (K)", seq(2, 20, by = 2), FALSE)
pn <- make_fpr_plot(sum_NA,  "na_rate", "NA Rate in Omic 1", NA_GRID, FALSE)
pp <- make_power_plot(sum_power)

# RankProd panels (Fig B; left panels carry the y title)
pr_rp      <- make_fpr_plot(sum_rho_rp,     "rho",     "Inter-Omics Correlation",
                            RHO_GRID, TRUE,  ax_rp)
pk_rp      <- make_fpr_plot(sum_K_rp,       "K",       "Number of Omics (K)",
                            seq(2, 20, by = 2), FALSE, ax_rp)
pn_rp      <- make_fpr_plot(sum_NA_rp,      "na_rate", "NA Rate in Omic 1",
                            NA_GRID, TRUE,  ax_rp)
pn_rp_skip <- make_fpr_plot(sum_NA_rp_skip, "na_rate", "NA Rate in Omic 1",
                            NA_GRID, FALSE, ax_rp)
pp_B       <- make_power_plot_both(sum_power, sum_power_rp)

# ============================================================================
# 8. Scatter panels
# ============================================================================

# --- Fig A row 1: Z-scatter, ORBIT-significant features highlighted ---------
build_scatter <- function(sim, panel_title = NULL) {
  direction <- setNames(rep(1, 2), names(sim$omics_list))
  omics_for_cor <- Map(function(o, nm) {
    o$significance <-
      if (!SCATTER_ORACLE_COR)            0L
    else if (!is.null(sim$is_signal_omic)) sim$is_signal_omic[o$feature, nm]
    else if (!is.null(sim$is_signal))      sim$is_signal[o$feature]
    else                                   0L
    o
  }, sim$omics_list, names(sim$omics_list))
  
  rho_est  <- suppressMessages(ORBIT_cor(omics_for_cor, direction, sig_mode = "column"))$rho
  res_Rank <- ORBIT_Rank(sim$omics_list, direction, rho = rho_est, seed = 42)
  res_Rank$fdr <- p.adjust(res_Rank$P, method = "BH")
  
  o1 <- sim$omics_list[[1]]; o2 <- sim$omics_list[[2]]
  common <- intersect(o1$feature, o2$feature)
  i1 <- match(common, o1$feature); i2 <- match(common, o2$feature)
  iR <- match(common, res_Rank$Feature)
  
  df <- data.frame(x   = o1$sign[i1] * qnorm(o1$stat[i1] / 2, lower.tail = FALSE),
                   y   = o2$sign[i2] * qnorm(o2$stat[i2] / 2, lower.tail = FALSE),
                   sig = !is.na(res_Rank$fdr[iR]) & res_Rank$fdr[iR] < 0.05)
  
  p <- ggplot(df, aes(x, y)) +
    geom_point(data = subset(df, !sig), alpha = 0.4, size = 1.0,
               color = "grey75", shape = 16) +
    geom_point(data = subset(df, sig), color = SIG_COLOR, alpha = 0.85,
               size = 2, shape = 16) +
    geom_hline(yintercept = 0, color = "grey25", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey25", linewidth = 0.4) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
    labs(x = "Omic 1: Z-statistic", y = "Omic 2: Z-statistic") +
    theme_panel() + theme(legend.position = "none")
  if (!is.null(panel_title)) p <- p + ggtitle(panel_title)
  p
}

# --- Fig B panel E: -log10 P_ORBIT vs -log10 P_RankProd on discordant data --
build_disc_p_scatter <- function(sim, seed = 42) {
  feats     <- sim$omics_list[[1]]$feature
  direction <- setNames(rep(1, length(sim$omics_list)), names(sim$omics_list))
  
  # ORBIT (oracle rho, as in the Fig A scatter panels)
  omics_for_cor <- lapply(sim$omics_list, function(o) {
    o$significance <- as.integer(sim$is_signal[o$feature]); o })
  rho_est <- suppressMessages(ORBIT_cor(omics_for_cor, direction, sig_mode = "column"))$rho
  res     <- suppressMessages(ORBIT_Rank(sim$omics_list, direction, rho = rho_est, seed = seed))
  p_orbit <- setNames(res$P, as.character(res$Feature))[feats]
  
  # RankProd (same method as the benchmark)
  p_rp <- rp_two_sided_p(omics_to_Z(sim$omics_list), seed = seed)[feats]
  
  df <- data.frame(nl_orbit = -log10(pmax(p_orbit, 1e-300)),
                   nl_rp    = -log10(pmax(p_rp,    1e-300)),
                   opposite = unname(sim$is_signal[feats]))
  df$type <- factor(ifelse(df$opposite, "Opposite-direction signal", "No signal"),
                    levels = c("No signal", "Opposite-direction signal"))
  
  fdr_o <- p.adjust(p_orbit, "BH"); fdr_r <- p.adjust(p_rp, "BH")
  message(sprintf("  [disc scatter] rho_est = %.3f | opposite features BH<0.05: ORBIT %d, RankProd %d (of %d)",
                  rho_est, sum(fdr_o[df$opposite] < ALPHA),
                  sum(fdr_r[df$opposite] < ALPHA), sum(df$opposite)))
  
  lim <- c(0, max(df$nl_orbit, df$nl_rp, -log10(ALPHA)) * 1.04)
  
  ggplot(df, aes(nl_orbit, nl_rp)) +
    geom_point(data = subset(df, !opposite), aes(color = type),
               alpha = 0.4, size = 1.0, shape = 16) +
    geom_point(data = subset(df, opposite), aes(color = type),
               alpha = 0.85, size = 1.8, shape = 16) +
    geom_abline(slope = 1, intercept = 0, color = "grey40",
                linetype = "dotted", linewidth = 0.5) +
    geom_hline(yintercept = -log10(ALPHA), color = "grey40",
               linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = -log10(ALPHA), color = "grey40",
               linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = c(`No signal` = COL_NULL_PT,
                                  `Opposite-direction signal` = COL_OPP_PT),
                       name = NULL) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = lim, ylim = lim) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 5, alpha = 1))) +
    labs(x = expression(ORBIT: -log[10] * " P"),
         y = expression(RankProd: -log[10] * " P")) +
    theme_panel() +
    theme(legend.title = element_blank(),
          legend.key   = element_blank(),
          axis.title   = element_text(size = BASE_SIZE + 12, family = BASE_FAMILY,
                                      color = "grey25"),
          axis.text    = element_text(size = BASE_SIZE + 3, family = BASE_FAMILY))
}

# --- data for the scatter panels ---------------------------------------------
message("\n=== Building scatter panels ===")
sim_null <- ORBIT_simulate_null(n = 10000, K = 2, rho = 0, seed = 42)
sim_two  <- local({
  s <- simulate_fpr_scenario("individual", K = 2, rho = 0, n = 10000,
                             signal_prop = 0.2, mu = 1.5, seed = 1)
  list(omics_list = s$omics_list, is_signal_omic = s$sig_mat)
})
sim_opp  <- simulate_different_signal(n = 10000, pi_null = 0.80,
                                      pi_concordant = 0, pi_discordant = 0.2,
                                      mu_signal = 1, rho_null = 0,
                                      rho_signal = 0.85, seed = 1)
sim_sig  <- ORBIT_simulate_signal(n = 10000, K = 2, pi_null = 0.8,
                                  pi_up = 0.1, pi_down = 0.1,
                                  mu_signal = 1, rho_null = 0, rho_signal = 0.85)

scatter_specs <- list(list(sim_null, "No Signal"),
                      list(sim_two,  "Two Independent Signals"),
                      list(sim_opp,  "Opposite Direction Signal"),
                      list(sim_sig,  "Same Direction Signal"))
pb_s <- progress_init(length(scatter_specs), "Scatter panels")
scatters <- lapply(scatter_specs, function(sp) {
  p <- build_scatter(sp[[1]], sp[[2]])
  progress_step(pb_s, sp[[2]])
  p
})

# Fig B panel E: same data as Fig A "Opposite Direction Signal"
p_disc <- build_disc_p_scatter(sim_opp)

# ============================================================================
# 9. Figure A : scatter + ORBIT
# ============================================================================
s1 <- scatters[[1]] + theme(axis.title.x = element_blank())
s2 <- scatters[[2]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
s3 <- scatters[[3]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
s4 <- scatters[[4]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

xtitle_scatter <- ggplot() +
  labs(title = "Omic 1: Z-statistic") +
  theme_void(base_family = BASE_FAMILY) +
  theme(plot.title = element_text(hjust = 0.5, color = "grey25", family = BASE_FAMILY,
                                  size = BASE_SIZE + 7, margin = margin(0, 0, 0, 0)),
        plot.margin = margin(0, 0, 0, 0))

scatter_legend <- get_leg(
  ggplot(data.frame(x = 1, y = 1, g = factor("sig")), aes(x, y, color = g)) +
    geom_point(size = 2.5, shape = 16, alpha = 0.85) +
    scale_color_manual(values = c(sig = SIG_COLOR), labels = c(sig = "ORBIT\nSignificant")) +
    labs(color = NULL) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    theme(legend.text = element_text(size = 23, family = BASE_FAMILY, face = "bold"),
          legend.key = element_blank(),
          legend.key.size = unit(0.7, "cm")),
  pos = "right")

power_legend <- get_leg(pp, pos = "right")

# bottom legend: colour (signal proportion) + linetype (scenario), horizontal
make_bottom_legend <- function(p) get_leg(
  p +
    guides(color    = guide_legend(nrow = 1, order = 1,
                                   override.aes = list(linewidth = 2, linetype = "solid")),
           fill     = "none",
           linetype = guide_legend(nrow = 1, order = 2,
                                   override.aes = list(color = "black", linewidth = 1))) +
    theme(legend.box = "horizontal", legend.spacing.x = unit(0.8, "cm")),
  pos = "bottom")
bottom_legend_A <- make_bottom_legend(pr)

noleg <- theme(legend.position = "none")

design_A <- c(
  area(1, 1), area(1, 2), area(1, 3), area(1, 4), area(1, 5),   # scatter x4 + legend
  area(2, 1, 2, 4),                                             # shared x title
  area(3, 1), area(3, 2), area(3, 3), area(3, 4),               # ORBIT row
  area(3, 5),                                                   # power legend
  area(4, 1, 4, 4)                                              # bottom legend
)

fig_A <- wrap_plots(
  s1, s2, s3, s4, scatter_legend,
  xtitle_scatter,
  pr + noleg + ggtitle("ORBIT"), pk + noleg, pn + noleg, pp + noleg,
  power_legend,
  bottom_legend_A,
  design  = design_A,
  widths  = c(1, 1, 1, 1, 0.5),
  heights = c(1, 0.04, 1, 0.12)
) & theme(plot.margin = margin(5, 5, 1.5, 1.5))

# ============================================================================
# 10. Figure B : 3 x 2 panel grid + legend column
#     row 1 : RankProd FPR vs rho | vs K
#     row 2 : RankProd FPR vs NA rate, na.rm = TRUE | na.rm = FALSE
#     row 3 : discordant-scenario P scatter (ORBIT vs RankProd) | power
#     Legends: FPR legend spans rows 1-2, power legend on row 3; panel E's
#     legend is drawn inside the panel (bottom-left).
# ============================================================================

# vertical FPR legend (colour block above linetype block)
make_right_legend <- function(p) get_leg(
  p +
    guides(color    = guide_legend(ncol = 1, order = 1,
                                   override.aes = list(linewidth = 2, linetype = "solid")),
           fill     = "none",
           linetype = guide_legend(ncol = 1, order = 2,
                                   override.aes = list(color = "black", linewidth = 1))) +
    theme(legend.box       = "vertical",
          legend.box.just  = "left",
          legend.spacing.y = unit(0.8, "cm"),
          legend.key.width = unit(1.4, "cm"),
          legend.title     = element_text(family = BASE_FAMILY, face = "bold", size = 22),
          legend.text      = element_text(family = BASE_FAMILY, size = 20)),
  pos = "right")

legend_fpr_B <- make_right_legend(pr_rp)
legend_power <- get_leg(pp_B, pos = "right")

na_title <- theme(plot.title = element_text(size = BASE_SIZE + 2))

pA <- pr_rp      + noleg + ggtitle("RankProd")
pB <- pk_rp      + noleg
pC <- pn_rp      + noleg + ggtitle(RP_MODE_LABEL["TRUE"])  + na_title
pD <- pn_rp_skip + noleg + ggtitle(RP_MODE_LABEL["FALSE"]) + na_title
pE <- legend_inside(p_disc) +
  theme(legend.direction  = "vertical",
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.margin     = margin(4, 8, 4, 6),
        legend.text       = element_text(family = BASE_FAMILY, size = 18),
        legend.key.size   = unit(0.6, "cm"))
pF <- pp_B + noleg

design_B <- c(
  area(1, 1), area(1, 2),                  # rho  | K
  area(2, 1), area(2, 2),                  # NA T | NA F
  area(1, 3, 2, 3),                        # FPR legend spanning rows 1-2
  area(3, 1), area(3, 2),                  # disc scatter | power
  area(3, 3)                               # power legend
)

fig_B <- wrap_plots(
  tag_it(pA, "A"), tag_it(pB, "B"),
  tag_it(pC, "C"), tag_it(pD, "D"),
  legend_fpr_B,
  tag_it(pE, "E"), tag_it(pF, "F"),
  legend_power,
  design  = design_B,
  widths  = c(1, 1, 0.6),
  heights = c(1, 1, 1)
) & theme(plot.margin = margin(5, 5, 1.5, 1.5))

# ============================================================================
# 11. Render + save
# ============================================================================
message("\n=== Rendering figures ===")
print(fig_A)
print(fig_B)

# portable PDF writer: quartz on macOS, cairo_pdf elsewhere
save_pdf <- function(p, file, width, height) {
  file <- file.path(OUT_DIR, file)
  if (capabilities("aqua")) {
    quartz(type = "pdf", file = file, width = width, height = height, family = BASE_FAMILY)
    print(p); dev.off()
  } else {
    ggsave(file, p, width = width, height = height, device = cairo_pdf)
  }
  message("[pdf] ", file)
}

save_pdf(fig_A, "Figure_2.pdf", width = 27, height = 12)
save_pdf(fig_B, "Sup Fig 3.pdf",      width = 16, height = 26)
