# ============================================================================
# Supplementary Figure 8 -- FPR under the null + power comparison
#
#   FPR  : ORBIT_P vs ActivePathways p-value merging (Fisher, Brown,
#          Stouffer_directional, DPM) on ORBIT_simulate_null() data.
#          Panels: FPR vs inter-omics rho (K = 2) | vs K (rho = 0) | vs NA rate in omic 1
#
#   Power: DPM vs ORBIT_P on ORBIT_simulate_signal() data.
#          Grid = K x signal proportion, mu = 3, rho_signal = 0.9, rho_null = 0,
#          BH-adjusted at ALPHA (same design as the main-text power simulation).
#          Color = signal proportion, linetype = method (DPM solid, ORBIT dashed).
#
#   Fisher: implemented directly, NA-aware (sum over observed omics, df = 2 k).
#   Stouffer_directional / Brown / DPM: ActivePathways::merge_p_values with
#   scores_direction = sign matrix, constraints = +1; it does not accept NA, so
#   NA is imputed for these three (p -> 1, sign -> 0), as the package recommends.
#   ORBIT_P: each omic is NA-filtered independently, rho estimated by ORBIT_cor.
#
#   Output (in ../Data):
#     FPR_methods_comparison.pdf        3 FPR panels
#     Power_DPM_vs_ORBIT.pdf            power panel
#     FPR_power_methods_comparison.pdf  FPR row + power row
#     + CSV caches when CACHE = TRUE
# ============================================================================

## ------------------------------------------------------------------------
## Portable setup: locate this script, read/write ../Data
## ------------------------------------------------------------------------
script_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  f <- gsub("~+~", " ", f, fixed = TRUE)
  if (length(f)) return(dirname(normalizePath(f[1])))
  for (e in rev(sys.frames())) if (!is.null(e$ofile)) return(dirname(normalizePath(e$ofile)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
})
setwd(script_dir)
DATA_DIR <- normalizePath(file.path(script_dir, "..", "Data/P_merge"), mustWork = FALSE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
message("Script dir: ", getwd(), "  |  data dir -> ", DATA_DIR)

for (pkg in c("ActivePathways", "ggplot2", "dplyr", "tidyr", "patchwork"))
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
if (!requireNamespace("ORBIT", quietly = TRUE))
  stop("Package 'ORBIT' not installed. Install it first, e.g. remotes::install_github('<user>/ORBIT').")
stopifnot(packageVersion("ActivePathways") >= "2.0.0")

suppressPackageStartupMessages({
  library(ORBIT)
  library(ActivePathways)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ============================================================================
# 0. Config
# ============================================================================
N_FEAT  <- 10000
N_REP   <- 20
ALPHA   <- 0.05
CACHE   <- TRUE          # reuse CSV in DATA_DIR if present
VERBOSE <- TRUE          # per-run progress bar
P_CLAMP <- 1e-15
FONT    <- "Arial"

# --- FPR grids ---
RHO_GRID <- seq(0, 0.9, by = 0.1)
K_GRID   <- 2:20
NA_GRID  <- seq(0, 0.9, by = 0.1)

METHODS <- c("Fisher", "Brown", "Stouffer", "DPM", "ORBIT_P")
pal <- c(Fisher = "#6E5E7B", Brown = "#4F6378", Stouffer = "#6B7A60",
         DPM = "#B07A2A", ORBIT_P = "#E63946")
method_labels <- c(Fisher = "Fisher", Brown = "Brown", Stouffer = "Stouffer",
                   DPM = "DPM", ORBIT_P = "ORBIT")

# --- Power design (mirrors the main-text power simulation) ---
POWER_PROPS      <- c(0.05, 0.10, 0.20, 0.50)   # add 0 for a no-signal baseline (flat zero line)
POWER_K_GRID     <- 2:10                        # power only; FPR K panel still uses K_GRID
POWER_Y_MIN      <- 0.7                         # lower bound of the power y-axis
POWER_MU         <- 3
POWER_RHO_SIGNAL <- 0.9
POWER_METHODS    <- c("DPM", "ORBIT_P")
POWER_LTY        <- c(DPM = "solid", ORBIT_P = "dashed")
# TRUE : DPM covariance and ORBIT rho estimated from true-null features only
#        (oracle; equivalent to rho = 0 in the main-text power script)
# FALSE: both estimated from all features (no oracle knowledge)
POWER_ORACLE_COR <- TRUE

prop_label  <- function(p) ifelse(p == 0, "No signal", paste0(p * 100, "%"))
PROP_LEVELS <- prop_label(POWER_PROPS)
PROP_PAL    <- setNames(c("grey50", "#4C72B0", "#55A868", "#DD8452", "#C44E52")[
  match(PROP_LEVELS, prop_label(c(0, 0.05, 0.10, 0.20, 0.50)))], PROP_LEVELS)

data_path <- function(f) file.path(DATA_DIR, f)

# ============================================================================
# 1. Progress helpers
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
  now     <- Sys.time()
  elapsed <- as.numeric(difftime(now, env$t0, units = "secs"))
  per     <- elapsed / env$i
  eta     <- per * (env$total - env$i)
  frac    <- env$i / env$total
  nbar    <- 24L
  bar     <- paste0(strrep("=", round(frac * nbar)),
                    strrep("-", nbar - round(frac * nbar)))
  line <- sprintf("\r  [%s] %4d/%-4d %5.1f%% | %-34s | %5.2fs/run | el %s | eta %s",
                  bar, env$i, env$total, frac * 100, detail, per,
                  fmt_hms(elapsed), fmt_hms(eta))
  cat(line, file = stderr())
  utils::flush.console()
  if (env$i == env$total) {
    cat("\n", file = stderr())
    message("  done in ", fmt_hms(elapsed), "  (", sprintf("%.2f", per), " s/run)")
  }
  invisible(NULL)
}

# warn if a single run exceeds hb seconds
run_timed <- function(expr, detail, hb = 120) {
  t <- Sys.time()
  res <- force(expr)
  el <- as.numeric(difftime(Sys.time(), t, units = "secs"))
  if (VERBOSE && el > hb) message("\n  [slow] ", detail, " took ", fmt_hms(el))
  res
}

# jobs: data.frame, one row per run; FUN(i) returns a data.frame
run_grid <- function(jobs, label, FUN) {
  pb  <- progress_init(nrow(jobs), label)
  out <- vector("list", nrow(jobs))
  for (i in seq_len(nrow(jobs))) {
    vals   <- vapply(jobs[i, , drop = FALSE], function(v) format(v, trim = TRUE), character(1))
    detail <- paste(sprintf("%s=%s", names(jobs), vals), collapse = " ")
    out[[i]] <- run_timed(FUN(i), detail)
    progress_step(pb, detail)
  }
  do.call(rbind, out)
}

run_or_load <- function(file, expr) {
  path <- data_path(file)
  if (CACHE && file.exists(path)) {
    message("[cache] ", path, " exists -- skipping simulation")
    return(read.csv(path, stringsAsFactors = FALSE))
  }
  res <- force(expr)
  write.csv(res, path, row.names = FALSE)
  message("[write] ", path, "  (", nrow(res), " rows)")
  res
}

# ============================================================================
# 2. Utilities
# ============================================================================
# omics_list -> feature x omic matrices for ActivePathways (NA kept as NA).
omics_to_PS <- function(omics_list) {
  feats <- omics_list[[1]]$feature
  P <- sapply(omics_list, function(df) df$stat[match(feats, df$feature)])
  S <- sapply(omics_list, function(df) df$sign[match(feats, df$feature)])
  dimnames(P) <- dimnames(S) <- list(feats, names(omics_list))
  P <- pmin(pmax(P, P_CLAMP), 1 - P_CLAMP)      # NA stays NA
  list(P = P, S = S)
}

# Fisher by hand, NA-aware: only observed omics enter the sum, df = 2 * n_obs
fisher_na <- function(P) {
  k    <- rowSums(!is.na(P))
  stat <- -2 * rowSums(log(P), na.rm = TRUE)
  p    <- pchisq(stat, df = 2 * k, lower.tail = FALSE)
  p[k == 0] <- NA_real_
  setNames(p, rownames(P))
}

# Stouffer_directional / Brown / DPM via ActivePathways::merge_p_values,
# which does not accept NA -> impute p -> 1, sign -> 0 (package recommendation).
ap_merge <- function(P, S, method) {
  if (method == "Fisher") return(fisher_na(P))
  cons <- rep(1, ncol(P))
  P[is.na(P)] <- 1 - P_CLAMP
  S[is.na(S)] <- 0
  p <- switch(method,
              Brown    = merge_p_values(P, method = "Brown"),
              Stouffer = merge_p_values(P, method = "Stouffer_directional",
                                        scores_direction = S, constraints_vector = cons),
              DPM      = merge_p_values(P, method = "DPM",
                                        scores_direction = S, constraints_vector = cons))
  setNames(as.numeric(p), rownames(P))
}

# mean and 95% t-CI of y_var within groups
agg_ci <- function(df, y_var, groups) {
  df %>%
    group_by(across(all_of(groups))) %>%
    summarise(mean  = mean(.data[[y_var]], na.rm = TRUE),
              sd    = sd(.data[[y_var]],   na.rm = TRUE),
              n     = sum(!is.na(.data[[y_var]])),
              se    = sd / sqrt(n),
              ci_lo = pmax(mean - qt(0.975, df = n - 1) * se, 0),
              ci_hi = mean + qt(0.975, df = n - 1) * se,
              .groups = "drop")
}

# portable PDF writer: quartz on macOS, cairo_pdf elsewhere
save_pdf <- function(p, file, width, height) {
  path <- data_path(file)
  if (capabilities("aqua")) {
    quartz(type = "pdf", file = path, width = width, height = height, family = FONT)
    print(p); dev.off()
  } else {
    ggsave(path, p, width = width, height = height, device = cairo_pdf)
  }
  message("[pdf] ", path)
}

# ============================================================================
# 3. FPR: one simulation -> FPR per method
# ============================================================================
run_one <- function(K = 2, rho = 0, na_rate = 0, seed = 42) {
  sim       <- ORBIT_simulate_null(n = N_FEAT, K = K, rho = rho, seed = seed)
  direction <- setNames(rep(1, K), names(sim$omics_list))
  
  if (na_rate > 0) {                               # NA injected into omic 1 only
    set.seed(seed + 10000L)
    n_na <- round(na_rate * N_FEAT)
    if (n_na > 0) {
      na_idx <- sample(seq_len(N_FEAT), n_na)
      sim$omics_list[[1]]$stat[na_idx] <- NA
      sim$omics_list[[1]]$sign[na_idx] <- NA
    }
  }
  
  ps   <- omics_to_PS(sim$omics_list)
  p_ap <- lapply(c("Fisher", "Brown", "Stouffer", "DPM"),
                 function(m) ap_merge(ps$P, ps$S, m))
  names(p_ap) <- c("Fisher", "Brown", "Stouffer", "DPM")
  
  # ORBIT_P: per-omic NA filtering, rho estimated from the data
  omics_clean   <- lapply(sim$omics_list, function(df)
    df[!is.na(df$stat) & !is.na(df$sign), , drop = FALSE])
  omics_for_cor <- lapply(omics_clean, function(df) { df$significance <- 0L; df })
  rho_est <- suppressMessages(ORBIT_cor(omics_for_cor, direction, sig_mode = "column"))$rho
  res_P   <- suppressMessages(ORBIT_P(omics_clean, direction, rho = rho_est))
  
  data.frame(method  = METHODS,
             rho_est = rho_est,
             fpr     = c(sapply(p_ap, function(p) mean(p < ALPHA, na.rm = TRUE)),
                         mean(res_P$P < ALPHA, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}

# ============================================================================
# 4. FPR grids (CSV cache, progress bar)
# ============================================================================
df_rho <- run_or_load("FPR_methods_rho.csv", {
  jobs <- expand.grid(rho = RHO_GRID, rep = seq_len(N_REP))
  run_grid(jobs, "FPR 1/3  varying rho", function(i) {
    res <- run_one(K = 2, rho = jobs$rho[i], na_rate = 0, seed = jobs$rep[i])
    res$rho <- jobs$rho[i]; res$rep <- jobs$rep[i]; res
  })
})

df_K <- run_or_load("FPR_methods_K.csv", {
  jobs <- expand.grid(K = K_GRID, rep = seq_len(N_REP))
  run_grid(jobs, "FPR 2/3  varying K", function(i) {
    res <- run_one(K = jobs$K[i], rho = 0, na_rate = 0, seed = jobs$rep[i])
    res$K <- jobs$K[i]; res$rep <- jobs$rep[i]; res
  })
})

df_NA <- run_or_load("FPR_methods_NA.csv", {
  jobs <- expand.grid(na_rate = NA_GRID, rep = seq_len(N_REP))
  run_grid(jobs, "FPR 3/3  varying NA rate", function(i) {
    res <- run_one(K = 2, rho = 0, na_rate = jobs$na_rate[i], seed = jobs$rep[i])
    res$na_rate <- jobs$na_rate[i]; res$rep <- jobs$rep[i]; res
  })
})

# ============================================================================
# 5. Power: DPM vs ORBIT_P on the same simulated data
# ============================================================================
# DPM for the power study. Unlike the FPR part, the covariance can be estimated
# from a chosen subset of features (cov_rows), so that the oracle setting
# (true-null features only) is possible. cov_rows = NULL -> all features,
# identical to merge_p_values(P, method = "DPM", ...).
dpm_power <- function(P, S, cov_rows = NULL) {
  cons <- rep(1, ncol(P))
  if (is.null(cov_rows))
    return(setNames(as.numeric(merge_p_values(P, method = "DPM",
                                              scores_direction = S,
                                              constraints_vector = cons)),
                    rownames(P)))
  cov_mat <- ActivePathways:::calculateCovariances(t(P[cov_rows, , drop = FALSE]))
  dimnames(cov_mat) <- list(colnames(P), colnames(P))
  p <- ActivePathways::DPM(P, cov_matrix = cov_mat,
                           scores_direction = S, constraints_vector = cons)
  setNames(as.numeric(p), rownames(P))
}

power_one <- function(K, signal_prop, seed) {
  set.seed(seed)
  sim <- ORBIT_simulate_signal(
    n = N_FEAT, K = K,
    pi_null = 1 - signal_prop,
    pi_up   = signal_prop / 2,
    pi_down = signal_prop / 2,
    mu_signal = POWER_MU, rho_null = 0, rho_signal = POWER_RHO_SIGNAL
  )
  direction <- setNames(rep(1, K), names(sim$omics_list))
  truth_all <- sim$is_signal
  if (is.null(names(truth_all))) names(truth_all) <- sim$omics_list[[1]]$feature
  
  ## --- DPM ---
  ps       <- omics_to_PS(sim$omics_list)          # no NA in power sims
  cov_rows <- if (POWER_ORACLE_COR) which(!truth_all[rownames(ps$P)]) else NULL
  p_dpm    <- dpm_power(ps$P, ps$S, cov_rows = cov_rows)
  
  ## --- ORBIT_P ---
  omics_for_cor <- lapply(sim$omics_list, function(df) {
    df$significance <- if (POWER_ORACLE_COR) as.integer(truth_all[df$feature]) else 0L
    df
  })
  rho_est <- suppressMessages(ORBIT_cor(omics_for_cor, direction, sig_mode = "column"))$rho
  res_P   <- suppressMessages(ORBIT_P(sim$omics_list, direction, rho = rho_est))
  stopifnot(all(c("Feature", "P") %in% names(res_P)))
  p_orbit <- setNames(res_P$P, as.character(res_P$Feature))
  
  ## --- evaluate: BH at ALPHA, power = TPR, fdp = FDP ---
  eval_power <- function(p, method, rho_est = NA_real_) {
    truth <- truth_all[names(p)]
    if (anyNA(truth)) stop(method, ": features absent from is_signal.")
    det <- p.adjust(p, method = "BH") < ALPHA
    det[is.na(det)] <- FALSE
    data.frame(method = method, K = K, signal_prop = signal_prop, rho_est = rho_est,
               power = sum(det & truth) / max(sum(truth), 1),
               fdp   = if (sum(det) > 0) sum(det & !truth) / sum(det) else 0,
               stringsAsFactors = FALSE)
  }
  rbind(eval_power(p_dpm,   "DPM"),
        eval_power(p_orbit, "ORBIT_P", rho_est))
}

df_power <- run_or_load("Power_DPM_ORBIT_K_signalprop.csv", {
  jobs <- expand.grid(K = POWER_K_GRID, signal_prop = POWER_PROPS, rep = seq_len(N_REP))
  run_grid(jobs, "Power  DPM vs ORBIT  varying K x signal prop", function(i) {
    r <- tryCatch(
      power_one(jobs$K[i], jobs$signal_prop[i], seed = 4000L + jobs$rep[i]),
      error = function(e) {
        message("\n  [FAIL] K=", jobs$K[i], " prop=", jobs$signal_prop[i],
                " rep=", jobs$rep[i], " : ", conditionMessage(e))
        data.frame(method = POWER_METHODS, K = jobs$K[i], signal_prop = jobs$signal_prop[i],
                   rho_est = NA_real_, power = NA_real_, fdp = NA_real_,
                   stringsAsFactors = FALSE)
      })
    r$rep <- jobs$rep[i]; r
  })
})

# ============================================================================
# 6. Aggregate
# ============================================================================
aggregate_fpr <- function(df, x_var) {
  agg_ci(df, "fpr", c(x_var, "method")) %>%
    mutate(method = factor(method, levels = METHODS))
}
sum_rho <- aggregate_fpr(df_rho, "rho")
sum_K   <- aggregate_fpr(df_K,   "K")
sum_NA  <- aggregate_fpr(df_NA,  "na_rate")

# keep only POWER_K_GRID (so an older cache with a wider K range still works)
df_power  <- df_power %>% filter(K %in% POWER_K_GRID)
sum_power <- agg_ci(df_power, "power", c("K", "signal_prop", "method")) %>%
  mutate(method       = factor(method, levels = POWER_METHODS),
         signal_label = factor(prop_label(signal_prop), levels = PROP_LEVELS))
sum_fdp   <- agg_ci(df_power, "fdp", c("K", "signal_prop", "method"))

# console tables
for (s in list(sum_rho, sum_K, sum_NA))
  print(s %>% mutate(label = sprintf("%.4f [%.4f, %.4f]", mean, ci_lo, ci_hi)) %>%
          select(1, method, label) %>% pivot_wider(names_from = method, values_from = label))
message("\nPower (mean [95% CI]) by K x signal proportion:")
print(sum_power %>% mutate(label = sprintf("%.3f [%.3f, %.3f]", mean, ci_lo, ci_hi)) %>%
        select(K, signal_prop, method, label) %>%
        pivot_wider(names_from = method, values_from = label), n = Inf)
message("\nMean FDP (should be <= ", ALPHA, "):")
print(sum_fdp %>% select(K, signal_prop, method, mean) %>%
        pivot_wider(names_from = method, values_from = mean), n = Inf)

# ============================================================================
# 7. FPR panels: one row per method (5 x 3 grid), fixed y-axis for all rows
# ============================================================================
FPR_Y_MAX    <- 0.15
FPR_Y_BREAKS <- seq(0, FPR_Y_MAX, by = 0.05)

make_fpr_panel <- function(d, x_var, x_label, x_breaks, m,
                           show_y_title = FALSE, show_x_title = FALSE, row_title = NULL) {
  ggplot(d, aes(x = .data[[x_var]], y = mean)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = pal[[m]], alpha = 0.18) +
    geom_line(color = pal[[m]], linewidth = 1.5) +
    geom_hline(yintercept = ALPHA, linetype = "dashed", color = "grey40") +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = FPR_Y_BREAKS, expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, FPR_Y_MAX)) +   # clip, don't drop, values above the cap
    labs(x = if (show_x_title) x_label else NULL,
         y = if (show_y_title) "False Positive Rate" else NULL,
         title = row_title) +
    theme_minimal(base_size = 15, base_family = FONT) +
    theme(panel.grid   = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          aspect.ratio = 1,
          axis.text    = element_text(size = 12),
          axis.title   = element_text(size = 20),
          plot.title   = element_text(size = 22, face = "bold", hjust = 0,
                                      color = pal[[m]], family = FONT))
}

fpr_panels <- unlist(lapply(METHODS, function(m) {
  bottom <- m == METHODS[length(METHODS)]
  list(
    make_fpr_panel(filter(sum_rho, method == m), "rho", "Inter-Omics Correlation",
                   RHO_GRID, m, show_y_title = TRUE, show_x_title = bottom,
                   row_title = method_labels[[m]]),
    make_fpr_panel(filter(sum_K,   method == m), "K", "Number of Omics (K)",
                   seq(2, 20, by = 2), m, show_x_title = bottom),
    make_fpr_panel(filter(sum_NA,  method == m), "na_rate", "NA Rate in Omic 1",
                   NA_GRID, m, show_x_title = bottom)
  )
}), recursive = FALSE)

fig_fpr <- wrap_plots(fpr_panels, ncol = 3, byrow = TRUE) &
  theme(plot.margin = margin(4, 8, 4, 4))

# ============================================================================
# 8. Power panel (color = signal proportion, linetype = method)
# ============================================================================
make_power_plot <- function(d) {
  ggplot(d, aes(x = K, y = mean,
                color = signal_label, fill = signal_label, linetype = method,
                group = interaction(method, signal_label))) +
    geom_ribbon(aes(ymin = pmax(ci_lo, 0), ymax = pmin(ci_hi, 1)),
                alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 1.6) +
    scale_x_continuous(breaks = POWER_K_GRID) +
    scale_y_continuous(breaks = seq(POWER_Y_MIN, 1, by = 0.1),
                       expand = expansion(mult = c(0, 0.02))) +
    # clip (not drop) anything below POWER_Y_MIN
    coord_cartesian(ylim = c(POWER_Y_MIN, 1), expand = TRUE) +
    scale_color_manual(values = PROP_PAL, limits = PROP_LEVELS, name = "Signal proportion") +
    scale_fill_manual( values = PROP_PAL, limits = PROP_LEVELS, name = "Signal proportion") +
    scale_linetype_manual(values = POWER_LTY, breaks = POWER_METHODS,
                          labels = method_labels[POWER_METHODS], name = "Method") +
    guides(color    = guide_legend(order = 1, ncol = 1,
                                   override.aes = list(linewidth = 2, linetype = "solid")),
           fill     = "none",
           linetype = guide_legend(order = 2, ncol = 1,
                                   override.aes = list(color = "black", linewidth = 1.3))) +
    labs(x = "Number of Omics (K)", y = "Power") +
    theme_minimal(base_size = 15, base_family = FONT) +
    theme(panel.grid        = element_blank(),
          panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.6),
          aspect.ratio      = 1,
          legend.position   = "right",
          legend.direction  = "vertical",
          legend.box        = "vertical",
          legend.title      = element_text(family = FONT, face = "bold", size = 20),
          legend.text       = element_text(family = FONT, size = 18),
          legend.key.width  = unit(1.6, "cm"),   # wide enough to show dashes
          legend.spacing.y  = unit(0.5, "cm"),
          axis.text         = element_text(size = 16),
          axis.title        = element_text(size = 24, color = "grey25"))
}

p_power <- make_power_plot(sum_power)

# combined: 5 FPR rows on top, power below
fig_all <- wrap_plots(fig_fpr, p_power, ncol = 1, heights = c(5, 1.7)) &
  theme(plot.margin = margin(5, 5, 1.5, 1.5))

# ============================================================================
# 9. Render / save
# ============================================================================
message("\n=== Rendering figures ===")
print(fig_fpr)
print(p_power)

save_pdf(fig_fpr, "FPR_methods_comparison.pdf",       width = 14, height = 24)
save_pdf(p_power, "Power_DPM_vs_ORBIT.pdf",           width = 9,  height = 7.5)
save_pdf(fig_all, "FPR_power_methods_comparison.pdf", width = 14, height = 33)
