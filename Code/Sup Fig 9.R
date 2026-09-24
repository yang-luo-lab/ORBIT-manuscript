###############################################################################
##  Supplemental Fig. 9
##  ORBIT vs VarianceGamma::pvg, arbitrated by Monte Carlo
##
##  MC samples the null directly (no special functions, no quadrature), so it
##  shares no failure mode with either method; its only error is sampling noise,
##  reported as a 95% CI.
##  x = Monte Carlo P, y = analytical P (ORBIT, pvg).
##  K = 2, rho = 0.3 -> r = 1.5385, sigma = 1.3 (rho = 0 would give r = 2, a
##  closed form in ORBIT that never touches the quadrature).
##  Deps: statmod (required), VarianceGamma (optional), ggplot2, scales
###############################################################################

## ------------------------------------------------------------------------
## Portable setup: locate this script, write PDFs to ../Data
## Layout: <repo>/Code/<this script>, <repo>/Data/
## Works with Rscript, source(), and RStudio Source/Run.
## ------------------------------------------------------------------------
script_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  f <- gsub("~+~", " ", f, fixed = TRUE)
  if (length(f)) return(dirname(normalizePath(f[1])))               # Rscript
  for (e in rev(sys.frames())) if (!is.null(e$ofile)) return(dirname(normalizePath(e$ofile)))  # source()
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path                  # RStudio
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
})
setwd(script_dir)
OUT_DIR <- normalizePath(file.path(script_dir, "..", "Data"), mustWork = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("Script dir: ", getwd(), "  |  outputs -> ", OUT_DIR)

## dependencies (installed only if missing)
for (pkg in c("statmod", "VarianceGamma", "ggplot2", "scales"))
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
if (!requireNamespace("ORBIT", quietly = TRUE))
  stop("Package 'ORBIT' not installed. Install it first, e.g. remotes::install_github('<user>/ORBIT').")

library(ORBIT)
library(VarianceGamma)
library(statmod)

set.seed(20260724)

K <- 2; RHO <- 0.3
LAMBDA  <- 1 + (K - 1) * RHO
R_SHAPE <- K / LAMBDA               # 1.53846
SIGMA   <- LAMBDA                   # 1.3

D_GRID  <- c(7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 44)  # P ~ 1e-2 .. 1e-14
N_PLAIN <- 5e8                      # crude MC, single pass
N_TILT  <- 4e6                      # tilted IS draws per point

###############################################################################
## 1. ORBIT evaluator (port of get_log_tail_prob() from ORBIT_Rank)
###############################################################################
GL <- statmod::gauss.quad(32, kind = "laguerre")
GLn <- GL$nodes; GLlw <- log(GL$weights)

orbit_log_tail <- function(abs_d, r, sigma) {
  if (!is.finite(abs_d) || r <= 0 || sigma <= 0) return(NA_real_)
  if (abs_d <= 1e-10) return(0)
  if (abs(r - 1) < 1e-3) return(-abs_d / sigma)
  if (abs(r - 2) < 1e-3) { z <- abs_d / sigma; return(log(0.5) + log(z + 2) - z) }
  nu <- r - 0.5
  lv <- vapply(seq_along(GLn), function(i) {
    u <- GLn[i]; x <- abs_d + sigma * u; z <- x / sigma
    k <- suppressWarnings(besselK(z, nu, expon.scaled = TRUE))
    if (!is.finite(k) || k <= 0) return(NA_real_)
    -log(sigma) - 0.5 * log(pi) - lgamma(r) +
      nu * (log(x) - log(2 * sigma)) + log(k) - z + u + GLlw[i]
  }, numeric(1))
  if (any(is.na(lv))) return(NA_real_)
  m <- max(lv)
  log(2) + log(sigma) + m + log(sum(exp(lv - m)))
}

###############################################################################
## 2. VarianceGamma::pvg
##    X = vgC + theta*G + sigma_p*sqrt(G)*Z, G ~ Gamma(1/nu, scale = nu)
##    -> vgC = 0, theta = 0, sigma_p = sigma*sqrt(2r), nu = 1/r
##    lower.tail = FALSE not implemented: upper tail = 1 - pvg(d)
###############################################################################
HAVE_VG <- requireNamespace("VarianceGamma", quietly = TRUE)
if (!HAVE_VG) message("VarianceGamma not installed; pvg columns will be NA.")

pkg_p <- function(d, r, sigma, tiny = 1e-10) {
  if (!HAVE_VG) return(NA_real_)
  o <- try(suppressWarnings(
    VarianceGamma::pvg(d, vgC = 0, sigma = sigma * sqrt(2 * r), theta = 0,
                       nu = 1 / r, lower.tail = TRUE, tiny = tiny)),
    silent = TRUE)
  if (inherits(o, "try-error") || !is.finite(o)) return(NA_real_)
  2 * (1 - o)
}

###############################################################################
## 3. Monte Carlo
##    D = signed sum of K iid Exp(1) = sigma*(G1 - G2), G_i ~ Gamma(r, 1)
##    (verified in 3c)
###############################################################################

## 3a. crude MC, all thresholds in one pass
mc_plain_all <- function(dvec, r, sigma, n, batch = 5e7) {
  cnt <- numeric(length(dvec)); left <- n
  while (left > 0) {
    m <- min(batch, left)
    A <- abs(sigma * (rgamma(m, shape = r, scale = 1) -
                        rgamma(m, shape = r, scale = 1)))
    cnt <- cnt + vapply(dvec, function(d) sum(A > d), numeric(1))
    left <- left - m
  }
  lo <- hi <- numeric(length(cnt))
  for (i in seq_along(cnt)) {
    ci <- if (cnt[i] > 0) poisson.test(cnt[i])$conf.int / n else c(0, 3.689 / n)
    lo[i] <- ci[1]; hi[i] <- ci[2]
  }
  data.frame(p = cnt / n, lo = lo, hi = hi, events = cnt)
}

## 3b. exponentially tilted IS, unbiased at any depth
##     saddlepoint t* = (sqrt(r^2 + a^2) - r)/a, a = d/sigma
mc_tilted <- function(d, r, sigma, n) {
  a  <- d / sigma
  t  <- min((sqrt(r^2 + a^2) - r) / a, 1 - 1e-13)
  th <- t / sigma
  D  <- sigma * (rgamma(n, shape = r, scale = 1 / (1 - t)) -
                   rgamma(n, shape = r, scale = 1 / (1 + t)))
  keep <- D > d
  if (!any(keep)) return(c(p = NA, lo = NA, hi = NA, rse = NA))
  lw <- -th * D[keep]; M <- max(lw); s <- exp(lw - M)
  ms <- sum(s) / n; ms2 <- sum(s^2) / n          # divide by n, not sum(keep)
  rse <- sqrt(max(ms2 - ms^2, 0) / n) / ms
  p <- exp(-r * log1p(-t^2) + M + log(ms) + log(2))
  c(p = p, lo = p * (1 - 1.96 * rse), hi = p * (1 + 1.96 * rse), rse = rse)
}

## 3c. sampler self-check: Gamma difference vs literal signed-Exp(1) sum
sampler_check <- function(n = 2e7) {
  r <- 2; sigma <- 1                       # integer r: both forms exact
  A <- abs(sigma * (rgamma(n, shape = r, scale = 1) -
                      rgamma(n, shape = r, scale = 1)))
  E <- matrix(rexp(n * r), ncol = r)
  S <- matrix(sample(c(-1, 1), n * r, replace = TRUE), ncol = r)
  B <- abs(sigma * rowSums(S * E))
  cat("sampler self-check (r = 2, sigma = 1):\n")
  for (d in c(3, 6, 9, 12))
    cat(sprintf("  d=%4.1f  Gamma-diff %.4e   signed-Exp sum %.4e   ratio %.4f\n",
                d, mean(A > d), mean(B > d), mean(A > d) / mean(B > d)))
}

###############################################################################
## 4. Run
###############################################################################
cat(sprintf("K = %d, rho = %.2f  ->  r = %.5f, sigma = %.2f\n",
            K, RHO, R_SHAPE, SIGMA))
cat(sprintf("pvg args: sigma = %.6f, nu = %.6f\n\n",
            SIGMA * sqrt(2 * R_SHAPE), 1 / R_SHAPE))
sampler_check()

cat(sprintf("\ncrude MC: %.0e draws (single pass)\n", N_PLAIN))
pl <- mc_plain_all(D_GRID, R_SHAPE, SIGMA, n = N_PLAIN)
ti <- t(vapply(D_GRID, mc_tilted, numeric(4), r = R_SHAPE, sigma = SIGMA, n = N_TILT))

res <- data.frame(
  D          = D_GRID,
  mc_p       = ti[, "p"],            # tilted IS = reference
  mc_lo      = ti[, "lo"],
  mc_hi      = ti[, "hi"],
  mc_rse     = ti[, "rse"],
  crude_p    = pl$p,                 # crude MC validates the IS estimator
  crude_ev   = pl$events,
  orbit_p    = exp(vapply(D_GRID, orbit_log_tail, numeric(1),
                          r = R_SHAPE, sigma = SIGMA)),
  pvg_p      = vapply(D_GRID, pkg_p, numeric(1), r = R_SHAPE, sigma = SIGMA),
  pvg_tuned  = vapply(D_GRID, pkg_p, numeric(1), r = R_SHAPE, sigma = SIGMA,
                      tiny = 1e-30))
res$orbit_over_mc <- res$orbit_p / res$mc_p
res$pvg_over_mc   <- ifelse(!is.na(res$pvg_p) & res$pvg_p > 0,
                            res$pvg_p / res$mc_p, NA)

cat("\n")
for (i in seq_len(nrow(res)))
  cat(sprintf("D=%5.1f | MC %10.3e (+-%.1f%%) | crude %10.3e (%s ev) | ORBIT %10.3e (x%.4f) | pvg %s\n",
              res$D[i], res$mc_p[i], 100 * 1.96 * res$mc_rse[i],
              res$crude_p[i], format(res$crude_ev[i]),
              res$orbit_p[i], res$orbit_over_mc[i],
              if (is.na(res$pvg_p[i])) "NA"
              else if (res$pvg_p[i] <= 0) "*** 0 ***"
              else sprintf("%10.3e (x%.4f)", res$pvg_p[i], res$pvg_over_mc[i])))

###############################################################################
## 5. Plots
##    A: -log10 computed P vs -log10 MC P
##    B: relative deviation (computed / MC) against the MC noise band
###############################################################################
library(ggplot2)
library(scales)

FONT <- "Arial"      # falls back to default sans if unavailable



mk <- function(nm, p) {
  zero <- !is.na(p) & p <= 0
  data.frame(
    method = nm,
    D      = res$D,
    x      = -log10(res$mc_p),
    xlo    = -log10(res$mc_hi),        # -log10 flips the interval
    xhi    = -log10(res$mc_lo),
    y      = ifelse(zero | is.na(p), NA_real_, -log10(p)),
    ratio  = ifelse(zero | is.na(p), NA_real_, p / res$mc_p),
    zero   = zero,
    stringsAsFactors = FALSE)
}
dd <- rbind(mk("ORBIT", res$orbit_p),
            mk("VarianceGamma::pvg", res$pvg_p))
dd$method <- factor(dd$method, levels = c("ORBIT", "VarianceGamma::pvg"))

band <- data.frame(x    = -log10(res$mc_p),
                   ymin = 1 - 1.96 * res$mc_rse,
                   ymax = 1 + 1.96 * res$mc_rse)

PAL <- c("ORBIT" = "#1f78b4", "VarianceGamma::pvg" = "#e31a1c")

base_theme <- theme_bw(base_size = 15, base_family = FONT) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
        axis.title       = element_text(size = 16),
        axis.text        = element_text(size = 14, colour = "black"),
        legend.position  = c(0.03, 0.97),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key       = element_blank(),
        legend.title     = element_blank(),
        legend.text      = element_text(size = 13),
        plot.title       = element_text(size = 15, face = "plain"),
        plot.tag         = element_text(size = 18, face = "bold"))

ttl <- sprintf("K = %d, rho = %.1f  (r = %.3f, sigma = %.1f)",
               K, RHO, R_SHAPE, SIGMA)

## Panel A
zt <- max(dd$x, na.rm = TRUE) + 1.2          # row for pvg zeros
zA <- subset(dd, zero)

pA <- ggplot(dd, aes(x, y, colour = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
  geom_errorbarh(aes(xmin = xlo, xmax = xhi), height = 0,
                 colour = "grey60", linewidth = 0.5, na.rm = TRUE) +
  geom_point(size = 3.2, na.rm = TRUE) +
  { if (nrow(zA))
    geom_point(data = zA, aes(x = x, y = zt), shape = 4, size = 3.4,
               stroke = 1.1, colour = PAL[2], inherit.aes = FALSE) } +
  { if (nrow(zA))
    annotate("text", x = mean(zA$x), y = zt + 1.0, label = "pvg returned 0",
             colour = PAL[2], size = 5, family = FONT) } +
  scale_colour_manual(values = PAL) +
  scale_shape_manual(values = c(16, 17)) +
  coord_equal(xlim = c(1.5, zt + 2.2), ylim = c(1.5, zt + 2.2)) +
  labs(x = expression(-log[10]~italic(P)[Monte~Carlo]),
       y = expression(-log[10]~italic(P)[computed]),
       title = ttl) +
  base_theme

quartz(file = file.path(OUT_DIR, "MC-VG.pdf"), type = "pdf", width = 6, height = 6)
print(pA)
dev.off()

## Panel B
yrng <- range(c(dd$ratio, band$ymin, band$ymax), na.rm = TRUE)
yrng <- yrng + c(-1, 1) * 0.04 * diff(yrng)
zB   <- subset(dd, zero)

pB <- ggplot() +
  geom_ribbon(data = band, aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey85") +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  geom_line(data = dd, aes(x, ratio, colour = method), linewidth = 0.6,
            na.rm = TRUE) +
  geom_point(data = dd, aes(x, ratio, colour = method, shape = method),
             size = 3.2, na.rm = TRUE) +
  { if (nrow(zB))
    geom_point(data = zB, aes(x = x, y = yrng[1]), shape = 4, size = 3.4,
               stroke = 1.1, colour = PAL[2], inherit.aes = FALSE) } +
  { if (nrow(zB))
    annotate("text", x = mean(zB$x), y = yrng[1] + 0.06 * diff(yrng),
             label = "pvg returned 0", colour = PAL[2], size = 5,
             family = FONT) } +
  scale_colour_manual(values = PAL) +
  scale_shape_manual(values = c(16, 17)) +
  scale_y_continuous(labels = function(v) percent(v - 1, accuracy = 0.1)) +
  coord_cartesian(xlim = c(1.5, max(dd$x) + 1.6), ylim = yrng) +
  labs(x = expression(-log[10]~italic(P)[Monte~Carlo]),
       y = "relative deviation from Monte Carlo",
       title = "grey band = Monte Carlo 95% CI") +
  base_theme +
  theme(legend.position = c(0.03, 0.75))

quartz(file = file.path(OUT_DIR, "MC-VG-2.pdf"), type = "pdf", width = 6, height = 6)
print(pB)
dev.off()