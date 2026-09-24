###############################################################################
##  Supplementary Figure 1 -- validity of the variance-gamma reference under
##  the continuous-rank null
##
##  ORBIT assumes u_j = r_j/(n_j+1) ~ U(0,1), x_j = -log(u_j) ~ Exp(1), and
##  refers D = sum_j s_j x_j (s_j = +-1) to VG(r = K, sigma = 1) at rho = 0.
##  Question: does the analytical tail as implemented in ORBIT_Rank (closed form
##  at K = 2, 32-node Gauss-Laguerre otherwise) equal the true tail of D under
##  this continuous model?  Rank discreteness is Supplementary Figure 2.
##
##  Panel: K = 2 / 3 / 5.  y = P_VG / P_MC vs x = -log10 P_VG, 1e-2 .. 1e-7.
##
##  Monte Carlo engines under the same model:
##    (i)  naive: draw u_j, s_j, form D, count exceedances.  Unbiased but needs
##         ~1e9 draws at P = 1e-7; used only for the cross-check.
##    (ii) conditional: s_j x_j is Laplace(0,1) = E - E', so D = G1 - G2 with
##         G1, G2 ~ Gamma(K, 1) and P(|D| > d) = 2 E[Q(K, d + G2)].  Only G2 is
##         simulated; relative SE ~ sqrt((4/3)^K / M) at any depth; uses
##         rgamma + pgamma only.  Drawn as lines with normal bands.
##
##  Console cross-checks:
##    (a) ORBIT_Rank quadrature vs terminating closed form (integer K)
##    (b) conditional vs naive Monte Carlo at K = 3
##
##  Deps: statmod, ggplot2
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

for (pkg in c("statmod", "ggplot2"))
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

set.seed(20260825)

K_PANEL <- c(2, 3, 5)
M_CMC   <- 1e6           # conditional MC draws per K
P_RANGE <- c(-2, -7)     # log10 P_VG range on the x axis
FONT    <- "Arial"

###############################################################################
## 1. analytical reference, copied verbatim from ORBIT_Rank
##    (rho = 0  =>  r_eff = K, sigma_eff = 1)
###############################################################################
GL       <- statmod::gauss.quad(32, kind = "laguerre")
GL_NODES <- GL$nodes
GL_LOG_W <- log(GL$weights)

log_besselK_safe <- function(z, nu) {
  if (!is.finite(z) || z <= 0) return(NA_real_)
  val <- suppressWarnings(besselK(z, nu, expon.scaled = TRUE))
  if (is.finite(val) && val > 0) return(log(val) - z)
  NA_real_
}
log_vg_density_local <- function(ax, r, sigma) {
  if (!is.finite(ax) || ax <= 0) return(NA_real_)
  nu <- r - 0.5
  z  <- ax / sigma
  log_pref  <- -log(sigma) - 0.5 * log(pi) - lgamma(r)
  log_power <- nu * (log(ax) - log(2 * sigma))
  logK <- log_besselK_safe(z, nu)
  if (!is.finite(logK)) return(NA_real_)
  log_pref + log_power + logK
}
get_log_tail_prob <- function(abs_d, r, sigma) {
  if (!is.finite(abs_d) || r <= 0 || sigma <= 0) return(NA_real_)
  if (abs_d <= 1e-10) return(0)
  if (abs(r - 1) < 1e-3) return(-abs_d / sigma)
  if (abs(r - 2) < 1e-3) {
    z <- abs_d / sigma
    return(log(0.5) + log(z + 2) - z)
  }
  log_vals <- vapply(seq_along(GL_NODES), function(i) {
    u <- GL_NODES[i]
    x <- abs_d + sigma * u
    lg <- log_vg_density_local(x, r, sigma)
    if (!is.finite(lg)) return(NA_real_)
    lg + u + GL_LOG_W[i]
  }, numeric(1))
  if (any(is.na(log_vals))) return(NA_real_)
  m <- max(log_vals)
  if (!is.finite(m)) return(NA_real_)
  log(2) + log(sigma) + m + log(sum(exp(log_vals - m)))
}
vg_P <- function(d, K) exp(get_log_tail_prob(d, K, 1))

## terminating closed form for integer K; deterministic check of the quadrature
vg_tail_closed <- function(d, K) {
  lt <- numeric(K); lc <- 0
  for (k in 0:(K - 1)) {
    lt[k + 1] <- lc + pgamma(d, shape = K - k, lower.tail = FALSE, log.p = TRUE)
    lc <- lc + log(K + k) - log(k + 1) - log(2)
  }
  m <- max(lt)
  exp((1 - K) * log(2) + m + log(sum(exp(lt - m))))
}

## |D| threshold at which the production P equals p
d_for_p <- function(p, K)
  uniroot(function(x) get_log_tail_prob(x, K, 1) - log(p),
          c(1e-8, 200), tol = 1e-12)$root

###############################################################################
## 2. conditional Monte Carlo: P(|D| > d) = 2 E[Q(K, d + G2)], G2 ~ Gamma(K)
##    one Gamma sample shared across thresholds (smooth curves)
###############################################################################
cmc_tail <- function(dvec, K, M) {
  G <- rgamma(M, shape = K, rate = 1)
  t(vapply(dvec, function(d) {
    h <- 2 * pgamma(d + G, shape = K, rate = 1, lower.tail = FALSE)
    c(p = mean(h), se = sd(h) / sqrt(M))
  }, numeric(2)))
}

###############################################################################
## 3. naive Monte Carlo of the literal generative process, one pass
###############################################################################
mc_counts <- function(dvec, K, M, batch = 2e6) {
  o <- order(dvec); ds <- dvec[o]
  nb <- length(ds) + 1L; tab <- integer(nb); left <- M
  while (left > 0) {
    m <- min(batch, left)
    D <- numeric(m)
    for (j in seq_len(K))
      D <- D + sample(c(-1, 1), m, replace = TRUE) * (-log(runif(m)))
    tab  <- tab + tabulate(findInterval(abs(D), ds) + 1L, nbins = nb)
    left <- left - m
  }
  cnt <- numeric(length(dvec))
  cnt[o] <- rev(cumsum(rev(tab)))[-1L]
  cnt
}

###############################################################################
## 4. cross-checks
###############################################################################
cat("cross-check (a): ORBIT_Rank quadrature vs terminating closed form\n")
for (K in K_PANEL) {
  dchk <- c(2, 5, 10, 15, 20, 25)
  rel  <- vapply(dchk, function(d) vg_P(d, K) / vg_tail_closed(d, K) - 1, numeric(1))
  cat(sprintf("  K = %d  max |rel dev| = %.2e   (d = %s)\n",
              K, max(abs(rel)), paste(dchk, collapse = ",")))
}

cat("\ncross-check (b): conditional vs naive Monte Carlo (K = 3, M_naive = 2e7)\n")
pchk <- 10^seq(-2, -5, by = -1)
dchk <- vapply(pchk, d_for_p, numeric(1), K = 3)
cchk <- mc_counts(dchk, 3, M = 2e7)
mchk <- cmc_tail(dchk, 3, M = 1e6)
for (i in seq_along(dchk))
  cat(sprintf("  P_VG=%.0e  d=%6.3f  conditional %.5e +- %.1e   naive %.5e  (%s events)\n",
              pchk[i], dchk[i], mchk[i, "p"], mchk[i, "se"],
              cchk[i] / 2e7, format(cchk[i])))

###############################################################################
## 5. panel data
###############################################################################
cat("\nbuilding panel (conditional MC, M = ", format(M_CMC, scientific = TRUE), " per K)\n", sep = "")
dat <- do.call(rbind, lapply(K_PANEL, function(K) {
  pv <- 10^seq(P_RANGE[1], P_RANGE[2], length.out = 60)
  dv <- vapply(pv, d_for_p, numeric(1), K = K)
  mc <- cmc_tail(dv, K, M_CMC)
  cat(sprintf("  K = %d done\n", K))
  data.frame(K = K, d = dv, vg_P = pv, mc_P = mc[, "p"], mc_se = mc[, "se"],
             ratio    = pv / mc[, "p"],
             ratio_lo = pv / (mc[, "p"] + 1.96 * mc[, "se"]),
             ratio_hi = pv / (mc[, "p"] - 1.96 * mc[, "se"]))
}))
dat$grp <- factor(dat$K, levels = K_PANEL, labels = paste("K =", K_PANEL))

###############################################################################
## 6. figure
###############################################################################
library(ggplot2)

th <- theme_bw(base_size = 15, base_family = FONT) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 14, colour = "black"),
        legend.position = c(0.03, 0.97),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key   = element_blank(),
        legend.title = element_text(size = 13),
        legend.text  = element_text(size = 13),
        plot.title   = element_text(size = 15, face = "plain"))

PAL <- c("#1f78b4", "#8073ac", "#c51b7d")

p1 <- ggplot(dat, aes(-log10(vg_P), ratio, colour = grp, group = grp)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  geom_ribbon(aes(ymin = ratio_lo, ymax = ratio_hi, fill = grp),
              alpha = 0.18, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  scale_y_log10(breaks = c(0.8, 0.9, 1, 1.1, 1.25)) +
  coord_cartesian(ylim = c(0.75, 1.33)) +
  scale_x_continuous(breaks = seq(2, 7, 1)) +
  scale_colour_manual(values = PAL, name = "omics layers") +
  scale_fill_manual(values = PAL) +
  labs(x = expression(-log[10]~italic(P)[VG]),
       y = expression(italic(P)[VG]/italic(P)[MC]),
       title = "Continuous-rank null") + th

quartz(file = file.path(OUT_DIR, "uniform-VG-1.pdf"), type = "pdf", width = 6, height = 6)
print(p1)
dev.off()