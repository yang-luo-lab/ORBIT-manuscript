###############################################################################
##  Supplementary Figure 2 -- discreteness bias of the rank transform
##
##  ORBIT assumes r/(n+1) ~ U(0,1), so x = -log(r/(n+1)) ~ Exp(1) and
##  D = sum_j s_j x_j is symmetric variance-gamma.  The rank transform lives on a
##  bounded lattice (|D| <= K log(n+1)), so the continuous reference is wider
##  than the exact null and P values are conservative.
##
##  Panel A: K = 2, n = 5,000 / 10,000 / 20,000, exact by enumeration.
##  Panel B: n = 10,000, K = 2 / 3 / 5; K = 2 exact, K = 3, 5 by permutation.
##  rho = 0 throughout: r = K is an integer, the continuous tail is a terminating
##  closed form, no quadrature anywhere.
##
##  K = 2 bounds:  min continuous P = (log(n+1)+1)/(n+1)^2,
##                 min rank P = 1/(2 n^2),  ratio = 2(log n + 1)
##
##  Deps: ggplot2
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

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

set.seed(20260724)

N_PANEL_A <- c(5000, 10000, 20000)
N_PANEL_B <- 10000
K_PANEL_B <- c(2, 3, 5)
M_PERM    <- 2e8
B_MINP    <- 3e-6
FONT      <- "Arial"

###############################################################################
## 1. continuous reference: D ~ VG_sym(r = K, sigma = 1), rho = 0
##    P(|D|>d) = 2^(1-K) sum_{k=0}^{K-1} c_k Q(K-k, d),  c_k = (K)_k/(k! 2^k)
###############################################################################
vg_tail <- function(d, K) {
  lt <- numeric(K); lc <- 0
  for (k in 0:(K - 1)) {
    lt[k + 1] <- lc + pgamma(d, shape = K - k, lower.tail = FALSE, log.p = TRUE)
    lc <- lc + log(K + k) - log(k + 1) - log(2)
  }
  m <- max(lt)
  exp((1 - K) * log(2) + m + log(sum(exp(lt - m))))
}
## K = 2 must reduce to 0.5*(d+2)*exp(-d)
stopifnot(abs(vg_tail(9, 2) / (0.5 * 11 * exp(-9)) - 1) < 1e-12)

d_ceiling <- function(n, K) K * log(n + 1)

d_for_p <- function(p, K, n) {
  cap <- d_ceiling(n, K)
  if (p <= vg_tail(cap, K)) return(cap * (1 - 1e-12))
  uniroot(function(x) log(vg_tail(x, K)) - log(p), c(1e-8, cap), tol = 1e-12)$root
}

###############################################################################
## 2. exact discrete null, K = 2
##    P(|D|>d) = 0.5*P(x1+x2 > d) + 0.5*P(|x1-x2| > d)
##      x1+x2 > d   <=>  r1*r2 < (n+1)^2 exp(-d)
##      |x1-x2| > d <=>  r2 > r1 exp(d)  or  r2 < r1 exp(-d)
###############################################################################
exact_k2 <- function(d, n) {
  r1  <- seq_len(n)
  tot <- as.numeric(n) * n
  C   <- exp(-d) * (n + 1)^2
  cs  <- sum(pmin(pmax(ceiling(C / r1) - 1, 0), n))
  cd  <- sum(pmax(n - floor(r1 * exp(d)), 0) +
               pmin(pmax(ceiling(r1 * exp(-d)) - 1, 0), n))
  0.5 * cs / tot + 0.5 * cd / tot
}

###############################################################################
## 3. permutation null, any K, one pass over all thresholds
###############################################################################
perm_counts <- function(dvec, n, K, M, batch = 2e6) {
  cnt <- numeric(length(dvec)); left <- M
  while (left > 0) {
    m <- min(batch, left)
    R <- matrix(sample.int(n, m * K, replace = TRUE), ncol = K)
    S <- matrix(sample(c(-1, 1), m * K, replace = TRUE), ncol = K)
    D <- abs(rowSums(S * (-log(R / (n + 1)))))
    cnt <- cnt + vapply(dvec, function(d) sum(D > d), numeric(1))
    left <- left - m
  }
  cnt
}

###############################################################################
## 4. cross-check: enumeration vs permutation, K = 2
###############################################################################
cat("cross-check, exact enumeration vs permutation (n = 10,000, K = 2)\n")
dchk <- c(4, 6, 8, 10, 12)
cchk <- perm_counts(dchk, 10000, 2, M = 2e7)
for (i in seq_along(dchk))
  cat(sprintf("  d=%5.1f  exact %.5e   permutation %.5e   ratio %.4f  (%s events)\n",
              dchk[i], exact_k2(dchk[i], 10000), cchk[i] / 2e7,
              (cchk[i] / 2e7) / exact_k2(dchk[i], 10000), format(cchk[i])))

###############################################################################
## 5. panel A data: K = 2, exact, down to the analytic floor
###############################################################################
cat("\nbuilding panel A (exact enumeration)\n")
datA <- do.call(rbind, lapply(N_PANEL_A, function(n) {
  cap   <- d_ceiling(n, 2)
  pfl   <- vg_tail(cap, 2)                       # smallest continuous P
  dv    <- vapply(10^seq(-2, log10(pfl), length.out = 60),
                  d_for_p, numeric(1), K = 2, n = n)
  data.frame(n = n,
             uniform_P = vapply(dv, vg_tail, numeric(1), K = 2),
             rank_P    = vapply(dv, exact_k2, numeric(1), n = n))
}))
datA$ratio <- datA$uniform_P / datA$rank_P
datA$grp   <- factor(datA$n, levels = N_PANEL_A,
                     labels = format(N_PANEL_A, big.mark = ",", trim = TRUE))

cat("\nanalytic endpoints, K = 2:\n")
for (n in N_PANEL_A) {
  cap <- d_ceiling(n, 2)
  cat(sprintf("  n = %6s | |D|max %.3f | min continuous P %.3e | min rank P %.3e | ratio %.2f (2(log n + 1) = %.2f)\n",
              format(n, big.mark = ","), cap, vg_tail(cap, 2), 1 / (2 * n^2),
              vg_tail(cap, 2) * 2 * n^2, 2 * (log(n) + 1)))
}

###############################################################################
## 6. panel B data: n = 10,000, K = 2 / 3 / 5
###############################################################################
cat(sprintf("\nbuilding panel B (permutation, M = %.0e per K)\n", M_PERM))
datB <- do.call(rbind, lapply(K_PANEL_B, function(K) {
  n   <- N_PANEL_B
  pfl <- max(vg_tail(d_ceiling(n, K), K), B_MINP)
  dv  <- vapply(10^seq(-2, log10(pfl), length.out = 25),
                d_for_p, numeric(1), K = K, n = n)
  up  <- vapply(dv, vg_tail, numeric(1), K = K)
  if (K == 2) {
    rp <- vapply(dv, exact_k2, numeric(1), n = n); ev <- NA_real_
    lo <- rp; hi <- rp
  } else {
    ev <- perm_counts(dv, n, K, M = M_PERM)
    rp <- ev / M_PERM
    lo <- qgamma(0.025, pmax(ev, 1e-9)) / M_PERM
    hi <- qgamma(0.975, ev + 1) / M_PERM
  }
  cat(sprintf("  K = %d done\n", K))
  data.frame(K = K, uniform_P = up, rank_P = rp, events = ev,
             ratio_lo = up / hi, ratio_hi = up / lo)
}))
datB$ratio <- datB$uniform_P / datB$rank_P
datB$grp   <- factor(datB$K, levels = K_PANEL_B,
                     labels = paste("K =", K_PANEL_B))

###############################################################################
## 7. figure
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

PAL_A <- c("#1f78b4", "#e08214", "#4d9221")
PAL_B <- c("#1f78b4", "#8073ac", "#c51b7d")

pA <- ggplot(datA, aes(-log10(uniform_P), ratio, colour = grp, group = grp)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  geom_line(linewidth = 0.9) +
  scale_y_log10(breaks = c(1, 1.5, 2, 3, 5, 10, 20)) +
  scale_x_continuous(breaks = seq(2, 8, 1)) +
  scale_colour_manual(values = PAL_A, name = "features per layer") +
  labs(x = expression(-log[10]~italic(P)[uniform]),
       y = expression(italic(P)[uniform]/italic(P)[rank]),
       title = "K = 2") + th

quartz(file = file.path(OUT_DIR, "rank-VG-1.pdf"), type = "pdf", width = 6, height = 6)
print(pA)
dev.off()

pB <- ggplot(datB, aes(-log10(uniform_P), ratio, colour = grp, group = grp)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45") +
  geom_ribbon(aes(ymin = ratio_lo, ymax = ratio_hi, fill = grp),
              alpha = 0.18, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.9) +
  scale_y_log10(breaks = c(1, 1.5, 2, 3, 5)) +
  scale_x_continuous(breaks = seq(2, 6, 1)) +
  scale_colour_manual(values = PAL_B, name = "omics layers") +
  scale_fill_manual(values = PAL_B) +
  labs(x = expression(-log[10]~italic(P)[uniform]),
       y = expression(italic(P)[uniform]/italic(P)[rank]),
       title = "10,000 features per layer") + th

quartz(file = file.path(OUT_DIR, "rank-VG-2.pdf"), type = "pdf", width = 6, height = 6)
print(pB)
dev.off()