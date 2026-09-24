# ============================================================================
# Figure 4 (DCM transcriptome meta-analysis) -- plotting
#   Input : Data/DCM/DCM_Transcriptome.rds
#             $limma_results  list by study: ID, t, P.Value, adj.P.Val
#             $orbit_full     Feature, P, padj, Direction, N
#             $orbit_loo      list by left-out study, same columns as orbit_full
#             $fisher         Feature, P
#             $gsea$Reactome  clusterProfiler gseaResult
#   Outputs: Data/DCM/*.pdf
#   4B signed-rank strips per study + ORBIT
#   4C leave-one-out Spearman heatmap (signed -log10 P)
#   4D leave-one-out significant-gene counts
#   4E Fisher->ORBIT rank-gap top20 bubble
#   4G ECM genes <-> Reactome pathways chord diagram
#   Sup Fig 5 ORA vs GSEA concordance on robust genes
#   Sup Fig 6 LOO heatmap (unsigned -log10 P) + counts for Fisher / Stouffer / DPM
#   ($fisher_loo etc. come from DCM_LOO_methods.R)
# ============================================================================

## ------------------------------------------------------------------------
## Portable setup: locate this script, read/write ../Data/DCM
## Layout: <repo>/Code/<this script>, <repo>/Data/DCM/
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
DATA_DIR <- normalizePath(file.path(script_dir, "..", "Data/DCM"), mustWork = FALSE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
message("Script dir: ", getwd(), "  |  data dir -> ", DATA_DIR)

for (pkg in c("dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork", "circlize", "msigdbr"))
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("ComplexHeatmap", "clusterProfiler"))
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  library(clusterProfiler)
  library(msigdbr)
})

FONT <- "Arial"
pdf_path <- function(f) file.path(DATA_DIR, f)

# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
DCM <- readRDS(pdf_path("DCM_Transcriptome.rds"))
lr   <- DCM$limma_results
orb  <- DCM$orbit_full
loo  <- DCM$orbit_loo
gsea_res <- as.data.frame(DCM$gsea$Reactome@result)

FDR <- 0.05
STUDY_ORDER <- c("Flam_2022", "Sweet_2018", "Hua_2019", "Spurrell_2022",
                 "vanHeesch_2019", "Hannenhalli_2006", "Kittleson_2005",
                 "Yang_2014", "Barth_2006", "Tarazon_2014")

# significance used for LOO panels: padj < 0.05 and supported by > 5 studies
is_sig <- function(d) !is.na(d$padj) & d$padj < FDR & d$N > 5

# ============================================================================
# 4B. Signed-rank score per study and for ORBIT; ORBIT top10 up/down traced
#   score = sign(stat) * (1 - rank(P)/(n+1))
# ============================================================================
TOP_LABEL <- 10
pal_dir   <- c("Up" = "#C0504D", "Down" = "#4A7BA7", "n.s." = "grey80")

signed_rank_score <- function(sign_stat, pval) {
  n <- sum(is.finite(pval))
  r <- rank(pval, ties.method = "average", na.last = "keep")
  sign(sign_stat) * (1 - r / (n + 1))
}

study_long <- bind_rows(lapply(names(lr), function(s)
  data.frame(Feature = lr[[s]]$ID, Study = s,
             score = signed_rank_score(lr[[s]]$t, lr[[s]]$P.Value),
             stringsAsFactors = FALSE)))

# group defined once from the ORBIT meta result, mapped onto every study
orbit_group <- data.frame(
  Feature = orb$Feature,
  group   = case_when(orb$padj < FDR & orb$Direction > 0 ~ "Up",
                      orb$padj < FDR & orb$Direction < 0 ~ "Down",
                      TRUE ~ "n.s."),
  stringsAsFactors = FALSE)

orbit_long <- data.frame(Feature = orb$Feature, Study = "ORBIT Transcriptome",
                         score = signed_rank_score(orb$Direction, orb$P),
                         stringsAsFactors = FALSE)
orbit_scored <- left_join(orbit_long, orbit_group, by = "Feature")

lvls <- c(STUDY_ORDER, "ORBIT Transcriptome")
dat_4B <- bind_rows(study_long, orbit_long) %>%
  left_join(orbit_group, by = "Feature") %>%
  mutate(group = ifelse(is.na(group), "n.s.", group),
         Study = factor(Study, levels = lvls),
         group = factor(group, levels = c("Up", "Down", "n.s."))) %>%
  arrange(group != "n.s.")                     # n.s. drawn first (underneath)

top_up   <- orbit_scored %>% filter(group == "Up") %>%
  slice_max(score, n = TOP_LABEL, with_ties = FALSE) %>% pull(Feature)
top_down <- orbit_scored %>% filter(group == "Down") %>%
  slice_min(score, n = TOP_LABEL, with_ties = FALSE) %>% pull(Feature)
top_genes <- c(top_up, top_down)

line_dat <- dat_4B %>% filter(Feature %in% top_genes) %>%
  mutate(dir = ifelse(Feature %in% top_up, "Up", "Down"))
lab_dat <- orbit_scored %>% filter(Feature %in% top_genes) %>%
  mutate(Study = factor("ORBIT Transcriptome", levels = lvls),
         dir   = ifelse(Feature %in% top_up, "Up", "Down"))
n_all <- length(lvls)

p_4B <- ggplot(dat_4B, aes(Study, score)) +
  geom_hline(yintercept = 0, color = "grey55") +
  geom_hline(yintercept = c(-0.5, 0.5), color = "grey88") +
  geom_vline(xintercept = length(STUDY_ORDER) + 0.5,
             linetype = "dashed", color = "grey55") +
  geom_point(aes(color = group), size = 0.35, alpha = 0.4) +
  geom_line(data = line_dat, aes(group = Feature, color = dir),
            alpha = 0.5, linewidth = 0.3, show.legend = FALSE) +
  geom_point(data = lab_dat, aes(color = dir), size = 1.6, show.legend = FALSE) +
  geom_text_repel(data = lab_dat, aes(label = Feature, color = dir),
                  family = FONT, size = 5,
                  hjust = 0, direction = "y", nudge_x = 0.35,
                  segment.size = 0.2, segment.alpha = 0.4,
                  box.padding = 0.1, max.overlaps = Inf,
                  xlim = c(n_all + 0.7, n_all + 2.8), show.legend = FALSE) +
  scale_color_manual(values = pal_dir, breaks = c("Up", "Down", "n.s."),
                     labels = c("Up regulate (Adjust P < 0.05)",
                                "Down Regulate (Adjust P < 0.05)",
                                "No Significant")) +
  scale_x_discrete(expand = expansion(add = c(0.6, 2))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Signed Rank Score", color = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
  theme_classic(base_size = 12, base_family = FONT) +
  theme(legend.position = "top",
        legend.text  = element_text(size = 11),
        axis.title.y = element_text(size = 14),
        axis.text.y  = element_text(size = 12),
        axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
        plot.margin  = margin(6, 10, 6, 6))

quartz(file = pdf_path("4B_signed_rank.pdf"), type = "pdf", width = 10, height = 6)
print(p_4B)
dev.off()

# ============================================================================
# Leave-one-out helpers (used for ORBIT in 4C/4D and for Fisher/Stouffer/DPM
# in Sup Fig 6). `loo` = list of per-left-out-study results, `full` = all-study
# result; both with columns Feature, N, P, Direction, padj.
# ============================================================================

# Spearman across LOO runs + full; upper triangle, hclust order, full-result
# tick label in black.
#   ORBIT                   : stat = Direction * -log10(P)  (ORBIT reports a direction)
#   Fisher / Stouffer / DPM : stat = -log10(P)              (ActivePathways does not)
loo_heatmap <- function(loo, full, full_label, title = NULL) {
  signed <- identical(full_label, "ORBIT")
  runs <- c(loo, setNames(list(full), full_label))
  long <- bind_rows(lapply(names(runs), function(s) {
    d <- runs[[s]]
    data.frame(Feature = d$Feature, run = s,
               stat = if (signed) d$Direction * -log10(d$P) else -log10(d$P),
               stringsAsFactors = FALSE)
  }))
  mat <- as.matrix(pivot_wider(long, names_from = run, values_from = stat)[, -1])
  M   <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
  ord <- hclust(as.dist(1 - M))$order
  M   <- M[ord, ord]
  lev <- rownames(M)
  M[lower.tri(M)] <- NA
  
  mdf <- as.data.frame(as.table(M)); colnames(mdf) <- c("run1", "run2", "rho")
  mdf <- mdf[!is.na(mdf$rho), ]
  xlev <- lev; ylev <- rev(lev)
  mdf$run1 <- factor(mdf$run1, levels = xlev)
  mdf$run2 <- factor(mdf$run2, levels = ylev)
  lim <- c(min(M, na.rm = TRUE), 1)
  mdf$txt <- ifelse(mdf$rho > mean(lim), "white", "grey20")
  xcol <- ifelse(xlev == full_label, "black", "grey30")
  ycol <- ifelse(ylev == full_label, "black", "grey30")
  
  ggplot(mdf, aes(run1, run2, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", rho), color = txt), size = 2.6, family = FONT) +
    scale_fill_gradientn(colours = c("#F7F0E1", "#D98E73", "#8C2D19"),
                         limits = lim, name = expression(Spearman~rho)) +
    scale_color_identity() +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 12, base_family = FONT) +
    theme(axis.text.x  = element_text(face = "bold", size = 11, family = FONT,
                                      colour = xcol, angle = 45, hjust = 1),
          axis.text.y  = element_text(face = "bold", size = 11, family = FONT, colour = ycol),
          axis.line    = element_blank(),
          panel.grid   = element_blank(),
          plot.title   = element_text(face = "bold", size = 14, hjust = 0),
          legend.title = element_text(family = FONT),
          legend.text  = element_text(family = FONT),
          legend.position      = c(0.80, 0.80),  # empty lower-right triangle
          legend.justification = c(0.5, 0.5),
          legend.background    = element_rect(fill = "white", colour = NA),
          plot.margin  = margin(10, 10, 10, 10))
}

# genes significant in every LOO run with a consistent direction
robust_genes <- function(loo) {
  sig_list <- lapply(loo, function(d) {
    s <- d[is_sig(d), c("Feature", "Direction")]
    s <- s[!duplicated(s$Feature), ]
    setNames(s$Direction, s$Feature)
  })
  shared <- Reduce(intersect, lapply(sig_list, names))
  if (length(shared) == 0) return(list(up = character(0), down = character(0)))
  dir_mat <- vapply(sig_list, function(v) unname(v[shared]), numeric(length(shared)))
  dir_mat <- matrix(dir_mat, nrow = length(shared), dimnames = list(shared, names(sig_list)))
  consistent <- shared[apply(dir_mat, 1, function(x) {
    x <- x[!is.na(x)]; length(x) > 0 && length(unique(x)) == 1
  })]
  gene_dir <- dir_mat[consistent, 1]
  list(up   = names(gene_dir)[gene_dir ==  1],
       down = names(gene_dir)[gene_dir == -1])
}

# significant-gene counts per LOO run (Up / Down), "Shared by all" pinned at
# the bottom, dashed line = full result
loo_counts <- function(loo, full, title = NULL, x_breaks = waiver()) {
  rg <- robust_genes(loo)
  cnt <- bind_rows(lapply(names(loo), function(s) {
    sig <- loo[[s]][is_sig(loo[[s]]), ]
    data.frame(group = s, Up = sum(sig$Direction == 1), Down = sum(sig$Direction == -1))
  }))
  dat <- bind_rows(cnt, data.frame(group = "Shared by all",
                                   Up = length(rg$up), Down = length(rg$down)))
  dat$Total <- dat$Up + dat$Down
  dat$group <- factor(dat$group, levels = c("Shared by all", cnt$group[order(cnt$Up + cnt$Down)]))
  long <- dat %>% pivot_longer(c(Up, Down), names_to = "dir", values_to = "n")
  long$dir <- factor(long$dir, levels = c("Up", "Down"))
  full_total <- sum(is_sig(full))
  
  ggplot(long, aes(x = n, y = group, fill = dir)) +
    geom_col(width = 0.72) +
    geom_text(data = dat, aes(x = Total, y = group, label = Total),
              inherit.aes = FALSE, hjust = -0.2, size = 3.2, family = FONT) +
    geom_vline(xintercept = full_total, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    annotate("text", x = full_total, y = 0.55, label = paste0("Full = ", full_total),
             hjust = 1.05, vjust = 0, size = 3.2, family = FONT, colour = "grey40") +
    scale_fill_manual(values = c(Up = "#B1453F", Down = "#3C6E96"),
                      labels = c(Up = "Up regulate", Down = "Down Regulate"), name = NULL) +
    scale_x_continuous(breaks = x_breaks, expand = expansion(mult = c(0, 0.12))) +
    labs(x = "Number of Significance", y = NULL, title = title) +
    theme_minimal(base_size = 12, base_family = FONT) +
    theme(axis.text.y  = element_text(size = 11, family = FONT, colour = "black"),
          axis.text.x  = element_text(size = 11, family = FONT, colour = "black"),
          axis.title.x = element_text(size = 12, family = FONT),
          panel.grid   = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
          plot.title   = element_text(face = "bold", size = 14, hjust = 0),
          legend.position = "top", legend.justification = "center",
          legend.text  = element_text(family = FONT),
          plot.margin  = margin(10, 20, 10, 10)) +
    coord_cartesian(clip = "off")
}

# ============================================================================
# 4C. ORBIT leave-one-out Spearman heatmap
# ============================================================================
p_4C <- loo_heatmap(loo, orb, full_label = "ORBIT")
quartz(file = pdf_path("4C_loo_spearman_heatmap.pdf"), type = "pdf",
       width = 8, height = 8, family = FONT)
print(p_4C)
dev.off()

# ============================================================================
# 4D. ORBIT leave-one-out significant-gene counts
# ============================================================================
rg_orbit   <- robust_genes(loo)
up_genes   <- rg_orbit$up
down_genes <- rg_orbit$down
cat("robust ORBIT genes:", length(up_genes) + length(down_genes),
    " up:", length(up_genes), " down:", length(down_genes), "\n")

p_4D <- loo_counts(loo, orb, x_breaks = c(0, 500, 1000))
quartz(file = pdf_path("4D_loo_sig_counts.pdf"), type = "pdf",
       width = 7, height = 8, family = FONT)
print(p_4D)
dev.off()

# ============================================================================
# 4E. Largest Fisher -> ORBIT rank drops (Fisher top 260), 20 genes;
#     per-study and ORBIT direction / significance bubble
# ============================================================================
fish <- DCM$fisher
fish$rank_fisher <- rank(fish$P, ties.method = "min")
orb$rank_orbit   <- rank(orb$P, ties.method = "min")
m_gap <- merge(fish[, c("Feature", "rank_fisher")], orb[, c("Feature", "rank_orbit")], by = "Feature")
m_gap <- m_gap[m_gap$rank_fisher <= 260, ]
m_gap$rank_diff <- m_gap$rank_orbit - m_gap$rank_fisher
gap_genes <- head(m_gap$Feature[order(-m_gap$rank_diff)], 20)

plot_gap_bubble <- function(genes, study_order = names(lr), fdr = FDR, cap = 15) {
  sd <- bind_rows(lapply(names(lr), function(s) {
    d <- lr[[s]]; d <- d[d$ID %in% genes, ]
    data.frame(Feature = d$ID, col = s, dir = sign(d$t), adjP = d$adj.P.Val,
               stringsAsFactors = FALSE)
  }))
  o  <- orb[orb$Feature %in% genes, ]
  od <- data.frame(Feature = o$Feature, col = "ORBIT", dir = o$Direction, adjP = o$padj,
                   stringsAsFactors = FALSE)
  dat <- rbind(sd, od) %>%
    mutate(neglogP = pmin(-log10(adjP), cap),
           sig     = !is.na(adjP) & adjP < fdr,
           cat     = paste0(ifelse(dir > 0, "Up", "Down"), ifelse(sig, "_sig", "_ns")),
           col     = factor(col, levels = c(study_order, "ORBIT")),
           Feature = factor(Feature, levels = rev(genes)))
  pal <- c(Up_sig = "#B2182B", Up_ns = "#F4A582", Down_sig = "#2166AC", Down_ns = "#92C5DE")
  
  ggplot(dat, aes(col, Feature, size = neglogP, color = cat)) +
    geom_vline(xintercept = length(study_order) + 0.5, linetype = "dashed", color = "grey60") +
    geom_point() +
    scale_color_manual(values = pal, name = NULL,
                       breaks = c("Up_sig", "Down_sig", "Up_ns", "Down_ns"),
                       labels = c("Significant Up", "Significant Down",
                                  "Not Significant Up", "Not Significant Down")) +
    scale_size_continuous(range = c(3, 7), name = expression(-log[10](adj~P))) +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 12, base_family = FONT) +
    theme(axis.text.x  = element_text(face = "bold", size = 12, family = FONT,
                                      angle = 45, hjust = 1),
          axis.text.y  = element_text(face = "bold", size = 14, colour = "grey30", family = FONT),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.line    = element_blank(),
          panel.grid   = element_blank(),
          legend.title = element_text(family = FONT),
          legend.text  = element_text(family = FONT),
          plot.margin  = margin(10, 20, 10, 10),
          legend.position = "right")
}

p_4E <- plot_gap_bubble(gap_genes)
quartz(file = pdf_path("4E_rankgap_top20_bubble.pdf"), type = "pdf",
       width = 8, height = 8, family = FONT)
print(p_4E)
dev.off()

# ============================================================================
# 4G. Chord diagram: ECM / collagen genes <-> Reactome pathways
#     edge = gene in pathway leading edge; sector colour = rank of -log10 padj
# ============================================================================
genesEcm <- c("UBB","SGCE",
              "LOXL1","EFEMP2","LAMA2","CAPN3","COL5A1","COL3A1","HTRA1","ASPN","LUM",
              "COL16A1","LTBP2","MFAP4","COL21A1","FMOD","JAM3",
              "COL14A1","COMP","COL4A5","LTBP1","PCOLCE2","HAPLN1","LTBP3","COLGALT2","CTSK",
              "SSPN","BMP4","P3H2","COL1A2","BGN","COL1A1","TLL2","COL8A1")

paths_ecm <- c(
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",
  "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION",
  "REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES",
  "REACTOME_CROSSLINKING_OF_COLLAGEN_FIBRILS",
  "REACTOME_COLLAGEN_FORMATION",
  "REACTOME_COLLAGEN_DEGRADATION",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
  "REACTOME_ECM_PROTEOGLYCANS",
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
  "REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS",
  "REACTOME_INFECTION_WITH_ENTEROBACTERIA",
  "REACTOME_ELASTIC_FIBRE_FORMATION",
  "REACTOME_MOLECULES_ASSOCIATED_WITH_ELASTIC_FIBRES")

short_ecm <- c(
  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION           = "ECM Organization",
  REACTOME_COLLAGEN_FORMATION                          = "Collagen Formation",
  REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES = "Collagen Biosynth. & Enzymes",
  REACTOME_COLLAGEN_CHAIN_TRIMERIZATION                = "Collagen Chain Trimerization",
  REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES = "Collagen Fibril Assembly",
  REACTOME_CROSSLINKING_OF_COLLAGEN_FIBRILS            = "Collagen Crosslinking",
  REACTOME_COLLAGEN_DEGRADATION                        = "Collagen Degradation",
  REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX     = "ECM Degradation",
  REACTOME_ECM_PROTEOGLYCANS                           = "ECM Proteoglycans",
  REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS          = "Integrin Cell-Surface",
  REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS      = "Non-integrin ECM",
  REACTOME_ELASTIC_FIBRE_FORMATION                     = "Elastic Fibre Formation",
  REACTOME_MOLECULES_ASSOCIATED_WITH_ELASTIC_FIBRES    = "Elastic-fibre Molecules",
  REACTOME_INFECTION_WITH_ENTEROBACTERIA               = "Infection with Enterobacteria")

res_ecm <- gsea_res[match(paths_ecm, gsea_res$ID), ]
core    <- setNames(strsplit(res_ecm$core_enrichment, "/"), paths_ecm)

edges <- do.call(rbind, lapply(paths_ecm, function(p) {
  hit <- intersect(genesEcm, core[[p]])
  if (length(hit) == 0) return(NULL)
  data.frame(gene = hit, pathway = short_ecm[[p]], value = 1, stringsAsFactors = FALSE)
}))
hits <- sapply(paths_ecm, function(p) length(intersect(genesEcm, core[[p]])))
print(data.frame(pathway = short_ecm[paths_ecm], n_hits = hits, row.names = NULL), row.names = FALSE)
paths_used <- paths_ecm[hits > 0]
genesEcm   <- genesEcm[genesEcm %in% edges$gene]      # drop genes with no edge

red_ramp <- colorRampPalette(c("#FEE0D2", "#FCBBA1", "#FC9272",
                               "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15"))
sig_to_col <- function(padj) {
  x <- -log10(padj); x[!is.finite(x)] <- max(x[is.finite(x)])
  red_ramp(100)[cut(rank(x, ties.method = "average"), 100, labels = FALSE)]
}
gene_padj <- setNames(orb$padj[match(genesEcm, orb$Feature)], genesEcm)
path_padj <- setNames(res_ecm$p.adjust[match(paths_used, paths_ecm)], short_ecm[paths_used])
grid.col  <- setNames(c(sig_to_col(gene_padj), sig_to_col(path_padj)),
                      c(genesEcm, short_ecm[paths_used]))
order_sectors <- c(genesEcm, short_ecm[paths_used])
zidx   <- setNames(seq_along(order_sectors), order_sectors)[edges$pathway]
n_gene <- length(genesEcm); n_path <- length(paths_used)

quartz(file = pdf_path("4G_ecm_pathways_chord.pdf"), type = "pdf",
       width = 8, height = 8, family = FONT)
circos.clear()
circos.par(start.degree = 100, clock.wise = FALSE,
           canvas.xlim = c(-0.95, 1.2), canvas.ylim = c(-0.92, 0.92),
           gap.after = c(rep(1, n_gene - 1), 30, rep(1, n_path - 1), 30))
chordDiagram(edges, order = order_sectors, grid.col = grid.col,
             link.lwd = 0.4, link.border = "grey85", link.sort = TRUE,
             link.decreasing = FALSE, link.zindex = zidx,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.20))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  s  <- get.cell.meta.data("sector.index")
  xl <- get.cell.meta.data("xlim")
  circos.text(mean(xl), get.cell.meta.data("ylim")[1] + 0.02, s,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
              cex = 0.9, family = FONT)
}, bg.border = NA)
circos.clear()

gene_x <- -log10(gene_padj); gene_x[!is.finite(gene_x)] <- max(gene_x[is.finite(gene_x)])
path_x <- -log10(path_padj); path_x[!is.finite(path_x)] <- max(path_x[is.finite(path_x)])
lg_gene <- Legend(col_fun = colorRamp2(seq(min(gene_x), max(gene_x), length = 5), red_ramp(5)),
                  title = "Gene  -log10(padj)", direction = "horizontal")
lg_path <- Legend(col_fun = colorRamp2(seq(min(path_x), max(path_x), length = 5), red_ramp(5)),
                  title = "Pathway  -log10(padj)", direction = "horizontal")
draw(packLegend(lg_gene, lg_path, direction = "vertical"),
     x = unit(1, "npc") - unit(4, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
dev.off()

# ============================================================================
# Sup Fig 6. ORA on robust up / down genes vs ORBIT GSEA (Reactome);
#            pathways significant in both
# ============================================================================
msig <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
  error = function(e)
    msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))
t2g <- msig[, c("gs_name", "gene_symbol")]
bg  <- intersect(unique(orb$Feature), unique(msig$gene_symbol))
cat("background genes:", length(bg), "\n")

run_ora <- function(genes)
  as.data.frame(enricher(gene = genes, TERM2GENE = t2g, universe = bg,
                         pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500))
up_df   <- run_ora(up_genes);   if (nrow(up_df))   up_df$ora_dir   <- "Up"
down_df <- run_ora(down_genes); if (nrow(down_df)) down_df$ora_dir <- "Down"
cat("ORA sig  up:", nrow(up_df), " down:", nrow(down_df), "\n")
ora_all <- bind_rows(up_df, down_df)

gsea_sig <- gsea_res[gsea_res$p.adjust < 0.05, ]
common   <- intersect(ora_all$ID, gsea_sig$ID)
cat("intersect (ORA sig / GSEA sig):", length(common), "\n")
o <- ora_all [match(common, ora_all$ID), ]
g <- gsea_sig[match(common, gsea_sig$ID), ]

df_S5 <- data.frame(
  ID        = common,
  ora_logP  = -log10(o$p.adjust),
  gsea_logP = -log10(g$p.adjust),
  count     = o$Count,
  NES       = g$NES,
  ora_dir   = o$ora_dir,
  stringsAsFactors = FALSE)
df_S5$dir <- ifelse(df_S5$NES > 0, "Up", "Down")

short_S5 <- c(
  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION                           = "ECM Organization",
  REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES                 = "Collagen Biosynthesis",
  REACTOME_COLLAGEN_FORMATION                                          = "Collagen Formation",
  REACTOME_ECM_PROTEOGLYCANS                                           = "ECM Proteoglycans",
  REACTOME_COLLAGEN_CHAIN_TRIMERIZATION                                = "Collagen Trimerization",
  REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES = "Collagen Fibril Assembly",
  REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM       = "GAG Metabolism Diseases",
  REACTOME_COLLAGEN_DEGRADATION                                        = "Collagen Degradation",
  REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS                          = "Integrin Interactions",
  REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX                     = "ECM Degradation")
df_S5$label <- unname(short_S5[df_S5$ID])
miss <- is.na(df_S5$label)
df_S5$label[miss] <- gsub("_", " ", sub("^REACTOME_", "", df_S5$ID[miss]))

df_S5 <- df_S5[order(-df_S5$gsea_logP), ]
df_S5$show_lab <- seq_len(nrow(df_S5)) <= 10
lo <- -log10(0.05)

p_S5 <- ggplot(df_S5, aes(x = ora_logP, y = gsea_logP)) +
  geom_hline(yintercept = lo, linetype = "dotted", colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = lo, linetype = "dotted", colour = "grey80", linewidth = 0.3) +
  geom_point(aes(size = count, fill = dir), shape = 21, colour = "black", stroke = 0.5, alpha = 0.9) +
  geom_text_repel(data = subset(df_S5, show_lab), aes(label = label),
                  size = 5, family = FONT, colour = "grey20",
                  nudge_x = max(df_S5$ora_logP) * 0.15, max.overlaps = 20,
                  segment.size = 0.2, segment.colour = "grey40",
                  box.padding = 0.4, min.segment.length = 0) +
  scale_fill_manual(values = c(Up = "#B1453F", Down = "#3C6E96"), name = "Direction") +
  scale_size_continuous(range = c(2, 10), name = "ORA gene count") +
  coord_cartesian(xlim = c(lo * 0.9, max(df_S5$ora_logP)  * 1.05),
                  ylim = c(lo * 0.9, max(df_S5$gsea_logP) * 1.05)) +
  labs(x = expression(ORA~-log[10]~(P[adj])),
       y = expression(GSEA~-log[10]~(P[adj]))) +
  theme_minimal(base_size = 12, base_family = FONT) +
  theme(panel.grid   = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
        axis.text    = element_text(colour = "black", family = FONT, size = 14),
        axis.title   = element_text(family = FONT, size = 15),
        legend.title = element_text(family = FONT, size = 14),
        legend.text  = element_text(family = FONT, size = 14),
        legend.position = "right",
        plot.margin  = margin(12, 12, 12, 12))

quartz(file = pdf_path("SupFig5_ora_gsea_concordance.pdf"), type = "pdf",
       width = 8, height = 10, family = FONT)
print(p_S5)
dev.off()

# ============================================================================
# Sup Fig /. Leave-one-out robustness of Fisher / Stouffer / DPM, drawn as in
#            4C (left) and 4D (right); one row per method
#   requires fisher_loo / stouffer_loo / dpm_loo (+ *_full) in the DCM RDS,
#   produced by DCM_LOO_methods.R
#   Fisher's method is non-directional and ActivePathways (Stouffer / DPM)
#   returns no direction, so bars show total significant counts without an
#   up/down split.
# ============================================================================
stopifnot(all(c("fisher_loo", "fisher_full", "stouffer_loo", "stouffer_full",
                "dpm_loo", "dpm_full") %in% names(DCM)))

## --- helpers ---------------------------------------------------------------
# collapse Direction to a single value so loo_counts() draws one bar per run;
# 1L (not 0L) in case loo_counts subsets by Direction > 0 / < 0
no_dir <- function(x) {
  x$full$Direction <- 1L
  x$loo <- lapply(x$loo, function(d) { d$Direction <- 1L; d })
  x
}
# neutral fill + no legend, layered on top of loo_counts()'s own scale
neutral_fill <- list(
  scale_fill_manual(values = c("grey55", "grey55"), guide = "none"),
  theme(legend.position = "none"))

## --- methods (P-value combination, non-directional display) ----------------
loo_pval_methods <- list(
  Fisher   = no_dir(list(loo = DCM$fisher_loo,   full = DCM$fisher_full)),
  Stouffer = no_dir(list(loo = DCM$stouffer_loo, full = DCM$stouffer_full)),
  DPM      = no_dir(list(loo = DCM$dpm_loo,      full = DCM$dpm_full)))

## --- panels ----------------------------------------------------------------
panels_S6 <- unlist(lapply(names(loo_pval_methods), function(m) {
  x <- loo_pval_methods[[m]]
  list(loo_heatmap(x$loo, x$full, full_label = m, title = m),
       loo_counts (x$loo, x$full, title = m) + neutral_fill)
}), recursive = FALSE)

p_S6 <- wrap_plots(panels_S6, ncol = 2, byrow = TRUE, widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", family = FONT))

## --- write -----------------------------------------------------------------
quartz(file = pdf_path("SupFig6_loo_pvalue_methods.pdf"), type = "pdf",
       width = 15, height = 22, family = FONT)
print(p_S6)
dev.off()

# ============================================================================
# Sup Fig /. Leave-one-out robustness of RankProd, drawn as in 4C / 4D
#   requires rankprod_loo / rankprod_full in the DCM RDS (DCM_LOO_RankProd.R,
#   RankProducts(..., na.rm = TRUE), i.e. gene-wise median imputation of
#   missing studies). RankProd returns separate up- and down-regulated
#   P values, so the bars keep the up/down split.
# ============================================================================
stopifnot(all(c("rankprod_loo", "rankprod_full") %in% names(DCM)))

## --- panels ----------------------------------------------------------------
p_S7 <- wrap_plots(
  loo_heatmap(DCM$rankprod_loo, DCM$rankprod_full,
              full_label = "RankProd", title = "RankProd"),
  loo_counts (DCM$rankprod_loo, DCM$rankprod_full, title = "RankProd"),
  ncol = 2, widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold", family = FONT))

## --- write -----------------------------------------------------------------
quartz(file = pdf_path("SupFig7_loo_rankprod.pdf"), type = "pdf",
       width = 15, height = 22 / 3, family = FONT)
print(p_S7)
dev.off()
