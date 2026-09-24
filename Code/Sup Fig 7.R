# ============================================================================
# Supplementary Figure 7 (DCM proteome meta-analysis) -- plotting
#   Input : Data/DCM/DCM_Prot_final_dat.rds
#             $orbit   proteome meta-analysis (Feature, ORBIT_Rank_P, ...)
#             $omics   list of proteome studies (gene, logFC, adj.P.Val)
#   Output: Data/DCM/SupFig7_DCM_prot_top10_bubble.pdf
# ============================================================================

## ------------------------------------------------------------------------
## Portable setup: locate this script, read/write ../Data/DCM
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
DATA_DIR <- normalizePath(file.path(script_dir, "..", "Data/DCM"), mustWork = FALSE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
message("Script dir: ", getwd(), "  |  data dir -> ", DATA_DIR)

for (pkg in c("dplyr", "purrr", "ggplot2"))
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
})

FONT <- "Arial"
pdf_path <- function(f) file.path(DATA_DIR, f)

# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
DCM_Prot <- readRDS(pdf_path("DCM_Prot_final_dat.rds"))

# ============================================================================
# Sup Fig 7. Proteome meta-analysis: ORBIT top-10 proteins across the proteome studies
# ============================================================================
top_prot   <- DCM_Prot$orbit %>% arrange(ORBIT_Rank_P) %>% slice_head(n = 10)
gene_order <- top_prot$Feature
study_lv   <- names(DCM_Prot$omics)

bubble_S7 <- imap_dfr(DCM_Prot$omics, function(df, study) {
  df %>% filter(gene %in% gene_order) %>%
    transmute(Feature = gene, study = study, logFC, adjP = adj.P.Val)
}) %>%
  mutate(neglog10_adjP = -log10(adjP),
         study   = factor(study, levels = study_lv),
         Feature = factor(Feature, levels = rev(gene_order)))
fc_max <- max(abs(bubble_S7$logFC), na.rm = TRUE)

p_S7 <- ggplot(bubble_S7, aes(x = study, y = Feature)) +
  geom_point(aes(size = neglog10_adjP, fill = logFC),
             shape = 21, colour = "black", stroke = 0.4) +
  scale_fill_gradient2(low = "#4C6FAD", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-fc_max, fc_max),
                       name = expression(log[2]~FC)) +
  scale_size_continuous(range = c(4, 10), breaks = c(0.1, 0.5, 1.0, 5.0),
                        name = expression(-log[10]~adj.P)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = FONT) +
  theme(axis.text.x  = element_text(face = "bold", size = 13, family = FONT, angle = 45, hjust = 1),
        axis.text.y  = element_text(face = "bold", size = 16, colour = "grey30", family = FONT),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.line    = element_blank(),
        legend.title = element_text(family = FONT),
        legend.text  = element_text(family = FONT),
        plot.margin  = margin(10, 20, 10, 10))

quartz(file = pdf_path("SupFig7_DCM_prot_top10_bubble.pdf"), type = "pdf", width = 6, height = 6)
print(p_S7)
dev.off()