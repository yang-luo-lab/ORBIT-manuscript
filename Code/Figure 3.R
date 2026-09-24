# ============================================================================
# Figure 3 (TI) -- plotting
#   Inputs : Data/CKD/TI_orbit.rds (master_table), Data/CKD/TI_gsea.rds (Reactome GSEA)
#   Outputs: Data/CKD/*.pdf
#   3A BioRender schematic (not generated here)
#   3B gene signed-rank lines | 3C ORBIT top10 bubble | 3D 3-set UpSet
#   3E pathway scatter | 3F TCA leading-edge bubble
#   3I 4-set UpSet (incl. Fisher) | 3J MECP2 leading-edge bubble
#   3K Spearman rho of ORBIT / DPM / Stouffer P vs single-omic P
#   Sup Fig 4 top10 bubbles for Fisher (A), DPM (B), Stouffer (C)
# ============================================================================

## ------------------------------------------------------------------------
## Portable setup: locate this script, read/write ../Data/CKD
## Layout: <repo>/Code/<this script>, <repo>/Data/CKD/
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
DATA_DIR <- normalizePath(file.path(script_dir, "..", "Data/CKD"), mustWork = FALSE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
message("Script dir: ", getwd(), "  |  data dir -> ", DATA_DIR)

for (pkg in c("dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork", "scales", "statebins"))
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

FONT <- "Arial"
pdf_path <- function(f) file.path(DATA_DIR, f)

# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
TI_orbit <- readRDS(pdf_path("TI_orbit.rds"))
stopifnot("master_table" %in% names(TI_orbit))
master_table <- TI_orbit$master_table
TI_gsea <- readRDS(pdf_path("TI_gsea.rds"))

alpha_cut <- 0.05
DB        <- "Reactome"
gsea_key  <- c(Transcriptome = "RNA", Proteome = "Protein",
               ORBIT = "ORBIT", Fisher = "Fisher")
sig_ids   <- function(g) g@result %>% filter(p.adjust < alpha_cut) %>% pull(ID)


# ============================================================================
# 3A. Schematic made in BioRender -- not generated here
# ============================================================================

# ============================================================================
# 3B. Cross-omics normalised signed rank (C3 / TIMP3 / MYH7)
#   nsr = sign * (1 - rank/(n+1)), rank 1 = strongest
# ============================================================================
omics_levels  <- c("Transcriptome", "Proteome", "ORBIT")
genes_to_plot <- c("C3", "TIMP3", "MYH7")
gene_cols <- c(C3    = "#C0392B",
               TIMP3 = "#2C6FAD",
               MYH7  = "#3E9B5F")

prepare_plot_df <- function(mt, omics_levels) {
  rna_layer <- mt %>%
    filter(!is.na(RNA_t), !is.na(RNA_logFC)) %>%
    transmute(Feature, omics = "Transcriptome",
              signal = abs(RNA_t), sgn = sign(RNA_logFC))
  prot_layer <- mt %>%
    filter(!is.na(Prot_P), !is.na(Prot_logFC)) %>%
    transmute(Feature, omics = "Proteome",
              signal = -log10(Prot_P), sgn = sign(Prot_logFC))
  orbit_layer <- mt %>%
    filter(!is.na(P), !is.na(Direction)) %>%
    transmute(Feature, omics = "ORBIT",
              signal = -log10(P), sgn = Direction)
  
  bind_rows(rna_layer, prot_layer, orbit_layer) %>%
    group_by(omics) %>%
    mutate(n    = n(),
           rank = rank(-signal, ties.method = "min"),
           nsr  = sgn * (1 - rank / (n + 1))) %>%
    ungroup() %>%
    mutate(omics = factor(omics, levels = omics_levels))
}

all_df  <- prepare_plot_df(master_table, omics_levels)
plot_df <- all_df %>%
  filter(Feature %in% genes_to_plot) %>%
  mutate(Feature = factor(Feature, levels = genes_to_plot))

p_3B <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(data = all_df, aes(x = omics, y = nsr),
             colour = "grey70", alpha = 0.4, size = 1) +
  geom_line(data = plot_df, aes(x = omics, y = nsr,
                                group = Feature, colour = Feature),
            linewidth = 1) +
  geom_point(data = plot_df, aes(x = omics, y = nsr, colour = Feature),
             size = 2.5, shape = 21, fill = "white", stroke = 1.4) +
  ggrepel::geom_text_repel(
    data = plot_df,
    aes(x = omics, y = nsr, label = Feature, colour = Feature),
    family = FONT,
    size = 6.5, fontface = "bold", show.legend = FALSE,
    bg.color = "white", bg.r = 0.15,
    box.padding = 0.6, point.padding = 0.4,
    min.segment.length = 0, segment.size = 0.4, segment.colour = "grey60",
    max.overlaps = Inf, seed = 42
  ) +
  scale_colour_manual(values = gene_cols) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("-1\n(Down)", "-0.5", "0", "0.5", "1\n(Up)")) +
  labs(x = NULL, y = "Normalised Signed Rank", colour = NULL) +
  theme_classic(base_size = 14, base_family = FONT) +
  theme(axis.text.x  = element_text(face = "bold", size = 16, family = FONT),
        axis.text.y  = element_text(family = FONT, size = 14),
        axis.title.y = element_text(size = 16, family = FONT),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(10, 30, 10, 10))

quartz(file = pdf_path("3B_gene_plot.pdf"), type = "pdf", width = 6, height = 6)
print(p_3B)
dev.off()

# ============================================================================
# 3C. ORBIT top-10 genes (N == 2) -- two-omics bubble plot
# ============================================================================
prepare_bubble_df <- function(mt) {
  mt_ranked <- mt %>% arrange(P) %>% mutate(Rank = row_number())
  top_genes <- mt_ranked %>% filter(N == 2) %>% slice_head(n = 10)
  gene_order <- top_genes$Feature
  
  top_genes %>%
    transmute(Feature, Rank,
              Transcriptome_logFC = RNA_logFC,
              Proteome_logFC      = Prot_logFC,
              Transcriptome_adjP  = RNA_adjP,
              Proteome_adjP       = Prot_adjP) %>%
    pivot_longer(cols = -c(Feature, Rank),
                 names_to = c("omics", ".value"), names_sep = "_") %>%
    mutate(neglog10_adjP = -log10(adjP),
           omics   = factor(omics,   levels = c("Transcriptome", "Proteome")),
           Feature = factor(Feature, levels = rev(gene_order)))
}

bubble_df <- prepare_bubble_df(master_table)
fc_max <- max(abs(bubble_df$logFC), na.rm = TRUE)
fc_lim <- c(-fc_max, fc_max)

rank_labels <- bubble_df %>%
  distinct(Feature, Rank) %>%
  arrange(Feature) %>%
  with(setNames(as.character(Rank), as.character(Feature)))

p_3C <- ggplot(bubble_df, aes(x = omics, y = Feature)) +
  geom_point(aes(size = neglog10_adjP, fill = logFC),
             shape = 21, colour = "black", stroke = 0.4) +
  scale_y_discrete(sec.axis = dup_axis(name = NULL,
                                       labels = function(x) rank_labels[x])) +
  scale_fill_gradient2(low = "#4C6FAD", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = fc_lim,
                       name = expression(log[2]~FC)) +
  scale_size_continuous(range = c(4, 10),
                        breaks = c(0.1, 0.5, 1.0, 5.0),
                        name = expression(-log[10]~adj.P)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14, base_family = FONT) +
  theme(axis.text.x       = element_text(face = "bold", size = 15, family = FONT),
        axis.text.y.left  = element_text(face = "bold", size = 16, colour = "grey30", family = FONT),
        axis.text.y.right = element_text(size = 15, colour = "grey30", family = FONT),
        axis.ticks.y.right = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.line    = element_blank(),
        legend.title = element_text(family = FONT),
        legend.text  = element_text(family = FONT),
        plot.margin  = margin(10, 20, 10, 10))

quartz(file = pdf_path("3C_Tub_ORBIT_top10_bubble.pdf"), type = "pdf", width = 6, height = 6)
print(p_3C)
dev.off()

# ============================================================================
# Shared helpers for the hand-built UpSet plots (3D, 3I)
# ============================================================================
build_upset <- function(set_names) {
  sig_list  <- lapply(set_names, function(s) sig_ids(TI_gsea[[ gsea_key[s] ]][[DB]]))
  names(sig_list) <- set_names
  all_paths <- sort(unique(unlist(sig_list)))
  upset_df  <- data.frame(pathway = all_paths, stringsAsFactors = FALSE)
  for (s in set_names) upset_df[[s]] <- upset_df$pathway %in% sig_list[[s]]
  message("Set totals: ",
          paste(set_names, sapply(set_names, function(s) sum(upset_df[[s]])),
                sep = "=", collapse = ", "))
  
  combo <- upset_df %>%
    mutate(combo = apply(across(all_of(set_names)), 1,
                         function(x) paste(as.integer(x), collapse = ""))) %>%
    count(combo, name = "size") %>%
    filter(combo != paste(rep(0, length(set_names)), collapse = "")) %>%
    arrange(desc(size)) %>%
    mutate(combo_id = factor(row_number(), levels = as.character(row_number())))
  
  combo_long <- combo %>%
    select(combo_id, combo) %>%
    tidyr::separate(combo, into = set_names, sep = "(?<=.)(?=.)", remove = FALSE) %>%
    pivot_longer(all_of(set_names), names_to = "set", values_to = "member") %>%
    mutate(set = factor(set, levels = rev(set_names)), member = member == "1")
  
  set_tot <- tibble(set  = factor(set_names, levels = rev(set_names)),
                    size = sapply(set_names, function(s) sum(upset_df[[s]])))
  
  list(upset_df = upset_df, combo = combo, combo_long = combo_long, set_tot = set_tot)
}

base_theme <- theme_classic(base_size = 14, base_family = FONT)
p_empty    <- ggplot() + theme_void()

# ============================================================================
# 3D. Reactome significant-pathway UpSet -- 3 methods, bars coloured by theme
# ============================================================================
set_names <- c("Transcriptome", "Proteome", "ORBIT")
set_cols  <- c(Transcriptome = "#3E9B5F",
               Proteome      = "#2C6FAD",
               ORBIT         = "#C0392B")
U <- build_upset(set_names)
combo <- U$combo; combo_long <- U$combo_long; set_tot <- U$set_tot

x_scale <- scale_x_discrete(limits = levels(combo$combo_id),
                            expand = expansion(add = 0.6))

# theme composition per intersection; counts taken from Cytoscape clusters (Fig. 3G/H)
theme_names <- c("Cell cycle & transcription",
                 "Inflammation & immune",
                 "ECM regulation",
                 "Mitochondria")
cat_levels  <- c(theme_names, "Other")
cat_cols <- c("Cell cycle & transcription" = "#9BDCEF",
              "Inflammation & immune"      = "#A99CF0",
              "ECM regulation"             = "#C77DC0",
              "Mitochondria"               = "#5B8FB0",
              "Other"                      = "grey60")

theme_counts <- tribble(
  ~Transcriptome, ~Proteome, ~ORBIT, ~`Cell cycle & transcription`, ~`Inflammation & immune`, ~`ECM regulation`, ~Mitochondria,
  TRUE,  FALSE, FALSE,  15, 12,  2,  1,   # transcriptome only
  FALSE, FALSE, TRUE,    4, 10,  4,  7,   # ORBIT only
  FALSE, TRUE,  FALSE,   4,  8,  1,  5,   # proteome only
  TRUE,  FALSE, TRUE,   74, 23,  4,  8,   # transcriptome + ORBIT
  FALSE, TRUE,  TRUE,    1, 13, 11,  3,   # proteome + ORBIT
  TRUE,  TRUE,  TRUE,    4, 20, 15, 12    # all three
) %>%
  mutate(combo = paste0(as.integer(Transcriptome),
                        as.integer(Proteome),
                        as.integer(ORBIT))) %>%
  select(-Transcriptome, -Proteome, -ORBIT)

combo_cat <- combo %>%
  select(combo, size, combo_id) %>%
  left_join(theme_counts, by = "combo") %>%
  mutate(Other = size - (`Cell cycle & transcription` + `Inflammation & immune` +
                           `ECM regulation` + Mitochondria)) %>%
  select(-size) %>%
  pivot_longer(all_of(cat_levels), names_to = "category", values_to = "count") %>%
  mutate(category = factor(category, levels = cat_levels))

p_bars <- ggplot(combo_cat, aes(combo_id, count, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(data = combo, aes(combo_id, size, label = size),
            inherit.aes = FALSE, vjust = -0.4, size = 5, family = FONT) +
  scale_fill_manual(values = cat_cols, name = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  x_scale +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank(),
        plot.margin  = margin(5, 5, 0, 0),
        legend.position        = "inside",
        legend.position.inside = c(0.99, 0.99),
        legend.justification   = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.6), colour = NA),
        legend.margin     = margin(2, 4, 2, 4),
        legend.key.size   = unit(1.3, "lines"),
        legend.spacing.y  = unit(1, "pt"),
        legend.text       = element_text(size = 13, family = FONT))

p_dots <- ggplot(combo_long, aes(combo_id, set)) +
  geom_point(aes(colour = member), size = 4) +
  geom_line(data = filter(combo_long, member),
            aes(group = combo_id), linewidth = 1, colour = "grey20") +
  scale_colour_manual(values = c(`TRUE` = "grey20", `FALSE` = "grey85"),
                      guide = "none") +
  x_scale +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line  = element_blank(),
        panel.background = element_rect(fill = "grey96", colour = NA),
        plot.margin = margin(0, 5, 5, 0))

p_title <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Pathways in intersection",
           angle = 90, family = FONT, size = 6) +
  theme_void()

p_names <- ggplot(set_tot, aes(x = 1, y = set)) +
  geom_text(aes(label = set), family = FONT, fontface = "bold",
            size = 5, hjust = 0.5) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  theme_void()

p_left <- ggplot(set_tot, aes(size, set, fill = as.character(set))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = size), hjust = 1.1, size = 4, family = FONT) +
  scale_fill_manual(values = set_cols) +
  guides(fill = "none") +
  scale_x_reverse(expand = expansion(mult = c(0.2, 0))) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  labs(x = "Total pathways", y = NULL) +
  base_theme +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        plot.margin  = margin(1, 2, 5, 5))

p_3D <- p_empty + p_title + p_bars +
  p_left  + p_names + p_dots +
  plot_layout(ncol = 3, nrow = 2,
              widths  = c(1.25, 1, 2.4),
              heights = c(3.1, 1))

quartz(file = pdf_path("3D_reactome_upset.pdf"), type = "pdf", width = 7, height = 7)
print(p_3D)
dev.off()

# ============================================================================
# 3E. Reactome pathway scatter: Transcriptome vs Proteome significance
# ============================================================================
grab <- function(g, tag) {
  g@result %>%
    transmute(ID,
              !!paste0(tag, "_sig")  := sign(NES) * -log10(p.adjust),
              !!paste0(tag, "_padj") := p.adjust)
}
rna   <- grab(TI_gsea$RNA[[DB]],     "RNA")
prot  <- grab(TI_gsea$Protein[[DB]], "Prot")
orbit <- grab(TI_gsea$ORBIT[[DB]],   "ORBIT")

scatter_df <- inner_join(rna, prot, by = "ID") %>%
  inner_join(orbit, by = "ID") %>%
  mutate(class = case_when(                         # five classes by priority
    ORBIT_padj < alpha_cut &
      RNA_padj  >= alpha_cut &
      Prot_padj >= alpha_cut          ~ "ORBIT Only",
    ORBIT_padj < alpha_cut            ~ "ORBIT with Others",
    RNA_padj   < alpha_cut            ~ "Transcriptome Only",
    Prot_padj  < alpha_cut            ~ "Proteome Only",
    TRUE                              ~ "Not Significant"
  )) %>%
  mutate(class = factor(class,
                        levels = c("ORBIT Only", "ORBIT with Others",
                                   "Transcriptome Only", "Proteome Only",
                                   "Not Significant")))

# labels: ORBIT-significant pathways, top 5 up + top 5 down by ORBIT_sig
orbit_sig_classes <- c("ORBIT Only", "ORBIT with Others")
lab_up   <- scatter_df %>% filter(class %in% orbit_sig_classes) %>% slice_max(ORBIT_sig, n = 5)
lab_down <- scatter_df %>% filter(class %in% orbit_sig_classes) %>% slice_min(ORBIT_sig, n = 5)

clean_label <- function(id) {
  x <- sub("^REACTOME_", "", id)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  x <- gsub("\\b(\\w)", "\\U\\1", x, perl = TRUE)
  abbr <- c("Slc"="SLC", "Ecm"="ECM", "Mhc"="MHC", "Tca"="TCA",
            "Dna"="DNA", "Rna"="RNA", "Atp"="ATP", "Gpcr"="GPCR",
            "Tcr"="TCR", "Bcr"="BCR", "Il"="IL", "Tnf"="TNF",
            "Nadh"="NADH", "Fgfr"="FGFR", "Egfr"="EGFR")
  for (a in names(abbr))
    x <- gsub(paste0("\\b", a, "\\b"), abbr[[a]], x)
  shorten <- c(
    "Immunoregulatory Interactions Between A Lymphoid And A Non Lymphoid Cell"
    = "Lymphoid\u2013Non-Lymphoid Interactions",
    "Aerobic Respiration And Respiratory Electron Transport"
    = "Aerobic Respiration & Electron Transport",
    "Metabolism Of Amino Acids And Derivatives"
    = "Amino Acid Metabolism"
  )
  if (x %in% names(shorten)) x <- shorten[[x]]
  x
}
lab_df <- bind_rows(lab_up, lab_down) %>%
  mutate(label = vapply(ID, clean_label, character(1)))

message("Classes: ",
        paste(names(table(scatter_df$class)), table(scatter_df$class),
              sep = "=", collapse = ", "))

class_cols <- c("ORBIT Only"         = "#C0392B",
                "ORBIT with Others"  = "#F1948A",
                "Transcriptome Only" = "#3E9B5F",
                "Proteome Only"      = "#2C6FAD",
                "Not Significant"    = "grey75")

p_3E <- ggplot(scatter_df, aes(RNA_sig, Prot_sig, colour = class)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(data = ~ filter(.x, class == "Not Significant"),
             size = 1.6, alpha = 0.6) +
  geom_point(data = ~ filter(.x, class != "Not Significant"),
             size = 2.2) +
  ggrepel::geom_text_repel(
    data = lab_df, aes(label = label),
    family = FONT, size = 4, colour = "black", fontface = "bold",
    bg.color = "white", bg.r = 0.12,
    box.padding = 0.5, point.padding = 0.3,
    min.segment.length = 0, segment.size = 0.3, segment.colour = "grey50",
    max.overlaps = Inf, seed = 42
  ) +
  scale_colour_manual(values = class_cols, name = NULL) +
  labs(x = "Transcriptome Pathways Significance",
       y = "Proteome Pathways Significance") +
  theme_classic(base_size = 14, base_family = FONT) +
  theme(legend.position    = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.direction   = "vertical",
        legend.background  = element_rect(fill = alpha("white", 0.7),
                                          colour = "grey70", linewidth = 0.3),
        legend.text  = element_text(size = 9, family = FONT),
        axis.title   = element_text(size = 15, family = FONT),
        axis.text    = element_text(size = 12, family = FONT),
        plot.margin  = margin(10, 15, 10, 10)) +
  guides(colour = guide_legend(override.aes = list(size = 3.5)))

quartz(file = pdf_path("3E_reactome_scatter.pdf"), type = "pdf", width = 6, height = 6)
print(p_3E)
dev.off()

# ============================================================================
# Shared helpers for leading-edge bubble plots (3F, 3J)
# ============================================================================
rea_rna    <- TI_gsea$RNA$Reactome@result
rea_prot   <- TI_gsea$Protein$Reactome@result
rea_orbit  <- TI_gsea$ORBIT$Reactome@result
rea_fisher <- TI_gsea$Fisher$Reactome@result

get_leading <- function(gsea_df, pid) {
  row <- gsea_df %>% filter(ID == pid)
  if (nrow(row) == 0) return(character(0))
  ce <- row$core_enrichment[1]
  if (is.na(ce) || ce == "") return(character(0))
  strsplit(ce, "/")[[1]]
}
get_padj <- function(gsea_df, pid) {
  v <- gsea_df %>% filter(ID == pid) %>% pull(p.adjust)
  if (length(v) == 0) NA else v[1]
}
fmt_p <- function(p) if (is.na(p)) "NA" else sprintf("%.3f", p)

# contiguous runs of y positions -> rectangle start/end
find_runs <- function(pos) {
  if (length(pos) == 0) return(data.frame(start = integer(), end = integer()))
  pos <- sort(pos)
  if (length(pos) == 1) return(data.frame(start = pos, end = pos))
  brk <- which(diff(pos) > 1)
  data.frame(start = pos[c(1, brk + 1)], end = pos[c(brk, length(pos))])
}
mk_runs <- function(le, gene_ypos, xmin, xmax, core) {
  pos <- unname(gene_ypos[names(gene_ypos) %in% le])
  find_runs(pos) %>% mutate(xmin = xmin, xmax = xmax, core = core)
}

# gene table (logFC / adjP from master_table) -> long format for bubbles
make_plot_long <- function(gene_table) {
  gene_table %>%
    dplyr::select(Gene, RNA_logFC, RNA_padj, Prot_logFC, Prot_padj) %>%
    pivot_longer(cols = -Gene, names_to = c("Omic", ".value"),
                 names_pattern = "(.+)_(.+)") %>%
    filter(!is.na(logFC) & !is.na(padj)) %>%
    mutate(neg_log_padj = -log10(pmax(padj, 1e-300)),
           Omic = factor(Omic, levels = c("RNA", "Prot"),
                         labels = c("RNA", "Protein")))
}

base_bubble <- function(plot_long, title, y_size = 15) {
  max_nlp     <- max(plot_long$neg_log_padj, na.rm = TRUE)
  size_breaks <- pretty(c(0, max_nlp), n = 4)
  size_breaks <- size_breaks[size_breaks <= max_nlp]
  ggplot(plot_long, aes(x = Omic, y = Gene)) +
    geom_point(aes(fill = logFC, size = neg_log_padj),
               shape = 21, color = "black", stroke = 0.6, alpha = 0.9) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = expression(log[2]~FC),
                         guide = guide_colorbar(order = 3)) +
    scale_size_continuous(range = c(3, 8), name = expression(-log[10]~adj.P),
                          breaks = size_breaks, limits = c(0, max_nlp),
                          guide = guide_legend(order = 2)) +
    scale_x_discrete(limits = c("RNA", "Protein"),
                     labels = c("RNA" = "Transcriptome", "Protein" = "Proteome")) +
    theme_bw(base_size = 13, base_family = FONT) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y = element_text(face = "bold", size = y_size),
          axis.text.x = element_text(size = 18, face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
          legend.position = "right",
          legend.title = element_text(size = 11, face = "bold")) +
    labs(title = title, x = NULL, y = NULL)
}

# ============================================================================
# 3F. ORBIT-only pathway demo: TCA cycle leading edge
# ============================================================================
sig_rna   <- rea_rna   %>% filter(p.adjust < alpha_cut) %>% pull(ID)
sig_prot  <- rea_prot  %>% filter(p.adjust < alpha_cut) %>% pull(ID)
sig_orbit <- rea_orbit %>% filter(p.adjust < alpha_cut) %>% pull(ID)
orbit_only_ids <- setdiff(sig_orbit, union(sig_rna, sig_prot))
cat("ORBIT significant:", length(sig_orbit),
    "| ORBIT-only:", length(orbit_only_ids), "\n\n")

orbit_only_view <- rea_orbit %>%
  filter(ID %in% orbit_only_ids) %>%
  transmute(ID,
            pathway    = sub("^REACTOME_", "", ID),
            ORBIT_padj = p.adjust,
            ORBIT_NES  = NES,
            ORBIT_size = setSize) %>%
  left_join(rea_rna  %>% transmute(ID, RNA_padj  = p.adjust), by = "ID") %>%
  left_join(rea_prot %>% transmute(ID, Prot_padj = p.adjust), by = "ID") %>%
  arrange(ORBIT_padj) %>%
  select(pathway, ORBIT_padj, ORBIT_NES, ORBIT_size, RNA_padj, Prot_padj)
print(as.data.frame(orbit_only_view))

target_id     <- "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE"
pathway_title <- "Citric Acid Cycle (TCA)"

le_rna   <- get_leading(rea_rna,   target_id)
le_prot  <- get_leading(rea_prot,  target_id)
le_orbit <- get_leading(rea_orbit, target_id)
all_genes <- unique(c(le_rna, le_prot, le_orbit))
cat("LE sizes -> RNA:", length(le_rna), "Prot:", length(le_prot),
    "ORBIT:", length(le_orbit), "| union:", length(all_genes), "\n")

gene_table <- data.frame(Gene = all_genes, stringsAsFactors = FALSE) %>%
  left_join(master_table %>% transmute(Gene = Feature,
                                       RNA_logFC = RNA_logFC, RNA_padj = RNA_adjP,
                                       Prot_logFC = Prot_logFC, Prot_padj = Prot_adjP),
            by = "Gene") %>%
  mutate(in_RNA = Gene %in% le_rna, in_Prot = Gene %in% le_prot,
         in_ORBIT = Gene %in% le_orbit,
         le_pattern = paste0(ifelse(in_RNA, "R", ""), ifelse(in_Prot, "P", ""),
                             ifelse(in_ORBIT, "O", "")),
         le_order = case_when(le_pattern == "R" ~ 1, le_pattern == "RO" ~ 2,
                              le_pattern == "RP" ~ 3, le_pattern == "RPO" ~ 4,
                              le_pattern == "PO" ~ 5, le_pattern == "P" ~ 6,
                              le_pattern == "O" ~ 7, TRUE ~ 8)) %>%
  arrange(le_order, desc(RNA_logFC))

gene_order <- gene_table$Gene
plot_long  <- make_plot_long(gene_table) %>%
  mutate(Gene = factor(Gene, levels = rev(gene_order)))

y_levels  <- levels(plot_long$Gene)
gene_ypos <- setNames(seq_along(y_levels), y_levels)
rect_data <- bind_rows(
  mk_runs(le_rna,   gene_ypos, 0.6, 1.4, "RNA core"),
  mk_runs(le_prot,  gene_ypos, 1.6, 2.4, "Protein core"),
  mk_runs(le_orbit, gene_ypos, 0.7, 2.3, "ORBIT core")
) %>%
  mutate(ymin = start - 0.5, ymax = end + 0.5,
         core = factor(core, levels = c("RNA core", "Protein core", "ORBIT core")))

padj_rna   <- get_padj(rea_rna,   target_id)
padj_prot  <- get_padj(rea_prot,  target_id)
padj_orbit <- get_padj(rea_orbit, target_id)

p_3F <- base_bubble(plot_long, pathway_title, y_size = 15) +
  statebins:::geom_rrect(
    data = rect_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = core),
    fill = NA, linetype = "dashed", linewidth = 1.5,
    radius = grid::unit(20, "pt"), inherit.aes = FALSE) +
  scale_color_manual(
    name = "Pathway core (adjusted P)",
    values = c("RNA core" = "#4DAF4A", "Protein core" = "#377EB8", "ORBIT core" = "#E41A1C"),
    labels = c("RNA core"     = sprintf("Transcriptome (%s)", fmt_p(padj_rna)),
               "Protein core" = sprintf("Proteome (%s)",      fmt_p(padj_prot)),
               "ORBIT core"   = sprintf("ORBIT (%s)",         fmt_p(padj_orbit))),
    breaks = c("RNA core", "Protein core", "ORBIT core"),
    guide  = guide_legend(order = 1,
                          override.aes = list(fill = NA, linetype = "dashed", linewidth = 2))) +
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text  = element_text(size = 13),
        legend.key.size = unit(1, "cm"),
        legend.spacing.y = unit(1, "cm"))

quartz(type = "pdf", file = pdf_path("3F_reactome_TCA_bubble.pdf"),
       width = 8, height = 8, family = FONT)
print(p_3F)
dev.off()

# ============================================================================
# 3I. Reactome significant-pathway UpSet -- 4 methods incl. Fisher
# ============================================================================
set_names <- c("Transcriptome", "Proteome", "ORBIT", "Fisher")
set_cols  <- c(Transcriptome = "#3E9B5F",
               Proteome      = "#2C6FAD",
               ORBIT         = "#C0392B",
               Fisher        = "#E8A200")
U <- build_upset(set_names)
combo <- U$combo; combo_long <- U$combo_long; set_tot <- U$set_tot

x_scale      <- scale_x_discrete(limits = levels(combo$combo_id),
                                 expand = expansion(add = 0.6))
y_expand_dot <- scale_y_discrete(expand = expansion(add = 0.5))

p_bars <- ggplot(combo, aes(combo_id, size)) +
  geom_col(width = 0.6, fill = "grey25") +
  geom_text(aes(label = size), vjust = -0.4, size = 4, family = FONT) +
  x_scale +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), plot.margin = margin(5, 5, 0, 2))

p_dots <- ggplot(combo_long, aes(combo_id, set)) +
  geom_point(aes(colour = as.character(member)), size = 3.5) +
  geom_line(data = filter(combo_long, member),
            aes(group = combo_id), linewidth = 0.9, colour = "grey20") +
  scale_colour_manual(values = c("TRUE" = "grey20", "FALSE" = "grey85"),
                      guide = "none") +
  x_scale + y_expand_dot +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(axis.text  = element_blank(), axis.ticks = element_blank(),
        axis.line  = element_blank(),
        panel.background = element_rect(fill = "grey96", colour = NA),
        plot.margin = margin(0, 5, 5, 2))

p_title <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Pathways in intersection",
           angle = 90, family = FONT, size = 5) +
  theme_void()

p_names <- ggplot(set_tot, aes(x = 1, y = set)) +
  geom_text(aes(label = set), family = FONT, fontface = "bold",
            size = 4.5, hjust = 0.5) +
  y_expand_dot +
  theme_void()

p_left <- ggplot(set_tot, aes(size, set, fill = as.character(set))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = size), hjust = 1.2, size = 4, family = FONT) +
  scale_fill_manual(values = set_cols, guide = "none") +
  scale_x_reverse(expand = expansion(mult = c(0.2, 0))) +
  y_expand_dot +
  labs(x = "Total pathways", y = NULL) +
  base_theme +
  theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(), plot.margin = margin(2, 2, 5, 5))

p_3I <- p_empty + p_title + p_bars +
  p_left  + p_names + p_dots +
  plot_layout(ncol = 3, nrow = 2,
              widths  = c(1.6, 1.1, 2.8),
              heights = c(5.2, 1.5))

quartz(file = pdf_path("3I_fisher_compare.pdf"), type = "pdf", width = 7, height = 7)
print(p_3I)
dev.off()

# ============================================================================
# 3J. Fisher-enriched pathway inspection: MECP2 leading edge, four methods
# ============================================================================
target_id     <- "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_MECP2"
pathway_title <- "Transcriptional Regulation by MECP2"

le_rna    <- get_leading(rea_rna,    target_id)
le_prot   <- get_leading(rea_prot,   target_id)
le_orbit  <- get_leading(rea_orbit,  target_id)
le_fisher <- get_leading(rea_fisher, target_id)
cat("Leading edge sizes -> RNA:", length(le_rna), "Prot:", length(le_prot),
    "ORBIT:", length(le_orbit), "Fisher:", length(le_fisher), "\n")
all_genes <- unique(c(le_rna, le_prot, le_orbit, le_fisher))
cat("Union:", length(all_genes), "\n")

gene_table <- data.frame(Gene = all_genes, stringsAsFactors = FALSE) %>%
  left_join(master_table %>% transmute(Gene = Feature,
                                       RNA_logFC = RNA_logFC, RNA_padj = RNA_adjP,
                                       Prot_logFC = Prot_logFC, Prot_padj = Prot_adjP),
            by = "Gene") %>%
  mutate(in_RNA    = Gene %in% le_rna,
         in_Prot   = Gene %in% le_prot,
         in_ORBIT  = Gene %in% le_orbit,
         in_Fisher = Gene %in% le_fisher,
         n_omic    = in_RNA + in_Prot + in_ORBIT + in_Fisher,
         pin_top   = Gene %in% c("PTEN", "MOV10")) %>%      # pinned to top
  arrange(desc(pin_top), desc(in_Fisher), desc(n_omic), RNA_padj)

gene_order <- gene_table$Gene
plot_long  <- make_plot_long(gene_table) %>%
  mutate(Gene = factor(Gene, levels = rev(gene_order)))

y_levels  <- levels(plot_long$Gene)
gene_ypos <- setNames(seq_along(y_levels), y_levels)
rect_data <- bind_rows(
  mk_runs(le_rna,    gene_ypos, 0.55, 1.45, "RNA core"),
  mk_runs(le_prot,   gene_ypos, 1.55, 2.45, "Protein core"),
  mk_runs(le_orbit,  gene_ypos, 0.65, 2.35, "ORBIT core"),
  mk_runs(le_fisher, gene_ypos, 0.75, 2.25, "Fisher core")
) %>%
  mutate(ymin = start - 0.5, ymax = end + 0.5,
         core = factor(core, levels = c("RNA core", "Protein core", "ORBIT core", "Fisher core")))

padj_rna    <- get_padj(rea_rna,    target_id)
padj_prot   <- get_padj(rea_prot,   target_id)
padj_orbit  <- get_padj(rea_orbit,  target_id)
padj_fisher <- get_padj(rea_fisher, target_id)

p_3J <- base_bubble(plot_long, pathway_title, y_size = 14) +
  statebins:::geom_rrect(
    data = rect_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = core),
    fill = NA, linetype = "dashed", linewidth = 1.3,
    radius = grid::unit(18, "pt"), inherit.aes = FALSE) +
  scale_color_manual(
    name = "Pathway core (adjusted P)",
    values = c("RNA core"     = "#3E9B5F",
               "Protein core" = "#2C6FAD",
               "ORBIT core"   = "#C0392B",
               "Fisher core"  = "#E8A200"),
    labels = c("RNA core"     = sprintf("Transcriptome (%s)", fmt_p(padj_rna)),
               "Protein core" = sprintf("Proteome (%s)",      fmt_p(padj_prot)),
               "ORBIT core"   = sprintf("ORBIT (%s)",         fmt_p(padj_orbit)),
               "Fisher core"  = sprintf("Fisher (%s)",        fmt_p(padj_fisher))),
    breaks = c("RNA core", "Protein core", "ORBIT core", "Fisher core"),
    guide  = guide_legend(order = 1,
                          override.aes = list(fill = NA, linetype = "dashed", linewidth = 1.3))) +
  theme(legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 12),
        legend.key.size = unit(0.9, "cm"))

quartz(type = "pdf", file = pdf_path("3J_fisher_MECP2_bubble.pdf"),
       width = 8, height = 8, family = FONT)
print(p_3J)
dev.off()

# ============================================================================
# 3K. Spearman correlation of combined P (ORBIT / DPM / Stouffer) with
#     single-omic P, two-omics genes only; bootstrap 95% CI
# ============================================================================
# merged_results lives in TI_orbit.rds;
if ("merged_results" %in% names(TI_orbit)) {
  merged_results <- as.data.frame(TI_orbit$merged_results)
} else {
  merged_results <- as.data.frame(readRDS(pdf_path("TI_ORBIT.rds"))$merged_results)
}

keep <- c(ORBIT = "ORBIT_P", DPM = "DPM_P", Stouffer = "Stouffer_dir_P",
          RNA = "RNA_P", Protein = "Prot_P")
X <- as.matrix(merged_results[merged_results$n_omics == 2, keep])
colnames(X) <- names(keep)
X <- X[complete.cases(X), ]
V <- colnames(X); n_X <- nrow(X)

set.seed(42)
r0  <- cor(X, method = "spearman")
B   <- 5000
arr <- array(NA_real_, c(B, 5, 5), dimnames = list(NULL, V, V))
for (b in seq_len(B))
  arr[b, , ] <- cor(X[sample.int(n_X, n_X, TRUE), ], method = "spearman")

meth <- c("ORBIT", "DPM", "Stouffer")
bar  <- expand.grid(method = meth, omic = c("RNA", "Protein"), stringsAsFactors = FALSE)
bar$rho <- mapply(function(m, o) r0[m, o], bar$method, bar$omic)
bar$lo  <- mapply(function(m, o) quantile(arr[, m, o], .025), bar$method, bar$omic)
bar$hi  <- mapply(function(m, o) quantile(arr[, m, o], .975), bar$method, bar$omic)
bar$method <- factor(bar$method, meth)
bar$omic   <- factor(bar$omic, c("RNA", "Protein"))

p_3K <- ggplot(bar, aes(method, rho, fill = omic)) +
  geom_col(position = position_dodge(.85), width = .85) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.85),
                width = .2, linewidth = .7, colour = "black") +
  geom_text(aes(label = sprintf("%.2f", rho)), position = position_dodge(.85),
            vjust = -1.4, size = 6, fontface = "bold",
            family = FONT, show.legend = FALSE) +
  scale_fill_manual(values = c(RNA = "#3E9B5F", Protein = "#2C6FAD"),
                    labels = c(RNA = "Transcriptome", Protein = "Proteome"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 0.8), expand = expansion(mult = c(0, .02))) +
  labs(x = NULL, y = "Spearman \u03c1") +
  theme_minimal(base_size = 16, base_family = FONT) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.line    = element_blank(),
        legend.position = "top", legend.justification = "center",
        legend.text  = element_text(size = 18),
        axis.text.x  = element_text(face = "bold", size = 18),
        axis.text.y  = element_text(size = 22),
        axis.title.y = element_text(size = 22))

quartz(file = pdf_path("3K_method_correlation_bars.pdf"), type = "pdf", width = 6, height = 8)
print(p_3K)
dev.off()

# ============================================================================
# Sup Fig 4. Top-10 genes (N == 2) by Fisher (A) / DPM (B) / Stouffer (C)
#   two-omics bubble; "*" = not in ORBIT top10
# ============================================================================
bubble_dat <- merge(
  merged_results[, c("Feature", "n_omics", "DPM_P", "Stouffer_dir_P", "ORBIT_P")],
  master_table[, c("Feature", "fisher_P", "RNA_logFC", "Prot_logFC", "RNA_adjP", "Prot_adjP")],
  by = "Feature")

orbit_top10 <- bubble_dat %>% filter(n_omics == 2) %>%
  arrange(ORBIT_P) %>% slice_head(n = 10) %>% pull(Feature)

prepare_top10_bubble <- function(df, rank_col, ref_genes) {
  top_genes <- df %>%
    filter(n_omics == 2, !is.na(.data[[rank_col]])) %>%
    arrange(.data[[rank_col]]) %>% slice_head(n = 10) %>%
    mutate(label = ifelse(Feature %in% ref_genes, Feature, paste0(Feature, " *")))
  lab_order <- top_genes$label
  message(rank_col, " top10: ", paste(top_genes$Feature, collapse = ", "),
          "  | starred: ",
          paste(top_genes$Feature[!top_genes$Feature %in% ref_genes], collapse = ", "))
  top_genes %>%
    transmute(label,
              Transcriptome_logFC = RNA_logFC, Proteome_logFC = Prot_logFC,
              Transcriptome_adjP  = RNA_adjP,  Proteome_adjP  = Prot_adjP) %>%
    pivot_longer(cols = -label,
                 names_to = c("omics", ".value"), names_sep = "_") %>%
    mutate(neglog10_adjP = -log10(adjP),
           omics = factor(omics, levels = c("Transcriptome", "Proteome")),
           label = factor(label, levels = rev(lab_order)))
}

make_top10_bubble <- function(bubble_df) {
  fc_max <- max(abs(bubble_df$logFC), na.rm = TRUE); fc_lim <- c(-fc_max, fc_max)
  ggplot(bubble_df, aes(x = omics, y = label)) +
    geom_point(aes(size = neglog10_adjP, fill = logFC),
               shape = 21, colour = "black", stroke = 0.4) +
    scale_fill_gradient2(low = "#4C6FAD", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = fc_lim, name = expression(log[2]~FC)) +
    scale_size_continuous(range = c(4, 10), breaks = c(0.1, 0.5, 1.0, 5.0),
                          name = expression(-log[10]~adj.P)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 14, base_family = FONT) +
    theme(axis.text.x  = element_text(face = "bold", size = 15, family = FONT),
          axis.text.y  = element_text(face = "bold", size = 14, colour = "grey30", family = FONT),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.line    = element_blank(),
          legend.title = element_text(family = FONT),
          legend.text  = element_text(family = FONT),
          plot.margin  = margin(10, 20, 10, 10))
}

p_S4A <- make_top10_bubble(prepare_top10_bubble(bubble_dat, "fisher_P",       orbit_top10))
p_S4B <- make_top10_bubble(prepare_top10_bubble(bubble_dat, "DPM_P",          orbit_top10))
p_S4C <- make_top10_bubble(prepare_top10_bubble(bubble_dat, "Stouffer_dir_P", orbit_top10))

quartz(file = pdf_path("SupFig4A_fisher_top10_bubble.pdf"), type = "pdf", width = 6, height = 6)
print(p_S4A); dev.off()
quartz(file = pdf_path("SupFig4B_DPM_top10_bubble.pdf"), type = "pdf", width = 6, height = 6)
print(p_S4B); dev.off()
quartz(file = pdf_path("SupFig4C_Stouffer_top10_bubble.pdf"), type = "pdf", width = 6, height = 6)
print(p_S4C); dev.off()




# ============================================================================
# Sup Fig 5. ORBIT vs RankProd at the gene level (two-omics genes only)
#   A  Transcriptome vs Proteome signed significance, highlighting
#        opposite-direction genes
#   B  -log10 P, ORBIT (x) vs RankProd (y), highlighting opposite-direction genes
#
#   signed significance : sign(logFC) * -log10(P) per omic
#   Ranking universe = UNION of all tested genes (master_table, any N):
#     method top  : smallest TOP_FRAC of ORBIT_P / RankProd_P over all genes
#                   (no longer plotted; only reported in the console message)
#     omic top    : smallest OMIC_TOP_FRAC of RNA_P (resp. Prot_P) over all genes
#                   measured in that omic
#   Only two-omics genes (N == 2) can be plotted, so the panels show the
#   two-omics subset of those sets.
#   "opposite direction": sign(RNA_logFC) != sign(Prot_logFC) AND omic-top in both
# ============================================================================
TOP_FRAC      <- 0.20    # method top fraction (ORBIT / RankProd, by P) - message only
OMIC_TOP_FRAC <- 0.20    # per-omic "strong" fraction for the opposite-direction set
P_FLOOR       <- 1e-300

COL_OPP   <- "#7B3294"
COL_OTHER <- "grey80"

pct  <- paste0(TOP_FRAC * 100, "%")
opct <- paste0(OMIC_TOP_FRAC * 100, "%")
LAB_OPP  <- sprintf("Opposite direction,\ntop %s in both omics", opct)
lvls_opp <- c("Other", LAB_OPP)

# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------
# TRUE for the smallest `frac` of non-NA p (ranked over everything non-NA)
rank_top <- function(p, frac) {
  r <- rank(p, ties.method = "first", na.last = "keep")
  !is.na(r) & r <= ceiling(frac * sum(!is.na(p)))
}

# union universe: every gene in master_table, any N
all_df <- master_table %>%
  select(Feature, N, ORBIT_P = P, RNA_P, Prot_P, RNA_logFC, Prot_logFC) %>%
  left_join(TI_orbit$RankProd %>% select(Feature, RankProd_P), by = "Feature") %>%
  mutate(orbit_top = rank_top(ORBIT_P,    TOP_FRAC),
         rp_top    = rank_top(RankProd_P, TOP_FRAC),
         rna_top   = rank_top(RNA_P,      OMIC_TOP_FRAC),
         prot_top  = rank_top(Prot_P,     OMIC_TOP_FRAC))

# plotted subset: two-omics genes with everything available
gene_df <- all_df %>%
  filter(N == 2) %>%
  filter(if_all(c(ORBIT_P, RankProd_P, RNA_P, Prot_P, RNA_logFC, Prot_logFC), ~ !is.na(.x))) %>%
  mutate(opp_sig  = sign(RNA_logFC) != sign(Prot_logFC) & rna_top & prot_top,
         RNA_sig  = sign(RNA_logFC)  * -log10(pmax(RNA_P,  P_FLOOR)),
         Prot_sig = sign(Prot_logFC) * -log10(pmax(Prot_P, P_FLOOR)),
         nl_orbit = -log10(pmax(ORBIT_P,    P_FLOOR)),
         nl_rp    = -log10(pmax(RankProd_P, P_FLOOR)),
         cat      = factor(ifelse(opp_sig, LAB_OPP, "Other"), levels = lvls_opp))

message(sprintf(paste0(
  "SupFig5 universe: ORBIT %d | RankProd %d | RNA %d | Prot %d genes (union ranking)\n",
  "  method top %s -> ORBIT %d (two-omics %d) | RankProd %d (two-omics %d)\n",
  "  plotted two-omics genes = %d | opposite & top %s in both omics = %d\n",
  "  in ORBIT top: opposite %d | in RankProd top: opposite %d\n",
  "  Spearman rho(ORBIT P, RankProd P) = %.3f"),
  sum(!is.na(all_df$ORBIT_P)), sum(!is.na(all_df$RankProd_P)),
  sum(!is.na(all_df$RNA_P)),   sum(!is.na(all_df$Prot_P)),
  pct, sum(all_df$orbit_top), sum(gene_df$orbit_top),
  sum(all_df$rp_top),         sum(gene_df$rp_top),
  nrow(gene_df), opct, sum(gene_df$opp_sig),
  sum(gene_df$orbit_top & gene_df$opp_sig), sum(gene_df$rp_top & gene_df$opp_sig),
  cor(gene_df$ORBIT_P, gene_df$RankProd_P, method = "spearman")))

# ---------------------------------------------------------------------------
# shared theme / scales
# ---------------------------------------------------------------------------
theme_sup5 <- function() {
  theme_classic(base_size = 14, base_family = FONT) +
    theme(legend.position      = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.direction     = "vertical",
          legend.background    = element_rect(fill = alpha("white", 0.7),
                                              colour = "grey70", linewidth = 0.3),
          legend.text       = element_text(size = 11, family = FONT),
          legend.key.height = unit(0.9, "cm"),
          axis.title        = element_text(size = 15, family = FONT),
          axis.text         = element_text(size = 12, family = FONT),
          panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.line         = element_blank(),
          aspect.ratio      = 1,
          plot.tag          = element_text(size = 20, face = "bold", family = FONT),
          plot.margin       = margin(10, 15, 10, 10))
}

scale_opp <- function() {
  scale_colour_manual(values = setNames(c(COL_OTHER, COL_OPP), lvls_opp),
                      breaks = rev(lvls_opp), name = NULL, drop = FALSE)
}

guide_opp <- function() {
  guides(colour = guide_legend(override.aes = list(size = 3.5, alpha = 1)))
}

# ---------------------------------------------------------------------------
# A : signed transcriptome vs proteome, opposite highlighted
# ---------------------------------------------------------------------------
lim_A <- max(abs(c(gene_df$RNA_sig, gene_df$Prot_sig)), na.rm = TRUE) * 1.04

p_S5A <- ggplot(gene_df, aes(RNA_sig, Prot_sig, colour = cat)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(data = ~ filter(.x, cat == "Other"), size = 1.2, alpha = 0.5) +
  geom_point(data = ~ filter(.x, cat == LAB_OPP), size = 1.8, alpha = 0.85) +
  scale_opp() +
  coord_cartesian(xlim = c(-lim_A, lim_A), ylim = c(-lim_A, lim_A)) +
  labs(x = expression(Transcriptome:~sign(log[2]*FC) %*% -log[10]*P),
       y = expression(Proteome:~sign(log[2]*FC) %*% -log[10]*P)) +
  theme_sup5() +
  guide_opp()

# ---------------------------------------------------------------------------
# B : -log10 P, ORBIT vs RankProd, opposite highlighted
# ---------------------------------------------------------------------------
lim_B <- max(c(gene_df$nl_orbit, gene_df$nl_rp, -log10(0.05)), na.rm = TRUE) * 1.04

p_S5B <- ggplot(gene_df, aes(nl_orbit, nl_rp, colour = cat)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_point(data = ~ filter(.x, cat == "Other"), size = 1.2, alpha = 0.5) +
  geom_point(data = ~ filter(.x, cat == LAB_OPP), size = 1.8, alpha = 0.85) +
  scale_opp() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, lim_B), ylim = c(0, lim_B)) +
  labs(x = expression(ORBIT:~-log[10]*P),
       y = expression(RankProd:~-log[10]*P)) +
  theme_sup5() +
  guide_opp()

# ---------------------------------------------------------------------------
# assemble
# ---------------------------------------------------------------------------
p_S5 <- (p_S5A | p_S5B) +
  plot_annotation(tag_levels = "A")

quartz(file = pdf_path("SupFig5_orbit_vs_rankprod_gene_scatter.pdf"), type = "pdf",
       width = 12, height = 6.5)
print(p_S5)
dev.off()

