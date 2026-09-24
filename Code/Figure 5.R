# ============================================================================
# Figure 5 (DCM three-omics integration) -- plotting
#   Input  : Data/DCM/DCM_3omics_ORBIT_result.rds
#              $omics$Transcriptome  Feature, P, padj, Direction
#              $omics$Translatome    gene, logFC, P.Value, adj.P.Val
#              $omics$Proteome       Feature, ORBIT_Rank_P, ORBIT_Rank_adjP, ORBIT_Rank_Direction
#              $integration$ORBIT_Rank  Feature, N, P, Direction
#              $pathway$<layer>$gsea_objects$Reactome  gseaResult
#   Outputs: Data/DCM/*.pdf
#   5B normalised signed rank across Transcriptome / Translatome / Proteome
#   5C Reactome UpSet, 3 omics + ORBIT, bars coloured by 9 themes
#   5D mitochondrial Ca2+ transport leading-edge bubble
#   5G mitochondrial genes <-> pathways chord   5H lipid genes <-> pathways chord
#   (proteome-meta top10 bubble is Sup Fig 7, separate script)
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

for (pkg in c("dplyr", "tidyr", "ggplot2", "patchwork", "scales", "circlize", "statebins"))
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap", ask = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

FONT <- "Arial"
pdf_path <- function(f) file.path(DATA_DIR, f)

# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
DCM_3    <- readRDS(pdf_path("DCM_3omics_ORBIT_result.rds"))

DB <- "Reactome"
gsea_src <- list(
  Transcriptome = DCM_3$pathway$Transcriptome$gsea_objects[[DB]],
  Translatome   = DCM_3$pathway$Translatome$gsea_objects[[DB]],
  Proteome      = DCM_3$pathway$Proteome$gsea_objects[[DB]],
  ORBIT         = DCM_3$pathway$Integrated_3omics$gsea_objects[[DB]])
orbit_rank <- DCM_3$integration$ORBIT_Rank

layer_cols <- c(Transcriptome = "#3E9B5F",
                Translatome   = "#8E44AD",
                Proteome      = "#2C6FAD",
                ORBIT         = "#C0392B")
dir_cols   <- c(Down = "#2C6FAD", Up = "#C0392B")

get_leading <- function(g, pid) {
  row <- g@result %>% filter(ID == pid)
  if (nrow(row) == 0) return(character(0))
  strsplit(row$core_enrichment[1], "/")[[1]]
}
get_padj <- function(g, pid) {
  v <- g@result %>% filter(ID == pid) %>% pull(p.adjust)
  if (length(v) == 0) NA else v[1]
}
fmt_p <- function(p) if (is.na(p)) "NA" else sprintf("%.3f", p)

# ============================================================================
# 5B. Normalised signed rank across the three layers; ORBIT (N = 3) top 5 up
#     and top 5 down traced, names listed at the right
# ============================================================================
omics_levels <- c("Transcriptome", "Translatome", "Proteome", "ORBIT")
x_map <- setNames(seq_along(omics_levels), omics_levels)

orbit_n3 <- orbit_rank %>% filter(N == 3)
up5   <- orbit_n3 %>% filter(Direction ==  1) %>% arrange(P) %>% slice_head(n = 5) %>% pull(Feature)
down5 <- orbit_n3 %>% filter(Direction == -1) %>% arrange(P) %>% slice_head(n = 5) %>% pull(Feature)
gene_dir <- orbit_n3 %>%
  filter(Feature %in% c(up5, down5)) %>%
  transmute(Feature, grp = if_else(Direction == 1, "Up", "Down"), orbit_P = P)

# per layer: signal + sign -> nsr = sign * (1 - rank/(n+1))
nsr_layer <- function(df, feature, signal, sgn, layer)
  transmute(df, Feature = {{ feature }}, omics = layer, signal = {{ signal }}, sgn = {{ sgn }}) %>%
  filter(!is.na(signal), !is.na(sgn))
all_5B <- bind_rows(
  nsr_layer(DCM_3$omics$Transcriptome, Feature, -log10(P), sign(Direction), "Transcriptome"),
  nsr_layer(DCM_3$omics$Translatome,   gene,    -log10(P.Value), sign(logFC), "Translatome"),
  nsr_layer(DCM_3$omics$Proteome,      Feature, -log10(ORBIT_Rank_P), sign(ORBIT_Rank_Direction), "Proteome"),
  nsr_layer(orbit_rank,                Feature, -log10(P), sign(Direction), "ORBIT")) %>%
  group_by(omics) %>%
  mutate(n = n(), rank = rank(-signal, ties.method = "min"),
         nsr = sgn * (1 - rank / (n + 1))) %>%
  ungroup() %>%
  mutate(x = x_map[omics])

layers3 <- c("Transcriptome", "Translatome", "Proteome")
bg_5B   <- all_5B %>% filter(omics %in% layers3)
line_5B <- all_5B %>%
  filter(Feature %in% gene_dir$Feature, omics %in% layers3) %>%
  left_join(select(gene_dir, Feature, grp), by = "Feature") %>%
  mutate(grp = factor(grp, levels = c("Down", "Up")))

# equally spaced labels left of the ORBIT tick, ordered by ORBIT P within group
label_x <- x_map["ORBIT"] - 0.1
lab_5B <- gene_dir %>%
  group_by(grp) %>% arrange(orbit_P, .by_group = TRUE) %>%
  mutate(k = row_number(), m = n()) %>% ungroup() %>%
  mutate(y_lab = if_else(grp == "Up",
                         0.95 - (k - 1) * (0.55 / pmax(m - 1, 1)),
                         -0.95 + (k - 1) * (0.55 / pmax(m - 1, 1))),
         grp = factor(grp, levels = c("Down", "Up")), x_lab = label_x)
seg_5B <- lab_5B %>%
  left_join(line_5B %>% filter(omics == "Proteome") %>% transmute(Feature, y_prot = nsr), by = "Feature") %>%
  filter(!is.na(y_prot)) %>%
  mutate(x_start = x_map["Proteome"], x_end = label_x)

p_5B <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(data = bg_5B, aes(x = x, y = nsr), colour = "grey70", alpha = 0.4, size = 1) +
  geom_line(data = line_5B, aes(x = x, y = nsr, group = Feature, colour = grp),
            linewidth = 0.9, alpha = 0.8) +
  geom_segment(data = seg_5B, aes(x = x_start, xend = x_end, y = y_prot, yend = y_lab, colour = grp),
               linewidth = 0.9, alpha = 0.8) +
  geom_point(data = line_5B, aes(x = x, y = nsr, colour = grp),
             size = 2.5, shape = 21, fill = "white", stroke = 1.4) +
  geom_text(data = lab_5B, aes(x = x_lab, y = y_lab, label = Feature, colour = grp),
            family = FONT, size = 4.5, fontface = "bold", hjust = 0, show.legend = FALSE) +
  scale_colour_manual(values = dir_cols, name = NULL) +
  scale_x_continuous(breaks = x_map, labels = omics_levels,
                     limits = c(0.7, x_map["ORBIT"] + 0.7), expand = expansion(add = 0)) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
                     labels = c("-1\n(Down)", "-0.5", "0", "0.5", "1\n(Up)")) +
  labs(x = NULL, y = "Normalised Signed Rank", colour = NULL) +
  theme_classic(base_size = 14, base_family = FONT) +
  theme(axis.text.x  = element_text(face = "bold", size = 12, family = FONT),
        axis.text.y  = element_text(family = FONT),
        axis.title.y = element_text(size = 16, family = FONT),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        plot.margin  = margin(10, 40, 10, 10))

quartz(file = pdf_path("5B_3omics_signed_rank.pdf"), type = "pdf", width = 6, height = 6)
print(p_5B)
dev.off()

# ============================================================================
# 5C. Reactome significant-pathway UpSet: 3 omics + ORBIT, bars by 9 themes
#     (theme counts per intersection entered manually from Cytoscape clusters)
# ============================================================================
alpha_cut <- 0.05
sig_ids   <- function(g) g@result %>% filter(p.adjust < alpha_cut) %>% pull(ID)
set_names <- names(layer_cols)

sig_list  <- setNames(lapply(set_names, function(s) sig_ids(gsea_src[[s]])), set_names)
all_paths <- sort(unique(unlist(sig_list)))
upset_df  <- tibble(pathway = all_paths)
for (s in set_names) upset_df[[s]] <- upset_df$pathway %in% sig_list[[s]]
message("Set totals: ", paste(set_names, sapply(set_names, function(s) sum(upset_df[[s]])),
                              sep = "=", collapse = ", "))

combo <- upset_df %>%
  mutate(combo = apply(across(all_of(set_names)), 1, function(x) paste(as.integer(x), collapse = ""))) %>%
  count(combo, name = "size") %>%
  filter(combo != strrep("0", length(set_names))) %>%
  arrange(desc(size)) %>%
  mutate(combo_id = factor(row_number(), levels = as.character(row_number())))
combo_long <- combo %>%
  select(combo_id, combo) %>%
  tidyr::separate(combo, into = set_names, sep = "(?<=.)(?=.)", remove = FALSE) %>%
  pivot_longer(all_of(set_names), names_to = "set", values_to = "member") %>%
  mutate(set = factor(set, levels = rev(set_names)), member = member == "1")
set_tot <- tibble(set  = factor(set_names, levels = rev(set_names)),
                  size = sapply(set_names, function(s) sum(upset_df[[s]])))

theme_names <- c("Cell Cycle & Junction", "Innate Immune", "ECM Remodeling", "Mitochondria",
                 "Translation", "Phagocytosis", "Lipid Metabolism", "Hemostasis", "Lipoprotein")
cat_levels <- c(theme_names, "Other")
cat_cols <- c("Mitochondria" = "#5B8FB0", "Innate Immune" = "#A99CF0",       # down: cool
              "Cell Cycle & Junction" = "#9BDCEF", "Lipid Metabolism" = "#6FD0BE",
              "ECM Remodeling" = "#C77DC0", "Phagocytosis" = "#F0946B",      # up: warm
              "Hemostasis" = "#B98A5E", "Translation" = "#F080B0", "Lipoprotein" = "#F5D877",
              "Other" = "grey60")

# combo bits = Transcriptome / Translatome / Proteome / ORBIT
theme_counts <- tribble(
  ~combo, ~`Cell Cycle & Junction`, ~`Innate Immune`, ~`ECM Remodeling`, ~Mitochondria, ~Translation, ~Phagocytosis, ~`Lipid Metabolism`, ~Hemostasis, ~Lipoprotein,
  "1000",   24, 4, 1, 0,  0,  0, 1,  0, 0,   # Transcriptome only
  "0100",    0, 0, 1, 0, 20,  0, 3,  0, 0,   # Translatome only
  "0010",    0, 0, 0, 0,  1, 11, 0, 10, 7,   # Proteome only
  "0001",    0, 1, 3, 3,  0,  0, 6,  0, 0,   # ORBIT only
  "1010",    0, 0, 0, 0,  0,  1, 0,  0, 0,   # Transcriptome + Proteome
  "1001",    0, 9, 2, 2,  0,  0, 0,  0, 0,   # Transcriptome + ORBIT
  "0101",    0, 0, 0, 1,  0,  0, 6,  0, 0,   # Translatome + ORBIT
  "0011",    0, 0, 2, 6,  1,  0, 0,  1, 0,   # Proteome + ORBIT
  "1101",    0, 1, 0, 1,  0,  0, 1,  0, 0,   # Transcriptome + Translatome + ORBIT
  "1011",    0, 0, 9, 7,  0,  0, 0,  0, 0,   # Transcriptome + Proteome + ORBIT
  "0111",    0, 0, 1, 0,  1,  0, 1,  0, 0,   # Translatome + Proteome + ORBIT
  "1111",    0, 0, 3, 1,  0,  0, 0,  0, 0)   # all four

combo_cat <- combo %>%
  select(combo, size, combo_id) %>%
  left_join(theme_counts, by = "combo") %>%
  mutate(across(all_of(theme_names), ~ replace_na(., 0)),
         theme_sum = rowSums(across(all_of(theme_names))),
         Other     = size - theme_sum)
print(as_tibble(select(combo_cat, combo, size, theme_sum, Other)), n = Inf)
bad <- filter(combo_cat, Other < 0)
if (nrow(bad)) warning("theme sum exceeds intersection size: ",
                       paste(sprintf("%s(size=%d, sum=%d)", bad$combo, bad$size, bad$theme_sum), collapse = "; "))
missing_combo <- base::setdiff(theme_counts$combo, combo$combo)
if (length(missing_combo)) warning("theme_counts combos absent from the UpSet: ",
                                   paste(missing_combo, collapse = ", "))

combo_cat_long <- combo_cat %>%
  select(combo_id, all_of(theme_names), Other) %>%
  pivot_longer(all_of(cat_levels), names_to = "category", values_to = "count") %>%
  mutate(category = factor(category, levels = cat_levels))

x_scale      <- scale_x_discrete(limits = levels(combo$combo_id), expand = expansion(add = 0.6))
y_expand_dot <- scale_y_discrete(expand = expansion(add = 0.5))
base_theme   <- theme_classic(base_size = 14, base_family = FONT)

p_bars <- ggplot(combo_cat_long, aes(combo_id, count, fill = category)) +
  geom_col(width = 0.6) +
  geom_text(data = combo, aes(combo_id, size, label = size),
            inherit.aes = FALSE, vjust = -0.4, size = 4, family = FONT) +
  scale_fill_manual(values = cat_cols, name = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  x_scale +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), plot.margin = margin(5, 5, 0, 2),
        legend.position = "inside", legend.position.inside = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.6), colour = NA),
        legend.margin = margin(2, 4, 2, 4), legend.key.size = unit(0.8, "lines"),
        legend.spacing.y = unit(1, "pt"), legend.text = element_text(size = 11, family = FONT))

p_dots <- ggplot(combo_long, aes(combo_id, set)) +
  geom_point(aes(colour = as.character(member)), size = 3.5) +
  geom_line(data = filter(combo_long, member), aes(group = combo_id), linewidth = 0.9, colour = "grey20") +
  scale_colour_manual(values = c("TRUE" = "grey20", "FALSE" = "grey85"), guide = "none") +
  x_scale + y_expand_dot +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        panel.background = element_rect(fill = "grey96", colour = NA),
        plot.margin = margin(0, 5, 5, 2))

p_title <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Pathways in intersection", angle = 90, family = FONT, size = 5) +
  theme_void()
p_names <- ggplot(set_tot, aes(x = 1, y = set)) +
  geom_text(aes(label = set), family = FONT, fontface = "bold", size = 4.5, hjust = 0.5) +
  y_expand_dot + theme_void()
p_left <- ggplot(set_tot, aes(size, set, fill = as.character(set))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = size), hjust = 1.2, size = 4, family = FONT) +
  scale_fill_manual(values = layer_cols, guide = "none") +
  scale_x_reverse(expand = expansion(mult = c(0.2, 0))) +
  y_expand_dot +
  labs(x = "Total pathways", y = NULL) +
  base_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), plot.margin = margin(0, 2, 5, 5))
p_empty <- ggplot() + theme_void()

p_5C <- p_empty + p_title + p_bars +
  p_left  + p_names + p_dots +
  plot_layout(ncol = 3, nrow = 2, widths = c(1.2, 1, 4), heights = c(2.5, 1))

quartz(file = pdf_path("5C_reactome_upset_3omics.pdf"), type = "pdf", width = 10, height = 5)
print(p_5C)
dev.off()

# ============================================================================
# 5D. Single-pathway leading-edge bubble: mitochondrial Ca2+ ion transport
#     size = -log10 adj.P, fill = direction, dashed rounded box = layer core
# ============================================================================
target_id     <- "REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT"
pathway_title <- "Mitochondrial Calcium Ion Transport"

le <- lapply(gsea_src, get_leading, pid = target_id)
all_genes <- unique(unlist(le))
cat("LE ->", paste(names(le), lengths(le), collapse = " "), "| union:", length(all_genes), "\n")

prep_layer <- function(df, gene_col, padj_col, dir_expr, layer)
  df %>% filter(.data[[gene_col]] %in% all_genes) %>%
  transmute(Gene = .data[[gene_col]], Omic = layer, dir = {{ dir_expr }},
            size = -log10(pmax(.data[[padj_col]], 1e-300)))
plot_5D <- bind_rows(
  prep_layer(DCM_3$omics$Transcriptome, "Feature", "padj",
             if_else(Direction > 0, "Up", "Down"), "Transcriptome"),
  prep_layer(DCM_3$omics$Translatome, "gene", "adj.P.Val",
             if_else(logFC > 0, "Up", "Down"), "Translatome"),
  prep_layer(DCM_3$omics$Proteome, "Feature", "ORBIT_Rank_adjP",
             if_else(ORBIT_Rank_Direction > 0, "Up", "Down"), "Proteome")) %>%
  filter(!is.na(size), !is.na(dir))

# gene order: pinned genes first, then by leading-edge membership pattern
gene_pattern <- tibble(Gene = all_genes) %>%
  mutate(in_R = Gene %in% le$Transcriptome, in_T = Gene %in% le$Translatome,
         in_P = Gene %in% le$Proteome, in_O = Gene %in% le$ORBIT,
         n_omic = in_R + in_T + in_P + in_O,
         key = paste0(as.integer(in_R), as.integer(in_T), as.integer(in_P), as.integer(in_O))) %>%
  arrange(desc(n_omic), key, Gene)
manual_top <- intersect(c("MICU3", "PMPCB", "PHB1", "MCU", "VDAC2", "AFG3L2"), all_genes)
gene_order <- c(manual_top, base::setdiff(gene_pattern$Gene, manual_top))
move_before <- function(order, gene, before) {
  if (!(gene %in% order) || !(before %in% order)) return(order)
  order <- order[order != gene]
  append(order, gene, after = which(order == before) - 1)
}
gene_order <- move_before(gene_order, "YME1L1", "MAIP1")

plot_5D <- plot_5D %>%
  mutate(Omic = factor(Omic, levels = layers3),
         Gene = factor(Gene, levels = rev(gene_order)),
         dir  = factor(dir, levels = c("Down", "Up")))
max_size    <- max(plot_5D$size, na.rm = TRUE)
size_breaks <- pretty(c(0, max_size), n = 4); size_breaks <- size_breaks[size_breaks <= max_size]

# dashed rounded boxes over contiguous leading-edge runs
y_levels  <- levels(plot_5D$Gene)
gene_ypos <- setNames(seq_along(y_levels), y_levels)
find_runs <- function(pos) {
  if (length(pos) == 0) return(data.frame(start = integer(), end = integer()))
  pos <- sort(pos)
  if (length(pos) == 1) return(data.frame(start = pos, end = pos))
  brk <- which(diff(pos) > 1)
  data.frame(start = pos[c(1, brk + 1)], end = pos[c(brk, length(pos))])
}
mk_runs <- function(genes, xmin, xmax, core)
  find_runs(unname(gene_ypos[names(gene_ypos) %in% genes])) %>%
  mutate(xmin = xmin, xmax = xmax, core = core)
rect_5D <- bind_rows(
  mk_runs(le$Transcriptome, 0.6,  1.4,  "Transcriptome"),
  mk_runs(le$Translatome,   1.6,  2.4,  "Translatome"),
  mk_runs(le$Proteome,      2.6,  3.4,  "Proteome"),
  mk_runs(le$ORBIT,         0.65, 3.35, "ORBIT")) %>%
  mutate(ymin = start - 0.5, ymax = end + 0.5, core = factor(core, levels = names(layer_cols)))
padj_5D <- sapply(gsea_src, get_padj, pid = target_id)

p_5D <- ggplot(plot_5D, aes(x = Omic, y = Gene)) +
  geom_point(aes(fill = dir, size = size), shape = 21, color = "black", stroke = 0.6, alpha = 0.9) +
  scale_fill_manual(values = dir_cols, name = "Direction",
                    guide = guide_legend(order = 2, ncol = 1, position = "right",
                                         override.aes = list(size = 5))) +
  scale_size_continuous(range = c(1, 14), name = expression(-log[10]~adj.P),
                        breaks = size_breaks, limits = c(0, max_size),
                        guide = guide_legend(order = 1, ncol = 1, position = "right")) +
  scale_x_discrete(limits = layers3) +
  statebins:::geom_rrect(
    data = rect_5D, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = core),
    fill = NA, linetype = "dashed", linewidth = 1.2,
    radius = grid::unit(16, "pt"), inherit.aes = FALSE) +
  scale_color_manual(
    name = "Pathway core (adj.P)", values = layer_cols, breaks = names(layer_cols),
    labels = sprintf("%s (%s)", names(layer_cols), sapply(padj_5D[names(layer_cols)], fmt_p)),
    guide = guide_legend(position = "top", nrow = 2, title.position = "top",
                         override.aes = list(fill = NA, linetype = "dashed", linewidth = 3))) +
  labs(title = pathway_title, x = NULL, y = NULL) +
  theme_bw(base_size = 13, base_family = FONT) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 16, face = "bold"),
        plot.title  = element_text(face = "bold", hjust = 0.5, size = 18),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text  = element_text(size = 11),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(2, 2, 2, 2))

quartz(type = "pdf", file = pdf_path("5D_reactome_mitoCa_bubble.pdf"),
       width = 7, height = 12, family = FONT)
print(p_5D)
dev.off()

# ============================================================================
# Chord diagram helper (5G, 5H): genes <-> ORBIT-layer Reactome pathways
#   edge = gene in pathway leading edge
#   gene sector colour = ORBIT rank (smaller rank = deeper)
#   pathway sector colour = rank of -log10 padj (deeper = more significant)
# ============================================================================
orbit_ranked <- orbit_rank %>% mutate(orbit_rank = rank(P, ties.method = "min"))
n_total      <- nrow(orbit_ranked)
orbit_res    <- as.data.frame(gsea_src$ORBIT@result)

draw_chord <- function(genes, paths, short_map, ramp_cols, file,
                       width = 9, height = 9, start_degree = 100, cex = 0.85) {
  res  <- orbit_res[match(paths, orbit_res$ID), ]
  core <- setNames(strsplit(res$core_enrichment, "/"), paths)
  
  edges <- do.call(rbind, lapply(paths, function(p) {
    hit <- intersect(genes, core[[p]])
    if (length(hit) == 0) return(NULL)
    data.frame(gene = hit, pathway = short_map[[p]], value = 1, stringsAsFactors = FALSE)
  }))
  hits <- sapply(paths, function(p) length(intersect(genes, core[[p]])))
  print(data.frame(pathway = short_map[paths], n_hits = hits, row.names = NULL), row.names = FALSE)
  paths_used <- paths[hits > 0]
  genes      <- genes[genes %in% edges$gene]           # drop genes with no edge
  
  ramp <- colorRampPalette(ramp_cols)
  rank_col <- function(x) ramp(100)[cut(rank(x, ties.method = "average"), 100, labels = FALSE)]
  gene_rank <- setNames(orbit_ranked$orbit_rank[match(genes, orbit_ranked$Feature)], genes)
  path_padj <- setNames(res$p.adjust[match(paths_used, paths)], short_map[paths_used])
  path_x    <- -log10(path_padj); path_x[!is.finite(path_x)] <- max(path_x[is.finite(path_x)])
  grid.col  <- setNames(c(rank_col(-gene_rank), rank_col(path_x)),
                        c(genes, short_map[paths_used]))
  order_sectors <- c(genes, short_map[paths_used])
  zidx   <- setNames(seq_along(order_sectors), order_sectors)[edges$pathway]
  n_gene <- length(genes); n_path <- length(paths_used)
  
  quartz(file = pdf_path(file), type = "pdf", width = width, height = height, family = FONT)
  circos.clear()
  circos.par(start.degree = start_degree, clock.wise = FALSE,
             canvas.xlim = c(-0.95, 1.2), canvas.ylim = c(-0.92, 0.92),
             gap.after = c(rep(1, n_gene - 1), 30, rep(1, n_path - 1), 30))
  chordDiagram(edges, order = order_sectors, grid.col = grid.col,
               link.lwd = 0.4, link.border = "grey85", link.sort = TRUE,
               link.decreasing = FALSE, link.zindex = zidx,
               annotationTrack = "grid", preAllocateTracks = list(track.height = 0.20))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    s <- get.cell.meta.data("sector.index"); xl <- get.cell.meta.data("xlim")
    circos.text(mean(xl), get.cell.meta.data("ylim")[1] + 0.02, s,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = cex, family = FONT)
  }, bg.border = NA)
  circos.clear()
  
  lg_path <- Legend(col_fun = colorRamp2(seq(min(path_x), max(path_x), length = 5), ramp(5)),
                    title = "Pathway  -log10(padj)", direction = "horizontal")
  rr <- range(gene_rank, na.rm = TRUE)
  lg_gene <- Legend(col_fun = colorRamp2(seq(-rr[2], -rr[1], length = 5), ramp(5)),
                    title = "Gene  ORBIT rank (top %)", direction = "horizontal",
                    at = seq(-rr[2], -rr[1], length = 5),
                    labels = sprintf("%.1f%%", 100 * seq(rr[2], rr[1], length = 5) / n_total))
  draw(packLegend(lg_gene, lg_path, direction = "vertical"),
       x = unit(1, "npc") - unit(4, "mm"), y = unit(4, "mm"), just = c("right", "bottom"))
  dev.off()
  invisible(edges)
}

# ============================================================================
# 5G. Mitochondrial genes (42, all down) <-> 16 mitochondrial pathways
# ============================================================================
genesMito <- c("LRPPRC","CRAT","TIMM44","GRPEL1","MRPS34","MRPS25","MRPS16","MRPS12","MRPS17",
               "MRPL42","MTIF2","MTIF3","HSCB","TTC19","UQCRFS1","NDUFAF2","NDUFAB1","NDUFA5",
               "NDUFV2","UQCRB","COQ10A","MT-ND2","COX5A","TMEM223","GOT1","IDH2","ACACB","MOCS1",
               "IDH1","LETM1","NAMPT","GPT","MOCS2","SLC5A6","LMBRD1","CD38","NAXE","PANK2",
               "FLAD1","BCO2","PDXK","COQ8A")
paths_mito <- c(
  "REACTOME_UBIQUINOL_BIOSYNTHESIS",
  "REACTOME_METABOLISM_OF_COFACTORS",
  "REACTOME_METABOLISM_OF_VITAMINS_AND_COFACTORS",
  "REACTOME_METABOLISM_OF_WATER_SOLUBLE_VITAMINS_AND_COFACTORS",
  "REACTOME_BIOTIN_TRANSPORT_AND_METABOLISM",
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE",
  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "REACTOME_COMPLEX_IV_ASSEMBLY",
  "REACTOME_COMPLEX_I_BIOGENESIS",
  "REACTOME_COMPLEX_III_ASSEMBLY",
  "REACTOME_MITOCHONDRIAL_TRANSLATION",
  "REACTOME_MITOCHONDRIAL_TRANSLATION_ELONGATION",
  "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT",
  "REACTOME_PROTEIN_LOCALIZATION",
  "REACTOME_MITOCHONDRIAL_RNA_DEGRADATION")
short_mito <- c(
  REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE                            = "TCA Cycle",
  REACTOME_RESPIRATORY_ELECTRON_TRANSPORT                         = "Respiratory ETC",
  REACTOME_COMPLEX_I_BIOGENESIS                                   = "Complex I Biogenesis",
  REACTOME_COMPLEX_III_ASSEMBLY                                   = "Complex III Assembly",
  REACTOME_MITOCHONDRIAL_TRANSLATION                              = "Mito Translation",
  REACTOME_COMPLEX_IV_ASSEMBLY                                    = "Complex IV Assembly",
  REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT = "Aerobic Respiration",
  REACTOME_MITOCHONDRIAL_RNA_DEGRADATION                          = "Mito RNA Degradation",
  REACTOME_BIOTIN_TRANSPORT_AND_METABOLISM                        = "Biotin Transport/Metab.",
  REACTOME_METABOLISM_OF_VITAMINS_AND_COFACTORS                   = "Vitamins & Cofactors",
  REACTOME_UBIQUINOL_BIOSYNTHESIS                                 = "Ubiquinol Biosynthesis",
  REACTOME_METABOLISM_OF_COFACTORS                                = "Cofactor Metabolism",
  REACTOME_METABOLISM_OF_WATER_SOLUBLE_VITAMINS_AND_COFACTORS     = "Water-sol. Vit. & Cofactors",
  REACTOME_MITOCHONDRIAL_TRANSLATION_ELONGATION                   = "Mito Translation Elong.",
  REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT                           = "Mito Protein Import",
  REACTOME_PROTEIN_LOCALIZATION                                   = "Protein Localization")

draw_chord(genesMito, paths_mito, short_mito,
           ramp_cols = c("#EAF1F6", "#CBDDE9", "#A3C4D8", "#5B8FB0", "#3E7196", "#2A567A", "#1A3A57"),
           file = "5G_mito_pathways_chord.pdf")

# ============================================================================
# 5H. Lipid genes (31, all down) <-> 14 lipid pathways
# ============================================================================
genesLipid <- c("PITPNM1","ETNPPL","PLA2G4F","DGAT2","GPD2","GPD1","ABHD3","CHPT1","PTDSS1",
                "LPCAT3","AGPAT3","THEM4","FADS1","NDUFAB1","PCTP","ACSL1","ACSF2","DECR1",
                "ACADVL","ALOX5AP","ALOX5","CRAT","MID1IP1","ACACB","RXRA","IDI1","DHCR7",
                "MBTPS2","SC5D","NSDHL","SREBF1")
paths_lipid <- c(
  "REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP",
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS_VIA_DESMOSTEROL_BLOCH_PATHWAY",
  "REACTOME_METABOLISM_OF_STEROIDS",
  "REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF",
  "REACTOME_CARNITINE_SHUTTLE",
  "REACTOME_FATTY_ACID_METABOLISM",
  "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
  "REACTOME_FATTY_ACYL_COA_BIOSYNTHESIS",
  "REACTOME_ALPHA_LINOLENIC_OMEGA3_AND_LINOLEIC_OMEGA6_ACID_METABOLISM",
  "REACTOME_SYNTHESIS_OF_PA",
  "REACTOME_TRIGLYCERIDE_METABOLISM",
  "REACTOME_GLYCEROPHOSPHOLIPID_BIOSYNTHESIS",
  "REACTOME_PHOSPHOLIPID_METABOLISM")
short_lipid <- c(
  REACTOME_ACTIVATION_OF_GENE_EXPRESSION_BY_SREBF_SREBP               = "SREBP Gene Activation",
  REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF      = "Cholesterol Reg. by SREBP",
  REACTOME_METABOLISM_OF_STEROIDS                                     = "Steroid Metabolism",
  REACTOME_CHOLESTEROL_BIOSYNTHESIS                                   = "Cholesterol Biosynthesis",
  REACTOME_CHOLESTEROL_BIOSYNTHESIS_VIA_DESMOSTEROL_BLOCH_PATHWAY     = "Cholesterol (BLOCH)",
  REACTOME_SYNTHESIS_OF_PA                                            = "PA Synthesis",
  REACTOME_PHOSPHOLIPID_METABOLISM                                    = "Phospholipid Metab.",
  REACTOME_GLYCEROPHOSPHOLIPID_BIOSYNTHESIS                           = "Glycerophospholipid Biosynth.",
  REACTOME_CARNITINE_SHUTTLE                                          = "Carnitine Shuttle",
  REACTOME_FATTY_ACYL_COA_BIOSYNTHESIS                                = "Fatty Acyl-CoA Biosynth.",
  REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION                    = "Mito FA \u03b2-Oxidation",
  REACTOME_ALPHA_LINOLENIC_OMEGA3_AND_LINOLEIC_OMEGA6_ACID_METABOLISM = "PUFA (\u03c93/\u03c96) Metab.",
  REACTOME_FATTY_ACID_METABOLISM                                      = "Fatty Acid Metabolism",
  REACTOME_TRIGLYCERIDE_METABOLISM                                    = "Triglyceride Metab.")

draw_chord(genesLipid, paths_lipid, short_lipid,
           ramp_cols = c("#E4DEF5", "#CDC2EE", "#B0A1E6", "#8B7AE8", "#6F5DD6", "#5544B8", "#382A8F"),
           file = "5H_lipid_pathways_chord.pdf")