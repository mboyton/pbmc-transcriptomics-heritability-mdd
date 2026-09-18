library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

# paths
mr_file_path <- "path/to/psoriasis_mr_results.txt.gz"

nk_lps_wgcna_path <- "path/to/nk_lps_wgcna_results.rda"
cd4t_mem_wgcna_path <- "path/to/cd4t_mem_cd3_cd28_wgcna_results.rda"
cd4t_naive_wgcna_path <- "path/to/cd4t_naive_cd3_cd28_wgcna_results.rda"
cd8t_mem_wgcna_path <- "path/to/cd8t_mem_cd3_cd28_wgcna_results.rda"
cd8t_naive_wgcna_path <- "path/to/cd8t_naive_cd3_cd28_wgcna_results.rda"
nk_cd3_cd28_wgcna_path <- "path/to/nk_cd3_cd28_wgcna_results.rda"

# settings
mr_method <- "Inverse variance weighted"
mr_pval_cutoff <- 0.05

kme_threshold <- 0.75
module_trait_fdr_cutoff <- 0.05


# load psoriasis MR results
mr <- fread(mr_file_path)

mr_ivw <- mr %>%
 filter(
  method == mr_method
 ) %>%
 mutate(
  b = as.numeric(b),
  se = as.numeric(se),
  pval = as.numeric(pval),
  z = b / se,
  abs_z = abs(z),
  OR = exp(b),
  OR_lower = exp(b - 1.96 * se),
  OR_upper = exp(b + 1.96 * se)
 ) %>%
 filter(
  !is.na(gene),
  !is.na(b),
  !is.na(se),
  !is.na(pval),
  !is.na(z)
 )


# empty container
all_core_mr <- list()


# NK rest vs LPS
load(nk_lps_wgcna_path)

cell_state <- "nk_rest_vs_lps"

core_wgcna_genes <- gene_module_table %>%
 mutate(
  kME = as.numeric(kME)
 ) %>%
 filter(
  !is.na(kME),
  kME >= kme_threshold,
  !is.na(module_trait_fdr),
  module_trait_fdr < module_trait_fdr_cutoff
 )

core_wgcna_mr <- mr_ivw %>%
 inner_join(
  core_wgcna_genes,
  by = "gene"
 ) %>%
 mutate(
  cell_state = cell_state,
  pval_fdr_within_cellstate =
   p.adjust(
    pval,
    method = "fdr"
   )
 )

all_core_mr[[cell_state]] <- core_wgcna_mr


# CD4T memory rest vs CD3/CD28
load(cd4t_mem_wgcna_path)

cell_state <- "cd4t_mem_rest_vs_cd3_cd28"

core_wgcna_genes <- gene_module_table %>%
 mutate(
  kME = as.numeric(kME)
 ) %>%
 filter(
  !is.na(kME),
  kME >= kme_threshold,
  !is.na(module_trait_fdr),
  module_trait_fdr < module_trait_fdr_cutoff
 )

core_wgcna_mr <- mr_ivw %>%
 inner_join(
  core_wgcna_genes,
  by = "gene"
 ) %>%
 mutate(
  cell_state = cell_state,
  pval_fdr_within_cellstate =
   p.adjust(
    pval,
    method = "fdr"
   )
 )

all_core_mr[[cell_state]] <- core_wgcna_mr


# CD4T naive rest vs CD3/CD28
load(cd4t_naive_wgcna_path)

cell_state <- "cd4t_naive_rest_vs_cd3_cd28"

core_wgcna_genes <- gene_module_table %>%
 mutate(
  kME = as.numeric(kME)
 ) %>%
 filter(
  !is.na(kME),
  kME >= kme_threshold,
  !is.na(module_trait_fdr),
  module_trait_fdr < module_trait_fdr_cutoff
 )

core_wgcna_mr <- mr_ivw %>%
 inner_join(
  core_wgcna_genes,
  by = "gene"
 ) %>%
 mutate(
  cell_state = cell_state,
  pval_fdr_within_cellstate =
   p.adjust(
    pval,
    method = "fdr"
   )
 )

all_core_mr[[cell_state]] <- core_wgcna_mr


# CD8T memory rest vs CD3/CD28
load(cd8t_mem_wgcna_path)

cell_state <- "cd8t_mem_rest_vs_cd3_cd28"

core_wgcna_genes <- gene_module_table %>%
 mutate(
  kME = as.numeric(kME)
 ) %>%
 filter(
  !is.na(kME),
  kME >= kme_threshold,
  !is.na(module_trait_fdr),
  module_trait_fdr < module_trait_fdr_cutoff
 )

core_wgcna_mr <- mr_ivw %>%
 inner_join(
  core_wgcna_genes,
  by = "gene"
 ) %>%
 mutate(
  cell_state = cell_state,
  pval_fdr_within_cellstate =
   p.adjust(
    pval,
    method = "fdr"
   )
 )

all_core_mr[[cell_state]] <- core_wgcna_mr


# CD8T naive rest vs CD3/CD28
load(cd8t_naive_wgcna_path)

cell_state <- "cd8t_naive_rest_vs_cd3_cd28"

core_wgcna_genes <- gene_module_table %>%
 mutate(
  kME = as.numeric(kME)
 ) %>%
 filter(
  !is.na(kME),
  kME >= kme_threshold,
  !is.na(module_trait_fdr),
  module_trait_fdr < module_trait_fdr_cutoff
 )

core_wgcna_mr <- mr_ivw %>%
 inner_join(
  core_wgcna_genes,
  by = "gene"
 ) %>%
 mutate(
  cell_state = cell_state,
  pval_fdr_within_cellstate =
   p.adjust(
    pval,
    method = "fdr"
   )
 )

all_core_mr[[cell_state]] <- core_wgcna_mr


# NK rest vs CD3/CD28
load(nk_cd3_cd28_wgcna_path)

cell_state <- "nk_rest_vs_cd3_cd28"

core_wgcna_genes <- gene_module_table %>%
 mutate(
  kME = as.numeric(kME)
 ) %>%
 filter(
  !is.na(kME),
  kME >= kme_threshold,
  !is.na(module_trait_fdr),
  module_trait_fdr < module_trait_fdr_cutoff
 )

core_wgcna_mr <- mr_ivw %>%
 inner_join(
  core_wgcna_genes,
  by = "gene"
 ) %>%
 mutate(
  cell_state = cell_state,
  pval_fdr_within_cellstate =
   p.adjust(
    pval,
    method = "fdr"
   )
 )

all_core_mr[[cell_state]] <- core_wgcna_mr


# combine cell states
all_core_mr_df <- bind_rows(all_core_mr)


# nominally significant MR results
psoriasis_mr_sig <- all_core_mr_df %>%
 filter(
  pval < mr_pval_cutoff
 )


# unique genes for MR forest plot
plot_genes <- psoriasis_mr_sig %>%
 group_by(
  gene,
  gene_symbol
 ) %>%
 summarise(
  z = first(z),
  abs_z = first(abs_z),
  pval = first(pval),
  OR = first(OR),
  OR_lower = first(OR_lower),
  OR_upper = first(OR_upper),
  .groups = "drop"
 ) %>%
 arrange(OR) %>%
 mutate(
  gene_label = factor(
   gene_symbol,
   levels = gene_symbol
  )
 )


# MR forest plot
p_psoriasis_mr <- ggplot(
 plot_genes,
 aes(
  x = gene_label,
  y = OR,
  ymin = OR_lower,
  ymax = OR_upper
 )
) +
 geom_hline(
  yintercept = 1,
  linetype = "dashed"
 ) +
 geom_pointrange() +
 coord_flip() +
 labs(
  x = NULL,
  y = "Odds ratio",
  title = "Psoriasis MR results",
  subtitle = paste0(
   "Nominal MR p < ",
   mr_pval_cutoff,
   "; WGCNA kME >= ",
   kme_threshold
  )
 ) +
 theme_bw() +
 theme(
  plot.title = element_text(
   face = "bold"
  )
 )


# clean cell-state labels
clean_cell_state <- function(x) {
 x %>%
  gsub("_", " ", .) %>%
  gsub("cd4t", "CD4T", ., ignore.case = TRUE) %>%
  gsub("cd8t", "CD8T", ., ignore.case = TRUE) %>%
  gsub("nk", "NK", ., ignore.case = TRUE) %>%
  gsub("mem", "Mem", ., ignore.case = TRUE) %>%
  gsub("naive", "Naive", ., ignore.case = TRUE) %>%
  gsub("rest vs ", "", ., ignore.case = TRUE) %>%
  gsub(" vs ", " + ", .) %>%
  gsub("cd3 cd28", "CD3/CD28", ., ignore.case = TRUE) %>%
  gsub("lps", "LPS", ., ignore.case = TRUE)
}


# transcriptional heatmap
gene_cell_heat <- psoriasis_mr_sig %>%
 semi_join(
  plot_genes,
  by = c(
   "gene",
   "gene_symbol"
  )
 ) %>%
 distinct(
  gene,
  gene_symbol,
  cell_state,
  mean_logcpm_diff_stim_minus_base
 ) %>%
 mutate(
  gene_label = factor(
   gene_symbol,
   levels = levels(
    plot_genes$gene_label
   )
  ),
  cell_state_clean =
   clean_cell_state(
    cell_state
   )
 )

p_psoriasis_heat <- ggplot(
 gene_cell_heat,
 aes(
  x = cell_state_clean,
  y = gene_label
 )
) +
 geom_tile(
  aes(
   fill = mean_logcpm_diff_stim_minus_base
  ),
  color = "white"
 ) +
 scale_fill_gradient2(
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0,
  name = "Δ logCPM\nStim - Resting"
 ) +
 labs(
  x = "Cell state",
  y = NULL
 ) +
 theme_bw() +
 theme(
  axis.text.x = element_text(
   angle = 45,
   hjust = 1
  ),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
 )


# combined plot
p_psoriasis_combined <- (
 p_psoriasis_mr +
  p_psoriasis_heat +
  plot_layout(
   widths = c(2.5, 1.3)
  )
)

p_psoriasis_combined