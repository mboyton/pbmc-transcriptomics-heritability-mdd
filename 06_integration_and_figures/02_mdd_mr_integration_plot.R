library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

# paths
mr_file_path <- "path/to/mdd_mr_results.txt.gz"
nk_rest_file_path <- "path/to/nk_resting_markers.rda"
nk_lps_wgcna_path <- "path/to/nk_lps_wgcna_results.rda"

# settings
mr_method <- "Inverse variance weighted"
mr_pval_cutoff <- 0.05

kme_threshold <- 0.75
module_trait_fdr_cutoff <- 0.05

nk_rest_fdr_cutoff <- 0.05
nk_rest_logfc_cutoff <- 1
nk_rest_logcpm_cutoff <- 5

top_n_genes <- 30


# min-max scale against full reference distribution
minmax_against_reference <- function(x, ref) {
 
 ref_min <- min(ref, na.rm = TRUE)
 ref_max <- max(ref, na.rm = TRUE)
 
 (x - ref_min) / (ref_max - ref_min)
}


# load MDD MR results
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


# resting NK programme
load(nk_rest_file_path)

nk_rest_program <- rest_vs_other_markers %>%
 mutate(
  logFC = as.numeric(logFC),
  logCPM = as.numeric(logCPM),
  FDR = as.numeric(FDR)
 ) %>%
 filter(
  !is.na(gene),
  !is.na(logFC),
  !is.na(logCPM),
  !is.na(FDR),
  FDR < nk_rest_fdr_cutoff,
  logFC > nk_rest_logfc_cutoff,
  logCPM > nk_rest_logcpm_cutoff
 )

rest_scaling_reference <- nk_rest_program$logFC

nk_rest_program <- nk_rest_program %>%
 mutate(
  relative_nk_enrichment =
   minmax_against_reference(
    logFC,
    rest_scaling_reference
   )
 )


# integrate resting NK programme with MR
nk_rest_mr_all <- mr_ivw %>%
 inner_join(
  nk_rest_program %>%
   select(
    gene,
    logFC,
    logCPM,
    FDR,
    relative_nk_enrichment
   ),
  by = "gene"
 ) %>%
 mutate(
  pval_fdr_within_program =
   p.adjust(
    pval,
    method = "fdr"
   )
 )

nk_rest_mr_sig <- nk_rest_mr_all %>%
 filter(
  pval < mr_pval_cutoff
 )

top_rest_genes <- nk_rest_mr_sig %>%
 arrange(
  desc(abs_z)
 ) %>%
 group_by(
  gene,
  gene_symbol
 ) %>%
 slice(1) %>%
 ungroup() %>%
 head(top_n_genes) %>%
 arrange(OR) %>%
 mutate(
  gene_label = factor(
   gene_symbol,
   levels = gene_symbol
  )
 )


# NK LPS WGCNA programme
load(nk_lps_wgcna_path)

core_wgcna_genes <- gene_module_table %>%
 mutate(
  kME = as.numeric(kME),
  module_trait_cor = as.numeric(module_trait_cor),
  module_trait_fdr = as.numeric(module_trait_fdr),
  mean_logcpm_diff_stim_minus_base =
   as.numeric(mean_logcpm_diff_stim_minus_base)
 ) %>%
 filter(
  !is.na(gene),
  !is.na(kME),
  kME >= kme_threshold,
  !is.na(module_trait_fdr),
  module_trait_fdr < module_trait_fdr_cutoff
 )


# integrate LPS programme with MR
nk_lps_mr_all <- mr_ivw %>%
 inner_join(
  core_wgcna_genes %>%
   select(
    gene,
    module,
    kME,
    module_trait_cor,
    module_trait_fdr,
    module_trait_direction,
    mean_logcpm_diff_stim_minus_base
   ),
  by = "gene"
 ) %>%
 mutate(
  stim_expression_change =
   mean_logcpm_diff_stim_minus_base,
  
  pval_fdr_within_program =
   p.adjust(
    pval,
    method = "fdr"
   )
 )

nk_lps_mr_sig <- nk_lps_mr_all %>%
 filter(
  pval < mr_pval_cutoff
 )

top_lps_genes <- nk_lps_mr_sig %>%
 arrange(
  desc(abs_z)
 ) %>%
 group_by(
  gene,
  gene_symbol
 ) %>%
 slice(1) %>%
 ungroup() %>%
 head(top_n_genes) %>%
 arrange(OR) %>%
 mutate(
  gene_label = factor(
   gene_symbol,
   levels = gene_symbol
  )
 )


# common MR OR axis
common_or_range <- range(
 c(
  top_rest_genes$OR_lower,
  top_rest_genes$OR_upper,
  top_lps_genes$OR_lower,
  top_lps_genes$OR_upper
 ),
 na.rm = TRUE
)

or_padding <- diff(common_or_range) * 0.05

common_or_limits <- c(
 common_or_range[1] - or_padding,
 common_or_range[2] + or_padding
)


# resting NK MR forest plot
p_mdd_rest_mr <- ggplot(
 top_rest_genes,
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
 scale_y_continuous(
  limits = common_or_limits
 ) +
 coord_flip() +
 labs(
  x = NULL,
  y = "Odds ratio",
  title = "MDD MR results in the resting NK programme",
  subtitle = paste0(
   "Nominal MR p < ",
   mr_pval_cutoff
  )
 ) +
 theme_bw() +
 theme(
  plot.title = element_text(
   face = "bold"
  )
 )


# resting NK relative enrichment heatmap
rest_heat_data <- top_rest_genes %>%
 select(
  gene,
  gene_symbol,
  logFC,
  relative_nk_enrichment
 ) %>%
 mutate(
  gene_label = factor(
   gene_symbol,
   levels = levels(
    top_rest_genes$gene_label
   )
  ),
  transcriptional_metric = "NK resting"
 )

p_mdd_rest_heat <- ggplot(
 rest_heat_data,
 aes(
  x = transcriptional_metric,
  y = gene_label
 )
) +
 geom_tile(
  aes(
   fill = relative_nk_enrichment
  ),
  color = "white"
 ) +
 scale_fill_gradient(
  low = "white",
  high = "red",
  limits = c(0, 1),
  breaks = c(
   0,
   0.25,
   0.50,
   0.75,
   1
  ),
  name = "Relative NK\nenrichment"
 ) +
 labs(
  x = NULL,
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

p_mdd_resting <- (
 p_mdd_rest_mr +
  p_mdd_rest_heat +
  plot_layout(
   widths = c(2.5, 1)
  )
)


# NK LPS MR forest plot
p_mdd_lps_mr <- ggplot(
 top_lps_genes,
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
 scale_y_continuous(
  limits = common_or_limits
 ) +
 coord_flip() +
 labs(
  x = NULL,
  y = "Odds ratio",
  title = "MDD MR results in the LPS-responsive NK programme",
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


# NK LPS raw expression-change heatmap
lps_heat_data <- top_lps_genes %>%
 select(
  gene,
  gene_symbol,
  stim_expression_change
 ) %>%
 mutate(
  gene_label = factor(
   gene_symbol,
   levels = levels(
    top_lps_genes$gene_label
   )
  ),
  transcriptional_metric = "NK LPS"
 )

p_mdd_lps_heat <- ggplot(
 lps_heat_data,
 aes(
  x = transcriptional_metric,
  y = gene_label
 )
) +
 geom_tile(
  aes(
   fill = stim_expression_change
  ),
  color = "white"
 ) +
 scale_fill_gradient2(
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0,
  name = expression(
   Delta ~ "logCPM" ~
    "(LPS - resting)"
  )
 ) +
 labs(
  x = NULL,
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

p_mdd_lps <- (
 p_mdd_lps_mr +
  p_mdd_lps_heat +
  plot_layout(
   widths = c(2.5, 1)
  )
)


# combined plot
p_mdd_combined <- (
 p_mdd_resting /
  p_mdd_lps
)

p_mdd_combined