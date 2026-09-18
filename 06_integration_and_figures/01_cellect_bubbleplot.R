library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

# --- Path ---
cellect_data <- fread("path/to/cellect_results.txt")

# --- Plot ---
bubble_df <- cellect_data %>%
 mutate(
  # remove prefix
  Name_clean = str_remove(Name, "^LawlorPBMC__"),
  
  # extract stimulation from suffix
  stimulation = case_when(
   str_detect(Name_clean, "_Baseline$")  ~ "Resting",
   str_detect(Name_clean, "_LPS$")       ~ "LPS",
   str_detect(Name_clean, "_CD3_CD28$")  ~ "CD3/CD28",
   str_detect(Name_clean, "_CD23_CD28$") ~ "CD23/CD28",
   TRUE ~ "Other"
  ),
  
  # remove stimulation suffix to get cell type
  cell_type = Name_clean %>%
   str_remove("_(Baseline|LPS|CD3_CD28|CD23_CD28)$") %>%
   str_replace_all("_", " "),
  
  significant = Coefficient_P_value < 0.05
 ) %>%
 mutate(
  stimulation = factor(
   stimulation,
   levels = c("Resting", "LPS", "CD3/CD28", "CD23/CD28", "Other")
  )
 )

ggplot(bubble_df, aes(x = stimulation, y = cell_type)) +
 geom_point(
  aes(size = abs(Coefficient), color = significant),
  alpha = 0.8
 ) +
 scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#DA2C43")) +
 scale_size_continuous(name = "|Coefficient|") +
 labs(
  x = "Stimulation",
  y = "Cell type",
  color = "p < 0.05"
 ) +
 theme_bw() +
 theme(
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1)
 )