---
title: "Association of Vitamin D Intake with Gut Microbiome Composition, Predicted Functional Potential and Disease Severity in Multiple Sclerosis Patients"
output: html_document
---

```{r}

###################### Alpha Diversity ######################

###################### Beta Diversity ######################

###################### Indicator Species Analysis ######################

###################### Differential Abundance Analysis ######################

###################### Functional Analysis ######################

library(tidyverse)
library(phyloseq)
library(ggpicrust2)

set.seed(421)

meta = readRDS('Datasets/phyloseq_taxonomy.rds') %>%
  .@sam_data |>
  data.frame() |>
  rownames_to_column('sample_name') |>
  filter(sample_name != c('X71802.0091','X71402.0259','X76402.0035')) |>
  filter(disease == 'MS')

metacyc = read.delim('Datasets/path_abun_unstrat.tsv')
metacyc_tidy = metacyc |>
  select('pathway',all_of(meta$sample_name)) |>
  column_to_rownames('pathway')

daa_results_df = pathway_daa(abundance = metacyc_tidy,
                             metadata = meta, 
                             group = "vitaminD_cat", 
                             daa_method = "Maaslin2", 
                             select = NULL, reference = NULL)

 daa_annotated_results_df = pathway_annotation(pathway = "MetaCyc",
                                               daa_results_df = daa_results_df)

source('ggpicrust2_errorbar_function_fixed.R')

peb = pathway_errorbar_fixed(abundance = metacyc_tidy, 
                     daa_results_df = daa_annotated_filt, 
                     Group = factor(meta$vitaminD_cat,
                                    levels = c('Supplement',
                                               'No Supplement')), 
                     wrap_label = T, wraplength=65,
                     fc_cutoff = 0.5, order_by_log = F,
                     p_values_threshold = 5e-4, 
                     order = "name", 
                     ko_to_kegg = FALSE, 
                     p_value_bar = TRUE, 
                     x_lab = "description")
peb

###################### Correlation Matrix ######################

library(tidyverse)
library(phyloseq)
library(pheatmap)

ps = readRDS('Datasets/phyloseq_taxonomy.rds') |> 
  tax_glom('Genus') |>
  subset_samples(disease == 'MS')

metacyc = read_tsv('Datasets/metacyc.tsv')

ps_rel = ps |>
  microbiome::transform('compositional')

# taxa

taxa = as.data.frame(tax_table(ps_rel)) |>
  filter(Genus == " g__Succinivibrio" |
           Genus == " g__Frisingicoccus") |>
  rownames()

ps_filt = prune_taxa(taxa,ps_rel)
otu = data.frame(ps_filt@otu_table)

# edss

meta = data.frame(sample_data(ps_filt)) |>
  select(edss) |>
  t()

# pathways

metacyc_select = metacyc |>
  filter(pathway == 'PWY-6263' |
           pathway == 'PWY-7371' |
           pathway == 'PWY-7374' |
           pathway == 'PWY66-409') |>
  column_to_rownames('pathway')

# merge

summary = rbind(otu, `EDSS` = meta, 
                metacyc_select)
rownames(summary) = c('Succinivibrio',
                      'Frisingicoccus',
                      'EDSS',
                      'superpathway of menaquinol-8 biosynthesis II',
                      '1,4-dihydroxy-6-naphthoate biosynthesis II',
                      '1,4-dihydroxy-6-naphthoate biosynthesis I',
                      'superpathway of purine nucleotide salvage')

# corr calculation

otu_cor = cor(t(summary), method = 'spearman')

corr_mat = pheatmap(otu_cor,
         clustering_method = "complete",
         color = colorRampPalette(c("#F2CF7F","white","#AAD9F4"))(50),
         breaks = seq(-1, 1, length.out = 51),
         main = "",
         display_numbers = T,
         angle_col = "315",
         fontsize_col = 7.5,
         fontsize_row = 8.5)


```
