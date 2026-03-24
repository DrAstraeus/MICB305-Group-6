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
library(igraph)

ps = readRDS('Datasets/phyloseq_taxonomy.rds') |> 
  tax_glom('Genus')

maaslin2 = read_tsv('Datasets/maaslin2.tsv')
metacyc = read_tsv('Datasets/metacyc.tsv')

ps_rel = ps |>
  microbiome::transform('compositional')

# taxa

taxa = as.data.frame(tax_table(ps_rel)) |>
  filter(Genus == " g__Succinivibrio" |
           Genus == " g__Frisingicoccus") |>
  rownames()

ps_filt = prune_taxa(indic,ps_rel)
otu = data.frame(ps_filt@otu_table)

# edss

meta = data.frame(sample_data(ps_filt)) |>
  rownames_to_column('sample.id')

match_edss = meta$edss[match(colnames(data.frame(ps_filt@otu_table)), meta$sample.id)]

# pathways

maaslin2 = maaslin2 |>
  mutate(pathway = feature) |>
  select(-feature)
metacyc = metacyc |>
  filter(pathway == 'PWY-6263' |
           pathway == 'PWY-7371' |
           pathway == 'PWY-7374' |
           pathway == 'PWY66-409')
pathway = metacyc |>
  column_to_rownames('pathway')
tpathway = as.data.frame(t(pathway)) |>
  rownames_to_column('sample.id')
match_pathway6263 = tpathway[match(colnames(otu), tpathway$sample.id), 'PWY-6263']
match_pathway7371 = tpathway[match(colnames(otu), tpathway$sample.id), 'PWY-7371']
match_pathway7374 = tpathway[match(colnames(otu), tpathway$sample.id), 'PWY-7374']
match_pathway409 = tpathway[match(colnames(otu), tpathway$sample.id), 'PWY66-409']

# merge

summary = rbind(otu, EDSS = match_edss, 
                   `superpathway of menaquinol-8 biosynthesis II` = match_pathway6263,
                   `1,4-dihydroxy-6-naphthoate biosynthesis II` = match_pathway7371,
                   `1,4-dihydroxy-6-naphthoate biosynthesis I` = match_pathway7374,
                   `superpathway of purine nucleotide salvage` = match_pathway409)
  t() |>
  as.data.frame() |>
  drop_na() |>
  t() |>
  as.data.frame()

# corr calculation

otu_cor = cor(t(summary), method = 'spearman')
cor_cutoff = 0.3

wmat = abs(otu_cor)
diag(wmat) = 0
wmat[wmat < cor_cutoff] = 0
g_w = graph_from_adjacency_matrix(wmat, mode = "undirected", weighted = TRUE, diag = FALSE)
E(g_w)$weight = 1 / pmax(E(g_w)$weight, .Machine$double.eps)
btw_w = betweenness(g_w, weights = E(g_w)$weight, normalized = TRUE)
btw_w %>% sort

corr_mat = pheatmap(otu_cor,
         clustering_method = "complete",
         color = colorRampPalette(c("blue","white","red"))(50),
         breaks = seq(-1, 1, length.out = 51),
         main = "Spearman's Correlation Heatmap")


```
