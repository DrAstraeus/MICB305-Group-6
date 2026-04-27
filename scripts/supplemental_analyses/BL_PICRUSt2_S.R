```{r}
                                   
library(tidyverse)
library(phyloseq)
library(ggpicrust2)

set.seed(421)

meta = readRDS('phyloseq_taxonomy.rds') %>%
  .@sam_data %>%
  data.frame() %>%
  rownames_to_column('sample_name') %>%
  filter(sample_name != c('X71802.0091','X71402.0259','X76402.0035')) %>%
  filter(disease == 'MS')

metacyc = read.delim('path_abun_unstrat.tsv')
metacyc_tidy = metacyc %>%
  select('pathway',all_of(meta$sample_name)) %>%
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
                     fc_cutoff = 0, order_by_log = F,
                     p_values_threshold = 0.05, 
                     order = "name", 
                     ko_to_kegg = FALSE, 
                     p_value_bar = TRUE, 
                     x_lab = "description")
peb

```
