```{r}
library(tidyverse)
library(phyloseq)
library(pheatmap)

ps = readRDS('phyloseq_taxonomy.rds') %>% 
  tax_glom('Genus') %>%
  subset_samples(disease == 'MS')

metacyc = read_tsv('metacyc.tsv')

ps_rel = ps %>%
  microbiome::transform('compositional')

# taxa

taxa = as.data.frame(tax_table(ps_rel)) %>%
  filter(Genus == " g__Succinivibrio" |
           Genus == " g__Frisingicoccus") %>%
  rownames()

ps_filt = prune_taxa(taxa,ps_rel)
otu = data.frame(ps_filt@otu_table)

# edss

meta = data.frame(sample_data(ps_filt)) %>%
  select(edss) %>%
  t()

# pathways

metacyc_select = metacyc %>%
  filter(pathway == 'PWY-6263' |
           pathway == 'PWY-7371' |
           pathway == 'PWY-7374' |
           pathway == 'PWY66-409') %>%
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
