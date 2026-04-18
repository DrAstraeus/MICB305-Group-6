---
title: "Association of Vitamin D Intake with Gut Microbiome Composition, Predicted Functional Potential and Disease Severity in Multiple Sclerosis Patients"
output: html_document
---

```{r}

###################### Alpha Diversity ######################

```{r}
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)
library(microbiome)
library(indicspecies)
library(writexl)
```

```{r}
metadata = read.delim("ms_metadata.tsv", row.names = 1, check.names = FALSE)
df = metadata %>% 
  rownames_to_column(var = "sample_id")
filtered_df = df %>%
  filter(sample_id != "X71802.0091") %>%
  filter(sample_id != "X71402.0259") %>%
  filter(sample_id != "X76402.0035") %>%
  filter(disease == "MS")
filtered_metadata = filtered_df %>%
  column_to_rownames(var = "sample_id")

taxonomy = read.delim("taxonomy.tsv", row.names = 1)

counts = read.delim("feature-table.txt", skip=1, row.names=1)

tree = read_tree("tree.nwk")
```

```{r}
taxonomy_formatted = taxonomy %>% 
  separate(col = Taxon, 
           into = c("Domain","Phylum","Class","Order",
                    "Family","Genus","Species"),
           sep=";", fill="right") %>% 
  select(-Confidence) %>% 
  as.matrix()

counts_formatted = counts %>% 
  as.matrix()
filtered_metadata = filtered_metadata %>%
  filter(!is.na(`vitamin D (IU)`)) %>%
  mutate(`vitaminD_status` = ifelse(`vitamin D (IU)` > 0, "Supplement", "No Supplement"))

ps <- phyloseq(sample_data(filtered_metadata), otu_table(counts_formatted, taxa_are_rows = T), tax_table(taxonomy_formatted), tree)

saveRDS(ps, "./phyloseq_taxonomy.rds")
```

```{r}
#Alpha Diversity Shannon
hist(sample_sums(ps)) 
rarecurve(t(data.frame(ps@otu_table)), step= 1000, label = FALSE)

set.seed(421)
psrare = ps %>% rarefy_even_depth(sample.size = 10378, rngseed = 421)
table(sample_sums(psrare))

comparisons = list(c("Supplement","No Supplement"))

p = plot_richness(psrare, x = "vitaminD_status", measures = "Shannon", color = "vitaminD_status")
pdata = p$data
sum(pdata$vitaminD_status == "Supplement")
sum(pdata$vitaminD_status == "No Supplement")

#Specific p-value
wtest = wilcox.test(value ~ vitaminD_status, data = pdata, conf.int = T)
wtest

#Making Plot
p_formatted_plot = pdata %>%
  ggplot(aes(vitaminD_status,value, fill = vitaminD_status)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height=0, width=0.2, alpha = 0.2) +
  theme_classic(base_size=10) +
  theme(axis.text.x = element_text(angle = 45,vjust=1, hjust=1, size = 12), text = element_text(size = 15)) +
  stat_compare_means(comparisons = comparisons, method="wilcox.test", label = "p.signif", size = 6) +
  labs(x = "Vitamin D Status", y = "Shannon Diversity")+
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  scale_fill_manual(name = "Vitamin D Status", values = c("Supplement" = "#ABDAF4", "No Supplement" = "#F3CF80")) 
p_formatted_plot
```

```{r}
set.seed(421)
table(sample_sums(psrare))

p2 = plot_richness(psrare, x = "vitaminD_status", 
                  measures = c("Observed","Shannon","Chao1","Simpson"),
                  color = 'vitaminD_status')
p2



pdata2 = p2$data

p_formatted_plot2 = pdata2 %>%
  ggplot(aes(vitaminD_status,value, fill = vitaminD_status)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(height=0, width=0.2, alpha = 0.2) +
  theme_classic(base_size=10) +
  facet_wrap(~variable,ncol=4,scales = 'free_y')+
  theme(axis.text.x = element_text(angle = 45,vjust=1, hjust=1, size = 12), text = element_text(size = 15)) +
  stat_compare_means(comparisons = comparisons, method="wilcox.test", label = "p.signif", size = 6) +
  labs(x = "Vitamin D Supplemented", y = "Alpha Diversity") +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  scale_fill_manual(name = "Vitamin D Status", values = c("Supplement" = "#ABDAF4", "No Supplement" = "#F3CF80")) 
p_formatted_plot2
```

###################### Beta Diversity ######################
```{r}
#Beta Diversity (Weighted Unifrac)

ps_bray = phyloseq::distance(psrare, method = "wunifrac")
View(as.matrix(ps_bray))

set.seed(421)
mds = metaMDS(ps_bray)

# Extract data
mds_data = mds$points %>% as.data.frame %>%
  merge(sample_data(psrare), by='row.names', sort=F)
head(mds_data)

p3 = mds_data %>%
  ggplot(aes(MDS1,MDS2, color = vitaminD_status)) +
  geom_point(alpha = 0.9) +
  stat_ellipse() +
  theme_classic() +
  scale_color_manual(name = "Vitamin D Status", values = c("Supplement" = "#ABDAF4", "No Supplement" = "#F3CF80")) +
  labs(color = "Vitamin D Status") 
p3

ps_metadata = sample_data(psrare) %>% data.frame() # to make the next lines easier

# Single variable
set.seed(421)
stats_univar = adonis2(ps_bray ~ vitaminD_status, data = ps_metadata)
stats_univar

stats = bind_rows('Univariate' = stats_univar %>% as.data.frame() %>% rownames_to_column('Variable'), .id = 'Model') 
stats
```


```{r}
#Beta Diversity (Bray Curtis)

ps_bray2 = phyloseq::distance(psrare, method = "bray")
View(as.matrix(ps_bray2))

set.seed(421)
mds2 = metaMDS(ps_bray2)

# Extract data
mds_data2 = mds2$points %>% as.data.frame %>%
  merge(sample_data(psrare), by='row.names', sort=F)
head(mds_data2)

p4 = mds_data2 %>%
  ggplot(aes(MDS1,MDS2, color = vitaminD_status)) +
  geom_point(alpha = 0.9) +
  stat_ellipse() +
  theme_classic() +
  scale_color_manual(name = "Vitamin D Status", values = c("Supplement" = "#ABDAF4", "No Supplement" = "#F3CF80")) +
  labs(color = "Vitamin D Status") 
p4

ps_metadata2 = sample_data(psrare) %>% data.frame()

# Single variable
set.seed(421)
stats_univar2 = adonis2(ps_bray2 ~ vitaminD_status, data = ps_metadata2)
stats_univar2

stats2 = bind_rows('Univariate' = stats_univar2 %>% as.data.frame() %>% 
                    rownames_to_column('Variable'),
                  .id = 'Model')
```

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
