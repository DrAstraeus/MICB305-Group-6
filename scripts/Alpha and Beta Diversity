---
title: "Alpha and Beta Diversity"
output: html_document
date: "2026-03-17"
---

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
  filter(sample_id != "X76402.0035")
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

ggsave('./Alpha Diversity.jpeg', plot = p_formatted_plot, height= 8, width = 6)
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

```{r}
#Beta Diversity

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
  geom_point(alpha = 0.5) +
  stat_ellipse() +
  theme_minimal() +
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
#Beta Diversity

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
  geom_point(alpha = 0.5) +
  stat_ellipse() +
  theme_minimal() +
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


