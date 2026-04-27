---
title: "Vitamin D supplementation is correlated with gut microbial diversity, abundance of specific taxa, and predicted functional potential of metabolic pathways in patients with multiple sclerosis"
output: html_document
---

```{r}

########## TABLE OF CONTENTS ##########
# 1. Alpha Diversity                  #
# 2. Beta Diversity                   #
# 3. Indicator Taxa Analysis          #
# 4. Differential Abundance Analysis  #
# 5. Functional Analysis              #
# 6. Correlation Matrix               #
# 7. Supplemental Figures             #
#######################################
```
################################################################## 1. Alpha Diversity ##################################################################

```{r}
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)

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

################################################################## 2. Beta Diversity ##################################################################
```{r}

# library
library(tidyverse)
library(phyloseq)
library(vegan)

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

################################################################## 3. Indicator Taxa Analysis ##################################################################
```{r}
library(tidyverse)
library(phyloseq)
library(indicspecies)
library(dplyr)
library(tibble)
library(ggplot2)

# Load Taxonomy
taxonomy <- read.delim("taxonomy.tsv", header = T, row.names = 1)

# Load Counts
counts_raw <- read.delim("feature-table.txt", skip = 1)
counts_formatted <- as.matrix(counts_raw)

# Make first column row names
rownames(counts_formatted) <- counts_formatted[, 1]

# Remove first column, only numbers left
counts_numeric_only <- counts_formatted[, -1]

counts_final <- matrix(as.numeric(counts_numeric_only), 
                       nrow = nrow(counts_numeric_only), 
                       ncol = ncol(counts_numeric_only))

# Rename rows and columns
rownames(counts_final) <- rownames(counts_numeric_only)
colnames(counts_final) <- colnames(counts_numeric_only)

# Load Tree
tree <- read_tree("tree.nwk")

# Taxonomy Wrangle
taxonomy_formatted = taxonomy %>% 
  rownames_to_column(var = "FeatureID") %>%
  separate(col = Taxon, 
           into = c('Domain','Phylum','Class','Order', 'Family','Genus','Species'),
           sep=';', fill='right') %>% 
  select(-Confidence) %>% 
  column_to_rownames(var = "FeatureID") %>%
  as.matrix()

# Metadata Wrangle
meta_names = read.delim("ms_metadata.tsv", row.names = 1) %>% names()
meta_data = read.delim("ms_metadata.tsv", row.names = 1, skip=2, header = F) 
names(meta_data) = meta_names

meta_data = meta_data %>%
  rownames_to_column(var = "sample_id") %>%
  filter(sample_id != "X71802.0091") %>%
  filter(sample_id != "X71402.0259") %>%
  filter(sample_id != "X76402.0035") %>%
  column_to_rownames(var = "sample_id")

# Phyloseq Object
ps = phyloseq(sample_data(meta_data),
              otu_table(counts_final, taxa_are_rows = T),
              tax_table(taxonomy_formatted),
              tree)

ps

### Genus-level plot.
#ps_genus = tax_glom(ps,'Genus')

## To make quicker
# Rebuild object without tree to speed up tax_glom
ps_notree <- phyloseq(
  otu_table(ps),
  sample_data(ps),
  tax_table(ps)
)

# group Genus-level
ps_genus = tax_glom(ps_notree, 'Genus')

# Transform to relative abundance
ps_relab = microbiome::transform(ps_genus, 'compositional')

# Filter out low-abundance taxa
ps_filt = filter_taxa(ps_relab, function(x) mean(x) > 0.001, TRUE)


# "Supplement" vs "No_Supplement"
sample_data(ps_filt)$vitamin.D..IU. <- ifelse(sample_data(ps_filt)$`vitamin.D..IU.` > 0, 
                                              "Supplement", "No_Supplement")

# Convert to factor
sample_data(ps_filt)$vitamin.D..IU. <- as.factor(sample_data(ps_filt)$vitamin.D..IU.)

# Remove NA
ps_filt <- subset_samples(ps_filt, !is.na(vitamin.D..IU.))

#Takes a VERY long time (10-20 mins ish?) bc I did nperm = 99,999
set.seed(421)
indval = multipatt(t(otu_table(ps_filt)), 
                   cluster = sample_data(ps_filt)$vitamin.D..IU.,
                   control = how(nperm = 99999))

# Convert to factor
sample_data(ps_filt)$vitamin.D..IU. <- as.factor(sample_data(ps_filt)$vitamin.D..IU.)

# Remove NA
ps_filt <- subset_samples(ps_filt, !is.na(vitamin.D..IU.))

#Takes a VERY long time (10-20 mins ish?) bc I did nperm = 99,999
set.seed(421)
indval = multipatt(t(otu_table(ps_filt)), 
                   cluster = sample_data(ps_filt)$vitamin.D..IU.,
                   control = how(nperm = 99999))

summary(indval, indvalcomp = TRUE)

indval_table = as.data.frame(indval$sign)

####Graph
phyla_to_plot = indval_table %>% 
  filter(s.No_Supplement == 1 & s.Supplement == 0) %>% 
  rownames()

df_of_taxa = prune_taxa(phyla_to_plot, ps_filt) %>% psmelt()

df_of_taxa %>% 
  ggplot(aes(vitamin.D..IU., Abundance, fill = vitamin.D..IU.)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 1, width = 0.2) +
  facet_wrap(~Genus, ncol = 3, scales = 'free')

# Use Relative Abundance to prevent negative values
ps_plot <- transform_sample_counts(ps_filt, function(x) x / sum(x))

# Only significant ASVs
significant_asvs_df <- indval_table %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV_ID") %>%
  filter(p.value < 0.05)

significant_asvs <- rownames(indval_table)[which(indval_table$p.value < 0.05)]
# This also filters out taxa only present in No supplements (only 1 taxa)

# Prune and Melt
df_to_plot <- ps_plot %>%
  prune_taxa(significant_asvs, .) %>%
  psmelt()

# Plot
ggplot(df_to_plot, aes(x = vitamin.D..IU., y = Abundance, fill = vitamin.D..IU.)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = 1.5) + 
  facet_wrap(~Genus, scales = "free_y") + 
  theme_bw() +
  labs(title = "Significant Indicator Taxa for Different Vitamin D Status",
       y = "Relative Abundance")


#### Make plot publication ready


# Remove g__
df_to_plot$Genus <- gsub("g__", "", df_to_plot$Genus)

# Plot with Log10 transformation to expand the bottom (change this code with what Avril suggested in the meeting notes)
ggplot(df_to_plot, aes(x = vitamin.D..IU., y = Abundance, fill = vitamin.D..IU.)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.4,
    alpha = 0.5,
    color = "black"
  ) +
  geom_jitter(
    height = 0,
    width = 0.1,
    alpha = 0.4,
    size = 1.5
  ) +
  scale_y_log10() +
  scale_x_discrete(
    labels = c(
      "No_Supplement" = "No Supplement",
      "Supplement" = "Supplement"
    )
  ) +
  scale_fill_manual(
    values = c(
      "No_Supplement" = "#E69F00",
      "Supplement" = "#56B4E9"
    ),
    labels = c(
      "No_Supplement" = "No supplement",
      "Supplement" = "Supplement"
    ),
    name = "Vitamin D Status"
  ) +
  geom_signif(
    comparisons = list(c("No_Supplement", "Supplement")),
    annotation = "*",
    vjust = 0.5,
    y_position = 0.6
  ) +
  facet_wrap(~Genus) +
  labs(
    x = "Vitamin D Status",
    y = "Relative Abundance"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
    strip.text = element_text(face = "italic", size = 14),
    strip.background = element_rect(fill = "gray95"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 11),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

####Indicator sig taxa plot

#Build a clean table
indval_df <- indval$sign %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV_ID") %>%
  filter(p.value < 0.05)

#Attach taxonomy (get Genus names)
tax_df <- as.data.frame(tax_table(ps_genus)) %>%
  rownames_to_column(var = "ASV_ID")

indval_df <- indval_df %>%
  left_join(tax_df, by = "ASV_ID")

indval_df$Genus <- gsub("g__", "", indval_df$Genus)

#Assign which group each taxon indicates
indval_df <- indval_df %>%
  mutate(
    Group = case_when(
      s.No_Supplement == 1 ~ "No supplement",
      s.Supplement == 1 ~ "Supplement"
    )
  )

#Plot IndVal values

ggplot(indval_df, aes(x = reorder(Genus, stat), y = stat, fill = Group)) +
  geom_col(width = 0.7, color = "black") +
  
  coord_flip() +
  
  scale_fill_manual(
    values = c("No supplement" = "#E69F00",
               "Supplement" = "#56B4E9")
  ) +
  
  labs(
    x = "Genus",
    y = "Indicator Value (IndVal)",
    fill = "Vitamin D Status"
  ) +
  
  theme_bw() +
  theme(
    axis.text.y = element_text(face = "italic"),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )
```
################################################################## 4. Differential Abundance Analysis ##################################################################
```{r}
#Load library
 # if(!requireNamespace("BiocManager", quietly = TRUE))
 #     install.packages("BiocManager")
 # BiocManager::install("Maaslin2")
library(Maaslin2)
library(tidyverse)
library(phyloseq)
library(ggplot2)

#clean data
clean_ms_metadata = read.delim('ms_metadata.tsv',row.names = 1) %>%
  filter(!is.na(vitamin.D..IU.)) %>%
  mutate(
    supp_status=factor(if_else (vitamin.D..IU.>0, "supplement", "no supplement"),
  levels = c("no supplement", "supplement")))

#save into table
write.table(clean_ms_metadata, 
            file = "clean_ms_metadata.tsv", 
            sep = "\t", 
            row.names = TRUE, 
            quote = FALSE)

# Load datasets
taxonomy = read.delim('taxonomy.tsv', row.names = 1)
tree = read_tree('tree.nwk')

counts = read.delim('feature-table.txt', skip=1, row.names=1) # First line is not data
metadata = read.delim('clean_ms_metadata.tsv', row.names = 1) # 1st col are names

# Wrangle Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Taxonomy
taxonomy_formatted = taxonomy %>% 
  separate(col = Taxon, 
           into = c('Domain','Phylum','Class','Order',
                    'Family','Genus','Species'),
           sep=';', fill='right') %>% 
  select(-Confidence) %>% 
  as.matrix()

# Counts
counts_formatted = counts %>% as.matrix()

# Metadata ~~~~~~~
View(metadata)
# Extract just the column names
meta_names = read.delim('clean_ms_metadata.tsv',row.names = 1) %>% names()

meta_data = read.delim('clean_ms_metadata.tsv',row.names = 1,skip=2,header = F) 
table(metadata[2,]== meta_data[1,]) 

# Set column names for metadata
names(meta_data) = meta_names

# Create the phyloseq object
ps = phyloseq(sample_data(meta_data),
              otu_table(counts_formatted, taxa_are_rows = T),
              tax_table(taxonomy_formatted),
              tree)

# Save as .rds or .Rdata object
saveRDS(ps,'phyloseq_taxonomy.rds')




# Load object, filter data in phyloseq object
ps = readRDS('phyloseq_taxonomy.rds')%>%
  subset_samples( sample_names(ps) != 'X71802.0091' & 
                   sample_names(ps) != 'X71402.0259' & 
                   sample_names(ps) != 'X76402.0035')%>%
  subset_samples(disease=="MS")


ps_glom = tax_glom(ps, 'Genus')

# Differential abundance with MaAsLin2

set.seed(421)
out = Maaslin2(
  input_data = data.frame(ps_glom@otu_table), 
  input_metadata = data.frame(ps_glom@sam_data), 
  output = 'to_delete', 
  fixed_effects = c('supp_status'),
  normalization = 'TSS',
  transform = 'AST',
  min_abundance = 0.001, 
  min_prevalence = 0.1,
  max_significance = 0.05,
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

statistical_table = out$results


# Filter statistical table
taxa_to_plot = statistical_table %>% 
  filter(qval < 0.05)
taxa_to_plot=taxa_to_plot %>%
  mutate(feature=str_replace(feature,"^X(?=[0-9])", "")) %>%
  mutate(feature=str_replace_all(feature, "\\.", "-"))


#calculate log 2 fold change

ps_relab = transform_sample_counts(ps_glom, function(x) x / sum(x))

ps_melt = psmelt(ps_relab) %>%
  filter(!is.na(Abundance)) %>%
  mutate(Genus=str_replace(Genus,"g__", ""))
head(ps_melt)

pseudocount_value = min(ps_melt$Abundance[ps_melt$Abundance > 0], na.rm=TRUE) / 2
pseudocount_value

ps_melt$Abundance =ps_melt$Abundance + pseudocount_value
min(ps_melt$Abundance)

avg_abundances = ps_melt %>%
  inner_join(taxa_to_plot, by = c("OTU"="feature")) %>%
  group_by(Genus, supp_status) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ungroup()

ratio = avg_abundances %>%
  pivot_wider(names_from = supp_status, values_from = Abundance) %>%
  mutate(ratio = supplement / `no supplement`) %>% 
  mutate(log2fc = log2(ratio))

head(ratio)

ratio_filt = ratio %>%
  filter(abs(log2fc)>=1)

ratio_filt



#plot taxa
lfc_plot= ratio_filt %>%
  ggplot(aes(x=log2fc,y=Genus, fill=log2fc>0)) +
  geom_col(width=0.4)+
  scale_fill_manual(values=c("TRUE"="#80B1D3","FALSE"="#FDB462"),
                    labels=c("TRUE"="Enriched with supplementation","FALSE"= "Depleted with supplementation"),
                    name="Abundance Trend")+
  theme(
  axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(size=6),
  axis.text.x = element_text(size = 6),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6))+
  labs(x= "Log2 Fold Change Abundance",
    y= "Genera")

lfc_plot

#save plot
ggsave("differential_abundance_plot.png", 
       plot = lfc_plot,    
       width = 6,             
       height = 2.5,           
       dpi = 300)
```
                                   
################################################################## 5. Functional Analysis ##################################################################
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

source('ggpicrust2_calculate_log2fc.R')                                
lfc_fixed = fix_lfc(abundance = metacyc_tidy,
                             metadata = meta, 
                             group = "vitaminD_cat",
                             reference = 'No Supplement')
                                   
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
                     p_values_threshold = 0.05, 
                     order = "name", 
                     ko_to_kegg = FALSE, 
                     p_value_bar = TRUE, 
                     x_lab = "description")
peb

```                                   
################################################################## 6. Correlation Matrix ##################################################################

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
                                   
################################################################## 7. Supplemental Figures ##################################################################


                                   
####### Alpha Diversity #######
```{r}

library(tidyverse)
library(phyloseq)
library(vegan)
library(ggpubr)

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

hist(sample_sums(ps)) 
rarecurve(t(data.frame(ps@otu_table)), step= 1000, label = FALSE)

set.seed(421)
psrare = ps %>% rarefy_even_depth(sample.size = 10378, rngseed = 421)
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

####### Differential Abundance ####### 
```{r}
library(Maaslin2)
library(tidyverse)
library(phyloseq)
library(ggplot2)

#clean data
clean_ms_metadata = read.delim('Datasets/ms_metadata.tsv',row.names = 1)%>%
  filter(!is.na(vitamin.D..IU.))%>%
  mutate(
    supp_status=factor(if_else (vitamin.D..IU.>0, "supplement", "no supplement"),
  levels = c("no supplement", "supplement")))

#save into table
write.table(clean_ms_metadata, 
            file = "Datasets/clean_ms_metadata.tsv", 
            sep = "\t", 
            row.names = TRUE, 
            quote = FALSE)

# Load datasets
taxonomy = read.delim('Datasets/taxonomy.tsv', row.names = 1)
tree = read_tree('Datasets/tree.nwk')

counts = read.delim('Datasets/feature-table.txt', skip=1, row.names=1) # First line is not data
metadata = read.delim('Datasets/clean_ms_metadata.tsv', row.names = 1) # 1st col are names

# Wrangle Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Taxonomy
taxonomy_formatted = taxonomy %>% 
  separate(col = Taxon, 
           into = c('Domain','Phylum','Class','Order',
                    'Family','Genus','Species'),
           sep=';', fill='right') %>% 
  select(-Confidence) %>% 
  as.matrix()

# Counts
counts_formatted = counts %>% as.matrix()

# Metadata ~~~~~~~
View(metadata)
# Extract just the column names
meta_names = read.delim('Datasets/clean_ms_metadata.tsv',row.names = 1) %>% names()

meta_data = read.delim('Datasets/clean_ms_metadata.tsv',row.names = 1,skip=2,header = F) 
table(metadata[2,]== meta_data[1,]) 

# Set column names for metadata
names(meta_data) = meta_names

# Create the phyloseq object
ps = phyloseq(sample_data(meta_data),
              otu_table(counts_formatted, taxa_are_rows = T),
              tax_table(taxonomy_formatted),
              tree)

# Save as .rds or .Rdata object
saveRDS(ps,'Datasets/phyloseq_taxonomy.rds')




# Load object, filter data in phyloseq object
ps = readRDS('Datasets/phyloseq_taxonomy.rds')%>%
  subset_samples( sample_names(ps) != 'X71802.0091' & 
                   sample_names(ps) != 'X71402.0259' & 
                   sample_names(ps) != 'X76402.0035')%>%
  subset_samples(disease=="MS")


ps_glom = tax_glom(ps, 'Genus')

# Differential abundance with MaAsLin2

set.seed(421)
out = Maaslin2(
  input_data = data.frame(ps_glom@otu_table), 
  input_metadata = data.frame(ps_glom@sam_data), 
  output = 'to_delete', 
  fixed_effects = c('supp_status'),
  normalization = 'TSS',
  transform = 'AST',
  min_abundance = 0.001, 
  min_prevalence = 0.1,
  max_significance = 0.05,
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

statistical_table = out$results

writexl::write_xlsx(statistical_table,'Maaslin2 Results.xlsx')


# Filter statistical table
taxa_to_plot = statistical_table %>% 
  filter(qval < 0.05)
taxa_to_plot=taxa_to_plot%>%
  mutate(feature=str_replace(feature,"^X(?=[0-9])", ""))%>%
  mutate(feature=str_replace_all(feature, "\\.", "-"))


#calculate log 2 fold change

ps_relab = transform_sample_counts(ps_glom, function(x) x / sum(x))

ps_melt = psmelt(ps_relab)%>%
  filter(!is.na(Abundance))%>%
  mutate(Genus=str_replace(Genus,"g__", ""))
head(ps_melt)

pseudocount_value = min(ps_melt$Abundance[ps_melt$Abundance > 0], na.rm=TRUE) / 2
pseudocount_value

ps_melt$Abundance =ps_melt$Abundance + pseudocount_value
min(ps_melt$Abundance)

avg_abundances = ps_melt %>%
  inner_join(taxa_to_plot, by = c("OTU"="feature"))%>%
  group_by(Genus, supp_status) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE)) %>% 
  ungroup()

ratio = avg_abundances %>%
  pivot_wider(names_from = supp_status, values_from = Abundance) %>%
  mutate(ratio = supplement / `no supplement`) %>% 
  mutate(log2fc = log2(ratio))

head(ratio)

ratio_filt = ratio %>%
  filter(abs(log2fc)>=1)

ratio_filt

                                   
lfc_plot= ratio%>%
  ggplot(aes(reorder(Genus, -log2fc), log2fc, fill=log2fc>0)) +
  geom_col(width=0.8, linewidth = 0.2, color="black", show.legend=FALSE)+
  scale_fill_manual(values=c("TRUE"="#AAD9F4","FALSE"="#F2CF7F"),
                    labels=c("TRUE"="Enriched with supplementation","FALSE"= "Depleted with supplementation"),
                    name="Abundance Trend")+
  theme(
  axis.text.y = element_text(size = 6),
  axis.title.y = element_text(size=8),
  axis.title.x = element_text(size=6),
  axis.text.x = element_text(size = 6),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6))+
     coord_flip() +
  labs(y= "Log2 Fold Change Abundance",
    x= "Genera")

lfc_plot                                   
```
####### Functional Analysis ####### 

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

####### Correlation Matrix ####### 
```{r}
library(tidyverse)
library(phyloseq)
library(pheatmap)
library(igraph)

ps = readRDS('Datasets/phyloseq_taxonomy.rds') %>% 
  tax_glom('Genus') %>%
  subset_samples(disease == 'MS')

metacyc = read_tsv('Datasets/metacyc.tsv')

ps_rel = ps %>%
  microbiome::transform('compositional')

# taxa

taxa = as.data.frame(tax_table(ps_rel)) %>%
  filter(Genus == " g__Succinivibrio" |
           Genus == " g__Frisingicoccus" |
           Genus == " g__Erysipelatoclostridium" |
           Genus == " g__Desulfovibrio") %>%
  rownames()

ps_filt = prune_taxa(taxa,ps_rel)
otu = data.frame(ps_filt@otu_table)

# edss

meta = data.frame(sample_data(ps_filt)) %>%
  select(edss, vitamin.D..IU.) %>%
  t()

# pathways

metacyc_select = metacyc %>%
  filter(pathway == 'PWY-6263' |
           pathway == 'PWY-7371' |
           pathway == 'PWY-7374' |
           pathway == 'PWY66-409') %>%
  column_to_rownames('pathway')

# merge

summary = rbind(otu, 
                meta, 
              metacyc_select)
rownames(summary) = c('Desulfovibrio',
                      'Succinivibrio',
                      'Erysipelatoclostridium',
                      'Frisingicoccus',
                      'EDSS',
                      'Vitamin D (IU)',
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
####### Vitamin D Intake, Zhou et al. ####### 
```{r}
library(tidyverse)
library(ggsignif)
graph = metadata %>%
  select(disease,`vitamin D (IU)`) %>%
  filter(!is.na(`vitamin D (IU)`))
sum(graph$disease == "MS")
sum(graph$disease == "Control")

MSvitD_summ = graph %>%
  ggplot(aes(x=disease,y=`vitamin D (IU)`, fill = disease)) +
   geom_violin() +
   geom_boxplot(
    outlier.shape = NA,
    width = 0.4,
    alpha = 0.5,
    color = "black") +
  geom_jitter(
    height = 0,
    width = 0.1,
    alpha = 0.4,
    size = 1.5) +
  scale_x_discrete(
    labels = c(
      "Control" = "HHC",
      "MS" = "MS")) +
  scale_fill_manual(
    values = c(
      "Control" = "#E69F00",
      "MS" = "#56B4E9"),
    labels = c(
      "Control" = "HHC",
      "MS" = "MS"),
    name = "Disease") +
  geom_signif(
    comparisons = list(c("Control", "MS")),
    annotation = "***",
    vjust = 0.5) +
  labs(x = NULL, y = 'Vitamin D (IU)') +
theme_bw()
MSvitD_summ                                   
```
