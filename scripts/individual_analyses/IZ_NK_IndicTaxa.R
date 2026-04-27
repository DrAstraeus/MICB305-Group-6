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
