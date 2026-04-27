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
