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
