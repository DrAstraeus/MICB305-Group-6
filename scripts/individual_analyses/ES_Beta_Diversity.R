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
