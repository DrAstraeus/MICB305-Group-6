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
