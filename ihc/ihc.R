library(tidyverse)
library(ggprism)
library(svglite)

#---- EGR1
egr1.poor.good <- read.csv("./egr1poorgood.csv")
ggplot(egr1.poor.good, aes(x = bev_resp, y = h_score)) +
  #scale_color_manual(values=c("#ffb464", "#126079")) +
  geom_violin(alpha = 0.5) +
  geom_dotplot(aes(fill = factor(animal_id)),
               binaxis= "y",
               stackdir = "center",
               dotsize = 0.5) +
  theme_prism() + 
  ggtitle("EGR1 Expression") +
  xlab("Bevacizumab Response Group") + ylab("Expression (H Score)")
ggsave("./egr1poorgood.svg", height = 6, width = 8)
wilcox.test((egr1.poor.good %>% dplyr::filter(bev_resp == "good"))$h_score, 
            (egr1.poor.good %>% dplyr::filter(bev_resp == "poor"))$h_score)

egr1.poor.placebo <- read.csv("./egr1poorplacebo.csv")
ggplot(egr1.poor.placebo %>% dplyr::filter(staining_type == "nuclear"), 
       aes(x = treatment, y = h_score)) +
  geom_violin(alpha = 0.5) +
  geom_dotplot(aes(fill = factor(animal_id)),
               binaxis= "y",
               stackdir = "center",
               dotsize = 0.5) +
  theme_prism() + 
  ggtitle("EGR1 Nuclear Expression") +
  xlab("Treatment Group") + ylab("Expression (H Score)")
ggsave("./egr1poorplacebonuc.svg", height = 6, width = 8)

ggplot(egr1.poor.placebo %>% dplyr::filter(staining_type == "cytoplasmic"), 
       aes(x = treatment, y = h_score)) +
  geom_violin(alpha = 0.5) +
  geom_dotplot(aes(fill = factor(animal_id)),
               binaxis= "y",
               stackdir = "center",
               dotsize = 0.5) +
  theme_prism() + 
  ggtitle("EGR1 Cytoplasmic Expression") +
  xlab("Treatment Group") + ylab("Expression (H Score)")
ggsave("./egr1poorplacebocyto.svg", height = 6, width = 8)
