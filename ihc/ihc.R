library(tidyverse)
library(ggprism)
library(ggpubr)

# EGR1
egr1.poor.good <- read.csv("./egr1poorgood.csv")
egr1.poor.good$fullname <- paste(egr1.poor.good$pdx_id, egr1.poor.good$animal_id)
ggplot(egr1.poor.good, aes(x = fct_rev(bev_resp), y = h_score)) +
  geom_boxplot(alpha = 0.5, fill = c("#ffb464", "#126079")) +
  geom_dotplot(aes(fill = factor(fullname)),
    binaxis = "y",
    stackdir = "center",
    dotsize = 0.5,
    binpositions = "all",
    stackgroups = TRUE) +
  theme_prism() +
  theme(
    legend.key.height = unit(10, "pt"),
    legend.title = element_text()) +
  guides(fill = guide_legend(title = "Animal ID")) +
  #stat_compare_means(label.x.npc = "left", label.y.npc = "bottom") +
  ggtitle("EGR1 Expression (Nuclear)") +
  xlab("Bevacizumab Response Group") + ylab("Expression (H Score)")
ggsave("../plots/ihc/egr1poorgood.png", height = 6, width = 6)
wilcox.test((egr1.poor.good %>% dplyr::filter(bev_resp == "good")) $h_score,
  (egr1.poor.good %>% dplyr::filter(bev_resp == "poor")) $h_score)

egr1.poor.placebo <- read.csv("./egr1poorplacebo.csv")
egr1.poor.placebo$fullname <- paste(egr1.poor.placebo$pdx_id, egr1.poor.placebo$animal_id)
ggplot(egr1.poor.placebo %>% dplyr::filter(staining_type == "nuclear"),
    aes(x = treatment, y = h_score)) +
  geom_boxplot(alpha = 0.5, fill = c("#ffb464", "#126079")) +
  geom_dotplot(aes(fill = factor(fullname)),
    binaxis = "y",
    stackdir = "center",
    dotsize = 0.5,
    binpositions = "all",
    stackgroups = TRUE) +
  theme_prism() +
  theme(
    legend.key.height = unit(10, "pt"),
    legend.title = element_text()) +
  guides(fill = guide_legend(title = "Animal ID")) +
  ggtitle("EGR1 Nuclear Expression") +
  xlab("Treatment Group") + ylab("Expression (H Score)")
ggsave("../plots/ihc/egr1poorplacebonuc.png", height = 6, width = 6)

ggplot(egr1.poor.placebo %>% dplyr::filter(staining_type == "cytoplasmic"),
    aes(x = treatment, y = h_score)) +
  geom_boxplot(alpha = 0.5, fill = c("#ffb464", "#126079")) +
  geom_dotplot(aes(fill = factor(fullname)),
    binaxis = "y",
    stackdir = "center",
    dotsize = 0.5,
    binpositions = "all",
    stackgroups = TRUE) +
  theme_prism() +
  theme(
    legend.key.height = unit(10, "pt"),
    legend.title = element_text()) +
  guides(fill = guide_legend(title = "Animal ID")) +
  #stat_compare_means(label.x.npc = "left", label.y.npc = "bottom") +
  ggtitle("EGR1 Cytoplasmic Expression") +
  xlab("Treatment Group") + ylab("Expression (H Score)")
ggsave("../plots/ihc/egr1poorplacebocyto.png", height = 6, width = 6)

# Chi-squared
rawihc <- read_csv("./ihcraw.csv")
ihc <- rawihc %>%
  group_by(animal_id, bev_resp) %>%
  summarise_at(vars(neg, wk, plus_one), list(mean = mean))

#write.csv(ihc, "./ihc.csv")
#sortedihc <- read_csv("./ihc.csv")
sig <- ihc %>%
  rowwise() %>%
  mutate(
    test_stat = chisq.test(c(neg_mean, wk_mean, plus_one_mean))$statistic,
    p_val = chisq.test(c(neg_mean, wk_mean, plus_one_mean))$p.value
  )
write.csv(sig, "ihc.csv")

ihc <- rawihc %>%
  group_by(bev_resp) %>%
  summarise_at(vars(neg, wk, plus_one), list(mean = mean))
chisq.test(c(0.103, 0.704, 0.187), c(0.0229, 0.38, 0.55))