# packages and data
pkgs <- c("UCSCXenaTools", "dplyr", "survival", "survminer", "ggbreak", "ggprism", "svglite")
#BiocManager::install(pkgs)
invisible(lapply(pkgs, function (x) suppressMessages(library(x, character.only = T))))

gbm_cohort = XenaData %>%
  filter(XenaHostNames == "tcgaHub") %>%
  XenaScan("TCGA Glioblastoma") # microarray dataset, CNVA dataset DNA level

#download clinical data-- --
cli_query = gbm_cohort %>%
  filter(DataSubtype == "phenotype") %>% # select clinical dataset
  XenaGenerate() %>% # generate a XenaHub object
  XenaQuery() %>%
  XenaDownload()
cli = XenaPrepare(cli_query)

# download gene expression data
ge = gbm_cohort %>%
  filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq")

ge = gbm_cohort %>%
  filter(DataSubtype == "protein expression RPPA", Label == "RPPA (replicate-base normalization)")
#try AFFYmetrix
EGR1 = fetch_dense_values(host = ge$XenaHosts,
  dataset = ge$XenaDatasets,
  identifiers = "EGR1",
  use_probeMap = TRUE) %>% .[1, ]
EGR3 = fetch_dense_values(host = ge$XenaHosts,
  dataset = ge$XenaDatasets,
  identifiers = "EGR3",
  use_probeMap = TRUE) %>% .[1, ]
SOX10 = fetch_dense_values(host = ge$XenaHosts,
  dataset = ge$XenaDatasets,
  identifiers = "SOX10",
  use_probeMap = TRUE) %>% .[1, ]
RAMP3 = fetch_dense_values(host = ge$XenaHosts,
  dataset = ge$XenaDatasets,
  identifiers = "RAMP3",
  use_probeMap = TRUE) %>% .[1, ]
CHRNA7 = fetch_dense_values(host = ge$XenaHosts,
  dataset = ge$XenaDatasets,
  identifiers = "CHRNA7",
  use_probeMap = TRUE) %>% .[1, ]

#SOX10-- --
merged_SOX10 = tibble(sample = names(SOX10),
    SOX10_expression = as.numeric(SOX10)) %>%
  left_join(cli$GBM_survival.txt, by = "sample") %>%
  select(sample, SOX10_expression, OS.time, OS) %>%
  rename(time = OS.time,
    status = OS)
ggplot(merged_SOX10, aes(x = SOX10_expression)) +
  geom_histogram(color="black", fill="white") +
  theme_prism()

fit_SOX10 = coxph(Surv(time, status) ~ SOX10_expression, data = merged_SOX10)
fit_SOX10

merged_SOX10 = merged_SOX10 %>%
  mutate(group = case_when(
    SOX10_expression > quantile(SOX10_expression, 0.9) ~'High',
    (SOX10_expression < quantile(SOX10_expression, 0.9) &
      SOX10_expression > quantile(SOX10_expression, 0.1)) ~'Normal',
    SOX10_expression < quantile(SOX10_expression, 0.1) ~'Low',
    TRUE~NA_character_
  )) %>%
  mutate(z = (SOX10_expression - mean(SOX10_expression)) / sd(SOX10_expression)) %>%
  mutate(group = case_when(
    z > 1.5~'High',
    z < -1.5~'Low',
    (z < 1.5) & (z > -1.5) ~'Normal',
    TRUE~NA_character_
  ))
fit_SOX10 = survfit(Surv(time, status) ~group,
  data = merged_SOX10 %>% dplyr::filter(group != "Low"))
SOX10_plot <- ggsurvplot(fit_SOX10,
  pval = TRUE,
  pval.coord = c(600, 0.5),
  #xlim = c(0, 1000),
  palette = c("#ffb464", "#126079"),
  #conf.int = TRUE,
  #pval = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs = c("High Expression", "Normal Expression"),
  surv.median.line = "hv",
  break.time.by = 250,
    ggtheme = theme_prism(),
    legend = c(0.7, 0.8),
    title = "SOX10-expression Stratisfied Survival Plot")
ggsave("../plots/survival/SOX10survival.png", plot = print(SOX10_plot), height = 6, width = 6)

#EGR3-- --
merged_EGR3 = tibble(sample = names(EGR3),
    EGR3_expression = as.numeric(EGR3)) %>%
  left_join(cli$GBM_survival.txt, by = "sample") %>%
  #filter(sample_type == "Primary Tumor") %>% # Keep only 'Primary Tumor'
  select(sample, EGR3_expression, OS.time, OS) %>%
  rename(time = OS.time,
    status = OS)

fit_EGR3 = coxph(Surv(time, status) ~EGR3_expression, data = merged_EGR3)
fit_EGR3

merged_EGR3 = merged_EGR3 %>%
  mutate(group = case_when(
    EGR3_expression > quantile(EGR3_expression, 0.9) ~'High',
    (EGR3_expression < quantile(EGR3_expression, 0.9) &
      EGR3_expression > quantile(EGR3_expression, 0.1)) ~'Normal',
    EGR3_expression < quantile(EGR3_expression, 0.1) ~'Low',
    TRUE~NA_character_
  )) %>%
  mutate(z = (EGR3_expression - mean(EGR3_expression)) / sd(EGR3_expression)) %>%
  mutate(group = case_when(
    z > 1.5~'High',
    z < -1.5~'Low',
    (z < 1.5) & (z > -1.5) ~'Normal',
    TRUE~NA_character_
  ))
fit_EGR3 = survfit(Surv(time, status) ~group,
  data = merged_EGR3 %>% dplyr::filter(group != "Low"))
EGR3_plot <- ggsurvplot(fit_EGR3,
  pval = TRUE,
  pval.coord = c(600, 0.5),
  #xlim = c(0, 1000),
  palette = c("#ffb464", "#126079"),
  #conf.int = TRUE,
  #pval = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs = c("High Expression", "Normal Expression"),
  surv.median.line = "hv",
  break.time.by = 250,
    ggtheme = theme_prism(),
    legend = c(0.7, 0.8),
    title = "EGR3-expression Stratisfied Survival Plot")
ggsave("../plots/survival/EGR3survival.png", plot = print(EGR3_plot), height = 6, width = 6)


#merge-- --
merged_EGR1 = tibble(sample = names(EGR1),
    EGR1_expression = as.numeric(EGR1)) %>%
  left_join(cli$GBM_survival.txt, by = "sample") %>%
  #filter(sample_type == "Primary Tumor") %>% # Keep only 'Primary Tumor'
  select(sample, EGR1_expression, OS.time, OS) %>%
  rename(time = OS.time,
    status = OS)

fit_EGR1 = coxph(Surv(time, status) ~EGR1_expression, data = merged_EGR1)
fit_EGR1

merged_EGR1 = merged_EGR1 %>%
  mutate(group = case_when(
    EGR1_expression > quantile(EGR1_expression, 0.9) ~'High',
    (EGR1_expression < quantile(EGR1_expression, 0.9) &
      EGR1_expression > quantile(EGR1_expression, 0.1)) ~'Normal',
    EGR1_expression < quantile(EGR1_expression, 0.1) ~'Low',
    TRUE~NA_character_
  )) %>%
  mutate(z = (EGR1_expression - mean(EGR1_expression)) / sd(EGR1_expression)) %>%
  mutate(group = case_when(
    z > 1.5~'High',
    z < -1.5~'Low',
    (z < 1.5) & (z > -1.5) ~'Normal',
    TRUE~NA_character_
  ))
fit_EGR1 = survfit(Surv(time, status) ~group,
  data = merged_EGR1 %>% dplyr::filter(group != "Low"))
EGR1_plot <- ggsurvplot(fit_EGR1,
  pval = TRUE,
  pval.coord = c(600, 0.5),
  #xlim = c(0, 1000),
  palette = c("#ffb464", "#126079"),
  #conf.int = TRUE,
  #pval = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs = c("High Expression", "Normal Expression"),
  surv.median.line = "hv",
  break.time.by = 250,
    ggtheme = theme_prism(),
    legend = c(0.7, 0.8),
    title = "EGR1-expression Stratisfied Survival Plot")
ggsave("../plots/survival/EGR1survival.svg", plot = print(EGR1_plot), height = 6, width = 6)

# RAMP3
merged_RAMP3 = tibble(sample = names(RAMP3),
    RAMP3_expression = as.numeric(RAMP3)) %>%
  left_join(cli$GBM_survival.txt, by = "sample") %>%
  #filter(sample_type == "Primary Tumor") %>% # Keep only 'Primary Tumor'
select(sample, RAMP3_expression, OS.time, OS) %>%
  rename(time = OS.time,
    status = OS)

fit = coxph(Surv(time, status) ~RAMP3_expression, data = merged_RAMP3)
fit

merged_RAMP3 = merged_RAMP3 %>%
  mutate(group = case_when(
    RAMP3_expression > quantile(RAMP3_expression, 0.9) ~'High',
    (RAMP3_expression < quantile(RAMP3_expression, 0.9) &
      RAMP3_expression > quantile(RAMP3_expression, 0.1)) ~'Normal',
    RAMP3_expression < quantile(RAMP3_expression, 0.1) ~'Low',
    TRUE~NA_character_
  )) %>%
  mutate(z = (RAMP3_expression - mean(RAMP3_expression)) / sd(RAMP3_expression)) %>%
  mutate(group = case_when(
    z > 1.5~'High',
    z < -1.5~'Low',
    (z < 1.5) & (z > -1.5) ~'Normal',
    TRUE~NA_character_
  ))
fit_RAMP3 = survfit(Surv(time, status) ~group,
  data = merged_RAMP3 %>% dplyr::filter(group != "Low"))
RAMP3_plot <- ggsurvplot(fit_RAMP3,
  pval = TRUE,
  pval.coord = c(600, 0.5),
  #xlim = c(0, 1000),
  palette = c("#ffb464", "#126079"),
  #conf.int = TRUE,
  #pval = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs = c("High Expression", "Normal Expression"),
  surv.median.line = "hv",
  break.time.by = 250,
    ggtheme = theme_prism(),
    legend = c(0.7, 0.8),
    title = "RAMP3-expression Stratisfied Survival Plot")

#CHRNA7-- --
merged_CHRNA7 = tibble(sample = names(CHRNA7),
    CHRNA7_expression = as.numeric(CHRNA7)) %>%
  left_join(cli$GBM_survival.txt, by = "sample") %>%
  #filter(sample_type == "Primary Tumor") %>% # Keep only 'Primary Tumor'
select(sample, CHRNA7_expression, OS.time, OS) %>%
  rename(time = OS.time,
    status = OS)

fit = coxph(Surv(time, status) ~CHRNA7_expression, data = merged_CHRNA7)
fit

merged_CHRNA7 = merged_CHRNA7 %>%
  mutate(group = case_when(
    CHRNA7_expression > quantile(CHRNA7_expression, 0.9) ~ 'High',
    (CHRNA7_expression < quantile(CHRNA7_expression, 0.9) &
      CHRNA7_expression > quantile(CHRNA7_expression, 0.1)) ~ 'Normal',
    CHRNA7_expression < quantile(CHRNA7_expression, 0.1) ~ 'Low',
    TRUE ~ NA_character_
  )) %>%
  mutate(z = (CHRNA7_expression - mean(CHRNA7_expression)) / sd(CHRNA7_expression)) %>%
  mutate(group = case_when(
    z > 1.5~'High',
    z < -1.5~'Low',
    (z < 1.5) & (z > -1.5) ~'Normal',
    TRUE~NA_character_
  ))

fit_CHRNA7 = survfit(Surv(time, status) ~group,
  data = merged_CHRNA7 %>% dplyr::filter(group != "Low"))
CHRNA7_plot <- ggsurvplot(fit_CHRNA7,
  pval = TRUE,
  pval.coord = c(600, 0.5),
  #xlim = c(0, 1000),
  palette = c("#ffb464", "#126079"),
  #conf.int = TRUE,
  #pval = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs = c("High Expression", "Normal Expression"),
  surv.median.line = "hv",
  break.time.by = 250,
    ggtheme = theme_prism(),
    legend = c(0.7, 0.8),
    title = "CHRNA7-expression Stratisfied Survival Plot")

# EGR1 and SOX10
EGR1_SOX10 <- merge(merged_EGR1, merged_SOX10, by = c("sample", "time", "status")) %>%
  select("sample", "time", "status", "z.x", "z.y") %>%
  mutate(group = case_when(
    (z.x > 1) & (z.y > 1) ~'High',
    (z.x < -1) & (z.y < -1) ~'Low',
    TRUE~'Normal'
  ))
EGR1_SOX10_corr <- ggplot(data = EGR1_SOX10, aes(x = z.x, y = z.y)) +
  geom_point() +
  theme_prism() +
  labs(title = "EGR1 and SOX10 expression correlation") +
  xlab("EGR1 expression") +
  ylab("SOX10 expression") +
  geom_smooth(method = "lm")
ggsave("../plots/survival/EGR1.SOX10.corr.svg", plot = EGR1_SOX10_corr, height = 6, width = 6)

fit_EGR1_SOX10 = survfit(Surv(time, status) ~group,
  data = EGR1_SOX10 %>% dplyr::filter(group != "Low"))
EGR1_SOX10_plot <- ggsurvplot(fit_EGR1_SOX10,
  pval = TRUE,
  pval.coord = c(600, 0.5),
  #xlim = c(0, 1000),
  palette = c("#ffb464", "#126079"),
  #conf.int = TRUE,
  #pval = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs = c("High Expression", "Normal Expression"),
  surv.median.line = "hv",
  break.time.by = 250,
    ggtheme = theme_prism(),
    legend = c(0.7, 0.8),
    title = "EGR1 and SOX10 expression Stratisfied Survival Plot")

# EGR1 and CHRNA7
EGR1_CHRNA7 <- merge(merged_EGR1, merged_CHRNA7, by = c("sample", "time", "status")) %>%
  select("sample", "time", "status", "z.x", "z.y") %>%
  mutate(group = case_when(
    (z.x > 1) & (z.y > 1) ~'High',
    (z.x < -1) & (z.y < -1) ~'Low',
    TRUE~'Normal'
  ))
EGR1_CHRNA7_corr <- ggplot(data = EGR1_CHRNA7, aes(x = z.x, y = z.y)) +
  geom_point() +
  theme_prism() +
  labs(title = "EGR1 and CHRNA7 expression correlation") +
  xlab("EGR1 expression") +
  ylab("CHRNA7 expression") +
  geom_smooth(method = "lm")
ggsave("../plots/survival/EGR1.CHRNA7.corr.svg", plot = EGR1_CHRNA7_corr, height = 6, width = 6)

fit_EGR1_CHRNA7 = survfit(Surv(time, status) ~group,
  data = EGR1_CHRNA7 %>% dplyr::filter(group != "Low"))
EGR1_CHRNA7_plot <- ggsurvplot(fit_EGR1_CHRNA7,
  pval = TRUE,
  pval.coord = c(600, 0.5),
  #xlim = c(0, 1000),
  palette = c("#ffb464", "#126079"),
  #conf.int = TRUE,
  #pval = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs = c("High Expression", "Normal Expression"),
  surv.median.line = "hv",
  break.time.by = 250,
    ggtheme = theme_prism(),
    legend = c(0.7, 0.8),
    title = "EGR1 and CHRNA7 expression Stratisfied Survival Plot")

# EGR1 and RAMP3
EGR1_RAMP3 <- merge(merged_EGR1, merged_RAMP3, by = c("sample", "time", "status")) %>%
  select("sample", "time", "status", "z.x", "z.y") %>%
  mutate(group = case_when(
    (z.x > 1.5) & (z.y > 1.5) ~'High',
    (z.x < -1.5) & (z.y < -1.5) ~'Low',
    TRUE~'Normal'
  ))

EGR1_RAMP3_corr <- ggplot(data = EGR1_RAMP3, aes(x = z.x, y = z.y)) +
  geom_point() +
  theme_prism() +
  labs(title = "EGR1 and CHRNA7 expression correlation") +
  xlab("EGR1 expression") +
  ylab("CHRNA7 expression") +
  geom_smooth(method = "lm")
ggsave("../plots/survival/EGR1.RAMP3.corr.svg", plot = EGR1_RAMP3_corr, height = 6, width = 6)

fit_EGR1_RAMP3 = survfit(Surv(time, status) ~group,
  data = EGR1_RAMP3 %>% dplyr::filter(group != "Low"))
EGR1_RAMP3_plot <- ggsurvplot(fit_EGR1_RAMP3,
  pval = TRUE,
  pval.coord = c(600, 0.5),
  #xlim = c(0, 1000),
  palette = c("#ffb464", "#126079"),
  #conf.int = TRUE,
  #pval = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.labs = c("High Expression", "Normal Expression"),
  surv.median.line = "hv",
  break.time.by = 250,
    ggtheme = theme_prism(),
    legend = c(0.7, 0.8),
    title = "EGR1 and RAMP3 expression Stratisfied Survival Plot")