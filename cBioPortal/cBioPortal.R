library(survival)
library(survminer)
library(lubridate)
library(ggprism)

columns <- c("case", "study", "number_at_risk", "status", "survival", "time")
stati <- 1:2
names(stati) <- c("censored", "deceased")

#chrna7----
chrna7_altered <- read.csv("./chrna7_altered.tsv", sep="\t" )
colnames(chrna7_altered) <- columns
chrna7_altered$status <- stati[chrna7_altered$status]
chrna7_unaltered <- read.csv("./chrna7_unaltered.tsv", sep="\t" )
colnames(chrna7_unaltered) <- columns
chrna7_unaltered$status <- stati[chrna7_unaltered$status]
chrna7_altered$group <- "altered"
chrna7_unaltered$group <- "unaltered"
chrna7 <- rbind(chrna7_altered, chrna7_unaltered)
chrna7_plot <- ggsurvplot(
  fit = survfit(Surv(time, status) ~ group, data = chrna7), 
  risk.table = TRUE,
  xlab = "Months", 
  ylab = "Overall Survival",
  title = "Median OS Under Altered CHRNA7 Expression", 
  palette = c("#126079", "#ffb464"),
  censor = FALSE,
  ggtheme = theme_prism(),
  legend = c(0.8,0.8),
  xlim = c(0, 50),
  pval = TRUE, 
  pval.method = TRUE, 
  surv.median.line = "hv",
  break.time.by = 5)
ggsave(filename = "./chrna7.png", height = 6, width=8)

#acvrl1----
acvrl1_altered <- read.csv("./acvrl1_altered.tsv", sep="\t" )
colnames(acvrl1_altered) <- columns
acvrl1_altered$status <- stati[acvrl1_altered$status]
acvrl1_unaltered <- read.csv("./acvrl1_unaltered.tsv", sep="\t" )
colnames(acvrl1_unaltered) <- columns
acvrl1_unaltered$status <- stati[acvrl1_unaltered$status]
acvrl1_altered$group <- "altered"
acvrl1_unaltered$group <- "unaltered"
acvrl1 <- rbind(acvrl1_altered, acvrl1_unaltered)
acvrl1_plot <- ggsurvplot(
  fit = survfit(Surv(time, status) ~ group, data = acvrl1), 
  risk.table = TRUE,
  xlab = "Months", 
  ylab = "Overall Survival",
  title = "Median OS Under Altered ACVRL1 Expression",
  palette = c("#126079", "#ffb464"),
  censor = FALSE,
  ggtheme = theme_prism(),
  legend = c(0.8,0.8),
  xlim = c(0, 50), 
  pval = TRUE, 
  pval.method = TRUE, 
  surv.median.line = "hv",
  break.time.by = 5)
ggsave(filename = "./acvrl1.png", height = 6, width=6)

#helper functions----
estimate_survival <- function(data, time) {
  summary(survfit(Surv(time, status) ~ group, data = data), times = time)
}
