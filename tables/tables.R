library(sjPlot)
library(dplyr)
library(tidyr)
library(readr)

design <- read_tsv("raw/studydesign.txt")
tab_df(design,
       alternate.rows = T,
       title = "Study design", 
       file="./tables/design.doc")

dge <- read_csv("deseq/deseq.csv")
dge <- dge %>% select(gene, log2FoldChange, pvalue, padj, baseMean) 
dge$log2FoldChange <- 2 ^ dge$log2FoldChange
colnames(dge) <- c("Gene", "Fold Change", "p-value", "Adjusted p-value", "Mean Expression")
tab_df(head(dge, 20),
       alternate.rows = T,
       title = "Differential Gene Expression", 
       file="./tables/DGE.doc")

gsea <- read_csv("gsea/h.csv")
gsea <- gsea %>% select(Description, NES, p.adjust, core_enrichment) 
colnames(gsea) <- c("Gene Set", "Normalized Enrichment", "Adjusted P-value", "Core Enrichment")
tab_df(head(gsea, 20),
       alternate.rows = T,
       title = "Gene Set Enrichment", 
       file="./tables/GSEA.doc")

angiogenesis <- read_csv("deseq/blood_vessel.csv")
angiogenesis <- angiogenesis %>% select(gene, log2FoldChange, pvalue, padj, baseMean) 
angiogenesis$log2FoldChange <- 2 ^ angiogenesis$log2FoldChange
colnames(angiogenesis) <- c("Gene", "Fold Change", "p-value", "Adjusted p-value", "Mean Expression")
tab_df(angiogenesis,
       alternate.rows = T,
       title = "Angiogenesis Differential Gene Expression", 
       file="./tables/angiogenesis.doc")
