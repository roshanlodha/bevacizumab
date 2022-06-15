#libraries----
setRepositories(ind = c(1:8))
for (package in c('biomaRt', 'tidyverse', 'tximport', 'ensembldb', 'EnsDb.Hsapiens.v86', 
                  'edgeR', 'matrixStats', 'cowplot', 'DESeq2', 'apeglm', 
                  'DT', 'plotly', 'gt',
                  'limma', 'ggrepel',
                  'GSEABase', 'Biobase', 'GSVA', 'gprofiler2', 
                  'clusterProfiler', 'msigdbr', 'enrichplot',
                  'ggthemes', 'ggprism',
                  'org.Hs.eg.db', 'annotate',
                  'clusterProfiler', 'enrichplot',
                  'gtsummary', 'ggforce', 'ggdendro',
                  'ggpubr')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

#prefiltered pca----
design <- read_tsv("./prefilterstudydesign.txt")
sampleLabels <- design$sample
group <- factor(design$group)

gbmexpr <- read_csv("./prefiltergbmexpr.csv")[2:13] #pre-filtered
gbmexpr.matrix <- as.matrix(gbmexpr[,-1])
rownames(gbmexpr.matrix) <- unlist(gbmexpr[,1])
myDGEList <- DGEList(gbmexpr.matrix)
myDGEList.filtered.norm <- calcNormFactors(myDGEList, method = "TMM") #normalize using TMM
log2.tpm.filtered.norm <- log2(gbmexpr.matrix+1)
log2.tpm.filtered.norm.df <- as_tibble(log2.tpm.filtered.norm, rownames = "geneID")

distance <- dist(t(log2.tpm.filtered.norm), method = "maximum")
clusters <- hclust(distance, method = "average")
den <- ggdendrogram(clusters) +
  labs(title="B. GBM PDX Dendrogram") +
  xlab("Distance") + 
  ylab("Sample") +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

ggsave(path = "./plots/", filename = "dendrogram.svg", 
       plot = den, height = 5, width=5)

sampleLabels <- substr(sampleLabels, 1, nchar(sampleLabels)-5)
pca.res <- prcomp(t(log2.tpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color = group) +
  geom_point(size=3) +
  geom_label_repel(aes(label=sampleLabels),hjust=0, vjust=0) +
  scale_color_manual(values=c("#ffb464", "#126079")) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="A. GBM PDX PCA Plot") +
  theme_prism() +
  theme(legend.position = c(0.9, 0.5))
ggsave(path = "./dim", filename = "prefilterpca.png", plot = pca.plot, height = 5, width = 7)
supplementalfigure1 <- ggarrange(pca.plot, den,
                               ncol = 2, nrow = 1)
ggexport(supplementalfigure1, filename = "supplementalfigure1.png", width = 1000, height = 500)

#preprocessing----
mart <- useMart("ENSEMBL_MART_ENSEMBL") #coding genes for cleaning
mart <- useDataset("hsapiens_gene_ensembl", mart)
coding_genes <- getBM(attributes = c( "hgnc_symbol"), 
                      filters = c("biotype"), 
                      values = list(biotype="protein_coding"), 
                      mart = mart)$hgnc_symbol

design <- read_tsv("./studydesign.txt")
sampleLabels <- design$sample
group <- factor(design$group)
mm <- model.matrix(~0 + group)

counts <- read.table(file = "././counts.tabular", header=TRUE, sep="\t")
colnames(counts) <- c("geneID", sampleLabels)
counts$gene <- getSYMBOL(as.character(counts$geneID), data='org.Hs.eg')
counts <- counts %>% select(gene, everything())
counts <- counts[rowSums(counts <= 0) <= 3, ] %>% drop_na() #filter: at most 3 zeros
write_csv(counts, "././counts.csv")
counts.matrix <- as.matrix(counts[3:11])
rownames(counts.matrix) <- counts$gene

#postfilter pca----
design <- read_tsv("./studydesign.txt")
sampleLabels <- design$sample
group <- factor(design$group)

pca.res <- prcomp(t(log(counts.matrix+1)), scale.=F, retx=T)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label = sampleLabels, color = group) +
  geom_point(size=3) +
  #geom_text(aes(label=sampleLabels),hjust=0, vjust=0) +
  #stat_ellipse() +
  scale_color_manual(values=c("#ffb464", "#126079")) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="GBM PCA plot") +
  #     caption=paste0("produced on ", Sys.time())) +
  stat_ellipse() +
  theme_prism()
ggsave(path = "./pca/", filename = "postfilterpca.png", plot = pca.plot, height = 5, width = 5)

#pca----
pca.res <- prcomp(t(log(counts.matrix+1)), scale.=F, retx=T)
pc.var <- pca.res$sdev^2
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label = sampleLabels, color = group) +
  geom_point(size=3) +
  #geom_text(aes(label=sampleLabels),hjust=0, vjust=0) +
  #stat_ellipse() +
  scale_color_manual(values=c("#ffb464", "#126079")) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="GBM PCA plot") +
  #     caption=paste0("produced on ", Sys.time())) +
  stat_ellipse() +
  theme_prism()
ggsave(path = "./pca/", filename = "pca.png", plot = pca.plot, height = 3.5, width = 6)

#deseq----
rownames(counts.matrix) <- counts$geneID
dds <- DESeqDataSetFromMatrix(countData = round(counts.matrix),
                              colData = design,
                              design = ~ group)
keep <- rowSums(counts(dds)) >= 350 #determined via hyperparameter exploration
dds <- dds[keep,]
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
deseqvoom <- voom(counts(dds, normalized=TRUE), mm, plot = T)
res <- results(dds)
res <- lfcShrink(dds, coef="group_poor_vs_good", type="ashr")
deseq <- as.data.frame(res) %>% drop_na() %>% arrange(padj) %>% arrange(desc(abs(log2FoldChange)))
deseq$gene <- getSYMBOL(rownames(deseq), data='org.Hs.eg')
deseq <- dplyr::filter(deseq, gene %in% coding_genes)
deseq <- deseq %>%
  mutate(enrichment = case_when(
    (padj < 0.05) & (log2FoldChange > 1) ~ "poor",
    (padj < 0.05) & (log2FoldChange < -1) ~ "good"))
deseq$enrichment[is.na(deseq$enrichment)] <- "none"
rownames(deseq) <- deseq$gene
write_csv(deseq, "./deseq/deseq.csv")

search <- function(x) {
  deseq %>% filter(grepl(x, gene))
}

#deseq plot----
deseqplot <- ggplot(deseq) +
  aes(y=-log10(padj), x=log2FoldChange, colour=enrichment) +
  scale_color_manual(values = c("good" = "#126079", "none"="grey", "poor"="#ffb464")) +
  geom_point(size=1.5, alpha=0.25) +
  #facet_zoom(xlim = c(9, 11), zoom.size = 1) +
  geom_rect(mapping=aes(xmin=9, xmax=11, ymin=2, ymax=54), alpha=0, color='black') +
  geom_text_repel(size=4, data=head(deseq, 20), aes(label=gene), max.overlaps = Inf, colour="black") +
  #geom_text_repel(size=3, data=subset(deseq, gene == 'EGR1'), aes(label=gene)) + 
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="black", size=.5) +
  geom_text(aes(-8, -log10(0.05), label = "p = 0.05", vjust = 1), colour="black") +
  labs(title="A. Differential RNA Expression") +
       #subtitle="Positive logFC indicates upregulation in poor Bevacizumab responders") +
  xlab("Fold Change (log2)") +
  ylab("Significance (log10)") +
  theme_prism() + 
  theme(legend.position = c(0.15, 0.8))
deseqplot
ggsave(filename = "./deseq/deseq.png", plot = deseqplot, height = 6, width=6)

kdr <- scan("././KDR_geneset.txt", character(), quote = "")
kdr <- deseq %>% filter((gene %in% kdr) 
                        & (enrichment!="none")
                        & (abs(log2FoldChange) > 4))

blood_vessel <- scan("././blood_vessel_geneset.txt", character(), quote = "")
blood_vessel <- deseq %>% filter((gene %in% blood_vessel) 
                                 & (enrichment!="none")
                                 & (abs(log2FoldChange) > 4))
write_csv(blood_vessel, "./deseq/blood_vessel.csv")

#kdr heatmap----
goi <- kdr$gene
kdr_heatmap <- counts %>% filter(gene %in% goi) %>% dplyr::select(-geneID)
rownames(kdr_heatmap) <- kdr_heatmap$gene
kdr_heatmap <- kdr_heatmap %>% dplyr::select(-gene)
kdr_heatmap <- as.data.frame(t(kdr_heatmap))
kdr_heatmap$sample <- rownames(kdr_heatmap)
rownames(kdr_heatmap) <- NULL
kdr_heatmap$group <- gsub("GBM[0-9]*_", "", kdr_heatmap$sample)
kdr_heatmap <- kdr_heatmap %>% pivot_longer(cols = -c("sample", "group"),
                                            names_to = "gene",
                                            values_to = "expression")

kdr_heatmap$log.expr <- log(kdr_heatmap$expression+1)
kdr_heatmap$sample <- substr(kdr_heatmap$sample,1,nchar(kdr_heatmap$sample)-5)
kdr_heatmap.plot <- ggplot(data = kdr_heatmap, mapping = aes(x = sample,
                                                             y = gene,
                                                             fill = log.expr)) +
  geom_tile() +
  scale_fill_gradient(high = "#ffb464", low = "#126079") +
  scale_colour_prism(palette = "colors") +
  xlab(label = "Patient Derived Xenograft") + # Add a nicer x-axis title
  ggtitle("B. RNA Expression Heatmap") +
  facet_grid(~ group, switch = "x", scales = "free_x", space = "free_x") + 
  #labs(color = "Your title here") +
  theme_prism() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.7),
        legend.position = "bottom")
kdr_heatmap.plot
ggsave(filename = "./deseq/kdr_heatmap.png", plot = kdr_heatmap.plot, height = 6, width=6)

#blood vessel heatmap----
goi <- blood_vessel$gene
b_v_heatmap <- counts %>% filter(gene %in% goi) %>% dplyr::select(-geneID)
rownames(b_v_heatmap) <- b_v_heatmap$gene
b_v_heatmap <- b_v_heatmap %>% dplyr::select(-gene)
b_v_heatmap <- as.data.frame(t(b_v_heatmap))
b_v_heatmap$sample <- rownames(b_v_heatmap)
rownames(b_v_heatmap) <- NULL
b_v_heatmap$group <- gsub("GBM[0-9]*_", "", b_v_heatmap$sample)
b_v_heatmap <- b_v_heatmap %>% pivot_longer(cols = -c("sample", "group"),
                                            names_to = "gene",
                                            values_to = "expression")

b_v_heatmap$log.expr <- log(b_v_heatmap$expression+1)
b_v_heatmap$sample <- substr(b_v_heatmap$sample,1,nchar(b_v_heatmap$sample)-5)

#blood vessel heatmap plot----
b_v_heatmap.plot <- ggplot(data = b_v_heatmap, mapping = aes(x = sample,
                                                             y = gene,
                                                             fill = log.expr)) +
  geom_tile() +
  scale_fill_gradient(high = "#ffb464", low = "#126079") +
  scale_colour_prism(palette = "colors") +
  xlab(label = "Patient Derived Xenograft") + # Add a nicer x-axis title
  ggtitle("C. Angiogenic Gene RNA Expression") +
  facet_grid(~ group, switch = "x", scales = "free_x", space = "free_x") + 
  #labs(color = "Your title here") +
  theme_prism() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 0.7),
        legend.position = "bottom")
b_v_heatmap.plot
ggsave(filename = "./deseq/b_v_heatmap.png", plot = b_v_heatmap.plot, height = 6, width=6)

#gsea----
hs_gsea <- msigdbr(species = "Homo sapiens")
hs_gsea %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
hs_gsea_h <- msigdbr(species = "Homo sapiens",
                      category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
hs_gsea_kegg <- msigdbr(species = "Homo sapiens",
                     category = "C2",
                     subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol)

deseq.GSEA.select <- dplyr::select(deseq, gene, log2FoldChange, padj)
deseq.gsea <- abs(deseq.GSEA.select$log2FoldChange)/deseq.GSEA.select$log2FoldChange * -log10(deseq.GSEA.select$padj)
names(deseq.gsea) <- as.character(deseq.GSEA.select$gene)
deseq.gsea <- sort(deseq.gsea, decreasing = TRUE)
deseq.gsea.res <- GSEA(deseq.gsea, pvalueCutoff = 1, TERM2GENE=hs_gsea_h, verbose=FALSE)
deseq.GSEA.df <- as_tibble(deseq.gsea.res@result)
deseq.GSEA.df <- deseq.GSEA.df %>%
  mutate(phenotype = case_when(
    (NES > 0) & (p.adjust<0.05) ~ "poor",
    (NES < 0) & (p.adjust<0.05) ~ "good"))
deseq.GSEA.df$phenotype[is.na(deseq.GSEA.df$phenotype)] <- "none"
deseq.GSEA.df$Description <- gsub("_", " ", deseq.GSEA.df$Description)
write_csv(deseq.GSEA.df, "./gsea/h.csv")

#gsea plot----
deseq.GSEA.plot <- ggplot(deseq.GSEA.df, aes(x=NES, y=-log10(p.adjust), color=phenotype)) + 
  geom_point(aes(size=setSize), alpha=0.5) +
  scale_color_manual(values = c("good" = "#126079", "none"="grey", "poor"="#ffb464")) +
  geom_text(aes(-1.5, -log10(0.05), label = "p = 0.05", vjust = -1)) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", size=.5) +
  geom_text_repel(size=4, data=(deseq.GSEA.df %>% dplyr::filter((NES > 0) & (p.adjust < 0.05)))[1,], 
                  aes(label=Description)) +
  labs(title="B. Hallmark Gene Set Enrichment") +
  #     subtitle="Positive NES indicates upregulation of gene set in poor Bevacizumab responders") +
  ylab("Significance (log10)") +
  xlab("Normalized Enrichment Score") +
  theme_prism() +
  theme(legend.position = "bottom")
deseq.GSEA.plot
ggsave(path = "./gsea/", filename = "h.png", plot = deseq.GSEA.plot, 
       width = 6, height = 6)

#figure 1----
figure1 <- ggarrange(ggarrange(deseqplot, deseq.GSEA.plot,
                               ncol = 1, nrow = 2),
                     b_v_heatmap.plot,
                     ncol = 2, nrow = 1)
ggexport(figure1, filename = "figure1.png", width = 1000, height = 1000)
