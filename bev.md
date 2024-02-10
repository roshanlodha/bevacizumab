Tumors that respond poorly to bevacizumab show upregulation of
angiogenesis genes.
================
Roshan Lodha
09 February, 2024

# Bulk RNA-sequencing analysis

## Loading required packages and tools

## Preprocessing

![](bev_files/figure-gfm/suppfig1-1.png)<!-- -->

## Postprocessing

## Bulk-RNA Analysis

### DGE Plot

![](bev_files/figure-gfm/suppfig2-1.png)<!-- -->

### GSEA Plot

    ## # A tibble: 23 × 2
    ##    gs_cat gs_subcat        
    ##    <chr>  <chr>            
    ##  1 C1     ""               
    ##  2 C2     "CGP"            
    ##  3 C2     "CP"             
    ##  4 C2     "CP:BIOCARTA"    
    ##  5 C2     "CP:KEGG"        
    ##  6 C2     "CP:PID"         
    ##  7 C2     "CP:REACTOME"    
    ##  8 C2     "CP:WIKIPATHWAYS"
    ##  9 C3     "MIR:MIRDB"      
    ## 10 C3     "MIR:MIR_Legacy" 
    ## # ℹ 13 more rows

![](bev_files/figure-gfm/suppfig3-1.png)<!-- -->

### Angiogenesis Gene Expression Plot

### Collagen Gene Expression Plot

### Figure 1

![](bev_files/figure-gfm/fig1-1.png)<!-- -->

# Protein Level Analysis

## EGR1 Expression

### Figure 2

![](bev_files/figure-gfm/fig2-1.png)<!-- -->

## Non-EGR1 Expression

### a7-nAChR Intensity

### SOX10 Expression

![](bev_files/figure-gfm/suppfig4-1.png)<!-- -->

## Correlation between a7-nAChR and EGR1

### Figure 3

![](bev_files/figure-gfm/fig3-1.png)<!-- -->

# Survival Analysis

## GBM Analysis

### TCGA RNA Analysis

#### EGR1-Stratified Angiogenesis Expression

#### EGR1-Stratified Collagen Expression

#### EGR1-Stratified Average Expression Boxplot

#### EGR1-Stratified Survival Plot

##### EGR1-Stratified Recurrent Tumor Survival Plot

![](bev_files/figure-gfm/EGR1-recurrent-survival-1.png)<!-- -->

![](bev_files/figure-gfm/suppfig6-1.png)<!-- -->

### Figure 5

![](bev_files/figure-gfm/fig5-1.png)<!-- -->

## Other Cancers

### Pan-Cancer Analysis

### Lung Cancer Analysis

## Other Genes

### CHRNA7 Survival Plot

### SOX10 Survival Plot

### EGR3 Survival Plot

### RAMP3 Survival Plot

### Other Genes Plot

![](bev_files/figure-gfm/suppfig8-1.png)<!-- -->
