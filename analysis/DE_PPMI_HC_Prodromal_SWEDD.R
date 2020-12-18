# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: June 4th, 2019

#Purpose: differential expression: healthy controls vs prodromal and swedd

#Imports:
# PPMI count table (filtered)

#Exports:
# 

## load libraries  ----------------------------------------------------------------------------------------------
library(tidyverse) #42; base for data/table manipulation
library(tidyr) #for pivot tables
library(edgeR) #expression normalization
library(limma) #more expression normalization
library(corrplot) #correlation plot matrix
library(beepr) #beep alerts when job is complete
library(ggpubr) #paper level visuals
library(seqsetvis) #ggplot2 venn
library(pheatmap) #generate pretty heatmaps
library(devtools) #load R scripts from github

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

outdir <- "analysis/DE/HC_prodromal_SWEDD/"
dir.create(outdir)

## functions to output DE contrast tables  ---------------------------------------------------------------------

resultsOut <- function(limmaFit, contr.matrix, outDir, logfc, pval){
  
  x <- 1
  while(x <= length(colnames(contr.matrix))){
    contrName <- colnames(contr.matrix)[x]
    print(contrName)
    
    dt <- decideTests(limmaFit)
    
    #output full table
    df <- rownames_to_column(topTable(efit, coef = x, n = Inf), var = "gene_id")
    df <- right_join(genes.anno, df)
    df.withMeans <- left_join(df, meanCPM)
    write_tsv(df.withMeans, paste(outDir, contrName,"_full.tsv", sep = ""))
    
    #output volcano plot
    png(paste(outDir, contrName, "_volcano.png", sep = ""))
    with(df, plot(logFC, -log10(P.Value), pch = 20, main = paste("Volcano Plot: ", contrName, sep = "")))
    with(subset(df, adj.P.Val < pval & abs(logFC) > logfc), points(logFC, -log10(P.Value), pch=20, col="red"))
    dev.off()
    
    #output MD plot
    png(paste(outDir, contrName, "_MDPlot.png", sep = ""))
    plotMD(limmaFit, column = x, status = dt[,x], main = contrName,
           xlim=c(-10,15))
    dev.off()
    
    #output cutoff table
    df.cut <- rownames_to_column(topTable(efit, coef = x, n = Inf, lfc = logfc, p.value = pval), var = "gene_id")
    df.cut <- right_join(genes.anno, df.cut)
    df.cut.withMeans <- left_join(df.cut, meanCPM)
    write_tsv(df.cut.withMeans, paste(outDir, contrName,"_cutoffs.tsv", sep = ""))
    
    
    x <- x + 1
  }
  
}

## load count table  ----------------------------------------------------------------------------------------------
countTable <- read_csv("data/rawCounts_poolPCAContributeRemoved.csv")
countTable <- column_to_rownames(countTable, var = "Ensembl_ID")

## biotypes and gene names ----------------------------------------------------------------------------------------
tmpG <- read_tsv("sourceFiles/GRCh38_GENCODE29_geneInfo.txt",
                 skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type")]
genes.anno <- unique(genes.anno)

## load metadata  ---------------------------------------------------------------------------------------------------
meta <- read_csv("data/megaMetaTable.csv")
dim(meta)

meta <- meta[meta$QCflagIR3 == "pass",] #filter out fails
dim(meta) #4756

meta <- meta[!meta$PoolAssign %in% c("PD_POOL", "HC_POOL"),] #filter out pools
dim(meta) #4650


#CAN START HERE AND LOAD SAMPLE DATA TABLE NOW THAT IT EXISTS
metaSub <- read_csv(paste(outdir, "differentialExpression_sampleData.csv", sep = ""))
#order to match
countsSub <- countTable[,names(countTable) %in% sort(metaSub$HudAlphaSampleName)]

i <- match(metaSub$HudAlphaSampleName, names(countsSub), nomatch = 0)
countsSub <- countsSub[,i]

identical(names(countsSub), metaSub$HudAlphaSampleName) #double checking

## sample design --------------------------------------------------------------------------

disease <- factor(gsub(" ", "_", metaSub$Disease_Status), levels = c("Healthy_Control",
                                                     "Prodromal",
                                                     "SWEDD"))
age <- factor(metaSub$ageBin, levels = c("under_55", "55_to_65", "over_65"))
sex <- factor(metaSub$SEX,
              levels = c("Female", "Male"))
plate <- factor(metaSub$Plate)
usableBases <- metaSub$PCT_USABLE_BASES

design <- model.matrix(~ 0 + disease + sex + plate + age + usableBases)
colnames(design) <- gsub("disease", "", colnames(design))
colnames(design)

contr.matrix <- makeContrasts(
  Prodromal_HC = Prodromal - Healthy_Control,
  SWEDD_HC = SWEDD - Healthy_Control,
  levels = colnames(design)
)

contr.matrix

## normalization and filtering in limma --------------------------------------------------------------------------

keep <- filterByExpr(countsSub, group = disease)
sum(keep) #20376 genes

dge <- DGEList(countsSub[keep,])
dge  <- calcNormFactors(dge)

CPM <- cpm(dge)

logCPM <- cpm(dge, log = TRUE, prior.count = 3)
write_tsv(rownames_to_column(as.data.frame(logCPM), var = "gene_id"), paste(outdir, "logCPM.tsv", sep = ""))


v <- voom(dge, design, plot=TRUE) #check voom plot for curve to see if we need to do more filtering
vfit <- lmFit(v, design)
vfit.con1 <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit.con1)
plotSA(efit)
summary(decideTests(efit))

DE.results <- as.data.frame(summary(decideTests(efit))) %>%
  pivot_wider(names_from = Var2, values_from = Freq)

write_tsv(DE.results, paste(outdir, "DE_result_summary.tsv", sep = ""))


## make table of mean CPM for various groups --------------------------------------------------------------------
#make table of mean expression for each group
pro.Group <- metaSub$HudAlphaSampleName[metaSub$Disease_Status == "Prodromal"]
swedd.Group <- metaSub$HudAlphaSampleName[metaSub$Disease_Status == "SWEDD"]
control.Group <- metaSub$HudAlphaSampleName[metaSub$Disease_Status == "Healthy Control"]

pro.CPM <- logCPM[,pro.Group]
swedd.CPM <- logCPM[,swedd.Group]
control.CPM <- logCPM[,control.Group]

meanCPM <- data.frame("gene_id" = row.names(logCPM),
                      "AvgExpr_Healthy_Control" = rowMeans(control.CPM),
                      "AvgExpr_Prodromal" = rowMeans(pro.CPM),
                      "AvgExpr_SWEDD" = rowMeans(swedd.CPM)
                      
)

#output DE results ------------------------------------------------------------------------------------------------

resultsOut(efit, contr.matrix, outdir, logfc = 0.1, pval = 0.05)
