# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: September 9th, 2019

#Purpose: differential expression

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


toolkitPath <- "~/Dropbox/repos/toolkit_ehutchins/" #for custom scripts and gene annotation files
source(paste(toolkitPath, "R-functions/createColorPalettes.R", sep = "")) #color palettes

outdir <- "analysis/DE/allPD_vs_allCon_sizeMatched/"
dir.create(outdir)

#resultsOut function ---------------------------------------------------------
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
tmpG <- read_tsv(paste0(toolkitPath, "sourceFiles/GRCh38_GENCODE29_geneInfo.txt"),
                 skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type")]
genes.anno <- unique(genes.anno)

## load metadata to match size matched data --------------------------------------------------------------------------------------
metaSub <- read_csv("data/differentialExpression_sampleData.csv")


#order to match
i <- match(metaSub$HudAlphaSampleName, names(countsSub), nomatch = 0)
countsSub <- countsSub[,i]

identical(names(countsSub), metaSub$HudAlphaSampleName) #double checking



## sample design --------------------------------------------------------------------------

visit <- factor(metaSub$CLINICAL_EVENT, levels = c("BL", "V02", "V04", "V06", "V08"))
subject <- metaSub$PATNO
disease <- factor(gsub(" ", "", metaSub$Disease_Status), levels = c("HealthyControl", "IdiopathicPD", "GeneticPD", "GeneticUnaffected"))
PD_status <- factor(disease, levels = c("Unaffected", "Affected"))
PD_status[disease %in% c("HealthyControl", "GeneticUnaffected")] <- "Unaffected"
PD_status[disease %in% c("IdiopathicPD", "GeneticPD")] <- "Affected"
group <- factor(paste(disease, visit, sep = "_"))
age <- factor(metaSub$ageBin, levels = c("under_55", "55_to_65", "over_65"))
sex <- factor(metaSub$SEX,
              levels = c("Female", "Male"))
plate <- factor(metaSub$Plate)
usableBases <- metaSub$PCT_USABLE_BASES

design <- model.matrix(~ 0 + PD_status + sex + plate + age + usableBases)
colnames(design) <- gsub("PD_status", "", colnames(design))
colnames(design)


contr.matrix <- makeContrasts(
  PD = Affected - Unaffected,
  levels = colnames(design)
)

contr.matrix

## normalization and filtering in limma --------------------------------------------------------------------------

keep <- filterByExpr(countsSub, group = PD_status)
sum(keep) #18996 genes

dge <- DGEList(countsSub[keep,])
dge  <- calcNormFactors(dge)


logCPM <- cpm(dge, log = TRUE, prior.count = 3)


#make table of mean expression for each group
PDGroup <- metaSub$HudAlphaSampleName[metaSub$Disease_Status %in% c("Idiopathic PD", "Genetic PD")]
conGroup <- metaSub$HudAlphaSampleName[metaSub$Disease_Status %in% c("Healthy Control", "Genetic Unaffected")]

PDCPM <- logCPM[,PDGroup]
conCPM <- logCPM[,conGroup]

meanCPM <- data.frame("gene_id" = row.names(logCPM),
                      "AvgExpr_allPD" = rowMeans(PDCPM),
                      "AvgExpr_allControls" = rowMeans(conCPM))


v <- voom(dge, design, plot=TRUE) #check voom plot for curve to see if we need to do more filtering
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
summary(decideTests(efit))


#output DE results ------------------------------------------------------------------------------------------------

resultsOut(efit, contr.matrix, outdir, logfc = 0.1, pval = 0.05)
