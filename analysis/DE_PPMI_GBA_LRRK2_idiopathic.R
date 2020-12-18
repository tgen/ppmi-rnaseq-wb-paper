# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: September 6th, 2019

#Purpose: differential expression

#Imports:
# PPMI count table (filtered)

#Exports:
# 


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

outdir <- "analysis/DE/GBA_LRRK2_idio_sizeMatched/"
dir.create(outdir)

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

## subset to group of interest; match count data with metadata -----------------------------------------------------
# metaSub <- meta[meta$Genetic_Status %in% c("GBA+", "LRRK2+", "LRRK2-/SNCA-/GBA-"),]
# unique(metaSub$Disease_Status)
# metaSub <- metaSub[metaSub$Disease_Status %in% c("Genetic Unaffected", "Genetic PD",
#                                                  "Healthy Control", "Idiopathic PD"),]
# metaSub$GenDiseaseStatus <- paste(metaSub$Genetic_Status, metaSub$Disease_Status, sep = "_")
# metaSub$GenDiseaseStatus <- gsub(" ", "_", metaSub$GenDiseaseStatus)
# metaSub$GenDiseaseStatus <- gsub("LRRK2-/SNCA-/GBA-_", "", metaSub$GenDiseaseStatus)
# metaSub$GenDiseaseStatus <- gsub("Genetic_", "", metaSub$GenDiseaseStatus)
# 
# 
# dim(metaSub)
# unique(metaSub$Genetic_Status)
# table(metaSub$Disease_Status, metaSub$CLINICAL_EVENT)
# rowSums(table(metaSub$Disease_Status, metaSub$CLINICAL_EVENT))
# colSums(table(metaSub$PATNO, metaSub$Disease_Status, metaSub$Genetic_Status) > 0) #302 patients; 126 PD, 176 unaffected, 197 HC, 177 idiopathic PD
# colSums(table(metaSub$PATNO, metaSub$GenDiseaseStatus) > 0) #302 patients; 126 PD, 176 unaffected, 197 HC, 177 idiopathic PD
# table(metaSub$GenDiseaseStatus) #samples
# 
# 
# set.seed(42)
# #downsample to same number of samples
# GBAPD <- subset(metaSub, GenDiseaseStatus == "GBA+_PD")
# GBAPD.IDs <- GBAPD$HudAlphaSampleName[sample(305,275)]
# GBAUn <- subset(metaSub, GenDiseaseStatus == "GBA+_Unaffected")
# GBAUn.IDs <- GBAUn$HudAlphaSampleName
# LRRK2PD <- subset(metaSub, GenDiseaseStatus == "LRRK2+_PD")
# LRRK2PD.IDs <- LRRK2PD$HudAlphaSampleName[sample(629,275)]
# LRRK2Un <- subset(metaSub, GenDiseaseStatus == "LRRK2+_Unaffected")
# LRRK2Un.IDs <- LRRK2Un$HudAlphaSampleName[sample(534,275)]
# HC <- subset(metaSub, GenDiseaseStatus == "Healthy_Control")
# HC.IDs <- HC$HudAlphaSampleName[sample(838,275)]
# idio <- subset(metaSub, GenDiseaseStatus == "Idiopathic_PD")
# idio.IDs <- idio$HudAlphaSampleName[sample(1445,275)]
# 
# 
# #sampleIDs <- metaSub$HudAlphaSampleName
# sampleIDs <- rbind(GBAPD.IDs, GBAUn.IDs)
# sampleIDs <- rbind(sampleIDs, LRRK2PD.IDs)
# sampleIDs <- rbind(sampleIDs, LRRK2Un.IDs)
# sampleIDs <- rbind(sampleIDs, HC.IDs)
# sampleIDs <- rbind(sampleIDs, idio.IDs)
# length(sampleIDs) #1650
# 
# countsSub <- countTable[,names(countTable) %in% sort(sampleIDs)]
# dim(countsSub)
# 
# metaSub <- metaSub[metaSub$HudAlphaSampleName %in% colnames(countsSub),]
# dim(metaSub)
# 
# rowSums(table(metaSub$GenDiseaseStatus, metaSub$CLINICAL_EVENT))
# colSums(table(metaSub$PATNO, metaSub$GenDiseaseStatus) > 0) #1035 patnos total
# 
# #order to match
# i <- match(metaSub$HudAlphaSampleName, names(countsSub), nomatch = 0)
# countsSub <- countsSub[,i]
# 
# identical(names(countsSub), metaSub$HudAlphaSampleName) #double checking
# 
# write_csv(metaSub, paste(outdir, "differentialExpression_sampleData.csv", sep = ""))


#CAN START HERE AND LOAD SAMPLE DATA TABLE NOW THAT IT EXISTS
metaSub <- read_csv(paste(outdir, "differentialExpression_sampleData.csv", sep = ""))
#order to match
countsSub <- countTable[,names(countTable) %in% sort(metaSub$HudAlphaSampleName)]

i <- match(metaSub$HudAlphaSampleName, names(countsSub), nomatch = 0)
countsSub <- countsSub[,i]

identical(names(countsSub), metaSub$HudAlphaSampleName) #double checking


## sample design --------------------------------------------------------------------------

visit <- factor(metaSub$CLINICAL_EVENT, levels = c("BL", "V02", "V04", "V06", "V08"))
subject <- metaSub$PATNO
disease <- factor(gsub("\\+", "", metaSub$GenDiseaseStatus), levels = c("GBA_PD", "GBA_Unaffected",
                                                                        "Healthy_Control", "Idiopathic_PD",
                                                                        "LRRK2_PD", "LRRK2_Unaffected"))
group <- factor(paste(disease, visit, sep = "_"))
age <- factor(metaSub$ageBin, levels = c("under_55", "55_to_65", "over_65"))
sex <- factor(metaSub$SEX,
              levels = c("Female", "Male"))
plate <- factor(metaSub$Plate)
usableBases <- metaSub$PCT_USABLE_BASES

design <- model.matrix(~ 0 + disease + sex + plate + age + usableBases)
colnames(design) <- gsub("disease", "", colnames(design))
colnames(design)


contr.matrix <- makeContrasts(
  GBAPD_GBAUn = GBA_PD - GBA_Unaffected,
  LRRK2PD_LRRK2Un = LRRK2_PD - LRRK2_Unaffected,
  GBAPD_LRRK2PD = GBA_PD - LRRK2_PD,
  GBAUn_LRRK2Un = GBA_Unaffected - LRRK2_Unaffected,
  GBAPD_IdioPD = GBA_PD - Idiopathic_PD,
  LRRK2PD_IdioPD = LRRK2_PD - Idiopathic_PD,
  GBAUn_HC = GBA_Unaffected - Healthy_Control,
  LRRK2Un_HC = LRRK2_Unaffected - Healthy_Control,
  IdioPD_HC = Idiopathic_PD - Healthy_Control,
  levels = colnames(design)
)

contr.matrix


## normalization and filtering in limma --------------------------------------------------------------------------

keep <- filterByExpr(countsSub, group = disease)
sum(keep) #26770 genes

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
GBA.PD.Group <- metaSub$HudAlphaSampleName[metaSub$GenDiseaseStatus == "GBA+_PD"]
GBA.Un.Group <- metaSub$HudAlphaSampleName[metaSub$GenDiseaseStatus == "GBA+_Unaffected"]
LRRK2.PD.Group <- metaSub$HudAlphaSampleName[metaSub$GenDiseaseStatus == "LRRK2+_Unaffected"]
LRRK2.Un.Group <- metaSub$HudAlphaSampleName[metaSub$GenDiseaseStatus == "LRRK2+_PD"]
idio.PD.Group <- metaSub$HudAlphaSampleName[metaSub$GenDiseaseStatus == "Idiopathic_PD"]
control.Group <- metaSub$HudAlphaSampleName[metaSub$GenDiseaseStatus == "Healthy_Control"]

GBA.PD.CPM <- logCPM[,GBA.PD.Group]
GBA.Un.CPM <- logCPM[,GBA.Un.Group]
LRRK2.PD.CPM <- logCPM[,LRRK2.PD.Group]
LRRK2.Un.CPM <- logCPM[,LRRK2.Un.Group]
idioCPM <- logCPM[,idio.PD.Group]
controlCPM <- logCPM[,control.Group]

meanCPM <- data.frame("gene_id" = row.names(logCPM),
                      "AvgExpr_Healthy_Control" = rowMeans(controlCPM),
                      "AvgExpr_Idiopathic_PD" = rowMeans(idioCPM),
                      "AvgExpr_GBA_Unaffected" = rowMeans(GBA.Un.CPM),
                      "AvgExpr_GBA_PD" = rowMeans(GBA.PD.CPM),
                      "AvgExpr_LRRK2_Unaffected" = rowMeans(LRRK2.Un.CPM),
                      "AvgExpr_LRRK2_PD" = rowMeans(LRRK2.PD.CPM)

)


#add genes to voom object
genes.anno <- genes.anno %>%
  mutate(ens_base = gsub("\\.(.*)", "", gene_id))
v$genes <- left_join(data.frame(gene_id = row.names(CPM)), genes.anno)




#output DE results ------------------------------------------------------------------------------------------------

resultsOut(efit, contr.matrix, outdir, logfc = 0.1, pval = 0.05)


