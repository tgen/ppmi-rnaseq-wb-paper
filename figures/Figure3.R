# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: May 14th, 2020

#Purpose: summarize DE results

#Imports:
# 


#Exports:
# Figure 3

# load packages -------------------------------------------------------------------------------------------------
library(tidyverse) #basis for data manipulation
library(ggpubr) #publication ready plots
library(UpSetR) #venn diagram alternative
library(ggrepel) #adds non-overlapping text to ggplot
library(EnhancedVolcano) #enhanced volcano plots
library(ComplexHeatmap) #for matrix function going into upset plot
library(devtools) #load R scripts from github

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

#color palettes
myPalette <- getPalette(21, "multi")
dev.off()
names(myPalette) <- NULL

diseasePalette <- myPalette[c(1:7)]
names(diseasePalette) <- c("Healthy Control", "Genetic Unaffected",
                           "Idiopathic PD", "Genetic PD",
                           "Prodromal", "SWEDD", "Other ND")

genDiseasePalette <- myPalette[c(1:15)]
names(genDiseasePalette) <- c(names(diseasePalette),
                              "GBA+ Affected", "GBA+ Unaffected",
                              "LRRK2+ Affected", "LRRK2+ Unaffected",
                              "Prodromal", "GBA+ Prodromal",
                              "SWEDD", "GBA+ SWEDD")


# import gene names -------------------------------------------------------------------------------------------------
tmpG <- read_tsv("sourceFiles/GRCh38_GENCODE29_geneInfo.txt",
                 skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type")]
genes.anno <- unique(genes.anno)


#load metadata  ----------------------------------------------------------------------------------------------
meta <- read_csv("data/megaMetaTable.csv")
meta <- subset(meta, QCflagIR3 == "pass")
names(meta)

meta$GenDisease_Status <- paste(meta$Genetic_Status, meta$Disease_Status, sep = " ")
meta$GenDisease_Status <- gsub("Genetic PD", "Affected", meta$GenDisease_Status)
meta$GenDisease_Status <- gsub("Genetic ", "", meta$GenDisease_Status)
meta$GenDisease_Status <- gsub("LRRK2-/SNCA-/GBA- ", "", meta$GenDisease_Status)
unique(meta$GenDisease_Status)

#load DE lists ------------------------------------------------------------------------------------------------
#grab file names
fileNames <- list.files(path = "analysis/DE/GBA_LRRK2_idio_sizeMatched", full.names = TRUE, recursive = TRUE, pattern = "_full.tsv$")
fileNames <- append(fileNames[c(1, 6, 8)], "analysis/DE/allPD_vs_allCon_sizeMatched/PD_full.tsv")
fileNames <- append(fileNames, "analysis/DE/HC_prodromal_SWEDD/Prodromal_HC_full.tsv")
fileNames <- append(fileNames, "analysis/DE/HC_prodromal_SWEDD/SWEDD_HC_full.tsv")
fileNames <- append(fileNames, "data/PDBP_case_vs_control_limmaResults.tsv")

#read in files
files <- lapply(fileNames, function(x) read_tsv(x))
files <- lapply(files, function(x) right_join(genes.anno, x))

head(files[[7]])

compName <- sapply(str_split(fileNames, "/"), '[', 4)
compName <- gsub("_full.tsv", "", compName)
compName <- gsub(".tsv", "", compName)
compName[4] <- "PPMI_All"
compName[7] <- "PDBP_All"




#make volcano plots ------------------------------------------------------------------------------------------------
#output volcano plot
plotTitles <- gsub("_", " versus ", compName)
plotTitles <- gsub("GBA", "GBA+ ", plotTitles)
plotTitles <- gsub("LRRK2", "LRRK2+ ", plotTitles)
plotTitles <- gsub("Un", "Unaffected", plotTitles)
plotTitles <- gsub("Idio", "Idiopathic ", plotTitles)
plotTitles <- gsub("HC", "Healthy Control", plotTitles)
plotTitles <- gsub("^PD$", "All PD versus All Unaffected", plotTitles)
plotTitles[7] <- "PDBP Cohort: Case versus Control"
plotTitles[4] <- "PPMI Cohort: Case versus Control"
plotTitles

createVolcano <- function(df, plotTitle, geneLabels) {
  EnhancedVolcano(df,
                  lab = df$gene_name,
                  selectLab = geneLabels,
                  boxedLabels = TRUE,
                  labFace = "bold",
                  drawConnectors = TRUE,
                  x = 'logFC',
                  y = 'adj.P.Val',
                  title = plotTitle,
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  legendLabels = c('NS', bquote(~Log[2]~ 'fold change'), 'Adjusted p-value',
                                   bquote('Adjusted p-value and ' ~Log[2]~ 'fold change')),
                  pCutoff = 0.05,
                  FCcutoff = 0.1,
                  pointSize = 2.5,
                  labSize = 3,
                  legendPosition = "bottom") +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_blank())
}



volPlots <- list()

dir.create("Figures/volcanoPlots")

x <- 1
while(x <= length(plotTitles)){
  
  df <- files[[x]] %>%
    left_join( myGenes) %>%
    mutate(category = factor(category, levels = c("other", "PD_risk")))
  
  df$category[is.na(df$category)] <- "other"
  
  df <- arrange(df, category)
  
  riskgenes <- df %>%
    filter(category == "PD_risk" & abs(logFC) > 0.1 & adj.P.Val < 0.05) %>%
    transmute(gene_name = gene_name)
    
  
  volPlot <- createVolcano(files[[x]], plotTitles[[x]], "")
  
  volPlots[[x]] <- volPlot
  
  ggsave(paste("Figures/volcanoPlots/", compName[[x]], "_volcano.svg", sep = ""),
         volPlot,
         dpi = 600, width = 8, height = 8)
  
  x <- x + 1
}



# volcano plots with the PD risk genes: Nails et al ---------------------------------------------------------------------------------
# create plots with PD risk genes labeled
riskPlot <- read_tsv("data/Nalls_2019_genes_manhattan.txt", col_names = "gene_name")
plotGenes <- unique(riskPlot$gene_name)
plotGenes <- append(plotGenes, c("GBAP1"))

createVolcanoCustom <- function(df, plotTitle, keyvals.colour, geneLabels) {
  EnhancedVolcano(df,
                  lab = df$gene_name,
                  selectLab = geneLabels,
                  boxedLabels = TRUE,
                  drawConnectors = TRUE,
                  x = 'logFC',
                  y = 'adj.P.Val',
                  title = plotTitle,
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  colCustom = keyvals.colour,
                  #shapeCustom = keyvals.shape,
                  colAlpha = 0.6,
                  pCutoff = 0.05,
                  FCcutoff = 0.1,
                  pointSize = 2.5,
                  labSize = 3,
                  legendPosition = "bottom") +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_blank())
}

volRiskPlots <- list()

myGenes <- data.frame(gene_name = plotGenes,
                      category = "PD_risk")

dir.create("Figures/PD_risk")

PPMI <- files[[4]] %>%
  left_join( myGenes) %>%
  mutate(category = factor(category, levels = c("other", "PD_risk")))
PPMI$category[is.na(PPMI$category)] <- "other"
PPMI <- arrange(PPMI, category)

riskgenes <- PPMI %>% filter(category == "PD_risk" & abs(logFC) > 0.1 & adj.P.Val < 0.05) %>%
   transmute(gene_name = gene_name)
geneLabels <- riskgenes$gene_name

x <- 1
while(x <= length(compName)){
  
  df <- files[[x]] %>%
    left_join( myGenes) %>%
    mutate(category = factor(category, levels = c("other", "PD_risk")))
  
  df$category[is.na(df$category)] <- "other"
  
  df <- arrange(df, category)
  
  riskgenes <- df %>%
    filter(category == "PD_risk" & abs(logFC) > 0.1 & adj.P.Val < 0.05) %>%
    transmute(gene_name = gene_name)

  #custom key values for enhanced volcano plot
  keyvals <- rep('gray60', nrow(df))
  names(keyvals) <- rep('other', nrow(df))
  
  # modify keyvals for genes in PD risk category
  keyvals[which(df$category == "PD_risk")] <- '#440154FF'
  names(keyvals)[which(df$category == "PD_risk")] <- 'PD_risk'
  
  unique(names(keyvals))
  
  volPlot <- createVolcanoCustom(df, plotTitles[[x]], keyvals, geneLabels)
  
  volPlot <- volPlot +
    geom_rug(data = df %>%
               filter(category == "PD_risk" & logFC > 0.1 & adj.P.Val < 0.05),
             aes(y = -log10(adj.P.Val), color = category),
             inherit.aes = FALSE,
             sides = "r", outside = FALSE) +
    geom_rug(data = df %>%
               filter(category == "PD_risk" & logFC < -0.1 & adj.P.Val < 0.05),
             aes(y = -log10(adj.P.Val), color = category),
             inherit.aes = FALSE,
             sides = "l", outside = FALSE)
  
  volRiskPlots[[x]] <- volPlot
  
  ggsave(paste("Figures/volcanoPlots/", compName[[x]], "_volcano_PDrisk.svg", sep = ""),
         volPlot,
         dpi = 600, width = 8, height = 8)
  
  x <- x + 1
}

