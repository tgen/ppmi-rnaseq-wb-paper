# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: April 1st, 2020

#Purpose: summarize DE results at the pathway level: neutrophil activation/response

#Imports:
# gene lists from limma

#Exports:
# 


# load packages -------------------------------------------------------------------------------------------------
library(tidyverse) #basis for data manipulation
library(limma) #for camera gene set testing
library(ggpubr) #publication ready plots
library(pheatmap) #make pretty heatmaps
library(dendsort) #sorts dendrograms
library(EnhancedVolcano) #enhanced volcano plots
library(corrr) #correlation functions
library(VennDiagram) #create venn diagrams
library(devtools) #load R scripts from github

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

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

genDiseaseGroups <- c("Healthy Control","GBA+ Unaffected", "LRRK2+ Unaffected",
                      "Idiopathic PD", "GBA+ Affected", "LRRK2+ Affected")

createVolcano <- function(df, plotTitle, keyvals.colour) {
  EnhancedVolcano(df,
                  lab = df$gene_id,
                  #selectLab = c("ENSG00000100985.7"),
                  selectLab = c(""),
                  x = 'logFC',
                  y = 'adj.P.Val',
                  title = plotTitle,
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  colCustom = keyvals.colour,
                  colAlpha = 0.7,
                  pCutoff = 0.05,
                  FCcutoff = 0.1,
                  pointSize = 4,
                  labSize = 3,
                  titleLabSize = 22,
                  legendPosition = "bottom") +
    theme(plot.title = element_text(size = 20),
          plot.subtitle = element_blank())
}

## load metadata table ----------------------------------------------------------------------------------------------

meta <- read_csv("data/megaMetaTable.csv")
meta <- subset(meta, QCflagIR3 == "pass")
names(meta)

meta <- meta %>%
  filter(is.na(PoolAssign))

meta$Disease_Type <- paste(meta$Genetic_Status, meta$Disease_Status, sep = " ")
meta$Disease_Type <- gsub("Genetic PD", "Affected", meta$Disease_Type)
meta$Disease_Type <- gsub("Genetic ", "", meta$Disease_Type)
meta$Disease_Type <- gsub("LRRK2-/SNCA-/GBA- ", "", meta$Disease_Type)
meta$Disease_Type <- gsub("NA Healthy Control", "Healthy Control", meta$Disease_Type)
meta$Disease_Type <- gsub("NA Idiopathic PD", "Idiopathic PD", meta$Disease_Type)
meta$Disease_Type <- gsub("LRRK2\\+ Other ND", "Other ND", meta$Disease_Type)
meta$Disease_Type <- gsub("GBA\\+ Other ND", "Other ND", meta$Disease_Type)
meta$Disease_Type <- gsub("GBA\\+ SWEDD", "SWEDD", meta$Disease_Type)
meta$Disease_Type <- gsub("NA SWEDD", "SWEDD", meta$Disease_Type)
meta$Disease_Type <- gsub("SNCA\\+ SWEDD", "SWEDD", meta$Disease_Type)
meta$Disease_Type <- gsub("GBA\\+ Prodromal", "Prodromal", meta$Disease_Type)
unique(meta$Disease_Type)

meta <- meta %>%
  filter(Disease_Type != "NA Unaffected")

metaOut <- meta %>%
  dplyr::select(PATNO, VISIT, Disease_Type)
write_tsv(metaOut, "Figures/sankey_input.tsv")

metaOutTally <- metaOut %>% group_by(VISIT, Disease_Type) %>%
  tally()
write_tsv(metaOutTally, "Figures/sankey_input_tallied.tsv")

metaByPATNO <- meta[,c("PATNO", "ageBin", "SEX", "age_at_consent", "Disease_Type")]
metaByPATNO <- unique(metaByPATNO)
metaByPATNO$Disease_Type <- factor(metaByPATNO$Disease_Type, levels = genDiseaseGroups)
metaByPATNO$ageBin <- factor(metaByPATNO$ageBin, levels = c("under_55", "55_to_65", "over_65"))

# blood count data plots --------------------------------------------
#add blood cell types
bloodchem <- read_csv("data/bloodChem_sequenced.csv") #this includes all samples from patnos that we sequenced
bloodchem$LSIRES <- as.numeric(bloodchem$LSIRES)
bloodchem$Disease_Status <- gsub("Parkinson's Disease", "Idiopathic PD", bloodchem$Disease_Status)
bloodchem.full <- bloodchem
bloodchem.full$VISIT <- sapply(strsplit(bloodchem.full$PATNO_VISIT, "_"), "[", 2)
bloodchem.full$visitTime <- bloodchem.full$VISIT
bloodchem.full$visitTime <- gsub("SC", "0", bloodchem.full$visitTime)
bloodchem.full$visitTime <- gsub("V04", "12", bloodchem.full$visitTime)
bloodchem.full$visitTime <- gsub("V06", "24", bloodchem.full$visitTime)
bloodchem.full$visitTime <- gsub("V08", "36", bloodchem.full$visitTime)
bloodchem.full$visitTime <- gsub("V10", "48", bloodchem.full$visitTime)
bloodchem.full$visitTime <- gsub("V12", "60", bloodchem.full$visitTime)
bloodchem.full$visitTime <- as.numeric(bloodchem.full$visitTime)
bloodchem.full <- left_join(bloodchem.full, metaByPATNO) %>% na.omit()
bloodchem.full$Disease_Type <- factor(bloodchem.full$Disease_Type, levels = c("Healthy Control", "Prodromal", "Idiopathic PD",
                                                                              "GBA+ Unaffected", "GBA+ Affected",
                                                                              "LRRK2+ Unaffected", "LRRK2+ Affected"))
bloodchem.full$Case <- "Control"
bloodchem.full$Case[bloodchem.full$Disease_Status == "Idiopathic PD"] <- "Case"
bloodchem.full$Case[bloodchem.full$Disease_Status == "Genetic PD"] <- "Case"
bloodchem.full$Case[bloodchem.full$Disease_Status == "Prodromal"] <- "Prodromal"
bloodchem.full$Case <- factor(bloodchem.full$Case, levels = c("Control", "Case", "Prodromal"))

my_comparisons = list(c("Healthy Control", "Prodromal"),
                      c("Prodromal", "Idiopathic PD"),
                      c("Healthy Control", "Idiopathic PD"),
                      c("GBA+ Unaffected", "GBA+ Affected"),
                      c("LRRK2+ Unaffected", "LRRK2+ Affected")
)
lymphPlot <- ggplot(bloodchem.full %>%
                      filter(LSIUNIT == "%") %>%
                      filter(LTSTNAME == "Lymphocytes (%)") %>%
                      filter(VISIT == "SC"),
                    aes(x = Disease_Type, y = LSIRES)) +
  geom_jitter(alpha = 0.3, width = 0.1,
              aes(color = Disease_Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0, lwd = 1) +
  ggtitle("Lymphocyte Percentage") +
  ylab("percent") +
  scale_color_manual(values = genDiseasePalette) +
  stat_compare_means(comparisons = my_comparisons) +
  xlab("") +
  scale_x_discrete(labels=c("Healthy Control" = "Healthy Control \n (n = 195)",
                            "Prodromal" = "Prodromal \n (n = 58)",
                            "Idiopathic PD" = "Idiopathic PD \n (n = 349)",
                            "GBA+ Unaffected" = "GBA+ Unaffected \n (n = 163)",
                            "GBA+ Affected" = "GBA+ Affected \n (n = 111)",
                            "LRRK2+ Unaffected" = "LRRK2+ Unffected \n (n = 223)",
                            "LRRK2+ Affected" = "LRRK2+ Affected \n n = 206")) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size = 16, face = "bold"),
        legend.position = "none")



neutPlot <- ggplot(bloodchem.full %>%
                     filter(LSIUNIT == "%") %>%
                     filter(LTSTNAME == "Neutrophils (%)") %>%
                     filter(VISIT == "SC"),
                   aes(x = Disease_Type, y = LSIRES)) +
  geom_jitter(alpha = 0.3, width = 0.1,
              aes(color = Disease_Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0, lwd = 1) +
  ggtitle("Neutrophil Percentage") +
  ylab("percent") +
  scale_color_manual(values = genDiseasePalette) +
  stat_compare_means(comparisons = my_comparisons) +
  xlab("") +
  scale_x_discrete(labels=c("Healthy Control" = "Healthy Control \n (n = 195)",
                            "Prodromal" = "Prodromal \n (n = 58)",
                            "Idiopathic PD" = "Idiopathic PD \n (n = 349)",
                            "GBA+ Unaffected" = "GBA+ Unaffected \n (n = 163)",
                            "GBA+ Affected" = "GBA+ Affected \n (n = 111)",
                            "LRRK2+ Unaffected" = "LRRK2+ Unffected \n (n = 223)",
                            "LRRK2+ Affected" = "LRRK2+ Affected \n n = 206")) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size = 16, face = "bold"),
        legend.position = "none")

ggarrange(neutPlot, lymphPlot,
          ncol = 2, nrow = 1,
          labels = "auto",
          font.label = list(size = 24))
ggsave("Figures/Figure5-AB.png", width = 14, height = 7)


# same plot, for all case vs all controls ---------------
allNeut <- ggplot(bloodchem.full %>%
                    filter(LSIUNIT == "%") %>%
                    filter(LTSTNAME == "Neutrophils (%)") %>%
                    filter(VISIT == "SC") %>%
                    filter(Case != "Prodromal"),
                  aes(x = Case, y = LSIRES)) +
  geom_jitter(alpha = 0.3, width = 0.1,
              aes(color = Case)) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0, lwd = 1) +
  ggtitle("Neutrophil Percentage") +
  ylab("percent") +
  scale_color_manual(values = c("#BE0032", "#0067A5")) +
  stat_compare_means(comparisons = list(c("Case", "Control"))) +
  xlab("") +
  scale_x_discrete(labels=c("Case" = "Case \n (n = 666)",
                            "Control" = "Control \n (n = 581)")) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size = 16, face = "bold"),
        legend.position = "none")

allLymph <- ggplot(bloodchem.full %>%
                     filter(LSIUNIT == "%") %>%
                     filter(LTSTNAME == "Lymphocytes (%)") %>%
                     filter(VISIT == "SC") %>%
                     filter(Case != "Prodromal"),
                   aes(x = Case, y = LSIRES)) +
  geom_jitter(alpha = 0.3, width = 0.1,
              aes(color = Case)) +
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0, lwd = 1) +
  ggtitle("Lymphocyte Percentage") +
  ylab("percent") +
  scale_color_manual(values = c("#BE0032", "#0067A5")) +
  stat_compare_means(comparisons = list(c("Case", "Control"))) +
  xlab("") +
  scale_x_discrete(labels=c("Case" = "Case \n (n = 666)",
                            "Control" = "Control \n (n = 581)")) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size = 16, face = "bold"),
        legend.position = "none")

ggarrange(allNeut, allLymph,
          ncol = 2, nrow = 1,
          labels = "auto",
          font.label = list(size = 24))
ggsave("Figures/Figure4-AB.png", width = 14, height = 7)


#load gene annotations ----------------------------------------------------------------------------------------------
tmpG <- read_tsv("sourceFiles/GRCh38_GENCODE29_geneInfo.txt",
                 skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type")]
genes.anno <- unique(genes.anno)
genes.anno$gene_IDBase <- gsub("\\.(.*)$", "", genes.anno$gene_id)

# load neutrophil enriched genes from the human protein atlas  ----------------------------------------------------
neutEnriched <- read_tsv("data/blood_cell_category_rna_neutrophil_Cell_enriched.tsv")

BcellEnriched <- read_tsv("data/blood_cell_category_rna_naive_bcell_enriched.tsv")
TcellEnriched <- read_tsv("data/blood_cell_category_rna_T-reg_Cell_enriched.tsv")
lymphEnriched <- rbind(BcellEnriched, TcellEnriched)

# create volcanos for neutrophil enhanced genes --------------------------------------------------------------------------------------------

#grab file names
fileNames <- list.files(path = "data/", full.names = TRUE, recursive = TRUE, pattern = "_full.tsv$")
fileNames <- append(fileNames[c(1, 6, 8)], "data/PD_full.tsv")
fileNames <- append(fileNames, "data/Prodromal_HC_full.tsv")
fileNames <- append(fileNames, "data/SWEDD_HC_full.tsv")
fileNames <- append(fileNames, "data/PDBP_case_vs_control_limmaResults.tsv")

#read in files
files <- lapply(fileNames, function(x) read_tsv(x))
files <- lapply(files, function(x) right_join(genes.anno, x))

head(files[[7]])

compName <- sapply(str_split(fileNames, "/"), '[', 4)
compName <- gsub("_full.tsv", "", compName)
compName <- gsub(".tsv", "", compName)
compName[7] <- "Case_Control"

neutEnrichedGenes <- data.frame(gene_name = as.character(neutEnriched$Gene),
                                category = "Neutrophil Enriched",
                                cellType = "Neutrophil",
                                group = "Enriched")
lymphEnrichedGenes <- data.frame(gene_name = as.character(lymphEnriched$Gene),
                                 category = "Lymphocyte Enriched",
                                 cellType = "Lymphocyte",
                                 group = "Enriched")
myGenes <- rbind(neutEnrichedGenes, lymphEnrichedGenes)

PPMI.myGenes <- myGenes[c(1,2)]  %>%
  mutate(cohort = "PPMI") %>%
  inner_join(files[[4]]) %>%
  mutate(AvgExpr_Control = AvgExpr_allControls,
         AvgExpr_Case = AvgExpr_allPD) %>%
  select(-AvgExpr_allPD, -AvgExpr_allControls)

PDBP.myGenes <- myGenes[c(1,2)]  %>%
  mutate(cohort = "PDBP") %>%
  inner_join(files[[7]])

myGenes.out <- rbind(PPMI.myGenes,
                     PDBP.myGenes)
write_tsv(myGenes.out, "Tables/Supplemental_Table3.tsv")


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

volPlots <- list() #initiate list to hold volcano plots
neutBars <- list()
lymphBars <- list()


dir.create("Figures/bloodCell")

x <- 1
while(x <= length(compName)){
  
  df <- files[[x]] %>%
    left_join( myGenes) %>%
    mutate(category = factor(category, levels = c("Other",
                                                  "Neutrophil Enhanced", "Neutrophil Enriched",
                                                  "Lymphocyte Enhanced", "Lymphocyte Enriched")),
           cellType = factor(cellType, levels = c("Other", "Neutrophil", "Lymphocyte")),
           group = factor(group, levels = c("Other", "Enhanced", "Enriched")))
  
  df$category[is.na(df$category)] <- "Other"
  df$cellType[is.na(df$cellType)] <- "Other"
  df$group[is.na(df$group)] <- "Other"
  
  df <- arrange(df, category)
  
  #custom key values for enhanced volcano plot
  keyvals <- rep('gray60', nrow(df))
  names(keyvals) <- rep('Other', nrow(df))
  
  # modify keyvals for genes in enriched and enhanced categories
  keyvals[which(df$category == "Neutrophil Enriched")] <- '#0067A5'
  names(keyvals)[which(df$category == "Neutrophil Enriched")] <- 'Neutrophil Enriched'
  
  keyvals[which(df$category == "Lymphocyte Enriched")] <- '#F3C300'
  names(keyvals)[which(df$category == "Lymphocyte Enriched")] <- 'Lymphocyte Enriched'
  
  unique(names(keyvals))
  
  
  volPlot <- createVolcano(df, plotTitles[[x]], keyvals)
  volPlot <- volPlot +
    geom_rug(data = df %>%
               filter(category == "Neutrophil Enriched" | category == "Lymphocyte Enriched") %>%
               filter(logFC > 0),
             aes(y = -log10(adj.P.Val), color = category),
             inherit.aes = FALSE,
             sides = "r", outside = FALSE) +
    geom_rug(data = df %>%
               filter(category == "Neutrophil Enriched" | category == "Lymphocyte Enriched") %>%
               filter(logFC < 0),
             aes(y = -log10(adj.P.Val), color = category),
             inherit.aes = FALSE,
             sides = "l", outside = FALSE) +
    ylim(c(0, ceiling(max(-log10(df$adj.P.Val)))))
  ggsave(paste("Figures/bloodCell/", compName[[x]], "_volcano.png", sep = ""),
         volPlot,
         dpi = 600, width = 8, height = 8)

  geneList <- list(neutrophil_enriched = neutEnriched$Ensembl,
                   lymphocyte_enriched = lymphEnriched$Ensembl)
  idx <- ids2indices(geneList, id=df$gene_IDBase)
  
  png(paste("Figures/bloodCell/", compName[[x]], "_neutBarcode.png", sep = ""),
         width = 800, height = 400)
  barcodeplot(df$t, index = idx$neutrophil_enriched, main = paste(plotTitles[[x]], "\n", "Neutrophil Enriched", sep = ""))
  par(oma = c(2,2,2,2))
  neutBar <- recordPlot()
  dev.off()

  png(paste("Figures/bloodCell/", compName[[x]], "_lymphBarcode.png", sep = ""),
      width = 800, height = 400)
  barcodeplot(df$t, index = idx$lymphocyte_enriched, main = paste(plotTitles[[x]], "\n", "Lymphocyte Enriched", sep = ""))
  par(oma = c(2,2,2,2))
  lymphBar <- recordPlot()
  dev.off()
  
  
  volPlots[[x]] <- volPlot
  neutBars[[x]] <- neutBar
  lymphBars[[x]] <- lymphBar
  
  
  x <- x + 1
}





#venn diagram -----------
plotTitles
de.list.neut <- list(PPMI = files[[4]] %>%
                       filter(adj.P.Val < 0.05) %>%
                       filter(abs(logFC) > 0.1) %>%
                       inner_join(neutEnriched, by = c("gene_IDBase" = "Ensembl")) %>%
                       transmute(gene_name),
                     PDBP = files[[7]] %>%
                       filter(adj.P.Val < 0.05) %>%
                       filter(abs(logFC) > 0.1) %>%
                       inner_join(neutEnriched, by = c("gene_IDBase" = "Ensembl")) %>%
                       transmute(gene_name)
                     )

de.list.lymph <- list(PPMI = files[[4]] %>%
                       filter(adj.P.Val < 0.05) %>%
                       filter(abs(logFC) > 0.1) %>%
                       inner_join(lymphEnriched, by = c("gene_IDBase" = "Ensembl")) %>%
                       transmute(gene_name),
                     PDBP = files[[7]] %>%
                       filter(adj.P.Val < 0.05) %>%
                       filter(abs(logFC) > 0.1) %>%
                       inner_join(lymphEnriched, by = c("gene_IDBase" = "Ensembl")) %>%
                       transmute(gene_name)
)



dev.off()
svg(filename = "Figures/Figure4i.svg", height = 6, width = 6)
draw.pairwise.venn(category = c("PPMI" , "PDBP"),
                   area1 = length(de.list.neut[[1]]$gene_name),
                   area2 = length(de.list.neut[[2]]$gene_name),
                   cross.area = sum(de.list.neut[[1]]$gene_name %in% de.list.neut[[2]]$gene_name),
                   fill = c("#482677FF", "#2D708EFF"),
                   alpha = c(0.4, 0.4), cex = 2,
                   cat.fontface = 2, cat.cex = 2, cat.pos = c(135,230), cat.dist = 0.05)
dev.off()


dev.off()
svg(filename = "Figures/Figure4j.svg", height = 6, width = 6)
draw.pairwise.venn(category = c("PPMI" , "PDBP"),
                   area1 = length(de.list.lymph[[1]]$gene_name),
                   area2 = length(de.list.lymph[[2]]$gene_name),
                   cross.area = sum(de.list.lymph[[1]]$gene_name %in% de.list.lymph[[2]]$gene_name),
                   fill = c("#482677FF", "#2D708EFF"),
                   alpha = c(0.4, 0.4), cex = 2,
                   cat.fontface = 2, cat.cex = 2, cat.pos = c(180,180))
dev.off()