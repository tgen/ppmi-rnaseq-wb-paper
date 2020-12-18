# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: October 9th, 2019

#Purpose: Create Figure 6

#Imports:
# PPMI count table (filtered)
# cell count data


#Exports:
# Figure 6


## load libraries  ----------------------------------------------------------------------------------------------
library(tidyverse) #42; base for data/table manipulation
library(edgeR) #expression normalization
library(ggpubr) #paper level visuals
library(reshape)
library(devtools) #load R scripts from github

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

`%notin%` <- Negate(`%in%`)


myPalette <- getPalette(21, "multi")
dev.off()
names(myPalette) <- NULL
myPalette <- myPalette[c(1,3,5)]
names(myPalette) <- c("Control", "Case", "Prodromal")


## load count table  ----------------------------------------------------------------------------------------------
countTable <- read_csv("data/PPMI_Subjects_GeneCountsB38_Filtered.csv")
countTable <- column_to_rownames(countTable, var = "Ensembl_ID")


## biotypes and gene names ----------------------------------------------------------------------------------------
tmpG <- read_tsv("sourceFiles/GRCh38_GENCODE29_geneInfo.txt",
                 skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type")]
genes.anno <- unique(genes.anno)

genes.anno$Ensembl <- gsub("\\.(.*)$", "", genes.anno$gene_id)


## load metadata  ---------------------------------------------------------------------------------------------------
meta <- read_csv("data/megaMetaTable.csv")
dim(meta)

meta <- meta[meta$QCflagIR3 == "pass",] #filter out fails
dim(meta) #4756

meta <- meta[!meta$PoolAssign %in% c("PD_POOL", "HC_POOL"),] #filter out pools
dim(meta) #4650

meta$Disease_Status <- factor(meta$Disease_Status,
                              levels = c("Healthy Control", "Genetic Unaffected",
                                         "Idiopathic PD", "Genetic PD",
                                         "Prodromal", "SWEDD", "Other ND")
)

meta$Genetic_Status <- factor(meta$Genetic_Status,
                              levels = c("LRRK2-/SNCA-/GBA-",
                                         "LRRK2+", "SNCA+", "GBA+"))

meta$GenDisease_Status <- paste(meta$Genetic_Status, meta$Disease_Status, sep = " ")
meta$GenDisease_Status <- gsub("Genetic PD", "Affected", meta$GenDisease_Status)
meta$GenDisease_Status <- gsub("Genetic ", "", meta$GenDisease_Status)
meta$GenDisease_Status <- gsub("LRRK2-/SNCA-/GBA- ", "", meta$GenDisease_Status)
unique(meta$GenDisease_Status)

meta.demo.byPATNO <- meta[,c("PATNO", "SEX", "ageBin", "Disease_Status", "Genetic_Status", "GenDisease_Status")] %>%
  unique()
meta.demo.byPATNO$PATNO <- as.character(meta.demo.byPATNO$PATNO)
meta.demo.byPATNO$group <- NA
meta.demo.byPATNO$group[meta.demo.byPATNO$Disease_Status == "Healthy Control"] <- "Control"
meta.demo.byPATNO$group[meta.demo.byPATNO$Disease_Status == "Genetic Unaffected"] <- "Control"
meta.demo.byPATNO$group[meta.demo.byPATNO$Disease_Status == "Genetic PD"] <- "Case"
meta.demo.byPATNO$group[meta.demo.byPATNO$Disease_Status == "Idiopathic PD"] <- "Case"
meta.demo.byPATNO$group[meta.demo.byPATNO$Disease_Status == "Prodromal"] <- "Prodromal"
meta.demo.byPATNO$group[meta.demo.byPATNO$Disease_Status == "SWEDD"] <- "SWEDD"

diseaseGroups <- c("Healthy Control", "Idiopathic PD",
                   "GBA+ Unaffected", "GBA+ Affected",
                   "LRRK2+ Unaffected", "LRRK2+ Affected",
                   "Prodromal", "SWEDD")


pro_to_pd <- meta %>%
  filter(Disease_Status == "Prodromal") %>%
  filter(Most_likely_primary_diagnosis == "Idiopathic PD") %>%
  transmute(PATNO = PATNO) %>%
  unique()

pro_to_pd.visits <- meta %>% filter(PATNO %in% pro_to_pd$PATNO) %>%
  transmute(PATNO = PATNO,
            PATNO_VISIT = PATNO_VISIT,
            PRIMDAG = Most_likely_primary_diagnosis)

pro.exclude <- c("14426_V02", "14426_V04", "14426_V06", "14426_V08", "14426_V10", "14426_V12",
                 "18567_V04", "18567_V06", "18567_V08", "18567_V10", "18567_V12",
                 "60013_V04", "60013_V06", "60013_V08", "60013_V10", "60013_V12",
                 "60024_V08", "60024_V10", "60024_V12",
                 "60065_V06", "60065_V08", "60065_V10", "60065_V12",
                 "90456_V08", "90456_V10", "90456_V12")

## normalization and filtering in edgeR --------------------------------------------------------------------------

# identify genes that pass expression cutoff
geneCounts <- countTable
isexpr <- rowSums(cpm(geneCounts)>0.5) >= 0.25 * ncol(geneCounts)

# create data structure with only expressed genes
gExpr <- DGEList(counts=geneCounts[isexpr,])

# Perform TMM normalization
gExpr <- calcNormFactors(gExpr)

logCPM <- cpm(gExpr, log = TRUE, prior.count = 3)
CPM <- cpm(gExpr)

## subset to neutrophil genes --------------------------------------------------------------------------
neutGenes <- read_tsv("data/blood_cell_category_rna_neutrophil_Cell_enriched.tsv")
neutGenes <- neutGenes[,c("Gene", "Ensembl", "RNA blood cell specific NX")]

neutGenesList <- genes.anno[genes.anno$Ensembl %in% neutGenes$Ensembl,]

genes.expr <- CPM[row.names(CPM) %in% neutGenesList$gene_id, ]
genes.expr <- left_join(rownames_to_column(as.data.frame(t(genes.expr)), var = "HudAlphaSampleName"), meta[,c("HudAlphaSampleName", "PATNO", "VISIT", "Disease_Status", "Genetic_Status", "GenDisease_Status")])
genes.expr <- genes.expr[genes.expr$Disease_Status %in% c("Healthy Control", "Genetic Unaffected",
                                                          "Idiopathic PD", "Genetic PD",
                                                          "Prodromal", "SWEDD"),]
genes.expr$PATNO <- as.character(genes.expr$PATNO)
genes.expr.melted <- melt(genes.expr)
genes.expr.melted <- left_join(genes.expr.melted, genes.anno, by = c("variable" = "gene_id")) %>%
  mutate(PATNO_VISIT = paste(PATNO, VISIT, sep = "_"))
genes.expr.melted <- genes.expr.melted[!is.na(genes.expr.melted$Genetic_Status),]
genes.expr.melted$Genetic_Status <- as.factor(genes.expr.melted$Genetic_Status)


# DE neut genes ---------
DE.genes <- read_tsv("analysis/DE/allPD_vs_allCon_sizeMatched/PD_full.tsv")
DE.genes.neut.up <- DE.genes %>%
  filter(logFC > 0.1) %>%
  filter(adj.P.Val < 0.05) %>%
  filter(gene_name %in% neutGenes$Gene)
genes.expr.melted.neutUP <- genes.expr.melted %>%
  filter(variable %in% DE.genes.neut.up$gene_id) %>%
  left_join(meta.demo.byPATNO) %>%
  filter(group %in% c("Case", "Control", "Prodromal")) %>%
  filter(ageBin != "under_55") %>%
  filter(PATNO_VISIT %notin% pro.exclude)
genes.expr.melted.neutUP$visitTime <- genes.expr.melted.neutUP$VISIT
genes.expr.melted.neutUP$visitTime <- gsub("BL", "0", genes.expr.melted.neutUP$visitTime)
genes.expr.melted.neutUP$visitTime <- gsub("V02", "0.5", genes.expr.melted.neutUP$visitTime)
genes.expr.melted.neutUP$visitTime <- gsub("V04", "1", genes.expr.melted.neutUP$visitTime)
genes.expr.melted.neutUP$visitTime <- gsub("V06", "2", genes.expr.melted.neutUP$visitTime)
genes.expr.melted.neutUP$visitTime <- gsub("V08", "3", genes.expr.melted.neutUP$visitTime)
genes.expr.melted.neutUP$visitTime <- as.numeric(genes.expr.melted.neutUP$visitTime)
genes.expr.melted.neutUP$logvalue <- log2(genes.expr.melted.neutUP$value + 3)
genes.expr.melted.neutUP$group <- factor(genes.expr.melted.neutUP$group, levels = c("Control", "Prodromal", "Case"))

#glm -----------
linModSum.patno <- genes.expr.melted.neutUP %>%
  group_by(PATNO) %>% 
  do({
    mod = glm(logvalue ~ visitTime, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  }) %>%
  left_join(meta.demo.byPATNO)
linModSum.patno$group <- factor(linModSum.patno$group, levels = c("Control", "Prodromal", "Case"))

linModSum.gene <- genes.expr.melted.neutUP %>%
  group_by(group, gene_name) %>% 
  do({
    mod = glm(logvalue ~ visitTime, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
linModSum.gene$group <- factor(linModSum.gene$group, levels = c("Control", "Prodromal", "Case"))

linModSum.gene.patno <- genes.expr.melted.neutUP %>%
  group_by(PATNO, gene_name) %>%
  do({
    mod = glm(logvalue ~ visitTime, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
linModSum.gene.patno <- left_join(linModSum.gene.patno, meta.demo.byPATNO)
linModSum.gene.patno$group <- factor(linModSum.gene.patno$group, levels = c("Control", "Prodromal", "Case"))

#glm -----------
my_comparisons <- list(c("Control", "Prodromal"),
                       c("Prodromal", "Case"),
                       c("Control", "Case"))

slope.neut.expr.byPATNO <- ggplot(linModSum.patno, aes(x = group, y = Slope, fill = group)) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  scale_fill_manual(values = myPalette) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = TRUE, vjust = 0.6) +
  theme_bw(base_size = 26) +
  ylim(c(-2, 3.5)) +
  ggtitle("") +
  labs(x = "", y = "slope of glm") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
slope.neut.expr.byPATNO$layers[[2]]$aes_params$textsize <- 12
slope.neut.expr.byPATNO


glm.expr.plot <- ggplot(genes.expr.melted.neutUP, aes(x = visitTime, y = log2(value + 1), color = group)) +
  geom_smooth(method = "glm", se = TRUE, alpha = 0.2) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 2, 3), labels = c("0" = 0, "0.5" = 0.5, "1" = 1, "2" = 2, "3" = 3)) +
  ggtitle("Neutrophil Expression") +
  scale_color_manual(values = myPalette) +
  theme_bw(base_size = 22) +
  labs(x = "Visit (years)", y = "log2(CPM + 1)") +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 18), axis.title = element_text(size = 18))
glm.expr.plot

# blood count data plots --------------------------------------------
#add blood cell types
bloodchem <- read_csv("data/bloodChem_sequenced.csv") #this includes all samples from patnos that we sequenced
bloodchem$LSIRES <- as.numeric(bloodchem$LSIRES)
bloodchem$Disease_Status <- gsub("Parkinson's Disease", "Idiopathic PD", bloodchem$Disease_Status)
bloodchem.full <- bloodchem
bloodchem.full$PATNO <- as.character(bloodchem.full$PATNO)
bloodchem.full$VISIT <- sapply(strsplit(bloodchem.full$PATNO_VISIT, "_"), "[", 2)
bloodchem.full$visitTime <- bloodchem.full$VISIT
bloodchem.full$visitTime <- gsub("SC", "0", bloodchem.full$visitTime)
bloodchem.full$visitTime <- gsub("V04", "1", bloodchem.full$visitTime)
bloodchem.full$visitTime <- gsub("V06", "2", bloodchem.full$visitTime)
bloodchem.full$visitTime <- gsub("V08", "3", bloodchem.full$visitTime)
bloodchem.full$visitTime <- gsub("V10", "4", bloodchem.full$visitTime)
bloodchem.full$visitTime <- gsub("V12", "5", bloodchem.full$visitTime)
bloodchem.full$visitTime <- as.numeric(bloodchem.full$visitTime)
bloodchem.full <- left_join(bloodchem.full, meta.demo.byPATNO) %>%
  na.omit() %>%
  filter(group %in% c("Case", "Control", "Prodromal"))
bloodchem.full$group <- factor(bloodchem.full$group, levels = c("Control", "Prodromal", "Case"))


bloodchem.full.filt <- bloodchem.full %>%
  filter(LSIUNIT == "%") %>%
  filter(LTSTNAME == "Neutrophils (%)") %>%
  filter(VISIT %in% c("SC", "V04", "V06", "V08", "V10", "V12")) %>%
  filter(ageBin != "under_55") %>%
  filter(PATNO_VISIT %notin% pro.exclude)




glm.cbc.plot <- ggplot(bloodchem.full.filt,
                   aes(x = visitTime, y = LSIRES, color = group)) +
  geom_smooth(method = "glm", se = TRUE, alpha = 0.2) +
  ggtitle("Neutrophil Percentage") +
  scale_color_manual(values = myPalette) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  xlab("") +
  theme_bw(base_size = 18) +
  labs(x = "Visit (years)", y = "percent") +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 18))
glm.cbc.plot


#glm for cell count data ------
linModSum.counts.patno <- bloodchem.full.filt %>%
  group_by(PATNO) %>% 
  do({
    mod = glm(LSIRES ~ visitTime, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  }) %>%
  left_join(meta.demo.byPATNO)
linModSum.counts.patno$group <- factor(linModSum.counts.patno$group, levels = c("Control", "Prodromal", "Case"))

slope.neut.counts.byPATNO <- ggplot(linModSum.counts.patno, aes(x = group, y = Slope, fill = group)) +
  geom_boxplot(alpha = 0.8, width = 0.5) +
  scale_fill_manual(values = myPalette) +
  #stat_compare_means(comparisons = my_comparisons) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = FALSE) +
  theme_bw(base_size = 26) +
  ggtitle("") +
  ylim(c(-15, 30)) +
  labs(x = "", y = "slope of glm") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
slope.neut.counts.byPATNO$layers[[2]]$aes_params$textsize <- 6
slope.neut.counts.byPATNO
ggsave("analysis/cellType_progression/prodromal_slope_overTime_cellCounts_byPATNO.png")


stable.neutCount <- linModSum.counts.patno %>%
  group_by(group) %>%
  summarise(Intercept = round(mean(Intercept, na.rm = TRUE), digits = 2),
            Slope = round(mean(Slope, na.rm = TRUE), digits = 2))
stable.neutCount


# individual gene expression --------------------------------------------

# Summary table plot
stable.FCGR2A <- linModSum.gene.patno %>%
  filter(gene_name == "FCGR2A") %>%
  group_by(group) %>%
  summarise(Intercept = round(mean(Intercept, na.rm = TRUE), digits = 2),
            Slope = round(mean(Slope, na.rm = TRUE), digits = 2))

compare_means(Intercept ~ group, linModSum.gene.patno %>% filter(gene_name == "FCGR2A"))
compare_means(Slope ~ group, linModSum.gene.patno %>% filter(gene_name == "FCGR2A"))

stable.FCGR2A$Intercept[2] <- "9.81 ****"
stable.FCGR2A$Slope[3] <- "0.13 **"

stable.FCGR2A.grob <- ggtexttable(stable.FCGR2A,
                                 rows = NULL, theme =  ttheme(base_size = 16))

FCGR2A <- ggplot(genes.expr.melted.neutUP %>% filter(gene_name == "FCGR2A"), aes(x = visitTime, y = logvalue)) +
  geom_boxplot(alpha = 0.8, aes(group = interaction(visitTime, group), fill = group)) +
  ggtitle("FCGR2A") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 2, 3), labels = c("0" = 0, "0.5" = 0.5, "1" = 1, "2" = 2, "3" = 3)) +
  scale_fill_manual(values = myPalette) +
  scale_color_manual(values = myPalette) +
  theme_bw(base_size = 18) +
  labs(x = "Visit (years)", y = "log2(CPM + 3)") +
  annotation_custom(ggplotGrob(stable.FCGR2A.grob),
                    xmin = 1.1, ymin = 14) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 18),
        plot.margin = margin(t = 65, r = 5.5, b = 5.5, l = 5.5, unit = "pt"))
FCGR2A



stable.LRRK2 <- linModSum.gene.patno %>%
  filter(gene_name == "LRRK2") %>%
  group_by(group) %>%
  summarise(Intercept = round(mean(Intercept, na.rm = TRUE), digits = 2),
            Slope = round(mean(Slope, na.rm = TRUE), digits = 2))

compare_means(Intercept ~ group, linModSum.gene.patno %>% filter(gene_name == "LRRK2"))
compare_means(Slope ~ group, linModSum.gene.patno %>% filter(gene_name == "LRRK2"))

stable.LRRK2$Intercept[2] <- "10.8 ***"
stable.LRRK2$Slope[3] <- "0.09 *"

stable.LRRK2.grob <- ggtexttable(stable.LRRK2,
                                 rows = NULL, theme =  ttheme(base_size = 16))

LRRK2 <- ggplot(genes.expr.melted.neutUP %>% filter(gene_name == "LRRK2"), aes(x = visitTime, y = logvalue)) +
  geom_boxplot(alpha = 0.8, aes(group = interaction(visitTime, group), fill = group)) +
  ggtitle("LRRK2") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 2, 3), labels = c("0" = 0, "0.5" = 0.5, "1" = 1, "2" = 2, "3" = 3)) +
  scale_fill_manual(values = myPalette) +
  scale_color_manual(values = myPalette) +
  theme_bw(base_size = 18) +
  labs(x = "Visit (years)", y = "log2(CPM + 3)") +
  annotation_custom(ggplotGrob(stable.LRRK2.grob),
                    xmin = 1.1, ymin = 15.5) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 18),
        plot.margin = margin(t = 65, r = 5.5, b = 5.5, l = 5.5, unit = "pt"))
LRRK2


# putting it all together --------------------------------------------
ggarrange(ggarrange(glm.cbc.plot, slope.neut.counts.byPATNO,
                    common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1,
                    font.label = list(size = 24),
                    widths = c(1.5, 1)),
          ggarrange(glm.expr.plot, slope.neut.expr.byPATNO,
                    common.legend = TRUE, legend = "bottom", ncol = 2, nrow = 1,
                    font.label = list(size = 24),
                    widths = c(1.5, 1)),
          LRRK2, FCGR2A,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          font.label = list(size = 24),
          legend = "bottom",
          labels = c("auto"))
ggsave("Figures/Figure6.png", width = 12, height = 12)




