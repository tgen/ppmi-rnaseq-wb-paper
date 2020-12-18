# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: June 15th, 2020

#Purpose: pathway analysis after differential expression

#Imports:
# gene lists from limma

#Exports:
# 


# load packages -------------------------------------------------------------------------------------------------
library(tidyverse)
library(fgsea) #functional gene set enrichment
library(GOfuncR)
library(viridis)
library(GO.db)
library(org.Hs.eg.db) #load homo sapiens database
library(ggpubr) #publication ready plots
library(devtools) #load R scripts from github

source_gist("c579c4ddc06fd2ceb121e690dd9df186") #color palettes

`%notin%` <- Negate(`%in%`)

myPalette <- getPalette(50, "sequential")
dev.off()
names(myPalette) <- NULL

#import MSigDB C5 file
gmtfile <- "data/c5.all.v7.1.symbols.gmt"
gmtfile.bp <- "data/c5.bp.v7.1.symbols.gmt"

c5.all <- gmtPathways(gmtfile)
c5.bp <- gmtPathways(gmtfile.bp)


# import gene names -------------------------------------------------------------------------------------------------
tmpG <- read_tsv(paste0(toolkitPath, "sourceFiles/GRCh38_GENCODE29_geneInfo.txt"),
                 skip = 1)
genes.anno <- tmpG[,c("gene_id", "gene_name", "gene_type")]
genes.anno <- unique(genes.anno)

#load gene lists ------------------------------------------------------------------------------------------------
#grab file names

PPMI.case.control.DE <- read_tsv("analysis/DE/allPD_vs_allCon_sizeMatched/PD_cutoffs.tsv")

compName <- "PPMI_Case_Control"

background <- read_tsv("analysis/DE/allPD_vs_allCon_sizeMatched/PD_full.tsv") #background list
backgroundGenes <- background$gene_name
length(backgroundGenes)


# GOfunc -----
candi_gene_ids = PPMI.case.control.DE$gene_name
bg_gene_ids = backgroundGenes
is_candidate = c(rep(1,length(candi_gene_ids)), rep(0,length(bg_gene_ids)))
input_hyper_bg = data.frame(gene_ids = c(candi_gene_ids, bg_gene_ids),
                            is_candidate)
head(input_hyper_bg)
res_hyper_bg = go_enrich(input_hyper_bg, n_randsets=100)
head(res_hyper_bg[[1]])

by(res_hyper_bg[[1]], res_hyper_bg[[1]][,'ontology'], head)

res_hyper_bp <- res_hyper_bg[[1]] %>%
  filter(ontology == "biological_process")

top_gos_hyper = res_hyper_bp[1:5, 'node_id']

plot_anno_scores(res_hyper_bg, res_hyper_bp$node_id)

plot_stats <- plot_anno_scores(res_hyper_bg, res_hyper_bp$node_id) %>%
  mutate(geneRatio = candi_genes/(candi_genes + bg_genes))

refined = refine(res_hyper_bg, fwer=0.1)
by(refined[[1]], refined[[1]][,'ontology'], head)

refined.bp <- as_tibble(refined) %>%
  filter(signif == TRUE) %>%
  filter(ontology == "biological_process") %>%
  left_join(plot_stats, by = c("node_id" = "go_id")) %>%
  mutate(full_GO = paste(node_id, node_name, sep = " "))



x <- 1
refined.bp$total_term_genes <- NA

while(x <= length(refined.bp$node_id)){
  go_id <- refined.bp$node_id[x]
  allegs <- get(go_id, org.Hs.egGO2ALLEGS)
  genes <- unlist(mget(allegs,org.Hs.egSYMBOL))
  ngenes <- length(unique(genes))
  refined.bp$total_term_genes[x] <- ngenes
  
  x <- x + 1
}


GOfunc.out <- refined.bp %>%
  mutate(inData_term_genes = candi_genes + bg_genes) %>%
  mutate(inGeneSummary = paste(candi_genes, inData_term_genes, sep = " / ")) %>%
  dplyr::select(node_id, ontology, node_name,
                total_term_genes, inData_term_genes, inGeneSummary,
                candi_genes, bg_genes, 
                root_candi_genes, root_bg_genes,
                odds_ratio, geneRatio, refined_p_overrep, full_GO)
write_tsv(GOfunc.out, "analysis/GOFunc_PPMI_PD_vs_Controls.tsv")



ggplot(GOfunc.out[order(GOfunc.out$refined_p_overrep, decreasing = FALSE),][c(1:20),], aes(reorder(full_GO, -refined_p_overrep), -log10(refined_p_overrep), label = inGeneSummary)) +
  geom_col(aes(fill = -log10(refined_p_overrep)), color = "black") +
  coord_flip() +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) + 
  scale_fill_gradient(low = "white", high = "#440154FF") +
  labs(x="GO Biological Process", y="-log10(adjusted p value)",
       title=" ") + 
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(fill = "-log10(adjusted p value)")
ggsave("Figures/Figure3b.png", height = 6, width = 10)
