# DESCRIPTION -------------------------------------------------------------------------------------------------
# Author: Elizabeth Hutchins
# Date: May 19th, 2020

#Purpose: overall data summary

#Imports:
# gene lists from limma

#Exports:
# 

# load packages -------------------------------------------------------------------------------------------------
library(tidyverse) #basis for data manipulation
library(ggpubr) #publication ready plots


#load metadata  ----------------------------------------------------------------------------------------------
meta <- read_csv("data/megaMetaTable.csv")
meta <- subset(meta, QCflagIR3 != "remove")
names(meta)


#PCA plot  ----------------------------------------------------------------------------------------------
#create new QC flag column for editing
meta$QCflagIR3Plot <- as.character(meta$QCflagIR3)

#set QCflagIR3Plot to "pool" for all pool samples
meta$QCflagIR3Plot[meta$DIAGNOSIS=="HCPOOL"] <- "HC pool"
meta$QCflagIR3Plot[meta$DIAGNOSIS=="PDPOOL"] <- "PD pool"
meta$QCflagIR3Plot[meta$SEX =="Male"] <- "Male"
meta$QCflagIR3Plot[meta$SEX =="Female"] <- "Female"
meta$QCflagIR3Plot[meta$QCflagIR3 =="fail"] <- "fail"
meta$QCflagIR3Plot <- factor(meta$QCflagIR3Plot, levels = c("fail", "Female", "Male", "HC pool", "PD pool"))

df <- arrange(meta, QCflagIR3Plot)
#plot QC coded PCA
PCA.byPass <- ggplot(df, aes(x=PC1, y=PC2, color=QCflagIR3Plot)) +
  geom_point(alpha = 0.6) +
  theme_bw(base_size = 18) +
  xlab("PC1: 25.4% variance") +
  ylab("PC2: 16.7% variance") +
  scale_color_manual(values = c("gray70", "#222222", "#2171b5", "#ef3b2c", "#EE9402")) +
  theme(legend.title = element_blank())
PCA.byPass

# input read histogram --------------------------------------------------------------------------------------
meta$millionPairs <- meta$total_reads/1000000

readHist.byQC <- gghistogram(meta, x = "millionPairs",
                        add = "mean", rug = TRUE,
                        color = "QCflagIR3Plot", fill = "QCflagIR3Plot",
                        palette = c("gray70", "#222222", "#2171b5", "#ef3b2c", "#EE9402"),
                        xlab = "input reads (million pairs)",
                        ylab = "frequency"
) +
  theme_bw(base_size = 18) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
readHist.byQC


meta %>%
  filter(QCflagIR3Plot != "fail") %>%
  summarize(mean = mean(total_reads),
            mediam = median(total_reads),
            IQR = IQR(total_reads))


#gene and transcript diversity ----------------------------------------------------------------------------
geneDetectedBioTable <- read_csv("data/genesDetected_byBiotype_passFilter.csv")
transcriptDetectedBioTable <- read_csv("data/transcriptsDetected_byBiotype_passFilter.csv")

genesDetected <- left_join(geneDetectedBioTable, meta[,c("HudAlphaSampleName", "QCflagIR3Plot")],
                           by = c("Sample" = "HudAlphaSampleName"))

genesDetected.longer <- genesDetected %>% 
  pivot_longer(-c(Sample, QCflagIR3Plot), names_to = "variable", values_to = "value") %>%
  filter(grepl(" > 0", variable)) %>% #only include > 0, not > 50
  mutate(biotype = sapply(str_split(variable, " "), '[', 1),
         cutoff = gsub("protein_coding|lncRNA|ncRNA|pseudogene|other", "count", variable),  #add biotype column
         type = "gene",
         group = gsub("(.*)ale", "Sample", QCflagIR3Plot))

transcriptsDetected <- left_join(transcriptDetectedBioTable, meta[,c("HudAlphaSampleName", "QCflagIR3Plot")],
                           by = c("Sample" = "HudAlphaSampleName"))

transcriptsDetected.longer <- transcriptsDetected %>% 
  pivot_longer(-c(Sample, QCflagIR3Plot), names_to = "variable", values_to = "value") %>%
  filter(grepl(" > 0", variable)) %>% #only include > 0, not > 50
  mutate(biotype = sapply(str_split(variable, " "), '[', 1),
         cutoff = gsub("protein_coding|lncRNA|ncRNA|pseudogene|other", "count", variable),  #add biotype column
         type = "transcript",
         group = gsub("(.*)ale", "Sample", QCflagIR3Plot))

diversity <- rbind(genesDetected.longer, transcriptsDetected.longer)

diversity$biotype <- factor(diversity$biotype,
                            levels = c("protein_coding",
                                       "lncRNA",
                                       "ncRNA",
                                       "pseudogene",
                                       "other"))
diversity$group <- factor(diversity$group,
                          levels = c("Sample", "HC pool", "PD pool"))

means <- genesDetected.longer %>%
  group_by(group, variable) %>%
  summarize(mean = mean(value),
            label = format(round(mean(value, na.rm=T), 0), big.mark=",")) 

fun_median <- function(x){
  return(data.frame(y = median(x), label = format(round(median(x, na.rm=T), 0), big.mark=",")
                    )
         )
  }

diversityPlot <- ggplot(diversity, aes(x = group, y = value, fill = group)
) +
  geom_boxplot(width = 0.75, alpha = 0.8, outlier.color = "grey30") +
  stat_summary(fun.data = fun_mean, geom = "text", vjust = -1) +
  stat_summary(fun = mean, geom = "point", colour = "black", size = 2) +
  # geom_text(aes(label = gene_no), position = position_stack(vjust = 0.5),
  #           colour = "white", size = 6) +
  #geom_text(data = means, aes(label = label, x = group, y = mean)) +
  xlab("") +
  facet_grid(type ~ biotype, scales = "free") +
  scale_fill_manual(values = c("gray40", "#ef3b2c", "#EE9402")) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ylab("number detected")
diversityPlot



genesDetectedTab <- genesDetected.longer %>%
  group_by(group, variable) %>%
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = variable, values_from = mean)

transcriptsDetectedTab <- transcriptsDetected.longer %>%
  group_by(group, variable) %>%
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = variable, values_from = mean)

write_csv(genesDetectedTab, "Figures/meanGenesDetected.tsv")
write_csv(transcriptsDetectedTab, "Figures/meanTranscriptsDetected.tsv")


genesDetected.longer %>%
  group_by(variable) %>%
  filter(QCflagIR3Plot == "Female" | QCflagIR3Plot == "Male") %>%
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = variable, values_from = mean)

transcriptsDetected.longer %>%
  group_by(variable) %>%
  filter(QCflagIR3Plot == "Female" | QCflagIR3Plot == "Male") %>%
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = variable, values_from = mean)


genesDetected.longer %>%
  group_by(variable) %>%
  filter(QCflagIR3Plot == "HC pool" | QCflagIR3Plot == "PD pool") %>%
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = variable, values_from = mean)

transcriptsDetected.longer %>%
  group_by(variable) %>%
  filter(QCflagIR3Plot == "HC pool" | QCflagIR3Plot == "PD pool") %>%
  summarize(mean = mean(value)) %>%
  pivot_wider(names_from = variable, values_from = mean)

# put it all together ---------------------

ggarrange(ggarrange(PCA.byPass, readHist.byQC, ncol = 2, nrow = 1,
                    common.legend = TRUE,
                    legend = "bottom",
                    labels = c("b", "c"),
                    font.label = list(size = 18)
                    ),
          diversityPlot,
          ncol = 1, nrow = 2,
          labels = c("", "d"),
          font.label = list(size = 18),
          heights = c(1,1))
ggsave("Figures/Figure2_BCD.png", height = 12, width = 12)

