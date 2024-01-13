### =====================================================================================
###   LOADING AND FORMATTING DATA
### =====================================================================================

# setwd
setwd("/home/pacificasommers/R/antarctica/S03/big_experiment/")

#load packages and functions
library(tidyverse)
library(vegan)
library(ggpubr)
library(compositions)
library(zCompositions)
library(phyloseq)
library(ANCOMBC)
library(reshape2)
`%notin%` <- function(x,y) !(x %in% y)

# #load data
# tax_table <- read.table('esvs_16s_171819.txt', header = TRUE)
# #set ESVID as rowname
# rownames(tax_table) <- tax_table$X.ESV_ID
# tax_table <- tax_table[,-1]
# #format taxonomy into multiple columns
# tax_table <-separate(tax_table,taxonomy,c("Kingdom","Phylum","Class","Order","Family","Genus",NA),sep = ";")
# tax_table <- tax_table %>% filter(Kingdom!="Eukaryota") %>% filter(Order!="Chloroplast") %>% filter(Family!="Mitochondria")
# 
# # And also remove contaminants. These were chosen by scrutinizing ESVtable.txt as a spreadsheet. Gotta re-do for fully combined.
# contaminants16s <- read.table("esvs_16s_171819_contaminants.txt")
# contaminants16s$V1 <- as.character(contaminants16s$V1)
# taxa_input_bac_nocontam <- tax_table[rownames(tax_table) %notin% contaminants16s$V1,]
# 
# pfx.map <- read.csv("PFX_metadata_withnats.csv")
# # subset only samples and ESVs of interest for this analysis
# colnames(taxa_input_bac_nocontam)[633:687] <- gsub("B","P",colnames(taxa_input_bac_nocontam)[633:687])
# pfx.bac.tax1 <- taxa_input_bac_nocontam[,832:837]
# pfx.bac <- taxa_input_bac_nocontam[,colnames(taxa_input_bac_nocontam) %in% pfx.map$SampleID]
# pfx.bac <- pfx.bac[rowSums(pfx.bac)>0,]
# pfx.bac.tax <- pfx.bac.tax1[rownames(pfx.bac.tax1) %in% rownames(pfx.bac),]
# pfx.bac.wt <- cbind(pfx.bac,pfx.bac.tax)
# 
# pfx.map <- pfx.map[order(pfx.map$SampleID),]
# 
# # First subset groups for only field experiment
# fieldonly <- subset(pfx.map,Field_or_lab=="Field")
# fieldonly <- fieldonly[fieldonly$SampleID != "P115",]
# natsonly <- subset(pfx.map,Field_or_lab=="Natural")
# # Add in previously sequenced mats
# lab <- subset(pfx.map,Field_or_lab=="Lab")
# lab0 <- subset(lab,Block==0)
# field <- rbind(fieldonly,natsonly,lab0)
# field$Order <- factor(field$Order, levels = c("BB","BO","AL","OB","OO","NO","Natural"))
# field$Timing <- factor(field$Timing, levels = c("0","1","2","Natural"))
# field <- field[order(field$SampleID),]
# 
# field.bacs <- pfx.bac[,colnames(pfx.bac) %in% field$SampleID]
# field.bacs <- field.bacs[,order(colnames(field.bacs))]
# field.bacs <- field.bacs[rowSums(field.bacs)>0,]
# 
# field.bacs.tax <- pfx.bac.tax[rownames(pfx.bac.tax) %in% rownames(field.bacs),]
# field.bacs.tax <- field.bacs.tax[order(rownames(field.bacs.tax)),]
# field.bacs <- field.bacs[order(rownames(field.bacs)),]
# field.bacs.wt <- cbind(field.bacs,field.bacs.tax)
# 
# for (i in 1:6){ field.bacs.tax[,i] <- gsub("NA","",field.bacs.tax[,i])}
# for (i in 1:nrow(field.bacs.tax)){
#   if (field.bacs.tax[i,2] == ""){
#     kingdom <- paste("Kingdom ", field.bacs.tax[i,1], sep = "")
#     field.bacs.tax[i, 2:6] <- kingdom
#   } else if (field.bacs.tax[i,3] == ""){
#     phylum <- paste("Phylum ", field.bacs.tax[i,2], sep = "")
#     field.bacs.tax[i, 3:6] <- phylum
#   } else if (field.bacs.tax[i,4] == ""){
#     class <- paste("Class ", field.bacs.tax[i,3], sep = "")
#     field.bacs.tax[i, 4:6] <- class
#   } else if (field.bacs.tax[i,5] == ""){
#     order <- paste("Order ", field.bacs.tax[i,4], sep = "")
#     field.bacs.tax[i, 5:6] <- order
#   } else if (field.bacs.tax[i,6] == ""){
#     family <- paste("Family ", field.bacs.tax[i,5], sep = "")
#     field.bacs.tax[i, 6] <- family
#   }
# }

# # Export CSVs of the ESV table and mapping file used for figures and analyses in this paper
# write.csv(field.bacs,"ESV_table_PFX_Science_2023.csv")
# write.csv(field.bacs.tax,"ESV_taxonomy_PFX_Science_2023.csv")
# write.csv(field,"Metadata_mapping_table_PFX_Science_2023.csv")

field.bacs <- read_csv("ESV_table_PFX_Science_2023.csv",header = TRUE)
field.bacs.tax <- read_csv("ESV_taxonomy_PFX_Science_2023.csv",header = TRUE)
field <- read_csv("Metadata_mapping_table_PFX_Science_2023.csv",header = TRUE)


### =====================================================================================
###  RANK ABUNDANCE PLOTS
### =====================================================================================
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
rank_cols = gg_color_hue(6)

# Make field$SampleID into rownames(field) so it can later be phyloseqized
rownames(field) <- field$SampleID
field$Order <- factor(field$Order)
field$Group <- as.factor(paste0(field$Order,field$Timing))
OTU <- otu_table(field.bacs,taxa_are_rows = TRUE)
field.bacs.tax.w7 <- field.bacs.tax
field.bacs.tax.w7$taxonomy7 <- rownames(field.bacs.tax.w7)
field.bacs.tax.fp <- field.bacs.tax.w7
field.bacs.tax.fp$Taxon <- paste0(field.bacs.tax.fp$Genus," (SV ",sub(".*_","",field.bacs.tax.fp$taxonomy7),")")
field.bacs.tax.fp$Taxon <- gsub("Tychonema_CCAP_1459-11B","Tychonema",field.bacs.tax.fp$Taxon)
field.bacs.tax.fp$Taxon <- gsub("Nostoc_PCC-73102","Nostoc",field.bacs.tax.fp$Taxon)
field.bacs.tax.fp$Taxon <- gsub("Leptolyngbya_ANT.L67.1","Leptolyngbya",field.bacs.tax.fp$Taxon)
field.bacs.tax.fp$Taxon <- gsub("Wilmottia_Ant-Ph58","Wilmottia",field.bacs.tax.fp$Taxon)
field.bacs.tax.fp$Taxon <- gsub("Aphanizomenon_NIES81","Aphanizomenon",field.bacs.tax.fp$Taxon)
field.bacs.tax.fp$Taxon <- gsub("Pseudanabaena_PCC-7429","Pseudanabaena",field.bacs.tax.fp$Taxon)
field.bacs.tax.fp$Taxon <- gsub("Clostridium_sensu_stricto_13","Clostridium",field.bacs.tax.fp$Taxon)


oo0 <- subset(field,Order=="OO" & Timing=="0")
bb0 <- subset(field,Order=="BB" & Timing=="0")
ob1 <- subset(field,Order=="OB" & Timing=="1")
ob2 <- subset(field,Order=="OB" & Timing=="2")
bo1 <- subset(field,Order=="BO" & Timing=="1")
bo2 <- subset(field,Order=="BO" & Timing=="2")

oo0.bacs <- field.bacs[,colnames(field.bacs) %in% oo0$SampleID]
bb0.bacs <- field.bacs[,colnames(field.bacs) %in% bb0$SampleID]
ob1.bacs <- field.bacs[,colnames(field.bacs) %in% ob1$SampleID]
bo1.bacs <- field.bacs[,colnames(field.bacs) %in% bo1$SampleID]
ob2.bacs <- field.bacs[,colnames(field.bacs) %in% ob2$SampleID]
bo2.bacs <- field.bacs[,colnames(field.bacs) %in% bo2$SampleID]

oo0.bacs.per <- 100*sweep(oo0.bacs,MARGIN=2,FUN="/",STATS=colSums(oo0.bacs))
oo0.bacs.per <- oo0.bacs.per[order(rowSums(oo0.bacs.per),decreasing = TRUE),]
oo0.bacs.per <- oo0.bacs.per[1:10,]
oo0.tax.10 <- field.bacs.tax.fp[field.bacs.tax.fp$taxonomy7 %in% rownames(oo0.bacs.per),]
oo0.tax.10 <- oo0.tax.10[match(rownames(oo0.bacs.per),rownames(oo0.tax.10)),]
oo0.bacs.per$SequenceVariant <- oo0.tax.10$Taxon
oo0.bacs.per$Phylum <- oo0.tax.10$Phylum
oo0.bacs.fp <- melt(oo0.bacs.per)
oo0.bacs.fp$SequenceVariant <- reorder(oo0.bacs.fp$SequenceVariant, -oo0.bacs.fp$value)
oo0.bacs.fp$SequenceVariant <- reorder(oo0.bacs.fp$SequenceVariant, -oo0.bacs.fp$value)
oo0.bacs.fp$Phylum <- factor(oo0.bacs.fp$Phylum,levels = c("Actinobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Proteobacteria","Verrucomicrobia"))
oo0.bacs.plot <- ggplot(oo0.bacs.fp, aes(x = SequenceVariant, y = value))  +
  scale_y_continuous(expand = c(0, 0),limits = c(0,50)) +
  geom_boxplot(aes(color=Phylum),outlier.colour = "white") +
  geom_point(aes(color=Phylum)) +
  scale_color_manual(limits=levels(oo0.bacs.fp$Phylum),values=rank_cols) +
  labs(x="",y="Percent abundance",
       title="Orange mat before experiment") +
  theme(
    plot.title = element_text(color="black",size=12) 
    ,plot.background = element_blank()
    ,panel.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,axis.line = element_line(colour="black", linewidth=1)
    ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
    ,axis.text.y = element_text(color="black", size=12)
    ,axis.title.x = element_text(colour="black", size=12)
    ,axis.title.y = element_text(color="black", size=12)
    #,legend.position = "none"
  )
oo0.bacs.plot

bb0.bacs.per <- 100*sweep(bb0.bacs,MARGIN=2,FUN="/",STATS=colSums(bb0.bacs))
bb0.bacs.per <- bb0.bacs.per[order(rowSums(bb0.bacs.per),decreasing = TRUE),]
bb0.bacs.per <- bb0.bacs.per[1:10,]
bb0.tax.10 <- field.bacs.tax.fp[field.bacs.tax.fp$taxonomy7 %in% rownames(bb0.bacs.per),]
bb0.tax.10 <- bb0.tax.10[match(rownames(bb0.bacs.per),rownames(bb0.tax.10)),]
bb0.bacs.per$SequenceVariant <- bb0.tax.10$Taxon
bb0.bacs.per$Phylum <- bb0.tax.10$Phylum
bb0.bacs.fp <- melt(bb0.bacs.per)
bb0.bacs.fp$SequenceVariant <- reorder(bb0.bacs.fp$SequenceVariant, -bb0.bacs.fp$value)
bb0.bacs.fp$SequenceVariant <- reorder(bb0.bacs.fp$SequenceVariant, -bb0.bacs.fp$value)
bb0.bacs.fp$Phylum <- factor(bb0.bacs.fp$Phylum,levels = c("Actinobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Proteobacteria","Verrucomicrobia"))
bb0.bacs.plot <- ggplot(bb0.bacs.fp, aes(x = SequenceVariant, y = value))  +
  scale_y_continuous(expand = c(0, 0),limits = c(0,50)) +
  geom_boxplot(aes(color=Phylum),outlier.colour = "white") +
  geom_point(aes(color=Phylum)) +
  scale_color_manual(limits=levels(bb0.bacs.fp$Phylum),values=rank_cols) +
  labs(x="",y="",
       title="Black mat before experiment") +
  theme(
    plot.title = element_text(color="black",size=12) 
    ,plot.background = element_blank()
    ,panel.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,axis.line = element_line(colour="black", linewidth=1)
    ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
    ,axis.text.y = element_text(color="black", size=12)
    ,axis.title.x = element_text(colour="black", size=12)
    ,axis.title.y = element_text(color="black", size=12)
    #,legend.position = "none"
  )
bb0.bacs.plot


ob1.bacs.per <- 100*sweep(ob1.bacs,MARGIN=2,FUN="/",STATS=colSums(ob1.bacs))
ob1.bacs.per <- ob1.bacs.per[order(rowSums(ob1.bacs.per),decreasing = TRUE),]
ob1.bacs.per <- ob1.bacs.per[1:10,]
ob1.tax.10 <- field.bacs.tax.fp[field.bacs.tax.fp$taxonomy7 %in% rownames(ob1.bacs.per),]
ob1.tax.10 <- ob1.tax.10[match(rownames(ob1.bacs.per),rownames(ob1.tax.10)),]
ob1.bacs.per$SequenceVariant <- ob1.tax.10$Taxon
ob1.bacs.per$Phylum <- ob1.tax.10$Phylum
ob1.bacs.fp <- melt(ob1.bacs.per)
ob1.bacs.fp$SequenceVariant <- reorder(ob1.bacs.fp$SequenceVariant, -ob1.bacs.fp$value)
ob1.bacs.fp$SequenceVariant <- reorder(ob1.bacs.fp$SequenceVariant, -ob1.bacs.fp$value)
ob1.bacs.fp$Phylum <- factor(ob1.bacs.fp$Phylum,levels = c("Actinobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Proteobacteria","Verrucomicrobia"))
ob1.bacs.plot <- ggplot(ob1.bacs.fp, aes(x = SequenceVariant, y = value))  +
  scale_y_continuous(expand = c(0, 0),limits = c(0,50)) +
  geom_boxplot(aes(color=Phylum),outlier.colour = "white") +
  geom_point(aes(color=Phylum)) +
  scale_color_manual(limits=levels(ob1.bacs.fp$Phylum),values=rank_cols) +
  labs(x="",y="Percent abundance",
       title="Orange then black mat after 1st season") +
  theme(
    plot.title = element_text(color="black",size=12) 
    ,plot.background = element_blank()
    ,panel.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,axis.line = element_line(colour="black", linewidth=1)
    ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
    ,axis.text.y = element_text(color="black", size=12)
    ,axis.title.x = element_text(colour="black", size=12)
    ,axis.title.y = element_text(color="black", size=12)
    #,legend.position = "none"
  )
ob1.bacs.plot

bo1.bacs.per <- 100*sweep(bo1.bacs,MARGIN=2,FUN="/",STATS=colSums(bo1.bacs))
bo1.bacs.per <- bo1.bacs.per[order(rowSums(bo1.bacs.per),decreasing = TRUE),]
bo1.bacs.per <- bo1.bacs.per[1:10,]
bo1.tax.10 <- field.bacs.tax.fp[field.bacs.tax.fp$taxonomy7 %in% rownames(bo1.bacs.per),]
bo1.tax.10 <- bo1.tax.10[match(rownames(bo1.bacs.per),rownames(bo1.tax.10)),]
bo1.bacs.per$SequenceVariant <- bo1.tax.10$Taxon
bo1.bacs.per$Phylum <- bo1.tax.10$Phylum
bo1.bacs.fp <- melt(bo1.bacs.per)
bo1.bacs.fp$SequenceVariant <- reorder(bo1.bacs.fp$SequenceVariant, -bo1.bacs.fp$value)
bo1.bacs.fp$SequenceVariant <- reorder(bo1.bacs.fp$SequenceVariant, -bo1.bacs.fp$value)
bo1.bacs.fp$Phylum <- factor(bo1.bacs.fp$Phylum,levels = c("Actinobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Proteobacteria","Verrucomicrobia"))
bo1.bacs.plot <- ggplot(bo1.bacs.fp, aes(x = SequenceVariant, y = value))  +
  scale_y_continuous(expand = c(0, 0),limits = c(0,50)) +
  geom_boxplot(aes(color=Phylum),outlier.colour = "white") +
  geom_point(aes(color=Phylum)) +
  scale_color_manual(limits=levels(bo1.bacs.fp$Phylum),values=rank_cols) +
  labs(x="",y="",
       title="Black then orange mat after 1st season") +
  theme(
    plot.title = element_text(color="black",size=12) 
    ,plot.background = element_blank()
    ,panel.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,axis.line = element_line(colour="black", linewidth=1)
    ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
    ,axis.text.y = element_text(color="black", size=12)
    ,axis.title.x = element_text(colour="black", size=12)
    ,axis.title.y = element_text(color="black", size=12)
    #,legend.position = "none"
  )
bo1.bacs.plot


ob2.bacs.per <- 100*sweep(ob2.bacs,MARGIN=2,FUN="/",STATS=colSums(ob2.bacs))
ob2.bacs.per <- ob2.bacs.per[order(rowSums(ob2.bacs.per),decreasing = TRUE),]
ob2.bacs.per <- ob2.bacs.per[1:10,]
ob2.tax.10 <- field.bacs.tax.fp[field.bacs.tax.fp$taxonomy7 %in% rownames(ob2.bacs.per),]
ob2.tax.10 <- ob2.tax.10[match(rownames(ob2.bacs.per),rownames(ob2.tax.10)),]
ob2.bacs.per$SequenceVariant <- ob2.tax.10$Taxon
ob2.bacs.per$Phylum <- ob2.tax.10$Phylum
ob2.bacs.fp <- melt(ob2.bacs.per)
ob2.bacs.fp$SequenceVariant <- reorder(ob2.bacs.fp$SequenceVariant, -ob2.bacs.fp$value)
ob2.bacs.fp$SequenceVariant <- reorder(ob2.bacs.fp$SequenceVariant, -ob2.bacs.fp$value)
ob2.bacs.fp$Phylum <- factor(ob2.bacs.fp$Phylum,levels = c("Actinobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Proteobacteria","Verrucomicrobia"))
ob2.bacs.plot <- ggplot(ob2.bacs.fp, aes(x = SequenceVariant, y = value))  +
  scale_y_continuous(expand = c(0, 0),limits = c(0,50)) +
  geom_boxplot(aes(color=Phylum),outlier.colour = "white") +
  geom_point(aes(color=Phylum)) +
  scale_color_manual(limits=levels(ob2.bacs.fp$Phylum),values=rank_cols) +
  labs(x="",y="Percent abundance",
       title="Orange then black mat after 2nd season") +
  theme(
    plot.title = element_text(color="black",size=12) 
    ,plot.background = element_blank()
    ,panel.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,axis.line = element_line(colour="black", linewidth=1)
    ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
    ,axis.text.y = element_text(color="black", size=12)
    ,axis.title.x = element_text(colour="black", size=12)
    ,axis.title.y = element_text(color="black", size=12)
    #,legend.position = "none"
  )
ob2.bacs.plot

bo2.bacs.per <- 100*sweep(bo2.bacs,MARGIN=2,FUN="/",STATS=colSums(bo2.bacs))
bo2.bacs.per <- bo2.bacs.per[order(rowSums(bo2.bacs.per),decreasing = TRUE),]
bo2.bacs.per <- bo2.bacs.per[1:10,]
bo2.tax.10 <- field.bacs.tax.fp[field.bacs.tax.fp$taxonomy7 %in% rownames(bo2.bacs.per),]
bo2.tax.10 <- bo2.tax.10[match(rownames(bo2.bacs.per),rownames(bo2.tax.10)),]
bo2.bacs.per$SequenceVariant <- bo2.tax.10$Taxon
bo2.bacs.per$Phylum <- bo2.tax.10$Phylum
bo2.bacs.fp <- melt(bo2.bacs.per)
bo2.bacs.fp$SequenceVariant <- reorder(bo2.bacs.fp$SequenceVariant, -bo2.bacs.fp$value)
bo2.bacs.fp$SequenceVariant <- reorder(bo2.bacs.fp$SequenceVariant, -bo2.bacs.fp$value)
bo2.bacs.fp$Phylum <- factor(bo2.bacs.fp$Phylum,levels = c("Actinobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Proteobacteria","Verrucomicrobia"))
bo2.bacs.plot <- ggplot(bo2.bacs.fp, aes(x = SequenceVariant, y = value))  +
  scale_y_continuous(expand = c(0, 0),limits = c(0,50)) +
  geom_boxplot(aes(color=Phylum),outlier.colour = "white") +
  geom_point(aes(color=Phylum)) +
  scale_color_manual(limits=levels(bo2.bacs.fp$Phylum),values=rank_cols) +
  labs(x="",y="",
       title="Black then orange mat after 2nd season") +
  theme(
    plot.title = element_text(color="black",size=12) 
    ,plot.background = element_blank()
    ,panel.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,axis.line = element_line(colour="black", linewidth=1)
    ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
    ,axis.text.y = element_text(color="black", size=12)
    ,axis.title.x = element_text(colour="black", size=12)
    ,axis.title.y = element_text(color="black", size=12)
    #,legend.position = "none"
  )
bo2.bacs.plot

rank.abund.plot <- ggarrange(oo0.bacs.plot, bb0.bacs.plot, 
                             ob1.bacs.plot,bo1.bacs.plot,
                             ob2.bacs.plot,bo2.bacs.plot,
                                  labels = c("A", "B", "C", "D", "E", "F"),
                                  ncol = 2, nrow = 3,
                                  align = "hv",
                                  common.legend = TRUE)

svg(filename="PFX_bac_rank.svg", 
    width=9, 
    height=15, 
    pointsize=12)
rank.abund.plot
dev.off()


### =====================================================================================
###  GENERA DRIVING DIFFERENCE AMONG COMMUNITIES
### =====================================================================================

# Use ANCOM BC to highlight sequences variants driving differences among groups.
# Specifically focus on cyanos for Fig 4 in main text
TAX <- tax_table(as.matrix(field.bacs.tax.w7))
# Now fix missing taxonomy info so you can make nice plots later
for (i in 1:6){ TAX[,i] <- as.character(TAX[,i])}
for (i in 1:6){ TAX[,i] <- gsub("NA","",TAX[,i])}
for (i in 1:nrow(TAX)){
  if (TAX[i,2] == ""){
    kingdom <- paste("Kingdom_", TAX[i,1], sep = "")
    TAX[i, 2:6] <- kingdom
  } else if (TAX[i,3] == ""){
    phylum <- paste("Phylum_", TAX[i,2], sep = "")
    TAX[i, 3:6] <- phylum
  } else if (TAX[i,4] == ""){
    class <- paste("Class_", TAX[i,3], sep = "")
    TAX[i, 4:6] <- class
  } else if (TAX[i,5] == ""){
    order <- paste("Order_", TAX[i,4], sep = "")
    TAX[i, 5:6] <- order
  } else if (TAX[i,6] == ""){
    family <- paste("Family_", TAX[i,5], sep = "")
    TAX[i, 6] <- family
  }
}

sampledata <- sample_data(field)
pseq <- phyloseq(OTU,TAX,sampledata)
# Now subset to be what you need
pseq.oobb <- subset_samples(pseq, Order=="OO" | Order=="BB")
pseq.oobb.0 <- subset_samples(pseq.oobb, Timing=="0")
pseq.obbo <- subset_samples(pseq, Order=="OB" | Order=="BO")
pseq.obbo.1 <- subset_samples(pseq.obbo, Timing=="1")
pseq.obbo.2 <- subset_samples(pseq.obbo, Timing=="2")


# --------------------- FIRST COMPARE OO & BB PRE-GLACIER COMMUNITIES -------------------
# Previously ran it on genus level, but let's try ESV level:
# set.seed(1802)
# oobb0_out = ancombc2(data = pseq.oobb.0,assay_name = "counts",
#                      tax_level = "Genus", fix_formula = "Order",
#                     p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE, prv_cut = 0.1, lib_cut = 0,
#                     s0_perc = 0.05,
#                     group = NULL, struc_zero = FALSE, neg_lb = FALSE,  alpha = 0.05,   n_cl = 1,
#                     verbose = FALSE,
#                     global = FALSE,
#                     pairwise = FALSE,
#                     dunnet = FALSE,
#                     trend = FALSE,
#                     iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
#                     em_control = list(tol = 1e-05, max_iter = 100),
#                     lme_control = lme4::lmerControl(),
#                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
#                     trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)
# )
# 
# 
# oobb0_df <- oobb0_out$res %>% as_tibble() %>%
#   mutate(Genus=rownames(oobb0_out$res)) #%>%
# #  filter(diff_OrderOO=="TRUE")
# oobb0_df <- oobb0_df %>% left_join(field.bacs.tax %>% as_tibble() %>%
#                                      filter(Genus %in% oobb0_df$Genus) %>%
#                                      dplyr::select(Genus,Phylum) %>%
#                                      unique())
# oobb0_df_cyanos <- oobb0_df %>% filter(Phylum=="Cyanobacteria")
# 
# 
# # Now invasion holes after one season
# set.seed(1455)
# obbo1_out = ancombc2(data = pseq.obbo.1,assay_name = "counts",
#                      tax_level = "Genus", fix_formula = "Order",
#                      p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE, prv_cut = 0.1, lib_cut = 0,
#                      s0_perc = 0.05,
#                      group = NULL, struc_zero = FALSE, neg_lb = FALSE,  alpha = 0.05,   n_cl = 1,
#                      verbose = FALSE,
#                      global = FALSE,
#                      pairwise = FALSE,
#                      dunnet = FALSE,
#                      trend = FALSE,
#                      iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
#                      em_control = list(tol = 1e-05, max_iter = 100),
#                      lme_control = lme4::lmerControl(),
#                      mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
#                      trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)
# )
# 
# obbo1_df <- obbo1_out$res %>% as_tibble() %>%
#   mutate(Genus=rownames(obbo1_out$res)) #%>%
# #  filter(diff_OrderOO=="TRUE")
# obbo1_df <- obbo1_df %>% left_join(field.bacs.tax %>% as_tibble() %>%
#                                      filter(Genus %in% obbo1_df$Genus) %>%
#                                      dplyr::select(Genus,Phylum) %>%
#                                      unique())
# obbo1_df_cyanos <- obbo1_df %>% filter(Phylum=="Cyanobacteria")
# 
# 
# 
# # And invasion holes after both seasons
# set.seed(1511)
# obbo2_out = ancombc2(data = pseq.obbo.2,assay_name = "counts",
#                      tax_level = "Genus", fix_formula = "Order",
#                      p_adj_method = "bonferroni", pseudo = 0, pseudo_sens = TRUE, prv_cut = 0.1, lib_cut = 0,
#                      s0_perc = 0.05,
#                      group = NULL, struc_zero = FALSE, neg_lb = FALSE,  alpha = 0.05,   n_cl = 1,
#                      verbose = FALSE,
#                      global = FALSE,
#                      pairwise = FALSE,
#                      dunnet = FALSE,
#                      trend = FALSE,
#                      iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
#                      em_control = list(tol = 1e-05, max_iter = 100),
#                      lme_control = lme4::lmerControl(),
#                      mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
#                      trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)
# )
# 
# obbo2_df <- obbo2_out$res %>% as_tibble() %>%
#   mutate(Genus=rownames(obbo2_out$res)) #%>%
# #  filter(diff_OrderOB=="TRUE")
# obbo2_df <- obbo2_df %>% left_join(field.bacs.tax %>% as_tibble() %>%
#               filter(Genus %in% obbo2_df$Genus) %>%
#               dplyr::select(Genus,Phylum) %>%
#               unique())
# obbo2_df_cyanos <- obbo2_df %>% filter(Phylum=="Cyanobacteria")
# 
# write.csv(oobb0_df,"ancombc_output_oobb0_bonferroni_20221222.csv")
# write.csv(obbo1_df,"ancombc_output_obbo1_bonferroni_20221222.csv")
# write.csv(obbo2_df,"ancombc_output_obbo2_bonferroni_20221222.csv")

oobb0_df <- read.csv("ancombc_output_oobb0_20221216.csv")
  oobb0_df_cyanos <- oobb0_df %>% filter(Phylum=="Cyanobacteria")
obbo1_df <- read.csv("ancombc_output_obbo1_20221216.csv")
  obbo1_df_cyanos <- obbo1_df %>% filter(Phylum=="Cyanobacteria")
obbo2_df <- read.csv("ancombc_output_obbo2_20221216.csv")
  obbo2_df_cyanos <- obbo2_df %>% filter(Phylum=="Cyanobacteria")


# Collect all the ESVs that show up in one of these three plots:
cyano.esvs <- unique(c(oobb0_df_cyanos$Genus,obbo1_df_cyanos$Genus,obbo2_df_cyanos$Genus))
oobb0_df_fp <- oobb0_df %>% filter(Genus %in% cyano.esvs) %>% dplyr::select(Genus,lfc_OrderOO,se_OrderOO,diff_OrderOO)
obbo1_df_fp <- obbo1_df %>% filter(Genus %in% cyano.esvs) %>% dplyr::select(Genus,lfc_OrderOB1=lfc_OrderOB,se_OrderOB1=se_OrderOB,diff_OrderOB1=diff_OrderOB)
obbo2_df_fp <- obbo2_df %>% filter(Genus %in% cyano.esvs) %>% dplyr::select(Genus,lfc_OrderOB2=lfc_OrderOB,se_OrderOB2=se_OrderOB,diff_OrderOB2=diff_OrderOB)
df_fp <- full_join(oobb0_df_fp,obbo1_df_fp,by=c("Genus"))
df_fp <- full_join(df_fp,obbo2_df_fp,by=c("Genus"))

col0 <- c(rep(0,nrow(df_fp)))
col1 <- c(rep(0,nrow(df_fp)))
col2 <- c(rep(0,nrow(df_fp)))
for(c in 1:nrow(df_fp)){
  if (is.na(df_fp$diff_OrderOO[c])){
    col0[c] <- "NA"
  } else if (df_fp$diff_OrderOO[c]=="FALSE") {
    col0[c] <- "gray"
  } else if (df_fp$diff_OrderOO[c]=="TRUE" & df_fp$lfc_OrderOO[c] > 0) {
      col0[c] <- "orange"
  } else if (df_fp$diff_OrderOO[c]=="TRUE" && df_fp$lfc_OrderOO[c] < 0) {
      col0[c] <- "black"
  }
    if (is.na(df_fp$diff_OrderOB1[c])){
      col1[c] <- "NA"
    } else if (df_fp$diff_OrderOB1[c]=="FALSE") {
      col1[c] <- "gray"
    } else if (df_fp$diff_OrderOB1[c]=="TRUE" & df_fp$lfc_OrderOB1[c] > 1) {
      col1[c] <- "orange"
    } else if (df_fp$diff_OrderOB1[c]=="TRUE" && df_fp$lfc_OrderOB1[c] < 1) {
      col1[c] <- "black"
    }
  if (is.na(df_fp$diff_OrderOB2[c])){
    col2[c] <- "NA"
  } else if (df_fp$diff_OrderOB2[c]=="FALSE") {
    col2[c] <- "gray"
  } else if (df_fp$diff_OrderOB2[c]=="TRUE" & df_fp$lfc_OrderOB2[c] > 2) {
    col2[c] <- "orange"
  } else if (df_fp$diff_OrderOB2[c]=="TRUE" && df_fp$lfc_OrderOB2[c] < 2) {
    col2[c] <- "black"
  }
}
df_col <- data.frame(col0,col1,col2)

df_fp <-bind_cols(df_fp,df_col)

oobb0_fig <- ggplot(df_fp,aes(x = Genus, y = lfc_OrderOO,colour = col0,fill=col0))  +
  geom_bar(stat="identity",alpha=0.5)+coord_flip()+
  geom_errorbar(aes(x=Genus, ymin=lfc_OrderOO-se_OrderOO, ymax=lfc_OrderOO+se_OrderOO,colour=col0), width=0.4, alpha=0.8, linewidth=.5) +
  scale_fill_manual(values=c("black","gray","red","orange"),aesthetics = c("colour", "fill")) +
  labs(x="Genus",y="Log-fold difference in relative abundance",
       title="Black vs. orange mat pre-experiment") +
  theme(
    plot.title = element_text(color="black",size=11,hjust=0.5)
    ,plot.background = element_blank()
    ,panel.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    #   ,axis.line = element_line(colour="black", linewidth=1)
    ,axis.text.x = element_text(color="black", size=11, angle = 70, hjust = 1)
    ,axis.text.y = element_text(color="black", size=11)
    ,axis.title.x = element_blank()
    ,plot.margin = margin(0,0,0,0, "cm")
    ,axis.title.y = element_text(color="black", size=11)
    ,legend.position = "none"
  )
oobb0_fig


obbo1_fig <- ggplot(df_fp,aes(x = Genus, y = lfc_OrderOB1,colour = col1,fill=col1))  +
  geom_bar(stat="identity",alpha=0.5)+coord_flip()+
  geom_errorbar(aes(x=Genus, ymin=lfc_OrderOB1-se_OrderOB1, ymax=lfc_OrderOB1+se_OrderOB1,colour=col1), width=0.4, alpha=0.8, linewidth=.5) +
  scale_fill_manual(values=c("gray","black","orange"),aesthetics = c("colour", "fill")) +
  labs(x="Genus",y="Log-fold difference in relative abundance",
       title="B-O vs. O-B after 1st season") +
  theme(
    plot.title = element_text(color="black",size=11,hjust=0.5)
    ,plot.background = element_blank()
    ,panel.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    #   ,axis.line = element_line(colour="black", linewidth=1)
    ,axis.text.x = element_text(color="black", size=11, angle = 70, hjust = 1)
    ,axis.text.y = element_blank()
    ,axis.title.x = element_text(colour="black", size=11)
    ,axis.title.y = element_blank()
    ,plot.margin = margin(0,0,0,0, "cm")
    ,legend.position = "none"
  )
obbo1_fig

obbo2_fig <- ggplot(df_fp,aes(x = Genus, y = lfc_OrderOB2,colour = col2,fill=col2))  +
  geom_bar(stat="identity",alpha=0.5)+coord_flip()+
  geom_errorbar(aes(x=Genus, ymin=lfc_OrderOB2-se_OrderOB2, ymax=lfc_OrderOB2+se_OrderOB2,colour=col2), width=0.4, alpha=0.8, linewidth=.5) +
 scale_fill_manual(values=c("gray","black","orange"),aesthetics = c("colour", "fill")) +
  labs(#x="Genus",y="Log-fold difference in relative abundance",
       title="B-O vs. O-B after 2nd season") +
  theme(
    plot.title = element_text(color="black",size=12,hjust=0.5)
    ,plot.background = element_blank()
    ,panel.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
#   ,axis.line = element_line(colour="black", linewidth=1)
    ,axis.text.x = element_text(color="black", size=12, angle = 70, hjust = 1)
    ,axis.text.y = element_blank()
    ,axis.title.x = element_blank()
    ,axis.title.y = element_blank()
    ,plot.margin = margin(0,0,0,0, "cm")
    ,legend.position = "none"
  )
obbo2_fig

da_plot <- ggarrange(oobb0_fig,obbo1_fig,obbo2_fig,
                     ncol = 3, nrow = 1,
                     align = "hv",
                     common.legend = FALSE)
da_plot



### calc some summary stats

# First, how many of the ESVs in the mats are cyanos, 
# what percent of the relative abundance do they make up, 
# and what proportion of differentially abundant ESVs do they make up?
oobb_map <- subset(pfx.map,Order=="OO" | Order=="BB")
oobb0_map <- subset(oobb_map,Timing=="0")
oobb0_bacs <- field.bacs[,colnames(field.bacs) %in% oobb0_map$SampleID]
oobb0_bacs <- oobb0_bacs[rowSums(oobb0_bacs)>0,]
oobb0_tax <- field.bacs.tax[rownames(field.bacs.tax) %in% rownames(oobb0_bacs),]
oobb0_cyano_tax <- subset(oobb0_tax,Phylum=="Cyanobacteria")
oobb0_cyano_abund <- oobb0_bacs[rownames(oobb0_bacs) %in% rownames(oobb0_cyano_tax),]
oobb0_cyano_percent <- round(100*colSums(oobb0_cyano_abund)/colSums(oobb0_bacs),2)
oobb0_cyano_per_df <- data.frame(colnames(oobb0_bacs),oobb0_cyano_percent)




### =====================================================================================
###  NMDS plot and pairwise permanova
### =====================================================================================

# Exclude any samples from natural cryoconite holes in the data set:
field.nonat <- subset(field,Order != "Natural")
field.bacs.nonat <- field.bacs[,colnames(field.bacs) %in% field.nonat$SampleID]

# Define groups
field.nonat.al0 <- which(field.nonat$Order=="AL" & field.nonat$Timing==0)
field.nonat.al1 <- which(field.nonat$Order=="AL" & field.nonat$Timing==1)
field.nonat.al2 <- which(field.nonat$Order=="AL" & field.nonat$Timing==2)
field.nonat.bo2 <- which(field.nonat$Order=="BO" & field.nonat$Timing==2)
field.nonat.ob2 <- which(field.nonat$Order=="OB" & field.nonat$Timing==2)
field.nonat.bo1 <- which(field.nonat$Order=="BO" & field.nonat$Timing==1)
field.nonat.ob1 <- which(field.nonat$Order=="OB" & field.nonat$Timing==1)
field.nonat.oo1 <- which(field.nonat$Order=="OO" & field.nonat$Timing==1)
field.nonat.bb1 <- which(field.nonat$Order=="BB" & field.nonat$Timing==1)
field.nonat.oo0 <- which(field.nonat$Order=="OO" & field.nonat$Timing==0)
field.nonat.bb0 <- which(field.nonat$Order=="BB" & field.nonat$Timing==0)
field.nonat.oo2 <- which(field.nonat$Order=="OO" & field.nonat$Timing==2)
field.nonat.bb2 <- which(field.nonat$Order=="BB" & field.nonat$Timing==2)
field.nonat.nat <- which(field.nonat$Order=="Natural")
field.nonat.no0 <- which(field.nonat$Order=="NO" & field.nonat$Timing==0)
field.nonat.no1 <- which(field.nonat$Order=="NO" & field.nonat$Timing==1)
field.nonat.no2 <- which(field.nonat$Order=="NO" & field.nonat$Timing==2)

# First get rid of some low-percentage taxa to reduce the zeros:
field.bac.count <- 5
field.bacs.nonat.5 <- data.frame(field.bacs.nonat[which(apply(field.bacs.nonat, 1, function(x){mean(x)}) > field.bac.count),], 
                                 check.names=F) # now 4841 -> 360 ESVs
per.com.5.nonat <- 100-round(100*(colSums(field.bacs.nonat) - colSums(field.bacs.nonat.5))/colSums(field.bacs.nonat),2)
field.nonat$per.com.5 <- per.com.5.nonat


# Plot how much of the community reads are maintained:

# Make a nice fig in ggplot
plot.per.com.5.nonat <-  ggplot(field.nonat, aes(x = Order, y = per.com.5, color = Timing)) +
  # scale_fill_manual(values=c("black","blue","green")) +
  scale_color_manual(values=c("gray","black","blue")) +
  geom_boxplot(outlier.color = 'white') +
  geom_point(position=position_dodge(w=.8)) +
  labs(x="Order",y="Percent community represented",
       title="Mean count > 5 (360 ESVs)") +
  theme(
    plot.title = element_text(face="bold",size=12) 
    ,plot.background = element_blank()
    ,panel.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
    ,axis.line = element_line(colour="black", size=1)
    ,axis.text.x = element_text(color="black", size=12)
    ,axis.text.y = element_text(color="black", size=12)
    ,axis.title.x = element_text(colour="black", size=12)
    ,axis.title.y = element_text(color="black", size=12)
    ,legend.position="none"
  )
plot.per.com.5.nonat


# to calc ILR (Aitchison dist), samples must be ROWS, so the data frame will be transposed using t()
field.bacs.nonat.czm <- cmultRepl(t(field.bacs.nonat.5),  label=0, method="CZM")
# REMEMBER: IF samples are columns, margin = 2, if samples are ROWS, margin = 1
field.bacs.nonat.ilr <- ilr(field.bacs.nonat.czm)

#Run NMDS withOUT natural holes
field.bacs.nonat.ilr.mat <- as.data.frame(ilr(field.bacs.nonat.czm))
set.seed(1523)
field.bacs.nonat.nmds <- metaMDS(field.bacs.nonat.ilr.mat, 
                                 distance = "euclidean",
                                 k = 2, try = 20, trymax = 20,
                                 engine = "monoMDS",
                                 autotransform = FALSE,
                                 noshare = FALSE,
                                 wascores = FALSE,
                                 expand = FALSE)
plot(field.bacs.nonat.nmds)
data.scores.nonat<-as.data.frame(scores(field.bacs.nonat.nmds))

dev.off()
# Use NMDS for the three-panel figure 2
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
}

par(mfrow = c(1,3))
par(mar=c(2,2,1.5,1.5),xpd=TRUE)
# First up, Pre-glacier:
plot(data.scores.nonat[,2] ~ data.scores.nonat[,1], 
     xlab="", 
     ylab="",
     #  main="A   Pre-experiment",
     cex=0.1, col="white",
     # xaxt='n',
     # yaxt='n'
)
title("A   Pre-experiment",line=0.5, adj=0)
# title(xlab=field.bacs.PC1, line=0.5)
# title(ylab=field.bacs.PC2, line=0.5)
legend("topleft",                   # Create legend outside of plot
       legend = c("Orange then orange mat",
                  "Black then black mat",
                  "Orange then black mat",
                  "Black then orange mat"),
       pch = c(1,1,19,19,19),
       col = c(rgb(red=.94, green = .52, blue = 0),
               rgb(red=0, green = 0, blue = 0),
               rgb(red=.94, green = .52, blue = 0),
               rgb(red=0, green = 0, blue = 0)),
       bty="n")
points(data.scores.nonat[field.nonat.oo0,2] ~ data.scores.nonat[field.nonat.oo0,1], pch=1,
       cex=1.5, col=rgb(red=.94, green = .52, blue = 0))
points(data.scores.nonat[field.nonat.bb0,2] ~ data.scores.nonat[field.nonat.bb0,1], pch=1,
       cex=1.5, col=rgb(red=0, green = 0, blue = 0))

plot(data.scores.nonat[,2] ~ data.scores.nonat[,1], 
     xlab="", 
     ylab="",
     # main="B   After first season",
     cex=0.1, col="white",
     # xaxt='n',
     # yaxt='n'
)
title("B   After first season",line=0.5, adj=0)
# title(xlab=field.bacs.PC1, line=0.5)
# title(ylab=field.bacs.PC2, line=0.5)


points(data.scores.nonat[field.nonat.oo1,2] ~ data.scores.nonat[field.nonat.oo1,1], pch=1,
       cex=1.5, col=rgb(red=.94, green = .52, blue = 0))
points(data.scores.nonat[field.nonat.bb1,2] ~ data.scores.nonat[field.nonat.bb1,1], pch=1,
       cex=1.5, col=rgb(red=0, green = 0, blue = 0))
# points(data.scores.nonat[field.nonat.al1,2] ~ data.scores.nonat[field.nonat.al1,1], pch=1,
#        cex=1.5, col=rgb(red=0, green = 0, blue = 1))
points(data.scores.nonat[field.nonat.ob1,2] ~ data.scores.nonat[field.nonat.ob1,1], pch=19,
       cex=1.5, col=rgb(red=.94, green = .52, blue = 0))
Plot_ConvexHull(xcoord = data.scores.nonat[field.nonat.ob1,1], ycoord = data.scores.nonat[field.nonat.ob1,2],
                lcolor = rgb(red=.94, green = .52, blue = 0))
points(data.scores.nonat[field.nonat.bo1,2] ~ data.scores.nonat[field.nonat.bo1,1], pch=19,
       cex=1.5, col=rgb(red=0, green = 0, blue = 0))
Plot_ConvexHull(xcoord = data.scores.nonat[field.nonat.bo1,1], ycoord = data.scores.nonat[field.nonat.bo1,2],
                lcolor = "black")

plot(data.scores.nonat[,2] ~ data.scores.nonat[,1], 
     xlab="", 
     ylab="",
     #  main="B   After second season",
     cex=0.1, col="white",
     # xaxt='n',
     # yaxt='n'
)
title("C   After second season",line=0.5, adj=0)
# title(xlab=field.bacs.PC1, line=0.5)
# title(ylab=field.bacs.PC2, line=0.5)

points(data.scores.nonat[field.nonat.oo2,2] ~ data.scores.nonat[field.nonat.oo2,1], pch=1,
       cex=1.5, col=rgb(red=.94, green = .52, blue = 0))
points(data.scores.nonat[field.nonat.bb2,2] ~ data.scores.nonat[field.nonat.bb2,1], pch=1,
       cex=1.5, col=rgb(red=0, green = 0, blue = 0))
# points(data.scores.nonat[field.nonat.al2,2] ~ data.scores.nonat[field.nonat.al2,1], pch=1,
#        cex=1.5, col=rgb(red=0, green = 0, blue = 1))
points(data.scores.nonat[field.nonat.ob2,2] ~ data.scores.nonat[field.nonat.ob2,1], pch=19,
       cex=1.5, col=rgb(red=.94, green = .52, blue = 0))
Plot_ConvexHull(xcoord = data.scores.nonat[field.nonat.ob2,1], ycoord = data.scores.nonat[field.nonat.ob2,2], 
                lcolor = rgb(red=.94, green = .52, blue = 0))
points(data.scores.nonat[field.nonat.bo2,2] ~ data.scores.nonat[field.nonat.bo2,1], pch=19,
       cex=1.5, col="black")
Plot_ConvexHull(xcoord = data.scores.nonat[field.nonat.bo2,1], ycoord = data.scores.nonat[field.nonat.bo2,2], 
                lcolor = "black")


# #Apply a statistical test to contrasts of all groups:


##start copy here for function pairwise.adonis()
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

ord.stats.table <- pairwise.adonis(field.bacs.nonat.ilr.mat, field.nonat$Group, sim.method = 'euclidean', p.adjust.m = 'fdr')
write.csv(ord.stats.table,"pairwise_permmanova_results.csv")

