### Metabarcoding as a Tool to Examine Cryptic Algae in the
### Diets of Two Common Grazing Surgeonfishes, Acanthurus triostegus
### and A. nigrofuscus ###
### Eileen M. Nalley ### 
### April 2021 ### 

#### LOADING PACKAGES FOR BIOINFORMATICS ####
rm(list=ls())
library("devtools")
devtools::install_github("benjjneb/dada2"); library("dada2")
install.packages("ape"); library("ape")
devtools::install_github("shaunpwilkinson/insect"); library("insect")
install.packages("aphid"); library("aphid")

#### LOADING PACKAGES FOR ANALYSIS ####
library(ape); library(bipartite); library(car); library(cluster);
library(data.table);
library(DESeq2); library(dplyr); library(effects); library(eulerr);
library(funrar);
library(ggmap); library(ggplot2); library(ggpubr); library(goeveg);
library(iNEXT); library(knitr); library(labdsv); library(lattice);
library(lme4);
library(lmerTest); library(MASS); library(metagMisc);
library(microbiome);
library(modelr); library(MuMIn); library(pbkrtest); library(phyloseq);
library(plotly);
library(prodlim);  library(randomcoloR);  library(RColorBrewer);
library(rcompanion); library(RInSp);
library(rsample); library(reshape2); library(scales); library(sjPlot);
library(tibble);
library(tidyverse); library(vegan); library(vsn);
library(wesanderson); library(psycho)

## setting path
## Raw sequences available on GitHub
path_all <- "/AllSequences_NoZero"
## checking that files are in the path
list.files(path_all)

#### FILTERING ####
## filter forward and reverse reads separately
fnFs_all <- sort(list.files(path_all,pattern="_R1.fastq.gz",full.names
                            = TRUE))
fnRs_all <- sort(list.files(path_all,pattern="_R2.fastq.gz",full.names
                            = TRUE))
length(fnFs_all)
length(fnRs_all)
saveRDS(fnFs_all, file = "/fnFs_all.rds")
saveRDS(fnRs_all, file = "/fnRs_all.rds")
## extract sample names
sample.names_all <- sapply(strsplit(basename(fnFs_all)))
saveRDS(sample.names_all, file = "/sample.names_all.rds")
## inspect read quality profiles by ploting read length vs. quality
score
QSforward_all <- plotQualityProfile(fnFs_all, n=2, aggregate = TRUE)
saveRDS(QSforward_all, file = "/QSforward_all.rds")
#
QSreverse_all <- plotQualityProfile(fnRs_all, n=2, aggregate = TRUE)
saveRDS(QSreverse_all, file = "/QSreverse_all.rds")
# plotting
plot(QSforward_all)
plot(QSreverse_all)

#### TRIM ####
# assign the filenames for the filtered forward and reverse files
filtFs_all <-
  file.path(path_all,"filtered",paste0(sample.names_all,"_F_filt.fastq.g
                                       z"))
saveRDS(filtFs_all, file = "/filtFs_all.rds")
#
filtRs_all <-
  file.path(path_all,"filtered",paste0(sample.names_all,"_R_filt.fastq.g
                                       z"))
saveRDS(filtRs_all, file = "/filtRs_all.rds")
# filter reads
out_all <-
  filterAndTrim(fnFs_all,filtFs_all,fnRs_all,filtRs_all,truncLen=c(275,2
                                                                   50),
                maxN=0,maxEE=c(2,2),truncQ=2,rm.phix = TRUE,
                compress=TRUE,multithread=TRUE, trimLeft =
                  c(20, 20))
out_all ## checking output
saveRDS(out_all, file = "/out_all.rds")
## creating a CSV file of the data for reference
write.csv(out_all, file = "/FilteredReads_AllSequences.csv")
#### DEREPLICATE ####
sum(file.size(filtFs_all))
#
derepFs_all <- derepFastq(filtFs_all,verbose=TRUE)
saveRDS(derepFs_all, file = "/derepFs_all.rds")
#
derepRs_all <- derepFastq(filtRs_all,verbose=TRUE)
saveRDS(derepRs_all, file = "/derepRs_all.rds")
#### LEARN ERROR RATES ####
# have to remove all zeros in first steps
errF_all <- learnErrors(derepFs_all, multithread=TRUE)
saveRDS(errF_all, file = "/errF_all.rds")
#
errR_all <- learnErrors(derepRs_all, multithread=TRUE)
saveRDS(errR_all, file = "/errR_all.rds")
#
errFplot_all <- plotErrors(errF_all, nominalQ=TRUE)
errRplot_all <- plotErrors(errR_all, nominalQ=TRUE)
#### PSUEDOPOOLING DEREPLICATED DATA ####
dadaFsPool_all <- dada(derepFs_all,err=errF_all,multithread=TRUE, pool
                       = "pseudo")
saveRDS(dadaFsPool_all, file = "/dadaFsPool_all.rds")
dadaRsPool_all <- dada(derepRs_all,err=errR_all,multithread=TRUE, pool
                       = "pseudo")
saveRDS(dadaRsPool_all, file = "/dadaRsPool_all.rds")


#### MERGE READS ####
mergersPool_all <-  mergePairs(dadaFsPool_all, derepFs_all,
                               dadaRsPool_all, derepRs_all, verbose=TRUE)
# check that they merged properly
head(mergersPool_all)
write.csv(mergersPool_all, file = "/mergersPool_all.csv")

#### CREATE TABLE OF SEQUENCES ####
seqtabPool_all <- makeSequenceTable(mergersPool_all)
saveRDS(seqtabPool_all, file = "/seqtabPool_all.rds")
# summarize table --> number of samples and number of sequences
dim(seqtabPool_all)
# table of sequence lengths
table(nchar(getSequences(seqtabPool_all)))
# remove chimeras
seqtab.nochimPool_all <- removeBimeraDenovo(seqtabPool_all,
                                            method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochimPool_all)
# assess percentage of data that included chimeras
sum(seqtab.nochimPool_all)/sum(seqtabPool_all)

# track your sequence loss thus far in the pipeline
getN <- function(x) sum(getUniques(x))
trackPool_all <- cbind(out_all, sapply(dadaFsPool_all, getN),
                       sapply(dadaRsPool_all, getN), sapply(mergersPool_all, getN),
                       rowSums(seqtab.nochimPool_all))
colnames(trackPool_all) <- c("input", "filtered", "denoisedF",
                             "denoisedR", "merged", "nonchim")
rownames(trackPool_all) <- sample.names_all
trackPool_all
write.table(trackPool_all, "/23S_dada2_tracked_data_Poolall.txt",
            sep="\t")


# get data ready to classify with INSECT
xPool_all <- char2dna(colnames(seqtab.nochimPool_all))
names(xPool_all) <- paste0("ASV", seq_along(xPool_all))
tree <- readRDS("/classifier.rds")
names(attributes(tree))
#### TAXONOMIC ASSIGNMENT USING INSECT ####
longDFPool_all <- classify(xPool_all, tree, cores = 4, tabulize =
                             FALSE)
saveRDS(longDFPool_all, file = "/longDFPool_all.rds")
# integrate the data from DADA2
longDFPool_all2 <- cbind(longDFPool_all, t(seqtab.nochimPool_all))
write.csv(longDFPool_all2, file = "/longDFPool_all2.csv")
saveRDS(longDFPool_all2, file = "/longDFPool_all2.rds")
# make a more user friendly version of taxonomic assignment results
taxaPool_all <- aggregate(longDFPool_all2[3:12],
                          longDFPool_all2["taxID"], head, 1)
countsPool_all <- aggregate(longDFPool_all2[13:ncol(longDFPool_all2)],
                            longDFPool_all2["taxID"], sum)
shortDFPool_all <- merge(taxaPool_all, countsPool_all, by = "taxID")
write.csv(shortDFPool_all, file = "/shortDFPool_all.csv")

### ANALYSIS OF DATA ###
#### Importing, Cleaning & Organizing Data ####
dietPool <- as.data.frame(read.csv("shortDFPool.csv", header=TRUE))
str(dietPool)
#
metadata <- as.data.frame(read.csv("SpcsSite.csv", header=TRUE))
str(metadata)
rownames(metadata) <- metadata[,c("Sample_orig")]
## getting data into the same order
dietPool2 <- aggregate(dietPool[,c(10:ncol(dietPool))],dietPool["taxon"],sum)
ord.dietPool <- dietPool2
#
rownames(ord.dietPool) <- ord.dietPool$taxon
ord.dietPool <- ord.dietPool[,-c(1)]
ord.dietPool.t <- as.data.frame(t(ord.dietPool))
ord.dietPool.md <- ord.dietPool.t[match(rownames(metadata),
                                        rownames(ord.dietPool.t)), ]

### renaming columns now to be the same as TechnicalReplicates in
metadata
## cannot have duplicate row names so need to transpose
metadata.tr <- t(metadata)
colnames(metadata.tr) <- metadata.tr[2,]
ord.dietPool.tr <- as.data.frame(t(ord.dietPool.md))
colnames(ord.dietPool.tr) <- colnames(metadata.tr)

### Adding together TechnicalReplicates ###
ord.dietPool.tr2 <- t(rowsum(t(ord.dietPool.tr), group =
                               rownames(t(ord.dietPool.tr))))
sum(ord.dietPool.tr2[,35])
total.reads <- as.data.frame(colSums(ord.dietPool.tr2))
summary(total.reads)


### Selecting samples with >10000 reads ###
ord.dietPool.tr.10000 = ord.dietPool.tr2[,colSums(ord.dietPool.tr2) >
                                           10000]

### selecting rows in metadata that are in asv table ###
metadata.tr.10000 <- as.data.frame(t(metadata.tr))
## and removing rows from metadata so that just 1 per technical rep
metadata.tr.10000 <- metadata.tr.10000[!duplicated(metadata.tr.
                                                   10000$TechnicalReplicates),]
rownames(metadata.tr.10000) <- metadata.tr.10000$TechnicalReplicates
## only selecting samples that are in diet data file
metadata.tr.10000 <- subset(metadata.tr.10000, TechnicalReplicates
                            %in% colnames(ord.dietPool.tr.10000))
str(metadata.tr.10000)

### dealing with biological reps ###
## 1st getting richness & diversity & adding to df
adiv <- specnumber(t(ord.dietPool.tr.10000))
simp <- diversity(t(ord.dietPool.tr.10000), index = "simpson")
adivmd <- cbind(metadata.tr.10000, adiv, simp)

### choose rows where unique bio rep has highest diversity value ###
biorep_summary <- adivmd %>%
  group_by(BiologicalReplicate) %>%
  summarise(maxsimp = max(simp),
            minsimp = min(simp),
            maxrich = max(adiv),
            minrich = min(adiv),
            reps = n_distinct(TechnicalReplicates))
biorep_hisimp <- adivmd %>% ## setting df
  group_by(BiologicalReplicate) %>% ## grouping by biological rep
  arrange(desc(simp)) %>% ## ordering bio rep groups by dec simp div
  top_n(1) ## selecting top row of each biorep group (i.e., row with
max simp)
biorep_hirich <- adivmd %>%
  group_by(BiologicalReplicate) %>%
  arrange(desc(adiv)) %>%
  top_n(1)
### now matching asv df to metadata by technical replicates ###
asv_biorep_hisimp <- as.data.frame(t(ord.dietPool.tr.10000))
asv_biorep_hisimp$TechnicalReplicates <- rownames(asv_biorep_hisimp)
asv_biorep_hisimp <- asv_biorep_hisimp %>%
  semi_join(biorep_hisimp, by = "TechnicalReplicates")
str(asv_biorep_hisimp)
asv_biorep_hirich <- as.data.frame(t(ord.dietPool.tr.10000))
asv_biorep_hirich$TechnicalReplicates <- rownames(asv_biorep_hirich)
asv_biorep_hirich <- asv_biorep_hirich %>%
  semi_join(biorep_hirich, by = "TechnicalReplicates")
str(asv_biorep_hirich)
spcsmn <- biorep_hirich %>%
  group_by(SpeciesCode) %>%
  summarise(rich.mn = mean(adiv),
            rich.sd = sd(adiv),
            rich.max = max(adiv),
            rich.min = min(adiv),
            simp.mn = mean(simp),
            simp.sd = sd(simp),
            simp.max = max(simp),
            simp.min = min(simp),
            samp.size = n_distinct(BiologicalReplicate))
spcsmn 
#### Phyloseq ####
## putting into phyloseq format
## taxa table
tax_data <- dietPool[,c(1:9)]
tax_data2 <- semi_join(tax_data, dietPool2, by = "taxon")
tax_data3 <- tax_data2[!duplicated(tax_data2$taxon),]
rownames(tax_data3) <- tax_data3$taxon
tax_data4 <- tax_data3[,c(3:ncol(tax_data3))]
colnames(tax_data4) <- c("Kingdom", "Phylum", "Class", "Order",
                         "Family", "Genus", "Species")
taxmat <- as.matrix(tax_data4)
TAX = tax_table(taxmat)

## asv table
rownames(asv_biorep_hirich) <- asv_biorep_hirich$TechnicalReplicates
asv_phyloseq <- as.data.frame(t(subset(asv_biorep_hirich, select=-
                                         c(TechnicalReplicates))))
str(asv_phyloseq) # 138 obs. of  215 variables
asv_table <- asv_phyloseq # 138 obs. of  215 variables
str(asv_table)
otumat <- as.matrix(asv_table)
OTU = otu_table(otumat, taxa_are_rows = TRUE)

## combining to make physeq object
physeq = phyloseq(OTU, TAX)
physeq

## metadata
sampledata <- biorep_hirich
rownames(sampledata) <- biorep_hirich$TechnicalReplicates
str(sampledata)
sampledata2 = sample_data(data.frame(sampledata))

## tree
random_tree = rtree(ntaxa(physeq), rooted=TRUE,
                    tip.label=taxa_names(physeq))
## merging
physeq1 = merge_phyloseq(physeq, sampledata2, random_tree)
physeq1

#### Rarefaction ####
physeq.rareT = rarefy_even_depth(physeq1, replace=TRUE, rngseed =
                                   TRUE)
physeq.rareF = rarefy_even_depth(physeq1, replace=FALSE, rngseed =
                                   TRUE)
write_phyloseq(physeq.rareT, type = "OTU", path = "/")
write_phyloseq(physeq.rareT, type = "TAXONOMY", path = "/")
write_phyloseq(physeq.rareT, type = "METADATA", path = "/")

#### Preprocessing ####
physeq.RRA <- filter_taxa(physeq.rareT, function(x) mean(x) > 1e-5,
                          TRUE)
physeq.RRA <- transform_sample_counts(physeq.RRA, function(x) x /
                                        sum(x) )
write_phyloseq(physeq.RRA, type = "OTU", path = "/")

#### Prepping Data ####
OTU.IS <- as.data.frame(read.csv("otu_table_headers.csv"))
str(OTU.IS)
rownames(OTU.IS) <- OTU.IS[,1]
OTU.IS <- OTU.IS[,-which(colnames(OTU.IS)=='ACLE007A')] # outlier
OTU.IS <- OTU.IS[,-which(colnames(OTU.IS)=='ACTR078')] # outlier
OTU.IS <- OTU.IS[,-which(colnames(OTU.IS)=='ACLE002C')] # outlier
str(OTU.IS)
metadata.IS <- as.data.frame(read.csv("metadata_table_headers.csv"))
rownames(metadata.IS) <- metadata.IS$TechnicalReplicates
metadata.IS$Site[metadata.IS$Site == "ChinaWalls_OAH"] <-
  "Maunalua_OAH"
metadata.IS <- metadata.IS[!(metadata.IS$Site=="UNKNOWN"),]
metadata.IS <- metadata.IS[-which(rownames(metadata.IS)=='ACTR078'),]
metadata.IS <- metadata.IS[-which(rownames(metadata.IS)=='ACLE002C'),]

### check that taxa data file is ordered correctly before importing!
taxa.IS <- as.data.frame(read.csv("taxonomy_table_headers.csv"))
rownames(taxa.IS) <- taxa.IS$taxon
str(taxa.IS)
OTU.IS2 <- OTU.IS[ order(match(OTU.IS$taxon, taxa.IS$taxon)), ]
OTU.IS2 <- OTU.IS2[,-1]
#
OTU.IS.sum <- as.data.frame(t(OTU.IS2))
OTU.IS.sum$TechnicalReplicates <- rownames(OTU.IS.sum)
summaryall <- full_join(OTU.IS.sum, metadata.IS, by =
                          "TechnicalReplicates")
rownames(summaryall) <- summaryall$TechnicalReplicates


summaryall2 <- summaryall %>%
  group_by(Site, SpeciesCode) %>%
  summarise(count = n_distinct(TechnicalReplicates))
summaryall2

## create a list of all the taxa that are found in each species ###
alldietitems <- summaryall[,c(137,1:131)] %>%
  group_by(SpeciesCode) %>%
  summarise_each(funs(sum))
alldietitems2 <- as.data.frame(alldietitems[,2:ncol(alldietitems)])
rownames(alldietitems2) <- alldietitems$SpeciesCode
t.alldietitems <- t(alldietitems2)
write.csv(t.alldietitems, file = "/AllDietItems.csv")

#### Examining Differences in size, richness, diversity in ACTR & ACNF ####
summaryall_acnfactr <- summaryall %>%  
  filter(SpeciesCode == "ACTR" | SpeciesCode == "ACNF")
#
#### ACTR #### 
# ACTR Specific Plots  (NMDS + Diet Composition) #
ACTR.IS <- as.data.frame(t(OTU.IS2[,grepl("ACTR", colnames(OTU.IS2))])) #89
ACTR.md <- as.data.frame(metadata.IS[ which(metadata.IS$SpeciesCode=='ACTR'), ])
all.ACTR <- cbind(ACTR.md, ACTR.IS) # 89 OBS OF 157 VARS

### adding in new variables for ACTR ###
ACTR.md$FultonsK <- ACTR.md$Wt_g/((ACTR.md$TL_cm)^3)
ACTR.md$Liver_Wt_g <- as.numeric(ACTR.md$Liver_Wt_g)
ACTR.md$HepatosomaticIndex <- 
  (ACTR.md$Liver_Wt_g)/(ACTR.md$Wt_g)
ACTR.md$StCtns_Wt_g <- as.numeric(ACTR.md$StCtns_Wt_g)
ACTR.md$GutFullness <- 
  (ACTR.md$StCtns_Wt_g)/(ACTR.md$Wt_g)

### create a list of all the taxa that are found in the populations by
site ###
testactr <- all.ACTR[,c(7,22:ncol(all.ACTR))] %>%
  group_by(Site) %>%
  summarise_each(funs(sum))
t.testactr <- t(testactr)
write.csv(t.testactr, file = "/ACTRPopDietList.csv")
# NMDS ACTR #
### without Salt Pond, Kauai ###
ACTR.IS2 <- ACTR.IS[c(1:84),]
ACTR.md2 <- ACTR.md[c(1:84),]
dimcheckMDS(ACTR.IS2, distance = "bray", k = 6, trymax = 20,
            autotransform = FALSE)
ACTR.mds.dietDiv <- metaMDS(ACTR.IS2, distance = "bray", k = 4,
                            trymax = 100, autotransform = FALSE)
ACTR.mds.dietDiv # k = 4, stress = 0.091
#
ACTR.mds.dietDiv.spp <- envfit(ACTR.mds.dietDiv, ACTR.IS2,
                               permutations = 999)
ACTR.mds.dietDiv.spp
#
ACTR.mds.MD.vars <- ACTR.md2[,c(7, 9:10, 20:21)]
ACTR.mds.MD.vars <- ACTR.mds.MD.vars %>%
  dplyr::rename(
    Weight = Wt_g,
    Length = TL_cm,
    'Diet Diversity' = simp,
    'Diet Richness' = adiv
  )
ACTR.mds.MD.vars <- envfit(ACTR.mds.dietDiv, ACTR.mds.MD.vars,
                           ACTR.mds.MD.vars
                           permutations = 999, na.rm = TRUE)

## generating site scores
ACTR.DDsite.scrs <- as.data.frame(scores(ACTR.mds.dietDiv, display =
                                           "sites"))
ACTR.DDsite.scrs <- cbind(ACTR.DDsite.scrs, Location = ACTR.md2$Site)
## species scores
ACTR.mds.dietDiv.spp.scrs <-
  as.data.frame(scores(ACTR.mds.dietDiv.spp, display = "vectors"))
ACTR.mds.dietDiv.spp.scrs <- cbind(ACTR.mds.dietDiv.spp.scrs, Taxa =
                                     rownames(ACTR.mds.dietDiv.spp.scrs))
ACTR.mds.dietDiv.spp.scrs <- cbind(ACTR.mds.dietDiv.spp.scrs, pval =
                                     ACTR.mds.dietDiv.spp$vectors$pvals)
sig.ACTR.mds.dietDiv.spp.scrs <- subset(ACTR.mds.dietDiv.spp.scrs,
                                        pval<=0.01)
head(ACTR.mds.dietDiv.spp.scrs)
## ecological variable scores
ACTR.env.scores <- as.data.frame(scores(ACTR.mds.MD.vars, display =
                                          "vectors"))
ACTR.env.scores <- cbind(ACTR.env.scores, variables =
                           rownames(ACTR.env.scores))
ACTR.env.scores <- cbind(ACTR.env.scores, pval =
                           ACTR.mds.MD.vars$vectors$pvals)
sig.ACTR.env.scores <- subset(ACTR.env.scores, pval<=0.05)
head(ACTR.env.scores)

### plot
ACTR.DDnmds.plot <- ggplot(ACTR.DDsite.scrs, aes(x=NMDS1, y=NMDS2,
                                                 colour = factor(ACTR.DDsite.scrs$Location))) +
  geom_point(aes(NMDS1, NMDS2, colour =
                   factor(ACTR.DDsite.scrs$Location)), size = 2)+
  coord_fixed()+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black",
                                        size = 1, linetype = "solid"))+
  labs(colour = element_blank(), shape = element_blank())+
  theme(legend.position = "none", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text =
          element_text(size = 10)) +
  stat_ellipse(linetype = 2) +
  scale_color_manual(labels = c("Electric Beach, Oahu", "Puena Point,
                                Oahu", "Kahana Bay, Oahu", "Kailua Bay, Oahu", "Kakaako, Oahu", "Kaneohe Bay, Oahu",
                                "Makaha, Oahu", "Maunalua, Oahu", "Waihee, Maui"), values = c("tomato4", "palevioletred1", "plum2", "orange2", "gold", "springgreen", "darkcyan", "darkblue", "purple"))
ACTR.DDnmds.plot

### with sig species plotting on top
ACTR.DDnmds.plot2 <- ggplot(ACTR.DDsite.scrs, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour =
                   factor(ACTR.DDsite.scrs$Location)), size = 2)+
  coord_fixed()+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black",
                                        size = 1, linetype = "solid"))+
  labs(colour = 'Species', shape = element_blank())+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text =
          element_text(size = 10)) +
  geom_segment(data = sig.ACTR.mds.dietDiv.spp.scrs,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               colour = "grey10", lwd=0.5) + ggrepel::geom_text_repel(data = sig.ACTR.mds.dietDiv.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Taxa), size = 3.5, cex = 3, direction = "both", segment.size = 0.25) + scale_color_manual(labels = c("Electric Beach", "Puena Point", "Kahana Bay", "Kailua Bay", "Kakaako", "Kaneohe Bay", "Makaha", "Maunalua", "Waihee"), values = c("tomato4", "palevioletred1", "plum2", "orange2", "gold", "springgreen", "darkcyan", "darkblue", "purple"))
ACTR.DDnmds.plot2
#
ggarrange(ACTR.DDnmds.plot, ACTR.DDnmds.plot2, ncol = 2, nrow = 1, align = "hv") 
### with both ecological & environmental variables
ggplot(ACTR.DDsite.scrs, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(NMDS1, NMDS2, colour =
                   factor(ACTR.DDsite.scrs$Location)), size = 2)+
  coord_fixed()+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black",
                                        size = 1, linetype = "solid"))+
  labs(colour = element_blank(), shape = element_blank())+
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text =
          element_text(size = 10))  +
  geom_segment(data = sig.ACTR.env.scores,
               aes(x = 0, xend=NMDS1, y=0, yend=NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), size = 4,
               colour = "grey10", lwd=0.5) +
  ggrepel::geom_text_repel(data = sig.ACTR.env.scores,
                           aes(x=NMDS1, y=NMDS2, label = variables),
                           cex = 4, direction = "both", segment.size = 0.25)+ scale_color_manual(labels = c("Electric Beach", "Puena Point",
                                                                                                            "Kahana Bay", "Kailua Bay",

ACTRall.rich.hist <- ggplot(all.ACTR, aes(all.ACTR$adiv)) +
  geom_histogram(breaks=seq(20, 65, by = 2), col = "darkcyan", fill =
                   "orange2", alpha = 0.2) +
  geom_vline(aes(xintercept = mean(adiv)), col = "darkcyan", size = 1, linetype = "dashed") + labs(x = "Diet Richness", y = "") +  theme_minimal() +
  theme(text = element_text(size=10))
ACTRall.rich.hist

ACTRall.div.hist <- ggplot(all.ACTR, aes(all.ACTR$simp)) + 
  geom_histogram(breaks=seq(0.4, 1, by = 0.0267), col = "sienna4", fill = "lightcoral",
                 alpha=0.2) + 
  geom_vline(aes(xintercept = mean(simp)), col = "sienna4", size = 1, linetype = "dashed") + labs(x = "Diet Diversity (Simpson's Diversity Index)", y = "Number of Individuals") +  
  theme_minimal() +
  theme(text = element_text(size=14))
ACTRall.div.hist


#### ACNF ####
# ACNF Specific Plots  (NMDS + Diet Composition) #
ACNF.IS <- as.data.frame(t(OTU.IS2[,grepl("ACNF", colnames(OTU.IS2))])) # 73
ACNF.md <- as.data.frame(metadata.IS[ which(metadata.IS$SpeciesCode=='ACNF'), ])
all.ACNF <- cbind(ACNF.md, ACNF.IS) # 73 obs of 152 vars

### adding in new variables for ACNF ###
ACNF.md$FultonsK <- ACNF.md$Wt_g/((ACNF.md$TL_cm)^3)
ACNF.md$Liver_Wt_g <- as.numeric(ACNF.md$Liver_Wt_g)
ACNF.md$HepatosomaticIndex <- 
  (ACNF.md$Liver_Wt_g)/(ACNF.md$Wt_g)
ACNF.md$StCtns_Wt_g <- as.numeric(ACNF.md$StCtns_Wt_g)
ACNF.md$GutFullness <- 
  (ACNF.md$StCtns_Wt_g)/(ACNF.md$Wt_g)

### COMPARING ACTR & ACNF
Wilcox.FultonsK <- wilcox.test(ACTR.md2$FultonsK, ACNF.md$FultonsK, paired=FALSE)
Wilcox.FultonsK 
Wilcox.Hepatosomatic <- wilcox.test(ACTR.md2$HepatosomaticIndex, ACNF.md$HepatosomaticIndex,
                                    paired = FALSE)
Wilcox.Hepatosomatic 
Wilcox.GutFullness <- wilcox.test(ACTR.md2$GutFullness, ACNF.md$GutFullness, paired=FALSE)
Wilcox.GutFullness
Wilcox.Wt <- wilcox.test(ACTR.md2$Wt_g, ACNF.md$Wt_g, paired=FALSE)
Wilcox.Wt

### COMPARING BETWEEN SITES WITHIN SPECIES
RichACTR <- pairwisePercentileTest(adiv ~ Site, 
                                   data = summaryall_acnfactr[which(summaryall_acnfactr$SpeciesCode == "ACTR"),], 
                                   test = "mean", r = 10000)
RichACTR 
cldList(p.adjust ~ Comparison, data = RichACTR) # significant differences
DivACTR <- pairwisePercentileTest(simp ~ Site,
                                  data = summaryall_acnfactr[which(summaryall_acnfactr$SpeciesCode == "ACTR"),], 
                                  test = "mean", r = 10000)
DivACTR
cldList(p.adjust ~ Comparison, data = DivACTR) # significant differences
RichACNF <- pairwisePercentileTest(adiv ~ Site, 
                                   data = summaryall_acnfactr[which(summaryall_acnfactr$SpeciesCode == "ACNF"),], 
                                   test = "mean", r = 10000)
RichACNF
cldList(p.adjust ~ Comparison, data = RichACNF) # no significant differences
DivACNF <- pairwisePercentileTest(simp ~ Site, 
                                  data = summaryall_acnfactr[which(summaryall_acnfactr$SpeciesCode == "ACNF"),], 
                                  test = "mean", r = 10000)
DivACNF
cldList(p.adjust ~ Comparison, data = DivACNF) # no significant differences


#### Generating Summary Stats for Species ####
summarystats <- summaryall_acnfactr %>%
  group_by(SpeciesCode) %>%
  summarise(MnWt = mean(Wt_g), MaxWt = max(Wt_g), MnTL = mean(TL_cm),
            max(adiv),
            max(simp),
            MaxTL = max(TL_cm), MnRich = mean(adiv), MaxRich =
              MinRich = min(adiv), MnSimp = mean(simp), MaxSimp =
              MinSimp = min(simp), count=n())


