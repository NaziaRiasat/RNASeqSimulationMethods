#intall
install.packages(c("factoextra", "FactoMineR"))
#intall
install.packages(c("factoextra", "FactoMineR"))

#load
library("factoextra")
library("FactoMineR")
#load
library("factoextra")
library("FactoMineR")
library(tidyverse)
library(edgeR)
library(readxl) 
library(SPsimSeq)
library(devtools)
library(ssizeRNA)
library(DESeq2)
library(VennDiagram)
library(limma)
library(edgeR)
library(qvalue)
library(ggplot2)
library(magrittr)
set.seed(25081988)
library(ggfortify)
library(tibble)
library(tidyr)

setwd("D:/Users/Nazia Khan/Documents/Fall 2021/Home_Work/Project")

# read the data file
Counts<- read.csv("project_read_counts.csv", row.names=1)
Annot <- read.csv("gene_annotations.csv")

##first 5 rows and first 5 columns of count matrix
Counts[1:5,1:24]

res1 <- t(Counts)

pca.res.orig <- prcomp(res1, scale. = FALSE,center = TRUE)
summary(pca.res.orig)

pca.res.orig1 <- prcomp(res1, scale. = TRUE,center = FALSE)
summary(pca.res.orig1)

#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
pca.pro <- ((pca.res.orig$sdev^2) / (sum(pca.res.orig$sdev^2)))*100

barplot(pca.pro, cex.names=1, 
        xlab=paste("Principal component (PC), 1-", 
        length(pca.res.orig$sdev)),
        ylab="Proportion of variation (%)", 
        main="Scree plot", ylim=c(0,100))


##gene IDs

dim(Counts)
dim(Annot)


##Q2-b 

##Get gene annotations
colnames(Annot) <- c("Ensembl_ID","ENTREZ_ID")
rownames(Annot) <- Annot$Ensembl_ID


##Determine which genes have defined annotations
Genes <- Annot[rownames(Counts),]

head(Genes)
dim(Genes)


##keep only expressed genes with defined annotations
hasannot <- rowSums(is.na(Genes))==0
table(hasannot)

##Get rid of genes with no defined annotations
Counts1 <- Counts[hasannot,13:24]
dim(Counts1)
head(Counts1)

Genes1 <- Annot[rownames(Counts1),]
dim(Genes1)

##create grouping variable 
Group <- factor(c(rep("Nox2KO_PBS",6), rep("Nox2KO_LPS",6))) #, rep("Nox2KO_PBS",6), rep("Nox2KO_LPS",6)))

##Q3-a

##e the DGEList function to create an object with the (filtered) read counts, 
##gene names, and group indicators


##Form a DGEList object combining read counts and associated annotation
y <- DGEList(counts=Counts1, genes=Genes1, group=Group)
y$counts[1:5,1:12]
y$counts1 <- y$counts[sample(nrow(y$counts), 25000), ]
dim(y$counts)
head(y$samples)
head(y$genes[1:5,1:2])

##Q3-b 

## filter out the genes with low read counts.

##Determine which genes are expressed in a worthwhile number of samples.
isexpr <- filterByExpr(Counts1, group=Group)
table(isexpr)

##Determine which genes have defined annotations
hasannot1 <- rowSums(is.na(y$genes))==0
table(hasannot1)


##Get rid of genes with no defined annotations
Counts1 <- Counts1[hasannot1,]
dim(Counts)

##Keep only expressed genes with defined annotation and 
##recompute library sizes
y <- y[isexpr & hasannot1, keep.lib.sizes=FALSE]
dim(y$counts)
head(y)


##  Create bar plot of library sizes for experimental units

barplot(y$samples$lib.size*1e-6, names=1:24,
        main="Bar plot of library sizes for Experimental Units",
        xlab = NULL, ylab="Library size (millions)", border="Blue")


##  calculating the normalizing constant 
##  for each experimental unit using TMM normalization


##Apply TMM normalization
y1 <- calcNormFactors(y)
head(y1)
head(y$samples)
colnames(y$samples)
row.names(y$samples)


##Estimate dispersion parameters
y <- estimateDisp(y, robust=TRUE)

levels(y$samples$group)
head(y$counts)
dim(y$counts)
names(y)

dy1 <- estimateCommonDisp(y)
dy <- estimateTagwiseDisp(y)
dispy <- dy$tagwise.dispersion
head(dispy)
length(dispy)
hist(dispy)

#######################PCA ANALYSIS Results without TMM Normalization##################
res <- t(y$counts)
head(res)
colnames(res)
row.names(res)

res[1:5,1:5]

pca.res <- prcomp(res, scale. = TRUE,center = TRUE)
pca.ress <- prcomp(res, scale. = TRUE,center = FALSE)
pca.res.real <- prcomp(res, scale. = FALSE,center = TRUE)
pca.res.reall <- prcomp(res, scale. = FALSE,center = FALSE)


summary(pca.res)
summary(pca.ress)
summary(pca.res.real)
summary(pca.res.reall)

#Results-1
pc_evalues.res <- pca.res$sdev^2

pc_evalues.res <- tibble(PC = factor(1:length(pc_evalues.res)), 
                         variance = pc_evalues.res) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.res

#screeplot 
pc_evalues.res %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca.res, addlabels = TRUE, ylim = c(0, 70))



# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- pca.res$x

#Data frame
pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores

#PC plot
pc_scores %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings <- pca.res$rotation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "Gene")

# print the result
pc_loadings

library(ggplot2)
ggplot(pc_scores, aes(x = PC1,y=PC2, col=sample)) #+ geom_point()
ggplot(pc_loadings, aes(x = PC1,y=PC2, col=Genes)) #+ geom_point()

library(tidyverse)
library(dplyr)

top_genes <- pc_loadings %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes

top_loadings <- pc_loadings %>% 
  filter(gene %in% top_genes)

loadings_plot <- ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot

library(ggfortify)
autoplot(pca.res, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca.res$rotation[pca.res$rotation[,"PC1"]==max(pca.res$rotation[,"PC1"]), "PC1", drop=F]



#Results-2

pc_evalues.ress <- pca.ress$sdev^2

pc_evalues.ress <- tibble(PC = factor(1:length(pc_evalues.ress)), 
                         variance = pc_evalues.ress) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.ress


#screeplot 
pc_evalues.ress %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca.ress, addlabels = TRUE, ylim = c(0, 70))


# The PC scores are stored in the "x" value of the prcomp object
pc_scoress <- pca.ress$x

#Data frame
pc_scoress <- pc_scoress %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scoress

#PC plot
pc_scoress %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadingss <- pca.ress$rotation

pc_loadingss <- pc_loadingss %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadingss


top_geness <- pc_loadingss %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_geness

top_loadingss <- pc_loadingss %>% 
  filter(gene %in% top_geness)

loadings_plott <- ggplot(data = top_loadingss) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot

pca.ress$rotation[pca.ress$rotation[,"PC1"]==max(pca.ress$rotation[,"PC1"]), "PC1", drop=F]


#Results-3

pc_evalues.real <- pca.res.real$sdev^2

pc_evalues.real <- tibble(PC = factor(1:length(pc_evalues.real)), 
                          variance = pc_evalues.real) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.real

#screeplot 
pc_evalues.real %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca.res.real, addlabels = TRUE, ylim = c(0, 80))

# The PC scores are stored in the "x" value of the prcomp object
pc_scoresl <- pca.res.real$x

#Data frame
pc_scoresl <- pc_scoresl %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scoresl

#PC plot
pc_scoresl %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadingsl <- pca.res.real$rotation

pc_loadingsl <- pc_loadingsl %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadingsl


top_genesl <- pc_loadingsl %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genesl

top_loadingsl <- pc_loadingsl %>% 
  filter(gene %in% top_genesl)

loadings_plotl <- ggplot(data = top_loadingsl) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plotl

pca.res.real$rotation[pca.res.real$rotation[,"PC1"]==max(pca.res.real$rotation[,"PC1"]), "PC1", drop=F]



#Results-4

pc_evalues.reall <- pca.res.reall$sdev^2

pc_evalues.reall <- tibble(PC = factor(1:length(pc_evalues.reall)), 
                          variance = pc_evalues.reall) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.reall


#screeplot 
pc_evalues.reall %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca.res.real1, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scoresll <- pca.res.reall$x

#Data frame
pc_scoresll <- pc_scoresll %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scoresll

#PC plot
pc_scoresll %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadingsll <- pca.res.reall$rotation

pc_loadingsll <- pc_loadingsll %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadingsll


top_genesll <- pc_loadingsll %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genesll

top_loadingsll <- pc_loadingsll %>% 
  filter(gene %in% top_genesll)

loadings_plotll <- ggplot(data = top_loadingsll) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plotll

pca.res.reall$rotation[pca.res.reall$rotation[,"PC1"]==max(pca.res.reall$rotation[,"PC1"]), "PC1", drop=F]


#######################PCA ANALYSIS Results with TMM Normalization##################
res1 <- t(y1$counts)

res1[1:5,1:5]

pca.res1 <- prcomp(res1, scale. = TRUE,center = TRUE)
pca.ress1 <- prcomp(res1, scale. = TRUE,center = FALSE)
pca.res.real1 <- prcomp(res1, scale. = FALSE,center = TRUE)
pca.res.reall1 <- prcomp(res1, scale. = FALSE,center = FALSE)


summary(pca.res1)
summary(pca.ress1)
summary(pca.res.real1)
summary(pca.res.reall1)

#Results-1
pc_evalues.res1 <- pca.res1$sdev^2

pc_evalues.res1 <- tibble(PC = factor(1:length(pc_evalues.res1)), 
                         variance = pc_evalues.res1) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.res1

#Results-2

pc_evalues.ress <- pca.ress$sdev^2

pc_evalues.ress <- tibble(PC = factor(1:length(pc_evalues.ress)), 
                          variance = pc_evalues.ress) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.ress


#Results-3

pc_evalues.real <- pca.res.real$sdev^2

pc_evalues.real <- tibble(PC = factor(1:length(pc_evalues.real)), 
                          variance = pc_evalues.real) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.real

#Results-4

pc_evalues.reall <- pca.res.reall$sdev^2

pc_evalues.reall <- tibble(PC = factor(1:length(pc_evalues.reall)), 
                           variance = pc_evalues.reall) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.reall



pca.res[1:5,1:5]
autoplot(pca.res, data = y$counts, colour = 'Group',
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

pca.res$rotation[1:5, 1:5]
dim(pca.res$rotation)

pca.res.real$rotation[1:5, 1:5]
dim(pca.res.real$rotation)

pc_loadings <- pca.res$rotation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings

top_genes <- pc_loadings %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes


top_loadings <- pc_loadings %>% 
  filter(gene %in% top_genes)

loadings_plot <- ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend =PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))


# get the PC scores from prcomp object
pca_res$x %>% 
  # convert it to a tibble
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table
  FIXME(sample_info, by = "FIXME") %>% 
  # make the plot
  ggplot(aes(x = PC1, y = PC2, 
             FIXME = factor(minute), shape = FIXME)) +
  geom_point()





library(patchwork)

# Adjust some aspects of each plot
pca_plot <- pca_plot + 
  coord_fixed(ratio = 0.4) + 
  labs(title = "PC scores") +
  theme(legend.position = "none")

loadings_plot <- loadings_plot + 
  coord_fixed(ratio = 0.4) + 
  labs(title = "PC loadings")

# Put them together
(pca_plot | loadings_plot) + plot_annotation(tag_levels = "A")



# The sum of the squared loadings is equal to 1 for each PC
sum(pca.res$rotation[,"PC1"]^2)
sum(pca_res$rotation[,"PC4"]^2)
sum(pca_res$rotation[,"PC12"]^2)

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca.res$rotation[pca.res$rotation[,"PC1"]==max(pca.res$rotation[,"PC1"]), "PC1", drop=F]

pca.res.real$rotation[pca.res.real$rotation[,"PC1"]==max(pca.res.real$rotation[,"PC1"]), "PC1", drop=F]


par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(pca.res$x[,1:6], col="orange", main="Principal components analysis bi-plot1\nPCs 1-6", pch=16)
pairs(pca.res$x[,7:12], col="orange", main="Principal components analysis bi-plot1\nPCs 7-12", pch=16)


par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(pca.res.real$x[,1:6], col="orange", main="Principal components analysis bi-plot1\nPCs 1-6", pch=16)
pairs(pca.res.real$x[,7:12], col="orange", main="Principal components analysis bi-plot1\nPCs 7-12", pch=16)





#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
pca.res.pro <- ((pca.res$sdev^2) / (sum(pca.res$sdev^2)))*100

barplot(pca.res.pro, cex.names=1, 
        xlab=paste("Principal component (PC), 1-", 
                   length(pca.res$sdev)),
        ylab="Proportion of variation (%)", 
        main="Scree plot", ylim=c(0,100))


###############################Exact Test#######################################
###############################################################################

et_WT <- exactTest(y)

names(et_WT)

genes <- row.names(et_WT$table)


head(et_WT$table)



topTags(et_WT, n=10)$table

#p-values
p.23 <- et_WT $table$PValue



head(p.23)
hist(p.23)
##histogram of the p-values between KO_PBS and KO_LPS
##  Create a histogram of p-values for testing for differential expression 
##  between the KO_PBS and KO_LPS 

et_WT$comparison

hist(et_WT$table$PValue, main="Histogram of P-values Between KO_PBS and KO_LPS",
     border="blue")


## Adjusted p-values using B&H procedure
p.bh <- p.adjust(p.23, method="BH")

qvalue <- qvalue(p.23)
hist(qvalue)
head(qvalue)

de.genes.ex <- genes[qvalue$qvalues<=0.05]
head(de.genes.ex)

ee.genes.ex <- genes[qvalue$qvalues>0.05]
head(ee.genes.ex)
length(ee.genes.ex)
length(de.genes.ex)



padj_WT <- topTags(et_WT, n=dim(y$counts)[1])$table[,-4,-7]
head(padj_WT)

dim(padj_WT)


## genes that are DDE when controlling FDR at 5% and abs(LFC > 1)
sum(qvalue$qvalues <= 0.05)


cangene.ex <- (qvalue$qvalues <= 0.05)
non.cangene.ex <- (qvalue$qvalues > 0.05)



head(cangene.ex)
head(non.cangene.ex)

sum(padj_WT$FDR<0.05 & abs(padj_WT$logFC)>1, na.rm=TRUE)


DDE_WT <- padj_WT[padj_WT$FDR<0.05 & abs(padj_WT$logFC)>1,]

name_WT<-rownames(DDE_WT)
head(name_WT)
length(name_WT)


###############################limma-trend#####################################
###############################################################################

##calculate logCPM values
##(edgeR function)

lcpm <- cpm(y, log=TRUE)
head(lcpm)
dim(lcpm)

##Create a design matrix
design <- model.matrix(~0+Group)
design

##Between NO_PBS and NO_LPS

##identify contrast of interest and create corresponding vector
contrasts_KO <- makeContrasts(Nox2KO_PBS.vs.Nox2KO_LPS =
                                GroupNox2KO_PBS - GroupNox2KO_LPS,
                              levels  = design)

##perform limma trend
##(limma functions)
fitt.K <- lmFit(lcpm, design)
fitt.K <- contrasts.fit(fitt.K, contrasts_KO)
efitt.K <- eBayes(fitt.K, trend=TRUE)

##p-values for testing difference in gene expression
p<- efitt.K$p.value

head(efitt.K$p.value)


# get q-value object
qobj <- qvalue(p)
qval <- qobj$qvalue

pio <- qobj$pi0

genes.limma <- row.names(efitt.K)
head(genes.limma)

de.cangene.limma <- genes.limma[qobj$qvalue <= 0.05]
ee.non.cangene.limma <- genes.limma[qobj$qvalue >0.05]

head(ee.cangene.limma)
head(de.non.cangene.limma)
length(ee.cangene.limma)
length(de.non.cangene.limma)
 
head(qval)
head(pio)

# get summary results from q-value object
summary(qobj, cuts=c(0.05))

plot(qobj)
hist(qobj)
hist(qval)


##histogram of the p-values between KO_PBS and KO_LPS
hist(efitt.K$p.value, main="Histogram of P-values Between KO_PBS and KO_LPS",
     border="Red")


##Get top 10 results
top.table<- topTable(efitt.K,n=10)
results.K <- decideTests(efitt.K)

p_K <- (efitt.K$p.value<=0.05)
head(p_K)
sum(p_K <= 0.05)

DDE3 <- p_K[,-1]
head(DDE3)
name_K <- rownames(DDE3)
length(name_K)
name_K[1:10,]

##Get top 10 results


rownames(topTable(efitt.K,n=10))

summary(results.K)

length(which(top.table$adj.P.Val < 0.05))



###############################limma-voom#####################################
###############################################################################


##DE between WT_LPS and NO_LPS

##identify contrast of interest and create corresponding vector
contrasts_KO <- makeContrasts(Nox2KO_PBS.vs.Nox2KO_LPS =
                                GroupNox2KO_PBS - GroupNox2KO_LPS,
                              levels  = design)

v <- voom(y, design, plot=TRUE)

##Between KO_PBS and KO_LPS


fitv <- lmFit(v, design)
fitv.W <- contrasts.fit(fitv, contrasts_KO)
efitv.W <- eBayes(fitv.W)

genes.voom <- row.names(efitv.W)
head(genes.voom)


##Get top 10 results
topTable(efitv.W,n=10)
resultsv.W <- decideTests(efitv.W)

##p-values for testing difference in gene expression
##between males and females
p.v <- efitv.W$p.value

head(efitv.W$p.value)

p_wild <- (p.v<=0.05)
head(p_wild)

# get q-value object
qobj.v <- qvalue(p.v)
qval.v <- qobj.v$qvalue

pio.v <- qobj.v$pi0

pio.v

head(qval.v)

de.cangene.voom <- genes.voom[qval.v <= 0.05]
ee.non.cangene.voom <- genes.voom[qval.v >0.05]

head(ee.non.cangene.voom)
head(de.cangene.voom)
length(de.cangene.voom)
length(ee.non.cangene.voom)


# get summary results from q-value object
summary(qobj.v, cuts=c(0.05))

plot(qobj.v)
hist(qobj.v)
hist(qval.v)


##histogram of the p-values between KO_PBS and KO_LPS
hist(p.v, main="Histogram of P-values Between KO_PBS and KO_LPS",
     border="Red")

##histogram of the q-values between KO_PBS and KO_LPS
hist(qval.v, main="Histogram of Q-values Between KO_PBS and KO_LPS",
     border="Red")


DDE5 <- p_wild[,-1]
head(DDE5)

name_wi <- (DDE5)
name_wi[1:10,]

##Get top 10 results


rownames(topTable(efitv.W,n=10))

summary(decideTests(efitv.W))


################################################################################
######################Power/sample size calculation:############################


## mean counts in control group for all genes
mu <- apply(y$counts[, y$samples$group == "Nox2KO_LPS"], 1, mean)
summary(mu)

mu[mu==0] <- 0.1

head(mu)
length(mu)

## dispersion for all genes
disp <- y$tagwise.dispersion
head(disp)
length(disp)
hist(dispy)


## 2-fold change for DE genes using Exact test
fc <- 2^(et_WT$table$logFC)
length(fc)

hist((et_WT$table$logFC))

lfc1 <- sort((et_WT$table$logFC))[1:2000]  
length(lfc1)
hist(lfc1)

lfc2 <- sort((et_WT$table$logFC), decreasing = TRUE)[1:520] 
length(lfc2)

fc1 <- 2^(c(lfc1,lfc2))
hist(fc1)
length(fc1)

## 2-fold change for DE genes using limma-trend test
x <- topTable(efitt.K,n=12600)
x1 <- topTable(efitt.K,n=2000)
x2 <- topTable(efitt.K,n=520)
head(x)
dim(x)

fc_v <- 2^(x$logFC)
length(fc_v)

hist((fc_v))

lfc1_v <- sort((x$logFC))[1:2000]  
length(lfc1_v)
hist(lfc1_v)

lfc2_v <- sort((x$logFC), decreasing = TRUE)[1:520] 
length(lfc2_v)
hist(lfc2_v)

fc1_v <- 2^(c(lfc1_v,lfc2_v))
hist(fc1_v)
length(fc1_v)

## 2-fold change for DE genes using limma-voom test
z <- topTable(efitt.K,n=12600)
z1 <- topTable(efitt.K,n=2000)
z2 <- topTable(efitt.K,n=520)
head(z)
dim(z)

fc_z <- 2^(z$logFC)
length(fc_z)

hist((fc_z))

lfc1_z <- sort((z$logFC))[1:2000]  
length(lfc1_z)
hist(lfc1_z)

lfc2_z <- sort((z$logFC), decreasing = TRUE)[1:520] 
length(lfc2_z)
hist(lfc2_z)

fc1_z <- 2^(c(lfc1_z,lfc2_z))
hist(fc1_z)
length(fc1_z)



##Average Power and True FDR Based on limma/voom RNAseq Analysis Pipeline

check.power(nGenes = 12600, pi0 = 0.8, m=200, mu=mu, disp=disp, fc=fc1, up = 0.5,
            replace = TRUE, fdr = 0.05, sims = 1)

check.power(nGenes = 12600, pi0 = 0.8, m=200, mu=mu, disp=disp, fc=fc1_v, up = 0.5,
            replace = TRUE, fdr = 0.05, sims = 1)
##pow_bh_ave----average power when controlling FDR by Benjamini and Hochberg (1995) method.
##fdr_bh_ave----true false discovery rate when controlling FDR by Benjamini and Hochberg (1995) method.
##pow_bh_ave----average power when controlling FDR by q-value procedure (Storey et al., 2004).
##fdr_bh_ave----true false discovery rate when controlling FDR by q-value procedure (Storey et
                                          

##RNA-seq Count Data Simulation from Negative-Binomial Distribution
Sim_counts<- sim.counts(nGenes = 12600, pi0 = 0.8, m=30, mu=mu, disp=disp, fc=fc1,
                        up = 0.5,replace = FALSE)
Sim_counts$counts ## count data matrix

##Value
##counts RNA-seq count data matrix.
##group treatment group vector.
##lambda0 mean counts in control group for each gene.
##phi0 dispersion parameter for each gene.
##de differentially expressed genes indicator: 0 for non-differentially expressed genes,1 for up-regulated genes, -1 for down-regulated genes.
##delta log2 fold change for each gene between treatment group and control group.

set.seed(2016)

size <- ssizeRNA_vary(nGenes = 12600, pi0 = 0.65, m = 200, mu=mu, disp=disp, fc=fc1,
                      up = 0.5, replace = FALSE, fdr = 0.05, power = 0.9, maxN = 100)

size$ssize

set.seed(2016)
size <- ssizeRNA_vary(nGenes = 12600, pi0 = 0.65, m = 200, mu=mu, disp=dispy, fc=fc1,
                      up = 0.5, replace = FALSE, fdr = 0.05, power = 0.9, maxN = 100)

size$ssize

size_v <- ssizeRNA_vary(nGenes = 12600, pi0 = 0.6, m = 200, mu=mu, disp=disp, fc=fc1_z,
                        up = 0.5, replace = FALSE, fdr = 0.05, power = 0.8, maxN = 100)

size_v$ssize

size_v <- ssizeRNA_vary(nGenes = 12600, pi0 = 0.8, m = 200, mu=mu, disp=dispy, fc=fc1_v,
                      up = 0.5, replace = FALSE, fdr = 0.05, power = 0.8, maxN = 100)

size_v$ssize## first sample size to reach desired power
size$power ## calculated power for each sample size
size$crit.vals ## calculated critical value for each sample size


#construct the polyomial matrix, using Horner's rule
X <- buildXmat(y$counts, 12)

#Calculates counts per millions of reads, possibly with log-transform
CPM <- calculateCPM(y$counts, 0.1, prior.count=1)

#Check for data validity
checkInputValidity(data,group=Group,batch = NULL,
                   llStat.thrld=0.5,w = 0.5, 
                   pDE=0.1,  prior.count=1, group.config=c(0.5, 0.5),
                   result.format="list")


#Select candidate genes

#This function can be used to independently select candidate genes from a given real RNA-srq data
#(bulk/single) for the SPsimSeq simulation. 

data <- as.matrix(y$AveLogCPM)
data <- as.matrix(y$counts)


obtCount(CPM, 0.5)


cand.DE.gene <- chooseCandGenes(CPM,  group=Group,  lfc.thrld=0.5,
                                llStat.thrld=0.5,t.thrld=0.5,w = NULL, 
                                max.frac.zeror.diff = Inf,
                                pDE=0.2,   n.genes=12600,  prior.count=1)



# SPsimSeq simulation

sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = y$counts,group=Group,
                          n.genes = 12600, batch.config = 1,llStat.thrld = 5,
                          group.config = c(0.5, 0.5), tot.samples = 72,
                          pDE = 0.2, lfc.thrld = 1.5, t.thrld=2.5,w=0.5,
                          model.zero.prob = FALSE, result.format = "list")


sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = y$counts,group=Group,n.genes = 12600,
                          batch.config = 1,group.config = c(0.5, 0.5),
                          tot.samples = 70 ,pDE = 0.2, lfc.thrld = 0.5,
                          result.format = "list")

sim.data.bulk1 <- SPsimSeq(n.sim = 1, s.data = y$counts,group=Group,
                          n.genes = 1260, batch.config = 1,llStat.thrld = 0.5,
                          group.config = c(0.5, 0.5), tot.samples = 80,
                          pDE = 0.1, lfc.thrld = 0.5, t.thrld=2.5,w=NULL,
                          cand.DE.genes = list(c(ee.cangene.limma),c(de.cangene.limma)),
                          model.zero.prob = FALSE, result.format = "list")

sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = y$counts,group=Group,
                          n.genes = 12600, batch.config = 1,llStat.thrld = 5,
                          group.config = c(0.5, 0.5), tot.samples = 100,
                          pDE = 0.2, lfc.thrld = 1.5, t.thrld=2.5,w=0.5,
                          model.zero.prob = FALSE, result.format = "list")


sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = y$counts,group=Group,
                          n.genes = 12600, batch.config = 1,llStat.thrld = 0,
                          group.config = c(0.5, 0.5), tot.samples = 12,
                          pDE = 0.2, lfc.thrld = 0.25, t.thrld=1,w=0.5,
                          model.zero.prob = FALSE, result.format = "list")



sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = y$counts,group=Group,
                          n.genes = 12600, batch.config = 1,llStat.thrld = 0.5,
                          group.config = c(0.5, 0.5), tot.samples = 140,
                          pDE = 0.1, lfc.thrld = 0.5, t.thrld=1,w=0.5,
                          model.zero.prob = FALSE, result.format = "list",
                          return.details = TRUE)



##data Exploration
sim.data.bulk1 <- sim.data.bulk[[1]]                              
head(sim.data.bulk1$counts[, seq_len(8)])  # count data
head(sim.data.bulk1$colData)        # sample info
head(sim.data.bulk1$rowData)        # gene info
summary(sim.data.bulk1)

sim.data.bulk2 <- sim.data.bulk[[1]]                              
head(sim.data.bulk2$counts[, seq_len(5)])  # count data
head(sim.data.bulk2$colData)        # sample info
head(sim.data.bulk2$rowData)        # gene info
summary(sim.data.bulk2)

sim.data.bulk3 <- sim.data.bulk[[1]]                              
head(sim.data.bulk3$counts[, seq_len(5)])  # count data
head(sim.data.bulk3$colData)        # sample info
head(sim.data.bulk3$rowData)        # gene info
summary(sim.data.bulk3)


sim.data.bulk4 <- sim.data.bulk[[1]]                              
head(sim.data.bulk4$counts[, seq_len(5)])  # count data
head(sim.data.bulk4$colData)        # sample info
head(sim.data.bulk4$rowData)        # gene info
summary(sim.data.bulk4)

sim.data.bulk5 <- sim.data.bulk[[1]]                              
head(sim.data.bulk5$counts[, seq_len(5)])  # count data
dim(sim.data.bulk1$counts)
head(sim.data.bulk5$colData)        # sample info
dim(sim.data.bulk5$colData)
head(sim.data.bulk5$rowData)  # gene info
dim(sim.data.bulk5$rowData)
summary(sim.data.bulk5)

counts.SPsimseq <- sim.data.bulk1$counts  # Simulated Count matrix from SpsimSeq
dim (counts.SPsimseq)

genes.samp <- sim.data.bulk1$rowData # Genes sampled from source matrix

samp.col <- sim.data.bulk1$colData# Columns sampled in SimSeq algorithm

###############################################################################
#####################PCA Analysis#############################################
##############################################################################

##Form a DGEList object combining read counts and associated annotation
Datta11 <- DGEList(counts=sim.data.bulk1$counts, genes=sim.data.bulk1$rowData, group=Group)
Datta11$counts[1:5,1:12]
dim(Genes1)
dim(Datta11)

## filter out the genes with low read counts.

##Determine which genes are expressed in a worthwhile number of samples.
isexpr1 <- filterByExpr(Datta11, group=Group)
table(isexpr1)

##Determine which genes have defined annotations
hasannot1w <- rowSums(is.na(Datta11$genes))==0
table(hasannot1w)


##Keep only expressed genes with defined annotation and 
##recompute library sizes
Datta11 <- Datta11[isexpr1 & hasannot1w, keep.lib.sizes=FALSE]
dim(Datta11$counts)
head(y)

##Apply TMM normalization
Datta112 <- calcNormFactors(Datta11)
head(Datta112)


######Data transpose############
a <- t(Datta11$counts)
head(a[1:5,1:5])

##simulation II##


##Form a DGEList object combining read counts and associated annotation
Datta22 <- DGEList(counts=sim.data.bulk2$counts, genes=sim.data.bulk2$rowData, group=Group)
Datta22$counts[1:5,1:12]
dim(Genes1)
dim(Datta22)

## filter out the genes with low read counts.

##Determine which genes are expressed in a worthwhile number of samples.
isexpr12 <- filterByExpr(Datta22, group=Group)
table(isexpr12)

##Determine which genes have defined annotations
hasannot2w <- rowSums(is.na(Datta22$genes))==0
table(hasannot2w)


##Keep only expressed genes with defined annotation and 
##recompute library sizes
Datta22 <- Datta22[isexpr12 & hasannot2w, keep.lib.sizes=FALSE]
dim(Datta22$counts)
head(y)

##Apply TMM normalization
Datta221 <- calcNormFactors(Datta22)
head(Datta221)


####Simulation-III###


##Form a DGEList object combining read counts and associated annotation
Datta33 <- DGEList(counts=sim.data.bulk3$counts, genes=sim.data.bulk3$rowData, group=Group)
Datta33$counts[1:5,1:12]
dim((Datta33))

## filter out the genes with low read counts.

##Determine which genes are expressed in a worthwhile number of samples.
isexpr13 <- filterByExpr(Datta33, group=Group)
table(isexpr13)

##Determine which genes have defined annotations
hasannot3w <- rowSums(is.na(Datta33$genes))==0
table(hasannot3w)


##Keep only expressed genes with defined annotation and 
##recompute library sizes
Datta33 <- Datta33[isexpr13 & hasannot3w, keep.lib.sizes=FALSE]
dim(Datta33$counts)
head(y)

##Apply TMM normalization
Datta331 <- calcNormFactors(Datta33)
head(Datta331)


####Simulation-IV###


##Form a DGEList object combining read counts and associated annotation
Datta44 <- DGEList(counts=sim.data.bulk4$counts, genes=sim.data.bulk4$rowData, group=Group)
Datta44$counts[1:5,1:12]
dim(Datta44)

## filter out the genes with low read counts.

##Determine which genes are expressed in a worthwhile number of samples.
isexpr14 <- filterByExpr(Datta44, group=Group)
table(isexpr14)

##Determine which genes have defined annotations
hasannot4w <- rowSums(is.na(Datta44$genes))==0
table(hasannot4w)


##Keep only expressed genes with defined annotation and 
##recompute library sizes
Datta44 <- Datta44[isexpr14 & hasannot4w, keep.lib.sizes=FALSE]
dim(Datta44$counts)
head(y)

##Apply TMM normalization
Datta441 <- calcNormFactors(Datta44)
head(Datta441)



####Simulation-V###


##Form a DGEList object combining read counts and associated annotation
Datta55 <- DGEList(counts=sim.data.bulk5$counts, genes=sim.data.bulk5$rowData, group=Group)
Datta55$counts[1:5,1:12]
dim(Datta55)

## filter out the genes with low read counts.

##Determine which genes are expressed in a worthwhile number of samples.
isexpr15 <- filterByExpr(Datta55, group=Group)
table(isexpr15)

##Determine which genes have defined annotations
hasannot5w <- rowSums(is.na(Datta55$genes))==0
table(hasannot5w)


##Keep only expressed genes with defined annotation and 
##recompute library sizes
Datta55 <- Datta55[isexpr15 & hasannot5w, keep.lib.sizes=FALSE]
dim(Datta55$counts)
head(y)

##Apply TMM normalization
Datta551 <- calcNormFactors(Datta55)
head(Datta551)


######Data transpose############
a <- t(Datta11$counts)
head(a[1:5,1:5])

aa <- t(Datta112$counts)
head(a[1:5,1:5])

b <- t(Datta22$counts)
head(a[1:5,1:5])

bb <- t(Datta221$counts)
head(a[1:5,1:5])

c <- t(Datta33$counts)
head(a[1:5,1:5])

cc <- t(Datta331$counts)
head(a[1:5,1:5])

d <- t(Datta44$counts)
head(a[1:5,1:5])

dd <- t(Datta441$counts)
head(a[1:5,1:5])

e <- t(Datta55$counts)
head(a[1:5,1:5])

ee <- t(Datta551$counts)
head(a[1:5,1:5])


Datta1 <- t(sim.data.bulk1$counts)
head(Datta1[1:5,1:5])

Datta2 <- t(sim.data.bulk2$counts)
head(Datta2[1:5,1:5])

Datta3 <- t(sim.data.bulk3$counts)
head(Datta3[1:5,1:5])

Datta4 <- t(sim.data.bulk4$counts)
head(Datta4[1:5,1:5])

Datta5 <- t(sim.data.bulk5$counts)
head(Datta5[1:5,1:5])

library(stats)
library(magrittr)
library(ggplot2)
###PCa Results on Simulated Data####


#Result-1



##Perform PCA and plot PC1 vs. PC2
pca_a <- prcomp(a, scale.=TRUE, center=TRUE)
summary(pca_a)

pc_evalues.a <- pca_a$sdev^2

pc_evalues.a <- tibble(PC = factor(1:length(pc_evalues.a)), 
                         variance = pc_evalues.a) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.a

#screeplot 
pc_evalues.a %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_a, addlabels = TRUE, ylim = c(0, 70))


ggplot(plot_data, aes(x = PC1,y=PC2, col=Status)) + geom_point()
# The PC scores are stored in the "x" value of the prcomp object
pc_scores.a <- pca_a$x

#Data frame
pc_scores.a <- pc_scores.a %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.a

#PC plot
pc_scores.a %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.a <- pca_a$rotation

pc_loadings.a <- pc_loadings.a %>% 
  as_tibble(rownames = "Gene")

# print the result
pc_loadings.a

library(tidyverse)

top_genes.a <- pc_loadings.a %>% 
  # select only the PCs we are interested in
  select(Gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(Gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.a

top_loadings.a <- pc_loadings.a %>% 
  filter(gene %in% top_genes.a)

loadings_plot.a <- ggplot(data = top_loadings.a) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.a

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_a$rotation[pca_a$rotation[,"PC1"]==max(pca_a$rotation[,"PC1"]), "PC1", drop=F]


##Results-2

pca_aa <- prcomp(a, scale.=FALSE, center=TRUE)
summary(pca_aa)

pc_evalues.aa <- pca_aa$sdev^2

pc_evalues.aa <- tibble(PC = factor(1:length(pc_evalues.aa)), 
                       variance = pc_evalues.aa) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.aa

#screeplot 
pc_evalues.aa %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_aa, addlabels = TRUE, ylim = c(0, 80))

# The PC scores are stored in the "x" value of the prcomp object
pc_scores.aa <- pca_aa$x

#Data frame
pc_scores.aa <- pc_scores.aa %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.aa

#PC plot
pc_scores.aa %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.aa <- pca_aa$rotation

pc_loadings.aa <- pc_loadings.aa %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.aa

library(tidyverse)

top_genes.aa <- pc_loadings.aa %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.aa

top_loadings.aa <- pc_loadings.aa %>% 
  filter(gene %in% top_genes.aa)

loadings_plot.aa <- ggplot(data = top_loadings.aa) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.aa

library(ggfortify)
autoplot(pca_aa, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_aa$rotation[pca_aa$rotation[,"PC1"]==max(pca_aa$rotation[,"PC1"]), "PC1", drop=F]

#results-3
pca_aaa <- prcomp(aa, center = FALSE, scale.=TRUE)
summary(pca_aaa)

pc_evalues.aaa <- pca_aaa$sdev^2

pc_evalues.aaa <- tibble(PC = factor(1:length(pc_evalues.aaa)), 
                       variance = pc_evalues.aaa) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.aaa

#screeplot 
pc_evalues.aaa %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_aaa, addlabels = TRUE, ylim = c(0, 80))



# The PC scores are stored in the "x" value of the prcomp object
pc_scores.aaa <- pca_aaa$x

#Data frame
pc_scores.aaa <- pc_scores.aaa %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.aaa

#PC plot
pc_scores.aaa %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.aaa <- pca_aaa$rotation

pc_loadings.aaa <- pc_loadings.aaa %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.aaa

library(tidyverse)

top_genes.aaa <- pc_loadings.aaa %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.aaa

top_loadings.aaa <- pc_loadings.aaa %>% 
  filter(gene %in% top_genes.aaa)

loadings_plot.aaa <- ggplot(data = top_loadings.aaa) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.aaa

library(ggfortify)
autoplot(pca_aaa, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_aaa$rotation[pca_aaa$rotation[,"PC1"]==max(pca_aaa$rotation[,"PC1"]), "PC1", drop=F]

#results-4

pca_aaa1 <- prcomp(aa, scale.=FALSE, center = FALSE)
summary(pca_aaa1)

pc_evalues.aaa1 <- pca_aaa1$sdev^2

pc_evalues.aaa1 <- tibble(PC = factor(1:length(pc_evalues.aaa1)), 
                       variance = pc_evalues.aaa1) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.aaa1

#screeplot 
pc_evalues.aaa1 %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_aaa1, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.aaa1 <- pca_aaa1$x

#Data frame
pc_scores.aaa1 <- pc_scores.aaa1 %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.aaa1

#PC plot
pc_scores.aaa1 %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.aaa1 <- pca_aaa1$rotation

pc_loadings.aaa1 <- pc_loadings.aaa1 %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.aaa1

library(tidyverse)

top_genes.aaa1 <- pc_loadings.aaa1 %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.aaa1

top_loadings.aaa1 <- pc_loadings.aaa1 %>% 
  filter(gene %in% top_genes.aaa1)

loadings_plot.aaa1 <- ggplot(data = top_loadings.aaa1) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.aaa1

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_aaa1$rotation[pca_aaa1$rotation[,"PC1"]==max(pca_aaa1$rotation[,"PC1"]), "PC1", drop=F]

##Simulation-2

#Reults-1

pca_b <- prcomp(b, scale.=FALSE, center=TRUE)
summary(pca_b)

pc_evalues.b <- pca_b$sdev^2

pc_evalues.b <- tibble(PC = factor(1:length(pc_evalues.b)), 
                       variance = pc_evalues.b) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.b

#screeplot 
pc_evalues.b %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_b, addlabels = TRUE, ylim = c(0, 80))

# The PC scores are stored in the "x" value of the prcomp object
pc_scores.b <- pca_b$x

#Data frame
pc_scores.b <- pc_scores.b %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.b

#PC plot
pc_scores.b %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.b <- pca_b$rotation

pc_loadings.b <- pc_loadings.b %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.b

library(tidyverse)

top_genes.b <- pc_loadings.b %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.b

top_loadings.b <- pc_loadings.b %>% 
  filter(gene %in% top_genes.b)

loadings_plot.b <- ggplot(data = top_loadings.b) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.b

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_b$rotation[pca_b$rotation[,"PC1"]==max(pca_b$rotation[,"PC1"]), "PC1", drop=F]

#Results-2

pca_bb <- prcomp(b, scale.=TRUE, center=TRUE)
summary(pca_bb)

pc_evalues.bb <- pca_bb$sdev^2

pc_evalues.bb <- tibble(PC = factor(1:length(pc_evalues.bb)), 
                       variance = pc_evalues.bb) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.bb

#screeplot 
pc_evalues.bb %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_aa, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.bb <- pca_bb$x

#Data frame
pc_scores.bb <- pc_scores.bb %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.bb

#PC plot
pc_scores.bb %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.bb <- pca_bb$rotation

pc_loadings.bb <- pc_loadings.bb %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.bb

library(tidyverse)

top_genes.bb <- pc_loadings.bb %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.bb

top_loadings.bb <- pc_loadings.bb %>% 
  filter(gene %in% top_genes.bb)

loadings_plot.bb <- ggplot(data = top_loadings.bb) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.bb

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_bb$rotation[pca_bb$rotation[,"PC1"]==max(pca_bb$rotation[,"PC1"]), "PC1", drop=F]

#Results-3

pca_bbb <- prcomp(bb, scale.=TRUE, center=FALSE)
summary(pca_bbb)

pc_evalues.bbb <- pca_bbb$sdev^2

pc_evalues.bbb <- tibble(PC = factor(1:length(pc_evalues.bbb)), 
                       variance = pc_evalues.bbb) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.bbb

#screeplot 
pc_evalues.bbb %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_bbb, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.bbb <- pca_bbb$x

#Data frame
pc_scores.bbb <- pc_scores.bbb %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.bbb

#PC plot
pc_scores.bbb %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.bbb <- pca_bbb$rotation

pc_loadings.bbb <- pc_loadings.bbb %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.bbb

library(tidyverse)

top_genes.bbb <- pc_loadings.bbb %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.bbb

top_loadings.bbb <- pc_loadings.bbb %>% 
  filter(gene %in% top_genes.bbb)

loadings_plot.bbb <- ggplot(data = top_loadings.bbb) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.bbb

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_bbb$rotation[pca_bbb$rotation[,"PC1"]==max(pca_bbb$rotation[,"PC1"]), "PC1", drop=F]


##Results-4

pca_bbbb <- prcomp(b, scale.=FALSE, center=FALSE)
summary(pca_bbbb)

pc_evalues.bbbb <- pca_bbbb$sdev^2

pc_evalues.bbbb <- tibble(PC = factor(1:length(pc_evalues.bbbb)), 
                       variance = pc_evalues.bbbb) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.bbbb

#screeplot 
pc_evalues.bbbb %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_bbbb, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.bbbb <- pca_bbbb$x

#Data frame
pc_scores.bbbb <- pc_scores.bbbb %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.bbbb

#PC plot
pc_scores.bbbb %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.bbbb <- pca_bbbb$rotation

pc_loadings.bbbb <- pc_loadings.bbbb %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.bbbb

library(tidyverse)

top_genes.bbbb <- pc_loadings.bbbb %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.bbbb

top_loadings.bbbb <- pc_loadings.bbbb %>% 
  filter(gene %in% top_genes.bbbb)

loadings_plot.bbbb <- ggplot(data = top_loadings.bbbb) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.bbbb

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_bbbb$rotation[pca_bbbb$rotation[,"PC1"]==max(pca_bbbb$rotation[,"PC1"]), "PC1", drop=F]


##Simulation-3

#Results-1

pca_c <- prcomp(c, scale.=FALSE, center=TRUE)
summary(pca_c)

pc_evalues.c <- pca_c$sdev^2

pc_evalues.c <- tibble(PC = factor(1:length(pc_evalues.c)), 
                       variance = pc_evalues.c) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.c

#screeplot 
pc_evalues.c %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_c, addlabels = TRUE, ylim = c(0, 80))



# The PC scores are stored in the "x" value of the prcomp object
pc_scores.c <- pca_c$x

#Data frame
pc_scores.c <- pc_scores.c %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.c

#PC plot
pc_scores.c %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.c <- pca_c$rotation

pc_loadings.c <- pc_loadings.c %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.c

library(tidyverse)

top_genes.c <- pc_loadings.c %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.c

top_loadings.c <- pc_loadings.c %>% 
  filter(gene %in% top_genes.c)

loadings_plot.c <- ggplot(data = top_loadings.c) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.c

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_c$rotation[pca_c$rotation[,"PC1"]==max(pca_c$rotation[,"PC1"]), "PC1", drop=F]

## Results-2

pca_cc <- prcomp(c, scale.=TRUE, center=TRUE)
summary(pca_cc)

pc_evalues.cc <- pca_cc$sdev^2

pc_evalues.cc <- tibble(PC = factor(1:length(pc_evalues.cc)), 
                       variance = pc_evalues.cc) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.cc

#screeplot 
pc_evalues.cc %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_cc, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.cc <- pca_cc$x

#Data frame
pc_scores.cc <- pc_scores.cc %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.cc

#PC plot
pc_scores.cc %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.cc <- pca_cc$rotation

pc_loadings.cc <- pc_loadings.cc %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.cc

library(tidyverse)

top_genes.cc <- pc_loadings.cc %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.cc

top_loadings.cc <- pc_loadings.%>% 
  filter(gene %in% top_genes.cc)

loadings_plot.cc <- ggplot(data = top_loadings.cc) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.cc

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_cc$rotation[pca_cc$rotation[,"PC1"]==max(pca_cc$rotation[,"PC1"]), "PC1", drop=F]

#Result-3

pca_ccc <- prcomp(cc, scale.=TRUE, center=FALSE)
summary(pca_ccc)

pc_evalues.ccc <- pca_ccc$sdev^2

pc_evalues.ccc <- tibble(PC = factor(1:length(pc_evalues.ccc)), 
                       variance = pc_evalues.ccc) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.ccc

#screeplot 
pc_evalues.ccc %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_ccc, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.ccc <- pca_ccc$x

#Data frame
pc_scores.ccc <- pc_scores.ccc %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.ccc

#PC plot
pc_scores.ccc %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.ccc <- pca_ccc$rotation

pc_loadings.ccc <- pc_loadings.ccc %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.ccc

library(tidyverse)

top_genes.ccc <- pc_loadings.ccc %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.ccc

top_loadings.ccc <- pc_loadings.ccc %>% 
  filter(gene %in% top_genes.ccc)

loadings_plot.ccc <- ggplot(data = top_loadings.ccc) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.ccc

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_ccc$rotation[pca_ccc$rotation[,"PC1"]==max(pca_ccc$rotation[,"PC1"]), "PC1", drop=F]

#Results-4

pca_cccc <- prcomp(cc, scale.=FALSE, center =FALSE)
summary(pca_cccc)

pc_evalues.cccc <- pca_cccc$sdev^2

pc_evalues.cccc <- tibble(PC = factor(1:length(pc_evalues.cccc)), 
                       variance = pc_evalues.cccc) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.cccc

#screeplot 
pc_evalues.cccc %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_cccc, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.cccc <- pca_cccc$x

#Data frame
pc_scores.cccc <- pc_scores.cccc %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.cccc

#PC plot
pc_scores.cccc %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.cccc <- pca_cccc$rotation

pc_loadings.cccc <- pc_loadings.cccc %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.cccc

library(tidyverse)

top_genes.cccc <- pc_loadings.cccc %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.cccc

top_loadings.cccc <- pc_loadings.cccc %>% 
  filter(gene %in% top_genes.cccc)

loadings_plot.cccc <- ggplot(data = top_loadings.cccc) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.cccc

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_cccc$rotation[pca_cccc$rotation[,"PC1"]==max(pca_cccc$rotation[,"PC1"]), "PC1", drop=F]



###Simulation-4

##Results-1

##Perform PCA and plot PC1 vs. PC2
pca_d <- prcomp(d, scale.=TRUE, center=TRUE)
summary(pca_d)

pc_evalues.d <- pca_d$sdev^2

pc_evalues.d <- tibble(PC = factor(1:length(pc_evalues.d)), 
                       variance = pc_evalues.d) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.d

#screeplot 
pc_evalues.d %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_d, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.d <- pca_d$x

#Data frame
pc_scores.d <- pc_scores.d %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.d

#PC plot
pc_scores.d %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.d <- pca_d$rotation

pc_loadings.d <- pc_loadings.d %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.d

library(tidyverse)

top_genes.d <- pc_loadings.d %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.d

top_loadings.d <- pc_loadings.d %>% 
  filter(gene %in% top_genes.d)

loadings_plot.d <- ggplot(data = top_loadings.d) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.d

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_d$rotation[pca_d$rotation[,"PC1"]==max(pca_d$rotation[,"PC1"]), "PC1", drop=F]


##Results-2

##Perform PCA and plot PC1 vs. PC2
pca_dd <- prcomp(d, scale.=FALSE, center=TRUE)
summary(pca_dd)

pc_evalues.dd <- pca_dd$sdev^2

pc_evalues.dd <- tibble(PC = factor(1:length(pc_evalues.dd)), 
                       variance = pc_evalues.dd) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.dd

#screeplot 
pc_evalues.dd %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_dd, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.dd <- pca_dd$x

#Data frame
pc_scores.dd <- pc_scores.dd %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.dd

#PC plot
pc_scores.dd %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.dd <- pca_dd$rotation

pc_loadings.dd <- pc_loadings.dd %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.dd

library(tidyverse)

top_genes.dd <- pc_loadings.dd %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.dd

top_loadings.dd <- pc_loadings.dd %>% 
  filter(gene %in% top_genes.dd)

loadings_plot.dd <- ggplot(data = top_loadings.dd) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.dd

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_dd$rotation[pca_dd$rotation[,"PC1"]==max(pca_dd$rotation[,"PC1"]), "PC1", drop=F]


#Results-3

##Perform PCA and plot PC1 vs. PC2
pca_ddd <- prcomp(d, scale.=TRUE, center=FALSE)
summary(pca_ddd)

pc_evalues.ddd <- pca_ddd$sdev^2

pc_evalues.ddd <- tibble(PC = factor(1:length(pc_evalues.ddd)), 
                       variance = pc_evalues.ddd) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.ddd

#screeplot 
pc_evalues.ddd %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_ddd, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.ddd <- pca_ddd$x

#Data frame
pc_scores.ddd <- pc_scores.ddd %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.ddd

#PC plot
pc_scores.ddd %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.ddd <- pca_ddd$rotation

pc_loadings.ddd <- pc_loadings.ddd %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.ddd

library(tidyverse)

top_genes.ddd <- pc_loadings.ddd %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.ddd

top_loadings.ddd <- pc_loadings.ddd %>% 
  filter(gene %in% top_genes.ddd)

loadings_plot.ddd <- ggplot(data = top_loadings.ddd) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.ddd

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_ddd$rotation[pca_ddd$rotation[,"PC1"]==max(pca_ddd$rotation[,"PC1"]), "PC1", drop=F]


##Results-5

##Perform PCA and plot PC1 vs. PC2
pca_dddd <- prcomp(d, scale.=FALSE, center=FALSE)
summary(pca_dddd)

pc_evalues.dddd <- pca_dddd$sdev^2

pc_evalues.dddd <- tibble(PC = factor(1:length(pc_evalues.dddd)), 
                       variance = pc_evalues.dddd) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.dddd

#screeplot 
pc_evalues.dddd %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_dddd, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.dddd <- pca_dddd$x

#Data frame
pc_scores.dddd <- pc_scores.dddd %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.dddd

#PC plot
pc_scores.dddd %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.dddd <- pca_dddd$rotation

pc_loadings.dddd <- pc_loadings.dddd %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.dddd

library(tidyverse)

top_genes.dddd <- pc_loadings.dddd %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.dddd

top_loadings.dddd <- pc_loadings.dddd %>% 
  filter(gene %in% top_genes.dddd)

loadings_plot.dddd <- ggplot(data = top_loadings.dddd) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.dddd

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_dddd$rotation[pca_dddd$rotation[,"PC1"]==max(pca_dddd$rotation[,"PC1"]), "PC1", drop=F]


#Simulation-5

#Results-1

##Perform PCA and plot PC1 vs. PC2
pca_e <- prcomp(e, scale.=TRUE, center=TRUE)
summary(pca_e)

pc_evalues.e <- pca_e$sdev^2

pc_evalues.e <- tibble(PC = factor(1:length(pc_evalues.e)), 
                       variance = pc_evalues.e) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.e

#screeplot 
pc_evalues.e %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_e, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.e <- pca_e$x

#Data frame
pc_scores.e <- pc_scores.e %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.e

#PC plot
pc_scores.e %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.e <- pca_e$rotation

pc_loadings.e <- pc_loadings.e %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.e

library(tidyverse)

top_genes.e <- pc_loadings.e %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.e

top_loadings.e <- pc_loadings.e %>% 
  filter(gene %in% top_genes.e)

loadings_plot.e <- ggplot(data = top_loadings.e) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.e

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_e$rotation[pca_e$rotation[,"PC1"]==max(pca_e$rotation[,"PC1"]), "PC1", drop=F]

#Results-2

##Perform PCA and plot PC1 vs. PC2
pca_ee <- prcomp(e, scale.=FALSE, center=TRUE)
summary(pca_ee)

pc_evalues.ee <- pca_ee$sdev^2

pc_evalues.ee <- tibble(PC = factor(1:length(pc_evalues.ee)), 
                       variance = pc_evalues.ee) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.ee

#screeplot 
pc_evalues.ee %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_ee, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.ee <- pca_ee$x

#Data frame
pc_scores.ee <- pc_scores.ee %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.ee

#PC plot
pc_scores.ee %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.ee <- pca_ee$rotation

pc_loadings.ee <- pc_loadings.ee %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.ee

library(tidyverse)

top_genes.ee <- pc_loadings.ee %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.ee

top_loadings.ee <- pc_loadings.ee %>% 
  filter(gene %in% top_genes.ee)

loadings_plot.ee <- ggplot(data = top_loadings.ee) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.ee

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_ee$rotation[pca_ee$rotation[,"PC1"]==max(pca_ee$rotation[,"PC1"]), "PC1", drop=F]


#Results-3

##Perform PCA and plot PC1 vs. PC2
pca_eee <- prcomp(e, scale.=TRUE, center=FALSE)
summary(pca_eee)

pc_evalues.eee <- pca_eee$sdev^2

pc_evalues.eee <- tibble(PC = factor(1:length(pc_evalues.eee)), 
                       variance = pc_evalues.eee) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.eee

#screeplot 
pc_evalues.eee %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_eee, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.eee <- pca_eee$x

#Data frame
pc_scores.eee <- pc_scores.eee %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.eee

#PC plot
pc_scores.eee %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.eee <- pca_eee$rotation

pc_loadings.eee <- pc_loadings.eee %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.eee

library(tidyverse)

top_genes.eee <- pc_loadings.eee %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.eee

top_loadings.eee <- pc_loadings.eee %>% 
  filter(gene %in% top_genes.eee)

loadings_plot.eee <- ggplot(data = top_loadings.eee) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.eee

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_eee$rotation[pca_eee$rotation[,"PC1"]==max(pca_eee$rotation[,"PC1"]), "PC1", drop=F]

##Results-5

##Perform PCA and plot PC1 vs. PC2
pca_eeee <- prcomp(e, scale.=FALSE, center=FALSE)
summary(pca_e)

pc_evalues.eeee <- pca_eeee$sdev^2

pc_evalues.eeee <- tibble(PC = factor(1:length(pc_evalues.eeee)), 
                       variance = pc_evalues.eeee) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_evalues.eeee

#screeplot 
pc_evalues.eeee %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

#visulaizing samples on PC space
fviz_eig(pca_eee, addlabels = TRUE, ylim = c(0, 80))


# The PC scores are stored in the "x" value of the prcomp object
pc_scores.eeee <- pca_eeee$x

#Data frame
pc_scores.eeee <- pc_scores.eeee %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")


# print the result
pc_scores.eeee

#PC plot
pc_scores.eeee %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Exploring correlation between genes and PCs

pc_loadings.eeee <- pca_eeee$rotation

pc_loadings.eeee <- pc_loadings.eeee %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings.eeee

library(tidyverse)

top_genes.eeee <- pc_loadings.eeee %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

top_genes.eeee

top_loadings.eeee <- pc_loadings.eeee %>% 
  filter(gene %in% top_genes.eeee)

loadings_plot.eeee <- ggplot(data = top_loadings.eeee) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot.eeee

library(ggfortify)
autoplot(pca_a, data = res, colour = "sample")#, shape = "sample")
plotPCA(res,intgroup=c("Group"))

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_eeee$rotation[pca_eeee$rotation[,"PC1"]==max(pca_eeee$rotation[,"PC1"]), "PC1", drop=F]






















########################################################

pca_sim1 <- prcomp(Datta1, scale. =TRUE, center = FALSE)
summary(pca_sim1)


pca_sim11 <- prcomp(Datta1, center = TRUE, scale. = TRUE)
summary(pca_sim11)


pca_sim2 <- prcomp(Datta2, scale = TRUE)
summary(pca_sim2)

pca_sim22 <- prcomp(Datta2, center = TRUE, scale. = FALSE)
summary(pca_sim22)

pca_sim3 <- prcomp(Datta3, scale=TRUE)
summary(pca_sim3)

pca_sim33 <- prcomp(Datta3, center = TRUE, scale. = FALSE)
summary(pca_sim33)

pca_sim4 <- prcomp(Datta4, scale=TRUE)
summary(pca_sim4)

pca_sim44 <- prcomp(Datta4, center = TRUE, scale. = FALSE)
summary(pca_sim44)

pca_sim5 <- prcomp(Datta5, scale=TRUE)
summary(pca_sim5)

pca_sim55 <- prcomp(Datta5, center = TRUE, scale. = FALSE)
summary(pca_sim55)


## pca.res is a list with 5 elements 
length(pca_sim1)
length(pca_sim2)
length(pca_sim3)
length(pca_sim4)
length(pca_sim5)


names(pca_sim1)
names(pca_sim2)
names(pca_sim3)
names(pca_sim4)
names(pca_sim5)

str(pca_sim1)





# The rotation element contains the loadings/weights/phi for each gene in each component
pca_sim1$rotation[1:5, 1:5]
dim(pca_sim1$rotation)


pca_sim2$rotation[1:5, 1:5]
dim(pca_sim2$rotation)


pca_sim3$rotation[1:5, 1:5]
dim(pca_sim3$rotation)


pca_sim4$rotation[1:5, 1:5]
dim(pca_sim4$rotation)


pca_sim5$rotation[1:5, 1:5]
dim(pca_sim5$rotation)





# The sum of the squared loadings is equal to 1 for each PC
sum(pca_sim4$rotation[,"PC1"]^2)
sum(pca_sim4$rotation[,"PC4"]^2)
sum(pca_sim4$rotation[,"PC12"]^2)

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_sim1$rotation[pca_sim1$rotation[,"PC1"]==max(pca_sim1$rotation[,"PC1"]), "PC1", drop=F]

pca_sim2$rotation[pca_sim2$rotation[,"PC1"]==max(pca_sim2$rotation[,"PC1"]), "PC1", drop=F]

pca_sim3$rotation[pca_sim3$rotation[,"PC1"]==max(pca_sim3$rotation[,"PC1"]), "PC1", drop=F]

pca_sim4$rotation[pca_sim4$rotation[,"PC1"]==max(pca_sim4$rotation[,"PC1"]), "PC1", drop=F]

pca_sim5$rotation[pca_sim4$rotation[,"PC1"]==max(pca_sim5$rotation[,"PC1"]), "PC1", drop=F]




pc_eigen1 <- pca_sim1$sdev^2
pc_eigen2 <- pca_sim2$sdev^2
pc_eigen3 <- pca_sim3$sdev^2
pc_eigen4 <- pca_sim4$sdev^2
pc_eigen5 <- pca_sim5$sdev^2


pc_eigen11 <- pca_a$sdev^2
pc_eigen22 <- pca_b$sdev^2
pc_eigen33 <- pca_c$sdev^2
pc_eigen44 <- pca_d$sdev^2
pc_eigen55 <- pca_e$sdev^2

length(pc_eigen1)

library(tidyverse)

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigen1 <- tibble(PC = factor(1:length(pc_eigen1)), 
                         variance = pc_eigen1) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigen1

pc_eigen1 %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")


## The summary() function will give us statistical information about each PC
summary(pca_sim1)$importance[,c("PC1","PC2","PC3","PC4","PC5")]
summary(pca_sim2)$importance[,c("PC1","PC2","PC3","PC4","PC5")]
summary(pca_sim3)$importance[,c("PC1","PC2","PC3","PC4","PC5")]
summary(pca_sim4)$importance[,c("PC1","PC2","PC3","PC4","PC5")]
summary(pca_sim5)$importance[,c("PC1","PC2","PC3","PC4","PC5")]

## NOTE: 25% of the variation in the data is explained by PC1
## the remaining components explain a lot less, lets keep this information in an object
varex=round(100*summary(pca_sim1)$importance[2,c("PC1","PC2")])

## The element named x, this tells us the 'projected' values for each sample, 
# eg what is the expression of each sample when projected into PC1
pca_sim1$x[1:5,1:5]
pca_sim2$x[1:5,1:5]
pca_sim3$x[1:5,1:5]
pca_sim4$x[1:5,1:5]
pca_sim5$x[1:5,1:5]


# notice the variance decreases with each PC, hence the data projected onto PC1 has highest variance
var(pca_sim1$x[,"PC1"]) # >3,000
var(pca_sim1$x[,"PC2"]) # ~900
var(pca_sim1$x[,"PC10"]) # ~200

# The looking at the max() weight/loading on PC1 can tell us which gene contributes most to the variation.
pca_sim1$rotation[pca_sim1$rotation[,"PC1"]==max(pca_sim1$rotation[,"PC1"]), "PC1", drop=F]


## Lets extract the values for each samples projected on PC1 and PC2
db.proj.sim1<-pca_sim1$x[,c("PC1","PC2")]
## and then add our phenotype data to db.proj
db.proj.sim1 <- as.data.frame(db.proj.sim1)
# make sure the rownames are correct order
head(db.proj.sim1)

## Now lets plot the PC data, without any phenotype information
## The x-axis will be PC1 and the y-axis will be PC2
library(ggplot2) # 
### Plot 1: Baseline By tissue
ggplot(db.proj.sim1, aes(x=PC1, y=PC2)) +
  geom_point(size=2) +
  labs(x=paste('PC-1 (',varex[1],'%)',sep=''),
       y=paste('PC-2 (',varex[2],'%)',sep=''), 
       title="") + 
  theme_bw()



ggsave(file='Plot1BLByTissuie.pdf', height=4, width=5)

fviz_eig(pca_sim1, addlabels = TRUE, ylim = c(0, 70))



#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
pca.propvar1 <- ((pca_sim1$sdev^2) / (sum(pca_sim1$sdev^2)))*100
pca.propvar2 <- ((pca_sim2$sdev^2) / (sum(pca_sim2$sdev^2)))*100
pca.propvar3 <- ((pca_sim3$sdev^2) / (sum(pca_sim3$sdev^2)))*100
pca.propvar4 <- ((pca_sim4$sdev^2) / (sum(pca_sim4$sdev^2)))*100
pca.propvar5 <- ((pca_sim5$sdev^2) / (sum(pca_sim5$sdev^2)))*100

##Screeplot####
q1 <- barplot(pca.propvar1, cex.names=1,
        xlab=paste("Principal component (PC), 1-",
                   length(pca_sim1$sdev)), ylab="Proportion of variation (%)", 
        main="Scree plot", ylim=c(0,100))

q2 <- barplot(pca.propvar2, cex.names=1,
        xlab=paste("Principal component (PC), 1-",
                   length(pca_sim2$sdev)), ylab="Proportion of variation (%)", 
        main="Scree plot", ylim=c(0,100))

q3 <- barplot(pca.propvar3, cex.names=1,
        xlab=paste("Principal component (PC), 1-",
                   length(pca_sim3$sdev)), ylab="Proportion of variation (%)", 
        main="Scree plot", ylim=c(0,100))

q4 <- barplot(pca.propvar4, cex.names=1,
        xlab=paste("Principal component (PC), 1-",
                   length(pca_sim4$sdev)), ylab="Proportion of variation (%)", 
        main="Scree plot", ylim=c(0,100))

q5 <- barplot(pca.propvar5, cex.names=1,
        xlab=paste("Principal component (PC), 1-",
                   length(pca_sim5$sdev)), ylab="Proportion of variation (%)", 
        main="Scree plot", ylim=c(0,100))



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(q1, q2, q3, q4, q5)


#####Pairs plots#####

par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(pca_sim1$x[,1:6], col="orange", main="Principal components analysis bi-plot1\nPCs 1-6", pch=16)
pairs(pca_sim1$x[,7:12], col="orange", main="Principal components analysis bi-plot1\nPCs 7-12", pch=16)


pairs(pca_sim2$x[,1:6], col="green", main="Principal components analysis bi-plo2\nPCs 1-6", pch=16)
pairs(pca_sim2$x[,7:12], col="green", main="Principal components analysis bi-plot2\nPCs 7-12", pch=16)


pairs(pca_sim3$x[,1:6], col="blue", main="Principal components analysis bi-plo3\nPCs 1-6", pch=16)
pairs(pca_sim3$x[,7:12], col="blue", main="Principal components analysis bi-plot3\nPCs 7-12", pch=16)


pairs(pca_sim4$x[,1:6], col="gray", main="Principal components analysis bi-plot4\nPCs 1-6", pch=16)
pairs(pca_sim4$x[,7:12], col="gray", main="Principal components analysis bi-plot4\nPCs 7-12", pch=16)


pairs(pca_sim5$x[,1:6], col="red", main="Principal components analysis bi-plot5\nPCs 1-6", pch=16)
pairs(pca_sim5$x[,7:12], col="red", main="Principal components analysis bi-plot5\nPCs 7-12", pch=16)


##Bi-plots

par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)

#Plots scatter plot for PC 1 and 2
plot(pca_sim1$x, type="n", main="Principal components analysis bi-plot1", xlab=paste("PC1, ", round(pca.propvar1[1], 2), "%"), ylab=paste("PC2, ", round(pca.propvar1[2], 2), "%"))
points(pca_sim1$x, col="red", pch=16, cex=1)


#Plots scatter plot for PC 1 and 3
plot(pca_sim1$x[,1], pca_sim1$x[,3], type="n", main="Principal components analysis bi-plot1", xlab=paste("PC1, ", round(pca.propvar1[1], 2), "%"), ylab=paste("PC3, ", round(pca.propvar1[3], 2), "%"))
points(pca_sim1$x[,1], pca_sim1$x[,3], col="blue", pch=16, cex=1)

#Plots scatter plot for PC 2 and 3
plot(pca_sim1$x[,2], pca_sim1$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC2, ", round(pca.propvar1[2], 2), "%"), ylab=paste("PC3, ", round(pca.propvar1[3], 2), "%"))
points(pca_sim1$x[,2], pca_sim1$x[,3], col="green", pch=16, cex=1)


#####Tri-plot

require(scatterplot3d)
par(mar=c(4,4,4,4), cex=1.0, cex.main=0.8, cex.axis=0.8)

scatterplot3d(pca_sim1$x[,1:3], angle=-40, main="", color="blue", pch=17, xlab=paste("PC1, ", round(pca.propvar1[1], 2), "%"), ylab=paste("PC2, ", round(pca.propvar1[2], 2), "%"), zlab=paste("PC3, ", round(pca.propvar1[3], 2), "%"), grid=FALSE, box=FALSE)
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(pca_sim1$x[,1:3], grid = c("xy", "xz", "yz"))

source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(pca_sim1$x[,1:3], grid = c("xy", "xz", "yz"))




library(BiocGenerics)
library(DESeq2)
plotPCA(Datta1) #, intgroup = "condition", ntop = 500, returnData = FALSE)


de.genes <- genes.samp$DE.ind   # DE genes sampled from source matrix
table(de.genes)
length(de.genes)

###############################################################################
simy <- DGEList(counts=counts.SPsimseq, genes=sim.data.bulk1$rowData, group=samp.col$Group)
simy <- calcNormFactors(simy)
simy <- estimateDisp(simy, robust=TRUE)

##Determine which genes are expressed in a worthwhile number of samples.
isexprr <- filterByExpr(counts.SPsimseq, group=samp.col$Group)
table(isexpr)

##Determine which genes have defined annotations
hasannot11 <- rowSums(is.na(simy$genes))==0
table(hasannot11)

summary(simy)
##Get rid of genes with no defined annotations
counts.SPsimseq <- counts.SPsimseq[hasannot11,]
dim(counts.SPsimseq)

##Keep only expressed genes with defined annotation and 
##recompute library sizes
simy <- simy[isexprr & hasannot11, keep.lib.sizes=FALSE]
dim(simy$counts)
head(simy$counts)
names(simy)
head(simy)
dim(simy$counts)

################################################################################
et_sim <- exactTest(simy)

names(et_sim)
head(et_sim$table)
topTags(et_sim, n=10)$table

#p-values
p.sim <- et_sim $table$PValue
head(p.sim)
hist(p.sim)



ee.genes <- genes.samp[ ! (genes.samp %in% de.genes) ] # EE genes sampled from source matrix
length(ee.genes)

de.genes.sim <- counts.SPsimseq$genes.subset %in% de.genes # logical vector giving which genes are DE in simulted matrix

###############################################################################

##evaluate Densities
geneDens = evaluateDensities(sim.data.bulk1, newData = rownames(y$counts)[1])

#This returns for every sample, the midpoints (mids) and associated densities (gy)







## ----comparison, warning=FALSE, fig.width=8, fig.height=4---------------------

# compare the distributions of the mean expressions, variability, 
# and fraction of zero counts per gene
library(LSD) # for generating heatmap plots

# normalize counts for comparison  
Y0E.log.cpm <- log2(edgeR::cpm(y$counts)+1)
Y1E.log.cpm <- log2(edgeR::cpm(sim.data.bulk4$counts)+1)
summary(Y0E.log.cpm)
summary(Y1E.log.cpm)


Y0L.log.cpm <- (limma::voom(y$counts))
Y1L.log.cpm <- (limma::voom(sim.data.bulk4$counts))
summary(Y0L.log.cpm)
summary(Y1L.log.cpm)


Y0E.log.cpm <- Y0E.log.cpm[rowMeans(Y0E.log.cpm>0)>=0.1, ]
Y1E.log.cpm <- Y1E.log.cpm[rowMeans(Y1E.log.cpm>0)>=0.1, ]

Y0L.log.cpm <- Y0L.log.cpm$E[rowMeans(Y0L.log.cpm$E>0)>=0.1, ]
Y1L.log.cpm <- Y1L.log.cpm$E[rowMeans(Y1L.log.cpm$E>0)>=0.1, ]


rowVars <- function(X){apply(X, 1, var, na.rm=TRUE)}
rowCVs <- function(X){apply(X, 1, function(x) sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))}
par(mfrow=c(1, 3))

boxplot(list(real.data=log(colSums(y$counts)), 
             simulated.data=log(sim.data.bulk4$colData$sim.Lib.Size)), 
        main="library size") 

boxplot(list(real.data=rowMeans(Y0E.log.cpm), 
             simulated.data=rowMeans(Y1E.log.cpm)), 
        main="mean expression of genes(EdgeR)") 

boxplot(list(real.data=rowMeans(Y0L.log.cpm), 
             simulated.data=rowMeans(Y1L.log.cpm)), 
        main="mean expression of genes(limma-voom)") 


boxplot(list(real.data=rowVars(Y0E.log.cpm), 
             simulated.data=rowVars(Y1E.log.cpm)), 
        main="variance of gene expressions(EdgeR)") 


boxplot(list(real.data=rowVars(Y0L.log.cpm), 
             simulated.data=rowVars(Y1L.log.cpm)), 
        main="variance of gene expressions(limma-voom)") 


# compare the relationship between the mean and variability
par(mfrow=c(1,3), mar=c(4,4,4,1))
heatscatter(rowMeans(Y0L.log.cpm), rowCVs(Y0L.log.cpm), ylim=c(0, 5), xlim=c(0, 16),
            colpal="bl2gr2rd", main="real data", xlab="mean log2-CPM", 
            ylab="coefficients of variation", cexplot=0.5, alpha = 60, cex.lab=1.25)

heatscatter(rowMeans(Y1L.log.cpm), rowCVs(Y1L.log.cpm), ylim=c(0, 6), xlim=c(0, 16),
            main="SPsimSeq", xlab="mean log2-CPM", ylab="coefficients of variation", 
            cexplot=0.5, alpha = 60, colpal="bl2gr2rd", cex.lab=1.25)

n.gride <- 1000
min.g   <- seq(0, 20, length.out = n.gride+1)[-n.gride]
max.g   <- seq(0, 20, length.out = n.gride+1)[-1] 
mid.g   <- (min.g+max.g)/2
f.real  <- vapply(seq_len(n.gride), FUN.VALUE = double(1), function(r){
  x <- Y0L.log.cpm[rowMeans(Y0L.log.cpm)<=max.g[r] & rowMeans(Y0L.log.cpm)>min.g[r],]
  y <- ifelse(!is.null(dim(x)), mean(rowCVs(x)), mean(sd(x)/mean(x))) 
  y
})
f.SPsim <- vapply(seq_len(n.gride), FUN.VALUE = double(1), function(r){
  x <- Y1L.log.cpm[rowMeans(Y1L.log.cpm)<=max.g[r] & rowMeans(Y1L.log.cpm)>min.g[r],]
  y <- ifelse(!is.null(dim(x)), mean(rowCVs(x)), mean(sd(x)/mean(x))) 
  y
})

sm1 <- loess(I(f.SPsim-f.real)~mid.g) 
plot(mid.g, f.SPsim-f.real, xlim=c(0, 14), col="lightskyblue", pch=20, cex.lab=1.25,
     cex.main=1.4, main="SPsimSeq - real data", ylab="difference", xlab="mean log2-CPM")
lines(mid.g,predict(sm1, newdata = mid.g), col="blue", lwd=3) 



# compare the correlation between genes and samples 
cor.mat.Y0E <- cor(t(Y0E.log.cpm))
cor.mat.Y1E <- cor(t(Y1E.log.cpm)) 
cor.vec.Y0E <- cor.mat.Y0E[upper.tri(cor.mat.Y0E)]
cor.vec.Y1E <- cor.mat.Y1E[upper.tri(cor.mat.Y1E)] 
par(mfrow=c(1,3), mar=c(4,4,3.5,1))
hist(cor.vec.Y0E, nclass = 30, probability = TRUE, 
     border="gray", col="steelblue1", main="real data", xlab="Genewise correlations", 
     ylim=c(0, 1.5), xlim=c(-1.5, 1.5), cex.lab=1.25)
hist(cor.vec.Y1E, nclass = 30, probability = TRUE, border="gray",
     col="steelblue1",  main="SPsimSeq", xlab="Genewise correlations",
     ylim=c(0, 1.5), xlim=c(-1.5, 1.5), cex.lab=1.25)
plot(seq(-1, 1, 0.1), seq(-1, 1, 0.1), type="n", xlab="quantile (real data)", 
     ylab="quantile (simulated data)",  main="correlation quantile-quantile plot")
abline(0, 1, col="gray")
points(quantile(cor.vec.Y0E, seq(0, 1, 0.001)), quantile(cor.vec.Y1E, seq(0, 1, 0.001)), 
       col="blue", pch=20, cex=1.5, cex.lab=1.25)  



# compare the correlation between genes and samples 
cor.mat.Y0L <- cor(t(Y0L.log.cpm))
cor.mat.Y1L <- cor(t(Y1L.log.cpm)) 
cor.vec.Y0L <- cor.mat.Y0L[upper.tri(cor.mat.Y0L)]
cor.vec.Y1L <- cor.mat.Y1L[upper.tri(cor.mat.Y1L)] 
par(mfrow=c(1,3), mar=c(4,4,3.5,1))
hist(cor.vec.Y0L, nclass = 30, probability = TRUE, 
     border="gray", col="steelblue1", main="real data", xlab="Genewise correlations", 
     ylim=c(0, 1.5), xlim=c(-1.5, 1.5), cex.lab=1.25)
hist(cor.vec.Y1L, nclass = 30, probability = TRUE, border="gray",
     col="steelblue1",  main="SPsimSeq", xlab="Genewise correlations",
     ylim=c(0, 1.5), xlim=c(-1.5, 1.5), cex.lab=1.25)
plot(seq(-1, 1, 0.1), seq(-1, 1, 0.1), type="n", xlab="quantile (real data)", 
     ylab="quantile (simulated data)",  main="correlation quantile-quantile plot")
abline(0, 1, col="gray")
points(quantile(cor.vec.Y0L, seq(0, 1, 0.001)), quantile(cor.vec.Y1L, seq(0, 1, 0.001)), 
       col="blue", pch=20, cex=1.5, cex.lab=1.25) 