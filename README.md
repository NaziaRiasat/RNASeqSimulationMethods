# RNASeqSimulationMethods

Semi-parametric simulation:

The SPsimSeq simulation procedure involves two major steps: (1) estimate the probability distributions for the observed gene expression levels in a real RNA-seq dataset (bulk or single-cell), and (2) simulate new gene expression levels from the estimated distributions.

We use SPsimSeq to simulate bulk RNA-seq data used in Shahraz et al. (2020). Briefly, 12 NOX-2 knock-out (KO) mice were used in this RNA-Seq experiment. The 12 mice in the genotype (KO) were randomly assigned to one of two treatment groups (PBS control or LPS challenge), with six mice in each group. 

Simulate bulk RNA-seq data with the following features:

- 12600 genes (n.genes = 12600)
- A total number of    samples (tot. samples =172)
- The samples are divided into two groups (group. config = c(0.47, 0.53)) with 10% of the genes are DE between the groups (pDE = 0.1). The group composition is similar to that of the source data.
- The DE genes have an LFC at least 0.5 (lfc.thrld=0.5) with the t-statistic threshold of 2.5 (t.thrld=2.5) and ll. threshold of 5 (ll.thrld=5)
- All the samples are generated in a single batch (batch. config = 1), similar to the source data
- Since zero inflation is not an issue in bulk RNA-seq dataset, we do not model the zeros separately (model.zero.prob = FALSE)
- The number of classes to construct the distributions of the gene expression levels is 50% of the sample size (n) (i.e. w=0.5).



