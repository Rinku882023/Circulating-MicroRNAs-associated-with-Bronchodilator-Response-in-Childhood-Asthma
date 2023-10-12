### miRNA target gene Identification

library(multiMiR)
multimir_results <- get_multimir(org     = 'hsa',
                                 mirna   = c("hsa-miR-1246","hsa-miR-200b-3p"),
                                 table   = 'validated',
                                 summary = TRUE)
Targets=multimir_results@data
write.csv(Targets,"Table_RepmiR_validatedTarget.csv")

### Enrichment Analysis

AllTargetGene <- read.csv("AllTargetGene.txt", sep="")
target200b_3p <- read.csv("target200b_3p.txt", sep="")
target1246 <- read.csv("target1246.txt", sep="")

library(clusterProfiler)
library(ReactomePA)
library(enrichplot)

Reactome <- enrichPathway(gene=AllTargetGene$Total_validatedTarget, pvalueCutoff = 0.05, readable=TRUE)
edox2 <- pairwise_termsim(Reactome)
p1 <- treeplot(edox2,hclust_method = "ward.D2")
p1
Reactome1246 <- enrichPathway(gene=target1246$target_entrez_1246, pvalueCutoff = 0.05, readable=TRUE)
edox1246 <- pairwise_termsim(Reactome1246)
p2 <- treeplot(edox1246,hclust_method = "ward.D2",nCluster=1)
p2
Reactome200b <- enrichPathway(gene=target200b_3p$target_entrez_200b_3p, pvalueCutoff = 0.05, readable=TRUE)
edox200b <- pairwise_termsim(Reactome200b)
p3<- treeplot(edox200b,hclust_method = "ward.D2")
p3
