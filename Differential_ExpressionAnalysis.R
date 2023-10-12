library(tidyverse)
library(DESeq2)
library(gPCA)
library(EnhancedVolcano)

RawCount_CRA <- read.delim("mir_raw.txt")

Pheno_CRA <- read.csv("pheno.csv")

xvar=c("S_SUBJECTID","age","gender","pctpred_fev1_pre_BD","Inhaled_Steroids","bdrbase")
Pheno_CRA2=Pheno_CRA[xvar]

Pheno_CRA2=Pheno_CRA2 %>% drop_na(bdrbase)
quantile(Pheno_CRA2$bdrbase)

Pheno_CRA2$Q[Pheno_CRA2$bdrbase<0.2526283]="Q1"
Pheno_CRA2$Q[Pheno_CRA2$bdrbase>8.7242359]="Q3"
Pheno_CRA3=Pheno_CRA2 %>% drop_na(Q)

Pheno_CRA3=Pheno_CRA3 %>% filter(
   Pheno_CRA3$S_SUBJECTID %in% RawCount_CRA$S_SUBJECTID)

RawCount_CRA3=RawCount_CRA %>% filter(
    RawCount_CRA$S_SUBJECTID %in% Pheno_CRA3$S_SUBJECTID)

Pheno_CRA3=Pheno_CRA3[order(Pheno_CRA3$S_SUBJECTID),]
RawCount_CRA3=RawCount_CRA3[order(RawCount_CRA3$S_SUBJECTID),]
identical(RawCount_CRA3$S_SUBJECTID,Pheno_CRA3$S_SUBJECTID)
#### Testing for batch effect
#Batch_CRA <- read.delim("Batch.txt")
#xvar2=c("Norgen_Batch","S_SUBJECTID")
 #Batch_CRA=Batch_CRA %>% filter(
    #Batch_CRA$Norgen_Batch !="exosomal")
 #Batch_CRA=Batch_CRA %>% filter(
    #Batch_CRA$S_SUBJECTID %in% Pheno_CRA3$S_SUBJECTID)
#Batch_CRA=Batch_CRA[order(Batch_CRA$S_SUBJECTID),]
#identical(Batch_CRA$S_SUBJECTID,RawCount_CRA3$S_SUBJECTID)

#out_cra<-gPCA.batchdetect(x=as.matrix((RawCount_CRA3[6:322])),batch=as.factor(Batch_CRA$Norgen_Batch),center=FALSE,nperm=1000)
#out_cra$delta ; out_cra$p.val
#gDist(out_cra)
Pheno_CRA3$gender=as.factor(Pheno_CRA3$gender)
Pheno_CRA3$Q=as.factor(Pheno_CRA3$Q)
Pheno_CRA3$Inhaled_Steroids=as.factor(Pheno_CRA3$Inhaled_Steroids)
dds2 <- DESeqDataSetFromMatrix(countData=round(t(RawCount_CRA3[6:322])), colData =Pheno_CRA3,design = ~Q+ gender+age+pctpred_fev1_pre_BD+Inhaled_Steroids)
DeSeq2_CRA=DESeq(dds2)
DeSeq2_Q_CRA=results(DeSeq2_CRA,name="Q_Q3_vs_Q1",pAdjustMethod = "fdr")
EnhancedVolcano(DeSeq2_Q_CRA,
                lab = rownames(DeSeq2_Q_CRA),
                x = 'log2FoldChange',
                y = 'padj',
                ylim=c(0,3),ylab=bquote(~Log[10]~ 'P (FDR)'),
                pCutoff = 0.1,
                FCcutoff = 0.2,
                pointSize = 3.0,
                labSize = 4.0,labFace = 'bold', drawConnectors = TRUE,
                widthConnectors = 0.75)
write.csv(DeSeq2_Q_CRA,"Result_CRA_Q3vQ1_age_gender_race_ppfev1.csv")
