library(tidyverse)
# library(readr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ReactomePA)
library(annotate)
library(seqinr)

bpd_aug_table = read.table("bpd_pheno_refined.txt")
bpd_aug = as.character(bpd_aug_table$V1)

mdd_aug_table = read.table("mdd_pheno_refined.txt")
mdd_aug = as.character(mdd_aug_table$V1)

scz_aug_table = read.table('scz_pheno_refined.txt')
scz_aug = as.character(scz_aug_table$V1)

bpd_tbl = read.table("diseaseSignatures/Bipolar_disorders_and_related_disorders.txt")
bipolar = as.character(bpd_tbl$V1)

mdd_tbl = read.table("diseaseSignatures/Depressive_disorders.txt")
depress = as.character(mdd_tbl$V1)

scz_tbl = read.table("diseaseSignatures/Schizophrenia_spectrum_and_other_psychotic_disorders.txt")
scz = as.character(scz_tbl$V1)

alcohol_tbl = read.table("diseaseSignatures/Alcohol_use_disorders.txt")
alcohol = as.character(alcohol_tbl$V1)

bipolar_tbl = read.table("diseaseSignatures/Bipolar_disorders_and_related_disorders.txt")
bipolar = as.character(bipolar_tbl$V1)

cannabis_tbl = read.table("diseaseSignatures/Cannabis_use_disorders.txt")
cannabis = as.character(cannabis_tbl$V1)

cocaine_tbl = read.table("diseaseSignatures/Cocaine_use_disorders.txt")
cocaine = as.character(cocaine_tbl$V1)

depress_tbl = read.table("diseaseSignatures/Depressive_disorders.txt")
depress = as.character(depress_tbl$V1)

drug_induced_psychosis_tbl = read.table("diseaseSignatures/Drug-induced_psychosis.txt")
drug_induced_psychosis = as.character(drug_induced_psychosis_tbl$V1)

schizo_tbl = read.table("diseaseSignatures/Schizophrenia_spectrum_and_other_psychotic_disorders.txt")
schizo = as.character(schizo_tbl$V1)

substance_induced_depress_tbl =
  read.table("diseaseSignatures/Substance_drug_induced_depressive_disorder.txt")
substance_induced_depress = as.character(substance_induced_depress_tbl$V1)

all_genes = c(alcohol, bipolar, cannabis, cocaine, depress, drug_induced_psychosis, schizo, substance_induced_depress)

all_diseases = list('BPD' = bipolar, 'BPD+' = bpd_aug, 'MDD' = depress, 
                    'MDD+' = mdd_aug, 'SCZ' = scz, 'SCZ+' = scz_aug)

disease_reactome <- compareCluster(geneCluster = all_diseases, pAdjustMethod = "BH", pvalueCutoff = 0.01, fun = "enrichPathway", readable = TRUE)
head(as.data.frame(disease_reactome))
dotplot(disease_reactome)  + theme_grey(base_size = 20) + theme_grey(base_size = 22)
#write.csv(disease_reactome, file='diseaseSignatures/Figures/disease_reactome.csv')
ggsave(paste0('disease_reactome_dot_phenoenriched.jpg'),width = 15, height = 8, units = c("in"))

