install.packages("BiocManager")
BiocManager::install("PharmacoGx")
BiocManager::install("biomaRt")
library(PharmacoGx)
library(biomaRt)
data <- readRDS("/data/gCSI2.rds")
data <- updateObject(data)
#needed - expression data 
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
assay <- t(assay(summarizeMolecularProfiles(data,mDataType = "Kallisto_0.46.1.rnaseq",fill.missing = F)))
genes <- gsub(rep="", x=colnames(assay), pat="\\.[0-9]+$")
G_list <- getBM(filters= "ensembl_gene_id", attributes= c( "ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
colnames(assay) <- genes
assay <- as.data.frame(assay)
assay <- assay[,G_list$ensembl_gene_id]
colnames(assay) <- G_list$hgnc_symbol
write.table(assay,file = "/scratch/CTRPv2/Expression.csv",sep = ",")
samples <- data@sample
samples <- samples[rownames(assay),]
cell_samples <- data.frame(line=rownames(samples),pharma_id=samples$PharmacoDB.id,site=samples$tissueid,histology=samples$Cellosaurus.Disease.Type)

write.table(cell_samples,"/scratch/CTRPv2/Cell_Lines_Details.csv",sep = ",")

mutation <- assay(data@molecularProfiles$mutation)
mutation <- mutation[,intersect(colnames(mutation),rownames(assay))]
mutation[mutation == "wt"] <- 0
mutation[is.na(mutation)] <- 0
mutation[mutation != 0] <- 1
mutation <- t(mutation)
mutation <- as.data.frame(mutation)
write.table(mutation,file = "/scratch/CTRPv2/Mutations.csv",sep = ",")
total_aac_df <- PharmacoGx::summarizeSensitivityProfiles(data,c("aac_recomputed"))
total_aac_df <- total_aac_df[,colSums(is.na(total_aac_df))<nrow(total_aac_df)]
total_ic50_df <- PharmacoGx::summarizeSensitivityProfiles(data,c("ic50_recomputed"))
total_ic50_df <- total_ic50_df[,colSums(is.na(total_ic50_df))<nrow(total_ic50_df)]
common_drugs <- intersect(rownames(total_aac_df),rownames(total_ic50_df))
common_cols <- intersect(colnames(total_aac_df),colnames(total_ic50_df))
new_drug_response_df <- data.frame(cell=c(),drug=c(),aac=c(),ic50=c())
for (i in common_cols) {
  cell_line = unlist(rep(list(c(i)), nrow(total_aac_df) + nrow(total_ic50_df)))
  new_drug_response_df <- rbind(new_drug_response_df,data.frame(cell=cell_line,drug=rownames(total_aac_df),aac=total_aac_df[,i],ic50=total_ic50_df[,i]))
}

new_drug_response_df <- new_drug_response_df[complete.cases(new_drug_response_df), ]
write.table(new_drug_response_df,"/scratch/CTRPv2/Drug_Response.csv",sep = ",")
