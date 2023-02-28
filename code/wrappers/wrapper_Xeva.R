
library(Xeva)
library(biomaRt)
data <- readRDS("/data/Xeva_PDXE.rds")

#needed - expression data 
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
assay <- getMolecularProfiles(data,data.type = "RNASeq")
assay <- assay@assayData$exprs

cell_line_names <- names(data@experiment)
model_name_list <-c()
drug_name_list <- c()
for (i in cell_line_names) {
  model_name_list <- append(model_name_list,i)
  drug_name_list <- append(drug_name_list,data@experiment[[i]]@drug$join.name)
}
drug_name_df <- data.frame(cell=model_name_list,drugs=drug_name_list)
#drug_name_df <- drug_name_df[drug_name_df$drugs=='paclitaxel',]
drug_name_df <- drug_name_df[drug_name_df$drugs=="erlotinib",]
interested_models <- drug_name_df$cell
samples <- data@model
samples <- samples[samples$tissue.name=="Non-small Cell Lung Carcinoma",]
samples <- samples[interested_models,]
samples <- samples[!is.na(samples$patient.id),]
interested_patients <- intersect(samples$patient.id,colnames(assay))
assay <- assay[,interested_patients]
samples = samples[!duplicated(samples$patient.id),]
rownames(samples) <- samples$patient.id
samples <- samples[common_models,]
models <- samples$model.id
samples <-samples[,c("tissue.name","patient.id")]
colnames(samples) <- c("Site","Line")
samples <- samples[!is.na(samples$Site),]
samples[samples$Site=="Non-small Cell Lung Carcinoma",]$Site = "PDX_Lung"
#cell_samples <- data.frame(line=rownames(samples),site=rep("BeatAML",nrow(samples)))
write.table(assay,file = "/results/PDXE/Expression.csv",sep = ",")

write.table(samples,"/results/PDXE/Cell_Lines_Details.csv",sep = ",")

mutation <- getMolecularProfiles(data,data.type = "mutation")
mutation_df <- mutation@assayData$exprs
mutation_df <- as.data.frame(mutation_df)
mutation_df[mutation_df != 0] = 1
mutation_df <- mutation_df[,interested_patients]
write.table(mutation,file = "/results/PDXE/Mutations.csv",sep = ",")
drug_response <- data@sensitivity$model
drug_response <- drug_response[interested_models,]
rownames(drug_response) <- interested_models
drug_response <- drug_response[,c("model.id","mRECIST","AUC")]
# total_auc_df <- PharmacoGx::summarizeSensitivityProfiles(data,c("aac_recomputed"))
# total_auc_df <- total_aac_df[,colSums(is.na(total_aac_df))<nrow(total_aac_df)]
# total_ic50_df <- PharmacoGx::summarizeSensitivityProfiles(data,c("ic50_recomputed"))
# total_ic50_df <- total_ic50_df[,colSums(is.na(total_ic50_df))<nrow(total_ic50_df)]
# common_drugs <- intersect(rownames(total_auc_df),rownames(total_ic50_df))
# common_cols <- intersect(colnames(total_auc_df),colnames(total_ic50_df))
# new_drug_response_df <- data.frame(cell=c(),drug=c(),aac=c(),ic50=c())
# for (i in common_cols) {
#   cell_line = unlist(rep(list(c(i)), nrow(total_auc_df) + nrow(total_ic50_df)))
#   new_drug_response_df <- rbind(new_drug_response_df,data.frame(cell=cell_line,drug=rownames(total_auc_df),aac=total_auc_df[,i],ic50=total_ic50_df[,i]))
# }
drug_response$model.id <- rownames(drug_response)
write.table(drug_response,"/results/PDXE/Drug_Response.csv",sep = ",")
