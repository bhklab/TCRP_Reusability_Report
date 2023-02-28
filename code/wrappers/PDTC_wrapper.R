install.packages("BiocManager")
BiocManager::install("PharmacoGx")
BiocManager::install("biomaRt")
library(PharmacoGx)
library(biomaRt)
data <- readRDS("/data/PDTXBreast.rds")
data <- updateObject(data)
#needed - expression data 
assay <- t(assay(summarizeMolecularProfiles(data,mDataType = "rna",fill.missing = F)))

assay <- as.data.frame(assay)
samples <- data@sample
samples <- samples[rownames(assay),]
samples <- samples[!duplicated(samples$MODEL),]
assay <- assay[, !duplicated(rownames(assay))]
assay <- assay[,complete.cases(assay)]
models <- samples$MODEL
assay <- assay[rownames(samples),]
rownames(assay) <- models
assay <- assay[, !duplicated(rownames(assay))]
rownames(samples) <- rownames(assay)
write.table(assay,file = "/results/PDTC/Expression.csv",sep = ",")


cell_samples <- data.frame(line=rownames(samples), site=c("PDTC"))

#cell_samples <- data.frame(line=rownames(samples),site=rep("BeatAML",nrow(samples)))
write.table(cell_samples,"/results/PDTC/Cell_Lines_Details.csv",sep = ",")

#total_aac_df <- PharmacoGx::summarizeSensitivityProfiles(data,c("aac_recomputed"))
#total_aac_df <- total_aac_df[,colSums(is.na(total_aac_df))<nrow(total_aac_df)]
#total_aac_df <- data.frame(response=total_aac_df["Paclitaxel",],row.names = colnames(total_aac_df))
total_ic50_df <- PharmacoGx::summarizeSensitivityProfiles(data,c("ic50_published"))
total_ic50_df <- total_ic50_df[,colSums(is.na(total_ic50_df))<nrow(total_ic50_df)]
total_ic50_df <- data.frame(response=total_ic50_df["Paclitaxel",],row.names = colnames(total_ic50_df))
model_samples <- data@sample[rownames(total_ic50_df),]
model_samples <- model_samples[!duplicated(model_samples$MODEL),]
total_ic50_df <- data.frame(ic50 = total_ic50_df[rownames(model_samples),],row.names = rownames(model_samples))
rownames(total_ic50_df) <- model_samples$MODEL

# total_ic50_df <- (total_ic50_df[!is.na(total_ic50_df$response),])
# total_aac_df <- (total_aac_df[!is.na(total_aac_df$response),])

# new_drug_response_df <- data.frame(cell=c(),drug=c(),aac=c(),ic50=c())
# for (i in common_cols) {
#   #cell_line = unlist(rep(list(c(i)), nrow(total_aac_df) + nrow(total_ic50_df)))
#   cell_line = unlist(rep(list(c(i)),nrow(total_ic50_df)))
#   new_drug_response_df <- rbind(new_drug_response_df,data.frame(cell=cell_line,drug=rownames(total_ic50_df),ic50=total_ic50_df[,i]))
# }
total_ic50_df$drug =rep(c("Paclitaxel"),nrow(total_ic50_df))
total_ic50_df <- total_ic50_df[complete.cases(total_ic50_df), ]
total_ic50_df$cell <- rownames(total_ic50_df)
#total_ic50_df <- new_drug_response_df[new_drug_response_df$drug=="Paclitaxel",]
write.table(total_ic50_df,"/results/PDTC/Drug_Response.csv",sep = ",")
