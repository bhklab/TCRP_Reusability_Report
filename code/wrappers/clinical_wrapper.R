library(stringr)
#install.packages("BiocManager")
library(GEOquery)
library(SummarizedExperiment)
library(qs)
library(affy)
library(biomaRt)
library(purrr)
setwd("/data/")
dataset <- "GSE25066"
new_experiment <- readRDS(paste("/data/",dataset,".rds",sep = ""))
# cel.files <- list.files("./CCLE_Expression.Arrays_2013-03-18/")
# sapply(file.path("./CCLE_Expression.Arrays_2013-03-18/", cel.files), GEOquery::gunzip)
# cels <- affy::list.celfiles(file.path("./CCLE_Expression.Arrays_2013-03-18/"), full.names = TRUE)
# rma.norm <- affy::justRMA(filenames = cels, verbose = TRUE, cdfname = NULL)
expression <- assay(new_experiment@ExperimentList$expr_gene_counts)

ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
genes <- rownames(expression)

mapping <- getBM(attributes = c("ensembl_gene_id_version",
                                "hgnc_symbol"),
                 values = genes, mart = ensembl,filters = "ensembl_gene_id_version")
mapping <- mapping[mapping$hgnc_symbol != "",]
mapping <- mapping[!duplicated(mapping$hgnc_symbol), ]

expression <- as.data.frame(expression)
expression <- expression[mapping$ensembl_gene_id_version,]
rownames(expression) <- mapping$hgnc_symbol
system(paste("mkdir /results/",dataset,sep = ""))

samples <- as.data.frame(new_experiment@colData)
sample_ids <- samples$patientid
rownames(samples) <- sample_ids
drug_response <- samples[c("patientid","treatmentid","response")]
#drug_response$characteristics_ch1.10 <- strsplit(drug_response$characteristics_ch1.10,": ") %>% map_chr(`[`, 2)
colnames(drug_response) <- c("line","drug","response")
rownames(drug_response) <- drug_response$line
samples <- samples[drug_response$line,]
expression <- expression[,drug_response$line]
drug_response[(drug_response$response=="NR"),]$response = 0
drug_response[(drug_response$response=="R"),]$response = 1
write.table(expression,paste("/results/",dataset,"/Expression.csv",sep = ""),sep=",")
write.table(drug_response,paste("/results/",dataset,"/Drug_Response.csv",sep = ""),sep=",")
cell_line_table <- samples[c("patientid","tissueid")]
colnames(cell_line_table) <- c("line","site")
cell_line_table$line <- rownames(samples)
cell_line_table$site <- paste(dataset,"_Breast",sep="")
rownames(cell_line_table) <- cell_line_table$line
write.table(cell_line_table,paste("/results/",dataset,"/Cell_Lines_Details.csv",sep = ""),sep = ",")
