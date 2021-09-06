library(DESeq2)
library(dplyr)
#Reading from 6 files for each accession
#df <- read.table("~/counts_gene_norm_1_rpkm.txt", header = FALSE) %>% select(2)
#y <- c("~/counts_gene_norm_2_rpkm.txt","~/counts_gene_norm_3_rpkm.txt","~/counts_gene_norm_4_rpkm.txt","~/counts_gene_norm_5_rpkm.txt", "~/counts_gene_norm_6_rpkm.txt")


#for(i in y){
#        df <- cbind(df, read.table(i,header = FALSE) %>% select(2))
#}
#Reading from one file with all accessions
#df <- read.table("~/quartile_gene_norm.txt", header = TRUE) %>% select(2,3,4,5,6,7)
#print(1)
#df <- format(df, scientific = F)
#print(2)
df_1 <- read.table("~/Downloads/heart_norm_DESeq2.txt", sep = " ", header = TRUE)
#Processing data

df_2 <- read.table("~/Downloads/fibroblasts_norm_DESeq2.txt", sep = " ", header = TRUE)
samples <- sample.int(ncol(df_2), ncol(df_1))
df_2 <- df_2 %>% select(samples)

#df <- cbind(as.double(df_1), as.double(df_2))
df <- cbind(df_1, df_2)
df <- format(df, scientific = F)
mat <- data.matrix(df)
#print(3)
print(ncol(mat))

df_2 <- read.table("~/Downloads/decitabine_counts.txt", sep = "\t", header = TRUE)
#print(ncol(df_1))
#print(ncol(df_2))
#samples <- sample.int(ncol(df_2), 610)
#df_2 <- df_2 %>% select(samples)
#samples <- sample.int(ncol(df_1), 610)
#df_1 <- df_1 %>% select(samples)
#df_1 <- df_1[ ,c(1:ncol(df_2)) ]
#df <- cbind(as.double(df_1), as.double(df_2))
#df <- cbind(df_1, df_2)
df <- df_2
genes = df$GeneID
#df <- format(df, scientific = F)
genes

df = within(df, rm(GeneID, Length))
rownames(df) = genes
View(df)
mat <- data.matrix(df)
condition <- factor(c(replicate(3, "Heart"), replicate(3, "fibroblasts")))
dds <- DESeqDataSetFromMatrix(mat, DataFrame(condition), ~condition)

# DE analyses
dds <- DESeq(dds, test = "Wald")
res <- results(dds)
#genes <- read.table("~/counts_gene_norm_1_deseq.txt", header = FALSE) %>% select(1)
#res <- cbind(genes, res)

write.table(res, "~/DE_decitabine_counts_deseq.txt", quote = FALSE, col.names = TRUE, row.names = TRUE)

