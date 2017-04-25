setwd("/Users/hegemb/mesothelioma/")
data.s1<-read.table("quant_S1.sf.txt",header=TRUE)
data.s2<-read.table("quant_S2.sf.txt",header=TRUE)
data.s3<-read.table("quant_S3.sf.txt",header=TRUE)
data.s4<-read.table("quant_S4.sf.txt",header=TRUE)
countData<-cbind(data.s1$TPM,data.s2$TPM,data.s3$TPM,data.s4$TPM)
row.names(countData)<-data.s1[,1]
colnames(countData)<-c("s1","s2","s3","s4")
countData<-round(countData)
colData<-matrix(c("after",rep("before",3)),ncol=1)
colnames(colData)<-"condition"
rownames(colData)<-colnames(countData)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData,colData=colData,design=~condition)
dds <- DESeq(dds)
res <- results(dds)
ix<-which(res$padj<0.1)
o<-order(res$padj)
res2<-res[o,]
#top20
v<-match(rownames(res2)[1:20],rownames(res))
countData[v,]
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#K<-length(ix)
K<-5
R<-NULL
for (k in 1:K){
	print(k)
	results <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
                 filters = "ensembl_transcript_id", values=substr(row.names(res2)[k],1,15),mart = ensembl)
	R<-rbind(R,results)    
}             
res3<-cbind(R,res2[1:K,],countData[o[1:K],])
write.table(res3,file="deseq-res.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)