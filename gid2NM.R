# the input file you gave me:
x1 = read.csv("~/Desktop/gene_symbol.csv", as.is=TRUE)
head(x1)
library(biomaRt)

ensembl = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
attrs=listAttributes(ensembl)
write.csv(attrs,"attrs.csv")
aa=c("refseq_mrna" , "chromosome_name" , "start_position" ,  "end_position" ,
     "entrezgene_id" ,      "hgnc_symbol",      "description" );

bmt_NM=getBM(attributes = aa, filters = "hgnc_symbol",values=x1$Gene_Symbol, mart = ensembl)

x1$refseq = unlist(lapply( 1:nrow(x1), function(i) {
  i = which(bmt_NM$hgnc_symbol == x1$Gene_Symbol)
  if(length(i) >0)
    return(sample( bmt_NM$refseq_mrna[i],1)) # picking one at random
  return("NA")
  
  }))
write.csv(x1,"~/Desktop/gene_symbol.mod.csv")

