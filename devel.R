#devel workspace  2020-06-03
# adding and testing functionality for better data input
# and better data output
# maybe non-default diagramming
# but of course the vignette needs the diagrams  ... 

install.packages("~/Workspace/rgsepd_1.29.1.tgz", repos=NULL)
library(rgsepd)
# following DESeq2's examples of using data from several sources.

library("pasilla")
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnno <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
#read.csv(pasCts, sep="\t",nrow=10)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnno, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

#-==--
library(rgsepd)

sampleMeta = data.frame(Sample=rownames(coldata),
                        Condition= coldata$condition )

G <- GSEPD_INIT(Output_Folder="OUTZi",
                finalCounts= cts,
                sampleMeta=sampleMeta,
                COLORS=c("green","black","red"))
G$LIMIT$LFC=0.1
G <- GSEPD_ChangeConditions( G, unique(coldata$condition)) #set testing groups  


G$GeneIDSystem <- list(Transcript="REFSEQ",
                           GeneID="ENTREZ",
                           GeneName="HGNC")

G$GeneIDSystem <- list(Transcript="ENSEMBL",
                           GeneID="ENSEMBL",
                           GeneName="ENSEMBL",
                           Species="DMel")

for(x in list.files("./R",pattern="*R", full.names=TRUE))
  source(x)



G <- GSEPD_Process( G ) #have to have processed results to plot them



X <- vector(c("q","r","s"), names=c("a","b","c"))
#gns=hash::values(HGNC[rownames(sfc)])



BiocManager::install("HSMMSingleCell")
library(HSMMSingleCell)

??HSMMSingleCell

?HSMMSingleCell

install.packages("monocle")
BiocManager::install("monocle")
library(monocle)

x <- load_HSMM()
