# updating for BIOC 3.11 , R. 4.0 on osx 10.12 needs this magic incantation to install:
#Sys.setenv(XML_CONFIG = "/usr/local/Cellar/libxml2/2.9.10/bin/xml2-config")
# after of course putting libxml there. 
# not necesary for osx 10.15.
#



if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "devel")
BiocManager::valid()   
BiocManager::install("BiocGenerics")
BiocManager::install("goseq")

BiocManager::install("geneLenDataBase")
BiocManager::install("biomaRt")
BiocManager::install("AnnotationDbi")

library(rgsepd)

# attempt to recreate the original pictures in the manuscript. that was a comparison from this dataset

fc = read.csv("~/Documents/OUT_pairwise_D2016/DESEQ.counts.D0x2.D1x2.csv", header=TRUE, as.is=TRUE, row.names=1)

head(fc)
colnames(fc)
library(rgsepd)

sampleMeta = data.frame(Sample=colnames(fc),
                        Condition=rep(c("D0","D1","D3","D5","D8","D14"),each=2))
sampleMeta$SHORTNAME = c("D0I12","D0I6","D1T3","D1T4","3I12","D3I6",
                         "D5T7","D5T8","8T10","D8T9","4T11","4T12")

G <- GSEPD_INIT(Output_Folder="OUTZ2",
                finalCounts= round(fc*100),
                sampleMeta=sampleMeta,
                COLORS=c("green","black","red"))
G$LIMIT$LFC=0.1
G <- GSEPD_ChangeConditions( G, c("D3","D5")) #set testing groups first!           
G <- GSEPD_Process( G ) #have to have processed results to plot them


GOFILE = read.csv("OUTZ2/GSEPD.RES.D3x2.D5x2.GO2.csv", header=TRUE, as.is=TRUE)
sG = subset(GOFILE, Term=="cardiac atrium morphogenesis")
which(sG$HGNC == "NOTCH1")
which(sG$HGNC == "PITX2")
pdf()
G$sampleMeta$SHORTNAME = c("D0I12","D0I6","D1T3","D1T4","3I12","D3I6",
                         "D5T7","D5T8","8T10","D8T9","4T11","4T12")

ExtractProjection(GSEPD=G, txids=sG$REFSEQ,
                  DRAWING=TRUE, PRINTING=FALSE, GN=c(16,21),
                  plotTitle = "cardiac atrium morphogenesis")
dev.off()




library(rgsepd)

data("IlluminaBodymap")
data("IlluminaBodymapMeta")
set.seed(1000) #fixed randomness
isoform_ids <- Name_to_RefSeq(c("HIF1A","EGFR","MYH7","CD33","BRCA2"))
rows_of_interest <- unique( c( isoform_ids ,
                               sample(rownames(IlluminaBodymap),
                                      size=5000,replace=FALSE)))

G <- GSEPD_INIT(Output_Folder="OUT",
                finalCounts=round(IlluminaBodymap[rows_of_interest , ]),
                sampleMeta=IlluminaBodymapMeta,
                COLORS=c("green","black","red"))

G <- GSEPD_ChangeConditions( G, c("A","B")) #set testing groups first!    
G <- GSEPD_Process( G ) #have to have processed results to plot them

# looking at genes 2 and 3 will show us a view in dimensions "EGFR" and "MYH7"
# and an axis through five dimensional space.
ExtractProjection(GSEPD=G, txids=isoform_ids,
                  DRAWING=TRUE, PRINTING=FALSE, GN=c(2,3))






data("IlluminaBodymap")
data("IlluminaBodymapMeta")
set.seed(1000) #fixed randomness
isoform_ids <- Name_to_RefSeq(c("HIF1A","EGFR","MYH7","CD33","BRCA2"))
rows_of_interest <- unique( c( isoform_ids ,
                               sample(rownames(IlluminaBodymap),
                                      size=20000,replace=FALSE)))
G <- GSEPD_INIT(Output_Folder="OUT",
                finalCounts=round(IlluminaBodymap[rows_of_interest , ]),
                sampleMeta=IlluminaBodymapMeta,
                COLORS=c("green","black","red"))

#loosen statistics to ensure this example gets some gene sets 'enriched'
#G$LIMIT$LFC=0.1; G$LIMIT$GO_PVAL=0.5        
G <- GSEPD_ChangeConditions( G, c("A","B")) #set testing groups first!    
G <- GSEPD_Process( G ) #have to have processed results to plot them
