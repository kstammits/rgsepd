
# 2022-10-10
# version 1.29.1 in development for R 4.2.1
# is broken again, thanks bioc
# I call out MYH7 and it has changed transcript ids. 
# that select has returned 1:1 message was a warning
#that we have experienced now, my rownames are no longer found in the H2R map



#  library(rgsepd)
data("IlluminaBodymap")
data("IlluminaBodymapMeta")
set.seed(1000) #fixed randomness
x <- Name_to_RefSeq(c("HIF1A","EGFR","MYH7","CD33","BRCA2"))
isoform_ids = intersect(x, rownames(IlluminaBodymap))
rows_of_interest <- unique( c( isoform_ids ,
                               sample(rownames(IlluminaBodymap),
                                      size=2000,replace=FALSE)))
G <- GSEPD_INIT(Output_Folder="OUT",
                finalCounts=round(IlluminaBodymap[rows_of_interest , ]),
                sampleMeta=IlluminaBodymapMeta,
                COLORS=c("green","black","red"))

G <- GSEPD_ChangeConditions( G, c("A","B")) #set testing groups first!    
G <- GSEPD_Process( G ) #have to have processed results to plot them

# looking at genes 2 and 3 will show us a view in dimensions "EGFR" and "MYH7"
# and an axis through five dimensional space.
ExtractProjection(GSEPD=G, txids=isoform_ids, 
                  DRAWING=TRUE, PRINTING=TRUE, GN=c(2,3))

isoform_ids %in% rownames(G$finalCounts)
# two of isoform_ids not present in finalCounts row names
# that's the root of the problem but how did that happen
Name_to_RefSeq("MYH7") %in% rownames(G$finalCounts)
Name_to_RefSeq("MYH7") %in% rownames(IlluminaBodymap)
# so I gather it has changed, since it surely was in that datafile earlier
# my datafile no longer has the entry myh7 is giving. 

# which is a hashmap from H2R

rownames(IlluminaBodymap)
egs = RefSeq_to_EG(rownames(IlluminaBodymap))
nms = EG_to_Name(egs)
which(nms=="MYH7")
nms[238]
egs[238]
rownames(IlluminaBodymap)[238]
# NM_000256
RefSeq_to_EG("NM_000256")
EG_to_Name("4607")
nd = RefSeq_to_Name(rownames(IlluminaBodymap))
which(nd=="MYH7")
rownames(IlluminaBodymap)[239]
# NM_000257 has been updated to another thing. 
# my old dataset is calling it nm 000257
# and the newer org hs eg etc db calls it 
Name_to_RefSeq("MYH7")
#   NM_001407004




