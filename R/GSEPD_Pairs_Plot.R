
GSEPD_Pairs_Plot <- function(GSEPD, IDs, GO_Term, GO_P, Seg_P){
  
  fc <- GSEPD$normCounts[IDs,]
  rownames(fc) <- DisplayName(IDs)
  
  Sanitized_GOName <- strsplit(GO_Term,":")[[1]] #assuming this string ends with GO:NNNN, we cut at the colon
  outfile=sprintf("%s/Pairs.%s.%s.GO%s.pdf", GSEPD$Output_SCGO , GSEPD$C2T[1], GSEPD$C2T[2],
                  Sanitized_GOName[length(Sanitized_GOName)] )
  
  cols=rep(1,ncol(fc))
  cols[as.character(GSEPD$sampleMeta$Condition)==GSEPD$Conditions[1]]=GSEPD$COLORS[1]
  cols[as.character(GSEPD$sampleMeta$Condition)==GSEPD$Conditions[2]]=GSEPD$COLORS[3]
  pchs=rep(3,ncol(fc));
  pchs[as.character(GSEPD$sampleMeta$Condition)==GSEPD$Conditions[1]]=1
  pchs[as.character(GSEPD$sampleMeta$Condition)==GSEPD$Conditions[2]]=2
  
  square_dimension = length(IDs)*1.55+2
  pdf(outfile, width = square_dimension, height = square_dimension)
  pairs( t(fc) , main=GO_Term,
        sub=sprintf("%d Genes / %d Samples / GOSeq p=%f / segregation p=%f",length(IDs), ncol(fc), GO_P, Seg_P),
        pch=pchs, col=cols)
  dev.off();
  
}

