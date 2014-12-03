
#plot selected values of your dataset
# rgsepd-0.3.3 June 2014 KStamm
#make sure your G object is finished 
#and that your specified genes are rownames of finalCounts
#version 0.99.7 adding some error handling November2014-KStamm


GSEPD_Heatmap <- function(G,genes,cap.range=3,cellnote="log10"){
  
  G<-GSEPD_CheckCounts(G) #have to have finished GSEPD_Process() before we can use it
  
  KeeperRows <- genes %in% rownames(G$normCounts)
  if(any(KeeperRows==FALSE))
    message(sprintf("Warning: %d gene ids are not found in count data",sum(KeeperRows==FALSE)))
  genes<-genes[KeeperRows];rm(KeeperRows)
  
  fc<-G$normCounts[genes,G$sampleMeta$Sample] #doublecheck the column order
  fc=fc[genes,] #subset rows
  hmData <- fc
  if(cellnote=="log10")
  cellnote=signif(fc,digits=2)
  
  Sample_Columns <- which(G$sampleMeta$Condition %in% G$Conditions)
  rowMeans <- apply(hmData[,Sample_Columns],1,mean)
  rowSD <- apply(hmData[,Sample_Columns],1,sd)+0.01
  zData <- cap( t(t(hmData-rowMeans)/rowSD),0-cap.range,cap.range)
  rownames(zData)<- mapply( function(x){
    j<-DisplayName(x); if(is.na(j) || nchar(j)<2) j<-x; j}, rownames(zData))
  
  ColumnLabels <- colnames(zData)
  for (j in 1:length(ColumnLabels)) 
    ColumnLabels[j] <- paste( G$sampleMeta$SHORTNAME[G$sampleMeta$Sample ==  ColumnLabels[j]],
                              G$sampleMeta$Condition[G$sampleMeta$Sample ==  ColumnLabels[j]])
  
  ColLabelColors <- rep(G$COLORS[2],ncol(fc));#defaults to mid
  ColLabelColors[G$sampleMeta$Condition == G$Conditions[1]] <- G$COLORS[1];
  ColLabelColors[G$sampleMeta$Condition == G$Conditions[2]] <- G$COLORS[3];
  heatmap.2(zData,  
            scale="none", trace="none", dendrogram="both", margins=c(10,10),
            cellnote=cellnote, notecol="black", notecex=0.75,labCol=ColumnLabels,
            ColSideColors=ColLabelColors, col=G$COLORFUNCTION)  
}

