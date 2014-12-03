
GSEPD_DEGHeatmap <- function(G){
  
  infile=DESEQ_AFile(G)
  
  data=read.csv(infile,as.is=TRUE,header=TRUE)
  data=subset(data,!duplicated(data$ENTREZ))
  sdata=subset(data, data$ENTREZ != "NA" & data$baseMean >= G$LIMIT$baseMean)
  DEG=ifelse(sdata$PADJ <= G$LIMIT$PADJ & abs(sdata$LOG2.X.Y.) >= G$LIMIT$LFC , 1,0)
  
  if(!(G$LIMIT$HARD)) { #if the limit is not hard, move it around
    ToLevel=round((0.0010*length(DEG)))
    PThresh = sort(sdata$PVAL,decreasing=FALSE)[ToLevel]
    if(sum(DEG)<(0.0010*length(DEG))){
      #not enough genes... lets move that PADJ threshold
      DEG=ifelse(sdata$PVAL <= PThresh , 1,0)
      message(sprintf("Not many genes found differentially expressed, (re)moving the significance threshold to raw p=%f so we can use %d genes in the heatmap.",PThresh,sum(DEG)))
    }
    if(G$MAX_Genes_for_Heatmap>nrow(sdata))
      G$MAX_Genes_for_Heatmap=nrow(sdata)
    if(sum(DEG)>G$MAX_Genes_for_Heatmap){
      PThresh = sort(sdata$PVAL,decreasing=FALSE)[G$MAX_Genes_for_Heatmap]
      DEG=ifelse(sdata$PVAL <= PThresh , 1,0)
      message(sprintf("Too many genes found differentially expressed, changing the  threshold to raw p=%f so we can use %d genes in the heatmap.",PThresh,sum(DEG)))
    }
  }
  ROI=sdata$REFSEQ[DEG==1]
  
  plotFile=sprintf("%s/HM.%s.%s.%d.pdf",G$Output_Folder, G$C2T[1],G$C2T[2],sum(DEG))
  Message_Generate(plotFile) #let the user know how we're busy
  if(length(ROI)<2){
    warning("Can't make heatmaps of less than 2 rows. Check your DESEQ RES Annote_Filter file.")
    return();
  }
  
  G<-GSEPD_CheckCounts(G)
  fc<-G$normCounts[,G$sampleMeta$Sample] #doublecheck the column order
  
  fc=fc[ROI,] #subset rows
  hmData <- fc
  cellnote=signif(fc,digits=2)
  
  Sample_Columns <- which(G$sampleMeta$Condition %in% G$Conditions)
  rowMeans <- apply(hmData[,Sample_Columns],1,mean)
  rowSD <- apply(hmData[,Sample_Columns],1,sd)+0.01
  zData <- cap( t(t(hmData-rowMeans)/rowSD),-3,3)
  rm(rowMeans,rowSD)
  rownames(zData)<- DisplayName(rownames(zData))############
  
  ColumnLabels <- colnames(zData)
  for (j in 1:length(ColumnLabels)) 
    ColumnLabels[j] <- paste( G$sampleMeta$SHORTNAME[G$sampleMeta$Sample ==  ColumnLabels[j]],
                              G$sampleMeta$Condition[G$sampleMeta$Sample ==  ColumnLabels[j]])
  
  ColLabelColors <- rep(G$COLORS[2],ncol(fc));
  ColLabelColors[G$sampleMeta$Condition == G$Conditions[1]] <- G$COLORS[1];
  ColLabelColors[G$sampleMeta$Condition == G$Conditions[2]] <- G$COLORS[3];
  
  
  #primary, everyone, and large
  pdf(plotFile,width=max(c(3.5,ncol(fc)*0.30))+1,height=max(c(2,sum(DEG)*0.2))+2)
  heatmap.2(zData,  
            scale="none", trace="none", dendrogram="both", margins=c(10,10),
            cellnote=cellnote, notecol="black", notecex=0.75,labCol=ColumnLabels,
            ColSideColors=ColLabelColors, col=G$COLORFUNCTION)
  dev.off()
  
  plotFile=sprintf("%s/HM-.%s.%s.%d.pdf",G$Output_Folder, G$C2T[1],G$C2T[2],sum(DEG))
  #and a small one:
  pdf(plotFile,width=max(c(4,ncol(fc)*0.22))+1,height=max(c(2,sum(DEG)*0.1))+2)
  heatmap.2(zData,  labRow="", key=FALSE,  margins=c(10,10),
            scale="none", trace="none",dendrogram="none",labCol=ColumnLabels,
            ColSideColors=ColLabelColors, col=G$COLORFUNCTION)
  dev.off()
  
  COI <- (ColLabelColors != G$COLORS[2])
  
  plotFile <- sprintf("%s/HMS.%s.%s.%d.pdf",G$Output_Folder, G$C2T[1],G$C2T[2],sum(DEG))
  pdf(plotFile,width=max(c(4,sum(COI)*0.30))+1,height=max(c(2,sum(DEG)*0.2))+2)
  heatmap.2(zData[,COI],  
            scale="none", trace="none", dendrogram="both", margins=c(10,10),
            cellnote=cellnote[,COI], notecol="black", notecex=0.75,labCol=ColumnLabels[COI],
            ColSideColors=ColLabelColors[COI], col=G$COLORFUNCTION)
  dev.off()
  
  plotFile=sprintf("%s/HMS-.%s.%s.%d.pdf",G$Output_Folder, G$C2T[1],G$C2T[2],sum(DEG))
  pdf(plotFile,width=max(c(4,sum(COI)*0.22))+1,height=max(c(2,sum(DEG)*0.1))+2)
  heatmap.2(zData[,COI],  labRow="", key=FALSE,  margins=c(10,10),
            scale="none", trace="none",dendrogram="none",labCol=ColumnLabels[COI],
            ColSideColors=ColLabelColors[COI], col=G$COLORFUNCTION)
  dev.off()

}
