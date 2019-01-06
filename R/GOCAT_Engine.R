###
## Summer 2013 - 2014
## rgsepd-0.3.3 KStamm 
## a tool to perform linear vector projection in GO-Term space
###

GSEPD_ProjectionProcessor <- function(GSEPD) {
  
  GSEPD <- GSEPD_CheckCounts(GSEPD) # make sure normCounts is calculated!
  finalCounts <- GSEPD$normCounts
  C2T <- GSEPD$C2T
  sampleMeta <- GSEPD$sampleMeta
  
  Type1=GSEPD$Conditions[1]
  Type2=GSEPD$Conditions[2]
  
  DefinedGroup <- subset(sampleMeta, sampleMeta$Condition %in% c(Type1, Type2))
  ArcheTypes <- as.character(DefinedGroup$Sample)
  G1=ArcheTypes[DefinedGroup$Condition == Type1];
  G2=ArcheTypes[DefinedGroup$Condition == Type2]
  cols=rep(1,ncol(finalCounts))
  pchs=rep(3,ncol(finalCounts))
  
  if(length(GSEPD$COLORS)!=3){
    stop("Vector of GSEPD$COLORS should be three elements (Group1, middle, Group2).")
  }#doesnt this really belong elsewhere. 
  
  cols[colnames(finalCounts) %in% G1] <- GSEPD$COLORS[1]
  cols[colnames(finalCounts) %in% G2] <- GSEPD$COLORS[3]
  pchs[colnames(finalCounts) %in% G1] <- 1
  pchs[colnames(finalCounts) %in% G2] <- 2
  
  #pull up the merged doc and grab the minor GO term
  M.data <- read.csv(GSEPD_MFile(GSEPD), as.is=TRUE, header=TRUE,row.names=1) 

  AssembleGOSTAT <- function(GOTERM){
    x=list();
    x$Term=GOTERM;
    DEG=subset(M.data,M.data$category==GOTERM)
    x$Name=DEG$Term.x[1]
    x$DataProjection <- ExtractProjection(GSEPD, DEG$REFSEQ ,
                                          PRINTING = FALSE, DRAWING = FALSE)
    x
  }
  
  tcats <- table(M.data$category)
  cats <- unique( subset( M.data, M.data$over_represented_padj.x < GSEPD$LIMIT$GO_PVAL)$category)
  #omit those with too few or too many
  cats = setdiff(cats, names(tcats)[tcats<GSEPD$MinGenesInSet | tcats>GSEPD$MaxGenesInSet]  ) 

  if(length(cats)>1) {
    outfile=paste(GSEPD$Output_Folder, "/SCA.GSEPD.",C2T[1],".",C2T[2],".pdf",sep="")
    Message_Generate(outfile)
    pdf(outfile)
    for(i in seq(1,length(cats)-1,2)){
      #print(paste(cats[i],cats[i+1]))
      DIM_1=AssembleGOSTAT(cats[i])
      DIM_2=AssembleGOSTAT(cats[i+1])
      xs=(DIM_1$DataProjection$alpha) #GOSet 1 axis score
      ys=(DIM_2$DataProjection$alpha) #GOSet 2 axis score
      plot(x=xs, y=ys, main=paste(ncol(finalCounts), "Samples by 2 GO:Terms"),
           xlab=paste(Type2,DIM_1$Name,Type1, sep=" - "), ylab=paste(Type2,DIM_2$Name,Type1, sep=" - "),
           pch=pchs, col=cols )
      legend("topleft",legend=c(GSEPD$Conditions[1],GSEPD$Conditions[2],"Other"),pch=c(1,2,3),col=c(GSEPD$COLORS[1],GSEPD$COLORS[3],1))
    }
    dev.off()
  }else{
    warning("Skipping generation of SCA.GSEPD.pdf because there aren't enough enriched GO Categories.")
  }
  
  if(length(cats)>0) {
    message("Calculating Projections and Segregation Significance")
    OM=matrix(nrow=length(cats), ncol=ncol(finalCounts))
    WhiteOut <- OM;  G1OM <- OM;  G2OM <- OM;
    Segregation_PV=rep(1,nrow(OM)) # see GSEPD$LIMIT$Seg_P later.
    Segregation_Val=rep(1,nrow(OM)) # see GSEPD$LIMIT$Seg_P later.
    
    colnames(WhiteOut)<-colnames(finalCounts)
    rownames(OM)<-cats
    rownames(WhiteOut)<-cats
    rownames(G1OM)<-cats
    rownames(G2OM)<-cats
    names(Segregation_PV)<- cats
    names(Segregation_Val)<- cats
    
    for(i in 1:nrow(OM)){
      ExDP <- AssembleGOSTAT(cats[i])$DataProjection
      OM[i,] <-ExDP$alpha
      WhiteOut[i,]<-ExDP$beta
      G1OM[i,] <-ExDP$gamma1
      G2OM[i,] <-ExDP$gamma2
      #Segregation_P[i] <- Resampled_Significance.t(GSEPD,data.frame(alpha=ExDP$alpha, type=GSEPD$sampleMeta$Condition ), GSEPD$Conditions)
      O <- Resampled_Significance.k(GSEPD, 
                      ROI=M.data$REFSEQ[M.data$category==cats[i]],
                      SOI=ArcheTypes)
      Segregation_PV[i] <- O$PV
      Segregation_Val[i] <- O$Validity
    }
    colnames(OM)<-colnames(finalCounts)
    colnames(WhiteOut)<-colnames(finalCounts)
    colnames(G1OM)<-colnames(finalCounts)
    colnames(G2OM)<-colnames(finalCounts)
    SIGFIG=3
    write.csv(signif(OM,SIGFIG), paste(GSEPD$Output_Folder, "/GSEPD.Alpha.", C2T[1],".",C2T[2],".csv",sep=""))
    write.csv(signif(WhiteOut,SIGFIG), paste(GSEPD$Output_Folder, "/GSEPD.Beta.", C2T[1],".",C2T[2],".csv",sep=""))
    write.csv(cbind(Segregation_Val, Segregation_PV), GSEPD_Seg_File(GSEPD))
    write.csv(signif(G1OM,SIGFIG), GSEPD_HMG1CSV_File(GSEPD))
    write.csv(signif(G2OM,SIGFIG), GSEPD_HMG2CSV_File(GSEPD))
    
    GOPs <- rep(0,length(cats))
    GONames <- rep("",length(cats))
    for(i in 1:length(cats)) { 
      roi <- (M.data$category == cats[i]) #find the entry in Mdata
      GONames[i] <- sprintf("%s %s",M.data$Term.x[roi][1], M.data$category[roi][1]); 
      GOPs[i] <- M.data$over_represented_padj.x[roi][1];
      #if this GOSet fails to separate the classes, we can drop it from further processing
      if( Segregation_PV[i] < GSEPD$LIMIT$Seg_P){
        #if it segregates let's make a pairs plot to go with
        deg_ngenes <- sum(roi)
        if(deg_ngenes >= 2 && deg_ngenes <= 10){# only on those sets with a decent number of genes
          GSEPD_Pairs_Plot(GSEPD, M.data$REFSEQ[roi], GONames[i], GOPs[i], Segregation_PV[i])
        }
        pdf(paste(GSEPD$Output_SCGO,"/Scatter.",C2T[1],".",C2T[2],".GO",substring(cats[i],4),".pdf",sep=""))
          GSEPD_PCA_Spec(GSEPD=GSEPD,GOT=cats[i],MDATA=M.data)
          title(sub=sprintf("%d Genes / %d Samples / GOSeq p=%f / segregation p=%f",deg_ngenes, ncol(finalCounts), GOPs[i], Segregation_PV[i]))
        dev.off()
      }
    }
    
    #we'll use the combination to choose great sets
    GO_SumPs <- GOPs + Segregation_PV
    
    Message_Generate(GSEPD_HMA_File(GSEPD))
   
    if(length(G1)>1 && nrow(OM)>1)
      OM.A1 <- apply(OM[,G1], 1, mean)-min(OM)+1
    else
      OM.A1 <-OM[,G1] - min(OM)+1
    
    if(length(G2)>1 && nrow(OM)>1)
      OM.A2 <- apply(OM[,G2], 1, mean)-min(OM)+1
    else
      OM.A2 <-OM[,G2] - min(OM)+1
    
    OM.min <- apply(cbind(OM.A1,OM.A2),1,min)#*0.90
    OM.max <- apply(cbind(OM.A1,OM.A2),1,max)#*1.10
    #Z-scored OutMatrix
    zOM <- OM-min(OM)+1;
    for(i in 1:nrow(zOM)){
      zOM[i, zOM[i,]< OM.min[i] ] <- OM.min[i]
      zOM[i, zOM[i,]> OM.max[i] ] <- OM.max[i]
      zOM[i,] <- (zOM[i,] - OM.min[i]) / (OM.max[i] - OM.min[i])
    }#goes pretty much 0-1;

    #generate the asterisks you see in the heatmap
    #first lets get some Z-scores on the projection distance
    zWhiteOut <- t(scale(t(WhiteOut)))
    cellnote=matrix(" ", nrow=nrow(zOM),ncol=ncol(zOM))
    for(j in 1:nrow(zOM)){
     cellnote[j,  zWhiteOut[j,] > GSEPD$VECTOR_DISTANCE_ZTHRESH_Moderate] <- "." 
     cellnote[j,  zWhiteOut[j,] > GSEPD$VECTOR_DISTANCE_ZTHRESH_Severe ] <- rawToChar(as.raw(149)) #a fat circle bullet character
    }
    
    ColumnLabels <- colnames(zOM);
    for(j in 1:length(ColumnLabels))
      ColumnLabels[j] <- paste(sampleMeta$SHORTNAME[sampleMeta$Sample==ColumnLabels[j]],
                               sampleMeta$Condition[sampleMeta$Sample==ColumnLabels[j]])
    
    ColLabelColors <- rep("black",ncol(OM));
    names(ColLabelColors) <- colnames(OM)
    ColLabelColors[G1] <- GSEPD$COLORS[1];
    ColLabelColors[G2] <- GSEPD$COLORS[3];
    
    #which rows are worth looking at??
    sr=1:length(GO_SumPs) 
    sr=sort(GO_SumPs, decreasing=FALSE, index.return=TRUE)$ix[1:min(GSEPD$MAX_GOs_for_Heatmap,length(cats))]
    #and remove those non-significant
    sr <- sr[ GO_SumPs[sr] <= (GSEPD$LIMIT$GO_PVAL + GSEPD$LIMIT$Seg_P)] ;
    
    if(length(sr) > 1) {
      RowLabelColors <- rep("black",length(sr));
      names(RowLabelColors)<-cats[sr]
      sMdata<-subset(M.data, !duplicated(M.data$category), select=c("category","GOSEQ_DEG_Type"))
      for(j in 1:length(sr)){
        cj<-cats[sr[j]];DEGType<-sMdata$GOSEQ_DEG_Type[sMdata$category == cj];
        RowLabelColors[j] <- ifelse(DEGType == "Up",GSEPD$COLORS[3],
                                    ifelse(DEGType == "Down",GSEPD$COLORS[1],
                                           GSEPD$COLORS[2]))
      } ; rm(sMdata)
      
      pdf(GSEPD_HMA_File(GSEPD), height=6+length(sr)/7, width=6+ncol(zOM)*0.33)
      heatmap.2(zOM[sr,], labRow=GONames[sr],
            scale="none", trace="none", margins=c(10,35),
            cexRow=1.25, labCol=ColumnLabels, cellnote=cellnote[sr,],
            notecex=3, notecol="white",
            ColSideColors=ColLabelColors, RowSideColors=RowLabelColors,
            col=GSEPD$COLORFUNCTION);
      dev.off()
      
      #then plot HMG#####
      #it's overlaying the zOM
      #but for the HMG file now I need to make a new 0-1 score from the G1OM/G2OM matrices.
      zOM[,] <- 0.5 ; #default all to indeterminate/black
      zOM[G1OM < 0.45] <- 0.35 #dark green
      zOM[G1OM < (0.45*(2/3)) ] <- 0 # bright green
      zOM[G2OM < 0.45] <- 0.65 # dark red
      zOM[G2OM < (0.45*(2/3)) ] <- 1 # bright red
                   
      pdf(GSEPD_HMG_File(GSEPD), height=6+length(sr)/7, width=6+ncol(zOM)*0.33)
      heatmap.2(zOM[sr,], labRow=GONames[sr],
                scale="none", trace="none", margins=c(10,35),
                cexRow=1.25,labCol=ColumnLabels,
                ColSideColors=ColLabelColors, RowSideColors=RowLabelColors,
                col=GSEPD$COLORFUNCTION);
      dev.off()
      
    }else{
      warning("Not generating HMA/HMG files:  <2 significant rows. See your GOSEQ and Segregation tables for details.")
    }
    
    #for the rows seen in the heatmap...
    for(j in sr){ #make lots of scatter plots, one for each GO term with some enrichment
      thisDEG=subset(M.data, M.data$category==cats[j])
      if( nrow(thisDEG) >= GSEPD$MinGenesInSet ) { 
#thought this would be taken care of above, but the GO term whitelist could create entries with invalid gene counts.
        pdf(paste(GSEPD$Output_SCGO,"/GSEPD.", C2T[1],".",C2T[2],".GO",substring(cats[j],4),".pdf",sep=""))
        nGenes=length(unique(thisDEG$REFSEQ))
        if(nGenes>1){
          for(i in seq(1,nGenes-1,2))
              ExtractProjection(GSEPD, txids = unique(thisDEG$REFSEQ),
                                DRAWING=TRUE, GN=c(i,i+1),PRINTING=FALSE,
                                plotTitle=thisDEG$Term.x[1])
        }
        dev.off();
      }
    }

  }else{ #length(cats)<1
    warning("GSEPD : Can't make output heatmap images when zero GO-Categories enriched.")
  }

}  # end GOCATPROCESSOR for Type1, Type2



