
GSEPD_PCA_Plot <- function(GSEPD){

  sampleMeta<- GSEPD$sampleMeta
  Condition<- GSEPD$Conditions
  if(is.null(Condition)) 
    stop("Must specify conditions to test with GSEPD_ChangeConditions before GSEPD_PCA_Plot")
   #if only due to the CheckCounts below.
  
  COLORS=GSEPD$COLORS
  if(!(all(Condition %in% sampleMeta$Condition))) 
    stop("specified conditions not found in sample metadata table")
  
  GSEPD<-GSEPD_CheckCounts(GSEPD) #ensure #normCounts exists....
  fc <- GSEPD$normCounts
  dropRows <- apply(fc,1,mean) < 0 #normCounts are DESeq2 VST, nearly log2 space
  fc<-fc[!dropRows,]
  sfc<-t(scale(t( (fc))))
  dropRows <- apply(sfc,1,function(x){any(is.nan(x))})
  sfc=sfc[!dropRows,]
  
  if(nrow(sfc) <= ncol(sfc)){
    #not enough DEG to bother PCAing...
    warning("fewer expressed genes than subjects: won't PCA")
    return();
  }
  
  sms=colnames(sfc)
  cols=rep(1,ncol(sfc))
  cols[as.character(sampleMeta$Condition)==Condition[1]]=COLORS[1]
  cols[as.character(sampleMeta$Condition)==Condition[2]]=COLORS[3]
  pchs=rep(3,ncol(sfc));
  pchs[as.character(sampleMeta$Condition)==Condition[1]]=1
  pchs[as.character(sampleMeta$Condition)==Condition[2]]=2
  
  pc.md<-prcomp(sfc)
  
  lbs=unlist(lapply(1:ncol(sfc),function(i){sampleMeta$SHORTNAME[sampleMeta$Sample==colnames(sfc)[i]]}))
  dims=c("PC1","PC2")
  outfilename=sprintf("%s/GSEPD.PCA_AG.%s.%s.pdf",GSEPD$Output_Folder,GSEPD$C2T[1],GSEPD$C2T[2])
  Message_Generate(outfilename)
  pdf(outfilename,width=6,height=6)
    plot(pc.md$rotation[,dims], xlab="", ylab="", col=cols, pch=pchs)
    x=pc.md$rotation[,dims[1]]
    y=pc.md$rotation[,dims[2]]
    xo=rep(0,ncol(sfc)) ; yo=xo;
    text(x=x+xo,y=y-0.015+yo,labels=lbs, xpd=TRUE, cex=0.75)
    ANNOTE_GENES=min(c(5,nrow(sfc))) #the most important genes are 
    PC1= names(sort( abs(pc.md$x[,dims[1]]) ,
                     decreasing=TRUE )[1:ANNOTE_GENES])
    gns=DisplayName(PC1);
    mtext(side=1,text=gns,at=seq(min(x),max(x),length.out=ANNOTE_GENES),line=3)
    PC2= names(sort( abs(pc.md$x[,dims[2]]) ,
                     decreasing=TRUE )[1:ANNOTE_GENES])
    gns=DisplayName(PC2);
    mtext(side=2,text=gns,at=seq(min(y),max(y),length.out=ANNOTE_GENES),line=3) 
    legend(GSEPD$PCA_LEGEND,legend=c(Condition,"Other"), pch=c(1,2,3), col=c(COLORS[c(1,3)],1))
    title(main=paste("PCA over",nrow(sfc),"transcripts"))
  dev.off();
  
  if(file.exists(DESEQ_RFile(GSEPD))){
    DE.data<-read.csv(DESEQ_RFile(GSEPD),as.is=TRUE,header=TRUE)
    DE.data <- subset(DE.data, !(DE.data$id %in% GSEPD$EXCLUDES))
    keepRows <- c();
    if(GSEPD$LIMIT$HARD){
      keepRows=DE.data$id[ DE.data$padj <= GSEPD$LIMIT$PADJ &
                             abs(DE.data$log2FoldChange) >= GSEPD$LIMIT$LFC &
                             DE.data$baseMean >= GSEPD$LIMIT$baseMean ]
    }else{ #soft limits:
      TargetCount = (ncol(sfc)+1)
      if(sum(DE.data$padj<0.05) > TargetCount){
        keepRows=DE.data$id[DE.data$padj<0.05]
      }else if(sum(DE.data$pval<0.01) > TargetCount) {
        keepRows=DE.data$id[DE.data$pval<0.01]
      }else if(sum(DE.data$pval<0.05) > TargetCount){
        keepRows=DE.data$id[DE.data$pval<0.05]
      }else {
        pthresh=sort(DE.data$pval,decreasing=FALSE)[TargetCount]
        keepRows=DE.data$id[DE.data$pval<pthresh]
      }
    }
    keepRows=keepRows[keepRows%in%rownames(sfc)]
    sfc=sfc[keepRows,]
    
    if(nrow(sfc) < 3){
      #not enough DEG to bother PCAing...
      warning("less than three filtered genes: can't meaningfully PCA for DEG")
      return();
    }

    pc.md<-prcomp(sfc)
    outfilename=sprintf("%s/GSEPD.PCA_DEG.%s.%s.pdf",GSEPD$Output_Folder,GSEPD$C2T[1],GSEPD$C2T[2])
    Message_Generate(outfilename)
    pdf(outfilename,width=6,height=6)
      plot(pc.md$rotation[,dims], xlab="", ylab="", col=cols, pch=pchs)
      x=pc.md$rotation[,dims[1]]
      y=pc.md$rotation[,dims[2]]
      xo=rep(0,ncol(sfc)) ; yo=xo;
      text(x=x+xo,y=y-0.015+yo,labels=lbs, xpd=TRUE, cex=0.75)
      ANNOTE_GENES=min(c(5,nrow(sfc))) #the most important genes are 
      PC1= names(sort( abs(pc.md$x[,dims[1]]) ,
                       decreasing=TRUE )[1:ANNOTE_GENES])
      gns=DisplayName(PC1)
      mtext(side=1,text=gns,at=seq(min(x),max(x),length.out=ANNOTE_GENES),line=3)
      PC2= names(sort( abs(pc.md$x[,dims[2]]) ,
                       decreasing=TRUE )[1:ANNOTE_GENES])
      gns=DisplayName(PC2)
      mtext(side=2,text=gns,at=seq(min(y),max(y),length.out=ANNOTE_GENES),line=3) 
      legend(GSEPD$PCA_LEGEND,legend=c(Condition,"Other"), pch=c(1,2,3), col=c(COLORS[c(1,3)],1))    
      title(main=paste("PCA over",nrow(sfc),"transcripts"))
    dev.off();
  }
}



GSEPD_PCA_Spec <- function(GSEPD, GOT, MDATA=NULL){
  
  sampleMeta<- GSEPD$sampleMeta
  Condition<- GSEPD$Conditions
  COLORS=GSEPD$COLORS
  if(!(all(Condition %in% sampleMeta$Condition))) 
    stop("specified conditions not found in sample metadata table")
  
  ## pick up Mdata
  #pull up the merged doc and grab the minor GO term
  if(missing(MDATA) || is.null(MDATA)) {
   if(file.exists(GSEPD_MFile(GSEPD)))
      M.data <- read.csv(GSEPD_MFile(GSEPD), as.is=TRUE, header=TRUE,row.names=1) 
    else
      stop("Cannot PCA_Spec the given GO term because the .MERGE file is not found, be sure to run a full GSEPD_Process() first.")
  }
  else
    M.data <- MDATA
  
  HGNC <- hash(M.data$REFSEQ, M.data$HGNC)
  roi <- (M.data$category==GOT) #find the entry in Mdata
  GOName<-sprintf("%s %s",M.data$Term.x[roi][1], M.data$category[roi][1]); 
  #select GOI
  GOI = M.data$REFSEQ[roi];
  fc <- GSEPD$normCounts[GOI,]
  NS=ncol(GSEPD$normCounts)
  NG=length(GOI)
  
  cols=rep(1,NS)
  cols[as.character(sampleMeta$Condition)==Condition[1]]=COLORS[1]
  cols[as.character(sampleMeta$Condition)==Condition[2]]=COLORS[3]
  pchs=rep(3,NS);
  pchs[as.character(sampleMeta$Condition)==Condition[1]]=1
  pchs[as.character(sampleMeta$Condition)==Condition[2]]=2
  lbs=unlist(lapply(1:NS,function(i){
    sampleMeta$SHORTNAME[sampleMeta$Sample==colnames(fc)[i]]}))
  
  if(NG<1){
    stop(sprintf("Cannot GSEPD_PCA_Spec ... no genes found in %s for category %s",GSEPD_MFile(GSEPD),GOT))
  }else if(NG<2){ ## == 1 ... histogram??
    hist(fc, main=sprintf("Normalized Counts at %s",DisplayName(GOI)), xlab="Reads")
    points(x=fc, y=runif(length(fc))*5+2, cex=2,pch=pchs,col=cols)    
  }else if(NG<3){ ## == 2 ... regular scatter plot
    sfc=fc
    gns=hash::values(HGNC[rownames(sfc)])
    plot(x=sfc[1,], y=sfc[2,],
         xlab=sprintf("%s Normalized Counts",gns[1]),
         ylab=sprintf("%s Normalized Counts",gns[2]),
         main=sprintf("Scatter for %s",GOName),
         pch=pchs, col=cols )    
    xo=rep(0,NS) ; yo=xo;
    text(x=sfc[1,]+xo,y=sfc[2,]-0.015+yo,labels=lbs, xpd=TRUE, cex=0.75)
  }else { ## prcomp when >2 genes
    sfc=t(scale(t(fc)))
    pc.md<-prcomp(sfc)
    dims=c("PC1","PC2")
    plot(pc.md$rotation[,dims], xlab="", ylab="", col=cols, pch=pchs, main=sprintf("PCA for %s",GOName))
    x=pc.md$rotation[,dims[1]]
    y=pc.md$rotation[,dims[2]]
    xo=rep(0,NS) ; yo=xo;
    text(x=x+xo,y=y-0.015+yo,labels=lbs, xpd=TRUE, cex=0.75)
    ANNOTE_GENES=min(c(5,nrow(sfc))) #the most important genes are 
    PC1= names(sort( abs(pc.md$x[,dims[1]]) ,
                     decreasing=TRUE )[1:ANNOTE_GENES])
    gns=hash::values(HGNC[PC1])
    mtext(side=1,text=gns,at=seq(min(x),max(x),length.out=ANNOTE_GENES),line=3)
    PC2= names(sort( abs(pc.md$x[,dims[2]]) ,
                     decreasing=TRUE )[1:ANNOTE_GENES])
    gns=hash::values(HGNC[PC2])
    mtext(side=2,text=gns,at=seq(min(y),max(y),length.out=ANNOTE_GENES),line=3) 
  }
  legend(GSEPD$PCA_LEGEND,legend=c(Condition,"Other"), pch=c(1,2,3), col=c(COLORS[c(1,3)],1))
}
