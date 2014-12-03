
GSEPD_ChangeConditions <- function(GSEPD, newConditions){
  if(length(newConditions)==2){
    samples_conds <-  GSEPD$sampleMeta$Condition[GSEPD$sampleMeta$Condition %in% newConditions]
    tconds=table(samples_conds)
    if(any(tconds)<1){#check that we're not missing any
      missingno <- names(tconds[tconds<1])[1]
      om <- sprintf("cannot set the comparison conditions to %s vs %s because sampleMeta$Condition has no one of type %s",
                    newConditions[1],newConditions[2],missingno)
      stop(om)
    }
    conds=paste(samples_conds,tconds[samples_conds],sep=GSEPD$C2T_Delimiter)
    GSEPD$Conditions <- sort(newConditions)
    GSEPD$C2T<- sort(unique( as.character( factor(conds))  ))
    return(GSEPD)
  }else{
    stop("must specify two conditions")
  }
}


GSEPD_ChangeOutput <- function(GSEPD, newFolder){
  if(length(newFolder)==1){
    
    GSEPD$Output_Folder <- newFolder
    
    if(!(file.exists(GSEPD$Output_Folder)))
      dir.create(GSEPD$Output_Folder)
    GSEPD$Output_SCGO <- paste(GSEPD$Output_Folder,"/SCGO",sep="")
    
    if(!(file.exists(GSEPD$Output_SCGO)))
      dir.create(GSEPD$Output_SCGO)
    
    return(GSEPD)
  }else{
    stop("must specify one output folder name")
  }
}

GSEPD_PullDEG <- function(GSEPD, PTHRESH){
  INF<-DESEQ_RFile(GSEPD)
  if(file.exists(INF)){
    DE.data<-read.csv(INF,as.is=TRUE,header=TRUE)
    return( unique(subset(DE.data, DE.data$padj<=PTHRESH)$id ))
  }else{
    om<-sprintf("Failure accessing %s to get the gene list for you.", INF)
    stop(om)
  }
}

