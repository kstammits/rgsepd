
GSEPD_INIT <- function( Output_Folder="OUT" ,finalCounts=NULL, sampleMeta=NULL,
                        DESeqDataSet = NULL,
                        COLORS = c("green", "gray", "red"),
                        C2T = "x" ){
  
  GSEPD=list(); #master object antipattern
  
  if(is.null(DESeqDataSet)){
  #old interface, expect FC and SM.
    if(is.null(finalCounts) || is.null(sampleMeta)){
      stop(" not enough data given, you must specify either DESeqDataSet or finalCounts and sampleMeta")
    }else{#this would be a nice place to sanity check
      sampleMeta$Sample <- as.character(sampleMeta$Sample)
      sampleMeta$Condition <- as.character(sampleMeta$Condition)
      GSEPD$finalCounts <- finalCounts
      GSEPD$sampleMeta <- sampleMeta
    }
  }else{
    #new interface, DDS given.
    if(is.null(finalCounts) || is.null(sampleMeta)){
      GSEPD$finalCounts <- counts(DESeqDataSet)
      GSEPD$sampleMeta <- DESeqDataSet@colData
    }else{
      stop(" too much data given, you must specify either DESeqDataSet or finalCounts and sampleMeta")
    }
  }
  
  KeepRows<-rowSums(GSEPD$finalCounts)>0  & !is.na(rowSums(GSEPD$finalCounts))
  if(sum(KeepRows) < nrow(GSEPD$finalCounts))
     message(sprintf("Keeping rows with counts (%d of %d)",sum(KeepRows),length(KeepRows))) 
  GSEPD$finalCounts <- GSEPD$finalCounts[KeepRows,] ;rm(KeepRows)
     
  GSEPD$Output_Folder = Output_Folder #new files are generated into this sub-folder
  
  if(!(file.exists(GSEPD$Output_Folder)))
    dir.create(GSEPD$Output_Folder)
  
  GSEPD$Output_SCGO <- paste(GSEPD$Output_Folder,"/SCGO",sep="")
  
  if(!(file.exists(GSEPD$Output_SCGO)))
    dir.create(GSEPD$Output_SCGO)
  
  GSEPD$MinGenesInSet <- 2  ## filters which GOTerm/GeneSets need to be displayed
  GSEPD$MaxGenesInSet <- 30 #preventing the analysis of enormous sets
  
  GSEPD$LIMIT <- list(HARD=TRUE,  # should I stick with the given thresholds, or let them vary, default: stick to it
                      baseMean=20,  # defines which genes are expressed
                      LFC=1, # defines which genes are differentially expressed for the purposes of gene set enrichment
                      PADJ=0.050,  # defines which genes are differentially expressed for the purposes of gene set enrichment
                      GO_PVAL=0.050, # defines which pathways are to be included in the heatmap figure (MERGE.over_represented_pvalue < x)                      
                      Seg_P=0.050 ) #only PCA/Scatter plot those GOSets with segregating ability.
  
  GSEPD$Segregation_Precision <- 0.01 #about how certain do we have to be in the permuted kmeans. We'll empirically calculate between 1 and 4 times the reciprocal of this number. See GroupSignificance.R for permutation loops.

  GSEPD$MAX_GOs_for_Heatmap <- 40 #caps the size of the final heatmap to only the N most significant sets (by GOP+SegP)
  
  GSEPD$MAX_Genes_for_Heatmap <- 50 #caps the size of the DEG heatmap to only the N most significant Genes
          #in the event you have more than that you'll have to see the .Annote table
          #not applicable if LIMIT$HARD == TRUE, we'll report all of them
  
  GSEPD$COLORS <- COLORS #the heatmap (AND OTHER PLACES) will use this scale

  GSEPD$C2T_Delimiter <- C2T # five items of A versus four items of B is labeled Ax5.Bx4 . choose a symbol you don't use in condition names!
  
  #force the legend to appear in a vacant corner of the plot.
  #but I don't have a good way of choosing the vacant corner yet
  GSEPD$PCA_LEGEND <- "topright"
  
  GSEPD$QUIET <- FALSE #when true, most status and progress messages are suppressed.
  
  #for the heatmap and vector projection's color definition
  #this is a scale between group 1 and group2's profiles.
  #by default I double the first and last entry, so the center is squished.
  GSEPD$COLORFUNCTION <- colorRampPalette( c(GSEPD$COLORS[1],GSEPD$COLORS[1],
                                             GSEPD$COLORS[2],
                                             GSEPD$COLORS[3],GSEPD$COLORS[3] ))
  
  #we set a cutoff for how divergent a sample is before being marked on the heatmap.
  GSEPD$VECTOR_DISTANCE_ZTHRESH_Moderate <- 2
  GSEPD$VECTOR_DISTANCE_ZTHRESH_Severe <- 3
  #value is the z-score of the distance of all samples to the line
  
  #BY default empty
  GSEPD$EXCLUDES<-c("NM_0000000"); #vector of NM ids as strings
  GSEPD$Force_GO_Include<-c("GO:0000000"); #vector of GO ids as strings
  
  #theoretically you can change the background used for goseq
  #but i've only ever used/tested Entrez IDs on hg19
  GSEPD$GOSEQ <- list(genome="hg19" ,
                      system="knownGene",
                      use_genes_without_cat=TRUE)
  
  GSEPD$GeneIDSystem <- list(Transcript="REFSEQ",
                             GeneID="ENTREZ",
                             GeneName="HGNC")
  
  
  return(GSEPD)
}

