##
## Some helper functions.



RefSeq_to_EG <- function(nm){
  if(length(nm)>1)
    return(unlist(lapply(nm,RefSeq_to_EG)))
  if(length(nm)==0)
    return(NA)
  #else; do one:
  #hash::values(get("R2E", pkg_globals)[nm])
  get("R2E", pkg_globals)[[nm]]
}

EG_to_Name <- function(eg){
  if(length(eg)>1)
    return(unlist(lapply(eg,EG_to_Name)))
  if(length(eg)==0)
    return(NA)
  #else; do one:
  #hash::values(get("E2H", pkg_globals)[eg])
  get("E2H", pkg_globals)[[eg]]
}

Name_to_RefSeq <- function(x){
  if(length(x)>1)
    return(unlist(lapply(x,Name_to_RefSeq)))
  if(length(x)==0)
    return(NA)
  #else; do one:
  #hash::values(get("H2R", pkg_globals)[x])
  get("H2R", pkg_globals)[[x]]
}

RefSeq_to_Name <- function(x){
  if(length(x)>1)
    return(unlist(lapply(x,RefSeq_to_Name)))
  if(length(x)==0)
    return(NA)
  #else; do one:
  #hash::values(get("H2R", pkg_globals)[x])
  y = get("R2H", pkg_globals)[[x]]
  if(length(y)==0)
    return("NA")
  y
}


DisplayName <- function(txid){
  if(length(txid)>1)
    return(unlist(lapply(txid,DisplayName)))
  if(length(txid)==0)
    return("NA")
  #else; do one:
  #hash::values(get("R2H", pkg_globals)[txid])
  get("R2H", pkg_globals)[[txid]]
}







## of course now these tests wont work until after
# .onLoad has set up the hashes. 
#so this is the wrong place for tests...
if(FALSE) {
  RefSeq_to_EG("NM_153215")
  EG_to_Name("132228")
  Name_to_RefSeq("LSMEM2")
  
  #test:::
  stopifnot("BRCA1"==EG_to_Name(RefSeq_to_EG(Name_to_RefSeq("BRCA1"))))
  stopifnot("CD33"==EG_to_Name(RefSeq_to_EG(Name_to_RefSeq("CD33"))))
  stopifnot("CD33"==DisplayName(Name_to_RefSeq("CD33")) )
}


build_Annotation_Matrix <- function(nms){
  #from a list of NM refseq IDs, (and NR), we want a matrix of three columns
  O<-matrix(NA,nrow=length(nms),ncol=3)
  colnames(O)<-c("REFSEQ","HGNC","ENTREZ")
  O[,"REFSEQ"]<-nms;
  Es<-RefSeq_to_EG(nms)
  O[,"ENTREZ"]<-Es
  #EG to name doesnt like NAs.
  ROI <- !is.na(O[,"ENTREZ"])
  O[ROI,"HGNC"]<-EG_to_Name(O[ROI,"ENTREZ"])
  return(O)
}
  
#

EnsureCleanIDs <- function(NM){
 #make sure all NMs are available in R2E and R2H before it causes a crash
  hR2H <- get("R2H", pkg_globals)
  missingKeys <- NM[!(NM %in% names(hR2H))]
  if(length(missingKeys)>0) {
    for(j in 1:length(missingKeys))
      hR2H[[missingKeys[j]]] <- missingKeys[j]
    assign("R2H", hR2H, pkg_globals)
  }
  
  hR2E <- get("R2E", pkg_globals)
  missingKeys <- NM[!(NM %in% names(hR2E))]
  if(length(missingKeys)>0) {
    for(j in 1:length(missingKeys))
       hR2E[[missingKeys[j]]] <- missingKeys[j]
    assign("R2E", hR2E, pkg_globals)
  }
  
  hE2H <- get("E2H", pkg_globals)
  allEG <- as.character(unlist(as.list(hR2E)) )

  missingKeys <- allEG[!(allEG %in% names(hE2H))]
  if(length(missingKeys)>0) {
    for(j in 1:length(missingKeys))
      hE2H[[missingKeys[j]]] <- missingKeys[j]
    # their assigned HGNC is equal to their EntrezGeneID
    assign("E2H", hE2H, pkg_globals)
  }
  
}


