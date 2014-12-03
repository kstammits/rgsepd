
## V-Measure for external cluster validation
## a conditional entropy-based evaluation measure
#a harmonic mean of the entropic homogeneity and completeness
# Rosenberg & Hirschberg 2007, 

Homogeneity<-function(A){
  N <- sum(A)
  A<-A[,colSums(A)>0] #filter that C=0 column
  #print(A)
  HCK <- 0- sum(mapply(function(K)
    sum(mapply(function(C)
      ifelse(A[K,C]==0, 0, 
             (A[K,C]/N)*log(A[K,C]/sum(A[K,]),2) )
               ,1:ncol(A)))
                       ,1:nrow(A)))
  
  if(HCK==0) return(1);
  HC <- 0- sum(mapply(function(C){
    sum(A[,C])/N * log( sum(A[,C])/N ,2)
  },1:ncol(A)))
  
  return(1-HCK/HC)
}
Completeness <- function(A){
  N <- sum(A)
  A<-A[,colSums(A)>0] #filter that C=0 column
  HCK <- 0- sum(mapply(function(C)
    sum(mapply(function(K)
      ifelse(A[K,C]==0,0,
      (A[K,C]/N)*log(A[K,C]/sum(A[,C]),2) )
               ,1:nrow(A)))
                       ,1:ncol(A)))
  
  if(HCK==0) return(1);
  HK <- 0- sum(mapply(function(K){ sum(A[K,])/N * log( sum(A[K,])/N ,2)  },1:nrow(A)))
  
  return(1-HCK/HK)
}
VMeasure <- function(A, beta=1){
  H<-Homogeneity(A)
  C<-Completeness(A)
  V<-((1+beta)*H*C)/(beta*H + C)
  return( ifelse(is.nan(V),0,V)  )
}


#ROI is REFSEQ IDs, rows of the counts table.
#SOI is subject names; columns of the counts table.

Resampled_Significance.k<- function(GSEPD, ROI, SOI) {
  D<-GSEPD$normCounts[ROI,SOI]
  
  #scale induces a NaN if var==0
  #omit that: it wont help the clustering.
  rowVariation <- apply(D,1,var)
  #i expect we've induced a new crash potential, when <2 rows left
  keepRows <- (rowVariation>0.050)
  if(sum(keepRows) < 2){
    #but im not sure how to handle it. 
    #for now, call it a losing cluster with no segregation.
    return(list(Validity=0, PV=1))
  }
  D<-D[rowVariation>0.050,]

  D<-t(scale(t(  D  )))
  ob<-kmeans(t(D),centers=2, nstart=4)
  
  S1 <- names(ob$cluster[ob$cluster==1])
  S2 <- names(ob$cluster[ob$cluster==2])
  A <- Contingency_Matrix(GSEPD,S1,S2)
  TS <- VMeasure(A)
  PV <- 0
  #to compute random labels we'll repeatedly ask for the number of samples
  lSOI <- length(SOI) #so doing that once can be faster
  sample_rate <- length(S1)/lSOI
  Stepcount = pmax(round((1/GSEPD$Segregation_Precision)),2)
  if(TS < 1){
    TS_Range <- unlist(lapply(1:Stepcount, function(i){
      RB <- runif(lSOI) < sample_rate
      VMeasure(Contingency_Matrix(GSEPD,RB,!RB))
    }))
    PV<- mean( TS_Range >= TS )
    if(PV>0 && PV<0.50){ #some were better? let's go for more precision:
      TSR2 <- unlist(lapply(1:Stepcount, function(i){
        RB<-runif(lSOI) < sample_rate
        VMeasure(Contingency_Matrix(GSEPD,RB,!RB))
      }))
      TS_Range=c(TS_Range,TSR2)
      PV<- mean( TS_Range >= TS )
      #a third refinement for good luck
      if(PV>0.01 && PV<0.10){ # let's go for more precision:
        TSR2 <- unlist(lapply(1:(2*Stepcount), function(i){
          RB<-runif(lSOI) < sample_rate
          VMeasure(Contingency_Matrix(GSEPD,RB,!RB))
        }))
        TS_Range=c(TS_Range,TSR2)
        PV<- mean( TS_Range >= TS )
      }
    }
  }
  return(list(Validity=TS, PV=PV))
}

Contingency_Matrix <- function(GSEPD, G1,G2){
  #Group1 and Group2 are Sample names for Conditions1 and 2.
  Con<-GSEPD$sampleMeta$Condition
  names(Con)<-GSEPD$sampleMeta$Sample
  T1 <- c(sum(Con[G1]==GSEPD$Conditions[1]), sum(Con[G1]==GSEPD$Conditions[2]))
  T2 <- c(sum(Con[G2]==GSEPD$Conditions[1]), sum(Con[G2]==GSEPD$Conditions[2]))
  rbind(T1,T2)
}
