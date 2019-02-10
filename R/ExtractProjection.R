


ExtractProjection <- function(GSEPD, txids, DRAWING=FALSE, GN=c(1,2), PRINTING=FALSE, plotTitle="") {
  
  
  # comparison group sample IDs:
  Group1Set = GSEPD$sampleMeta$Sample[which(GSEPD$sampleMeta$Condition == GSEPD$Conditions[1])]
  Group2Set = GSEPD$sampleMeta$Sample[which(GSEPD$sampleMeta$Condition == GSEPD$Conditions[2])]
  
  #z-score normalization by gene of the samples' normalized counts:
  myPoints <- t(scale( t( GSEPD$normCounts[txids,])))
  
  cols=rep(1,ncol(myPoints))
  pchs=rep(3,ncol(myPoints))
  cols[colnames(myPoints) %in% Group1Set] <- GSEPD$COLORS[1]
  cols[colnames(myPoints) %in% Group2Set] <- GSEPD$COLORS[3]
  pchs[colnames(myPoints) %in% Group1Set] <- 1
  pchs[colnames(myPoints) %in% Group2Set] <- 2
  
  if(PRINTING) {
    print(myPoints)
    print(GN)
    print(Group1Set)
  }
  if(length(Group1Set) > 1)
    AA <- apply(myPoints[,Group1Set],1,mean)#take the mean
  else
    AA <- myPoints[,Group1Set]#take the mean
  if(length(Group2Set)>1)
    AS <- apply(myPoints[,Group2Set],1,mean)
  else
    AS <-  myPoints[,Group2Set]#take the mean
  
  VectorBetween <- (AA - AS)
  uVectorBetween <- VectorBetween/ sqrt(sum(VectorBetween^2)) #unit
  uVectorBetween = as.vector(uVectorBetween)
  #what is the scalar value of each archetype?
  AA.s <- AA %*% uVectorBetween
  AS.s <- AS %*% uVectorBetween 
  
  #then we can easily get the distance to each centroid
  #AA is mean of gene expressions. genes are on the row of myPoints
  Gamma1 <- apply(myPoints , 2, function(x) sqrt(sum(  (x-AA)^2  )))
  Gamma2 <- apply(myPoints , 2, function(x) sqrt(sum(  (x-AS)^2  )))
  #theyre NOT zero-one yet.
  #scale by the distance between them 
  #so it's in units of (AA to AS)
  Gamma1 <- Gamma1 / (sqrt(sum(VectorBetween^2))) ; 
  Gamma2 <- Gamma2 / (sqrt(sum(VectorBetween^2))) ; 
  
  xlim=c(min(AA[GN[1]],min(myPoints[GN[1],]),AS[GN[1]]),
         max(AA[GN[1]],max(myPoints[GN[1],]),AS[GN[1]]))+c(-0.5,0.5)
  ylim=c(min(AA[GN[2]],min(myPoints[GN[2],]),AS[GN[2]]),
         max(AA[GN[2]],max(myPoints[GN[2],]),AS[GN[2]]))+c(-0.5,0.5)
  #each one needs to get projected onto VectorBetween
  alpha=rep(NA,ncol(myPoints)) # a 0-100% scale along the line
  distance_to_line=rep(NA,ncol(myPoints)) #a sample's divergence from the mean/mean line
  
  if(DRAWING==TRUE){
    genenames = c(  DisplayName(txids[GN[1]]) , DisplayName(txids[GN[2]]) )
    plot(cbind(AA, AS)[GN[1],] , cbind(AA, AS)[GN[2],] , type="l", lwd=2,
         xlab=paste(genenames[1],"Z-Score Norm Counts"),
         ylab=paste(genenames[2],"Z-Score Norm Counts"), 
         xlim=xlim, ylim=ylim, main=plotTitle ) 
    text( x=cbind(AA, AS)[GN[1],] , y=cbind(AA, AS)[GN[2],], pos=1, labels=GSEPD$C2T)
    legend("topleft",legend=c(GSEPD$Conditions[1],GSEPD$Conditions[2],"Other"),
           pch=c(1,2,3),col=c(GSEPD$COLORS[1],GSEPD$COLORS[3],1))
  }
  for(i in 1:ncol(myPoints)){ #across people
    scalarProjection <- c((myPoints[,i]-AA) %*% (uVectorBetween))
    scalarProjection <- as.vector(scalarProjection)# new in R 3.5.0
    projected <- scalarProjection * uVectorBetween +AA; 
    distance_to_line[i] <- sqrt(sum( (myPoints[,i] - projected)^2))
    if(DRAWING==TRUE){
      #dots along the line are pch=16
      points(projected[GN[1]],projected[GN[2]], col="black", pch=16);
      points(myPoints[GN[1],i],myPoints[GN[2],i],col=cols[i], pch=pchs[i],cex=2)
      thisGuy <- GSEPD$sampleMeta$SHORTNAME[GSEPD$sampleMeta$Sample==colnames(myPoints)[i]]
      text( x=myPoints[GN[1],i] , y=myPoints[GN[2],i], pos=1, labels=thisGuy,cex=0.5)
    }
    #alpha is the scalar value..
    #now with all these scalars
    #shift their space down by tail and divide out the distance for a 0-100% scale along the line
    alpha[i] <- (scalarProjection - AA.s) / (AS.s-AA.s)
    #theoretically that vector is all zero,zero,zero, one,one,one
  }
  
  #unknown bug has alpha not scaled 0-1. let's do it again.
  alpha_G1 <- mean(alpha[which(colnames(myPoints) %in% Group1Set)]) ;
  alpha=alpha-alpha_G1
  alpha_G2 <- mean(alpha[which(colnames(myPoints) %in% Group2Set)]) ;
  alpha=alpha/alpha_G2
  
  if(DRAWING==TRUE){
    zD <- scale( distance_to_line ) #calculate line styles the same as the white dots on the heatmap
    CF=GSEPD$COLORFUNCTION(100)
    for(i in 1:ncol(myPoints)){ #across people
      scalarProjection <- (myPoints[,i]-AA) %*% (uVectorBetween)
      scalarProjection <- as.vector(scalarProjection)
      projected <- scalarProjection * uVectorBetween + AA; 
      lty=ifelse(zD[i] < GSEPD$VECTOR_DISTANCE_ZTHRESH_Moderate, 1 , 
                 ifelse(zD[i] < GSEPD$VECTOR_DISTANCE_ZTHRESH_Severe, 2 , 3))
      col=CF[cap(round( alpha[i]*100), 1,100 )]
      lines(x=c(projected[GN[1]],myPoints[GN[1],i]), y=c(projected[GN[2]],myPoints[GN[2],i]),lty=lty, col=col)
    }
  }
  
  ValidityScore = Resampled_Significance.k(GSEPD,
                                           ROI=txids,
                                           SOI=c(Group1Set,Group2Set))
  
  x=list();
  x$alpha=alpha
  x$beta=distance_to_line / nrow(myPoints)
  #adding a "gamma" so we can measure the distance to each centroid
  x$gamma1=Gamma1
  x$gamma2=Gamma2
  x$Validity.Score=ValidityScore$Validity
  x$Validity.P=ValidityScore$PV
  x
}


