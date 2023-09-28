likFunc.munic<-function(par){
  
  extTranMatDat.tmp$pars$beta<-par[1]
  extTranMatDat.tmp$pars$gamma<-par[2]
  
  ##### Put this back to fit more
  nInfecLoc<-rep(1,9)
  nInfecLoc[1:8]<-exp(par[3:10])
  nInfecLoc<-nInfecLoc/sum(nInfecLoc)
  ################################
  
  
  # tmp.pHome.MUNIC<-pop2019.town
  # tmp.pHome.MUNIC<-tmp.pHome.MUNIC/sum(tmp.pHome.MUNIC)
  tmp.pHome.PROV<-extTranMatDat.tmp$popbyCell
  tmp.pHome.PROV<-tmp.pHome.PROV/sum(tmp.pHome.PROV)
  
  
  
  #### gravity model parameter fit
  popsize<-pop2019.town
  # beta=1.5
  # beta1=log(extTranMatDat.tmp$pars$beta/(1+extTranMatDat.tmp$pars$beta))
  beta=exp(extTranMatDat.tmp$pars$beta)
  
  ### reciprocal logit
  # beta=1/(1+exp(-extTranMatDat.tmp$pars$beta))
  # beta <- min.range+beta1*(max.range-min.range)
  
  
  # tmppar1 <- exp(beta_tmp)/(1+exp(beta_tmp))
  # gamma=2
  gamma=exp(extTranMatDat.tmp$pars$gamma)
  dists<-pairwise_geodist.town
  diag(dists)<-10
  n<-dim(dists)[1]
  probMove<-matrix(NaN,n,n)
  for(i in 1:n){
    for(j in 1:n){
      # probMove[i,j]<-popsize[i]*popsize[j]^beta/dists[i,j]^gamma
      probMove[i,j]<-log(popsize[j])^beta/dists[i,j]^gamma
      
    }
  }
  rowtot<-rowSums(probMove)
  probMove_preadj <- apply(probMove, 2, function(x) x/rowtot  )
  
  
  
  ###  method to adjust for infectious period
  timeWindow<-35
  probStay<-1-(diag(probMove_preadj))^timeWindow
  # probStay[which(probStay>0.99999)]<-0.99999
  tmp<-probMove_preadj
  diag(tmp)<-0
  tmp1<-sweep(tmp,1,rowSums(tmp),"/")
  tmp2<-sweep(tmp1,1,(1-probStay)/(1-diag(tmp1)),"*")
  diag(tmp2)<-probStay
  tmpbase<-tmp2
  
  
  ####### Collapse mobility matrix to 9X9 here
  nloc.munic=nrow(tmpbase)
  
  splitnames <- matrix(nrow=nloc.munic,ncol=2)
  for (us in 1:nloc.munic) {
    for (prov in 1:2) { 
      splitnames[us,prov] <- strsplit(names(pop2019.town),"_")[[us]][prov]
    }
  }
  
  tmpbase.1 <- apply(tmpbase, 1, function(x){
    split(x,splitnames[,2]) %>% sapply(sum)
  })  %>% transpose()

  
  tmpbase2 <- apply(tmpbase.1,2,function(x) {
    mapply( weighted.mean
            , x = split(x,splitnames[,2])
            , w = split(pop2019.town,splitnames[,2])
    )
  })
  
  move3<-tcrossprod(tmpbase2,tmpbase2)
  move4<-sweep(move3,2,tmp.pHome.PROV,"*")
  TranMat.tmp2<-sweep(move4,1,rowSums(move4),"/")
  
  
  ######### WEIGHT BY DESTINATION POPULATION BEFORE COLLAPSING
  
  # move3<-tcrossprod(tmp.sick,tmp)
  # move4<-sweep(move3,2,tmp.pHome,"*")
  # TranMat.tmp<-sweep(move4,1,rowSums(move4),"/")

  
#  # #### COLLAPSE BEFORE MATRIX MULTIPLICATION AFTER WEIGHTING DESTINATION
#   nloc.munic=nrow(TranMat.tmp)
#   ### Create index for municipality to province
#   splitnames <- matrix(nrow=nloc.munic,ncol=2)
#   for (us in 1:nloc.munic) {
#     for (prov in 1:2) { 
#       splitnames[us,prov] <- strsplit(names(pop2019.town),"_")[[us]][prov]
#     }
#   }
#   
#   TranMat.tmp.1 <- apply(TranMat.tmp, 1, function(x){
#     split(x,splitnames[,2]) %>% sapply(sum)
#   })  %>% transpose()
  
  # # ### Collapse rows to province level by weighting by the population sizes
  # TranMat.tmp2 <- apply(TranMat.tmp.1,2,function(x) {
  #   mapply( weighted.mean
  #     , x = split(x,splitnames[,2])
  #     , w = split(pop2019.town,splitnames[,2])
  #   )
  # })
  
######RUN MATRIX MULTIPLICATION
  gensA<-gensB<-(1:maxGen)
  maxgen.tmp<-maxGen
  TranMatArray<-array(NA,c(nlocs,nlocs,maxgen.tmp))
  TranMatArray[,,1]<-TranMat.tmp2
  for (j in 2:maxgen.tmp){TranMatArray[,,j]<-eigenMapMatMult(TranMatArray[,,j-1], TranMat.tmp2)
  }
  
  
# #### Test Applying
#   ix = c(1,1,2)
#   xx = array(1:27, dim=c(3,3,3))
#   apply(xx, c(1,3), function(x){
#     split(x,ix) %>% sapply(sum)
#   })  %>%
#   aperm(c(2,1,3))
  
### True Condense from 234X234XmaxGen into 234X9XmaxGen#### 
  # TranMatArray<-apply(TranMatArray.1, c(1,3), function(x){
  #   split(x,splitnames[,2]) %>% sapply(sum)
  # })  %>%
  #   aperm(c(2,1,3))
  
  
  
  mrcaVec<-extTranMatDat.tmp$popbyCell
  # mrcaVec<-pop2019.town
  
  llIndPair<-function(ii){
    
    probByGenA<-probByGen.tmp[[dat.in2$GPSC[ii]]][dat.in2[ii,1],dat.in2[ii,2],1:maxGen]
    probByGenB<-probByGen.tmp[[dat.in2$GPSC[ii]]][dat.in2[ii,2],dat.in2[ii,1],1:maxGen]
    
    minProb=1e-4
    
    gensA<-(1:(maxGen))[which(probByGenA>minProb)]
    gensB<-(1:(maxGen))[which(probByGenB>minProb)]
    
    
    # nInfecLoc<-mrcaVec*extTranMatDat.tmp$popbyCell
    
    ndetect<-NoDetectByYearGPSC[,dat.in2$GPSC[ii],dat.in2$YearREF1[ii]]
    # ndetect<-NoDetectByYear[,dat.in2$YearREF1[ii]]
    probByloc1<-ndetect/nInfecLoc
    
    ndetect<-NoDetectByYearGPSC[,dat.in2$GPSC[ii],dat.in2$YearREF2[ii]]
    # ndetect<-NoDetectByYear[,dat.in2$YearREF2[ii]]
    probByloc2<-ndetect/nInfecLoc
    
    # probByloc1<-probByloc2<-probByloc
    
    # probByloc1[]<-probByloc2[]<-1
    
    TranMatArrayA<-TranMatArray[,,gensA]
    # TranMatArrayA<-sweep(TranMatArrayA,3,probByGenA[gensA],"*")
    if (length(gensA)>1) { TranMatArrayA<-sweep(TranMatArrayA,3,probByGenA[gensA],"*")} else {TranMatArrayA<-TranMatArrayA*probByGenA[gensA]}
    TranMatArrayA<-sweep(TranMatArrayA,2,probByloc1,"*")
    
    TranMatArrayB<-TranMatArray[,,gensB]
    if (length(gensB)>1) { TranMatArrayB<-sweep(TranMatArrayB,3,probByGenB[gensB],"*")} else {TranMatArrayB<-TranMatArrayB*probByGenB[gensB]}
    # TranMatArrayB<-sweep(TranMatArrayB,3,probByGenB[gensB],"*")
    TranMatArrayB<-sweep(TranMatArrayB,2,probByloc2,"*")
    TranMatArrayB.2<-sweep(TranMatArrayB,1,mrcaVec,"*")
    
    TranMatArrayA2<-sweep(TranMatArrayA,1,mrcaVec,"*")
    TranMatArrayA3<-matrix(TranMatArrayA2,nlocs,nlocs*length(gensA))
    # TranMatArrayA3<-matrix(TranMatArrayA2,nloc.munic,nlocs*length(gensA))
    TranMatArrayB2<-matrix(TranMatArrayB,nlocs*length(gensB),nlocs,byrow=T)
    # TranMatArrayB2<-matrix(TranMatArrayB,nlocs*length(gensB),nloc.munic,byrow=T)
    probAllPrs<-eigenMapMatMult(TranMatArrayB2,TranMatArrayA3)
    den<-sum(probAllPrs)
    
    
    TranMatArrayA3<-matrix(TranMatArrayA2[,dat.in2$loc1[ii],],nlocs,length(gensA))
    TranMatArrayB2<-matrix(TranMatArrayB[,dat.in2$loc2[ii],],length(gensB),nlocs,byrow=T)
    num<-sum(eigenMapMatMult(TranMatArrayB2,TranMatArrayA3))
    
    lik=log(num)-log(den)
    # if(lik<(-6))lik<-(-6)
    
    if(calcAllProbs){
      outlist<-list(all.probs, lik)
      return(outlist)
    }else{
      return(lik)
    }
  }
  
  n.per.core=ceiling(npairs/ncore)
  ntot=n.per.core*ncore
  
  lines.mat<-matrix(1:ntot,nc=ncore,nr=n.per.core)
  
  func<-function(lines){
    if(calcAllProbs==F){
      out<-rep(NaN,length(lines))
      for (iii in 1:length(lines)){out[iii]<-llIndPair(ii=lines[iii])}
      return(out)}else{
        outprobs<-array(NaN,c(length(lines),nlocs,nlocs))
        out2<-rep(NaN,length(lines))
        for (iii in 1:length(lines)){
          a<-llIndPair(lines[iii])
          outprobs[iii,,]<-a[[1]]
          out2[iii]<-a[[2]]
        }
        outlist<-list(outprobs,out2)
        return(outlist)
      }
  }
  
  if(ncore>1){
    aa<-foreach (jj=1:ncore)%dopar%func(lines=lines.mat[,jj])}else{
      aa<-list(func(lines=lines.mat[,1]))
    }
  
  aa2<-do.call("c",aa)
  LL<-sum(aa2,na.rm=T)
  out=sum(LL)
  #LL+sum(log(dnorm(extTranMatDat.tmp$pars$beta)
  
  
  print(out)
  return(out)
}
