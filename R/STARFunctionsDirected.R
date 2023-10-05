require("Rcpp")
require("RcppArmadillo")
#require("MASS")
require("matrixcalc")
require("Matrix")
require("mvtnorm")
#require("RcppEigen")
require("truncnorm")
require("MCMCpack")

# G means -----------------------------------------------------------------

.Gmean=function(A,ind,decay=NULL){
  # input: n*n SPARSE matrix A, type of mean structure ind
  # output: n*n SPARSE matrix, unless ind %in% c(1,2,9)
  
  #if (class(A) != 'dgCMatrix') as(A, 'dgCMatrix')
  
  if(ind==1 | ind== "outDeg"){
    return(matrix(rowSums(A),dim(A)[1],dim(A)[1]))
  }
  if(ind==2 | ind== "inDeg"){
    return(matrix(colSums(A),dim(A)[1],dim(A)[1],byrow=T))
  }
  if(ind==3 | ind== "stability"){
    return(as(A, "dgCMatrix"))
  }
  if(ind==4 | ind== "reciprocity"){
    return(as(t(A), "dgCMatrix"))
  }
  if(ind==5 | ind== "trans1"){
    return(as(A%*%A, "dgCMatrix"))
  }
  if(ind==6 | ind== "trans2"){
    return(as(tcrossprod(A), "dgCMatrix"))
  }
  if(ind==7 | ind== "trans3"){
    return(as(crossprod(A), "dgCMatrix"))
  }
  if(ind==8 | ind== "3cycle"){
    return(as(t(A)%*%t(A), "dgCMatrix"))
  }
  if(ind==9 | ind== "gwdeg"){
    #Snijders et al. 2006 eq. (17)
    temp0 = 1L*(A+t(A)>0)
    temp1 = rowSums(temp0)
    temp2 = matrix(temp1,nrow(A),nrow(A)) #- temp0
    temp2[which(temp0==1)] = temp2[which(temp0==1)]-1
    return(-(1-exp(-decay))*(exp(-decay*temp2)+exp(-decay*t(temp2))))
  }
  if(ind==10 | ind=="isolatesSendSym"){
    temp0 = 1L*(A+t(A)>0)
    temp1 = rowSums(temp0)
    temp2 = Matrix(0L,nrow(A),ncol(A))
    try({temp2[which(temp1==0),]=1L},silent=T)
    return(temp2)
  }
  if(ind==11 | ind=="isolatesRecSym"){
    temp0 = 1L*(A+t(A)>0)
    temp1 = rowSums(temp0)
    temp2 = Matrix(0L,nrow(A),ncol(A))
    try({temp2[,which(temp1==0)]=1L},silent=T)
    return(temp2)
  }
  if(ind==12 | ind=="isolatesSend"){
    temp1 = rowSums(A)
    temp2 = Matrix(0L,nrow(A),ncol(A))
    try({temp2[which(temp1==0),]=1L},silent=T)
    return(temp2)
  }
  if(ind==13 | ind=="isolatesRec"){
    temp1 = colSums(A)
    temp2 = Matrix(0L,nrow(A),ncol(A))
    try({temp2[,which(temp1==0)]=1L},silent=T)
    return(temp2)
  }
}



# Trace -------------------------------------------------------------------

.Tr = function(X){
  sum(diag(X))
}



# Main VB Function with random effects --------------------------------------------------------------

.STAR.RE = function(A,X,lagFunctions=lagFunctions, LL,
                    maxIter=1000,eps=1e-5,
                    sigB,sigTh,
                    #mu_sr,
                    s2R,#aRec0,bRec0,
                    aS0, bS0, aR0, bR0,
                    aOm0, BOm0,
                    Xvarnames=NULL){
  
  
  # Create/initalize objects ----------------------------------------------------------
  
  #Dimensions:
  n=dim(A[[1]])[1] 
  TT=length(X) 
  pp=length(X[[1]])
  
  if(is.null(Xvarnames)){
    if(class(lagFunctions$functions)=="character"){
      Xvarnames=c(1:pp,
                  paste(rep(lagFunctions$functions,LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }else{
      temp = c("outDeg","inDeg","stability","reciprocity","trans1","trans2",
               "trans3","3cycle","gwdeg","isolatesSendSym","isolatesRecSym",
               "isolatesSend","isolatesRec")
      Xvarnames=c(1:pp,
                  paste(rep(temp[lagFunctions$functions],LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }
  }else{
    if(class(lagFunctions$functions)=="character"){
      Xvarnames=c(Xvarnames,
                  paste(rep(lagFunctions$functions,LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }else{
      temp = c("outDeg","inDeg","stability","reciprocity","trans1","trans2",
               "trans3","3cycle","gwdeg","isolatesSendSym","isolatesRecSym",
               "isolatesSend","isolatesRec")
      Xvarnames=c(Xvarnames,
                  paste(rep(temp[lagFunctions$functions],LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }
  }
  qq=length(lagFunctions$functions)
  
  #Misc
  Imin = matrix(1,n,n)-diag(n)
  ones = matrix(1,n,1)
  crit=numeric(maxIter)
  
  #Astar in matrix and vector form
  MAt.vec = matrix(0,n*(n-1),TT)
  for(tt in 1:TT) MAt.vec[,tt] = cVecSp(A[[tt+LL]])
  Aind= which(MAt.vec!=0)
  nAind=length(Aind)
  
  
  Hr.Inv = Hs.Inv = Hr= Hs = array(0,c(n,n,TT))
  for(tt in 1:TT){
    diag(A[[tt+LL-1]])=1
    Din = diag(1/sqrt(colSums(A[[tt+LL-1]])))
    Dout = diag(1/sqrt(rowSums(A[[tt+LL-1]])))
    
    Hs[,,tt] = cRoundMatrix(as.matrix(tcrossprod(Dout%*%A[[tt+LL-1]])))
    Hr[,,tt] = cRoundMatrix(as.matrix(crossprod(A[[tt+LL-1]]%*%Din)))

    if(!is.positive.definite(Hs[,,tt])){
      Hs[,,tt] = Hs[,,tt] + diag(n)*1e-5 
    }
    if(!is.positive.definite(Hr[,,tt])){
      Hr[,,tt] = Hr[,,tt] + diag(n)*1e-5
    }
    diag(A[[tt+LL-1]])=0
  }
  
  for(tt in 1:TT){
    Hs.Inv[,,tt] <- cRoundMatrix(chol2inv(chol(Hs[,,tt])))
    Hr.Inv[,,tt] <- cRoundMatrix(chol2inv(chol(Hr[,,tt])))
  }
  
  
  #covariates
  X.vec = array(0,c(n*(n-1),pp+qq*LL,TT))
  for(tt in 1:TT){
    for(ell in 1:pp){
      X.vec[,ell,tt] = cVec(X[[tt]][[ell]])
    }
    for(LAG in 1:LL){
      for(ell in 1:qq){
        if(lagFunctions$functions[ell] %in% c(1,2,9) |
           lagFunctions$functions[ell] %in% c("outDeg","inDeg","gwdeg")){
          X.vec[,pp+ell+qq*(LAG-1),tt] = 
            cVec(.Gmean(A[[tt+LL-LAG]],lagFunctions$functions[ell],lagFunctions$decay))
        }else{
          X.vec[,pp+ell+qq*(LAG-1),tt] = 
            cVecSp(.Gmean(A[[tt+LL-LAG]],lagFunctions$functions[ell],lagFunctions$decay))
        }
      }
    }
  }
  
  #beta and theta
  SigM.Inv = matrix(0,pp+qq*LL,pp+qq*LL)
  for(tt in 1:TT){
    SigM.Inv = SigM.Inv + crossprod(X.vec[,,tt])
  }
  diag(SigM.Inv) = diag(SigM.Inv) + rep(1/c(sigB,sigTh),c(pp,qq*LL))
  SigM = cRoundMatrix(chol2inv(chol(SigM.Inv)))
  
  temp = .STAR.noRE(A=A,X=X,lagFunctions=lagFunctions,LL=LL,
                        maxIter=50,eps=eps,
                        sigB=sigB,sigTh=sigTh,
                        Xvarnames=NULL,showPB=FALSE)
  muM = c(temp$BetaMean,temp$ThetaMean)
  MAt.vec = temp$AstarMean
  
  #Variance components
  aS = aS0+0.5*n*TT
  bS = bS0
  aR = aR0+0.5*n*TT
  bR = bR0
  aOm = aOm0 + n*TT
  BOm = BOm0
  BOm.Inv = cRoundMatrix(chol2inv(chol(BOm)))
  
  #random effects
  musr = matrix(0,4*n,TT)
  Sigsr4 = Sigsr3 = Sigsr = array(0,c(4*n,4*n,TT))
  Sigsr5 = matrix(0, 4*n,4*n)
  
  Sigsr1 = Sigsr2 = matrix(0,4*n,4*n)
  for(i in 1:4){
    if(i%%2==1){
      Sigsr1[cbind(rep(n*(i-1)+1:n,2),c(1:n,2*n+1:n))]=n-1
      Sigsr2[n*(i-1)+1:n,n+1:n]=Imin
      Sigsr2[n*(i-1)+1:n,3*n+1:n]=Imin
      
    }else{
      Sigsr1[cbind(rep(n*(i-1)+1:n,2),c(n+1:n,3*n+1:n))]=n-1
      Sigsr2[n*(i-1)+1:n,1:n]=Imin
      Sigsr2[n*(i-1)+1:n,2*n+1:n]=Imin
    }
  }
  
  for(tt in 1:TT) Sigsr3[1:n,1:n,tt] = Hs.Inv[,,tt]
  for(tt in 1:TT) Sigsr4[n+1:n,n+1:n,tt] = Hr.Inv[,,tt]
  Sigsr5[2*n +1:(2*n), 2*n +1:(2*n)] = aOm*BOm.Inv%x%diag(n)
  for(tt in 1:TT){
    Sigsr[,,tt] = Sigsr1 + 
      Sigsr2 + 
      aS0/bS0*Sigsr3[,,tt] +
      aR0/bR0*Sigsr4[,,tt] +
      Sigsr5
    Sigsr[,,tt] = cRoundMatrix(chol2inv(chol(Sigsr[,,tt])))
  }
  
  Atilde = matrix(0,n,n)
  XBTh = matrix(0,n*(n-1),TT)
  musr.vec = matrix(0,n*(n-1),TT)
  for(tt in 1:TT){
    XBTh[,tt] = X.vec[,,tt]%*%muM
    Atilde = cUnvec(MAt.vec[,tt]-XBTh[,tt])
    temp1= rowSums(Atilde);temp2=colSums(Atilde)
    musr[,tt] = Sigsr[,,tt]%*%rep(c(temp1,temp2),2)
    musr.vec[,tt] = cVec(tcrossprod(musr[1:n,tt]+musr[2*n+1:n,tt], ones)) + 
                     cVec(tcrossprod(ones, musr[n+1:n,tt]+musr[3*n+1:n,tt]))
  }
  
  # Reciprocity
  SigR = s2R/(1+2*s2R)
  
  MRt.vec=RecA.vec=matrix(0,n*(n-1),TT)
  for(tt in 1:TT){
   RecA.vec[,tt] = MAt.vec[,tt] - XBTh[,tt] - musr.vec[,tt]
   MRt.vec[,tt] = SigR*cSumXtX(RecA.vec[,tt])
  }
  
  # Begin Variational EM Algorithm ------------------------------------------
  
  pb=txtProgressBar(min=1,max=maxIter,style=3)
  for(it in 1:maxIter){
    
    #Astar (Result 3)
    MAt.vec=XBTh 
    for(tt in 1:TT){
      MAt.vec[,tt] = MAt.vec[,tt] + musr.vec[,tt] + MRt.vec[,tt]
    }
    MAt.vec[Aind] = etruncnorm(a=0,mean=MAt.vec[Aind])
    MAt.vec[-Aind] = etruncnorm(b=0,mean=MAt.vec[-Aind])
    
    #beta and theta (Result 1)
    muM.Old = muM
    muM[]=0
    for(tt in 1:TT){
      muM = muM + crossprod(X.vec[,,tt], 
                            MAt.vec[,tt] - musr.vec[,tt] - MRt.vec[,tt])            
    }
    muM = SigM%*%muM
    
    for(tt in 1:TT){
      XBTh[,tt] = X.vec[,,tt]%*%muM
    }
    
    #variance components (Result 2)
    bS = bS0
    bR = bR0
    for(tt in 1:TT){
      bS = bS + 0.5*( .Tr(Sigsr[1:n,1:n,tt]%*%Hs.Inv[,,tt]) + 
                        drop( crossprod(musr[1:n,tt], Hs.Inv[,,tt])%*%musr[1:n,tt] ) )
      bR = bR + 0.5*( .Tr(Sigsr[n+1:n,n+1:n,tt]%*%Hr.Inv[,,tt]) + 
                        drop(crossprod(musr[n+1:n,tt], Hr.Inv[,,tt])%*%musr[n+1:n,tt]))
    }
    
    BOm = cBOmega1(BOm0,Sigsr[c(2*n+1:(2*n)),c(c(2*n+1:(2*n))),],musr[2*n+1:(2*n),])
    BOm.Inv = cRoundMatrix(chol2inv(chol(BOm)))
    
    #Random Effects (Result 4)
    Sigsr5 = matrix(c(0,0,0,1),2,2)%x%(aOm*BOm.Inv%x%diag(n))
    for(tt in 1:TT){
      Sigsr[,,tt] = Sigsr1 + Sigsr2 + aS/bS*Sigsr3[,,tt] +
        aR/bR*Sigsr4[,,tt] +Sigsr5
      Atilde = cUnvec(MAt.vec[,tt]-XBTh[,tt]-MRt.vec[,tt])
      temp1= rowSums(Atilde);temp2=colSums(Atilde)
      Sigsr[,,tt] = cRoundMatrix(chol2inv(chol(Sigsr[,,tt])))
      musr[,tt] = Sigsr[,,tt]%*%rep(c(temp1,temp2),2)
      musr.vec[,tt] = cVec(tcrossprod(musr[1:n,tt]+musr[2*n+1:n,tt], ones)) + 
                        cVec(tcrossprod(ones, musr[n+1:n,tt]+musr[3*n+1:n,tt]))
    }
    
    
    # Reciprocity (Result 5)
    for(tt in 1:TT){
      RecA.vec[,tt] = MAt.vec[,tt] - XBTh[,tt] - musr.vec[,tt]
      MRt.vec[,tt] = SigR*cSumXtX(RecA.vec[,tt])
    }
    
    #Check if we can break out of the loop
    crit[it] = max(na.omit(abs(muM-muM.Old)/abs(muM.Old)))
    if(crit[it]<eps){
      crit = crit[1:it]
      break
    }
    
    setTxtProgressBar(pb,it)
  }
  
  names(muM)=Xvarnames
  rownames(SigM) = colnames(SigM) = Xvarnames
  
  retObj = list(BetaMean=muM[1:pp],
                ThetaMean=muM[pp+1:(qq*LL)],
                BetaThetaCov=SigM,
                sShape=aS,
                sScale=bS,
                rShape=aR,
                rScale=bR,
                OmegaDOF=aOm,
                OmegaScaleMat=BOm,
                s2Means=musr[1:n,],
                s1Means=musr[2*n+1:n,],
                r2Means=musr[n+1:n,],
                r1Means=musr[3*n+1:n,],
                srCov=Sigsr,
                RecMean=MRt.vec,
                RecVar=SigR,
                numIter=it,
                convergence=crit,
                RE=TRUE,method="VB")
  class(retObj) = "STAR"
  return(retObj)
  
}



# Main VB Function without random effects --------------------------------------------------------------

.STAR.noRE = function(A,X,lagFunctions=lagFunctions,LL,
                      maxIter=1000,eps=1e-5,
                      sigB,sigTh,
                      Xvarnames=NULL,showPB=TRUE){
  
  
  # Create/initalize objects ----------------------------------------------------------
  
  #Dimensions:
  n = dim(A[[1]])[1]
  TT = length(X) 
  pp = length(X[[1]]) 
  
  if(is.null(Xvarnames)){
    if(class(lagFunctions$functions)=="character"){
      Xvarnames=c(1:pp,
                  paste(rep(lagFunctions$functions,LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }else{
      temp = c("outDeg","inDeg","stability","reciprocity","trans1","trans2",
               "trans3","3cycle","gwdeg","isolatesSendSym","isolatesRecSym",
               "isolatesSend","isolatesRec")
      Xvarnames=c(1:pp,
                  paste(rep(temp[lagFunctions$functions],LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }
  }else{
    if(class(lagFunctions$functions)=="character"){
      Xvarnames=c(Xvarnames,
                  paste(rep(lagFunctions$functions,LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }else{
      temp = c("outDeg","inDeg","stability","reciprocity","trans1","trans2",
               "trans3","3cycle","gwdeg","isolatesSendSym","isolatesRecSym",
               "isolatesSend","isolatesRec")
      Xvarnames=c(Xvarnames,
                  paste(rep(temp[lagFunctions$functions],LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }
  }
  qq=length(lagFunctions$functions)
  
  #Misc
  Imin = matrix(1,n,n)-diag(n)
  ones = matrix(1,n,1)
  crit=numeric(maxIter)
  
  #Astar in matrix and vector form
  MAt.vec = matrix(0,n*(n-1),TT)
  for(tt in 1:TT) MAt.vec[,tt] = cVecSp(A[[tt+LL]])
  Aind= which(MAt.vec!=0)
  nAind=length(Aind)
  MAt.vec[Aind] = etruncnorm(a=0)
  MAt.vec[-Aind] = etruncnorm(b=0) 
  
  
  
  #covariates
  X.vec = array(0,c(n*(n-1),pp+qq*LL,TT))
  for(tt in 1:TT){
    for(ell in 1:pp){
      X.vec[,ell,tt] = cVec(X[[tt]][[ell]])
    }
    for(LAG in 1:LL){
      for(ell in 1:qq){
        if(lagFunctions$functions[ell] %in% c(1,2,9) |
           lagFunctions$functions[ell] %in% c("outDeg","inDeg","gwdeg")){
          X.vec[,pp+ell+qq*(LAG-1),tt] = 
            cVec(.Gmean(A[[tt+LL-LAG]],lagFunctions$functions[ell],lagFunctions$decay))
        }else{
          X.vec[,pp+ell+qq*(LAG-1),tt] = 
            cVecSp(.Gmean(A[[tt+LL-LAG]],lagFunctions$functions[ell],lagFunctions$decay))
        }
      }
    }
  }
  
  #beta and theta
  SigM.Inv = matrix(0,pp+qq*LL,pp+qq*LL)
  for(tt in 1:TT){
    SigM.Inv = SigM.Inv + crossprod(X.vec[,,tt])
  }
  diag(SigM.Inv) = diag(SigM.Inv) + rep(1/c(sigB,sigTh),c(pp,qq*LL))
  SigM = cRoundMatrix(chol2inv(chol(SigM.Inv)))
  
  muM = numeric(pp+qq*LL)
  for(tt in 1:TT){
    muM= muM + crossprod(X.vec[,,tt], MAt.vec[,tt]) 
  }
  muM= SigM%*%muM 
  
  XBTh = matrix(0,n*(n-1),TT)
  for(tt in 1:TT){
    XBTh[,tt] = X.vec[,,tt]%*%muM
  }
  

  
  # Begin Variational EM Algorithm ------------------------------------------
  
  if(showPB)pb=txtProgressBar(min=1,max=maxIter,style=3)
  for(it in 1:maxIter){
    
    #Astar
    MAt.vec=XBTh 
    MAt.vec[Aind] = etruncnorm(a=0,mean=MAt.vec[Aind])
    MAt.vec[-Aind] = etruncnorm(b=0,mean=MAt.vec[-Aind])
    
    #beta and theta
    muM.Old = muM
    muM[]=0
    for(tt in 1:TT){
      muM = muM + crossprod(X.vec[,,tt], MAt.vec[,tt]) 
    }
    muM = SigM%*%muM
    
    for(tt in 1:TT){
      XBTh[,tt] = X.vec[,,tt]%*%muM
    }
    
    #Check if we can break out of the loop
    crit[it] = max(na.omit(abs(muM-muM.Old)/abs(muM.Old)))
    if(crit[it]<eps){
      crit = crit[1:it]
      break
    }
    
    if(showPB) setTxtProgressBar(pb,it)
  }
  
  names(muM)=Xvarnames
  rownames(SigM) = colnames(SigM) = Xvarnames
  
  retObj = 
    list(AstarMean=MAt.vec,
       BetaMean=muM[1:pp],
       ThetaMean=muM[pp+1:(qq*LL)],
       BetaThetaCov=SigM,
       numIter=it,
       convergence=crit,
       RE=FALSE,method="VB")
  class(retObj) = "STAR"
  return(retObj)
}


# Main Gibbs Function with random effects --------------------------------------------------------------

.STAR.RE.Gibbs = function(A,X,lagFunctions=NULL,LL,
                    nSims=1000,
                    sigB,sigTh,
                    #mu_sr,
                    aRec0,bRec0,
                    aS0, bS0, aR0, bR0,
                    aOm0, BOm0,
                    Xvarnames=NULL){
  
  
  # Create/initalize objects ----------------------------------------------------------
  
  nSims = nSims + 1
  
  #Dimensions:
  n=dim(A[[1]])[1] 
  TT=length(X) 
  pp=length(X[[1]])
  
  if(is.null(Xvarnames)){
    if(class(lagFunctions$functions)=="character"){
      Xvarnames=c(1:pp,
                  paste(rep(lagFunctions$functions,LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }else{
      temp = c("outDeg","inDeg","stability","reciprocity","trans1","trans2",
               "trans3","3cycle","gwdeg","isolatesSendSym","isolatesRecSym",
               "isolatesSend","isolatesRec")
      Xvarnames=c(1:pp,
                  paste(rep(temp[lagFunctions$functions],LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }
  }else{
    if(class(lagFunctions$functions)=="character"){
      Xvarnames=c(Xvarnames,
                  paste(rep(lagFunctions$functions,LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }else{
      temp = c("outDeg","inDeg","stability","reciprocity","trans1","trans2",
               "trans3","3cycle","gwdeg","isolatesSendSym","isolatesRecSym",
               "isolatesSend","isolatesRec")
      Xvarnames=c(Xvarnames,
                  paste(rep(temp[lagFunctions$functions],LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }
  }
  qq=length(lagFunctions$functions)
  
  #MCMC objects
  s2R = tauS1 = tauR1 = numeric(nSims)
  betaThetas = matrix(0.0,pp+qq*LL,nSims)
  rEff2 = rEff1 = sEff2 = sEff1 = array(0.0,c(n,TT,nSims))
  Omegas = array(0.0,c(2,2,nSims))
  
  tauS1[1] = bS0/(aS0-1)
  tauR1[1] = bR0/(aR0-1)
  Omegas[,,1] = BOm0/(aOm0-3)
  
  s2R[1] = (aRec0/bRec0)/(1+2*aRec0/bRec0)
  
  #Misc
  Imin = matrix(1,n,n)-diag(n)
  ones = matrix(1,n,1)
  
  #Astar in matrix and vector form
  MAt.vec = matrix(0,n*(n-1),TT)
  for(tt in 1:TT) MAt.vec[,tt] = cVecSp(A[[tt+LL]])
  Aind= which(MAt.vec!=0)
  nAind=length(Aind)
  
  
  Hr.Inv = Hs.Inv = Hr= Hs = array(0,c(n,n,TT))
  for(tt in 1:TT){
    diag(A[[tt+LL-1]])=1
    Din = diag(1/sqrt(colSums(A[[tt+LL-1]])))
    Dout = diag(1/sqrt(rowSums(A[[tt+LL-1]])))
    
    Hs[,,tt] = cRoundMatrix(as.matrix(tcrossprod(Dout%*%A[[tt+LL-1]])))
    Hr[,,tt] = cRoundMatrix(as.matrix(crossprod(A[[tt+LL-1]]%*%Din)))
    
    if(!is.positive.definite(Hs[,,tt])){
      Hs[,,tt] = Hs[,,tt] + diag(n)*1e-5 
    }
    if(!is.positive.definite(Hr[,,tt])){
      Hr[,,tt] = Hr[,,tt] + diag(n)*1e-5
    }
    diag(A[[tt+LL-1]])=0
  }
  
  for(tt in 1:TT){
    Hs.Inv[,,tt] <- cRoundMatrix(chol2inv(chol(Hs[,,tt])))
    Hr.Inv[,,tt] <- cRoundMatrix(chol2inv(chol(Hr[,,tt])))
  }
  
  
  #covariates
  X.vec = array(0,c(n*(n-1),pp+qq*LL,TT))
  for(tt in 1:TT){
    for(ell in 1:pp){
      X.vec[,ell,tt] = cVec(X[[tt]][[ell]])
    }
    for(LAG in 1:LL){
      for(ell in 1:qq){
        if(lagFunctions$functions[ell] %in% c(1,2,9) |
             lagFunctions$functions[ell] %in% c("outDeg","inDeg","gwdeg")){
          X.vec[,pp+ell+qq*(LAG-1),tt] = 
            cVec(.Gmean(A[[tt+LL-LAG]],lagFunctions$functions[ell],lagFunctions$decay))
        }else{
          X.vec[,pp+ell+qq*(LAG-1),tt] = 
            cVecSp(.Gmean(A[[tt+LL-LAG]],lagFunctions$functions[ell],lagFunctions$decay))
        }
      }
    }
  }
  
  #beta and theta
  SigM = matrix(0,pp+qq*LL,pp+qq*LL)
  for(tt in 1:TT){
    SigM = SigM + crossprod(X.vec[,,tt])
  }
  diag(SigM) = diag(SigM) + rep(1/c(sigB,sigTh),c(pp,qq*LL))
  SigM = cRoundMatrix(chol2inv(chol(SigM)))
  eg = eigen(SigM)
  SigMsqrt = tcrossprod(eg$vectors%*%diag(sqrt(eg$values)),eg$vectors)
  
  
  temp = .STAR.noRE(A=A,X=X,lagFunctions=lagFunctions,LL=LL,
                    maxIter=50,eps=1e-5,
                    sigB=sigB,sigTh=sigTh,
                    Xvarnames=NULL,showPB=FALSE)
  betaThetas[,1] = c(temp$BetaMean,temp$ThetaMean)
  MAt.vec = temp$AstarMean
  
  
  #Variance components
  aS = aS0+0.5*n*TT
  bS = bS0
  aR = aR0+0.5*n*TT
  bR = bR0
  aOm = aOm0 + n*TT
  BOm = BOm0
  
  #random effects
  musr = numeric(4*n)
  Sigsr4 = Sigsr3 = list()
  Sigsr = Matrix(0,4*n,4*n)
  for(tt in 1:TT){
    Sigsr4[[tt]] = Sigsr3[[tt]] = Matrix(0,4*n,4*n)
  }
  Sigsr5 = Matrix(0, 4*n,4*n)
  
  Sigsr1 = Sigsr2 = Matrix(0,4*n,4*n)
  for(i in 1:4){
    if(i%%2==1){
      Sigsr1[cbind(rep(n*(i-1)+1:n,2),c(1:n,2*n+1:n))]=n-1
      Sigsr2[n*(i-1)+1:n,n+1:n]=Imin
      Sigsr2[n*(i-1)+1:n,3*n+1:n]=Imin
      
    }else{
      Sigsr1[cbind(rep(n*(i-1)+1:n,2),c(n+1:n,3*n+1:n))]=n-1
      Sigsr2[n*(i-1)+1:n,1:n]=Imin
      Sigsr2[n*(i-1)+1:n,2*n+1:n]=Imin
    }
  }
  
  for(tt in 1:TT) Sigsr3[[tt]][1:n,1:n] = Hs.Inv[,,tt]
  for(tt in 1:TT) Sigsr4[[tt]][n+1:n,n+1:n] = Hr.Inv[,,tt]
  Sigsr5[2*n +1:(2*n), 2*n +1:(2*n)] = aOm*solve(BOm)%x%diag(n)
  
  Atilde = matrix(0,n,n)
  XBTh = matrix(0,n*(n-1),TT)
  sr.vec = matrix(0,n*(n-1),TT)
  for(tt in 1:TT){
   XBTh[,tt] = X.vec[,,tt]%*%betaThetas[,1]
   Atilde = cUnvec(MAt.vec[,tt]-XBTh[,tt])
   temp1= rowSums(Atilde);temp2=colSums(Atilde)
   Sigsr = Sigsr1 + 
     Sigsr2 + 
     aS0/bS0*Sigsr3[[tt]] +
     aR0/bR0*Sigsr4[[tt]] +
     Sigsr5
   Sigsr = cRoundMatrix(chol2inv(chol(as.matrix(Sigsr))))
   musr = Sigsr%*%rep(c(temp1,temp2),2)
   sEff1[,tt,1] <- musr[1:n]
   rEff1[,tt,1] <- musr[n+1:n]
   sEff2[,tt,1] <- musr[2*n+1:n]
   rEff2[,tt,1] <- musr[3*n+1:n]
   sr.vec[,tt] = cVec(tcrossprod(sEff1[,tt,1]+sEff2[,tt,1],ones)+
                       tcrossprod(ones,rEff1[,tt,1]+rEff2[,tt,1]))
  }
  
  # Reciprocity
  aRec = aRec0+TT*n*(n-1)/4
  bRec = bRec0
  SigR = bRec0/(aRec0+2*bRec0)
  
  MRt.vec=RecA.vec=matrix(0,n*(n-1),TT)
  for(tt in 1:TT){
   RecA.vec[,tt] = MAt.vec[,tt] - XBTh[,tt] - sr.vec[,tt]
   MRt.vec[,tt] = SigR/(1+2*SigR)*cSumXtX(RecA.vec[,tt])
  }
  
  # Begin Gibbs Algorithm ------------------------------------------
  
  pb=txtProgressBar(min=2,max=nSims,style=3)
  for(it in 2:nSims){
    
    #beta and theta (Result 1)
    for(tt in 1:TT){
      betaThetas[,it] = betaThetas[,it] + 
        crossprod(X.vec[,,tt],MAt.vec[,tt] - sr.vec[,tt] - MRt.vec[,tt])          
    }
    betaThetas[,it] = SigM%*%betaThetas[,it]
    betaThetas[,it] = betaThetas[,it] + SigMsqrt%*%rnorm(pp+qq*LL)
    
    for(tt in 1:TT){
     XBTh[,tt] = X.vec[,,tt]%*%betaThetas[,it]
    }
    
    
    #Astar (Result 3)
    MAt.vec[]=0
    for(tt in 1:TT){
      MAt.vec[,tt] = XBTh[,tt] + sr.vec[,tt] + MRt.vec[,tt]
    }
    MAt.vec[Aind] = rtruncnorm(length(Aind),a=0,mean=MAt.vec[Aind])
    MAt.vec[-Aind] = rtruncnorm(n*(n-1)*TT-length(Aind),b=0,mean=MAt.vec[-Aind])

    
    #Random Effects (Result 4)
    Sigsr5[2*n+1:(2*n),2*n+1:(2*n)] = chol2inv(chol(Omegas[,,it-1]))%x%diag(n)
    for(tt in 1:TT){
      Sigsr = Sigsr1 + Sigsr2 + Sigsr3[[tt]]/tauS1[it-1] +
        Sigsr4[[tt]]/tauR1[it-1] +Sigsr5
      Atilde = cUnvec(MAt.vec[,tt]-XBTh[,tt]-MRt.vec[,tt])
      temp1= rowSums(Atilde);temp2=colSums(Atilde)
      Sigsr = cRoundMatrix(chol2inv(chol(as.matrix(Sigsr))))
      musr = Sigsr%*%rep(c(temp1,temp2),2)
      musr = drop(rmvnorm(1,mean=musr,sigma=Sigsr))
      sEff1[,tt,it] <- musr[1:n]
      rEff1[,tt,it] <- musr[n+1:n]
      sEff2[,tt,it] <- musr[2*n+1:n]
      rEff2[,tt,it] <- musr[3*n+1:n]
      sr.vec[,tt] = cVec(tcrossprod(sEff1[,tt,it]+sEff2[,tt,it],ones)+
                           tcrossprod(ones,rEff1[,tt,it]+rEff2[,tt,it]))
    }
    
    
    #variance components (Result 2)
    bS = bS0
    bR = bR0
    for(tt in 1:TT){
      bS = bS + 0.5*drop(crossprod(sEff1[,tt,it], Hs.Inv[,,tt]%*%sEff1[,tt,it]))
      bR = bR + 0.5*drop(crossprod(rEff1[,tt,it], Hr.Inv[,,tt]%*%rEff1[,tt,it]))
    }
    
    tauS1[it] <- rinvgamma(1,aS,bS)
    tauR1[it] <- rinvgamma(1,aR,bR)
    
    BOm = cBOmegaGibbs(BOm0,sEff2[,,it],rEff2[,,it])
    Omegas[,,it] <- riwish(aOm,BOm)
    
    
    
    #Recip Variance (Result 6)
    bRec=bRec0
    for(tt in 1:TT){
      bRec = bRec + 0.5*cUpperMRt(MRt.vec[,tt])
    }
    
    s2R[it]=rinvgamma(1,aRec,bRec)
    
    # Reciprocity (Result 5)
    SigR = s2R[it]/(1+2*s2R[it])
    for(tt in 1:TT){
      RecA.vec[,tt] = MAt.vec[,tt] - XBTh[,tt] - sr.vec[,tt]
      MRt.vec[,tt] = cSumXtXAndDraw(RecA.vec[,tt],SigR)
    }
    
    
    setTxtProgressBar(pb,it)
  }
  
  rownames(betaThetas) = Xvarnames
  
  retObj = 
    list(Beta=betaThetas[1:pp,-1],
         Theta=betaThetas[pp+1:(qq*LL),-1],
         Omega=Omegas[,,-1],
         tauS=tauS1[-1],
         tauR=tauR1[-1],
         s1=sEff2[,,-1],
         s2=sEff1[,,-1],
         r1=rEff2[,,-1],
         r2=rEff1[,,-1],
         s2R = s2R[-1],
         RE=TRUE,method="Gibbs")
  class(retObj) = "STAR"
  return(retObj)
  
}


# Extract Posterior Means from Gibbs sampler ------------------------------

extractPostMeans = function(STARObj,burnin=0,thinning=1){
  if(STARObj$method!='Gibbs')stop("Method must be 'Gibbs'")
  nSims = ncol(STARObj$Beta)
  if(burnin>nSims)stop('burn-in cannot be larger than the number of simulations')
  N = seq(from=burnin+1,to=nSims,by=thinning)
  ret = list()
  ret$BetaMean = rowMeans(STARObj$Beta[,N])
  ret$ThetaMean = rowMeans(STARObj$Theta[,N])
  if(STARObj$RE){
    ret$OmegaMean = matrix(0,2,2)
    for(i in 1:2)for(j in i:2) ret$OmegaMean[i,j] = mean(STARObj$Omega[i,j,N])
    ret$OmegaMean[2,1] = ret$OmegaMean[1,2]
    ret$tauS = mean(STARObj$tauS[N])
    ret$tauR = mean(STARObj$tauR[N])
    ret$s2R = mean(STARObj$s2R) 
  }
  
  STARObj$PostMean = ret
  return(STARObj)
}


# summary -----------------------------------------------------------------

summary.STAR <- function(STARObj,alpha=0.05,burnin=0,thinning=1,...){
  sObj = list(); class(sObj) = 'summary.STAR'
  
  if(STARObj$method=='VB'){
    
    sObj$method = 'VB'
    sObj$RE = STARObj$RE
    pp = NROW(STARObj$BetaMean)
    qq = NROW(STARObj$ThetaMean)
    tempSD = sqrt(diag(STARObj$BetaThetaCov))
    sObj$coefficients = 
      data.frame(#Variable=c(names(STARObj$BetaMean),names(STARObj$ThetaMean)),
                 Estimate=c(STARObj$BetaMean,STARObj$ThetaMean),
                 sigma=tempSD,
                 c(STARObj$BetaMean,STARObj$ThetaMean)+qnorm(alpha/2)*tempSD,
                 c(STARObj$BetaMean,STARObj$ThetaMean)+qnorm(1-alpha/2)*tempSD)
    colnames(sObj$coefficients)[3:4] = paste(100*c(alpha/2,1-alpha/2),'%',sep='')
    
    print(sObj,...)
    
  }else{
    
    sObj$method = 'Gibbs'
    sObj$RE = STARObj$RE
    pp = nrow(STARObj$Beta)
    qq = nrow(STARObj$Theta)
    nSims = ncol(STARObj$Beta)
    N = seq(from = burnin + 1, to = nSims, by = thinning)
    
    if(is.null(STARObj$PostMean)){
      warning("No posterior mean previously computed.  Reduce redundant computation by running 'extractPostMeans()' first.")
      STARObj = extractPostMeans(STARObj,burnin=burnin,thinning=thinning)
    }
    
    tempQuant1 = apply(STARObj$Beta[,N],1,quantile,probs=c(alpha/2,1-alpha/2))
    tempQuant2 = apply(STARObj$Theta[,N],1,quantile,probs=c(alpha/2,1-alpha/2))
    sObj$coefficients = 
      data.frame(Estimate = c(STARObj$PostMean$BetaMean,STARObj$PostMean$ThetaMean),
                 c(tempQuant1[1,],tempQuant2[1,]),
                 c(tempQuant1[2,],tempQuant2[2,]))
    colnames(sObj$coefficients)[2:3] = paste(100*c(alpha/2,1-alpha/2),'%',sep='')
    
    print(sObj,...)
    
  }
  
  invisible(sObj)
}


print.summary.STAR = function(sObj,digits = max(3, getOption("digits")-3),...){
  if(sObj$method=='VB'){
    
    cat('Method: Variational Bayes \n')
    
    if(sObj$RE){cat('Simultaneous dependence: Included \n')
    }else{cat('Simultaneous dependence: Ignored \n')}
    
    cat('Coefficients: \n')
    print(round(sObj$coef,digits))
    
  }else{
    
    cat('Method: Gibbs sampler \n')
    
    if(sObj$RE){cat('Simultaneous dependence: Included \n')
    }else{cat('Simultaneous dependence: Ignored \n')}
    
    cat('Coefficients: \n')
    print(round(sObj$coef,digits))
    
  }
}

# Main Gibbs Function without random effects --------------------------------------------------------------

.STAR.noRE.Gibbs = function(A,X,lagFunctions=lagFunctions,LL,
                            nSims,sigB,sigTh,
                            Xvarnames=NULL){
  
  # Create/initalize objects ----------------------------------------------------------
  
  nSims = nSims + 1
  
  #Dimensions:
  n = dim(A[[1]])[1]
  TT = length(X) 
  pp = length(X[[1]]) 
  
  if(is.null(Xvarnames)){
    if(class(lagFunctions$functions)=="character"){
      Xvarnames=c(1:pp,
                  paste(rep(lagFunctions$functions,LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }else{
      temp = c("outDeg","inDeg","stability","reciprocity","trans1","trans2",
               "trans3","3cycle","gwdeg","isolatesSendSym","isolatesRecSym",
               "isolatesSend","isolatesRec")
      Xvarnames=c(1:pp,
                  paste(rep(temp[lagFunctions$functions],LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }
  }else{
    if(class(lagFunctions$functions)=="character"){
      Xvarnames=c(Xvarnames,
                  paste(rep(lagFunctions$functions,LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }else{
      temp = c("outDeg","inDeg","stability","reciprocity","trans1","trans2",
               "trans3","3cycle","gwdeg","isolatesSendSym","isolatesRecSym",
               "isolatesSend","isolatesRec")
      Xvarnames=c(Xvarnames,
                  paste(rep(temp[lagFunctions$functions],LL),"_t-",
                        rep(1:LL,each=length(lagFunctions$functions)),sep=""))
    }
  }
  qq=length(lagFunctions$functions)
  
  #MCMC objects
  betaThetas = matrix(0.0,pp+qq*LL,nSims)
  
  #Misc
  Imin = matrix(1,n,n)-diag(n)
  ones = matrix(1,n,1)
  
  #Astar in matrix and vector form
  MAt.vec = matrix(0,n*(n-1),TT)
  for(tt in 1:TT) MAt.vec[,tt] = cVecSp(A[[tt+LL]])
  Aind= which(MAt.vec!=0)
  nAind=length(Aind)
  MAt.vec[Aind] = etruncnorm(a=0)
  MAt.vec[-Aind] = etruncnorm(b=0) 
  
  
  
  #covariates
  X.vec = array(0,c(n*(n-1),pp+qq*LL,TT))
  for(tt in 1:TT){
    for(ell in 1:pp){
      X.vec[,ell,tt] = cVec(X[[tt]][[ell]])
    }
    for(LAG in 1:LL){
      for(ell in 1:qq){
        if(lagFunctions$functions[ell] %in% c(1,2,9) |
           lagFunctions$functions[ell] %in% c("outDeg","inDeg","gwdeg")){
          X.vec[,pp+ell+qq*(LAG-1),tt] = 
            cVec(.Gmean(A[[tt+LL-LAG]],lagFunctions$functions[ell],lagFunctions$decay))
        }else{
          X.vec[,pp+ell+qq*(LAG-1),tt] = 
            cVecSp(.Gmean(A[[tt+LL-LAG]],lagFunctions$functions[ell],lagFunctions$decay))
        }
      }
    }
  }
  
  #beta and theta
  SigM = matrix(0,pp+qq*LL,pp+qq*LL)
  for(tt in 1:TT){
    SigM = SigM + crossprod(X.vec[,,tt])
  }
  diag(SigM) = diag(SigM) + rep(1/c(sigB,sigTh),c(pp,qq*LL))
  SigM = cRoundMatrix(chol2inv(chol(SigM)))
  
  XBTh = matrix(0,n*(n-1),TT)
  
  
  # Begin Gibbs Algorithm ------------------------------------------
  
  pb=txtProgressBar(min=2,max=nSims,style=3)
  for(it in 2:nSims){
    
    #beta and theta (Result 1)
    for(tt in 1:TT){
      betaThetas[,it] = betaThetas[,it] + 
                        crossprod(X.vec[,,tt],MAt.vec[,tt])
    }
    betaThetas[,it] = SigM%*%betaThetas[,it]
    betaThetas[,it] = drop(rmvnorm(1,mean=betaThetas[,it],sigma=SigM))
    
    for(tt in 1:TT){
      XBTh[,tt] = X.vec[,,tt]%*%betaThetas[,it]
    }
    
    #Astar (Result 3)
    MAt.vec[]=0
    for(tt in 1:TT){
      MAt.vec[,tt] = XBTh[,tt]
    }
    MAt.vec[Aind] = rtruncnorm(length(Aind),a=0,mean=MAt.vec[Aind])
    MAt.vec[-Aind] = rtruncnorm(n*(n-1)*TT-length(Aind),b=0,mean=MAt.vec[-Aind])
  
  
    setTxtProgressBar(pb,it)
  }
  
  rownames(betaThetas)=Xvarnames
  
  retObj = list(Beta=betaThetas[1:pp,],
              Theta=betaThetas[pp+1:(qq*LL),],
              RE=FALSE,method="Gibbs")
  class(retObj) = "STAR"
  return(retObj)
}


# STAR Main Function ------------------------------------------------------

STARnet = function(A,X,
                lagFunctions=list(functions = 1:8,decay=NULL),
                RandEff=TRUE,
                method=c("VB","Gibbs"),
                maxIter=1000,eps=1e-5,nSims=1000,
                hyperparms=NULL,
                Xvarnames=NULL){
  
  if(length(method)>1)method=method[1]
  
  # Set hyperparameters
  if(is.null(hyperparms$mu_s)){
    mu_s=0.5
  }else{mu_s=hyperparms$mu_s}
  if(is.null(hyperparms$mu_r)){
    mu_r=0.5
  }else{mu_r=hyperparms$mu_r}
  if(is.null(hyperparms$delta)){
    delta=0.05
  }else{delta=hyperparms$delta}
  if(is.null(hyperparms$mu_sr)){
    mu_sr=matrix(c(1,0.5,0.5,1),2,2)
  }else{mu_sr=hyperparms$mu_sr}
  if(is.null(hyperparms$sigB)){
    sigB=100
  }else{sigB=hyperparms$sigB}
  if(is.null(hyperparms$sigTh)){
    sigTh=100
  }else{sigTh=hyperparms$sigTh}
  if(is.null(hyperparms$bRec0)){
    bRec0=1+delta
  }else{bRec0=hyperparms$bRec0}
  if(is.null(hyperparms$aRec0)){
    aRec0=2+delta
  }else{aRec0=hyperparms$aRec0}
  if(is.null(hyperparms$s2R)){
    s2R=3
  }else{s2R=hyperparms$s2R}
  aS0=2+delta; bS0=mu_s*(aS0-1)
  aR0=2+delta; bR0=mu_r*(aR0-1)
  LL=length(A)-length(X)
  aOm0=5+delta; 
  BOm0=mu_sr*(aOm0-3)
  rm(mu_s,mu_r,mu_sr,delta)
  
  # Test that A is list of sparse matrices
  test <- array(NA, length(A))
  for (tt in 1:length(A)) test[tt] <- (class(A[[tt]])=="dgCMatrix")
  if (!all(test)){
    for (tt in 1:length(A)) A[[tt]] <- as(as.matrix(A[[tt]]), "dgCMatrix")
  }
  rm(test)
  
  # Diagonals for A and X should be 0's
  for (tt in 1:length(A)){
    if ( !all(diag(A[[tt]]) == 0) ){
      diag(A[[tt]]) <- 0
      warning("Diagonals of A_t should be zeros.  Continuing with diagonals coerced to zero.")
    }
    A[[tt]] <- as(as.matrix(A[[tt]]), "dgCMatrix")
  }
  for (tt in 1:length(X)){
    for (ell in 1:length(X[[1]])){
      if ( !all(diag(X[[tt]][[ell]]) == 0) ){
        diag(X[[tt]][[ell]]) <- 0
        warning("Diagonals of X_{tl} should be zeros.  Continuing with diagonals coerced to zero.")
      }
    }
  }
  
  if (method == "VB"){
    if(RandEff){
       return(.STAR.RE(A=A,X=X,lagFunctions=lagFunctions,LL=LL,
                        maxIter=maxIter,eps=eps,
                        sigB=sigB, sigTh=sigTh,
                        s2R,
                        aS0=aS0,bS0=bS0,
                        aR0=aR0, bR0=bR0,
                        aOm0=aOm0, BOm0=BOm0,
                        Xvarnames=Xvarnames))
    }else{
      return(.STAR.noRE(A=A,X=X,lagFunctions=lagFunctions,LL=LL,
                        maxIter=maxIter,eps=eps,
                        sigB=sigB,
                        sigTh=sigTh,
                        Xvarnames=Xvarnames))
    }
  }else{
    if(RandEff){
        return(.STAR.RE.Gibbs(A=A,X=X,lagFunctions=lagFunctions,LL=LL,
                              nSims=nSims,
                              sigB=sigB,sigTh=sigTh,
                              aRec0=aRec0,bRec0=bRec0,
                              aS0=aS0, bS0=bS0, aR0=aR0, bR0=bR0,
                              aOm0=aOm0, BOm0=BOm0,
                              Xvarnames=Xvarnames))
    }else{
      return(.STAR.noRE.Gibbs(A=A,X=X,lagFunctions=lagFunctions,LL=LL,
                              nSims=nSims,
                              sigB=sigB,
                              sigTh=sigTh,
                              Xvarnames=Xvarnames))
    }
  }
}


# Prob mass of RE's in sphere ---------------------------------------------

.STAR.RETest = function(STARObj,burnin=0,N=5e3,propVar=0.5,Titles=TRUE,savePlots=FALSE,
                        prefix=""){
  if(STARObj$method == "VB"){
    s2R=var(c(STARObj$RecMean))
    Means = rbind(STARObj$s1Means,STARObj$r1Means,STARObj$s2Means,STARObj$r2Means)
    n=dim(STARObj$s1Means)[1]
    TT=dim(STARObj$s1Means)[2]
    samp = matrix(0,N,TT)
    VV = propVar*(1+s2R)/(4*(1-propVar))
    for(tt in 1:TT){
      temp = rmvnorm(N,mean=Means[,tt],sigma=STARObj$srCov[,,tt])
      samp[,tt] = apply(temp,1,function(x) sqrt(crossprod(x)))
    }
    
    temp = uniroot(f=function(x){ pgamma(x^2,shape=2*n,scale=2*max(VV))-0.99},interval=c(0,100))$root
    XL = c(0,max(quantile(c(samp),1-1/N),temp))
    if(savePlots){
      jpeg(filename=paste(prefix,"REEvidence",".jpg",sep=""),
           height=800,width=800)}
    par(mar=c(5,6,1,1))
    if(Titles){
      plot(1,type="n",ylim=c(0,1),xlim=XL,xlab=expression(epsilon),
           ylab=expression("P"(REs%in%B[epsilon])),cex=2,cex.axis=2,cex.lab=2)
      for(tt in 1:TT)lines(ecdf(samp[,tt]),lwd=2)
    }else{
      plot(1,type="n",ylim=c(0,1),xlim=XL,cex=2,cex.axis=2,xlab="",ylab="")
      for(tt in 1:TT)lines(ecdf(samp[,tt]),lwd=2)
    }
    radius = seq(XL[1],XL[2],length.out=1000)
    for(vv in VV){
      lines(pgamma(radius^2,shape=2*n,scale=2*vv)~radius,lty=2,col="blue",lwd=2)
    }
    if(savePlots){dev.off()}
  }
  
  if(STARObj$method=="Gibbs"){

    n = nrow(STARObj$s1)
    TT = ncol(STARObj$s1)
    burnin = burnin +1
    MM = burnin:dim(STARObj$s1)[3]
    samp = matrix(0,length(MM),TT)
    s2R = mean(STARObj$s2R[MM])
    VV = propVar*(1+s2R)/(4*(1-propVar))
    
    for(tt in 1:TT){
      temp = rbind(STARObj$s1[,tt,MM],STARObj$r1[,tt,MM],
                   STARObj$s2[,tt,MM],STARObj$r2[,tt,MM])
      samp[,tt] = apply(temp,2,function(x) sqrt(crossprod(x)))
    }
    
    temp = uniroot(f=function(x){ pgamma(x^2,shape=2*n,scale=2*max(VV))-0.99},interval=c(0,100))$root
    XL = c(0,max(quantile(c(samp),1-1/length(MM)),temp))
    if(savePlots){
      jpeg(filename=paste(prefix,"REEvidence",".jpg",sep=""),
           height=800,width=800)}
    par(mar=c(5,6,1,1))
    if(Titles){
      plot(1,type="n",ylim=c(0,1),xlim=XL,xlab=expression(epsilon),
           ylab=expression("P"(REs%in%B[epsilon])),cex=2,cex.axis=2,cex.lab=2)
    }else{
      plot(1,type="n",ylim=c(0,1),xlim=XL,cex=2,cex.axis=2,xlab="",ylab="")
    }
    for(tt in 1:TT)lines(ecdf(samp[,tt]),lwd=2)
    radius = seq(XL[1],XL[2],length.out=1000)
    for(vv in VV){
      lines(pgamma(radius^2,shape=2*n,scale=2*vv)~radius,lty=2,col="blue",lwd=2)
    }
    if(savePlots){dev.off()}
  }
  
}


# Plot STAR object --------------------------------------------------------

.normCurvePlot = function(Mu,Sig,alpha=0.05){
  XLIM=Mu+c(-1,1)*4*Sig
  curve(dnorm(x,mean=Mu,sd=Sig),from=XLIM[1],to=XLIM[2],
        main="",xlab="",ylab="",cex.axis=2,lwd=2,yaxt="n")
  xcoords= seq(from=Mu+Sig*qnorm(alpha/2),to=Mu+Sig*qnorm(1-alpha/2),length.out=250)
  ycoords = c(0,dnorm(xcoords,Mu,Sig),0)
  xcoords = c(xcoords[1],xcoords,xcoords[length(xcoords)])
  polygon(xcoords,ycoords,col="skyblue")
}
.normCurvePlot2 = function(Mu,Sig){
  XLIM=c(min(Mu-4*Sig),max(Mu+4*Sig))
  curve(dnorm(x,mean=Mu[1],sd=Sig[1]),from=XLIM[1],to=XLIM[2],n=100000,
        main="",xlab="",ylab="",cex.axis=2,lwd=2,yaxt="n")
  par(new=TRUE)
  curve(dnorm(x,mean=Mu[2],sd=Sig[2]),from=XLIM[1],to=XLIM[2],n=100000,
        main="",xlab="",ylab="",cex.axis=2,lwd=2,yaxt="n",lty=2)
  xcoords= seq(from=XLIM[1],to=XLIM[2],length.out=500)
  ycoords1 = c(0,dnorm(xcoords,Mu[1],Sig[1]),0)
  ycoords2 = c(0,dnorm(xcoords,Mu[2],Sig[2]),0)
  ycoords = apply(cbind(ycoords1,ycoords2),1,min)
  xcoords = c(xcoords[1],xcoords,xcoords[length(xcoords)])
  polygon(xcoords,ycoords,col="skyblue")
}

.STAR.Compare = function(STARObj1,STARObj2,burnin=0,Titles=TRUE,savePlots=FALSE,prefix=""){
  if(STARObj1$method!=STARObj2$method){stop('Estimation methods of the two models are not the same')}
  
  if(STARObj1$method=='VB'){
    pp=length(STARObj1$BetaMean)
    qq=length(STARObj1$ThetaMean)
    for(ell in 1:pp){
      if(savePlots){
        jpeg(filename=paste(prefix,"betaComp",names(STARObj1$BetaMean)[ell],".jpg",sep=""),
             height=800,width=800)}
      if(Titles){par(mar=c(5,1,4,1))}else{par(mar=c(5,1,1,1))}
      .normCurvePlot2(c(STARObj1$BetaMean[ell],STARObj2$BetaMean[ell]),
                      c(sqrt(STARObj1$BetaThetaCov[ell,ell]),sqrt(STARObj2$BetaThetaCov[ell,ell])))
      if(Titles)title(main=paste("Comparing the effect of",names(STARObj1$BetaMean)[ell]),cex.main=2)
      if(savePlots){dev.off()}else{cat("Hit 'Enter' for next plot");readline()}
    }
    for(ell in 1:qq){
      if(savePlots){
        jpeg(filename=paste(prefix,"thetaComp",names(STARObj1$ThetaMean)[ell],".jpg",sep=""),
             height=800,width=800)
      }
      if(Titles){par(mar=c(5,1,4,1))}else{par(mar=c(5,1,1,1))}
      .normCurvePlot2(c(STARObj1$ThetaMean[ell],STARObj2$ThetaMean[ell]),
                      c(sqrt(STARObj1$BetaThetaCov[pp+ell,pp+ell]),
                        sqrt(STARObj2$BetaThetaCov[pp+ell,pp+ell])))
      if(Titles)title(main=paste("Comparing the effect of",names(STARObj1$ThetaMean)[ell]),cex.main=2)
      if(savePlots){dev.off()}else{if(ell<qq){cat("Hit 'Enter' for next plot");readline()}}
    }
  }else{
    pp=nrow(STARObj1$Beta)
    MM = ncol(STARObj1$Beta)
    qq=nrow(STARObj1$Theta)
    burnin = burnin+1
    for(ell in 1:pp){
      if(savePlots){
        jpeg(filename=paste(prefix,"betaComp",rownames(STARObj1$Beta)[ell],".jpg",sep=""),
             height=800,width=800)}
      if(Titles){par(mar=c(5,1,4,1))}else{par(mar=c(5,1,1,1))}
      d1 = density(STARObj1$Beta[ell,burnin:MM],adjust=3)
      d2 = density(STARObj2$Beta[ell,burnin:MM],adjust=3)
      xl = range(c(d1$x,d2$x))
      d1 = density(STARObj1$Beta[ell,burnin:MM],adjust=3,from=xl[1],to=xl[2])
      d2 = density(STARObj2$Beta[ell,burnin:MM],adjust=3,from=xl[1],to=xl[2])
      yl = range(c(d1$y,d2$y))
      plot(d1,ylim=yl,main="",xlab="",ylab="",cex.axis=2,lwd=2,yaxt="n")
      lines(d2,lwd=2,lty=2)
      yover = c(0,pmin(d1$y,d2$y),0)
      xcoords = c(d1$x[1],d1$x,d1$x[length(d1$x)])
      polygon(xcoords,yover,col="skyblue")
      if(Titles)title(main=paste("Comparing the effect of",rownames(STARObj1$Beta)[ell]),cex.main=2)
      if(savePlots){dev.off()}else{cat("Hit 'Enter' for next plot");readline()}
    }
    for(ell in 1:qq){
      if(savePlots){
        jpeg(filename=paste(prefix,"thetaComp",rownames(STARObj1$Theta)[ell],".jpg",sep=""),
             height=800,width=800)
      }
      if(Titles){par(mar=c(5,1,4,1))}else{par(mar=c(5,1,1,1))}
      d1 = density(STARObj1$Theta[ell,burnin:MM],adjust=3)
      d2 = density(STARObj2$Theta[ell,burnin:MM],adjust=3)
      xl = range(c(d1$x,d2$x))
      d1 = density(STARObj1$Theta[ell,burnin:MM],adjust=3,from=xl[1],to=xl[2])
      d2 = density(STARObj2$Theta[ell,burnin:MM],adjust=3,from=xl[1],to=xl[2])
      yl = range(c(d1$y,d2$y))
      plot(d1,ylim=yl,main="",xlab="",ylab="",cex.axis=2,lwd=2,yaxt="n")
      lines(d2,lwd=2,lty=2)
      yover = c(0,pmin(d1$y,d2$y),0)
      xcoords = c(d1$x[1],d1$x,d1$x[length(d1$x)])
      polygon(xcoords,yover,col="skyblue")
      if(Titles)title(main=paste("Comparing the effect of",rownames(STARObj1$Theta)[ell]),cex.main=2)
      if(savePlots){dev.off()}else{if(ell<qq){cat("Hit 'Enter' for next plot");readline()}}
    }
  }
}

plot.STAR = function(STARObj1,STARObj2=NULL,burnin=0,plotType=c("trace","epsBall","Comparison"),
                     Titles=TRUE,savePlots=FALSE,prefix="",
                     N=5e3,propVar=0.5){
  
  if("trace" %in% plotType){
    if(STARObj1$method != 'Gibbs'){warning('Cannot create trace plot for VB')}else{
      pp=nrow(STARObj1$Beta)
      MM = ncol(STARObj1$Beta)
      qq=nrow(STARObj1$Theta)
      for(ell in 1:pp){
        if(savePlots){
          jpeg(filename=paste(prefix,"betaTrace",rownames(STARObj1$Beta)[ell],".jpg",sep=""),
               height=800,width=800)}
        if(Titles){par(mar=c(5,4,4,1))}else{par(mar=c(5,4,1,1))}
        plot(STARObj1$Beta[ell,],type="l",lwd=2,cex.axis=2,ylab="",xlab="",main="")
        if(Titles)title(main=paste("Trace plot for",rownames(STARObj1$Beta)[ell]),cex.main=2)
        if(savePlots){dev.off()}else{cat("Hit 'Enter' for next plot");readline()}
      }
      for(ell in 1:qq){
        if(savePlots){
          jpeg(filename=paste(prefix,"thetaTrace",rownames(STARObj1$Theta)[ell],".jpg",sep=""),
               height=800,width=800)}
        if(Titles){par(mar=c(5,4,4,1))}else{par(mar=c(5,4,1,1))}
        plot(STARObj1$Theta[ell,],type="l",lwd=2,cex.axis=2,ylab="",xlab="",main="")
        if(Titles)title(main=paste("Trace plot for",rownames(STARObj1$Theta)[ell]),cex.main=2)
        if(savePlots){dev.off()}else{cat("Hit 'Enter' for next plot");readline()}
      }
    }
  }
  
  if("Comparison" %in% plotType){
    if(is.null(STARObj1) | is.null(STARObj2) ){
      warning("Must have both STARObj1 and STARObj2 to do plotType = 'Comparison'")
    }else{
      .STAR.Compare(STARObj1,STARObj2,burnin,Titles,savePlots,prefix)
    }
  }
  
  if("epsBall" %in% plotType){
    if(!STARObj1$RE){
      warning("STAR object must include random effects to do plotType = 'epsBall'")
    }else{
      .STAR.RETest(STARObj1,burnin,N,propVar,Titles,savePlots,prefix)
    }
  }
  
}



