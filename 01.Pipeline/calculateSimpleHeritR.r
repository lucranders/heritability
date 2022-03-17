library(dplyr)
library(foreach)
# Read binaries
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  close(BinFile)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
  
}
# Functions to estimate heritability
# By maximum likelihood - Multiple random variables can be considered
maximumLikelihood = function(listParams_){
  
  X = listParams_[['X']]
  y = listParams_[['y']]
  givenMatrixes = listParams_[['givenMatrixes']]
  sigmas2Init = listParams_[['sigmas2Init']]
  conv = listParams_[['conv']]
  maxIt = listParams_[['maxIt']]
  
  dimMatrixes = nrow(givenMatrixes[[1]])
  givenMatrixes[['Residual']] = diag(1,dimMatrixes)
  numberMatrixes = length(givenMatrixes)
  
  betas = list()
  logLik = c()
  vars2 = list()
  namesComponents = names(givenMatrixes)
  
  V = matrix(0 , nrow = dimMatrixes , ncol = dimMatrixes)
  for(i in 1:numberMatrixes){
    
    V = V + (givenMatrixes[[i]] * sigmas2Init[i])
  }
  vInverse = solve(V)
  beta_ = solve((t(X) %*% vInverse %*% X)) %*% t(X) %*% vInverse %*% y
  difMean = y- (X%*%beta_)
  likelihood_ = determinant(V)$modulus[[1]] +  t(difMean) %*% vInverse %*% difMean
  
  likelihood_ = -.5 * likelihood_
  logLik = c(likelihood_)
  betas[[1]] = beta_    
  
  
  k = 1
  while(k < maxIt){
    
    s = c()
    for(i in 1:numberMatrixes){
      
      # Score vector element
      x = sum(diag( vInverse %*% givenMatrixes[[i]] ) ) - t(difMean) %*% vInverse %*% givenMatrixes[[i]] %*% vInverse %*% difMean  
      s = c(s, -0.5*x)
      
    }
    
    
    Hessian = matrix(0, nrow = numberMatrixes, ncol = numberMatrixes)
    for(i in 1:numberMatrixes){
      for(j in 1:numberMatrixes){
        
        x = sum(diag( vInverse %*% givenMatrixes[[i]] %*% vInverse %*% givenMatrixes[[j]] ) )
        Hessian[i,j] = -.5*x
        
      }
    }
    
    invHessian = solve(Hessian)
    
    thetaNew = sigmas2Init - invHessian %*% s

    V = matrix(0 , nrow = dimMatrixes , ncol = dimMatrixes)    
    for(i in 1:numberMatrixes){
      sigmas2Init[i] = thetaNew[i]
      vars2[[namesComponents[i]]] = c(vars2[[namesComponents[i]]],thetaNew[i])
      V = V + (givenMatrixes[[i]] * thetaNew[i])
    }

    
    vInverse = solve(V)
    beta_ = solve((t(X) %*% vInverse %*% X)) %*% t(X) %*% vInverse %*% y
    betas[[i]] = beta_
    difMean = y- (X%*%beta_)
    
    likelihood_ = determinant(V)$modulus[[1]] +  t(difMean) %*% vInverse %*% difMean
    likelihood_ = -.5 * likelihood_
    
    logLik = c(logLik , likelihood_)
    
    k = k+1
    if( abs(logLik[k]-logLik[k-1]) < conv){
      break
    }
    
    
    
  }
  return_ = list()
  return_[['sigmas2']] = vars2
  # return_[['logLikelihood']] = logLik
  return_[['logLikelihood']] = likelihood_
  return_[['betas']] = betas
  return(return_)
}

# By classical REML - described in the literature
REML = function(listParams_){
  
  
  X = listParams_[['X']]
  y = listParams_[['y']]
  givenMatrixes = listParams_[['givenMatrixes']]
  sigmas2Init = listParams_[['sigmas2Init']]
  conv = listParams_[['conv']]
  maxIt = listParams_[['maxIt']]
  
  dimMatrixes = nrow(givenMatrixes[[1]])
  givenMatrixes[['Residual']] = diag(1,dimMatrixes)
  numberMatrixes = length(givenMatrixes)
  
  betas = list()
  logLik = c()
  vars2 = list()
  namesComponents = names(givenMatrixes)
  
  V = matrix(0 , nrow = dimMatrixes , ncol = dimMatrixes)
  for(i in 1:numberMatrixes){
    
    V = V + (givenMatrixes[[i]] * sigmas2Init[i])

  }
  
  vInverse = solve(V)
  beta_ = solve((t(X) %*% vInverse %*% X)) %*% t(X) %*% vInverse %*% y
  difMean = y- (X%*%beta_)
  P = vInverse - ( vInverse %*% X %*% solve( t(X) %*% vInverse %*% X)  %*% t(X) %*% vInverse )
  likelihood_ = determinant(V)$modulus[[1]] + 
    determinant(t(X) %*% solve(V) %*% X)$modulus[[1]] + 
    t(difMean) %*% vInverse %*% difMean
  
  likelihood_ = -.5 * likelihood_
  logLik = c(likelihood_)
  betas[[1]] = beta_   
  
  
  k = 1
  while(k < maxIt){

    s = c()
    for(i in 1:numberMatrixes){
      
      # Score vector element
      x = sum(diag( P %*% givenMatrixes[[i]] ) ) - t(difMean) %*% vInverse %*% givenMatrixes[[i]] %*% vInverse %*% difMean  
      s = c(s, -0.5*x)
      
    }
    
    Hessian = matrix(0, nrow = numberMatrixes, ncol = numberMatrixes)
    for(i in 1:numberMatrixes){
      for(j in 1:numberMatrixes){
        
        x = sum(diag( P %*% givenMatrixes[[i]] %*% P %*% givenMatrixes[[j]] ) )
        Hessian[i,j] = -.5*x
        
      }
    }
    
    invHessian = solve(Hessian)
    thetaNew = sigmas2Init - invHessian %*% s
    V = matrix(0 , nrow = dimMatrixes , ncol = dimMatrixes)
    for(i in 1:numberMatrixes){
      sigmas2Init[i] = thetaNew[i]
      vars2[[namesComponents[i]]] = c(vars2[[namesComponents[i]]],thetaNew[i])
      V = V + (givenMatrixes[[i]] * thetaNew[i])
    }
        
    vInverse = solve(V)
    beta_ = solve((t(X) %*% vInverse %*% X)) %*% t(X) %*% vInverse %*% y
    
    betas[[k]] = beta_
    difMean = y- (X%*%beta_)
    P = vInverse - ( vInverse %*% X %*% solve( t(X) %*% vInverse %*% X)  %*% t(X) %*% vInverse )
    likelihood_ = determinant(V)$modulus[[1]] + 
    determinant(t(X) %*% solve(V) %*% X)$modulus[[1]] + 
    t(difMean) %*% vInverse %*% difMean
  
    likelihood_ = -.5 * likelihood_
    logLik = c(logLik , likelihood_)
    
    k = k+1
    if( abs(logLik[k]-logLik[k-1]) < conv){
      break
    }

    
  }
  return_ = list()
  return_[['sigmas2']] = vars2
  # return_[['logLikelihood']] = logLik
  return_[['logLikelihood']] = likelihood_
  return_[['betas']] = betas
  return(return_)
}


# By REML - as literaly described in the literature
REML2 = function(listParams_){
  
  X = listParams_[['X']]
  y = listParams_[['y']]
  givenMatrixes = listParams_[['givenMatrixes']]
  sigmas2Init = listParams_[['sigmas2Init']]
  conv = listParams_[['conv']]
  maxIt = listParams_[['maxIt']]
  samp_ = listParams_[['samp_']]
  
  dimMatrixes = nrow(givenMatrixes[[1]])
  givenMatrixes[['Residual']] = diag(1,dimMatrixes)
  numberMatrixes = length(givenMatrixes)
  
  
  # transf = diag(1,nrow = nrow(X)) - X%*% solve(t(X) %*% X) %*% t(X)
  ACt = diag(1 , nrow(X)) - (X %*% solve(t(X) %*% X) %*% t(X) )
  ASample = ACt[samp_,]
  
  newX = ASample %*% X
  newY = ASample %*% y
  
  betas = list()
  logLik = c()
  vars2 = list()
  namesComponents = names(givenMatrixes)
  
  V = matrix(0 , nrow = dimMatrixes , ncol = dimMatrixes)
  for(i in 1:numberMatrixes){
    
    V = V + (givenMatrixes[[i]] * sigmas2Init[i])
  }
  
  newV =  ASample %*% V %*% t(ASample)
  newVInverse = solve(newV)
  likelihood_ = determinant(newV)$modulus[[1]] +  t(newY) %*% newVInverse %*% newY
  
  likelihood_ = -.5 * likelihood_
  logLik = c(likelihood_)
  vInverse = solve(V)
  beta_ = solve((t(X) %*% vInverse %*% X)) %*% t(X) %*% vInverse %*% y
  betas[[1]] = beta_
  
  k = 1
  while(k < maxIt){
    
    s = c()
    for(i in 1:numberMatrixes){
      
      # Derivative from transformation
      x0 = ASample %*% givenMatrixes[[i]] %*% t(ASample)
      # Element from score vector
      x1 = sum(diag( newVInverse %*% x0 ) ) - t(newY) %*% newVInverse %*% x0 %*% newVInverse %*% newY  
      # Score vector element
      s = c(s, -0.5*x1)
      
    }
    
    Hessian = matrix(0, nrow = numberMatrixes, ncol = numberMatrixes)
    for(i in 1:numberMatrixes){
      for(j in 1:numberMatrixes){
        
        x01 = ASample %*% givenMatrixes[[i]] %*% t(ASample)
        x02 = ASample %*% givenMatrixes[[j]] %*% t(ASample)
        x = sum(diag( newVInverse %*% x01 %*% newVInverse %*% x02 ) )
        Hessian[i,j] = -.5*x
        
      }
    }
    
    invHessian = solve(Hessian)
    
    thetaNew = sigmas2Init - invHessian %*% s
    V = matrix(0 , nrow = dimMatrixes , ncol = dimMatrixes)    
    for(i in 1:numberMatrixes){
      sigmas2Init[i] = thetaNew[i]
      vars2[[namesComponents[i]]] = c(vars2[[namesComponents[i]]] , thetaNew[i])
      V = V + (givenMatrixes[[i]] * thetaNew[i])
    }
    
    newV =  ASample %*% V %*% t(ASample)
    newVInverse = solve(newV)
    vInverse = solve(V)

    beta_ = solve((t(X) %*% vInverse %*% X)) %*% t(X) %*% vInverse %*% y
    betas[[k]] = beta_

    likelihood_ = determinant(newV)$modulus[[1]] +  t(newY) %*% newVInverse %*% newY 
    likelihood_ = -.5 * likelihood_
    logLik = c(logLik,likelihood_)
    
    k = k+1
    if( abs(logLik[k]-logLik[k-1]) < conv){
      break
    }
    
    
  }
  return_ = list()
  return_[['sigmas2']] = vars2
  # return_[['logLikelihood']] = logLik
  return_[['logLikelihood']] = likelihood_
  return_[['betas']] = betas
  return(return_)
}

# Useful to choose not LD observations (necessary do REML2)
notLD = function(X,seed_,tol_){
  
  set.seed(seed_)
  ACt = diag(1 , nrow(X)) - (X %*% solve(t(X) %*% X) %*% t(X) )
  samp_ = sample(1:nrow(X) , replace = F , size = (nrow(X) - (ncol(X) ) ) )
  ASample = ACt[samp_,]
  teste_ = svd(ASample)
  while(sum(teste_$d < tol_) > 0){
    samp_ = sample(1:nrow(X) , replace = F , size = (nrow(X) - (ncol(X) ) ) )
    ASample = ACt[samp_,]
    teste_ = svd(ASample)
    print(sum(teste_$d < tol_))
  }
  return(samp_)
  
}

# considering populations as different groups (i.e, independent samples)
#### TO DO: rewrite this function
likelihoodGroups = function(listXs,listYs,listAs,sigmaA,sigmaE,conv,maxIt){
  
  # transf = diag(1,nrow = nrow(X)) - X%*% solve(t(X) %*% X) %*% t(X)
  # newX = 
  valsA = c()
  valsE = c()
  betas = list()
  logLik = c()
  
  i=1
  
  Vs = list()
  vInverses = list()
  pt1 = matrix(0,ncol = ncol(listXs[[1]]), nrow = ncol(listXs[[1]]))
  pt2 = matrix(0,ncol = 1, nrow = ncol(listXs[[1]]))
  logliks = 0
  cols_ = 1
  for(name_ in names(listAs)){
    
    Xk = listXs[[name_]]
    Ak = listAs[[name_]]
    yk = listYs[[name_]]
    rows_ = dim(Ak)[[1]]
    E_ = diag(1,rows_)* sigmaE
    V_ = (Ak*sigmaA) + (E_)
    Vs[[name_]] = V_
    VInverse_ = solve(V_)
    vInverses[[name_]] = VInverse_
    
    pt1 = pt1 + t(Xk) %*% VInverse_ %*% Xk
    pt2 = pt2 + t(Xk) %*% VInverse_ %*% yk
    
  }
  
  betaMean = solve(pt1) %*% pt2
  betas[[1]] = betaMean
  likelihood1_ = 0
  likelihood2_ = 0
  likelihood3_ = 0
  s1 = 0
  s2 = 0
  H11 = 0
  H12 = 0
  H21 = 0
  H22 = 0
  for(name_ in names(listAs)){
    
    V_ = Vs[[name_]]
    vInverse = solve(V_)
    Xk = listXs[[name_]]
    yk = listYs[[name_]]
    difMean =  yk - (Xk %*% betaMean)
    likelihood1_ = likelihood1_ + determinant(V_)$modulus[[1]]
    likelihood2_ = likelihood2_ + determinant( t(Xk) %*% vInverse %*% Xk )$modulus[[1]]
    likelihood3_ = likelihood3_ + t(difMean) %*% vInverse %*% difMean
    derivSig1Pop = listAs[[name_]]
    derivSig2Pop = diag(1,nrow = nrow(listAs[[name_]]))
    
    
    s1 = s1 + sum(diag( vInverse %*%  derivSig1Pop ) ) - t(difMean) %*% vInverse %*% derivSig1Pop %*% vInverse %*% difMean  
    s2 = s2 + sum(diag( vInverse %*%  derivSig2Pop ) ) - t(difMean) %*% vInverse %*% derivSig2Pop %*% vInverse %*% difMean
    
    
    H11 = H11 + sum(diag( vInverse %*% derivSig1Pop %*% vInverse %*% derivSig1Pop ) )
    H12 = H12 + sum(diag( vInverse %*% derivSig1Pop %*% vInverse %*% derivSig2Pop ) )
    H21 = H21 + sum(diag( vInverse %*% derivSig2Pop %*% vInverse %*% derivSig1Pop ) )
    H22 = H22 + sum(diag( vInverse %*% derivSig2Pop %*% vInverse %*% derivSig2Pop ) )
    
  }
  
  s = -.5 * c(s1,s2)
  Hessian = -.5*matrix(c(H11,H12,H21,H22),byrow = T,nrow=2)
  invHessian = solve(Hessian)
  
  thetaNew = c(sigmaA,sigmaE) - invHessian %*% s
  
  
  sigmaA = thetaNew[1]
  sigmaE = thetaNew[2]
  valsA = c(valsA,sigmaA)
  valsE = c(valsE,sigmaE)
  
  likelihood_ = -.5 * (likelihood1_ + likelihood2_ + likelihood3_)
  logLik = c(logLik , likelihood_)
  
  
  
  while(i < maxIt){
    
    i = i+1
    logliks = 0
    cols_ = 1
    pt1 = matrix(0,ncol = ncol(listXs[[1]]), nrow = ncol(listXs[[1]]))
    pt2 = matrix(0,ncol = 1, nrow = ncol(listXs[[1]]))
    for(name_ in names(listAs)){
      
      Xk = listXs[[name_]]
      Ak = listAs[[name_]]
      yk = listYs[[name_]]
      rows_ = dim(Ak)[[1]]
      E_ = diag(1,rows_)* sigmaE
      V_ = (Ak*sigmaA) + (E_)
      Vs[[name_]] = V_
      VInverse_ = solve(V_)
      vInverses[[name_]] = VInverse_
      
      
      pt1 = pt1 + t(Xk) %*% VInverse_ %*% Xk
      pt2 = pt2 + t(Xk) %*% VInverse_ %*% yk
      
      
    }
    betaMean = solve(pt1) %*% pt2
    betas[[i]] = betaMean
    likelihood1_ = 0
    likelihood2_ = 0
    likelihood3_ = 0
    s1 = 0
    s2 = 0
    H11 = 0
    H12 = 0
    H21 = 0
    H22 = 0
    for(name_ in names(listAs)){
      
      V_ = Vs[[name_]]
      vInverse = solve(V_)
      Xk = listXs[[name_]]
      yk = listYs[[name_]]
      difMean =  yk - (Xk %*% betaMean)
      likelihood1_ = likelihood1_ + determinant(V_)$modulus[[1]]
      likelihood2_ = likelihood2_ + determinant( t(Xk) %*% vInverse %*% Xk )$modulus[[1]]
      likelihood3_ = likelihood3_ + t(difMean) %*% vInverse %*% difMean
      derivSig1Pop = listAs[[name_]]
      derivSig2Pop = diag(1,nrow = nrow(listAs[[name_]]))
      
      
      s1 = s1 + sum(diag( vInverse %*%  derivSig1Pop ) ) - t(difMean) %*% vInverse %*% derivSig1Pop %*% vInverse %*% difMean  
      s2 = s2 + sum(diag( vInverse %*%  derivSig2Pop ) ) - t(difMean) %*% vInverse %*% derivSig2Pop %*% vInverse %*% difMean
      
      
      H11 = H11 + sum(diag( vInverse %*% derivSig1Pop %*% vInverse %*% derivSig1Pop ) )
      H12 = H12 + sum(diag( vInverse %*% derivSig1Pop %*% vInverse %*% derivSig2Pop ) )
      H21 = H21 + sum(diag( vInverse %*% derivSig2Pop %*% vInverse %*% derivSig1Pop ) )
      H22 = H22 + sum(diag( vInverse %*% derivSig2Pop %*% vInverse %*% derivSig2Pop ) )
      
    }
    
    s = -.5 * c(s1,s2)
    Hessian = -.5*matrix(c(H11,H12,H21,H22),byrow = T,nrow=2)
    invHessian = solve(Hessian)
    
    thetaNew = c(sigmaA,sigmaE) - invHessian %*% s
    
    
    sigmaA = thetaNew[1]
    sigmaE = thetaNew[2]
    valsA = c(valsA,sigmaA)
    valsE = c(valsE,sigmaE)
    
    likelihood_ = -.5 * (likelihood1_ + likelihood2_ + likelihood3_)
    logLik = c(logLik , likelihood_)
    if( abs(logLik[i]-logLik[i-1]) < conv){
      break
    }
    
    
    
    
  }
  return_ = list()
  return_[['sigmas2A']] = valsA
  return_[['sigmas2E']] = valsE
  # return_[['logLikelihood']] = logLik
  return_[['logLikelihood']] = likelihood_
  return_[['betas']] = betas
  return(return_)
}


bestStart = function(model_,starts_,cl_,listParams_,nameInterest_){
  
  cl <- parallel::makeCluster(cl_)
  doParallel::registerDoParallel(cl)
  func_ = get(model_)
  
  run_ <- foreach(sigma_ = starts_,.combine = 'rbind') %dopar% {
    listParams_[['sigmas2Init']] = rep(sigma_,(length(listParams_[['givenMatrixes']]) + 1) )
    exec_ = func_(listParams_)
    
    loglik = exec_[['logLikelihood']]
    sigmas2 = exec_[['sigmas2']]
    sigmaInterest = sigmas2[[nameInterest_]][[length(sigmas2[[nameInterest_]])]]
    sigma2Tot = 0
    for(name_ in names(sigmas2)){
      sigma2Tot = sigma2Tot + sigmas2[[name_]]
    }
    sigma2TotFinal = sigma2Tot[[length(sigmas2[[nameInterest_]])]]
    data.frame(start = sigma_,likelihood = loglik,heritability = sigmaInterest/sigma2TotFinal , sigma2A = sigmaInterest , sigma2Tot = sigma2TotFinal)
    
  }
  
  parallel::stopCluster(cl)
  return(run_)
}

returnBestStart = function(df_,model_,listParams_){
  bestOpt = df_ %>% summarise(b_ = min(likelihood))
  filt_ = df_ %>% filter(likelihood == bestOpt)
  conf_ = filt_[1,'start']
  func_ = get(model_)
  lastExec = func_(listParams_)
  return(lastExec)
}

pathTmp_ =  commandArgs(TRUE)[1]
gene_ = commandArgs(TRUE)[2]
matrixName_ = commandArgs(TRUE)[3]
method_ = commandArgs(TRUE)[4]
fileSave_ = commandArgs(TRUE)[5]

dbXName = paste0(pathTmp_,'/X_',gsub('-','',gene_),'.txt')
dbYName = paste0(pathTmp_,'/Y_',gsub('-','',gene_),'.txt')
GRMName = paste0(pathTmp_,'/',matrixName_)

X_covs = read.table(dbXName,sep = ' ',header = F)
X_covs$V1 = NULL
X_covs$V2 = NULL
Y_ = read.table(dbYName,sep = ' ',header = F)
Y_$V1 = NULL
Y_$V2 = NULL

GRM = ReadGRMBin(GRMName)
M <- matrix(0, dim(GRM[[3]])[1], dim(GRM[[3]])[1])
M[upper.tri(M, diag = FALSE)] = GRM[[2]]
M = t(M)
M[upper.tri(M, diag = FALSE)] = GRM[[2]]
diag(M) = GRM[[1]]
(M[1:4,1:4])

listParams_ = list()

listParams_[['X']] = model.matrix(~. , X_covs)
listParams_[['y']] = Y_[,1]
listParams_[['A']] = M
listParams_[['givenMatrixes']] = list(Genes = M)
val_ = 30
listParams_[['sigmas2Init']] = rep(val_ , (length(listParams_[['givenMatrixes']]) + 1) )
listParams_[['conv']] = 1e-4
listParams_[['maxIt']] = 100
# listParams_[['samp_']] = notLD(listParams_[['X']],9297791,0.01)

df_ = bestStart(model_ = method_,starts_ = c(sd(listParams_[['y']]),sd(listParams_[['y']])/2,var(listParams_[['y']]),var(listParams_[['y']])/2),cl_ = 4,listParams_ = listParams_,'Genes')
mod_ = returnBestStart(df_ = df_,model_ = method_,listParams_ = listParams_)

sigmas2 = mod_[['sigmas2']]
sigma2Tot = 0
for(name_ in names(sigmas2)){
  sigma2Tot = sigma2Tot + sigmas2[[name_]]
}
heritability = sigmas2[['Genes']]/sigma2Tot

finalSaveFile = paste0(pathTmp_,'/',fileSave_,'_',gsub('-','',gene_),'.txt')
heritFinal = heritability[[length(heritability)]]
write.table(paste0("V(G)/Vp ",heritFinal),finalSaveFile,append = T, row.names = F, quote = F, sep = '\t')
sigma2GeneFinal = sigmas2[['Genes']][[length(heritability)]]
write.table(paste0("V(G)  ",sigma2GeneFinal),finalSaveFile,append = T, row.names = F, quote = F, sep = '\t')
sigma2TotFinal = sigma2Tot[[length(heritability)]]
write.table(paste0("Vp  ",sigma2TotFinal),finalSaveFile,append = T, row.names = F, quote = F, sep = '\t')