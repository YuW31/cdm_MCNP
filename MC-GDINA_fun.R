# # functions for MC-DINA model
# dat = N0.response
# Qc = mc.q
# model = "MCDINA"
# key = key.all
# conv.crit = .001
# maxitr=2000
# conv.type="pr"
# SE=FALSE

MCmodel_Yu <- function(dat, Qc, model = "MCDINA", key = NULL, H,conv.crit = .001,maxitr=2000,conv.type="pr",SE=FALSE){
  Q <- Qc[, -c(1:2)]
  N <- nrow(dat)
  J <- ncol(dat)
  unique.code <- apply(dat, 2, unique)
  ### code to deal with the conflict: H-1 is the ID
  missing.item = which(sapply(unique.code,
                             function(x){
                               length(x)<H
                             }))
  ### DONE ###
  if (is.list(unique.code)) {
    C <- sapply(unique.code, length)
    
    for (j in 1:length(missing.item)){
      miss.j = missing.item[j]
      diff.c = H - length(unlist(unique.code[miss.j]))
      miss.c = (1:H)[-unlist(unique.code[miss.j])]
      for (c in 1:diff.c){
        
        reorder.c = sub.mcQ(Qc,miss.j)[,2,drop = FALSE]
        reorder.c[reorder.c > miss.c[c]] = reorder.c[reorder.c > miss.c[c]] - diff.c
        Qc[Qc[,1]==miss.j,2] = reorder.c
        
        dat[dat[,miss.j] > miss.c[c],miss.j] = dat[dat[,miss.j] > miss.c[c],miss.j]-diff.c
        
        key[miss.j] = ifelse(key[miss.j] > miss.c[c],key[miss.j]-diff.c,key[miss.j])
        
      }
      
      
    }
    
    
  } else if (is.matrix(unique.code)) {
    C <- rep(nrow(unique.code), ncol(unique.code))
  }
  # C <- rep(H,J)
  S0 <- sum(C)
  item.no <- rep(1:J, times = C - 1)
  item.no.0 <- rep(1:J, times = C)
  K <- ncol(Q)
  L <- 2^K
  
  
  
  wg.break <- tryCatch(
    
    {init <- eta.mcm(Qc, model = model, no.options = C, key=key) # C to H
    param <- pl(init, model = model)
    eta2 <- init$eta
    logprior <- log(rep(1/L, L))
    dif <- LL.2 <- 1
    itr <- 0
    
    while (itr < maxitr) {
  
    p0 <- param
    p <- pl2pm(param, eta2)
    
    if (any(p < 0) || any(p > 1))
      stop("some item success probabilities are not between 0 and 1.",
           call. = FALSE)
    
    likepost <- Lik_DTM(as.matrix(p), as.matrix(dat - 1),
                        C - 1, logprior)
    
    
    R1 <- Rljs_DTM(likepost$logpost, dat - 1, C - 1)
    for (j in seq_len(J)) {
      param[[j]] <- ColNormalize(aggregateCol(R1[which(item.no.0 == j), ], eta2[j, ]))
    }
    
    LL.1 <- -2 * likepost$LL
    if (tolower(conv.type) == "LL") {
      dif <- abs(LL.1 - LL.2)
    } else if (tolower(conv.type) == "pr") {
      dif <- max(abs(unlist(param) - unlist(p0)))
    }
    itr <- itr + 1
    cat("\rIter =", itr, " Max. abs. change =", formatC(dif, digits = 5, format = "f"), " Deviance  =", formatC(LL.1, digits = 3, format = "f"), "                 ")
    
    if (dif < conv.crit)
      break
    LL.2 <- LL.1
    logprior <- c(likepost$logprior)
  }},error = function(e){
    wg.break =TRUE
    return(wg.break)
  })
  
  patt <- attributepattern(K)
  
  if (is.null(wg.break)){
    p <- pl2pm(param,eta2)
    est_item = vector(mode = "list",length = J) #likelihood
    for (j in 1:J){
      if (j==1){
        est_item[[j]] = p[1:C[j],]
      } else {
        est_item[[j]] = p[(sum(C[1:(j-1)])+1):sum(C[1:j]),]
      }
      
      
      
    }
    
    likepost <- Lik_DTM(as.matrix(p),as.matrix(dat-1),C-1,logprior)
    att <- list(EAP = ((exp(likepost$logpost) %*% patt) > 0.5) * 1,
                MAP = patt[apply(likepost$logpost, 1, which.max.randomtie), ],
                MLE = patt[apply(likepost$loglik, 1, which.max.randomtie), ])
    par.att <- att$MLE
  } else {
    par.att=matrix(9,nrow=N,ncol=K)
    est_item = NULL
    
  }

  postP = matrix(exp(likepost$logprior),nrow = 2^K,ncol = 1) #posterior of each latent class
  
  return(list("person" = par.att,"item" = est_item,"posterior.prob" = postP))
}

MCmodel_Yu_C <- function(dat, Qc, model = "MCDINA", key = NULL, H,conv.crit = .001,maxitr=2000,conv.type="pr",SE=FALSE){
  Q <- Qc[, -c(1:2)]
  N <- nrow(dat)
  J <- ncol(dat)
  unique.code <- apply(dat, 2, unique)
  # if (is.list(unique.code)) {
  #   C <- sapply(unique.code, length)
  # } else if (is.matrix(unique.code)) {
  #   C <- rep(nrow(unique.code), ncol(unique.code))
  # }
  C <- rep(H,J)
  S0 <- sum(C)
  item.no <- rep(1:J, times = C - 1)
  item.no.0 <- rep(1:J, times = C)
  K <- ncol(Q)
  L <- 2^K
  
  init <- eta.mcm(Qc, model = model, no.options = C,key=key)
  param <- pl(init, model = model)
  eta2 <- init$eta
  logprior <- log(rep(1/L, L))
  dif <- LL.2 <- 1
  itr <- 0
  
  
  wg.break <- tryCatch( while (itr < maxitr) {
    p0 <- param
    p <- pl2pm(param, eta2)
    
    if (any(p < 0) || any(p > 1))
      stop("some item success probabilities are not between 0 and 1.",
           call. = FALSE)
    
    likepost <- Lik_DTM(as.matrix(p), as.matrix(dat - 1),
                        C - 1, logprior)
    
    
    R1 <- Rljs_DTM(likepost$logpost, dat - 1, C - 1)
    for (j in seq_len(J)) {
      param[[j]] <- ColNormalize(aggregateCol(R1[which(item.no.0 == j), ], eta2[j, ]))
    }
    
    LL.1 <- -2 * likepost$LL
    if (tolower(conv.type) == "LL") {
      dif <- abs(LL.1 - LL.2)
    } else if (tolower(conv.type) == "pr") {
      dif <- max(abs(unlist(param) - unlist(p0)))
    }
    itr <- itr + 1
    cat("\rIter =", itr, " Max. abs. change =", formatC(dif, digits = 5, format = "f"), " Deviance  =", formatC(LL.1, digits = 3, format = "f"), "                 ")
    
    if (dif < conv.crit)
      break
    LL.2 <- LL.1
    logprior <- c(likepost$logprior)
  },error = function(e){
    wg.break =TRUE
    return(wg.break)
  })
  
  patt <- attributepattern(K)
  
  if (is.null(wg.break)){
    p <- pl2pm(param,eta2)
    est_item = vector(mode = "list",length = J) #likelihood
    for (j in 1:J){
      if (j==1){
        est_item[[j]] = p[1:C[j],]
      } else {
        est_item[[j]] = p[(sum(C[1:(j-1)])+1):sum(C[1:j]),]
      }
      
      
      
    }
    
    likepost <- Lik_DTM(as.matrix(p),as.matrix(dat-1),C-1,logprior)
    att <- list(EAP = ((exp(likepost$logpost) %*% patt) > 0.5) * 1,
                MAP = patt[apply(likepost$logpost, 1, which.max.randomtie), ],
                MLE = patt[apply(likepost$loglik, 1, which.max.randomtie), ])
    par.att <- att$MLE
  } else {
    par.att=matrix(9,nrow=N,ncol=K)
    est_item = NULL
    
  }
  
  postP = matrix(exp(likepost$logprior),nrow = 2^K,ncol = 1) #posterior of each latent class
  
  return(list("person" = par.att,"item" = est_item,"posterior.prob" = postP))
}


eta <- function(Q, AlphaPattern = NULL) {
  .Call('_GDINA_eta', PACKAGE = 'GDINA', Q, AlphaPattern)
}   # get possible class

eta.mcm <- function(Qc, model = "MCDINA",no.options = C, key=key){
  
  K = ncol(Qc)-2
  Q <- Qc[,3:(2+K),drop = FALSE]
  
  item.no <- Qc[,1]
  coded.cat.no <- Qc[,2]
  
  J <- length(unique(item.no))
  
  Cj <- tabulate(item.no)
  
  K <- ncol(Q)
  
  if(length(no.options)==1){
    no.options <- rep(no.options, J)
  }else if(length(no.options)!=J){
    stop("no.options must be of length 1 or the number of items.",call. = FALSE)
  }
  
  et <- matrix(0,J,2^K)
  et.label <- list()
  m <- list()
  if(model=="MCDINA"){
    for(j in 1:J){
      
      Qj <- Q[which(item.no==j),,drop=FALSE]
      Kj <- rowSums(Qj)
      kj.order <- order(Kj,decreasing = TRUE)
      coded.cat.no.j <- coded.cat.no[which(item.no==j)]
      if(any(duplicated(coded.cat.no.j))||any(Kj==0))
        stop("Q-matrix for item",j,"is not correctly specified.",call. = FALSE)
      
      if (Cj[j]>1){
        if(!is.null(key)){
          if(key[j]%in%coded.cat.no.j){
            key.loc <- which(coded.cat.no.j==key[j])
            if(key.loc!=kj.order[1]){
              kj.order <- c(key.loc,setdiff(kj.order,key.loc))  # make the answer goes first
            }
          }else{
            stop("Option key for item",j,"is not given in the Q-matrix.",call. = FALSE)
          }
        }
        Qj <- Qj[kj.order,]
        et.label[[j]] <- c(apply(Qj,1,paste,collapse=""),paste(rep("*",ncol(Qj)),collapse = ""))
        m[[j]] <- matrix(0,no.options[j],length(et.label[[j]]))
        m[[j]][matrix(c(coded.cat.no.j[kj.order],seq_len(Cj[j])),nrow = Cj[j])] <- 1
        # m class?
        tmp <- eta(Qj)
        eta.j <- apply(rbind(0,(tmp==apply(tmp,1,max))*1),2,which.max)  #get the first max
        # apply(tmp,1,max): the most possible class for each choice
        max.j <- max(eta.j)
        eta.j <- eta.j - 1
        eta.j[eta.j==0] <- max.j
        et[j,] <- eta.j
        # ideal response of class: 1: correct; 2,3,4 remaining choices; the last number means random
        
      }else{
        et.label[[j]] <- c(paste(Qj,collapse=""),paste(rep("*",ncol(Qj)),collapse = ""))
        eta.j <- eta(Qj)
        et[j,] <- 2 - (eta.j==max(eta.j))
        m[[j]] <- matrix(0,no.options[j],length(et.label[[j]]))
        m[[j]][matrix(c(coded.cat.no.j[kj.order],seq_len(Cj[j])),nrow = Cj[j])] <- 1
      }
      
    }
  }else if(model=="GNDM"){
    new.Q <- unrestrQ(Qc)[which(!duplicated(item.no)),-c(1:2)]
    et <- eta(new.Q)
    Kj <- rowSums(new.Q)
    for(j in 1:J){
      Qj <- Q[which(item.no==j),which(new.Q[j,]==1),drop=FALSE]
      Kjc <- rowSums(Qj)
      coded.cat.no.j <- coded.cat.no[which(item.no==j)]
      if(any(duplicated(coded.cat.no.j))||any(Kjc==0))
        stop("Q-matrix for item",j,"is not correctly specified.",call. = FALSE)
      att <- attributepattern(Kj[j])
      et.label[[j]] <- apply(att,1,paste,collapse="")
      m[[j]] <- matrix(0,no.options[j],length(et.label[[j]]))
      if(is.null(key)){
        loc.key <- which.max(Kjc)
      }else{
        if(!(key[j]%in%coded.cat.no.j)){
          stop("Option key for item",j,"is not given in the Q-matrix.",call. = FALSE)
        }else{
          loc.key <- which(coded.cat.no.j==key[j])
        }
      }
      tmp <- c(att%*%Qj[loc.key,])
      key.col <- which(tmp==max(tmp))
      m[[j]][coded.cat.no.j[loc.key],key.col] <- 1
      if(length(coded.cat.no.j)>1){
        for(jj in seq_len(length(coded.cat.no.j))){
          if(jj!=loc.key){
            tmp <- c(att%*%Qj[jj,])
            new.col <- setdiff(which(tmp==max(tmp)),key.col)
            m[[j]][coded.cat.no.j[jj],new.col] <- 1
            key.col <- c(new.col,key.col)
          }
        }
      }
    }
  }
  
  colnames(et) <- apply(attributepattern(K),1,paste,collapse="")
  rownames(et) <- names(et.label) <- paste("Item",seq_len(J))
  
  return(list(m=m,eta=et,label=et.label))
  
}

# m: eta for each option

pl <- function(eta.obj,model = "MCDINA"){
  if(model=="MCDINA"){
    p <- eta.obj$m
    for(j in seq_len(length(p))){
      #p[[j]][p[[j]]==1] <- .8
      #p[[j]][p[[j]]==0] <- .2/(nrow(p[[j]])-1)
      rp = runif(1,0.5,0.9999999999)
      p[[j]][p[[j]]==1] <- rp
      p[[j]][p[[j]]==0] <- (1-rp)/(nrow(p[[j]])-1)
      p[[j]][,ncol(p[[j]])] <- 1/nrow(p[[j]])
      
    }
  }else if(model=="GNDM"){
    p <- eta.obj$m
    for(j in seq_len(length(p))){
      cs <- colSums(p[[j]])
      p[[j]][,which(cs==0)] <- 1/nrow(p[[j]])
      p[[j]][p[[j]]==1] <- .8
      p[[j]][p[[j]]==0] <- .2/(nrow(p[[j]])-1)
    }
  }
  p
}

pl2pm <- function(plist,eta){
  S0 <- sum(sapply(plist,nrow))
  L <- ncol(eta)
  J <- nrow(eta)
  p <- matrix(0,S0,L)

  for(l in seq_len(L)){
    loc <- 1
    for(j in seq_len(J)){
      pj <- plist[[j]][,eta[j,l]]
      p[loc:(loc+length(pj)-1),l] <- pj
      loc <- loc + length(pj)
    }
  }
  p
}

Lik_DTM <- function(mP, mX, vC, vlogPrior) {
  .Call('_GDINA_Lik_DTM', PACKAGE = 'GDINA', mP, mX, vC, vlogPrior)
}

# rewrite Lik_DTM in R:
# likelihood = matrix(nrow = N, ncol = 2^K)
# for (i in 1:N){
#   p.list = vector(mode = "list",length = 2^K)
#   j_P_2k = matrix(nrow = J,ncol = 2^K)
#   for (m in 1:2^K){
#     p.list[[m]] = matrix(p[,m],nrow = 4,ncol = J,byrow = FALSE)
#     
#     j_P <- c()
#     for (j in 1:J){
#       j_P <- c(j_P,p.list[[m]][observed.response[i,j],j])
#     }
#     j_P_2k[,m] = j_P
#   }
#   
#   likelihood[i,] = apply(j_P_2k,2,prod)
# }
# 
# log(likelihood)[1,]
# likepost$loglik[1,]
# 
# mlogpost = log(likelihood) - matrix(logprior, nrow = N, ncol = 2^K, byrow = TRUE)
# logsumPost = log(rowSums(exp(mlogpost)))
# mlogpost = mlogpost-matrix(logsumPost,nrow = N, ncol = 8, byrow = FALSE)
# mlogpost[1,]
# likepost$logpost[1,]


Rljs_DTM <- function(mlogPost, mX, vC) {
  .Call('_GDINA_Rljs_DTM', PACKAGE = 'GDINA', mlogPost, mX, vC)
}

ColNormalize <- function(X) {
  .Call('_GDINA_ColNormalize', PACKAGE = 'GDINA', X)
}

score.mcm.j <- function(Xj,
                        parloc.j,
                        catprob.j,
                        logpost){

  Xj[is.na(Xj)] <- -1 # remove NA from the data
  post <- exp(logpost) # posterior N x 2^K
  score.p <- vector("list",nrow(catprob.j)-1)
  for(r in 1:length(score.p)){
    score.p[[r]] <- aggregateCol(post,parloc.j)*
      (outer((Xj==r),catprob.j[r,],"/")-outer((Xj==nrow(catprob.j)),(catprob.j[nrow(catprob.j),]),"/"))
  }
  return(score.p)
}

inverse_crossprod <- function(x) {
  if(!is.null(x))  MASS::ginv(crossprod(x))
}

which.max.randomtie <- function(x,na.rm=TRUE){
  loc <- which(x==max(x,na.rm = na.rm))
  if(length(loc)>1){
    loc <- sample(loc,1)
  }
  return(loc)
}


uP <- function(mloc, mpar) {
  .Call('_GDINA_uP', PACKAGE = 'GDINA', mloc, mpar)
}

ObsLogLik <- function(mpar, mX, vlogPrior, vgroup, mloc, weights) {
  .Call('_GDINA_ObsLogLik', PACKAGE = 'GDINA', mpar, mX, vlogPrior, vgroup, mloc, weights)
}

LikNR <- function(mpar, mX, vlogPrior, vgroup, mloc, weights, simplify = TRUE) {
  .Call('_GDINA_LikNR', PACKAGE = 'GDINA', mpar, mX, vlogPrior, vgroup, mloc, weights, simplify)
}

LikNR_LC <- function(mP, mX, vlogPrior, vgroup, weights, simplify = 1L) {
  .Call('_GDINA_LikNR_LC', PACKAGE = 'GDINA', mP, mX, vlogPrior, vgroup, weights, simplify)
}

Lik_DTM <- function(mP, mX, vC, vlogPrior) {
  .Call('_GDINA_Lik_DTM', PACKAGE = 'GDINA', mP, mX, vC, vlogPrior)
}

fast_GDINA_EM <- function(mloc, mpar, mX, vlogPrior, model_numeric, maxitr, lP, uP, smallNcorrection, vbeta, prior = FALSE, crit = 0.0001) {
  .Call('_GDINA_fast_GDINA_EM', PACKAGE = 'GDINA', mloc, mpar, mX, vlogPrior, model_numeric, maxitr, lP, uP, smallNcorrection, vbeta, prior, crit)
}

LouisC <- function(mX, np, mlogPost, itemparmLC, parloc, weight, SEtype) {
  .Call('_GDINA_LouisC', PACKAGE = 'GDINA', mX, np, mlogPost, itemparmLC, parloc, weight, SEtype)
}

Mord <- function(item_no, LCprob, prior) {
  .Call('_GDINA_Mord', PACKAGE = 'GDINA', item_no, LCprob, prior)
}

Calc_Pj <- function(par, designMj, linkfunc, boundary = FALSE, eps = 1e-16) {
  .Call('_GDINA_Calc_Pj', PACKAGE = 'GDINA', par, designMj, linkfunc, boundary, eps)
}

Calc_Dj <- function(par, designMj, linkfunc, boundary = FALSE, eps = 1e-16) {
  .Call('_GDINA_Calc_Dj', PACKAGE = 'GDINA', par, designMj, linkfunc, boundary, eps)
}

Calc_Pj_jac <- function(par, designMj, linkfunc, boundary = FALSE, eps = 1e-16) {
  .Call('_GDINA_Calc_Pj_jac', PACKAGE = 'GDINA', par, designMj, linkfunc, boundary, eps)
}

Mstep_obj_fn <- function(par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix = NULL, eps = 1e-16, ConstrType = 0L, greaterthan0 = TRUE) {
  .Call('_GDINA_Mstep_obj_fn', PACKAGE = 'GDINA', par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix, eps, ConstrType, greaterthan0)
}

Mstep_obj_fn_prior <- function(par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix = NULL, eps = 1e-16, ConstrType = 0L, greaterthan0 = TRUE, m = 0, sd = 5) {
  .Call('_GDINA_Mstep_obj_fn_prior', PACKAGE = 'GDINA', par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix, eps, ConstrType, greaterthan0, m, sd)
}

Mstep_obj_fn_max <- function(par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix = NULL, eps = 1e-16, ConstrType = 0L, greaterthan0 = TRUE) {
  .Call('_GDINA_Mstep_obj_fn_max', PACKAGE = 'GDINA', par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix, eps, ConstrType, greaterthan0)
}

Mstep_obj_gr <- function(par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix = NULL, eps = 1e-16, ConstrType = 0L, greaterthan0 = TRUE) {
  .Call('_GDINA_Mstep_obj_gr', PACKAGE = 'GDINA', par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix, eps, ConstrType, greaterthan0)
}

Mstep_ineq_fn <- function(par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix = NULL, eps = 1e-16, ConstrType = 0L, greaterthan0 = TRUE) {
  .Call('_GDINA_Mstep_ineq_fn', PACKAGE = 'GDINA', par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix, eps, ConstrType, greaterthan0)
}

Mstep_ineq_jac <- function(par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix = NULL, eps = 1e-16, ConstrType = 0L, greaterthan0 = TRUE) {
  .Call('_GDINA_Mstep_ineq_jac', PACKAGE = 'GDINA', par, Nj, Rj, designMj, uPj, lPj, linkfunc, ConstrMatrix, eps, ConstrType, greaterthan0)
}

NgRg <- function(mlogPost, mX, mloc, weights) {
  .Call('_GDINA_NgRg', PACKAGE = 'GDINA', mlogPost, mX, mloc, weights)
}

Rljs_DTM <- function(mlogPost, mX, vC) {
  .Call('_GDINA_Rljs_DTM', PACKAGE = 'GDINA', mlogPost, mX, vC)
}

SE <- function(mX, mlogPost, itmpar, parloc, model, mIndmiss, SE_type) {
  .Call('_GDINA_SE', PACKAGE = 'GDINA', mX, mlogPost, itmpar, parloc, model, mIndmiss, SE_type)
}

aggregateCol <- function(mX, ind) {
  .Call('_GDINA_aggregateCol', PACKAGE = 'GDINA', mX, ind)
}

fitstats <- function(mX, Xfit, cor = TRUE) {
  .Call('_GDINA_fitstats', PACKAGE = 'GDINA', mX, Xfit, cor)
}

scorefun <- function(mX, mlogPost, itmpar, parloc, model) {
  .Call('_GDINA_scorefun', PACKAGE = 'GDINA', mX, mlogPost, itmpar, parloc, model)
}

sequP <- function(mloc, mpar, vC) {
  .Call('_GDINA_sequP', PACKAGE = 'GDINA', mloc, mpar, vC)
}

combnCpp <- function(n, k) {
  .Call('_GDINA_combnCpp', PACKAGE = 'GDINA', n, k)
}

rowProd <- function(m, v) {
  .Call('_GDINA_rowProd', PACKAGE = 'GDINA', m, v)
}

whichrow_AinB <- function(A, B) {
  .Call('_GDINA_whichrow_AinB', PACKAGE = 'GDINA', A, B)
}

whichcol_AinB <- function(A, B) {
  .Call('_GDINA_whichcol_AinB', PACKAGE = 'GDINA', A, B)
}

unique_rows <- function(A) {
  .Call('_GDINA_unique_rows', PACKAGE = 'GDINA', A)
}

alpha2 <- function(K) {
  .Call('_GDINA_alpha2', PACKAGE = 'GDINA', K)
}

alphap <- function(maxlevel) {
  .Call('_GDINA_alphap', PACKAGE = 'GDINA', maxlevel)
}

ColNormalize <- function(X) {
  .Call('_GDINA_ColNormalize', PACKAGE = 'GDINA', X)
}

RowNormalize <- function(X) {
  .Call('_GDINA_RowNormalize', PACKAGE = 'GDINA', X)
}

Pr_2PL <- function(theta, a, b) {
  .Call('_GDINA_Pr_2PL', PACKAGE = 'GDINA', theta, a, b)
}

Pr_2PL_vec <- function(theta, a, b, minvalue = 1e-16, maxvalue = 1 - 1e-16) {
  .Call('_GDINA_Pr_2PL_vec', PACKAGE = 'GDINA', theta, a, b, minvalue, maxvalue)
}

logLikPattern <- function(AlphaPattern, theta, a, b) {
  .Call('_GDINA_logLikPattern', PACKAGE = 'GDINA', AlphaPattern, theta, a, b)
}

PostTheta <- function(AlphaPattern, theta, f_theta, a, b) {
  .Call('_GDINA_PostTheta', PACKAGE = 'GDINA', AlphaPattern, theta, f_theta, a, b)
}

expectedNR <- function(AlphaPattern, nc, theta, f_theta, a, b) {
  .Call('_GDINA_expectedNR', PACKAGE = 'GDINA', AlphaPattern, nc, theta, f_theta, a, b)
}

logP_AlphaPattern <- function(AlphaPattern, theta, f_theta, a, b) {
  .Call('_GDINA_logP_AlphaPattern', PACKAGE = 'GDINA', AlphaPattern, theta, f_theta, a, b)
}

HoIRTlogLik <- function(AlphaPattern, ns, theta, f_theta, a, b) {
  .Call('_GDINA_HoIRTlogLik', PACKAGE = 'GDINA', AlphaPattern, ns, theta, f_theta, a, b)
}

HoIRTlogLik3 <- function(ns, mX, theta, f_theta, a, b) {
  .Call('_GDINA_HoIRTlogLik3', PACKAGE = 'GDINA', ns, mX, theta, f_theta, a, b)
}

incomplogL <- function(a, b, logL, AlphaPattern, theta, f_theta) {
  .Call('_GDINA_incomplogL', PACKAGE = 'GDINA', a, b, logL, AlphaPattern, theta, f_theta)
}

designM <- function(Kj, rule, AlphaPattern = NULL) {
  .Call('_GDINA_designM', PACKAGE = 'GDINA', Kj, rule, AlphaPattern)
}

matchMatrix <- function(A, B) {
  .Call('_GDINA_matchMatrix', PACKAGE = 'GDINA', A, B)
}


item_latent_group <- function(Q, AlphaPattern = NULL) {
  .Call('_GDINA_item_latent_group', PACKAGE = 'GDINA', Q, AlphaPattern)
}

varsigma <- function(mloc, mP, vw) {
  .Call('_GDINA_varsigma', PACKAGE = 'GDINA', mloc, mP, vw)
}
